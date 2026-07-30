//! Five source-closed, non-wire base E2E projections.
//!
//! These projections exercise real public constructors and validators in
//! process. They do not create or execute the downstream-owned shell scripts,
//! emit lifecycle records, publish bundles, or mint authority.

#![allow(
    deprecated,
    reason = "this module constructs and verifies the retained legacy coverage and root compatibility projections"
)]

use crate::budget::{
    AdmittedRunnerBudgetsV2, RunnerBudgetExpectationV2, RunnerBudgetFieldV2, RunnerBudgetUnitV2,
    RunnerBudgetValueV2, RunnerBudgetViolationKindV2, RunnerBudgetViolationV2,
    RunnerBudgetsCandidateV2, RunnerBudgetsV2,
};
use crate::canonical::CanonicalFrameV1;
use crate::capability::{
    NarrowedPolicyViewV2, OverlapPolicyRegistrationV2, RootCapabilityPolicyV2,
    RootPolicyRegistryProjectionV2, expected_rights,
};
use crate::catalog::{
    ArtifactDispositionV2, ArtifactRoleV2, DecisionDetailNamespaceRegistryV2,
    DestinationAdmissionModeV2, DiagnosticCodeV2, DigestRoleV2, LifecycleRecordKindV2,
    LogicalExtentAxisV2, LogicalUnitV2, NotRunCauseCodeV2, OverlapPolicyRelationV2,
    PlatformPathProfileV2, ProofExitV2, PublicationProtocolV2, RefusedReasonV2, RepairActionKindV2,
    RetryabilityV2, RootCapabilityAccessV2, RootCapabilityRightV2, RootClassV2, RunProfileV2,
    RunnerApiGeneration, RunnerCommandV2, RunnerWireVersion, StateBearingRecordRoleV2,
    TypedOptionTagV1, TypedValueTagV2, WirePredecessorPolicyV1,
};
use crate::command::{
    CommandIntentV2, CommandSelectionV2, CommandSelectorCardinalityV2, CommandSelectorFieldV2,
    CommandSelectorPresenceV2, validate_command_selector_presence_v2,
};
use crate::construction::{
    ConstructionClosedSemanticV2, ConstructionErrorKindV2, ConstructionErrorV2,
    ConstructionFixedObservationV2, ConstructionObservedDataClassV2, ConstructionObservedV2,
};
use crate::coverage::{
    BaseCoverageCaseDeclarationV1, BaseCoverageCheckedReportV1,
    BaseCoverageCloseCapabilityProfileV1, BaseCoverageCloseDecisionV1,
    BaseCoverageCloseDeferredEvidenceEnvelopeV1, BaseCoverageCloseDeferredObservationContractV1,
    BaseCoverageCloseDownstreamContributionV2, BaseCoverageCloseExecutionScopeV1,
    BaseCoverageCloseFacetV1, BaseCoverageCloseGroupV1, BaseCoverageCloseManifestCellV1,
    BaseCoverageCloseManifestV1, BaseCoverageClosePartitionV1, BaseCoverageCloseReasonCodeV1,
    BaseCoverageCloseRetainedRelativeArtifactPolicyV1, BaseCoverageManifestClassV1,
    BaseCoverageManifestV1, BaseCoveragePresentedOutcomeV1, BaseCoveragePresentedResultV1,
    CanonicalSchemaImpactDispositionV1, CompatibleSourceSnapshotRootV1, DeferredReasonV1,
    RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1, RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1,
    RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1, RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1,
    RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1, RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
    RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1, RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
    RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1, RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1,
    compatible_source_snapshot_root_from_exact_frame_v1,
    frozen_downstream_close_contribution_budgets_v1,
};
use crate::dependency::{
    CURRENT_DIRECT_DEPENDENCY_DECLARATION_DOMAIN_V1, current_direct_dependency_declaration_root_v1,
};
use crate::diagnostic::{
    ActionableDiagnosticV2, DiagnosticCodeRefV2, DiagnosticEnvelopeGrantsV2, DiagnosticValueV2,
    RegisteredDecisionDetailProjectionV2, RepairActionV2,
};
use crate::extension::{
    ArtifactCodecIdV2, BaseExtensionRegistryProjectionV2, RegisteredLogicalUnitDescriptorV2,
};
use crate::identity::{
    BuildIdentityRootV2, CancelledStopRootV2, CaseManifestRootV2, DrainedInternalErrorRootV2,
    IdentityError, NoClaimScopeRootV1, SourceIdentityRootV2, TimedOutStopRootV2,
    ToolchainIdentityRootV2,
};
use crate::limits::{
    ArtifactStorageProjectionV2, PublicationStorageProjectionV2, RunnerFamilyLimitRequirementsV2,
    RunnerLimitExpectationV2, RunnerLimitFieldV2, RunnerLimitTightenabilityV2, RunnerLimitUnitV2,
    RunnerLimitValueV2, RunnerLimitsV2, RunnerLimitsViolationKindV2, RunnerLimitsViolationV2,
    SystemObjectStorageProjectionV2, SystemPublicationObjectRoleV2,
};
use crate::logging::{
    BaseE2eLogEventV1, BaseE2eLogFieldV1, BaseE2eLogKindV1, BaseE2eLogV1, BaseE2eOutcomeV1,
    SymbolicReproductionArgV1,
};
use crate::path::{
    ContentStoreObjectKeyV1, LogicalBundlePathV1, PathError, PathSetAdjudicationV1,
    adjudicate_logical_bundle_path_set,
};
use crate::publication::{
    PublicationSelectionV2, PublicationTargetV2, SymbolicCommandResultPlanV2,
};
use crate::schema_impact::{
    RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1, RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1,
    RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1, SchemaImpactManifestRelationV1,
    SchemaImpactManifestV1, runner_v2_base_schema_impact_manifest_v1,
};
use crate::state::{
    NotRunBasisErrorV2, NotRunBasisV2, NotRunCauseV2, PresentedDrainRootKindV2,
    StateValidationErrorV2, StateValidationInputV2, validate_state_v2,
};
use crate::value::{OpaqueBytesV2, RationalV2, StableTokenV2, TypedValueV2, ValueError};
use core::fmt::Write as _;
use fs_blake3::{ContentHash, hash_domain};

/// Typed witness of the exact existing compiled source-snapshot root.
///
/// Only this source-closure module can mint the production witness. Consumers
/// receive it from [`RunnerV2BaseSourceClosureV1::compatible_snapshot`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CompatibleSourceSnapshotV1 {
    root: CompatibleSourceSnapshotRootV1,
}

impl CompatibleSourceSnapshotV1 {
    fn from_exact_frame(frame: &CanonicalFrameV1) -> Result<Self, ConstructionErrorV2> {
        Ok(Self {
            root: compatible_source_snapshot_root_from_exact_frame_v1(frame)?,
        })
    }

    /// Exact nominal snapshot root.
    #[must_use]
    pub const fn root(self) -> CompatibleSourceSnapshotRootV1 {
        self.root
    }
}

/// Crate-private proof that one logical source path is an exact member of the
/// compatible compiled source closure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CompatibleSourceMemberV1 {
    path: LogicalBundlePathV1,
    snapshot: CompatibleSourceSnapshotV1,
}

impl CompatibleSourceMemberV1 {
    pub(crate) const fn path(&self) -> &LogicalBundlePathV1 {
        &self.path
    }

    pub(crate) const fn snapshot(&self) -> CompatibleSourceSnapshotV1 {
        self.snapshot
    }
}

/// Overall non-wire projection root domain.
pub const BASE_E2E_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-projection.v1";
/// Per-journey non-wire projection root domain.
pub const BASE_E2E_JOURNEY_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-journey-projection.v1";
/// Domain for one closed base E2E semantic-row descriptor.
pub const BASE_E2E_SEMANTIC_ROW_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-semantic-row.v1";
/// Domain for one journey-specific mapping of a semantic row.
pub const BASE_E2E_JOURNEY_ROW_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-journey-row.v1";
/// Domain for one checked row-result projection.
pub const BASE_E2E_ROW_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-row-result.v1";
/// Domain for one checked row result carrying a private execution witness.
pub const BASE_E2E_EXECUTED_ROW_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-executed-row-result.v1";
/// Domain for one caller-presented row observation before manifest joining.
pub const BASE_E2E_PRESENTED_ROW_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-presented-row-result.v1";
/// Domain for one public, comparison-only journey report.
pub const BASE_E2E_JOURNEY_COMPARISON_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-journey-comparison.v1";
/// Domain for one exact expected or observed refusal/unsupported cell.
pub const BASE_E2E_DECISION_DETAIL_CELL_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-decision-detail-cell.v1";
/// Domain for one ordered per-case decision-detail manifest.
pub const BASE_E2E_DECISION_DETAIL_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-decision-detail-manifest.v1";
/// Maximum canonical bytes retained while checking one caller-presented
/// decision-detail cell.
pub const BASE_E2E_DETAIL_CELL_MAX_ENCODED_BYTES_V1: usize = 4 * 1024;
/// Closed stable ID for a red row contract whose typed detail cells are exact.
pub const BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1: &str = "row.contract";
/// Legacy closed failure ID formerly used when comparison and execution
/// authority shared one report type.
#[deprecated(
    since = "0.1.0",
    note = "comparison reports no longer claim execution; use comparison_root and exact_match"
)]
pub const BASE_E2E_EXECUTION_OBSERVATION_REQUIRED_ID_V1: &str = "execution.observation";
/// Legacy private observation domain retained only as a source-compatible
/// constant.
#[deprecated(
    since = "0.1.0",
    note = "use BASE_E2E_IN_PROCESS_ROW_EXECUTION_WITNESS_DOMAIN_V1"
)]
pub const BASE_E2E_EXECUTION_OBSERVATION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-execution-observation.v1";
/// Domain for one private, in-process, row-specific execution witness.
pub const BASE_E2E_IN_PROCESS_ROW_EXECUTION_WITNESS_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-in-process-row-execution-witness.v1";
/// First semantic ordinal used when command applicability/setup itself fails.
const COMMAND_APPLICABILITY_SETUP_SEMANTIC_ORDINAL_V1: u32 = 1;
/// First semantic ordinal used when the base budget cannot be constructed.
const BUDGET_BASE_CONSTRUCTION_SEMANTIC_ORDINAL_V1: u32 = 1;
/// First semantic ordinal used when the capability registry cannot be constructed.
const CAPABILITY_REGISTRY_SETUP_SEMANTIC_ORDINAL_V1: u32 = 1;
/// Canonical identity domain for one handwritten semantic oracle manifest.
pub const BASE_E2E_ORACLE_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-oracle-manifest.v1";
/// Domain for the opaque, non-authoritative case-conformance detail fixture
/// referenced by the base diagnostic projection.
pub const BASE_E2E_REGISTERED_DETAIL_FIXTURE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-registered-detail-fixture.v1";
/// Domain for one journey execution under exact presented build context.
pub const BASE_E2E_JOURNEY_EXECUTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-journey-execution.v1";
/// Domain for the ordered five-journey aggregate execution under one harness.
pub const BASE_E2E_PROJECTION_EXECUTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-projection-execution.v1";
/// Domain for an exact canonical feature set.
pub const BASE_E2E_FEATURE_SET_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-feature-set.v1";
/// Domain for an exact target token.
pub const BASE_E2E_TARGET_DOMAIN_V1: &str = "org.frankensim.fs-evidence-runner.base-e2e-target.v1";
/// Domain for the complete presented harness context.
pub const BASE_E2E_HARNESS_CONTEXT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-harness-context.v1";
/// Domain for one exact embedded source file's raw bytes.
pub const BASE_SOURCE_FILE_CONTENT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-file-content.v1";
/// Domain for one declarative, non-live-proof source identity.
pub const BASE_EXPECTED_SOURCE_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-expected-source-identity.v1";
/// Domain for one path-, length-, and content-bound source entry.
pub const BASE_SOURCE_FILE_ENTRY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-file-entry.v1";
/// Domain for the common exact compile-time source snapshot.
pub const BASE_SOURCE_SNAPSHOT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-snapshot.v1";
/// Domain for the exact ordered base-schema source closure.
pub const BASE_SOURCE_CLOSURE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-closure.v1";
/// Domain for the immutable, result-free coverage-source inventory.
pub const BASE_COVERAGE_INVENTORY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-inventory.v1";

/// Domain for one V2 row projected from the source-frozen schema-impact DAG.
pub const BASE_SCHEMA_IMPACT_PROJECTION_ROW_DOMAIN_V2: &str =
    "org.frankensim.fs-evidence-runner.base-schema-impact-projection-row.v2";
/// Domain for the exact V2 projection of the base schema-impact manifest.
pub const BASE_SCHEMA_IMPACT_PROJECTION_DOMAIN_V2: &str =
    "org.frankensim.fs-evidence-runner.base-schema-impact-projection.v2";
/// Domain for the immutable Phase-1 contract payload declaration.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_PAYLOAD_DOMAIN_V2: &str =
    "org.frankensim.fs-evidence-runner.phase-one-contract-payload.v2";
/// Domain for the downstream-owned Phase-1 case-manifest declaration.
pub const RUNNER_V2_PHASE_ONE_CASE_MANIFEST_CONTRACT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-case-manifest-contract.v1";
/// Domain for the exact Phase-1 expected-partition declaration.
pub const RUNNER_V2_PHASE_ONE_EXPECTED_PARTITIONS_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-expected-partitions.v1";
/// Domain for the Phase-1 portable target-route declaration.
pub const RUNNER_V2_PHASE_ONE_TARGET_CONTRACT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-target-contract.v1";
/// Domain for the Phase-1 feature-set declaration.
pub const RUNNER_V2_PHASE_ONE_FEATURE_CONTRACT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-feature-contract.v1";
/// Domain for the Phase-1 build-input declaration.
pub const RUNNER_V2_PHASE_ONE_BUILD_CONTRACT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-build-contract.v1";
/// Domain for the Phase-1 toolchain-input declaration.
pub const RUNNER_V2_PHASE_ONE_TOOLCHAIN_CONTRACT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-toolchain-contract.v1";
/// Domain for the complete contribution plus separate Deferred envelope join.
pub const RUNNER_V2_PHASE_ONE_CONTRIBUTION_PROJECTION_DOMAIN_V2: &str =
    "org.frankensim.fs-evidence-runner.phase-one-contribution-projection.v2";

/// Exact base schema-impact projection row count.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2: usize = 3;
/// Exact number of source paths represented by the base schema-impact rows.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SOURCE_PATH_COUNT_V2: usize = 1;
/// Exact number of owned rows in the base schema-impact projection.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_OWNED_ROW_COUNT_V2: usize = 3;
/// Exact number of consumed rows in the base schema-impact projection.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_CONSUMED_ROW_COUNT_V2: usize = 0;

const RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SCHEMA_IDS_V2: [&str;
    RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2] = [
    "canonical-schema-field-descriptor",
    "canonical-schema-frame-descriptor",
    "canonical-schema-version-slot-descriptor",
];

const RUNNER_V2_PHASE_ONE_EXPECTED_PARTITION_NAMES_V1: [&str; 6] = [
    "positive",
    "expected-refusal",
    "expected-failure",
    "mutation",
    "unsupported",
    "inapplicable",
];

const RUNNER_V2_PHASE_ONE_SEMANTIC_JOURNEY_V1: &str = "runner-v2-phase1-contract";
const RUNNER_V2_PHASE_ONE_LEDGER_ROW_NO_CLAIM_V1: &str =
    "phase1-ledger-row-is-expected-contract-not-execution-actual-roots-require-real-observation";
const RUNNER_V2_PHASE_ONE_PARTITION_VOCABULARY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-partition-vocabulary.v1";
const RUNNER_V2_PHASE_ONE_EXPECTED_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.phase-one-expected-result.v1";

/// Sole declaration owner of one compiled base source input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseSourceOwnerV1 {
    /// Shared workspace build and governance input.
    FrankensimWorkspaceGovernance = 1,
    /// Runner V2 base-schema declaration owned by this leaf.
    RunnerV2BaseSchema = 2,
}

impl BaseSourceOwnerV1 {
    /// Stable owner token retained in human-readable evidence.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::FrankensimWorkspaceGovernance => "frankensim-workspace-governance",
            Self::RunnerV2BaseSchema => "runner-v2-base-schema",
        }
    }
}

/// Closed compile-time route by which one base source input is consumed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseSourceRouteV1 {
    /// Workspace Cargo configuration.
    WorkspaceCargoConfig = 1,
    /// Workspace Cargo dependency lock.
    WorkspaceLockfile = 2,
    /// Workspace Cargo manifest.
    WorkspaceManifest = 3,
    /// Workspace constellation lock.
    WorkspaceConstellationLock = 4,
    /// Runner V2 crate contract.
    CrateContract = 5,
    /// Runner V2 crate manifest.
    CrateManifest = 6,
    /// Runner V2 Rust module.
    CrateModule = 7,
    /// Workspace Rust toolchain declaration.
    WorkspaceToolchain = 8,
}

impl BaseSourceRouteV1 {
    /// Stable route token retained in human-readable evidence.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::WorkspaceCargoConfig => "workspace-cargo-config",
            Self::WorkspaceLockfile => "workspace-lockfile",
            Self::WorkspaceManifest => "workspace-manifest",
            Self::WorkspaceConstellationLock => "workspace-constellation-lock",
            Self::CrateContract => "crate-contract",
            Self::CrateManifest => "crate-manifest",
            Self::CrateModule => "crate-module",
            Self::WorkspaceToolchain => "workspace-toolchain",
        }
    }
}

/// Closed snapshot policy for every compiled base source entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseSourceSnapshotPolicyV1 {
    /// The entry must carry the one exact common compiled-snapshot root.
    ExactCommonCompiledSnapshot = 1,
}

impl BaseSourceSnapshotPolicyV1 {
    /// Stable policy token retained in human-readable evidence.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ExactCommonCompiledSnapshot => "exact-common-compiled-snapshot-v1",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct EmbeddedSourceFileV1 {
    path: &'static str,
    bytes: &'static [u8],
    owner: BaseSourceOwnerV1,
    source_route: BaseSourceRouteV1,
    expected_source_identity: &'static str,
    snapshot_policy: BaseSourceSnapshotPolicyV1,
}

// This is the exact bytewise-lexicographic source set owned or consumed by the
// base leaf. `include_bytes!` makes the compiled projection move whenever any
// source, contract, manifest, or lock input changes.
const EMBEDDED_SOURCE_FILES_V1: [EmbeddedSourceFileV1; 27] = [
    EmbeddedSourceFileV1 {
        path: ".cargo/config.toml",
        bytes: include_bytes!("../../../.cargo/config.toml"),
        owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
        source_route: BaseSourceRouteV1::WorkspaceCargoConfig,
        expected_source_identity: "frankensim.workspace.cargo-config.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "Cargo.lock",
        bytes: include_bytes!("../../../Cargo.lock"),
        owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
        source_route: BaseSourceRouteV1::WorkspaceLockfile,
        expected_source_identity: "frankensim.workspace.cargo-lock.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "Cargo.toml",
        bytes: include_bytes!("../../../Cargo.toml"),
        owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
        source_route: BaseSourceRouteV1::WorkspaceManifest,
        expected_source_identity: "frankensim.workspace.manifest.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "constellation.lock",
        bytes: include_bytes!("../../../constellation.lock"),
        owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
        source_route: BaseSourceRouteV1::WorkspaceConstellationLock,
        expected_source_identity: "frankensim.workspace.constellation-lock.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/CONTRACT.md",
        bytes: include_bytes!("../CONTRACT.md"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateContract,
        expected_source_identity: "frankensim.fs-evidence-runner.contract.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/Cargo.toml",
        bytes: include_bytes!("../Cargo.toml"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateManifest,
        expected_source_identity: "frankensim.fs-evidence-runner.manifest.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/budget.rs",
        bytes: include_bytes!("budget.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.budget.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/canonical.rs",
        bytes: include_bytes!("canonical.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.canonical.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/capability.rs",
        bytes: include_bytes!("capability.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.capability.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/catalog.rs",
        bytes: include_bytes!("catalog.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.catalog.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/command.rs",
        bytes: include_bytes!("command.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.command.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/construction.rs",
        bytes: include_bytes!("construction.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.construction.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/coverage.rs",
        bytes: include_bytes!("coverage.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.coverage.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/dependency.rs",
        bytes: include_bytes!("dependency.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.dependency.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/diagnostic.rs",
        bytes: include_bytes!("diagnostic.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.diagnostic.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/extension.rs",
        bytes: include_bytes!("extension.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.extension.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/identity.rs",
        bytes: include_bytes!("identity.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.identity.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/lib.rs",
        bytes: include_bytes!("lib.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.lib.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/limits.rs",
        bytes: include_bytes!("limits.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.limits.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/logging.rs",
        bytes: include_bytes!("logging.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.logging.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/path.rs",
        bytes: include_bytes!("path.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.path.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/projection.rs",
        bytes: include_bytes!("projection.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.projection.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/publication.rs",
        bytes: include_bytes!("publication.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.publication.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/schema_impact.rs",
        bytes: include_bytes!("schema_impact.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.schema-impact.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/state.rs",
        bytes: include_bytes!("state.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.state.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/value.rs",
        bytes: include_bytes!("value.rs"),
        owner: BaseSourceOwnerV1::RunnerV2BaseSchema,
        source_route: BaseSourceRouteV1::CrateModule,
        expected_source_identity: "frankensim.fs-evidence-runner.src.value.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
    EmbeddedSourceFileV1 {
        path: "rust-toolchain.toml",
        bytes: include_bytes!("../../../rust-toolchain.toml"),
        owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
        source_route: BaseSourceRouteV1::WorkspaceToolchain,
        expected_source_identity: "frankensim.workspace.rust-toolchain.v1",
        snapshot_policy: BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot,
    },
];

/// Raw, untrusted input to exact source-closure reconstruction.
///
/// This type intentionally makes no canonicality or proof claim. Only
/// [`RunnerV2BaseSourceClosureV1::reconstruct`] can admit the exact compiled
/// source set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseSourceClosureInputV1 {
    path: String,
    owner_code: u16,
    source_route_code: u16,
    expected_source_identity_root: ContentHash,
    snapshot_policy_code: u16,
    snapshot_root: ContentHash,
    encoded_bytes: u64,
    content_root: ContentHash,
    bytes: Vec<u8>,
}

impl BaseSourceClosureInputV1 {
    /// Presents bytes using the exact registered metadata when the path is
    /// known.
    ///
    /// An unknown path receives deliberately unregistered metadata and will
    /// refuse during exact-set reconstruction. Use
    /// [`Self::presented_with_metadata`] to exercise individual metadata
    /// mutations.
    #[must_use]
    pub fn presented(path: impl Into<String>, bytes: impl Into<Vec<u8>>) -> Self {
        let path = path.into();
        let bytes = bytes.into();
        if let Some(declaration) = source_declaration(&path)
            && let Ok(snapshot_root) = compiled_source_snapshot_root()
            && let Ok(input) = Self::from_declaration(declaration, bytes.clone(), snapshot_root)
        {
            return input;
        }
        let encoded_bytes = u64::try_from(bytes.len()).unwrap_or(u64::MAX);
        let content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, &bytes);
        Self {
            expected_source_identity_root: hash_domain(
                BASE_EXPECTED_SOURCE_IDENTITY_DOMAIN_V1,
                path.as_bytes(),
            ),
            snapshot_root: hash_domain(
                BASE_SOURCE_SNAPSHOT_DOMAIN_V1,
                b"unregistered-source-input",
            ),
            path,
            owner_code: 0,
            source_route_code: 0,
            snapshot_policy_code: 0,
            encoded_bytes,
            content_root,
            bytes,
        }
    }

    /// Presents one exact registered source input.
    pub fn exact_presented(
        path: impl Into<String>,
        bytes: impl Into<Vec<u8>>,
    ) -> Result<Self, ConstructionErrorV2> {
        let path = path.into();
        let declaration = source_declaration(&path).ok_or_else(|| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Unexpected,
                "base_source_closure.path",
                "a member of the exact embedded source set",
                ConstructionObservedDataClassV2::CallerControlledText,
            )
        })?;
        let snapshot_root = compiled_source_snapshot_root()?;
        Self::from_declaration(declaration, bytes.into(), snapshot_root)
    }

    /// Constructs a completely caller-presented metadata input.
    ///
    /// No field is trusted by this constructor. Exact reconstruction compares
    /// all fields with the independent compile-time declaration table.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn presented_with_metadata(
        path: impl Into<String>,
        bytes: impl Into<Vec<u8>>,
        owner_code: u16,
        source_route_code: u16,
        expected_source_identity_root: ContentHash,
        snapshot_policy_code: u16,
        snapshot_root: ContentHash,
        encoded_bytes: u64,
        content_root: ContentHash,
    ) -> Self {
        Self {
            path: path.into(),
            owner_code,
            source_route_code,
            expected_source_identity_root,
            snapshot_policy_code,
            snapshot_root,
            encoded_bytes,
            content_root,
            bytes: bytes.into(),
        }
    }

    fn from_declaration(
        declaration: &EmbeddedSourceFileV1,
        bytes: Vec<u8>,
        snapshot_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        let encoded_bytes = u64::try_from(bytes.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_source_closure.encoded_bytes",
                "a u64 source byte length",
                bytes.len(),
            )
        })?;
        let content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, &bytes);
        Ok(Self {
            path: declaration.path.to_owned(),
            owner_code: declaration.owner as u16,
            source_route_code: declaration.source_route as u16,
            expected_source_identity_root: expected_source_identity_root(declaration),
            snapshot_policy_code: declaration.snapshot_policy as u16,
            snapshot_root,
            encoded_bytes,
            content_root,
            bytes,
        })
    }

    /// Presented relative source path.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Presented source bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Presented owner code.
    #[must_use]
    pub const fn owner_code(&self) -> u16 {
        self.owner_code
    }

    /// Presented source-route code.
    #[must_use]
    pub const fn source_route_code(&self) -> u16 {
        self.source_route_code
    }

    /// Presented declarative expected-source-identity root.
    #[must_use]
    pub const fn expected_source_identity_root(&self) -> ContentHash {
        self.expected_source_identity_root
    }

    /// Presented snapshot-policy code.
    #[must_use]
    pub const fn snapshot_policy_code(&self) -> u16 {
        self.snapshot_policy_code
    }

    /// Presented common compiled-snapshot root.
    #[must_use]
    pub const fn snapshot_root(&self) -> ContentHash {
        self.snapshot_root
    }

    /// Presented encoded source length.
    #[must_use]
    pub const fn encoded_bytes(&self) -> u64 {
        self.encoded_bytes
    }

    /// Presented content root.
    #[must_use]
    pub const fn content_root(&self) -> ContentHash {
        self.content_root
    }
}

/// One exact metadata-, path-, length-, content-, and snapshot-bound entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseSourceClosureEntryV1 {
    path: &'static str,
    owner: BaseSourceOwnerV1,
    source_route: BaseSourceRouteV1,
    expected_source_identity: &'static str,
    expected_source_identity_root: ContentHash,
    snapshot_policy: BaseSourceSnapshotPolicyV1,
    encoded_bytes: u64,
    content_root: ContentHash,
    snapshot_root: ContentHash,
    entry_root: ContentHash,
}

impl BaseSourceClosureEntryV1 {
    /// Exact workspace-relative source path.
    #[must_use]
    pub const fn path(self) -> &'static str {
        self.path
    }

    /// Sole declaration owner for this compile-time input.
    #[must_use]
    pub const fn owner(self) -> BaseSourceOwnerV1 {
        self.owner
    }

    /// Closed source route used to obtain the input.
    #[must_use]
    pub const fn source_route(self) -> BaseSourceRouteV1 {
        self.source_route
    }

    /// Unique declarative source-identity token.
    ///
    /// This token identifies the expected compiled input and does not attest
    /// any ambient or live-tree source.
    #[must_use]
    pub const fn expected_source_identity(self) -> &'static str {
        self.expected_source_identity
    }

    /// Domain-separated declarative expected-source-identity root.
    #[must_use]
    pub const fn expected_source_identity_root(self) -> ContentHash {
        self.expected_source_identity_root
    }

    /// Exact common-snapshot admission policy.
    #[must_use]
    pub const fn snapshot_policy(self) -> BaseSourceSnapshotPolicyV1 {
        self.snapshot_policy
    }

    /// Exact included source byte length.
    #[must_use]
    pub const fn encoded_bytes(self) -> u64 {
        self.encoded_bytes
    }

    /// Domain-separated root of the exact raw file bytes.
    #[must_use]
    pub const fn content_root(self) -> ContentHash {
        self.content_root
    }

    /// Common exact compile-time snapshot bound into every entry.
    #[must_use]
    pub const fn snapshot_root(self) -> ContentHash {
        self.snapshot_root
    }

    /// Domain-separated root binding all exact entry metadata.
    #[must_use]
    pub const fn entry_root(self) -> ContentHash {
        self.entry_root
    }
}

/// Exact immutable source closure for the Runner V2 base-schema leaf.
///
/// ```
/// use fs_evidence_runner::projection::RunnerV2BaseSourceClosureV1;
///
/// let closure = RunnerV2BaseSourceClosureV1::frozen()?;
/// assert!(!closure.entries().is_empty());
/// # Ok::<(), fs_evidence_runner::ConstructionErrorV2>(())
/// ```
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::projection::RunnerV2BaseSourceClosureV1;
///
/// let mut closure = RunnerV2BaseSourceClosureV1::frozen().unwrap();
/// closure.entries.clear();
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseSourceClosureV1 {
    entries: Box<[BaseSourceClosureEntryV1]>,
    compatible_snapshot: CompatibleSourceSnapshotV1,
    dependency_declaration_root: ContentHash,
    root: ContentHash,
}

impl RunnerV2BaseSourceClosureV1 {
    /// Reconstructs the closure from the exact compile-time embedded inputs.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        static FROZEN_SOURCE_CLOSURE_V1: std::sync::OnceLock<
            Result<RunnerV2BaseSourceClosureV1, ConstructionErrorV2>,
        > = std::sync::OnceLock::new();
        FROZEN_SOURCE_CLOSURE_V1
            .get_or_init(Self::build_frozen)
            .clone()
    }

    fn build_frozen() -> Result<Self, ConstructionErrorV2> {
        let inputs = EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| BaseSourceClosureInputV1::exact_presented(file.path, file.bytes.to_vec()))
            .collect::<Result<Vec<_>, _>>()?;
        Self::reconstruct(&inputs)
    }

    /// Checks and reconstructs the one exact ordered source closure.
    ///
    /// Duplicate, missing, extra, reordered, or byte-mutated inputs refuse
    /// before a closure root is returned.
    pub fn reconstruct(inputs: &[BaseSourceClosureInputV1]) -> Result<Self, ConstructionErrorV2> {
        Self::reconstruct_with_dependency_declaration(
            inputs,
            current_direct_dependency_declaration_root_v1(),
        )
    }

    /// Reconstruct with caller-presented declaration-time dependency identity.
    ///
    /// The dependency root is static policy data, not live supply-chain
    /// evidence. A stale or substituted root refuses before source admission.
    pub fn reconstruct_with_dependency_declaration(
        inputs: &[BaseSourceClosureInputV1],
        dependency_declaration_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        let expected_dependency_root = current_direct_dependency_declaration_root_v1();
        if dependency_declaration_root != expected_dependency_root {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_source_closure.dependency_declaration_root",
                "the exact current declaration-time direct-dependency root",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        validate_source_input_set(inputs)?;
        let expected_paths = EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| file.path)
            .collect::<std::collections::BTreeSet<_>>();
        let compatible_snapshot = compiled_source_snapshot()?;
        let snapshot_root = compatible_snapshot.root().content_hash();
        for (ordinal, (input, expected)) in inputs
            .iter()
            .zip(EMBEDDED_SOURCE_FILES_V1.iter())
            .enumerate()
        {
            validate_source_input(input, expected, ordinal, &expected_paths, snapshot_root)?;
        }
        let entries = EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|expected| source_closure_entry(expected, snapshot_root))
            .collect::<Result<Vec<_>, _>>()?;
        let root = source_closure_root(&entries, snapshot_root)?;
        Ok(Self {
            entries: entries.into_boxed_slice(),
            compatible_snapshot,
            dependency_declaration_root,
            root,
        })
    }

    /// Entries in exact canonical source order.
    #[must_use]
    pub fn entries(&self) -> &[BaseSourceClosureEntryV1] {
        &self.entries
    }

    /// Domain-separated root of the complete exact source closure.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Common content-derived compile-time snapshot identity.
    #[must_use]
    pub const fn snapshot_root(&self) -> ContentHash {
        self.compatible_snapshot.root().content_hash()
    }

    /// Nominally typed witness of the exact existing snapshot frame.
    ///
    /// This is the same content hash returned by [`Self::snapshot_root`].
    /// No wrapper frame or digest-of-digest transformation is applied.
    #[must_use]
    pub const fn compatible_snapshot(&self) -> CompatibleSourceSnapshotV1 {
        self.compatible_snapshot
    }

    /// Resolve one exact embedded source member for crate-owned declaration
    /// admission.
    pub(crate) fn compatible_source_member(
        &self,
        path: &'static str,
    ) -> Result<CompatibleSourceMemberV1, ConstructionErrorV2> {
        if !self.entries.iter().any(|entry| entry.path() == path) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Missing,
                "source_closure.compatible_member",
                "one exact path in the compiled source closure",
                ConstructionObservedDataClassV2::PhysicalLocator,
            ));
        }
        let path = LogicalBundlePathV1::new(path).map_err(|_| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "source_closure.compatible_member",
                "one valid logical path in the compiled source closure",
                ConstructionObservedDataClassV2::PhysicalLocator,
            )
        })?;
        Ok(CompatibleSourceMemberV1 {
            path,
            snapshot: self.compatible_snapshot,
        })
    }

    /// Exact declaration-time direct-dependency policy root.
    ///
    /// This is static policy data, not live Cargo, lockfile, constellation, or
    /// checkout proof.
    #[must_use]
    pub const fn dependency_declaration_root(&self) -> ContentHash {
        self.dependency_declaration_root
    }
}

/// One exact row in the additive V2 projection of the schema-impact manifest.
///
/// The row is a projection identity only. It does not parse schema bytes,
/// perform a migration, validate an execution, or mint schema authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseSchemaImpactProjectionRowV2 {
    local_ordinal: u32,
    relation: SchemaImpactManifestRelationV1,
    schema_id: StableTokenV2,
    source_path: LogicalBundlePathV1,
    disposition: CanonicalSchemaImpactDispositionV1,
    schema_impact_row_root: ContentHash,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    root: ContentHash,
}

impl RunnerV2BaseSchemaImpactProjectionRowV2 {
    /// Derived one-based manifest-local ordinal.
    #[must_use]
    pub const fn local_ordinal(&self) -> u32 {
        self.local_ordinal
    }

    /// Exact Owned/Consumed relationship to the base leaf.
    #[must_use]
    pub const fn relation(&self) -> SchemaImpactManifestRelationV1 {
        self.relation
    }

    /// Exact stable schema ID.
    #[must_use]
    pub const fn schema_id(&self) -> &StableTokenV2 {
        &self.schema_id
    }

    /// Exact source-closed declaration path.
    #[must_use]
    pub const fn source_path(&self) -> &LogicalBundlePathV1 {
        &self.source_path
    }

    /// Exact source-frozen schema-impact disposition.
    #[must_use]
    pub const fn disposition(&self) -> CanonicalSchemaImpactDispositionV1 {
        self.disposition
    }

    /// Exact admitted schema-impact row root.
    #[must_use]
    pub const fn schema_impact_row_root(&self) -> ContentHash {
        self.schema_impact_row_root
    }

    /// Exact compatible source-snapshot root.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Domain-separated V2 projection-row root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Exact V2 projection wrapper over the base source-frozen schema-impact DAG.
///
/// This is intentionally a small declaration-only bridge. It does not imply
/// that the larger V2 semantic-row, journey, result, execution, report, or log
/// schema family has been implemented.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseSchemaImpactProjectionV2 {
    manifest: SchemaImpactManifestV1,
    rows: Box<[RunnerV2BaseSchemaImpactProjectionRowV2]>,
    source_path_count: u32,
    owned_row_count: u32,
    consumed_row_count: u32,
    root: ContentHash,
}

impl RunnerV2BaseSchemaImpactProjectionV2 {
    /// Reconstruct the exact source-frozen base projection.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        static FROZEN_SCHEMA_IMPACT_PROJECTION_V2: std::sync::OnceLock<
            Result<RunnerV2BaseSchemaImpactProjectionV2, ConstructionErrorV2>,
        > = std::sync::OnceLock::new();
        FROZEN_SCHEMA_IMPACT_PROJECTION_V2
            .get_or_init(|| {
                let manifest = runner_v2_base_schema_impact_manifest_v1()?;
                Self::build_exact(manifest)
            })
            .clone()
    }

    /// Reconstruct only when the presented manifest is byte-for-byte the
    /// source-owned base manifest.
    pub fn reconstruct_exact(
        manifest: &SchemaImpactManifestV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let expected = runner_v2_base_schema_impact_manifest_v1()?;
        if manifest != &expected {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.schema_impact.manifest",
                "the exact source-frozen base schema-impact manifest",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        Self::build_exact(manifest.clone())
    }

    fn build_exact(manifest: SchemaImpactManifestV1) -> Result<Self, ConstructionErrorV2> {
        validate_schema_impact_projection_manifest_v2(&manifest)?;
        let rows = manifest
            .entries()
            .iter()
            .map(|entry| {
                let schema_id = StableTokenV2::new(entry.row().schema_id().as_str().to_owned())
                    .map_err(|_| {
                        ConstructionErrorV2::new_redacted(
                            ConstructionErrorKindV2::Incompatible,
                            "projection.schema_impact.schema_id",
                            "one exact bounded stable schema ID",
                            ConstructionObservedDataClassV2::CallerControlledText,
                        )
                    })?;
                let source_path = entry.row().source_path().clone();
                let schema_impact_row_root = entry.row().root().content_hash();
                let compatible_source_snapshot_root = entry.row().compatible_source_snapshot_root();
                let root = schema_impact_projection_row_root_v2(
                    entry.local_ordinal(),
                    entry.relation(),
                    &schema_id,
                    &source_path,
                    entry.row().disposition(),
                    schema_impact_row_root,
                    compatible_source_snapshot_root,
                )?;
                Ok(RunnerV2BaseSchemaImpactProjectionRowV2 {
                    local_ordinal: entry.local_ordinal(),
                    relation: entry.relation(),
                    schema_id,
                    source_path,
                    disposition: entry.row().disposition(),
                    schema_impact_row_root,
                    compatible_source_snapshot_root,
                    root,
                })
            })
            .collect::<Result<Vec<_>, ConstructionErrorV2>>()?;
        let source_path_count = checked_u32(
            rows.iter()
                .map(|row| row.source_path.as_str())
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            "projection.schema_impact.source_path_count",
        )?;
        let owned_row_count = checked_u32(
            rows.iter()
                .filter(|row| row.relation == SchemaImpactManifestRelationV1::Owned)
                .count(),
            "projection.schema_impact.owned_row_count",
        )?;
        let consumed_row_count = checked_u32(
            rows.iter()
                .filter(|row| row.relation == SchemaImpactManifestRelationV1::Consumed)
                .count(),
            "projection.schema_impact.consumed_row_count",
        )?;
        let root = schema_impact_projection_root_v2(
            &manifest,
            &rows,
            source_path_count,
            owned_row_count,
            consumed_row_count,
        )?;
        Ok(Self {
            manifest,
            rows: rows.into_boxed_slice(),
            source_path_count,
            owned_row_count,
            consumed_row_count,
            root,
        })
    }

    /// Exact source-frozen schema-impact manifest.
    #[must_use]
    pub const fn manifest(&self) -> &SchemaImpactManifestV1 {
        &self.manifest
    }

    /// Complete rows in exact manifest-local topological order.
    #[must_use]
    pub fn rows(&self) -> &[RunnerV2BaseSchemaImpactProjectionRowV2] {
        &self.rows
    }

    /// Number of distinct exact source paths across the projected rows.
    #[must_use]
    pub const fn source_path_count(&self) -> u32 {
        self.source_path_count
    }

    /// Exact number of base-leaf-owned rows.
    #[must_use]
    pub const fn owned_row_count(&self) -> u32 {
        self.owned_row_count
    }

    /// Exact number of externally owned rows consumed by this base wrapper.
    #[must_use]
    pub const fn consumed_row_count(&self) -> u32 {
        self.consumed_row_count
    }

    /// Domain-separated V2 schema-impact projection root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// This declaration-only projection contains no runtime observation.
    #[must_use]
    pub const fn execution_observation_count(&self) -> u32 {
        0
    }
}

/// Immutable Phase-1 contribution plus its distinct Deferred envelope.
///
/// This wrapper provides upstream declaration evidence only. It does not read
/// or execute either downstream wrapper, does not read the downstream-owned
/// case-manifest file, and cannot report that the designated route passed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2PhaseOneContractContributionV2 {
    schema_impact_projection: RunnerV2BaseSchemaImpactProjectionV2,
    source_closure_root: ContentHash,
    base_coverage_manifest_root: ContentHash,
    case_manifest_contract_root: ContentHash,
    expected_ledger: RunnerV2PhaseOneExpectedLedgerV1,
    expected_partitions_root: ContentHash,
    payload_root: ContentHash,
    contribution: BaseCoverageCloseDownstreamContributionV2,
    deferred_envelope: BaseCoverageCloseDeferredEvidenceEnvelopeV1,
    root: ContentHash,
}

impl RunnerV2PhaseOneContractContributionV2 {
    /// Reconstruct the one exact source-owned Phase-1 contribution.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        static FROZEN_PHASE_ONE_CONTRIBUTION_V2: std::sync::OnceLock<
            Result<RunnerV2PhaseOneContractContributionV2, ConstructionErrorV2>,
        > = std::sync::OnceLock::new();
        FROZEN_PHASE_ONE_CONTRIBUTION_V2
            .get_or_init(Self::build_frozen)
            .clone()
    }

    fn build_frozen() -> Result<Self, ConstructionErrorV2> {
        let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
        let base_coverage_manifest = BaseCoverageManifestV1::frozen()?;
        let base_coverage_close_manifest = BaseCoverageCloseManifestV1::frozen()?;
        let schema_impact_projection = RunnerV2BaseSchemaImpactProjectionV2::frozen()?;
        if schema_impact_projection
            .manifest()
            .compatible_source_snapshot_root()
            != source_closure.compatible_snapshot().root()
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.phase_one.compatible_source_snapshot",
                "one exact compatible source snapshot across contribution inputs",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let expected_ledger = RunnerV2PhaseOneExpectedLedgerV1::build_from(
            &base_coverage_close_manifest,
            &source_closure,
        )?;
        let expected_partitions_root = expected_ledger.root();
        let case_manifest_contract_root = phase_one_case_manifest_contract_root_v1(
            &base_coverage_close_manifest,
            expected_partitions_root,
        )?;
        let target_root = phase_one_target_contract_root_v1()?;
        let feature_set_root = phase_one_feature_contract_root_v1()?;
        let build_contract_root =
            phase_one_build_contract_root_v1(&source_closure, &schema_impact_projection)?;
        let toolchain_contract_root = phase_one_toolchain_contract_root_v1(&source_closure)?;
        let source_root = parse_phase_one_source_identity(source_closure.root())?;
        let build_root = parse_phase_one_build_identity(build_contract_root)?;
        let toolchain_root = parse_phase_one_toolchain_identity(toolchain_contract_root)?;
        let schema_impact_log_schema_root = crate::logging::schema_impact_log_schema_root_v1()?;
        let payload_root = phase_one_payload_root_v2(
            &source_closure,
            base_coverage_manifest.root(),
            &schema_impact_projection,
            expected_partitions_root,
            target_root,
            feature_set_root,
            schema_impact_log_schema_root,
        )?;
        let observer_contract = BaseCoverageCloseDeferredObservationContractV1::new(
            RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
            BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
            token(RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1)?,
            RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
            RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
            case_manifest_contract_root,
            DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
            token(RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1)?,
        )?;
        validate_exact_phase_one_observer_contract_v1(&observer_contract)?;
        let contribution = BaseCoverageCloseDownstreamContributionV2::new(
            RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1,
            payload_root,
            expected_partitions_root,
            frozen_downstream_close_contribution_budgets_v1(1, 1)?,
            frozen_downstream_close_contribution_budgets_v1(4, 2)?,
            schema_impact_projection.root(),
            schema_impact_log_schema_root,
            source_root,
            build_root,
            toolchain_root,
            target_root,
            feature_set_root,
            BaseCoverageCloseRetainedRelativeArtifactPolicyV1::OwnerEnvelopeRelativePathsOnly,
            observer_contract,
            token(RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1)?,
        )?;
        let deferred_envelope = BaseCoverageCloseDeferredEvidenceEnvelopeV1::new(
            &contribution,
            token(RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1)?,
        )?;
        let root = phase_one_contribution_projection_root_v2(
            &schema_impact_projection,
            source_closure.root(),
            base_coverage_manifest.root(),
            case_manifest_contract_root,
            expected_partitions_root,
            payload_root,
            &contribution,
            &deferred_envelope,
        )?;
        Ok(Self {
            schema_impact_projection,
            source_closure_root: source_closure.root(),
            base_coverage_manifest_root: base_coverage_manifest.root(),
            case_manifest_contract_root,
            expected_ledger,
            expected_partitions_root,
            payload_root,
            contribution,
            deferred_envelope,
            root,
        })
    }

    /// Exact declaration-only schema-impact projection.
    #[must_use]
    pub const fn schema_impact_projection(&self) -> &RunnerV2BaseSchemaImpactProjectionV2 {
        &self.schema_impact_projection
    }

    /// Exact compiled source-closure root.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Exact base source-coverage manifest root.
    #[must_use]
    pub const fn base_coverage_manifest_root(&self) -> ContentHash {
        self.base_coverage_manifest_root
    }

    /// Declaration root for the downstream-owned case manifest.
    #[must_use]
    pub const fn case_manifest_contract_root(&self) -> ContentHash {
        self.case_manifest_contract_root
    }

    /// Complete source-authoritative result-free expected ledger.
    #[must_use]
    pub const fn expected_ledger(&self) -> &RunnerV2PhaseOneExpectedLedgerV1 {
        &self.expected_ledger
    }

    /// Exact declaration-side expected partition root.
    #[must_use]
    pub const fn expected_partitions_root(&self) -> ContentHash {
        self.expected_partitions_root
    }

    /// Exact immutable Phase-1 payload root.
    #[must_use]
    pub const fn payload_root(&self) -> ContentHash {
        self.payload_root
    }

    /// Exact result-free V2 contribution.
    #[must_use]
    pub const fn contribution(&self) -> &BaseCoverageCloseDownstreamContributionV2 {
        &self.contribution
    }

    /// Separate exact Deferred evidence envelope.
    #[must_use]
    pub const fn deferred_envelope(&self) -> &BaseCoverageCloseDeferredEvidenceEnvelopeV1 {
        &self.deferred_envelope
    }

    /// Domain-separated contribution/envelope join root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// This upstream declaration performs no downstream execution.
    #[must_use]
    pub const fn downstream_execution_count(&self) -> u32 {
        0
    }
}

fn validate_schema_impact_projection_manifest_v2(
    manifest: &SchemaImpactManifestV1,
) -> Result<(), ConstructionErrorV2> {
    if manifest.issuer_leaf_id().as_str() != RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1 {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.schema_impact.issuer",
            "the exact base schema-impact owner leaf",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    if manifest.entries().len() != RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2 {
        return Err(ConstructionErrorV2::new(
            if manifest.entries().len() < RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2 {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "projection.schema_impact.row_count",
            "the exact three-row base schema-impact manifest",
            manifest.entries().len(),
        ));
    }
    if manifest.graph_edge_count() != 0
        || !manifest.leaf_extension_registries().is_empty()
        || manifest.base_partition_nominal_role_count() != 44
        || manifest.frozen_base_nominal_role_count() != 47
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.schema_impact.manifest_shape",
            "zero edges, zero extensions, 44 base roles, and 47 frozen roles",
            ConstructionObservedV2::unsigned_pair(
                u64::from(manifest.graph_edge_count()),
                manifest.leaf_extension_registries().len() as u64,
            ),
        ));
    }
    for (index, (entry, expected_schema_id)) in manifest
        .entries()
        .iter()
        .zip(RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SCHEMA_IDS_V2)
        .enumerate()
    {
        let expected_ordinal = checked_u32(index + 1, "projection.schema_impact.local_ordinal")?;
        if entry.local_ordinal() != expected_ordinal
            || entry.relation() != SchemaImpactManifestRelationV1::Owned
            || entry.row().schema_id().as_str() != expected_schema_id
            || entry.row().source_path().as_str() != RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1
            || entry.row().disposition() != CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor
            || entry.row().no_claim().as_str() != RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1
            || entry.row().compatible_source_snapshot_root()
                != manifest.compatible_source_snapshot_root()
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "projection.schema_impact.row",
                "the exact ordered owned row ID, path, disposition, snapshot, and no-claim",
                ConstructionObservedV2::indexed_redacted(
                    index,
                    ConstructionObservedDataClassV2::CallerControlledText,
                ),
            ));
        }
    }
    Ok(())
}

#[allow(
    clippy::too_many_arguments,
    reason = "the row projection root visibly binds every exact schema-impact coordinate"
)]
fn schema_impact_projection_row_root_v2(
    local_ordinal: u32,
    relation: SchemaImpactManifestRelationV1,
    schema_id: &StableTokenV2,
    source_path: &LogicalBundlePathV1,
    disposition: CanonicalSchemaImpactDispositionV1,
    schema_impact_row_root: ContentHash,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    if local_ordinal == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "projection.schema_impact.local_ordinal",
            "a nonzero one-based manifest-local ordinal",
            local_ordinal,
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSBASESCHEMAIMPACTPROJECTIONROW\x02", 4 * 1024)?;
    frame.push_u16("projection.api_generation", 2)?;
    frame.push_u16("projection.wire_version", 1)?;
    frame.push_u16("projection.wire_predecessor_policy", 1)?;
    frame.push_u32("projection.local_ordinal", local_ordinal)?;
    frame.push_u16("projection.relation", relation.code())?;
    frame.push_str("projection.schema_id", schema_id.as_str())?;
    frame.push_str("projection.source_path", source_path.as_str())?;
    frame.push_u16("projection.disposition", disposition.code())?;
    frame.push_bytes(
        "projection.schema_impact_row_root",
        schema_impact_row_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.compatible_source_snapshot_root",
        compatible_source_snapshot_root.content_hash().as_bytes(),
    )?;
    frame.push_str(
        "projection.no_claim",
        "schema-impact-projection-row-proves-declaration-binding-not-migration-execution-or-authority",
    )?;
    Ok(frame.root(BASE_SCHEMA_IMPACT_PROJECTION_ROW_DOMAIN_V2))
}

fn schema_impact_projection_root_v2(
    manifest: &SchemaImpactManifestV1,
    rows: &[RunnerV2BaseSchemaImpactProjectionRowV2],
    source_path_count: u32,
    owned_row_count: u32,
    consumed_row_count: u32,
) -> Result<ContentHash, ConstructionErrorV2> {
    let exact_row_count = checked_u32(
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2,
        "projection.schema_impact.row_count",
    )?;
    let observed_row_count = checked_u32(rows.len(), "projection.schema_impact.row_count")?;
    if observed_row_count != exact_row_count {
        return Err(ConstructionErrorV2::new(
            if observed_row_count < exact_row_count {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "projection.schema_impact.row_count",
            "the exact three projected rows",
            observed_row_count,
        ));
    }
    if source_path_count != RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SOURCE_PATH_COUNT_V2 as u32
        || owned_row_count != RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_OWNED_ROW_COUNT_V2 as u32
        || consumed_row_count
            != RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_CONSUMED_ROW_COUNT_V2 as u32
        || owned_row_count
            .checked_add(consumed_row_count)
            .ok_or_else(sequence_overflow)?
            != observed_row_count
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.schema_impact.counts",
            "one source path, three owned rows, zero consumed rows, and an exact partition sum",
            ConstructionObservedV2::unsigned_pair(
                u64::from(owned_row_count),
                u64::from(consumed_row_count),
            ),
        ));
    }
    for (index, (row, manifest_entry)) in rows.iter().zip(manifest.entries()).enumerate() {
        let expected_ordinal = checked_u32(index + 1, "projection.schema_impact.local_ordinal")?;
        let expected_row_root = schema_impact_projection_row_root_v2(
            row.local_ordinal,
            row.relation,
            &row.schema_id,
            &row.source_path,
            row.disposition,
            row.schema_impact_row_root,
            row.compatible_source_snapshot_root,
        )?;
        if row.local_ordinal != expected_ordinal
            || row.relation != SchemaImpactManifestRelationV1::Owned
            || row.schema_id.as_str()
                != RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SCHEMA_IDS_V2[index]
            || row.source_path.as_str() != RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1
            || row.disposition != CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor
            || row.schema_impact_row_root != manifest_entry.row().root().content_hash()
            || row.compatible_source_snapshot_root != manifest.compatible_source_snapshot_root()
            || row.root != expected_row_root
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "projection.schema_impact.rows",
                "the exact three-row source-frozen order and common snapshot",
                index,
            ));
        }
    }
    let mut frame = CanonicalFrameV1::new(b"FSBASESCHEMAIMPACTPROJECTION\x02", 16 * 1024)?;
    frame.push_u16("projection.api_generation", 2)?;
    frame.push_u16("projection.wire_version", 1)?;
    frame.push_u16("projection.wire_predecessor_policy", 1)?;
    frame.push_str(
        "projection.issuer_leaf_id",
        manifest.issuer_leaf_id().as_str(),
    )?;
    frame.push_bytes(
        "projection.compatible_source_snapshot_root",
        manifest
            .compatible_source_snapshot_root()
            .content_hash()
            .as_bytes(),
    )?;
    frame.push_bytes(
        "projection.schema_impact_manifest_root",
        manifest.root().content_hash().as_bytes(),
    )?;
    frame.push_u32(
        "projection.base_partition_nominal_role_count",
        manifest.base_partition_nominal_role_count(),
    )?;
    frame.push_u32(
        "projection.frozen_base_nominal_role_count",
        manifest.frozen_base_nominal_role_count(),
    )?;
    frame.push_u32(
        "projection.extension_fragment_count",
        checked_u32(
            manifest.leaf_extension_registries().len(),
            "projection.schema_impact.extension_fragment_count",
        )?,
    )?;
    frame.push_u32("projection.graph_edge_count", manifest.graph_edge_count())?;
    frame.push_u32("projection.source_path_count", source_path_count)?;
    frame.push_u32("projection.owned_row_count", owned_row_count)?;
    frame.push_u32("projection.consumed_row_count", consumed_row_count)?;
    frame.push_u32("projection.row_count", observed_row_count)?;
    for row in rows {
        frame.push_u32("projection.row_local_ordinal", row.local_ordinal)?;
        frame.push_u16("projection.row_relation", row.relation.code())?;
        frame.push_str("projection.row_schema_id", row.schema_id.as_str())?;
        frame.push_str("projection.row_source_path", row.source_path.as_str())?;
        frame.push_bytes(
            "projection.row_schema_impact_root",
            row.schema_impact_row_root.as_bytes(),
        )?;
        frame.push_bytes("projection.row_projection_root", row.root.as_bytes())?;
    }
    frame.push_str(
        "projection.no_claim",
        "schema-impact-v2-projection-proves-exact-declaration-rows-not-parser-migration-execution-or-authority",
    )?;
    Ok(frame.root(BASE_SCHEMA_IMPACT_PROJECTION_DOMAIN_V2))
}

fn validate_exact_phase_one_observer_contract_v1(
    observer: &BaseCoverageCloseDeferredObservationContractV1,
) -> Result<(), ConstructionErrorV2> {
    let exact_owner = RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1;
    if observer.route_id().as_str() != RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1
        || observer.semantic_consumer() != exact_owner
        || observer.execution_owner() != exact_owner
        || observer.driver_owner() != exact_owner
        || observer.posix_wrapper_owner() != exact_owner
        || observer.windows_wrapper_owner() != exact_owner
        || observer.capability_profile() != BaseCoverageCloseCapabilityProfileV1::ReleaseControl
        || observer.driver().as_str() != RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1
        || observer.posix_route() != RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1
        || observer.windows_route() != RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1
        || observer.case_manifest_path() != RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1
        || observer.deferred_reason()
            != DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        || observer.no_claim().as_str() != RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.phase_one.observer_contract",
            "the exact owner, capability, driver, POSIX, Windows, manifest, and Deferred route",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    Ok(())
}

/// One immutable, result-free row in the source-authoritative Phase-1 ledger.
///
/// Every field is expected contract data. No member records a matched result,
/// runtime observation, retained execution artifact, or downstream success.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2PhaseOneExpectedLedgerCellV1 {
    ledger_ordinal: u32,
    source_ordinal: u32,
    cell_id: StableTokenV2,
    source_case_id: Box<str>,
    source_class: BaseCoverageManifestClassV1,
    source_path: Box<str>,
    group: BaseCoverageCloseGroupV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_result_root: Option<ContentHash>,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    semantic_journey: StableTokenV2,
    route_id: StableTokenV2,
    route_owner: Box<str>,
    compatible_source_snapshot_root: ContentHash,
    case_manifest_root: ContentHash,
    source_closure_root: ContentHash,
    downstream_contribution_root: Option<ContentHash>,
    five_explicits_root: ContentHash,
    cell_root: ContentHash,
    no_claim: StableTokenV2,
}

impl RunnerV2PhaseOneExpectedLedgerCellV1 {
    fn from_close_cell(
        ledger_ordinal: u32,
        cell: &BaseCoverageCloseManifestCellV1,
        close_manifest: &BaseCoverageCloseManifestV1,
        source_closure: &RunnerV2BaseSourceClosureV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let cell_id = token(&format!("phase1-ledger-cell-{ledger_ordinal}"))?;
        let expected_result_root = if cell.expected_reason().is_none() {
            Some(phase_one_expected_result_root_v1(
                ledger_ordinal,
                &cell_id,
                cell.source_case_id(),
                cell.partition(),
                cell.expected_decision(),
                cell.root(),
                source_closure.snapshot_root(),
                close_manifest.root(),
                source_closure.root(),
            )?)
        } else {
            None
        };
        Ok(Self {
            ledger_ordinal,
            source_ordinal: cell.source_ordinal(),
            cell_id,
            source_case_id: cell.source_case_id().into(),
            source_class: cell.source_class(),
            source_path: cell.source_path().into(),
            group: cell.group(),
            facet: cell.facet(),
            execution_scope: cell.execution_scope(),
            partition: cell.partition(),
            expected_decision: cell.expected_decision(),
            expected_result_root,
            expected_reason: cell.expected_reason(),
            semantic_journey: token(RUNNER_V2_PHASE_ONE_SEMANTIC_JOURNEY_V1)?,
            route_id: token(RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1)?,
            route_owner: RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1.into(),
            compatible_source_snapshot_root: source_closure.snapshot_root(),
            case_manifest_root: close_manifest.root(),
            source_closure_root: source_closure.root(),
            downstream_contribution_root: cell
                .downstream_contribution()
                .map(crate::coverage::BaseCoverageCloseDownstreamContributionV1::root),
            five_explicits_root: cell.five_explicits().root(),
            cell_root: cell.root(),
            no_claim: token(RUNNER_V2_PHASE_ONE_LEDGER_ROW_NO_CLAIM_V1)?,
        })
    }

    /// One-based ordinal in the complete canonical Phase-1 ledger.
    #[must_use]
    pub const fn ledger_ordinal(&self) -> u32 {
        self.ledger_ordinal
    }

    /// One-based ordinal of the source close-manifest cell.
    #[must_use]
    pub const fn source_ordinal(&self) -> u32 {
        self.source_ordinal
    }

    /// Stable Phase-1 ledger-cell identifier.
    #[must_use]
    pub const fn cell_id(&self) -> &StableTokenV2 {
        &self.cell_id
    }

    /// Exact source-coverage case identifier.
    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    /// Exact frozen source-coverage class.
    #[must_use]
    pub const fn source_class(&self) -> BaseCoverageManifestClassV1 {
        self.source_class
    }

    /// Exact slash-separated source path.
    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
    }

    /// Exact AC53 close group.
    #[must_use]
    pub const fn group(&self) -> BaseCoverageCloseGroupV1 {
        self.group
    }

    /// Exact AC53 close facet.
    #[must_use]
    pub const fn facet(&self) -> BaseCoverageCloseFacetV1 {
        self.facet
    }

    /// Exact declaration-side execution scope.
    #[must_use]
    pub const fn execution_scope(&self) -> BaseCoverageCloseExecutionScopeV1 {
        self.execution_scope
    }

    /// Exact expected close partition.
    #[must_use]
    pub const fn partition(&self) -> BaseCoverageClosePartitionV1 {
        self.partition
    }

    /// Exact expected closed decision.
    #[must_use]
    pub const fn expected_decision(&self) -> BaseCoverageCloseDecisionV1 {
        self.expected_decision
    }

    /// Source-declared expected-decision contract root when the row has no reason.
    ///
    /// This is not an observed/output root. A downstream executor must derive
    /// its actual root from the real outcome and must never copy this value.
    #[must_use]
    pub const fn expected_result_root(&self) -> Option<ContentHash> {
        self.expected_result_root
    }

    /// Exact registered reason when an expected-decision root is inapplicable.
    #[must_use]
    pub const fn expected_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.expected_reason
    }

    /// Exact semantic-journey identifier.
    #[must_use]
    pub const fn semantic_journey(&self) -> &StableTokenV2 {
        &self.semantic_journey
    }

    /// Exact downstream Phase-1 route identifier.
    #[must_use]
    pub const fn route_id(&self) -> &StableTokenV2 {
        &self.route_id
    }

    /// Exact designated downstream route owner.
    #[must_use]
    pub fn route_owner(&self) -> &str {
        &self.route_owner
    }

    /// Exact compatible source-snapshot root.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> ContentHash {
        self.compatible_source_snapshot_root
    }

    /// Exact immutable source close-manifest root.
    #[must_use]
    pub const fn case_manifest_root(&self) -> ContentHash {
        self.case_manifest_root
    }

    /// Exact compiled source-closure root.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Compatibility V1 downstream-contribution root, when the source cell has one.
    #[must_use]
    pub const fn downstream_contribution_root(&self) -> Option<ContentHash> {
        self.downstream_contribution_root
    }

    /// Exact Five Explicits declaration root.
    #[must_use]
    pub const fn five_explicits_root(&self) -> ContentHash {
        self.five_explicits_root
    }

    /// Exact source close-manifest cell root.
    #[must_use]
    pub const fn cell_root(&self) -> ContentHash {
        self.cell_root
    }

    /// Exact result-free ledger-row no-claim.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "the expected-result identity visibly binds every non-reason Phase-1 ledger coordinate"
)]
fn phase_one_expected_result_root_v1(
    ledger_ordinal: u32,
    cell_id: &StableTokenV2,
    source_case_id: &str,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    close_cell_root: ContentHash,
    compatible_source_snapshot_root: ContentHash,
    case_manifest_root: ContentHash,
    source_closure_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    if ledger_ordinal == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "projection.phase_one.expected_result.ledger_ordinal",
            "one nonzero Phase-1 ledger ordinal",
            ledger_ordinal,
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1EXPECTEDRESULT\x01", 2 * 1024)?;
    frame.push_u32("expected_result.ledger_ordinal", ledger_ordinal)?;
    frame.push_str("expected_result.cell_id", cell_id.as_str())?;
    frame.push_str("expected_result.source_case_id", source_case_id)?;
    frame.push_u16("expected_result.partition", partition.code())?;
    frame.push_u16(
        "expected_result.expected_decision",
        expected_decision.code(),
    )?;
    frame.push_bytes(
        "expected_result.close_cell_root",
        close_cell_root.as_bytes(),
    )?;
    frame.push_bytes(
        "expected_result.compatible_source_snapshot_root",
        compatible_source_snapshot_root.as_bytes(),
    )?;
    frame.push_bytes(
        "expected_result.case_manifest_root",
        case_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "expected_result.source_closure_root",
        source_closure_root.as_bytes(),
    )?;
    frame.push_str(
        "expected_result.no_claim",
        "expected-result-root-is-a-source-declared-expected-decision-contract-not-an-observed-output-downstream-actual-roots-must-be-derived-from-real-outcomes-not-copied",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_EXPECTED_RESULT_DOMAIN_V1))
}

fn phase_one_partition_vocabulary_root_v1() -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1PARTITIONVOCAB\x01", 4 * 1024)?;
    frame.push_u32(
        "vocabulary.partition_count",
        checked_u32(
            RUNNER_V2_PHASE_ONE_EXPECTED_PARTITION_NAMES_V1.len(),
            "projection.phase_one.partition_vocabulary_count",
        )?,
    )?;
    for (index, name) in RUNNER_V2_PHASE_ONE_EXPECTED_PARTITION_NAMES_V1
        .iter()
        .enumerate()
    {
        let partition = BaseCoverageClosePartitionV1::ALL[index];
        if partition.stable_name() != *name {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.phase_one.partition_vocabulary",
                "the exact six-name closed partition vocabulary in code order",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        frame.push_u16("vocabulary.partition_code", partition.code())?;
        frame.push_str("vocabulary.partition_name", name)?;
    }
    frame.push_str(
        "vocabulary.no_claim",
        "partition-vocabulary-root-proves-closed-names-and-codes-not-ledger-membership-or-execution",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_PARTITION_VOCABULARY_DOMAIN_V1))
}

/// Source-authoritative, immutable, result-free Phase-1 expected ledger.
///
/// The downstream owner can inspect every exact row through this public
/// surface. The ledger contains no matched counts, observations, or execution
/// artifacts and cannot claim that the owner ran the designated route.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2PhaseOneExpectedLedgerV1 {
    source_manifest_root: ContentHash,
    reason_registry_root: ContentHash,
    close_manifest_root: ContentHash,
    compatible_source_snapshot_root: ContentHash,
    source_closure_root: ContentHash,
    partition_vocabulary_root: ContentHash,
    cells: Box<[RunnerV2PhaseOneExpectedLedgerCellV1]>,
    partition_cell_counts: [u32; 6],
    root: ContentHash,
}

impl RunnerV2PhaseOneExpectedLedgerV1 {
    /// Reconstruct the exact frozen ledger without performing execution.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let close_manifest = BaseCoverageCloseManifestV1::frozen()?;
        let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
        Self::build_from(&close_manifest, &source_closure)
    }

    fn build_from(
        close_manifest: &BaseCoverageCloseManifestV1,
        source_closure: &RunnerV2BaseSourceClosureV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let cells = phase_one_close_ledger_cells_for_source_v1(close_manifest, source_closure)?;
        let root = phase_one_expected_partitions_root_from_ledger_v1(
            close_manifest,
            source_closure,
            &cells,
        )?;
        let mut partition_cell_counts = [0_u32; 6];
        for cell in &cells {
            let index = usize::from(cell.partition.code() - 1);
            partition_cell_counts[index] = partition_cell_counts[index]
                .checked_add(1)
                .ok_or_else(sequence_overflow)?;
        }
        Ok(Self {
            source_manifest_root: close_manifest.source_manifest_root(),
            reason_registry_root: close_manifest.reason_registry_root(),
            close_manifest_root: close_manifest.root(),
            compatible_source_snapshot_root: source_closure.snapshot_root(),
            source_closure_root: source_closure.root(),
            partition_vocabulary_root: phase_one_partition_vocabulary_root_v1()?,
            cells: cells.into_boxed_slice(),
            partition_cell_counts,
            root,
        })
    }

    /// Exact source-coverage manifest root.
    #[must_use]
    pub const fn source_manifest_root(&self) -> ContentHash {
        self.source_manifest_root
    }

    /// Exact closed reason-registry root.
    #[must_use]
    pub const fn reason_registry_root(&self) -> ContentHash {
        self.reason_registry_root
    }

    /// Exact full-set close-manifest root.
    #[must_use]
    pub const fn close_manifest_root(&self) -> ContentHash {
        self.close_manifest_root
    }

    /// Exact compatible source-snapshot root shared by every row.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> ContentHash {
        self.compatible_source_snapshot_root
    }

    /// Exact compiled source-closure root shared by every row.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Exact independent six-name/code partition-vocabulary root.
    #[must_use]
    pub const fn partition_vocabulary_root(&self) -> ContentHash {
        self.partition_vocabulary_root
    }

    /// Complete ledger rows in canonical close-manifest order.
    #[must_use]
    pub fn cells(&self) -> &[RunnerV2PhaseOneExpectedLedgerCellV1] {
        &self.cells
    }

    /// Independently checked exact cell count for one closed partition.
    #[must_use]
    pub const fn partition_cell_count(&self, partition: BaseCoverageClosePartitionV1) -> u32 {
        self.partition_cell_counts[partition as usize - 1]
    }

    /// Canonical root binding all rows, vocabulary, and partition ID sets.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Runtime observations embedded in this result-free declaration.
    #[must_use]
    pub const fn execution_observation_count(&self) -> u32 {
        0
    }
}

fn phase_one_case_manifest_contract_root_v1(
    close_manifest: &BaseCoverageCloseManifestV1,
    expected_ledger_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    phase_one_case_manifest_contract_root_from_components_v1(
        close_manifest.source_manifest_root(),
        close_manifest.reason_registry_root(),
        close_manifest.root(),
        checked_u32(
            close_manifest.cells().len(),
            "projection.phase_one.close_manifest_cell_count",
        )?,
        expected_ledger_root,
    )
}

fn phase_one_case_manifest_contract_root_from_components_v1(
    source_manifest_root: ContentHash,
    reason_registry_root: ContentHash,
    close_manifest_root: ContentHash,
    close_manifest_cell_count: u32,
    expected_ledger_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    if close_manifest_cell_count == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "projection.phase_one.close_manifest_cell_count",
            "the nonzero exact full-set close-manifest cell count",
            close_manifest_cell_count,
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1CASEMANIFEST\x01", 4 * 1024)?;
    frame.push_str(
        "manifest.route_id",
        RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
    )?;
    frame.push_str(
        "manifest.execution_owner",
        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
    )?;
    frame.push_str("manifest.driver", RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1)?;
    frame.push_str(
        "manifest.posix_route",
        RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
    )?;
    frame.push_str(
        "manifest.windows_route",
        RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
    )?;
    frame.push_str(
        "manifest.path",
        RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
    )?;
    frame.push_bytes(
        "manifest.source_manifest_root",
        source_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "manifest.reason_registry_root",
        reason_registry_root.as_bytes(),
    )?;
    frame.push_bytes(
        "manifest.close_manifest_root",
        close_manifest_root.as_bytes(),
    )?;
    frame.push_u32(
        "manifest.close_manifest_cell_count",
        close_manifest_cell_count,
    )?;
    frame.push_bytes(
        "manifest.expected_ledger_root",
        expected_ledger_root.as_bytes(),
    )?;
    frame.push_u16("manifest.declaration_only", 1)?;
    frame.push_str(
        "manifest.no_claim",
        "source-exact-full-set-ledger-and-downstream-path-contract-are-not-file-content-or-execution-proof",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_CASE_MANIFEST_CONTRACT_DOMAIN_V1))
}

fn phase_one_close_ledger_cells_v1(
    close_manifest: &BaseCoverageCloseManifestV1,
) -> Result<Vec<RunnerV2PhaseOneExpectedLedgerCellV1>, ConstructionErrorV2> {
    let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
    phase_one_close_ledger_cells_for_source_v1(close_manifest, &source_closure)
}

fn phase_one_close_ledger_cells_for_source_v1(
    close_manifest: &BaseCoverageCloseManifestV1,
    source_closure: &RunnerV2BaseSourceClosureV1,
) -> Result<Vec<RunnerV2PhaseOneExpectedLedgerCellV1>, ConstructionErrorV2> {
    close_manifest
        .cells()
        .iter()
        .enumerate()
        .map(|(index, cell)| {
            RunnerV2PhaseOneExpectedLedgerCellV1::from_close_cell(
                checked_u32(index + 1, "projection.phase_one.close_ledger_ordinal")?,
                cell,
                close_manifest,
                source_closure,
            )
        })
        .collect()
}

fn phase_one_expected_partitions_root_from_ledger_v1(
    close_manifest: &BaseCoverageCloseManifestV1,
    source_closure: &RunnerV2BaseSourceClosureV1,
    cells: &[RunnerV2PhaseOneExpectedLedgerCellV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let expected_cells = close_manifest.cells();
    if cells.len() != expected_cells.len() {
        return Err(ConstructionErrorV2::new(
            if cells.len() < expected_cells.len() {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "projection.phase_one.close_ledger_cell_count",
            "one ledger row for every exact full-set close-manifest cell",
            cells.len(),
        ));
    }
    for (index, (presented, expected)) in cells.iter().zip(expected_cells).enumerate() {
        let expected_ordinal = checked_u32(index + 1, "projection.phase_one.close_ledger_ordinal")?;
        if presented.ledger_ordinal != expected_ordinal
            || presented.source_ordinal != expected_ordinal
            || presented.source_ordinal != expected.source_ordinal()
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "projection.phase_one.close_ledger_ordinal",
                "the exact contiguous source-manifest ordinal in canonical order",
                presented.ledger_ordinal,
            ));
        }
        let exact = RunnerV2PhaseOneExpectedLedgerCellV1::from_close_cell(
            expected_ordinal,
            expected,
            close_manifest,
            source_closure,
        )?;
        if presented != &exact {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.phase_one.close_ledger_cell",
                "the exact IDs, source, classification, result-or-reason, route, snapshot, manifests, contribution, Five Explicits, cell root, and no-claim",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }

    let mut presented_partition_counts = [0_u32; 6];
    for cell in cells {
        let index = usize::from(cell.partition.code() - 1);
        presented_partition_counts[index] = presented_partition_counts[index]
            .checked_add(1)
            .ok_or_else(sequence_overflow)?;
    }
    let mut independent_manifest_partition_counts = [0_u32; 6];
    for cell in expected_cells {
        let index = usize::from(cell.partition().code() - 1);
        independent_manifest_partition_counts[index] = independent_manifest_partition_counts[index]
            .checked_add(1)
            .ok_or_else(sequence_overflow)?;
    }
    if presented_partition_counts != independent_manifest_partition_counts {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.phase_one.close_ledger_partition_counts",
            "the independently recomputed exact per-partition cardinalities",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    for partition in BaseCoverageClosePartitionV1::ALL {
        let presented_ids = cells
            .iter()
            .filter(|cell| cell.partition == partition)
            .map(|cell| {
                (
                    cell.ledger_ordinal,
                    cell.cell_id.as_str().to_owned(),
                    cell.source_case_id.to_string(),
                )
            })
            .collect::<Vec<_>>();
        let independent_manifest_ids = expected_cells
            .iter()
            .enumerate()
            .filter(|(_, cell)| cell.partition() == partition)
            .map(|(index, cell)| {
                let ordinal =
                    checked_u32(index + 1, "projection.phase_one.partition_id_set_ordinal")?;
                Ok((
                    ordinal,
                    format!("phase1-ledger-cell-{ordinal}"),
                    cell.source_case_id().to_owned(),
                ))
            })
            .collect::<Result<Vec<_>, ConstructionErrorV2>>()?;
        if presented_ids != independent_manifest_ids {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.phase_one.partition_id_sets",
                "the independently recomputed exact ordered cell-ID set for each partition",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    let counted_total = presented_partition_counts
        .iter()
        .try_fold(0_u32, |total, count| {
            total.checked_add(*count).ok_or_else(sequence_overflow)
        })?;
    let exact_total = checked_u32(
        expected_cells.len(),
        "projection.phase_one.close_ledger_cell_count",
    )?;
    if counted_total != exact_total {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.phase_one.close_ledger_partition_sum",
            "an exact partition sum equal to the full close-manifest cell count",
            ConstructionObservedV2::unsigned_pair(u64::from(counted_total), u64::from(exact_total)),
        ));
    }

    let mut frame = CanonicalFrameV1::new(b"FSPHASE1PARTITIONS\x01", 2 * 1024 * 1024)?;
    frame.push_bytes(
        "partitions.source_manifest_root",
        close_manifest.source_manifest_root().as_bytes(),
    )?;
    frame.push_bytes(
        "partitions.reason_registry_root",
        close_manifest.reason_registry_root().as_bytes(),
    )?;
    frame.push_bytes(
        "partitions.close_manifest_root",
        close_manifest.root().as_bytes(),
    )?;
    frame.push_bytes(
        "partitions.compatible_source_snapshot_root",
        source_closure.snapshot_root().as_bytes(),
    )?;
    frame.push_bytes(
        "partitions.source_closure_root",
        source_closure.root().as_bytes(),
    )?;
    frame.push_bytes(
        "partitions.vocabulary_root",
        phase_one_partition_vocabulary_root_v1()?.as_bytes(),
    )?;
    frame.push_u32(
        "partitions.count",
        checked_u32(
            RUNNER_V2_PHASE_ONE_EXPECTED_PARTITION_NAMES_V1.len(),
            "projection.phase_one.partition_count",
        )?,
    )?;
    for (index, name) in RUNNER_V2_PHASE_ONE_EXPECTED_PARTITION_NAMES_V1
        .iter()
        .enumerate()
    {
        let partition = BaseCoverageClosePartitionV1::ALL[index];
        if partition.stable_name() != *name {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.phase_one.partition_vocabulary",
                "the exact six-name closed partition vocabulary in code order",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        frame.push_u16("partitions.code", partition.code())?;
        frame.push_str("partitions.name", name)?;
        frame.push_u32(
            "partitions.expected_cell_count",
            presented_partition_counts[index],
        )?;
        frame.push_u16("partitions.eventual_owner_requires_exact_match", 1)?;
    }
    frame.push_u32("partitions.cell_count", exact_total)?;
    for cell in cells {
        frame.push_u32("partitions.cell.ledger_ordinal", cell.ledger_ordinal)?;
        frame.push_u32("partitions.cell.ordinal", cell.source_ordinal)?;
        frame.push_str("partitions.cell.cell_id", cell.cell_id.as_str())?;
        frame.push_str("partitions.cell.stable_id", &cell.source_case_id)?;
        frame.push_u16("partitions.cell.source_class", cell.source_class as u16)?;
        frame.push_str("partitions.cell.source_path", &cell.source_path)?;
        frame.push_u16("partitions.cell.group", cell.group.code())?;
        frame.push_u16("partitions.cell.facet", cell.facet.code())?;
        frame.push_u16(
            "partitions.cell.execution_scope",
            cell.execution_scope.code(),
        )?;
        frame.push_u16("partitions.cell.partition", cell.partition.code())?;
        frame.push_u16(
            "partitions.cell.expected_decision",
            cell.expected_decision.code(),
        )?;
        frame.push_presence(
            "partitions.cell.expected_result_present",
            cell.expected_result_root.is_some(),
        )?;
        if let Some(root) = cell.expected_result_root {
            frame.push_bytes("partitions.cell.expected_result_root", root.as_bytes())?;
        }
        frame.push_presence(
            "partitions.cell.expected_reason_present",
            cell.expected_reason.is_some(),
        )?;
        if let Some(reason) = cell.expected_reason {
            frame.push_u16("partitions.cell.expected_reason", reason.code())?;
        }
        frame.push_str(
            "partitions.cell.semantic_journey",
            cell.semantic_journey.as_str(),
        )?;
        frame.push_str("partitions.cell.route_id", cell.route_id.as_str())?;
        frame.push_str("partitions.cell.route_owner", &cell.route_owner)?;
        frame.push_bytes(
            "partitions.cell.compatible_source_snapshot_root",
            cell.compatible_source_snapshot_root.as_bytes(),
        )?;
        frame.push_bytes(
            "partitions.cell.case_manifest_root",
            cell.case_manifest_root.as_bytes(),
        )?;
        frame.push_bytes(
            "partitions.cell.source_closure_root",
            cell.source_closure_root.as_bytes(),
        )?;
        frame.push_presence(
            "partitions.cell.downstream_contribution_present",
            cell.downstream_contribution_root.is_some(),
        )?;
        if let Some(root) = cell.downstream_contribution_root {
            frame.push_bytes(
                "partitions.cell.downstream_contribution_root",
                root.as_bytes(),
            )?;
        }
        frame.push_bytes(
            "partitions.cell.five_explicits_root",
            cell.five_explicits_root.as_bytes(),
        )?;
        frame.push_bytes("partitions.cell.root", cell.cell_root.as_bytes())?;
        frame.push_str("partitions.cell.no_claim", cell.no_claim.as_str())?;
    }
    frame.push_u32(
        "partitions.partition_id_set_count",
        checked_u32(
            BaseCoverageClosePartitionV1::ALL.len(),
            "projection.phase_one.partition_id_set_count",
        )?,
    )?;
    for partition in BaseCoverageClosePartitionV1::ALL {
        let members = cells
            .iter()
            .filter(|cell| cell.partition == partition)
            .collect::<Vec<_>>();
        frame.push_u16("partitions.id_set.partition", partition.code())?;
        frame.push_u32(
            "partitions.id_set.cell_count",
            checked_u32(
                members.len(),
                "projection.phase_one.partition_id_set_cell_count",
            )?,
        )?;
        for member in members {
            frame.push_u32("partitions.id_set.ledger_ordinal", member.ledger_ordinal)?;
            frame.push_str("partitions.id_set.cell_id", member.cell_id.as_str())?;
            frame.push_str("partitions.id_set.source_case_id", &member.source_case_id)?;
        }
    }
    frame.push_u32("partitions.required_unexpected_mismatches", 0)?;
    frame.push_u32("partitions.required_execution_failures", 0)?;
    frame.push_u32("partitions.required_unexplained_skips", 0)?;
    frame.push_str(
        "partitions.no_claim",
        "expected-partitions-do-not-prove-owner-execution-or-matched-results",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_EXPECTED_PARTITIONS_DOMAIN_V1))
}

fn phase_one_target_contract_root_v1() -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1TARGET\x01", 2 * 1024)?;
    frame.push_str(
        "target.posix_route",
        RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
    )?;
    frame.push_str(
        "target.windows_route",
        RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
    )?;
    frame.push_str(
        "target.no_claim",
        "portable-route-declaration-does-not-prove-platform-execution",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_TARGET_CONTRACT_DOMAIN_V1))
}

fn phase_one_feature_contract_root_v1() -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1FEATURES\x01", 1_024)?;
    frame.push_u32("features.required_count", 0)?;
    frame.push_str(
        "features.no_claim",
        "empty-required-feature-declaration-does-not-prove-owner-build-configuration",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_FEATURE_CONTRACT_DOMAIN_V1))
}

fn phase_one_build_contract_root_v1(
    source_closure: &RunnerV2BaseSourceClosureV1,
    schema_impact_projection: &RunnerV2BaseSchemaImpactProjectionV2,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1BUILD\x01", 2 * 1024)?;
    frame.push_bytes(
        "build.source_closure_root",
        source_closure.root().as_bytes(),
    )?;
    frame.push_bytes(
        "build.compatible_source_snapshot_root",
        source_closure.snapshot_root().as_bytes(),
    )?;
    frame.push_bytes(
        "build.schema_impact_projection_root",
        schema_impact_projection.root().as_bytes(),
    )?;
    frame.push_str(
        "build.no_claim",
        "build-input-declaration-does-not-prove-a-build-or-release-artifact",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_BUILD_CONTRACT_DOMAIN_V1))
}

fn phase_one_toolchain_contract_root_v1(
    source_closure: &RunnerV2BaseSourceClosureV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let toolchain_entry = source_closure
        .entries()
        .iter()
        .find(|entry| entry.path() == "rust-toolchain.toml")
        .ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "projection.phase_one.toolchain_source",
                "the exact rust-toolchain.toml source-closure entry",
                source_closure.entries().len(),
            )
        })?;
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1TOOLCHAIN\x01", 2 * 1024)?;
    frame.push_str("toolchain.source_path", toolchain_entry.path())?;
    frame.push_bytes(
        "toolchain.source_entry_root",
        toolchain_entry.entry_root().as_bytes(),
    )?;
    frame.push_bytes(
        "toolchain.compatible_source_snapshot_root",
        source_closure.snapshot_root().as_bytes(),
    )?;
    frame.push_str(
        "toolchain.no_claim",
        "toolchain-input-declaration-does-not-prove-owner-toolchain-execution",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_TOOLCHAIN_CONTRACT_DOMAIN_V1))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the immutable payload explicitly binds every currently implemented Phase-1 declaration input"
)]
fn phase_one_payload_root_v2(
    source_closure: &RunnerV2BaseSourceClosureV1,
    base_coverage_manifest_root: ContentHash,
    schema_impact_projection: &RunnerV2BaseSchemaImpactProjectionV2,
    expected_partitions_root: ContentHash,
    target_root: ContentHash,
    feature_set_root: ContentHash,
    schema_impact_log_schema_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1PAYLOAD\x02", 4 * 1024)?;
    frame.push_u16("payload.api_generation", 2)?;
    frame.push_u16("payload.wire_version", 1)?;
    frame.push_str("payload.route_id", RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1)?;
    frame.push_str(
        "payload.source_owner",
        RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1,
    )?;
    frame.push_bytes(
        "payload.source_closure_root",
        source_closure.root().as_bytes(),
    )?;
    frame.push_bytes(
        "payload.compatible_source_snapshot_root",
        source_closure.snapshot_root().as_bytes(),
    )?;
    frame.push_bytes(
        "payload.dependency_declaration_root",
        source_closure.dependency_declaration_root().as_bytes(),
    )?;
    frame.push_bytes(
        "payload.base_coverage_manifest_root",
        base_coverage_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "payload.schema_impact_projection_root",
        schema_impact_projection.root().as_bytes(),
    )?;
    frame.push_bytes(
        "payload.expected_partitions_root",
        expected_partitions_root.as_bytes(),
    )?;
    frame.push_bytes("payload.target_root", target_root.as_bytes())?;
    frame.push_bytes("payload.feature_set_root", feature_set_root.as_bytes())?;
    frame.push_bytes(
        "payload.schema_impact_log_schema_root",
        schema_impact_log_schema_root.as_bytes(),
    )?;
    frame.push_str(
        "payload.no_claim",
        "phase1-payload-is-result-free-and-does-not-complete-the-v2-execution-result-or-log-family",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_CONTRACT_PAYLOAD_DOMAIN_V2))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the aggregate root binds the contribution and separate envelope without introducing a reverse edge"
)]
fn phase_one_contribution_projection_root_v2(
    schema_impact_projection: &RunnerV2BaseSchemaImpactProjectionV2,
    source_closure_root: ContentHash,
    base_coverage_manifest_root: ContentHash,
    case_manifest_contract_root: ContentHash,
    expected_partitions_root: ContentHash,
    payload_root: ContentHash,
    contribution: &BaseCoverageCloseDownstreamContributionV2,
    deferred_envelope: &BaseCoverageCloseDeferredEvidenceEnvelopeV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    validate_exact_phase_one_observer_contract_v1(contribution.observer_contract())?;
    if contribution.source_owner() != RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1
        || contribution.payload_root() != payload_root
        || contribution.expected_partitions_root() != expected_partitions_root
        || contribution.schema_root() != schema_impact_projection.root()
        || contribution.log_schema_root() != crate::logging::schema_impact_log_schema_root_v1()?
        || contribution.observer_contract().case_manifest_root() != case_manifest_contract_root
        || contribution.retained_artifact_policy()
            != BaseCoverageCloseRetainedRelativeArtifactPolicyV1::OwnerEnvelopeRelativePathsOnly
        || contribution.no_claim().as_str() != RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1
        || deferred_envelope.contribution_root() != contribution.root()
        || deferred_envelope.observer_contract_root() != contribution.observer_contract().root()
        || deferred_envelope.evidence_kind()
            != crate::coverage::BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
        || deferred_envelope.disposition()
            != crate::coverage::RuntimeObservationDispositionV1::Deferred
        || deferred_envelope.deferred_reason()
            != DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        || deferred_envelope.retained_artifact_policy()
            != BaseCoverageCloseRetainedRelativeArtifactPolicyV1::OwnerEnvelopeRelativePathsOnly
        || deferred_envelope.no_claim().as_str()
            != RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1
        || contribution.execution_actual_field_count() != 0
        || deferred_envelope.execution_actual_field_count() != 0
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.phase_one.contribution_envelope_join",
            "the exact cycle-free contribution, observer, payload, manifest, and Deferred envelope join",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSPHASE1CONTRIBUTIONPROJECTION\x02", 4 * 1024)?;
    frame.push_u16("projection.api_generation", 2)?;
    frame.push_u16("projection.wire_version", 1)?;
    frame.push_str(
        "projection.route_id",
        RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
    )?;
    frame.push_bytes(
        "projection.schema_impact_projection_root",
        schema_impact_projection.root().as_bytes(),
    )?;
    frame.push_bytes(
        "projection.source_closure_root",
        source_closure_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.base_coverage_manifest_root",
        base_coverage_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.case_manifest_contract_root",
        case_manifest_contract_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.expected_partitions_root",
        expected_partitions_root.as_bytes(),
    )?;
    frame.push_bytes("projection.payload_root", payload_root.as_bytes())?;
    frame.push_bytes(
        "projection.contribution_root",
        contribution.root().as_bytes(),
    )?;
    frame.push_bytes(
        "projection.deferred_envelope_root",
        deferred_envelope.root().content_hash().as_bytes(),
    )?;
    frame.push_u32("projection.downstream_execution_count", 0)?;
    frame.push_str(
        "projection.no_claim",
        "immutable-phase1-contribution-projection-does-not-prove-downstream-execution-or-success",
    )?;
    Ok(frame.root(RUNNER_V2_PHASE_ONE_CONTRIBUTION_PROJECTION_DOMAIN_V2))
}

fn parse_phase_one_source_identity(
    root: ContentHash,
) -> Result<SourceIdentityRootV2, ConstructionErrorV2> {
    SourceIdentityRootV2::parse_presented(
        SourceIdentityRootV2::DESCRIPTOR.role(),
        SourceIdentityRootV2::DESCRIPTOR.domain(),
        &root.to_hex(),
    )
    .map_err(|_| phase_one_identity_error("projection.phase_one.source_root"))
}

fn parse_phase_one_build_identity(
    root: ContentHash,
) -> Result<BuildIdentityRootV2, ConstructionErrorV2> {
    BuildIdentityRootV2::parse_presented(
        BuildIdentityRootV2::DESCRIPTOR.role(),
        BuildIdentityRootV2::DESCRIPTOR.domain(),
        &root.to_hex(),
    )
    .map_err(|_| phase_one_identity_error("projection.phase_one.build_root"))
}

fn parse_phase_one_toolchain_identity(
    root: ContentHash,
) -> Result<ToolchainIdentityRootV2, ConstructionErrorV2> {
    ToolchainIdentityRootV2::parse_presented(
        ToolchainIdentityRootV2::DESCRIPTOR.role(),
        ToolchainIdentityRootV2::DESCRIPTOR.domain(),
        &root.to_hex(),
    )
    .map_err(|_| phase_one_identity_error("projection.phase_one.toolchain_root"))
}

fn phase_one_identity_error(field: &'static str) -> ConstructionErrorV2 {
    ConstructionErrorV2::new_redacted(
        ConstructionErrorKindV2::Incompatible,
        field,
        "one exact role- and domain-checked declaration identity",
        ConstructionObservedDataClassV2::CallerControlledText,
    )
}

fn checked_u32(value: usize, field: &'static str) -> Result<u32, ConstructionErrorV2> {
    u32::try_from(value).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            field,
            "a bounded u32 count",
            value,
        )
    })
}

/// Result-free source-coverage classes retained by the immutable manifest.
#[deprecated(
    note = "use BaseCoverageManifestClassV1 from the source-authoritative coverage manifest"
)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageClassV1 {
    /// Focused unit source cases.
    Unit = 1,
    /// Zero/below/exact/above and extrema source cases.
    Boundary = 2,
    /// Property and metamorphic source cases.
    PropertyMetamorphic = 3,
    /// API privacy, compile-fail, and doctest source cases.
    CompileFailDoctest = 4,
    /// Literal schema and descriptor source cases.
    SchemaDescriptor = 5,
    /// One-field and malformed-input mutation source cases.
    Mutation = 6,
    /// No-mock public-constructor integration source cases.
    Integration = 7,
    /// Five-journey mapped in-process E2E source cases.
    ProjectionE2e = 8,
    /// Deterministic structured-logging source cases.
    Logging = 9,
    /// Exact-set source-closure reconstruction source cases.
    SourceClosure = 10,
}

impl BaseCoverageClassV1 {
    /// Every coverage class in frozen order.
    pub const ALL: [Self; 10] = [
        Self::Unit,
        Self::Boundary,
        Self::PropertyMetamorphic,
        Self::CompileFailDoctest,
        Self::SchemaDescriptor,
        Self::Mutation,
        Self::Integration,
        Self::ProjectionE2e,
        Self::Logging,
        Self::SourceClosure,
    ];

    const fn code(self) -> u16 {
        self as u16
    }
}

#[allow(deprecated)]
impl ConstructionClosedSemanticV2 for BaseCoverageClassV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::Unit => "unit",
            Self::Boundary => "boundary",
            Self::PropertyMetamorphic => "property-metamorphic",
            Self::CompileFailDoctest => "compile-fail-doctest",
            Self::SchemaDescriptor => "schema-descriptor",
            Self::Mutation => "mutation",
            Self::Integration => "integration",
            Self::ProjectionE2e => "projection-e2e",
            Self::Logging => "logging",
            Self::SourceClosure => "source-closure",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct CoverageCaseTemplateV1 {
    class: BaseCoverageClassV1,
    id: &'static str,
    source_path: &'static str,
}

const COVERAGE_CASE_TEMPLATES_V1: [CoverageCaseTemplateV1; 44] = [
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-catalog",
        source_path: "crates/fs-evidence-runner/src/catalog.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-construction",
        source_path: "crates/fs-evidence-runner/src/construction.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-dependency",
        source_path: "crates/fs-evidence-runner/src/dependency.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-limits",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-logging",
        source_path: "crates/fs-evidence-runner/src/logging.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-projection",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-publication",
        source_path: "crates/fs-evidence-runner/src/publication.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Unit,
        id: "unit-value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Boundary,
        id: "boundary-value-extrema",
        source_path: "crates/fs-evidence-runner/src/value.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Boundary,
        id: "boundary-path-limits",
        source_path: "crates/fs-evidence-runner/src/path.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Boundary,
        id: "boundary-runner-limits",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Boundary,
        id: "boundary-runner-budgets",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Boundary,
        id: "boundary-diagnostic-envelope",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::PropertyMetamorphic,
        id: "property-rational-equivalence",
        source_path: "crates/fs-evidence-runner/src/value.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::PropertyMetamorphic,
        id: "property-path-permutation",
        source_path: "crates/fs-evidence-runner/src/path.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::PropertyMetamorphic,
        id: "property-identity-movement",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::PropertyMetamorphic,
        id: "property-limit-tightening",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::CompileFailDoctest,
        id: "compile-fail-sealed-domain",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::CompileFailDoctest,
        id: "compile-fail-immutable-closure",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SchemaDescriptor,
        id: "schema-catalog-tags",
        source_path: "crates/fs-evidence-runner/src/catalog.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SchemaDescriptor,
        id: "schema-limit-fields",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SchemaDescriptor,
        id: "schema-budget-fields",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SchemaDescriptor,
        id: "schema-identity-domains",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-unknown-catalog",
        source_path: "crates/fs-evidence-runner/src/catalog.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-noncanonical-value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-path-alias",
        source_path: "crates/fs-evidence-runner/src/path.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-identity-field",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-budget-field",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Mutation,
        id: "mutation-source-entry",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Integration,
        id: "integration-command-public-api",
        source_path: "crates/fs-evidence-runner/src/command.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Integration,
        id: "integration-publication-capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Integration,
        id: "integration-state-diagnostic",
        source_path: "crates/fs-evidence-runner/src/state.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Integration,
        id: "integration-budget-storage",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Logging,
        id: "logging-canonical-order",
        source_path: "crates/fs-evidence-runner/src/logging.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Logging,
        id: "logging-no-ambient-fields",
        source_path: "crates/fs-evidence-runner/src/logging.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::Logging,
        id: "logging-symbolic-reproduction",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SourceClosure,
        id: "source-closure-exact",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
    CoverageCaseTemplateV1 {
        class: BaseCoverageClassV1::SourceClosure,
        id: "source-closure-rejections",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
    },
];

/// Deprecated auxiliary compatibility view of one immutable source case.
///
/// This result-free view records no execution or pass/fail result and is not
/// the AC38 coverage authority. [`BaseCoverageManifestV1`] is the sole
/// source-authoritative coverage manifest.
#[deprecated(
    note = "use BaseCoverageManifestCaseV1 from the source-authoritative coverage manifest"
)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageSourceCaseV1 {
    class: BaseCoverageClassV1,
    ordinal: u32,
    id: StableTokenV2,
    source_path: LogicalBundlePathV1,
}

impl BaseCoverageSourceCaseV1 {
    /// Coverage class.
    #[must_use]
    pub const fn class(&self) -> BaseCoverageClassV1 {
        self.class
    }

    /// One-based ordinal within the coverage class.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Stable source-case ID.
    #[must_use]
    pub const fn id(&self) -> &StableTokenV2 {
        &self.id
    }

    /// Relative source file containing the case.
    #[must_use]
    pub const fn source_path(&self) -> &LogicalBundlePathV1 {
        &self.source_path
    }
}

/// Deprecated auxiliary compatibility inventory with no execution results.
///
/// This view is retained only for existing callers. It must not be used for
/// AC38 closure, test discovery, or pass/fail claims; [`BaseCoverageManifestV1`]
/// is the sole source of truth for those decisions.
#[deprecated(note = "use BaseCoverageManifestV1 as the source-authoritative coverage inventory")]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageInventoryV1 {
    cases: Box<[BaseCoverageSourceCaseV1]>,
    root: ContentHash,
}

impl BaseCoverageInventoryV1 {
    /// Source cases in exact class and within-class ordinal order.
    #[must_use]
    pub fn cases(&self) -> &[BaseCoverageSourceCaseV1] {
        &self.cases
    }

    /// Number of immutable source cases in one coverage class.
    #[must_use]
    pub fn source_case_count(&self, class: BaseCoverageClassV1) -> usize {
        self.cases
            .iter()
            .filter(|source_case| source_case.class == class)
            .count()
    }

    /// Domain-separated root of the exact result-free source inventory.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Exact downstream journey keys and sole-owner script mappings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseE2eJourneyV1 {
    /// Publication typestate and durability journey.
    PublicationState = 1,
    /// Cross-backend publication journey.
    PublicationV2 = 2,
    /// Independent verifier journey.
    VerifierV2 = 3,
    /// Canonical controller/process journey.
    CanonicalRunnerV2 = 4,
    /// Independent Runner-to-rjoq handoff journey.
    RjoqHandoffV1 = 5,
}

impl BaseE2eJourneyV1 {
    /// Every journey in frozen order.
    pub const ALL: [Self; 5] = [
        Self::PublicationState,
        Self::PublicationV2,
        Self::VerifierV2,
        Self::CanonicalRunnerV2,
        Self::RjoqHandoffV1,
    ];

    /// Exact non-wire journey tag.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable journey key.
    #[must_use]
    pub const fn key(self) -> &'static str {
        match self {
            Self::PublicationState => "publication-state-v2",
            Self::PublicationV2 => "publication-v2",
            Self::VerifierV2 => "verifier-v2",
            Self::CanonicalRunnerV2 => "canonical-runner-v2",
            Self::RjoqHandoffV1 => "rjoq-handoff-v1",
        }
    }

    /// Exact downstream-owned script path.
    #[must_use]
    pub const fn script_path(self) -> &'static str {
        match self {
            Self::PublicationState => "scripts/ci/e2e_evidence_runner_publication_state_v2.sh",
            Self::PublicationV2 => "scripts/ci/e2e_evidence_runner_publication_v2.sh",
            Self::VerifierV2 => "scripts/ci/e2e_evidence_verifier_v2.sh",
            Self::CanonicalRunnerV2 => "scripts/ci/canonical_evidence_runner_v2.sh",
            Self::RjoqHandoffV1 => "scripts/ci/verify_runner_rjoq_handoff_v1.sh",
        }
    }

    /// Exact downstream-owned release-driver name.
    #[must_use]
    pub const fn driver_name(self) -> &'static str {
        match self {
            Self::PublicationState => "e2e-evidence-runner-publication-state-v2-driver",
            Self::PublicationV2 => "e2e-evidence-runner-publication-v2-driver",
            Self::VerifierV2 => "e2e-evidence-verifier-v2-driver",
            Self::CanonicalRunnerV2 => "canonical-evidence-runner-v2-e2e-driver",
            Self::RjoqHandoffV1 => "verify-runner-rjoq-handoff-v1-driver",
        }
    }

    /// Exact downstream-owned immutable suite-manifest path.
    #[must_use]
    pub const fn manifest_path(self) -> &'static str {
        match self {
            Self::PublicationState => {
                "scripts/ci/manifests/evidence_runner_publication_state_v2_cases.v1.json"
            }
            Self::PublicationV2 => {
                "scripts/ci/manifests/evidence_runner_publication_v2_cases.v1.json"
            }
            Self::VerifierV2 => "scripts/ci/manifests/evidence_verifier_v2_cases.v1.json",
            Self::CanonicalRunnerV2 => {
                "scripts/ci/manifests/canonical_evidence_runner_v2_cases.v1.json"
            }
            Self::RjoqHandoffV1 => {
                "scripts/ci/manifests/runner_rjoq_handoff_verifier_v1_cases.v1.json"
            }
        }
    }

    /// Sole downstream Bead that owns the mapped release script.
    #[must_use]
    pub const fn downstream_owner(self) -> &'static str {
        match self {
            Self::PublicationState => "frankensim-epic-foundations-huq.24.2.2.2",
            Self::PublicationV2 => "frankensim-epic-foundations-huq.24.2.2.3.3",
            Self::VerifierV2 => "frankensim-epic-foundations-huq.24.3.3.3.3",
            Self::CanonicalRunnerV2 => "frankensim-epic-foundations-huq.24.4.1.4",
            Self::RjoqHandoffV1 => "frankensim-epic-foundations-huq.24.5.3.1",
        }
    }

    /// Explicit reason this downstream journey consumes its selected rows.
    #[must_use]
    pub const fn consumption_rationale(self) -> &'static str {
        match self {
            Self::PublicationState => "publication-state-contract-input",
            Self::PublicationV2 => "cross-backend-publication-contract-input",
            Self::VerifierV2 => "independent-verifier-contract-input",
            Self::CanonicalRunnerV2 => "canonical-controller-contract-input",
            Self::RjoqHandoffV1 => "runner-rjoq-handoff-contract-input",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eJourneyV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.key()
    }
}

/// Expected pure validation decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseE2eExpectedDecisionV1 {
    /// Constructor or validator must accept.
    Accept = 1,
    /// Constructor or validator must refuse.
    Refuse = 2,
    /// Platform-owned adjudication is explicitly unsupported locally.
    Unsupported = 3,
}

impl BaseE2eExpectedDecisionV1 {
    const fn code(self) -> u16 {
        self as u16
    }

    const fn name(self) -> &'static str {
        match self {
            Self::Accept => "accept",
            Self::Refuse => "refuse",
            Self::Unsupported => "unsupported",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eExpectedDecisionV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.name()
    }
}

/// Closed semantic case vocabulary shared by the five base journey manifests.
///
/// Each variant identifies one independently expected decision and a fixed
/// semantic-cell count; it does not carry an execution result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseE2eCaseKindV1 {
    /// Every frozen literal catalog and registered outer-tag round trip.
    CatalogLiterals = 1,
    /// Unknown base catalog code refusal.
    UnknownCatalogCode = 2,
    /// Canonical rational equivalence.
    CanonicalRational = 3,
    /// Stable-token upper-bound refusal.
    OverlongStableToken = 4,
    /// Portable logical-path admission.
    LogicalPath = 5,
    /// Reserved ContentStore key-prefix refusal.
    ReservedContentStorePrefix = 6,
    /// Explicit Windows non-ASCII alias unsupported cell.
    WindowsUnicodeAlias = 7,
    /// Complete limit descriptor/profile/refusal matrix.
    LimitCatalog = 8,
    /// Complete budget field/unit/profile matrix.
    BudgetAdmission = 9,
    /// Invalid child parallelism relation refusal.
    BudgetChildRelation = 10,
    /// Exact publication protocol/profile/mode matrix.
    PublicationSelection = 11,
    /// Cross-profile publication-selection refusal.
    PublicationCrossCell = 12,
    /// Exact least-privilege capability policy matrix.
    CapabilityLeastPrivilege = 13,
    /// One-right capability policy mutant matrix.
    CapabilityExtraRight = 14,
    /// Exhaustive lifecycle state and NotRun matrix.
    StatePass = 15,
    /// Usage-in-lifecycle refusal.
    StateUsageInLifecycle = 16,
    /// Diagnostic code/retryability/repair-kind matrix.
    Diagnostic = 17,
    /// Noncontiguous diagnostic repair-rank refusal.
    DiagnosticRankGap = 18,
    /// Nominal source/build/toolchain identity mutation matrix.
    IdentityMutation = 19,
    /// No-claim nominality and mutation matrix.
    NoClaimNominality = 20,
    /// Valid atomic command-result projection.
    AtomicResult = 21,
    /// Invalid durable-result presence projection.
    AtomicResultPresence = 22,
    /// Whole-publication stored-byte accounting.
    PublicationStorage = 23,
    /// Exact command selection and disposition matrix.
    CommandList = 24,
}

impl BaseE2eCaseKindV1 {
    /// Every case kind in exact tag order.
    pub const ALL: [Self; 24] = [
        Self::CatalogLiterals,
        Self::UnknownCatalogCode,
        Self::CanonicalRational,
        Self::OverlongStableToken,
        Self::LogicalPath,
        Self::ReservedContentStorePrefix,
        Self::WindowsUnicodeAlias,
        Self::LimitCatalog,
        Self::BudgetAdmission,
        Self::BudgetChildRelation,
        Self::PublicationSelection,
        Self::PublicationCrossCell,
        Self::CapabilityLeastPrivilege,
        Self::CapabilityExtraRight,
        Self::StatePass,
        Self::StateUsageInLifecycle,
        Self::Diagnostic,
        Self::DiagnosticRankGap,
        Self::IdentityMutation,
        Self::NoClaimNominality,
        Self::AtomicResult,
        Self::AtomicResultPresence,
        Self::PublicationStorage,
        Self::CommandList,
    ];

    /// Exact non-wire case-kind tag.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable semantic case name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::CatalogLiterals => "catalog-literals",
            Self::UnknownCatalogCode => "unknown-catalog-code",
            Self::CanonicalRational => "canonical-rational",
            Self::OverlongStableToken => "overlong-stable-token",
            Self::LogicalPath => "logical-path",
            Self::ReservedContentStorePrefix => "reserved-content-store-prefix",
            Self::WindowsUnicodeAlias => "windows-unicode-alias",
            Self::LimitCatalog => "limit-catalog",
            Self::BudgetAdmission => "budget-admission",
            Self::BudgetChildRelation => "budget-child-relation",
            Self::PublicationSelection => "publication-selection",
            Self::PublicationCrossCell => "publication-cross-cell",
            Self::CapabilityLeastPrivilege => "capability-least-privilege",
            Self::CapabilityExtraRight => "capability-extra-right",
            Self::StatePass => "state-pass",
            Self::StateUsageInLifecycle => "state-usage-in-lifecycle",
            Self::Diagnostic => "diagnostic",
            Self::DiagnosticRankGap => "diagnostic-rank-gap",
            Self::IdentityMutation => "identity-mutation",
            Self::NoClaimNominality => "no-claim-nominality",
            Self::AtomicResult => "atomic-result",
            Self::AtomicResultPresence => "atomic-result-presence",
            Self::PublicationStorage => "publication-storage",
            Self::CommandList => "command-list",
        }
    }

    /// Exact number of bounded semantic cells executed by this row.
    #[must_use]
    pub const fn semantic_cell_count(self) -> u32 {
        match self {
            Self::CatalogLiterals => 186,
            Self::LimitCatalog => 284,
            Self::BudgetAdmission => 44,
            Self::PublicationSelection => 6,
            Self::CapabilityLeastPrivilege => 12,
            Self::CapabilityExtraRight => 390,
            Self::StatePass => 32_460,
            Self::Diagnostic => 30,
            Self::IdentityMutation => 105,
            Self::NoClaimNominality => 5,
            Self::CommandList => 11,
            Self::UnknownCatalogCode
            | Self::CanonicalRational
            | Self::OverlongStableToken
            | Self::LogicalPath
            | Self::ReservedContentStorePrefix
            | Self::WindowsUnicodeAlias
            | Self::BudgetChildRelation
            | Self::PublicationCrossCell
            | Self::StateUsageInLifecycle
            | Self::DiagnosticRankGap
            | Self::AtomicResult
            | Self::AtomicResultPresence
            | Self::PublicationStorage => 1,
        }
    }

    /// Independently frozen semantic cells expected to accept.
    #[must_use]
    pub const fn positive_cell_count(self) -> u32 {
        match self {
            Self::CatalogLiterals => 186,
            Self::CanonicalRational
            | Self::LogicalPath
            | Self::AtomicResult
            | Self::PublicationStorage => 1,
            Self::LimitCatalog => 142,
            Self::BudgetAdmission => 36,
            Self::PublicationSelection => 6,
            Self::CapabilityLeastPrivilege => 12,
            Self::StatePass => 69,
            Self::Diagnostic => 30,
            Self::IdentityMutation => 96,
            Self::NoClaimNominality => 2,
            Self::CommandList => 11,
            Self::UnknownCatalogCode
            | Self::OverlongStableToken
            | Self::ReservedContentStorePrefix
            | Self::WindowsUnicodeAlias
            | Self::BudgetChildRelation
            | Self::PublicationCrossCell
            | Self::CapabilityExtraRight
            | Self::StateUsageInLifecycle
            | Self::DiagnosticRankGap
            | Self::AtomicResultPresence => 0,
        }
    }

    /// Independently frozen semantic cells expected to refuse.
    #[must_use]
    pub const fn expected_refusal_cell_count(self) -> u32 {
        match self {
            Self::UnknownCatalogCode
            | Self::OverlongStableToken
            | Self::ReservedContentStorePrefix
            | Self::BudgetChildRelation
            | Self::PublicationCrossCell
            | Self::StateUsageInLifecycle
            | Self::DiagnosticRankGap
            | Self::AtomicResultPresence => 1,
            Self::LimitCatalog => 142,
            Self::BudgetAdmission => 8,
            Self::CapabilityExtraRight => 390,
            Self::StatePass => 32_391,
            Self::IdentityMutation => 9,
            Self::NoClaimNominality => 3,
            Self::CatalogLiterals
            | Self::CanonicalRational
            | Self::LogicalPath
            | Self::WindowsUnicodeAlias
            | Self::PublicationSelection
            | Self::CapabilityLeastPrivilege
            | Self::Diagnostic
            | Self::AtomicResult
            | Self::PublicationStorage
            | Self::CommandList => 0,
        }
    }

    /// Independently frozen semantic cells expected to be unsupported.
    #[must_use]
    pub const fn unsupported_cell_count(self) -> u32 {
        match self {
            Self::WindowsUnicodeAlias => 1,
            _ => 0,
        }
    }

    /// Logical unit used for the row's checked-cell count or primary bound.
    #[must_use]
    pub const fn unit(self) -> LogicalUnitV2 {
        match self {
            Self::StatePass | Self::StateUsageInLifecycle => LogicalUnitV2::Records,
            Self::PublicationStorage | Self::AtomicResult | Self::AtomicResultPresence => {
                LogicalUnitV2::StoredBytes
            }
            _ => LogicalUnitV2::Count,
        }
    }

    /// Stable exact refusal/unsupported class, when the row does not accept.
    #[must_use]
    pub const fn expected_detail_code(self) -> Option<&'static str> {
        match self {
            Self::UnknownCatalogCode => Some("unknown-code"),
            Self::OverlongStableToken => Some("token-too-long"),
            Self::ReservedContentStorePrefix => Some("reserved-prefix"),
            Self::WindowsUnicodeAlias => Some("windows-unicode-alias-unsupported"),
            Self::BudgetChildRelation => Some("parallel-child-relation"),
            Self::PublicationCrossCell => Some("publication-cross-cell"),
            Self::CapabilityExtraRight => Some("least-privilege-rights"),
            Self::StateUsageInLifecycle => Some("usage-in-lifecycle"),
            Self::DiagnosticRankGap => Some("repair-rank-gap"),
            Self::AtomicResultPresence => Some("durable-result-presence"),
            _ => None,
        }
    }

    /// Closed fixture or subcase-manifest reference for this semantic row.
    #[must_use]
    pub const fn fixture_reference(self) -> &'static str {
        match self {
            Self::CatalogLiterals | Self::UnknownCatalogCode => "catalog-literal-oracle-v1",
            Self::CanonicalRational | Self::OverlongStableToken => "typed-value-oracle-v1",
            Self::LogicalPath | Self::ReservedContentStorePrefix | Self::WindowsUnicodeAlias => {
                "logical-path-oracle-v1"
            }
            Self::LimitCatalog => "limit-71-field-oracle-v1",
            Self::BudgetAdmission | Self::BudgetChildRelation => {
                "budget-18-field-16-unit-2-profile-8-refusal-oracle-v1"
            }
            Self::PublicationSelection | Self::PublicationCrossCell => {
                "publication-selection-oracle-v1"
            }
            Self::CapabilityLeastPrivilege | Self::CapabilityExtraRight => {
                "capability-12-cell-oracle-v1"
            }
            Self::StatePass | Self::StateUsageInLifecycle => "state-cartesian-oracle-v1",
            Self::Diagnostic | Self::DiagnosticRankGap => "diagnostic-repair-oracle-v1",
            Self::IdentityMutation => "identity-one-field-mutation-oracle-v1",
            Self::NoClaimNominality => "no-claim-nominality-oracle-v1",
            Self::AtomicResult | Self::AtomicResultPresence => "atomic-result-oracle-v1",
            Self::PublicationStorage => "publication-storage-oracle-v1",
            Self::CommandList => "command-table-oracle-v1",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eCaseKindV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.name()
    }
}

#[derive(Debug, Clone, Copy)]
struct BaseCaseTemplateV1 {
    id: &'static str,
    kind: BaseE2eCaseKindV1,
    expected: BaseE2eExpectedDecisionV1,
}

const BASE_CASE_TEMPLATES_V1: [BaseCaseTemplateV1; 24] = [
    BaseCaseTemplateV1 {
        id: "catalog-literals",
        kind: BaseE2eCaseKindV1::CatalogLiterals,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "unknown-catalog-code",
        kind: BaseE2eCaseKindV1::UnknownCatalogCode,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "canonical-rational",
        kind: BaseE2eCaseKindV1::CanonicalRational,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "overlong-stable-token",
        kind: BaseE2eCaseKindV1::OverlongStableToken,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "logical-path",
        kind: BaseE2eCaseKindV1::LogicalPath,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "reserved-content-store-prefix",
        kind: BaseE2eCaseKindV1::ReservedContentStorePrefix,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "windows-unicode-alias",
        kind: BaseE2eCaseKindV1::WindowsUnicodeAlias,
        expected: BaseE2eExpectedDecisionV1::Unsupported,
    },
    BaseCaseTemplateV1 {
        id: "limit-catalog",
        kind: BaseE2eCaseKindV1::LimitCatalog,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "budget-admission",
        kind: BaseE2eCaseKindV1::BudgetAdmission,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "budget-child-relation",
        kind: BaseE2eCaseKindV1::BudgetChildRelation,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "publication-selection",
        kind: BaseE2eCaseKindV1::PublicationSelection,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "publication-cross-cell",
        kind: BaseE2eCaseKindV1::PublicationCrossCell,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "capability-least-privilege",
        kind: BaseE2eCaseKindV1::CapabilityLeastPrivilege,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "capability-extra-right",
        kind: BaseE2eCaseKindV1::CapabilityExtraRight,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "state-pass",
        kind: BaseE2eCaseKindV1::StatePass,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "state-usage-in-lifecycle",
        kind: BaseE2eCaseKindV1::StateUsageInLifecycle,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "diagnostic",
        kind: BaseE2eCaseKindV1::Diagnostic,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "diagnostic-rank-gap",
        kind: BaseE2eCaseKindV1::DiagnosticRankGap,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "identity-mutation",
        kind: BaseE2eCaseKindV1::IdentityMutation,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "no-claim-nominality",
        kind: BaseE2eCaseKindV1::NoClaimNominality,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "atomic-result",
        kind: BaseE2eCaseKindV1::AtomicResult,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "atomic-result-presence",
        kind: BaseE2eCaseKindV1::AtomicResultPresence,
        expected: BaseE2eExpectedDecisionV1::Refuse,
    },
    BaseCaseTemplateV1 {
        id: "publication-storage",
        kind: BaseE2eCaseKindV1::PublicationStorage,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
    BaseCaseTemplateV1 {
        id: "command-list",
        kind: BaseE2eCaseKindV1::CommandList,
        expected: BaseE2eExpectedDecisionV1::Accept,
    },
];

const PUBLICATION_STATE_CASES_V1: &[BaseE2eCaseKindV1] = &[
    BaseE2eCaseKindV1::CatalogLiterals,
    BaseE2eCaseKindV1::UnknownCatalogCode,
    BaseE2eCaseKindV1::LogicalPath,
    BaseE2eCaseKindV1::ReservedContentStorePrefix,
    BaseE2eCaseKindV1::LimitCatalog,
    BaseE2eCaseKindV1::BudgetAdmission,
    BaseE2eCaseKindV1::BudgetChildRelation,
    BaseE2eCaseKindV1::PublicationSelection,
    BaseE2eCaseKindV1::PublicationCrossCell,
    BaseE2eCaseKindV1::CapabilityLeastPrivilege,
    BaseE2eCaseKindV1::CapabilityExtraRight,
    BaseE2eCaseKindV1::StatePass,
    BaseE2eCaseKindV1::StateUsageInLifecycle,
    BaseE2eCaseKindV1::Diagnostic,
    BaseE2eCaseKindV1::DiagnosticRankGap,
    BaseE2eCaseKindV1::IdentityMutation,
    BaseE2eCaseKindV1::NoClaimNominality,
    BaseE2eCaseKindV1::AtomicResult,
    BaseE2eCaseKindV1::AtomicResultPresence,
    BaseE2eCaseKindV1::PublicationStorage,
    BaseE2eCaseKindV1::CommandList,
];

const PUBLICATION_V2_CASES_V1: &[BaseE2eCaseKindV1] = &[
    BaseE2eCaseKindV1::CatalogLiterals,
    BaseE2eCaseKindV1::UnknownCatalogCode,
    BaseE2eCaseKindV1::LogicalPath,
    BaseE2eCaseKindV1::ReservedContentStorePrefix,
    BaseE2eCaseKindV1::WindowsUnicodeAlias,
    BaseE2eCaseKindV1::LimitCatalog,
    BaseE2eCaseKindV1::BudgetAdmission,
    BaseE2eCaseKindV1::PublicationSelection,
    BaseE2eCaseKindV1::PublicationCrossCell,
    BaseE2eCaseKindV1::CapabilityLeastPrivilege,
    BaseE2eCaseKindV1::CapabilityExtraRight,
    BaseE2eCaseKindV1::Diagnostic,
    BaseE2eCaseKindV1::DiagnosticRankGap,
    BaseE2eCaseKindV1::IdentityMutation,
    BaseE2eCaseKindV1::NoClaimNominality,
    BaseE2eCaseKindV1::AtomicResult,
    BaseE2eCaseKindV1::PublicationStorage,
    BaseE2eCaseKindV1::CommandList,
];

const VERIFIER_V2_CASES_V1: &[BaseE2eCaseKindV1] = &[
    BaseE2eCaseKindV1::CatalogLiterals,
    BaseE2eCaseKindV1::UnknownCatalogCode,
    BaseE2eCaseKindV1::CanonicalRational,
    BaseE2eCaseKindV1::OverlongStableToken,
    BaseE2eCaseKindV1::LogicalPath,
    BaseE2eCaseKindV1::ReservedContentStorePrefix,
    BaseE2eCaseKindV1::WindowsUnicodeAlias,
    BaseE2eCaseKindV1::LimitCatalog,
    BaseE2eCaseKindV1::BudgetAdmission,
    BaseE2eCaseKindV1::StatePass,
    BaseE2eCaseKindV1::StateUsageInLifecycle,
    BaseE2eCaseKindV1::Diagnostic,
    BaseE2eCaseKindV1::DiagnosticRankGap,
    BaseE2eCaseKindV1::IdentityMutation,
    BaseE2eCaseKindV1::NoClaimNominality,
    BaseE2eCaseKindV1::AtomicResult,
    BaseE2eCaseKindV1::AtomicResultPresence,
    BaseE2eCaseKindV1::PublicationStorage,
    BaseE2eCaseKindV1::CommandList,
];

const CANONICAL_RUNNER_V2_CASES_V1: &[BaseE2eCaseKindV1] = &BaseE2eCaseKindV1::ALL;

const RJOQ_HANDOFF_V1_CASES_V1: &[BaseE2eCaseKindV1] = &[
    BaseE2eCaseKindV1::CatalogLiterals,
    BaseE2eCaseKindV1::UnknownCatalogCode,
    BaseE2eCaseKindV1::CanonicalRational,
    BaseE2eCaseKindV1::OverlongStableToken,
    BaseE2eCaseKindV1::LimitCatalog,
    BaseE2eCaseKindV1::BudgetAdmission,
    BaseE2eCaseKindV1::BudgetChildRelation,
    BaseE2eCaseKindV1::PublicationSelection,
    BaseE2eCaseKindV1::CapabilityLeastPrivilege,
    BaseE2eCaseKindV1::StatePass,
    BaseE2eCaseKindV1::Diagnostic,
    BaseE2eCaseKindV1::IdentityMutation,
    BaseE2eCaseKindV1::NoClaimNominality,
    BaseE2eCaseKindV1::AtomicResult,
    BaseE2eCaseKindV1::AtomicResultPresence,
    BaseE2eCaseKindV1::CommandList,
];

// Independent literal row-ID oracles. These intentionally do not call
// `BaseE2eCaseKindV1::name`; drift between the semantic case arrays and these
// downstream-consumed IDs refuses during frozen projection construction.
const PUBLICATION_STATE_ROW_IDS_V1: &[&str] = &[
    "catalog-literals",
    "unknown-catalog-code",
    "logical-path",
    "reserved-content-store-prefix",
    "limit-catalog",
    "budget-admission",
    "budget-child-relation",
    "publication-selection",
    "publication-cross-cell",
    "capability-least-privilege",
    "capability-extra-right",
    "state-pass",
    "state-usage-in-lifecycle",
    "diagnostic",
    "diagnostic-rank-gap",
    "identity-mutation",
    "no-claim-nominality",
    "atomic-result",
    "atomic-result-presence",
    "publication-storage",
    "command-list",
];

const PUBLICATION_V2_ROW_IDS_V1: &[&str] = &[
    "catalog-literals",
    "unknown-catalog-code",
    "logical-path",
    "reserved-content-store-prefix",
    "windows-unicode-alias",
    "limit-catalog",
    "budget-admission",
    "publication-selection",
    "publication-cross-cell",
    "capability-least-privilege",
    "capability-extra-right",
    "diagnostic",
    "diagnostic-rank-gap",
    "identity-mutation",
    "no-claim-nominality",
    "atomic-result",
    "publication-storage",
    "command-list",
];

const VERIFIER_V2_ROW_IDS_V1: &[&str] = &[
    "catalog-literals",
    "unknown-catalog-code",
    "canonical-rational",
    "overlong-stable-token",
    "logical-path",
    "reserved-content-store-prefix",
    "windows-unicode-alias",
    "limit-catalog",
    "budget-admission",
    "state-pass",
    "state-usage-in-lifecycle",
    "diagnostic",
    "diagnostic-rank-gap",
    "identity-mutation",
    "no-claim-nominality",
    "atomic-result",
    "atomic-result-presence",
    "publication-storage",
    "command-list",
];

const CANONICAL_RUNNER_V2_ROW_IDS_V1: &[&str] = &[
    "catalog-literals",
    "unknown-catalog-code",
    "canonical-rational",
    "overlong-stable-token",
    "logical-path",
    "reserved-content-store-prefix",
    "windows-unicode-alias",
    "limit-catalog",
    "budget-admission",
    "budget-child-relation",
    "publication-selection",
    "publication-cross-cell",
    "capability-least-privilege",
    "capability-extra-right",
    "state-pass",
    "state-usage-in-lifecycle",
    "diagnostic",
    "diagnostic-rank-gap",
    "identity-mutation",
    "no-claim-nominality",
    "atomic-result",
    "atomic-result-presence",
    "publication-storage",
    "command-list",
];

const RJOQ_HANDOFF_V1_ROW_IDS_V1: &[&str] = &[
    "catalog-literals",
    "unknown-catalog-code",
    "canonical-rational",
    "overlong-stable-token",
    "limit-catalog",
    "budget-admission",
    "budget-child-relation",
    "publication-selection",
    "capability-least-privilege",
    "state-pass",
    "diagnostic",
    "identity-mutation",
    "no-claim-nominality",
    "atomic-result",
    "atomic-result-presence",
    "command-list",
];

const fn journey_case_kinds(journey: BaseE2eJourneyV1) -> &'static [BaseE2eCaseKindV1] {
    match journey {
        BaseE2eJourneyV1::PublicationState => PUBLICATION_STATE_CASES_V1,
        BaseE2eJourneyV1::PublicationV2 => PUBLICATION_V2_CASES_V1,
        BaseE2eJourneyV1::VerifierV2 => VERIFIER_V2_CASES_V1,
        BaseE2eJourneyV1::CanonicalRunnerV2 => CANONICAL_RUNNER_V2_CASES_V1,
        BaseE2eJourneyV1::RjoqHandoffV1 => RJOQ_HANDOFF_V1_CASES_V1,
    }
}

const fn journey_row_id_oracle(journey: BaseE2eJourneyV1) -> &'static [&'static str] {
    match journey {
        BaseE2eJourneyV1::PublicationState => PUBLICATION_STATE_ROW_IDS_V1,
        BaseE2eJourneyV1::PublicationV2 => PUBLICATION_V2_ROW_IDS_V1,
        BaseE2eJourneyV1::VerifierV2 => VERIFIER_V2_ROW_IDS_V1,
        BaseE2eJourneyV1::CanonicalRunnerV2 => CANONICAL_RUNNER_V2_ROW_IDS_V1,
        BaseE2eJourneyV1::RjoqHandoffV1 => RJOQ_HANDOFF_V1_ROW_IDS_V1,
    }
}

fn case_template(kind: BaseE2eCaseKindV1) -> &'static BaseCaseTemplateV1 {
    &BASE_CASE_TEMPLATES_V1[usize::from(kind.code() - 1)]
}

/// One immutable source-closed projection row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eProjectionRowV1 {
    id: StableTokenV2,
    kind: BaseE2eCaseKindV1,
    journey: BaseE2eJourneyV1,
    downstream_owner: StableTokenV2,
    downstream_driver: StableTokenV2,
    downstream_script: LogicalBundlePathV1,
    downstream_manifest: LogicalBundlePathV1,
    consumption_rationale: StableTokenV2,
    fixture_reference: StableTokenV2,
    expected: BaseE2eExpectedDecisionV1,
    expected_detail: Option<StableTokenV2>,
    semantic_cell_count: u32,
    positive_cell_count: u32,
    expected_refusal_cell_count: u32,
    unsupported_cell_count: u32,
    expected_detail_cell_count: u32,
    expected_detail_manifest_root: ContentHash,
    oracle_manifest_root: ContentHash,
    registered_decision_detail: Option<RegisteredDecisionDetailProjectionV2>,
    semantic_manifest_root: ContentHash,
    unit: LogicalUnitV2,
    no_claim_scope: StableTokenV2,
    source_closure_root: ContentHash,
    log_schema_root: ContentHash,
    mapping_root: ContentHash,
}

impl BaseE2eProjectionRowV1 {
    /// Stable row ID.
    #[must_use]
    pub const fn id(&self) -> &StableTokenV2 {
        &self.id
    }

    /// Expected pure decision.
    #[must_use]
    pub const fn expected(&self) -> BaseE2eExpectedDecisionV1 {
        self.expected
    }

    /// Public semantic case kind executed by the scoped runner.
    #[must_use]
    pub const fn kind(&self) -> BaseE2eCaseKindV1 {
        self.kind
    }

    /// Journey that consumes this exact mapping.
    #[must_use]
    pub const fn journey(&self) -> BaseE2eJourneyV1 {
        self.journey
    }

    /// Sole downstream Bead owner.
    #[must_use]
    pub const fn downstream_owner(&self) -> &StableTokenV2 {
        &self.downstream_owner
    }

    /// Exact downstream release-driver name.
    #[must_use]
    pub const fn downstream_driver(&self) -> &StableTokenV2 {
        &self.downstream_driver
    }

    /// Sole downstream script mapping.
    #[must_use]
    pub const fn downstream_script(&self) -> &LogicalBundlePathV1 {
        &self.downstream_script
    }

    /// Exact immutable downstream suite-manifest mapping.
    #[must_use]
    pub const fn downstream_manifest(&self) -> &LogicalBundlePathV1 {
        &self.downstream_manifest
    }

    /// Explicit journey-specific consumption rationale.
    #[must_use]
    pub const fn consumption_rationale(&self) -> &StableTokenV2 {
        &self.consumption_rationale
    }

    /// Closed fixture or subcase-manifest reference.
    #[must_use]
    pub const fn fixture_reference(&self) -> &StableTokenV2 {
        &self.fixture_reference
    }

    /// Exact refusal or unsupported detail code, when applicable.
    #[must_use]
    pub const fn expected_detail(&self) -> Option<&StableTokenV2> {
        self.expected_detail.as_ref()
    }

    /// Exact bounded semantic subcase count.
    #[must_use]
    pub const fn semantic_cell_count(&self) -> u32 {
        self.semantic_cell_count
    }

    /// Frozen expected-accept subcase count.
    #[must_use]
    pub const fn positive_cell_count(&self) -> u32 {
        self.positive_cell_count
    }

    /// Frozen expected-refusal subcase count.
    #[must_use]
    pub const fn expected_refusal_cell_count(&self) -> u32 {
        self.expected_refusal_cell_count
    }

    /// Frozen explicitly unsupported subcase count.
    #[must_use]
    pub const fn unsupported_cell_count(&self) -> u32 {
        self.unsupported_cell_count
    }

    /// Exact expected refusal plus unsupported detail-cell count.
    #[must_use]
    pub const fn expected_detail_cell_count(&self) -> u32 {
        self.expected_detail_cell_count
    }

    /// Ordered exact expected refusal/unsupported detail-manifest root.
    #[must_use]
    pub const fn expected_detail_manifest_root(&self) -> ContentHash {
        self.expected_detail_manifest_root
    }

    /// Domain-separated identity of the handwritten literal oracle consumed by
    /// this semantic row.
    #[must_use]
    pub const fn oracle_manifest_root(&self) -> ContentHash {
        self.oracle_manifest_root
    }

    /// Closed expected detail-manifest descriptor.
    #[must_use]
    pub const fn expected_detail_manifest(&self) -> BaseE2eDecisionDetailManifestV1 {
        BaseE2eDecisionDetailManifestV1 {
            cell_count: self.expected_detail_cell_count,
            root: self.expected_detail_manifest_root,
        }
    }

    /// Ordered, publicly inspectable independent detail-cell descriptors.
    #[must_use]
    pub fn expected_detail_cells(&self) -> &'static [BaseE2eDetailCellV1] {
        expected_detail_cells(self.kind)
    }

    /// Optional bounded registered-family detail reference carried by this
    /// containing row frame.
    ///
    /// The current base projection uses this only for the lane-neutral
    /// `case.conformance_mismatch` diagnostic. The reference is opaque and
    /// non-authoritative; the downstream family remains the sole owner of its
    /// content schema and lane semantics.
    #[must_use]
    pub const fn registered_decision_detail(
        &self,
    ) -> Option<&RegisteredDecisionDetailProjectionV2> {
        self.registered_decision_detail.as_ref()
    }

    /// Closed semantic-row manifest root.
    #[must_use]
    pub const fn semantic_manifest_root(&self) -> ContentHash {
        self.semantic_manifest_root
    }

    /// Logical unit associated with checked cells or the primary bound.
    #[must_use]
    pub const fn unit(&self) -> LogicalUnitV2 {
        self.unit
    }

    /// Stable no-claim classification for this pure validation row.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &StableTokenV2 {
        &self.no_claim_scope
    }

    /// Compiled source closure bound into this journey mapping.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Closed deterministic logging-schema root bound into this mapping.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    /// Domain-separated journey-specific row mapping root.
    #[must_use]
    pub const fn mapping_root(&self) -> ContentHash {
        self.mapping_root
    }
}

fn projection_row(
    journey: BaseE2eJourneyV1,
    downstream_owner: &StableTokenV2,
    downstream_driver: &StableTokenV2,
    downstream_script: &LogicalBundlePathV1,
    downstream_manifest: &LogicalBundlePathV1,
    source_closure_root: ContentHash,
    log_schema_root: ContentHash,
    kind: BaseE2eCaseKindV1,
) -> Result<BaseE2eProjectionRowV1, ConstructionErrorV2> {
    let template = case_template(kind);
    if template.kind != kind || template.id != kind.name() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.case_template",
            "the exact independently tagged case name and kind",
            ConstructionObservedV2::closed(&kind),
        ));
    }
    let id = token(template.id)?;
    let expected_detail = kind.expected_detail_code().map(token).transpose()?;
    let no_claim_scope = token("pure-base-validation-no-authority")?;
    let consumption_rationale = token(journey.consumption_rationale())?;
    let fixture_reference = token(kind.fixture_reference())?;
    let semantic_cell_count = kind.semantic_cell_count();
    let positive_cell_count = kind.positive_cell_count();
    let expected_refusal_cell_count = kind.expected_refusal_cell_count();
    let unsupported_cell_count = kind.unsupported_cell_count();
    let expected_detail_manifest = expected_detail_manifest(kind);
    let oracle_manifest_root = case_oracle_manifest_root(kind);
    let registered_decision_detail = registered_decision_detail_for_case(kind)?;
    let expected_detail_cell_count = expected_refusal_cell_count
        .checked_add(unsupported_cell_count)
        .ok_or_else(sequence_overflow)?;
    if expected_detail_manifest.cell_count != expected_detail_cell_count {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.expected_detail_manifest",
            "one exact detail for every expected-refusal or unsupported cell",
            ConstructionObservedV2::unsigned_pair(
                u64::from(expected_detail_manifest.cell_count),
                u64::from(expected_detail_cell_count),
            ),
        ));
    }
    let partition_count = positive_cell_count
        .checked_add(expected_refusal_cell_count)
        .and_then(|count| count.checked_add(unsupported_cell_count))
        .ok_or_else(sequence_overflow)?;
    if partition_count != semantic_cell_count {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.semantic_partition",
            "positive + expected-refusal + unsupported equals semantic cell count",
            ConstructionObservedV2::unsigned_pair(
                u64::from(partition_count),
                u64::from(semantic_cell_count),
            ),
        ));
    }
    let unit = kind.unit();
    let semantic_manifest_root = semantic_row_root(
        kind,
        template.expected,
        expected_detail.as_ref(),
        semantic_cell_count,
        positive_cell_count,
        expected_refusal_cell_count,
        unsupported_cell_count,
        expected_detail_cell_count,
        expected_detail_manifest.root,
        oracle_manifest_root,
        registered_decision_detail.as_ref(),
        unit,
        &no_claim_scope,
    )?;
    let mapping_root = journey_row_root(
        journey,
        downstream_owner,
        downstream_driver,
        downstream_script,
        downstream_manifest,
        &consumption_rationale,
        &fixture_reference,
        source_closure_root,
        log_schema_root,
        semantic_manifest_root,
    )?;
    Ok(BaseE2eProjectionRowV1 {
        id,
        kind,
        journey,
        downstream_owner: downstream_owner.clone(),
        downstream_driver: downstream_driver.clone(),
        downstream_script: downstream_script.clone(),
        downstream_manifest: downstream_manifest.clone(),
        consumption_rationale,
        fixture_reference,
        expected: template.expected,
        expected_detail,
        semantic_cell_count,
        positive_cell_count,
        expected_refusal_cell_count,
        unsupported_cell_count,
        expected_detail_cell_count,
        expected_detail_manifest_root: expected_detail_manifest.root,
        oracle_manifest_root,
        registered_decision_detail,
        semantic_manifest_root,
        unit,
        no_claim_scope,
        source_closure_root,
        log_schema_root,
        mapping_root,
    })
}

#[allow(
    clippy::too_many_arguments,
    reason = "the semantic row root intentionally binds every independently frozen row field"
)]
fn semantic_row_root(
    kind: BaseE2eCaseKindV1,
    expected: BaseE2eExpectedDecisionV1,
    expected_detail: Option<&StableTokenV2>,
    semantic_cell_count: u32,
    positive_cell_count: u32,
    expected_refusal_cell_count: u32,
    unsupported_cell_count: u32,
    expected_detail_cell_count: u32,
    expected_detail_manifest_root: ContentHash,
    oracle_manifest_root: ContentHash,
    registered_decision_detail: Option<&RegisteredDecisionDetailProjectionV2>,
    unit: LogicalUnitV2,
    no_claim_scope: &StableTokenV2,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASESEMANTICROW\x01", 2048)?;
    frame.push_u16("semantic.algorithm_version", 1)?;
    frame.push_u16("semantic.case_kind", kind.code())?;
    frame.push_str("semantic.case_name", kind.name())?;
    frame.push_u16("semantic.expected", expected.code())?;
    frame.push_u16(
        "semantic.expected_detail_presence",
        u16::from(expected_detail.is_some()),
    )?;
    if let Some(expected_detail) = expected_detail {
        frame.push_str("semantic.expected_detail", expected_detail.as_str())?;
    }
    frame.push_u32("semantic.cell_count", semantic_cell_count)?;
    frame.push_u32("semantic.positive_cell_count", positive_cell_count)?;
    frame.push_u32(
        "semantic.expected_refusal_cell_count",
        expected_refusal_cell_count,
    )?;
    frame.push_u32("semantic.unsupported_cell_count", unsupported_cell_count)?;
    frame.push_u32(
        "semantic.expected_detail_cell_count",
        expected_detail_cell_count,
    )?;
    frame.push_bytes(
        "semantic.expected_detail_manifest_root",
        expected_detail_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "semantic.oracle_manifest_root",
        oracle_manifest_root.as_bytes(),
    )?;
    frame.push_u16(
        "semantic.registered_decision_detail_presence",
        u16::from(registered_decision_detail.is_some()),
    )?;
    if let Some(detail) = registered_decision_detail {
        frame.push_bytes(
            "semantic.registered_decision_detail_root",
            detail.root().as_bytes(),
        )?;
        frame.push_bytes(
            "semantic.registered_decision_detail_registry_root",
            detail.registry_root().as_bytes(),
        )?;
        frame.push_u16(
            "semantic.registered_decision_detail_namespace",
            detail.namespace().code(),
        )?;
        frame.push_u16(
            "semantic.registered_decision_detail_code",
            detail.detail_code(),
        )?;
        frame.push_bytes(
            "semantic.registered_decision_detail_content_root",
            detail.content_root().as_bytes(),
        )?;
        frame.push_u32(
            "semantic.registered_decision_detail_encoded_length",
            detail.encoded_length(),
        )?;
    }
    frame.push_u16("semantic.unit_tag", unit.tag())?;
    if let Some(registered_id) = unit.registered_id() {
        frame.push_u16("semantic.unit_registered_id", registered_id)?;
    }
    frame.push_str("semantic.no_claim_scope", no_claim_scope.as_str())?;
    Ok(frame.root(BASE_E2E_SEMANTIC_ROW_DOMAIN_V1))
}

fn registered_decision_detail_for_case(
    kind: BaseE2eCaseKindV1,
) -> Result<Option<RegisteredDecisionDetailProjectionV2>, ConstructionErrorV2> {
    const NAMESPACE: u16 = 7;
    const DETAIL_CODE: u16 = 1;
    const OPAQUE_FIXTURE: &[u8] = b"lane-neutral-case-conformance-detail-reference-v1";

    if kind != BaseE2eCaseKindV1::Diagnostic {
        return Ok(None);
    }
    let registry = DecisionDetailNamespaceRegistryV2::frozen();
    let descriptor = registry.lookup_registered_family(NAMESPACE)?;
    if descriptor.stable_name() != "case-conformance-detail"
        || descriptor.owner() != "frankensim-epic-foundations-huq.24.1.1.3.1"
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.registered_decision_detail",
            "namespace 7 case-conformance-detail owned by 24.1.1.3.1",
            ConstructionObservedV2::closed(&descriptor),
        ));
    }
    let encoded_length = u32::try_from(OPAQUE_FIXTURE.len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_projection.registered_decision_detail",
            "a bounded u32 opaque detail reference",
            OPAQUE_FIXTURE.len(),
        )
    })?;
    RegisteredDecisionDetailProjectionV2::new(
        &registry,
        NAMESPACE,
        DETAIL_CODE,
        hash_domain(BASE_E2E_REGISTERED_DETAIL_FIXTURE_DOMAIN_V1, OPAQUE_FIXTURE),
        encoded_length,
    )
    .map(Some)
}

fn case_oracle_manifest_root_from_table_root(
    kind: BaseE2eCaseKindV1,
    owner: &str,
    literal_table_root: ContentHash,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(256);
    bytes.extend_from_slice(&kind.code().to_be_bytes());
    detail_push_str(&mut bytes, kind.name());
    detail_push_str(&mut bytes, owner);
    bytes.extend_from_slice(literal_table_root.as_bytes());
    hash_domain(BASE_E2E_ORACLE_MANIFEST_DOMAIN_V1, &bytes)
}

fn case_oracle_manifest_root(kind: BaseE2eCaseKindV1) -> ContentHash {
    let (owner, literal_table_root) = match kind {
        BaseE2eCaseKindV1::CatalogLiterals | BaseE2eCaseKindV1::UnknownCatalogCode => {
            ("catalog", catalog_literal_oracle_root())
        }
        BaseE2eCaseKindV1::LimitCatalog => ("limits", limit_oracle_table_root(limit_oracle_rows())),
        BaseE2eCaseKindV1::BudgetAdmission | BaseE2eCaseKindV1::BudgetChildRelation => {
            ("budgets", budget_oracle_table_root())
        }
        BaseE2eCaseKindV1::Diagnostic | BaseE2eCaseKindV1::DiagnosticRankGap => {
            ("diagnostics", diagnostic_oracle_table_root())
        }
        BaseE2eCaseKindV1::CommandList => ("commands", command_oracle_table_root()),
        _ => (
            "fixture-reference",
            hash_domain(
                "org.frankensim.fs-evidence-runner.base-e2e-fixture-oracle-table.v1",
                kind.fixture_reference().as_bytes(),
            ),
        ),
    };
    case_oracle_manifest_root_from_table_root(kind, owner, literal_table_root)
}

#[allow(
    clippy::too_many_arguments,
    reason = "the journey row root intentionally binds every downstream mapping and source-closure field"
)]
fn journey_row_root(
    journey: BaseE2eJourneyV1,
    downstream_owner: &StableTokenV2,
    downstream_driver: &StableTokenV2,
    downstream_script: &LogicalBundlePathV1,
    downstream_manifest: &LogicalBundlePathV1,
    consumption_rationale: &StableTokenV2,
    fixture_reference: &StableTokenV2,
    source_closure_root: ContentHash,
    log_schema_root: ContentHash,
    semantic_manifest_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEJOURNEYROW\x01", 4096)?;
    frame.push_u16("mapping.journey", journey.code())?;
    frame.push_str("mapping.downstream_owner", downstream_owner.as_str())?;
    frame.push_str("mapping.downstream_driver", downstream_driver.as_str())?;
    frame.push_str("mapping.downstream_script", downstream_script.as_str())?;
    frame.push_str("mapping.downstream_manifest", downstream_manifest.as_str())?;
    frame.push_str(
        "mapping.consumption_rationale",
        consumption_rationale.as_str(),
    )?;
    frame.push_str("mapping.fixture_reference", fixture_reference.as_str())?;
    frame.push_bytes(
        "mapping.source_closure_root",
        source_closure_root.as_bytes(),
    )?;
    frame.push_bytes("mapping.log_schema_root", log_schema_root.as_bytes())?;
    frame.push_bytes(
        "mapping.semantic_manifest_root",
        semantic_manifest_root.as_bytes(),
    )?;
    Ok(frame.root(BASE_E2E_JOURNEY_ROW_DOMAIN_V1))
}

/// One journey-keyed immutable projection and root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eJourneyProjectionV1 {
    journey: BaseE2eJourneyV1,
    downstream_owner: StableTokenV2,
    driver_name: StableTokenV2,
    script_path: LogicalBundlePathV1,
    manifest_path: LogicalBundlePathV1,
    rows: Box<[BaseE2eProjectionRowV1]>,
    source_closure_root: ContentHash,
    log_schema_root: ContentHash,
    root: ContentHash,
}

impl BaseE2eJourneyProjectionV1 {
    /// Journey.
    #[must_use]
    pub const fn journey(&self) -> BaseE2eJourneyV1 {
        self.journey
    }

    /// Sole-owner downstream script path.
    #[must_use]
    pub const fn script_path(&self) -> &LogicalBundlePathV1 {
        &self.script_path
    }

    /// Sole downstream Bead that owns the mapped script.
    #[must_use]
    pub const fn downstream_owner(&self) -> &StableTokenV2 {
        &self.downstream_owner
    }

    /// Exact downstream-owned release-driver name.
    #[must_use]
    pub const fn driver_name(&self) -> &StableTokenV2 {
        &self.driver_name
    }

    /// Exact downstream-owned immutable suite-manifest path.
    #[must_use]
    pub const fn manifest_path(&self) -> &LogicalBundlePathV1 {
        &self.manifest_path
    }

    /// Exact row set.
    #[must_use]
    pub fn rows(&self) -> &[BaseE2eProjectionRowV1] {
        &self.rows
    }

    /// Immutable journey projection root.
    #[deprecated(note = "use manifest_root to distinguish manifests from execution results")]
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Explicitly named immutable journey-manifest root.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.root
    }

    /// Compile-time base source closure bound into this journey.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Closed deterministic logging-schema root.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }
}

/// Complete five-journey, non-wire projection manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseE2eProjectionV1 {
    journeys: Box<[BaseE2eJourneyProjectionV1]>,
    source_closure: RunnerV2BaseSourceClosureV1,
    coverage_inventory: BaseCoverageInventoryV1,
    coverage_manifest: BaseCoverageManifestV1,
    close_manifest: BaseCoverageCloseManifestV1,
    log_schema_root: ContentHash,
    root: ContentHash,
}

impl RunnerV2BaseE2eProjectionV1 {
    /// Construct the exact five journey projections and roots.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
        let log_schema_root = crate::logging::base_e2e_log_schema_root_v1()?;
        let mut journeys = Vec::with_capacity(BaseE2eJourneyV1::ALL.len());
        for journey in BaseE2eJourneyV1::ALL {
            let downstream_owner = token(journey.downstream_owner())?;
            let driver_name = token(journey.driver_name())?;
            let script_path =
                LogicalBundlePathV1::new(journey.script_path()).map_err(|_error| {
                    ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_projection.script_path",
                        "the frozen logical relative script path",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    )
                })?;
            let manifest_path =
                LogicalBundlePathV1::new(journey.manifest_path()).map_err(|_error| {
                    ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_projection.manifest_path",
                        "the frozen logical relative immutable suite-manifest path",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    )
                })?;
            let rows = journey_case_kinds(journey)
                .iter()
                .copied()
                .map(|kind| {
                    projection_row(
                        journey,
                        &downstream_owner,
                        &driver_name,
                        &script_path,
                        &manifest_path,
                        source_closure.root(),
                        log_schema_root,
                        kind,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            let row_id_oracle = journey_row_id_oracle(journey);
            if rows.len() != row_id_oracle.len()
                || rows
                    .iter()
                    .zip(row_id_oracle)
                    .any(|(row, expected_id)| row.id().as_str() != *expected_id)
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_projection.journey_row_oracle",
                    "the exact independent literal row-ID sequence",
                    ConstructionObservedV2::closed(&journey),
                ));
            }
            let root = journey_root(
                journey,
                &downstream_owner,
                &driver_name,
                &script_path,
                &manifest_path,
                &rows,
                source_closure.root(),
                log_schema_root,
            )?;
            journeys.push(BaseE2eJourneyProjectionV1 {
                journey,
                downstream_owner,
                driver_name,
                script_path,
                manifest_path,
                rows: rows.into_boxed_slice(),
                source_closure_root: source_closure.root(),
                log_schema_root,
                root,
            });
        }
        let coverage_inventory = coverage_inventory(&journeys)?;
        let coverage_manifest = exact_coverage_manifest(&journeys)?;
        let close_manifest = BaseCoverageCloseManifestV1::frozen()?;
        if close_manifest.source_manifest_root() != coverage_manifest.root() {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_projection.close_source_manifest_root",
                "the exact source-authoritative full coverage manifest root",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let root = projection_root(
            &journeys,
            source_closure.root(),
            coverage_inventory.root(),
            coverage_manifest.root(),
            close_manifest.root(),
            log_schema_root,
        )?;
        Ok(Self {
            journeys: journeys.into_boxed_slice(),
            source_closure,
            coverage_inventory,
            coverage_manifest,
            close_manifest,
            log_schema_root,
            root,
        })
    }

    /// Five projections in frozen journey order.
    #[must_use]
    pub fn journeys(&self) -> &[BaseE2eJourneyProjectionV1] {
        &self.journeys
    }

    /// Exact compile-time source closure bound into this projection.
    #[must_use]
    pub const fn source_closure(&self) -> &RunnerV2BaseSourceClosureV1 {
        &self.source_closure
    }

    /// Immutable result-free source coverage inventory.
    ///
    /// This deprecated auxiliary view is not the AC38 coverage authority; use
    /// [`Self::coverage_manifest`] for source-authoritative coverage.
    #[deprecated(note = "use coverage_manifest as the AC38 source of truth")]
    #[must_use]
    pub const fn coverage_inventory(&self) -> &BaseCoverageInventoryV1 {
        &self.coverage_inventory
    }

    /// Source-authoritative two-stage coverage manifest.
    #[must_use]
    pub const fn coverage_manifest(&self) -> &BaseCoverageManifestV1 {
        &self.coverage_manifest
    }

    /// Source-authoritative, result-free, full-set-only AC53 close manifest.
    ///
    /// This projection binds the manifest but does not synthesize a close
    /// report: unit, compile-fail, and downstream-owned execution remain
    /// separately presented proof obligations.
    #[must_use]
    pub const fn close_manifest(&self) -> &BaseCoverageCloseManifestV1 {
        &self.close_manifest
    }

    /// Closed deterministic logging-schema root bound by every journey.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    /// Complete immutable projection root.
    #[deprecated(note = "use manifest_root to distinguish manifests from execution results")]
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Explicitly named complete projection-manifest root.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.root
    }
}

/// Semantic build inputs retained in every deterministic projection log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eHarnessIdentityV1 {
    source: SourceIdentityRootV2,
    build: BuildIdentityRootV2,
    toolchain: ToolchainIdentityRootV2,
    target: StableTokenV2,
    target_root: ContentHash,
    features: Box<[StableTokenV2]>,
    feature_set_root: ContentHash,
    no_claim_scope: NoClaimScopeRootV1,
    context_root: ContentHash,
}

impl BaseE2eHarnessIdentityV1 {
    /// Validate a duplicate-free canonical feature set.
    pub fn new(
        source: SourceIdentityRootV2,
        build: BuildIdentityRootV2,
        toolchain: ToolchainIdentityRootV2,
        target: StableTokenV2,
        mut features: Vec<StableTokenV2>,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let mut seen = std::collections::BTreeSet::new();
        for feature in &features {
            if !seen.insert(feature.as_str()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_harness.features",
                    "unique stable feature tokens",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        features.sort();
        let target_root = base_e2e_target_root(&target)?;
        let feature_set_root = base_e2e_feature_set_root(&features)?;
        let context_root = base_e2e_harness_context_root(
            &source,
            &build,
            &toolchain,
            &target,
            target_root,
            &features,
            feature_set_root,
            &no_claim_scope,
        )?;
        Ok(Self {
            source,
            build,
            toolchain,
            target,
            target_root,
            features: features.into_boxed_slice(),
            feature_set_root,
            no_claim_scope,
            context_root,
        })
    }

    /// Caller-presented source identity; this is context, not live-tree proof.
    #[must_use]
    pub const fn source(&self) -> &SourceIdentityRootV2 {
        &self.source
    }

    /// Caller-presented build identity; this is context, not derived proof.
    #[must_use]
    pub const fn build(&self) -> &BuildIdentityRootV2 {
        &self.build
    }

    /// Caller-presented toolchain identity.
    #[must_use]
    pub const fn toolchain(&self) -> &ToolchainIdentityRootV2 {
        &self.toolchain
    }

    /// Exact target token.
    #[must_use]
    pub const fn target(&self) -> &StableTokenV2 {
        &self.target
    }

    /// Domain-separated target root.
    #[must_use]
    pub const fn target_root(&self) -> ContentHash {
        self.target_root
    }

    /// Exact canonical feature tokens.
    #[must_use]
    pub fn features(&self) -> &[StableTokenV2] {
        &self.features
    }

    /// Domain-separated root of the exact feature set.
    #[must_use]
    pub const fn feature_set_root(&self) -> ContentHash {
        self.feature_set_root
    }

    /// Presented no-claim scope.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Root binding every presented harness-context field.
    #[must_use]
    pub const fn context_root(&self) -> ContentHash {
        self.context_root
    }
}

fn base_e2e_target_root(target: &StableTokenV2) -> Result<ContentHash, ConstructionErrorV2> {
    crate::logging::base_e2e_target_root_v1(target)
}

fn base_e2e_feature_set_root(
    features: &[StableTokenV2],
) -> Result<ContentHash, ConstructionErrorV2> {
    crate::logging::base_e2e_feature_set_root_v1(features)
}

#[allow(
    clippy::too_many_arguments,
    reason = "the frame deliberately binds every Five-Explicits build-context component"
)]
fn base_e2e_harness_context_root(
    source: &SourceIdentityRootV2,
    build: &BuildIdentityRootV2,
    toolchain: &ToolchainIdentityRootV2,
    target: &StableTokenV2,
    target_root: ContentHash,
    features: &[StableTokenV2],
    feature_set_root: ContentHash,
    no_claim_scope: &NoClaimScopeRootV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEHARNESSCONTEXT\x01", 16 * 1024)?;
    frame.push_bytes("harness.presented_source", source.bytes())?;
    frame.push_bytes("harness.presented_build", build.bytes())?;
    frame.push_bytes("harness.presented_toolchain", toolchain.bytes())?;
    frame.push_str("harness.target", target.as_str())?;
    frame.push_bytes("harness.target_root", target_root.as_bytes())?;
    frame.push_u32(
        "harness.feature_count",
        u32::try_from(features.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "harness.feature_count",
                "a u32 canonical feature count",
                features.len(),
            )
        })?,
    )?;
    for feature in features {
        frame.push_str("harness.feature", feature.as_str())?;
    }
    frame.push_bytes("harness.feature_set_root", feature_set_root.as_bytes())?;
    frame.push_bytes("harness.no_claim_scope", no_claim_scope.bytes())?;
    Ok(frame.root(BASE_E2E_HARNESS_CONTEXT_DOMAIN_V1))
}

/// Closed retained-artifact claim for the pure in-process base projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseE2eRetainedArtifactClaimV1 {
    /// No artifact is retained; source-closure paths are source inventory, not
    /// execution artifacts.
    Absent,
}

impl BaseE2eRetainedArtifactClaimV1 {
    /// Whether the typed claim is explicit absence.
    #[must_use]
    pub const fn is_absent(self) -> bool {
        matches!(self, Self::Absent)
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eRetainedArtifactClaimV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::Absent => "absent",
        }
    }
}

/// Exact aggregate projection execution counts and typed log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eProjectionReportV1 {
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    projection_rows_checked: u32,
    projection_e2e_checked: u32,
    logging_events_checked: u32,
    source_closure_positive_eligible: u32,
    source_closure_positive_matched: u32,
    source_closure_expected_refusals: u32,
    source_closure_expected_refusals_matched: u32,
    source_closure_unexpected_mismatches: u32,
    projection_root: ContentHash,
    source_closure_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    source_closure_paths: Box<[LogicalBundlePathV1]>,
    journey_executions: Box<[BaseE2eJourneyExecutionReportV1]>,
    retained_artifact_claim: BaseE2eRetainedArtifactClaimV1,
    execution_root: ContentHash,
    log: BaseE2eLogV1,
    coverage_report: BaseCoverageCheckedReportV1,
}

impl BaseE2eProjectionReportV1 {
    /// Rows independently expected to accept.
    #[must_use]
    pub const fn positive_eligible(&self) -> u32 {
        self.positive_eligible
    }

    /// Expected-accept rows whose observed decision and cell count matched.
    #[must_use]
    pub const fn positive_matched(&self) -> u32 {
        self.positive_matched
    }

    /// Rows independently expected to refuse.
    #[must_use]
    pub const fn expected_refusals(&self) -> u32 {
        self.expected_refusals
    }

    /// Expected-refusal rows whose observed decision and cell count matched.
    #[must_use]
    pub const fn expected_refusals_matched(&self) -> u32 {
        self.expected_refusals_matched
    }

    /// Exactly matched, explicitly unsupported platform-owned rows.
    #[must_use]
    pub const fn unsupported(&self) -> u32 {
        self.unsupported
    }

    /// Any unexpected decision or semantic-cell-count mismatch.
    #[must_use]
    pub const fn unexpected_mismatches(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Exact manifest rows joined to checked results.
    #[must_use]
    pub const fn projection_rows_checked(&self) -> u32 {
        self.projection_rows_checked
    }

    /// Projection-E2E subcases actually checked across all five journeys.
    #[must_use]
    pub const fn projection_e2e_checked(&self) -> u32 {
        self.projection_e2e_checked
    }

    /// Deterministic structured log events validated and retained.
    #[must_use]
    pub const fn logging_events_checked(&self) -> u32 {
        self.logging_events_checked
    }

    /// Positive source-closure checks executed by this report.
    #[must_use]
    pub const fn source_closure_positive_eligible(&self) -> u32 {
        self.source_closure_positive_eligible
    }

    /// Positive source-closure checks whose observed decision matched.
    #[must_use]
    pub const fn source_closure_positive_matched(&self) -> u32 {
        self.source_closure_positive_matched
    }

    /// Expected source-closure refusals executed by this report.
    #[must_use]
    pub const fn source_closure_expected_refusals(&self) -> u32 {
        self.source_closure_expected_refusals
    }

    /// Expected source-closure refusals whose observed decision matched.
    #[must_use]
    pub const fn source_closure_expected_refusals_matched(&self) -> u32 {
        self.source_closure_expected_refusals_matched
    }

    /// Unexpected source-closure mismatches.
    #[must_use]
    pub const fn source_closure_unexpected_mismatches(&self) -> u32 {
        self.source_closure_unexpected_mismatches
    }

    /// Backward-compatible total locally adjudicable rows.
    #[deprecated(
        note = "use positive_eligible and expected_refusals to preserve the exact partitions"
    )]
    #[must_use]
    pub const fn eligible(&self) -> u32 {
        self.positive_eligible + self.expected_refusals
    }

    /// Backward-compatible total matched locally adjudicable rows.
    #[deprecated(
        note = "use positive_matched and expected_refusals_matched to preserve the exact partitions"
    )]
    #[must_use]
    pub const fn passed(&self) -> u32 {
        self.positive_matched + self.expected_refusals_matched
    }

    /// Backward-compatible alias for unexpected mismatches.
    #[deprecated(note = "use unexpected_mismatches")]
    #[must_use]
    pub const fn failed(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Backward-compatible total locally adjudicable source-closure checks.
    #[deprecated(
        note = "use source_closure_positive_eligible and source_closure_expected_refusals"
    )]
    #[must_use]
    pub const fn source_closure_eligible(&self) -> u32 {
        self.source_closure_positive_eligible + self.source_closure_expected_refusals
    }

    /// Backward-compatible total matched source-closure checks.
    #[deprecated(
        note = "use source_closure_positive_matched and source_closure_expected_refusals_matched"
    )]
    #[must_use]
    pub const fn source_closure_passed(&self) -> u32 {
        self.source_closure_positive_matched + self.source_closure_expected_refusals_matched
    }

    /// Backward-compatible alias for unexpected source-closure mismatches.
    #[deprecated(note = "use source_closure_unexpected_mismatches")]
    #[must_use]
    pub const fn source_closure_failed(&self) -> u32 {
        self.source_closure_unexpected_mismatches
    }

    /// Projection root executed.
    #[deprecated(note = "use manifest_root to distinguish manifests from execution results")]
    #[must_use]
    pub const fn projection_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Explicitly named immutable projection-manifest root.
    ///
    /// This is the canonical name for [`Self::projection_root`], retained as a
    /// compatibility alias.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Exact source-closure root reconstructed before row execution.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    /// Ordered five typed journey executions retained by the aggregate.
    #[must_use]
    pub fn journey_executions(&self) -> &[BaseE2eJourneyExecutionReportV1] {
        &self.journey_executions
    }

    /// Explicit typed absent-retained-artifact claim.
    #[must_use]
    pub const fn retained_artifact_claim(&self) -> BaseE2eRetainedArtifactClaimV1 {
        self.retained_artifact_claim
    }

    /// Compatibility-shaped artifact accessor. The pure projection never
    /// retains an execution artifact, so this is always typed absence.
    #[deprecated(note = "use retained_artifact_claim for the typed absence contract")]
    #[must_use]
    pub const fn retained_artifact(&self) -> Option<&LogicalBundlePathV1> {
        None
    }

    /// Context-bound aggregate execution root over all five ordered journey
    /// execution roots and the typed absent-artifact claim.
    #[must_use]
    pub const fn execution_root(&self) -> ContentHash {
        self.execution_root
    }

    /// Presented source identity retained by the harness.
    #[must_use]
    pub const fn source_root(&self) -> &SourceIdentityRootV2 {
        &self.source_root
    }

    /// Presented build identity retained by the harness.
    #[must_use]
    pub const fn build_root(&self) -> &BuildIdentityRootV2 {
        &self.build_root
    }

    /// Exact relative source paths bound by the compile-time closure.
    #[must_use]
    pub fn source_closure_paths(&self) -> &[LogicalBundlePathV1] {
        &self.source_closure_paths
    }

    /// Deterministic detailed log.
    #[must_use]
    pub const fn log(&self) -> &BaseE2eLogV1 {
        &self.log
    }

    /// Exact-joined in-process coverage results for the cells this run owns.
    #[must_use]
    pub const fn coverage_report(&self) -> &BaseCoverageCheckedReportV1 {
        &self.coverage_report
    }

    /// Whether every exact local and source-closure partition matched and both
    /// independently reconciled subordinate reports are green.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        self.positive_matched == self.positive_eligible
            && self.expected_refusals_matched == self.expected_refusals
            && self.unexpected_mismatches == 0
            && self.source_closure_positive_matched == self.source_closure_positive_eligible
            && self.source_closure_expected_refusals_matched
                == self.source_closure_expected_refusals
            && self.source_closure_unexpected_mismatches == 0
            && self.log.is_green()
            && self.coverage_report.is_green()
    }
}

/// Closed descriptor for one ordered exact refusal/unsupported detail
/// manifest.
///
/// The root binds the case kind, exact detail-cell count, semantic ordinals,
/// cell identities, expected or observed decisions, and typed detail payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseE2eDecisionDetailManifestV1 {
    cell_count: u32,
    root: ContentHash,
}

/// Concise public name for the closed per-case decision-detail manifest.
///
/// The longer original name remains available for source compatibility.
pub type BaseE2eDetailManifestV1 = BaseE2eDecisionDetailManifestV1;

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eDetailExecutionV1 {
    expected: BaseE2eDecisionDetailManifestV1,
    observed: BaseE2eDecisionDetailManifestV1,
    expected_cells: Option<Box<[BaseE2eDetailCellV1]>>,
    observed_cells: Option<Box<[BaseE2eDetailCellV1]>>,
    matched_cells: u32,
    first_divergent_cell: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eCaseExecutionV1 {
    decision: BaseE2eExpectedDecisionV1,
    checked_cells: u32,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    first_failed_cell: Option<String>,
    detail: BaseE2eDetailExecutionV1,
}

impl BaseE2eCaseExecutionV1 {
    fn accepted(checked_cells: u32, detail: BaseE2eDetailExecutionV1) -> Self {
        assert!(
            checked_cells > 0,
            "an accepted case execution must observe at least one semantic cell"
        );
        Self {
            decision: BaseE2eExpectedDecisionV1::Accept,
            checked_cells,
            positive_eligible: checked_cells,
            positive_matched: checked_cells,
            expected_refusals: 0,
            expected_refusals_matched: 0,
            unsupported: 0,
            unexpected_mismatches: 0,
            first_failed_cell: None,
            detail,
        }
    }

    fn refused(checked_cells: u32, detail: BaseE2eDetailExecutionV1) -> Self {
        assert!(
            checked_cells > 0,
            "a refused case execution must observe at least one semantic cell"
        );
        Self {
            decision: BaseE2eExpectedDecisionV1::Refuse,
            checked_cells,
            positive_eligible: 0,
            positive_matched: 0,
            expected_refusals: checked_cells,
            expected_refusals_matched: checked_cells,
            unsupported: 0,
            unexpected_mismatches: 0,
            first_failed_cell: None,
            detail,
        }
    }

    fn unsupported(checked_cells: u32, detail: BaseE2eDetailExecutionV1) -> Self {
        assert!(
            checked_cells > 0,
            "an unsupported case execution must observe at least one semantic cell"
        );
        Self {
            decision: BaseE2eExpectedDecisionV1::Unsupported,
            checked_cells,
            positive_eligible: 0,
            positive_matched: 0,
            expected_refusals: 0,
            expected_refusals_matched: 0,
            unsupported: checked_cells,
            unexpected_mismatches: 0,
            first_failed_cell: None,
            detail,
        }
    }

    fn mixed(
        positive_eligible: u32,
        expected_refusals: u32,
        unsupported: u32,
        detail: BaseE2eDetailExecutionV1,
    ) -> Self {
        let checked_cells = positive_eligible
            .checked_add(expected_refusals)
            .and_then(|count| count.checked_add(unsupported))
            .expect("the frozen mixed semantic-cell inventory fits u32");
        assert!(
            checked_cells > 0,
            "a mixed case execution must observe at least one semantic cell"
        );
        Self {
            decision: BaseE2eExpectedDecisionV1::Accept,
            checked_cells,
            positive_eligible,
            positive_matched: positive_eligible,
            expected_refusals,
            expected_refusals_matched: expected_refusals,
            unsupported,
            unexpected_mismatches: 0,
            first_failed_cell: None,
            detail,
        }
    }

    fn with_failure(
        decision: BaseE2eExpectedDecisionV1,
        checked_cells: u32,
        first_failed_cell: impl Into<String>,
        detail: BaseE2eDetailExecutionV1,
    ) -> Self {
        assert!(
            checked_cells > 0,
            "a failed case execution must identify a real observed semantic cell"
        );
        Self {
            decision,
            checked_cells,
            positive_eligible: checked_cells,
            positive_matched: checked_cells - 1,
            expected_refusals: 0,
            expected_refusals_matched: 0,
            unsupported: 0,
            unexpected_mismatches: 1,
            first_failed_cell: Some(first_failed_cell.into()),
            detail,
        }
    }
}

/// One expected/matched partition in a presented semantic-cell observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseE2eMatchedPartitionV1 {
    eligible: u32,
    matched: u32,
}

impl BaseE2eMatchedPartitionV1 {
    /// Construct one bounded partition.
    pub fn new(eligible: u32, matched: u32) -> Result<Self, ConstructionErrorV2> {
        if matched > eligible {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_result.partition",
                "matched cells no greater than eligible cells",
                ConstructionObservedV2::unsigned_pair(u64::from(matched), u64::from(eligible)),
            ));
        }
        Ok(Self { eligible, matched })
    }

    /// Eligible cells in this partition.
    #[must_use]
    pub const fn eligible(self) -> u32 {
        self.eligible
    }

    /// Cells whose observations matched.
    #[must_use]
    pub const fn matched(self) -> u32 {
        self.matched
    }
}

/// Intrinsically valid observed semantic-cell partitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseE2eObservedCountsV1 {
    positive: BaseE2eMatchedPartitionV1,
    expected_refusals: BaseE2eMatchedPartitionV1,
    unsupported: u32,
    unexpected_mismatches: u32,
    checked_cells: u32,
}

impl BaseE2eObservedCountsV1 {
    /// Construct exact partitions and reconcile their mismatch count.
    pub fn new(
        positive: BaseE2eMatchedPartitionV1,
        expected_refusals: BaseE2eMatchedPartitionV1,
        unsupported: u32,
        unexpected_mismatches: u32,
    ) -> Result<Self, ConstructionErrorV2> {
        let positive_gap = positive
            .eligible
            .checked_sub(positive.matched)
            .ok_or_else(sequence_overflow)?;
        let refusal_gap = expected_refusals
            .eligible
            .checked_sub(expected_refusals.matched)
            .ok_or_else(sequence_overflow)?;
        let expected_mismatches = positive_gap
            .checked_add(refusal_gap)
            .ok_or_else(sequence_overflow)?;
        if unexpected_mismatches != expected_mismatches {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_result.unexpected_mismatches",
                "the exact positive and expected-refusal partition gaps",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(unexpected_mismatches),
                    u64::from(expected_mismatches),
                ),
            ));
        }
        let checked_cells = positive
            .eligible
            .checked_add(expected_refusals.eligible)
            .and_then(|count| count.checked_add(unsupported))
            .ok_or_else(sequence_overflow)?;
        if checked_cells == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "base_e2e_result.checked_cells",
                "at least one checked semantic cell",
                0,
            ));
        }
        Ok(Self {
            positive,
            expected_refusals,
            unsupported,
            unexpected_mismatches,
            checked_cells,
        })
    }

    /// Expected-accept partition.
    #[must_use]
    pub const fn positive(self) -> BaseE2eMatchedPartitionV1 {
        self.positive
    }

    /// Expected-refusal partition.
    #[must_use]
    pub const fn expected_refusals(self) -> BaseE2eMatchedPartitionV1 {
        self.expected_refusals
    }

    /// Exactly matched typed unsupported cells.
    #[must_use]
    pub const fn unsupported(self) -> u32 {
        self.unsupported
    }

    /// Exact mismatch gap across both locally adjudicable partitions.
    #[must_use]
    pub const fn unexpected_mismatches(self) -> u32 {
        self.unexpected_mismatches
    }

    /// Reconciled checked-cell count.
    #[must_use]
    pub const fn checked_cells(self) -> u32 {
        self.checked_cells
    }
}

/// One caller-presented row observation before exact manifest joining.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2ePresentedRowResultV1 {
    journey: BaseE2eJourneyV1,
    row_id: StableTokenV2,
    semantic_manifest_root: ContentHash,
    observed: BaseE2eExpectedDecisionV1,
    counts: BaseE2eObservedCountsV1,
    typed_detail_cells_presented: bool,
    observed_detail_manifest_root: ContentHash,
    observed_detail_cell_count: u32,
    detail_cells_matched: u32,
    first_unexpected_cell: Option<StableTokenV2>,
    first_observed_detail_cell: Option<BaseE2eDetailCellV1>,
    root: ContentHash,
}

impl BaseE2ePresentedRowResultV1 {
    /// Construct one intrinsically valid, non-authoritative compatibility
    /// observation without typed detail cells.
    #[deprecated(
        since = "0.1.0",
        note = "use new_with_observed_detail_cells so kind, order, count, root, and the first typed divergence are checked"
    )]
    pub fn new(
        journey: BaseE2eJourneyV1,
        row_id: StableTokenV2,
        semantic_manifest_root: ContentHash,
        observed: BaseE2eExpectedDecisionV1,
        counts: BaseE2eObservedCountsV1,
        first_unexpected_cell: Option<StableTokenV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new_with_detail_manifest(
            journey,
            row_id,
            semantic_manifest_root,
            observed,
            counts,
            BaseE2eDecisionDetailManifestV1::empty(BaseE2eCaseKindV1::CatalogLiterals).root,
            0,
            0,
            first_unexpected_cell,
        )
    }

    /// Construct one observation carrying its ordered exact decision-detail
    /// manifest. The legacy [`Self::new`] constructor deliberately carries an
    /// empty, unmatched manifest and therefore cannot satisfy a row that owns
    /// refusal or unsupported detail cells.
    #[deprecated(
        since = "0.1.0",
        note = "use new_with_observed_detail_cells; this compatibility constructor cannot validate or retain typed observed cells"
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_detail_manifest(
        journey: BaseE2eJourneyV1,
        row_id: StableTokenV2,
        semantic_manifest_root: ContentHash,
        observed: BaseE2eExpectedDecisionV1,
        counts: BaseE2eObservedCountsV1,
        observed_detail_manifest_root: ContentHash,
        observed_detail_cell_count: u32,
        detail_cells_matched: u32,
        first_unexpected_cell: Option<StableTokenV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new_with_detail_manifest_and_first_observed_cell(
            journey,
            row_id,
            semantic_manifest_root,
            observed,
            counts,
            observed_detail_manifest_root,
            observed_detail_cell_count,
            detail_cells_matched,
            first_unexpected_cell,
            None,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_with_detail_manifest_and_first_observed_cell(
        journey: BaseE2eJourneyV1,
        row_id: StableTokenV2,
        semantic_manifest_root: ContentHash,
        observed: BaseE2eExpectedDecisionV1,
        counts: BaseE2eObservedCountsV1,
        observed_detail_manifest_root: ContentHash,
        observed_detail_cell_count: u32,
        detail_cells_matched: u32,
        first_unexpected_cell: Option<StableTokenV2>,
        first_observed_detail_cell: Option<BaseE2eDetailCellV1>,
        typed_detail_cells_presented: bool,
    ) -> Result<Self, ConstructionErrorV2> {
        if detail_cells_matched > observed_detail_cell_count {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_result.detail_cells",
                "matched detail cells no greater than observed detail cells",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(detail_cells_matched),
                    u64::from(observed_detail_cell_count),
                ),
            ));
        }
        if first_unexpected_cell.is_some() != (counts.unexpected_mismatches() > 0) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_result.first_unexpected_cell",
                "presence exactly when unexpected_mismatches is nonzero",
                first_unexpected_cell.is_some(),
            ));
        }
        let root = presented_row_result_root(
            journey,
            &row_id,
            semantic_manifest_root,
            observed,
            counts,
            observed_detail_manifest_root,
            observed_detail_cell_count,
            detail_cells_matched,
            first_unexpected_cell.as_ref(),
            first_observed_detail_cell.as_ref(),
            typed_detail_cells_presented,
        )?;
        Ok(Self {
            journey,
            row_id,
            semantic_manifest_root,
            observed,
            counts,
            typed_detail_cells_presented,
            observed_detail_manifest_root,
            observed_detail_cell_count,
            detail_cells_matched,
            first_unexpected_cell,
            first_observed_detail_cell,
            root,
        })
    }

    /// Construct an observation from one already assembled closed detail
    /// descriptor.
    ///
    /// This is the migration path for callers that previously used
    /// [`Self::new`]: the caller must still present the observed descriptor and
    /// exact matched-cell count, so this convenience does not infer or mint
    /// refusal evidence from a row identifier.
    #[deprecated(
        since = "0.1.0",
        note = "use new_with_observed_detail_cells; descriptor-only input cannot prove cell kind, order, or first typed divergence"
    )]
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_detail_descriptor(
        journey: BaseE2eJourneyV1,
        row_id: StableTokenV2,
        semantic_manifest_root: ContentHash,
        observed: BaseE2eExpectedDecisionV1,
        counts: BaseE2eObservedCountsV1,
        observed_detail_manifest: BaseE2eDetailManifestV1,
        detail_cells_matched: u32,
        first_unexpected_cell: Option<StableTokenV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new_with_detail_manifest(
            journey,
            row_id,
            semantic_manifest_root,
            observed,
            counts,
            observed_detail_manifest.root(),
            observed_detail_manifest.cell_count(),
            detail_cells_matched,
            first_unexpected_cell,
        )
    }

    /// Construct one checked observation from an exact row and the complete
    /// ordered caller-observed detail-cell slice.
    ///
    /// The supplied descriptor must be the exact descriptor reconstructed
    /// from `observed_detail_cells`. Every cell must name the row's case kind,
    /// use a one-based ordinal within that case's semantic matrix, have a
    /// unique stable ID and root, and appear in strictly increasing ordinal
    /// order. The exact matched-cell count and first typed observed divergence
    /// are derived here against the row's independent expected cells; callers
    /// cannot present either value directly. When the typed slices diverge,
    /// `first_unexpected_cell` is mandatory and must exactly equal the first
    /// observed cell ID, or the first missing expected cell ID when the
    /// observed slice ends early. When the typed slices are exact but the
    /// semantic partitions are red, the only admitted ID is the closed
    /// [`BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1`] sentinel.
    ///
    /// This remains a comparison-only, non-authoritative observation. It
    /// proves only that the typed input is internally closed. An exact public
    /// comparison may report equality for copied expected cells and counts,
    /// but it has no path to an execution report, witness, or execution root.
    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        reason = "the public typed constructor performs one indivisible fail-closed validation of every observed detail-cell invariant"
    )]
    pub fn new_with_observed_detail_cells(
        row: &BaseE2eProjectionRowV1,
        observed: BaseE2eExpectedDecisionV1,
        counts: BaseE2eObservedCountsV1,
        observed_detail_manifest: BaseE2eDetailManifestV1,
        observed_detail_cells: &[BaseE2eDetailCellV1],
        first_unexpected_cell: Option<StableTokenV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        let observed_count = u32::try_from(observed_detail_cells.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_result.observed_detail_cells",
                "a u32-bounded detail-cell slice",
                observed_detail_cells.len(),
            )
        })?;
        if observed_detail_manifest.cell_count() > observed_count {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_result.observed_detail_cells",
                "every cell declared by the observed detail manifest",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(observed_count),
                    u64::from(observed_detail_manifest.cell_count()),
                ),
            ));
        }
        if observed_detail_manifest.cell_count() < observed_count {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "base_e2e_result.observed_detail_cells",
                "no cell beyond the observed detail manifest count",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(observed_count),
                    u64::from(observed_detail_manifest.cell_count()),
                ),
            ));
        }
        let reconstructed =
            BaseE2eDecisionDetailManifestV1::from_cells(row.kind(), observed_detail_cells)?;
        if reconstructed.root() != observed_detail_manifest.root() {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_result.observed_detail_manifest_root",
                "the exact root reconstructed from the ordered observed detail cells",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let expected_cells = row.expected_detail_cells();
        let detail_cells_matched = u32::try_from(
            expected_cells
                .iter()
                .zip(observed_detail_cells)
                .filter(|(expected, observed)| expected == observed)
                .count(),
        )
        .expect("a case-bounded matched detail-cell count fits u32");
        let first_divergent_index = (0..expected_cells.len().max(observed_detail_cells.len()))
            .find(|&index| expected_cells.get(index) != observed_detail_cells.get(index));
        if let Some(index) = first_divergent_index {
            let required_id = observed_detail_cells
                .get(index)
                .or_else(|| expected_cells.get(index))
                .expect("a divergent index names an expected or observed cell")
                .stable_id();
            let Some(presented_id) = first_unexpected_cell.as_ref() else {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_e2e_result.first_unexpected_cell",
                    "the exact first typed detail-divergence cell ID",
                    ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
                ));
            };
            if presented_id.as_str() != required_id {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_result.first_unexpected_cell",
                    "the exact first typed detail-divergence cell ID",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        } else if counts.unexpected_mismatches() > 0 {
            let Some(presented_id) = first_unexpected_cell.as_ref() else {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_e2e_result.first_unexpected_cell",
                    "the closed row.contract divergence ID when typed detail cells are exact",
                    ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
                ));
            };
            if presented_id.as_str() != BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1 {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_result.first_unexpected_cell",
                    "the closed row.contract divergence ID when typed detail cells are exact",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        let first_observed_detail_cell = first_divergent_index
            .and_then(|index| observed_detail_cells.get(index))
            .cloned();

        Self::new_with_detail_manifest_and_first_observed_cell(
            row.journey(),
            row.id().clone(),
            row.semantic_manifest_root(),
            observed,
            counts,
            reconstructed.root(),
            reconstructed.cell_count(),
            detail_cells_matched,
            first_unexpected_cell,
            first_observed_detail_cell,
            true,
        )
    }

    /// Journey claimed by the observation.
    #[must_use]
    pub const fn journey(&self) -> BaseE2eJourneyV1 {
        self.journey
    }

    /// Stable row ID claimed by the observation.
    #[must_use]
    pub const fn row_id(&self) -> &StableTokenV2 {
        &self.row_id
    }

    /// Semantic manifest root claimed by the observation.
    #[must_use]
    pub const fn semantic_manifest_root(&self) -> ContentHash {
        self.semantic_manifest_root
    }

    /// Observed aggregate row decision.
    #[must_use]
    pub const fn observed(&self) -> BaseE2eExpectedDecisionV1 {
        self.observed
    }

    /// Exact observed semantic-cell partitions.
    #[must_use]
    pub const fn counts(&self) -> BaseE2eObservedCountsV1 {
        self.counts
    }

    /// Whether the caller supplied the complete ordered typed detail-cell
    /// slice and the constructor verified it against the presented
    /// descriptor.
    ///
    /// Deprecated root/count-only compatibility constructors always return
    /// `false`; an opaque descriptor alone cannot establish membership or
    /// order and therefore can never produce a green exact join.
    #[must_use]
    pub const fn typed_detail_cells_presented(&self) -> bool {
        self.typed_detail_cells_presented
    }

    /// Legacy compatibility accessor.
    ///
    /// Caller-presented rows are now comparison-only by type, so this always
    /// returns `false`. Execution witnesses exist only on checked rows returned
    /// by [`run_base_e2e_journey_v1`].
    #[deprecated(
        since = "0.1.0",
        note = "caller-presented rows are comparison-only; inspect execution_witness_root on executed checked rows"
    )]
    #[must_use]
    pub const fn has_internal_execution_observation(&self) -> bool {
        false
    }

    /// Ordered observed refusal/unsupported detail-manifest root.
    #[must_use]
    pub const fn observed_detail_manifest_root(&self) -> ContentHash {
        self.observed_detail_manifest_root
    }

    /// Closed observed detail-manifest descriptor.
    #[must_use]
    pub const fn observed_detail_manifest(&self) -> BaseE2eDecisionDetailManifestV1 {
        BaseE2eDecisionDetailManifestV1 {
            cell_count: self.observed_detail_cell_count,
            root: self.observed_detail_manifest_root,
        }
    }

    /// Observed detail cells in the ordered manifest.
    #[must_use]
    pub const fn observed_detail_cell_count(&self) -> u32 {
        self.observed_detail_cell_count
    }

    /// Observed detail cells exactly equal to their independent oracle.
    #[must_use]
    pub const fn detail_cells_matched(&self) -> u32 {
        self.detail_cells_matched
    }

    /// First unexpected cell, if any.
    #[must_use]
    pub const fn first_unexpected_cell(&self) -> Option<&StableTokenV2> {
        self.first_unexpected_cell.as_ref()
    }

    /// Bounded typed first observed detail cell supplied by an in-process
    /// execution, when one was available.
    #[must_use]
    pub const fn first_observed_detail_cell(&self) -> Option<&BaseE2eDetailCellV1> {
        self.first_observed_detail_cell.as_ref()
    }

    /// Domain-separated presented-observation root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Checked result for one public semantic projection row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eRowResultV1 {
    kind: BaseE2eCaseKindV1,
    row_id: StableTokenV2,
    expected: BaseE2eExpectedDecisionV1,
    observed: BaseE2eExpectedDecisionV1,
    checked_cells: u32,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    expected_detail_manifest_root: ContentHash,
    observed_detail_manifest_root: ContentHash,
    expected_detail_cell_count: u32,
    observed_detail_cell_count: u32,
    detail_cells_matched: u32,
    observed_detail_cells_verified: bool,
    execution_witness_root: Option<ContentHash>,
    unexpected_mismatches: u32,
    matched: bool,
    first_unexpected_cell: Option<String>,
    first_observed_detail_divergence: Option<BaseE2eDetailDivergenceV1>,
    first_divergence_root: Option<ContentHash>,
    root: ContentHash,
}

impl BaseE2eRowResultV1 {
    /// Semantic case kind that owns this checked result.
    #[must_use]
    pub const fn kind(&self) -> BaseE2eCaseKindV1 {
        self.kind
    }

    /// Stable manifest row ID.
    #[must_use]
    pub const fn row_id(&self) -> &StableTokenV2 {
        &self.row_id
    }

    /// Expected three-way pure decision.
    #[must_use]
    pub const fn expected(&self) -> BaseE2eExpectedDecisionV1 {
        self.expected
    }

    /// Observed three-way pure decision.
    #[must_use]
    pub const fn observed(&self) -> BaseE2eExpectedDecisionV1 {
        self.observed
    }

    /// Exact semantic cells checked for this row.
    #[must_use]
    pub const fn checked_cells(&self) -> u32 {
        self.checked_cells
    }

    /// Semantic cells independently expected to accept.
    #[must_use]
    pub const fn positive_eligible(&self) -> u32 {
        self.positive_eligible
    }

    /// Expected-accept semantic cells whose observations matched.
    #[must_use]
    pub const fn positive_matched(&self) -> u32 {
        self.positive_matched
    }

    /// Semantic cells independently expected to refuse.
    #[must_use]
    pub const fn expected_refusals(&self) -> u32 {
        self.expected_refusals
    }

    /// Expected-refusal semantic cells whose observations matched.
    #[must_use]
    pub const fn expected_refusals_matched(&self) -> u32 {
        self.expected_refusals_matched
    }

    /// Exactly matched, explicitly unsupported semantic cells.
    #[must_use]
    pub const fn unsupported(&self) -> u32 {
        self.unsupported
    }

    /// Ordered independent expected decision-detail manifest.
    #[must_use]
    pub const fn expected_detail_manifest_root(&self) -> ContentHash {
        self.expected_detail_manifest_root
    }

    /// Closed expected detail-manifest descriptor.
    #[must_use]
    pub const fn expected_detail_manifest(&self) -> BaseE2eDecisionDetailManifestV1 {
        BaseE2eDecisionDetailManifestV1 {
            cell_count: self.expected_detail_cell_count,
            root: self.expected_detail_manifest_root,
        }
    }

    /// Ordered observed decision-detail manifest.
    #[must_use]
    pub const fn observed_detail_manifest_root(&self) -> ContentHash {
        self.observed_detail_manifest_root
    }

    /// Closed observed detail-manifest descriptor.
    #[must_use]
    pub const fn observed_detail_manifest(&self) -> BaseE2eDecisionDetailManifestV1 {
        BaseE2eDecisionDetailManifestV1 {
            cell_count: self.observed_detail_cell_count,
            root: self.observed_detail_manifest_root,
        }
    }

    /// Expected refusal plus unsupported detail cells.
    #[must_use]
    pub const fn expected_detail_cell_count(&self) -> u32 {
        self.expected_detail_cell_count
    }

    /// Observed detail cells.
    #[must_use]
    pub const fn observed_detail_cell_count(&self) -> u32 {
        self.observed_detail_cell_count
    }

    /// Detail cells that exactly matched their independent oracle.
    #[must_use]
    pub const fn detail_cells_matched(&self) -> u32 {
        self.detail_cells_matched
    }

    /// Complete typed observed detail cells when exact equality was verified.
    ///
    /// No duplicate cell vector is retained. Exact observations borrow the
    /// immutable cached oracle only after typed equality was established.
    /// Descriptor-only and mismatched typed observations return `None` and
    /// expose only
    /// [`Self::first_observed_detail_divergence`].
    #[must_use]
    pub fn observed_detail_cells(&self) -> Option<&'static [BaseE2eDetailCellV1]> {
        self.observed_detail_cells_verified
            .then(|| expected_detail_cells(self.kind))
    }

    /// Exact private in-process execution witness, when this checked row came
    /// from [`run_base_e2e_journey_v1`].
    ///
    /// Rows in [`BaseE2eJourneyComparisonReportV1`] always return `None`.
    #[must_use]
    pub const fn execution_witness_root(&self) -> Option<ContentHash> {
        self.execution_witness_root
    }

    /// Legacy boolean view of [`Self::execution_witness_root`].
    #[deprecated(
        since = "0.1.0",
        note = "use execution_witness_root to retain and inspect the exact witness identity"
    )]
    #[must_use]
    pub const fn execution_observation_verified(&self) -> bool {
        self.execution_witness_root.is_some()
    }

    /// First bounded observed-detail divergence, when the complete observed
    /// manifest cannot be represented by the cached expected slice.
    #[must_use]
    pub const fn first_observed_detail_divergence(&self) -> Option<&BaseE2eDetailDivergenceV1> {
        self.first_observed_detail_divergence.as_ref()
    }

    /// Typed bounded root of the first failed detail cell or row-contract
    /// divergence.
    #[must_use]
    pub const fn first_divergence_root(&self) -> Option<ContentHash> {
        self.first_divergence_root
    }

    /// Unexpected semantic-cell mismatches.
    #[must_use]
    pub const fn unexpected_mismatches(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Whether expected and observed decisions and subcase counts agreed.
    #[must_use]
    pub const fn matched(&self) -> bool {
        self.matched
    }

    /// First unexpected semantic cell, if the row disagreed.
    #[must_use]
    pub fn first_unexpected_cell(&self) -> Option<&str> {
        self.first_unexpected_cell.as_deref()
    }

    /// Domain-separated checked row-result root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// One journey's public, comparison-only exact manifest join.
///
/// This report can state whether caller-presented values equal the immutable
/// oracle. It carries no execution witness, harness binding, execution report,
/// or execution-root API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eJourneyComparisonReportV1 {
    journey: BaseE2eJourneyV1,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    checked_cells: u32,
    results: Box<[BaseE2eRowResultV1]>,
    projection_root: ContentHash,
    log_schema_root: ContentHash,
    root: ContentHash,
}

impl BaseE2eJourneyComparisonReportV1 {
    /// Journey whose immutable manifest was used as the comparison oracle.
    #[must_use]
    pub const fn journey(&self) -> BaseE2eJourneyV1 {
        self.journey
    }

    /// Rows expected to accept.
    #[must_use]
    pub const fn positive_eligible(&self) -> u32 {
        self.positive_eligible
    }

    /// Expected-accept rows whose caller-presented values matched.
    #[must_use]
    pub const fn positive_matched(&self) -> u32 {
        self.positive_matched
    }

    /// Rows expected to refuse.
    #[must_use]
    pub const fn expected_refusals(&self) -> u32 {
        self.expected_refusals
    }

    /// Expected-refusal rows whose caller-presented values matched.
    #[must_use]
    pub const fn expected_refusals_matched(&self) -> u32 {
        self.expected_refusals_matched
    }

    /// Exactly matched typed unsupported rows.
    #[must_use]
    pub const fn unsupported(&self) -> u32 {
        self.unsupported
    }

    /// Any caller-presented decision, partition, or detail disagreement.
    #[must_use]
    pub const fn unexpected_mismatches(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Total bounded semantic cells represented by the comparison.
    #[must_use]
    pub const fn checked_cells(&self) -> u32 {
        self.checked_cells
    }

    /// Exact comparison rows in manifest order.
    #[must_use]
    pub fn results(&self) -> &[BaseE2eRowResultV1] {
        &self.results
    }

    /// Whether every comparison partition exactly matched.
    #[must_use]
    pub const fn exact_match(&self) -> bool {
        self.positive_matched == self.positive_eligible
            && self.expected_refusals_matched == self.expected_refusals
            && self.unexpected_mismatches == 0
    }

    /// Immutable journey-manifest root used as the comparison oracle.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Compatibility alias for [`Self::manifest_root`].
    #[deprecated(note = "use manifest_root")]
    #[must_use]
    pub const fn projection_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Closed logging schema associated with the immutable journey manifest.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    /// Domain-separated comparison-only root.
    #[must_use]
    pub const fn comparison_root(&self) -> ContentHash {
        self.root
    }

    /// Compatibility alias for [`Self::comparison_root`].
    #[deprecated(note = "use comparison_root; this report has no execution root")]
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// One journey's exact in-process execution under a presented harness context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eJourneyExecutionReportV1 {
    journey: BaseE2eJourneyV1,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    checked_cells: u32,
    results: Box<[BaseE2eRowResultV1]>,
    projection_root: ContentHash,
    harness_context_root: ContentHash,
    log_schema_root: ContentHash,
    root: ContentHash,
}

impl BaseE2eJourneyExecutionReportV1 {
    /// Journey whose exact manifest was executed.
    #[must_use]
    pub const fn journey(&self) -> BaseE2eJourneyV1 {
        self.journey
    }

    /// Rows expected to accept.
    #[must_use]
    pub const fn positive_eligible(&self) -> u32 {
        self.positive_eligible
    }

    /// Expected-accept rows that matched.
    #[must_use]
    pub const fn positive_matched(&self) -> u32 {
        self.positive_matched
    }

    /// Rows expected to refuse.
    #[must_use]
    pub const fn expected_refusals(&self) -> u32 {
        self.expected_refusals
    }

    /// Expected-refusal rows that matched.
    #[must_use]
    pub const fn expected_refusals_matched(&self) -> u32 {
        self.expected_refusals_matched
    }

    /// Exactly matched typed unsupported rows.
    #[must_use]
    pub const fn unsupported(&self) -> u32 {
        self.unsupported
    }

    /// Any row decision or semantic-cell-count disagreement.
    #[must_use]
    pub const fn unexpected_mismatches(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Total bounded semantic cells actually checked.
    #[must_use]
    pub const fn checked_cells(&self) -> u32 {
        self.checked_cells
    }

    /// Exact row results in manifest order.
    #[must_use]
    pub fn results(&self) -> &[BaseE2eRowResultV1] {
        &self.results
    }

    /// Immutable journey projection root that was consumed.
    #[deprecated(note = "use manifest_root to distinguish manifests from execution results")]
    #[must_use]
    pub const fn projection_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Explicitly named immutable journey-manifest root that was consumed.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Presented harness-context root under which execution occurred.
    #[must_use]
    pub const fn harness_context_root(&self) -> ContentHash {
        self.harness_context_root
    }

    /// Closed logging schema under which the journey is reportable.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    /// Domain-separated checked journey-execution root.
    #[deprecated(note = "use execution_root to distinguish execution from manifest identity")]
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Explicitly named context-bound journey-execution root.
    #[must_use]
    pub const fn execution_root(&self) -> ContentHash {
        self.root
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eCheckedObservationV1 {
    observed: BaseE2eExpectedDecisionV1,
    counts: BaseE2eObservedCountsV1,
    typed_detail_cells_presented: bool,
    observed_detail_manifest: BaseE2eDecisionDetailManifestV1,
    detail_cells_matched: u32,
    first_unexpected_cell: Option<String>,
    first_observed_detail_cell: Option<BaseE2eDetailCellV1>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eExecutedRowV1 {
    row_ordinal: u32,
    observation: BaseE2eCheckedObservationV1,
    comparison_root: ContentHash,
    witness_root: ContentHash,
    result: BaseE2eRowResultV1,
}

fn observation_from_presented(
    presented: &BaseE2ePresentedRowResultV1,
) -> BaseE2eCheckedObservationV1 {
    BaseE2eCheckedObservationV1 {
        observed: presented.observed(),
        counts: presented.counts(),
        typed_detail_cells_presented: presented.typed_detail_cells_presented(),
        observed_detail_manifest: presented.observed_detail_manifest(),
        detail_cells_matched: presented.detail_cells_matched(),
        first_unexpected_cell: presented
            .first_unexpected_cell()
            .map(|cell| cell.as_str().to_owned()),
        first_observed_detail_cell: presented.first_observed_detail_cell().cloned(),
    }
}

#[allow(
    clippy::too_many_lines,
    reason = "the private finalizer derives and reconciles every redundant execution field in one fail-closed transaction"
)]
fn observation_from_execution(
    row: &BaseE2eProjectionRowV1,
    execution: BaseE2eCaseExecutionV1,
) -> Result<BaseE2eCheckedObservationV1, ConstructionErrorV2> {
    let expected_cells = execution.detail.expected_cells.as_deref().ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_execution.expected_detail_cells",
            "the complete in-process expected detail-cell slice",
            ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
        )
    })?;
    if expected_cells != row.expected_detail_cells()
        || execution.detail.expected != row.expected_detail_manifest()
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.expected_detail_manifest",
            "the exact independent row detail oracle",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let observed_cells = execution.detail.observed_cells.as_deref().ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_execution.observed_detail_cells",
            "the complete in-process observed detail-cell slice",
            ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
        )
    })?;
    let observed_detail_manifest =
        BaseE2eDecisionDetailManifestV1::from_cells(row.kind(), observed_cells)?;
    if observed_detail_manifest != execution.detail.observed {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.observed_detail_manifest",
            "the exact descriptor reconstructed from in-process observed cells",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let detail_cells_matched = u32::try_from(
        row.expected_detail_cells()
            .iter()
            .zip(observed_cells)
            .filter(|(expected, observed)| expected == observed)
            .count(),
    )
    .expect("a case-bounded matched detail-cell count fits u32");
    let first_divergent_index = (0..row.expected_detail_cells().len().max(observed_cells.len()))
        .find(|&index| row.expected_detail_cells().get(index) != observed_cells.get(index));
    let first_divergent_cell = first_divergent_index.map(|index| {
        observed_cells
            .get(index)
            .or_else(|| row.expected_detail_cells().get(index))
            .expect("a divergent index names an expected or observed cell")
            .stable_id()
            .to_owned()
    });
    if execution.detail.matched_cells != detail_cells_matched
        || execution.detail.first_divergent_cell != first_divergent_cell
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.detail_reconciliation",
            "the exact matched count and first divergence derived from in-process cells",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let counts = BaseE2eObservedCountsV1::new(
        BaseE2eMatchedPartitionV1::new(execution.positive_eligible, execution.positive_matched)?,
        BaseE2eMatchedPartitionV1::new(
            execution.expected_refusals,
            execution.expected_refusals_matched,
        )?,
        execution.unsupported,
        execution.unexpected_mismatches,
    )?;
    if counts.checked_cells() != execution.checked_cells {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.checked_cells",
            "the exact sum of the in-process observed partitions",
            ConstructionObservedV2::unsigned_pair(
                u64::from(execution.checked_cells),
                u64::from(counts.checked_cells()),
            ),
        ));
    }
    if execution.first_failed_cell.is_some() != (counts.unexpected_mismatches() > 0) {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.first_failed_cell",
            "presence exactly when the in-process mismatch count is nonzero",
            execution.first_failed_cell.is_some(),
        ));
    }
    if let Some(detail_divergence) = first_divergent_cell.as_deref()
        && execution.first_failed_cell.as_deref() != Some(detail_divergence)
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution.first_failed_cell",
            "the exact first typed detail divergence",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let first_observed_detail_cell = first_divergent_index
        .and_then(|index| observed_cells.get(index))
        .cloned();

    Ok(BaseE2eCheckedObservationV1 {
        observed: execution.decision,
        counts,
        typed_detail_cells_presented: true,
        observed_detail_manifest,
        detail_cells_matched,
        first_unexpected_cell: execution.first_failed_cell,
        first_observed_detail_cell,
    })
}

#[allow(
    clippy::too_many_lines,
    reason = "the checked-row constructor derives every comparison partition, divergence, and root field together"
)]
fn checked_row_result(
    row: &BaseE2eProjectionRowV1,
    observation: &BaseE2eCheckedObservationV1,
    execution_witness_root: Option<ContentHash>,
) -> Result<BaseE2eRowResultV1, ConstructionErrorV2> {
    let counts = observation.counts;
    let positive = counts.positive();
    let refusals = counts.expected_refusals();
    let observed_detail_cells_verified = observation.typed_detail_cells_presented
        && observation.observed_detail_manifest.root() == row.expected_detail_manifest_root()
        && observation.observed_detail_manifest.cell_count() == row.expected_detail_cell_count()
        && observation.detail_cells_matched == row.expected_detail_cell_count();
    let matched = observation.observed == row.expected()
        && counts.checked_cells() == row.semantic_cell_count()
        && positive.eligible() == row.positive_cell_count()
        && refusals.eligible() == row.expected_refusal_cell_count()
        && counts.unsupported() == row.unsupported_cell_count()
        && positive.matched() == positive.eligible()
        && refusals.matched() == refusals.eligible()
        && observed_detail_cells_verified
        && counts.unexpected_mismatches() == 0;
    let unexpected_mismatches = if !matched && counts.unexpected_mismatches() == 0 {
        1
    } else {
        counts.unexpected_mismatches()
    };
    let first_unexpected_cell = observation.first_unexpected_cell.clone().or_else(|| {
        (!matched).then(|| {
            if observed_detail_cells_verified {
                BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1.to_owned()
            } else {
                "detail.manifest".to_owned()
            }
        })
    });
    let first_observed_detail_divergence = (!observed_detail_cells_verified).then(|| {
        detail_divergence(
            row,
            observation.observed_detail_manifest.root(),
            observation.observed_detail_manifest.cell_count(),
            observation.detail_cells_matched,
            first_unexpected_cell.as_deref(),
            observation.first_observed_detail_cell.as_ref(),
        )
    });
    let first_divergence_root = if matched {
        None
    } else if let Some(divergence) = first_observed_detail_divergence.as_ref() {
        Some(divergence.root())
    } else {
        Some(row_contract_divergence_root(
            row,
            observation,
            first_unexpected_cell
                .as_deref()
                .unwrap_or(BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1),
        )?)
    };
    let execution = BaseE2eCaseExecutionV1 {
        decision: observation.observed,
        checked_cells: counts.checked_cells(),
        positive_eligible: positive.eligible(),
        positive_matched: positive.matched(),
        expected_refusals: refusals.eligible(),
        expected_refusals_matched: refusals.matched(),
        unsupported: counts.unsupported(),
        unexpected_mismatches,
        first_failed_cell: first_unexpected_cell.clone(),
        detail: BaseE2eDetailExecutionV1 {
            expected: row.expected_detail_manifest(),
            observed: observation.observed_detail_manifest,
            expected_cells: None,
            observed_cells: None,
            matched_cells: observation.detail_cells_matched,
            first_divergent_cell: (!observed_detail_cells_verified)
                .then(|| "detail.manifest".to_owned()),
        },
    };
    let root = row_result_root(
        row,
        &execution,
        matched,
        observed_detail_cells_verified,
        execution_witness_root,
        first_observed_detail_divergence.as_ref(),
        first_divergence_root,
    )?;
    Ok(BaseE2eRowResultV1 {
        kind: row.kind(),
        row_id: row.id().clone(),
        expected: row.expected(),
        observed: observation.observed,
        checked_cells: counts.checked_cells(),
        positive_eligible: positive.eligible(),
        positive_matched: positive.matched(),
        expected_refusals: refusals.eligible(),
        expected_refusals_matched: refusals.matched(),
        unsupported: counts.unsupported(),
        expected_detail_manifest_root: row.expected_detail_manifest_root(),
        observed_detail_manifest_root: observation.observed_detail_manifest.root(),
        expected_detail_cell_count: row.expected_detail_cell_count(),
        observed_detail_cell_count: observation.observed_detail_manifest.cell_count(),
        detail_cells_matched: observation.detail_cells_matched,
        observed_detail_cells_verified,
        execution_witness_root,
        unexpected_mismatches,
        matched,
        first_unexpected_cell,
        first_observed_detail_divergence,
        first_divergence_root,
        root,
    })
}

fn in_process_row_execution_witness_root(
    manifest: &BaseE2eJourneyProjectionV1,
    row_ordinal: u32,
    row: &BaseE2eProjectionRowV1,
    harness: &BaseE2eHarnessIdentityV1,
    observation: &BaseE2eCheckedObservationV1,
    comparison: &BaseE2eRowResultV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    if comparison.execution_witness_root().is_some() {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution_witness.comparison_row",
            "a comparison row with explicit witness absence",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let counts = observation.counts;
    let mut frame = CanonicalFrameV1::new(b"FSBASEINPROCESSROWWITNESS\x01", 4096)?;
    frame.push_str(
        "witness.execution_class",
        "source-closed-in-process-case-execution",
    )?;
    frame.push_bytes(
        "witness.journey_manifest_root",
        manifest.manifest_root().as_bytes(),
    )?;
    frame.push_u16("witness.journey", manifest.journey().code())?;
    frame.push_u32("witness.row_ordinal", row_ordinal)?;
    frame.push_str("witness.row_id", row.id().as_str())?;
    frame.push_u16("witness.row_kind", row.kind().code())?;
    frame.push_bytes(
        "witness.semantic_manifest_root",
        row.semantic_manifest_root().as_bytes(),
    )?;
    frame.push_bytes("witness.mapping_root", row.mapping_root().as_bytes())?;
    frame.push_bytes(
        "witness.harness_context_root",
        harness.context_root().as_bytes(),
    )?;
    frame.push_u16("witness.observed", observation.observed.code())?;
    frame.push_u32("witness.positive_eligible", counts.positive().eligible())?;
    frame.push_u32("witness.positive_matched", counts.positive().matched())?;
    frame.push_u32(
        "witness.expected_refusals",
        counts.expected_refusals().eligible(),
    )?;
    frame.push_u32(
        "witness.expected_refusals_matched",
        counts.expected_refusals().matched(),
    )?;
    frame.push_u32("witness.unsupported", counts.unsupported())?;
    frame.push_u32(
        "witness.unexpected_mismatches",
        counts.unexpected_mismatches(),
    )?;
    frame.push_u32("witness.checked_cells", counts.checked_cells())?;
    frame.push_bytes(
        "witness.observed_detail_manifest_root",
        observation.observed_detail_manifest.root().as_bytes(),
    )?;
    frame.push_u32(
        "witness.observed_detail_cell_count",
        observation.observed_detail_manifest.cell_count(),
    )?;
    frame.push_u32(
        "witness.detail_cells_matched",
        observation.detail_cells_matched,
    )?;
    frame.push_u16(
        "witness.first_unexpected_presence",
        u16::from(comparison.first_unexpected_cell().is_some()),
    )?;
    if let Some(first_unexpected_cell) = comparison.first_unexpected_cell() {
        frame.push_str("witness.first_unexpected_cell", first_unexpected_cell)?;
    }
    frame.push_u16(
        "witness.divergence_presence",
        u16::from(comparison.first_divergence_root().is_some()),
    )?;
    if let Some(first_divergence_root) = comparison.first_divergence_root() {
        frame.push_bytes(
            "witness.first_divergence_root",
            first_divergence_root.as_bytes(),
        )?;
    }
    frame.push_bytes("witness.comparison_row_root", comparison.root().as_bytes())?;
    Ok(frame.root(BASE_E2E_IN_PROCESS_ROW_EXECUTION_WITNESS_DOMAIN_V1))
}

fn finalize_executed_row(
    manifest: &BaseE2eJourneyProjectionV1,
    row_ordinal: u32,
    row: &BaseE2eProjectionRowV1,
    harness: &BaseE2eHarnessIdentityV1,
    execution: BaseE2eCaseExecutionV1,
) -> Result<BaseE2eExecutedRowV1, ConstructionErrorV2> {
    let row_index = row_ordinal.checked_sub(1).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "base_e2e_execution_witness.row_ordinal",
            "a one-based row ordinal",
            row_ordinal,
        )
    })?;
    let expected_row = manifest
        .rows()
        .get(usize::try_from(row_index).map_err(|_| sequence_overflow())?)
        .ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_e2e_execution_witness.row_ordinal",
                "a one-based ordinal within the exact journey manifest",
                row_ordinal,
            )
        })?;
    if expected_row != row {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_execution_witness.row_binding",
            "the exact row at the presented journey ordinal",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let observation = observation_from_execution(row, execution)?;
    let comparison = checked_row_result(row, &observation, None)?;
    let comparison_root = comparison.root();
    let witness_root = in_process_row_execution_witness_root(
        manifest,
        row_ordinal,
        row,
        harness,
        &observation,
        &comparison,
    )?;
    let result = checked_row_result(row, &observation, Some(witness_root))?;
    Ok(BaseE2eExecutedRowV1 {
        row_ordinal,
        observation,
        comparison_root,
        witness_root,
        result,
    })
}

fn validate_executed_row(
    manifest: &BaseE2eJourneyProjectionV1,
    expected_ordinal: u32,
    row: &BaseE2eProjectionRowV1,
    harness: &BaseE2eHarnessIdentityV1,
    executed: &BaseE2eExecutedRowV1,
) -> Result<(), ConstructionErrorV2> {
    if executed.row_ordinal != expected_ordinal {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfOrder,
            "base_e2e_journey_execution.row_ordinal",
            "the exact one-based journey row ordinal",
            executed.row_ordinal,
        ));
    }
    let comparison = checked_row_result(row, &executed.observation, None)?;
    if comparison.root() != executed.comparison_root {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_journey_execution.comparison_root",
            "the independently reconstructed comparison row root",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let witness_root = in_process_row_execution_witness_root(
        manifest,
        expected_ordinal,
        row,
        harness,
        &executed.observation,
        &comparison,
    )?;
    if witness_root != executed.witness_root {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_journey_execution.witness_root",
            "the independently reconstructed in-process row witness",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let result = checked_row_result(row, &executed.observation, Some(witness_root))?;
    if result != executed.result {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_journey_execution.row_result",
            "the exact witness-bound checked row result",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    Ok(())
}

fn finalize_journey_execution(
    manifest: &BaseE2eJourneyProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
    executed_rows: Vec<BaseE2eExecutedRowV1>,
) -> Result<BaseE2eJourneyExecutionReportV1, ConstructionErrorV2> {
    if executed_rows.len() < manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_journey_execution.rows",
            "one ordered in-process witness for every journey row",
            executed_rows.len(),
        ));
    }
    if executed_rows.len() > manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_journey_execution.rows",
            "no witness beyond the exact journey row set",
            executed_rows.len(),
        ));
    }

    let mut positive_eligible = 0_u32;
    let mut positive_matched = 0_u32;
    let mut expected_refusals = 0_u32;
    let mut expected_refusals_matched = 0_u32;
    let mut unsupported = 0_u32;
    let mut unexpected_mismatches = 0_u32;
    let mut checked_cells = 0_u32;
    let mut witness_roots = std::collections::BTreeSet::new();
    let mut results = Vec::with_capacity(executed_rows.len());
    for (index, (row, executed)) in manifest.rows().iter().zip(executed_rows).enumerate() {
        let row_ordinal = u32::try_from(index + 1).map_err(|_| sequence_overflow())?;
        validate_executed_row(manifest, row_ordinal, row, harness, &executed)?;
        if !witness_roots.insert(executed.witness_root) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_journey_execution.witness_root",
                "one distinct in-process witness per ordered journey row",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let result = executed.result;
        positive_eligible = positive_eligible
            .checked_add(result.positive_eligible())
            .ok_or_else(sequence_overflow)?;
        positive_matched = positive_matched
            .checked_add(result.positive_matched())
            .ok_or_else(sequence_overflow)?;
        expected_refusals = expected_refusals
            .checked_add(result.expected_refusals())
            .ok_or_else(sequence_overflow)?;
        expected_refusals_matched = expected_refusals_matched
            .checked_add(result.expected_refusals_matched())
            .ok_or_else(sequence_overflow)?;
        unsupported = unsupported
            .checked_add(result.unsupported())
            .ok_or_else(sequence_overflow)?;
        unexpected_mismatches = unexpected_mismatches
            .checked_add(result.unexpected_mismatches())
            .ok_or_else(sequence_overflow)?;
        checked_cells = checked_cells
            .checked_add(result.checked_cells())
            .ok_or_else(sequence_overflow)?;
        results.push(result);
    }

    let root = journey_execution_root(
        manifest,
        harness,
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unsupported,
        unexpected_mismatches,
        checked_cells,
        &results,
    )?;
    Ok(BaseE2eJourneyExecutionReportV1 {
        journey: manifest.journey(),
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unsupported,
        unexpected_mismatches,
        checked_cells,
        results: results.into_boxed_slice(),
        projection_root: manifest.manifest_root(),
        harness_context_root: harness.context_root(),
        log_schema_root: manifest.log_schema_root(),
        root,
    })
}

/// Execute one exact journey-specific manifest through real public
/// constructors and validators.
///
/// The function accepts only a member of the immutable frozen projection and
/// binds the checked results to the presented source/build/toolchain/target/
/// feature context. It performs no filesystem, process, publication, or
/// authority effect.
pub fn run_base_e2e_journey_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    journey: BaseE2eJourneyV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<BaseE2eJourneyExecutionReportV1, ConstructionErrorV2> {
    let expected = RunnerV2BaseE2eProjectionV1::frozen()?;
    if projection != &expected {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_journey.projection",
            "the exact frozen projection",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    run_base_e2e_journey_rows_v1(projection, journey, harness)
}

fn run_base_e2e_journey_rows_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    journey: BaseE2eJourneyV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<BaseE2eJourneyExecutionReportV1, ConstructionErrorV2> {
    let manifest = projection
        .journeys()
        .iter()
        .find(|candidate| candidate.journey() == journey)
        .ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_journey.manifest",
                "one exact journey manifest",
                ConstructionObservedV2::closed(&journey),
            )
        })?;

    let mut executed_rows = Vec::with_capacity(manifest.rows().len());
    for (index, row) in manifest.rows().iter().enumerate() {
        // Every manifest row is an execution claim. Re-run the public semantic
        // case for every row instead of reusing an observation from another
        // row or journey.
        let execution = execute_case(row.kind(), harness);
        let row_ordinal = u32::try_from(index + 1).map_err(|_| sequence_overflow())?;
        executed_rows.push(finalize_executed_row(
            manifest,
            row_ordinal,
            row,
            harness,
            execution,
        )?);
    }
    finalize_journey_execution(manifest, harness, executed_rows)
}

/// Compare caller-presented row observations with one frozen journey.
///
/// Missing, extra, duplicate, reordered, stale-root, unmapped, and
/// cross-journey observations refuse before any summary is returned. Exact
/// equality is comparison data only and cannot produce an execution report,
/// witness, or execution root.
#[allow(
    clippy::too_many_lines,
    reason = "the exact journey comparison is one fail-closed validation transaction whose ordered checks define its refusal precedence"
)]
pub fn compare_base_e2e_journey_results_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    journey: BaseE2eJourneyV1,
    presented: &[BaseE2ePresentedRowResultV1],
) -> Result<BaseE2eJourneyComparisonReportV1, ConstructionErrorV2> {
    let expected_projection = RunnerV2BaseE2eProjectionV1::frozen()?;
    if projection != &expected_projection {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_join.projection",
            "the exact frozen projection",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let manifest = projection
        .journeys()
        .iter()
        .find(|candidate| candidate.journey() == journey)
        .ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_join.manifest",
                "one exact journey manifest",
                ConstructionObservedV2::closed(&journey),
            )
        })?;

    let mut seen_ids = std::collections::BTreeSet::new();
    let mut seen_roots = std::collections::BTreeSet::new();
    for result in presented {
        if !seen_ids.insert(result.row_id().as_str()) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_join.row_id",
                "one result per exact journey row",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        if !seen_roots.insert(result.root()) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_join.presented_root",
                "one distinct result root per exact journey row",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    if presented.len() < manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_join.results",
            "one ordered result for every exact manifest row",
            ConstructionObservedV2::unsigned_pair(
                presented.len() as u64,
                manifest.rows().len() as u64,
            ),
        ));
    }
    if presented.len() > manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_join.results",
            "no result beyond the exact manifest row set",
            ConstructionObservedV2::unsigned_pair(
                presented.len() as u64,
                manifest.rows().len() as u64,
            ),
        ));
    }

    let target_ids = manifest
        .rows()
        .iter()
        .map(|row| row.id().as_str())
        .collect::<std::collections::BTreeSet<_>>();
    for result in presented {
        if result.journey() != journey {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_join.journey",
                "the selected journey on every result",
                ConstructionObservedV2::closed(&result.journey()),
            ));
        }
        if !target_ids.contains(result.row_id().as_str()) {
            let mapped_elsewhere = projection.journeys().iter().any(|candidate| {
                candidate.journey() != journey
                    && candidate
                        .rows()
                        .iter()
                        .any(|row| row.id() == result.row_id())
            });
            return Err(ConstructionErrorV2::new_redacted(
                if mapped_elsewhere {
                    ConstructionErrorKindV2::Incompatible
                } else {
                    ConstructionErrorKindV2::Unexpected
                },
                "base_e2e_join.row_id",
                "an exact row mapped to the selected journey",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }

    let ordered_ids_match = manifest
        .rows()
        .iter()
        .zip(presented)
        .all(|(row, result)| row.id() == result.row_id());
    if !ordered_ids_match {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfOrder,
            "base_e2e_join.results",
            "exact manifest row order",
            ConstructionObservedV2::fixed(ConstructionFixedObservationV2::ExactSetDifferentOrder),
        ));
    }

    let mut positive_eligible = 0_u32;
    let mut positive_matched = 0_u32;
    let mut expected_refusals = 0_u32;
    let mut expected_refusals_matched = 0_u32;
    let mut unsupported = 0_u32;
    let mut unexpected_mismatches = 0_u32;
    let mut checked_cells = 0_u32;
    let mut results = Vec::with_capacity(manifest.rows().len());

    for (row, presented_result) in manifest.rows().iter().zip(presented) {
        if presented_result.semantic_manifest_root() != row.semantic_manifest_root() {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_join.semantic_manifest_root",
                "the exact source-closed row manifest root",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let observation = observation_from_presented(presented_result);
        let result = checked_row_result(row, &observation, None)?;
        positive_eligible = positive_eligible
            .checked_add(result.positive_eligible())
            .ok_or_else(sequence_overflow)?;
        positive_matched = positive_matched
            .checked_add(result.positive_matched())
            .ok_or_else(sequence_overflow)?;
        expected_refusals = expected_refusals
            .checked_add(result.expected_refusals())
            .ok_or_else(sequence_overflow)?;
        expected_refusals_matched = expected_refusals_matched
            .checked_add(result.expected_refusals_matched())
            .ok_or_else(sequence_overflow)?;
        unsupported = unsupported
            .checked_add(result.unsupported())
            .ok_or_else(sequence_overflow)?;
        unexpected_mismatches = unexpected_mismatches
            .checked_add(result.unexpected_mismatches())
            .ok_or_else(sequence_overflow)?;
        checked_cells = checked_cells
            .checked_add(result.checked_cells())
            .ok_or_else(sequence_overflow)?;
        results.push(result);
    }

    let root = journey_comparison_root(
        manifest,
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unsupported,
        unexpected_mismatches,
        checked_cells,
        &results,
    )?;
    Ok(BaseE2eJourneyComparisonReportV1 {
        journey,
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unsupported,
        unexpected_mismatches,
        checked_cells,
        results: results.into_boxed_slice(),
        projection_root: manifest.manifest_root(),
        log_schema_root: manifest.log_schema_root(),
        root,
    })
}

/// Deprecated compatibility name for
/// [`compare_base_e2e_journey_results_v1`].
///
/// The harness parameter is retained for source migration only. It does not
/// turn comparison data into execution evidence, and the returned type has no
/// execution-root API.
#[deprecated(
    since = "0.1.0",
    note = "use compare_base_e2e_journey_results_v1; public joins now return comparison-only reports"
)]
pub fn join_base_e2e_journey_results_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    journey: BaseE2eJourneyV1,
    _harness: &BaseE2eHarnessIdentityV1,
    presented: &[BaseE2ePresentedRowResultV1],
) -> Result<BaseE2eJourneyComparisonReportV1, ConstructionErrorV2> {
    compare_base_e2e_journey_results_v1(projection, journey, presented)
}

#[allow(
    clippy::too_many_lines,
    reason = "the canonical row-result encoder keeps every authority-bearing field in one auditable domain-separated frame"
)]
fn row_result_root(
    row: &BaseE2eProjectionRowV1,
    execution: &BaseE2eCaseExecutionV1,
    matched: bool,
    observed_detail_cells_verified: bool,
    execution_witness_root: Option<ContentHash>,
    first_observed_detail_divergence: Option<&BaseE2eDetailDivergenceV1>,
    first_divergence_root: Option<ContentHash>,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEROWRESULT\x01", 4096)?;
    frame.push_str("result.row_id", row.id().as_str())?;
    frame.push_bytes(
        "result.semantic_manifest_root",
        row.semantic_manifest_root().as_bytes(),
    )?;
    frame.push_u16("result.expected", row.expected().code())?;
    frame.push_u16("result.observed", execution.decision.code())?;
    frame.push_u32("result.expected_cell_count", row.semantic_cell_count())?;
    frame.push_u32("result.checked_cell_count", execution.checked_cells)?;
    frame.push_u32("result.expected_positive_cells", row.positive_cell_count())?;
    frame.push_u32(
        "result.observed_positive_eligible",
        execution.positive_eligible,
    )?;
    frame.push_u32(
        "result.observed_positive_matched",
        execution.positive_matched,
    )?;
    frame.push_u32(
        "result.expected_refusal_cells",
        row.expected_refusal_cell_count(),
    )?;
    frame.push_u32(
        "result.observed_expected_refusals",
        execution.expected_refusals,
    )?;
    frame.push_u32(
        "result.observed_expected_refusals_matched",
        execution.expected_refusals_matched,
    )?;
    frame.push_u32(
        "result.expected_unsupported_cells",
        row.unsupported_cell_count(),
    )?;
    frame.push_u32("result.observed_unsupported", execution.unsupported)?;
    frame.push_u32(
        "result.expected_detail_cell_count",
        row.expected_detail_cell_count(),
    )?;
    frame.push_bytes(
        "result.expected_detail_manifest_root",
        row.expected_detail_manifest_root().as_bytes(),
    )?;
    frame.push_u32(
        "result.observed_detail_cell_count",
        execution.detail.observed.cell_count,
    )?;
    frame.push_u32(
        "result.detail_cells_matched",
        execution.detail.matched_cells,
    )?;
    frame.push_bytes(
        "result.observed_detail_manifest_root",
        execution.detail.observed.root.as_bytes(),
    )?;
    frame.push_u16(
        "result.observed_detail_cells_verified",
        u16::from(observed_detail_cells_verified),
    )?;
    frame.push_u16(
        "result.execution_witness_presence",
        u16::from(execution_witness_root.is_some()),
    )?;
    if let Some(execution_witness_root) = execution_witness_root {
        frame.push_bytes(
            "result.execution_witness_root",
            execution_witness_root.as_bytes(),
        )?;
    }
    frame.push_u32(
        "result.unexpected_mismatches",
        execution.unexpected_mismatches,
    )?;
    frame.push_u16("result.matched", u16::from(matched))?;
    frame.push_u16(
        "result.first_unexpected_presence",
        u16::from(execution.first_failed_cell.is_some()),
    )?;
    if let Some(first_failed_cell) = &execution.first_failed_cell {
        frame.push_str("result.first_unexpected_cell", first_failed_cell)?;
    }
    frame.push_u16(
        "result.detail_divergence_presence",
        u16::from(first_observed_detail_divergence.is_some()),
    )?;
    if let Some(divergence) = first_observed_detail_divergence {
        frame.push_bytes(
            "result.detail_divergence_root",
            divergence.root().as_bytes(),
        )?;
    }
    frame.push_u16(
        "result.first_divergence_root_presence",
        u16::from(first_divergence_root.is_some()),
    )?;
    if let Some(first_divergence_root) = first_divergence_root {
        frame.push_bytes(
            "result.first_divergence_root",
            first_divergence_root.as_bytes(),
        )?;
    }
    Ok(frame.root(if execution_witness_root.is_some() {
        BASE_E2E_EXECUTED_ROW_RESULT_DOMAIN_V1
    } else {
        BASE_E2E_ROW_RESULT_DOMAIN_V1
    }))
}

fn row_contract_divergence_root(
    row: &BaseE2eProjectionRowV1,
    observed: &BaseE2eCheckedObservationV1,
    failed_cell: &str,
) -> Result<ContentHash, ConstructionErrorV2> {
    let counts = observed.counts;
    let mut frame = CanonicalFrameV1::new(b"FSBASEROWCONTRACTDIVERGENCE\x01", 4096)?;
    frame.push_str("divergence.failed_cell", failed_cell)?;
    frame.push_u16("divergence.kind", row.kind().code())?;
    frame.push_bytes(
        "divergence.expected_semantic_manifest_root",
        row.semantic_manifest_root().as_bytes(),
    )?;
    frame.push_bytes(
        "divergence.observed_semantic_manifest_root",
        row.semantic_manifest_root().as_bytes(),
    )?;
    frame.push_u16("divergence.expected_decision", row.expected().code())?;
    frame.push_u16("divergence.observed_decision", observed.observed.code())?;
    frame.push_u32(
        "divergence.expected_semantic_cells",
        row.semantic_cell_count(),
    )?;
    frame.push_u32("divergence.observed_semantic_cells", counts.checked_cells())?;
    frame.push_u32(
        "divergence.expected_positive_cells",
        row.positive_cell_count(),
    )?;
    frame.push_u32(
        "divergence.observed_positive_eligible",
        counts.positive().eligible(),
    )?;
    frame.push_u32(
        "divergence.observed_positive_matched",
        counts.positive().matched(),
    )?;
    frame.push_u32(
        "divergence.expected_refusal_cells",
        row.expected_refusal_cell_count(),
    )?;
    frame.push_u32(
        "divergence.observed_refusal_eligible",
        counts.expected_refusals().eligible(),
    )?;
    frame.push_u32(
        "divergence.observed_refusal_matched",
        counts.expected_refusals().matched(),
    )?;
    frame.push_bytes(
        "divergence.expected_detail_manifest_root",
        row.expected_detail_manifest_root().as_bytes(),
    )?;
    frame.push_bytes(
        "divergence.observed_detail_manifest_root",
        observed.observed_detail_manifest.root().as_bytes(),
    )?;
    frame.push_u32(
        "divergence.expected_detail_cells",
        row.expected_detail_cell_count(),
    )?;
    frame.push_u32(
        "divergence.observed_detail_cells",
        observed.observed_detail_manifest.cell_count(),
    )?;
    Ok(frame.root("org.frankensim.fs-evidence-runner.base-e2e-row-contract-divergence.v1"))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the presented-result root intentionally binds every caller-observed partition and detail-cell field"
)]
fn presented_row_result_root(
    journey: BaseE2eJourneyV1,
    row_id: &StableTokenV2,
    semantic_manifest_root: ContentHash,
    observed: BaseE2eExpectedDecisionV1,
    counts: BaseE2eObservedCountsV1,
    observed_detail_manifest_root: ContentHash,
    observed_detail_cell_count: u32,
    detail_cells_matched: u32,
    first_unexpected_cell: Option<&StableTokenV2>,
    first_observed_detail_cell: Option<&BaseE2eDetailCellV1>,
    typed_detail_cells_presented: bool,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEPRESENTEDRESULT\x01", 4096)?;
    frame.push_u16("presented.journey", journey.code())?;
    frame.push_str("presented.row_id", row_id.as_str())?;
    frame.push_bytes(
        "presented.semantic_manifest_root",
        semantic_manifest_root.as_bytes(),
    )?;
    frame.push_u16("presented.observed", observed.code())?;
    frame.push_u32("presented.positive_eligible", counts.positive().eligible())?;
    frame.push_u32("presented.positive_matched", counts.positive().matched())?;
    frame.push_u32(
        "presented.expected_refusals",
        counts.expected_refusals().eligible(),
    )?;
    frame.push_u32(
        "presented.expected_refusals_matched",
        counts.expected_refusals().matched(),
    )?;
    frame.push_u32("presented.unsupported", counts.unsupported())?;
    frame.push_u32(
        "presented.unexpected_mismatches",
        counts.unexpected_mismatches(),
    )?;
    frame.push_u32("presented.checked_cells", counts.checked_cells())?;
    frame.push_u32(
        "presented.observed_detail_cell_count",
        observed_detail_cell_count,
    )?;
    frame.push_u32("presented.detail_cells_matched", detail_cells_matched)?;
    frame.push_bytes(
        "presented.observed_detail_manifest_root",
        observed_detail_manifest_root.as_bytes(),
    )?;
    frame.push_u16(
        "presented.typed_detail_cells_presented",
        u16::from(typed_detail_cells_presented),
    )?;
    frame.push_u16(
        "presented.first_unexpected_presence",
        u16::from(first_unexpected_cell.is_some()),
    )?;
    if let Some(first_unexpected_cell) = first_unexpected_cell {
        frame.push_str(
            "presented.first_unexpected_cell",
            first_unexpected_cell.as_str(),
        )?;
    }
    frame.push_u16(
        "presented.first_observed_detail_cell_presence",
        u16::from(first_observed_detail_cell.is_some()),
    )?;
    if let Some(cell) = first_observed_detail_cell {
        frame.push_bytes(
            "presented.first_observed_detail_cell_root",
            cell.root().as_bytes(),
        )?;
    }
    Ok(frame.root(BASE_E2E_PRESENTED_ROW_RESULT_DOMAIN_V1))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the comparison root deliberately binds every exact comparison partition without execution context"
)]
fn journey_comparison_root(
    manifest: &BaseE2eJourneyProjectionV1,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    checked_cells: u32,
    results: &[BaseE2eRowResultV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    if results.len() != manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            if results.len() < manifest.rows().len() {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "base_e2e_journey_comparison.results",
            "one ordered comparison result for every manifest row",
            results.len(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSBASEJOURNEYCOMPARISON\x01", 64 * 1024)?;
    frame.push_bytes(
        "comparison.journey_manifest_root",
        manifest.manifest_root().as_bytes(),
    )?;
    frame.push_u16("comparison.journey", manifest.journey().code())?;
    frame.push_bytes(
        "comparison.source_closure_root",
        manifest.source_closure_root().as_bytes(),
    )?;
    frame.push_bytes(
        "comparison.log_schema_root",
        manifest.log_schema_root().as_bytes(),
    )?;
    frame.push_u32("comparison.positive_eligible", positive_eligible)?;
    frame.push_u32("comparison.positive_matched", positive_matched)?;
    frame.push_u32("comparison.expected_refusals", expected_refusals)?;
    frame.push_u32(
        "comparison.expected_refusals_matched",
        expected_refusals_matched,
    )?;
    frame.push_u32("comparison.unsupported", unsupported)?;
    frame.push_u32("comparison.unexpected_mismatches", unexpected_mismatches)?;
    frame.push_u32("comparison.checked_cells", checked_cells)?;
    frame.push_u32(
        "comparison.result_count",
        u32::try_from(results.len()).map_err(|_| sequence_overflow())?,
    )?;
    for (index, (row, result)) in manifest.rows().iter().zip(results).enumerate() {
        if result.row_id() != row.id() || result.kind() != row.kind() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "base_e2e_journey_comparison.row",
                "the exact journey row order",
                index,
            ));
        }
        if result.execution_witness_root().is_some() {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_journey_comparison.execution_witness",
                "explicit witness absence on every comparison row",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        frame.push_bytes("comparison.result_root", result.root().as_bytes())?;
    }
    Ok(frame.root(BASE_E2E_JOURNEY_COMPARISON_DOMAIN_V1))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the root deliberately binds every exact proof-vocabulary count"
)]
fn journey_execution_root(
    manifest: &BaseE2eJourneyProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    checked_cells: u32,
    results: &[BaseE2eRowResultV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    if results.len() != manifest.rows().len() {
        return Err(ConstructionErrorV2::new(
            if results.len() < manifest.rows().len() {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "base_e2e_journey_execution.results",
            "one witness-bound result for every ordered journey row",
            results.len(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSBASEJOURNEYEXECUTION\x01", 64 * 1024)?;
    frame.push_bytes("execution.projection_root", manifest.root().as_bytes())?;
    frame.push_bytes(
        "execution.source_closure_root",
        manifest.source_closure_root().as_bytes(),
    )?;
    frame.push_bytes(
        "execution.log_schema_root",
        manifest.log_schema_root().as_bytes(),
    )?;
    frame.push_bytes(
        "execution.harness_context_root",
        harness.context_root().as_bytes(),
    )?;
    frame.push_bytes("execution.source_root", harness.source().bytes())?;
    frame.push_bytes("execution.build_root", harness.build().bytes())?;
    frame.push_bytes("execution.toolchain_root", harness.toolchain().bytes())?;
    frame.push_bytes("execution.target_root", harness.target_root().as_bytes())?;
    frame.push_bytes(
        "execution.feature_set_root",
        harness.feature_set_root().as_bytes(),
    )?;
    frame.push_u32("execution.positive_eligible", positive_eligible)?;
    frame.push_u32("execution.positive_matched", positive_matched)?;
    frame.push_u32("execution.expected_refusals", expected_refusals)?;
    frame.push_u32(
        "execution.expected_refusals_matched",
        expected_refusals_matched,
    )?;
    frame.push_u32("execution.unsupported", unsupported)?;
    frame.push_u32("execution.unexpected_mismatches", unexpected_mismatches)?;
    frame.push_u32("execution.checked_cells", checked_cells)?;
    frame.push_u32(
        "execution.result_count",
        u32::try_from(results.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "execution.result_count",
                "a u32 result count",
                results.len(),
            )
        })?,
    )?;
    for (index, (row, result)) in manifest.rows().iter().zip(results).enumerate() {
        if result.row_id() != row.id() || result.kind() != row.kind() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "base_e2e_journey_execution.row",
                "the exact journey row order",
                index,
            ));
        }
        let witness_root = result.execution_witness_root().ok_or_else(|| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Missing,
                "base_e2e_journey_execution.witness_root",
                "one private in-process witness for every result",
                ConstructionObservedDataClassV2::CallerControlledText,
            )
        })?;
        frame.push_u32(
            "execution.row_ordinal",
            u32::try_from(index + 1).map_err(|_| sequence_overflow())?,
        )?;
        frame.push_bytes("execution.row_witness_root", witness_root.as_bytes())?;
        frame.push_bytes("execution.result_root", result.root().as_bytes())?;
    }
    Ok(frame.root(BASE_E2E_JOURNEY_EXECUTION_DOMAIN_V1))
}

fn projection_execution_root(
    projection: &RunnerV2BaseE2eProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
    journey_executions: &[BaseE2eJourneyExecutionReportV1],
    retained_artifact_claim: BaseE2eRetainedArtifactClaimV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    if journey_executions.len() < projection.journeys().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_projection_execution.journeys",
            "one execution for every ordered journey manifest",
            journey_executions.len(),
        ));
    }
    if journey_executions.len() > projection.journeys().len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_projection_execution.journeys",
            "no execution beyond the five ordered journey manifests",
            journey_executions.len(),
        ));
    }

    let mut frame = CanonicalFrameV1::new(b"FSBASEPROJECTIONEXECUTION\x01", 4096)?;
    frame.push_bytes(
        "aggregate.projection_manifest_root",
        projection.manifest_root().as_bytes(),
    )?;
    frame.push_bytes(
        "aggregate.harness_context_root",
        harness.context_root().as_bytes(),
    )?;
    frame.push_u16(
        "aggregate.retained_artifact_presence",
        match retained_artifact_claim {
            BaseE2eRetainedArtifactClaimV1::Absent => 0,
        },
    )?;
    frame.push_u32(
        "aggregate.journey_count",
        u32::try_from(journey_executions.len()).map_err(|_| sequence_overflow())?,
    )?;
    for (index, (manifest, execution)) in projection
        .journeys()
        .iter()
        .zip(journey_executions)
        .enumerate()
    {
        if execution.journey() != manifest.journey() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "base_e2e_projection_execution.journey",
                "the exact frozen journey order",
                ConstructionObservedV2::closed_and_usize(&execution.journey(), index),
            ));
        }
        if execution.manifest_root() != manifest.manifest_root()
            || execution.harness_context_root() != harness.context_root()
            || execution.log_schema_root() != manifest.log_schema_root()
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_projection_execution.binding",
                "the exact journey manifest, harness context, and log schema roots",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let reconstructed_execution_root = journey_execution_root(
            manifest,
            harness,
            execution.positive_eligible,
            execution.positive_matched,
            execution.expected_refusals,
            execution.expected_refusals_matched,
            execution.unsupported,
            execution.unexpected_mismatches,
            execution.checked_cells,
            &execution.results,
        )?;
        if execution.execution_root() != reconstructed_execution_root {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_projection_execution.journey_execution_root",
                "the independently reconstructed context-bound journey execution root",
                ConstructionObservedV2::indexed_redacted(
                    index,
                    ConstructionObservedDataClassV2::CallerControlledText,
                ),
            ));
        }
        frame.push_u16("aggregate.journey", execution.journey().code())?;
        frame.push_bytes(
            "aggregate.journey_manifest_root",
            manifest.manifest_root().as_bytes(),
        )?;
        frame.push_bytes(
            "aggregate.journey_execution_root",
            execution.execution_root().as_bytes(),
        )?;
    }
    Ok(frame.root(BASE_E2E_PROJECTION_EXECUTION_DOMAIN_V1))
}

/// Run every frozen row through real in-process public constructors and
/// validators with deterministic detailed logging.
#[allow(
    clippy::too_many_lines,
    reason = "the aggregate runner keeps its five-journey accounting, exact coverage join, and deterministic log sequence in one auditable transaction"
)]
pub fn run_base_e2e_projection_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<BaseE2eProjectionReportV1, ConstructionErrorV2> {
    let expected = RunnerV2BaseE2eProjectionV1::frozen()?;
    if projection != &expected {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.root",
            "the exact immutable frozen projection",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let source_closure_report = run_source_closure_checks(&projection.source_closure);
    let mut sequence = 0_u32;
    let mut positive_eligible = 0_u32;
    let mut positive_matched = 0_u32;
    let mut expected_refusals = 0_u32;
    let mut expected_refusals_matched = 0_u32;
    let mut unsupported = 0_u32;
    let mut unexpected_mismatches = 0_u32;
    let mut legacy_eligible_rows = 0_u32;
    let mut legacy_passed_rows = 0_u32;
    let mut legacy_failed_rows = 0_u32;
    let mut projection_rows_checked = 0_u32;
    let mut projection_e2e_checked = 0_u32;
    let mut events = Vec::new();
    let mut coverage_observations = std::collections::BTreeMap::new();
    let mut journey_executions = Vec::with_capacity(projection.journeys.len());
    for journey in &projection.journeys {
        let journey_report = run_base_e2e_journey_v1(projection, journey.journey(), harness)?;
        events.push(log_event(
            sequence,
            journey,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            harness,
            None,
            vec![field(
                "expected-row-count",
                TypedValueV2::U32(
                    u32::try_from(journey.rows().len()).map_err(|_| sequence_overflow())?,
                ),
            )?],
        )?);
        sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
        positive_eligible = positive_eligible
            .checked_add(journey_report.positive_eligible())
            .ok_or_else(sequence_overflow)?;
        positive_matched = positive_matched
            .checked_add(journey_report.positive_matched())
            .ok_or_else(sequence_overflow)?;
        expected_refusals = expected_refusals
            .checked_add(journey_report.expected_refusals())
            .ok_or_else(sequence_overflow)?;
        expected_refusals_matched = expected_refusals_matched
            .checked_add(journey_report.expected_refusals_matched())
            .ok_or_else(sequence_overflow)?;
        unsupported = unsupported
            .checked_add(journey_report.unsupported())
            .ok_or_else(sequence_overflow)?;
        unexpected_mismatches = unexpected_mismatches
            .checked_add(journey_report.unexpected_mismatches())
            .ok_or_else(sequence_overflow)?;
        projection_rows_checked = projection_rows_checked
            .checked_add(
                u32::try_from(journey_report.results().len()).map_err(|_| sequence_overflow())?,
            )
            .ok_or_else(sequence_overflow)?;
        projection_e2e_checked = projection_e2e_checked
            .checked_add(journey_report.checked_cells())
            .ok_or_else(sequence_overflow)?;

        let mut journey_legacy_eligible_rows = 0_u32;
        let mut journey_legacy_passed_rows = 0_u32;
        let mut journey_legacy_failed_rows = 0_u32;
        for (row, result) in journey.rows.iter().zip(journey_report.results()) {
            let outcome = match (row.expected(), result.matched()) {
                (BaseE2eExpectedDecisionV1::Unsupported, true) => BaseE2eOutcomeV1::Unsupported,
                (_, true) => BaseE2eOutcomeV1::Passed,
                (_, false) => BaseE2eOutcomeV1::Failed,
            };
            let coverage_id = format!(
                "projection-e2e:{}:{}",
                journey.journey().key(),
                row.id().as_str()
            );
            let coverage_outcome = match (row.expected(), result.matched()) {
                (_, false) => BaseCoveragePresentedOutcomeV1::UnexpectedMismatch,
                (BaseE2eExpectedDecisionV1::Accept, true) => {
                    BaseCoveragePresentedOutcomeV1::PositiveMatched
                }
                (BaseE2eExpectedDecisionV1::Refuse, true) => {
                    BaseCoveragePresentedOutcomeV1::ExpectedRefusalMatched
                }
                (BaseE2eExpectedDecisionV1::Unsupported, true) => {
                    BaseCoveragePresentedOutcomeV1::ExpectedUnsupportedMatched
                }
            };
            if coverage_observations
                .insert(coverage_id.clone(), (coverage_outcome, result.root()))
                .is_some()
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.runtime_result",
                    "one result per exact coverage ID",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            match outcome {
                BaseE2eOutcomeV1::Passed => {
                    journey_legacy_eligible_rows = journey_legacy_eligible_rows
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                    journey_legacy_passed_rows = journey_legacy_passed_rows
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                }
                BaseE2eOutcomeV1::Failed => {
                    journey_legacy_eligible_rows = journey_legacy_eligible_rows
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                    journey_legacy_failed_rows = journey_legacy_failed_rows
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                }
                BaseE2eOutcomeV1::Unsupported => {}
                BaseE2eOutcomeV1::NotApplicable => {
                    unreachable!("terminal rows always have a terminal outcome")
                }
            }
            events.push(case_terminal_log_event(
                sequence, journey, row, result, outcome, harness,
            )?);
            sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
        }
        legacy_eligible_rows = legacy_eligible_rows
            .checked_add(journey_legacy_eligible_rows)
            .ok_or_else(sequence_overflow)?;
        legacy_passed_rows = legacy_passed_rows
            .checked_add(journey_legacy_passed_rows)
            .ok_or_else(sequence_overflow)?;
        legacy_failed_rows = legacy_failed_rows
            .checked_add(journey_legacy_failed_rows)
            .ok_or_else(sequence_overflow)?;
        let mut journey_summary_fields = count_fields(
            journey_legacy_eligible_rows,
            journey_legacy_passed_rows,
            journey_legacy_failed_rows,
            journey_report.unsupported(),
        )?;
        journey_summary_fields.extend(partition_fields(
            journey_report.positive_eligible(),
            journey_report.positive_matched(),
            journey_report.expected_refusals(),
            journey_report.expected_refusals_matched(),
            journey_report.unexpected_mismatches(),
        )?);
        journey_summary_fields.extend([
            field(
                "row-count",
                TypedValueV2::U32(
                    u32::try_from(journey.rows().len()).map_err(|_| sequence_overflow())?,
                ),
            )?,
            field(
                "result-count",
                TypedValueV2::U32(
                    u32::try_from(journey_report.results().len())
                        .map_err(|_| sequence_overflow())?,
                ),
            )?,
            field(
                "checked-cells",
                TypedValueV2::U32(journey_report.checked_cells()),
            )?,
        ]);
        events.push(log_event(
            sequence,
            journey,
            None,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eOutcomeV1::NotApplicable,
            harness,
            Some(journey_report.execution_root()),
            journey_summary_fields,
        )?);
        sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
        journey_executions.push(journey_report);
    }
    let retained_artifact_claim = BaseE2eRetainedArtifactClaimV1::Absent;
    let execution_root = projection_execution_root(
        projection,
        harness,
        &journey_executions,
        retained_artifact_claim,
    )?;
    let logging_events_checked = u32::try_from(events.len() + 1).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_projection.logging_events_checked",
            "a u32 deterministic event count",
            events.len() + 1,
        )
    })?;
    let mut summary_fields = count_fields(
        legacy_eligible_rows,
        legacy_passed_rows,
        legacy_failed_rows,
        unsupported,
    )?;
    summary_fields.extend(partition_fields(
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unexpected_mismatches,
    )?);
    summary_fields.extend([
        field(
            "journey-count",
            TypedValueV2::U32(
                u32::try_from(projection.journeys.len()).map_err(|_| sequence_overflow())?,
            ),
        )?,
        field("row-count", TypedValueV2::U32(projection_rows_checked))?,
        field("result-count", TypedValueV2::U32(projection_rows_checked))?,
        field(
            "coverage-source-cases",
            TypedValueV2::U32(
                u32::try_from(projection.coverage_manifest.cases().len())
                    .expect("the frozen inventory is bounded"),
            ),
        )?,
        field(
            "logging-events-checked",
            TypedValueV2::U32(logging_events_checked),
        )?,
        field(
            "projection-e2e-checked",
            TypedValueV2::U32(projection_e2e_checked),
        )?,
        field(
            "source-closure-eligible",
            TypedValueV2::U32(
                source_closure_report
                    .positive_eligible
                    .checked_add(source_closure_report.expected_refusals)
                    .ok_or_else(sequence_overflow)?,
            ),
        )?,
        field(
            "source-closure-failed",
            TypedValueV2::U32(source_closure_report.unexpected_mismatches),
        )?,
        field(
            "source-closure-passed",
            TypedValueV2::U32(
                source_closure_report
                    .positive_matched
                    .checked_add(source_closure_report.expected_refusals_matched)
                    .ok_or_else(sequence_overflow)?,
            ),
        )?,
        field(
            "source-closure-root",
            opaque_root(projection.source_closure.root())?,
        )?,
    ]);
    events.push(projection_summary_log_event(
        sequence,
        projection,
        harness,
        execution_root,
        summary_fields,
    )?);
    let log = BaseE2eLogV1::new(events)?;
    if let Some(event) = log
        .events()
        .iter()
        .find(|event| event.relative_artifact().is_some())
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_projection.retained_artifact",
            "typed absence for every pure local projection event",
            event.logical_sequence(),
        ));
    }
    coverage_observations.insert(
        RUNTIME_LOG_COVERAGE_ID_V1.to_owned(),
        (BaseCoveragePresentedOutcomeV1::PositiveMatched, log.root()),
    );
    for (index, id) in SOURCE_CLOSURE_COVERAGE_IDS_V1.into_iter().enumerate() {
        let matched = source_closure_report.matched_cases[index];
        let expected_outcome = if index == 0 {
            BaseCoveragePresentedOutcomeV1::PositiveMatched
        } else {
            BaseCoveragePresentedOutcomeV1::ExpectedRefusalMatched
        };
        coverage_observations.insert(
            id.to_owned(),
            (
                if matched {
                    expected_outcome
                } else {
                    BaseCoveragePresentedOutcomeV1::UnexpectedMismatch
                },
                source_coverage_evidence_root(id, projection.source_closure.root()),
            ),
        );
    }
    let coverage_report = reconstruct_exact_local_coverage_report(
        &projection.coverage_manifest,
        &coverage_observations,
    )?;
    let source_closure_paths = projection
        .source_closure
        .entries()
        .iter()
        .map(|entry| {
            LogicalBundlePathV1::new(entry.path()).map_err(|_error| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_projection.source_closure_path",
                    "a source-closure-relative path",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_boxed_slice();
    Ok(BaseE2eProjectionReportV1 {
        positive_eligible,
        positive_matched,
        expected_refusals,
        expected_refusals_matched,
        unsupported,
        unexpected_mismatches,
        projection_rows_checked,
        projection_e2e_checked,
        logging_events_checked,
        source_closure_positive_eligible: source_closure_report.positive_eligible,
        source_closure_positive_matched: source_closure_report.positive_matched,
        source_closure_expected_refusals: source_closure_report.expected_refusals,
        source_closure_expected_refusals_matched: source_closure_report.expected_refusals_matched,
        source_closure_unexpected_mismatches: source_closure_report.unexpected_mismatches,
        projection_root: projection.root,
        source_closure_root: projection.source_closure.root(),
        source_root: harness.source.clone(),
        build_root: harness.build.clone(),
        source_closure_paths,
        journey_executions: journey_executions.into_boxed_slice(),
        retained_artifact_claim,
        execution_root,
        log,
        coverage_report,
    })
}

fn source_coverage_evidence_root(id: &str, closure_root: ContentHash) -> ContentHash {
    let mut bytes = Vec::with_capacity(id.len() + closure_root.as_bytes().len());
    bytes.extend_from_slice(id.as_bytes());
    bytes.extend_from_slice(closure_root.as_bytes());
    hash_domain(
        "org.frankensim.fs-evidence-runner.source-coverage-evidence.v1",
        &bytes,
    )
}

/// Exact constructor stage at which a capability mutant refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseE2eCapabilityRefusalStageV1 {
    /// The rights set was not intrinsically legal for any destination mode.
    IntrinsicPolicy = 1,
    /// The rights set was intrinsically legal but not for the selected cell.
    ContextualNarrowing = 2,
}

impl BaseE2eCapabilityRefusalStageV1 {
    /// Stable non-wire stage code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eCapabilityRefusalStageV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::IntrinsicPolicy => "intrinsic-policy",
            Self::ContextualNarrowing => "contextual-narrowing",
        }
    }
}

/// Exact normalized path-set adjudication retained by a detail cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseE2ePathAdjudicationDetailV1 {
    /// The path set was exact.
    Exact,
    /// One exact duplicate path was present.
    Duplicate(String),
    /// One path was a strict segment prefix of another.
    StrictSegmentPrefix {
        /// Prefix path.
        prefix: String,
        /// Strict descendant path.
        descendant: String,
    },
    /// Two paths were ASCII aliases under the Windows profile.
    WindowsAsciiAlias {
        /// First path in canonical input order.
        first: String,
        /// Second aliasing path.
        second: String,
    },
    /// Non-ASCII Windows alias adjudication is explicitly unsupported.
    UnsupportedWindowsNonAsciiAlias(String),
}

/// Closed, typed expected or observed payload carried by a detail cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseE2eDetailPayloadV1 {
    /// Unknown closed-catalog code.
    UnknownCatalog {
        /// Catalog name.
        catalog: &'static str,
        /// Refused numeric code.
        code: u16,
    },
    /// Typed-value constructor refusal.
    Value(ValueError),
    /// Logical-path or object-key refusal.
    Path(PathError),
    /// Complete path-set adjudication.
    PathAdjudication(BaseE2ePathAdjudicationDetailV1),
    /// Exact limit refusal plus its stable diagnostic/repair metadata.
    Limit {
        /// Refusal kind.
        kind: RunnerLimitsViolationKindV2,
        /// Refused limit field.
        field: RunnerLimitFieldV2,
        /// Semantic unit.
        unit: RunnerLimitUnitV2,
        /// Exact admission predicate.
        expected: RunnerLimitExpectationV2,
        /// Exact violating value.
        observed: RunnerLimitValueV2,
        /// Stable diagnostic owner.
        owner: &'static str,
        /// One-based repair rank.
        repair_rank: u8,
        /// Closed non-executable repair kind.
        repair_kind: RepairActionKindV2,
        /// Stable structured repair target.
        repair_target: &'static str,
    },
    /// Exact budget refusal plus its stable diagnostic owner.
    Budget {
        /// Refusal kind.
        kind: RunnerBudgetViolationKindV2,
        /// Refused budget field.
        field: RunnerBudgetFieldV2,
        /// Semantic unit.
        unit: RunnerBudgetUnitV2,
        /// Exact admission predicate.
        expected: RunnerBudgetExpectationV2,
        /// Exact violating value.
        observed: RunnerBudgetValueV2,
        /// Stable diagnostic owner.
        owner: &'static str,
        /// One-based repair rank.
        repair_rank: u8,
        /// Closed non-executable repair kind.
        repair_kind: RepairActionKindV2,
        /// Stable structured repair target.
        repair_target: &'static str,
    },
    /// Bounded generic construction refusal.
    Construction(ConstructionErrorV2),
    /// Exact capability refusal stage and payload.
    Capability {
        /// Stage that refused.
        stage: BaseE2eCapabilityRefusalStageV1,
        /// Exact construction refusal.
        error: ConstructionErrorV2,
    },
    /// Exact lifecycle-state refusal.
    State(StateValidationErrorV2),
    /// Exact NotRun ordinal refusal.
    NotRun(NotRunBasisErrorV2),
    /// Exact nominal-identity refusal.
    Identity(IdentityError),
    /// Sentinel proving that a cell unexpectedly accepted.
    AcceptedInstead,
}

impl BaseE2eDetailPayloadV1 {
    /// Stable descriptive payload-kind name.
    #[must_use]
    pub const fn kind_name(&self) -> &'static str {
        match self {
            Self::UnknownCatalog { .. } => "unknown-catalog",
            Self::Value(_) => "value",
            Self::Path(_) => "path",
            Self::PathAdjudication(_) => "path-adjudication",
            Self::Limit { .. } => "limit",
            Self::Budget { .. } => "budget",
            Self::Construction(_) => "construction",
            Self::Capability { .. } => "capability",
            Self::State(_) => "state",
            Self::NotRun(_) => "not-run",
            Self::Identity(_) => "identity",
            Self::AcceptedInstead => "accepted-instead",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseE2eDetailPayloadV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.kind_name()
    }
}

/// One publicly inspectable expected or observed decision-detail cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eDetailCellV1 {
    kind: BaseE2eCaseKindV1,
    semantic_ordinal: u32,
    stable_id: String,
    decision: BaseE2eExpectedDecisionV1,
    payload: BaseE2eDetailPayloadV1,
    root: ContentHash,
}

impl BaseE2eDetailCellV1 {
    /// Construct one bounded typed observed detail cell.
    ///
    /// The ordinal is one-based and cannot exceed the containing case's
    /// semantic-cell count. `Accept` is represented only by the
    /// `AcceptedInstead` sentinel; refusal and unsupported cells must carry a
    /// typed non-sentinel payload. Limit and budget repair metadata must be
    /// present. Canonical cell encoding is bounded before the root is frozen.
    pub fn new(
        kind: BaseE2eCaseKindV1,
        semantic_ordinal: u32,
        stable_id: &StableTokenV2,
        decision: BaseE2eExpectedDecisionV1,
        payload: BaseE2eDetailPayloadV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if semantic_ordinal == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "base_e2e_detail_cell.semantic_ordinal",
                "a one-based semantic ordinal",
                semantic_ordinal,
            ));
        }
        if semantic_ordinal > kind.semantic_cell_count() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_e2e_detail_cell.semantic_ordinal",
                "an ordinal within the containing case semantic matrix",
                semantic_ordinal,
            ));
        }
        let accepted_sentinel = matches!(&payload, BaseE2eDetailPayloadV1::AcceptedInstead);
        let unsupported_adjudication = matches!(
            &payload,
            BaseE2eDetailPayloadV1::PathAdjudication(
                BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(_)
            )
        );
        let decision_payload_compatible = match decision {
            BaseE2eExpectedDecisionV1::Accept => accepted_sentinel,
            BaseE2eExpectedDecisionV1::Refuse => !accepted_sentinel && !unsupported_adjudication,
            BaseE2eExpectedDecisionV1::Unsupported => unsupported_adjudication,
        };
        if !decision_payload_compatible {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_detail_cell.payload",
                "Accept with AcceptedInstead, Unsupported with the non-ASCII Windows adjudication, or Refuse with another typed refusal payload",
                ConstructionObservedV2::closed_pair(&decision, &payload),
            ));
        }
        match &payload {
            BaseE2eDetailPayloadV1::UnknownCatalog { catalog: "", .. } => {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_e2e_detail_cell.catalog",
                    "a nonempty stable catalog name",
                    ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Absent),
                ));
            }
            BaseE2eDetailPayloadV1::Limit {
                owner,
                repair_rank,
                repair_target,
                ..
            }
            | BaseE2eDetailPayloadV1::Budget {
                owner,
                repair_rank,
                repair_target,
                ..
            } => {
                if !(1..=16).contains(repair_rank) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::OutOfRange,
                        "base_e2e_detail_cell.repair_rank",
                        "an inclusive repair rank from 1 through 16",
                        *repair_rank,
                    ));
                }
                for (field, value) in [
                    ("base_e2e_detail_cell.owner", *owner),
                    ("base_e2e_detail_cell.repair_target", *repair_target),
                ] {
                    StableTokenV2::new(value).map_err(|_error| {
                        ConstructionErrorV2::new_redacted(
                            ConstructionErrorKindV2::Incompatible,
                            field,
                            "a bounded lowercase stable token",
                            ConstructionObservedDataClassV2::CallerControlledText,
                        )
                    })?;
                }
            }
            _ => {}
        }

        let stable_id = stable_id.as_str().to_owned();
        let root =
            checked_detail_cell_root(kind, semantic_ordinal, &stable_id, decision, &payload)?;
        Ok(Self {
            kind,
            semantic_ordinal,
            stable_id,
            decision,
            payload,
            root,
        })
    }

    /// Containing semantic case kind.
    #[must_use]
    pub const fn kind(&self) -> BaseE2eCaseKindV1 {
        self.kind
    }

    /// One-based ordinal in the complete semantic matrix.
    #[must_use]
    pub const fn semantic_ordinal(&self) -> u32 {
        self.semantic_ordinal
    }

    /// Stable human- and agent-inspectable cell identifier.
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    /// Expected or observed three-way decision.
    #[must_use]
    pub const fn decision(&self) -> BaseE2eExpectedDecisionV1 {
        self.decision
    }

    /// Closed typed payload, including field/code/owner/expected/observed/unit
    /// and repair metadata where those concepts apply.
    #[must_use]
    pub const fn payload(&self) -> &BaseE2eDetailPayloadV1 {
        &self.payload
    }

    /// Domain-separated cell root.
    #[must_use]
    pub const fn cell_root(&self) -> ContentHash {
        self.root
    }

    /// Backward-compatible concise root accessor.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.cell_root()
    }
}

/// One bounded descriptor for the first observed detail-manifest divergence.
///
/// A descriptor-only compatibility observation cannot reveal an arbitrary
/// caller-owned typed cell from its root. In that case `observed_cell` is
/// explicitly absent while the exact observed manifest root/count and the
/// first independently expected cell remain bound. This avoids retaining an
/// unbounded mismatched cell vector.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eDetailDivergenceV1 {
    stable_id: String,
    semantic_ordinal: Option<u32>,
    expected_cell: Option<BaseE2eDetailCellV1>,
    observed_cell: Option<BaseE2eDetailCellV1>,
    expected_manifest_root: ContentHash,
    observed_manifest_root: ContentHash,
    expected_cell_count: u32,
    observed_cell_count: u32,
    root: ContentHash,
}

impl BaseE2eDetailDivergenceV1 {
    /// Stable first-divergence identity.
    #[must_use]
    pub fn stable_id(&self) -> &str {
        &self.stable_id
    }

    /// Semantic ordinal when an independently expected cell exists.
    #[must_use]
    pub const fn semantic_ordinal(&self) -> Option<u32> {
        self.semantic_ordinal
    }

    /// Independently expected typed cell at the first unmatched ordinal.
    #[must_use]
    pub const fn expected_cell(&self) -> Option<&BaseE2eDetailCellV1> {
        self.expected_cell.as_ref()
    }

    /// Caller-presented typed observed cell, when available.
    ///
    /// Root/count-only compatibility observations return explicit absence.
    #[must_use]
    pub const fn observed_cell(&self) -> Option<&BaseE2eDetailCellV1> {
        self.observed_cell.as_ref()
    }

    /// Exact independent expected detail-manifest root.
    #[must_use]
    pub const fn expected_manifest_root(&self) -> ContentHash {
        self.expected_manifest_root
    }

    /// Exact caller-presented observed detail-manifest root.
    #[must_use]
    pub const fn observed_manifest_root(&self) -> ContentHash {
        self.observed_manifest_root
    }

    /// Expected detail-cell count.
    #[must_use]
    pub const fn expected_cell_count(&self) -> u32 {
        self.expected_cell_count
    }

    /// Observed detail-cell count.
    #[must_use]
    pub const fn observed_cell_count(&self) -> u32 {
        self.observed_cell_count
    }

    /// Domain-separated bounded divergence-descriptor root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

struct BaseE2eDetailPairBuilderV1 {
    kind: BaseE2eCaseKindV1,
    expected_cells: Vec<BaseE2eDetailCellV1>,
    observed_cells: Vec<BaseE2eDetailCellV1>,
    matched_cells: u32,
    first_divergent_cell: Option<String>,
}

impl BaseE2eDecisionDetailManifestV1 {
    /// Reconstruct one bounded exact descriptor from a complete ordered cell
    /// slice.
    ///
    /// This validates case membership, one-based bounded ordinals, strict
    /// ordinal order, stable-ID uniqueness, cell-root uniqueness, and every
    /// cell's own canonical root before freezing the manifest root.
    pub fn from_cells(
        kind: BaseE2eCaseKindV1,
        cells: &[BaseE2eDetailCellV1],
    ) -> Result<Self, ConstructionErrorV2> {
        validate_detail_cell_sequence(kind, cells)?;
        Ok(detail_manifest_from_cells(kind, cells))
    }

    /// Exact refusal plus unsupported detail-cell count.
    #[must_use]
    pub const fn cell_count(self) -> u32 {
        self.cell_count
    }

    /// Domain-separated ordered detail-manifest root.
    #[must_use]
    pub const fn root(self) -> ContentHash {
        self.root
    }

    /// Explicitly named domain-separated manifest root.
    #[must_use]
    pub const fn manifest_root(self) -> ContentHash {
        self.root
    }

    fn empty(kind: BaseE2eCaseKindV1) -> Self {
        detail_manifest_from_cells(kind, &[])
    }
}

impl BaseE2eDetailPairBuilderV1 {
    fn new(kind: BaseE2eCaseKindV1) -> Self {
        Self {
            kind,
            expected_cells: Vec::new(),
            observed_cells: Vec::new(),
            matched_cells: 0,
            first_divergent_cell: None,
        }
    }

    fn push(
        &mut self,
        semantic_ordinal: u32,
        cell_id: String,
        expected_decision: BaseE2eExpectedDecisionV1,
        expected_payload: BaseE2eDetailPayloadV1,
        observed_decision: BaseE2eExpectedDecisionV1,
        observed_payload: BaseE2eDetailPayloadV1,
    ) {
        let expected = detail_cell(
            self.kind,
            semantic_ordinal,
            cell_id.clone(),
            expected_decision,
            expected_payload,
        );
        let observed = detail_cell(
            self.kind,
            semantic_ordinal,
            cell_id.clone(),
            observed_decision,
            observed_payload,
        );
        if expected == observed {
            self.matched_cells += 1;
        } else if self.first_divergent_cell.is_none() {
            self.first_divergent_cell = Some(cell_id);
        }
        self.expected_cells.push(expected);
        self.observed_cells.push(observed);
    }

    fn finish(self) -> BaseE2eDetailExecutionV1 {
        validate_detail_cell_sequence(self.kind, &self.expected_cells)
            .expect("the internal expected-detail builder emits one canonical global sequence");
        validate_detail_cell_sequence(self.kind, &self.observed_cells)
            .expect("the internal observed-detail builder emits one canonical global sequence");
        let expected = detail_manifest_from_cells(self.kind, &self.expected_cells);
        let observed = detail_manifest_from_cells(self.kind, &self.observed_cells);
        BaseE2eDetailExecutionV1 {
            expected,
            observed,
            expected_cells: Some(self.expected_cells.into_boxed_slice()),
            observed_cells: Some(self.observed_cells.into_boxed_slice()),
            matched_cells: self.matched_cells,
            first_divergent_cell: self.first_divergent_cell,
        }
    }
}

fn detail_manifest_from_cells(
    kind: BaseE2eCaseKindV1,
    cells: &[BaseE2eDetailCellV1],
) -> BaseE2eDecisionDetailManifestV1 {
    let mut bytes = Vec::with_capacity(6 + cells.len() * 32);
    bytes.extend_from_slice(&kind.code().to_be_bytes());
    bytes.extend_from_slice(
        &u32::try_from(cells.len())
            .expect("base detail count is bounded by the semantic manifest")
            .to_be_bytes(),
    );
    for cell in cells {
        bytes.extend_from_slice(cell.root.as_bytes());
    }
    BaseE2eDecisionDetailManifestV1 {
        cell_count: u32::try_from(cells.len())
            .expect("base detail count is bounded by the semantic manifest"),
        root: hash_domain(BASE_E2E_DECISION_DETAIL_MANIFEST_DOMAIN_V1, &bytes),
    }
}

fn validate_detail_cell_sequence(
    kind: BaseE2eCaseKindV1,
    cells: &[BaseE2eDetailCellV1],
) -> Result<(), ConstructionErrorV2> {
    let cell_count = u32::try_from(cells.len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_detail_manifest.cells",
            "a u32-bounded detail-cell slice",
            cells.len(),
        )
    })?;
    if cell_count > kind.semantic_cell_count() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_detail_manifest.cells",
            "no more cells than the containing semantic matrix",
            ConstructionObservedV2::unsigned_pair(
                u64::from(cell_count),
                u64::from(kind.semantic_cell_count()),
            ),
        ));
    }
    let mut stable_ids = std::collections::BTreeSet::new();
    let mut roots = std::collections::BTreeSet::new();
    let mut previous_ordinal = None;
    for cell in cells {
        if cell.kind() != kind {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_detail_manifest.cell_kind",
                "the containing detail-manifest case kind",
                ConstructionObservedV2::closed(&cell.kind()),
            ));
        }
        if cell.semantic_ordinal() == 0 || cell.semantic_ordinal() > kind.semantic_cell_count() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_e2e_detail_manifest.semantic_ordinal",
                "a one-based ordinal within the containing semantic matrix",
                cell.semantic_ordinal(),
            ));
        }
        if let Some(previous) = previous_ordinal {
            if cell.semantic_ordinal() == previous {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_detail_manifest.semantic_ordinal",
                    "one cell per semantic ordinal",
                    cell.semantic_ordinal(),
                ));
            }
            if cell.semantic_ordinal() < previous {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "base_e2e_detail_manifest.semantic_ordinal",
                    "strictly increasing semantic ordinals",
                    ConstructionObservedV2::unsigned_pair(
                        u64::from(cell.semantic_ordinal()),
                        u64::from(previous),
                    ),
                ));
            }
        }
        previous_ordinal = Some(cell.semantic_ordinal());
        if !stable_ids.insert(cell.stable_id()) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_detail_manifest.stable_id",
                "one unique stable ID per detail cell",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        if !roots.insert(cell.root()) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_detail_manifest.cell_root",
                "one unique canonical root per detail cell",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let recomputed = checked_detail_cell_root(
            cell.kind(),
            cell.semantic_ordinal(),
            cell.stable_id(),
            cell.decision(),
            cell.payload(),
        )?;
        if recomputed != cell.root() {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_detail_manifest.cell_root",
                "the exact root reconstructed from the typed detail cell",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    Ok(())
}

fn checked_detail_cell_root(
    kind: BaseE2eCaseKindV1,
    semantic_ordinal: u32,
    stable_id: &str,
    decision: BaseE2eExpectedDecisionV1,
    payload: &BaseE2eDetailPayloadV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut bytes = Vec::with_capacity(256);
    bytes.extend_from_slice(&kind.code().to_be_bytes());
    bytes.extend_from_slice(&semantic_ordinal.to_be_bytes());
    detail_push_str(&mut bytes, stable_id);
    bytes.extend_from_slice(&decision.code().to_be_bytes());
    encode_detail_payload(&mut bytes, payload);
    if bytes.len() > BASE_E2E_DETAIL_CELL_MAX_ENCODED_BYTES_V1 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_detail_cell.encoded_bytes",
            "a bounded canonical detail-cell encoding",
            bytes.len(),
        ));
    }
    Ok(hash_domain(BASE_E2E_DECISION_DETAIL_CELL_DOMAIN_V1, &bytes))
}

fn detail_cell(
    kind: BaseE2eCaseKindV1,
    semantic_ordinal: u32,
    stable_id: String,
    decision: BaseE2eExpectedDecisionV1,
    payload: BaseE2eDetailPayloadV1,
) -> BaseE2eDetailCellV1 {
    let root = checked_detail_cell_root(kind, semantic_ordinal, &stable_id, decision, &payload)
        .expect("the frozen bounded internal detail-cell table is valid");
    BaseE2eDetailCellV1 {
        kind,
        semantic_ordinal,
        stable_id,
        decision,
        payload,
        root,
    }
}

fn detail_divergence(
    row: &BaseE2eProjectionRowV1,
    observed_manifest_root: ContentHash,
    observed_cell_count: u32,
    matched_cells: u32,
    first_unexpected_cell: Option<&str>,
    first_observed_detail_cell: Option<&BaseE2eDetailCellV1>,
) -> BaseE2eDetailDivergenceV1 {
    let expected_cell = first_unexpected_cell
        .and_then(|stable_id| {
            row.expected_detail_cells()
                .iter()
                .find(|cell| cell.stable_id() == stable_id)
        })
        .or_else(|| {
            first_observed_detail_cell.and_then(|observed| {
                row.expected_detail_cells()
                    .iter()
                    .find(|expected| expected.semantic_ordinal() == observed.semantic_ordinal())
            })
        })
        .cloned();
    let stable_id = expected_cell
        .as_ref()
        .map(|cell| cell.stable_id().to_owned())
        .or_else(|| first_unexpected_cell.map(str::to_owned))
        .unwrap_or_else(|| "detail.manifest".to_owned());
    let semantic_ordinal = expected_cell
        .as_ref()
        .map(BaseE2eDetailCellV1::semantic_ordinal);
    let observed_cell = first_observed_detail_cell.cloned();
    let mut bytes = Vec::with_capacity(256);
    detail_push_str(&mut bytes, &stable_id);
    bytes.extend_from_slice(&semantic_ordinal.unwrap_or(0).to_be_bytes());
    bytes.extend_from_slice(row.expected_detail_manifest_root().as_bytes());
    bytes.extend_from_slice(observed_manifest_root.as_bytes());
    bytes.extend_from_slice(&row.expected_detail_cell_count().to_be_bytes());
    bytes.extend_from_slice(&observed_cell_count.to_be_bytes());
    bytes.extend_from_slice(&matched_cells.to_be_bytes());
    bytes.extend_from_slice(
        expected_cell
            .as_ref()
            .map_or_else(
                || hash_domain(BASE_E2E_DECISION_DETAIL_CELL_DOMAIN_V1, b"absent"),
                BaseE2eDetailCellV1::root,
            )
            .as_bytes(),
    );
    bytes.extend_from_slice(
        observed_cell
            .as_ref()
            .map_or_else(
                || {
                    hash_domain(
                        BASE_E2E_DECISION_DETAIL_CELL_DOMAIN_V1,
                        b"observed-cell-unavailable",
                    )
                },
                BaseE2eDetailCellV1::root,
            )
            .as_bytes(),
    );
    let root = hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-detail-divergence.v1",
        &bytes,
    );
    BaseE2eDetailDivergenceV1 {
        stable_id,
        semantic_ordinal,
        expected_cell,
        observed_cell,
        expected_manifest_root: row.expected_detail_manifest_root(),
        observed_manifest_root,
        expected_cell_count: row.expected_detail_cell_count(),
        observed_cell_count,
        root,
    }
}

fn detail_push_str(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(
        &u32::try_from(value.len())
            .expect("base detail strings are bounded fixture data")
            .to_be_bytes(),
    );
    bytes.extend_from_slice(value.as_bytes());
}

fn detail_push_bool(bytes: &mut Vec<u8>, value: bool) {
    bytes.push(u8::from(value));
}

fn encode_detail_payload(bytes: &mut Vec<u8>, payload: &BaseE2eDetailPayloadV1) {
    match payload {
        BaseE2eDetailPayloadV1::UnknownCatalog { catalog, code } => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
            detail_push_str(bytes, catalog);
            bytes.extend_from_slice(&code.to_be_bytes());
        }
        BaseE2eDetailPayloadV1::Value(error) => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            encode_value_error(bytes, error);
        }
        BaseE2eDetailPayloadV1::Path(error) => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            encode_path_error(bytes, error);
        }
        BaseE2eDetailPayloadV1::PathAdjudication(detail) => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            encode_path_adjudication(bytes, detail);
        }
        BaseE2eDetailPayloadV1::Limit {
            kind,
            field,
            unit,
            expected,
            observed,
            owner,
            repair_rank,
            repair_kind,
            repair_target,
        } => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            bytes.extend_from_slice(&limit_violation_kind_code(*kind).to_be_bytes());
            bytes.extend_from_slice(&field.ordinal().to_be_bytes());
            bytes.extend_from_slice(&limit_unit_code(*unit).to_be_bytes());
            encode_limit_expectation(bytes, *expected);
            encode_limit_value(bytes, *observed);
            detail_push_str(bytes, owner);
            bytes.push(*repair_rank);
            bytes.extend_from_slice(&repair_kind.code().to_be_bytes());
            detail_push_str(bytes, repair_target);
        }
        BaseE2eDetailPayloadV1::Budget {
            kind,
            field,
            unit,
            expected,
            observed,
            owner,
            repair_rank,
            repair_kind,
            repair_target,
        } => {
            bytes.extend_from_slice(&6_u16.to_be_bytes());
            bytes.extend_from_slice(&budget_violation_kind_code(*kind).to_be_bytes());
            bytes.extend_from_slice(&field.ordinal().to_be_bytes());
            bytes.extend_from_slice(&budget_unit_code(*unit).to_be_bytes());
            encode_budget_expectation(bytes, *expected);
            encode_budget_value(bytes, *observed);
            detail_push_str(bytes, owner);
            bytes.push(*repair_rank);
            bytes.extend_from_slice(&repair_kind.code().to_be_bytes());
            detail_push_str(bytes, repair_target);
        }
        BaseE2eDetailPayloadV1::Construction(error) => {
            bytes.extend_from_slice(&7_u16.to_be_bytes());
            encode_construction_error(bytes, error);
        }
        BaseE2eDetailPayloadV1::Capability { stage, error } => {
            bytes.extend_from_slice(&8_u16.to_be_bytes());
            bytes.extend_from_slice(&(*stage as u16).to_be_bytes());
            encode_construction_error(bytes, error);
        }
        BaseE2eDetailPayloadV1::State(error) => {
            bytes.extend_from_slice(&9_u16.to_be_bytes());
            encode_state_error(bytes, *error);
        }
        BaseE2eDetailPayloadV1::NotRun(error) => {
            bytes.extend_from_slice(&10_u16.to_be_bytes());
            encode_not_run_error(bytes, *error);
        }
        BaseE2eDetailPayloadV1::Identity(error) => {
            bytes.extend_from_slice(&11_u16.to_be_bytes());
            encode_identity_error(bytes, error);
        }
        BaseE2eDetailPayloadV1::AcceptedInstead => {
            bytes.extend_from_slice(&12_u16.to_be_bytes());
        }
    }
}

fn encode_value_error(bytes: &mut Vec<u8>, error: &ValueError) {
    match error {
        ValueError::ZeroRationalDenominator => bytes.extend_from_slice(&1_u16.to_be_bytes()),
        ValueError::NonCanonicalRational => bytes.extend_from_slice(&2_u16.to_be_bytes()),
        ValueError::DecimalScaleOutOfRange { observed } => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(&observed.to_be_bytes());
        }
        ValueError::DecimalNormalizationScaleOutOfRange => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
        }
        ValueError::NonCanonicalDecimal => bytes.extend_from_slice(&5_u16.to_be_bytes()),
        ValueError::UnitScaleNotPositive => bytes.extend_from_slice(&6_u16.to_be_bytes()),
        ValueError::StableTokenEmpty => bytes.extend_from_slice(&7_u16.to_be_bytes()),
        ValueError::StableTokenTooLong { observed, maximum } => {
            bytes.extend_from_slice(&8_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*maximum)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        ValueError::StableTokenInvalidByte { index, byte } => {
            bytes.extend_from_slice(&9_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*index)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.push(*byte);
        }
        ValueError::StableTokenEmptySegment { index } => {
            bytes.extend_from_slice(&10_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*index)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        ValueError::TextTooLong { observed, maximum } => {
            bytes.extend_from_slice(&11_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*maximum)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        ValueError::OpaqueBytesTooLong { observed, maximum } => {
            bytes.extend_from_slice(&12_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*maximum)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
    }
}

fn encode_path_error(bytes: &mut Vec<u8>, error: &PathError) {
    match error {
        PathError::Empty => bytes.extend_from_slice(&1_u16.to_be_bytes()),
        PathError::Absolute => bytes.extend_from_slice(&2_u16.to_be_bytes()),
        PathError::DriveDesignator => bytes.extend_from_slice(&3_u16.to_be_bytes()),
        PathError::Backslash { index } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*index)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::Nul { index } => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*index)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::EmptySegment { segment } => {
            bytes.extend_from_slice(&6_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*segment)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::DotSegment { segment } => {
            bytes.extend_from_slice(&7_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*segment)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::DotDotSegment { segment } => {
            bytes.extend_from_slice(&8_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*segment)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::TooManyBytes { observed, maximum } => {
            bytes.extend_from_slice(&9_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*maximum)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::TooManySegments { observed, maximum } => {
            bytes.extend_from_slice(&10_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*maximum)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        PathError::ReservedContentStorePrefix { prefix } => {
            bytes.extend_from_slice(&11_u16.to_be_bytes());
            detail_push_str(bytes, prefix);
        }
    }
}

fn encode_path_adjudication(bytes: &mut Vec<u8>, detail: &BaseE2ePathAdjudicationDetailV1) {
    match detail {
        BaseE2ePathAdjudicationDetailV1::Exact => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
        }
        BaseE2ePathAdjudicationDetailV1::Duplicate(path) => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            detail_push_str(bytes, path);
        }
        BaseE2ePathAdjudicationDetailV1::StrictSegmentPrefix { prefix, descendant } => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            detail_push_str(bytes, prefix);
            detail_push_str(bytes, descendant);
        }
        BaseE2ePathAdjudicationDetailV1::WindowsAsciiAlias { first, second } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            detail_push_str(bytes, first);
            detail_push_str(bytes, second);
        }
        BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(path) => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            detail_push_str(bytes, path);
        }
    }
}

const fn construction_error_kind_code(kind: ConstructionErrorKindV2) -> u16 {
    match kind {
        ConstructionErrorKindV2::Missing => 1,
        ConstructionErrorKindV2::Unexpected => 2,
        ConstructionErrorKindV2::UnknownCode => 3,
        ConstructionErrorKindV2::Zero => 4,
        ConstructionErrorKindV2::Duplicate => 5,
        ConstructionErrorKindV2::OutOfOrder => 6,
        ConstructionErrorKindV2::OutOfRange => 7,
        ConstructionErrorKindV2::ArithmeticOverflow => 8,
        ConstructionErrorKindV2::Incompatible => 9,
        ConstructionErrorKindV2::TooLarge => 10,
        ConstructionErrorKindV2::Unsupported => 11,
    }
}

fn encode_construction_error(bytes: &mut Vec<u8>, error: &ConstructionErrorV2) {
    bytes.extend_from_slice(&construction_error_kind_code(error.kind()).to_be_bytes());
    detail_push_str(bytes, error.field());
    detail_push_str(bytes, error.expected());
    detail_push_str(bytes, error.observed());
}

const fn limit_violation_kind_code(kind: RunnerLimitsViolationKindV2) -> u16 {
    match kind {
        RunnerLimitsViolationKindV2::WrongWidth => 1,
        RunnerLimitsViolationKindV2::ExceedsBaseCeiling => 2,
        RunnerLimitsViolationKindV2::FixedFieldChanged => 3,
        RunnerLimitsViolationKindV2::BelowStructuralMinimum => 4,
        RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder => 5,
        RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum => 6,
        RunnerLimitsViolationKindV2::DeclaredMinimumUnmet => 7,
        RunnerLimitsViolationKindV2::ExecutableCaseSetEmpty => 8,
        RunnerLimitsViolationKindV2::NonExecutableCaseSetPresent => 9,
        RunnerLimitsViolationKindV2::CaseCountExceeded => 10,
        RunnerLimitsViolationKindV2::FamilyRowsExceeded => 11,
        RunnerLimitsViolationKindV2::ArithmeticOverflow => 12,
        RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient => 13,
        RunnerLimitsViolationKindV2::JointFeasibilityViolation => 14,
        RunnerLimitsViolationKindV2::ProtocolStoredLengthMismatch => 15,
        RunnerLimitsViolationKindV2::EnvelopeOverheadExceeded => 16,
        RunnerLimitsViolationKindV2::ArtifactCountExceeded => 17,
        RunnerLimitsViolationKindV2::SystemObjectSetMismatch => 18,
        RunnerLimitsViolationKindV2::AggregateMismatch => 19,
    }
}

const fn limit_unit_code(unit: RunnerLimitUnitV2) -> u16 {
    match unit {
        RunnerLimitUnitV2::Count => 1,
        RunnerLimitUnitV2::Records => 2,
        RunnerLimitUnitV2::Rows => 3,
        RunnerLimitUnitV2::EncodedBytes => 4,
        RunnerLimitUnitV2::ExpandedBytes => 5,
        RunnerLimitUnitV2::StoredBytes => 6,
        RunnerLimitUnitV2::LogicalBytes => 7,
        RunnerLimitUnitV2::Depth => 8,
        RunnerLimitUnitV2::Nodes => 9,
        RunnerLimitUnitV2::Digits => 10,
        RunnerLimitUnitV2::Segments => 11,
        RunnerLimitUnitV2::Diagnostics => 12,
        RunnerLimitUnitV2::Prerequisites => 13,
        RunnerLimitUnitV2::Repairs => 14,
        RunnerLimitUnitV2::Artifacts => 15,
        RunnerLimitUnitV2::Namespaces => 16,
        RunnerLimitUnitV2::Classes => 17,
        RunnerLimitUnitV2::Visits => 18,
        RunnerLimitUnitV2::DecimalScale => 19,
    }
}

fn encode_limit_value(bytes: &mut Vec<u8>, value: RunnerLimitValueV2) {
    match value {
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

fn encode_limit_expectation(bytes: &mut Vec<u8>, expectation: RunnerLimitExpectationV2) {
    match expectation {
        RunnerLimitExpectationV2::Width(width) => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
            bytes.extend_from_slice(
                &(match width {
                    crate::limits::RunnerLimitWidthV2::U32 => 1_u16,
                    crate::limits::RunnerLimitWidthV2::U64 => 2_u16,
                })
                .to_be_bytes(),
            );
        }
        RunnerLimitExpectationV2::AtMost(value) => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            encode_limit_value(bytes, value);
        }
        RunnerLimitExpectationV2::AtLeast(value) => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            encode_limit_value(bytes, value);
        }
        RunnerLimitExpectationV2::Exactly(value) => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            encode_limit_value(bytes, value);
        }
        RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
        }
    }
}

const fn budget_violation_kind_code(kind: RunnerBudgetViolationKindV2) -> u16 {
    match kind {
        RunnerBudgetViolationKindV2::Zero => 1,
        RunnerBudgetViolationKindV2::ParallelChildrenExceedTotal => 2,
        RunnerBudgetViolationKindV2::TimeoutSumOverflow => 3,
        RunnerBudgetViolationKindV2::TimeoutSumExceedsWall => 4,
        RunnerBudgetViolationKindV2::ProfileCeilingExceeded => 5,
        RunnerBudgetViolationKindV2::LimitExceeded => 6,
        RunnerBudgetViolationKindV2::ContextualZeroRequired => 7,
        RunnerBudgetViolationKindV2::ContextualNonZeroRequired => 8,
        RunnerBudgetViolationKindV2::CommandResultCannotContainLifecycle => 9,
        RunnerBudgetViolationKindV2::ArtifactStoredBelowEncoded => 10,
        RunnerBudgetViolationKindV2::PublicationSumOverflow => 11,
        RunnerBudgetViolationKindV2::PublicationEquationMismatch => 12,
        RunnerBudgetViolationKindV2::UnregisteredLogicalWorkUnit => 13,
    }
}

const fn budget_unit_code(unit: RunnerBudgetUnitV2) -> u16 {
    match unit {
        RunnerBudgetUnitV2::Nanoseconds => 1,
        RunnerBudgetUnitV2::LogicalBytes => 2,
        RunnerBudgetUnitV2::Count => 3,
        RunnerBudgetUnitV2::LogicalWork => 4,
        RunnerBudgetUnitV2::LogicalWorkUnit => 5,
        RunnerBudgetUnitV2::EncodedBytes => 6,
        RunnerBudgetUnitV2::StoredBytes => 7,
        RunnerBudgetUnitV2::ExpandedBytes => 8,
    }
}

fn encode_budget_value(bytes: &mut Vec<u8>, value: RunnerBudgetValueV2) {
    match value {
        RunnerBudgetValueV2::U32(value) => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::U64(value) => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::U128(value) => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::LogicalUnit { tag, registered_id } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            bytes.extend_from_slice(&tag.to_be_bytes());
            detail_push_bool(bytes, registered_id.is_some());
            if let Some(registered_id) = registered_id {
                bytes.extend_from_slice(&registered_id.to_be_bytes());
            }
        }
    }
}

fn encode_budget_expectation(bytes: &mut Vec<u8>, expectation: RunnerBudgetExpectationV2) {
    match expectation {
        RunnerBudgetExpectationV2::NonZero => bytes.extend_from_slice(&1_u16.to_be_bytes()),
        RunnerBudgetExpectationV2::Zero => bytes.extend_from_slice(&2_u16.to_be_bytes()),
        RunnerBudgetExpectationV2::AtMost(value) => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            encode_budget_value(bytes, value);
        }
        RunnerBudgetExpectationV2::AtLeast(value) => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            encode_budget_value(bytes, value);
        }
        RunnerBudgetExpectationV2::Exactly(value) => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            encode_budget_value(bytes, value);
        }
        RunnerBudgetExpectationV2::RegisteredInExtensionRegistry => {
            bytes.extend_from_slice(&6_u16.to_be_bytes());
        }
    }
}

fn encode_state_error(bytes: &mut Vec<u8>, error: StateValidationErrorV2) {
    match error {
        StateValidationErrorV2::StateNotAllowedForRole { role, state } => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
            bytes.extend_from_slice(&role.code().to_be_bytes());
            bytes.extend_from_slice(&state.code().to_be_bytes());
        }
        StateValidationErrorV2::MissingRefusedReason => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
        }
        StateValidationErrorV2::UnexpectedRefusedReason { state, observed } => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(&state.code().to_be_bytes());
            bytes.extend_from_slice(&observed.code().to_be_bytes());
        }
        StateValidationErrorV2::MissingDiagnostic { expected } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            bytes.extend_from_slice(&expected.code().to_be_bytes());
        }
        StateValidationErrorV2::UnexpectedDiagnostic { observed } => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            bytes.extend_from_slice(&observed.code().to_be_bytes());
        }
        StateValidationErrorV2::WrongDiagnostic { expected, observed } => {
            bytes.extend_from_slice(&6_u16.to_be_bytes());
            bytes.extend_from_slice(&expected.code().to_be_bytes());
            bytes.extend_from_slice(&observed.code().to_be_bytes());
        }
        StateValidationErrorV2::MissingDrainBasis { expected } => {
            bytes.extend_from_slice(&7_u16.to_be_bytes());
            bytes.extend_from_slice(&drain_code(expected).to_be_bytes());
        }
        StateValidationErrorV2::UnexpectedDrainBasis { state, observed } => {
            bytes.extend_from_slice(&8_u16.to_be_bytes());
            bytes.extend_from_slice(&state.code().to_be_bytes());
            bytes.extend_from_slice(&drain_code(observed).to_be_bytes());
        }
        StateValidationErrorV2::WrongDrainBasis { expected, observed } => {
            bytes.extend_from_slice(&9_u16.to_be_bytes());
            bytes.extend_from_slice(&drain_code(expected).to_be_bytes());
            bytes.extend_from_slice(&drain_code(observed).to_be_bytes());
        }
    }
}

fn encode_not_run_error(bytes: &mut Vec<u8>, error: NotRunBasisErrorV2) {
    match error {
        NotRunBasisErrorV2::EmptyManifest => bytes.extend_from_slice(&1_u16.to_be_bytes()),
        NotRunBasisErrorV2::ManifestCaseCountExceedsMaximum { observed, maximum } => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            bytes.extend_from_slice(&observed.to_be_bytes());
            bytes.extend_from_slice(&maximum.to_be_bytes());
        }
        NotRunBasisErrorV2::LowestRemainingOrdinalOutOfRange {
            observed,
            ordered_case_count,
        } => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(&observed.to_be_bytes());
            bytes.extend_from_slice(&ordered_case_count.to_be_bytes());
        }
    }
}

fn encode_identity_error(bytes: &mut Vec<u8>, error: &IdentityError) {
    match error {
        IdentityError::InvalidDomain => bytes.extend_from_slice(&1_u16.to_be_bytes()),
        IdentityError::WrongDigestLength { observed, expected } => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*expected)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        IdentityError::WrongLowerHexLength { observed, expected } => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*observed)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.extend_from_slice(
                &u64::try_from(*expected)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
        }
        IdentityError::NonCanonicalLowerHex { index, byte } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            bytes.extend_from_slice(
                &u64::try_from(*index)
                    .expect("usize fits u64 on supported targets")
                    .to_be_bytes(),
            );
            bytes.push(*byte);
        }
        IdentityError::WrongRole { expected, observed } => {
            bytes.extend_from_slice(&5_u16.to_be_bytes());
            bytes.extend_from_slice(&expected.code().to_be_bytes());
            bytes.extend_from_slice(&observed.code().to_be_bytes());
        }
        IdentityError::WrongDomain { expected, observed } => {
            bytes.extend_from_slice(&6_u16.to_be_bytes());
            detail_push_str(bytes, expected);
            detail_push_str(bytes, observed);
        }
    }
}

fn normalize_limit_violation(error: RunnerLimitsViolationV2) -> BaseE2eDetailPayloadV1 {
    BaseE2eDetailPayloadV1::Limit {
        kind: error.kind(),
        field: error.field(),
        unit: error.unit(),
        expected: error.expected(),
        observed: error.observed(),
        owner: error.owner(),
        repair_rank: error.repair_rank(),
        repair_kind: error.repair_kind(),
        repair_target: error.repair_target(),
    }
}

fn normalize_budget_violation(error: RunnerBudgetViolationV2) -> BaseE2eDetailPayloadV1 {
    BaseE2eDetailPayloadV1::Budget {
        kind: error.kind(),
        field: error.field(),
        unit: error.unit(),
        expected: error.expected(),
        observed: error.observed(),
        owner: error.owner(),
        repair_rank: error.repair_rank(),
        repair_kind: error.repair_kind(),
        repair_target: error.repair_target(),
    }
}

fn normalize_path_adjudication(
    adjudication: PathSetAdjudicationV1<'_>,
) -> BaseE2ePathAdjudicationDetailV1 {
    match adjudication {
        PathSetAdjudicationV1::Exact => BaseE2ePathAdjudicationDetailV1::Exact,
        PathSetAdjudicationV1::Duplicate { path } => {
            BaseE2ePathAdjudicationDetailV1::Duplicate(path.to_owned())
        }
        PathSetAdjudicationV1::StrictSegmentPrefix { prefix, descendant } => {
            BaseE2ePathAdjudicationDetailV1::StrictSegmentPrefix {
                prefix: prefix.to_owned(),
                descendant: descendant.to_owned(),
            }
        }
        PathSetAdjudicationV1::WindowsAsciiAlias { first, second } => {
            BaseE2ePathAdjudicationDetailV1::WindowsAsciiAlias {
                first: first.to_owned(),
                second: second.to_owned(),
            }
        }
        PathSetAdjudicationV1::UnsupportedWindowsNonAsciiAlias { path } => {
            BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(path.to_owned())
        }
    }
}

struct BaseE2eExpectedDetailCatalogV1 {
    manifest: BaseE2eDecisionDetailManifestV1,
    cells: Box<[BaseE2eDetailCellV1]>,
}

static EXPECTED_DETAIL_CATALOGS_V1: std::sync::OnceLock<[BaseE2eExpectedDetailCatalogV1; 24]> =
    std::sync::OnceLock::new();

const EXPECTED_DETAIL_CASE_KINDS_ORACLE_V1: [BaseE2eCaseKindV1; 24] = [
    BaseE2eCaseKindV1::CatalogLiterals,
    BaseE2eCaseKindV1::UnknownCatalogCode,
    BaseE2eCaseKindV1::CanonicalRational,
    BaseE2eCaseKindV1::OverlongStableToken,
    BaseE2eCaseKindV1::LogicalPath,
    BaseE2eCaseKindV1::ReservedContentStorePrefix,
    BaseE2eCaseKindV1::WindowsUnicodeAlias,
    BaseE2eCaseKindV1::LimitCatalog,
    BaseE2eCaseKindV1::BudgetAdmission,
    BaseE2eCaseKindV1::BudgetChildRelation,
    BaseE2eCaseKindV1::PublicationSelection,
    BaseE2eCaseKindV1::PublicationCrossCell,
    BaseE2eCaseKindV1::CapabilityLeastPrivilege,
    BaseE2eCaseKindV1::CapabilityExtraRight,
    BaseE2eCaseKindV1::StatePass,
    BaseE2eCaseKindV1::StateUsageInLifecycle,
    BaseE2eCaseKindV1::Diagnostic,
    BaseE2eCaseKindV1::DiagnosticRankGap,
    BaseE2eCaseKindV1::IdentityMutation,
    BaseE2eCaseKindV1::NoClaimNominality,
    BaseE2eCaseKindV1::AtomicResult,
    BaseE2eCaseKindV1::AtomicResultPresence,
    BaseE2eCaseKindV1::PublicationStorage,
    BaseE2eCaseKindV1::CommandList,
];

fn expected_detail_catalog(kind: BaseE2eCaseKindV1) -> &'static BaseE2eExpectedDetailCatalogV1 {
    let catalogs = EXPECTED_DETAIL_CATALOGS_V1.get_or_init(|| {
        EXPECTED_DETAIL_CASE_KINDS_ORACLE_V1.map(|case_kind| {
            let mut detail = expected_detail_execution(case_kind);
            BaseE2eExpectedDetailCatalogV1 {
                manifest: detail.expected,
                cells: detail
                    .expected_cells
                    .take()
                    .expect("oracle construction always retains expected cells"),
            }
        })
    });
    &catalogs[usize::from(kind.code() - 1)]
}

fn expected_detail_manifest(kind: BaseE2eCaseKindV1) -> BaseE2eDecisionDetailManifestV1 {
    expected_detail_catalog(kind).manifest
}

fn expected_detail_cells(kind: BaseE2eCaseKindV1) -> &'static [BaseE2eDetailCellV1] {
    &expected_detail_catalog(kind).cells
}

fn execute_detail_manifest(
    kind: BaseE2eCaseKindV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> BaseE2eDetailExecutionV1 {
    decision_detail_execution(kind, Some(harness))
}

fn push_detail_observation(
    builder: &mut BaseE2eDetailPairBuilderV1,
    semantic_ordinal: u32,
    cell_id: String,
    expected_decision: BaseE2eExpectedDecisionV1,
    expected_payload: BaseE2eDetailPayloadV1,
    observed: Option<(BaseE2eExpectedDecisionV1, BaseE2eDetailPayloadV1)>,
) {
    let (observed_decision, observed_payload) =
        observed.unwrap_or_else(|| (expected_decision, expected_payload.clone()));
    builder.push(
        semantic_ordinal,
        cell_id,
        expected_decision,
        expected_payload,
        observed_decision,
        observed_payload,
    );
}

fn expected_detail_execution(kind: BaseE2eCaseKindV1) -> BaseE2eDetailExecutionV1 {
    let mut builder = BaseE2eDetailPairBuilderV1::new(kind);
    match kind {
        BaseE2eCaseKindV1::CatalogLiterals
        | BaseE2eCaseKindV1::CanonicalRational
        | BaseE2eCaseKindV1::LogicalPath
        | BaseE2eCaseKindV1::PublicationSelection
        | BaseE2eCaseKindV1::CapabilityLeastPrivilege
        | BaseE2eCaseKindV1::Diagnostic
        | BaseE2eCaseKindV1::AtomicResult
        | BaseE2eCaseKindV1::PublicationStorage
        | BaseE2eCaseKindV1::CommandList => {}
        BaseE2eCaseKindV1::UnknownCatalogCode => push_detail_observation(
            &mut builder,
            1,
            "catalog.unknown-code".to_owned(),
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::UnknownCatalog {
                catalog: "ProofExitV2",
                code: 1,
            },
            None,
        ),
        BaseE2eCaseKindV1::OverlongStableToken => push_detail_observation(
            &mut builder,
            1,
            "value.overlong-token".to_owned(),
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::Value(ValueError::StableTokenTooLong {
                observed: 129,
                maximum: 128,
            }),
            None,
        ),
        BaseE2eCaseKindV1::ReservedContentStorePrefix => push_detail_observation(
            &mut builder,
            1,
            "path.reserved-prefix".to_owned(),
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::Path(PathError::ReservedContentStorePrefix {
                prefix: "__runner_",
            }),
            None,
        ),
        BaseE2eCaseKindV1::WindowsUnicodeAlias => push_detail_observation(
            &mut builder,
            1,
            "path.windows-unicode-alias".to_owned(),
            BaseE2eExpectedDecisionV1::Unsupported,
            BaseE2eDetailPayloadV1::PathAdjudication(
                BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(
                    "résumé/a".to_owned(),
                ),
            ),
            None,
        ),
        BaseE2eCaseKindV1::LimitCatalog => append_limit_detail_cells(&mut builder, false),
        BaseE2eCaseKindV1::BudgetAdmission => {
            append_budget_admission_detail_cells(&mut builder, false);
        }
        BaseE2eCaseKindV1::BudgetChildRelation => {
            append_budget_child_detail(&mut builder, false);
        }
        BaseE2eCaseKindV1::PublicationCrossCell => {
            append_publication_cross_detail(&mut builder, false);
        }
        BaseE2eCaseKindV1::CapabilityExtraRight => {
            append_capability_detail_cells(&mut builder, None);
        }
        BaseE2eCaseKindV1::StatePass => append_state_detail_cells(&mut builder, false),
        BaseE2eCaseKindV1::StateUsageInLifecycle => {
            append_state_usage_detail(&mut builder, false);
        }
        BaseE2eCaseKindV1::DiagnosticRankGap => {
            append_diagnostic_rank_detail(&mut builder, None);
        }
        BaseE2eCaseKindV1::IdentityMutation => {
            append_identity_detail_cells(&mut builder, None);
        }
        BaseE2eCaseKindV1::NoClaimNominality => {
            append_no_claim_detail_cells(&mut builder, None);
        }
        BaseE2eCaseKindV1::AtomicResultPresence => {
            append_atomic_presence_detail(&mut builder, false);
        }
    }
    builder.finish()
}

#[allow(
    clippy::too_many_lines,
    reason = "the decision-detail executor mirrors the closed case catalog in one exhaustive match so omissions remain compiler-visible"
)]
fn decision_detail_execution(
    kind: BaseE2eCaseKindV1,
    harness: Option<&BaseE2eHarnessIdentityV1>,
) -> BaseE2eDetailExecutionV1 {
    let mut builder = BaseE2eDetailPairBuilderV1::new(kind);
    match kind {
        BaseE2eCaseKindV1::CatalogLiterals
        | BaseE2eCaseKindV1::CanonicalRational
        | BaseE2eCaseKindV1::LogicalPath
        | BaseE2eCaseKindV1::PublicationSelection
        | BaseE2eCaseKindV1::CapabilityLeastPrivilege
        | BaseE2eCaseKindV1::Diagnostic
        | BaseE2eCaseKindV1::AtomicResult
        | BaseE2eCaseKindV1::PublicationStorage
        | BaseE2eCaseKindV1::CommandList => {}
        BaseE2eCaseKindV1::UnknownCatalogCode => {
            let expected = BaseE2eDetailPayloadV1::UnknownCatalog {
                catalog: "ProofExitV2",
                code: 1,
            };
            let observed = harness.map(|_| match ProofExitV2::from_code(1) {
                Err(error) => (
                    BaseE2eExpectedDecisionV1::Refuse,
                    BaseE2eDetailPayloadV1::UnknownCatalog {
                        catalog: error.catalog(),
                        code: error.code(),
                    },
                ),
                Ok(_) => (
                    BaseE2eExpectedDecisionV1::Accept,
                    BaseE2eDetailPayloadV1::AcceptedInstead,
                ),
            });
            push_detail_observation(
                &mut builder,
                1,
                "catalog.unknown-code".to_owned(),
                BaseE2eExpectedDecisionV1::Refuse,
                expected,
                observed,
            );
        }
        BaseE2eCaseKindV1::OverlongStableToken => {
            let expected = BaseE2eDetailPayloadV1::Value(ValueError::StableTokenTooLong {
                observed: 129,
                maximum: 128,
            });
            let observed = harness.map(|_| match StableTokenV2::new("a".repeat(129)) {
                Err(error) => (
                    BaseE2eExpectedDecisionV1::Refuse,
                    BaseE2eDetailPayloadV1::Value(error),
                ),
                Ok(_) => (
                    BaseE2eExpectedDecisionV1::Accept,
                    BaseE2eDetailPayloadV1::AcceptedInstead,
                ),
            });
            push_detail_observation(
                &mut builder,
                1,
                "value.overlong-token".to_owned(),
                BaseE2eExpectedDecisionV1::Refuse,
                expected,
                observed,
            );
        }
        BaseE2eCaseKindV1::ReservedContentStorePrefix => {
            let expected = BaseE2eDetailPayloadV1::Path(PathError::ReservedContentStorePrefix {
                prefix: "__runner_",
            });
            let observed =
                harness.map(
                    |_| match ContentStoreObjectKeyV1::new("__runner_private/object") {
                        Err(error) => (
                            BaseE2eExpectedDecisionV1::Refuse,
                            BaseE2eDetailPayloadV1::Path(error),
                        ),
                        Ok(_) => (
                            BaseE2eExpectedDecisionV1::Accept,
                            BaseE2eDetailPayloadV1::AcceptedInstead,
                        ),
                    },
                );
            push_detail_observation(
                &mut builder,
                1,
                "path.reserved-prefix".to_owned(),
                BaseE2eExpectedDecisionV1::Refuse,
                expected,
                observed,
            );
        }
        BaseE2eCaseKindV1::WindowsUnicodeAlias => {
            let expected = BaseE2eDetailPayloadV1::PathAdjudication(
                BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(
                    "résumé/a".to_owned(),
                ),
            );
            let observed = harness.map(|_| {
                let paths =
                    [LogicalBundlePathV1::new("résumé/a").expect("valid UTF-8 path fixture")];
                let detail = normalize_path_adjudication(adjudicate_logical_bundle_path_set(
                    PlatformPathProfileV2::WindowsHandleRelativeV1,
                    &paths,
                ));
                (
                    if matches!(
                        detail,
                        BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(_)
                    ) {
                        BaseE2eExpectedDecisionV1::Unsupported
                    } else {
                        BaseE2eExpectedDecisionV1::Refuse
                    },
                    BaseE2eDetailPayloadV1::PathAdjudication(detail),
                )
            });
            push_detail_observation(
                &mut builder,
                1,
                "path.windows-unicode-alias".to_owned(),
                BaseE2eExpectedDecisionV1::Unsupported,
                expected,
                observed,
            );
        }
        BaseE2eCaseKindV1::LimitCatalog => {
            append_limit_detail_cells(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::BudgetAdmission => {
            append_budget_admission_detail_cells(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::BudgetChildRelation => {
            append_budget_child_detail(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::PublicationCrossCell => {
            append_publication_cross_detail(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::CapabilityExtraRight => {
            append_capability_detail_cells(&mut builder, harness);
        }
        BaseE2eCaseKindV1::StatePass => {
            append_state_detail_cells(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::StateUsageInLifecycle => {
            append_state_usage_detail(&mut builder, harness.is_some());
        }
        BaseE2eCaseKindV1::DiagnosticRankGap => {
            append_diagnostic_rank_detail(&mut builder, harness);
        }
        BaseE2eCaseKindV1::IdentityMutation => {
            append_identity_detail_cells(&mut builder, harness);
        }
        BaseE2eCaseKindV1::NoClaimNominality => {
            append_no_claim_detail_cells(&mut builder, harness);
        }
        BaseE2eCaseKindV1::AtomicResultPresence => {
            append_atomic_presence_detail(&mut builder, harness.is_some());
        }
    }
    builder.finish()
}

fn append_limit_detail_cells(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    for profile in [RunProfileV2::Smoke, RunProfileV2::Full] {
        let (profile_name, profile_semantic_offset) = match profile {
            RunProfileV2::Smoke => ("smoke", 0_u32),
            RunProfileV2::Full => ("full", 142_u32),
        };
        for &(field, ordinal, name, unit, _width, tightenability, smoke_value, full_value) in
            limit_oracle_rows()
        {
            let base_value = match profile {
                RunProfileV2::Smoke => smoke_value,
                RunProfileV2::Full => full_value,
            };
            let observed_value = match base_value {
                RunnerLimitValueV2::U32(value) => RunnerLimitValueV2::U32(
                    value
                        .checked_add(1)
                        .expect("every frozen u32 limit oracle value is below u32::MAX"),
                ),
                RunnerLimitValueV2::U64(value) => RunnerLimitValueV2::U64(
                    value
                        .checked_add(1)
                        .expect("every frozen u64 limit oracle value is below u64::MAX"),
                ),
            };
            let expected_kind = match tightenability {
                RunnerLimitTightenabilityV2::Fixed => {
                    RunnerLimitsViolationKindV2::FixedFieldChanged
                }
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling
                }
            };
            let expected_expectation = match tightenability {
                RunnerLimitTightenabilityV2::Fixed => RunnerLimitExpectationV2::Exactly(base_value),
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitExpectationV2::AtMost(base_value)
                }
            };
            let expected = BaseE2eDetailPayloadV1::Limit {
                kind: expected_kind,
                field,
                unit,
                expected: expected_expectation,
                observed: observed_value,
                owner: "fs-evidence-runner.runner-limits",
                repair_rank: 1,
                repair_kind: match expected_expectation {
                    RunnerLimitExpectationV2::AtMost(_) => RepairActionKindV2::ReduceResourceDemand,
                    RunnerLimitExpectationV2::Width(_)
                    | RunnerLimitExpectationV2::AtLeast(_)
                    | RunnerLimitExpectationV2::Exactly(_)
                    | RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
                        RepairActionKindV2::UpdatePolicyOrCapability
                    }
                },
                repair_target: name,
            };
            let observed = observe.then(|| {
                let observed_field = RunnerLimitFieldV2::from_ordinal(ordinal)
                    .expect("observed schema resolves every literal oracle ordinal");
                let mut one_over = RunnerLimitsV2::base(profile).to_candidate();
                one_over
                    .set_value(observed_field, observed_value)
                    .expect("one-over fixture preserves the observed field width");
                match RunnerLimitsV2::admit_family(
                    profile,
                    one_over,
                    RunnerFamilyLimitRequirementsV2::NONE,
                ) {
                    Err(error) => (
                        BaseE2eExpectedDecisionV1::Refuse,
                        normalize_limit_violation(error),
                    ),
                    Ok(_) => (
                        BaseE2eExpectedDecisionV1::Accept,
                        BaseE2eDetailPayloadV1::AcceptedInstead,
                    ),
                }
            });
            push_detail_observation(
                builder,
                profile_semantic_offset + u32::from(ordinal) * 2,
                format!("limit.{profile_name}.{name}.one-over"),
                BaseE2eExpectedDecisionV1::Refuse,
                expected,
                observed,
            );
        }
    }
}

fn append_budget_admission_detail_cells(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    for &(
        semantic_ordinal,
        profile,
        _profile_code,
        profile_name,
        field,
        _field_ordinal,
        field_name,
        unit,
        ceiling,
        one_over,
    ) in budget_profile_refusal_oracle_rows()
    {
        let expected = BaseE2eDetailPayloadV1::Budget {
            kind: RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            field,
            unit,
            expected: RunnerBudgetExpectationV2::AtMost(ceiling),
            observed: one_over,
            owner: "fs-evidence-runner.runner-budgets",
            repair_rank: 1,
            repair_kind: RepairActionKindV2::ReduceResourceDemand,
            repair_target: field_name,
        };
        let observed = observe.then(|| observe_budget_profile_refusal(profile, field, one_over));
        push_detail_observation(
            builder,
            semantic_ordinal,
            format!("budget.profile.{profile_name}.{field_name}.one-over"),
            BaseE2eExpectedDecisionV1::Refuse,
            expected,
            observed,
        );
    }
}

fn observe_budget_profile_refusal(
    profile: RunProfileV2,
    field: RunnerBudgetFieldV2,
    one_over: RunnerBudgetValueV2,
) -> (BaseE2eExpectedDecisionV1, BaseE2eDetailPayloadV1) {
    let mut candidate = durable_budget_candidate();
    match (field, one_over) {
        (RunnerBudgetFieldV2::WallTimeNs, RunnerBudgetValueV2::U64(value)) => {
            candidate.wall_time_ns = value;
        }
        (RunnerBudgetFieldV2::MaxResidentBytes, RunnerBudgetValueV2::U64(value)) => {
            candidate.max_resident_bytes = value;
        }
        (RunnerBudgetFieldV2::MaxParallelChildren, RunnerBudgetValueV2::U32(value)) => {
            let total_children = budget_profile_oracle_rows()
                .iter()
                .find(|&&(candidate_profile, ..)| candidate_profile == profile)
                .map(|row| row.6)
                .expect("every observed profile has a literal ceiling row");
            candidate.max_child_processes = total_children;
            candidate.max_parallel_children = value;
        }
        (RunnerBudgetFieldV2::MaxChildProcesses, RunnerBudgetValueV2::U32(value)) => {
            candidate.max_child_processes = value;
        }
        _ => {
            return (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            );
        }
    }
    match RunnerBudgetsV2::try_new(candidate).and_then(|budgets| {
        budgets.admit(
            profile,
            ArtifactDispositionV2::DurableBundleRequired,
            &RunnerLimitsV2::base(profile),
        )
    }) {
        Err(error) => (
            BaseE2eExpectedDecisionV1::Refuse,
            normalize_budget_violation(error),
        ),
        Ok(_) => (
            BaseE2eExpectedDecisionV1::Accept,
            BaseE2eDetailPayloadV1::AcceptedInstead,
        ),
    }
}

fn append_budget_child_detail(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    let expected = BaseE2eDetailPayloadV1::Budget {
        kind: RunnerBudgetViolationKindV2::ParallelChildrenExceedTotal,
        field: RunnerBudgetFieldV2::MaxParallelChildren,
        unit: RunnerBudgetUnitV2::Count,
        expected: RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(8)),
        observed: RunnerBudgetValueV2::U32(9),
        owner: "fs-evidence-runner.runner-budgets",
        repair_rank: 1,
        repair_kind: RepairActionKindV2::ReduceResourceDemand,
        repair_target: "max_parallel_children",
    };
    let observed = observe.then(|| {
        let mut candidate = durable_budget_candidate();
        candidate.max_parallel_children = candidate.max_child_processes + 1;
        match RunnerBudgetsV2::try_new(candidate) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                normalize_budget_violation(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        1,
        "budget.parallel-children".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        expected,
        observed,
    );
}

fn publication_cross_error() -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Incompatible,
        "publication.profile_protocol_target",
        "the exact frozen profile/protocol/target cell",
        ConstructionObservedV2::closed_triple(
            &PlatformPathProfileV2::WindowsHandleRelativeV1,
            &PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            &PlatformPathProfileV2::WindowsHandleRelativeV1,
        ),
    )
}

fn append_publication_cross_detail(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    let observed = observe.then(|| {
        let path = LogicalBundlePathV1::new("runner/seal").expect("fixture path");
        match PublicationSelectionV2::new(
            PlatformPathProfileV2::WindowsHandleRelativeV1,
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            DestinationAdmissionModeV2::Absent,
            PublicationTargetV2::WindowsRelative(path),
        ) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Construction(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        1,
        "publication.cross-cell".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Construction(publication_cross_error()),
        observed,
    );
}

fn capability_oracle_rights(
    profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    mode: DestinationAdmissionModeV2,
) -> Vec<RootCapabilityRightV2> {
    use RootCapabilityRightV2 as Right;
    match (profile, access, mode) {
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            _,
        ) => vec![Right::Traverse, Right::ReadObject, Right::Enumerate],
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
        ) => vec![
            Right::Traverse,
            Right::Enumerate,
            Right::CreateObject,
            Right::SyncObject,
            Right::SyncContainer,
        ],
        (
            PlatformPathProfileV2::PosixDescriptorRelativeV1
            | PlatformPathProfileV2::WindowsHandleRelativeV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ) => vec![
            Right::Traverse,
            Right::Enumerate,
            Right::CreateObject,
            Right::PopulateEmptyDestination,
            Right::SyncObject,
            Right::SyncContainer,
        ],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::ReadOnlyInput,
            _,
        ) => vec![Right::ReadObject, Right::Enumerate, Right::QueryGeneration],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::Absent,
        ) => vec![
            Right::CreateObject,
            Right::QueryGeneration,
            Right::CommitCompareAndSwap,
        ],
        (
            PlatformPathProfileV2::ContentStoreObjectKeyV1,
            RootCapabilityAccessV2::DurableOutput,
            DestinationAdmissionModeV2::PreExistingEmpty,
        ) => vec![
            Right::Enumerate,
            Right::CreateObject,
            Right::AcquireExclusiveLease,
            Right::QueryGeneration,
            Right::CommitCompareAndSwap,
        ],
    }
}

fn expected_capability_refusal(
    profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    rights: &[RootCapabilityRightV2],
) -> (BaseE2eCapabilityRefusalStageV1, ConstructionErrorV2) {
    let intrinsically_legal = DestinationAdmissionModeV2::ALL
        .into_iter()
        .any(|mode| capability_oracle_rights(profile, access, mode).as_slice() == rights);
    let stage = if intrinsically_legal {
        BaseE2eCapabilityRefusalStageV1::ContextualNarrowing
    } else {
        BaseE2eCapabilityRefusalStageV1::IntrinsicPolicy
    };
    let expected = match stage {
        BaseE2eCapabilityRefusalStageV1::IntrinsicPolicy => {
            "the exact rights of at least one legal profile/access/mode cell"
        }
        BaseE2eCapabilityRefusalStageV1::ContextualNarrowing => match access {
            RootCapabilityAccessV2::ReadOnlyInput => "the exact least-privilege input cell",
            RootCapabilityAccessV2::DurableOutput => "the exact least-privilege output cell",
        },
    };
    (
        stage,
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "root_capability_policy.rights",
            expected,
            ConstructionObservedV2::closed_sequence(rights.iter()),
        ),
    )
}

fn observed_capability_refusal(
    profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    mode: DestinationAdmissionModeV2,
    rights: Vec<RootCapabilityRightV2>,
    no_claim_scope: NoClaimScopeRootV1,
) -> (BaseE2eExpectedDecisionV1, BaseE2eDetailPayloadV1) {
    let policy = match RootCapabilityPolicyV2::new(
        root_class_for_access(access),
        profile,
        access,
        rights,
        1,
        1,
        1,
        no_claim_scope,
    ) {
        Ok(policy) => policy,
        Err(error) => {
            return (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Capability {
                    stage: BaseE2eCapabilityRefusalStageV1::IntrinsicPolicy,
                    error,
                },
            );
        }
    };
    let result = match access {
        RootCapabilityAccessV2::ReadOnlyInput => {
            NarrowedPolicyViewV2::for_read_only(&policy).map(|_| ())
        }
        RootCapabilityAccessV2::DurableOutput => selection_for_profile(profile, mode)
            .and_then(|selection| NarrowedPolicyViewV2::for_publication(&policy, &selection))
            .map(|_| ()),
    };
    match result {
        Err(error) => (
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::Capability {
                stage: BaseE2eCapabilityRefusalStageV1::ContextualNarrowing,
                error,
            },
        ),
        Ok(()) => (
            BaseE2eExpectedDecisionV1::Accept,
            BaseE2eDetailPayloadV1::AcceptedInstead,
        ),
    }
}

fn append_capability_detail_cells(
    builder: &mut BaseE2eDetailPairBuilderV1,
    harness: Option<&BaseE2eHarnessIdentityV1>,
) {
    let mut semantic_ordinal = 0_u32;
    for profile in PlatformPathProfileV2::ALL {
        for access in RootCapabilityAccessV2::ALL {
            for mode in DestinationAdmissionModeV2::ALL {
                let exact = capability_oracle_rights(profile, access, mode);
                for right in RootCapabilityRightV2::ALL {
                    let mut mutant = exact.clone();
                    if let Some(index) = mutant.iter().position(|candidate| *candidate == right) {
                        mutant.remove(index);
                    } else {
                        mutant.push(right);
                    }
                    mutant.sort_unstable_by_key(|candidate| candidate.code());
                    semantic_ordinal += 1;
                    let (stage, expected_error) =
                        expected_capability_refusal(profile, access, &mutant);
                    let observed = harness.map(|harness| {
                        observed_capability_refusal(
                            profile,
                            access,
                            mode,
                            mutant.clone(),
                            harness.no_claim_scope.clone(),
                        )
                    });
                    push_detail_observation(
                        builder,
                        semantic_ordinal,
                        format!(
                            "capability.mutant.{}.{}.{}.{}",
                            profile.code(),
                            access.code(),
                            mode.code(),
                            right.code()
                        ),
                        BaseE2eExpectedDecisionV1::Refuse,
                        BaseE2eDetailPayloadV1::Capability {
                            stage,
                            error: expected_error,
                        },
                        observed,
                    );
                }
                for removed in exact.iter().copied() {
                    for replacement in RootCapabilityRightV2::ALL
                        .into_iter()
                        .filter(|candidate| !exact.contains(candidate))
                    {
                        let mut mutant = exact.clone();
                        let index = mutant
                            .iter()
                            .position(|candidate| *candidate == removed)
                            .expect("removed right is in the exact literal oracle cell");
                        mutant[index] = replacement;
                        mutant.sort_unstable_by_key(|candidate| candidate.code());
                        semantic_ordinal += 1;
                        let (stage, expected_error) =
                            expected_capability_refusal(profile, access, &mutant);
                        let observed = harness.map(|harness| {
                            observed_capability_refusal(
                                profile,
                                access,
                                mode,
                                mutant.clone(),
                                harness.no_claim_scope.clone(),
                            )
                        });
                        push_detail_observation(
                            builder,
                            semantic_ordinal,
                            format!(
                                "capability.substitution.{}.{}.{}.{}.{}",
                                profile.code(),
                                access.code(),
                                mode.code(),
                                removed.code(),
                                replacement.code()
                            ),
                            BaseE2eExpectedDecisionV1::Refuse,
                            BaseE2eDetailPayloadV1::Capability {
                                stage,
                                error: expected_error,
                            },
                            observed,
                        );
                    }
                }
            }
        }
    }
    debug_assert_eq!(semantic_ordinal, 390);
}

fn expected_state_validation(
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
    reason: Option<RefusedReasonV2>,
    diagnostic: Option<DiagnosticCodeV2>,
    drain: Option<PresentedDrainRootKindV2>,
) -> Result<(), StateValidationErrorV2> {
    let state_allowed = match role {
        StateBearingRecordRoleV2::PreRunDiagnostic => matches!(
            state,
            ProofExitV2::Usage
                | ProofExitV2::Refused
                | ProofExitV2::NoData
                | ProofExitV2::Stale
                | ProofExitV2::EnvironmentInvalid
                | ProofExitV2::Blocked
                | ProofExitV2::Unsupported
                | ProofExitV2::Cancelled
                | ProofExitV2::TimedOut
                | ProofExitV2::InternalError
        ),
        StateBearingRecordRoleV2::ExecutedCaseTerminal | StateBearingRecordRoleV2::RunTerminal => {
            !matches!(state, ProofExitV2::Usage | ProofExitV2::NotRun)
        }
        StateBearingRecordRoleV2::SuppressedCaseTerminal => state == ProofExitV2::NotRun,
    };
    if !state_allowed {
        return Err(StateValidationErrorV2::StateNotAllowedForRole { role, state });
    }
    match (state, reason) {
        (ProofExitV2::Refused, None) => {
            return Err(StateValidationErrorV2::MissingRefusedReason);
        }
        (ProofExitV2::Refused, Some(_)) | (_, None) => {}
        (_, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedRefusedReason { state, observed });
        }
    }
    match (expected_diagnostic(state), diagnostic) {
        (None, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedDiagnostic { observed });
        }
        (Some(expected), None) => {
            return Err(StateValidationErrorV2::MissingDiagnostic { expected });
        }
        (Some(expected), Some(observed)) if expected != observed => {
            return Err(StateValidationErrorV2::WrongDiagnostic { expected, observed });
        }
        (None, None) | (Some(_), Some(_)) => {}
    }
    match (expected_drain(role, state), drain) {
        (None, Some(observed)) => {
            return Err(StateValidationErrorV2::UnexpectedDrainBasis { state, observed });
        }
        (Some(expected), None) => {
            return Err(StateValidationErrorV2::MissingDrainBasis { expected });
        }
        (Some(expected), Some(observed)) if expected != observed => {
            return Err(StateValidationErrorV2::WrongDrainBasis { expected, observed });
        }
        (None, None) | (Some(_), Some(_)) => {}
    }
    Ok(())
}

#[allow(
    clippy::too_many_lines,
    reason = "the state-detail builder exhaustively enumerates the complete lifecycle cross-product in canonical semantic order"
)]
fn append_state_detail_cells(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    let reasons = core::iter::once(None)
        .chain(RefusedReasonV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let diagnostics = core::iter::once(None)
        .chain(DiagnosticCodeV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let drains = [
        None,
        Some(PresentedDrainRootKindV2::CancelledStopRoot),
        Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
    ];
    let mut semantic_ordinal = 0_u32;
    for role in StateBearingRecordRoleV2::ALL {
        for state in ProofExitV2::ALL {
            for reason in &reasons {
                for diagnostic in &diagnostics {
                    for drain in drains {
                        semantic_ordinal += 1;
                        let expected =
                            expected_state_validation(role, state, *reason, *diagnostic, drain);
                        let Err(expected_error) = expected else {
                            continue;
                        };
                        let observed = observe.then(|| {
                            match validate_state_v2(StateValidationInputV2::new(
                                role,
                                state,
                                *reason,
                                *diagnostic,
                                drain,
                            )) {
                                Err(error) => (
                                    BaseE2eExpectedDecisionV1::Refuse,
                                    BaseE2eDetailPayloadV1::State(error),
                                ),
                                Ok(_) => (
                                    BaseE2eExpectedDecisionV1::Accept,
                                    BaseE2eDetailPayloadV1::AcceptedInstead,
                                ),
                            }
                        });
                        push_detail_observation(
                            builder,
                            semantic_ordinal,
                            format!(
                                "state.matrix.{}.{}.{}.{}.{}",
                                role.code(),
                                state.code(),
                                reason.map_or(0, RefusedReasonV2::code),
                                diagnostic.map_or(0, DiagnosticCodeV2::code),
                                drain.map_or(0, drain_code)
                            ),
                            BaseE2eExpectedDecisionV1::Refuse,
                            BaseE2eDetailPayloadV1::State(expected_error),
                            observed,
                        );
                    }
                }
            }
        }
    }

    let (cancelled, timed_out, internal_error) =
        presented_stop_fixture_roots().expect("nominal stop-root fixtures");
    let causes = [
        NotRunCauseV2::PriorCancelled(cancelled),
        NotRunCauseV2::PriorTimedOut(timed_out),
        NotRunCauseV2::PriorControlledInternalError(internal_error),
    ];
    for cause in causes {
        let code = cause.code();
        semantic_ordinal += 2;
        semantic_ordinal += 1;
        let observed = observe.then(|| match NotRunBasisV2::new(cause.clone(), 256, 256) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::NotRun(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        });
        push_detail_observation(
            builder,
            semantic_ordinal,
            format!("not-run.{code}.one-over"),
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::NotRun(NotRunBasisErrorV2::LowestRemainingOrdinalOutOfRange {
                observed: 256,
                ordered_case_count: 256,
            }),
            observed,
        );

        semantic_ordinal += 1;
        let observed = observe.then(|| match NotRunBasisV2::new(cause, 0, 0) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::NotRun(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        });
        push_detail_observation(
            builder,
            semantic_ordinal,
            format!("not-run.{code}.empty"),
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::NotRun(NotRunBasisErrorV2::EmptyManifest),
            observed,
        );
    }
    debug_assert_eq!(semantic_ordinal, 32_460);
}

fn append_state_usage_detail(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    let expected = StateValidationErrorV2::StateNotAllowedForRole {
        role: StateBearingRecordRoleV2::ExecutedCaseTerminal,
        state: ProofExitV2::Usage,
    };
    let observed = observe.then(|| {
        match validate_state_v2(StateValidationInputV2::new(
            StateBearingRecordRoleV2::ExecutedCaseTerminal,
            ProofExitV2::Usage,
            None,
            Some(DiagnosticCodeV2::RunnerUsage),
            None,
        )) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::State(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        1,
        "state.usage-in-lifecycle".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::State(expected),
        observed,
    );
}

fn diagnostic_rank_error() -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::OutOfOrder,
        "diagnostic.repair_rank",
        "contiguous ranks beginning at one",
        2,
    )
}

fn append_diagnostic_rank_detail(
    builder: &mut BaseE2eDetailPairBuilderV1,
    harness: Option<&BaseE2eHarnessIdentityV1>,
) {
    let observed = harness.map(
        |harness| match diagnostic(harness.no_claim_scope.clone(), 2) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Construction(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        },
    );
    push_detail_observation(
        builder,
        1,
        "diagnostic.rank-gap".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Construction(diagnostic_rank_error()),
        observed,
    );
}

#[allow(
    clippy::too_many_lines,
    reason = "the identity-detail builder keeps all three nominal identity mutation matrices adjacent and identically ordered"
)]
fn append_identity_detail_cells(
    builder: &mut BaseE2eDetailPairBuilderV1,
    harness: Option<&BaseE2eHarnessIdentityV1>,
) {
    let mut semantic_ordinal = 0_u32;
    macro_rules! append_identity {
        (
            $type:ty,
            $label:literal,
            $expected_role:expr,
            $wrong_role:expr,
            $expected_domain:literal
        ) => {{
            semantic_ordinal += 32;
            semantic_ordinal += 1;
            let expected = IdentityError::WrongRole {
                expected: $expected_role,
                observed: $wrong_role,
            };
            let observed = harness.map(|_| {
                match <$type>::parse_presented($wrong_role, $expected_domain, &"00".repeat(32)) {
                    Err(error) => (
                        BaseE2eExpectedDecisionV1::Refuse,
                        BaseE2eDetailPayloadV1::Identity(error),
                    ),
                    Ok(_) => (
                        BaseE2eExpectedDecisionV1::Accept,
                        BaseE2eDetailPayloadV1::AcceptedInstead,
                    ),
                }
            });
            push_detail_observation(
                builder,
                semantic_ordinal,
                format!("identity.{}.wrong-role", $label),
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(expected),
                observed,
            );

            semantic_ordinal += 1;
            let wrong_domain = "org.frankensim.fs-evidence-runner.wrong-domain.v1";
            let expected = IdentityError::WrongDomain {
                expected: $expected_domain,
                observed: wrong_domain.to_owned(),
            };
            let observed = harness.map(|_| {
                match <$type>::parse_presented($expected_role, wrong_domain, &"00".repeat(32)) {
                    Err(error) => (
                        BaseE2eExpectedDecisionV1::Refuse,
                        BaseE2eDetailPayloadV1::Identity(error),
                    ),
                    Ok(_) => (
                        BaseE2eExpectedDecisionV1::Accept,
                        BaseE2eDetailPayloadV1::AcceptedInstead,
                    ),
                }
            });
            push_detail_observation(
                builder,
                semantic_ordinal,
                format!("identity.{}.wrong-domain", $label),
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(expected),
                observed,
            );

            semantic_ordinal += 1;
            let expected = IdentityError::WrongLowerHexLength {
                observed: 2,
                expected: 64,
            };
            let observed = harness.map(|_| {
                match <$type>::parse_presented($expected_role, $expected_domain, "00") {
                    Err(error) => (
                        BaseE2eExpectedDecisionV1::Refuse,
                        BaseE2eDetailPayloadV1::Identity(error),
                    ),
                    Ok(_) => (
                        BaseE2eExpectedDecisionV1::Accept,
                        BaseE2eDetailPayloadV1::AcceptedInstead,
                    ),
                }
            });
            push_detail_observation(
                builder,
                semantic_ordinal,
                format!("identity.{}.wrong-length", $label),
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(expected),
                observed,
            );
        }};
    }
    append_identity!(
        SourceIdentityRootV2,
        "source",
        DigestRoleV2::Source,
        DigestRoleV2::Build,
        "org.frankensim.fs-evidence-runner.source-identity.v1"
    );
    append_identity!(
        BuildIdentityRootV2,
        "build",
        DigestRoleV2::Build,
        DigestRoleV2::Toolchain,
        "org.frankensim.fs-evidence-runner.build-identity.v1"
    );
    append_identity!(
        ToolchainIdentityRootV2,
        "toolchain",
        DigestRoleV2::Toolchain,
        DigestRoleV2::Source,
        "org.frankensim.fs-evidence-runner.toolchain-identity.v1"
    );
    debug_assert_eq!(semantic_ordinal, 105);
}

fn append_no_claim_detail_cells(
    builder: &mut BaseE2eDetailPairBuilderV1,
    harness: Option<&BaseE2eHarnessIdentityV1>,
) {
    const DOMAIN: &str = "org.frankensim.fs-evidence-runner.no-claim-scope.v1";
    let hex = "00".repeat(32);
    let mut semantic_ordinal = 2_u32;

    semantic_ordinal += 1;
    let observed = harness.map(|_| {
        match NoClaimScopeRootV1::parse_presented(DigestRoleV2::Policy, DOMAIN, &hex) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        semantic_ordinal,
        "no-claim.wrong-role".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Identity(IdentityError::WrongRole {
            expected: DigestRoleV2::ClaimScope,
            observed: DigestRoleV2::Policy,
        }),
        observed,
    );

    semantic_ordinal += 1;
    let wrong_domain = "org.frankensim.fs-evidence-runner.wrong-no-claim.v1";
    let observed = harness.map(|_| {
        match NoClaimScopeRootV1::parse_presented(DigestRoleV2::ClaimScope, wrong_domain, &hex) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        semantic_ordinal,
        "no-claim.wrong-domain".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Identity(IdentityError::WrongDomain {
            expected: DOMAIN,
            observed: wrong_domain.to_owned(),
        }),
        observed,
    );

    semantic_ordinal += 1;
    let observed = harness.map(|_| {
        match NoClaimScopeRootV1::parse_presented(DigestRoleV2::ClaimScope, DOMAIN, "00") {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Identity(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        semantic_ordinal,
        "no-claim.wrong-length".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Identity(IdentityError::WrongLowerHexLength {
            observed: 2,
            expected: 64,
        }),
        observed,
    );
    debug_assert_eq!(semantic_ordinal, 5);
}

fn atomic_presence_error() -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Unexpected,
        "result.catalog",
        "zero bytes for the command result variant",
        1,
    )
}

fn append_atomic_presence_detail(builder: &mut BaseE2eDetailPairBuilderV1, observe: bool) {
    let observed = observe.then(|| {
        match SymbolicCommandResultPlanV2::new(RunnerCommandV2::Run, 32, 128, 1, 128, 1024, 1024) {
            Err(error) => (
                BaseE2eExpectedDecisionV1::Refuse,
                BaseE2eDetailPayloadV1::Construction(error),
            ),
            Ok(_) => (
                BaseE2eExpectedDecisionV1::Accept,
                BaseE2eDetailPayloadV1::AcceptedInstead,
            ),
        }
    });
    push_detail_observation(
        builder,
        1,
        "result.atomic-presence".to_owned(),
        BaseE2eExpectedDecisionV1::Refuse,
        BaseE2eDetailPayloadV1::Construction(atomic_presence_error()),
        observed,
    );
}

#[allow(
    clippy::too_many_lines,
    reason = "the case executor is the exhaustive compiler-checked dispatch over the closed base E2E case catalog"
)]
fn execute_case(
    kind: BaseE2eCaseKindV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> BaseE2eCaseExecutionV1 {
    let detail = execute_detail_manifest(kind, harness);
    match kind {
        BaseE2eCaseKindV1::CatalogLiterals => aggregate_accept(catalog_literal_matrix(), detail),
        BaseE2eCaseKindV1::UnknownCatalogCode => refuse_if(
            ProofExitV2::from_code(1)
                .is_err_and(|error| error.catalog() == "ProofExitV2" && error.code() == 1),
            "catalog.unknown-code",
            detail,
        ),
        BaseE2eCaseKindV1::CanonicalRational => accept_if(
            RationalV2::new(6, 8).ok() == RationalV2::new(3, 4).ok(),
            "value.rational-equivalence",
            detail,
        ),
        BaseE2eCaseKindV1::OverlongStableToken => refuse_if(
            StableTokenV2::new("a".repeat(129)).is_err_and(|error| {
                error
                    == ValueError::StableTokenTooLong {
                        observed: 129,
                        maximum: 128,
                    }
            }),
            "value.overlong-token",
            detail,
        ),
        BaseE2eCaseKindV1::LogicalPath => accept_if(
            LogicalBundlePathV1::new("runner/seal").is_ok(),
            "path.logical",
            detail,
        ),
        BaseE2eCaseKindV1::ReservedContentStorePrefix => refuse_if(
            ContentStoreObjectKeyV1::new("__runner_private/object").is_err_and(|error| {
                error
                    == PathError::ReservedContentStorePrefix {
                        prefix: "__runner_",
                    }
            }),
            "path.reserved-prefix",
            detail,
        ),
        BaseE2eCaseKindV1::WindowsUnicodeAlias => {
            let paths = [LogicalBundlePathV1::new("résumé/a").expect("valid UTF-8 path")];
            match adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &paths,
            ) {
                PathSetAdjudicationV1::UnsupportedWindowsNonAsciiAlias { .. } => {
                    BaseE2eCaseExecutionV1::unsupported(1, detail)
                }
                _ => BaseE2eCaseExecutionV1::with_failure(
                    BaseE2eExpectedDecisionV1::Refuse,
                    1,
                    "path.windows-unicode-alias",
                    detail,
                ),
            }
        }
        BaseE2eCaseKindV1::LimitCatalog => aggregate_mixed(
            mixed_progress_from_first_failure(limit_matrix(), 284, limit_matrix_partition),
            142,
            142,
            0,
            detail,
        ),
        BaseE2eCaseKindV1::BudgetAdmission => aggregate_mixed(
            mixed_progress_from_first_failure(budget_matrix(), 44, budget_matrix_partition),
            36,
            8,
            0,
            detail,
        ),
        BaseE2eCaseKindV1::BudgetChildRelation => {
            let mut candidate = durable_budget_candidate();
            candidate.max_parallel_children = candidate.max_child_processes + 1;
            refuse_if(
                RunnerBudgetsV2::try_new(candidate).is_err_and(|error| {
                    normalize_budget_violation(error)
                        == BaseE2eDetailPayloadV1::Budget {
                            kind: RunnerBudgetViolationKindV2::ParallelChildrenExceedTotal,
                            field: RunnerBudgetFieldV2::MaxParallelChildren,
                            unit: RunnerBudgetUnitV2::Count,
                            expected: RunnerBudgetExpectationV2::AtMost(RunnerBudgetValueV2::U32(
                                8,
                            )),
                            observed: RunnerBudgetValueV2::U32(9),
                            owner: "fs-evidence-runner.runner-budgets",
                            repair_rank: 1,
                            repair_kind: RepairActionKindV2::ReduceResourceDemand,
                            repair_target: "max_parallel_children",
                        }
                }),
                "budget.parallel-children",
                detail,
            )
        }
        BaseE2eCaseKindV1::PublicationSelection => {
            aggregate_accept(publication_selection_matrix(), detail)
        }
        BaseE2eCaseKindV1::PublicationCrossCell => {
            let path = LogicalBundlePathV1::new("runner/seal").expect("fixture path");
            refuse_if(
                PublicationSelectionV2::new(
                    PlatformPathProfileV2::WindowsHandleRelativeV1,
                    PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
                    DestinationAdmissionModeV2::Absent,
                    PublicationTargetV2::WindowsRelative(path),
                )
                .is_err_and(|error| error == publication_cross_error()),
                "publication.cross-cell",
                detail,
            )
        }
        BaseE2eCaseKindV1::CapabilityLeastPrivilege => {
            aggregate_accept(capability_valid_matrix(&harness.no_claim_scope), detail)
        }
        BaseE2eCaseKindV1::CapabilityExtraRight => {
            aggregate_refusal(capability_invalid_matrix(&harness.no_claim_scope), detail)
        }
        BaseE2eCaseKindV1::StatePass => {
            let partitions = state_matrix_partitions();
            aggregate_mixed(
                mixed_progress_from_first_failure(state_and_not_run_matrix(), 32_460, |ordinal| {
                    partitions
                        [usize::try_from(ordinal - 1).expect("state matrix ordinal fits usize")]
                }),
                69,
                32_391,
                0,
                detail,
            )
        }
        BaseE2eCaseKindV1::StateUsageInLifecycle => refuse_if(
            validate_state_v2(StateValidationInputV2::new(
                StateBearingRecordRoleV2::ExecutedCaseTerminal,
                ProofExitV2::Usage,
                None,
                Some(DiagnosticCodeV2::RunnerUsage),
                None,
            ))
            .is_err_and(|error| {
                error
                    == StateValidationErrorV2::StateNotAllowedForRole {
                        role: StateBearingRecordRoleV2::ExecutedCaseTerminal,
                        state: ProofExitV2::Usage,
                    }
            }),
            "state.usage-in-lifecycle",
            detail,
        ),
        BaseE2eCaseKindV1::Diagnostic => {
            aggregate_accept(diagnostic_matrix(&harness.no_claim_scope), detail)
        }
        BaseE2eCaseKindV1::DiagnosticRankGap => refuse_if(
            diagnostic(harness.no_claim_scope.clone(), 2)
                .is_err_and(|error| error == diagnostic_rank_error()),
            "diagnostic.rank-gap",
            detail,
        ),
        BaseE2eCaseKindV1::IdentityMutation => aggregate_mixed(
            mixed_progress_from_first_failure(
                identity_mutation_matrix(harness),
                105,
                identity_matrix_partition,
            ),
            96,
            9,
            0,
            detail,
        ),
        BaseE2eCaseKindV1::NoClaimNominality => aggregate_mixed(
            mixed_progress_from_first_failure(
                no_claim_matrix(&harness.no_claim_scope),
                5,
                no_claim_matrix_partition,
            ),
            2,
            3,
            0,
            detail,
        ),
        BaseE2eCaseKindV1::AtomicResult => accept_if(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::List, 32, 0, 128, 0, 1024, 1024)
                .is_ok(),
            "result.atomic",
            detail,
        ),
        BaseE2eCaseKindV1::AtomicResultPresence => refuse_if(
            SymbolicCommandResultPlanV2::new(RunnerCommandV2::Run, 32, 128, 1, 128, 1024, 1024)
                .is_err_and(|error| error == atomic_presence_error()),
            "result.atomic-presence",
            detail,
        ),
        BaseE2eCaseKindV1::PublicationStorage => {
            accept_if(publication_storage().is_ok(), "publication.storage", detail)
        }
        BaseE2eCaseKindV1::CommandList => aggregate_accept(command_matrix(), detail),
    }
}

fn detail_is_exact(detail: &BaseE2eDetailExecutionV1) -> bool {
    detail.expected == detail.observed
        && detail.matched_cells == detail.expected.cell_count
        && detail.first_divergent_cell.is_none()
}

fn accept_if(
    condition: bool,
    failed_cell: &'static str,
    detail: BaseE2eDetailExecutionV1,
) -> BaseE2eCaseExecutionV1 {
    if condition && detail_is_exact(&detail) {
        BaseE2eCaseExecutionV1::accepted(1, detail)
    } else {
        let detail_cell = detail.first_divergent_cell.clone();
        BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Refuse,
            1,
            detail_cell.as_deref().unwrap_or(failed_cell),
            detail,
        )
    }
}

fn refuse_if(
    condition: bool,
    failed_cell: &'static str,
    detail: BaseE2eDetailExecutionV1,
) -> BaseE2eCaseExecutionV1 {
    if condition && detail_is_exact(&detail) {
        BaseE2eCaseExecutionV1::refused(1, detail)
    } else {
        let detail_cell = detail.first_divergent_cell.clone();
        BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Accept,
            1,
            detail_cell.as_deref().unwrap_or(failed_cell),
            detail,
        )
    }
}

fn aggregate_accept(
    result: Result<u32, (u32, String)>,
    detail: BaseE2eDetailExecutionV1,
) -> BaseE2eCaseExecutionV1 {
    match result {
        Ok(checked_cells) if detail_is_exact(&detail) => {
            BaseE2eCaseExecutionV1::accepted(checked_cells, detail)
        }
        Ok(checked_cells) => {
            let detail_cell = detail.first_divergent_cell.clone();
            BaseE2eCaseExecutionV1::with_failure(
                BaseE2eExpectedDecisionV1::Refuse,
                checked_cells,
                detail_cell.as_deref().unwrap_or("detail.manifest"),
                detail,
            )
        }
        Err((checked_cells, failed_cell)) => BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Refuse,
            checked_cells,
            failed_cell,
            detail,
        ),
    }
}

fn aggregate_refusal(
    result: Result<u32, (u32, String)>,
    detail: BaseE2eDetailExecutionV1,
) -> BaseE2eCaseExecutionV1 {
    match result {
        Ok(checked_cells) if detail_is_exact(&detail) => {
            BaseE2eCaseExecutionV1::refused(checked_cells, detail)
        }
        Ok(checked_cells) => {
            let detail_cell = detail.first_divergent_cell.clone();
            BaseE2eCaseExecutionV1::with_failure(
                BaseE2eExpectedDecisionV1::Accept,
                checked_cells,
                detail_cell.as_deref().unwrap_or("detail.manifest"),
                detail,
            )
        }
        Err((checked_cells, failed_cell)) => BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Accept,
            checked_cells,
            failed_cell,
            detail,
        ),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BaseE2eMixedPartitionV1 {
    Positive,
    ExpectedRefusal,
    Unsupported,
}

impl ConstructionClosedSemanticV2 for BaseE2eMixedPartitionV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::ExpectedRefusal => "expected-refusal",
            Self::Unsupported => "unsupported",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eMixedProgressV1 {
    checked_cells: u32,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unsupported: u32,
    unexpected_mismatches: u32,
    first_failed_cell: Option<String>,
    first_failed_partition: Option<BaseE2eMixedPartitionV1>,
    last_partition: Option<BaseE2eMixedPartitionV1>,
}

impl BaseE2eMixedProgressV1 {
    const fn new() -> Self {
        Self {
            checked_cells: 0,
            positive_eligible: 0,
            positive_matched: 0,
            expected_refusals: 0,
            expected_refusals_matched: 0,
            unsupported: 0,
            unexpected_mismatches: 0,
            first_failed_cell: None,
            first_failed_partition: None,
            last_partition: None,
        }
    }

    fn record(
        &mut self,
        partition: BaseE2eMixedPartitionV1,
        matched: bool,
        stable_id: impl Into<String>,
    ) -> bool {
        self.checked_cells = self
            .checked_cells
            .checked_add(1)
            .expect("the frozen semantic matrices fit u32");
        self.last_partition = Some(partition);
        match partition {
            BaseE2eMixedPartitionV1::Positive => {
                self.positive_eligible += 1;
                self.positive_matched += u32::from(matched);
            }
            BaseE2eMixedPartitionV1::ExpectedRefusal => {
                self.expected_refusals += 1;
                self.expected_refusals_matched += u32::from(matched);
            }
            BaseE2eMixedPartitionV1::Unsupported => {
                self.unsupported += u32::from(matched);
            }
        }
        if !matched {
            self.unexpected_mismatches += 1;
            if self.first_failed_cell.is_none() {
                self.first_failed_cell = Some(stable_id.into());
                self.first_failed_partition = Some(partition);
            }
        }
        matched
    }

    fn invalidate_last(&mut self, stable_id: impl Into<String>) {
        if self.unexpected_mismatches > 0 {
            return;
        }
        match self
            .last_partition
            .expect("a total-count check follows at least one semantic cell")
        {
            BaseE2eMixedPartitionV1::Positive => self.positive_matched -= 1,
            BaseE2eMixedPartitionV1::ExpectedRefusal => self.expected_refusals_matched -= 1,
            BaseE2eMixedPartitionV1::Unsupported => self.unsupported -= 1,
        }
        self.unexpected_mismatches = 1;
        self.first_failed_cell = Some(stable_id.into());
        self.first_failed_partition = self.last_partition;
    }

    const fn is_green(&self) -> bool {
        self.positive_eligible == self.positive_matched
            && self.expected_refusals == self.expected_refusals_matched
            && self.unexpected_mismatches == 0
    }
}

fn mixed_progress_from_first_failure(
    result: Result<u32, (u32, String)>,
    expected_total: u32,
    partition_for_ordinal: impl Fn(u32) -> BaseE2eMixedPartitionV1,
) -> BaseE2eMixedProgressV1 {
    assert!(
        expected_total > 0,
        "a mixed matrix must declare at least one expected semantic cell"
    );
    let mut progress = BaseE2eMixedProgressV1::new();
    match result {
        Ok(checked) => {
            assert!(
                checked > 0,
                "a successful mixed matrix must report real observed semantic-cell progress"
            );
            for ordinal in 1..=checked {
                let partition = partition_for_ordinal(ordinal);
                progress.record(partition, true, "matched");
            }
            if checked != expected_total {
                progress.invalidate_last("matrix.partition-count");
            }
        }
        Err((checked, failed_cell)) => {
            assert!(
                checked > 0,
                "a mixed-matrix failure must identify a real observed semantic ordinal"
            );
            let failed_ordinal = checked;
            for ordinal in 1..=failed_ordinal {
                let partition = partition_for_ordinal(ordinal);
                if !progress.record(partition, ordinal != failed_ordinal, failed_cell.clone()) {
                    break;
                }
            }
        }
    }
    progress
}

const fn limit_matrix_partition(ordinal: u32) -> BaseE2eMixedPartitionV1 {
    if ordinal % 2 == 1 {
        BaseE2eMixedPartitionV1::Positive
    } else {
        BaseE2eMixedPartitionV1::ExpectedRefusal
    }
}

const fn budget_matrix_partition(ordinal: u32) -> BaseE2eMixedPartitionV1 {
    if ordinal <= 36 {
        BaseE2eMixedPartitionV1::Positive
    } else {
        BaseE2eMixedPartitionV1::ExpectedRefusal
    }
}

const fn identity_matrix_partition(ordinal: u32) -> BaseE2eMixedPartitionV1 {
    if (ordinal - 1) % 35 < 32 {
        BaseE2eMixedPartitionV1::Positive
    } else {
        BaseE2eMixedPartitionV1::ExpectedRefusal
    }
}

const fn no_claim_matrix_partition(ordinal: u32) -> BaseE2eMixedPartitionV1 {
    if ordinal <= 2 {
        BaseE2eMixedPartitionV1::Positive
    } else {
        BaseE2eMixedPartitionV1::ExpectedRefusal
    }
}

fn state_matrix_partitions() -> Vec<BaseE2eMixedPartitionV1> {
    let reasons = core::iter::once(None)
        .chain(RefusedReasonV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let diagnostics = core::iter::once(None)
        .chain(DiagnosticCodeV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let drains = [
        None,
        Some(PresentedDrainRootKindV2::CancelledStopRoot),
        Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
    ];
    let mut partitions = Vec::with_capacity(32_460);
    for role in StateBearingRecordRoleV2::ALL {
        for state in ProofExitV2::ALL {
            for reason in &reasons {
                for diagnostic in &diagnostics {
                    for drain in drains {
                        partitions.push(
                            if expected_state_cell(role, state, *reason, *diagnostic, drain) {
                                BaseE2eMixedPartitionV1::Positive
                            } else {
                                BaseE2eMixedPartitionV1::ExpectedRefusal
                            },
                        );
                    }
                }
            }
        }
    }
    for _ in 0..3 {
        partitions.extend([
            BaseE2eMixedPartitionV1::Positive,
            BaseE2eMixedPartitionV1::Positive,
            BaseE2eMixedPartitionV1::ExpectedRefusal,
            BaseE2eMixedPartitionV1::ExpectedRefusal,
        ]);
    }
    debug_assert_eq!(partitions.len(), 32_460);
    partitions
}

fn aggregate_mixed(
    mut progress: BaseE2eMixedProgressV1,
    positive_eligible: u32,
    expected_refusals: u32,
    unsupported: u32,
    detail: BaseE2eDetailExecutionV1,
) -> BaseE2eCaseExecutionV1 {
    let expected_checked_cells = positive_eligible
        .checked_add(expected_refusals)
        .and_then(|count| count.checked_add(unsupported))
        .expect("the frozen mixed semantic-cell inventory fits u32");
    assert!(
        expected_checked_cells > 0,
        "a mixed aggregate must declare at least one expected semantic cell"
    );
    assert!(
        progress.checked_cells > 0,
        "a mixed aggregate must retain real observed semantic-cell progress"
    );
    let totals_exact = progress.positive_eligible == positive_eligible
        && progress.expected_refusals == expected_refusals
        && progress.unsupported == unsupported
        && progress.checked_cells == expected_checked_cells;
    if progress.is_green() && totals_exact && detail_is_exact(&detail) {
        return BaseE2eCaseExecutionV1::mixed(
            positive_eligible,
            expected_refusals,
            unsupported,
            detail,
        );
    }
    if progress.is_green() && !totals_exact {
        progress.invalidate_last("matrix.partition-count");
    } else if progress.is_green() {
        let detail_cell = detail
            .first_divergent_cell
            .clone()
            .unwrap_or_else(|| "detail.manifest".to_owned());
        let expected_partition = detail
            .expected_cells
            .as_deref()
            .and_then(|cells| cells.iter().find(|cell| cell.stable_id() == detail_cell))
            .map_or(
                BaseE2eMixedPartitionV1::ExpectedRefusal,
                |cell| match cell.decision() {
                    BaseE2eExpectedDecisionV1::Accept => BaseE2eMixedPartitionV1::Positive,
                    BaseE2eExpectedDecisionV1::Refuse => BaseE2eMixedPartitionV1::ExpectedRefusal,
                    BaseE2eExpectedDecisionV1::Unsupported => BaseE2eMixedPartitionV1::Unsupported,
                },
            );
        match expected_partition {
            BaseE2eMixedPartitionV1::Positive if progress.positive_matched > 0 => {
                progress.positive_matched -= 1;
            }
            BaseE2eMixedPartitionV1::ExpectedRefusal if progress.expected_refusals_matched > 0 => {
                progress.expected_refusals_matched -= 1;
            }
            BaseE2eMixedPartitionV1::Unsupported if progress.unsupported > 0 => {
                progress.unsupported -= 1;
            }
            _ => progress.invalidate_last("detail.manifest"),
        }
        progress.unexpected_mismatches = 1;
        progress.first_failed_cell = Some(detail_cell);
        progress.first_failed_partition = Some(expected_partition);
    }

    assert!(
        progress.unexpected_mismatches > 0,
        "a red mixed aggregate must retain an observed mismatch"
    );
    let first_failed_partition = progress
        .first_failed_partition
        .expect("a red mixed aggregate must retain the failed semantic-cell partition");
    let first_failed_cell = progress
        .first_failed_cell
        .expect("a red mixed aggregate must retain the failed semantic-cell ID");
    BaseE2eCaseExecutionV1 {
        decision: match first_failed_partition {
            BaseE2eMixedPartitionV1::ExpectedRefusal => BaseE2eExpectedDecisionV1::Accept,
            BaseE2eMixedPartitionV1::Unsupported | BaseE2eMixedPartitionV1::Positive => {
                BaseE2eExpectedDecisionV1::Refuse
            }
        },
        checked_cells: progress.checked_cells,
        positive_eligible: progress.positive_eligible,
        positive_matched: progress.positive_matched,
        expected_refusals: progress.expected_refusals,
        expected_refusals_matched: progress.expected_refusals_matched,
        unsupported: progress.unsupported,
        unexpected_mismatches: progress.unexpected_mismatches,
        first_failed_cell: Some(first_failed_cell),
        detail,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SourceClosureExecutionReportV1 {
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unexpected_mismatches: u32,
    matched_cases: [bool; 15],
}

fn source_closure_refusal_matches(
    result: &Result<RunnerV2BaseSourceClosureV1, ConstructionErrorV2>,
    expected: ConstructionErrorV2,
) -> bool {
    result == &Err(expected)
}

#[allow(
    clippy::too_many_lines,
    reason = "the source-closure audit is an exhaustive ordered matrix whose fifteen mutations share one exact-input fixture"
)]
fn run_source_closure_checks(
    closure: &RunnerV2BaseSourceClosureV1,
) -> SourceClosureExecutionReportV1 {
    let exact_inputs = || {
        EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| BaseSourceClosureInputV1::presented(file.path, file.bytes.to_vec()))
            .collect::<Vec<_>>()
    };
    let mut matched_cases = [false; 15];

    let exact = exact_inputs();
    if RunnerV2BaseSourceClosureV1::reconstruct(&exact).as_ref() == Ok(closure) {
        matched_cases[0] = true;
    }

    let mut missing = exact_inputs();
    missing.pop();
    matched_cases[1] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&missing),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_source_closure.entries",
            "all exact embedded source entries",
            missing.len(),
        ),
    );

    let mut extra = exact_inputs();
    extra.push(BaseSourceClosureInputV1::presented(
        "crates/fs-evidence-runner/src/unowned-ambient.rs",
        b"unowned ambient source".to_vec(),
    ));
    matched_cases[2] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&extra),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_source_closure.entries",
            "no entries beyond the exact embedded source set",
            extra.len(),
        ),
    );

    let mut duplicate = exact_inputs();
    duplicate[1] = duplicate[0].clone();
    matched_cases[3] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&duplicate),
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Duplicate,
            "base_source_closure.path",
            "one unique path per exact source entry",
            ConstructionObservedDataClassV2::CallerControlledText,
        ),
    );

    let mut reordered = exact_inputs();
    reordered.swap(0, 1);
    matched_cases[4] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&reordered),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfOrder,
            "base_source_closure.path",
            "the exact bytewise-lexicographic source order",
            ConstructionObservedV2::indexed_redacted(
                0,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ),
    );

    let mut stale = exact_inputs();
    stale[0].bytes[0] ^= 1;
    stale[0].content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, &stale[0].bytes);
    matched_cases[5] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&stale),
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.bytes",
            "the exact compile-time included source bytes",
            ConstructionObservedDataClassV2::BulkPayload,
        ),
    );

    let mut owner = exact_inputs();
    owner[0].owner_code ^= u16::MAX;
    matched_cases[6] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&owner),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.owner",
            "the exact sole declaration owner",
            ConstructionObservedV2::unsigned_pair(0, u64::from(owner[0].owner_code())),
        ),
    );

    let mut source_route = exact_inputs();
    source_route[0].source_route_code ^= u16::MAX;
    matched_cases[7] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&source_route),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.source_route",
            "the exact closed compile-time source route",
            ConstructionObservedV2::unsigned_pair(
                0,
                u64::from(source_route[0].source_route_code()),
            ),
        ),
    );

    let mut source_identity = exact_inputs();
    source_identity[0].expected_source_identity_root = hash_domain(
        BASE_EXPECTED_SOURCE_IDENTITY_DOMAIN_V1,
        b"source-closure-wrong-expected-identity",
    );
    matched_cases[8] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&source_identity),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.expected_source_identity_root",
            "the exact declarative expected-source-identity root",
            ConstructionObservedV2::indexed_redacted(
                0,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ),
    );

    let mut snapshot_policy = exact_inputs();
    snapshot_policy[0].snapshot_policy_code ^= u16::MAX;
    matched_cases[9] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&snapshot_policy),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.snapshot_policy",
            "the exact common compiled-snapshot policy",
            ConstructionObservedV2::unsigned_pair(
                0,
                u64::from(snapshot_policy[0].snapshot_policy_code()),
            ),
        ),
    );

    let mut encoded_length = exact_inputs();
    encoded_length[0].encoded_bytes += 1;
    matched_cases[10] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&encoded_length),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.encoded_bytes",
            "the exact presented source byte length",
            ConstructionObservedV2::unsigned_pair(0, encoded_length[0].encoded_bytes()),
        ),
    );

    let mut content_root = exact_inputs();
    content_root[0].content_root = hash_domain(
        BASE_SOURCE_FILE_CONTENT_DOMAIN_V1,
        b"source-closure-wrong-content-root",
    );
    matched_cases[11] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&content_root),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.content_root",
            "the root of the exact presented source bytes",
            ConstructionObservedV2::indexed_redacted(
                0,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ),
    );

    let wrong_snapshot = hash_domain(
        BASE_SOURCE_SNAPSHOT_DOMAIN_V1,
        b"source-closure-wrong-common-snapshot",
    );
    let mut mixed_snapshot = exact_inputs();
    mixed_snapshot[0].snapshot_root = wrong_snapshot;
    matched_cases[12] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&mixed_snapshot),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.snapshot_root",
            "the one exact common compiled-snapshot root",
            ConstructionObservedV2::indexed_redacted(
                0,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ),
    );

    let mut wrong_common_snapshot = exact_inputs();
    for input in &mut wrong_common_snapshot {
        input.snapshot_root = wrong_snapshot;
    }
    matched_cases[13] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct(&wrong_common_snapshot),
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.snapshot_root",
            "the one exact common compiled-snapshot root",
            ConstructionObservedV2::indexed_redacted(
                0,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ),
    );

    let wrong_dependency_root = hash_domain(
        CURRENT_DIRECT_DEPENDENCY_DECLARATION_DOMAIN_V1,
        b"source-closure-wrong-dependency-declaration",
    );
    matched_cases[14] = source_closure_refusal_matches(
        &RunnerV2BaseSourceClosureV1::reconstruct_with_dependency_declaration(
            &exact_inputs(),
            wrong_dependency_root,
        ),
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.dependency_declaration_root",
            "the exact current declaration-time direct-dependency root",
            ConstructionObservedDataClassV2::CallerControlledText,
        ),
    );

    let positive_matched = u32::from(matched_cases[0]);
    let expected_refusals_matched = matched_cases[1..]
        .iter()
        .map(|matched| u32::from(*matched))
        .sum();
    SourceClosureExecutionReportV1 {
        positive_eligible: 1,
        positive_matched,
        expected_refusals: 14,
        expected_refusals_matched,
        unexpected_mismatches: (1 - positive_matched) + (14 - expected_refusals_matched),
        matched_cases,
    }
}

fn encode_lower_hex(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to a String cannot fail");
    }
    encoded
}

fn identity_mutation_matrix(harness: &BaseE2eHarnessIdentityV1) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    let mut positive_cells = 0_u32;
    let mut expected_refusal_cells = 0_u32;

    macro_rules! check_identity {
        ($root:expr, $type:ty, $wrong_role:expr, $label:literal) => {{
            for byte_index in 0..32 {
                checked += 1;
                positive_cells += 1;
                let mut bytes = *$root.bytes();
                bytes[byte_index] ^= 1;
                let lower_hex = encode_lower_hex(&bytes);
                let mutated = <$type>::parse_presented(
                    <$type>::DESCRIPTOR.role(),
                    <$type>::DESCRIPTOR.domain(),
                    &lower_hex,
                )
                .map_err(|_| (checked, format!("identity.{}.byte.{byte_index}", $label)))?;
                if mutated.bytes() == $root.bytes() {
                    return Err((checked, format!("identity.{}.byte.{byte_index}", $label)));
                }
            }

            checked += 1;
            expected_refusal_cells += 1;
            if <$type>::parse_presented($wrong_role, <$type>::DESCRIPTOR.domain(), &"00".repeat(32))
                != Err(IdentityError::WrongRole {
                    expected: <$type>::DESCRIPTOR.role(),
                    observed: $wrong_role,
                })
            {
                return Err((checked, format!("identity.{}.wrong-role", $label)));
            }

            checked += 1;
            expected_refusal_cells += 1;
            let wrong_domain = "org.frankensim.fs-evidence-runner.wrong-domain.v1";
            if <$type>::parse_presented(<$type>::DESCRIPTOR.role(), wrong_domain, &"00".repeat(32))
                != Err(IdentityError::WrongDomain {
                    expected: <$type>::DESCRIPTOR.domain(),
                    observed: wrong_domain.to_owned(),
                })
            {
                return Err((checked, format!("identity.{}.wrong-domain", $label)));
            }

            checked += 1;
            expected_refusal_cells += 1;
            if <$type>::parse_presented(
                <$type>::DESCRIPTOR.role(),
                <$type>::DESCRIPTOR.domain(),
                "00",
            ) != Err(IdentityError::WrongLowerHexLength {
                observed: 2,
                expected: 64,
            }) {
                return Err((checked, format!("identity.{}.wrong-length", $label)));
            }
        }};
    }

    check_identity!(
        &harness.source,
        SourceIdentityRootV2,
        DigestRoleV2::Build,
        "source"
    );
    check_identity!(
        &harness.build,
        BuildIdentityRootV2,
        DigestRoleV2::Toolchain,
        "build"
    );
    check_identity!(
        &harness.toolchain,
        ToolchainIdentityRootV2,
        DigestRoleV2::Source,
        "toolchain"
    );

    if checked != 105 || positive_cells != 96 || expected_refusal_cells != 9 {
        return Err((checked, "identity.total-count".to_owned()));
    }
    Ok(checked)
}

fn no_claim_matrix(scope: &NoClaimScopeRootV1) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    let mut positive_cells = 0_u32;
    let mut expected_refusal_cells = 0_u32;
    let lower_hex = encode_lower_hex(scope.bytes());

    checked += 1;
    positive_cells += 1;
    let exact = NoClaimScopeRootV1::parse_presented(
        NoClaimScopeRootV1::DESCRIPTOR.role(),
        NoClaimScopeRootV1::DESCRIPTOR.domain(),
        &lower_hex,
    )
    .map_err(|_| (checked, "no-claim.exact".to_owned()))?;
    if exact != *scope {
        return Err((checked, "no-claim.exact".to_owned()));
    }

    checked += 1;
    positive_cells += 1;
    let mut mutated_bytes = *scope.bytes();
    mutated_bytes[0] ^= 1;
    let mutated_hex = encode_lower_hex(&mutated_bytes);
    let mutated = NoClaimScopeRootV1::parse_presented(
        NoClaimScopeRootV1::DESCRIPTOR.role(),
        NoClaimScopeRootV1::DESCRIPTOR.domain(),
        &mutated_hex,
    )
    .map_err(|_| (checked, "no-claim.mutation".to_owned()))?;
    if mutated == *scope {
        return Err((checked, "no-claim.mutation".to_owned()));
    }

    checked += 1;
    expected_refusal_cells += 1;
    if NoClaimScopeRootV1::parse_presented(
        DigestRoleV2::Policy,
        NoClaimScopeRootV1::DESCRIPTOR.domain(),
        &lower_hex,
    ) != Err(IdentityError::WrongRole {
        expected: NoClaimScopeRootV1::DESCRIPTOR.role(),
        observed: DigestRoleV2::Policy,
    }) {
        return Err((checked, "no-claim.wrong-role".to_owned()));
    }

    checked += 1;
    expected_refusal_cells += 1;
    let wrong_domain = "org.frankensim.fs-evidence-runner.wrong-no-claim.v1";
    if NoClaimScopeRootV1::parse_presented(
        NoClaimScopeRootV1::DESCRIPTOR.role(),
        wrong_domain,
        &lower_hex,
    ) != Err(IdentityError::WrongDomain {
        expected: NoClaimScopeRootV1::DESCRIPTOR.domain(),
        observed: wrong_domain.to_owned(),
    }) {
        return Err((checked, "no-claim.wrong-domain".to_owned()));
    }

    checked += 1;
    expected_refusal_cells += 1;
    if NoClaimScopeRootV1::parse_presented(
        NoClaimScopeRootV1::DESCRIPTOR.role(),
        NoClaimScopeRootV1::DESCRIPTOR.domain(),
        "00",
    ) != Err(IdentityError::WrongLowerHexLength {
        observed: 2,
        expected: 64,
    }) {
        return Err((checked, "no-claim.wrong-length".to_owned()));
    }

    if positive_cells != 2 || expected_refusal_cells != 3 || checked != 5 {
        return Err((checked, "no-claim.total-count".to_owned()));
    }
    Ok(checked)
}

type BaseCommandOracleRowV1 = (
    RunnerCommandV2,
    u16,
    &'static str,
    RunProfileV2,
    u16,
    &'static str,
    ArtifactDispositionV2,
    u16,
    &'static str,
);

type BaseCommandApplicabilityOracleRowV1 = (
    RunnerCommandV2,
    u16,
    &'static str,
    [CommandSelectorCardinalityV2; 5],
);

const fn command_list_oracle() -> (RunnerCommandV2, u16, &'static str) {
    (RunnerCommandV2::List, 0, "list")
}

#[allow(
    clippy::too_many_lines,
    reason = "the command oracle is a frozen literal table whose complete rows must remain locally auditable"
)]
fn command_intent_oracle_rows() -> &'static [BaseCommandOracleRowV1] {
    const INTENT_ORACLE: [BaseCommandOracleRowV1; 10] = [
        (
            RunnerCommandV2::Check,
            1,
            "check",
            RunProfileV2::Smoke,
            1,
            "smoke",
            ArtifactDispositionV2::LifecycleOnlyNoBundle,
            1,
            "lifecycle-only-no-bundle",
        ),
        (
            RunnerCommandV2::Check,
            1,
            "check",
            RunProfileV2::Full,
            2,
            "full",
            ArtifactDispositionV2::LifecycleOnlyNoBundle,
            1,
            "lifecycle-only-no-bundle",
        ),
        (
            RunnerCommandV2::SelfTest,
            2,
            "self-test",
            RunProfileV2::Smoke,
            1,
            "smoke",
            ArtifactDispositionV2::LifecycleOnlyNoBundle,
            1,
            "lifecycle-only-no-bundle",
        ),
        (
            RunnerCommandV2::SelfTest,
            2,
            "self-test",
            RunProfileV2::Full,
            2,
            "full",
            ArtifactDispositionV2::LifecycleOnlyNoBundle,
            1,
            "lifecycle-only-no-bundle",
        ),
        (
            RunnerCommandV2::Run,
            3,
            "run",
            RunProfileV2::Smoke,
            1,
            "smoke",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
        (
            RunnerCommandV2::Run,
            3,
            "run",
            RunProfileV2::Full,
            2,
            "full",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
        (
            RunnerCommandV2::Negative,
            4,
            "negative",
            RunProfileV2::Smoke,
            1,
            "smoke",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
        (
            RunnerCommandV2::Negative,
            4,
            "negative",
            RunProfileV2::Full,
            2,
            "full",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
        (
            RunnerCommandV2::Replay,
            5,
            "replay",
            RunProfileV2::Smoke,
            1,
            "smoke",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
        (
            RunnerCommandV2::Replay,
            5,
            "replay",
            RunProfileV2::Full,
            2,
            "full",
            ArtifactDispositionV2::DurableBundleRequired,
            2,
            "durable-bundle-required",
        ),
    ];
    &INTENT_ORACLE
}

fn command_applicability_oracle_rows() -> &'static [BaseCommandApplicabilityOracleRowV1] {
    use CommandSelectorCardinalityV2::{Absent, Singular};
    const APPLICABILITY_ORACLE: [BaseCommandApplicabilityOracleRowV1; 6] = [
        (RunnerCommandV2::List, 0, "list", [Absent; 5]),
        (RunnerCommandV2::Check, 1, "check", [Absent; 5]),
        (RunnerCommandV2::SelfTest, 2, "self-test", [Absent; 5]),
        (
            RunnerCommandV2::Run,
            3,
            "run",
            [Singular, Singular, Singular, Absent, Absent],
        ),
        (
            RunnerCommandV2::Negative,
            4,
            "negative",
            [Absent, Absent, Absent, Singular, Absent],
        ),
        (
            RunnerCommandV2::Replay,
            5,
            "replay",
            [Absent, Absent, Absent, Absent, Singular],
        ),
    ];
    &APPLICABILITY_ORACLE
}

const COMMAND_SELECTOR_FIELD_ORACLE_V1: [CommandSelectorFieldV2; 5] = [
    CommandSelectorFieldV2::Family,
    CommandSelectorFieldV2::Mode,
    CommandSelectorFieldV2::Profile,
    CommandSelectorFieldV2::NegativeCase,
    CommandSelectorFieldV2::ReplaySource,
];

const fn command_selector_requirement_oracle_code(
    requirement: CommandSelectorCardinalityV2,
) -> u16 {
    match requirement {
        CommandSelectorCardinalityV2::Absent => 0,
        CommandSelectorCardinalityV2::Singular => 1,
        CommandSelectorCardinalityV2::Duplicate => 2,
        CommandSelectorCardinalityV2::Ambiguous => 3,
    }
}

const fn command_oracle_identity_code(command: RunnerCommandV2) -> u16 {
    match command {
        RunnerCommandV2::List => 0,
        RunnerCommandV2::Check => 1,
        RunnerCommandV2::SelfTest => 2,
        RunnerCommandV2::Run => 3,
        RunnerCommandV2::Negative => 4,
        RunnerCommandV2::Replay => 5,
    }
}

const fn command_profile_oracle_identity_code(profile: RunProfileV2) -> u16 {
    match profile {
        RunProfileV2::Smoke => 1,
        RunProfileV2::Full => 2,
    }
}

const fn command_disposition_oracle_identity_code(disposition: ArtifactDispositionV2) -> u16 {
    match disposition {
        ArtifactDispositionV2::LifecycleOnlyNoBundle => 1,
        ArtifactDispositionV2::DurableBundleRequired => 2,
    }
}

fn command_oracle_table_root_from_rows(
    list: (RunnerCommandV2, u16, &'static str),
    rows: &[BaseCommandOracleRowV1],
    applicability: &[BaseCommandApplicabilityOracleRowV1],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(2048);
    let (list_command, list_code, list_name) = list;
    bytes.extend_from_slice(&command_oracle_identity_code(list_command).to_be_bytes());
    bytes.extend_from_slice(&list_code.to_be_bytes());
    detail_push_str(&mut bytes, list_name);
    bytes.extend_from_slice(
        &u32::try_from(rows.len())
            .expect("the literal command oracle count fits u32")
            .to_be_bytes(),
    );
    for &(
        command,
        command_code,
        command_name,
        profile,
        profile_code,
        profile_name,
        disposition,
        disposition_code,
        disposition_name,
    ) in rows
    {
        bytes.extend_from_slice(&command_oracle_identity_code(command).to_be_bytes());
        bytes.extend_from_slice(&command_code.to_be_bytes());
        detail_push_str(&mut bytes, command_name);
        bytes.extend_from_slice(&command_profile_oracle_identity_code(profile).to_be_bytes());
        bytes.extend_from_slice(&profile_code.to_be_bytes());
        detail_push_str(&mut bytes, profile_name);
        bytes.extend_from_slice(
            &command_disposition_oracle_identity_code(disposition).to_be_bytes(),
        );
        bytes.extend_from_slice(&disposition_code.to_be_bytes());
        detail_push_str(&mut bytes, disposition_name);
    }
    bytes.extend_from_slice(
        &u32::try_from(applicability.len())
            .expect("the literal command applicability count fits u32")
            .to_be_bytes(),
    );
    for &(command, command_code, command_name, requirements) in applicability {
        bytes.extend_from_slice(&command_oracle_identity_code(command).to_be_bytes());
        bytes.extend_from_slice(&command_code.to_be_bytes());
        detail_push_str(&mut bytes, command_name);
        for requirement in requirements {
            bytes.extend_from_slice(
                &command_selector_requirement_oracle_code(requirement).to_be_bytes(),
            );
        }
    }
    hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-command-literal-oracle.v1",
        &bytes,
    )
}

fn command_oracle_table_root() -> ContentHash {
    command_oracle_table_root_from_rows(
        command_list_oracle(),
        command_intent_oracle_rows(),
        command_applicability_oracle_rows(),
    )
}

fn command_matrix() -> Result<u32, (u32, String)> {
    for &(command, command_code, command_name, requirements) in command_applicability_oracle_rows()
    {
        if command.code() != command_code || command.name() != command_name {
            return Err((
                COMMAND_APPLICABILITY_SETUP_SEMANTIC_ORDINAL_V1,
                format!("command.{command_name}.applicability-identity"),
            ));
        }
        let observed = CommandSelectorPresenceV2::exact_for(command);
        if COMMAND_SELECTOR_FIELD_ORACLE_V1
            .iter()
            .zip(requirements)
            .any(|(&field, expected)| observed.cardinality(field) != expected)
            || validate_command_selector_presence_v2(command, observed).is_err()
        {
            return Err((
                COMMAND_APPLICABILITY_SETUP_SEMANTIC_ORDINAL_V1,
                format!("command.{command_name}.applicability"),
            ));
        }
    }

    let (list_command, list_code, list_name) = command_list_oracle();
    let list = CommandIntentV2::list();
    if list.command() != list_command
        || list.command().code() != list_code
        || list.command().name() != list_name
        || list.selection().is_some()
        || list.budgets().is_some()
        || list.disposition().is_some()
        || list.publication_selection().is_some()
    {
        return Err((1, "command.list".to_owned()));
    }
    let mut checked = 1_u32;

    for &(
        command,
        command_code,
        command_name,
        profile,
        profile_code,
        profile_name,
        disposition,
        disposition_code,
        disposition_name,
    ) in command_intent_oracle_rows()
    {
        checked += 1;
        let selection = command_selection(command, profile)
            .map_err(|_| (checked, format!("command.{}.selection", command.name())))?;
        let budgets = command_budgets(profile, disposition)
            .map_err(|_| (checked, format!("command.{}.budgets", command.name())))?;
        let publication = if disposition == ArtifactDispositionV2::DurableBundleRequired {
            Some(
                publication_selection("results/bundle")
                    .map_err(|_| (checked, format!("command.{}.publication", command.name())))?,
            )
        } else {
            None
        };
        let intent = CommandIntentV2::new(command, selection, budgets, publication)
            .map_err(|_| (checked, format!("command.{}.intent", command.name())))?;
        if command.code() != command_code
            || command.name() != command_name
            || profile.code() != profile_code
            || profile.name() != profile_name
            || disposition.code() != disposition_code
            || disposition.name() != disposition_name
            || intent.command() != command
            || intent
                .selection()
                .is_none_or(|selection| selection.profile() != profile)
            || intent.disposition() != Some(disposition)
            || intent.publication_selection().is_some()
                != (disposition == ArtifactDispositionV2::DurableBundleRequired)
        {
            return Err((checked, format!("command.{}.cell", command.name())));
        }
    }

    if checked != 11 {
        return Err((checked, "command.total-count".to_owned()));
    }
    Ok(checked)
}

fn command_selection(
    command: RunnerCommandV2,
    profile: RunProfileV2,
) -> Result<CommandSelectionV2, ConstructionErrorV2> {
    let family = token("family.fixture")?;
    let mode = token("mode.fixture")?;
    let manifest = |byte: u8| {
        CaseManifestRootV2::parse_presented(
            CaseManifestRootV2::DESCRIPTOR.role(),
            CaseManifestRootV2::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .map_err(|_| presented_identity_error("command.case_manifest"))
    };
    match command {
        RunnerCommandV2::Check => Ok(CommandSelectionV2::fixed_preflight(
            manifest(1)?,
            family,
            mode,
            profile,
        )),
        RunnerCommandV2::SelfTest => Ok(CommandSelectionV2::fixed_self_test(
            manifest(2)?,
            family,
            mode,
            profile,
        )),
        RunnerCommandV2::Run => Ok(CommandSelectionV2::caller_run(family, mode, profile)),
        RunnerCommandV2::Negative => Ok(CommandSelectionV2::sealed_negative(
            manifest(3)?,
            family,
            mode,
            profile,
        )),
        RunnerCommandV2::Replay => {
            let source = SourceIdentityRootV2::parse_presented(
                SourceIdentityRootV2::DESCRIPTOR.role(),
                SourceIdentityRootV2::DESCRIPTOR.domain(),
                &"04".repeat(32),
            )
            .map_err(|_| presented_identity_error("command.replay_source"))?;
            Ok(CommandSelectionV2::sealed_replay(
                source, family, mode, profile,
            ))
        }
        RunnerCommandV2::List => Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "command.selection",
            "typed absence for List",
            ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Present),
        )),
    }
}

fn command_budgets(
    profile: RunProfileV2,
    disposition: ArtifactDispositionV2,
) -> Result<AdmittedRunnerBudgetsV2, RunnerBudgetViolationV2> {
    let mut candidate = durable_budget_candidate();
    if disposition == ArtifactDispositionV2::LifecycleOnlyNoBundle {
        candidate.max_child_processes = 0;
        candidate.max_parallel_children = 0;
        candidate.combined_child_stdout_bytes = 0;
        candidate.combined_child_stderr_bytes = 0;
        candidate.artifact_encoded_bytes = 0;
        candidate.artifact_stored_bytes = 0;
        candidate.artifact_expanded_bytes = 0;
        candidate.system_publication_stored_bytes = 0;
        candidate.publication_stored_bytes = 0;
    }
    RunnerBudgetsV2::try_new(candidate)?.admit(profile, disposition, &RunnerLimitsV2::base(profile))
}

fn presented_identity_error(field: &'static str) -> ConstructionErrorV2 {
    ConstructionErrorV2::new_redacted(
        ConstructionErrorKindV2::Incompatible,
        field,
        "an exact nominal presented identity fixture",
        ConstructionObservedDataClassV2::CallerControlledText,
    )
}

fn opaque_root(root: ContentHash) -> Result<TypedValueV2, ConstructionErrorV2> {
    OpaqueBytesV2::new(root.as_bytes().to_vec())
        .map(TypedValueV2::OpaqueBytes)
        .map_err(|_error| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.opaque_root",
                "an exact 32-byte domain-separated root",
                ConstructionObservedDataClassV2::BulkPayload,
            )
        })
}

/// Return frozen expected-inventory metadata for one case kind.
///
/// These catalog-size fields describe the closed matrix the case intends to
/// execute. They are not observed progress counters: a red execution may
/// report a shorter nonzero `checked_cells` prefix without changing this
/// expected inventory.
#[allow(
    clippy::too_many_lines,
    reason = "the log-detail inventory is an exhaustive compiler-checked match over the closed base case catalog"
)]
fn case_detail_fields(
    kind: BaseE2eCaseKindV1,
) -> Result<Vec<BaseE2eLogFieldV1>, ConstructionErrorV2> {
    let fields = match kind {
        BaseE2eCaseKindV1::CatalogLiterals => {
            vec![field("catalog-literal-cells", TypedValueV2::U32(186))?]
        }
        BaseE2eCaseKindV1::LimitCatalog => vec![
            field("limit-field-count", TypedValueV2::U32(71))?,
            field("limit-profile-cells", TypedValueV2::U32(284))?,
        ],
        BaseE2eCaseKindV1::BudgetAdmission => vec![
            field("budget-field-count", TypedValueV2::U32(18))?,
            field("logical-unit-count", TypedValueV2::U32(16))?,
        ],
        BaseE2eCaseKindV1::CapabilityLeastPrivilege | BaseE2eCaseKindV1::CapabilityExtraRight => {
            vec![
                field("capability-valid-cells", TypedValueV2::U32(12))?,
                field("capability-mutant-cells", TypedValueV2::U32(390))?,
                field("capability-right-count", TypedValueV2::U32(10))?,
            ]
        }
        BaseE2eCaseKindV1::StatePass | BaseE2eCaseKindV1::StateUsageInLifecycle => {
            let (cancelled, timed_out, internal_error) = presented_stop_fixture_roots()?;
            vec![
                field(
                    "cancelled-causal-root",
                    TypedValueV2::Digest(cancelled.digest().clone()),
                )?,
                field(
                    "diagnostic-code-count",
                    TypedValueV2::U32(
                        u32::try_from(DiagnosticCodeV2::ALL.len())
                            .expect("closed diagnostic count fits u32"),
                    ),
                )?,
                field(
                    "internal-error-causal-root",
                    TypedValueV2::Digest(internal_error.digest().clone()),
                )?,
                field("lowest-manifest-ordinal", TypedValueV2::U32(0))?,
                field("maximum-manifest-ordinal", TypedValueV2::U32(255))?,
                field(
                    "record-role-count",
                    TypedValueV2::U32(
                        u32::try_from(StateBearingRecordRoleV2::ALL.len())
                            .expect("closed role count fits u32"),
                    ),
                )?,
                field(
                    "refused-reason-count",
                    TypedValueV2::U32(
                        u32::try_from(RefusedReasonV2::ALL.len())
                            .expect("closed reason count fits u32"),
                    ),
                )?,
                field("state-matrix-cells", TypedValueV2::U32(32_448))?,
                field(
                    "timed-out-causal-root",
                    TypedValueV2::Digest(timed_out.digest().clone()),
                )?,
            ]
        }
        BaseE2eCaseKindV1::Diagnostic | BaseE2eCaseKindV1::DiagnosticRankGap => vec![
            field("diagnostic-code-count", TypedValueV2::U32(12))?,
            field("diagnostic-expected", TypedValueV2::U64(4))?,
            field("diagnostic-observed", TypedValueV2::U64(5))?,
            field(
                "diagnostic-owner",
                TypedValueV2::Token(token("runner.owner")?),
            )?,
            field("diagnostic-prerequisite-count", TypedValueV2::U32(1))?,
            field("diagnostic-repair-count", TypedValueV2::U32(1))?,
            field(
                "diagnostic-retryability-count",
                TypedValueV2::U32(
                    u32::try_from(RetryabilityV2::ALL.len())
                        .expect("closed retryability count fits u32"),
                ),
            )?,
            field(
                "repair-kind-count",
                TypedValueV2::U32(
                    u32::try_from(RepairActionKindV2::ALL.len())
                        .expect("closed repair-kind count fits u32"),
                ),
            )?,
        ],
        BaseE2eCaseKindV1::IdentityMutation => {
            vec![field("identity-mutation-cells", TypedValueV2::U32(105))?]
        }
        BaseE2eCaseKindV1::PublicationStorage => {
            let storage = publication_storage()?;
            vec![
                field("artifact-stored-bytes", TypedValueV2::U64(storage.artifact))?,
                field(
                    "system-publication-stored-bytes",
                    TypedValueV2::U64(storage.system_publication),
                )?,
                field(
                    "publication-stored-bytes",
                    TypedValueV2::U64(storage.publication),
                )?,
                field(
                    "stored-byte-unit",
                    TypedValueV2::Token(token("stored-bytes")?),
                )?,
            ]
        }
        _ => Vec::new(),
    };
    Ok(fields)
}

type BaseCatalogOracleRowV1 = (u8, &'static str, u16, &'static str, bool);

#[allow(
    clippy::too_many_lines,
    reason = "the literal matrix keeps every independent catalog oracle and its observed lookup adjacent in canonical order"
)]
fn catalog_literal_matrix_and_rows(
    observe: bool,
) -> Result<(u32, Vec<BaseCatalogOracleRowV1>), (u32, String)> {
    let mut checked = 0_u32;
    let mut oracle_rows = Vec::with_capacity(186);

    macro_rules! check_closed_catalog {
        ($catalog:ty, $label:literal, $rows:expr) => {
            for &(code, name) in $rows {
                checked += 1;
                oracle_rows.push((1, $label, code, name, false));
                if observe {
                    let Ok(value) = <$catalog>::from_code(code) else {
                        return Err((checked, format!("catalog.{}.{}", $label, code)));
                    };
                    if value.code() != code || value.name() != name {
                        return Err((checked, format!("catalog.{}.{}", $label, code)));
                    }
                }
            }
        };
    }

    macro_rules! check_registered_catalog {
        ($catalog:ty, $label:literal, $rows:expr) => {
            for &(tag, name, requires_registered_id) in $rows {
                checked += 1;
                oracle_rows.push((2, $label, tag, name, requires_registered_id));
                if observe {
                    let registered_id = requires_registered_id.then_some(7);
                    let Ok(value) = <$catalog>::from_tag(tag, registered_id) else {
                        return Err((checked, format!("catalog.{}.{}", $label, tag)));
                    };
                    if value.tag() != tag
                        || value.name() != name
                        || value.registered_id() != registered_id
                    {
                        return Err((checked, format!("catalog.{}.{}", $label, tag)));
                    }
                }
            }
        };
    }

    check_closed_catalog!(
        RunnerApiGeneration,
        "api-generation",
        &[(2, "RunnerSpecV2")]
    );
    check_closed_catalog!(RunnerWireVersion, "wire-version", &[(1, "runner-wire-v1")]);
    check_closed_catalog!(
        ArtifactCodecIdV2,
        "artifact-codec",
        &[(0, "identity"), (1, "zstd-frame-v1")]
    );
    checked += 1;
    oracle_rows.push((3, "wire-predecessor", 1, "no-predecessor", false));
    if observe
        && (WirePredecessorPolicyV1::NoPredecessor.name() != "no-predecessor"
            || WirePredecessorPolicyV1::NoPredecessor
                .predecessor()
                .is_some())
    {
        return Err((checked, "catalog.wire-predecessor.1".to_owned()));
    }
    check_closed_catalog!(
        ProofExitV2,
        "proof-exit",
        &[
            (0, "pass"),
            (10, "failed"),
            (11, "refused"),
            (12, "no-data"),
            (13, "stale"),
            (14, "environment-invalid"),
            (15, "blocked"),
            (16, "unsupported"),
            (17, "not-run"),
            (18, "cancelled"),
            (19, "timed-out"),
            (64, "usage"),
            (70, "internal-error"),
        ]
    );
    check_closed_catalog!(
        RefusedReasonV2,
        "refused-reason",
        &[
            (1, "invalid-evidence"),
            (2, "non-canonical-evidence"),
            (3, "evidence-identity-mismatch"),
            (4, "evidence-tampered"),
            (5, "limit-exceeded"),
            (6, "unsafe-artifact-placement"),
            (7, "artifact-collision"),
            (8, "lifecycle-violation"),
            (9, "policy-refused"),
            (10, "authority-boundary-violation"),
            (11, "migration-refused"),
        ]
    );
    check_closed_catalog!(
        RunnerCommandV2,
        "runner-command",
        &[
            (0, "list"),
            (1, "check"),
            (2, "self-test"),
            (3, "run"),
            (4, "negative"),
            (5, "replay"),
        ]
    );
    check_closed_catalog!(RunProfileV2, "run-profile", &[(1, "smoke"), (2, "full")]);
    check_closed_catalog!(
        ArtifactDispositionV2,
        "artifact-disposition",
        &[
            (1, "lifecycle-only-no-bundle"),
            (2, "durable-bundle-required"),
        ]
    );
    check_closed_catalog!(
        PlatformPathProfileV2,
        "path-profile",
        &[
            (1, "posix-descriptor-relative-v1"),
            (2, "windows-handle-relative-v1"),
            (3, "content-store-object-key-v1"),
        ]
    );
    check_closed_catalog!(
        LifecycleRecordKindV2,
        "record-kind",
        &[
            (1, "run-start"),
            (2, "case-start"),
            (3, "family-row"),
            (4, "case-terminal"),
            (5, "run-summary"),
            (6, "run-terminal"),
        ]
    );
    check_closed_catalog!(
        StateBearingRecordRoleV2,
        "record-role",
        &[
            (1, "pre-run-diagnostic"),
            (2, "executed-case-terminal"),
            (3, "suppressed-case-terminal"),
            (4, "run-terminal"),
        ]
    );
    check_closed_catalog!(
        DiagnosticCodeV2,
        "diagnostic-code",
        &[
            (1, "case.conformance_mismatch"),
            (2, "runner.not_run"),
            (3, "runner.refused"),
            (4, "runner.no_data"),
            (5, "runner.stale"),
            (6, "runner.environment_invalid"),
            (7, "runner.blocked"),
            (8, "runner.unsupported"),
            (9, "runner.cancelled"),
            (10, "runner.timed_out"),
            (11, "runner.usage"),
            (12, "runner.internal_error"),
        ]
    );
    check_closed_catalog!(
        RetryabilityV2,
        "retryability",
        &[
            (0, "never"),
            (1, "same-invocation"),
            (2, "after-input-change"),
            (3, "after-environment-change"),
            (4, "after-prerequisite-change"),
        ]
    );
    check_closed_catalog!(
        RepairActionKindV2,
        "repair-kind",
        &[
            (1, "change-arguments"),
            (2, "supply-evidence"),
            (3, "regenerate-canonical-evidence"),
            (4, "refresh-evidence"),
            (5, "reduce-resource-demand"),
            (6, "choose-safe-artifact-destination"),
            (7, "restore-lifecycle"),
            (8, "update-policy-or-capability"),
            (9, "register-migration"),
            (10, "retry-same-invocation"),
            (11, "contact-owner"),
            (12, "inspect-retained-artifact"),
        ]
    );
    check_closed_catalog!(
        NotRunCauseCodeV2,
        "not-run-cause",
        &[
            (1, "prior-cancelled"),
            (2, "prior-timed-out"),
            (3, "prior-controlled-internal-error"),
        ]
    );
    check_closed_catalog!(
        TypedValueTagV2,
        "typed-value",
        &[
            (1, "i8"),
            (2, "i16"),
            (3, "i32"),
            (4, "i64"),
            (5, "i128"),
            (6, "u8"),
            (7, "u16"),
            (8, "u32"),
            (9, "u64"),
            (10, "u128"),
            (11, "rational"),
            (12, "decimal"),
            (13, "f32-bits"),
            (14, "f64-bits"),
            (15, "digest"),
            (16, "quantity"),
            (17, "token"),
            (18, "text"),
            (19, "relative-path"),
            (20, "opaque-bytes"),
        ]
    );
    check_closed_catalog!(
        TypedOptionTagV1,
        "typed-option",
        &[(0, "absent"), (1, "present")]
    );
    check_closed_catalog!(
        DigestRoleV2,
        "digest-role",
        &[
            (1, "spec"),
            (2, "invocation"),
            (3, "run"),
            (4, "source"),
            (5, "build"),
            (6, "toolchain"),
            (7, "case-manifest"),
            (8, "artifact-encoded"),
            (9, "artifact-content"),
            (10, "stored-object"),
            (11, "artifact-inventory"),
            (12, "lifecycle-log"),
            (13, "run-summary"),
            (14, "run-terminal"),
            (15, "bundle-manifest"),
            (16, "durable-publication"),
            (17, "seal"),
            (18, "published-bundle-receipt"),
            (19, "policy"),
            (20, "candidate-bytes"),
            (21, "candidate-schema"),
            (22, "source-closure"),
            (23, "claim-scope"),
            (24, "producer-manifest"),
            (25, "registered-family-domain"),
        ]
    );
    check_closed_catalog!(
        PublicationProtocolV2,
        "publication-protocol",
        &[
            (1, "posix-descriptor-rename-and-directory-sync-v1"),
            (2, "windows-handle-replace-and-directory-flush-v1"),
            (3, "content-store-atomic-commit-v1"),
        ]
    );
    check_closed_catalog!(
        DestinationAdmissionModeV2,
        "destination-mode",
        &[(1, "absent"), (2, "pre-existing-empty")]
    );
    check_closed_catalog!(
        RootCapabilityAccessV2,
        "capability-access",
        &[(1, "read-only-input"), (2, "durable-output")]
    );
    check_closed_catalog!(
        RootCapabilityRightV2,
        "capability-right",
        &[
            (1, "traverse"),
            (2, "read-object"),
            (3, "enumerate"),
            (4, "create-object"),
            (5, "populate-empty-destination"),
            (6, "sync-object"),
            (7, "sync-container"),
            (8, "acquire-exclusive-lease"),
            (9, "query-generation"),
            (10, "commit-compare-and-swap"),
        ]
    );
    check_closed_catalog!(
        OverlapPolicyRelationV2,
        "overlap-relation",
        &[(1, "require-input-output-disjoint")]
    );
    check_registered_catalog!(
        RootClassV2,
        "root-class",
        &[
            (1, "input-artifact-root", false),
            (2, "output-artifact-root", false),
            (3, "other", true),
        ]
    );
    check_registered_catalog!(
        LogicalUnitV2,
        "logical-unit",
        &[
            (1, "encoded-bytes", false),
            (2, "expanded-bytes", false),
            (3, "stored-bytes", false),
            (4, "logical-bytes", false),
            (5, "count", false),
            (6, "records", false),
            (7, "rows", false),
            (8, "elements", false),
            (9, "samples", false),
            (10, "iterations", false),
            (11, "operations", false),
            (12, "cycles", false),
            (13, "nanoseconds", false),
            (14, "seconds", false),
            (15, "dimensionless", false),
            (16, "registered-unit", true),
        ]
    );
    check_registered_catalog!(
        ArtifactRoleV2,
        "artifact-role",
        &[
            (1, "observation", false),
            (2, "comparison-detail", false),
            (3, "effect-detail", false),
            (4, "diagnostic-log", false),
            (5, "family-evidence", false),
            (6, "performance-evidence", false),
            (7, "replay-support", false),
            (8, "registered-family-role", true),
        ]
    );
    check_registered_catalog!(
        LogicalExtentAxisV2,
        "logical-axis",
        &[
            (1, "payload", false),
            (2, "records", false),
            (3, "rows", false),
            (4, "elements", false),
            (5, "samples", false),
            (6, "iterations", false),
            (7, "operations", false),
            (8, "cycles", false),
            (9, "duration", false),
            (10, "registered-axis", true),
        ]
    );

    if checked != 186 {
        return Err((checked, "catalog.total-count".to_owned()));
    }
    debug_assert_eq!(usize::try_from(checked).ok(), Some(oracle_rows.len()));
    Ok((checked, oracle_rows))
}

fn catalog_oracle_table_root(rows: &[BaseCatalogOracleRowV1]) -> ContentHash {
    let mut bytes = Vec::with_capacity(4 + rows.len() * 64);
    bytes.extend_from_slice(
        &u32::try_from(rows.len())
            .expect("the literal catalog oracle count fits u32")
            .to_be_bytes(),
    );
    for &(family, label, code, name, requires_registered_id) in rows {
        bytes.push(family);
        detail_push_str(&mut bytes, label);
        bytes.extend_from_slice(&code.to_be_bytes());
        detail_push_str(&mut bytes, name);
        detail_push_bool(&mut bytes, requires_registered_id);
    }
    hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-catalog-literal-oracle.v1",
        &bytes,
    )
}

fn catalog_literal_matrix_and_root(observe: bool) -> Result<(u32, ContentHash), (u32, String)> {
    let (checked, rows) = catalog_literal_matrix_and_rows(observe)?;
    Ok((checked, catalog_oracle_table_root(&rows)))
}

#[cfg(test)]
fn catalog_literal_oracle_rows() -> Vec<BaseCatalogOracleRowV1> {
    catalog_literal_matrix_and_rows(false)
        .expect("literal-only catalog oracle assembly cannot refuse")
        .1
}

fn catalog_literal_oracle_root() -> ContentHash {
    catalog_literal_matrix_and_root(false)
        .expect("literal-only catalog oracle encoding cannot refuse")
        .1
}

fn catalog_literal_matrix() -> Result<u32, (u32, String)> {
    catalog_literal_matrix_and_root(true).map(|(checked, _)| checked)
}

type BaseLimitOracleRowV1 = (
    RunnerLimitFieldV2,
    u16,
    &'static str,
    RunnerLimitUnitV2,
    crate::limits::RunnerLimitWidthV2,
    RunnerLimitTightenabilityV2,
    RunnerLimitValueV2,
    RunnerLimitValueV2,
);

const fn limit_field_oracle_identity_code(field: RunnerLimitFieldV2) -> u16 {
    match field {
        RunnerLimitFieldV2::ArgvTokens => 1,
        RunnerLimitFieldV2::ArgvTokenBytes => 2,
        RunnerLimitFieldV2::ArgvAggregateBytes => 3,
        RunnerLimitFieldV2::LifecycleRecordEncodedBytes => 4,
        RunnerLimitFieldV2::CaseLifecycleRecords => 5,
        RunnerLimitFieldV2::CaseLifecycleEncodedBytes => 6,
        RunnerLimitFieldV2::FamilyRowsPerCase => 7,
        RunnerLimitFieldV2::InvocationCases => 8,
        RunnerLimitFieldV2::LifecycleDocumentRecords => 9,
        RunnerLimitFieldV2::LifecycleDocumentEncodedBytes => 10,
        RunnerLimitFieldV2::CommandResultStdoutBytes => 11,
        RunnerLimitFieldV2::ChildStdoutBytes => 12,
        RunnerLimitFieldV2::CombinedChildStdoutBytes => 13,
        RunnerLimitFieldV2::ChildStderrBytes => 14,
        RunnerLimitFieldV2::CombinedChildStderrBytes => 15,
        RunnerLimitFieldV2::ManifestEncodedBytes => 16,
        RunnerLimitFieldV2::NestingDepth => 17,
        RunnerLimitFieldV2::ComparisonNodes => 18,
        RunnerLimitFieldV2::EffectNodes => 19,
        RunnerLimitFieldV2::TextBytes => 20,
        RunnerLimitFieldV2::StableTokenBytes => 21,
        RunnerLimitFieldV2::BundleRelativePathBytes => 22,
        RunnerLimitFieldV2::DiagnosticsPerCase => 23,
        RunnerLimitFieldV2::DiagnosticsPerRun => 24,
        RunnerLimitFieldV2::PrerequisitesPerDiagnostic => 25,
        RunnerLimitFieldV2::RepairsPerDiagnostic => 26,
        RunnerLimitFieldV2::Artifacts => 27,
        RunnerLimitFieldV2::ArtifactEncodedBytes => 28,
        RunnerLimitFieldV2::ArtifactExpandedBytes => 29,
        RunnerLimitFieldV2::ArtifactStoredBytes => 30,
        RunnerLimitFieldV2::BundleEncodedBytes => 31,
        RunnerLimitFieldV2::BundleExpandedBytes => 32,
        RunnerLimitFieldV2::ArtifactStoredAggregateBytes => 33,
        RunnerLimitFieldV2::SystemPublicationStoredBytes => 34,
        RunnerLimitFieldV2::PublicationStoredBytes => 35,
        RunnerLimitFieldV2::ChildStreamDiscardBytes => 36,
        RunnerLimitFieldV2::ModesPerFamily => 37,
        RunnerLimitFieldV2::ExtensionDiagnosticsPerFamily => 38,
        RunnerLimitFieldV2::ArtifactRolesPerFamily => 39,
        RunnerLimitFieldV2::RootPoliciesPerFamily => 40,
        RunnerLimitFieldV2::RegisteredUnitsPerFamily => 41,
        RunnerLimitFieldV2::DigestDomainsPerFamily => 42,
        RunnerLimitFieldV2::ExtensionSchemasPerFamily => 43,
        RunnerLimitFieldV2::ExecutableDescriptorsPerFamily => 44,
        RunnerLimitFieldV2::MapEntries => 45,
        RunnerLimitFieldV2::GenericArrayItems => 46,
        RunnerLimitFieldV2::PathSegments => 47,
        RunnerLimitFieldV2::IntegerDigits => 48,
        RunnerLimitFieldV2::RationalComponentBytes => 49,
        RunnerLimitFieldV2::DecimalCoefficientBytes => 50,
        RunnerLimitFieldV2::DecimalAbsoluteScale => 51,
        RunnerLimitFieldV2::LogicalExtentsPerArtifact => 52,
        RunnerLimitFieldV2::ObservationKeysPerCase => 53,
        RunnerLimitFieldV2::DecisionDetailNamespaces => 54,
        RunnerLimitFieldV2::OutputClasses => 55,
        RunnerLimitFieldV2::OpaqueValueBytes => 56,
        RunnerLimitFieldV2::RetainedUnknownExtensionBytes => 57,
        RunnerLimitFieldV2::ExpressionEdges => 58,
        RunnerLimitFieldV2::MemoizedEvaluationVisits => 59,
        RunnerLimitFieldV2::RepairActionEncodedBytes => 60,
        RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes => 61,
        RunnerLimitFieldV2::FailureStderrEncodedBytes => 62,
        RunnerLimitFieldV2::RunnerCatalogEncodedBytes => 63,
        RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes => 64,
        RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes => 65,
        RunnerLimitFieldV2::RegisteredExtentAxesPerFamily => 66,
        RunnerLimitFieldV2::RegisteredObservationKeysPerFamily => 67,
        RunnerLimitFieldV2::RegisteredAuthorityScopesPerFamily => 68,
        RunnerLimitFieldV2::RegisteredExternalRootClassesPerFamily => 69,
        RunnerLimitFieldV2::RegisteredEvaluationUnitsPerFamily => 70,
        RunnerLimitFieldV2::RegisteredResourceIdentitiesPerFamily => 71,
    }
}

const fn limit_unit_oracle_code(unit: RunnerLimitUnitV2) -> u16 {
    match unit {
        RunnerLimitUnitV2::Count => 1,
        RunnerLimitUnitV2::Records => 2,
        RunnerLimitUnitV2::Rows => 3,
        RunnerLimitUnitV2::EncodedBytes => 4,
        RunnerLimitUnitV2::ExpandedBytes => 5,
        RunnerLimitUnitV2::StoredBytes => 6,
        RunnerLimitUnitV2::LogicalBytes => 7,
        RunnerLimitUnitV2::Depth => 8,
        RunnerLimitUnitV2::Nodes => 9,
        RunnerLimitUnitV2::Digits => 10,
        RunnerLimitUnitV2::Segments => 11,
        RunnerLimitUnitV2::Diagnostics => 12,
        RunnerLimitUnitV2::Prerequisites => 13,
        RunnerLimitUnitV2::Repairs => 14,
        RunnerLimitUnitV2::Artifacts => 15,
        RunnerLimitUnitV2::Namespaces => 16,
        RunnerLimitUnitV2::Classes => 17,
        RunnerLimitUnitV2::Visits => 18,
        RunnerLimitUnitV2::DecimalScale => 19,
    }
}

const fn limit_width_oracle_code(width: crate::limits::RunnerLimitWidthV2) -> u16 {
    match width {
        crate::limits::RunnerLimitWidthV2::U32 => 1,
        crate::limits::RunnerLimitWidthV2::U64 => 2,
    }
}

const fn limit_tightenability_oracle_code(tightenability: RunnerLimitTightenabilityV2) -> u16 {
    match tightenability {
        RunnerLimitTightenabilityV2::Tightenable => 1,
        RunnerLimitTightenabilityV2::Fixed => 2,
    }
}

fn push_limit_oracle_value(bytes: &mut Vec<u8>, value: RunnerLimitValueV2) {
    match value {
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

fn limit_oracle_table_root(rows: &[BaseLimitOracleRowV1]) -> ContentHash {
    let mut bytes = Vec::with_capacity(rows.len() * 96);
    bytes.extend_from_slice(
        &u32::try_from(rows.len())
            .expect("the literal limit oracle count fits u32")
            .to_be_bytes(),
    );
    for &(field, ordinal, name, unit, width, tightenability, smoke, full) in rows {
        bytes.extend_from_slice(&limit_field_oracle_identity_code(field).to_be_bytes());
        bytes.extend_from_slice(&ordinal.to_be_bytes());
        detail_push_str(&mut bytes, name);
        bytes.extend_from_slice(&limit_unit_oracle_code(unit).to_be_bytes());
        bytes.extend_from_slice(&limit_width_oracle_code(width).to_be_bytes());
        bytes.extend_from_slice(&limit_tightenability_oracle_code(tightenability).to_be_bytes());
        push_limit_oracle_value(&mut bytes, smoke);
        push_limit_oracle_value(&mut bytes, full);
    }
    hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-limit-literal-oracle.v1",
        &bytes,
    )
}

#[allow(
    clippy::too_many_lines,
    reason = "the 71-row limit oracle is deliberately handwritten and contiguous so every literal remains directly auditable"
)]
fn limit_oracle_rows() -> &'static [BaseLimitOracleRowV1] {
    use crate::limits::RunnerLimitWidthV2::{U32 as WidthU32, U64 as WidthU64};
    use RunnerLimitTightenabilityV2::{Fixed, Tightenable};
    use RunnerLimitUnitV2::{
        Artifacts, Classes, Count, DecimalScale, Depth, Diagnostics, Digits, EncodedBytes,
        ExpandedBytes, LogicalBytes, Namespaces, Nodes, Prerequisites, Records, Repairs, Rows,
        Segments, StoredBytes, Visits,
    };
    use RunnerLimitValueV2::{U32, U64};

    const ORACLE: [BaseLimitOracleRowV1; 71] = [
        (
            RunnerLimitFieldV2::ArgvTokens,
            1,
            "argv_tokens",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::ArgvTokenBytes,
            2,
            "argv_token_bytes",
            LogicalBytes,
            WidthU64,
            Tightenable,
            U64(8_192),
            U64(8_192),
        ),
        (
            RunnerLimitFieldV2::ArgvAggregateBytes,
            3,
            "argv_aggregate_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(65_536),
            U64(65_536),
        ),
        (
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            4,
            "lifecycle_record_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(16_384),
            U64(16_384),
        ),
        (
            RunnerLimitFieldV2::CaseLifecycleRecords,
            5,
            "case_lifecycle_records",
            Records,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
            6,
            "case_lifecycle_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(262_144),
            U64(262_144),
        ),
        (
            RunnerLimitFieldV2::FamilyRowsPerCase,
            7,
            "family_rows_per_case",
            Rows,
            WidthU32,
            Tightenable,
            U32(254),
            U32(254),
        ),
        (
            RunnerLimitFieldV2::InvocationCases,
            8,
            "invocation_cases",
            Count,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            9,
            "lifecycle_document_records",
            Records,
            WidthU32,
            Tightenable,
            U32(4_096),
            U32(4_096),
        ),
        (
            RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
            10,
            "lifecycle_document_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(4_194_304),
            U64(4_194_304),
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            11,
            "command_result_stdout_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(5_242_880),
            U64(5_242_880),
        ),
        (
            RunnerLimitFieldV2::ChildStdoutBytes,
            12,
            "child_stdout_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(4_194_304),
            U64(4_194_304),
        ),
        (
            RunnerLimitFieldV2::CombinedChildStdoutBytes,
            13,
            "combined_child_stdout_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(16_777_216),
            U64(134_217_728),
        ),
        (
            RunnerLimitFieldV2::ChildStderrBytes,
            14,
            "child_stderr_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(65_536),
            U64(65_536),
        ),
        (
            RunnerLimitFieldV2::CombinedChildStderrBytes,
            15,
            "combined_child_stderr_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(262_144),
            U64(262_144),
        ),
        (
            RunnerLimitFieldV2::ManifestEncodedBytes,
            16,
            "manifest_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(1_048_576),
            U64(1_048_576),
        ),
        (
            RunnerLimitFieldV2::NestingDepth,
            17,
            "nesting_depth",
            Depth,
            WidthU32,
            Tightenable,
            U32(32),
            U32(32),
        ),
        (
            RunnerLimitFieldV2::ComparisonNodes,
            18,
            "comparison_nodes",
            Nodes,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::EffectNodes,
            19,
            "effect_nodes",
            Nodes,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::TextBytes,
            20,
            "text_bytes",
            LogicalBytes,
            WidthU64,
            Tightenable,
            U64(8_192),
            U64(8_192),
        ),
        (
            RunnerLimitFieldV2::StableTokenBytes,
            21,
            "stable_token_bytes",
            LogicalBytes,
            WidthU64,
            Tightenable,
            U64(128),
            U64(128),
        ),
        (
            RunnerLimitFieldV2::BundleRelativePathBytes,
            22,
            "bundle_relative_path_bytes",
            LogicalBytes,
            WidthU64,
            Tightenable,
            U64(240),
            U64(240),
        ),
        (
            RunnerLimitFieldV2::DiagnosticsPerCase,
            23,
            "diagnostics_per_case",
            Diagnostics,
            WidthU32,
            Tightenable,
            U32(32),
            U32(32),
        ),
        (
            RunnerLimitFieldV2::DiagnosticsPerRun,
            24,
            "diagnostics_per_run",
            Diagnostics,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::PrerequisitesPerDiagnostic,
            25,
            "prerequisites_per_diagnostic",
            Prerequisites,
            WidthU32,
            Tightenable,
            U32(16),
            U32(16),
        ),
        (
            RunnerLimitFieldV2::RepairsPerDiagnostic,
            26,
            "repairs_per_diagnostic",
            Repairs,
            WidthU32,
            Tightenable,
            U32(16),
            U32(16),
        ),
        (
            RunnerLimitFieldV2::Artifacts,
            27,
            "artifacts",
            Artifacts,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::ArtifactEncodedBytes,
            28,
            "artifact_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(67_108_864),
            U64(67_108_864),
        ),
        (
            RunnerLimitFieldV2::ArtifactExpandedBytes,
            29,
            "artifact_expanded_bytes",
            ExpandedBytes,
            WidthU64,
            Tightenable,
            U64(67_108_864),
            U64(67_108_864),
        ),
        (
            RunnerLimitFieldV2::ArtifactStoredBytes,
            30,
            "artifact_stored_bytes",
            StoredBytes,
            WidthU64,
            Tightenable,
            U64(67_112_960),
            U64(67_112_960),
        ),
        (
            RunnerLimitFieldV2::BundleEncodedBytes,
            31,
            "bundle_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(67_108_864),
            U64(536_870_912),
        ),
        (
            RunnerLimitFieldV2::BundleExpandedBytes,
            32,
            "bundle_expanded_bytes",
            ExpandedBytes,
            WidthU64,
            Tightenable,
            U64(67_108_864),
            U64(536_870_912),
        ),
        (
            RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            33,
            "artifact_stored_aggregate_bytes",
            StoredBytes,
            WidthU64,
            Tightenable,
            U64(68_157_440),
            U64(537_919_488),
        ),
        (
            RunnerLimitFieldV2::SystemPublicationStoredBytes,
            34,
            "system_publication_stored_bytes",
            StoredBytes,
            WidthU64,
            Tightenable,
            U64(8_388_608),
            U64(8_388_608),
        ),
        (
            RunnerLimitFieldV2::PublicationStoredBytes,
            35,
            "publication_stored_bytes",
            StoredBytes,
            WidthU64,
            Tightenable,
            U64(76_546_048),
            U64(546_308_096),
        ),
        (
            RunnerLimitFieldV2::ChildStreamDiscardBytes,
            36,
            "child_stream_discard_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(1_048_576),
            U64(1_048_576),
        ),
        (
            RunnerLimitFieldV2::ModesPerFamily,
            37,
            "modes_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::ExtensionDiagnosticsPerFamily,
            38,
            "extension_diagnostics_per_family",
            Diagnostics,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::ArtifactRolesPerFamily,
            39,
            "artifact_roles_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RootPoliciesPerFamily,
            40,
            "root_policies_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RegisteredUnitsPerFamily,
            41,
            "registered_units_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::DigestDomainsPerFamily,
            42,
            "digest_domains_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::ExtensionSchemasPerFamily,
            43,
            "extension_schemas_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::ExecutableDescriptorsPerFamily,
            44,
            "executable_descriptors_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::MapEntries,
            45,
            "map_entries",
            Count,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::GenericArrayItems,
            46,
            "generic_array_items",
            Count,
            WidthU32,
            Tightenable,
            U32(4_096),
            U32(4_096),
        ),
        (
            RunnerLimitFieldV2::PathSegments,
            47,
            "path_segments",
            Segments,
            WidthU32,
            Tightenable,
            U32(32),
            U32(32),
        ),
        (
            RunnerLimitFieldV2::IntegerDigits,
            48,
            "integer_digits",
            Digits,
            WidthU32,
            Fixed,
            U32(39),
            U32(39),
        ),
        (
            RunnerLimitFieldV2::RationalComponentBytes,
            49,
            "rational_component_bytes",
            EncodedBytes,
            WidthU64,
            Fixed,
            U64(16),
            U64(16),
        ),
        (
            RunnerLimitFieldV2::DecimalCoefficientBytes,
            50,
            "decimal_coefficient_bytes",
            EncodedBytes,
            WidthU64,
            Fixed,
            U64(16),
            U64(16),
        ),
        (
            RunnerLimitFieldV2::DecimalAbsoluteScale,
            51,
            "decimal_absolute_scale",
            DecimalScale,
            WidthU32,
            Fixed,
            U32(6_144),
            U32(6_144),
        ),
        (
            RunnerLimitFieldV2::LogicalExtentsPerArtifact,
            52,
            "logical_extents_per_artifact",
            Count,
            WidthU32,
            Tightenable,
            U32(16),
            U32(16),
        ),
        (
            RunnerLimitFieldV2::ObservationKeysPerCase,
            53,
            "observation_keys_per_case",
            Count,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
        (
            RunnerLimitFieldV2::DecisionDetailNamespaces,
            54,
            "decision_detail_namespaces",
            Namespaces,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::OutputClasses,
            55,
            "output_classes",
            Classes,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::OpaqueValueBytes,
            56,
            "opaque_value_bytes",
            LogicalBytes,
            WidthU64,
            Tightenable,
            U64(8_192),
            U64(8_192),
        ),
        (
            RunnerLimitFieldV2::RetainedUnknownExtensionBytes,
            57,
            "retained_unknown_extension_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(65_536),
            U64(65_536),
        ),
        (
            RunnerLimitFieldV2::ExpressionEdges,
            58,
            "expression_edges",
            Count,
            WidthU32,
            Tightenable,
            U32(512),
            U32(512),
        ),
        (
            RunnerLimitFieldV2::MemoizedEvaluationVisits,
            59,
            "memoized_evaluation_visits",
            Visits,
            WidthU32,
            Tightenable,
            U32(4_096),
            U32(4_096),
        ),
        (
            RunnerLimitFieldV2::RepairActionEncodedBytes,
            60,
            "repair_action_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(1_024),
            U64(1_024),
        ),
        (
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            61,
            "actionable_diagnostic_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(8_192),
            U64(8_192),
        ),
        (
            RunnerLimitFieldV2::FailureStderrEncodedBytes,
            62,
            "failure_stderr_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(16_384),
            U64(16_384),
        ),
        (
            RunnerLimitFieldV2::RunnerCatalogEncodedBytes,
            63,
            "runner_catalog_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(1_048_576),
            U64(1_048_576),
        ),
        (
            RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes,
            64,
            "published_bundle_receipt_encoded_bytes",
            EncodedBytes,
            WidthU64,
            Tightenable,
            U64(1_048_576),
            U64(1_048_576),
        ),
        (
            RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes,
            65,
            "content_store_envelope_non_payload_bytes",
            StoredBytes,
            WidthU64,
            Tightenable,
            U64(4_096),
            U64(4_096),
        ),
        (
            RunnerLimitFieldV2::RegisteredExtentAxesPerFamily,
            66,
            "registered_extent_axes_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RegisteredObservationKeysPerFamily,
            67,
            "registered_observation_keys_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(4_096),
            U32(4_096),
        ),
        (
            RunnerLimitFieldV2::RegisteredAuthorityScopesPerFamily,
            68,
            "registered_authority_scopes_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RegisteredExternalRootClassesPerFamily,
            69,
            "registered_external_root_classes_per_family",
            Classes,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RegisteredEvaluationUnitsPerFamily,
            70,
            "registered_evaluation_units_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(64),
            U32(64),
        ),
        (
            RunnerLimitFieldV2::RegisteredResourceIdentitiesPerFamily,
            71,
            "registered_resource_identities_per_family",
            Count,
            WidthU32,
            Tightenable,
            U32(256),
            U32(256),
        ),
    ];
    &ORACLE
}

fn limit_matrix() -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    for profile in [RunProfileV2::Smoke, RunProfileV2::Full] {
        let admitted = RunnerLimitsV2::base(profile);
        for &(field, ordinal, name, unit, width, tightenability, smoke_value, full_value) in
            limit_oracle_rows()
        {
            checked += 1;
            let observed_field = RunnerLimitFieldV2::from_ordinal(ordinal)
                .ok_or_else(|| (checked, format!("limit.{ordinal}.unknown-field")))?;
            let expected_value = match profile {
                RunProfileV2::Smoke => smoke_value,
                RunProfileV2::Full => full_value,
            };
            let descriptor = observed_field.descriptor();
            let base_value = admitted.value(observed_field);
            let descriptor_failure = || format!("limit.{}.{}.descriptor", profile.name(), name);
            if descriptor.field != field
                || descriptor.ordinal != ordinal
                || descriptor.name != name
                || descriptor.unit != unit
                || descriptor.width != width
                || descriptor.tightenability != tightenability
                || field.ordinal() != ordinal
                || base_value != expected_value
            {
                return Err((checked, descriptor_failure()));
            }

            checked += 1;
            let mutation_failure = || format!("limit.{}.{}.one-over", profile.name(), name);
            let mut one_over = admitted.to_candidate();
            let one_over_value = match expected_value {
                RunnerLimitValueV2::U32(value) => RunnerLimitValueV2::U32(
                    value
                        .checked_add(1)
                        .expect("every frozen u32 limit is below u32::MAX"),
                ),
                RunnerLimitValueV2::U64(value) => RunnerLimitValueV2::U64(
                    value
                        .checked_add(1)
                        .expect("every frozen u64 limit is below u64::MAX"),
                ),
            };
            if one_over.set_value(observed_field, one_over_value).is_err() {
                return Err((checked, mutation_failure()));
            }
            let Err(error) = RunnerLimitsV2::admit_family(
                profile,
                one_over,
                RunnerFamilyLimitRequirementsV2::NONE,
            ) else {
                return Err((checked, mutation_failure()));
            };
            let expected_kind = match tightenability {
                RunnerLimitTightenabilityV2::Fixed => {
                    RunnerLimitsViolationKindV2::FixedFieldChanged
                }
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling
                }
            };
            let expected_expectation = match tightenability {
                RunnerLimitTightenabilityV2::Fixed => {
                    RunnerLimitExpectationV2::Exactly(expected_value)
                }
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitExpectationV2::AtMost(expected_value)
                }
            };
            let expected_repair_kind = match tightenability {
                RunnerLimitTightenabilityV2::Fixed => RepairActionKindV2::UpdatePolicyOrCapability,
                RunnerLimitTightenabilityV2::Tightenable => {
                    RepairActionKindV2::ReduceResourceDemand
                }
            };
            if error.kind() != expected_kind
                || error.field() != field
                || error.unit() != unit
                || error.expected() != expected_expectation
                || error.observed() != one_over_value
                || error.owner() != "fs-evidence-runner.runner-limits"
                || error.repair_rank() != 1
                || error.repair_kind() != expected_repair_kind
                || error.repair_target() != name
            {
                return Err((checked, mutation_failure()));
            }
        }
    }
    if checked != 284 {
        return Err((checked, "limit.total-count".to_owned()));
    }
    Ok(checked)
}

type BaseBudgetFieldOracleRowV1 = (
    RunnerBudgetFieldV2,
    u16,
    &'static str,
    RunnerBudgetUnitV2,
    crate::budget::RunnerBudgetWidthV2,
    RunnerBudgetValueV2,
);

#[allow(
    clippy::too_many_lines,
    reason = "the 18-row budget oracle is deliberately handwritten and contiguous so order, width, unit, and fixture values remain auditable"
)]
fn budget_field_oracle_rows() -> &'static [BaseBudgetFieldOracleRowV1] {
    use crate::budget::RunnerBudgetWidthV2::{
        LogicalUnitTaggedSum, U32 as WidthU32, U64 as WidthU64, U128 as WidthU128,
    };
    use RunnerBudgetUnitV2::{
        Count, EncodedBytes, ExpandedBytes, LogicalBytes, LogicalWork, LogicalWorkUnit,
        Nanoseconds, StoredBytes,
    };
    use RunnerBudgetValueV2::{LogicalUnit, U32, U64, U128};
    const FIELD_ORACLE: [BaseBudgetFieldOracleRowV1; 18] = [
        (
            RunnerBudgetFieldV2::WallTimeNs,
            1,
            "wall_time_ns",
            Nanoseconds,
            WidthU64,
            U64(100_000_000_000),
        ),
        (
            RunnerBudgetFieldV2::MaxResidentBytes,
            2,
            "max_resident_bytes",
            LogicalBytes,
            WidthU64,
            U64(1_073_741_824),
        ),
        (
            RunnerBudgetFieldV2::MaxChildProcesses,
            3,
            "max_child_processes",
            Count,
            WidthU32,
            U32(8),
        ),
        (
            RunnerBudgetFieldV2::MaxParallelChildren,
            4,
            "max_parallel_children",
            Count,
            WidthU32,
            U32(4),
        ),
        (
            RunnerBudgetFieldV2::LogicalWorkLimit,
            5,
            "logical_work_limit",
            LogicalWork,
            WidthU128,
            U128(1_000),
        ),
        (
            RunnerBudgetFieldV2::LogicalWorkUnit,
            6,
            "logical_work_unit",
            LogicalWorkUnit,
            LogicalUnitTaggedSum,
            LogicalUnit {
                tag: 11,
                registered_id: None,
            },
        ),
        (
            RunnerBudgetFieldV2::LifecycleEncodedBytes,
            7,
            "lifecycle_encoded_bytes",
            EncodedBytes,
            WidthU64,
            U64(1_000),
        ),
        (
            RunnerBudgetFieldV2::CommandResultStdoutBytes,
            8,
            "command_result_stdout_bytes",
            EncodedBytes,
            WidthU64,
            U64(4_000),
        ),
        (
            RunnerBudgetFieldV2::CombinedChildStdoutBytes,
            9,
            "combined_child_stdout_bytes",
            EncodedBytes,
            WidthU64,
            U64(2_000),
        ),
        (
            RunnerBudgetFieldV2::CombinedChildStderrBytes,
            10,
            "combined_child_stderr_bytes",
            EncodedBytes,
            WidthU64,
            U64(1_000),
        ),
        (
            RunnerBudgetFieldV2::ArtifactEncodedBytes,
            11,
            "artifact_encoded_bytes",
            EncodedBytes,
            WidthU64,
            U64(100),
        ),
        (
            RunnerBudgetFieldV2::ArtifactStoredBytes,
            12,
            "artifact_stored_bytes",
            StoredBytes,
            WidthU64,
            U64(104),
        ),
        (
            RunnerBudgetFieldV2::ArtifactExpandedBytes,
            13,
            "artifact_expanded_bytes",
            ExpandedBytes,
            WidthU64,
            U64(200),
        ),
        (
            RunnerBudgetFieldV2::SystemPublicationStoredBytes,
            14,
            "system_publication_stored_bytes",
            StoredBytes,
            WidthU64,
            U64(72),
        ),
        (
            RunnerBudgetFieldV2::PublicationStoredBytes,
            15,
            "publication_stored_bytes",
            StoredBytes,
            WidthU64,
            U64(176),
        ),
        (
            RunnerBudgetFieldV2::StopObservationNs,
            16,
            "stop_observation_ns",
            Nanoseconds,
            WidthU64,
            U64(1_000_000_000),
        ),
        (
            RunnerBudgetFieldV2::DrainNs,
            17,
            "drain_ns",
            Nanoseconds,
            WidthU64,
            U64(1_000_000_000),
        ),
        (
            RunnerBudgetFieldV2::FinalizeNs,
            18,
            "finalize_ns",
            Nanoseconds,
            WidthU64,
            U64(1_000_000_000),
        ),
    ];
    &FIELD_ORACLE
}

fn budget_logical_unit_oracle_rows() -> &'static [(u16, &'static str, bool)] {
    const LOGICAL_UNIT_ORACLE: [(u16, &str, bool); 16] = [
        (1, "encoded-bytes", false),
        (2, "expanded-bytes", false),
        (3, "stored-bytes", false),
        (4, "logical-bytes", false),
        (5, "count", false),
        (6, "records", false),
        (7, "rows", false),
        (8, "elements", false),
        (9, "samples", false),
        (10, "iterations", false),
        (11, "operations", false),
        (12, "cycles", false),
        (13, "nanoseconds", false),
        (14, "seconds", false),
        (15, "dimensionless", false),
        (16, "registered-unit", true),
    ];
    &LOGICAL_UNIT_ORACLE
}

type BaseBudgetProfileOracleRowV1 = (RunProfileV2, u16, &'static str, u64, u64, u32, u32);

type BaseBudgetProfileRefusalOracleRowV1 = (
    u32,
    RunProfileV2,
    u16,
    &'static str,
    RunnerBudgetFieldV2,
    u16,
    &'static str,
    RunnerBudgetUnitV2,
    RunnerBudgetValueV2,
    RunnerBudgetValueV2,
);

fn budget_profile_oracle_rows() -> &'static [BaseBudgetProfileOracleRowV1] {
    const PROFILE_ORACLE: [BaseBudgetProfileOracleRowV1; 2] = [
        (
            RunProfileV2::Smoke,
            1,
            "smoke",
            900_000_000_000,
            17_179_869_184,
            32,
            256,
        ),
        (
            RunProfileV2::Full,
            2,
            "full",
            86_400_000_000_000,
            137_438_953_472,
            64,
            256,
        ),
    ];
    &PROFILE_ORACLE
}

fn budget_profile_refusal_oracle_rows() -> &'static [BaseBudgetProfileRefusalOracleRowV1] {
    const REFUSAL_ORACLE: [BaseBudgetProfileRefusalOracleRowV1; 8] = [
        (
            37,
            RunProfileV2::Smoke,
            1,
            "smoke",
            RunnerBudgetFieldV2::WallTimeNs,
            1,
            "wall_time_ns",
            RunnerBudgetUnitV2::Nanoseconds,
            RunnerBudgetValueV2::U64(900_000_000_000),
            RunnerBudgetValueV2::U64(900_000_000_001),
        ),
        (
            38,
            RunProfileV2::Smoke,
            1,
            "smoke",
            RunnerBudgetFieldV2::MaxResidentBytes,
            2,
            "max_resident_bytes",
            RunnerBudgetUnitV2::LogicalBytes,
            RunnerBudgetValueV2::U64(17_179_869_184),
            RunnerBudgetValueV2::U64(17_179_869_185),
        ),
        (
            39,
            RunProfileV2::Smoke,
            1,
            "smoke",
            RunnerBudgetFieldV2::MaxParallelChildren,
            4,
            "max_parallel_children",
            RunnerBudgetUnitV2::Count,
            RunnerBudgetValueV2::U32(32),
            RunnerBudgetValueV2::U32(33),
        ),
        (
            40,
            RunProfileV2::Smoke,
            1,
            "smoke",
            RunnerBudgetFieldV2::MaxChildProcesses,
            3,
            "max_child_processes",
            RunnerBudgetUnitV2::Count,
            RunnerBudgetValueV2::U32(256),
            RunnerBudgetValueV2::U32(257),
        ),
        (
            41,
            RunProfileV2::Full,
            2,
            "full",
            RunnerBudgetFieldV2::WallTimeNs,
            1,
            "wall_time_ns",
            RunnerBudgetUnitV2::Nanoseconds,
            RunnerBudgetValueV2::U64(86_400_000_000_000),
            RunnerBudgetValueV2::U64(86_400_000_000_001),
        ),
        (
            42,
            RunProfileV2::Full,
            2,
            "full",
            RunnerBudgetFieldV2::MaxResidentBytes,
            2,
            "max_resident_bytes",
            RunnerBudgetUnitV2::LogicalBytes,
            RunnerBudgetValueV2::U64(137_438_953_472),
            RunnerBudgetValueV2::U64(137_438_953_473),
        ),
        (
            43,
            RunProfileV2::Full,
            2,
            "full",
            RunnerBudgetFieldV2::MaxParallelChildren,
            4,
            "max_parallel_children",
            RunnerBudgetUnitV2::Count,
            RunnerBudgetValueV2::U32(64),
            RunnerBudgetValueV2::U32(65),
        ),
        (
            44,
            RunProfileV2::Full,
            2,
            "full",
            RunnerBudgetFieldV2::MaxChildProcesses,
            3,
            "max_child_processes",
            RunnerBudgetUnitV2::Count,
            RunnerBudgetValueV2::U32(256),
            RunnerBudgetValueV2::U32(257),
        ),
    ];
    &REFUSAL_ORACLE
}

const fn budget_unit_oracle_code(unit: RunnerBudgetUnitV2) -> u16 {
    match unit {
        RunnerBudgetUnitV2::Nanoseconds => 1,
        RunnerBudgetUnitV2::LogicalBytes => 2,
        RunnerBudgetUnitV2::Count => 3,
        RunnerBudgetUnitV2::LogicalWork => 4,
        RunnerBudgetUnitV2::LogicalWorkUnit => 5,
        RunnerBudgetUnitV2::EncodedBytes => 6,
        RunnerBudgetUnitV2::StoredBytes => 7,
        RunnerBudgetUnitV2::ExpandedBytes => 8,
    }
}

const fn budget_field_oracle_identity_code(field: RunnerBudgetFieldV2) -> u16 {
    match field {
        RunnerBudgetFieldV2::WallTimeNs => 1,
        RunnerBudgetFieldV2::MaxResidentBytes => 2,
        RunnerBudgetFieldV2::MaxChildProcesses => 3,
        RunnerBudgetFieldV2::MaxParallelChildren => 4,
        RunnerBudgetFieldV2::LogicalWorkLimit => 5,
        RunnerBudgetFieldV2::LogicalWorkUnit => 6,
        RunnerBudgetFieldV2::LifecycleEncodedBytes => 7,
        RunnerBudgetFieldV2::CommandResultStdoutBytes => 8,
        RunnerBudgetFieldV2::CombinedChildStdoutBytes => 9,
        RunnerBudgetFieldV2::CombinedChildStderrBytes => 10,
        RunnerBudgetFieldV2::ArtifactEncodedBytes => 11,
        RunnerBudgetFieldV2::ArtifactStoredBytes => 12,
        RunnerBudgetFieldV2::ArtifactExpandedBytes => 13,
        RunnerBudgetFieldV2::SystemPublicationStoredBytes => 14,
        RunnerBudgetFieldV2::PublicationStoredBytes => 15,
        RunnerBudgetFieldV2::StopObservationNs => 16,
        RunnerBudgetFieldV2::DrainNs => 17,
        RunnerBudgetFieldV2::FinalizeNs => 18,
    }
}

const fn budget_profile_oracle_identity_code(profile: RunProfileV2) -> u16 {
    match profile {
        RunProfileV2::Smoke => 1,
        RunProfileV2::Full => 2,
    }
}

const fn budget_width_oracle_code(width: crate::budget::RunnerBudgetWidthV2) -> u16 {
    match width {
        crate::budget::RunnerBudgetWidthV2::U32 => 1,
        crate::budget::RunnerBudgetWidthV2::U64 => 2,
        crate::budget::RunnerBudgetWidthV2::U128 => 3,
        crate::budget::RunnerBudgetWidthV2::LogicalUnitTaggedSum => 4,
    }
}

fn push_budget_oracle_value(bytes: &mut Vec<u8>, value: RunnerBudgetValueV2) {
    match value {
        RunnerBudgetValueV2::U32(value) => {
            bytes.extend_from_slice(&1_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::U64(value) => {
            bytes.extend_from_slice(&2_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::U128(value) => {
            bytes.extend_from_slice(&3_u16.to_be_bytes());
            bytes.extend_from_slice(&value.to_be_bytes());
        }
        RunnerBudgetValueV2::LogicalUnit { tag, registered_id } => {
            bytes.extend_from_slice(&4_u16.to_be_bytes());
            bytes.extend_from_slice(&tag.to_be_bytes());
            match registered_id {
                Some(registered_id) => {
                    bytes.push(1);
                    bytes.extend_from_slice(&registered_id.to_be_bytes());
                }
                None => bytes.push(0),
            }
        }
    }
}

fn budget_oracle_table_root_from_rows(
    fields: &[BaseBudgetFieldOracleRowV1],
    logical_units: &[(u16, &'static str, bool)],
    profiles: &[BaseBudgetProfileOracleRowV1],
    refusals: &[BaseBudgetProfileRefusalOracleRowV1],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(4096);
    bytes.extend_from_slice(
        &u32::try_from(fields.len())
            .expect("the literal budget field count fits u32")
            .to_be_bytes(),
    );
    for &(field, ordinal, name, unit, width, value) in fields {
        bytes.extend_from_slice(&budget_field_oracle_identity_code(field).to_be_bytes());
        bytes.extend_from_slice(&ordinal.to_be_bytes());
        detail_push_str(&mut bytes, name);
        bytes.extend_from_slice(&budget_unit_oracle_code(unit).to_be_bytes());
        bytes.extend_from_slice(&budget_width_oracle_code(width).to_be_bytes());
        push_budget_oracle_value(&mut bytes, value);
    }
    bytes.extend_from_slice(
        &u32::try_from(logical_units.len())
            .expect("the literal logical-unit count fits u32")
            .to_be_bytes(),
    );
    for &(tag, name, requires_registered_id) in logical_units {
        bytes.extend_from_slice(&tag.to_be_bytes());
        detail_push_str(&mut bytes, name);
        detail_push_bool(&mut bytes, requires_registered_id);
    }
    bytes.extend_from_slice(
        &u32::try_from(profiles.len())
            .expect("the literal budget profile count fits u32")
            .to_be_bytes(),
    );
    for &(
        profile,
        code,
        name,
        wall_time_ns,
        max_resident_bytes,
        max_parallel_children,
        max_child_processes,
    ) in profiles
    {
        bytes.extend_from_slice(&budget_profile_oracle_identity_code(profile).to_be_bytes());
        bytes.extend_from_slice(&code.to_be_bytes());
        detail_push_str(&mut bytes, name);
        bytes.extend_from_slice(&wall_time_ns.to_be_bytes());
        bytes.extend_from_slice(&max_resident_bytes.to_be_bytes());
        bytes.extend_from_slice(&max_parallel_children.to_be_bytes());
        bytes.extend_from_slice(&max_child_processes.to_be_bytes());
    }
    bytes.extend_from_slice(
        &u32::try_from(refusals.len())
            .expect("the literal budget profile-refusal count fits u32")
            .to_be_bytes(),
    );
    for &(
        semantic_ordinal,
        profile,
        profile_code,
        profile_name,
        field,
        field_ordinal,
        field_name,
        unit,
        ceiling,
        one_over,
    ) in refusals
    {
        bytes.extend_from_slice(&semantic_ordinal.to_be_bytes());
        bytes.extend_from_slice(&budget_profile_oracle_identity_code(profile).to_be_bytes());
        bytes.extend_from_slice(&profile_code.to_be_bytes());
        detail_push_str(&mut bytes, profile_name);
        bytes.extend_from_slice(&budget_field_oracle_identity_code(field).to_be_bytes());
        bytes.extend_from_slice(&field_ordinal.to_be_bytes());
        detail_push_str(&mut bytes, field_name);
        bytes.extend_from_slice(&budget_unit_oracle_code(unit).to_be_bytes());
        push_budget_oracle_value(&mut bytes, ceiling);
        push_budget_oracle_value(&mut bytes, one_over);
    }
    hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-budget-literal-oracle.v1",
        &bytes,
    )
}

fn budget_oracle_table_root() -> ContentHash {
    budget_oracle_table_root_from_rows(
        budget_field_oracle_rows(),
        budget_logical_unit_oracle_rows(),
        budget_profile_oracle_rows(),
        budget_profile_refusal_oracle_rows(),
    )
}

fn budget_registered_unit_registry()
-> Result<BaseExtensionRegistryProjectionV2, ConstructionErrorV2> {
    let registered_unit = LogicalUnitV2::from_tag(16, Some(7)).map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "budget.registry.logical_unit",
            "the exact registered-unit fixture ID 7",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })?;
    let no_claim_scope = NoClaimScopeRootV1::parse_presented(
        NoClaimScopeRootV1::DESCRIPTOR.role(),
        NoClaimScopeRootV1::DESCRIPTOR.domain(),
        &"00".repeat(32),
    )
    .map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "budget.registry.no_claim_scope",
            "the exact all-zero-present nominal no-claim fixture",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })?;
    let descriptor = RegisteredLogicalUnitDescriptorV2::new(
        registered_unit,
        StableTokenV2::new("org.frankensim.fixture.logical-work-unit").map_err(|_error| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "budget.registry.logical_unit.name",
                "the exact namespaced registered-unit fixture name",
                ConstructionObservedDataClassV2::CallerControlledText,
            )
        })?,
        StableTokenV2::new("org.frankensim.fixture.owner").map_err(|_error| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "budget.registry.logical_unit.owner",
                "the exact registered-unit fixture owner",
                ConstructionObservedDataClassV2::CallerControlledText,
            )
        })?,
        no_claim_scope,
    )?;
    BaseExtensionRegistryProjectionV2::try_new(
        &RunnerLimitsV2::base(RunProfileV2::Smoke),
        &[],
        &[descriptor],
        &[],
    )
}

#[allow(
    clippy::too_many_lines,
    reason = "the budget matrix exhaustively reconciles the frozen field, unit, profile, and one-over oracle partitions"
)]
fn budget_matrix() -> Result<u32, (u32, String)> {
    let base = RunnerBudgetsV2::try_new(durable_budget_candidate()).map_err(|_| {
        (
            BUDGET_BASE_CONSTRUCTION_SEMANTIC_ORDINAL_V1,
            "budget.base-construction".to_owned(),
        )
    })?;
    let base_root = base.semantic_root();
    let mut checked = 0_u32;

    for &(field, ordinal, name, unit, width, expected_value) in budget_field_oracle_rows() {
        checked += 1;
        let observed_field = RunnerBudgetFieldV2::from_ordinal(ordinal)
            .ok_or_else(|| (checked, format!("budget.field.{ordinal}")))?;
        let descriptor = observed_field.descriptor();
        let failure = || format!("budget.field.{name}");
        if descriptor.field != field
            || descriptor.ordinal != ordinal
            || descriptor.name != name
            || descriptor.unit != unit
            || descriptor.width != width
            || field.ordinal() != ordinal
            || base.value(field) != expected_value
        {
            return Err((checked, failure()));
        }
        let mutated = RunnerBudgetsV2::try_new(mutated_budget_candidate(observed_field))
            .map_err(|_| (checked, failure()))?;
        if mutated.value(field) == base.value(field)
            || mutated.semantic_root().bytes() == base_root.bytes()
        {
            return Err((checked, failure()));
        }
    }

    for &(tag, name, requires_registered_id) in budget_logical_unit_oracle_rows() {
        checked += 1;
        let unit = LogicalUnitV2::from_tag(tag, requires_registered_id.then_some(7))
            .map_err(|_| (checked, format!("budget.logical-unit.{tag}")))?;
        if unit.tag() != tag
            || unit.name() != name
            || unit.registered_id() != requires_registered_id.then_some(7)
        {
            return Err((checked, format!("budget.logical-unit.{tag}")));
        }
        let mut candidate = durable_budget_candidate();
        candidate.logical_work_unit = unit;
        let unit_matches = if requires_registered_id {
            let registry = budget_registered_unit_registry()
                .map_err(|_| (checked, format!("budget.logical-unit.{tag}.registry")))?;
            RunnerBudgetsV2::try_new_with_extension_registry(candidate, &registry)
                .map(|bound| bound.budgets().logical_work_unit() == unit)
        } else {
            RunnerBudgetsV2::try_new(candidate).map(|budgets| budgets.logical_work_unit() == unit)
        }
        .map_err(|_| (checked, format!("budget.logical-unit.{tag}")))?;
        if !unit_matches {
            return Err((checked, format!("budget.logical-unit.{tag}")));
        }
    }

    for &(
        profile,
        profile_code,
        profile_name,
        wall_time_ns,
        max_resident_bytes,
        max_parallel_children,
        max_child_processes,
    ) in budget_profile_oracle_rows()
    {
        checked += 1;
        let mut candidate = durable_budget_candidate();
        candidate.wall_time_ns = wall_time_ns;
        candidate.max_resident_bytes = max_resident_bytes;
        candidate.max_parallel_children = max_parallel_children;
        candidate.max_child_processes = max_child_processes;
        let exact = RunnerBudgetsV2::try_new(candidate)
            .map_err(|_| {
                (
                    checked,
                    format!("budget.admission.{profile_name}.intrinsic"),
                )
            })?
            .admit(
                profile,
                ArtifactDispositionV2::DurableBundleRequired,
                &RunnerLimitsV2::base(profile),
            )
            .map_err(|_| (checked, format!("budget.admission.{profile_name}.exact")))?;
        if profile.code() != profile_code
            || profile.name() != profile_name
            || exact.budgets().wall_time_ns() != wall_time_ns
            || exact.budgets().max_resident_bytes() != max_resident_bytes
            || exact.budgets().max_parallel_children() != max_parallel_children
            || exact.budgets().max_child_processes() != max_child_processes
        {
            return Err((checked, format!("budget.admission.{profile_name}")));
        }
    }

    for &(
        _semantic_ordinal,
        profile,
        _profile_code,
        profile_name,
        field,
        _field_ordinal,
        field_name,
        unit,
        ceiling,
        one_over,
    ) in budget_profile_refusal_oracle_rows()
    {
        checked += 1;
        let (decision, observed) = observe_budget_profile_refusal(profile, field, one_over);
        let expected = BaseE2eDetailPayloadV1::Budget {
            kind: RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
            field,
            unit,
            expected: RunnerBudgetExpectationV2::AtMost(ceiling),
            observed: one_over,
            owner: "fs-evidence-runner.runner-budgets",
            repair_rank: 1,
            repair_kind: RepairActionKindV2::ReduceResourceDemand,
            repair_target: field_name,
        };
        if decision != BaseE2eExpectedDecisionV1::Refuse || observed != expected {
            return Err((
                checked,
                format!("budget.profile.{profile_name}.{field_name}.one-over"),
            ));
        }
    }

    if checked != 44 {
        return Err((checked, "budget.total-count".to_owned()));
    }
    Ok(checked)
}

fn mutated_budget_candidate(field: RunnerBudgetFieldV2) -> RunnerBudgetsCandidateV2 {
    let mut candidate = durable_budget_candidate();
    match field {
        RunnerBudgetFieldV2::WallTimeNs => candidate.wall_time_ns += 1,
        RunnerBudgetFieldV2::MaxResidentBytes => candidate.max_resident_bytes += 1,
        RunnerBudgetFieldV2::MaxChildProcesses => candidate.max_child_processes += 1,
        RunnerBudgetFieldV2::MaxParallelChildren => candidate.max_parallel_children += 1,
        RunnerBudgetFieldV2::LogicalWorkLimit => candidate.logical_work_limit += 1,
        RunnerBudgetFieldV2::LogicalWorkUnit => {
            candidate.logical_work_unit = LogicalUnitV2::Cycles;
        }
        RunnerBudgetFieldV2::LifecycleEncodedBytes => candidate.lifecycle_encoded_bytes += 1,
        RunnerBudgetFieldV2::CommandResultStdoutBytes => {
            candidate.command_result_stdout_bytes += 1;
        }
        RunnerBudgetFieldV2::CombinedChildStdoutBytes => {
            candidate.combined_child_stdout_bytes += 1;
        }
        RunnerBudgetFieldV2::CombinedChildStderrBytes => {
            candidate.combined_child_stderr_bytes += 1;
        }
        RunnerBudgetFieldV2::ArtifactEncodedBytes => candidate.artifact_encoded_bytes += 1,
        RunnerBudgetFieldV2::ArtifactStoredBytes => candidate.artifact_stored_bytes += 1,
        RunnerBudgetFieldV2::ArtifactExpandedBytes => candidate.artifact_expanded_bytes += 1,
        RunnerBudgetFieldV2::SystemPublicationStoredBytes => {
            candidate.system_publication_stored_bytes += 1;
        }
        RunnerBudgetFieldV2::PublicationStoredBytes => candidate.publication_stored_bytes += 1,
        RunnerBudgetFieldV2::StopObservationNs => candidate.stop_observation_ns += 1,
        RunnerBudgetFieldV2::DrainNs => candidate.drain_ns += 1,
        RunnerBudgetFieldV2::FinalizeNs => candidate.finalize_ns += 1,
    }
    candidate
}

fn durable_budget_candidate() -> RunnerBudgetsCandidateV2 {
    RunnerBudgetsCandidateV2 {
        wall_time_ns: 100_000_000_000,
        max_resident_bytes: 1024 * 1024 * 1024,
        max_child_processes: 8,
        max_parallel_children: 4,
        logical_work_limit: 1000,
        logical_work_unit: crate::catalog::LogicalUnitV2::Operations,
        lifecycle_encoded_bytes: 1000,
        command_result_stdout_bytes: 4000,
        combined_child_stdout_bytes: 2000,
        combined_child_stderr_bytes: 1000,
        artifact_encoded_bytes: 100,
        artifact_stored_bytes: 104,
        artifact_expanded_bytes: 200,
        system_publication_stored_bytes: 72,
        publication_stored_bytes: 176,
        stop_observation_ns: 1_000_000_000,
        drain_ns: 1_000_000_000,
        finalize_ns: 1_000_000_000,
    }
}

fn publication_selection_matrix() -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    let mut roots = Vec::new();
    for profile in PlatformPathProfileV2::ALL {
        for mode in DestinationAdmissionModeV2::ALL {
            checked += 1;
            let selection = selection_for_profile(profile, mode).map_err(|_| {
                (
                    checked,
                    format!("publication.{}.{}", profile.name(), mode.name()),
                )
            })?;
            if selection.path_profile() != profile || selection.destination_mode() != mode {
                return Err((
                    checked,
                    format!("publication.{}.{}", profile.name(), mode.name()),
                ));
            }
            let root = selection.semantic_projection_root();
            if roots.contains(&root) {
                return Err((
                    checked,
                    format!("publication.{}.{}.root", profile.name(), mode.name()),
                ));
            }
            roots.push(root);
        }
    }
    if checked != 6 {
        return Err((checked, "publication.total-count".to_owned()));
    }
    Ok(checked)
}

fn selection_for_profile(
    profile: PlatformPathProfileV2,
    mode: DestinationAdmissionModeV2,
) -> Result<PublicationSelectionV2, ConstructionErrorV2> {
    let logical = LogicalBundlePathV1::new("runner/seal").map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.publication_path",
            "a valid logical path",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })?;
    let (protocol, target) = match profile {
        PlatformPathProfileV2::PosixDescriptorRelativeV1 => (
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            PublicationTargetV2::PosixRelative(logical),
        ),
        PlatformPathProfileV2::WindowsHandleRelativeV1 => (
            PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1,
            PublicationTargetV2::WindowsRelative(logical),
        ),
        PlatformPathProfileV2::ContentStoreObjectKeyV1 => (
            PublicationProtocolV2::ContentStoreAtomicCommitV1,
            PublicationTargetV2::ContentStoreLogicalKey(
                ContentStoreObjectKeyV1::new("objects/seal").map_err(|_error| {
                    ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Incompatible,
                        "projection.content_store_key",
                        "a valid logical object key",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    )
                })?,
            ),
        ),
    };
    PublicationSelectionV2::new(profile, protocol, mode, target)
}

fn publication_selection(path: &str) -> Result<PublicationSelectionV2, ConstructionErrorV2> {
    let path = LogicalBundlePathV1::new(path).map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.publication_path",
            "a valid logical path",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })?;
    PublicationSelectionV2::new(
        PlatformPathProfileV2::PosixDescriptorRelativeV1,
        PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
        DestinationAdmissionModeV2::Absent,
        PublicationTargetV2::PosixRelative(path),
    )
}

fn capability_registry() -> Result<RootPolicyRegistryProjectionV2, ConstructionErrorV2> {
    RootPolicyRegistryProjectionV2::new(
        vec![1],
        vec![1],
        vec![OverlapPolicyRegistrationV2::new(
            1,
            OverlapPolicyRelationV2::RequireInputOutputDisjoint,
        )?],
    )
}

fn capability_valid_matrix(no_claim_scope: &NoClaimScopeRootV1) -> Result<u32, (u32, String)> {
    let registry = capability_registry().map_err(|_| {
        (
            CAPABILITY_REGISTRY_SETUP_SEMANTIC_ORDINAL_V1,
            "capability.registry".to_owned(),
        )
    })?;
    let mut checked = 0_u32;
    for profile in PlatformPathProfileV2::ALL {
        for access in RootCapabilityAccessV2::ALL {
            for mode in DestinationAdmissionModeV2::ALL {
                checked += 1;
                let failure = || {
                    format!(
                        "capability.valid.{}.{}.{}",
                        profile.name(),
                        access.name(),
                        mode.name()
                    )
                };
                let oracle_rights = capability_oracle_rights(profile, access, mode);
                if expected_rights(profile, access, mode) != oracle_rights {
                    return Err((checked, format!("{}.oracle", failure())));
                }
                let policy = RootCapabilityPolicyV2::new(
                    root_class_for_access(access),
                    profile,
                    access,
                    oracle_rights,
                    1,
                    1,
                    1,
                    no_claim_scope.clone(),
                )
                .map_err(|_| (checked, failure()))?;
                policy
                    .validate_registration(&registry)
                    .map_err(|_| (checked, failure()))?;
                match access {
                    RootCapabilityAccessV2::ReadOnlyInput => {
                        let view = NarrowedPolicyViewV2::for_read_only(&policy)
                            .map_err(|_| (checked, failure()))?;
                        if view.rights() != policy.rights()
                            || view.destination_mode().is_some()
                            || view.path_profile() != profile
                        {
                            return Err((checked, failure()));
                        }
                    }
                    RootCapabilityAccessV2::DurableOutput => {
                        let selection = selection_for_profile(profile, mode)
                            .map_err(|_| (checked, failure()))?;
                        let view = NarrowedPolicyViewV2::for_publication(&policy, &selection)
                            .map_err(|_| (checked, failure()))?;
                        if view.rights() != policy.rights()
                            || view.destination_mode() != Some(mode)
                            || view.path_profile() != profile
                        {
                            return Err((checked, failure()));
                        }
                    }
                }
            }
        }
    }
    if checked != 12 {
        return Err((checked, "capability.valid.total-count".to_owned()));
    }
    Ok(checked)
}

fn capability_invalid_matrix(no_claim_scope: &NoClaimScopeRootV1) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    for profile in PlatformPathProfileV2::ALL {
        for access in RootCapabilityAccessV2::ALL {
            for mode in DestinationAdmissionModeV2::ALL {
                let exact = capability_oracle_rights(profile, access, mode);
                if expected_rights(profile, access, mode) != exact {
                    let failed_ordinal = checked
                        .checked_add(1)
                        .expect("the frozen capability matrix fits u32");
                    return Err((
                        failed_ordinal,
                        format!(
                            "capability.oracle.{}.{}.{}",
                            profile.name(),
                            access.name(),
                            mode.name()
                        ),
                    ));
                }
                for right in RootCapabilityRightV2::ALL {
                    let mut mutant = exact.clone();
                    if let Some(index) = mutant.iter().position(|candidate| *candidate == right) {
                        mutant.remove(index);
                    } else {
                        mutant.push(right);
                    }
                    mutant.sort_unstable_by_key(|candidate| candidate.code());
                    checked += 1;
                    if !capability_mutant_refuses(
                        profile,
                        access,
                        mode,
                        mutant,
                        no_claim_scope.clone(),
                    ) {
                        return Err((
                            checked,
                            format!(
                                "capability.mutant.{}.{}.{}.{}",
                                profile.name(),
                                access.name(),
                                mode.name(),
                                right.code()
                            ),
                        ));
                    }
                }
                for removed in exact.iter().copied() {
                    for replacement in RootCapabilityRightV2::ALL
                        .into_iter()
                        .filter(|candidate| !exact.contains(candidate))
                    {
                        let mut mutant = exact.clone();
                        let index = mutant
                            .iter()
                            .position(|candidate| *candidate == removed)
                            .expect("the removed right comes from the exact set");
                        mutant[index] = replacement;
                        mutant.sort_unstable_by_key(|candidate| candidate.code());
                        checked += 1;
                        if !capability_mutant_refuses(
                            profile,
                            access,
                            mode,
                            mutant,
                            no_claim_scope.clone(),
                        ) {
                            return Err((
                                checked,
                                format!(
                                    "capability.substitution.{}.{}.{}.{}.{}",
                                    profile.name(),
                                    access.name(),
                                    mode.name(),
                                    removed.code(),
                                    replacement.code()
                                ),
                            ));
                        }
                    }
                }
            }
        }
    }
    if checked != 390 {
        return Err((checked, "capability.invalid.total-count".to_owned()));
    }
    Ok(checked)
}

fn capability_mutant_refuses(
    profile: PlatformPathProfileV2,
    access: RootCapabilityAccessV2,
    mode: DestinationAdmissionModeV2,
    rights: Vec<RootCapabilityRightV2>,
    no_claim_scope: NoClaimScopeRootV1,
) -> bool {
    let (stage, error) = expected_capability_refusal(profile, access, &rights);
    observed_capability_refusal(profile, access, mode, rights, no_claim_scope)
        == (
            BaseE2eExpectedDecisionV1::Refuse,
            BaseE2eDetailPayloadV1::Capability { stage, error },
        )
}

const fn root_class_for_access(access: RootCapabilityAccessV2) -> RootClassV2 {
    match access {
        RootCapabilityAccessV2::ReadOnlyInput => RootClassV2::InputArtifactRoot,
        RootCapabilityAccessV2::DurableOutput => RootClassV2::OutputArtifactRoot,
    }
}

#[allow(
    clippy::too_many_lines,
    reason = "the state and NotRun matrix exhaustively enumerates the closed role-state-reason-diagnostic-drain cross-product"
)]
fn state_and_not_run_matrix() -> Result<u32, (u32, String)> {
    let reasons = core::iter::once(None)
        .chain(RefusedReasonV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let diagnostics = core::iter::once(None)
        .chain(DiagnosticCodeV2::ALL.into_iter().map(Some))
        .collect::<Vec<_>>();
    let drains = [
        None,
        Some(PresentedDrainRootKindV2::CancelledStopRoot),
        Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
    ];
    let mut checked = 0_u32;
    let mut positive_cells = 0_u32;
    let mut expected_refusal_cells = 0_u32;
    for role in StateBearingRecordRoleV2::ALL {
        for state in ProofExitV2::ALL {
            for reason in &reasons {
                for diagnostic in &diagnostics {
                    for drain in drains {
                        checked += 1;
                        let observed = validate_state_v2(StateValidationInputV2::new(
                            role,
                            state,
                            *reason,
                            *diagnostic,
                            drain,
                        ))
                        .map(|_| ());
                        let expected =
                            expected_state_validation(role, state, *reason, *diagnostic, drain);
                        debug_assert_eq!(
                            expected.is_ok(),
                            expected_state_cell(role, state, *reason, *diagnostic, drain)
                        );
                        if expected.is_ok() {
                            positive_cells += 1;
                        } else {
                            expected_refusal_cells += 1;
                        }
                        if observed != expected {
                            return Err((
                                checked,
                                format!(
                                    "state.matrix.{}.{}.{}.{}.{}",
                                    role.code(),
                                    state.code(),
                                    reason.map_or(0, RefusedReasonV2::code),
                                    diagnostic.map_or(0, DiagnosticCodeV2::code),
                                    drain.map_or(0, drain_code)
                                ),
                            ));
                        }
                    }
                }
            }
        }
    }

    let (cancelled, timed_out, internal_error) = presented_stop_fixture_roots().map_err(|_| {
        (
            checked
                .checked_add(1)
                .expect("the not-run fixture cell fits u32"),
            "not-run.fixture-roots".to_owned(),
        )
    })?;
    let causes = [
        NotRunCauseV2::PriorCancelled(cancelled),
        NotRunCauseV2::PriorTimedOut(timed_out),
        NotRunCauseV2::PriorControlledInternalError(internal_error),
    ];
    for cause in causes {
        let code = cause.code();
        checked += 1;
        positive_cells += 1;
        let first = NotRunBasisV2::new(cause.clone(), 0, 1)
            .map_err(|_| (checked, format!("not-run.{code}.first")))?;
        if first.remaining_case_count(1) != Ok(1)
            || first.diagnostic() != DiagnosticCodeV2::RunnerNotRun
            || first.state() != ProofExitV2::NotRun
        {
            return Err((checked, format!("not-run.{code}.first")));
        }

        checked += 1;
        positive_cells += 1;
        let last = NotRunBasisV2::new(cause.clone(), 255, 256)
            .map_err(|_| (checked, format!("not-run.{code}.last")))?;
        if last.remaining_case_count(256) != Ok(1) {
            return Err((checked, format!("not-run.{code}.last")));
        }

        checked += 1;
        expected_refusal_cells += 1;
        if NotRunBasisV2::new(cause.clone(), 256, 256)
            != Err(NotRunBasisErrorV2::LowestRemainingOrdinalOutOfRange {
                observed: 256,
                ordered_case_count: 256,
            })
        {
            return Err((checked, format!("not-run.{code}.one-over")));
        }

        checked += 1;
        expected_refusal_cells += 1;
        if NotRunBasisV2::new(cause, 0, 0) != Err(NotRunBasisErrorV2::EmptyManifest) {
            return Err((checked, format!("not-run.{code}.empty")));
        }
    }

    if checked != 32_460 || positive_cells != 69 || expected_refusal_cells != 32_391 {
        return Err((checked, "state.total-count".to_owned()));
    }
    Ok(checked)
}

fn presented_stop_fixture_roots() -> Result<
    (
        CancelledStopRootV2,
        TimedOutStopRootV2,
        DrainedInternalErrorRootV2,
    ),
    ConstructionErrorV2,
> {
    let parse_error = |field: &'static str| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            field,
            "a nominal presented stop fixture root",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    };
    let cancelled_hex = hash_domain(
        "org.frankensim.fs-evidence-runner.fixture.cancelled-stop.v1",
        b"cancelled-stop",
    )
    .to_hex();
    let timed_out_hex = hash_domain(
        "org.frankensim.fs-evidence-runner.fixture.timed-out-stop.v1",
        b"timed-out-stop",
    )
    .to_hex();
    let internal_error_hex = hash_domain(
        "org.frankensim.fs-evidence-runner.fixture.internal-error-stop.v1",
        b"internal-error-stop",
    )
    .to_hex();
    Ok((
        CancelledStopRootV2::parse_presented(
            CancelledStopRootV2::DESCRIPTOR.role(),
            CancelledStopRootV2::DESCRIPTOR.domain(),
            &cancelled_hex,
        )
        .map_err(|_| parse_error("fixture.cancelled_stop_root"))?,
        TimedOutStopRootV2::parse_presented(
            TimedOutStopRootV2::DESCRIPTOR.role(),
            TimedOutStopRootV2::DESCRIPTOR.domain(),
            &timed_out_hex,
        )
        .map_err(|_| parse_error("fixture.timed_out_stop_root"))?,
        DrainedInternalErrorRootV2::parse_presented(
            DrainedInternalErrorRootV2::DESCRIPTOR.role(),
            DrainedInternalErrorRootV2::DESCRIPTOR.domain(),
            &internal_error_hex,
        )
        .map_err(|_| parse_error("fixture.internal_error_root"))?,
    ))
}

fn expected_state_cell(
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
    reason: Option<RefusedReasonV2>,
    diagnostic: Option<DiagnosticCodeV2>,
    drain: Option<PresentedDrainRootKindV2>,
) -> bool {
    let state_allowed = match role {
        StateBearingRecordRoleV2::PreRunDiagnostic => matches!(
            state,
            ProofExitV2::Usage
                | ProofExitV2::Refused
                | ProofExitV2::NoData
                | ProofExitV2::Stale
                | ProofExitV2::EnvironmentInvalid
                | ProofExitV2::Blocked
                | ProofExitV2::Unsupported
                | ProofExitV2::Cancelled
                | ProofExitV2::TimedOut
                | ProofExitV2::InternalError
        ),
        StateBearingRecordRoleV2::ExecutedCaseTerminal | StateBearingRecordRoleV2::RunTerminal => {
            !matches!(state, ProofExitV2::Usage | ProofExitV2::NotRun)
        }
        StateBearingRecordRoleV2::SuppressedCaseTerminal => state == ProofExitV2::NotRun,
    };
    if !state_allowed || (state == ProofExitV2::Refused) != reason.is_some() {
        return false;
    }
    if diagnostic != expected_diagnostic(state) {
        return false;
    }
    drain == expected_drain(role, state)
}

const fn expected_diagnostic(state: ProofExitV2) -> Option<DiagnosticCodeV2> {
    match state {
        ProofExitV2::Pass => None,
        ProofExitV2::Failed => Some(DiagnosticCodeV2::CaseConformanceMismatch),
        ProofExitV2::Refused => Some(DiagnosticCodeV2::RunnerRefused),
        ProofExitV2::NoData => Some(DiagnosticCodeV2::RunnerNoData),
        ProofExitV2::Stale => Some(DiagnosticCodeV2::RunnerStale),
        ProofExitV2::EnvironmentInvalid => Some(DiagnosticCodeV2::RunnerEnvironmentInvalid),
        ProofExitV2::Blocked => Some(DiagnosticCodeV2::RunnerBlocked),
        ProofExitV2::Unsupported => Some(DiagnosticCodeV2::RunnerUnsupported),
        ProofExitV2::NotRun => Some(DiagnosticCodeV2::RunnerNotRun),
        ProofExitV2::Cancelled => Some(DiagnosticCodeV2::RunnerCancelled),
        ProofExitV2::TimedOut => Some(DiagnosticCodeV2::RunnerTimedOut),
        ProofExitV2::Usage => Some(DiagnosticCodeV2::RunnerUsage),
        ProofExitV2::InternalError => Some(DiagnosticCodeV2::RunnerInternalError),
    }
}

const fn expected_drain(
    role: StateBearingRecordRoleV2,
    state: ProofExitV2,
) -> Option<PresentedDrainRootKindV2> {
    if matches!(role, StateBearingRecordRoleV2::PreRunDiagnostic) {
        return None;
    }
    match state {
        ProofExitV2::Cancelled => Some(PresentedDrainRootKindV2::CancelledStopRoot),
        ProofExitV2::TimedOut => Some(PresentedDrainRootKindV2::TimedOutStopRoot),
        ProofExitV2::InternalError => Some(PresentedDrainRootKindV2::DrainedInternalErrorRoot),
        _ => None,
    }
}

const fn drain_code(drain: PresentedDrainRootKindV2) -> u16 {
    match drain {
        PresentedDrainRootKindV2::CancelledStopRoot => 1,
        PresentedDrainRootKindV2::TimedOutStopRoot => 2,
        PresentedDrainRootKindV2::DrainedInternalErrorRoot => 3,
    }
}

fn diagnostic_code_oracle_rows() -> &'static [(u16, &'static str)] {
    const CODE_ORACLE: [(u16, &str); 12] = [
        (1, "case.conformance_mismatch"),
        (2, "runner.not_run"),
        (3, "runner.refused"),
        (4, "runner.no_data"),
        (5, "runner.stale"),
        (6, "runner.environment_invalid"),
        (7, "runner.blocked"),
        (8, "runner.unsupported"),
        (9, "runner.cancelled"),
        (10, "runner.timed_out"),
        (11, "runner.usage"),
        (12, "runner.internal_error"),
    ];
    &CODE_ORACLE
}

const fn diagnostic_registered_code_oracle() -> (u16, u16) {
    (7, 9)
}

fn diagnostic_retryability_oracle_rows() -> &'static [(u16, &'static str)] {
    const RETRYABILITY_ORACLE: [(u16, &str); 5] = [
        (0, "never"),
        (1, "same-invocation"),
        (2, "after-input-change"),
        (3, "after-environment-change"),
        (4, "after-prerequisite-change"),
    ];
    &RETRYABILITY_ORACLE
}

fn diagnostic_repair_kind_oracle_rows() -> &'static [(u16, &'static str)] {
    const REPAIR_KIND_ORACLE: [(u16, &str); 12] = [
        (1, "change-arguments"),
        (2, "supply-evidence"),
        (3, "regenerate-canonical-evidence"),
        (4, "refresh-evidence"),
        (5, "reduce-resource-demand"),
        (6, "choose-safe-artifact-destination"),
        (7, "restore-lifecycle"),
        (8, "update-policy-or-capability"),
        (9, "register-migration"),
        (10, "retry-same-invocation"),
        (11, "contact-owner"),
        (12, "inspect-retained-artifact"),
    ];
    &REPAIR_KIND_ORACLE
}

fn diagnostic_oracle_table_root_from_rows(
    codes: &[(u16, &'static str)],
    registered: (u16, u16),
    retryability: &[(u16, &'static str)],
    repairs: &[(u16, &'static str)],
) -> ContentHash {
    let mut bytes = Vec::with_capacity(1024);
    bytes.extend_from_slice(
        &u32::try_from(codes.len())
            .expect("the literal diagnostic code count fits u32")
            .to_be_bytes(),
    );
    for &(code, name) in codes {
        bytes.extend_from_slice(&code.to_be_bytes());
        detail_push_str(&mut bytes, name);
    }
    let (registered_namespace, registered_code) = registered;
    bytes.extend_from_slice(&registered_namespace.to_be_bytes());
    bytes.extend_from_slice(&registered_code.to_be_bytes());
    bytes.extend_from_slice(
        &u32::try_from(retryability.len())
            .expect("the literal retryability count fits u32")
            .to_be_bytes(),
    );
    for &(code, name) in retryability {
        bytes.extend_from_slice(&code.to_be_bytes());
        detail_push_str(&mut bytes, name);
    }
    bytes.extend_from_slice(
        &u32::try_from(repairs.len())
            .expect("the literal repair-kind count fits u32")
            .to_be_bytes(),
    );
    for &(code, name) in repairs {
        bytes.extend_from_slice(&code.to_be_bytes());
        detail_push_str(&mut bytes, name);
    }
    hash_domain(
        "org.frankensim.fs-evidence-runner.base-e2e-diagnostic-literal-oracle.v1",
        &bytes,
    )
}

fn diagnostic_oracle_table_root() -> ContentHash {
    diagnostic_oracle_table_root_from_rows(
        diagnostic_code_oracle_rows(),
        diagnostic_registered_code_oracle(),
        diagnostic_retryability_oracle_rows(),
        diagnostic_repair_kind_oracle_rows(),
    )
}

fn diagnostic_matrix(no_claim_scope: &NoClaimScopeRootV1) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    for &(code, name) in diagnostic_code_oracle_rows() {
        checked += 1;
        let observed_code = DiagnosticCodeV2::from_code(code)
            .map_err(|_| (checked, format!("diagnostic.code.{code}")))?;
        if observed_code.code() != code || observed_code.name() != name {
            return Err((checked, format!("diagnostic.code.{code}")));
        }
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(observed_code),
            RetryabilityV2::AfterInputChange,
            RepairActionKindV2::ChangeArguments,
            1,
        )
        .map_err(|_| (checked, format!("diagnostic.code.{code}")))?;
        if value.code().code() != code {
            return Err((checked, format!("diagnostic.code.{code}")));
        }
    }

    checked += 1;
    let (registered_namespace, registered_code) = diagnostic_registered_code_oracle();
    let registered = DiagnosticCodeRefV2::registered(registered_namespace, registered_code)
        .map_err(|_| (checked, "diagnostic.registered".to_owned()))?;
    let registered_value = diagnostic_fixture(
        no_claim_scope.clone(),
        registered,
        RetryabilityV2::AfterPrerequisiteChange,
        RepairActionKindV2::ContactOwner,
        1,
    )
    .map_err(|_| (checked, "diagnostic.registered".to_owned()))?;
    if registered_value.code().registered_namespace() != Some(registered_namespace)
        || registered_value.code().code() != registered_code
    {
        return Err((checked, "diagnostic.registered".to_owned()));
    }

    for &(code, name) in diagnostic_retryability_oracle_rows() {
        checked += 1;
        let retryability = RetryabilityV2::from_code(code)
            .map_err(|_| (checked, format!("diagnostic.retryability.{code}")))?;
        if retryability.code() != code || retryability.name() != name {
            return Err((checked, format!("diagnostic.retryability.{code}")));
        }
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            retryability,
            RepairActionKindV2::ChangeArguments,
            1,
        )
        .map_err(|_| (checked, format!("diagnostic.retryability.{code}")))?;
        if value.retryability() != retryability {
            return Err((checked, format!("diagnostic.retryability.{code}")));
        }
    }

    for &(code, name) in diagnostic_repair_kind_oracle_rows() {
        checked += 1;
        let kind = RepairActionKindV2::from_code(code)
            .map_err(|_| (checked, format!("diagnostic.repair-kind.{code}")))?;
        if kind.code() != code || kind.name() != name {
            return Err((checked, format!("diagnostic.repair-kind.{code}")));
        }
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            RetryabilityV2::AfterInputChange,
            kind,
            1,
        )
        .map_err(|_| (checked, format!("diagnostic.repair-kind.{code}")))?;
        if value.repairs()[0].kind() != kind {
            return Err((checked, format!("diagnostic.repair-kind.{code}")));
        }
    }

    if checked != 30 {
        return Err((checked, "diagnostic.total-count".to_owned()));
    }
    Ok(checked)
}

fn diagnostic_fixture(
    no_claim_scope: NoClaimScopeRootV1,
    code: DiagnosticCodeRefV2,
    retryability: RetryabilityV2,
    kind: RepairActionKindV2,
    rank: u8,
) -> Result<ActionableDiagnosticV2, ConstructionErrorV2> {
    let repair = RepairActionV2::new(
        rank,
        kind,
        token("runner.arguments")?,
        Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
        Some(DiagnosticValueV2::Inline(TypedValueV2::U8(2))),
        token("runner.owner")?,
        Some("supply one canonical argument".to_owned()),
    )?;
    ActionableDiagnosticV2::new(
        code,
        retryability,
        Some(DiagnosticValueV2::Inline(TypedValueV2::U64(4))),
        Some(DiagnosticValueV2::Inline(TypedValueV2::U64(5))),
        token("runner.owner")?,
        vec![token("runner.arguments")?],
        no_claim_scope,
        vec![repair],
        DiagnosticEnvelopeGrantsV2::base_maxima(),
    )
}

fn diagnostic(
    no_claim_scope: NoClaimScopeRootV1,
    rank: u8,
) -> Result<ActionableDiagnosticV2, ConstructionErrorV2> {
    let repair = RepairActionV2::new(
        rank,
        RepairActionKindV2::ChangeArguments,
        token("runner.arguments")?,
        Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
        Some(DiagnosticValueV2::Inline(TypedValueV2::U8(2))),
        token("runner.owner")?,
        Some("supply one canonical argument".to_owned()),
    )?;
    ActionableDiagnosticV2::new(
        DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
        RetryabilityV2::AfterInputChange,
        None,
        None,
        token("runner.owner")?,
        vec![token("runner.arguments")?],
        no_claim_scope,
        vec![repair],
        DiagnosticEnvelopeGrantsV2::base_maxima(),
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BaseE2ePublicationStorageTotalsV1 {
    artifact: u64,
    system_publication: u64,
    publication: u64,
}

fn publication_storage() -> Result<BaseE2ePublicationStorageTotalsV1, ConstructionErrorV2> {
    let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
    let artifacts = [ArtifactStorageProjectionV2 {
        protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
        encoded_bytes: 1,
        stored_bytes: 1,
        envelope_non_payload_bytes: 0,
    }];
    let system_objects =
        SystemPublicationObjectRoleV2::ALL.map(|role| SystemObjectStorageProjectionV2 {
            role,
            protocol: PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            encoded_bytes: 1,
            stored_bytes: 1,
            envelope_non_payload_bytes: 0,
        });
    let artifact_encoded_bytes = artifacts
        .iter()
        .try_fold(0_u64, |total, artifact| {
            total.checked_add(artifact.encoded_bytes)
        })
        .expect("the frozen one-byte artifact encoded total cannot overflow");
    let artifact_stored_bytes = artifacts
        .iter()
        .try_fold(0_u64, |total, artifact| {
            total.checked_add(artifact.stored_bytes)
        })
        .expect("the frozen one-byte artifact stored total cannot overflow");
    let system_publication_stored_bytes = system_objects
        .iter()
        .try_fold(0_u64, |total, object| {
            total.checked_add(object.stored_bytes)
        })
        .expect("the frozen system-object stored total cannot overflow");
    let publication_stored_bytes = artifact_stored_bytes
        .checked_add(system_publication_stored_bytes)
        .expect("the frozen publication stored total cannot overflow");
    limits
        .validate_publication_storage(PublicationStorageProjectionV2 {
            artifacts: &artifacts,
            system_objects: &system_objects,
            artifact_encoded_bytes,
            artifact_stored_bytes,
            system_publication_stored_bytes,
            publication_stored_bytes,
        })
        .map_err(|_error| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "publication.storage",
                "the independently reconstructed base publication-storage projection",
                ConstructionObservedDataClassV2::BulkPayload,
            )
        })?;
    Ok(BaseE2ePublicationStorageTotalsV1 {
        artifact: artifact_stored_bytes,
        system_publication: system_publication_stored_bytes,
        publication: publication_stored_bytes,
    })
}

fn validate_source_input_set(
    inputs: &[BaseSourceClosureInputV1],
) -> Result<(), ConstructionErrorV2> {
    let mut paths = std::collections::BTreeSet::new();
    for input in inputs {
        if !paths.insert(input.path()) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_source_closure.path",
                "one unique path per exact source entry",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    if inputs.len() < EMBEDDED_SOURCE_FILES_V1.len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_source_closure.entries",
            "all exact embedded source entries",
            inputs.len(),
        ));
    }
    if inputs.len() > EMBEDDED_SOURCE_FILES_V1.len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_source_closure.entries",
            "no entries beyond the exact embedded source set",
            inputs.len(),
        ));
    }
    Ok(())
}

fn validate_source_input(
    input: &BaseSourceClosureInputV1,
    expected: &EmbeddedSourceFileV1,
    ordinal: usize,
    expected_paths: &std::collections::BTreeSet<&str>,
    snapshot_root: ContentHash,
) -> Result<(), ConstructionErrorV2> {
    validate_source_input_position(input, expected, ordinal, expected_paths)?;
    validate_source_input_metadata(input, expected, ordinal)?;
    validate_source_input_content(input, expected, ordinal, snapshot_root)
}

fn validate_source_input_position(
    input: &BaseSourceClosureInputV1,
    expected: &EmbeddedSourceFileV1,
    ordinal: usize,
    expected_paths: &std::collections::BTreeSet<&str>,
) -> Result<(), ConstructionErrorV2> {
    if input.path() == expected.path {
        return Ok(());
    }
    let (kind, expectation) = if expected_paths.contains(input.path()) {
        (
            ConstructionErrorKindV2::OutOfOrder,
            "the exact bytewise-lexicographic source order",
        )
    } else {
        (
            ConstructionErrorKindV2::Unexpected,
            "a member of the exact embedded source set",
        )
    };
    Err(ConstructionErrorV2::new(
        kind,
        "base_source_closure.path",
        expectation,
        ConstructionObservedV2::indexed_redacted(
            ordinal,
            ConstructionObservedDataClassV2::CallerControlledText,
        ),
    ))
}

fn validate_source_input_metadata(
    input: &BaseSourceClosureInputV1,
    expected: &EmbeddedSourceFileV1,
    ordinal: usize,
) -> Result<(), ConstructionErrorV2> {
    if input.owner_code() != expected.owner as u16 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.owner",
            "the exact sole declaration owner",
            ConstructionObservedV2::unsigned_pair(ordinal as u64, u64::from(input.owner_code())),
        ));
    }
    if input.source_route_code() != expected.source_route as u16 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.source_route",
            "the exact closed compile-time source route",
            ConstructionObservedV2::unsigned_pair(
                ordinal as u64,
                u64::from(input.source_route_code()),
            ),
        ));
    }
    if input.expected_source_identity_root() != expected_source_identity_root(expected) {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.expected_source_identity_root",
            "the exact declarative expected-source-identity root",
            ConstructionObservedV2::indexed_redacted(
                ordinal,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ));
    }
    if input.snapshot_policy_code() != expected.snapshot_policy as u16 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.snapshot_policy",
            "the exact common compiled-snapshot policy",
            ConstructionObservedV2::unsigned_pair(
                ordinal as u64,
                u64::from(input.snapshot_policy_code()),
            ),
        ));
    }
    Ok(())
}

fn validate_source_input_content(
    input: &BaseSourceClosureInputV1,
    expected: &EmbeddedSourceFileV1,
    ordinal: usize,
    snapshot_root: ContentHash,
) -> Result<(), ConstructionErrorV2> {
    let encoded_bytes = u64::try_from(input.bytes().len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_source_closure.encoded_bytes",
            "a u64 source byte length",
            input.bytes().len(),
        )
    })?;
    if input.encoded_bytes() != encoded_bytes {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.encoded_bytes",
            "the exact presented source byte length",
            ConstructionObservedV2::unsigned_pair(ordinal as u64, input.encoded_bytes()),
        ));
    }
    let content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, input.bytes());
    if input.content_root() != content_root {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.content_root",
            "the root of the exact presented source bytes",
            ConstructionObservedV2::indexed_redacted(
                ordinal,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ));
    }
    if input.bytes() != expected.bytes {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.bytes",
            "the exact compile-time included source bytes",
            ConstructionObservedDataClassV2::BulkPayload,
        ));
    }
    if input.snapshot_root() != snapshot_root {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_source_closure.snapshot_root",
            "the one exact common compiled-snapshot root",
            ConstructionObservedV2::indexed_redacted(
                ordinal,
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ));
    }
    Ok(())
}

fn source_closure_entry(
    declaration: &EmbeddedSourceFileV1,
    snapshot_root: ContentHash,
) -> Result<BaseSourceClosureEntryV1, ConstructionErrorV2> {
    let encoded_bytes = u64::try_from(declaration.bytes.len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_source_closure.encoded_bytes",
            "a u64 source byte length",
            declaration.bytes.len(),
        )
    })?;
    let content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, declaration.bytes);
    let source_identity_root = expected_source_identity_root(declaration);
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCEENTRY\x01", 2048)?;
    frame.push_str("source.path", declaration.path)?;
    frame.push_u16("source.owner", declaration.owner as u16)?;
    frame.push_u16("source.route", declaration.source_route as u16)?;
    frame.push_str(
        "source.expected_source_identity",
        declaration.expected_source_identity,
    )?;
    frame.push_bytes(
        "source.expected_source_identity_root",
        source_identity_root.as_bytes(),
    )?;
    frame.push_u16("source.snapshot_policy", declaration.snapshot_policy as u16)?;
    frame.push_u64("source.encoded_bytes", encoded_bytes)?;
    frame.push_bytes("source.content_root", content_root.as_bytes())?;
    frame.push_bytes("source.snapshot_root", snapshot_root.as_bytes())?;
    let entry_root = frame.root(BASE_SOURCE_FILE_ENTRY_DOMAIN_V1);
    Ok(BaseSourceClosureEntryV1 {
        path: declaration.path,
        owner: declaration.owner,
        source_route: declaration.source_route,
        expected_source_identity: declaration.expected_source_identity,
        expected_source_identity_root: source_identity_root,
        snapshot_policy: declaration.snapshot_policy,
        encoded_bytes,
        content_root,
        snapshot_root,
        entry_root,
    })
}

fn source_declaration(path: &str) -> Option<&'static EmbeddedSourceFileV1> {
    EMBEDDED_SOURCE_FILES_V1
        .binary_search_by(|declaration| declaration.path.as_bytes().cmp(path.as_bytes()))
        .ok()
        .map(|index| &EMBEDDED_SOURCE_FILES_V1[index])
}

fn expected_source_identity_root(declaration: &EmbeddedSourceFileV1) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"FSBASEEXPECTEDSOURCEIDENTITY\x01");
    for value in [declaration.path, declaration.expected_source_identity] {
        bytes.extend_from_slice(
            &u32::try_from(value.len())
                .expect("the static source-identity component length fits u32")
                .to_be_bytes(),
        );
        bytes.extend_from_slice(value.as_bytes());
    }
    bytes.extend_from_slice(&(declaration.owner as u16).to_be_bytes());
    bytes.extend_from_slice(&(declaration.source_route as u16).to_be_bytes());
    bytes.extend_from_slice(&(declaration.snapshot_policy as u16).to_be_bytes());
    hash_domain(BASE_EXPECTED_SOURCE_IDENTITY_DOMAIN_V1, &bytes)
}

fn compiled_source_snapshot() -> Result<CompatibleSourceSnapshotV1, ConstructionErrorV2> {
    static COMPILED_SOURCE_SNAPSHOT_V1: std::sync::OnceLock<
        Result<CompatibleSourceSnapshotV1, ConstructionErrorV2>,
    > = std::sync::OnceLock::new();
    COMPILED_SOURCE_SNAPSHOT_V1
        .get_or_init(|| {
            let frame = source_snapshot_frame(&EMBEDDED_SOURCE_FILES_V1)?;
            CompatibleSourceSnapshotV1::from_exact_frame(&frame)
        })
        .clone()
}

#[cfg(test)]
pub(crate) fn compatible_source_test_fixture_v1(
    discriminator: u8,
) -> Result<(CompatibleSourceSnapshotV1, CompatibleSourceMemberV1), ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCESNAPSHOT\x01", 1_024)?;
    frame.push_u8("snapshot.test_fixture_discriminator", discriminator)?;
    let snapshot = CompatibleSourceSnapshotV1::from_exact_frame(&frame)?;
    let path = LogicalBundlePathV1::new("crates/fs-evidence-runner/src/schema_impact.rs").map_err(
        |_| {
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "source_closure.test_fixture_member",
                "the exact schema-impact source path",
                ConstructionObservedDataClassV2::PhysicalLocator,
            )
        },
    )?;
    Ok((snapshot, CompatibleSourceMemberV1 { path, snapshot }))
}

fn compiled_source_snapshot_root() -> Result<ContentHash, ConstructionErrorV2> {
    Ok(compiled_source_snapshot()?.root().content_hash())
}

fn source_snapshot_frame(
    files: &[EmbeddedSourceFileV1],
) -> Result<CanonicalFrameV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCESNAPSHOT\x01", 64 * 1024)?;
    frame.push_bytes(
        "snapshot.current_direct_dependency_declaration_root",
        current_direct_dependency_declaration_root_v1().as_bytes(),
    )?;
    frame.push_u32(
        "snapshot.entry_count",
        u32::try_from(files.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "snapshot.entry_count",
                "a u32 compile-time input count",
                files.len(),
            )
        })?,
    )?;
    for file in files {
        frame.push_str("snapshot.path", file.path)?;
        frame.push_u16("snapshot.owner", file.owner as u16)?;
        frame.push_u16("snapshot.route", file.source_route as u16)?;
        frame.push_str(
            "snapshot.expected_source_identity",
            file.expected_source_identity,
        )?;
        frame.push_bytes(
            "snapshot.expected_source_identity_root",
            expected_source_identity_root(file).as_bytes(),
        )?;
        frame.push_u16("snapshot.snapshot_policy", file.snapshot_policy as u16)?;
        frame.push_u64(
            "snapshot.encoded_bytes",
            u64::try_from(file.bytes.len()).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "snapshot.encoded_bytes",
                    "a u64 compile-time input length",
                    file.bytes.len(),
                )
            })?,
        )?;
        frame.push_bytes(
            "snapshot.content_root",
            hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, file.bytes).as_bytes(),
        )?;
    }
    Ok(frame)
}

fn source_closure_root(
    entries: &[BaseSourceClosureEntryV1],
    snapshot_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCECLOSURE\x01", 16 * 1024)?;
    frame.push_bytes("source.snapshot_root", snapshot_root.as_bytes())?;
    frame.push_bytes(
        "source.current_direct_dependency_declaration_root",
        current_direct_dependency_declaration_root_v1().as_bytes(),
    )?;
    frame.push_u32(
        "source.entry_count",
        u32::try_from(entries.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "source.entry_count",
                "a u32 source entry count",
                entries.len(),
            )
        })?,
    )?;
    for entry in entries {
        frame.push_str("source.path", entry.path())?;
        frame.push_u16("source.owner", entry.owner() as u16)?;
        frame.push_u16("source.route", entry.source_route() as u16)?;
        frame.push_str(
            "source.expected_source_identity",
            entry.expected_source_identity(),
        )?;
        frame.push_bytes(
            "source.expected_source_identity_root",
            entry.expected_source_identity_root().as_bytes(),
        )?;
        frame.push_u16("source.snapshot_policy", entry.snapshot_policy() as u16)?;
        frame.push_u64("source.encoded_bytes", entry.encoded_bytes())?;
        frame.push_bytes("source.content_root", entry.content_root().as_bytes())?;
        frame.push_bytes(
            "source.entry_snapshot_root",
            entry.snapshot_root().as_bytes(),
        )?;
        frame.push_bytes("source.entry_root", entry.entry_root().as_bytes())?;
    }
    Ok(frame.root(BASE_SOURCE_CLOSURE_DOMAIN_V1))
}

const RUNTIME_LOG_COVERAGE_ID_V1: &str = "runtime-logging:aggregate-closed-log";
const LOCALLY_EXECUTED_COVERAGE_CLASSES_V1: [BaseCoverageManifestClassV1; 3] = [
    BaseCoverageManifestClassV1::ProjectionE2e,
    BaseCoverageManifestClassV1::RuntimeLogging,
    BaseCoverageManifestClassV1::SourceClosure,
];
const SOURCE_CLOSURE_COVERAGE_IDS_V1: [&str; 15] = [
    "source-closure:exact-positive",
    "source-closure:missing-entry-refusal",
    "source-closure:extra-entry-refusal",
    "source-closure:duplicate-entry-refusal",
    "source-closure:reordered-entry-refusal",
    "source-closure:stale-bytes-refusal",
    "source-closure:owner-refusal",
    "source-closure:source-route-refusal",
    "source-closure:expected-source-identity-refusal",
    "source-closure:snapshot-policy-refusal",
    "source-closure:encoded-length-refusal",
    "source-closure:content-root-refusal",
    "source-closure:mixed-snapshot-refusal",
    "source-closure:wrong-common-snapshot-refusal",
    "source-closure:dependency-declaration-refusal",
];

fn reconstruct_exact_local_coverage_report(
    manifest: &BaseCoverageManifestV1,
    observations: &std::collections::BTreeMap<
        String,
        (BaseCoveragePresentedOutcomeV1, ContentHash),
    >,
) -> Result<BaseCoverageCheckedReportV1, ConstructionErrorV2> {
    let selected_ids = manifest
        .cases()
        .iter()
        .filter(|case| LOCALLY_EXECUTED_COVERAGE_CLASSES_V1.contains(&case.class()))
        .map(crate::coverage::BaseCoverageManifestCaseV1::id)
        .collect::<Vec<_>>();
    let selected_id_set = selected_ids
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();

    if observations
        .keys()
        .any(|id| !selected_id_set.contains(id.as_str()))
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Unexpected,
            "coverage.runtime_result.source_case_id",
            "only independently declared locally executed coverage IDs",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    if selected_ids
        .iter()
        .any(|id| !observations.contains_key(*id))
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Missing,
            "coverage.runtime_result.source_case_id",
            "one observation for every independently declared locally executed coverage ID",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let selection = manifest.select_exact(&selected_ids)?;
    let results = selected_ids
        .iter()
        .map(|id| {
            let (outcome, evidence_root) = observations.get(*id).ok_or_else(|| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Missing,
                    "coverage.runtime_result.source_case_id",
                    "one observation for every independently declared locally executed coverage ID",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })?;
            BaseCoveragePresentedResultV1::new(manifest.root(), *id, *outcome, *evidence_root)
        })
        .collect::<Result<Vec<_>, _>>()?;
    BaseCoverageCheckedReportV1::reconstruct(manifest, &selection, &results)
}

#[allow(
    clippy::too_many_lines,
    reason = "the exact coverage manifest materializes every independently owned local and external class in one deterministic order"
)]
fn exact_coverage_manifest(
    journeys: &[BaseE2eJourneyProjectionV1],
) -> Result<BaseCoverageManifestV1, ConstructionErrorV2> {
    let mut extensions = Vec::new();
    for journey in [
        BaseE2eJourneyV1::CanonicalRunnerV2,
        BaseE2eJourneyV1::PublicationState,
        BaseE2eJourneyV1::PublicationV2,
        BaseE2eJourneyV1::RjoqHandoffV1,
        BaseE2eJourneyV1::VerifierV2,
    ] {
        for row_id in journey_row_id_oracle(journey) {
            extensions.push(BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                format!("projection-e2e:{}:{row_id}", journey.key()),
                "crates/fs-evidence-runner/src/projection.rs",
            )?);
        }
    }
    extensions.push(BaseCoverageCaseDeclarationV1::new(
        BaseCoverageManifestClassV1::RuntimeLogging,
        RUNTIME_LOG_COVERAGE_ID_V1,
        "crates/fs-evidence-runner/src/logging.rs",
    )?);
    for id in SOURCE_CLOSURE_COVERAGE_IDS_V1 {
        extensions.push(BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::SourceClosure,
            id,
            "crates/fs-evidence-runner/src/projection.rs",
        )?);
    }
    for (id, path) in [
        (
            "external-e2e:canonical-runner-v2",
            "scripts/ci/canonical_evidence_runner_v2.sh",
        ),
        (
            "external-e2e:publication-state-v2",
            "scripts/ci/e2e_evidence_runner_publication_state_v2.sh",
        ),
        (
            "external-e2e:publication-v2",
            "scripts/ci/e2e_evidence_runner_publication_v2.sh",
        ),
        (
            "external-e2e:rjoq-handoff-v1",
            "scripts/ci/verify_runner_rjoq_handoff_v1.sh",
        ),
        (
            "external-e2e:verifier-v2",
            "scripts/ci/e2e_evidence_verifier_v2.sh",
        ),
    ] {
        extensions.push(BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::ExternalE2eScript,
            id,
            path,
        )?);
    }
    extensions.push(BaseCoverageCaseDeclarationV1::new(
        BaseCoverageManifestClassV1::ExternalMutation,
        "external-mutation:base-contract-exact-result-join",
        "crates/fs-evidence-runner/src/projection.rs",
    )?);
    extensions.push(BaseCoverageCaseDeclarationV1::new(
        BaseCoverageManifestClassV1::ExternalGovernance,
        "external-governance:live-source-dependency-closure",
        "crates/fs-evidence-runner/src/dependency.rs",
    )?);
    extensions.sort_by(|left, right| {
        (left.class(), left.id(), left.source_path()).cmp(&(
            right.class(),
            right.id(),
            right.source_path(),
        ))
    });
    let manifest = BaseCoverageManifestV1::with_exact_extensions(&extensions)?;
    if manifest.case_count(BaseCoverageManifestClassV1::ProjectionE2e) != 98
        || manifest.case_count(BaseCoverageManifestClassV1::RuntimeLogging) != 1
        || manifest.case_count(BaseCoverageManifestClassV1::SourceClosure) != 15
        || manifest.case_count(BaseCoverageManifestClassV1::ExternalE2eScript) != 5
        || manifest.case_count(BaseCoverageManifestClassV1::ExternalMutation) != 1
        || manifest.case_count(BaseCoverageManifestClassV1::ExternalGovernance) != 1
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "coverage.extension_counts",
            "98 projection, 1 logging, 15 source, 5 external E2E, 1 mutation, and 1 governance case",
            extensions.len(),
        ));
    }
    for journey in journeys {
        let literal_ids = journey_row_id_oracle(journey.journey());
        if literal_ids.len() != journey.rows().len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "coverage.journey_row_count",
                "the exact independent literal row-ID count",
                ConstructionObservedV2::closed(&journey.journey()),
            ));
        }
        for (literal_id, row) in literal_ids.iter().zip(journey.rows()) {
            let coverage_id = format!("projection-e2e:{}:{literal_id}", journey.journey().key());
            if row.id().as_str() != *literal_id || manifest.case(&coverage_id).is_none() {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Missing,
                    "coverage.projection_e2e",
                    "one independently declared coverage case per journey row",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
    }
    Ok(manifest)
}

fn coverage_inventory(
    journeys: &[BaseE2eJourneyProjectionV1],
) -> Result<BaseCoverageInventoryV1, ConstructionErrorV2> {
    let mut cases = Vec::new();
    for class in BaseCoverageClassV1::ALL {
        let mut ordinal = 0_u32;
        if class == BaseCoverageClassV1::ProjectionE2e {
            for journey in journeys {
                for row in journey.rows() {
                    ordinal = ordinal.checked_add(1).ok_or_else(sequence_overflow)?;
                    cases.push(coverage_source_case(
                        class,
                        ordinal,
                        &format!("{}.{}", journey.journey().key(), row.id().as_str()),
                        "crates/fs-evidence-runner/src/projection.rs",
                    )?);
                }
            }
        } else {
            for template in COVERAGE_CASE_TEMPLATES_V1
                .iter()
                .filter(|template| template.class == class)
            {
                ordinal = ordinal.checked_add(1).ok_or_else(sequence_overflow)?;
                cases.push(coverage_source_case(
                    class,
                    ordinal,
                    template.id,
                    template.source_path,
                )?);
            }
        }
        if ordinal == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "coverage.source_cases",
                "at least one immutable source case per coverage class",
                ConstructionObservedV2::closed(&class),
            ));
        }
    }

    let expected_e2e_rows = BaseE2eJourneyV1::ALL
        .iter()
        .try_fold(0_usize, |total, journey| {
            total
                .checked_add(journey_case_kinds(*journey).len())
                .ok_or_else(sequence_overflow)
        })?;
    let actual_e2e_rows = cases
        .iter()
        .filter(|source_case| source_case.class == BaseCoverageClassV1::ProjectionE2e)
        .count();
    if actual_e2e_rows != expected_e2e_rows || actual_e2e_rows != 98 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "coverage.projection_e2e_source_cases",
            "the exact five journey-specific manifests, or 98 source cases",
            actual_e2e_rows,
        ));
    }

    let mut seen = std::collections::BTreeSet::new();
    for source_case in &cases {
        if !seen.insert((source_case.class.code(), source_case.id.as_str())) {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "coverage.source_case_id",
                "unique IDs within every coverage class",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    let root = coverage_inventory_root(&cases)?;
    Ok(BaseCoverageInventoryV1 {
        cases: cases.into_boxed_slice(),
        root,
    })
}

fn coverage_source_case(
    class: BaseCoverageClassV1,
    ordinal: u32,
    id: &str,
    source_path: &str,
) -> Result<BaseCoverageSourceCaseV1, ConstructionErrorV2> {
    if !EMBEDDED_SOURCE_FILES_V1
        .iter()
        .any(|source| source.path == source_path)
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Unexpected,
            "coverage.source_path",
            "a member of the exact embedded source closure",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let id = token(id)?;
    let source_path = LogicalBundlePathV1::new(source_path).map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "coverage.source_path",
            "a validated relative source path",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })?;
    Ok(BaseCoverageSourceCaseV1 {
        class,
        ordinal,
        id,
        source_path,
    })
}

fn coverage_inventory_root(
    cases: &[BaseCoverageSourceCaseV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASECOVERAGE\x01", 64 * 1024)?;
    frame.push_u32(
        "coverage.source_case_count",
        u32::try_from(cases.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "coverage.source_case_count",
                "a u32 source-case count",
                cases.len(),
            )
        })?,
    )?;
    for source_case in cases {
        frame.push_u16("coverage.class", source_case.class.code())?;
        frame.push_u32("coverage.ordinal", source_case.ordinal)?;
        frame.push_str("coverage.id", source_case.id.as_str())?;
        frame.push_str("coverage.source_path", source_case.source_path.as_str())?;
    }
    Ok(frame.root(BASE_COVERAGE_INVENTORY_DOMAIN_V1))
}

fn journey_root(
    journey: BaseE2eJourneyV1,
    downstream_owner: &StableTokenV2,
    driver_name: &StableTokenV2,
    script_path: &LogicalBundlePathV1,
    manifest_path: &LogicalBundlePathV1,
    rows: &[BaseE2eProjectionRowV1],
    source_closure_root: ContentHash,
    log_schema_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEJOURNEY\x01", 64 * 1024)?;
    frame.push_u16("projection.api_generation", 2)?;
    frame.push_u16("projection.wire_version", 1)?;
    frame.push_u16("projection.journey", journey.code())?;
    frame.push_str("projection.downstream_owner", downstream_owner.as_str())?;
    frame.push_str("projection.driver_name", driver_name.as_str())?;
    frame.push_str("projection.script_path", script_path.as_str())?;
    frame.push_str("projection.manifest_path", manifest_path.as_str())?;
    frame.push_bytes(
        "projection.source_closure_root",
        source_closure_root.as_bytes(),
    )?;
    frame.push_bytes("projection.log_schema_root", log_schema_root.as_bytes())?;
    frame.push_u32(
        "projection.row_count",
        u32::try_from(rows.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "projection.row_count",
                "a u32 row count",
                rows.len(),
            )
        })?,
    )?;
    for row in rows {
        if row.journey() != journey
            || row.downstream_owner() != downstream_owner
            || row.downstream_driver() != driver_name
            || row.downstream_script() != script_path
            || row.downstream_manifest() != manifest_path
            || row.source_closure_root() != source_closure_root
            || row.log_schema_root() != log_schema_root
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "projection.journey_row_context",
                "the exact journey, owner, driver, script, manifest, and source closure",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        frame.push_str("projection.row_id", row.id.as_str())?;
        frame.push_u16("projection.case_kind", row.kind.code())?;
        frame.push_u16("projection.expected", row.expected.code())?;
        frame.push_u32("projection.semantic_cell_count", row.semantic_cell_count)?;
        frame.push_u32("projection.positive_cell_count", row.positive_cell_count)?;
        frame.push_u32(
            "projection.expected_refusal_cell_count",
            row.expected_refusal_cell_count,
        )?;
        frame.push_u32(
            "projection.unsupported_cell_count",
            row.unsupported_cell_count,
        )?;
        frame.push_bytes(
            "projection.semantic_manifest_root",
            row.semantic_manifest_root.as_bytes(),
        )?;
        frame.push_bytes(
            "projection.journey_mapping_root",
            row.mapping_root.as_bytes(),
        )?;
        frame.push_str(
            "projection.consumption_rationale",
            row.consumption_rationale.as_str(),
        )?;
        frame.push_str(
            "projection.fixture_reference",
            row.fixture_reference.as_str(),
        )?;
        frame.push_u16("projection.unit_tag", row.unit.tag())?;
        if let Some(registered_id) = row.unit.registered_id() {
            frame.push_u16("projection.unit_registered_id", registered_id)?;
        }
        frame.push_str("projection.no_claim_scope", row.no_claim_scope.as_str())?;
    }
    Ok(frame.root(BASE_E2E_JOURNEY_PROJECTION_DOMAIN_V1))
}

fn projection_root(
    journeys: &[BaseE2eJourneyProjectionV1],
    source_closure_root: ContentHash,
    coverage_inventory_root: ContentHash,
    coverage_manifest_root: ContentHash,
    close_manifest_root: ContentHash,
    log_schema_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEPROJECTION\x01", 4096)?;
    frame.push_u16("projection.api_generation", 2)?;
    frame.push_u16("projection.wire_version", 1)?;
    frame.push_u16(
        "projection.journey_count",
        u16::try_from(journeys.len()).expect("exactly five journeys"),
    )?;
    for journey in journeys {
        frame.push_u16("projection.journey", journey.journey.code())?;
        frame.push_bytes("projection.journey_root", journey.root.as_bytes())?;
    }
    frame.push_bytes(
        "projection.source_closure_root",
        source_closure_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.coverage_inventory_root",
        coverage_inventory_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.coverage_manifest_root",
        coverage_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "projection.close_manifest_root",
        close_manifest_root.as_bytes(),
    )?;
    frame.push_bytes("projection.log_schema_root", log_schema_root.as_bytes())?;
    Ok(frame.root(BASE_E2E_PROJECTION_DOMAIN_V1))
}

fn case_terminal_log_event(
    sequence: u32,
    journey: &BaseE2eJourneyProjectionV1,
    row: &BaseE2eProjectionRowV1,
    result: &BaseE2eRowResultV1,
    outcome: BaseE2eOutcomeV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
    let mut case_fields = vec![
        field("checked-cells", TypedValueV2::U32(result.checked_cells()))?,
        field("expected", TypedValueV2::Token(token(row.expected.name())?))?,
        field(
            "observed",
            TypedValueV2::Token(token(result.observed().name())?),
        )?,
        field(
            "semantic-cell-count",
            TypedValueV2::U32(row.semantic_cell_count()),
        )?,
        field(
            "semantic-manifest-root",
            opaque_root(row.semantic_manifest_root())?,
        )?,
        field("row-result-root", opaque_root(result.root())?)?,
        field(
            "expected-detail-manifest-root",
            opaque_root(result.expected_detail_manifest_root())?,
        )?,
        field(
            "observed-detail-manifest-root",
            opaque_root(result.observed_detail_manifest_root())?,
        )?,
        field(
            "expected-detail-cells",
            TypedValueV2::U32(result.expected_detail_cell_count()),
        )?,
        field(
            "observed-detail-cells",
            TypedValueV2::U32(result.observed_detail_cell_count()),
        )?,
        field(
            "detail-cells-matched",
            TypedValueV2::U32(result.detail_cells_matched()),
        )?,
        field(
            "logical-unit",
            TypedValueV2::Token(token(row.unit().name())?),
        )?,
        field(
            "no-claim-scope",
            TypedValueV2::Digest(harness.no_claim_scope().digest().clone()),
        )?,
    ];
    case_fields.extend(partition_fields(
        result.positive_eligible(),
        result.positive_matched(),
        result.expected_refusals(),
        result.expected_refusals_matched(),
        result.unexpected_mismatches(),
    )?);
    case_fields.push(field(
        "unsupported",
        TypedValueV2::U32(result.unsupported()),
    )?);
    if let Some(expected_detail) = row.expected_detail() {
        case_fields.push(field(
            "expected-detail",
            TypedValueV2::Token(expected_detail.clone()),
        )?);
    }
    if let Some(first_failed_cell) = result.first_unexpected_cell() {
        case_fields.push(field(
            "first-failed-cell",
            TypedValueV2::Token(token(first_failed_cell)?),
        )?);
        case_fields.push(field(
            "first-detail-divergence-root",
            opaque_root(result.first_divergence_root().ok_or_else(|| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Missing,
                    "base_e2e_projection.first_divergence_root",
                    "one typed bounded detail-or-row-contract divergence root for every failed terminal",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })?)?,
        )?);
    }
    case_fields.extend(case_detail_fields(row.kind())?);
    log_event(
        sequence,
        journey,
        Some(row),
        BaseE2eLogKindV1::CaseTerminal,
        outcome,
        harness,
        None,
        case_fields,
    )
}

#[allow(
    clippy::too_many_arguments,
    reason = "the closed log-event constructor binds every event coordinate explicitly before adding the caller-specific field set"
)]
fn log_event(
    sequence: u32,
    journey: &BaseE2eJourneyProjectionV1,
    row: Option<&BaseE2eProjectionRowV1>,
    kind: BaseE2eLogKindV1,
    outcome: BaseE2eOutcomeV1,
    harness: &BaseE2eHarnessIdentityV1,
    execution_root: Option<ContentHash>,
    mut fields: Vec<BaseE2eLogFieldV1>,
) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
    fields.extend(environment_log_fields(harness)?);
    fields.extend([
        field("projection-root", opaque_root(journey.manifest_root())?)?,
        field("manifest-root", opaque_root(journey.manifest_root())?)?,
        field(
            "downstream-script-mapping",
            TypedValueV2::RelativePath(journey.script_path().clone()),
        )?,
    ]);
    if let Some(execution_root) = execution_root {
        fields.push(field("execution-root", opaque_root(execution_root)?)?);
    }
    BaseE2eLogEventV1::new(
        sequence,
        token(journey.journey.key())?,
        row.map(|row| row.id.clone()),
        kind,
        outcome,
        fields,
        None,
        vec![
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(token(journey.journey.key())?),
        ],
    )
}

fn projection_summary_log_event(
    sequence: u32,
    projection: &RunnerV2BaseE2eProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
    execution_root: ContentHash,
    mut fields: Vec<BaseE2eLogFieldV1>,
) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
    fields.extend(environment_log_fields(harness)?);
    fields.extend([
        field("projection-root", opaque_root(projection.manifest_root())?)?,
        field("manifest-root", opaque_root(projection.manifest_root())?)?,
        field("execution-root", opaque_root(execution_root)?)?,
    ]);
    BaseE2eLogEventV1::new(
        sequence,
        token("all")?,
        None,
        BaseE2eLogKindV1::ProjectionSummary,
        BaseE2eOutcomeV1::NotApplicable,
        fields,
        None,
        vec![
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(token("all")?),
        ],
    )
}

fn environment_log_fields(
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<Vec<BaseE2eLogFieldV1>, ConstructionErrorV2> {
    Ok(vec![
        field("api-generation", TypedValueV2::U16(2))?,
        field("wire-version", TypedValueV2::U16(1))?,
        field(
            "source-root",
            TypedValueV2::Digest(harness.source.digest().clone()),
        )?,
        field(
            "build-root",
            TypedValueV2::Digest(harness.build.digest().clone()),
        )?,
        field(
            "toolchain-root",
            TypedValueV2::Digest(harness.toolchain.digest().clone()),
        )?,
        field("target", TypedValueV2::Token(harness.target.clone()))?,
        field(
            "feature-count",
            TypedValueV2::U32(
                u32::try_from(harness.features.len()).expect("bounded feature fixture"),
            ),
        )?,
        field("feature-set-root", opaque_root(harness.feature_set_root())?)?,
        field("target-root", opaque_root(harness.target_root())?)?,
    ])
}

fn count_fields(
    eligible: u32,
    passed: u32,
    failed: u32,
    unsupported: u32,
) -> Result<Vec<BaseE2eLogFieldV1>, ConstructionErrorV2> {
    Ok(vec![
        field("eligible", TypedValueV2::U32(eligible))?,
        field("passed", TypedValueV2::U32(passed))?,
        field("failed", TypedValueV2::U32(failed))?,
        field("unsupported", TypedValueV2::U32(unsupported))?,
    ])
}

fn partition_fields(
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unexpected_mismatches: u32,
) -> Result<Vec<BaseE2eLogFieldV1>, ConstructionErrorV2> {
    Ok(vec![
        field("positive-eligible", TypedValueV2::U32(positive_eligible))?,
        field("positive-matched", TypedValueV2::U32(positive_matched))?,
        field("expected-refusals", TypedValueV2::U32(expected_refusals))?,
        field(
            "expected-refusals-matched",
            TypedValueV2::U32(expected_refusals_matched),
        )?,
        field(
            "unexpected-mismatches",
            TypedValueV2::U32(unexpected_mismatches),
        )?,
    ])
}

fn field(
    name: &'static str,
    value: TypedValueV2,
) -> Result<BaseE2eLogFieldV1, ConstructionErrorV2> {
    Ok(BaseE2eLogFieldV1::new(token(name)?, value))
}

fn token(value: &str) -> Result<StableTokenV2, ConstructionErrorV2> {
    StableTokenV2::new(value).map_err(|_error| {
        ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "projection.stable_token",
            "a valid stable token",
            ConstructionObservedDataClassV2::CallerControlledText,
        )
    })
}

fn sequence_overflow() -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::ArithmeticOverflow,
        "projection.logical_sequence",
        "a checked u32 sequence",
        ConstructionObservedV2::fixed(ConstructionFixedObservationV2::Overflow),
    )
}

#[cfg(test)]
mod tests {
    use super::{
        BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, BaseE2eCaseKindV1, BaseE2eHarnessIdentityV1,
        BaseE2eJourneyV1, BaseE2eMatchedPartitionV1, BaseE2eObservedCountsV1,
        BaseE2ePresentedRowResultV1, BaseSourceClosureInputV1, BaseSourceOwnerV1,
        BaseSourceRouteV1, BaseSourceSnapshotPolicyV1, EMBEDDED_SOURCE_FILES_V1,
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_CONSUMED_ROW_COUNT_V2,
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_OWNED_ROW_COUNT_V2,
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2,
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SCHEMA_IDS_V2,
        RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SOURCE_PATH_COUNT_V2,
        RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1, RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1,
        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1, RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1,
        RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1, RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
        RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1,
        RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
        RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1,
        RUNNER_V2_PHASE_ONE_EXPECTED_PARTITIONS_DOMAIN_V1,
        RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1, RunnerV2BaseE2eProjectionV1,
        RunnerV2BaseSchemaImpactProjectionV2, RunnerV2BaseSourceClosureV1,
        RunnerV2PhaseOneContractContributionV2, compare_base_e2e_journey_results_v1,
        encode_construction_error, execute_case, expected_source_identity_root, journey_root,
        journey_row_id_oracle, journey_row_root, phase_one_build_contract_root_v1,
        phase_one_case_manifest_contract_root_from_components_v1, phase_one_close_ledger_cells_v1,
        phase_one_contribution_projection_root_v2,
        phase_one_expected_partitions_root_from_ledger_v1, phase_one_payload_root_v2,
        reconstruct_exact_local_coverage_report, run_base_e2e_journey_v1,
        run_base_e2e_projection_v1, schema_impact_projection_root_v2,
        schema_impact_projection_row_root_v2, source_closure_entry,
        validate_exact_phase_one_observer_contract_v1,
    };
    use crate::budget::RunnerBudgetViolationKindV2;
    use crate::catalog::{DigestRoleV2, RepairActionKindV2};
    use crate::construction::ConstructionErrorKindV2;
    use crate::coverage::{
        BaseCoverageCloseCapabilityProfileV1, BaseCoverageCloseContributionBudgetsV1,
        BaseCoverageCloseDecisionV1, BaseCoverageCloseDeferredEvidenceEnvelopeV1,
        BaseCoverageCloseDeferredObservationContractV1, BaseCoverageCloseDownstreamContributionV2,
        BaseCoverageCloseEvidenceKindV1, BaseCoverageCloseManifestV1, BaseCoverageClosePartitionV1,
        BaseCoverageCloseReasonCodeV1, BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
        BaseCoverageManifestClassV1, CanonicalSchemaImpactDispositionV1, DeferredReasonV1,
        RuntimeObservationDispositionV1, compatible_source_snapshot_root_from_exact_frame_v1,
        frozen_downstream_close_contribution_budgets_v1,
    };
    use crate::dependency::current_direct_dependency_declaration_root_v1;
    use crate::identity::{
        BuildIdentityRootV2, NoClaimScopeRootV1, SourceIdentityRootV2, ToolchainIdentityRootV2,
    };
    use crate::limits::RunnerLimitValueV2;
    use crate::path::LogicalBundlePathV1;
    use crate::schema_impact::{
        RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1, RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1,
        SchemaImpactManifestRelationV1, SchemaImpactManifestRowSourceV1, SchemaImpactNoClaimV1,
        runner_v2_base_schema_impact_field_name_mutant_manifest_v1,
        source_frozen_schema_impact_manifest_v1,
    };
    use crate::value::StableTokenV2;
    use fs_blake3::hash_domain;

    const EXPECTED_SOURCE_PATHS_V1: [&str; 27] = [
        ".cargo/config.toml",
        "Cargo.lock",
        "Cargo.toml",
        "constellation.lock",
        "crates/fs-evidence-runner/CONTRACT.md",
        "crates/fs-evidence-runner/Cargo.toml",
        "crates/fs-evidence-runner/src/budget.rs",
        "crates/fs-evidence-runner/src/canonical.rs",
        "crates/fs-evidence-runner/src/capability.rs",
        "crates/fs-evidence-runner/src/catalog.rs",
        "crates/fs-evidence-runner/src/command.rs",
        "crates/fs-evidence-runner/src/construction.rs",
        "crates/fs-evidence-runner/src/coverage.rs",
        "crates/fs-evidence-runner/src/dependency.rs",
        "crates/fs-evidence-runner/src/diagnostic.rs",
        "crates/fs-evidence-runner/src/extension.rs",
        "crates/fs-evidence-runner/src/identity.rs",
        "crates/fs-evidence-runner/src/lib.rs",
        "crates/fs-evidence-runner/src/limits.rs",
        "crates/fs-evidence-runner/src/logging.rs",
        "crates/fs-evidence-runner/src/path.rs",
        "crates/fs-evidence-runner/src/projection.rs",
        "crates/fs-evidence-runner/src/publication.rs",
        "crates/fs-evidence-runner/src/schema_impact.rs",
        "crates/fs-evidence-runner/src/state.rs",
        "crates/fs-evidence-runner/src/value.rs",
        "rust-toolchain.toml",
    ];

    fn presented<T>(
        role: DigestRoleV2,
        domain: &str,
        parser: impl FnOnce(DigestRoleV2, &str, &str) -> Result<T, crate::identity::IdentityError>,
    ) -> T {
        parser(role, domain, &"00".repeat(32)).expect("presented fixture")
    }

    fn presented_byte<T>(
        role: DigestRoleV2,
        domain: &str,
        byte: u8,
        parser: impl FnOnce(DigestRoleV2, &str, &str) -> Result<T, crate::identity::IdentityError>,
    ) -> T {
        parser(role, domain, &format!("{byte:02x}").repeat(32)).expect("presented fixture")
    }

    fn harness() -> BaseE2eHarnessIdentityV1 {
        BaseE2eHarnessIdentityV1::new(
            presented(
                DigestRoleV2::Source,
                SourceIdentityRootV2::DESCRIPTOR.domain(),
                SourceIdentityRootV2::parse_presented,
            ),
            presented(
                DigestRoleV2::Build,
                BuildIdentityRootV2::DESCRIPTOR.domain(),
                BuildIdentityRootV2::parse_presented,
            ),
            presented(
                DigestRoleV2::Toolchain,
                ToolchainIdentityRootV2::DESCRIPTOR.domain(),
                ToolchainIdentityRootV2::parse_presented,
            ),
            StableTokenV2::new("aarch64-apple-darwin").expect("target"),
            vec![StableTokenV2::new("default").expect("feature")],
            presented(
                DigestRoleV2::ClaimScope,
                NoClaimScopeRootV1::DESCRIPTOR.domain(),
                NoClaimScopeRootV1::parse_presented,
            ),
        )
        .expect("harness identity")
    }

    #[allow(clippy::too_many_arguments)]
    fn harness_variant(
        source_byte: u8,
        build_byte: u8,
        toolchain_byte: u8,
        target: &str,
        features: &[&str],
        no_claim_byte: u8,
    ) -> BaseE2eHarnessIdentityV1 {
        BaseE2eHarnessIdentityV1::new(
            presented_byte(
                DigestRoleV2::Source,
                SourceIdentityRootV2::DESCRIPTOR.domain(),
                source_byte,
                SourceIdentityRootV2::parse_presented,
            ),
            presented_byte(
                DigestRoleV2::Build,
                BuildIdentityRootV2::DESCRIPTOR.domain(),
                build_byte,
                BuildIdentityRootV2::parse_presented,
            ),
            presented_byte(
                DigestRoleV2::Toolchain,
                ToolchainIdentityRootV2::DESCRIPTOR.domain(),
                toolchain_byte,
                ToolchainIdentityRootV2::parse_presented,
            ),
            StableTokenV2::new(target).expect("target"),
            features
                .iter()
                .map(|feature| StableTokenV2::new(*feature).expect("feature"))
                .collect(),
            presented_byte(
                DigestRoleV2::ClaimScope,
                NoClaimScopeRootV1::DESCRIPTOR.domain(),
                no_claim_byte,
                NoClaimScopeRootV1::parse_presented,
            ),
        )
        .expect("harness variant")
    }

    fn frozen_source_inputs() -> Vec<BaseSourceClosureInputV1> {
        EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| BaseSourceClosureInputV1::presented(file.path, file.bytes.to_vec()))
            .collect()
    }

    #[derive(Clone)]
    struct PhaseOneContributionFixture {
        source_owner: String,
        payload_root: fs_blake3::ContentHash,
        expected_partitions_root: fs_blake3::ContentHash,
        per_cell_budgets: BaseCoverageCloseContributionBudgetsV1,
        shard_budgets: BaseCoverageCloseContributionBudgetsV1,
        schema_root: fs_blake3::ContentHash,
        log_schema_root: fs_blake3::ContentHash,
        source_root: SourceIdentityRootV2,
        build_root: BuildIdentityRootV2,
        toolchain_root: ToolchainIdentityRootV2,
        target_root: fs_blake3::ContentHash,
        feature_set_root: fs_blake3::ContentHash,
        retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
        observer_contract: BaseCoverageCloseDeferredObservationContractV1,
        no_claim: StableTokenV2,
    }

    impl PhaseOneContributionFixture {
        fn from_contribution(contribution: &BaseCoverageCloseDownstreamContributionV2) -> Self {
            Self {
                source_owner: contribution.source_owner().to_owned(),
                payload_root: contribution.payload_root(),
                expected_partitions_root: contribution.expected_partitions_root(),
                per_cell_budgets: contribution.per_cell_budgets(),
                shard_budgets: contribution.shard_budgets(),
                schema_root: contribution.schema_root(),
                log_schema_root: contribution.log_schema_root(),
                source_root: contribution.source_root().clone(),
                build_root: contribution.build_root().clone(),
                toolchain_root: contribution.toolchain_root().clone(),
                target_root: contribution.target_root(),
                feature_set_root: contribution.feature_set_root(),
                retained_artifact_policy: contribution.retained_artifact_policy(),
                observer_contract: contribution.observer_contract().clone(),
                no_claim: contribution.no_claim().clone(),
            }
        }

        fn build(
            &self,
        ) -> Result<BaseCoverageCloseDownstreamContributionV2, crate::ConstructionErrorV2> {
            BaseCoverageCloseDownstreamContributionV2::new(
                self.source_owner.clone(),
                self.payload_root,
                self.expected_partitions_root,
                self.per_cell_budgets,
                self.shard_budgets,
                self.schema_root,
                self.log_schema_root,
                self.source_root.clone(),
                self.build_root.clone(),
                self.toolchain_root.clone(),
                self.target_root,
                self.feature_set_root,
                self.retained_artifact_policy,
                self.observer_contract.clone(),
                self.no_claim.clone(),
            )
        }
    }

    fn vocabulary_only_phase_one_partition_root_v1() -> fs_blake3::ContentHash {
        let literals = [
            (1_u16, "positive"),
            (2_u16, "expected-refusal"),
            (3_u16, "expected-failure"),
            (4_u16, "mutation"),
            (5_u16, "unsupported"),
            (6_u16, "inapplicable"),
        ];
        let mut frame =
            crate::canonical::CanonicalFrameV1::new(b"FSPHASE1PARTITIONVOCAB\x01", 4 * 1024)
                .expect("independent vocabulary-only frame");
        frame
            .push_u32("vocabulary.partition_count", 6)
            .expect("partition count");
        for (code, name) in literals {
            frame
                .push_u16("vocabulary.partition_code", code)
                .expect("partition code");
            frame
                .push_str("vocabulary.partition_name", name)
                .expect("partition name");
        }
        frame
            .push_str(
                "vocabulary.no_claim",
                "partition-vocabulary-root-proves-closed-names-and-codes-not-ledger-membership-or-execution",
            )
            .expect("vocabulary no-claim");
        frame.root("org.frankensim.fs-evidence-runner.phase-one-partition-vocabulary.v1")
    }

    fn count_only_phase_one_partition_root_v1(
        close_manifest: &BaseCoverageCloseManifestV1,
        source_closure: &RunnerV2BaseSourceClosureV1,
    ) -> fs_blake3::ContentHash {
        let mut frame =
            crate::canonical::CanonicalFrameV1::new(b"FSPHASE1PARTITIONS\x01", 64 * 1024)
                .expect("count-only frame");
        frame
            .push_bytes(
                "partitions.source_manifest_root",
                close_manifest.source_manifest_root().as_bytes(),
            )
            .expect("source-manifest root");
        frame
            .push_bytes(
                "partitions.reason_registry_root",
                close_manifest.reason_registry_root().as_bytes(),
            )
            .expect("reason-registry root");
        frame
            .push_bytes(
                "partitions.close_manifest_root",
                close_manifest.root().as_bytes(),
            )
            .expect("close-manifest root");
        frame
            .push_bytes(
                "partitions.compatible_source_snapshot_root",
                source_closure.snapshot_root().as_bytes(),
            )
            .expect("compatible snapshot root");
        frame
            .push_bytes(
                "partitions.source_closure_root",
                source_closure.root().as_bytes(),
            )
            .expect("source-closure root");
        frame
            .push_bytes(
                "partitions.vocabulary_root",
                vocabulary_only_phase_one_partition_root_v1().as_bytes(),
            )
            .expect("independent vocabulary root");
        frame
            .push_u32("partitions.count", 6)
            .expect("partition count");
        for partition in BaseCoverageClosePartitionV1::ALL {
            let count = u32::try_from(
                close_manifest
                    .cells()
                    .iter()
                    .filter(|cell| cell.partition() == partition)
                    .count(),
            )
            .expect("bounded partition cardinality");
            frame
                .push_u16("partitions.code", partition.code())
                .expect("partition code");
            frame
                .push_str("partitions.name", partition.stable_name())
                .expect("partition name");
            frame
                .push_u32("partitions.expected_cell_count", count)
                .expect("partition count");
            frame
                .push_u16("partitions.eventual_owner_requires_exact_match", 1)
                .expect("partition exact-match flag");
        }
        frame
            .push_u32(
                "partitions.cell_count",
                u32::try_from(close_manifest.cells().len()).expect("bounded close cell count"),
            )
            .expect("close cell count");
        frame
            .push_u32("partitions.required_unexpected_mismatches", 0)
            .expect("zero mismatches");
        frame
            .push_u32("partitions.required_execution_failures", 0)
            .expect("zero failures");
        frame
            .push_u32("partitions.required_unexplained_skips", 0)
            .expect("zero skips");
        frame
            .push_str(
                "partitions.no_claim",
                "expected-partitions-do-not-prove-owner-execution-or-matched-results",
            )
            .expect("partition no-claim");
        frame.root(RUNNER_V2_PHASE_ONE_EXPECTED_PARTITIONS_DOMAIN_V1)
    }

    fn independent_phase_one_expected_result_root_v1(
        ledger_ordinal: u32,
        cell_id: &str,
        cell: &crate::coverage::BaseCoverageCloseManifestCellV1,
        close_manifest: &BaseCoverageCloseManifestV1,
        source_closure: &RunnerV2BaseSourceClosureV1,
    ) -> fs_blake3::ContentHash {
        let mut frame =
            crate::canonical::CanonicalFrameV1::new(b"FSPHASE1EXPECTEDRESULT\x01", 2 * 1024)
                .expect("independent expected-result frame");
        frame
            .push_u32("expected_result.ledger_ordinal", ledger_ordinal)
            .expect("ledger ordinal");
        frame
            .push_str("expected_result.cell_id", cell_id)
            .expect("cell ID");
        frame
            .push_str("expected_result.source_case_id", cell.source_case_id())
            .expect("source case ID");
        frame
            .push_u16("expected_result.partition", cell.partition().code())
            .expect("partition");
        frame
            .push_u16(
                "expected_result.expected_decision",
                cell.expected_decision().code(),
            )
            .expect("expected decision");
        frame
            .push_bytes("expected_result.close_cell_root", cell.root().as_bytes())
            .expect("close-cell root");
        frame
            .push_bytes(
                "expected_result.compatible_source_snapshot_root",
                source_closure.snapshot_root().as_bytes(),
            )
            .expect("snapshot root");
        frame
            .push_bytes(
                "expected_result.case_manifest_root",
                close_manifest.root().as_bytes(),
            )
            .expect("case-manifest root");
        frame
            .push_bytes(
                "expected_result.source_closure_root",
                source_closure.root().as_bytes(),
            )
            .expect("source-closure root");
        frame
            .push_str(
                "expected_result.no_claim",
                "expected-result-root-is-a-source-declared-expected-decision-contract-not-an-observed-output-downstream-actual-roots-must-be-derived-from-real-outcomes-not-copied",
            )
            .expect("expected-result no-claim");
        frame.root("org.frankensim.fs-evidence-runner.phase-one-expected-result.v1")
    }

    #[allow(
        clippy::too_many_lines,
        reason = "the independent handwritten oracle intentionally repeats the complete AC64 frame instead of calling the production encoder"
    )]
    fn independent_phase_one_full_ledger_root_v1(
        close_manifest: &BaseCoverageCloseManifestV1,
        source_closure: &RunnerV2BaseSourceClosureV1,
    ) -> fs_blake3::ContentHash {
        let partition_literals = [
            (BaseCoverageClosePartitionV1::Positive, 1_u16, "positive"),
            (
                BaseCoverageClosePartitionV1::ExpectedRefusal,
                2_u16,
                "expected-refusal",
            ),
            (
                BaseCoverageClosePartitionV1::ExpectedFailure,
                3_u16,
                "expected-failure",
            ),
            (BaseCoverageClosePartitionV1::Mutation, 4_u16, "mutation"),
            (
                BaseCoverageClosePartitionV1::Unsupported,
                5_u16,
                "unsupported",
            ),
            (
                BaseCoverageClosePartitionV1::Inapplicable,
                6_u16,
                "inapplicable",
            ),
        ];
        let mut counts = [0_u32; 6];
        for cell in close_manifest.cells() {
            let index = usize::from(cell.partition().code() - 1);
            counts[index] = counts[index].checked_add(1).expect("bounded count");
        }
        let mut frame =
            crate::canonical::CanonicalFrameV1::new(b"FSPHASE1PARTITIONS\x01", 2 * 1024 * 1024)
                .expect("independent full-ledger frame");
        frame
            .push_bytes(
                "partitions.source_manifest_root",
                close_manifest.source_manifest_root().as_bytes(),
            )
            .expect("source-manifest root");
        frame
            .push_bytes(
                "partitions.reason_registry_root",
                close_manifest.reason_registry_root().as_bytes(),
            )
            .expect("reason-registry root");
        frame
            .push_bytes(
                "partitions.close_manifest_root",
                close_manifest.root().as_bytes(),
            )
            .expect("close-manifest root");
        frame
            .push_bytes(
                "partitions.compatible_source_snapshot_root",
                source_closure.snapshot_root().as_bytes(),
            )
            .expect("snapshot root");
        frame
            .push_bytes(
                "partitions.source_closure_root",
                source_closure.root().as_bytes(),
            )
            .expect("source-closure root");
        frame
            .push_bytes(
                "partitions.vocabulary_root",
                vocabulary_only_phase_one_partition_root_v1().as_bytes(),
            )
            .expect("independent vocabulary root");
        frame
            .push_u32("partitions.count", 6)
            .expect("partition vocabulary count");
        for (index, (_, code, name)) in partition_literals.iter().copied().enumerate() {
            frame
                .push_u16("partitions.code", code)
                .expect("partition code");
            frame
                .push_str("partitions.name", name)
                .expect("partition name");
            frame
                .push_u32("partitions.expected_cell_count", counts[index])
                .expect("partition cell count");
            frame
                .push_u16("partitions.eventual_owner_requires_exact_match", 1)
                .expect("exact-match flag");
        }
        frame
            .push_u32(
                "partitions.cell_count",
                u32::try_from(close_manifest.cells().len()).expect("bounded cell count"),
            )
            .expect("cell count");
        for (index, cell) in close_manifest.cells().iter().enumerate() {
            let ordinal = u32::try_from(index + 1).expect("bounded ledger ordinal");
            let cell_id = format!("phase1-ledger-cell-{ordinal}");
            frame
                .push_u32("partitions.cell.ledger_ordinal", ordinal)
                .expect("ledger ordinal");
            frame
                .push_u32("partitions.cell.ordinal", cell.source_ordinal())
                .expect("source ordinal");
            frame
                .push_str("partitions.cell.cell_id", &cell_id)
                .expect("cell ID");
            frame
                .push_str("partitions.cell.stable_id", cell.source_case_id())
                .expect("source case ID");
            frame
                .push_u16("partitions.cell.source_class", cell.source_class() as u16)
                .expect("source class");
            frame
                .push_str("partitions.cell.source_path", cell.source_path())
                .expect("source path");
            frame
                .push_u16("partitions.cell.group", cell.group().code())
                .expect("group");
            frame
                .push_u16("partitions.cell.facet", cell.facet().code())
                .expect("facet");
            frame
                .push_u16(
                    "partitions.cell.execution_scope",
                    cell.execution_scope().code(),
                )
                .expect("execution scope");
            frame
                .push_u16("partitions.cell.partition", cell.partition().code())
                .expect("partition");
            frame
                .push_u16(
                    "partitions.cell.expected_decision",
                    cell.expected_decision().code(),
                )
                .expect("expected decision");
            let expected_result_root = cell.expected_reason().is_none().then(|| {
                independent_phase_one_expected_result_root_v1(
                    ordinal,
                    &cell_id,
                    cell,
                    close_manifest,
                    source_closure,
                )
            });
            frame
                .push_presence(
                    "partitions.cell.expected_result_present",
                    expected_result_root.is_some(),
                )
                .expect("expected-result presence");
            if let Some(root) = expected_result_root {
                frame
                    .push_bytes("partitions.cell.expected_result_root", root.as_bytes())
                    .expect("expected-result root");
            }
            frame
                .push_presence(
                    "partitions.cell.expected_reason_present",
                    cell.expected_reason().is_some(),
                )
                .expect("reason presence");
            if let Some(reason) = cell.expected_reason() {
                frame
                    .push_u16("partitions.cell.expected_reason", reason.code())
                    .expect("reason");
            }
            frame
                .push_str(
                    "partitions.cell.semantic_journey",
                    "runner-v2-phase1-contract",
                )
                .expect("semantic journey");
            frame
                .push_str(
                    "partitions.cell.route_id",
                    RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
                )
                .expect("route ID");
            frame
                .push_str(
                    "partitions.cell.route_owner",
                    RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                )
                .expect("route owner");
            frame
                .push_bytes(
                    "partitions.cell.compatible_source_snapshot_root",
                    source_closure.snapshot_root().as_bytes(),
                )
                .expect("row snapshot root");
            frame
                .push_bytes(
                    "partitions.cell.case_manifest_root",
                    close_manifest.root().as_bytes(),
                )
                .expect("row case-manifest root");
            frame
                .push_bytes(
                    "partitions.cell.source_closure_root",
                    source_closure.root().as_bytes(),
                )
                .expect("row source-closure root");
            let downstream_root = cell
                .downstream_contribution()
                .map(crate::coverage::BaseCoverageCloseDownstreamContributionV1::root);
            frame
                .push_presence(
                    "partitions.cell.downstream_contribution_present",
                    downstream_root.is_some(),
                )
                .expect("downstream contribution presence");
            if let Some(root) = downstream_root {
                frame
                    .push_bytes(
                        "partitions.cell.downstream_contribution_root",
                        root.as_bytes(),
                    )
                    .expect("downstream contribution root");
            }
            frame
                .push_bytes(
                    "partitions.cell.five_explicits_root",
                    cell.five_explicits().root().as_bytes(),
                )
                .expect("Five Explicits root");
            frame
                .push_bytes("partitions.cell.root", cell.root().as_bytes())
                .expect("cell root");
            frame
                .push_str(
                    "partitions.cell.no_claim",
                    "phase1-ledger-row-is-expected-contract-not-execution-actual-roots-require-real-observation",
                )
                .expect("row no-claim");
        }
        frame
            .push_u32("partitions.partition_id_set_count", 6)
            .expect("partition ID-set count");
        for (partition, _, _) in partition_literals {
            let members = close_manifest
                .cells()
                .iter()
                .enumerate()
                .filter(|(_, cell)| cell.partition() == partition)
                .collect::<Vec<_>>();
            frame
                .push_u16("partitions.id_set.partition", partition.code())
                .expect("ID-set partition");
            frame
                .push_u32(
                    "partitions.id_set.cell_count",
                    u32::try_from(members.len()).expect("bounded ID-set count"),
                )
                .expect("ID-set cell count");
            for (index, cell) in members {
                let ordinal = u32::try_from(index + 1).expect("bounded member ordinal");
                frame
                    .push_u32("partitions.id_set.ledger_ordinal", ordinal)
                    .expect("ID-set ordinal");
                frame
                    .push_str(
                        "partitions.id_set.cell_id",
                        &format!("phase1-ledger-cell-{ordinal}"),
                    )
                    .expect("ID-set cell ID");
                frame
                    .push_str("partitions.id_set.source_case_id", cell.source_case_id())
                    .expect("ID-set source case ID");
            }
        }
        frame
            .push_u32("partitions.required_unexpected_mismatches", 0)
            .expect("zero mismatches");
        frame
            .push_u32("partitions.required_execution_failures", 0)
            .expect("zero failures");
        frame
            .push_u32("partitions.required_unexplained_skips", 0)
            .expect("zero skips");
        frame
            .push_str(
                "partitions.no_claim",
                "expected-partitions-do-not-prove-owner-execution-or-matched-results",
            )
            .expect("ledger no-claim");
        frame.root("org.frankensim.fs-evidence-runner.phase-one-expected-partitions.v1")
    }

    #[test]
    fn schema_impact_v2_projection_has_exact_rows_paths_counts_and_canonical_roots() {
        let projection =
            RunnerV2BaseSchemaImpactProjectionV2::frozen().expect("exact V2 projection");
        let reconstructed =
            RunnerV2BaseSchemaImpactProjectionV2::reconstruct_exact(projection.manifest())
                .expect("exact manifest reconstructs");

        assert_eq!(reconstructed, projection);
        assert_eq!(
            projection.rows().len(),
            RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_ROW_COUNT_V2
        );
        assert_eq!(
            projection.source_path_count(),
            RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SOURCE_PATH_COUNT_V2 as u32
        );
        assert_eq!(
            projection.owned_row_count(),
            RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_OWNED_ROW_COUNT_V2 as u32
        );
        assert_eq!(
            projection.consumed_row_count(),
            RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_CONSUMED_ROW_COUNT_V2 as u32
        );
        assert_eq!(projection.execution_observation_count(), 0);
        assert_eq!(
            projection.manifest().no_claim().as_str(),
            RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1
        );

        let mut projection_roots = std::collections::BTreeSet::new();
        for (index, (row, manifest_entry)) in projection
            .rows()
            .iter()
            .zip(projection.manifest().entries())
            .enumerate()
        {
            assert_eq!(
                row.local_ordinal(),
                u32::try_from(index + 1).expect("bounded")
            );
            assert_eq!(row.relation(), SchemaImpactManifestRelationV1::Owned);
            assert_eq!(
                row.schema_id().as_str(),
                RUNNER_V2_BASE_SCHEMA_IMPACT_PROJECTION_SCHEMA_IDS_V2[index]
            );
            assert_eq!(
                row.source_path().as_str(),
                RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1
            );
            assert_eq!(
                row.disposition(),
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor
            );
            assert_eq!(
                row.schema_impact_row_root(),
                manifest_entry.row().root().content_hash()
            );
            assert_eq!(
                row.compatible_source_snapshot_root(),
                projection.manifest().compatible_source_snapshot_root()
            );
            assert_eq!(
                row.root(),
                schema_impact_projection_row_root_v2(
                    row.local_ordinal(),
                    row.relation(),
                    row.schema_id(),
                    row.source_path(),
                    row.disposition(),
                    row.schema_impact_row_root(),
                    row.compatible_source_snapshot_root(),
                )
                .expect("exact row root reconstructs")
            );
            assert!(projection_roots.insert(row.root()));
        }
        assert_eq!(projection_roots.len(), projection.rows().len());
        assert_eq!(
            projection.root(),
            schema_impact_projection_root_v2(
                projection.manifest(),
                projection.rows(),
                projection.source_path_count(),
                projection.owned_row_count(),
                projection.consumed_row_count(),
            )
            .expect("exact aggregate root reconstructs")
        );
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "one adversarial matrix checks every schema-impact projection coordinate and exact aggregate boundary"
    )]
    fn schema_impact_v2_projection_rejects_missing_extra_reordered_path_count_and_root_mutants() {
        let projection =
            RunnerV2BaseSchemaImpactProjectionV2::frozen().expect("exact V2 projection");
        let manifest = projection.manifest();
        let source_closure = RunnerV2BaseSourceClosureV1::frozen().expect("source closure");
        let alternate_manifest = source_frozen_schema_impact_manifest_v1(
            manifest.issuer_leaf_id().clone(),
            source_closure.compatible_snapshot(),
            manifest.frozen_base_registry(),
            Vec::new(),
            manifest
                .entries()
                .iter()
                .map(|entry| SchemaImpactManifestRowSourceV1 {
                    relation: entry.relation(),
                    row: entry.row().clone(),
                })
                .collect(),
            SchemaImpactNoClaimV1::new("alternate-valid-no-claim-boundary")
                .expect("alternate no-claim"),
        )
        .expect("alternate manifest remains internally valid");
        assert_eq!(
            RunnerV2BaseSchemaImpactProjectionV2::reconstruct_exact(&alternate_manifest)
                .expect_err("any non-exact manifest must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        assert_eq!(
            schema_impact_projection_root_v2(manifest, &projection.rows()[..2], 1, 2, 0)
                .expect_err("missing row must refuse")
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut extra = projection.rows().to_vec();
        extra.push(projection.rows()[2].clone());
        assert_eq!(
            schema_impact_projection_root_v2(manifest, &extra, 1, 4, 0)
                .expect_err("extra row must refuse")
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let mut reordered = projection.rows().to_vec();
        reordered.swap(0, 1);
        assert_eq!(
            schema_impact_projection_root_v2(manifest, &reordered, 1, 3, 0)
                .expect_err("reordered rows must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let mut wrong_path = projection.rows().to_vec();
        wrong_path[0].source_path =
            LogicalBundlePathV1::new("crates/fs-evidence-runner/src/projection.rs")
                .expect("mutated path");
        assert_eq!(
            schema_impact_projection_root_v2(manifest, &wrong_path, 2, 3, 0)
                .expect_err("path substitution must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_schema_row_root = projection.rows().to_vec();
        wrong_schema_row_root[0].schema_impact_row_root =
            hash_domain("projection-test-schema-impact-row-mutant.v1", b"row");
        assert_eq!(
            schema_impact_projection_root_v2(manifest, &wrong_schema_row_root, 1, 3, 0)
                .expect_err("schema row-root substitution must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let mut wrong_projection_row_root = projection.rows().to_vec();
        wrong_projection_row_root[0].root =
            hash_domain("projection-test-schema-projection-row-mutant.v1", b"row");
        assert_eq!(
            schema_impact_projection_root_v2(manifest, &wrong_projection_row_root, 1, 3, 0)
                .expect_err("projection row-root substitution must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(
            schema_impact_projection_root_v2(manifest, projection.rows(), 2, 3, 0)
                .expect_err("source-path count inflation must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            schema_impact_projection_root_v2(manifest, projection.rows(), 1, 2, 1)
                .expect_err("owned/consumed count substitution must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let row = &projection.rows()[0];
        assert_eq!(
            schema_impact_projection_row_root_v2(
                0,
                row.relation(),
                row.schema_id(),
                row.source_path(),
                row.disposition(),
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect_err("zero ordinal must refuse")
            .kind(),
            ConstructionErrorKindV2::Zero
        );
        let mut alternate_snapshot_frame =
            crate::canonical::CanonicalFrameV1::new(b"FSBASESOURCESNAPSHOT\x01", 256)
                .expect("alternate snapshot frame");
        alternate_snapshot_frame
            .push_str("snapshot.mutant", "one")
            .expect("alternate snapshot field");
        let alternate_snapshot =
            compatible_source_snapshot_root_from_exact_frame_v1(&alternate_snapshot_frame)
                .expect("alternate nominal snapshot");
        let alternate_schema_id = StableTokenV2::new("canonical-schema-field-descriptor-mutant")
            .expect("alternate schema ID");
        let alternate_source_path =
            LogicalBundlePathV1::new("crates/fs-evidence-runner/src/projection.rs")
                .expect("alternate source path");
        let variants = [
            schema_impact_projection_row_root_v2(
                row.local_ordinal() + 1,
                row.relation(),
                row.schema_id(),
                row.source_path(),
                row.disposition(),
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect("ordinal variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                SchemaImpactManifestRelationV1::Consumed,
                row.schema_id(),
                row.source_path(),
                row.disposition(),
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect("relation variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                row.relation(),
                &alternate_schema_id,
                row.source_path(),
                row.disposition(),
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect("schema-ID variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                row.relation(),
                row.schema_id(),
                &alternate_source_path,
                row.disposition(),
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect("source-path variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                row.relation(),
                row.schema_id(),
                row.source_path(),
                CanonicalSchemaImpactDispositionV1::UnchangedV1,
                row.schema_impact_row_root(),
                row.compatible_source_snapshot_root(),
            )
            .expect("disposition variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                row.relation(),
                row.schema_id(),
                row.source_path(),
                row.disposition(),
                hash_domain("projection-test-schema-impact-root-axis.v1", b"mutant"),
                row.compatible_source_snapshot_root(),
            )
            .expect("schema-root variant"),
            schema_impact_projection_row_root_v2(
                row.local_ordinal(),
                row.relation(),
                row.schema_id(),
                row.source_path(),
                row.disposition(),
                row.schema_impact_row_root(),
                alternate_snapshot,
            )
            .expect("snapshot variant"),
        ];
        assert!(variants.iter().all(|root| *root != row.root()));
        assert_eq!(
            variants
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            variants.len()
        );
    }

    #[test]
    fn phase_one_contribution_is_exact_result_free_deferred_and_no_execution() {
        let projection =
            RunnerV2PhaseOneContractContributionV2::frozen().expect("phase-one contribution");
        let contribution = projection.contribution();
        let observer = contribution.observer_contract();
        let envelope = projection.deferred_envelope();

        assert_eq!(
            RunnerV2PhaseOneContractContributionV2::frozen().expect("deterministic reconstruction"),
            projection
        );
        assert_eq!(
            contribution.source_owner(),
            RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1
        );
        assert_eq!(
            observer.route_id().as_str(),
            RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1
        );
        assert_eq!(
            observer.semantic_consumer(),
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
        );
        assert_eq!(
            observer.execution_owner(),
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
        );
        assert_eq!(
            observer.driver_owner(),
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
        );
        assert_eq!(
            observer.posix_wrapper_owner(),
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
        );
        assert_eq!(
            observer.windows_wrapper_owner(),
            RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
        );
        assert_eq!(
            observer.capability_profile(),
            BaseCoverageCloseCapabilityProfileV1::ReleaseControl
        );
        assert_eq!(
            observer.capability_profile().stable_id(),
            "fs-evidence-runner.close-capability.release-control.v1"
        );
        assert_eq!(
            observer.driver().as_str(),
            RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1
        );
        assert_eq!(
            observer.posix_route(),
            RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1
        );
        assert_eq!(
            observer.windows_route(),
            RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1
        );
        assert_eq!(
            observer.case_manifest_path(),
            RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1
        );
        assert_eq!(
            observer.case_manifest_root(),
            projection.case_manifest_contract_root()
        );
        assert_eq!(
            observer.deferred_reason(),
            DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        );
        assert_eq!(
            observer.deferred_reason().descriptor().name(),
            "immutable-contribution-awaits-designated-release-owner"
        );
        assert_eq!(
            observer.no_claim().as_str(),
            RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1
        );
        assert_eq!(contribution.payload_root(), projection.payload_root());
        assert_eq!(
            contribution.expected_partitions_root(),
            projection.expected_partitions_root()
        );
        assert_eq!(
            contribution.schema_root(),
            projection.schema_impact_projection().root()
        );
        let schema_impact_log_schema_root =
            crate::logging::schema_impact_log_schema_root_v1().expect("schema-impact log schema");
        assert_eq!(
            contribution.log_schema_root(),
            schema_impact_log_schema_root
        );
        let source_closure = RunnerV2BaseSourceClosureV1::frozen().expect("exact source closure");
        assert_eq!(
            projection.payload_root(),
            phase_one_payload_root_v2(
                &source_closure,
                projection.base_coverage_manifest_root(),
                projection.schema_impact_projection(),
                projection.expected_partitions_root(),
                contribution.target_root(),
                contribution.feature_set_root(),
                schema_impact_log_schema_root,
            )
            .expect("independent exact Phase-1 payload oracle")
        );
        assert_ne!(
            projection.payload_root(),
            phase_one_payload_root_v2(
                &source_closure,
                projection.base_coverage_manifest_root(),
                projection.schema_impact_projection(),
                projection.expected_partitions_root(),
                contribution.target_root(),
                contribution.feature_set_root(),
                crate::logging::base_e2e_log_schema_root_v1()
                    .expect("compatibility base-E2E log schema"),
            )
            .expect("compatibility-root mutant payload")
        );
        assert_eq!(contribution.per_cell_budgets().max_child_processes(), 1);
        assert_eq!(contribution.per_cell_budgets().max_parallel_children(), 1);
        assert_eq!(contribution.shard_budgets().max_child_processes(), 4);
        assert_eq!(contribution.shard_budgets().max_parallel_children(), 2);
        assert_eq!(
            contribution.retained_artifact_policy(),
            BaseCoverageCloseRetainedRelativeArtifactPolicyV1::OwnerEnvelopeRelativePathsOnly
        );
        assert_eq!(
            contribution.no_claim().as_str(),
            RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1
        );
        assert_eq!(
            envelope.evidence_kind(),
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
        );
        assert_eq!(
            envelope.disposition(),
            RuntimeObservationDispositionV1::Deferred
        );
        assert_eq!(envelope.contribution_root(), contribution.root());
        assert_eq!(envelope.observer_contract_root(), observer.root());
        assert_eq!(
            envelope.deferred_reason_registry_root(),
            observer.deferred_reason_registry_root()
        );
        assert_eq!(envelope.deferred_reason(), observer.deferred_reason());
        assert_eq!(
            envelope.retained_artifact_policy(),
            contribution.retained_artifact_policy()
        );
        assert_eq!(
            envelope.no_claim().as_str(),
            RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1
        );
        assert_eq!(contribution.execution_actual_field_count(), 0);
        assert_eq!(envelope.execution_actual_field_count(), 0);
        assert_eq!(projection.downstream_execution_count(), 0);
        assert_ne!(contribution.root(), envelope.root().content_hash());
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "one exhaustive mutation matrix proves every contribution identity axis and boundary is bound or refused"
    )]
    fn phase_one_contribution_roots_move_and_budget_owner_cycle_boundaries_fail_closed() {
        let projection =
            RunnerV2PhaseOneContractContributionV2::frozen().expect("phase-one contribution");
        let base = PhaseOneContributionFixture::from_contribution(projection.contribution());
        let original_root = projection.contribution().root();
        let mut variants = Vec::new();

        for axis in 0_u8..14 {
            let mut fixture = base.clone();
            match axis {
                0 => fixture.source_owner = "frankensim-test-mutated-source-owner".to_owned(),
                1 => {
                    fixture.payload_root =
                        hash_domain("projection-test-phase-one-payload-axis.v1", b"mutant")
                }
                2 => {
                    fixture.expected_partitions_root =
                        hash_domain("projection-test-phase-one-partitions-axis.v1", b"mutant")
                }
                3 => {
                    fixture.per_cell_budgets = frozen_downstream_close_contribution_budgets_v1(2, 1)
                        .expect("mutated per-cell budget")
                }
                4 => {
                    fixture.shard_budgets = frozen_downstream_close_contribution_budgets_v1(5, 2)
                        .expect("mutated shard budget")
                }
                5 => {
                    fixture.schema_root =
                        hash_domain("projection-test-phase-one-schema-axis.v1", b"mutant")
                }
                6 => {
                    fixture.log_schema_root =
                        hash_domain("projection-test-phase-one-log-axis.v1", b"mutant")
                }
                7 => {
                    fixture.source_root = super::parse_phase_one_source_identity(hash_domain(
                        "projection-test-phase-one-source-axis.v1",
                        b"mutant",
                    ))
                    .expect("mutated source identity")
                }
                8 => {
                    fixture.build_root = super::parse_phase_one_build_identity(hash_domain(
                        "projection-test-phase-one-build-axis.v1",
                        b"mutant",
                    ))
                    .expect("mutated build identity")
                }
                9 => {
                    fixture.toolchain_root = super::parse_phase_one_toolchain_identity(hash_domain(
                        "projection-test-phase-one-toolchain-axis.v1",
                        b"mutant",
                    ))
                    .expect("mutated toolchain identity")
                }
                10 => {
                    fixture.target_root =
                        hash_domain("projection-test-phase-one-target-axis.v1", b"mutant")
                }
                11 => {
                    fixture.feature_set_root =
                        hash_domain("projection-test-phase-one-feature-axis.v1", b"mutant")
                }
                12 => {
                    fixture.observer_contract = BaseCoverageCloseDeferredObservationContractV1::new(
                        RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                        BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
                        StableTokenV2::new(RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1).expect("driver"),
                        RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
                        RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
                        projection.case_manifest_contract_root(),
                        DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
                        StableTokenV2::new("mutated-phase1-observer-no-claim")
                            .expect("mutated observer no-claim"),
                    )
                    .expect("structurally valid mutated observer")
                }
                13 => {
                    fixture.no_claim = StableTokenV2::new("mutated-phase1-contribution-no-claim")
                        .expect("mutated contribution no-claim")
                }
                _ => unreachable!("closed mutation axis"),
            }
            let variant = fixture.build().expect("individually valid mutation");
            assert_ne!(
                variant.root(),
                original_root,
                "axis {axis} must move the root"
            );
            variants.push(variant.root());
        }
        assert_eq!(
            variants
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            variants.len()
        );

        let close_manifest =
            BaseCoverageCloseManifestV1::frozen().expect("exact full-set close manifest");
        let source_closure = RunnerV2BaseSourceClosureV1::frozen().expect("exact source closure");
        let exact_ledger =
            phase_one_close_ledger_cells_v1(&close_manifest).expect("exact Phase-1 close ledger");
        assert_eq!(exact_ledger.len(), close_manifest.cells().len());
        let mut exact_cell_ids = std::collections::BTreeSet::new();
        let mut exact_case_ids = std::collections::BTreeSet::new();
        for (index, cell) in exact_ledger.iter().enumerate() {
            let ordinal = u32::try_from(index + 1).expect("bounded ledger ordinal");
            assert_eq!(cell.ledger_ordinal, ordinal);
            assert_eq!(cell.source_ordinal, ordinal);
            assert_eq!(
                cell.cell_id.as_str(),
                format!("phase1-ledger-cell-{ordinal}")
            );
            assert!(exact_cell_ids.insert(cell.cell_id.as_str()));
            assert!(exact_case_ids.insert(cell.source_case_id.as_ref()));
            assert_ne!(
                cell.expected_result_root.is_some(),
                cell.expected_reason.is_some(),
                "every row carries exactly one expected result root or registered reason"
            );
            assert_eq!(
                cell.semantic_journey.as_str(),
                super::RUNNER_V2_PHASE_ONE_SEMANTIC_JOURNEY_V1
            );
            assert_eq!(
                cell.route_id.as_str(),
                RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1
            );
            assert_eq!(
                cell.route_owner.as_ref(),
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1
            );
            assert_eq!(
                cell.compatible_source_snapshot_root,
                source_closure.snapshot_root()
            );
            assert_eq!(cell.case_manifest_root, close_manifest.root());
            assert_eq!(cell.source_closure_root, source_closure.root());
            assert_eq!(
                cell.no_claim.as_str(),
                super::RUNNER_V2_PHASE_ONE_LEDGER_ROW_NO_CLAIM_V1
            );
        }
        assert_eq!(exact_cell_ids.len(), exact_ledger.len());
        assert_eq!(exact_case_ids.len(), exact_ledger.len());
        assert_eq!(projection.expected_ledger().cells(), exact_ledger);
        assert_eq!(
            projection.expected_ledger(),
            &super::RunnerV2PhaseOneExpectedLedgerV1::frozen()
                .expect("public frozen Phase-1 expected ledger")
        );
        assert_eq!(
            projection.expected_ledger().source_manifest_root(),
            close_manifest.source_manifest_root()
        );
        assert_eq!(
            projection.expected_ledger().reason_registry_root(),
            close_manifest.reason_registry_root()
        );
        assert_eq!(
            projection.expected_ledger().close_manifest_root(),
            close_manifest.root()
        );
        assert_eq!(
            projection
                .expected_ledger()
                .compatible_source_snapshot_root(),
            source_closure.snapshot_root()
        );
        assert_eq!(
            projection.expected_ledger().source_closure_root(),
            source_closure.root()
        );
        assert_eq!(
            projection.expected_ledger().partition_vocabulary_root(),
            vocabulary_only_phase_one_partition_root_v1(),
            "the production vocabulary root must match an independent literal oracle"
        );
        assert_eq!(
            projection.expected_ledger().execution_observation_count(),
            0
        );
        for partition in BaseCoverageClosePartitionV1::ALL {
            assert_eq!(
                usize::try_from(projection.expected_ledger().partition_cell_count(partition))
                    .expect("bounded partition count"),
                exact_ledger
                    .iter()
                    .filter(|cell| cell.partition == partition)
                    .count()
            );
        }
        assert_eq!(
            projection.case_manifest_contract_root(),
            phase_one_case_manifest_contract_root_from_components_v1(
                close_manifest.source_manifest_root(),
                close_manifest.reason_registry_root(),
                close_manifest.root(),
                u32::try_from(close_manifest.cells().len()).expect("bounded close cell count"),
                projection.expected_partitions_root(),
            )
            .expect("exact case-manifest contract root")
        );
        assert_eq!(
            projection.expected_partitions_root(),
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &exact_ledger,
            )
            .expect("exact full-ledger partition root")
        );
        assert_eq!(
            projection.expected_partitions_root(),
            independent_phase_one_full_ledger_root_v1(&close_manifest, &source_closure),
            "the production ledger root must match the handwritten literal encoder"
        );
        assert_ne!(
            projection.expected_partitions_root(),
            vocabulary_only_phase_one_partition_root_v1(),
            "the six-name vocabulary alone cannot identify the full ledger"
        );
        assert_ne!(
            projection.expected_partitions_root(),
            count_only_phase_one_partition_root_v1(&close_manifest, &source_closure),
            "exact aggregate counts alone cannot identify ordered ledger rows"
        );

        let exact_cell_count =
            u32::try_from(close_manifest.cells().len()).expect("bounded close cell count");
        let case_manifest_root_variants = [
            phase_one_case_manifest_contract_root_from_components_v1(
                hash_domain(
                    "projection-test-phase-one-source-manifest-axis.v1",
                    b"mutant",
                ),
                close_manifest.reason_registry_root(),
                close_manifest.root(),
                exact_cell_count,
                projection.expected_partitions_root(),
            )
            .expect("source-manifest-root variant"),
            phase_one_case_manifest_contract_root_from_components_v1(
                close_manifest.source_manifest_root(),
                hash_domain(
                    "projection-test-phase-one-reason-registry-axis.v1",
                    b"mutant",
                ),
                close_manifest.root(),
                exact_cell_count,
                projection.expected_partitions_root(),
            )
            .expect("reason-registry-root variant"),
            phase_one_case_manifest_contract_root_from_components_v1(
                close_manifest.source_manifest_root(),
                close_manifest.reason_registry_root(),
                hash_domain(
                    "projection-test-phase-one-close-manifest-axis.v1",
                    b"mutant",
                ),
                exact_cell_count,
                projection.expected_partitions_root(),
            )
            .expect("close-manifest-root variant"),
            phase_one_case_manifest_contract_root_from_components_v1(
                close_manifest.source_manifest_root(),
                close_manifest.reason_registry_root(),
                close_manifest.root(),
                exact_cell_count
                    .checked_add(1)
                    .expect("bounded count variant"),
                projection.expected_partitions_root(),
            )
            .expect("close-manifest-count variant"),
            phase_one_case_manifest_contract_root_from_components_v1(
                close_manifest.source_manifest_root(),
                close_manifest.reason_registry_root(),
                close_manifest.root(),
                exact_cell_count,
                hash_domain(
                    "projection-test-phase-one-expected-ledger-axis.v1",
                    b"mutant",
                ),
            )
            .expect("expected-ledger-root variant"),
        ];
        assert!(
            case_manifest_root_variants
                .iter()
                .all(|root| *root != projection.case_manifest_contract_root())
        );
        assert_eq!(
            case_manifest_root_variants
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            case_manifest_root_variants.len()
        );

        let mut missing_ledger_cell = exact_ledger.clone();
        missing_ledger_cell.pop();
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &missing_ledger_cell,
            )
            .expect_err("a missing full-ledger row must refuse")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut reordered_ledger = exact_ledger.clone();
        reordered_ledger.swap(0, 1);
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &reordered_ledger,
            )
            .expect_err("reordered full-ledger rows must refuse")
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let mut reclassified_ledger = exact_ledger.clone();
        reclassified_ledger[0].source_class =
            if reclassified_ledger[0].source_class == BaseCoverageManifestClassV1::Unit {
                BaseCoverageManifestClassV1::Boundary
            } else {
                BaseCoverageManifestClassV1::Unit
            };
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &reclassified_ledger,
            )
            .expect_err("source reclassification must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_cell_root_ledger = exact_ledger.clone();
        wrong_cell_root_ledger[0].cell_root =
            hash_domain("projection-test-phase-one-cell-root-axis.v1", b"mutant");
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_cell_root_ledger,
            )
            .expect_err("close-cell root substitution must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let expected_result_index = exact_ledger
            .iter()
            .position(|cell| cell.expected_result_root.is_some())
            .expect("at least one exact expected-result row");
        let mut wrong_result_ledger = exact_ledger.clone();
        wrong_result_ledger[expected_result_index].expected_result_root = Some(hash_domain(
            "projection-test-phase-one-expected-result-axis.v1",
            b"mutant",
        ));
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_result_ledger,
            )
            .expect_err("wrong expected-result root must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut decision_root_mismatch_ledger = exact_ledger.clone();
        let original_expected_decision =
            decision_root_mismatch_ledger[expected_result_index].expected_decision;
        decision_root_mismatch_ledger[expected_result_index].expected_decision =
            BaseCoverageCloseDecisionV1::ALL
                .into_iter()
                .find(|decision| *decision != original_expected_decision)
                .expect("more than one closed decision");
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &decision_root_mismatch_ledger,
            )
            .expect_err("decision/root mismatch must refuse before close reporting")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let reason_index = exact_ledger
            .iter()
            .position(|cell| cell.expected_reason.is_some())
            .expect("at least one registered-reason row");
        let original_reason = exact_ledger[reason_index]
            .expected_reason
            .expect("registered reason");
        let alternate_reason = BaseCoverageCloseReasonCodeV1::ALL
            .into_iter()
            .find(|reason| *reason != original_reason)
            .expect("more than one registered reason");
        let mut wrong_reason_ledger = exact_ledger.clone();
        wrong_reason_ledger[reason_index].expected_reason = Some(alternate_reason);
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_reason_ledger,
            )
            .expect_err("wrong registered reason must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mut wrong_owner_ledger = exact_ledger.clone();
        wrong_owner_ledger[0].route_owner = "frankensim-test-wrong-route-owner".into();
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_owner_ledger,
            )
            .expect_err("wrong route owner must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_route_ledger = exact_ledger.clone();
        wrong_route_ledger[0].route_id =
            StableTokenV2::new("runner-v2.route.mutant.v1").expect("valid mutant route");
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_route_ledger,
            )
            .expect_err("wrong route must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut stale_snapshot_ledger = exact_ledger.clone();
        stale_snapshot_ledger[0].compatible_source_snapshot_root = hash_domain(
            "projection-test-phase-one-stale-snapshot-axis.v1",
            b"mutant",
        );
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &stale_snapshot_ledger,
            )
            .expect_err("mixed compatible source snapshot must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_source_closure_ledger = exact_ledger.clone();
        wrong_source_closure_ledger[0].source_closure_root = hash_domain(
            "projection-test-phase-one-source-closure-axis.v1",
            b"mutant",
        );
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_source_closure_ledger,
            )
            .expect_err("wrong source-closure root must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_row_no_claim_ledger = exact_ledger.clone();
        wrong_row_no_claim_ledger[0].no_claim =
            StableTokenV2::new("mutant-phase1-ledger-row-no-claim")
                .expect("valid mutant row no-claim");
        assert_eq!(
            phase_one_expected_partitions_root_from_ledger_v1(
                &close_manifest,
                &source_closure,
                &wrong_row_no_claim_ledger,
            )
            .expect_err("wrong row no-claim must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mut every_row_root_mutant = exact_ledger.clone();
        for index in 0..every_row_root_mutant.len() {
            let original = every_row_root_mutant[index].cell_root;
            every_row_root_mutant[index].cell_root = hash_domain(
                "projection-test-phase-one-every-row-root-axis.v1",
                &u32::try_from(index)
                    .expect("bounded row index")
                    .to_be_bytes(),
            );
            assert_eq!(
                phase_one_expected_partitions_root_from_ledger_v1(
                    &close_manifest,
                    &source_closure,
                    &every_row_root_mutant,
                )
                .expect_err("every individual close-row root mutant must refuse")
                .kind(),
                ConstructionErrorKindV2::Incompatible,
                "row index {index}"
            );
            every_row_root_mutant[index].cell_root = original;
        }

        assert_eq!(
            frozen_downstream_close_contribution_budgets_v1(0, 1)
                .expect_err("zero child budget must refuse")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            frozen_downstream_close_contribution_budgets_v1(1, 2)
                .expect_err("parallel greater than total must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            frozen_downstream_close_contribution_budgets_v1(257, 1)
                .expect_err("child budget above ceiling must refuse")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut inverted_budget_scope = base.clone();
        inverted_budget_scope.per_cell_budgets =
            frozen_downstream_close_contribution_budgets_v1(5, 3)
                .expect("valid standalone per-cell budget");
        assert_eq!(
            inverted_budget_scope
                .build()
                .expect_err("per-cell scope wider than shard must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut colliding_roots = base.clone();
        colliding_roots.feature_set_root = colliding_roots.target_root;
        assert_eq!(
            colliding_roots
                .build()
                .expect_err("semantic-root collision must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut owner_alias = base.clone();
        owner_alias.source_owner = RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1.to_owned();
        assert_eq!(
            owner_alias
                .build()
                .expect_err("source/execution owner alias must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseDeferredObservationContractV1::new(
                RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
                StableTokenV2::new(RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1).expect("driver"),
                RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
                projection.case_manifest_contract_root(),
                DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
                StableTokenV2::new(RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1)
                    .expect("observer no-claim"),
            )
            .expect_err("upstream semantic consumer must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseDeferredObservationContractV1::new(
                RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
                StableTokenV2::new(RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1).expect("driver"),
                RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
                projection.case_manifest_contract_root(),
                DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
                StableTokenV2::new(RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1)
                    .expect("observer no-claim"),
            )
            .expect_err("duplicate platform routes must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let alternate_observer_contribution = {
            let mut fixture = base.clone();
            fixture.observer_contract = BaseCoverageCloseDeferredObservationContractV1::new(
                RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1,
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
                StableTokenV2::new(RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1).expect("driver"),
                RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1,
                RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1,
                projection.case_manifest_contract_root(),
                DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
                StableTokenV2::new("mutated-phase1-observer-no-claim")
                    .expect("mutated observer no-claim"),
            )
            .expect("structurally valid alternate observer");
            fixture.build().expect("alternate observer contribution")
        };
        assert_eq!(
            validate_exact_phase_one_observer_contract_v1(
                alternate_observer_contribution.observer_contract()
            )
            .expect_err("exact Phase-1 validator must reject observer substitution")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let alternate_observer_envelope = BaseCoverageCloseDeferredEvidenceEnvelopeV1::new(
            &alternate_observer_contribution,
            StableTokenV2::new(RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1)
                .expect("envelope no-claim"),
        )
        .expect("structurally valid alternate-observer envelope");
        assert_eq!(
            phase_one_contribution_projection_root_v2(
                projection.schema_impact_projection(),
                projection.source_closure_root(),
                projection.base_coverage_manifest_root(),
                projection.case_manifest_contract_root(),
                projection.expected_partitions_root(),
                projection.payload_root(),
                &alternate_observer_contribution,
                &alternate_observer_envelope,
            )
            .expect_err("aggregate must reject observer substitution")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            phase_one_contribution_projection_root_v2(
                projection.schema_impact_projection(),
                projection.source_closure_root(),
                projection.base_coverage_manifest_root(),
                projection.case_manifest_contract_root(),
                projection.expected_partitions_root(),
                hash_domain("projection-test-phase-one-join-mutant.v1", b"payload"),
                projection.contribution(),
                projection.deferred_envelope(),
            )
            .expect_err("aggregate must reject a payload/contribution mismatch")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let wrong_log_contribution = {
            let mut fixture = base.clone();
            fixture.log_schema_root = crate::logging::base_e2e_log_schema_root_v1()
                .expect("compatibility base-E2E log schema");
            fixture
                .build()
                .expect("structurally valid wrong-log contribution")
        };
        let wrong_log_envelope = BaseCoverageCloseDeferredEvidenceEnvelopeV1::new(
            &wrong_log_contribution,
            StableTokenV2::new(RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1)
                .expect("envelope no-claim"),
        )
        .expect("wrong-log Deferred envelope");
        assert_eq!(
            phase_one_contribution_projection_root_v2(
                projection.schema_impact_projection(),
                projection.source_closure_root(),
                projection.base_coverage_manifest_root(),
                projection.case_manifest_contract_root(),
                projection.expected_partitions_root(),
                projection.payload_root(),
                &wrong_log_contribution,
                &wrong_log_envelope,
            )
            .expect_err("aggregate must reject the compatibility log schema root")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mutant_manifest = runner_v2_base_schema_impact_field_name_mutant_manifest_v1()
            .expect("production-shape field-name mutant manifest");
        let exact_manifest = projection.schema_impact_projection().manifest();
        assert_eq!(
            mutant_manifest.issuer_leaf_id(),
            exact_manifest.issuer_leaf_id()
        );
        assert_eq!(
            mutant_manifest.compatible_source_snapshot_root(),
            exact_manifest.compatible_source_snapshot_root()
        );
        assert_eq!(
            mutant_manifest.frozen_base_registry(),
            exact_manifest.frozen_base_registry()
        );
        assert_eq!(
            mutant_manifest.leaf_extension_registries(),
            exact_manifest.leaf_extension_registries()
        );
        assert_eq!(
            mutant_manifest.graph_edge_count(),
            exact_manifest.graph_edge_count()
        );
        assert_eq!(mutant_manifest.no_claim(), exact_manifest.no_claim());
        assert_ne!(mutant_manifest.root(), exact_manifest.root());
        assert_eq!(
            &mutant_manifest.entries()[1..],
            &exact_manifest.entries()[1..],
            "the frame-descriptor and version-slot meta rows remain exact siblings"
        );

        let exact_field_row = exact_manifest.entries()[0].row();
        let mutant_field_row = mutant_manifest.entries()[0].row();
        assert_eq!(mutant_field_row.schema_id(), exact_field_row.schema_id());
        assert_eq!(
            mutant_field_row.disposition(),
            exact_field_row.disposition()
        );
        assert_eq!(
            mutant_field_row.migration_policy(),
            exact_field_row.migration_policy()
        );
        assert_eq!(
            mutant_field_row.prior_frame(),
            exact_field_row.prior_frame()
        );
        assert_eq!(
            mutant_field_row.legacy_container(),
            exact_field_row.legacy_container()
        );
        assert_eq!(
            mutant_field_row.owner_leaf_id(),
            exact_field_row.owner_leaf_id()
        );
        assert_eq!(
            mutant_field_row.source_path(),
            exact_field_row.source_path()
        );
        assert_eq!(
            mutant_field_row.authority_surfaces(),
            exact_field_row.authority_surfaces()
        );
        assert_eq!(
            mutant_field_row.construction_predecessors(),
            exact_field_row.construction_predecessors()
        );
        assert_eq!(
            mutant_field_row.legal_parent_slots(),
            exact_field_row.legal_parent_slots()
        );
        assert_eq!(
            mutant_field_row.legal_child_slots(),
            exact_field_row.legal_child_slots()
        );
        assert_eq!(
            mutant_field_row.compatible_source_snapshot_root(),
            exact_field_row.compatible_source_snapshot_root()
        );
        assert_eq!(mutant_field_row.no_claim(), exact_field_row.no_claim());
        assert_ne!(mutant_field_row.root(), exact_field_row.root());

        let exact_field_binding = exact_field_row
            .authoritative_frame()
            .expect("exact field meta frame");
        let mutant_field_binding = mutant_field_row
            .authoritative_frame()
            .expect("mutant field meta frame");
        assert_eq!(
            mutant_field_binding.authority_state(),
            exact_field_binding.authority_state()
        );
        let exact_field_frame = exact_field_binding.descriptor();
        let mutant_field_frame = mutant_field_binding.descriptor();
        assert_eq!(
            mutant_field_frame.rust_schema_name(),
            exact_field_frame.rust_schema_name()
        );
        assert_eq!(
            mutant_field_frame.frame_version(),
            exact_field_frame.frame_version()
        );
        assert_eq!(mutant_field_frame.domain(), exact_field_frame.domain());
        assert_eq!(mutant_field_frame.magic(), exact_field_frame.magic());
        assert_eq!(
            mutant_field_frame.nominal_role(),
            exact_field_frame.nominal_role()
        );
        assert_eq!(
            mutant_field_frame.fields().len(),
            exact_field_frame.fields().len()
        );
        assert_eq!(
            &mutant_field_frame.fields()[1..],
            &exact_field_frame.fields()[1..],
            "only the first field descriptor changes"
        );
        let exact_first_field = &exact_field_frame.fields()[0];
        let mutant_first_field = &mutant_field_frame.fields()[0];
        assert_eq!(exact_first_field.field_name().as_str(), "ordinal");
        assert_eq!(mutant_first_field.field_name().as_str(), "ordinal-mutant");
        assert_eq!(mutant_first_field.ordinal(), exact_first_field.ordinal());
        assert_eq!(
            mutant_first_field.field_code(),
            exact_first_field.field_code()
        );
        assert_eq!(
            mutant_first_field.semantic_type_id(),
            exact_first_field.semantic_type_id()
        );
        assert_eq!(
            mutant_first_field.wire_kind(),
            exact_first_field.wire_kind()
        );
        assert_eq!(mutant_first_field.layout(), exact_first_field.layout());
        assert_eq!(
            mutant_first_field.related_field_code(),
            exact_first_field.related_field_code()
        );
        assert_eq!(
            mutant_first_field.version_slot_code(),
            exact_first_field.version_slot_code()
        );
        assert_ne!(
            mutant_first_field.descriptor_identity(),
            exact_first_field.descriptor_identity()
        );
        assert_ne!(
            mutant_field_frame.descriptor_identity(),
            exact_field_frame.descriptor_identity()
        );

        let mutant_schema_projection =
            RunnerV2BaseSchemaImpactProjectionV2::build_exact(mutant_manifest)
                .expect("production-shape mutant projects through exact V2 shape validation");
        assert_eq!(
            mutant_schema_projection.source_path_count(),
            projection.schema_impact_projection().source_path_count()
        );
        assert_eq!(
            mutant_schema_projection.owned_row_count(),
            projection.schema_impact_projection().owned_row_count()
        );
        assert_eq!(
            mutant_schema_projection.consumed_row_count(),
            projection.schema_impact_projection().consumed_row_count()
        );
        assert_ne!(
            mutant_schema_projection.rows()[0].schema_impact_row_root(),
            projection.schema_impact_projection().rows()[0].schema_impact_row_root()
        );
        assert_ne!(
            mutant_schema_projection.rows()[0].root(),
            projection.schema_impact_projection().rows()[0].root()
        );
        assert_eq!(
            &mutant_schema_projection.rows()[1..],
            &projection.schema_impact_projection().rows()[1..]
        );
        assert_ne!(
            mutant_schema_projection.root(),
            projection.schema_impact_projection().root()
        );
        assert_eq!(mutant_schema_projection.execution_observation_count(), 0);

        let mutant_payload_root = phase_one_payload_root_v2(
            &source_closure,
            projection.base_coverage_manifest_root(),
            &mutant_schema_projection,
            projection.expected_partitions_root(),
            projection.contribution().target_root(),
            projection.contribution().feature_set_root(),
            projection.contribution().log_schema_root(),
        )
        .expect("mutant projection propagates into the Phase-1 payload");
        assert_ne!(mutant_payload_root, projection.payload_root());
        let mutant_build_identity = super::parse_phase_one_build_identity(
            phase_one_build_contract_root_v1(&source_closure, &mutant_schema_projection)
                .expect("mutant projection propagates into the build contract"),
        )
        .expect("mutant build identity");
        assert_ne!(
            &mutant_build_identity,
            projection.contribution().build_root()
        );
        let mutant_contribution = {
            let mut fixture = base.clone();
            fixture.payload_root = mutant_payload_root;
            fixture.schema_root = mutant_schema_projection.root();
            fixture.build_root = mutant_build_identity;
            fixture
                .build()
                .expect("mutant meta-field propagates into a result-free contribution")
        };
        assert_ne!(mutant_contribution.root(), projection.contribution().root());
        assert_eq!(mutant_contribution.payload_root(), mutant_payload_root);
        assert_eq!(
            mutant_contribution.schema_root(),
            mutant_schema_projection.root()
        );
        assert_eq!(mutant_contribution.execution_actual_field_count(), 0);
        let mutant_deferred_envelope = BaseCoverageCloseDeferredEvidenceEnvelopeV1::new(
            &mutant_contribution,
            StableTokenV2::new(RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1)
                .expect("exact Deferred no-claim"),
        )
        .expect("result-free mutant Deferred envelope");
        assert_eq!(
            mutant_deferred_envelope.disposition(),
            RuntimeObservationDispositionV1::Deferred
        );
        assert_eq!(
            mutant_deferred_envelope.evidence_kind(),
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
        );
        assert_eq!(
            mutant_deferred_envelope.deferred_reason(),
            DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        );
        assert_eq!(mutant_deferred_envelope.execution_actual_field_count(), 0);
        assert_ne!(
            mutant_deferred_envelope.root(),
            projection.deferred_envelope().root()
        );
        let mutant_aggregate_root = phase_one_contribution_projection_root_v2(
            &mutant_schema_projection,
            projection.source_closure_root(),
            projection.base_coverage_manifest_root(),
            projection.case_manifest_contract_root(),
            projection.expected_partitions_root(),
            mutant_payload_root,
            &mutant_contribution,
            &mutant_deferred_envelope,
        )
        .expect("mutant contribution remains a valid result-free Deferred aggregate");
        assert_ne!(mutant_aggregate_root, projection.root());

        let alternate_envelope = BaseCoverageCloseDeferredEvidenceEnvelopeV1::new(
            projection.contribution(),
            StableTokenV2::new("alternate-deferred-envelope-no-claim")
                .expect("alternate envelope no-claim"),
        )
        .expect("alternate declaration-only envelope");
        assert_eq!(projection.contribution().root(), original_root);
        assert_ne!(
            alternate_envelope.root(),
            projection.deferred_envelope().root()
        );
        assert_eq!(
            alternate_envelope.contribution_root(),
            projection.contribution().root()
        );
    }

    fn presented_results(
        projection: &RunnerV2BaseE2eProjectionV1,
        journey: BaseE2eJourneyV1,
        harness: &BaseE2eHarnessIdentityV1,
    ) -> Vec<BaseE2ePresentedRowResultV1> {
        projection
            .journeys()
            .iter()
            .find(|manifest| manifest.journey() == journey)
            .expect("journey manifest")
            .rows()
            .iter()
            .map(|row| {
                let execution = execute_case(row.kind(), harness);
                let counts = BaseE2eObservedCountsV1::new(
                    BaseE2eMatchedPartitionV1::new(
                        execution.positive_eligible,
                        execution.positive_matched,
                    )
                    .expect("valid positive partition"),
                    BaseE2eMatchedPartitionV1::new(
                        execution.expected_refusals,
                        execution.expected_refusals_matched,
                    )
                    .expect("valid refusal partition"),
                    execution.unsupported,
                    execution.unexpected_mismatches,
                )
                .expect("reconciled observed counts");
                let first_unexpected_cell = execution
                    .first_failed_cell
                    .as_deref()
                    .map(|id| StableTokenV2::new(id).expect("stable failed-cell ID"));
                BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                    row,
                    execution.decision,
                    counts,
                    execution.detail.observed,
                    execution
                        .detail
                        .observed_cells
                        .as_deref()
                        .expect("in-process execution retains observed detail cells"),
                    first_unexpected_cell,
                )
                .expect("comparison-only public presentation")
            })
            .collect()
    }

    fn reconstruct_result(
        source: &BaseE2ePresentedRowResultV1,
        journey: BaseE2eJourneyV1,
        row_id: StableTokenV2,
        semantic_manifest_root: fs_blake3::ContentHash,
    ) -> BaseE2ePresentedRowResultV1 {
        BaseE2ePresentedRowResultV1::new_with_detail_descriptor(
            journey,
            row_id,
            semantic_manifest_root,
            source.observed(),
            source.counts(),
            source.observed_detail_manifest(),
            source.detail_cells_matched(),
            source.first_unexpected_cell().cloned(),
        )
        .expect("mutated presented result")
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "this exhaustive contract test audits every frozen journey row, semantic root, and independent oracle mutation in one matrix"
    )]
    fn manifest_exactly_maps_five_scripts_and_all_base_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        assert_eq!(projection.journeys().len(), BaseE2eJourneyV1::ALL.len());
        let expected_row_counts = [21_usize, 18, 19, 24, 16];
        let expected_routes = [
            (
                BaseE2eJourneyV1::PublicationState,
                "publication-state-v2",
                "frankensim-epic-foundations-huq.24.2.2.2",
                "e2e-evidence-runner-publication-state-v2-driver",
                "scripts/ci/e2e_evidence_runner_publication_state_v2.sh",
                "scripts/ci/manifests/evidence_runner_publication_state_v2_cases.v1.json",
            ),
            (
                BaseE2eJourneyV1::PublicationV2,
                "publication-v2",
                "frankensim-epic-foundations-huq.24.2.2.3.3",
                "e2e-evidence-runner-publication-v2-driver",
                "scripts/ci/e2e_evidence_runner_publication_v2.sh",
                "scripts/ci/manifests/evidence_runner_publication_v2_cases.v1.json",
            ),
            (
                BaseE2eJourneyV1::VerifierV2,
                "verifier-v2",
                "frankensim-epic-foundations-huq.24.3.3.3.3",
                "e2e-evidence-verifier-v2-driver",
                "scripts/ci/e2e_evidence_verifier_v2.sh",
                "scripts/ci/manifests/evidence_verifier_v2_cases.v1.json",
            ),
            (
                BaseE2eJourneyV1::CanonicalRunnerV2,
                "canonical-runner-v2",
                "frankensim-epic-foundations-huq.24.4.1.4",
                "canonical-evidence-runner-v2-e2e-driver",
                "scripts/ci/canonical_evidence_runner_v2.sh",
                "scripts/ci/manifests/canonical_evidence_runner_v2_cases.v1.json",
            ),
            (
                BaseE2eJourneyV1::RjoqHandoffV1,
                "rjoq-handoff-v1",
                "frankensim-epic-foundations-huq.24.5.3.1",
                "verify-runner-rjoq-handoff-v1-driver",
                "scripts/ci/verify_runner_rjoq_handoff_v1.sh",
                "scripts/ci/manifests/runner_rjoq_handoff_verifier_v1_cases.v1.json",
            ),
        ];
        let mut row_count = 0_usize;
        let mut detail_manifests = std::collections::BTreeMap::new();
        for (index, journey) in projection.journeys().iter().enumerate() {
            let (expected_journey, key, owner, driver, script, manifest) = expected_routes[index];
            assert_eq!(journey.journey(), expected_journey);
            assert_eq!(journey.journey(), BaseE2eJourneyV1::ALL[index]);
            assert_eq!(journey.journey().key(), key);
            assert_eq!(journey.downstream_owner().as_str(), owner);
            assert_eq!(journey.driver_name().as_str(), driver);
            assert_eq!(journey.script_path().as_str(), script);
            assert_eq!(journey.manifest_path().as_str(), manifest);
            assert_eq!(journey.rows().len(), expected_row_counts[index]);
            assert_eq!(
                journey
                    .rows()
                    .iter()
                    .map(|row| row.id().as_str())
                    .collect::<Vec<_>>(),
                journey_row_id_oracle(journey.journey())
            );
            assert!(journey.rows().iter().all(|row| {
                row.journey() == journey.journey()
                    && row.downstream_owner() == journey.downstream_owner()
                    && row.downstream_driver() == journey.driver_name()
                    && row.downstream_script() == journey.script_path()
                    && row.downstream_manifest() == journey.manifest_path()
                    && row.source_closure_root() == projection.source_closure().root()
                    && row.log_schema_root() == projection.log_schema_root()
                    && !row.consumption_rationale().as_str().is_empty()
                    && !row.fixture_reference().as_str().is_empty()
            }));
            for row in journey.rows() {
                let descriptor = row.expected_detail_manifest();
                let cells = row.expected_detail_cells();
                assert_eq!(
                    descriptor.cell_count(),
                    row.expected_refusal_cell_count() + row.unsupported_cell_count()
                );
                assert_eq!(descriptor.cell_count(), row.expected_detail_cell_count());
                assert_eq!(descriptor.root(), row.expected_detail_manifest_root());
                assert_eq!(
                    usize::try_from(descriptor.cell_count()).expect("bounded detail count"),
                    cells.len()
                );
                assert_eq!(
                    super::detail_manifest_from_cells(row.kind(), cells),
                    descriptor
                );
                assert!(
                    cells
                        .windows(2)
                        .all(|pair| { pair[0].semantic_ordinal() < pair[1].semantic_ordinal() })
                );
                assert!(cells.iter().all(|cell| {
                    cell.kind() == row.kind()
                        && !cell.stable_id().is_empty()
                        && cell.cell_root() == cell.root()
                }));
                assert_eq!(
                    cells
                        .iter()
                        .map(super::BaseE2eDetailCellV1::stable_id)
                        .collect::<std::collections::BTreeSet<_>>()
                        .len(),
                    cells.len()
                );
                assert_eq!(
                    cells
                        .iter()
                        .map(|cell| cell.root().to_hex())
                        .collect::<std::collections::BTreeSet<_>>()
                        .len(),
                    cells.len()
                );
                let substituted_oracle_root = hash_domain(
                    "projection-test-substituted-oracle-manifest.v1",
                    row.id().as_str().as_bytes(),
                );
                let oracle_substitution = super::semantic_row_root(
                    row.kind(),
                    row.expected(),
                    row.expected_detail(),
                    row.semantic_cell_count(),
                    row.positive_cell_count(),
                    row.expected_refusal_cell_count(),
                    row.unsupported_cell_count(),
                    row.expected_detail_cell_count(),
                    row.expected_detail_manifest_root(),
                    substituted_oracle_root,
                    row.registered_decision_detail(),
                    row.unit(),
                    row.no_claim_scope(),
                )
                .expect("oracle-root substitution remains structurally encodable");
                assert_ne!(oracle_substitution, row.semantic_manifest_root());
                match row.kind() {
                    BaseE2eCaseKindV1::LimitCatalog => {
                        assert_eq!(cells.len(), 142);
                        assert_eq!(
                            cells
                                .iter()
                                .map(super::BaseE2eDetailCellV1::semantic_ordinal)
                                .collect::<Vec<_>>(),
                            (1_u32..=142).map(|index| index * 2).collect::<Vec<_>>(),
                            "limit refusals follow the exact global smoke-then-full matrix order"
                        );
                        assert_eq!(
                            cells.first().map(super::BaseE2eDetailCellV1::stable_id),
                            Some("limit.smoke.argv_tokens.one-over")
                        );
                        assert_eq!(
                            cells.get(70).map(super::BaseE2eDetailCellV1::stable_id),
                            Some("limit.smoke.registered_resource_identities_per_family.one-over")
                        );
                        assert_eq!(
                            cells.get(71).map(super::BaseE2eDetailCellV1::stable_id),
                            Some("limit.full.argv_tokens.one-over")
                        );
                        assert_eq!(
                            cells.last().map(super::BaseE2eDetailCellV1::stable_id),
                            Some("limit.full.registered_resource_identities_per_family.one-over")
                        );
                        assert!(cells.iter().all(|cell| matches!(
                            cell.payload(),
                            super::BaseE2eDetailPayloadV1::Limit {
                                owner: "fs-evidence-runner.runner-limits",
                                repair_rank: 1,
                                repair_target,
                                ..
                            } if !repair_target.is_empty()
                        )));
                    }
                    BaseE2eCaseKindV1::BudgetAdmission => {
                        assert_eq!(cells.len(), 8);
                        assert!(cells.iter().all(|cell| matches!(
                            cell.payload(),
                            super::BaseE2eDetailPayloadV1::Budget {
                                kind: RunnerBudgetViolationKindV2::ProfileCeilingExceeded,
                                owner: "fs-evidence-runner.runner-budgets",
                                repair_rank: 1,
                                repair_kind: RepairActionKindV2::ReduceResourceDemand,
                                repair_target,
                                ..
                            } if !repair_target.is_empty()
                        )));
                    }
                    BaseE2eCaseKindV1::BudgetChildRelation => {
                        assert!(!cells.is_empty());
                        assert!(cells.iter().all(|cell| matches!(
                            cell.payload(),
                            super::BaseE2eDetailPayloadV1::Budget {
                                owner: "fs-evidence-runner.runner-budgets",
                                repair_rank: 1,
                                repair_kind: RepairActionKindV2::ReduceResourceDemand,
                                repair_target: "max_parallel_children",
                                ..
                            }
                        )));
                    }
                    BaseE2eCaseKindV1::Diagnostic => {
                        let registered = row
                            .registered_decision_detail()
                            .expect("diagnostic carries the bounded downstream detail reference");
                        assert_eq!(registered.namespace().code(), 7);
                        assert_eq!(registered.detail_code(), 1);
                        assert!(registered.encoded_length() > 0);
                        assert_eq!(
                            registered.registry_root(),
                            crate::catalog::DecisionDetailNamespaceRegistryV2::frozen().root()
                        );
                        let without_registered = super::semantic_row_root(
                            row.kind(),
                            row.expected(),
                            row.expected_detail(),
                            row.semantic_cell_count(),
                            row.positive_cell_count(),
                            row.expected_refusal_cell_count(),
                            row.unsupported_cell_count(),
                            row.expected_detail_cell_count(),
                            row.expected_detail_manifest_root(),
                            row.oracle_manifest_root(),
                            None,
                            row.unit(),
                            row.no_claim_scope(),
                        )
                        .expect("semantic root without registered reference");
                        assert_ne!(without_registered, row.semantic_manifest_root());
                    }
                    _ => assert!(row.registered_decision_detail().is_none()),
                }
                if let Some(previous) = detail_manifests.insert(row.kind(), descriptor) {
                    assert_eq!(
                        previous, descriptor,
                        "one case kind must have one journey-independent detail manifest"
                    );
                }
            }
            row_count += journey.rows().len();
        }
        assert_eq!(row_count, 98);
        assert_eq!(detail_manifests.len(), BaseE2eCaseKindV1::ALL.len());
        let distinct_detail_roots = detail_manifests
            .values()
            .map(|descriptor| descriptor.root().to_hex())
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(
            distinct_detail_roots.len(),
            BaseE2eCaseKindV1::ALL.len(),
            "the manifest root binds the case kind even for empty manifests"
        );

        let assert_oracle_mutation =
            |kind: BaseE2eCaseKindV1,
             owner: &str,
             original_table_root: fs_blake3::ContentHash,
             mutated_table_root: fs_blake3::ContentHash| {
                assert_ne!(
                    mutated_table_root,
                    original_table_root,
                    "changing one actual literal must move the table root for {}",
                    kind.name()
                );
                let mutated_oracle_root = super::case_oracle_manifest_root_from_table_root(
                    kind,
                    owner,
                    mutated_table_root,
                );
                assert_ne!(
                    mutated_oracle_root,
                    super::case_oracle_manifest_root(kind),
                    "changing one actual literal must move the containing oracle root for {}",
                    kind.name()
                );
                let row = projection
                    .journeys()
                    .iter()
                    .flat_map(super::BaseE2eJourneyProjectionV1::rows)
                    .find(|row| row.kind() == kind)
                    .expect("every oracle-backed case kind has a projection row");
                let mutated_semantic_root = super::semantic_row_root(
                    row.kind(),
                    row.expected(),
                    row.expected_detail(),
                    row.semantic_cell_count(),
                    row.positive_cell_count(),
                    row.expected_refusal_cell_count(),
                    row.unsupported_cell_count(),
                    row.expected_detail_cell_count(),
                    row.expected_detail_manifest_root(),
                    mutated_oracle_root,
                    row.registered_decision_detail(),
                    row.unit(),
                    row.no_claim_scope(),
                )
                .expect("one literal mutation remains structurally encodable");
                assert_ne!(
                    mutated_semantic_root,
                    row.semantic_manifest_root(),
                    "the semantic row root must bind the mutated oracle root for {}",
                    kind.name()
                );
            };

        macro_rules! assert_catalog_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::catalog_literal_oracle_rows();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::CatalogLiterals,
                    "catalog",
                    super::catalog_oracle_table_root(&super::catalog_literal_oracle_rows()),
                    super::catalog_oracle_table_root(&rows),
                );
            }};
        }
        assert_catalog_column!(0, 9_u8);
        assert_catalog_column!(1, "api-generation-mutant");
        assert_catalog_column!(2, 3_u16);
        assert_catalog_column!(3, "RunnerSpecV2-mutant-column");
        assert_catalog_column!(4, true);

        macro_rules! assert_limit_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::limit_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::LimitCatalog,
                    "limits",
                    super::limit_oracle_table_root(super::limit_oracle_rows()),
                    super::limit_oracle_table_root(&rows),
                );
            }};
        }
        assert_limit_column!(0, crate::limits::RunnerLimitFieldV2::ArgvTokenBytes);
        assert_limit_column!(1, 2_u16);
        assert_limit_column!(2, "argv_tokens_mutant");
        assert_limit_column!(3, crate::limits::RunnerLimitUnitV2::Rows);
        assert_limit_column!(4, crate::limits::RunnerLimitWidthV2::U64);
        assert_limit_column!(5, crate::limits::RunnerLimitTightenabilityV2::Fixed);
        assert_limit_column!(6, RunnerLimitValueV2::U32(65));
        assert_limit_column!(7, RunnerLimitValueV2::U32(66));

        let mut mutated_catalog_literals = super::catalog_literal_oracle_rows();
        mutated_catalog_literals[0].3 = "RunnerSpecV2-mutant";
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CatalogLiterals,
            "catalog",
            super::catalog_oracle_table_root(&super::catalog_literal_oracle_rows()),
            super::catalog_oracle_table_root(&mutated_catalog_literals),
        );

        let mut mutated_limit_literals = super::limit_oracle_rows().to_vec();
        mutated_limit_literals[0].6 = RunnerLimitValueV2::U32(65);
        assert_oracle_mutation(
            BaseE2eCaseKindV1::LimitCatalog,
            "limits",
            super::limit_oracle_table_root(super::limit_oracle_rows()),
            super::limit_oracle_table_root(&mutated_limit_literals),
        );
        let mut mutated_limit_identity = super::limit_oracle_rows().to_vec();
        mutated_limit_identity[0].0 = crate::limits::RunnerLimitFieldV2::ArgvTokenBytes;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::LimitCatalog,
            "limits",
            super::limit_oracle_table_root(super::limit_oracle_rows()),
            super::limit_oracle_table_root(&mutated_limit_identity),
        );

        macro_rules! assert_budget_field_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::budget_field_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::BudgetAdmission,
                    "budgets",
                    super::budget_oracle_table_root(),
                    super::budget_oracle_table_root_from_rows(
                        &rows,
                        super::budget_logical_unit_oracle_rows(),
                        super::budget_profile_oracle_rows(),
                        super::budget_profile_refusal_oracle_rows(),
                    ),
                );
            }};
        }
        assert_budget_field_column!(0, crate::budget::RunnerBudgetFieldV2::MaxResidentBytes);
        assert_budget_field_column!(1, 2_u16);
        assert_budget_field_column!(2, "wall_time_ns_mutant");
        assert_budget_field_column!(3, crate::budget::RunnerBudgetUnitV2::LogicalBytes);
        assert_budget_field_column!(4, crate::budget::RunnerBudgetWidthV2::U32);
        assert_budget_field_column!(5, crate::budget::RunnerBudgetValueV2::U64(100_000_000_001));

        macro_rules! assert_budget_unit_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::budget_logical_unit_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::BudgetAdmission,
                    "budgets",
                    super::budget_oracle_table_root(),
                    super::budget_oracle_table_root_from_rows(
                        super::budget_field_oracle_rows(),
                        &rows,
                        super::budget_profile_oracle_rows(),
                        super::budget_profile_refusal_oracle_rows(),
                    ),
                );
            }};
        }
        assert_budget_unit_column!(0, 2_u16);
        assert_budget_unit_column!(1, "encoded-bytes-mutant-column");
        assert_budget_unit_column!(2, true);

        macro_rules! assert_budget_profile_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::budget_profile_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::BudgetAdmission,
                    "budgets",
                    super::budget_oracle_table_root(),
                    super::budget_oracle_table_root_from_rows(
                        super::budget_field_oracle_rows(),
                        super::budget_logical_unit_oracle_rows(),
                        &rows,
                        super::budget_profile_refusal_oracle_rows(),
                    ),
                );
            }};
        }
        assert_budget_profile_column!(0, crate::catalog::RunProfileV2::Full);
        assert_budget_profile_column!(1, 2_u16);
        assert_budget_profile_column!(2, "smoke-mutant");
        assert_budget_profile_column!(3, 900_000_000_001_u64);
        assert_budget_profile_column!(4, 17_179_869_185_u64);
        assert_budget_profile_column!(5, 33_u32);
        assert_budget_profile_column!(6, 257_u32);

        macro_rules! assert_budget_refusal_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::budget_profile_refusal_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::BudgetAdmission,
                    "budgets",
                    super::budget_oracle_table_root(),
                    super::budget_oracle_table_root_from_rows(
                        super::budget_field_oracle_rows(),
                        super::budget_logical_unit_oracle_rows(),
                        super::budget_profile_oracle_rows(),
                        &rows,
                    ),
                );
            }};
        }
        assert_budget_refusal_column!(0, 38_u32);
        assert_budget_refusal_column!(1, crate::catalog::RunProfileV2::Full);
        assert_budget_refusal_column!(2, 2_u16);
        assert_budget_refusal_column!(3, "smoke-mutant-refusal");
        assert_budget_refusal_column!(4, crate::budget::RunnerBudgetFieldV2::MaxResidentBytes);
        assert_budget_refusal_column!(5, 2_u16);
        assert_budget_refusal_column!(6, "wall_time_ns_mutant");
        assert_budget_refusal_column!(7, crate::budget::RunnerBudgetUnitV2::LogicalBytes);
        assert_budget_refusal_column!(8, crate::budget::RunnerBudgetValueV2::U64(900_000_000_002));
        assert_budget_refusal_column!(9, crate::budget::RunnerBudgetValueV2::U64(900_000_000_003));

        let mut mutated_budget_literals = super::budget_field_oracle_rows().to_vec();
        mutated_budget_literals[0].5 = crate::budget::RunnerBudgetValueV2::U64(100_000_000_001);
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                &mutated_budget_literals,
                super::budget_logical_unit_oracle_rows(),
                super::budget_profile_oracle_rows(),
                super::budget_profile_refusal_oracle_rows(),
            ),
        );
        let mut mutated_budget_field_identity = super::budget_field_oracle_rows().to_vec();
        mutated_budget_field_identity[0].0 = crate::budget::RunnerBudgetFieldV2::MaxResidentBytes;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                &mutated_budget_field_identity,
                super::budget_logical_unit_oracle_rows(),
                super::budget_profile_oracle_rows(),
                super::budget_profile_refusal_oracle_rows(),
            ),
        );
        let mut mutated_budget_units = super::budget_logical_unit_oracle_rows().to_vec();
        mutated_budget_units[0].1 = "encoded-bytes-mutant";
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                super::budget_field_oracle_rows(),
                &mutated_budget_units,
                super::budget_profile_oracle_rows(),
                super::budget_profile_refusal_oracle_rows(),
            ),
        );
        let mut mutated_budget_profiles = super::budget_profile_oracle_rows().to_vec();
        mutated_budget_profiles[0].0 = crate::catalog::RunProfileV2::Full;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                super::budget_field_oracle_rows(),
                super::budget_logical_unit_oracle_rows(),
                &mutated_budget_profiles,
                super::budget_profile_refusal_oracle_rows(),
            ),
        );
        let mut mutated_budget_refusal_profile =
            super::budget_profile_refusal_oracle_rows().to_vec();
        mutated_budget_refusal_profile[0].1 = crate::catalog::RunProfileV2::Full;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                super::budget_field_oracle_rows(),
                super::budget_logical_unit_oracle_rows(),
                super::budget_profile_oracle_rows(),
                &mutated_budget_refusal_profile,
            ),
        );
        let mut mutated_budget_refusal_field = super::budget_profile_refusal_oracle_rows().to_vec();
        mutated_budget_refusal_field[0].4 = crate::budget::RunnerBudgetFieldV2::MaxResidentBytes;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::BudgetAdmission,
            "budgets",
            super::budget_oracle_table_root(),
            super::budget_oracle_table_root_from_rows(
                super::budget_field_oracle_rows(),
                super::budget_logical_unit_oracle_rows(),
                super::budget_profile_oracle_rows(),
                &mutated_budget_refusal_field,
            ),
        );

        macro_rules! assert_diagnostic_code_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::diagnostic_code_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::Diagnostic,
                    "diagnostics",
                    super::diagnostic_oracle_table_root(),
                    super::diagnostic_oracle_table_root_from_rows(
                        &rows,
                        super::diagnostic_registered_code_oracle(),
                        super::diagnostic_retryability_oracle_rows(),
                        super::diagnostic_repair_kind_oracle_rows(),
                    ),
                );
            }};
        }
        assert_diagnostic_code_column!(0, 2_u16);
        assert_diagnostic_code_column!(1, "case.conformance_mismatch-mutant-column");

        for registered in [(8_u16, 9_u16), (7_u16, 10_u16)] {
            assert_oracle_mutation(
                BaseE2eCaseKindV1::Diagnostic,
                "diagnostics",
                super::diagnostic_oracle_table_root(),
                super::diagnostic_oracle_table_root_from_rows(
                    super::diagnostic_code_oracle_rows(),
                    registered,
                    super::diagnostic_retryability_oracle_rows(),
                    super::diagnostic_repair_kind_oracle_rows(),
                ),
            );
        }

        macro_rules! assert_diagnostic_retryability_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::diagnostic_retryability_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::Diagnostic,
                    "diagnostics",
                    super::diagnostic_oracle_table_root(),
                    super::diagnostic_oracle_table_root_from_rows(
                        super::diagnostic_code_oracle_rows(),
                        super::diagnostic_registered_code_oracle(),
                        &rows,
                        super::diagnostic_repair_kind_oracle_rows(),
                    ),
                );
            }};
        }
        assert_diagnostic_retryability_column!(0, 1_u16);
        assert_diagnostic_retryability_column!(1, "never-mutant");

        macro_rules! assert_diagnostic_repair_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::diagnostic_repair_kind_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::Diagnostic,
                    "diagnostics",
                    super::diagnostic_oracle_table_root(),
                    super::diagnostic_oracle_table_root_from_rows(
                        super::diagnostic_code_oracle_rows(),
                        super::diagnostic_registered_code_oracle(),
                        super::diagnostic_retryability_oracle_rows(),
                        &rows,
                    ),
                );
            }};
        }
        assert_diagnostic_repair_column!(0, 2_u16);
        assert_diagnostic_repair_column!(1, "change-arguments-mutant");

        let mut mutated_diagnostic_literals = super::diagnostic_code_oracle_rows().to_vec();
        mutated_diagnostic_literals[0].1 = "case.conformance_mismatch-mutant";
        assert_oracle_mutation(
            BaseE2eCaseKindV1::Diagnostic,
            "diagnostics",
            super::diagnostic_oracle_table_root(),
            super::diagnostic_oracle_table_root_from_rows(
                &mutated_diagnostic_literals,
                super::diagnostic_registered_code_oracle(),
                super::diagnostic_retryability_oracle_rows(),
                super::diagnostic_repair_kind_oracle_rows(),
            ),
        );

        macro_rules! assert_command_list_column {
            ($column:tt, $value:expr) => {{
                let mut list = super::command_list_oracle();
                list.$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::CommandList,
                    "commands",
                    super::command_oracle_table_root(),
                    super::command_oracle_table_root_from_rows(
                        list,
                        super::command_intent_oracle_rows(),
                        super::command_applicability_oracle_rows(),
                    ),
                );
            }};
        }
        assert_command_list_column!(0, crate::catalog::RunnerCommandV2::Check);
        assert_command_list_column!(1, 1_u16);
        assert_command_list_column!(2, "list-mutant-column");

        macro_rules! assert_command_intent_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::command_intent_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::CommandList,
                    "commands",
                    super::command_oracle_table_root(),
                    super::command_oracle_table_root_from_rows(
                        super::command_list_oracle(),
                        &rows,
                        super::command_applicability_oracle_rows(),
                    ),
                );
            }};
        }
        assert_command_intent_column!(0, crate::catalog::RunnerCommandV2::Run);
        assert_command_intent_column!(1, 2_u16);
        assert_command_intent_column!(2, "check-mutant-column");
        assert_command_intent_column!(3, crate::catalog::RunProfileV2::Full);
        assert_command_intent_column!(4, 2_u16);
        assert_command_intent_column!(5, "smoke-mutant-column");
        assert_command_intent_column!(
            6,
            crate::catalog::ArtifactDispositionV2::DurableBundleRequired
        );
        assert_command_intent_column!(7, 2_u16);
        assert_command_intent_column!(8, "lifecycle-only-mutant-column");

        macro_rules! assert_command_applicability_column {
            ($column:tt, $value:expr) => {{
                let mut rows = super::command_applicability_oracle_rows().to_vec();
                rows[0].$column = $value;
                assert_oracle_mutation(
                    BaseE2eCaseKindV1::CommandList,
                    "commands",
                    super::command_oracle_table_root(),
                    super::command_oracle_table_root_from_rows(
                        super::command_list_oracle(),
                        super::command_intent_oracle_rows(),
                        &rows,
                    ),
                );
            }};
        }
        assert_command_applicability_column!(0, crate::catalog::RunnerCommandV2::Run);
        assert_command_applicability_column!(1, 1_u16);
        assert_command_applicability_column!(2, "list-mutant-applicability");
        for requirement_index in 0..super::COMMAND_SELECTOR_FIELD_ORACLE_V1.len() {
            let mut rows = super::command_applicability_oracle_rows().to_vec();
            rows[0].3[requirement_index] = crate::command::CommandSelectorCardinalityV2::Singular;
            assert_oracle_mutation(
                BaseE2eCaseKindV1::CommandList,
                "commands",
                super::command_oracle_table_root(),
                super::command_oracle_table_root_from_rows(
                    super::command_list_oracle(),
                    super::command_intent_oracle_rows(),
                    &rows,
                ),
            );
        }

        let mut mutated_command_literals = super::command_intent_oracle_rows().to_vec();
        mutated_command_literals[0].2 = "check-mutant";
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                super::command_list_oracle(),
                &mutated_command_literals,
                super::command_applicability_oracle_rows(),
            ),
        );
        let mut mutated_command_identity = super::command_intent_oracle_rows().to_vec();
        mutated_command_identity[0].0 = crate::catalog::RunnerCommandV2::Run;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                super::command_list_oracle(),
                &mutated_command_identity,
                super::command_applicability_oracle_rows(),
            ),
        );
        let mut mutated_command_profile = super::command_intent_oracle_rows().to_vec();
        mutated_command_profile[0].3 = crate::catalog::RunProfileV2::Full;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                super::command_list_oracle(),
                &mutated_command_profile,
                super::command_applicability_oracle_rows(),
            ),
        );
        let mut mutated_command_disposition = super::command_intent_oracle_rows().to_vec();
        mutated_command_disposition[0].6 =
            crate::catalog::ArtifactDispositionV2::DurableBundleRequired;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                super::command_list_oracle(),
                &mutated_command_disposition,
                super::command_applicability_oracle_rows(),
            ),
        );
        let mut mutated_command_applicability = super::command_applicability_oracle_rows().to_vec();
        mutated_command_applicability[0].0 = crate::catalog::RunnerCommandV2::Run;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                super::command_list_oracle(),
                super::command_intent_oracle_rows(),
                &mutated_command_applicability,
            ),
        );
        let mut mutated_list_identity = super::command_list_oracle();
        mutated_list_identity.0 = crate::catalog::RunnerCommandV2::Check;
        assert_oracle_mutation(
            BaseE2eCaseKindV1::CommandList,
            "commands",
            super::command_oracle_table_root(),
            super::command_oracle_table_root_from_rows(
                mutated_list_identity,
                super::command_intent_oracle_rows(),
                super::command_applicability_oracle_rows(),
            ),
        );
    }

    #[test]
    fn journey_routes_reject_wrong_owner_driver_script_or_manifest() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let wrong_owner = StableTokenV2::new("wrong-owner").expect("wrong owner");
        let wrong_driver = StableTokenV2::new("wrong-driver").expect("wrong driver");
        let wrong_script =
            LogicalBundlePathV1::new("scripts/ci/wrong_route.sh").expect("wrong script");
        let wrong_manifest = LogicalBundlePathV1::new("scripts/ci/manifests/wrong_route.v1.json")
            .expect("wrong manifest");

        for journey in projection.journeys() {
            let exact = journey_root(
                journey.journey(),
                journey.downstream_owner(),
                journey.driver_name(),
                journey.script_path(),
                journey.manifest_path(),
                journey.rows(),
                journey.source_closure_root(),
                journey.log_schema_root(),
            )
            .expect("exact route");
            assert_eq!(exact, journey.manifest_root());

            for (mutation, result) in [
                (
                    "owner",
                    journey_root(
                        journey.journey(),
                        &wrong_owner,
                        journey.driver_name(),
                        journey.script_path(),
                        journey.manifest_path(),
                        journey.rows(),
                        journey.source_closure_root(),
                        journey.log_schema_root(),
                    ),
                ),
                (
                    "driver",
                    journey_root(
                        journey.journey(),
                        journey.downstream_owner(),
                        &wrong_driver,
                        journey.script_path(),
                        journey.manifest_path(),
                        journey.rows(),
                        journey.source_closure_root(),
                        journey.log_schema_root(),
                    ),
                ),
                (
                    "script",
                    journey_root(
                        journey.journey(),
                        journey.downstream_owner(),
                        journey.driver_name(),
                        &wrong_script,
                        journey.manifest_path(),
                        journey.rows(),
                        journey.source_closure_root(),
                        journey.log_schema_root(),
                    ),
                ),
                (
                    "manifest",
                    journey_root(
                        journey.journey(),
                        journey.downstream_owner(),
                        journey.driver_name(),
                        journey.script_path(),
                        &wrong_manifest,
                        journey.rows(),
                        journey.source_closure_root(),
                        journey.log_schema_root(),
                    ),
                ),
            ] {
                let error = result.expect_err(mutation);
                assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
                assert_eq!(error.field(), "projection.journey_row_context");
            }
        }
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "this end-to-end test validates every aggregate counter, root binding, and deterministic logging invariant together"
    )]
    fn all_real_constructor_rows_agree_and_logs_are_deterministic() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let first = run_base_e2e_projection_v1(&projection, &harness()).expect("projection run");
        let second = run_base_e2e_projection_v1(&projection, &harness()).expect("projection run");
        assert_eq!(first, second);
        assert!(first.is_green());
        assert_eq!(first.unexpected_mismatches(), 0);
        assert_eq!(
            [
                super::COMMAND_APPLICABILITY_SETUP_SEMANTIC_ORDINAL_V1,
                super::BUDGET_BASE_CONSTRUCTION_SEMANTIC_ORDINAL_V1,
                super::CAPABILITY_REGISTRY_SETUP_SEMANTIC_ORDINAL_V1,
            ],
            [1, 1, 1],
            "every pre-matrix setup refusal owns the first real semantic ordinal"
        );

        let catalog_inventory = super::case_detail_fields(BaseE2eCaseKindV1::CatalogLiterals)
            .expect("catalog expected-inventory fields");
        let catalog_inventory_cells = catalog_inventory
            .iter()
            .find(|field| field.name().as_str() == "catalog-literal-cells")
            .expect("catalog inventory exposes its expected cell count");
        assert_eq!(
            catalog_inventory_cells.value(),
            &crate::value::TypedValueV2::U32(186)
        );
        let observed_catalog_prefix = super::aggregate_accept(
            Err((
                1,
                "catalog.api-generation.synthetic-prefix-failure".to_owned(),
            )),
            super::execute_detail_manifest(BaseE2eCaseKindV1::CatalogLiterals, &harness()),
        );
        assert_eq!(observed_catalog_prefix.checked_cells, 1);
        assert_eq!(observed_catalog_prefix.positive_eligible, 1);
        assert_eq!(observed_catalog_prefix.positive_matched, 0);
        assert_ne!(
            observed_catalog_prefix.checked_cells, 186,
            "observed progress must not be replaced by expected inventory"
        );

        assert!(
            std::panic::catch_unwind(|| super::aggregate_accept(
                Err((0, "catalog.fabricated-zero-progress".to_owned())),
                super::expected_detail_execution(BaseE2eCaseKindV1::CatalogLiterals),
            ))
            .is_err(),
            "the case execution helper must flag a zero-cell failure"
        );
        assert!(
            std::panic::catch_unwind(|| super::mixed_progress_from_first_failure(
                Ok(0),
                44,
                super::budget_matrix_partition,
            ))
            .is_err(),
            "the mixed-progress helper must flag a zero-cell success"
        );
        assert!(
            std::panic::catch_unwind(|| super::mixed_progress_from_first_failure(
                Err((0, "budget.fabricated-zero-progress".to_owned())),
                44,
                super::budget_matrix_partition,
            ))
            .is_err(),
            "the mixed-progress helper must flag a zero-cell failure"
        );
        assert!(
            std::panic::catch_unwind(|| super::aggregate_mixed(
                super::BaseE2eMixedProgressV1::new(),
                36,
                8,
                0,
                super::expected_detail_execution(BaseE2eCaseKindV1::BudgetAdmission),
            ))
            .is_err(),
            "the mixed aggregate must flag absent observed progress"
        );

        let forced_positive = super::aggregate_mixed(
            super::mixed_progress_from_first_failure(
                Err((1, "budget.field.wall_time_ns".to_owned())),
                44,
                super::budget_matrix_partition,
            ),
            36,
            8,
            0,
            super::execute_detail_manifest(BaseE2eCaseKindV1::BudgetAdmission, &harness()),
        );
        assert_eq!(forced_positive.checked_cells, 1);
        assert_eq!(forced_positive.positive_eligible, 1);
        assert_eq!(forced_positive.positive_matched, 0);
        assert_eq!(forced_positive.expected_refusals, 0);
        assert_eq!(forced_positive.expected_refusals_matched, 0);
        assert_eq!(forced_positive.unexpected_mismatches, 1);
        assert_eq!(
            forced_positive.first_failed_cell.as_deref(),
            Some("budget.field.wall_time_ns")
        );
        let forced_refusal = super::aggregate_mixed(
            super::mixed_progress_from_first_failure(
                Err((37, "budget.profile.smoke.wall_time_ns.one-over".to_owned())),
                44,
                super::budget_matrix_partition,
            ),
            36,
            8,
            0,
            super::execute_detail_manifest(BaseE2eCaseKindV1::BudgetAdmission, &harness()),
        );
        assert_eq!(forced_refusal.checked_cells, 37);
        assert_eq!(forced_refusal.positive_eligible, 36);
        assert_eq!(forced_refusal.positive_matched, 36);
        assert_eq!(forced_refusal.expected_refusals, 1);
        assert_eq!(forced_refusal.expected_refusals_matched, 0);
        assert_eq!(forced_refusal.unexpected_mismatches, 1);
        assert_eq!(
            forced_refusal.first_failed_cell.as_deref(),
            Some("budget.profile.smoke.wall_time_ns.one-over")
        );
        for (total, expected_positive, expected_refusals, classifier) in [
            (
                284,
                142,
                142,
                super::limit_matrix_partition as fn(u32) -> super::BaseE2eMixedPartitionV1,
            ),
            (
                44,
                36,
                8,
                super::budget_matrix_partition as fn(u32) -> super::BaseE2eMixedPartitionV1,
            ),
            (
                105,
                96,
                9,
                super::identity_matrix_partition as fn(u32) -> super::BaseE2eMixedPartitionV1,
            ),
            (
                5,
                2,
                3,
                super::no_claim_matrix_partition as fn(u32) -> super::BaseE2eMixedPartitionV1,
            ),
        ] {
            let progress = super::mixed_progress_from_first_failure(Ok(total), total, classifier);
            assert!(progress.is_green());
            assert_eq!(progress.positive_eligible, expected_positive);
            assert_eq!(progress.expected_refusals, expected_refusals);
        }
        let state_partitions = super::state_matrix_partitions();
        let state_progress =
            super::mixed_progress_from_first_failure(Ok(32_460), 32_460, |ordinal| {
                state_partitions
                    [usize::try_from(ordinal - 1).expect("state matrix ordinal fits usize")]
            });
        assert!(state_progress.is_green());
        assert_eq!(state_progress.positive_eligible, 69);
        assert_eq!(state_progress.expected_refusals, 32_391);
        assert!(first.positive_eligible() + first.expected_refusals() > 0);
        assert!(first.positive_matched() + first.expected_refusals_matched() > 0);
        assert!(first.unsupported() > 0);
        assert!(!first.log().events().is_empty());
        assert_eq!(first.manifest_root(), projection.manifest_root());
        assert_ne!(first.execution_root(), first.manifest_root());
        assert_eq!(first.journey_executions().len(), 5);
        assert!(first.retained_artifact_claim().is_absent());
        assert!(first.retained_artifact().is_none());
        assert!(
            first
                .log()
                .events()
                .iter()
                .all(|event| event.relative_artifact().is_none())
        );
        assert_eq!(
            first
                .journey_executions()
                .iter()
                .map(super::BaseE2eJourneyExecutionReportV1::journey)
                .collect::<Vec<_>>(),
            BaseE2eJourneyV1::ALL
        );
        let journey_manifest_roots = first
            .journey_executions()
            .iter()
            .map(super::BaseE2eJourneyExecutionReportV1::manifest_root)
            .collect::<std::collections::BTreeSet<_>>();
        let journey_execution_roots = first
            .journey_executions()
            .iter()
            .map(super::BaseE2eJourneyExecutionReportV1::execution_root)
            .collect::<std::collections::BTreeSet<_>>();
        let row_execution_witnesses = first
            .journey_executions()
            .iter()
            .flat_map(super::BaseE2eJourneyExecutionReportV1::results)
            .map(|result| {
                result
                    .execution_witness_root()
                    .expect("aggregate retains only witness-bound execution rows")
            })
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(journey_manifest_roots.len(), BaseE2eJourneyV1::ALL.len());
        assert_eq!(journey_execution_roots.len(), BaseE2eJourneyV1::ALL.len());
        assert_eq!(
            row_execution_witnesses.len(),
            usize::try_from(first.projection_rows_checked()).expect("bounded row count")
        );
        assert!(journey_manifest_roots.is_disjoint(&journey_execution_roots));
        assert!(!journey_manifest_roots.contains(&first.execution_root()));
        assert!(!journey_execution_roots.contains(&first.execution_root()));
        assert!(!journey_manifest_roots.contains(&first.manifest_root()));
        assert!(!journey_execution_roots.contains(&first.manifest_root()));

        let exact_executions = first.journey_executions().to_vec();
        assert_eq!(
            super::projection_execution_root(
                &projection,
                &harness(),
                &exact_executions,
                first.retained_artifact_claim(),
            )
            .expect("exact aggregate execution root"),
            first.execution_root()
        );
        let missing = super::projection_execution_root(
            &projection,
            &harness(),
            &exact_executions[..exact_executions.len() - 1],
            first.retained_artifact_claim(),
        )
        .expect_err("missing journey execution must refuse");
        assert_eq!(missing.kind(), ConstructionErrorKindV2::Missing);
        let mut extra = exact_executions.clone();
        extra.push(exact_executions[0].clone());
        let extra_error = super::projection_execution_root(
            &projection,
            &harness(),
            &extra,
            first.retained_artifact_claim(),
        )
        .expect_err("extra journey execution must refuse");
        assert_eq!(extra_error.kind(), ConstructionErrorKindV2::Unexpected);
        let mut reordered = exact_executions.clone();
        reordered.swap(0, 1);
        let reorder_error = super::projection_execution_root(
            &projection,
            &harness(),
            &reordered,
            first.retained_artifact_claim(),
        )
        .expect_err("reordered journey executions must refuse");
        assert_eq!(reorder_error.kind(), ConstructionErrorKindV2::OutOfOrder);
        let context_error = super::projection_execution_root(
            &projection,
            &harness_variant(1, 0, 0, "aarch64-apple-darwin", &["default"], 0),
            &exact_executions,
            first.retained_artifact_claim(),
        )
        .expect_err("cross-context journey executions must refuse");
        assert_eq!(context_error.kind(), ConstructionErrorKindV2::Incompatible);
        let mut substituted_root = exact_executions.clone();
        substituted_root[0].root = projection.manifest_root();
        let substitution_error = super::projection_execution_root(
            &projection,
            &harness(),
            &substituted_root,
            first.retained_artifact_claim(),
        )
        .expect_err("substituted journey execution root must refuse");
        assert_eq!(
            substitution_error.field(),
            "base_e2e_projection_execution.journey_execution_root"
        );
        let comparison_journey = exact_executions[0].journey();
        let comparison_presented = presented_results(&projection, comparison_journey, &harness());
        let comparison = compare_base_e2e_journey_results_v1(
            &projection,
            comparison_journey,
            &comparison_presented,
        )
        .expect("exact public comparison");
        let mut comparison_substitution = exact_executions.clone();
        comparison_substitution[0].results = comparison.results.clone();
        let comparison_substitution_error = super::projection_execution_root(
            &projection,
            &harness(),
            &comparison_substitution,
            first.retained_artifact_claim(),
        )
        .expect_err("comparison rows cannot substitute for witness-bound execution rows");
        assert_eq!(
            comparison_substitution_error.field(),
            "base_e2e_journey_execution.witness_root"
        );

        let storage = super::publication_storage().expect("publication storage projection");
        assert_eq!(storage.artifact, 1);
        assert_eq!(storage.system_publication, 6);
        assert_eq!(storage.publication, 7);
        let publication_event = first
            .log()
            .events()
            .iter()
            .find(|event| {
                event
                    .case()
                    .is_some_and(|case| case.as_str() == "publication-storage")
            })
            .expect("publication-storage terminal event");
        for (name, expected) in [
            ("artifact-stored-bytes", 1_u64),
            ("system-publication-stored-bytes", 6_u64),
            ("publication-stored-bytes", 7_u64),
        ] {
            let value = publication_event
                .fields()
                .iter()
                .find(|field| field.name().as_str() == name)
                .expect("storage field");
            assert_eq!(value.value(), &crate::value::TypedValueV2::U64(expected));
        }
        assert_eq!(
            publication_event
                .fields()
                .iter()
                .find(|field| field.name().as_str() == "stored-byte-unit")
                .expect("stored byte unit")
                .value(),
            &crate::value::TypedValueV2::Token(
                StableTokenV2::new("stored-bytes").expect("stored byte unit token")
            )
        );

        let mut journey_summary_index = 0_usize;
        for event in first.log().events() {
            let manifest = event
                .fields()
                .iter()
                .find(|field| field.name().as_str() == "manifest-root")
                .expect("manifest root on every event")
                .value();
            let legacy = event
                .fields()
                .iter()
                .find(|field| field.name().as_str() == "projection-root")
                .expect("legacy projection root on every event")
                .value();
            assert_eq!(manifest, legacy);
            let execution = event
                .fields()
                .iter()
                .find(|field| field.name().as_str() == "execution-root");
            match event.kind() {
                crate::logging::BaseE2eLogKindV1::JourneySummary => {
                    let execution = execution.expect("journey summary execution root");
                    assert_ne!(execution.value(), manifest);
                    assert_eq!(
                        execution.value(),
                        &super::opaque_root(
                            first.journey_executions()[journey_summary_index].execution_root()
                        )
                        .expect("typed journey execution root")
                    );
                    journey_summary_index += 1;
                }
                crate::logging::BaseE2eLogKindV1::ProjectionSummary => {
                    let execution = execution.expect("projection summary execution root");
                    assert_ne!(execution.value(), manifest);
                    assert_eq!(
                        execution.value(),
                        &super::opaque_root(first.execution_root())
                            .expect("typed aggregate execution root")
                    );
                }
                crate::logging::BaseE2eLogKindV1::JourneyStart
                | crate::logging::BaseE2eLogKindV1::CaseTerminal => {
                    assert!(execution.is_none());
                }
            }
        }
        assert_eq!(
            journey_summary_index,
            first.journey_executions().len(),
            "every retained journey execution has one ordered summary root"
        );
    }

    #[test]
    fn one_field_projection_mutation_moves_the_root() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let journey = &projection.journeys()[0];
        let row = &journey.rows()[0];
        let exact = journey_row_root(
            row.journey(),
            row.downstream_owner(),
            row.downstream_driver(),
            row.downstream_script(),
            row.downstream_manifest(),
            row.consumption_rationale(),
            row.fixture_reference(),
            row.source_closure_root(),
            row.log_schema_root(),
            row.semantic_manifest_root(),
        )
        .expect("exact mapping root");
        assert_eq!(exact, row.mapping_root());

        let other_owner = StableTokenV2::new("other-owner").expect("owner");
        let other_driver = StableTokenV2::new("other-driver").expect("driver");
        let other_script = LogicalBundlePathV1::new("scripts/ci/other.sh").expect("script");
        let other_manifest =
            LogicalBundlePathV1::new("scripts/ci/manifests/other.v1.json").expect("manifest");
        let other_rationale = StableTokenV2::new("other-rationale").expect("rationale");
        let other_fixture = StableTokenV2::new("other-fixture").expect("fixture");
        let other_root = hash_domain("projection-test-field-mutation.v1", b"other");
        let mutations = [
            journey_row_root(
                row.journey(),
                &other_owner,
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("owner mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                &other_driver,
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("driver mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                &other_script,
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("script mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                &other_manifest,
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("manifest mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                &other_rationale,
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("rationale mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                &other_fixture,
                row.source_closure_root(),
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("fixture mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                other_root,
                row.log_schema_root(),
                row.semantic_manifest_root(),
            )
            .expect("source mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                other_root,
                row.semantic_manifest_root(),
            )
            .expect("logging mutation"),
            journey_row_root(
                row.journey(),
                row.downstream_owner(),
                row.downstream_driver(),
                row.downstream_script(),
                row.downstream_manifest(),
                row.consumption_rationale(),
                row.fixture_reference(),
                row.source_closure_root(),
                row.log_schema_root(),
                other_root,
            )
            .expect("semantic mutation"),
        ];
        assert!(mutations.into_iter().all(|root| root != exact));
    }

    #[test]
    fn every_harness_context_field_moves_exactly_one_context_root() {
        let exact = harness_variant(0, 0, 0, "aarch64-apple-darwin", &["default"], 0);
        let roots = [
            harness_variant(1, 0, 0, "aarch64-apple-darwin", &["default"], 0).context_root(),
            harness_variant(0, 1, 0, "aarch64-apple-darwin", &["default"], 0).context_root(),
            harness_variant(0, 0, 1, "aarch64-apple-darwin", &["default"], 0).context_root(),
            harness_variant(0, 0, 0, "x86_64-unknown-linux-gnu", &["default"], 0).context_root(),
            harness_variant(0, 0, 0, "aarch64-apple-darwin", &["default", "strict"], 0)
                .context_root(),
            harness_variant(0, 0, 0, "aarch64-apple-darwin", &["default"], 1).context_root(),
        ];
        assert!(roots.into_iter().all(|root| root != exact.context_root()));
        assert_eq!(
            roots
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            roots.len()
        );
    }

    #[test]
    fn source_closure_membership_and_order_are_exact_and_content_bound() {
        let closure = RunnerV2BaseSourceClosureV1::frozen().expect("frozen source closure");
        assert_eq!(closure.entries().len(), EXPECTED_SOURCE_PATHS_V1.len());
        assert_eq!(
            closure.dependency_declaration_root(),
            current_direct_dependency_declaration_root_v1()
        );
        assert_eq!(
            closure.compatible_snapshot().root().content_hash(),
            closure.snapshot_root(),
            "the typed compatibility witness must be the existing snapshot hash, not a wrapper hash"
        );
        let mut identities = std::collections::BTreeSet::new();
        let mut entry_roots = std::collections::BTreeSet::new();
        for ((entry, expected_path), embedded) in closure
            .entries()
            .iter()
            .zip(EXPECTED_SOURCE_PATHS_V1)
            .zip(EMBEDDED_SOURCE_FILES_V1)
        {
            assert_eq!(entry.path(), expected_path);
            assert_eq!(embedded.path, expected_path);
            assert_eq!(
                entry.encoded_bytes(),
                u64::try_from(embedded.bytes.len()).expect("source length fits u64")
            );
            assert_eq!(
                entry.content_root(),
                hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, embedded.bytes)
            );
            assert_eq!(entry.owner(), embedded.owner);
            assert_eq!(entry.source_route(), embedded.source_route);
            assert_eq!(
                entry.expected_source_identity(),
                embedded.expected_source_identity
            );
            assert_eq!(
                entry.expected_source_identity_root(),
                expected_source_identity_root(&embedded)
            );
            assert_eq!(
                entry.snapshot_policy(),
                BaseSourceSnapshotPolicyV1::ExactCommonCompiledSnapshot
            );
            assert_eq!(entry.snapshot_root(), closure.snapshot_root());
            assert!(identities.insert(entry.expected_source_identity()));
            assert!(entry_roots.insert(entry.entry_root()));
        }
        assert_eq!(
            closure
                .entries()
                .iter()
                .filter(|entry| {
                    entry.owner() == BaseSourceOwnerV1::FrankensimWorkspaceGovernance
                })
                .count(),
            5
        );
        assert_eq!(
            closure
                .entries()
                .iter()
                .filter(|entry| entry.owner() == BaseSourceOwnerV1::RunnerV2BaseSchema)
                .count(),
            22
        );
        assert_eq!(
            closure
                .entries()
                .iter()
                .filter(|entry| entry.source_route() == BaseSourceRouteV1::CrateModule)
                .count(),
            20
        );
        for route in [
            BaseSourceRouteV1::WorkspaceCargoConfig,
            BaseSourceRouteV1::WorkspaceLockfile,
            BaseSourceRouteV1::WorkspaceManifest,
            BaseSourceRouteV1::WorkspaceConstellationLock,
            BaseSourceRouteV1::CrateContract,
            BaseSourceRouteV1::CrateManifest,
            BaseSourceRouteV1::WorkspaceToolchain,
        ] {
            assert_eq!(
                closure
                    .entries()
                    .iter()
                    .filter(|entry| entry.source_route() == route)
                    .count(),
                1
            );
        }
        assert!(
            EXPECTED_SOURCE_PATHS_V1
                .windows(2)
                .all(|pair| pair[0].as_bytes() < pair[1].as_bytes())
        );

        let sentinel = "caller-controlled-no-echo-sentinel";
        let error = BaseSourceClosureInputV1::exact_presented(sentinel, b"rejected".to_vec())
            .expect_err("an unregistered caller-presented source path must refuse");
        assert_eq!(error.observed(), "<redacted:caller-controlled-text>");
        let mut encoded_error = Vec::new();
        encode_construction_error(&mut encoded_error, &error);
        assert!(
            !encoded_error
                .windows(sentinel.len())
                .any(|window| window == sentinel.as_bytes()),
            "canonical error detail must not retain rejected sentinel bytes"
        );
        assert!(
            encoded_error
                .windows(b"<redacted:caller-controlled-text>".len())
                .any(|window| window == b"<redacted:caller-controlled-text>")
        );
    }

    #[test]
    fn exact_source_reconstruction_is_deterministic() {
        let frozen = RunnerV2BaseSourceClosureV1::frozen().expect("frozen source closure");
        let reconstructed = RunnerV2BaseSourceClosureV1::reconstruct(&frozen_source_inputs())
            .expect("exact reconstruction");
        assert_eq!(reconstructed, frozen);
        assert_eq!(reconstructed.root(), frozen.root());
    }

    #[test]
    fn source_reconstruction_rejects_missing_entry() {
        let mut inputs = frozen_source_inputs();
        inputs.pop();
        let error = RunnerV2BaseSourceClosureV1::reconstruct(&inputs)
            .expect_err("missing source must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
    }

    #[test]
    fn source_reconstruction_rejects_extra_entry() {
        let mut inputs = frozen_source_inputs();
        inputs.push(BaseSourceClosureInputV1::presented(
            "crates/fs-evidence-runner/src/ambient.rs",
            b"ambient".to_vec(),
        ));
        let error = RunnerV2BaseSourceClosureV1::reconstruct(&inputs)
            .expect_err("extra source must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
    }

    #[test]
    fn source_reconstruction_rejects_duplicate_entry() {
        let mut inputs = frozen_source_inputs();
        inputs[1] = inputs[0].clone();
        let error = RunnerV2BaseSourceClosureV1::reconstruct(&inputs)
            .expect_err("duplicate source must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Duplicate);
    }

    #[test]
    fn source_reconstruction_rejects_reordered_entries() {
        let mut inputs = frozen_source_inputs();
        inputs.swap(0, 1);
        let error = RunnerV2BaseSourceClosureV1::reconstruct(&inputs)
            .expect_err("reordered source must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::OutOfOrder);
    }

    #[test]
    fn source_reconstruction_rejects_mutated_bytes() {
        let mut inputs = frozen_source_inputs();
        inputs[0].bytes[0] ^= 1;
        let error = RunnerV2BaseSourceClosureV1::reconstruct(&inputs)
            .expect_err("mutated source must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
    }

    type SourceInputMutationV1 = (&'static str, fn(&mut BaseSourceClosureInputV1));

    #[test]
    fn source_reconstruction_rejects_owner_route_identity_and_policy_mutations() {
        let exact = frozen_source_inputs();
        let mutations: [SourceInputMutationV1; 4] = [
            (
                "base_source_closure.owner",
                |input: &mut BaseSourceClosureInputV1| input.owner_code += 1,
            ),
            (
                "base_source_closure.source_route",
                |input: &mut BaseSourceClosureInputV1| input.source_route_code += 1,
            ),
            (
                "base_source_closure.expected_source_identity_root",
                |input: &mut BaseSourceClosureInputV1| {
                    input.expected_source_identity_root =
                        hash_domain("source-test-wrong-identity.v1", b"wrong");
                },
            ),
            (
                "base_source_closure.snapshot_policy",
                |input: &mut BaseSourceClosureInputV1| input.snapshot_policy_code += 1,
            ),
        ];
        for (field, mutate) in mutations {
            let mut mutant = exact.clone();
            mutate(&mut mutant[0]);
            let error = RunnerV2BaseSourceClosureV1::reconstruct(&mutant)
                .expect_err("metadata mutation must refuse");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), field);
        }
    }

    #[test]
    fn source_reconstruction_rejects_length_content_and_snapshot_mutations() {
        let exact = frozen_source_inputs();

        let mut wrong_length = exact.clone();
        wrong_length[0].encoded_bytes += 1;
        assert_eq!(
            RunnerV2BaseSourceClosureV1::reconstruct(&wrong_length)
                .expect_err("length mutation")
                .field(),
            "base_source_closure.encoded_bytes"
        );

        let mut wrong_content = exact.clone();
        wrong_content[0].content_root = hash_domain("source-test-wrong-content.v1", b"wrong");
        assert_eq!(
            RunnerV2BaseSourceClosureV1::reconstruct(&wrong_content)
                .expect_err("content-root mutation")
                .field(),
            "base_source_closure.content_root"
        );

        let wrong_snapshot = hash_domain("source-test-wrong-snapshot.v1", b"wrong");
        let mut mixed_snapshot = exact.clone();
        mixed_snapshot[0].snapshot_root = wrong_snapshot;
        assert_eq!(
            RunnerV2BaseSourceClosureV1::reconstruct(&mixed_snapshot)
                .expect_err("mixed snapshot")
                .field(),
            "base_source_closure.snapshot_root"
        );

        let mut common_wrong_snapshot = exact;
        for input in &mut common_wrong_snapshot {
            input.snapshot_root = wrong_snapshot;
        }
        assert_eq!(
            RunnerV2BaseSourceClosureV1::reconstruct(&common_wrong_snapshot)
                .expect_err("common but wrong snapshot")
                .field(),
            "base_source_closure.snapshot_root"
        );
    }

    #[test]
    fn source_reconstruction_rejects_resealed_bytes_and_dependency_identity() {
        let mut resealed = frozen_source_inputs();
        resealed[0].bytes[0] ^= 1;
        resealed[0].encoded_bytes = u64::try_from(resealed[0].bytes.len()).expect("fixture length");
        resealed[0].content_root =
            hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, &resealed[0].bytes);
        let error =
            RunnerV2BaseSourceClosureV1::reconstruct(&resealed).expect_err("resealed stale bytes");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "base_source_closure.bytes");

        let wrong_dependency = hash_domain("source-test-wrong-dependency-declaration.v1", b"wrong");
        let error = RunnerV2BaseSourceClosureV1::reconstruct_with_dependency_declaration(
            &frozen_source_inputs(),
            wrong_dependency,
        )
        .expect_err("dependency declaration mutation");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(
            error.field(),
            "base_source_closure.dependency_declaration_root"
        );
    }

    #[test]
    fn source_entry_root_moves_for_each_admitted_metadata_field() {
        let exact = EMBEDDED_SOURCE_FILES_V1[6];
        let snapshot = RunnerV2BaseSourceClosureV1::frozen()
            .expect("closure")
            .snapshot_root();
        let root = source_closure_entry(&exact, snapshot)
            .expect("exact entry")
            .entry_root();
        let mutations = [
            super::EmbeddedSourceFileV1 {
                path: "crates/fs-evidence-runner/src/alternate.rs",
                ..exact
            },
            super::EmbeddedSourceFileV1 {
                owner: BaseSourceOwnerV1::FrankensimWorkspaceGovernance,
                ..exact
            },
            super::EmbeddedSourceFileV1 {
                source_route: BaseSourceRouteV1::CrateContract,
                ..exact
            },
            super::EmbeddedSourceFileV1 {
                expected_source_identity: "frankensim.fs-evidence-runner.src.alternate.v1",
                ..exact
            },
            super::EmbeddedSourceFileV1 {
                bytes: b"alternate source bytes",
                ..exact
            },
        ];
        assert!(mutations.into_iter().all(|mutation| {
            source_closure_entry(&mutation, snapshot)
                .expect("mutant entry")
                .entry_root()
                != root
        }));
        assert_ne!(
            source_closure_entry(
                &exact,
                hash_domain("source-test-alternate-snapshot.v1", b"alternate")
            )
            .expect("snapshot mutant")
            .entry_root(),
            root
        );
    }

    #[test]
    fn ac38_coverage_manifest_source_of_truth_and_checked_report_are_exact() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let expected_counts = [
            (BaseCoverageManifestClassV1::Unit, 12),
            (BaseCoverageManifestClassV1::CompileFailDoctest, 78),
            (BaseCoverageManifestClassV1::ManifestContract, 29),
            (BaseCoverageManifestClassV1::ProjectionE2e, 98),
            (BaseCoverageManifestClassV1::RuntimeLogging, 1),
            (BaseCoverageManifestClassV1::SourceClosure, 15),
            (BaseCoverageManifestClassV1::ExternalE2eScript, 5),
            (BaseCoverageManifestClassV1::ExternalMutation, 1),
            (BaseCoverageManifestClassV1::ExternalGovernance, 1),
            (BaseCoverageManifestClassV1::Boundary, 53),
            (BaseCoverageManifestClassV1::PropertyMetamorphic, 24),
            (BaseCoverageManifestClassV1::SchemaDescriptor, 50),
            (BaseCoverageManifestClassV1::Mutation, 61),
            (BaseCoverageManifestClassV1::NoMockIntegration, 17),
        ];
        let manifest = projection.coverage_manifest();
        assert_eq!(manifest.cases().len(), 445);
        for (class, expected_count) in expected_counts {
            assert_eq!(manifest.case_count(class), expected_count);
            assert!(manifest.case_count(class) > 0);
            let ordinals = manifest
                .cases()
                .iter()
                .filter(|source_case| source_case.class() == class)
                .map(crate::coverage::BaseCoverageManifestCaseV1::ordinal)
                .collect::<Vec<_>>();
            assert!(ordinals.windows(2).all(|pair| pair[0] < pair[1]));
        }
        let close_manifest = projection.close_manifest();
        assert_eq!(close_manifest.source_manifest_root(), manifest.root());
        assert_eq!(close_manifest.cells().len(), manifest.cases().len());
        assert_eq!(close_manifest.cells().len(), 445);
        for group in crate::coverage::BaseCoverageCloseGroupV1::ALL {
            assert!(
                close_manifest.group_count(group) > 0,
                "group {} must have an exact source cell",
                group.stable_name()
            );
        }
        for facet in crate::coverage::BaseCoverageCloseFacetV1::ALL {
            assert!(
                close_manifest.facet_count(facet) > 0,
                "facet {} must be applicable or exactly reasoned",
                facet.stable_name()
            );
        }

        let report =
            run_base_e2e_projection_v1(&projection, &harness()).expect("aggregate coverage report");
        assert!(report.coverage_report().is_green());
        assert_eq!(report.coverage_report().manifest_root(), manifest.root());
        assert_eq!(report.coverage_report().results().len(), 114);
        let exact_observations = report
            .coverage_report()
            .results()
            .iter()
            .map(|result| {
                (
                    result.source_case_id().to_owned(),
                    (result.outcome(), result.evidence_root()),
                )
            })
            .collect::<std::collections::BTreeMap<_, _>>();
        assert_eq!(
            &reconstruct_exact_local_coverage_report(manifest, &exact_observations)
                .expect("the independently declared local selection must reconstruct"),
            report.coverage_report()
        );

        let mut missing = exact_observations.clone();
        let missing_id = report.coverage_report().results()[0].source_case_id();
        assert!(missing.remove(missing_id).is_some());
        let error = reconstruct_exact_local_coverage_report(manifest, &missing)
            .expect_err("one omitted local observation must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
        assert_eq!(error.field(), "coverage.runtime_result.source_case_id");

        let mut extra = exact_observations;
        let extra_id = manifest
            .cases()
            .iter()
            .find(|case| case.class() == BaseCoverageManifestClassV1::Unit)
            .expect("the frozen manifest has unit cases")
            .id()
            .to_owned();
        assert!(
            extra
                .insert(
                    extra_id,
                    (
                        crate::coverage::BaseCoveragePresentedOutcomeV1::PositiveMatched,
                        hash_domain(
                            "org.frankensim.fs-evidence-runner.coverage-extra-test.v1",
                            b"extra",
                        ),
                    ),
                )
                .is_none()
        );
        let error = reconstruct_exact_local_coverage_report(manifest, &extra)
            .expect_err("one non-local observation must refuse");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
        assert_eq!(error.field(), "coverage.runtime_result.source_case_id");
        assert_eq!(report.projection_rows_checked(), 98);
        assert_eq!(report.projection_e2e_checked(), 134_455);
        assert_eq!(report.source_closure_positive_eligible(), 1);
        assert_eq!(report.source_closure_expected_refusals(), 14);
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "this exhaustive boundary test keeps comparison, private execution-witness, and typed-detail mutations adjacent"
    )]
    fn exact_result_join_reconstructs_the_frozen_journey() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let harness = harness();
        let journey = BaseE2eJourneyV1::CanonicalRunnerV2;
        let presented = presented_results(&projection, journey, &harness);
        let first = compare_base_e2e_journey_results_v1(&projection, journey, &presented)
            .expect("exact comparison");
        let second = compare_base_e2e_journey_results_v1(&projection, journey, &presented)
            .expect("deterministic exact comparison");
        assert_eq!(first, second);
        assert!(first.exact_match());
        assert_eq!(first.results().len(), 24);
        assert_eq!(first.unexpected_mismatches(), 0);
        assert_eq!(first.positive_eligible(), first.positive_matched());
        assert_eq!(first.expected_refusals(), first.expected_refusals_matched());
        let manifest = projection
            .journeys()
            .iter()
            .find(|candidate| candidate.journey() == journey)
            .expect("canonical journey");
        for (row, result) in manifest.rows().iter().zip(first.results()) {
            assert_eq!(
                result.expected_detail_manifest(),
                row.expected_detail_manifest()
            );
            assert_eq!(
                result.observed_detail_manifest(),
                row.expected_detail_manifest()
            );
            assert_eq!(
                result.detail_cells_matched(),
                row.expected_detail_cell_count()
            );
            assert!(result.execution_witness_root().is_none());
            assert_eq!(
                result
                    .observed_detail_cells()
                    .expect("exact typed comparisons borrow the cached detail cells"),
                row.expected_detail_cells()
            );
            assert!(result.first_observed_detail_divergence().is_none());
            assert!(result.first_divergence_root().is_none());
        }

        let execution_report = run_base_e2e_journey_v1(&projection, journey, &harness)
            .expect("source-closed in-process execution");
        assert_eq!(execution_report.unexpected_mismatches(), 0);
        assert_ne!(execution_report.execution_root(), first.comparison_root());
        assert_eq!(execution_report.results().len(), manifest.rows().len());
        let execution_witnesses = execution_report
            .results()
            .iter()
            .map(|result| {
                result
                    .execution_witness_root()
                    .expect("every executed row retains its exact witness")
            })
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(execution_witnesses.len(), manifest.rows().len());
        let second_journey = BaseE2eJourneyV1::PublicationState;
        let second_execution = run_base_e2e_journey_v1(&projection, second_journey, &harness)
            .expect("second source-closed journey execution");
        let second_witnesses = second_execution
            .results()
            .iter()
            .map(|result| {
                result
                    .execution_witness_root()
                    .expect("every second-journey row retains a witness")
            })
            .collect::<std::collections::BTreeSet<_>>();
        assert!(
            execution_witnesses.is_disjoint(&second_witnesses),
            "journey, manifest, mapping, and ordinal binding keep witnesses distinct"
        );

        let detail_index = presented
            .iter()
            .position(|result| result.observed_detail_cell_count() > 2)
            .expect("canonical journey contains a multi-cell refusal detail");
        let exact = &presented[detail_index];
        let detail_row = &manifest.rows()[detail_index];
        let legacy_descriptor_only = BaseE2ePresentedRowResultV1::new_with_detail_descriptor(
            journey,
            exact.row_id().clone(),
            exact.semantic_manifest_root(),
            exact.observed(),
            exact.counts(),
            exact.observed_detail_manifest(),
            exact.detail_cells_matched(),
            exact.first_unexpected_cell().cloned(),
        )
        .expect("the descriptor-only compatibility observation remains constructible");
        assert!(!legacy_descriptor_only.typed_detail_cells_presented());
        let mut legacy_descriptor_results = presented.clone();
        legacy_descriptor_results[detail_index] = legacy_descriptor_only;
        let legacy_descriptor_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &legacy_descriptor_results)
                .expect("descriptor-only compatibility input remains inspectable as a mismatch");
        assert_eq!(legacy_descriptor_report.unexpected_mismatches(), 1);
        let legacy_descriptor_result = &legacy_descriptor_report.results()[detail_index];
        assert!(!legacy_descriptor_result.matched());
        assert!(
            legacy_descriptor_result.observed_detail_cells().is_none(),
            "an opaque descriptor must never borrow the independent expected cells as observed"
        );
        assert!(
            legacy_descriptor_result
                .first_observed_detail_divergence()
                .is_some(),
            "the unverified descriptor boundary must remain explicit"
        );
        let checked_exact = BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
            detail_row,
            exact.observed(),
            exact.counts(),
            detail_row.expected_detail_manifest(),
            detail_row.expected_detail_cells(),
            None,
        )
        .expect("the public typed observation path accepts the exact closed manifest");
        let mut checked_exact_presented = presented.clone();
        checked_exact_presented[detail_index] = checked_exact;
        let checked_exact_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &checked_exact_presented)
                .expect("the public typed comparison reports exact equality");
        assert!(checked_exact_report.exact_match());
        assert_eq!(checked_exact_report.unexpected_mismatches(), 0);
        let copied_oracle_result = &checked_exact_report.results()[detail_index];
        assert!(copied_oracle_result.matched());
        assert!(copied_oracle_result.execution_witness_root().is_none());
        assert_eq!(
            copied_oracle_result
                .observed_detail_cells()
                .expect("typed copied cells remain inspectable comparison data"),
            detail_row.expected_detail_cells()
        );
        assert!(copied_oracle_result.first_unexpected_cell().is_none());

        let positive_only_index = manifest
            .rows()
            .iter()
            .position(|row| {
                row.expected_detail_cell_count() == 0
                    && row.positive_cell_count() > 0
                    && row.expected() == super::BaseE2eExpectedDecisionV1::Accept
            })
            .expect("canonical journey contains a positive-only empty-detail row");
        let positive_only_row = &manifest.rows()[positive_only_index];
        let positive_only_exact = &presented[positive_only_index];
        let copied_empty_detail = BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
            positive_only_row,
            positive_only_exact.observed(),
            positive_only_exact.counts(),
            positive_only_row.expected_detail_manifest(),
            &[],
            None,
        )
        .expect("empty typed comparison remains structurally inspectable");
        let mut copied_empty_results = presented.clone();
        copied_empty_results[positive_only_index] = copied_empty_detail;
        let copied_empty_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &copied_empty_results)
                .expect("positive-only empty-detail comparison reports equality");
        assert!(copied_empty_report.exact_match());
        assert_eq!(copied_empty_report.unexpected_mismatches(), 0);
        let copied_empty_result = &copied_empty_report.results()[positive_only_index];
        assert!(copied_empty_result.matched());
        assert!(copied_empty_result.execution_witness_root().is_none());
        assert_eq!(
            copied_empty_result.observed_detail_cells(),
            Some(positive_only_row.expected_detail_cells())
        );
        assert!(copied_empty_result.first_unexpected_cell().is_none());
        let executed_positive = &execution_report.results()[positive_only_index];
        assert!(executed_positive.matched());
        assert!(executed_positive.execution_witness_root().is_some());
        assert_ne!(executed_positive.root(), copied_empty_result.root());

        let executed_rows = manifest
            .rows()
            .iter()
            .enumerate()
            .map(|(index, row)| {
                super::finalize_executed_row(
                    manifest,
                    u32::try_from(index + 1).expect("bounded row ordinal"),
                    row,
                    &harness,
                    execute_case(row.kind(), &harness),
                )
            })
            .collect::<Result<Vec<_>, _>>()
            .expect("one private executed row per exact manifest row");
        assert_eq!(
            super::finalize_journey_execution(manifest, &harness, executed_rows.clone())
                .expect("the private witness finalizer reconstructs the execution report"),
            execution_report
        );
        let mut wrong_ordinal = executed_rows.clone();
        wrong_ordinal[0].row_ordinal = 2;
        assert_eq!(
            super::finalize_journey_execution(manifest, &harness, wrong_ordinal)
                .expect_err("a repeated or shifted row ordinal must refuse")
                .field(),
            "base_e2e_journey_execution.row_ordinal"
        );
        let mut wrong_witness = executed_rows.clone();
        wrong_witness[0].witness_root =
            hash_domain("projection-test-wrong-execution-witness.v1", b"wrong");
        assert_eq!(
            super::finalize_journey_execution(manifest, &harness, wrong_witness)
                .expect_err("a substituted private witness must refuse")
                .field(),
            "base_e2e_journey_execution.witness_root"
        );
        assert_eq!(
            super::finalize_journey_execution(
                manifest,
                &harness_variant(1, 0, 0, "aarch64-apple-darwin", &["default"], 0),
                executed_rows.clone(),
            )
            .expect_err("a witness cannot move to a different harness context")
            .field(),
            "base_e2e_journey_execution.witness_root"
        );
        assert_eq!(
            super::finalize_executed_row(
                manifest,
                2,
                &manifest.rows()[0],
                &harness,
                execute_case(manifest.rows()[0].kind(), &harness),
            )
            .expect_err("an execution cannot bind the right row to a wrong ordinal")
            .field(),
            "base_e2e_execution_witness.row_binding"
        );
        let mut wrong_mapping_row = manifest.rows()[0].clone();
        wrong_mapping_row.mapping_root =
            hash_domain("projection-test-wrong-row-mapping.v1", b"wrong");
        assert_eq!(
            super::finalize_executed_row(
                manifest,
                1,
                &wrong_mapping_row,
                &harness,
                execute_case(wrong_mapping_row.kind(), &harness),
            )
            .expect_err("a substituted downstream mapping root must refuse")
            .field(),
            "base_e2e_execution_witness.row_binding"
        );

        let exact_counts = exact.counts();
        let row_contract_counts = BaseE2eObservedCountsV1::new(
            BaseE2eMatchedPartitionV1::new(
                exact_counts.positive().eligible(),
                exact_counts.positive().matched() - 1,
            )
            .expect("one positive row-contract mismatch"),
            exact_counts.expected_refusals(),
            exact_counts.unsupported(),
            1,
        )
        .expect("exact row-contract red partition");
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                row_contract_counts,
                detail_row.expected_detail_manifest(),
                detail_row.expected_detail_cells(),
                None,
            )
            .expect_err("a red row contract requires its closed sentinel")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                row_contract_counts,
                detail_row.expected_detail_manifest(),
                detail_row.expected_detail_cells(),
                Some(StableTokenV2::new("arbitrary.row-gap").expect("valid wrong ID")),
            )
            .expect_err("an arbitrary red row-contract ID must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let row_contract_presented = BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
            detail_row,
            exact.observed(),
            row_contract_counts,
            detail_row.expected_detail_manifest(),
            detail_row.expected_detail_cells(),
            Some(
                StableTokenV2::new(super::BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1)
                    .expect("closed row-contract ID"),
            ),
        )
        .expect("the exact row.contract sentinel is admitted");
        let row_contract_observation = super::observation_from_presented(&row_contract_presented);
        let expected_row_contract_root = super::row_contract_divergence_root(
            detail_row,
            &row_contract_observation,
            super::BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1,
        )
        .expect("closed row-contract divergence root");
        let mut row_contract_results = presented.clone();
        row_contract_results[detail_index] = row_contract_presented;
        let row_contract_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &row_contract_results)
                .expect("red row contract compares as a checked mismatch");
        let row_contract_result = &row_contract_report.results()[detail_index];
        assert_eq!(
            row_contract_result.first_unexpected_cell(),
            Some(super::BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1)
        );
        assert!(
            row_contract_result
                .first_observed_detail_divergence()
                .is_none()
        );
        assert_eq!(
            row_contract_result.first_divergence_root(),
            Some(expected_row_contract_root)
        );

        assert_eq!(
            super::BaseE2eDetailCellV1::new(
                detail_row.kind(),
                0,
                &StableTokenV2::new("detail.zero-ordinal").expect("valid token"),
                super::BaseE2eExpectedDecisionV1::Accept,
                super::BaseE2eDetailPayloadV1::AcceptedInstead,
            )
            .expect_err("zero detail ordinal must refuse")
            .kind(),
            ConstructionErrorKindV2::Zero
        );
        let mut missing_cells = detail_row.expected_detail_cells().to_vec();
        missing_cells.pop();
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                exact.counts(),
                detail_row.expected_detail_manifest(),
                &missing_cells,
                None,
            )
            .expect_err("a declared observed cell cannot be missing")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut reordered_cells = detail_row.expected_detail_cells().to_vec();
        reordered_cells.swap(0, 1);
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                exact.counts(),
                detail_row.expected_detail_manifest(),
                &reordered_cells,
                None,
            )
            .expect_err("typed observed detail cells must remain in canonical order")
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let budget_row = manifest
            .rows()
            .iter()
            .find(|row| row.kind() == BaseE2eCaseKindV1::BudgetAdmission)
            .expect("canonical journey owns budget admission");
        let budget_template = budget_row.expected_detail_cells()[0].clone();
        let budget_template_payload = budget_template.payload().clone();
        let budget_payload = |repair_rank: u8, owner: &'static str, repair_target: &'static str| {
            let super::BaseE2eDetailPayloadV1::Budget {
                kind,
                field,
                unit,
                expected,
                observed,
                repair_kind,
                ..
            } = budget_template_payload.clone()
            else {
                panic!("budget template payload")
            };
            super::BaseE2eDetailPayloadV1::Budget {
                kind,
                field,
                unit,
                expected,
                observed,
                owner,
                repair_rank,
                repair_kind,
                repair_target,
            }
        };
        assert!(
            super::BaseE2eDetailCellV1::new(
                budget_row.kind(),
                budget_template.semantic_ordinal(),
                &StableTokenV2::new("budget.rank-16").expect("valid ID"),
                super::BaseE2eExpectedDecisionV1::Refuse,
                budget_payload(
                    16,
                    "fs-evidence-runner.runner-budgets",
                    "max-resident-bytes"
                ),
            )
            .is_ok(),
            "repair rank 16 is the inclusive public bound"
        );
        assert_eq!(
            super::BaseE2eDetailCellV1::new(
                budget_row.kind(),
                budget_template.semantic_ordinal(),
                &StableTokenV2::new("budget.rank-17").expect("valid ID"),
                super::BaseE2eExpectedDecisionV1::Refuse,
                budget_payload(
                    17,
                    "fs-evidence-runner.runner-budgets",
                    "max-resident-bytes"
                ),
            )
            .expect_err("repair rank above the public bound must refuse")
            .kind(),
            ConstructionErrorKindV2::OutOfRange
        );
        for (owner, target) in [
            ("runner owner", "max-resident-bytes"),
            (
                "fs-evidence-runner.runner-budgets",
                "run-$(unexpected-command)",
            ),
        ] {
            assert_eq!(
                super::BaseE2eDetailCellV1::new(
                    budget_row.kind(),
                    budget_template.semantic_ordinal(),
                    &StableTokenV2::new("budget.invalid-repair-token").expect("valid ID"),
                    super::BaseE2eExpectedDecisionV1::Refuse,
                    budget_payload(1, owner, target),
                )
                .expect_err("repair owner and target must be bounded stable tokens")
                .kind(),
                ConstructionErrorKindV2::Incompatible
            );
        }
        assert_eq!(
            super::BaseE2eDetailCellV1::new(
                budget_row.kind(),
                budget_template.semantic_ordinal(),
                &StableTokenV2::new("budget.unsupported-cross-substitution").expect("valid ID"),
                super::BaseE2eExpectedDecisionV1::Unsupported,
                budget_payload(1, "fs-evidence-runner.runner-budgets", "max-resident-bytes"),
            )
            .expect_err("Unsupported cannot carry a Budget refusal payload")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            super::BaseE2eDetailCellV1::new(
                BaseE2eCaseKindV1::WindowsUnicodeAlias,
                1,
                &StableTokenV2::new("path.unsupported-cross-substitution").expect("valid ID"),
                super::BaseE2eExpectedDecisionV1::Refuse,
                super::BaseE2eDetailPayloadV1::PathAdjudication(
                    super::BaseE2ePathAdjudicationDetailV1::UnsupportedWindowsNonAsciiAlias(
                        "alias".to_owned(),
                    ),
                ),
            )
            .expect_err("Refuse cannot carry the explicit Unsupported adjudication")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mut wrong_root = presented.clone();
        wrong_root[detail_index] = BaseE2ePresentedRowResultV1::new_with_detail_manifest(
            journey,
            exact.row_id().clone(),
            exact.semantic_manifest_root(),
            exact.observed(),
            exact.counts(),
            hash_domain("projection-test-wrong-detail-root.v1", b"wrong"),
            exact.observed_detail_cell_count(),
            exact.detail_cells_matched(),
            exact.first_unexpected_cell().cloned(),
        )
        .expect("intrinsically valid wrong detail root");
        let wrong_root_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &wrong_root)
                .expect("detail disagreement is a checked mismatch");
        assert_eq!(wrong_root_report.unexpected_mismatches(), 1);
        assert!(!wrong_root_report.results()[detail_index].matched());
        assert_eq!(
            wrong_root_report.results()[detail_index].first_unexpected_cell(),
            Some("detail.manifest")
        );
        let wrong_root_result = &wrong_root_report.results()[detail_index];
        assert!(wrong_root_result.observed_detail_cells().is_none());
        let divergence = wrong_root_result
            .first_observed_detail_divergence()
            .expect("bounded mismatch descriptor");
        assert_eq!(
            divergence.expected_manifest_root(),
            manifest.rows()[detail_index].expected_detail_manifest_root()
        );
        assert_eq!(
            divergence.observed_manifest_root(),
            wrong_root[detail_index].observed_detail_manifest_root()
        );
        assert!(divergence.expected_cell().is_none());
        assert!(divergence.observed_cell().is_none());
        assert_eq!(
            wrong_root_result.first_divergence_root(),
            Some(divergence.root())
        );

        let expected_non_prefix = manifest.rows()[detail_index].expected_detail_cells()[2].clone();
        let observed_non_prefix = super::BaseE2eDetailCellV1::new(
            expected_non_prefix.kind(),
            expected_non_prefix.semantic_ordinal(),
            &StableTokenV2::new(expected_non_prefix.stable_id()).expect("stable expected cell ID"),
            super::BaseE2eExpectedDecisionV1::Accept,
            super::BaseE2eDetailPayloadV1::AcceptedInstead,
        )
        .expect("public bounded detail-cell constructor");
        let mut observed_non_prefix_cells = detail_row.expected_detail_cells().to_vec();
        observed_non_prefix_cells[2] = observed_non_prefix.clone();
        let observed_non_prefix_manifest = super::BaseE2eDecisionDetailManifestV1::from_cells(
            detail_row.kind(),
            &observed_non_prefix_cells,
        )
        .expect("closed observed detail manifest");
        let exact_counts = exact.counts();
        let red_counts = BaseE2eObservedCountsV1::new(
            exact_counts.positive(),
            BaseE2eMatchedPartitionV1::new(
                exact_counts.expected_refusals().eligible(),
                exact_counts.expected_refusals().matched() - 1,
            )
            .expect("one expected-refusal mismatch"),
            exact_counts.unsupported(),
            1,
        )
        .expect("exact red partition");
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                red_counts,
                observed_non_prefix_manifest,
                &observed_non_prefix_cells,
                None,
            )
            .expect_err("a typed divergence requires its first stable ID")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                red_counts,
                observed_non_prefix_manifest,
                &observed_non_prefix_cells,
                Some(StableTokenV2::new("detail.wrong-id").expect("valid wrong ID")),
            )
            .expect_err("a caller cannot redirect typed divergence selection")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
                detail_row,
                exact.observed(),
                exact.counts(),
                observed_non_prefix_manifest,
                detail_row.expected_detail_cells(),
                None,
            )
            .expect_err("the descriptor root must equal the supplied exact cell slice")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut non_prefix = presented.clone();
        non_prefix[detail_index] = BaseE2ePresentedRowResultV1::new_with_observed_detail_cells(
            detail_row,
            exact.observed(),
            red_counts,
            observed_non_prefix_manifest,
            &observed_non_prefix_cells,
            Some(
                StableTokenV2::new(observed_non_prefix.stable_id())
                    .expect("exact first divergent cell ID"),
            ),
        )
        .expect("bounded non-prefix detail observation");
        let non_prefix_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &non_prefix)
                .expect("non-prefix detail disagreement is a checked mismatch");
        let non_prefix_divergence = non_prefix_report.results()[detail_index]
            .first_observed_detail_divergence()
            .expect("non-prefix divergence descriptor");
        assert_eq!(
            non_prefix_divergence
                .expected_cell()
                .expect("expected cell at observed ordinal"),
            &expected_non_prefix
        );
        assert_eq!(
            non_prefix_divergence
                .observed_cell()
                .expect("bounded typed observed cell"),
            &observed_non_prefix
        );
        assert_eq!(
            non_prefix_report.results()[detail_index].first_divergence_root(),
            Some(non_prefix_divergence.root())
        );
        let red_result = &non_prefix_report.results()[detail_index];
        let red_terminal = super::case_terminal_log_event(
            0,
            manifest,
            &manifest.rows()[detail_index],
            red_result,
            crate::logging::BaseE2eOutcomeV1::Failed,
            &harness,
        )
        .expect("production CaseTerminal emitter retains a reconciled divergence");
        assert_eq!(
            red_terminal.kind(),
            crate::logging::BaseE2eLogKindV1::CaseTerminal
        );
        assert_eq!(
            red_terminal.outcome(),
            crate::logging::BaseE2eOutcomeV1::Failed
        );
        assert!(red_terminal.fields().iter().any(|field| {
            field.field_code() == Some(crate::logging::BaseE2eLogFieldCodeV1::FirstFailedCell)
        }));
        let logged_divergence = red_terminal
            .fields()
            .iter()
            .find(|field| {
                field.field_code()
                    == Some(crate::logging::BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot)
            })
            .expect("field78 typed divergence root");
        assert_eq!(
            logged_divergence.value(),
            &super::opaque_root(
                red_result
                    .first_divergence_root()
                    .expect("red row has a typed divergence root")
            )
            .expect("typed opaque divergence root")
        );

        let mut wrong_count = presented.clone();
        wrong_count[detail_index] = BaseE2ePresentedRowResultV1::new_with_detail_manifest(
            journey,
            exact.row_id().clone(),
            exact.semantic_manifest_root(),
            exact.observed(),
            exact.counts(),
            exact.observed_detail_manifest_root(),
            exact.observed_detail_cell_count() + 1,
            exact.detail_cells_matched(),
            exact.first_unexpected_cell().cloned(),
        )
        .expect("intrinsically valid wrong detail count");
        let wrong_count_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &wrong_count)
                .expect("detail count disagreement is a checked mismatch");
        assert_eq!(wrong_count_report.unexpected_mismatches(), 1);
        assert!(!wrong_count_report.results()[detail_index].matched());

        let mut wrong_matched = presented.clone();
        wrong_matched[detail_index] = BaseE2ePresentedRowResultV1::new_with_detail_manifest(
            journey,
            exact.row_id().clone(),
            exact.semantic_manifest_root(),
            exact.observed(),
            exact.counts(),
            exact.observed_detail_manifest_root(),
            exact.observed_detail_cell_count(),
            exact.detail_cells_matched() - 1,
            exact.first_unexpected_cell().cloned(),
        )
        .expect("intrinsically valid wrong matched-detail count");
        let wrong_matched_report =
            compare_base_e2e_journey_results_v1(&projection, journey, &wrong_matched)
                .expect("matched-detail disagreement is a checked mismatch");
        assert_eq!(wrong_matched_report.unexpected_mismatches(), 1);
        assert!(!wrong_matched_report.results()[detail_index].matched());
    }

    #[test]
    fn result_join_rejects_missing_first_middle_and_last_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let harness = harness();
        let journey = BaseE2eJourneyV1::CanonicalRunnerV2;
        let complete = presented_results(&projection, journey, &harness);
        for index in [0, complete.len() / 2, complete.len() - 1] {
            let mut missing = complete.clone();
            missing.remove(index);
            let error = compare_base_e2e_journey_results_v1(&projection, journey, &missing)
                .expect_err("missing row must refuse");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
        }
    }

    #[test]
    fn result_join_rejects_extra_duplicate_and_reordered_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let harness = harness();
        let journey = BaseE2eJourneyV1::CanonicalRunnerV2;
        let complete = presented_results(&projection, journey, &harness);

        let mut extra = complete.clone();
        extra.push(reconstruct_result(
            &complete[0],
            journey,
            StableTokenV2::new("unmapped-extra-row").expect("extra ID"),
            complete[0].semantic_manifest_root(),
        ));
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &extra)
                .expect_err("extra row must refuse")
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let mut duplicate = complete.clone();
        duplicate.push(complete[0].clone());
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &duplicate)
                .expect_err("duplicate row must refuse")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut reordered = complete.clone();
        reordered.swap(0, 1);
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &reordered)
                .expect_err("reordered rows must refuse")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
    }

    #[test]
    fn result_join_rejects_stale_unmapped_and_cross_journey_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let harness = harness();
        let journey = BaseE2eJourneyV1::CanonicalRunnerV2;
        let complete = presented_results(&projection, journey, &harness);

        let mut stale = complete.clone();
        stale[0] = reconstruct_result(
            &complete[0],
            journey,
            complete[0].row_id().clone(),
            hash_domain("projection-test-stale-manifest.v1", b"stale"),
        );
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &stale)
                .expect_err("stale semantic root must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mut unmapped = complete.clone();
        unmapped[0] = reconstruct_result(
            &complete[0],
            journey,
            StableTokenV2::new("unmapped-row").expect("unmapped ID"),
            complete[0].semantic_manifest_root(),
        );
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &unmapped)
                .expect_err("unmapped row must refuse")
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let mut cross_journey = complete.clone();
        cross_journey[0] = reconstruct_result(
            &complete[0],
            BaseE2eJourneyV1::VerifierV2,
            complete[0].row_id().clone(),
            complete[0].semantic_manifest_root(),
        );
        assert_eq!(
            compare_base_e2e_journey_results_v1(&projection, journey, &cross_journey)
                .expect_err("cross-journey row must refuse")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }
}
