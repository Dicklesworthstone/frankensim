//! Frozen, bounded Runner V2 base contracts.
//!
//! This TOOL crate owns declarations and pure validation only. It does not
//! execute cases, emit lifecycle records, parse hostile wire bytes, access a
//! filesystem, spawn a process, publish a bundle, verify scientific evidence,
//! or mint admission authority.

#![deny(unsafe_code)]

mod canonical;

pub mod budget;
pub mod capability;
pub mod catalog;
pub mod command;
pub mod construction;
pub mod dependency;
pub mod diagnostic;
pub mod identity;
pub mod limits;
pub mod logging;
pub mod path;
pub mod projection;
pub mod publication;
pub mod state;
pub mod value;

pub use budget::{AdmittedRunnerBudgetsV2, RunnerBudgetsV2};
pub use capability::{
    NarrowedPolicyViewV2, RootCapabilityPolicySetV2, RootCapabilityPolicyV2,
    RootPolicyRegistryProjectionV2,
};
pub use catalog::{
    ArtifactDispositionV2, ArtifactRoleV2, DestinationAdmissionModeV2, DiagnosticCodeV2,
    DigestRoleV2, LifecycleRecordKindV2, LogicalExtentAxisV2, LogicalUnitV2, NotRunCauseCodeV2,
    PlatformPathProfileV2, ProofExitV2, PublicationProtocolV2, RUNNER_SPEC_V2_API_GENERATION,
    RUNNER_V2_PREDECESSOR_POLICY, RUNNER_V2_WIRE_VERSION, RefusedReasonV2, RepairActionKindV2,
    RetryabilityV2, RootCapabilityAccessV2, RootCapabilityRightV2, RootClassV2, RunProfileV2,
    RunnerApiGeneration, RunnerCommandV2, RunnerWireVersion, StateBearingRecordRoleV2,
    TypedOptionTagV1, TypedValueTagV2, WirePredecessorPolicyV1,
};
pub use command::{CommandIntentV2, CommandSelectionV2};
pub use construction::{ConstructionErrorKindV2, ConstructionErrorV2};
pub use dependency::{
    CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1, DependencyOwnerPhaseV1, DependencyPolicyRowV1,
    DependencySourceRouteV1, EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1,
    PresentedDependencyRouteV1, validate_current_direct_dependencies_v1,
    validate_eventual_direct_dependencies_v1,
};
pub use diagnostic::{
    ActionableDiagnosticV2, DiagnosticCodeRefV2, DiagnosticOverflowRefV2, RepairActionV2,
};
pub use identity::{
    CancelledStopRootV2, DigestValueV2, DrainedInternalErrorRootV2, NoClaimScopeRootV1,
    TimedOutStopRootV2,
};
pub use limits::{RunnerLimitsV2, RunnerLimitsViolationV2};
pub use path::{
    ContentStoreObjectKeyV1, LogicalBundlePathV1, PathSetAdjudicationV1,
    adjudicate_content_store_object_key_set, adjudicate_logical_bundle_path_set,
};
pub use projection::{
    BaseCoverageClassV1, BaseCoverageInventoryV1, BaseCoverageSourceCaseV1,
    BaseSourceClosureEntryV1, BaseSourceClosureInputV1, RunnerV2BaseE2eProjectionV1,
    RunnerV2BaseSourceClosureV1, run_base_e2e_projection_v1,
};
pub use publication::{PublicationSelectionV2, PublicationTargetV2, SymbolicCommandResultPlanV2};
pub use state::{
    NotRunBasisErrorV2, NotRunBasisV2, NotRunCauseV2, StateValidationInputV2, validate_state_v2,
};
pub use value::{
    DecimalV2, NumericValueV2, RationalV2, StableTokenV2, TextV2, TypedOptionV1, TypedValueV2,
    UnitV2,
};
