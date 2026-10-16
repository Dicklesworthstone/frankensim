//! Frozen, bounded Runner V2 base contracts.
//!
//! This TOOL crate owns declarations and pure validation only. It does not
//! execute cases, emit lifecycle records, parse hostile wire bytes, access a
//! filesystem, spawn a process, publish a bundle, verify scientific evidence,
//! or mint admission authority.
//!
//! The primary declaration surfaces are:
//!
//! - [`identity`] for presented roots, the exact 43-row wrapper inventory,
//!   constructor-owner handoff, and root-free evaluator-member guard;
//! - [`extension`] for fixed codec IDs, logical extents, registered
//!   role/unit/axis descriptors, and exact conversions;
//! - [`limits`] and [`budget`] for bounded profile limits and registry-aware
//!   budget admission;
//! - [`coverage`] for source-authoritative coverage and leaf-close manifests;
//! - [`projection`] for independent E2E/source-closure projections;
//! - [`schema_impact`] for source-frozen canonical-frame evolution graphs; and
//! - [`logging`] for deterministic base-E2E and leaf-close evidence logs.
//!
//! These modules expose declaration and pure-validation contracts only. A
//! successful reconstruction proves exact agreement with the relevant frozen
//! data; it never proves that a process ran or that an artifact, resource,
//! durability event, scientific claim, or authority grant exists.

#![deny(unsafe_code)]

mod canonical;

pub mod budget;
pub mod capability;
pub mod catalog;
pub mod command;
pub mod construction;
pub mod coverage;
pub mod dependency;
pub mod diagnostic;
pub mod extension;
pub mod identity;
pub mod limits;
pub mod logging;
pub mod path;
pub mod projection;
pub mod publication;
pub mod runner_v2;
pub mod schema_impact;
pub mod state;
pub mod value;

pub use budget::{
    AdmittedRunnerBudgetsV2, RegistryBoundAdmittedRunnerBudgetsV2, RegistryBoundRunnerBudgetsV2,
    RunnerBudgetsV2,
};
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
pub use command::{
    CommandIntentV2, CommandSelectionV2, CommandSelectorCardinalityV2,
    CommandSelectorExpectationV2, CommandSelectorFieldV2, CommandSelectorPresenceV2,
    CommandSelectorUsageKindV2, CommandSelectorUsageV2, validate_command_selector_presence_v2,
};
pub use construction::{
    ConstructionErrorKindV2, ConstructionErrorV2, ConstructionObservedDataClassV2,
};
pub use coverage::{
    BaseCoverageCaseDeclarationV1, BaseCoverageCheckedReportV1, BaseCoverageExecutableSubsetV1,
    BaseCoverageManifestCaseV1, BaseCoverageManifestClassV1, BaseCoverageManifestV1,
    BaseCoveragePresentedOutcomeV1, BaseCoveragePresentedResultV1,
};
pub use dependency::{
    CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1, DependencyOwnerPhaseV1, DependencyPolicyRowV1,
    DependencyRouteQualifierV1, DependencySourceIdentityV1, DependencySourceRouteV1,
    EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1, PresentedDependencyRouteV1,
    current_direct_dependency_declaration_root_v1, validate_current_direct_dependencies_v1,
    validate_eventual_direct_dependencies_v1,
};
pub use diagnostic::{
    ActionableDiagnosticV2, DiagnosticCodeRefV2, DiagnosticOverflowRefV2,
    RegisteredDecisionDetailProjectionV2, RepairActionV2,
};
pub use extension::{
    ArtifactCodecIdV2, BaseExtensionRegistryProjectionV2, LogicalExtentFieldV1, LogicalExtentV2,
    LogicalUnitScaleToCanonicalV2, RegisteredArtifactRoleDescriptorV2,
    RegisteredLogicalExtentAxisDescriptorV2, RegisteredLogicalUnitDescriptorV2,
    convert_rational_quantity_v2, normalized_unit_scale_ratio_v2,
};
pub use identity::{
    ALL_PRESENTED_IDENTITY_DESCRIPTORS_V1, CancelledStopRootV2, ConstructorOwnerHandoffEntryV1,
    ConstructorOwnerHandoffProjectionV1, ConstructorOwnerV1, DigestValueV2,
    DrainedInternalErrorRootV2, FROZEN_CONSTRUCTOR_OWNER_HANDOFF_ENTRIES_V1, NoClaimScopeRootV1,
    PresentedConstructorOwnerHandoffEntryV1, PresentedIdentityDescriptorV1,
    PresentedRootFreeEvaluatorMemberGuardV1, RootFreeEvaluatorMemberGuardDescriptorV1,
    RootFreeEvaluatorMemberGuardProjectionV1, RootFreeEvaluatorMemberV1, TimedOutStopRootV2,
};
pub use limits::{RunnerLimitsV2, RunnerLimitsValidationReportV2, RunnerLimitsViolationV2};
pub use logging::{
    BaseE2eLogEventV1, BaseE2eLogFieldCodeV1, BaseE2eLogFieldV1, BaseE2eLogKindV1, BaseE2eLogV1,
    BaseE2eOutcomeV1, SymbolicReproductionArgV1,
};
pub use path::{
    ContentStoreObjectKeyV1, LogicalBundlePathV1, PathSetAdjudicationV1,
    adjudicate_content_store_object_key_set, adjudicate_logical_bundle_path_set,
};
#[allow(
    deprecated,
    reason = "retain deprecated coverage compatibility reexports while directing new callers to coverage::*"
)]
pub use projection::{
    BASE_E2E_DETAIL_CELL_MAX_ENCODED_BYTES_V1, BASE_E2E_ROW_CONTRACT_DIVERGENCE_ID_V1,
    BaseCoverageClassV1, BaseCoverageInventoryV1, BaseCoverageSourceCaseV1,
    BaseE2eCapabilityRefusalStageV1, BaseE2eCaseKindV1, BaseE2eDecisionDetailManifestV1,
    BaseE2eDetailCellV1, BaseE2eDetailDivergenceV1, BaseE2eDetailManifestV1,
    BaseE2eDetailPayloadV1, BaseE2eExpectedDecisionV1, BaseE2eHarnessIdentityV1,
    BaseE2eJourneyComparisonReportV1, BaseE2eJourneyExecutionReportV1, BaseE2eJourneyProjectionV1,
    BaseE2eJourneyV1, BaseE2eMatchedPartitionV1, BaseE2eObservedCountsV1,
    BaseE2ePathAdjudicationDetailV1, BaseE2ePresentedRowResultV1, BaseE2eProjectionReportV1,
    BaseE2eProjectionRowV1, BaseE2eRetainedArtifactClaimV1, BaseE2eRowResultV1,
    BaseSourceClosureEntryV1, BaseSourceClosureInputV1, BaseSourceOwnerV1, BaseSourceRouteV1,
    BaseSourceSnapshotPolicyV1, RunnerV2BaseE2eProjectionV1, RunnerV2BaseSourceClosureV1,
    compare_base_e2e_journey_results_v1, join_base_e2e_journey_results_v1, run_base_e2e_journey_v1,
    run_base_e2e_projection_v1,
};
pub use publication::{PublicationSelectionV2, PublicationTargetV2, SymbolicCommandResultPlanV2};
pub use schema_impact::{
    CanonicalFieldCodeV1, CanonicalFieldLayoutV1, CanonicalFieldNameV1, CanonicalFieldWireKindV1,
    CanonicalFrameVersionV1, CanonicalNominalRootRoleIdV1, CanonicalRustSchemaNameV1,
    CanonicalSchemaAuthorityStateV1, CanonicalSchemaAuthoritySurfaceV1, CanonicalSchemaDomainV1,
    CanonicalSchemaFieldDescriptorV1, CanonicalSchemaFrameBindingV1,
    CanonicalSchemaFrameDescriptorV1, CanonicalSchemaIdV1, CanonicalSchemaMagicV1,
    CanonicalSchemaSlotUseV1, CanonicalSchemaVersionSlotDescriptorV1, CanonicalSemanticTypeIdV1,
    CanonicalSlotIdV1, CanonicalVersionSlotCodeV1, CompatibleSourceSnapshotV1,
    FrozenBaseNominalRootRegistryFragmentV1, LeafExtensionNominalRootRegistryFragmentV1,
    LegacyNestedContainerRefV1, NominalRootRegistryIdV1, NominalRootRegistryKindV1,
    NominalRootRoleRefV1, RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1,
    RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1, RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1,
    SchemaImpactLeafIdV1, SchemaImpactManifestEntryV1, SchemaImpactManifestRelationV1,
    SchemaImpactManifestV1, SchemaImpactNoClaimV1, SchemaImpactRowV1,
    runner_v2_base_schema_impact_manifest_v1,
};
pub use state::{
    NotRunBasisErrorV2, NotRunBasisV2, NotRunCauseV2, StateValidationInputV2, validate_state_v2,
};
pub use value::{
    CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1, CASE_SEED_DERIVATION_NO_CLAIM_V1,
    CaseSeedDerivationDomainRegistryRootV1, CaseSeedDerivationDomainRegistryV1,
    CaseSeedDerivationDomainRootV1, CaseSeedPolicyRootV2, CaseSeedPolicyV2,
    CaseSeedProvenanceRootV2, CaseSeedResolutionRootV2, CaseSeedResolutionV2, DecimalV2,
    FixedManifestSeedBindingV2, InvocationDerivedSeedBindingV2, InvocationSeedSelectionRootV2,
    InvocationSeedSelectionV2, NoRandomnessSeedBindingV2, NumericValueV2, RationalV2,
    RegisteredCaseSeedDerivationDomainV1, SEED_CLI_FLAG_V2, SEED_CLI_PREFIX_V2,
    SEED_MATERIAL_BYTES_V2, SEED_MATERIAL_LOWER_HEX_BYTES_V2, STABLE_CASE_ID_MAX_BYTES_V2,
    SeedErrorV2, SeedGeneratorVersionV1, SeedInapplicableCodeV1, SeedMaterialRootV2,
    SeedMaterialV2, SeedMinimizerVersionV1, SeedSchemaDescriptorV1, StableCaseIdentityRootV2,
    StableCaseIdentityV2, StableTokenV2, TextV2, TypedOptionV1, TypedValueV2, UnitV2,
};
