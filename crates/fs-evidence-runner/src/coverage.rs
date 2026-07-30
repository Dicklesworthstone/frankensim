//! Exact, result-free coverage declarations and checked result accounting.
//!
//! This module deliberately separates three things that are easy to conflate:
//!
//! 1. the immutable source-case manifest;
//! 2. a caller-selected, manifest-ordered executable subset; and
//! 3. presented result records joined back to that exact subset.
//!
//! The frozen base enumerates the current 217 non-manifest Rust tests and all 78
//! `compile_fail` contracts. The historical aggregate is recorded as 130
//! ratified cases plus an eighty-seven-case delta, but this source does not pretend
//! to recover a per-case historical membership label that was never retained.
//! Every current Rust test has a handwritten evidence-class assignment:
//! focused unit, boundary, property/metamorphic, schema/descriptor, mutation,
//! or no-mock integration. Coverage-manifest contract tests are enumerated
//! separately.
//!
//! [`BaseCoverageManifestV1`] is the sole source-authoritative, result-free
//! AC38 inventory. `projection::BaseCoverageInventoryV1` is an older,
//! compatibility-only auxiliary projection: it cannot replace, widen, or
//! prove this manifest.
//!
//! No type in this module discovers tests, reads source files, executes a test,
//! interprets scientific evidence, or certifies that a presented result is
//! true. A green checked report proves only exact manifest/result accounting
//! for the selected IDs. The owning test runner or external harness remains
//! responsible for execution and evidence semantics.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::{
    DiagnosticCodeV2, LogicalUnitV2, RUNNER_SPEC_V2_API_GENERATION, RUNNER_V2_PREDECESSOR_POLICY,
    RUNNER_V2_WIRE_VERSION, RunnerApiGeneration, RunnerWireVersion, WirePredecessorPolicyV1,
};
use crate::construction::{
    ConstructionErrorKindV2, ConstructionErrorV2, ConstructionObservedDataClassV2,
};
use crate::identity::{BuildIdentityRootV2, SourceIdentityRootV2, ToolchainIdentityRootV2};
use crate::value::{NumericValueV2, SeedInapplicableCodeV1, SeedMaterialV2, StableTokenV2, UnitV2};
use core::num::NonZeroU16;
use fs_blake3::ContentHash;
use std::collections::{BTreeMap, BTreeSet};

/// Domain for the exact result-free coverage manifest.
pub const BASE_COVERAGE_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-manifest.v1";

/// Domain for a caller-selected exact executable subset.
pub const BASE_COVERAGE_SELECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-selection.v1";

/// Domain for one presented coverage result.
pub const BASE_COVERAGE_PRESENTED_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-presented-result.v1";

/// Domain for an exact checked coverage report.
pub const BASE_COVERAGE_CHECKED_REPORT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-checked-report.v1";

/// Domain for the AC53 source-authoritative, full-set-only close manifest.
pub const BASE_COVERAGE_CLOSE_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-manifest.v1";
/// Domain for all Five Explicits bound into one close cell.
pub const BASE_COVERAGE_CLOSE_FIVE_EXPLICITS_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits.v1";
/// Domain for one independently rooted semantic numeric profile.
pub const BASE_COVERAGE_CLOSE_NUMERIC_PROFILE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-numeric-profile.v1";

/// Domain for one AC53 result-free close cell.
pub const BASE_COVERAGE_CLOSE_CELL_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-cell.v1";

/// Domain for the closed AC53 inapplicability and Unsupported reason registry.
pub const BASE_COVERAGE_CLOSE_REASON_REGISTRY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-reason-registry.v1";

/// Domain for one caller-presented AC53 full-set result.
pub const BASE_COVERAGE_CLOSE_PRESENTED_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-presented-result.v1";

/// Domain for the checked AC53 full-set close report.
pub const BASE_COVERAGE_CLOSE_REPORT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-report.v1";

/// Domain for one immutable downstream-contribution declaration.
pub const BASE_COVERAGE_CLOSE_DOWNSTREAM_CONTRIBUTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-downstream-contribution.v1";

/// Domain for the additive, result-free downstream contribution V2 frame.
pub const BASE_COVERAGE_CLOSE_DOWNSTREAM_CONTRIBUTION_DOMAIN_V2: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-downstream-contribution.v2";

/// Exact route ID for the base leaf's Phase-1 contract contribution.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1: &str =
    "runner-v2.route.24-1-1-1.to.phase1-contract.v1";
/// Source leaf that semantically owns the immutable contribution.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1: &str =
    "frankensim-epic-foundations-huq.24.1.1.1";
/// Sole downstream owner that may execute the Phase-1 contract route.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_EXECUTION_OWNER_V1: &str =
    "frankensim-epic-foundations-huq.24.1.3.2";
/// Exact release driver owned by the downstream Phase-1 contract leaf.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_DRIVER_V1: &str = "runner-v2-phase1-contract-e2e-driver";
/// Exact downstream-owned POSIX wrapper route.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_POSIX_ROUTE_V1: &str =
    "scripts/ci/e2e_runner_v2_phase1_contract.sh";
/// Exact downstream-owned native-Windows PowerShell route.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_WINDOWS_ROUTE_V1: &str =
    "scripts/ci/e2e_runner_v2_phase1_contract.ps1";
/// Exact downstream-owned immutable case-manifest path.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_CASE_MANIFEST_PATH_V1: &str =
    "scripts/ci/manifests/runner_v2_phase1_contract_cases.v1.json";
/// Exact no-execution boundary for the immutable Phase-1 contribution.
pub const RUNNER_V2_PHASE_ONE_CONTRACT_NO_CLAIM_V1: &str =
    "phase1-contract-contribution-is-result-free-and-does-not-prove-owner-execution";
/// Exact no-execution boundary for the designated Phase-1 observer contract.
pub const RUNNER_V2_PHASE_ONE_OBSERVER_NO_CLAIM_V1: &str =
    "phase1-observer-contract-proves-no-owner-execution";
/// Exact no-execution boundary for its separate Deferred evidence envelope.
pub const RUNNER_V2_PHASE_ONE_DEFERRED_ENVELOPE_NO_CLAIM_V1: &str =
    "phase1-contract-deferred-envelope-proves-no-downstream-execution-or-success";

/// Domain for typed facet-applicability evidence identities.
pub const BASE_COVERAGE_CLOSE_APPLICABILITY_EVIDENCE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-close-applicability-evidence.v1";

/// Frozen descriptor shared by the AC54 locally derived nominal-root family.
///
/// A descriptor fixes schema identity and no-claim semantics only. It does not
/// construct a root, validate bytes, establish execution, or mint authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseNominalRootDescriptorV1 {
    schema_name: &'static str,
    domain: &'static str,
    api_generation: RunnerApiGeneration,
    wire_version: RunnerWireVersion,
    predecessor_policy: WirePredecessorPolicyV1,
    no_claim: &'static str,
}

impl BaseCoverageCloseNominalRootDescriptorV1 {
    const fn frozen(
        schema_name: &'static str,
        domain: &'static str,
        no_claim: &'static str,
    ) -> Self {
        Self {
            schema_name,
            domain,
            api_generation: RUNNER_SPEC_V2_API_GENERATION,
            wire_version: RUNNER_V2_WIRE_VERSION,
            predecessor_policy: RUNNER_V2_PREDECESSOR_POLICY,
            no_claim,
        }
    }

    /// Exact lowercase kebab-case schema name.
    #[must_use]
    pub const fn schema_name(self) -> &'static str {
        self.schema_name
    }

    /// Exact version-matched domain owned by this schema.
    ///
    /// The immutable 47-row FrozenBase inventory is entirely `.v1`; a
    /// source-frozen leaf-extension descriptor may instead name its exact
    /// `.v2` role without weakening that frozen base.
    #[must_use]
    pub const fn domain(self) -> &'static str {
        self.domain
    }

    /// Public Runner API generation, exactly two.
    #[must_use]
    pub const fn api_generation(self) -> RunnerApiGeneration {
        self.api_generation
    }

    /// Frozen wire version, exactly one.
    #[must_use]
    pub const fn wire_version(self) -> RunnerWireVersion {
        self.wire_version
    }

    /// Frozen wire-predecessor policy, exactly no predecessor.
    #[must_use]
    pub const fn predecessor_policy(self) -> WirePredecessorPolicyV1 {
        self.predecessor_policy
    }

    /// Exact no-claim boundary for this nominal role.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        self.no_claim
    }
}

macro_rules! define_base_coverage_close_nominal_root_v1 {
    ($name:ident, $schema_name:literal, $domain:literal, $no_claim:literal) => {
        #[doc = concat!("Nominal AC54 content identity for `", $schema_name, "`.")]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            content_hash: ContentHash,
        }

        impl $name {
            /// Frozen schema/domain/version/no-claim descriptor.
            pub const DESCRIPTOR: BaseCoverageCloseNominalRootDescriptorV1 =
                BaseCoverageCloseNominalRootDescriptorV1::frozen($schema_name, $domain, $no_claim);

            #[allow(
                dead_code,
                reason = "AC54 consumer migration follows the nominal-role freeze"
            )]
            const fn from_content_hash(content_hash: ContentHash) -> Self {
                Self { content_hash }
            }

            /// Read the role-bound content hash without changing nominal type.
            #[must_use]
            pub const fn content_hash(self) -> ContentHash {
                self.content_hash
            }
        }
    };
}

define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseNumericInputsRootV1,
    "numeric-inputs",
    "org.frankensim.fs-evidence-runner.base-coverage-close-numeric-inputs.v1",
    "numeric-inputs-root-does-not-prove-runtime-observation"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseNumericGrantsRootV1,
    "numeric-grants",
    "org.frankensim.fs-evidence-runner.base-coverage-close-numeric-grants.v1",
    "numeric-grants-root-does-not-prove-resource-acquisition"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseExpectedNumericObservationsRootV1,
    "expected-numeric-observations",
    "org.frankensim.fs-evidence-runner.base-coverage-close-expected-numeric-observations.v1",
    "expected-numeric-observations-root-does-not-prove-runtime-observation"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseActualNumericObservationsRootV1,
    "actual-numeric-observations",
    "org.frankensim.fs-evidence-runner.base-coverage-close-actual-numeric-observations.v1",
    "actual-numeric-observations-root-does-not-prove-scientific-validity"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseSeedRequirementRootV1,
    "seed-requirement",
    "org.frankensim.fs-evidence-runner.base-coverage-close-seed-requirement.v1",
    "seed-requirement-root-does-not-prove-seed-resolution"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseObservedSeedDispositionRootV1,
    "observed-seed-disposition",
    "org.frankensim.fs-evidence-runner.base-coverage-close-observed-seed-disposition.v1",
    "observed-seed-disposition-root-does-not-prove-randomness-quality"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseBudgetSetRootV1,
    "budget-set",
    "org.frankensim.fs-evidence-runner.base-coverage-close-budget-set.v1",
    "budget-set-root-does-not-prove-resource-enforcement"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseVersionRequirementsRootV1,
    "version-requirements",
    "org.frankensim.fs-evidence-runner.base-coverage-close-version-requirements.v1",
    "version-requirements-root-does-not-prove-runtime-version-match"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseObservedVersionsRootV1,
    "observed-versions",
    "org.frankensim.fs-evidence-runner.base-coverage-close-observed-versions.v1",
    "observed-versions-root-does-not-prove-source-or-build-trust"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseCapabilityDescriptorRootV1,
    "capability-descriptor",
    "org.frankensim.fs-evidence-runner.base-coverage-close-capability-descriptor.v1",
    "capability-descriptor-root-proves-declared-semantics-not-acquisition-use-return-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseCapabilityRegistryRootV1,
    "capability-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-capability-registry.v1",
    "capability-registry-root-proves-declared-membership-not-acquisition-use-return-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseCapabilityProfileRegistryRootV1,
    "capability-profile-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-capability-profile-registry.v1",
    "capability-profile-registry-root-proves-declared-profile-membership-not-acquisition-use-return-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseCapabilityContractRootV1,
    "capability-contract",
    "org.frankensim.fs-evidence-runner.base-coverage-close-capability-contract.v1",
    "capability-contract-root-proves-declared-bounds-not-acquisition-use-return-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseObservedCapabilitySetsRootV1,
    "observed-capability-sets",
    "org.frankensim.fs-evidence-runner.base-coverage-close-observed-capability-sets.v1",
    "observed-capability-sets-root-proves-structural-reconciliation-not-resource-return-effect-success-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseNoClaimRootV1,
    "no-claim",
    "org.frankensim.fs-evidence-runner.base-coverage-close-no-claim.v1",
    "no-claim-root-does-not-mint-scientific-admission-or-publication-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseFiveExplicitsProfileRegistryRootV1,
    "five-explicits-profile-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-profile-registry.v1",
    "five-explicits-profile-registry-root-does-not-prove-cell-coverage"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseFiveExplicitsCellOracleRootV1,
    "five-explicits-cell-oracle",
    "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-cell-oracle.v1",
    "five-explicits-cell-oracle-root-does-not-prove-execution"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseFiveExplicitsDeclarationRootV1,
    "five-explicits-declaration",
    "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-declaration.v1",
    "five-explicits-declaration-root-does-not-prove-runtime-observation"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRuntimeExplicitsRootV1,
    "runtime-explicits",
    "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-explicits.v1",
    "runtime-explicits-root-does-not-prove-result-correctness"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseLogExplicitJoinRootV1,
    "log-explicit-join",
    "org.frankensim.fs-evidence-runner.base-coverage-close-log-explicit-join.v1",
    "log-explicit-join-root-does-not-prove-log-retention"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseReproductionExplicitJoinRootV1,
    "reproduction-explicit-join",
    "org.frankensim.fs-evidence-runner.base-coverage-close-reproduction-explicit-join.v1",
    "reproduction-explicit-join-root-does-not-prove-reproduction-success"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseVersionSchemaRootV1,
    "version-schema",
    "org.frankensim.fs-evidence-runner.base-coverage-close-version-schema.v1",
    "version-schema-root-does-not-prove-runtime-version-match"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseTargetRequirementRootV1,
    "target-requirement",
    "org.frankensim.fs-evidence-runner.base-coverage-close-target-requirement.v1",
    "target-requirement-root-does-not-prove-target-execution"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageClosePlatformMatrixRootV1,
    "platform-matrix",
    "org.frankensim.fs-evidence-runner.base-coverage-close-platform-matrix.v1",
    "platform-matrix-root-does-not-prove-platform-execution"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseConfigurationProfileRootV1,
    "configuration-profile",
    "org.frankensim.fs-evidence-runner.base-coverage-close-configuration-profile.v1",
    "configuration-profile-root-does-not-prove-configuration-application"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseFeatureSetRootV1,
    "feature-set",
    "org.frankensim.fs-evidence-runner.base-coverage-close-feature-set.v1",
    "feature-set-root-does-not-prove-feature-activation"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseCapabilityPolicyRequirementsRootV1,
    "capability-policy-requirements",
    "org.frankensim.fs-evidence-runner.base-coverage-close-capability-policy-requirements.v1",
    "capability-policy-requirements-root-does-not-grant-or-acquire-capabilities"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseDeferredObservationContractRootV1,
    "deferred-observation-contract",
    "org.frankensim.fs-evidence-runner.base-coverage-close-deferred-observation-contract.v1",
    "deferred-observation-contract-root-does-not-prove-downstream-execution"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRuntimeObservationDispositionRootV1,
    "runtime-observation-disposition",
    "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation-disposition.v1",
    "runtime-observation-disposition-root-does-not-prove-observation-validity"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRuntimeObservationAggregateRootV1,
    "runtime-observation-aggregate",
    "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation-aggregate.v1",
    "runtime-observation-aggregate-root-does-not-prove-result-correctness"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseBudgetReconciliationRootV1,
    "budget-reconciliation",
    "org.frankensim.fs-evidence-runner.base-coverage-close-budget-reconciliation.v1",
    "budget-reconciliation-root-proves-structural-accounting-not-enforcement-return-or-effect-success"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRuntimeObservationRootV1,
    "runtime-observation",
    "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation.v1",
    "runtime-observation-root-proves-reported-observation-not-scientific-validity-or-effect-success"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseEvidenceEnvelopeRootV1,
    "evidence-envelope",
    "org.frankensim.fs-evidence-runner.base-coverage-close-evidence-envelope.v1",
    "evidence-envelope-root-proves-structural-evidence-binding-not-evidence-validity-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRetainedArtifactStateRootV1,
    "retained-artifact-state",
    "org.frankensim.fs-evidence-runner.base-coverage-close-retained-artifact-state.v1",
    "retained-artifact-state-root-proves-reported-retention-state-not-durability-completeness-or-validity"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseSafePartialEvidenceRootV1,
    "safe-partial-evidence",
    "org.frankensim.fs-evidence-runner.base-coverage-close-safe-partial-evidence.v1",
    "safe-partial-evidence-root-proves-explicit-partial-state-not-completeness-success-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseResourceReconciliationRootV1,
    "resource-reconciliation",
    "org.frankensim.fs-evidence-runner.base-coverage-close-resource-reconciliation.v1",
    "resource-reconciliation-root-proves-structural-accounting-not-leak-freedom-drain-success-or-effect-success"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseExecutionCompletenessRootV1,
    "execution-completeness",
    "org.frankensim.fs-evidence-runner.base-coverage-close-execution-completeness.v1",
    "execution-completeness-root-proves-declared-terminal-coverage-not-result-correctness-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseNotObservedReasonRegistryRootV1,
    "not-observed-reason-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-not-observed-reason-registry.v1",
    "not-observed-reason-registry-root-proves-frozen-reason-membership-not-runtime-cause-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseDeferredReasonRegistryRootV1,
    "deferred-reason-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-deferred-reason-registry.v1",
    "deferred-reason-registry-root-proves-frozen-reason-membership-not-downstream-execution-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseAttemptIdentityRootV1,
    "attempt-identity",
    "org.frankensim.fs-evidence-runner.base-coverage-close-attempt-identity.v1",
    "attempt-identity-root-proves-structural-attempt-binding-not-process-execution-success-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    CompatibleSourceSnapshotRootV1,
    "compatible-source-snapshot",
    "org.frankensim.fs-evidence-runner.base-source-snapshot.v1",
    "compatible-source-snapshot-root-proves-exact-source-closure-identity-not-build-execution-or-schema-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseNominalRootRegistryRootV1,
    "nominal-root-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-nominal-root-registry.v1",
    "nominal-root-registry-root-proves-frozen-role-descriptors-not-root-construction-validity-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    SchemaImpactRowRootV1,
    "schema-impact-row",
    "org.frankensim.fs-evidence-runner.schema-impact-row.v1",
    "schema-impact-row-root-proves-checked-schema-declaration-not-parser-safety-migration-success-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    SchemaImpactManifestRootV1,
    "schema-impact-manifest",
    "org.frankensim.fs-evidence-runner.schema-impact-manifest.v1",
    "schema-impact-manifest-root-proves-exact-schema-inventory-and-dag-not-implementation-correctness-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1,
    "registered-extension-capability-descriptor",
    "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-descriptor.v1",
    "registered-extension-capability-descriptor-root-proves-declared-extension-semantics-not-acquisition-use-return-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1,
    "registered-extension-capability-registry",
    "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-registry.v1",
    "registered-extension-capability-registry-root-proves-declared-extension-membership-not-base-membership-acquisition-or-authority"
);
define_base_coverage_close_nominal_root_v1!(
    BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1,
    "registered-extension-capability-set",
    "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-set.v1",
    "registered-extension-capability-set-root-proves-structural-extension-membership-not-base-membership-acquisition-or-authority"
);

/// Exact number of base nominal-root roles before registered extensions.
pub const BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1: usize = 44;
/// Exact number of locally derived AC54-AC61 nominal-root roles.
pub const BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1: usize = 47;

/// AC54-AC61 nominal-root descriptors in frozen semantic-role order.
pub const BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTORS_V1:
    [BaseCoverageCloseNominalRootDescriptorV1; BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1] = [
    BaseCoverageCloseNumericInputsRootV1::DESCRIPTOR,
    BaseCoverageCloseNumericGrantsRootV1::DESCRIPTOR,
    BaseCoverageCloseExpectedNumericObservationsRootV1::DESCRIPTOR,
    BaseCoverageCloseActualNumericObservationsRootV1::DESCRIPTOR,
    BaseCoverageCloseSeedRequirementRootV1::DESCRIPTOR,
    BaseCoverageCloseObservedSeedDispositionRootV1::DESCRIPTOR,
    BaseCoverageCloseBudgetSetRootV1::DESCRIPTOR,
    BaseCoverageCloseVersionRequirementsRootV1::DESCRIPTOR,
    BaseCoverageCloseObservedVersionsRootV1::DESCRIPTOR,
    BaseCoverageCloseCapabilityDescriptorRootV1::DESCRIPTOR,
    BaseCoverageCloseCapabilityRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseCapabilityProfileRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseCapabilityContractRootV1::DESCRIPTOR,
    BaseCoverageCloseObservedCapabilitySetsRootV1::DESCRIPTOR,
    BaseCoverageCloseNoClaimRootV1::DESCRIPTOR,
    BaseCoverageCloseFiveExplicitsProfileRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseFiveExplicitsCellOracleRootV1::DESCRIPTOR,
    BaseCoverageCloseFiveExplicitsDeclarationRootV1::DESCRIPTOR,
    BaseCoverageCloseRuntimeExplicitsRootV1::DESCRIPTOR,
    BaseCoverageCloseLogExplicitJoinRootV1::DESCRIPTOR,
    BaseCoverageCloseReproductionExplicitJoinRootV1::DESCRIPTOR,
    BaseCoverageCloseVersionSchemaRootV1::DESCRIPTOR,
    BaseCoverageCloseTargetRequirementRootV1::DESCRIPTOR,
    BaseCoverageClosePlatformMatrixRootV1::DESCRIPTOR,
    BaseCoverageCloseConfigurationProfileRootV1::DESCRIPTOR,
    BaseCoverageCloseFeatureSetRootV1::DESCRIPTOR,
    BaseCoverageCloseCapabilityPolicyRequirementsRootV1::DESCRIPTOR,
    BaseCoverageCloseDeferredObservationContractRootV1::DESCRIPTOR,
    BaseCoverageCloseRuntimeObservationDispositionRootV1::DESCRIPTOR,
    BaseCoverageCloseRuntimeObservationAggregateRootV1::DESCRIPTOR,
    BaseCoverageCloseBudgetReconciliationRootV1::DESCRIPTOR,
    BaseCoverageCloseRuntimeObservationRootV1::DESCRIPTOR,
    BaseCoverageCloseEvidenceEnvelopeRootV1::DESCRIPTOR,
    BaseCoverageCloseRetainedArtifactStateRootV1::DESCRIPTOR,
    BaseCoverageCloseSafePartialEvidenceRootV1::DESCRIPTOR,
    BaseCoverageCloseResourceReconciliationRootV1::DESCRIPTOR,
    BaseCoverageCloseExecutionCompletenessRootV1::DESCRIPTOR,
    BaseCoverageCloseNotObservedReasonRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseDeferredReasonRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseAttemptIdentityRootV1::DESCRIPTOR,
    CompatibleSourceSnapshotRootV1::DESCRIPTOR,
    BaseCoverageCloseNominalRootRegistryRootV1::DESCRIPTOR,
    SchemaImpactRowRootV1::DESCRIPTOR,
    SchemaImpactManifestRootV1::DESCRIPTOR,
    BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1::DESCRIPTOR,
    BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1::DESCRIPTOR,
    BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1::DESCRIPTOR,
];

/// AC54-AC61 nominal-root descriptors in exact frozen order.
#[must_use]
pub const fn base_coverage_close_nominal_root_descriptors_v1()
-> &'static [BaseCoverageCloseNominalRootDescriptorV1] {
    &BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTORS_V1
}

/// Construct one crate-owned, source-frozen nominal-role descriptor for an
/// AC60 leaf-extension registry fragment.
///
/// This is deliberately crate-private: downstream leaves in this crate may
/// declare literals, while public callers cannot assemble an open descriptor
/// vector or mint registry membership.
pub(crate) const fn source_frozen_nominal_root_descriptor_v1(
    schema_name: &'static str,
    domain: &'static str,
    no_claim: &'static str,
) -> BaseCoverageCloseNominalRootDescriptorV1 {
    BaseCoverageCloseNominalRootDescriptorV1::frozen(schema_name, domain, no_claim)
}

fn exact_nominal_root_frame_hash_v1(
    frame: &CanonicalFrameV1,
    magic: &'static [u8],
    domain: &'static str,
    field: &'static str,
) -> Result<ContentHash, ConstructionErrorV2> {
    if !frame.as_bytes().starts_with(magic) {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            field,
            "the exact role-specific canonical-frame magic",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    Ok(frame.root(domain))
}

/// Seal the existing exact source-snapshot frame under its nominal role
/// without wrapping or hashing its digest a second time.
pub(crate) fn compatible_source_snapshot_root_from_exact_frame_v1(
    frame: &CanonicalFrameV1,
) -> Result<CompatibleSourceSnapshotRootV1, ConstructionErrorV2> {
    let content_hash = exact_nominal_root_frame_hash_v1(
        frame,
        b"FSBASESOURCESNAPSHOT\x01",
        CompatibleSourceSnapshotRootV1::DESCRIPTOR.domain(),
        "coverage.close.compatible_source_snapshot.frame",
    )?;
    Ok(CompatibleSourceSnapshotRootV1::from_content_hash(
        content_hash,
    ))
}

/// Seal one checked tagged nominal-registry fragment.
pub(crate) fn nominal_root_registry_root_from_exact_frame_v1(
    frame: &CanonicalFrameV1,
) -> Result<BaseCoverageCloseNominalRootRegistryRootV1, ConstructionErrorV2> {
    let content_hash = exact_nominal_root_frame_hash_v1(
        frame,
        b"FSCLOSENOMINALREG\x01",
        BaseCoverageCloseNominalRootRegistryRootV1::DESCRIPTOR.domain(),
        "coverage.close.nominal_root_registry.frame",
    )?;
    Ok(BaseCoverageCloseNominalRootRegistryRootV1::from_content_hash(content_hash))
}

/// Seal one source-frozen, validated schema-impact row.
pub(crate) fn schema_impact_row_root_from_exact_frame_v1(
    frame: &CanonicalFrameV1,
) -> Result<SchemaImpactRowRootV1, ConstructionErrorV2> {
    let content_hash = exact_nominal_root_frame_hash_v1(
        frame,
        b"FSSCHEMAIMPACTROW\x01",
        SchemaImpactRowRootV1::DESCRIPTOR.domain(),
        "coverage.close.schema_impact.row_frame",
    )?;
    Ok(SchemaImpactRowRootV1::from_content_hash(content_hash))
}

/// Seal one source-frozen, validated schema-impact manifest.
pub(crate) fn schema_impact_manifest_root_from_exact_frame_v1(
    frame: &CanonicalFrameV1,
) -> Result<SchemaImpactManifestRootV1, ConstructionErrorV2> {
    let content_hash = exact_nominal_root_frame_hash_v1(
        frame,
        b"FSSCHEMAIMPACTMANIFEST\x01",
        SchemaImpactManifestRootV1::DESCRIPTOR.domain(),
        "coverage.close.schema_impact.manifest_frame",
    )?;
    Ok(SchemaImpactManifestRootV1::from_content_hash(content_hash))
}

/// Exact number of Rust tests in the ratified pre-manifest base inventory.
pub const BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1: usize = 130;

/// Exact Rust-test delta added by the same implementation train.
pub const BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1: usize = 87;

/// Exact current non-manifest Rust-test total frozen by this source inventory.
pub const BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1: usize =
    BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1
        + BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1;

/// Historical aggregate name retained for callers of the initial contract.
///
/// This is the total across all six Rust-test evidence classes, not the count
/// of [`BaseCoverageManifestClassV1::Unit`] alone.
#[deprecated(
    note = "use BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1 for the aggregate or BASE_COVERAGE_UNIT_CLASS_CASE_COUNT_V1 for the Unit class"
)]
pub const BASE_COVERAGE_UNIT_CASE_COUNT_V1: usize = BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1;

/// Exact number of focused unit-behavior cases.
pub const BASE_COVERAGE_UNIT_CLASS_CASE_COUNT_V1: usize = 12;

/// Exact number of boundary and checked-arithmetic cases.
pub const BASE_COVERAGE_BOUNDARY_CASE_COUNT_V1: usize = 53;

/// Exact number of property and metamorphic cases.
pub const BASE_COVERAGE_PROPERTY_METAMORPHIC_CASE_COUNT_V1: usize = 24;

/// Exact number of schema and descriptor cases.
pub const BASE_COVERAGE_SCHEMA_DESCRIPTOR_CASE_COUNT_V1: usize = 50;

/// Exact number of mutation and malformed-presentation cases.
pub const BASE_COVERAGE_MUTATION_CASE_COUNT_V1: usize = 61;

/// Exact number of no-mock, in-process public-API integration cases.
pub const BASE_COVERAGE_NO_MOCK_INTEGRATION_CASE_COUNT_V1: usize = 17;

/// Exact number of compile-fail contracts in the ratified base inventory.
pub const BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1: usize = 78;

/// Exact number of unit tests that protect this manifest contract itself.
pub const BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1: usize = 29;

/// Exact maximum number of named semantic numeric rows in either partition.
pub const BASE_COVERAGE_CLOSE_NUMERIC_EXPLICIT_MAX_V1: usize = 64;

/// Exact number of independently typed hard/soft budget axes.
pub const BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1: usize = 7;

/// Bounded canonical-frame grant for one maximum-shape Five Explicits value.
pub const BASE_COVERAGE_CLOSE_FIVE_EXPLICITS_FRAME_MAX_BYTES_V1: usize = 64 * 1024;

/// Bounded canonical-frame grant for one maximum-shape numeric profile.
pub const BASE_COVERAGE_CLOSE_NUMERIC_PROFILE_FRAME_MAX_BYTES_V1: usize = 16 * 1024;

/// Exact number of controlling AC53 category groups.
pub const BASE_COVERAGE_CLOSE_GROUP_COUNT_V1: usize = 9;

/// Exact number of controlling AC53 slash-separated facets.
pub const BASE_COVERAGE_CLOSE_FACET_COUNT_V1: usize = 22;

/// Exact number of registered AC53 Unsupported/inapplicability reasons.
pub const BASE_COVERAGE_CLOSE_REASON_COUNT_V1: usize = 5;

/// Maximum caller-declared extension cases admitted into one manifest.
pub const BASE_COVERAGE_EXTENSION_CASES_MAX_V1: usize = 1_024;

/// Maximum UTF-8 bytes in one stable coverage case ID.
pub const BASE_COVERAGE_CASE_ID_MAX_BYTES_V1: usize = 160;

/// Maximum UTF-8 bytes in one relative coverage source path.
pub const BASE_COVERAGE_SOURCE_PATH_MAX_BYTES_V1: usize = 240;

/// Closed source-coverage classes.
///
/// [`BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES`] partitions every
/// non-manifest Rust test through an explicit source declaration;
/// classification is never inferred from a test name. The original class codes
/// 1 through 9 remain stable, and the five new local evidence classes are
/// append-only codes 10 through 14. Compile-fail and manifest-contract classes
/// are separately frozen. Projection, logging, source-closure, and the final
/// three externally executed classes are exact extension lanes. A class label
/// declares the intended evidence shape only; it does not prove execution,
/// correctness, integration with an external system, or scientific authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageManifestClassV1 {
    /// Focused local behavior not primarily proving another evidence class.
    Unit = 1,
    /// The 78 source-frozen Rustdoc compile-fail contracts.
    CompileFailDoctest = 2,
    /// Unit tests protecting this manifest and join implementation.
    ManifestContract = 3,
    /// Real in-process projection journey rows declared by the projection owner.
    ProjectionE2e = 4,
    /// Runtime structured-log cases declared by the logging owner.
    RuntimeLogging = 5,
    /// Runtime exact source-closure cases declared by the governance owner.
    SourceClosure = 6,
    /// Canonical external end-to-end script cases.
    ExternalE2eScript = 7,
    /// External mutation-harness cases.
    ExternalMutation = 8,
    /// Live-tree governance and dependency-source checks.
    ExternalGovernance = 9,
    /// Exact zero/minimum/maximum/one-over and checked-arithmetic boundaries.
    Boundary = 10,
    /// Algebraic, exhaustive, round-trip, invariance, or metamorphic evidence.
    PropertyMetamorphic = 11,
    /// Literal catalogs, wire shapes, descriptor inventories, and closed tables.
    SchemaDescriptor = 12,
    /// One-field, malformed-input, missing/extra, collision, or stale mutants.
    Mutation = 13,
    /// Real public constructors and validators composed in process without mocks.
    NoMockIntegration = 14,
}

impl BaseCoverageManifestClassV1 {
    /// Every class in canonical order.
    pub const ALL: [Self; 14] = [
        Self::Unit,
        Self::CompileFailDoctest,
        Self::ManifestContract,
        Self::ProjectionE2e,
        Self::RuntimeLogging,
        Self::SourceClosure,
        Self::ExternalE2eScript,
        Self::ExternalMutation,
        Self::ExternalGovernance,
        Self::Boundary,
        Self::PropertyMetamorphic,
        Self::SchemaDescriptor,
        Self::Mutation,
        Self::NoMockIntegration,
    ];

    /// The six exact evidence classes partitioning the frozen Rust-test corpus.
    pub const RUST_TEST_EVIDENCE_CLASSES: [Self; 6] = [
        Self::Unit,
        Self::Boundary,
        Self::PropertyMetamorphic,
        Self::SchemaDescriptor,
        Self::Mutation,
        Self::NoMockIntegration,
    ];

    /// Externally executed classes. These cannot be satisfied by an in-process
    /// projection or by merely importing this crate.
    pub const EXTERNALLY_OWNED: [Self; 3] = [
        Self::ExternalE2eScript,
        Self::ExternalMutation,
        Self::ExternalGovernance,
    ];

    const fn code(self) -> u16 {
        self as u16
    }

    const fn stable_prefix(self) -> &'static str {
        match self {
            Self::Unit => "unit:",
            Self::Boundary => "boundary:",
            Self::PropertyMetamorphic => "property-metamorphic:",
            Self::SchemaDescriptor => "schema-descriptor:",
            Self::Mutation => "mutation:",
            Self::NoMockIntegration => "no-mock-integration:",
            Self::CompileFailDoctest => "compile-fail:",
            Self::ManifestContract => "manifest-contract:",
            Self::ProjectionE2e => "projection-e2e:",
            Self::RuntimeLogging => "runtime-logging:",
            Self::SourceClosure => "source-closure:",
            Self::ExternalE2eScript => "external-e2e:",
            Self::ExternalMutation => "external-mutation:",
            Self::ExternalGovernance => "external-governance:",
        }
    }

    const fn is_extension(self) -> bool {
        matches!(
            self,
            Self::ProjectionE2e
                | Self::RuntimeLogging
                | Self::SourceClosure
                | Self::ExternalE2eScript
                | Self::ExternalMutation
                | Self::ExternalGovernance
        )
    }
}

/// Exact authoritative partition of the frozen non-manifest Rust-test corpus.
pub const BASE_COVERAGE_RUST_TEST_CLASS_COUNTS_V1: [(BaseCoverageManifestClassV1, usize); 6] = [
    (
        BaseCoverageManifestClassV1::Unit,
        BASE_COVERAGE_UNIT_CLASS_CASE_COUNT_V1,
    ),
    (
        BaseCoverageManifestClassV1::Boundary,
        BASE_COVERAGE_BOUNDARY_CASE_COUNT_V1,
    ),
    (
        BaseCoverageManifestClassV1::PropertyMetamorphic,
        BASE_COVERAGE_PROPERTY_METAMORPHIC_CASE_COUNT_V1,
    ),
    (
        BaseCoverageManifestClassV1::SchemaDescriptor,
        BASE_COVERAGE_SCHEMA_DESCRIPTOR_CASE_COUNT_V1,
    ),
    (
        BaseCoverageManifestClassV1::Mutation,
        BASE_COVERAGE_MUTATION_CASE_COUNT_V1,
    ),
    (
        BaseCoverageManifestClassV1::NoMockIntegration,
        BASE_COVERAGE_NO_MOCK_INTEGRATION_CASE_COUNT_V1,
    ),
];

/// Independently declared, result-free coverage source case.
///
/// `source_path` names the workspace-relative source owner or designated
/// external harness for the case. A designated downstream harness need not
/// exist in the current leaf, and this declaration does not claim that it ran.
/// Construction validates only bounded stable-ID and relative-path grammar;
/// admission into a manifest additionally validates exact class ownership,
/// ordering, and uniqueness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCaseDeclarationV1 {
    class: BaseCoverageManifestClassV1,
    id: Box<str>,
    source_path: Box<str>,
}

impl BaseCoverageCaseDeclarationV1 {
    /// Construct one result-free case declaration.
    pub fn new(
        class: BaseCoverageManifestClassV1,
        id: impl Into<String>,
        source_path: impl Into<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        let id = id.into();
        let source_path = source_path.into();
        validate_case_id(class, &id)?;
        validate_source_path(&source_path)?;
        Ok(Self {
            class,
            id: id.into_boxed_str(),
            source_path: source_path.into_boxed_str(),
        })
    }

    /// Closed coverage class and execution owner.
    #[must_use]
    pub const fn class(&self) -> BaseCoverageManifestClassV1 {
        self.class
    }

    /// Globally unique stable source-case ID.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Exact workspace-relative source owner or designated harness path.
    ///
    /// This is a result-free ownership/execution mapping, not an existence or
    /// execution claim.
    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
    }
}

/// One admitted immutable manifest case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageManifestCaseV1 {
    ordinal: u32,
    declaration: BaseCoverageCaseDeclarationV1,
}

impl BaseCoverageManifestCaseV1 {
    /// One-based global ordinal in the exact manifest.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Closed coverage class and execution owner.
    #[must_use]
    pub const fn class(&self) -> BaseCoverageManifestClassV1 {
        self.declaration.class
    }

    /// Globally unique stable source-case ID.
    #[must_use]
    pub fn id(&self) -> &str {
        self.declaration.id()
    }

    /// Exact workspace-relative source owner or designated harness path.
    ///
    /// This is a result-free ownership/execution mapping, not an existence or
    /// execution claim.
    #[must_use]
    pub fn source_path(&self) -> &str {
        self.declaration.source_path()
    }

    /// Original result-free declaration.
    #[must_use]
    pub const fn declaration(&self) -> &BaseCoverageCaseDeclarationV1 {
        &self.declaration
    }
}

/// Sole source-authoritative, immutable, result-free AC38 coverage manifest.
///
/// Projection-local compatibility inventories are auxiliary data and cannot
/// replace, widen, or prove this exact manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageManifestV1 {
    cases: Box<[BaseCoverageManifestCaseV1]>,
    root: ContentHash,
}

impl BaseCoverageManifestV1 {
    /// Construct the frozen 217-Rust-test, 78-compile-fail, and
    /// manifest-contract source inventory without extension rows.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let declarations = frozen_base_declarations()?;
        manifest_from_declarations(declarations)
    }

    /// Reconstruct the frozen base from a caller-presented exact declaration
    /// sequence.
    ///
    /// Missing, extra, duplicate, reordered, or semantically changed
    /// declarations refuse before a manifest root is returned.
    pub fn reconstruct_exact_base(
        presented: &[BaseCoverageCaseDeclarationV1],
    ) -> Result<Self, ConstructionErrorV2> {
        let expected = frozen_base_declarations()?;
        validate_exact_declaration_sequence(&expected, presented)?;
        manifest_from_declarations(expected)
    }

    /// Construct the frozen base plus independently declared exact extension
    /// cases.
    ///
    /// Extensions must use only extension classes and must already be ordered
    /// by `(class code, stable ID, source path)`. This method never derives
    /// expected extension IDs from runtime results. The declaring owner must
    /// retain an independent literal oracle for the supplied declarations.
    pub fn with_exact_extensions(
        extensions: &[BaseCoverageCaseDeclarationV1],
    ) -> Result<Self, ConstructionErrorV2> {
        validate_extensions(extensions)?;
        let mut declarations = frozen_base_declarations()?;
        declarations.extend_from_slice(extensions);
        manifest_from_declarations(declarations)
    }

    /// Manifest cases in exact canonical order.
    #[must_use]
    pub fn cases(&self) -> &[BaseCoverageManifestCaseV1] {
        &self.cases
    }

    /// Number of source cases in one closed class.
    #[must_use]
    pub fn case_count(&self, class: BaseCoverageManifestClassV1) -> usize {
        self.cases
            .iter()
            .filter(|case| case.class() == class)
            .count()
    }

    /// Find a manifest case by its globally unique stable ID.
    #[must_use]
    pub fn case(&self, id: &str) -> Option<&BaseCoverageManifestCaseV1> {
        self.cases.iter().find(|case| case.id() == id)
    }

    /// Domain-separated root of the exact result-free manifest.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Select an exact executable subset in manifest order.
    ///
    /// Skipping cases is allowed and represents caller scope, not success.
    /// Unknown, duplicate, or reordered IDs refuse. An empty subset is valid
    /// and produces an exact zero-result accounting report.
    pub fn select_exact(
        &self,
        selected_ids: &[&str],
    ) -> Result<BaseCoverageExecutableSubsetV1, ConstructionErrorV2> {
        let positions = self
            .cases
            .iter()
            .enumerate()
            .map(|(index, case)| (case.id(), index))
            .collect::<BTreeMap<_, _>>();
        let mut seen = BTreeSet::new();
        let mut previous_position = None;
        let mut ids = Vec::with_capacity(selected_ids.len());
        for (ordinal, id) in selected_ids.iter().copied().enumerate() {
            validate_untyped_case_id(id)?;
            if !seen.insert(id) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.selection.source_case_id",
                    "one occurrence of every selected source-case ID",
                    id,
                ));
            }
            let Some(position) = positions.get(id).copied() else {
                return Err(refusal(
                    ConstructionErrorKindV2::UnknownCode,
                    "coverage.selection.source_case_id",
                    "an ID mapped by the exact manifest",
                    id,
                ));
            };
            if previous_position.is_some_and(|previous| position <= previous) {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.selection.source_case_id",
                    "strict manifest order",
                    format_args!("{ordinal}:{id}"),
                ));
            }
            previous_position = Some(position);
            ids.push(id.to_owned().into_boxed_str());
        }
        let root = selection_root(self.root, &ids)?;
        Ok(BaseCoverageExecutableSubsetV1 {
            manifest_root: self.root,
            source_case_ids: ids.into_boxed_slice(),
            root,
        })
    }
}

/// Manifest-bound caller selection for one executable run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageExecutableSubsetV1 {
    manifest_root: ContentHash,
    source_case_ids: Box<[Box<str>]>,
    root: ContentHash,
}

impl BaseCoverageExecutableSubsetV1 {
    /// Manifest root against which this subset was selected.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.manifest_root
    }

    /// Exact selected source-case IDs in manifest order.
    #[must_use]
    pub fn source_case_ids(&self) -> &[Box<str>] {
        &self.source_case_ids
    }

    /// Domain-separated root binding the manifest and exact selection.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Closed presented outcome partition.
///
/// These values are accounting labels supplied by the execution owner. This
/// module does not independently verify that the underlying evidence warrants
/// the label.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoveragePresentedOutcomeV1 {
    /// A positive case matched its independent oracle.
    PositiveMatched = 1,
    /// A deliberate invalid case produced its exact expected refusal.
    ExpectedRefusalMatched = 2,
    /// An explicitly unsupported case matched its exact unsupported oracle.
    ExpectedUnsupportedMatched = 3,
    /// At least one expected semantic cell did not match.
    UnexpectedMismatch = 4,
}

impl BaseCoveragePresentedOutcomeV1 {
    const fn code(self) -> u16 {
        self as u16
    }
}

/// One caller-presented result bound to a manifest and source-case ID.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoveragePresentedResultV1 {
    manifest_root: ContentHash,
    source_case_id: Box<str>,
    outcome: BaseCoveragePresentedOutcomeV1,
    evidence_root: ContentHash,
    root: ContentHash,
}

impl BaseCoveragePresentedResultV1 {
    /// Construct one bounded, manifest-bound presented result.
    pub fn new(
        manifest_root: ContentHash,
        source_case_id: impl Into<String>,
        outcome: BaseCoveragePresentedOutcomeV1,
        evidence_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_case_id = source_case_id.into();
        validate_untyped_case_id(&source_case_id)?;
        let root = presented_result_root(manifest_root, &source_case_id, outcome, evidence_root)?;
        Ok(Self {
            manifest_root,
            source_case_id: source_case_id.into_boxed_str(),
            outcome,
            evidence_root,
            root,
        })
    }

    /// Manifest root presented by the execution owner.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.manifest_root
    }

    /// Stable source-case ID presented by the execution owner.
    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    /// Closed presented outcome.
    #[must_use]
    pub const fn outcome(&self) -> BaseCoveragePresentedOutcomeV1 {
        self.outcome
    }

    /// Opaque evidence root. Its meaning belongs to the execution owner.
    #[must_use]
    pub const fn evidence_root(&self) -> ContentHash {
        self.evidence_root
    }

    /// Domain-separated root binding every presented field.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Exact manifest/subset/result accounting report.
///
/// A green report means there was exactly one manifest-bound presented result
/// for every selected ID, in order, and none was labelled
/// [`BaseCoveragePresentedOutcomeV1::UnexpectedMismatch`]. It does not prove
/// that a test binary ran or that an evidence root is trustworthy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCheckedReportV1 {
    manifest_root: ContentHash,
    selection_root: ContentHash,
    results: Box<[BaseCoveragePresentedResultV1]>,
    positive_matched: u32,
    expected_refusals_matched: u32,
    expected_unsupported_matched: u32,
    unexpected_mismatches: u32,
    root: ContentHash,
}

impl BaseCoverageCheckedReportV1 {
    /// Exact-join caller-presented results to one caller-selected subset.
    ///
    /// The reconstructor rejects stale manifest roots, globally unmapped IDs,
    /// IDs outside the selected subset, multiply reported IDs, missing
    /// results, extra results, and reordered results before returning a report.
    pub fn reconstruct(
        manifest: &BaseCoverageManifestV1,
        selection: &BaseCoverageExecutableSubsetV1,
        presented: &[BaseCoveragePresentedResultV1],
    ) -> Result<Self, ConstructionErrorV2> {
        validate_selection_against_manifest(manifest, selection)?;
        validate_presented_results(manifest, selection, presented)?;

        let mut positive_matched = 0_u32;
        let mut expected_refusals_matched = 0_u32;
        let mut expected_unsupported_matched = 0_u32;
        let mut unexpected_mismatches = 0_u32;
        for result in presented {
            let counter = match result.outcome {
                BaseCoveragePresentedOutcomeV1::PositiveMatched => &mut positive_matched,
                BaseCoveragePresentedOutcomeV1::ExpectedRefusalMatched => {
                    &mut expected_refusals_matched
                }
                BaseCoveragePresentedOutcomeV1::ExpectedUnsupportedMatched => {
                    &mut expected_unsupported_matched
                }
                BaseCoveragePresentedOutcomeV1::UnexpectedMismatch => &mut unexpected_mismatches,
            };
            *counter = counter.checked_add(1).ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "coverage.report.outcome_count",
                    "checked u32 outcome counts",
                    presented.len(),
                )
            })?;
        }
        let root = checked_report_root(selection.root, presented)?;
        Ok(Self {
            manifest_root: manifest.root,
            selection_root: selection.root,
            results: presented.to_vec().into_boxed_slice(),
            positive_matched,
            expected_refusals_matched,
            expected_unsupported_matched,
            unexpected_mismatches,
            root,
        })
    }

    /// Exact manifest root used by the join.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.manifest_root
    }

    /// Exact selection root used by the join.
    #[must_use]
    pub const fn selection_root(&self) -> ContentHash {
        self.selection_root
    }

    /// Presented results in exact selected order.
    #[must_use]
    pub fn results(&self) -> &[BaseCoveragePresentedResultV1] {
        &self.results
    }

    /// Positive cases labelled as independently matched.
    #[must_use]
    pub const fn positive_matched(&self) -> u32 {
        self.positive_matched
    }

    /// Expected refusals labelled as exactly matched.
    #[must_use]
    pub const fn expected_refusals_matched(&self) -> u32 {
        self.expected_refusals_matched
    }

    /// Expected unsupported cases labelled as exactly matched.
    #[must_use]
    pub const fn expected_unsupported_matched(&self) -> u32 {
        self.expected_unsupported_matched
    }

    /// Cases labelled as unexpected semantic mismatches.
    #[must_use]
    pub const fn unexpected_mismatches(&self) -> u32 {
        self.unexpected_mismatches
    }

    /// Whether no result was labelled as an unexpected mismatch.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        self.unexpected_mismatches == 0
    }

    /// Domain-separated root of the exact joined report.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// The nine preserved AC53 close-manifest category groups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseGroupV1 {
    /// Literal, unit, and boundary evidence.
    LiteralUnitBoundary = 1,
    /// Property and metamorphic evidence.
    PropertyMetamorphic = 2,
    /// State, typestate, model, and race evidence.
    StateTypestateModelRace = 3,
    /// Mutation and fuzz evidence.
    MutationFuzz = 4,
    /// API, compile-fail, and trait evidence.
    ApiCompileFailTrait = 5,
    /// Fault, resource, and cancellation evidence.
    FaultResourceCancellation = 6,
    /// Release-built execution or immutable E2E contribution evidence.
    ReleaseBuiltOrImmutableE2e = 7,
    /// Detailed deterministic logging and redaction evidence.
    DetailedLoggingRedaction = 8,
    /// Exact source-closure evidence.
    SourceClosure = 9,
}

impl BaseCoverageCloseGroupV1 {
    /// Every group in exact code order.
    pub const ALL: [Self; BASE_COVERAGE_CLOSE_GROUP_COUNT_V1] = [
        Self::LiteralUnitBoundary,
        Self::PropertyMetamorphic,
        Self::StateTypestateModelRace,
        Self::MutationFuzz,
        Self::ApiCompileFailTrait,
        Self::FaultResourceCancellation,
        Self::ReleaseBuiltOrImmutableE2e,
        Self::DetailedLoggingRedaction,
        Self::SourceClosure,
    ];

    /// Stable numeric code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable closed name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::LiteralUnitBoundary => "literal-unit-boundary",
            Self::PropertyMetamorphic => "property-metamorphic",
            Self::StateTypestateModelRace => "state-typestate-model-race",
            Self::MutationFuzz => "mutation-fuzz",
            Self::ApiCompileFailTrait => "api-compile-fail-trait",
            Self::FaultResourceCancellation => "fault-resource-cancellation",
            Self::ReleaseBuiltOrImmutableE2e => "release-built-or-immutable-e2e",
            Self::DetailedLoggingRedaction => "detailed-logging-redaction",
            Self::SourceClosure => "source-closure",
        }
    }
}

/// The twenty-two preserved AC53 slash-separated facets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseFacetV1 {
    Literal = 1,
    Unit = 2,
    Boundary = 3,
    Property = 4,
    Metamorphic = 5,
    State = 6,
    Typestate = 7,
    Model = 8,
    Race = 9,
    Mutation = 10,
    Fuzz = 11,
    Api = 12,
    CompileFail = 13,
    Trait = 14,
    Fault = 15,
    Resource = 16,
    Cancellation = 17,
    ReleaseBuiltNoMockE2e = 18,
    ImmutableE2eContribution = 19,
    DetailedDeterministicLogging = 20,
    Redaction = 21,
    SourceClosure = 22,
}

impl BaseCoverageCloseFacetV1 {
    /// Every facet in exact code order.
    pub const ALL: [Self; BASE_COVERAGE_CLOSE_FACET_COUNT_V1] = [
        Self::Literal,
        Self::Unit,
        Self::Boundary,
        Self::Property,
        Self::Metamorphic,
        Self::State,
        Self::Typestate,
        Self::Model,
        Self::Race,
        Self::Mutation,
        Self::Fuzz,
        Self::Api,
        Self::CompileFail,
        Self::Trait,
        Self::Fault,
        Self::Resource,
        Self::Cancellation,
        Self::ReleaseBuiltNoMockE2e,
        Self::ImmutableE2eContribution,
        Self::DetailedDeterministicLogging,
        Self::Redaction,
        Self::SourceClosure,
    ];

    /// Stable numeric code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable closed name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Literal => "literal",
            Self::Unit => "unit",
            Self::Boundary => "boundary",
            Self::Property => "property",
            Self::Metamorphic => "metamorphic",
            Self::State => "state",
            Self::Typestate => "typestate",
            Self::Model => "model",
            Self::Race => "race",
            Self::Mutation => "mutation",
            Self::Fuzz => "fuzz",
            Self::Api => "api",
            Self::CompileFail => "compile-fail",
            Self::Trait => "trait",
            Self::Fault => "fault",
            Self::Resource => "resource",
            Self::Cancellation => "cancellation",
            Self::ReleaseBuiltNoMockE2e => "release-built-no-mock-e2e",
            Self::ImmutableE2eContribution => "immutable-e2e-contribution",
            Self::DetailedDeterministicLogging => "detailed-deterministic-logging",
            Self::Redaction => "redaction",
            Self::SourceClosure => "source-closure",
        }
    }

    /// Sole owning category group.
    #[must_use]
    pub const fn group(self) -> BaseCoverageCloseGroupV1 {
        match self {
            Self::Literal | Self::Unit | Self::Boundary => {
                BaseCoverageCloseGroupV1::LiteralUnitBoundary
            }
            Self::Property | Self::Metamorphic => BaseCoverageCloseGroupV1::PropertyMetamorphic,
            Self::State | Self::Typestate | Self::Model | Self::Race => {
                BaseCoverageCloseGroupV1::StateTypestateModelRace
            }
            Self::Mutation | Self::Fuzz => BaseCoverageCloseGroupV1::MutationFuzz,
            Self::Api | Self::CompileFail | Self::Trait => {
                BaseCoverageCloseGroupV1::ApiCompileFailTrait
            }
            Self::Fault | Self::Resource | Self::Cancellation => {
                BaseCoverageCloseGroupV1::FaultResourceCancellation
            }
            Self::ReleaseBuiltNoMockE2e | Self::ImmutableE2eContribution => {
                BaseCoverageCloseGroupV1::ReleaseBuiltOrImmutableE2e
            }
            Self::DetailedDeterministicLogging | Self::Redaction => {
                BaseCoverageCloseGroupV1::DetailedLoggingRedaction
            }
            Self::SourceClosure => BaseCoverageCloseGroupV1::SourceClosure,
        }
    }
}

/// Exact execution ownership for one AC53 close cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseExecutionScopeV1 {
    CrateTest = 1,
    CompileFailDoctest = 2,
    InProcessProjection = 3,
    ImmutableDownstreamContribution = 4,
    FacetApplicabilityDeclaration = 5,
}

impl BaseCoverageCloseExecutionScopeV1 {
    pub const ALL: [Self; 5] = [
        Self::CrateTest,
        Self::CompileFailDoctest,
        Self::InProcessProjection,
        Self::ImmutableDownstreamContribution,
        Self::FacetApplicabilityDeclaration,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::CrateTest => "crate-test",
            Self::CompileFailDoctest => "compile-fail-doctest",
            Self::InProcessProjection => "in-process-projection",
            Self::ImmutableDownstreamContribution => "immutable-downstream-contribution",
            Self::FacetApplicabilityDeclaration => "facet-applicability-declaration",
        }
    }
}

/// Exact expected or observed close-cell decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseDecisionV1 {
    Accept = 1,
    Refuse = 2,
    Fail = 3,
    Unsupported = 4,
    Inapplicable = 5,
}

impl BaseCoverageCloseDecisionV1 {
    pub const ALL: [Self; 5] = [
        Self::Accept,
        Self::Refuse,
        Self::Fail,
        Self::Unsupported,
        Self::Inapplicable,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Accept => "accept",
            Self::Refuse => "refuse",
            Self::Fail => "fail",
            Self::Unsupported => "unsupported",
            Self::Inapplicable => "inapplicable",
        }
    }
}

/// Exact AC53 accounting partition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageClosePartitionV1 {
    Positive = 1,
    ExpectedRefusal = 2,
    ExpectedFailure = 3,
    Mutation = 4,
    Unsupported = 5,
    Inapplicable = 6,
}

impl BaseCoverageClosePartitionV1 {
    pub const ALL: [Self; 6] = [
        Self::Positive,
        Self::ExpectedRefusal,
        Self::ExpectedFailure,
        Self::Mutation,
        Self::Unsupported,
        Self::Inapplicable,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

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

    const fn is_adversarial(self) -> bool {
        matches!(
            self,
            Self::ExpectedRefusal | Self::ExpectedFailure | Self::Mutation
        )
    }
}

/// Closed typed reasons for Unsupported and genuinely inapplicable cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseReasonCodeV1 {
    RaceNotApplicablePureSingleThreadedValidator = 1,
    TraitNotApplicableNoPublicTraitContract = 2,
    CancellationNotApplicablePureBoundedValidator = 3,
    ReleaseExecutionDownstreamOwned = 4,
    WindowsNonasciiAliasLocallyUnadjudicable = 5,
}

impl BaseCoverageCloseReasonCodeV1 {
    pub const ALL: [Self; BASE_COVERAGE_CLOSE_REASON_COUNT_V1] = [
        Self::RaceNotApplicablePureSingleThreadedValidator,
        Self::TraitNotApplicableNoPublicTraitContract,
        Self::CancellationNotApplicablePureBoundedValidator,
        Self::ReleaseExecutionDownstreamOwned,
        Self::WindowsNonasciiAliasLocallyUnadjudicable,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn descriptor(self) -> &'static BaseCoverageCloseReasonDescriptorV1 {
        &BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1[self as usize - 1]
    }
}

/// One immutable typed reason descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCoverageCloseReasonDescriptorV1 {
    code: BaseCoverageCloseReasonCodeV1,
    name: &'static str,
    owner: &'static str,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    prerequisite: &'static str,
    no_claim: &'static str,
}

impl BaseCoverageCloseReasonDescriptorV1 {
    #[must_use]
    pub const fn code(self) -> BaseCoverageCloseReasonCodeV1 {
        self.code
    }

    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    #[must_use]
    pub const fn owner(self) -> &'static str {
        self.owner
    }

    #[must_use]
    pub const fn execution_scope(self) -> BaseCoverageCloseExecutionScopeV1 {
        self.execution_scope
    }

    #[must_use]
    pub const fn prerequisite(self) -> &'static str {
        self.prerequisite
    }

    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        self.no_claim
    }
}

/// Exact reason registry in stable code order.
pub const BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1: [BaseCoverageCloseReasonDescriptorV1;
    BASE_COVERAGE_CLOSE_REASON_COUNT_V1] = [
    BaseCoverageCloseReasonDescriptorV1 {
        code: BaseCoverageCloseReasonCodeV1::RaceNotApplicablePureSingleThreadedValidator,
        name: "race-not-applicable-pure-single-threaded-validator",
        owner: "fs-evidence-runner/coverage",
        execution_scope: BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration,
        prerequisite: "a-concurrent-stateful-execution-owner",
        no_claim: "no-concurrent-race-execution-claim",
    },
    BaseCoverageCloseReasonDescriptorV1 {
        code: BaseCoverageCloseReasonCodeV1::TraitNotApplicableNoPublicTraitContract,
        name: "trait-not-applicable-no-public-trait-contract",
        owner: "fs-evidence-runner/coverage",
        execution_scope: BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration,
        prerequisite: "a-public-trait-contract-owned-by-this-leaf",
        no_claim: "no-public-trait-conformance-claim",
    },
    BaseCoverageCloseReasonDescriptorV1 {
        code: BaseCoverageCloseReasonCodeV1::CancellationNotApplicablePureBoundedValidator,
        name: "cancellation-not-applicable-pure-bounded-validator",
        owner: "fs-evidence-runner/coverage",
        execution_scope: BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration,
        prerequisite: "an-asynchronous-cancellable-effect-owner",
        no_claim: "no-cancellation-drain-or-resource-return-claim",
    },
    BaseCoverageCloseReasonDescriptorV1 {
        code: BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned,
        name: "release-execution-downstream-owned",
        owner: "runner-v2-downstream-script-owners",
        execution_scope: BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution,
        prerequisite: "release-built-no-mock-downstream-script-execution",
        no_claim: "immutable-contribution-is-not-release-execution-proof",
    },
    BaseCoverageCloseReasonDescriptorV1 {
        code: BaseCoverageCloseReasonCodeV1::WindowsNonasciiAliasLocallyUnadjudicable,
        name: "windows-nonascii-alias-locally-unadjudicable",
        owner: "platform-path-collision-owner",
        execution_scope: BaseCoverageCloseExecutionScopeV1::InProcessProjection,
        prerequisite: "a-proved-windows-unicode-collision-key",
        no_claim: "no-windows-nonascii-alias-equivalence-claim",
    },
];

/// The reason descriptors in exact code order.
#[must_use]
pub const fn base_coverage_close_reason_descriptors_v1()
-> &'static [BaseCoverageCloseReasonDescriptorV1] {
    &BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1
}

/// Closed AC57 disposition for the runtime-observation side of one close cell.
///
/// This is deliberately distinct from semantic result status. A matched
/// positive, expected refusal, expected failure, mutation, or executed
/// Unsupported result can all be `Observed` when their runtime actuals are
/// complete.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RuntimeObservationDispositionV1 {
    /// An execution-owned attempt completed with every required actual.
    Observed = 1,
    /// An attempt exists, but complete actuals could not be established.
    NotObserved = 2,
    /// An immutable contribution awaits its designated execution owner.
    Deferred = 3,
    /// A source-registered applicability declaration excludes execution.
    Inapplicable = 4,
}

impl RuntimeObservationDispositionV1 {
    /// Every admitted disposition in exact wire-code order.
    pub const ALL: [Self; 4] = [
        Self::Observed,
        Self::NotObserved,
        Self::Deferred,
        Self::Inapplicable,
    ];

    /// Exact unsigned 16-bit wire code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Observed => "observed",
            Self::NotObserved => "not-observed",
            Self::Deferred => "deferred",
            Self::Inapplicable => "inapplicable",
        }
    }

    /// Parse one exact nonzero closed code. Zero and unknown codes refuse.
    pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::Observed),
            2 => Ok(Self::NotObserved),
            3 => Ok(Self::Deferred),
            4 => Ok(Self::Inapplicable),
            _ => Err(refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close.runtime_observation.disposition",
                "one exact AC57 runtime-observation disposition code in 1..=4",
                code,
            )),
        }
    }

    /// Nominal identity for this exact disposition value.
    pub fn root(
        self,
    ) -> Result<BaseCoverageCloseRuntimeObservationDispositionRootV1, ConstructionErrorV2> {
        close_runtime_observation_disposition_root(self)
    }
}

/// Exact AC57 reason for an attempted cell without complete observations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum NotObservedReasonV1 {
    /// The execution attempt failed before all actuals became complete.
    ExecutionFailedBeforeCompleteness = 1,
    /// The observation channel failed before all actuals became complete.
    ObservationChannelFailedBeforeCompleteness = 2,
    /// A cell expected to run was unexpectedly unstarted or skipped.
    UnexpectedUnstartedOrSkipped = 3,
}

impl NotObservedReasonV1 {
    /// Every admitted reason in exact wire-code order.
    pub const ALL: [Self; 3] = [
        Self::ExecutionFailedBeforeCompleteness,
        Self::ObservationChannelFailedBeforeCompleteness,
        Self::UnexpectedUnstartedOrSkipped,
    ];

    /// Exact unsigned 16-bit wire code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact immutable descriptor row.
    #[must_use]
    pub const fn descriptor(self) -> &'static NotObservedReasonDescriptorV1 {
        &NOT_OBSERVED_REASON_DESCRIPTORS_V1[self as usize - 1]
    }

    /// Parse one exact nonzero closed code. Zero and unknown codes refuse.
    pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::ExecutionFailedBeforeCompleteness),
            2 => Ok(Self::ObservationChannelFailedBeforeCompleteness),
            3 => Ok(Self::UnexpectedUnstartedOrSkipped),
            _ => Err(refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close.runtime_observation.not_observed_reason",
                "one exact AC57 NotObserved reason code in 1..=3",
                code,
            )),
        }
    }
}

/// One immutable source-owned NotObserved reason descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NotObservedReasonDescriptorV1 {
    reason: NotObservedReasonV1,
    name: &'static str,
    owner: &'static str,
    scope: &'static str,
    prerequisite: &'static str,
    diagnostic: DiagnosticCodeV2,
    no_claim: &'static str,
}

impl NotObservedReasonDescriptorV1 {
    /// Exact reason code.
    #[must_use]
    pub const fn reason(self) -> NotObservedReasonV1 {
        self.reason
    }

    /// Exact stable reason name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// Exact source owner.
    #[must_use]
    pub const fn owner(self) -> &'static str {
        self.owner
    }

    /// Exact closed semantic scope.
    #[must_use]
    pub const fn scope(self) -> &'static str {
        self.scope
    }

    /// Exact prerequisite for resolving this incomplete state.
    #[must_use]
    pub const fn prerequisite(self) -> &'static str {
        self.prerequisite
    }

    /// Exact actionable base diagnostic.
    #[must_use]
    pub const fn diagnostic(self) -> DiagnosticCodeV2 {
        self.diagnostic
    }

    /// Exact no-complete-observation boundary.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        self.no_claim
    }
}

/// Exact source-owned NotObserved reason rows in code order.
pub const NOT_OBSERVED_REASON_DESCRIPTORS_V1: [NotObservedReasonDescriptorV1; 3] = [
    NotObservedReasonDescriptorV1 {
        reason: NotObservedReasonV1::ExecutionFailedBeforeCompleteness,
        name: "execution-failed-before-completeness",
        owner: "fs-evidence-runner.coverage",
        scope: "execution-owned-attempt",
        prerequisite: "a-complete-terminal-attempt-and-runtime-observation",
        diagnostic: DiagnosticCodeV2::RunnerInternalError,
        no_claim: "not-observed-execution-failure-proves-no-complete-observation-or-success",
    },
    NotObservedReasonDescriptorV1 {
        reason: NotObservedReasonV1::ObservationChannelFailedBeforeCompleteness,
        name: "observation-channel-failed-before-completeness",
        owner: "fs-evidence-runner.coverage",
        scope: "execution-owned-observation-channel",
        prerequisite: "a-complete-redacted-runtime-observation-channel",
        diagnostic: DiagnosticCodeV2::RunnerNoData,
        no_claim: "not-observed-channel-failure-proves-no-complete-observation-or-success",
    },
    NotObservedReasonDescriptorV1 {
        reason: NotObservedReasonV1::UnexpectedUnstartedOrSkipped,
        name: "unexpected-unstarted-or-skipped",
        owner: "fs-evidence-runner.coverage",
        scope: "execution-owned-dispatch",
        prerequisite: "a-complete-execution-owned-attempt",
        diagnostic: DiagnosticCodeV2::RunnerNotRun,
        no_claim: "not-observed-unstarted-or-skipped-proves-no-execution-observation-or-success",
    },
];

/// Exact source-owned NotObserved reason registry and its nominal identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NotObservedReasonRegistryV1 {
    root: BaseCoverageCloseNotObservedReasonRegistryRootV1,
}

impl NotObservedReasonRegistryV1 {
    /// Reconstruct the sole frozen registry.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        Ok(Self {
            root: close_not_observed_reason_registry_root(&NOT_OBSERVED_REASON_DESCRIPTORS_V1)?,
        })
    }

    /// Exact rows in code order.
    #[must_use]
    pub const fn descriptors(&self) -> &'static [NotObservedReasonDescriptorV1; 3] {
        &NOT_OBSERVED_REASON_DESCRIPTORS_V1
    }

    /// Look up one exact registered reason.
    #[must_use]
    pub const fn descriptor(
        &self,
        reason: NotObservedReasonV1,
    ) -> &'static NotObservedReasonDescriptorV1 {
        reason.descriptor()
    }

    /// Nominal root of the complete exact registry.
    #[must_use]
    pub const fn root(self) -> BaseCoverageCloseNotObservedReasonRegistryRootV1 {
        self.root
    }
}

/// Exact AC57 deferred reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum DeferredReasonV1 {
    /// Immutable result-free contribution awaiting its designated owner.
    ImmutableContributionAwaitsDesignatedReleaseOwner = 1,
}

impl DeferredReasonV1 {
    /// Every admitted reason in exact wire-code order.
    pub const ALL: [Self; 1] = [Self::ImmutableContributionAwaitsDesignatedReleaseOwner];

    /// Exact unsigned 16-bit wire code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact immutable descriptor row.
    #[must_use]
    pub const fn descriptor(self) -> &'static DeferredReasonDescriptorV1 {
        &DEFERRED_REASON_DESCRIPTORS_V1[self as usize - 1]
    }

    /// Parse the sole exact nonzero code. Zero and unknown codes refuse.
    pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::ImmutableContributionAwaitsDesignatedReleaseOwner),
            _ => Err(refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close.runtime_observation.deferred_reason",
                "the exact AC57 Deferred reason code 1",
                code,
            )),
        }
    }
}

/// One immutable source-owned Deferred reason descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeferredReasonDescriptorV1 {
    reason: DeferredReasonV1,
    name: &'static str,
    owner: &'static str,
    scope: &'static str,
    prerequisite: &'static str,
    diagnostic: DiagnosticCodeV2,
    no_claim: &'static str,
}

impl DeferredReasonDescriptorV1 {
    /// Exact reason code.
    #[must_use]
    pub const fn reason(self) -> DeferredReasonV1 {
        self.reason
    }

    /// Exact stable reason name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// Exact source owner.
    #[must_use]
    pub const fn owner(self) -> &'static str {
        self.owner
    }

    /// Exact closed semantic scope.
    #[must_use]
    pub const fn scope(self) -> &'static str {
        self.scope
    }

    /// Exact prerequisite for replacing the Deferred envelope.
    #[must_use]
    pub const fn prerequisite(self) -> &'static str {
        self.prerequisite
    }

    /// Exact actionable base diagnostic.
    #[must_use]
    pub const fn diagnostic(self) -> DiagnosticCodeV2 {
        self.diagnostic
    }

    /// Exact no-execution-proof boundary.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        self.no_claim
    }
}

/// Exact source-owned Deferred reason rows in code order.
pub const DEFERRED_REASON_DESCRIPTORS_V1: [DeferredReasonDescriptorV1; 1] =
    [DeferredReasonDescriptorV1 {
        reason: DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner,
        name: "immutable-contribution-awaits-designated-release-owner",
        owner: "fs-evidence-runner.coverage",
        scope: "immutable-downstream-contribution",
        prerequisite: "release-built-no-mock-designated-owner-execution",
        diagnostic: DiagnosticCodeV2::RunnerNotRun,
        no_claim: "deferred-contribution-proves-no-designated-owner-execution-or-success",
    }];

/// Exact source-owned Deferred reason registry and its nominal identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeferredReasonRegistryV1 {
    root: BaseCoverageCloseDeferredReasonRegistryRootV1,
}

impl DeferredReasonRegistryV1 {
    /// Reconstruct the sole frozen registry.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        Ok(Self {
            root: close_deferred_reason_registry_root(&DEFERRED_REASON_DESCRIPTORS_V1)?,
        })
    }

    /// Exact rows in code order.
    #[must_use]
    pub const fn descriptors(&self) -> &'static [DeferredReasonDescriptorV1; 1] {
        &DEFERRED_REASON_DESCRIPTORS_V1
    }

    /// Look up the exact registered reason.
    #[must_use]
    pub const fn descriptor(
        &self,
        reason: DeferredReasonV1,
    ) -> &'static DeferredReasonDescriptorV1 {
        reason.descriptor()
    }

    /// Nominal root of the complete exact registry.
    #[must_use]
    pub const fn root(self) -> BaseCoverageCloseDeferredReasonRegistryRootV1 {
        self.root
    }
}

/// Closed AC58 classification of one canonical-schema impact row.
///
/// These codes describe canonical-schema migration. They do not replace or
/// reinterpret Runner API generation, Runner wire version, or Runner wire
/// predecessor policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalSchemaImpactDispositionV1 {
    /// A genuinely new V1 component with no canonical predecessor.
    NewV1NoPredecessor = 1,
    /// An existing V1 frame whose bytes and meaning remain unchanged.
    UnchangedV1 = 2,
    /// An authoritative V1 frame migrated coherently to a distinct V2 frame.
    MigratedV1ToV2 = 3,
    /// A legacy V1 frame retained only as decode/compatibility evidence.
    DecodeOnlyLegacyV1 = 4,
    /// A legacy V1 frame explicitly rejected by the authoritative path.
    RetiredV1 = 5,
    /// A source-owned item with no canonical frame to version or migrate.
    InapplicableNoCanonicalFrame = 6,
}

impl CanonicalSchemaImpactDispositionV1 {
    /// Every admitted disposition in exact wire-code order.
    pub const ALL: [Self; 6] = [
        Self::NewV1NoPredecessor,
        Self::UnchangedV1,
        Self::MigratedV1ToV2,
        Self::DecodeOnlyLegacyV1,
        Self::RetiredV1,
        Self::InapplicableNoCanonicalFrame,
    ];

    /// Exact unsigned 16-bit code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::NewV1NoPredecessor => "new-v1-no-predecessor",
            Self::UnchangedV1 => "unchanged-v1",
            Self::MigratedV1ToV2 => "migrated-v1-to-v2",
            Self::DecodeOnlyLegacyV1 => "decode-only-legacy-v1",
            Self::RetiredV1 => "retired-v1",
            Self::InapplicableNoCanonicalFrame => "inapplicable-no-canonical-frame",
        }
    }

    /// Parse one exact nonzero closed code. Zero and unknown codes refuse.
    pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::NewV1NoPredecessor),
            2 => Ok(Self::UnchangedV1),
            3 => Ok(Self::MigratedV1ToV2),
            4 => Ok(Self::DecodeOnlyLegacyV1),
            5 => Ok(Self::RetiredV1),
            6 => Ok(Self::InapplicableNoCanonicalFrame),
            _ => Err(refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close.schema_impact.disposition",
                "one exact AC58 canonical-schema impact code in 1..=6",
                code,
            )),
        }
    }
}

/// Closed AC58 policy for the canonical predecessor of one schema row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalSchemaMigrationPolicyV1 {
    /// The schema has no canonical predecessor.
    NoSchemaPredecessor = 1,
    /// V1 bytes remain decode-only compatibility evidence.
    V1DecodeOnlyCompatibilityEvidence = 2,
    /// V1 bytes are explicitly retired from authoritative construction.
    V1Retired = 3,
}

impl CanonicalSchemaMigrationPolicyV1 {
    /// Every admitted policy in exact code order.
    pub const ALL: [Self; 3] = [
        Self::NoSchemaPredecessor,
        Self::V1DecodeOnlyCompatibilityEvidence,
        Self::V1Retired,
    ];

    /// Exact unsigned 16-bit code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::NoSchemaPredecessor => "no-schema-predecessor",
            Self::V1DecodeOnlyCompatibilityEvidence => "v1-decode-only-compatibility-evidence",
            Self::V1Retired => "v1-retired",
        }
    }

    /// Parse one exact nonzero closed code. Zero and unknown codes refuse.
    pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::NoSchemaPredecessor),
            2 => Ok(Self::V1DecodeOnlyCompatibilityEvidence),
            3 => Ok(Self::V1Retired),
            _ => Err(refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close.schema_impact.migration_policy",
                "one exact AC58 canonical-schema migration policy code in 1..=3",
                code,
            )),
        }
    }
}

/// Exact target binding for one stable close cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCoverageCloseTargetV1 {
    /// Pure validation is independent of a machine target.
    TargetIndependentPureValidation,
    /// Compile-fail proof is bound to the declared Rust target.
    DeclaredRustTarget,
    /// An in-process projection is bound to the declared host target.
    DeclaredHostTarget,
    /// A result-free contribution is bound to its downstream platform matrix.
    DownstreamPlatformMatrix,
}

impl BaseCoverageCloseTargetV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::TargetIndependentPureValidation => 1,
            Self::DeclaredRustTarget => 2,
            Self::DeclaredHostTarget => 3,
            Self::DownstreamPlatformMatrix => 4,
        }
    }
}

/// Exact execution-profile binding for one stable close cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCoverageCloseProfileV1 {
    /// Ordinary crate test under the frozen test profile.
    CrateTest,
    /// Compile-fail or doctest compile profile.
    CompileFailDoctest,
    /// Pure in-process journey projection.
    InProcessProjection,
    /// Downstream release-built contribution contract.
    DownstreamRelease,
    /// Result-free facet-applicability declaration.
    ApplicabilityDeclaration,
}

impl BaseCoverageCloseProfileV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::CrateTest => 1,
            Self::CompileFailDoctest => 2,
            Self::InProcessProjection => 3,
            Self::DownstreamRelease => 4,
            Self::ApplicabilityDeclaration => 5,
        }
    }
}

/// One exact logical-unit reference.
///
/// Fixed units carry no registry identity. A registered unit is incomplete
/// without the exact extension-registry identity that gives its namespace-local
/// identifier meaning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCoverageCloseLogicalUnitReferenceV1 {
    unit: LogicalUnitV2,
    registry_identity: Option<ContentHash>,
}

impl BaseCoverageCloseLogicalUnitReferenceV1 {
    /// Construct a fixed or registry-bound logical-unit reference.
    pub fn new(
        unit: LogicalUnitV2,
        registry_identity: Option<ContentHash>,
    ) -> Result<Self, ConstructionErrorV2> {
        match (unit.registered_id(), registry_identity) {
            (None, None) | (Some(_), Some(_)) => Ok(Self {
                unit,
                registry_identity,
            }),
            (Some(_), None) => Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.logical_unit.registry_identity",
                "an exact extension-registry identity for a registered logical unit",
                unit.tag(),
            )),
            (None, Some(_)) => Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close.logical_unit.registry_identity",
                "no extension-registry identity for a fixed logical unit",
                unit.tag(),
            )),
        }
    }

    /// Construct one exact fixed logical-unit reference.
    pub fn fixed(unit: LogicalUnitV2) -> Result<Self, ConstructionErrorV2> {
        Self::new(unit, None)
    }

    /// The exact closed or registered logical unit.
    #[must_use]
    pub const fn unit(self) -> LogicalUnitV2 {
        self.unit
    }

    /// Registry identity required exactly for a registered logical unit.
    #[must_use]
    pub const fn registry_identity(self) -> Option<ContentHash> {
        self.registry_identity
    }
}

/// The exact seven hard/soft budget axes in canonical order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseBudgetAxisV1 {
    /// Wall-clock allowance in nanoseconds.
    Time = 1,
    /// Resident-memory allowance in logical bytes.
    Memory = 2,
    /// Logical work in one exact logical unit.
    LogicalWork = 3,
    /// Aggregate process launches, distinct from child-process shape.
    Processes = 4,
    /// Artifact output in canonical encoded bytes.
    Artifacts = 5,
    /// Command and child output in canonical encoded bytes.
    Output = 6,
    /// Detailed deterministic log output in canonical encoded bytes.
    Logs = 7,
}

impl BaseCoverageCloseBudgetAxisV1 {
    /// Exact canonical axis order.
    pub const ALL: [Self; BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1] = [
        Self::Time,
        Self::Memory,
        Self::LogicalWork,
        Self::Processes,
        Self::Artifacts,
        Self::Output,
        Self::Logs,
    ];

    /// Frozen one-based axis code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable source-facing axis name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Time => "time",
            Self::Memory => "memory",
            Self::LogicalWork => "logical-work",
            Self::Processes => "processes",
            Self::Artifacts => "artifacts",
            Self::Output => "output",
            Self::Logs => "logs",
        }
    }

    /// Frozen primitive width for both hard and soft values.
    #[must_use]
    pub const fn width(self) -> BaseCoverageCloseBudgetWidthV1 {
        match self {
            Self::Time | Self::Memory | Self::Artifacts | Self::Output | Self::Logs => {
                BaseCoverageCloseBudgetWidthV1::U64
            }
            Self::LogicalWork => BaseCoverageCloseBudgetWidthV1::U128,
            Self::Processes => BaseCoverageCloseBudgetWidthV1::U32,
        }
    }

    /// Fixed unit required by this axis, or `None` for source-declared work.
    #[must_use]
    pub const fn fixed_unit(self) -> Option<LogicalUnitV2> {
        match self {
            Self::Time => Some(LogicalUnitV2::Nanoseconds),
            Self::Memory => Some(LogicalUnitV2::LogicalBytes),
            Self::LogicalWork => None,
            Self::Processes => Some(LogicalUnitV2::Count),
            Self::Artifacts | Self::Output | Self::Logs => Some(LogicalUnitV2::EncodedBytes),
        }
    }
}

/// Exact primitive width of one budget axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseCoverageCloseBudgetWidthV1 {
    /// Unsigned 32-bit process count, tag 1.
    U32,
    /// Unsigned 64-bit time/byte quantity, tag 2.
    U64,
    /// Unsigned 128-bit logical-work quantity, tag 3.
    U128,
}

impl BaseCoverageCloseBudgetWidthV1 {
    /// Frozen width tag.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::U32 => 1,
            Self::U64 => 2,
            Self::U128 => 3,
        }
    }
}

/// Width-preserving budget quantity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseCoverageCloseBudgetValueV1 {
    /// Process count.
    U32(u32),
    /// Time or bytes.
    U64(u64),
    /// Logical work.
    U128(u128),
}

impl BaseCoverageCloseBudgetValueV1 {
    /// Exact primitive width.
    #[must_use]
    pub const fn width(self) -> BaseCoverageCloseBudgetWidthV1 {
        match self {
            Self::U32(_) => BaseCoverageCloseBudgetWidthV1::U32,
            Self::U64(_) => BaseCoverageCloseBudgetWidthV1::U64,
            Self::U128(_) => BaseCoverageCloseBudgetWidthV1::U128,
        }
    }

    #[must_use]
    const fn is_zero(self) -> bool {
        match self {
            Self::U32(value) => value == 0,
            Self::U64(value) => value == 0,
            Self::U128(value) => value == 0,
        }
    }

    fn exceeds(self, other: Self) -> Result<bool, ConstructionErrorV2> {
        match (self, other) {
            (Self::U32(left), Self::U32(right)) => Ok(left > right),
            (Self::U64(left), Self::U64(right)) => Ok(left > right),
            (Self::U128(left), Self::U128(right)) => Ok(left > right),
            _ => Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.budget.width",
                "matching frozen primitive widths",
                self.width().code(),
            )),
        }
    }
}

/// Source-owned profile identity for one resolved seven-axis budget set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseCoverageCloseBudgetProfileV1 {
    /// Shared profile for source-owned local pure validators.
    LocalSourceValidation,
    /// Shared profile for immutable downstream release contributions.
    DownstreamSourceContribution,
}

impl BaseCoverageCloseBudgetProfileV1 {
    /// Frozen profile code.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::LocalSourceValidation => 1,
            Self::DownstreamSourceContribution => 2,
        }
    }

    /// Exact source-declared profile name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::LocalSourceValidation => "base-close-local-source-validation-v1",
            Self::DownstreamSourceContribution => "base-close-downstream-source-contribution-v1",
        }
    }

    fn hard_ceiling(self, axis: BaseCoverageCloseBudgetAxisV1) -> BaseCoverageCloseBudgetValueV1 {
        const MIB: u64 = 1024 * 1024;
        const GIB: u64 = 1024 * 1024 * 1024;
        const SECOND_NS: u64 = 1_000_000_000;
        match (self, axis) {
            (Self::LocalSourceValidation, BaseCoverageCloseBudgetAxisV1::Time) => {
                BaseCoverageCloseBudgetValueV1::U64(900 * SECOND_NS)
            }
            (Self::DownstreamSourceContribution, BaseCoverageCloseBudgetAxisV1::Time) => {
                BaseCoverageCloseBudgetValueV1::U64(86_400 * SECOND_NS)
            }
            (Self::LocalSourceValidation, BaseCoverageCloseBudgetAxisV1::Memory) => {
                BaseCoverageCloseBudgetValueV1::U64(16 * GIB)
            }
            (Self::DownstreamSourceContribution, BaseCoverageCloseBudgetAxisV1::Memory) => {
                BaseCoverageCloseBudgetValueV1::U64(128 * GIB)
            }
            (_, BaseCoverageCloseBudgetAxisV1::LogicalWork) => {
                BaseCoverageCloseBudgetValueV1::U128(u128::MAX)
            }
            (_, BaseCoverageCloseBudgetAxisV1::Processes) => {
                BaseCoverageCloseBudgetValueV1::U32(256)
            }
            (_, BaseCoverageCloseBudgetAxisV1::Artifacts) => {
                BaseCoverageCloseBudgetValueV1::U64(64 * MIB)
            }
            (_, BaseCoverageCloseBudgetAxisV1::Output) => {
                BaseCoverageCloseBudgetValueV1::U64(5 * MIB)
            }
            (Self::LocalSourceValidation, BaseCoverageCloseBudgetAxisV1::Logs) => {
                BaseCoverageCloseBudgetValueV1::U64(64 * MIB)
            }
            (Self::DownstreamSourceContribution, BaseCoverageCloseBudgetAxisV1::Logs) => {
                BaseCoverageCloseBudgetValueV1::U64(512 * MIB)
            }
        }
    }
}

/// One explicit hard/soft budget row with exact axis, width, and unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCoverageCloseTypedBudgetV1 {
    axis: BaseCoverageCloseBudgetAxisV1,
    hard: BaseCoverageCloseBudgetValueV1,
    soft: BaseCoverageCloseBudgetValueV1,
    unit: BaseCoverageCloseLogicalUnitReferenceV1,
}

impl BaseCoverageCloseTypedBudgetV1 {
    /// Construct one typed row. Soft may be zero but is never inferred.
    pub fn new(
        axis: BaseCoverageCloseBudgetAxisV1,
        hard: BaseCoverageCloseBudgetValueV1,
        soft: BaseCoverageCloseBudgetValueV1,
        unit: BaseCoverageCloseLogicalUnitReferenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if hard.width() != axis.width() || soft.width() != axis.width() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.budget.width",
                "the exact primitive width declared by the budget axis",
                hard.width().code(),
            ));
        }
        if hard.is_zero() {
            return Err(refusal(
                ConstructionErrorKindV2::Zero,
                "coverage.close.five_explicits.budget.hard",
                "one explicit nonzero hard budget",
                axis.name(),
            ));
        }
        if soft.exceeds(hard)? {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.budget.soft",
                "a source-declared soft budget no greater than its hard budget",
                axis.name(),
            ));
        }
        if let Some(expected) = axis.fixed_unit()
            && unit != BaseCoverageCloseLogicalUnitReferenceV1::fixed(expected)?
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.budget.unit",
                "the exact fixed logical unit declared by the budget axis",
                unit.unit().tag(),
            ));
        }
        Ok(Self {
            axis,
            hard,
            soft,
            unit,
        })
    }

    #[must_use]
    pub const fn axis(self) -> BaseCoverageCloseBudgetAxisV1 {
        self.axis
    }

    #[must_use]
    pub const fn hard(self) -> BaseCoverageCloseBudgetValueV1 {
        self.hard
    }

    #[must_use]
    pub const fn soft(self) -> BaseCoverageCloseBudgetValueV1 {
        self.soft
    }

    #[must_use]
    pub const fn unit(self) -> BaseCoverageCloseLogicalUnitReferenceV1 {
        self.unit
    }
}

/// The exact ordered seven-axis budget profile resolved for one source cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCoverageCloseBudgetSetV1 {
    profile: BaseCoverageCloseBudgetProfileV1,
    rows: [BaseCoverageCloseTypedBudgetV1; BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1],
}

impl BaseCoverageCloseBudgetSetV1 {
    /// Construct an exact profile. Missing, extra, duplicate, reordered,
    /// wrong-width, wrong-unit, and over-ceiling rows refuse.
    pub fn new(
        profile: BaseCoverageCloseBudgetProfileV1,
        rows: Vec<BaseCoverageCloseTypedBudgetV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if rows.len() < BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1 {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.five_explicits.budget.rows",
                "all seven exact budget-axis rows",
                rows.len(),
            ));
        }
        if rows.len() > BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1 {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close.five_explicits.budget.rows",
                "no row beyond the seven exact budget axes",
                rows.len(),
            ));
        }
        let mut seen = BTreeSet::new();
        for row in &rows {
            if !seen.insert(row.axis) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.five_explicits.budget.axis",
                    "each exact budget axis exactly once",
                    row.axis.name(),
                ));
            }
        }
        for (expected, row) in BaseCoverageCloseBudgetAxisV1::ALL.iter().zip(&rows) {
            if row.axis != *expected {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.close.five_explicits.budget.axis",
                    "the exact seven-axis canonical order",
                    row.axis.name(),
                ));
            }
            if row.hard.exceeds(profile.hard_ceiling(*expected))? {
                return Err(refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "coverage.close.five_explicits.budget.hard",
                    "a hard budget within the source profile's governing Runner ceiling",
                    row.axis.name(),
                ));
            }
        }
        let rows: [BaseCoverageCloseTypedBudgetV1; BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1] =
            rows.try_into().map_err(|rows: Vec<_>| {
                refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.five_explicits.budget.rows",
                    "the exact fixed-size seven-axis row vector",
                    rows.len(),
                )
            })?;
        Ok(Self { profile, rows })
    }

    #[must_use]
    pub const fn profile(self) -> BaseCoverageCloseBudgetProfileV1 {
        self.profile
    }

    #[must_use]
    pub const fn rows(
        &self,
    ) -> &[BaseCoverageCloseTypedBudgetV1; BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1] {
        &self.rows
    }

    #[must_use]
    pub const fn row(&self, axis: BaseCoverageCloseBudgetAxisV1) -> BaseCoverageCloseTypedBudgetV1 {
        self.rows[(axis.code() - 1) as usize]
    }

    #[must_use]
    pub const fn time(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Time)
    }

    #[must_use]
    pub const fn memory(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Memory)
    }

    #[must_use]
    pub const fn logical_work(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::LogicalWork)
    }

    #[must_use]
    pub const fn processes(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Processes)
    }

    #[must_use]
    pub const fn artifacts(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Artifacts)
    }

    #[must_use]
    pub const fn output(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Output)
    }

    #[must_use]
    pub const fn logs(&self) -> BaseCoverageCloseTypedBudgetV1 {
        self.row(BaseCoverageCloseBudgetAxisV1::Logs)
    }
}

/// One semantic seed explicit or one exact registered inapplicability reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseCoverageCloseSeedExplicitV1 {
    /// One exact semantic workload seed and exact generator/minimizer versions.
    Applicable {
        material: SeedMaterialV2,
        generator_version: StableTokenV2,
        minimizer_version: StableTokenV2,
    },
    /// The cell consumes no semantic workload randomness.
    Inapplicable { reason: SeedInapplicableCodeV1 },
}

impl BaseCoverageCloseSeedExplicitV1 {
    #[must_use]
    pub const fn material(&self) -> Option<&SeedMaterialV2> {
        match self {
            Self::Applicable { material, .. } => Some(material),
            Self::Inapplicable { .. } => None,
        }
    }

    #[must_use]
    pub const fn generator_version(&self) -> Option<&StableTokenV2> {
        match self {
            Self::Applicable {
                generator_version, ..
            } => Some(generator_version),
            Self::Inapplicable { .. } => None,
        }
    }

    #[must_use]
    pub const fn minimizer_version(&self) -> Option<&StableTokenV2> {
        match self {
            Self::Applicable {
                minimizer_version, ..
            } => Some(minimizer_version),
            Self::Inapplicable { .. } => None,
        }
    }

    #[must_use]
    pub const fn inapplicable_reason(&self) -> Option<SeedInapplicableCodeV1> {
        match self {
            Self::Applicable { .. } => None,
            Self::Inapplicable { reason } => Some(*reason),
        }
    }
}

/// Exact physical or logical unit domain for one semantic numeric value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCoverageCloseNumericUnitV1 {
    /// Canonical physical unit with exact positive rational scale and seven
    /// ordered SI base-dimension exponents.
    Physical(UnitV2),
    /// Closed or registry-bound logical unit.
    Logical(BaseCoverageCloseLogicalUnitReferenceV1),
}

impl BaseCoverageCloseNumericUnitV1 {
    /// Construct one exact physical-unit reference.
    #[must_use]
    pub const fn physical(unit: UnitV2) -> Self {
        Self::Physical(unit)
    }

    /// Construct one exact fixed logical-unit reference.
    pub fn logical(unit: LogicalUnitV2) -> Result<Self, ConstructionErrorV2> {
        BaseCoverageCloseLogicalUnitReferenceV1::fixed(unit).map(Self::Logical)
    }

    /// Construct one exact registry-bound logical-unit reference.
    pub fn registered_logical(
        unit: LogicalUnitV2,
        registry_identity: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        BaseCoverageCloseLogicalUnitReferenceV1::new(unit, Some(registry_identity))
            .map(Self::Logical)
    }
}

/// The three independently rooted semantic numeric profiles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseNumericPartitionV1 {
    /// Values consumed as semantic case inputs.
    Inputs = 1,
    /// Values granted as explicit numeric limits or tolerances.
    Grants = 2,
    /// Values expected as semantic observations.
    Observations = 3,
}

impl BaseCoverageCloseNumericPartitionV1 {
    /// Frozen partition code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact source-facing profile name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Inputs => "semantic-numeric-inputs",
            Self::Grants => "semantic-numeric-grants",
            Self::Observations => "expected-numeric-observations",
        }
    }
}

/// One named semantic numeric input, grant, or expected observation.
///
/// Source classification fields and incidental counters do not enter this
/// surface. An exact-empty input, grant, or observation profile is represented
/// by an empty, present vector with its own component root in
/// [`BaseCoverageCloseFiveExplicitsV1`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseNumericExplicitV1 {
    name: StableTokenV2,
    value: NumericValueV2,
    unit: BaseCoverageCloseNumericUnitV1,
}

impl BaseCoverageCloseNumericExplicitV1 {
    /// Construct one complete closed-union numeric value and exact unit.
    #[must_use]
    pub fn new(
        name: StableTokenV2,
        value: NumericValueV2,
        unit: BaseCoverageCloseNumericUnitV1,
    ) -> Self {
        Self { name, value, unit }
    }

    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    #[must_use]
    pub const fn value(&self) -> &NumericValueV2 {
        &self.value
    }

    #[must_use]
    pub const fn unit(&self) -> BaseCoverageCloseNumericUnitV1 {
        self.unit
    }
}

/// Exact API/wire/schema/source/build/toolchain/target/profile/feature versions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseVersionSetV1 {
    api_generation: RunnerApiGeneration,
    wire_version: RunnerWireVersion,
    schema_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    toolchain_root: ToolchainIdentityRootV2,
    target: BaseCoverageCloseTargetV1,
    profile: BaseCoverageCloseProfileV1,
    feature_set_root: ContentHash,
}

impl BaseCoverageCloseVersionSetV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "every required version dimension remains separately typed"
    )]
    pub const fn new(
        api_generation: RunnerApiGeneration,
        wire_version: RunnerWireVersion,
        schema_root: ContentHash,
        source_root: SourceIdentityRootV2,
        build_root: BuildIdentityRootV2,
        toolchain_root: ToolchainIdentityRootV2,
        target: BaseCoverageCloseTargetV1,
        profile: BaseCoverageCloseProfileV1,
        feature_set_root: ContentHash,
    ) -> Self {
        Self {
            api_generation,
            wire_version,
            schema_root,
            source_root,
            build_root,
            toolchain_root,
            target,
            profile,
            feature_set_root,
        }
    }

    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        self.api_generation
    }

    #[must_use]
    pub const fn wire_version(&self) -> RunnerWireVersion {
        self.wire_version
    }

    #[must_use]
    pub const fn schema_root(&self) -> ContentHash {
        self.schema_root
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
    pub const fn toolchain_root(&self) -> &ToolchainIdentityRootV2 {
        &self.toolchain_root
    }

    #[must_use]
    pub const fn target(&self) -> BaseCoverageCloseTargetV1 {
        self.target
    }

    #[must_use]
    pub const fn profile(&self) -> BaseCoverageCloseProfileV1 {
        self.profile
    }

    #[must_use]
    pub const fn feature_set_root(&self) -> ContentHash {
        self.feature_set_root
    }
}

/// Exact number of source-declared close capability descriptors.
pub const BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1: usize = 5;
/// Maximum number of base semantic capability IDs accepted by a base set.
///
/// Registered-extension capabilities use a distinct nominal registry and
/// retain their independent 64-row ceiling; they cannot enter this base set.
pub const BASE_COVERAGE_CLOSE_CAPABILITY_SET_MAX_V1: usize =
    BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1;
/// Sole source owner for every exact close capability descriptor.
pub const BASE_COVERAGE_CLOSE_CAPABILITY_OWNER_V1: &str = "fs-evidence-runner.coverage";
/// Exact no-claim shared by all declaration-side capability contracts.
pub const BASE_COVERAGE_CLOSE_CAPABILITY_CONTRACT_NO_CLAIM_V1: &str =
    "capability-contract-proves-no-acquisition-effect-success-or-authority";

#[derive(Debug, Clone, Copy)]
struct BaseCoverageCloseCapabilityDefinitionV1 {
    stable_id: &'static str,
    policy: BaseCoverageCloseCapabilityPolicyV1,
    no_claim: &'static str,
}

const EXACT_CLOSE_CAPABILITY_DEFINITIONS_V1: [BaseCoverageCloseCapabilityDefinitionV1;
    BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1] = [
    BaseCoverageCloseCapabilityDefinitionV1 {
        stable_id: "fs-evidence-runner.close.control-input.read",
        policy: BaseCoverageCloseCapabilityPolicyV1::DeclaredControlInputRead,
        no_claim: "declared-control-input-read-proves-no-content-authenticity-or-execution",
    },
    BaseCoverageCloseCapabilityDefinitionV1 {
        stable_id: "fs-evidence-runner.close.release-process.control",
        policy: BaseCoverageCloseCapabilityPolicyV1::VersionBoundReleaseProcessControl,
        no_claim: "declared-process-control-proves-no-launch-drain-success-or-version-match",
    },
    BaseCoverageCloseCapabilityDefinitionV1 {
        stable_id: "fs-evidence-runner.close.retained-evidence.write",
        policy: BaseCoverageCloseCapabilityPolicyV1::AttemptConfinedRetainedEvidenceWrite,
        no_claim: "declared-retention-write-proves-no-completeness-durability-or-validity",
    },
    BaseCoverageCloseCapabilityDefinitionV1 {
        stable_id: "fs-evidence-runner.close.evidence-input.read",
        policy: BaseCoverageCloseCapabilityPolicyV1::DeclaredEvidenceInputRead,
        no_claim: "declared-evidence-input-read-proves-no-evidence-validity-or-verification",
    },
    BaseCoverageCloseCapabilityDefinitionV1 {
        stable_id: "fs-evidence-runner.close.publication-output.commit",
        policy: BaseCoverageCloseCapabilityPolicyV1::SelectionBoundPublicationOutputCommit,
        no_claim: "declared-publication-output-proves-no-write-durability-receipt-or-authority",
    },
];

/// Nonzero semantic capability identifier.
///
/// The private field prevents raw integers, capability-policy roots, owners,
/// routes, or physical resource identifiers from substituting for a semantic
/// capability ID. Registry-aware constructors still reject nonzero unknown
/// codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseCapabilityIdV1 {
    code: NonZeroU16,
}

impl BaseCoverageCloseCapabilityIdV1 {
    /// Construct one nonzero ID for subsequent registry-aware validation.
    pub fn new(code: u16) -> Result<Self, ConstructionErrorV2> {
        let code = NonZeroU16::new(code).ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::Zero,
                "coverage.close.capability.id",
                "a nonzero u16 semantic capability ID",
                code,
            )
        })?;
        Ok(Self { code })
    }

    /// Exact unsigned 16-bit code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.code.get()
    }
}

/// Closed semantic scope for one declared close capability.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseCapabilityPolicyV1 {
    /// Read exact source-declared control inputs.
    DeclaredControlInputRead = 1,
    /// Control only the exact version-bound release process graph.
    VersionBoundReleaseProcessControl = 2,
    /// Write evidence only inside the capability-confined retained output.
    AttemptConfinedRetainedEvidenceWrite = 3,
    /// Read exact evidence inputs for the registered verification route.
    DeclaredEvidenceInputRead = 4,
    /// Commit only the selection-bound publication output.
    SelectionBoundPublicationOutputCommit = 5,
}

impl BaseCoverageCloseCapabilityPolicyV1 {
    /// Every admitted semantic scope in registry order.
    pub const ALL: [Self; BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1] = [
        Self::DeclaredControlInputRead,
        Self::VersionBoundReleaseProcessControl,
        Self::AttemptConfinedRetainedEvidenceWrite,
        Self::DeclaredEvidenceInputRead,
        Self::SelectionBoundPublicationOutputCommit,
    ];

    /// Exact unsigned 16-bit code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact stable semantic-scope name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::DeclaredControlInputRead => "declared-control-input-read",
            Self::VersionBoundReleaseProcessControl => "version-bound-release-process-control",
            Self::AttemptConfinedRetainedEvidenceWrite => {
                "attempt-confined-retained-evidence-write"
            }
            Self::DeclaredEvidenceInputRead => "declared-evidence-input-read",
            Self::SelectionBoundPublicationOutputCommit => {
                "selection-bound-publication-output-commit"
            }
        }
    }
}

/// One exact source-declared semantic capability descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilityDescriptorV1 {
    id: BaseCoverageCloseCapabilityIdV1,
    stable_id: StableTokenV2,
    owner: StableTokenV2,
    policy: BaseCoverageCloseCapabilityPolicyV1,
    no_claim: StableTokenV2,
    root: BaseCoverageCloseCapabilityDescriptorRootV1,
}

impl BaseCoverageCloseCapabilityDescriptorV1 {
    fn frozen(
        code: u16,
        stable_id: &'static str,
        policy: BaseCoverageCloseCapabilityPolicyV1,
        no_claim: &'static str,
    ) -> Result<Self, ConstructionErrorV2> {
        let id = BaseCoverageCloseCapabilityIdV1::new(code)?;
        let stable_id = StableTokenV2::new(stable_id).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability.stable_id",
                "one exact bounded stable capability ID",
                code,
            )
        })?;
        let owner = StableTokenV2::new(BASE_COVERAGE_CLOSE_CAPABILITY_OWNER_V1).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability.owner",
                "the exact source-declared capability owner",
                code,
            )
        })?;
        let no_claim = StableTokenV2::new(no_claim).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability.no_claim",
                "one exact bounded capability no-claim",
                code,
            )
        })?;
        let root = close_capability_descriptor_root(id, &stable_id, &owner, policy, &no_claim)?;
        Ok(Self {
            id,
            stable_id,
            owner,
            policy,
            no_claim,
            root,
        })
    }

    /// Exact nonzero registry ID.
    #[must_use]
    pub const fn id(&self) -> BaseCoverageCloseCapabilityIdV1 {
        self.id
    }

    /// Exact stable semantic capability ID.
    #[must_use]
    pub const fn stable_id(&self) -> &StableTokenV2 {
        &self.stable_id
    }

    /// Sole source owner, never a downstream Bead or route owner.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Closed semantic scope.
    #[must_use]
    pub const fn policy(&self) -> BaseCoverageCloseCapabilityPolicyV1 {
        self.policy
    }

    /// Exact descriptor-specific no-claim.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal root of the complete descriptor row.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseCapabilityDescriptorRootV1 {
        self.root
    }
}

/// Exact ordered source registry of semantic close capabilities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilityRegistryV1 {
    rows: Box<[BaseCoverageCloseCapabilityDescriptorV1]>,
    root: BaseCoverageCloseCapabilityRegistryRootV1,
}

impl BaseCoverageCloseCapabilityRegistryV1 {
    /// Construct the sole exact five-row source registry.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let rows = EXACT_CLOSE_CAPABILITY_DEFINITIONS_V1
            .iter()
            .enumerate()
            .map(|(index, definition)| {
                BaseCoverageCloseCapabilityDescriptorV1::frozen(
                    u16::try_from(index + 1).map_err(|_| {
                        refusal(
                            ConstructionErrorKindV2::TooLarge,
                            "coverage.close.capability_registry.id",
                            "one exact u16 base capability ID",
                            index + 1,
                        )
                    })?,
                    definition.stable_id,
                    definition.policy,
                    definition.no_claim,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::try_from_rows(rows)
    }

    fn try_from_rows(
        rows: Vec<BaseCoverageCloseCapabilityDescriptorV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        validate_close_capability_registry_rows(&rows)?;
        let root = close_capability_registry_root(&rows)?;
        Ok(Self {
            rows: rows.into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct the complete registry against exact source rows and root.
    pub fn reconstruct_exact(
        &self,
        presented_rows: &[BaseCoverageCloseCapabilityDescriptorV1],
        presented_root: BaseCoverageCloseCapabilityRegistryRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if presented_rows.len() < self.rows.len() {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.capability_registry.rows",
                "the complete exact five-row source capability registry",
                presented_rows.len(),
            ));
        }
        if presented_rows.len() > self.rows.len() {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close.capability_registry.rows",
                "no row beyond the complete exact source capability registry",
                presented_rows.len(),
            ));
        }
        let candidate = Self::try_from_rows(presented_rows.to_vec())?;
        if candidate.rows.as_ref() != self.rows.as_ref() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_registry.rows",
                "the exact source descriptor rows without mutation",
                candidate.rows.len(),
            ));
        }
        if candidate.root != presented_root || candidate.root != self.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_registry.root",
                "the exact nominal root of the source capability registry",
                candidate.rows.len(),
            ));
        }
        Ok(candidate)
    }

    /// Exact ordered source rows.
    #[must_use]
    pub fn rows(&self) -> &[BaseCoverageCloseCapabilityDescriptorV1] {
        &self.rows
    }

    /// Find one registered descriptor by nominal ID.
    #[must_use]
    pub fn descriptor(
        &self,
        id: BaseCoverageCloseCapabilityIdV1,
    ) -> Option<&BaseCoverageCloseCapabilityDescriptorV1> {
        self.rows
            .binary_search_by_key(&id, BaseCoverageCloseCapabilityDescriptorV1::id)
            .ok()
            .map(|index| &self.rows[index])
    }

    /// Find one frozen descriptor by its exact stable semantic ID.
    ///
    /// Owners, routes, scripts, paths, and policy-root values are not accepted
    /// by this lookup and therefore cannot be promoted into capabilities.
    #[must_use]
    pub fn descriptor_by_stable_id(
        &self,
        stable_id: &str,
    ) -> Option<&BaseCoverageCloseCapabilityDescriptorV1> {
        self.rows
            .iter()
            .find(|descriptor| descriptor.stable_id.as_str() == stable_id)
    }

    /// Find one frozen descriptor by its closed semantic policy.
    #[must_use]
    pub fn descriptor_by_policy(
        &self,
        policy: BaseCoverageCloseCapabilityPolicyV1,
    ) -> Option<&BaseCoverageCloseCapabilityDescriptorV1> {
        self.rows
            .iter()
            .find(|descriptor| descriptor.policy == policy)
    }

    /// Nominal root of the exact ordered registry.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseCapabilityRegistryRootV1 {
        self.root
    }
}

/// Exact declaration-side close capability contract profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseCapabilityProfileV1 {
    None = 1,
    ReleaseControl = 2,
    ReleasePublication = 3,
    ReleaseVerification = 4,
    ReleaseCanonical = 5,
}

impl BaseCoverageCloseCapabilityProfileV1 {
    /// Every exact profile in contract-registry order.
    pub const ALL: [Self; 5] = [
        Self::None,
        Self::ReleaseControl,
        Self::ReleasePublication,
        Self::ReleaseVerification,
        Self::ReleaseCanonical,
    ];

    /// Exact unsigned 16-bit code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact stable profile ID.
    #[must_use]
    pub const fn stable_id(self) -> &'static str {
        match self {
            Self::None => "fs-evidence-runner.close-capability.none.v1",
            Self::ReleaseControl => "fs-evidence-runner.close-capability.release-control.v1",
            Self::ReleasePublication => {
                "fs-evidence-runner.close-capability.release-publication.v1"
            }
            Self::ReleaseVerification => {
                "fs-evidence-runner.close-capability.release-verification.v1"
            }
            Self::ReleaseCanonical => "fs-evidence-runner.close-capability.release-canonical.v1",
        }
    }

    fn required_codes(self) -> &'static [u16] {
        match self {
            Self::None => &[],
            Self::ReleaseControl => &[1, 2, 3],
            Self::ReleasePublication => &[1, 2, 3, 5],
            Self::ReleaseVerification => &[1, 2, 3, 4],
            Self::ReleaseCanonical => &[1, 2, 3, 4, 5],
        }
    }
}

/// One exact row in the source-owned capability-profile registry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilityProfileDescriptorV1 {
    profile: BaseCoverageCloseCapabilityProfileV1,
    stable_id: StableTokenV2,
    required: Box<[BaseCoverageCloseCapabilityIdV1]>,
    permitted: Box<[BaseCoverageCloseCapabilityIdV1]>,
    no_claim: StableTokenV2,
}

impl BaseCoverageCloseCapabilityProfileDescriptorV1 {
    fn frozen(
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        profile: BaseCoverageCloseCapabilityProfileV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let stable_id = StableTokenV2::new(profile.stable_id()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_profile.stable_id",
                "one exact bounded stable capability-profile ID",
                profile.code(),
            )
        })?;
        let required = profile
            .required_codes()
            .iter()
            .copied()
            .map(BaseCoverageCloseCapabilityIdV1::new)
            .collect::<Result<Vec<_>, _>>()?;
        validate_close_capability_id_set(
            "coverage.close.capability_profile.required",
            registry,
            &required,
        )?;
        let no_claim = StableTokenV2::new(BASE_COVERAGE_CLOSE_CAPABILITY_CONTRACT_NO_CLAIM_V1)
            .map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.capability_profile.no_claim",
                    "the exact frozen capability-contract no-claim",
                    profile.code(),
                )
            })?;
        Ok(Self {
            profile,
            stable_id,
            permitted: required.clone().into_boxed_slice(),
            required: required.into_boxed_slice(),
            no_claim,
        })
    }

    #[must_use]
    pub const fn profile(&self) -> BaseCoverageCloseCapabilityProfileV1 {
        self.profile
    }

    #[must_use]
    pub const fn stable_id(&self) -> &StableTokenV2 {
        &self.stable_id
    }

    #[must_use]
    pub fn required(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.required
    }

    #[must_use]
    pub fn permitted(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.permitted
    }

    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// Exact independently rooted ordered registry of capability contract profiles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilityProfileRegistryV1 {
    capability_registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    rows: Box<[BaseCoverageCloseCapabilityProfileDescriptorV1]>,
    root: BaseCoverageCloseCapabilityProfileRegistryRootV1,
}

impl BaseCoverageCloseCapabilityProfileRegistryV1 {
    /// Construct the sole exact five-profile source registry.
    pub fn frozen(
        capability_registry: &BaseCoverageCloseCapabilityRegistryV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let rows = BaseCoverageCloseCapabilityProfileV1::ALL
            .into_iter()
            .map(|profile| {
                BaseCoverageCloseCapabilityProfileDescriptorV1::frozen(capability_registry, profile)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::try_from_rows(capability_registry, rows)
    }

    fn try_from_rows(
        capability_registry: &BaseCoverageCloseCapabilityRegistryV1,
        rows: Vec<BaseCoverageCloseCapabilityProfileDescriptorV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        validate_close_capability_profile_rows(capability_registry, &rows)?;
        let root = close_capability_profile_registry_root(capability_registry.root(), &rows)?;
        Ok(Self {
            capability_registry_root: capability_registry.root(),
            rows: rows.into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct every exact profile row and the independent registry root.
    pub fn reconstruct_exact(
        &self,
        capability_registry: &BaseCoverageCloseCapabilityRegistryV1,
        presented_rows: &[BaseCoverageCloseCapabilityProfileDescriptorV1],
        presented_root: BaseCoverageCloseCapabilityProfileRegistryRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if presented_rows.len() < self.rows.len() {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.capability_profile_registry.rows",
                "the complete exact five-profile source registry",
                presented_rows.len(),
            ));
        }
        if presented_rows.len() > self.rows.len() {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close.capability_profile_registry.rows",
                "no row beyond the complete exact capability-profile registry",
                presented_rows.len(),
            ));
        }
        let candidate = Self::try_from_rows(capability_registry, presented_rows.to_vec())?;
        if candidate.rows.as_ref() != self.rows.as_ref()
            || candidate.capability_registry_root != self.capability_registry_root
            || candidate.root != presented_root
            || candidate.root != self.root
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_profile_registry.reconstruction",
                "the exact ordered profiles and independently rooted registry",
                presented_rows.len(),
            ));
        }
        Ok(candidate)
    }

    #[must_use]
    pub const fn capability_registry_root(&self) -> BaseCoverageCloseCapabilityRegistryRootV1 {
        self.capability_registry_root
    }

    #[must_use]
    pub fn rows(&self) -> &[BaseCoverageCloseCapabilityProfileDescriptorV1] {
        &self.rows
    }

    #[must_use]
    pub fn descriptor(
        &self,
        profile: BaseCoverageCloseCapabilityProfileV1,
    ) -> Option<&BaseCoverageCloseCapabilityProfileDescriptorV1> {
        self.rows
            .binary_search_by_key(
                &profile,
                BaseCoverageCloseCapabilityProfileDescriptorV1::profile,
            )
            .ok()
            .map(|index| &self.rows[index])
    }

    /// Find one exact contract profile by stable profile ID.
    #[must_use]
    pub fn descriptor_by_stable_id(
        &self,
        stable_id: &str,
    ) -> Option<&BaseCoverageCloseCapabilityProfileDescriptorV1> {
        self.rows
            .iter()
            .find(|descriptor| descriptor.stable_id.as_str() == stable_id)
    }

    /// Resolve one exact contract profile from its stable profile ID.
    #[must_use]
    pub fn profile_by_stable_id(
        &self,
        stable_id: &str,
    ) -> Option<BaseCoverageCloseCapabilityProfileV1> {
        self.descriptor_by_stable_id(stable_id)
            .map(BaseCoverageCloseCapabilityProfileDescriptorV1::profile)
    }

    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseCapabilityProfileRegistryRootV1 {
        self.root
    }
}

/// One exact declaration-side required/permitted capability contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilityContractV1 {
    registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    profile_registry_root: BaseCoverageCloseCapabilityProfileRegistryRootV1,
    profile: BaseCoverageCloseCapabilityProfileV1,
    required: Box<[BaseCoverageCloseCapabilityIdV1]>,
    permitted: Box<[BaseCoverageCloseCapabilityIdV1]>,
    no_claim: StableTokenV2,
    root: BaseCoverageCloseCapabilityContractRootV1,
}

impl BaseCoverageCloseCapabilityContractV1 {
    /// Resolve one exact profile through the frozen capability registry.
    pub fn for_profile(
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        profile_registry: &BaseCoverageCloseCapabilityProfileRegistryV1,
        profile: BaseCoverageCloseCapabilityProfileV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let descriptor = profile_registry.descriptor(profile).ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::UnknownCode,
                "coverage.close.capability_contract.profile",
                "one profile registered by the exact capability-profile registry",
                profile.code(),
            )
        })?;
        Self::try_from_parts(
            registry,
            profile_registry,
            profile,
            descriptor.required.to_vec(),
            descriptor.permitted.to_vec(),
        )
    }

    fn try_from_parts(
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        profile_registry: &BaseCoverageCloseCapabilityProfileRegistryV1,
        profile: BaseCoverageCloseCapabilityProfileV1,
        required: Vec<BaseCoverageCloseCapabilityIdV1>,
        permitted: Vec<BaseCoverageCloseCapabilityIdV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if profile_registry.capability_registry_root != registry.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_contract.profile_registry",
                "the profile registry rooted against the exact capability registry",
                profile.code(),
            ));
        }
        validate_close_capability_id_set(
            "coverage.close.capability_contract.required",
            registry,
            &required,
        )?;
        validate_close_capability_id_set(
            "coverage.close.capability_contract.permitted",
            registry,
            &permitted,
        )?;
        if !close_capability_set_is_subset(&required, &permitted) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_contract.required",
                "required capabilities contained by permitted capabilities",
                required.len(),
            ));
        }
        let descriptor = profile_registry.descriptor(profile).ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::UnknownCode,
                "coverage.close.capability_contract.profile",
                "one profile registered by the exact capability-profile registry",
                profile.code(),
            )
        })?;
        if required.as_slice() != descriptor.required()
            || permitted.as_slice() != descriptor.permitted()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_contract.profile",
                "the exact required-equals-permitted set for the selected profile",
                profile.code(),
            ));
        }
        let no_claim = StableTokenV2::new(BASE_COVERAGE_CLOSE_CAPABILITY_CONTRACT_NO_CLAIM_V1)
            .map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.capability_contract.no_claim",
                    "the exact frozen capability-contract no-claim",
                    profile.code(),
                )
            })?;
        let root = close_capability_contract_root(
            registry.root(),
            profile_registry.root(),
            profile,
            &required,
            &permitted,
            &no_claim,
        )?;
        Ok(Self {
            registry_root: registry.root(),
            profile_registry_root: profile_registry.root(),
            profile,
            required: required.into_boxed_slice(),
            permitted: permitted.into_boxed_slice(),
            no_claim,
            root,
        })
    }

    /// Reconstruct one exact profile contract and its presented nominal root.
    pub fn reconstruct_exact(
        &self,
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        profile_registry: &BaseCoverageCloseCapabilityProfileRegistryV1,
        presented_profile: BaseCoverageCloseCapabilityProfileV1,
        presented_required: &[BaseCoverageCloseCapabilityIdV1],
        presented_permitted: &[BaseCoverageCloseCapabilityIdV1],
        presented_root: BaseCoverageCloseCapabilityContractRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let candidate = Self::try_from_parts(
            registry,
            profile_registry,
            presented_profile,
            presented_required.to_vec(),
            presented_permitted.to_vec(),
        )?;
        if &candidate != self || candidate.root != presented_root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_contract.reconstruction",
                "the exact source profile, sets, registry root, no-claim, and contract root",
                presented_profile.code(),
            ));
        }
        Ok(candidate)
    }

    /// Exact registry identity that resolves every capability ID.
    #[must_use]
    pub const fn registry_root(&self) -> BaseCoverageCloseCapabilityRegistryRootV1 {
        self.registry_root
    }

    /// Exact source-owned capability-profile registry used for resolution.
    #[must_use]
    pub const fn profile_registry_root(&self) -> BaseCoverageCloseCapabilityProfileRegistryRootV1 {
        self.profile_registry_root
    }

    /// Exact selected contract profile.
    #[must_use]
    pub const fn profile(&self) -> BaseCoverageCloseCapabilityProfileV1 {
        self.profile
    }

    /// Exact required capability IDs.
    #[must_use]
    pub fn required(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.required
    }

    /// Exact permitted capability IDs, equal to required for every V1 profile.
    #[must_use]
    pub fn permitted(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.permitted
    }

    /// Exact declaration-only no-claim.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal root of the complete contract.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseCapabilityContractRootV1 {
        self.root
    }
}

/// Literal source-case to declaration-side capability-profile oracle.
///
/// No owner, Bead ID, driver, script, path, root policy, or narrowed policy
/// view participates in this decision.
pub fn base_coverage_close_capability_profile_for_source_case_v1(
    source_class: BaseCoverageManifestClassV1,
    source_case_id: &str,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
) -> Result<BaseCoverageCloseCapabilityProfileV1, ConstructionErrorV2> {
    use BaseCoverageCloseCapabilityProfileV1::{
        None, ReleaseCanonical, ReleaseControl, ReleasePublication, ReleaseVerification,
    };
    use BaseCoverageCloseExecutionScopeV1::{
        CompileFailDoctest, CrateTest, FacetApplicabilityDeclaration,
        ImmutableDownstreamContribution, InProcessProjection,
    };
    use BaseCoverageManifestClassV1::{ExternalE2eScript, ExternalGovernance, ExternalMutation};

    validate_untyped_case_id(source_case_id)?;
    match (source_class, execution_scope, source_case_id) {
        (
            ExternalE2eScript,
            ImmutableDownstreamContribution,
            "external-e2e:publication-state-v2",
        ) => Ok(ReleasePublication),
        (ExternalE2eScript, ImmutableDownstreamContribution, "external-e2e:publication-v2")
        | (
            ExternalGovernance,
            ImmutableDownstreamContribution,
            "external-governance:live-source-dependency-closure",
        ) => Ok(ReleaseControl),
        (
            ExternalE2eScript,
            ImmutableDownstreamContribution,
            "external-e2e:verifier-v2" | "external-e2e:rjoq-handoff-v1",
        ) => Ok(ReleaseVerification),
        (
            ExternalE2eScript,
            ImmutableDownstreamContribution,
            "external-e2e:canonical-runner-v2",
        )
        | (
            ExternalMutation,
            ImmutableDownstreamContribution,
            "external-mutation:base-contract-exact-result-join",
        ) => Ok(ReleaseCanonical),
        (
            ExternalE2eScript | ExternalMutation | ExternalGovernance,
            ImmutableDownstreamContribution,
            _,
        ) => Err(refusal(
            ConstructionErrorKindV2::UnknownCode,
            "coverage.close.capability_contract.source_case_id",
            "one exact literal downstream capability-profile mapping",
            source_case_id,
        )),
        (
            class,
            CrateTest | CompileFailDoctest | InProcessProjection | FacetApplicabilityDeclaration,
            _,
        ) if !matches!(
            class,
            ExternalE2eScript | ExternalMutation | ExternalGovernance
        ) =>
        {
            Ok(None)
        }
        _ => Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.capability_contract.source_scope",
            "one exact source-class, source-case, and execution-scope mapping",
            source_case_id,
        )),
    }
}

/// Runtime-reported semantic capability sets reconciled to one declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseObservedCapabilitySetsV1 {
    registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    contract_root: BaseCoverageCloseCapabilityContractRootV1,
    required: Box<[BaseCoverageCloseCapabilityIdV1]>,
    granted: Box<[BaseCoverageCloseCapabilityIdV1]>,
    observed: Box<[BaseCoverageCloseCapabilityIdV1]>,
    returned: Box<[BaseCoverageCloseCapabilityIdV1]>,
    revoked: Box<[BaseCoverageCloseCapabilityIdV1]>,
    root: BaseCoverageCloseObservedCapabilitySetsRootV1,
}

impl BaseCoverageCloseObservedCapabilitySetsV1 {
    /// Check every AC56 containment and return/revoke law.
    #[allow(
        clippy::too_many_arguments,
        reason = "the five independently inspectable semantic sets are the AC56 contract"
    )]
    pub fn new(
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        contract: &BaseCoverageCloseCapabilityContractV1,
        required: Vec<BaseCoverageCloseCapabilityIdV1>,
        granted: Vec<BaseCoverageCloseCapabilityIdV1>,
        observed: Vec<BaseCoverageCloseCapabilityIdV1>,
        returned: Vec<BaseCoverageCloseCapabilityIdV1>,
        revoked: Vec<BaseCoverageCloseCapabilityIdV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if contract.registry_root != registry.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.registry_root",
                "the exact registry bound by the declaration-side contract",
                required.len(),
            ));
        }
        for (field, values) in [
            (
                "coverage.close.observed_capabilities.required",
                required.as_slice(),
            ),
            (
                "coverage.close.observed_capabilities.granted",
                granted.as_slice(),
            ),
            (
                "coverage.close.observed_capabilities.observed",
                observed.as_slice(),
            ),
            (
                "coverage.close.observed_capabilities.returned",
                returned.as_slice(),
            ),
            (
                "coverage.close.observed_capabilities.revoked",
                revoked.as_slice(),
            ),
        ] {
            validate_close_capability_id_set(field, registry, values)?;
        }
        if required.as_slice() != contract.required() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.required",
                "actual required exactly equal to declared required",
                required.len(),
            ));
        }
        if !close_capability_set_is_subset(&required, &granted) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.granted",
                "every required capability present in granted",
                granted.len(),
            ));
        }
        if !close_capability_set_is_subset(&granted, contract.permitted()) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.permitted",
                "every granted capability declared permitted",
                granted.len(),
            ));
        }
        if !close_capability_set_is_subset(&observed, &granted) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.observed",
                "every observed capability present in granted",
                observed.len(),
            ));
        }
        if !close_capability_set_is_subset(&returned, &granted)
            || !close_capability_set_is_subset(&revoked, &granted)
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.terminal_sets",
                "returned and revoked are subsets of granted",
                returned.len() + revoked.len(),
            ));
        }
        if returned.iter().any(|id| revoked.binary_search(id).is_ok()) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.terminal_sets",
                "returned and revoked are disjoint",
                returned.len() + revoked.len(),
            ));
        }
        if granted
            .iter()
            .any(|id| returned.binary_search(id).is_err() && revoked.binary_search(id).is_err())
        {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.observed_capabilities.terminal_sets",
                "returned union revoked exactly equals granted",
                granted.len(),
            ));
        }
        let root = close_observed_capability_sets_root(
            registry.root,
            contract.root,
            &required,
            &granted,
            &observed,
            &returned,
            &revoked,
        )?;
        Ok(Self {
            registry_root: registry.root,
            contract_root: contract.root,
            required: required.into_boxed_slice(),
            granted: granted.into_boxed_slice(),
            observed: observed.into_boxed_slice(),
            returned: returned.into_boxed_slice(),
            revoked: revoked.into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct checked sets against one caller-presented nominal root.
    #[allow(
        clippy::too_many_arguments,
        reason = "the five independently inspectable semantic sets are the AC56 contract"
    )]
    pub fn reconstruct_exact(
        registry: &BaseCoverageCloseCapabilityRegistryV1,
        contract: &BaseCoverageCloseCapabilityContractV1,
        required: &[BaseCoverageCloseCapabilityIdV1],
        granted: &[BaseCoverageCloseCapabilityIdV1],
        observed: &[BaseCoverageCloseCapabilityIdV1],
        returned: &[BaseCoverageCloseCapabilityIdV1],
        revoked: &[BaseCoverageCloseCapabilityIdV1],
        presented_root: BaseCoverageCloseObservedCapabilitySetsRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let candidate = Self::new(
            registry,
            contract,
            required.to_vec(),
            granted.to_vec(),
            observed.to_vec(),
            returned.to_vec(),
            revoked.to_vec(),
        )?;
        if candidate.root != presented_root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.observed_capabilities.root",
                "the exact nominal root of the reconciled capability sets",
                required.len(),
            ));
        }
        Ok(candidate)
    }

    #[must_use]
    pub const fn registry_root(&self) -> BaseCoverageCloseCapabilityRegistryRootV1 {
        self.registry_root
    }

    #[must_use]
    pub const fn contract_root(&self) -> BaseCoverageCloseCapabilityContractRootV1 {
        self.contract_root
    }

    #[must_use]
    pub fn required(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.required
    }

    #[must_use]
    pub fn granted(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.granted
    }

    #[must_use]
    pub fn observed(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.observed
    }

    #[must_use]
    pub fn returned(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.returned
    }

    #[must_use]
    pub fn revoked(&self) -> &[BaseCoverageCloseCapabilityIdV1] {
        &self.revoked
    }

    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseObservedCapabilitySetsRootV1 {
        self.root
    }
}

/// Independent ceiling for registered-extension capability rows and sets.
pub const BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1: usize = 64;

/// Nonzero registry-local identifier for one extension capability.
///
/// This is nominally distinct from [`BaseCoverageCloseCapabilityIdV1`].
/// There is intentionally no conversion in either direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityIdV1 {
    code: NonZeroU16,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityIdV1 {
    /// Construct one bounded nonzero extension-registry code.
    pub fn new(code: u16) -> Result<Self, ConstructionErrorV2> {
        let code = NonZeroU16::new(code).ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::Zero,
                "coverage.close.extension_capability.id",
                "a nonzero extension capability code in 1..=64",
                code,
            )
        })?;
        if usize::from(code.get()) > BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1 {
            return Err(refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.extension_capability.id",
                "an extension capability code in 1..=64",
                code.get(),
            ));
        }
        Ok(Self { code })
    }

    /// Exact unsigned 16-bit registry-local code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.code.get()
    }
}

/// Stable semantic ID for one registered-extension capability.
///
/// This nominal wrapper prevents an owner, scope, or no-claim token from being
/// substituted even when the underlying text happens to satisfy the same token
/// grammar.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1 {
    value: StableTokenV2,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1 {
    /// Validate one fully namespaced ID disjoint from every base capability.
    pub fn new(value: StableTokenV2) -> Result<Self, ConstructionErrorV2> {
        if !value.as_str().contains('.') || is_base_close_capability_stable_id(value.as_str()) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.extension_capability.stable_id",
                "a fully namespaced extension capability ID disjoint from every base capability",
                value.as_str(),
            ));
        }
        Ok(Self { value })
    }

    /// Exact validated stable ID.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.value.as_str()
    }
}

/// Source-owner identity for one registered-extension capability.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1 {
    value: StableTokenV2,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1 {
    /// Validate one fully namespaced source-owner identity.
    pub fn new(value: StableTokenV2) -> Result<Self, ConstructionErrorV2> {
        if !value.as_str().contains('.') {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.extension_capability.owner",
                "a fully namespaced registered-extension capability owner",
                value.as_str(),
            ));
        }
        Ok(Self { value })
    }

    /// Exact validated source owner.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.value.as_str()
    }
}

/// Closed semantic scope for one registered-extension capability.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityScopeV1 {
    value: StableTokenV2,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityScopeV1 {
    /// Validate one fully namespaced closed semantic scope.
    pub fn new(value: StableTokenV2) -> Result<Self, ConstructionErrorV2> {
        if !value.as_str().contains('.') {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.extension_capability.scope",
                "a fully namespaced registered-extension capability scope",
                value.as_str(),
            ));
        }
        Ok(Self { value })
    }

    /// Exact validated closed semantic scope.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.value.as_str()
    }
}

/// Explicit no-acquisition/no-authority boundary for one extension capability.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1 {
    value: StableTokenV2,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1 {
    /// Preserve one grammar-checked no-claim under its nominal role.
    pub const fn new(value: StableTokenV2) -> Self {
        Self { value }
    }

    /// Exact validated no-claim token.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.value.as_str()
    }
}

/// One source-owned registered-extension capability descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1 {
    id: BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
    stable_id: BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1,
    owner: BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1,
    scope: BaseCoverageCloseRegisteredExtensionCapabilityScopeV1,
    no_claim: BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1,
    root: BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1 {
    /// Construct one bounded, fully namespaced extension descriptor.
    pub fn new(
        id: BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
        stable_id: BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1,
        owner: BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1,
        scope: BaseCoverageCloseRegisteredExtensionCapabilityScopeV1,
        no_claim: BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if stable_id.as_str() == owner.as_str()
            || stable_id.as_str() == scope.as_str()
            || owner.as_str() == scope.as_str()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.extension_capability.semantic_roles",
                "distinct stable capability, owner, and closed-scope identities",
                id.code(),
            ));
        }
        let root = close_registered_extension_capability_descriptor_root(
            id, &stable_id, &owner, &scope, &no_claim,
        )?;
        Ok(Self {
            id,
            stable_id,
            owner,
            scope,
            no_claim,
            root,
        })
    }

    /// Exact registry-local ID.
    #[must_use]
    pub const fn id(&self) -> BaseCoverageCloseRegisteredExtensionCapabilityIdV1 {
        self.id
    }

    /// Exact fully namespaced semantic capability ID.
    #[must_use]
    pub const fn stable_id(&self) -> &BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1 {
        &self.stable_id
    }

    /// Exact source owner.
    #[must_use]
    pub const fn owner(&self) -> &BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1 {
        &self.owner
    }

    /// Exact closed semantic scope.
    #[must_use]
    pub const fn scope(&self) -> &BaseCoverageCloseRegisteredExtensionCapabilityScopeV1 {
        &self.scope
    }

    /// Exact no-acquisition/no-authority boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1 {
        &self.no_claim
    }

    /// Nominal root of the complete descriptor.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1 {
        self.root
    }
}

/// Ordered source-owned registry of extension capabilities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1 {
    rows: Box<[BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1]>,
    root: BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1,
}

impl BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1 {
    /// Construct an exact zero-through-64-row registry in contiguous code
    /// order. Missing codes, duplicate namespaces, and reordered rows refuse.
    pub fn new(
        rows: Vec<BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if rows.len() > BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1 {
            return Err(refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.extension_capability_registry.rows",
                "at most 64 registered-extension capability rows",
                rows.len(),
            ));
        }
        let mut ids = BTreeSet::new();
        for row in &rows {
            if !ids.insert(row.id()) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.extension_capability_registry.code",
                    "globally unique extension capability codes",
                    row.id().code(),
                ));
            }
        }
        let mut stable_ids = BTreeSet::new();
        for row in &rows {
            if !stable_ids.insert(row.stable_id().as_str()) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.extension_capability_registry.stable_id",
                    "globally unique extension capability namespaces",
                    row.stable_id().as_str(),
                ));
            }
        }
        for expected_code in 1..=rows.len() {
            let expected_code = u16::try_from(expected_code).map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "coverage.close.extension_capability_registry.code",
                    "a contiguous u16 code in exact row order",
                    expected_code,
                )
            })?;
            let expected_id =
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(expected_code)?;
            if !ids.contains(&expected_id) {
                return Err(refusal(
                    ConstructionErrorKindV2::Missing,
                    "coverage.close.extension_capability_registry.code",
                    "every contiguous extension capability code starting at one",
                    expected_code,
                ));
            }
        }
        for (index, row) in rows.iter().enumerate() {
            let expected_code = u16::try_from(index + 1).map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "coverage.close.extension_capability_registry.code",
                    "a contiguous u16 code in exact row order",
                    index + 1,
                )
            })?;
            if row.id().code() != expected_code {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.close.extension_capability_registry.code",
                    "contiguous extension capability codes in exact ascending order",
                    row.id().code(),
                ));
            }
        }
        let root = close_registered_extension_capability_registry_root(&rows)?;
        Ok(Self {
            rows: rows.into_boxed_slice(),
            root,
        })
    }

    /// Exact rows in registry order.
    #[must_use]
    pub fn rows(&self) -> &[BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1] {
        &self.rows
    }

    /// Look up one exact registered extension ID.
    #[must_use]
    pub fn descriptor(
        &self,
        id: BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
    ) -> Option<&BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1> {
        self.rows
            .get(usize::from(id.code()) - 1)
            .filter(|row| row.id() == id)
    }

    /// Nominal root of the complete registry.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1 {
        self.root
    }
}

/// Ordered exact set of IDs from one registered-extension capability registry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseRegisteredExtensionCapabilitySetV1 {
    registry_root: BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1,
    values: Box<[BaseCoverageCloseRegisteredExtensionCapabilityIdV1]>,
    root: BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1,
}

impl BaseCoverageCloseRegisteredExtensionCapabilitySetV1 {
    /// Construct an exact zero-through-64-member set in registry order.
    pub fn new(
        registry: &BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1,
        values: Vec<BaseCoverageCloseRegisteredExtensionCapabilityIdV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if values.len() > BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1 {
            return Err(refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.extension_capability_set.values",
                "at most 64 exact registered-extension capability IDs",
                values.len(),
            ));
        }
        for value in &values {
            if registry.descriptor(*value).is_none() {
                return Err(refusal(
                    ConstructionErrorKindV2::UnknownCode,
                    "coverage.close.extension_capability_set.value",
                    "an ID in the exact presented extension registry",
                    value.code(),
                ));
            }
        }
        let mut seen = BTreeSet::new();
        for value in &values {
            if !seen.insert(*value) {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.extension_capability_set.value",
                    "globally unique extension capability IDs in registry order",
                    value.code(),
                ));
            }
        }
        for pair in values.windows(2) {
            if pair[0] > pair[1] {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.close.extension_capability_set.value",
                    "extension capability IDs in registry order",
                    pair[1].code(),
                ));
            }
        }
        let root = close_registered_extension_capability_set_root(registry.root(), &values)?;
        Ok(Self {
            registry_root: registry.root(),
            values: values.into_boxed_slice(),
            root,
        })
    }

    /// Exact registry identity.
    #[must_use]
    pub const fn registry_root(
        &self,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1 {
        self.registry_root
    }

    /// Exact IDs in registry order.
    #[must_use]
    pub fn values(&self) -> &[BaseCoverageCloseRegisteredExtensionCapabilityIdV1] {
        &self.values
    }

    /// Nominal root of the complete set.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1 {
        self.root
    }
}

/// Exact required/granted/observed/returned/revoked capability-ID sets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCapabilitySetsV1 {
    required: Box<[StableTokenV2]>,
    granted: Box<[StableTokenV2]>,
    observed: Box<[StableTokenV2]>,
    returned: Box<[StableTokenV2]>,
    revoked: Box<[StableTokenV2]>,
}

impl BaseCoverageCloseCapabilitySetsV1 {
    pub fn new(
        required: Vec<StableTokenV2>,
        granted: Vec<StableTokenV2>,
        observed: Vec<StableTokenV2>,
        returned: Vec<StableTokenV2>,
        revoked: Vec<StableTokenV2>,
    ) -> Result<Self, ConstructionErrorV2> {
        for (field, values) in [
            ("required", required.as_slice()),
            ("granted", granted.as_slice()),
            ("observed", observed.as_slice()),
            ("returned", returned.as_slice()),
            ("revoked", revoked.as_slice()),
        ] {
            validate_close_capability_set(field, values)?;
        }
        Ok(Self {
            required: required.into_boxed_slice(),
            granted: granted.into_boxed_slice(),
            observed: observed.into_boxed_slice(),
            returned: returned.into_boxed_slice(),
            revoked: revoked.into_boxed_slice(),
        })
    }

    #[must_use]
    pub fn required(&self) -> &[StableTokenV2] {
        &self.required
    }

    #[must_use]
    pub fn granted(&self) -> &[StableTokenV2] {
        &self.granted
    }

    #[must_use]
    pub fn observed(&self) -> &[StableTokenV2] {
        &self.observed
    }

    #[must_use]
    pub fn returned(&self) -> &[StableTokenV2] {
        &self.returned
    }

    #[must_use]
    pub fn revoked(&self) -> &[StableTokenV2] {
        &self.revoked
    }
}

/// All Five Explicits bound into one stable close cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseFiveExplicitsV1 {
    numeric_inputs: Box<[BaseCoverageCloseNumericExplicitV1]>,
    numeric_grants: Box<[BaseCoverageCloseNumericExplicitV1]>,
    numeric_observations: Box<[BaseCoverageCloseNumericExplicitV1]>,
    numeric_inputs_root: ContentHash,
    numeric_grants_root: ContentHash,
    numeric_observations_root: ContentHash,
    seed: BaseCoverageCloseSeedExplicitV1,
    budgets: BaseCoverageCloseBudgetSetV1,
    versions: BaseCoverageCloseVersionSetV1,
    capabilities: BaseCoverageCloseCapabilitySetsV1,
    no_claim: StableTokenV2,
    root: ContentHash,
}

impl BaseCoverageCloseFiveExplicitsV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "all Five Explicits are mandatory and separately inspectable"
    )]
    pub fn new(
        numeric_inputs: Vec<BaseCoverageCloseNumericExplicitV1>,
        numeric_grants: Vec<BaseCoverageCloseNumericExplicitV1>,
        numeric_observations: Vec<BaseCoverageCloseNumericExplicitV1>,
        seed: BaseCoverageCloseSeedExplicitV1,
        budgets: BaseCoverageCloseBudgetSetV1,
        versions: BaseCoverageCloseVersionSetV1,
        capabilities: BaseCoverageCloseCapabilitySetsV1,
        no_claim: StableTokenV2,
    ) -> Result<Self, ConstructionErrorV2> {
        validate_close_numeric_explicit_sequence(
            "coverage.close.five_explicits.numeric_inputs",
            &numeric_inputs,
        )?;
        validate_close_numeric_explicit_sequence(
            "coverage.close.five_explicits.numeric_grants",
            &numeric_grants,
        )?;
        validate_close_numeric_explicit_sequence(
            "coverage.close.five_explicits.numeric_observations",
            &numeric_observations,
        )?;
        let numeric_inputs_root = close_numeric_profile_root(
            BaseCoverageCloseNumericPartitionV1::Inputs,
            &numeric_inputs,
        )?;
        let numeric_grants_root = close_numeric_profile_root(
            BaseCoverageCloseNumericPartitionV1::Grants,
            &numeric_grants,
        )?;
        let numeric_observations_root = close_numeric_profile_root(
            BaseCoverageCloseNumericPartitionV1::Observations,
            &numeric_observations,
        )?;
        let root = close_five_explicits_root(
            numeric_inputs_root,
            numeric_grants_root,
            numeric_observations_root,
            &seed,
            budgets,
            &versions,
            &capabilities,
            &no_claim,
        )?;
        Ok(Self {
            numeric_inputs: numeric_inputs.into_boxed_slice(),
            numeric_grants: numeric_grants.into_boxed_slice(),
            numeric_observations: numeric_observations.into_boxed_slice(),
            numeric_inputs_root,
            numeric_grants_root,
            numeric_observations_root,
            seed,
            budgets,
            versions,
            capabilities,
            no_claim,
            root,
        })
    }

    #[must_use]
    pub fn numeric_inputs(&self) -> &[BaseCoverageCloseNumericExplicitV1] {
        &self.numeric_inputs
    }

    #[must_use]
    pub fn numeric_grants(&self) -> &[BaseCoverageCloseNumericExplicitV1] {
        &self.numeric_grants
    }

    #[must_use]
    pub fn numeric_observations(&self) -> &[BaseCoverageCloseNumericExplicitV1] {
        &self.numeric_observations
    }

    #[must_use]
    pub const fn numeric_inputs_root(&self) -> ContentHash {
        self.numeric_inputs_root
    }

    #[must_use]
    pub const fn numeric_grants_root(&self) -> ContentHash {
        self.numeric_grants_root
    }

    #[must_use]
    pub const fn numeric_observations_root(&self) -> ContentHash {
        self.numeric_observations_root
    }

    #[must_use]
    pub const fn seed(&self) -> &BaseCoverageCloseSeedExplicitV1 {
        &self.seed
    }

    #[must_use]
    pub const fn budgets(&self) -> BaseCoverageCloseBudgetSetV1 {
        self.budgets
    }

    #[must_use]
    pub const fn versions(&self) -> &BaseCoverageCloseVersionSetV1 {
        &self.versions
    }

    #[must_use]
    pub const fn capabilities(&self) -> &BaseCoverageCloseCapabilitySetsV1 {
        &self.capabilities
    }

    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Explicit bounded budgets attached to an immutable downstream contribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseCoverageCloseContributionBudgetsV1 {
    resolved: BaseCoverageCloseBudgetSetV1,
    max_child_processes: u32,
    max_parallel_children: u32,
}

impl BaseCoverageCloseContributionBudgetsV1 {
    /// Construct one downstream source profile and its independent process
    /// shape. The process hard/soft row is an aggregate grant; total-child and
    /// parallel-child limits remain separate structural constraints.
    pub fn new(
        resolved: BaseCoverageCloseBudgetSetV1,
        max_child_processes: u32,
        max_parallel_children: u32,
    ) -> Result<Self, ConstructionErrorV2> {
        if resolved.profile != BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_budgets.profile",
                "the exact downstream source-contribution budget profile",
                resolved.profile.name(),
            ));
        }
        if max_child_processes == 0 || max_parallel_children == 0 {
            return Err(refusal(
                ConstructionErrorKindV2::Zero,
                "coverage.close.contribution_budgets.process_shape",
                "explicit nonzero total-child and parallel-child shape constraints",
                0,
            ));
        }
        if max_parallel_children > max_child_processes {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_budgets.max_parallel_children",
                "a nonzero parallel-child cap no greater than total child processes",
                max_parallel_children,
            ));
        }
        if max_child_processes > 256 || max_parallel_children > 64 {
            return Err(refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.contribution_budgets.process_shape",
                "at most 256 total and 64 parallel children for the downstream profile",
                max_child_processes,
            ));
        }
        let process_hard = match resolved.processes().hard() {
            BaseCoverageCloseBudgetValueV1::U32(value) => value,
            other => {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.contribution_budgets.process_budget",
                    "one exact u32 aggregate-process hard budget",
                    other.width().code(),
                ));
            }
        };
        if max_child_processes > process_hard {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_budgets.max_child_processes",
                "a total-child shape no greater than the independent aggregate-process hard budget",
                max_child_processes,
            ));
        }
        Ok(Self {
            resolved,
            max_child_processes,
            max_parallel_children,
        })
    }

    #[must_use]
    pub const fn resolved(self) -> BaseCoverageCloseBudgetSetV1 {
        self.resolved
    }

    #[must_use]
    pub const fn max_child_processes(self) -> u32 {
        self.max_child_processes
    }

    #[must_use]
    pub const fn max_parallel_children(self) -> u32 {
        self.max_parallel_children
    }
}

/// Result-free payload required for every downstream-owned close cell.
///
/// These are declarative identities and budgets, never claims that the
/// downstream owner executed, retained, or verified anything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseDownstreamContributionV1 {
    literal_expectation_oracle_root: ContentHash,
    semantic_input_root: ContentHash,
    budgets: BaseCoverageCloseContributionBudgetsV1,
    schema_root: ContentHash,
    log_schema_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    downstream_owner: Box<str>,
    downstream_driver: StableTokenV2,
    downstream_script: Box<str>,
    downstream_manifest_path: Box<str>,
    downstream_manifest_root: ContentHash,
    no_claim: Box<str>,
    root: ContentHash,
}

impl BaseCoverageCloseDownstreamContributionV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "the immutable contribution retains every AC53 binding explicitly"
    )]
    pub fn new(
        literal_expectation_oracle_root: ContentHash,
        semantic_input_root: ContentHash,
        budgets: BaseCoverageCloseContributionBudgetsV1,
        schema_root: ContentHash,
        log_schema_root: ContentHash,
        source_root: SourceIdentityRootV2,
        build_root: BuildIdentityRootV2,
        downstream_owner: impl Into<String>,
        downstream_driver: StableTokenV2,
        downstream_script: impl Into<String>,
        downstream_manifest_path: impl Into<String>,
        downstream_manifest_root: ContentHash,
        no_claim: impl Into<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        let downstream_owner = downstream_owner.into();
        let downstream_script = downstream_script.into();
        let downstream_manifest_path = downstream_manifest_path.into();
        let no_claim = no_claim.into();
        validate_untyped_case_id(&downstream_owner)?;
        validate_source_path(&downstream_script)?;
        validate_source_path(&downstream_manifest_path)?;
        validate_untyped_case_id(&no_claim)?;
        let root = close_downstream_contribution_root(
            literal_expectation_oracle_root,
            semantic_input_root,
            budgets,
            schema_root,
            log_schema_root,
            &source_root,
            &build_root,
            &downstream_owner,
            &downstream_driver,
            &downstream_script,
            &downstream_manifest_path,
            downstream_manifest_root,
            &no_claim,
        )?;
        Ok(Self {
            literal_expectation_oracle_root,
            semantic_input_root,
            budgets,
            schema_root,
            log_schema_root,
            source_root,
            build_root,
            downstream_owner: downstream_owner.into_boxed_str(),
            downstream_driver,
            downstream_script: downstream_script.into_boxed_str(),
            downstream_manifest_path: downstream_manifest_path.into_boxed_str(),
            downstream_manifest_root,
            no_claim: no_claim.into_boxed_str(),
            root,
        })
    }

    #[must_use]
    pub const fn literal_expectation_oracle_root(&self) -> ContentHash {
        self.literal_expectation_oracle_root
    }

    #[must_use]
    pub const fn semantic_input_root(&self) -> ContentHash {
        self.semantic_input_root
    }

    #[must_use]
    pub const fn budgets(&self) -> BaseCoverageCloseContributionBudgetsV1 {
        self.budgets
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
    pub const fn source_root(&self) -> &SourceIdentityRootV2 {
        &self.source_root
    }

    #[must_use]
    pub const fn build_root(&self) -> &BuildIdentityRootV2 {
        &self.build_root
    }

    #[must_use]
    pub fn downstream_owner(&self) -> &str {
        &self.downstream_owner
    }

    #[must_use]
    pub const fn downstream_driver(&self) -> &StableTokenV2 {
        &self.downstream_driver
    }

    #[must_use]
    pub fn downstream_script(&self) -> &str {
        &self.downstream_script
    }

    #[must_use]
    pub fn downstream_manifest_path(&self) -> &str {
        &self.downstream_manifest_path
    }

    #[must_use]
    pub const fn downstream_manifest_root(&self) -> ContentHash {
        self.downstream_manifest_root
    }

    #[must_use]
    pub fn no_claim(&self) -> &str {
        &self.no_claim
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Declaration-side retained-artifact policy for an immutable contribution.
///
/// This policy permits only owner-reported relative paths in a later,
/// execution-owned envelope. It does not assert that an artifact exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseRetainedRelativeArtifactPolicyV1 {
    /// A later owner envelope may name only validated relative artifact paths.
    OwnerEnvelopeRelativePathsOnly = 1,
}

impl BaseCoverageCloseRetainedRelativeArtifactPolicyV1 {
    /// Exact nonzero policy code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Exact source-owned policy name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::OwnerEnvelopeRelativePathsOnly => "owner-envelope-relative-paths-only",
        }
    }
}

/// Exact designated observer for a result-free downstream contribution.
///
/// This value names who may observe a future execution and how that execution
/// must be routed. It contains no attempt, observation, result, artifact, or
/// execution-envelope root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseDeferredObservationContractV1 {
    route_id: StableTokenV2,
    semantic_consumer: Box<str>,
    execution_owner: Box<str>,
    driver_owner: Box<str>,
    posix_wrapper_owner: Box<str>,
    windows_wrapper_owner: Box<str>,
    capability_profile: BaseCoverageCloseCapabilityProfileV1,
    driver: StableTokenV2,
    posix_route: Box<str>,
    windows_route: Box<str>,
    case_manifest_path: Box<str>,
    case_manifest_root: ContentHash,
    deferred_reason_registry_root: BaseCoverageCloseDeferredReasonRegistryRootV1,
    deferred_reason: DeferredReasonV1,
    no_claim: StableTokenV2,
    root: BaseCoverageCloseDeferredObservationContractRootV1,
}

impl BaseCoverageCloseDeferredObservationContractV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "the observer contract deliberately preserves every distinct owner and route role"
    )]
    pub(crate) fn new(
        route_id: impl Into<String>,
        semantic_consumer: impl Into<String>,
        execution_owner: impl Into<String>,
        driver_owner: impl Into<String>,
        posix_wrapper_owner: impl Into<String>,
        windows_wrapper_owner: impl Into<String>,
        capability_profile: BaseCoverageCloseCapabilityProfileV1,
        driver: StableTokenV2,
        posix_route: impl Into<String>,
        windows_route: impl Into<String>,
        case_manifest_path: impl Into<String>,
        case_manifest_root: ContentHash,
        deferred_reason: DeferredReasonV1,
        no_claim: StableTokenV2,
    ) -> Result<Self, ConstructionErrorV2> {
        let route_id = StableTokenV2::new(route_id.into()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.deferred_observer.route_id",
                "one exact bounded stable route ID",
                0,
            )
        })?;
        let semantic_consumer = semantic_consumer.into();
        let execution_owner = execution_owner.into();
        let driver_owner = driver_owner.into();
        let posix_wrapper_owner = posix_wrapper_owner.into();
        let windows_wrapper_owner = windows_wrapper_owner.into();
        for owner in [
            semantic_consumer.as_str(),
            execution_owner.as_str(),
            driver_owner.as_str(),
            posix_wrapper_owner.as_str(),
            windows_wrapper_owner.as_str(),
        ] {
            validate_untyped_case_id(owner)?;
        }
        let posix_route = posix_route.into();
        let windows_route = windows_route.into();
        let case_manifest_path = case_manifest_path.into();
        validate_source_path(&posix_route)?;
        validate_source_path(&windows_route)?;
        validate_source_path(&case_manifest_path)?;
        if posix_route == windows_route
            || posix_route == case_manifest_path
            || windows_route == case_manifest_path
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.deferred_observer.route_paths",
                "distinct POSIX, native-Windows, and immutable case-manifest paths",
                route_id.as_str(),
            ));
        }
        if route_id.as_str() == RUNNER_V2_PHASE_ONE_CONTRACT_ROUTE_ID_V1
            && (semantic_consumer == RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1
                || execution_owner == RUNNER_V2_PHASE_ONE_CONTRACT_SOURCE_OWNER_V1)
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.deferred_observer.owner_roles",
                "the downstream semantic consumer and execution owner, never the upstream source owner",
                route_id.as_str(),
            ));
        }
        let deferred_registry = DeferredReasonRegistryV1::frozen()?;
        let root = close_deferred_observation_contract_root_v1(
            &route_id,
            &semantic_consumer,
            &execution_owner,
            &driver_owner,
            &posix_wrapper_owner,
            &windows_wrapper_owner,
            capability_profile,
            &driver,
            &posix_route,
            &windows_route,
            &case_manifest_path,
            case_manifest_root,
            deferred_registry.root(),
            deferred_reason,
            &no_claim,
        )?;
        Ok(Self {
            route_id,
            semantic_consumer: semantic_consumer.into_boxed_str(),
            execution_owner: execution_owner.into_boxed_str(),
            driver_owner: driver_owner.into_boxed_str(),
            posix_wrapper_owner: posix_wrapper_owner.into_boxed_str(),
            windows_wrapper_owner: windows_wrapper_owner.into_boxed_str(),
            capability_profile,
            driver,
            posix_route: posix_route.into_boxed_str(),
            windows_route: windows_route.into_boxed_str(),
            case_manifest_path: case_manifest_path.into_boxed_str(),
            case_manifest_root,
            deferred_reason_registry_root: deferred_registry.root(),
            deferred_reason,
            no_claim,
            root,
        })
    }

    /// Exact route obligation ID.
    #[must_use]
    pub const fn route_id(&self) -> &StableTokenV2 {
        &self.route_id
    }

    /// Leaf that semantically consumes the immutable payload.
    #[must_use]
    pub fn semantic_consumer(&self) -> &str {
        &self.semantic_consumer
    }

    /// Sole downstream execution owner.
    #[must_use]
    pub fn execution_owner(&self) -> &str {
        &self.execution_owner
    }

    /// Owner of the exact release driver.
    #[must_use]
    pub fn driver_owner(&self) -> &str {
        &self.driver_owner
    }

    /// Owner of the exact POSIX wrapper.
    #[must_use]
    pub fn posix_wrapper_owner(&self) -> &str {
        &self.posix_wrapper_owner
    }

    /// Owner of the exact native-Windows wrapper.
    #[must_use]
    pub fn windows_wrapper_owner(&self) -> &str {
        &self.windows_wrapper_owner
    }

    /// Declared downstream capability profile; this grants no capability.
    #[must_use]
    pub const fn capability_profile(&self) -> BaseCoverageCloseCapabilityProfileV1 {
        self.capability_profile
    }

    /// Exact downstream-owned release driver.
    #[must_use]
    pub const fn driver(&self) -> &StableTokenV2 {
        &self.driver
    }

    /// Exact downstream-owned POSIX wrapper path.
    #[must_use]
    pub fn posix_route(&self) -> &str {
        &self.posix_route
    }

    /// Exact downstream-owned native-Windows wrapper path.
    #[must_use]
    pub fn windows_route(&self) -> &str {
        &self.windows_route
    }

    /// Exact downstream-owned immutable case-manifest path.
    #[must_use]
    pub fn case_manifest_path(&self) -> &str {
        &self.case_manifest_path
    }

    /// Declaration root for the downstream-owned case manifest.
    ///
    /// The upstream leaf does not claim to have read or executed that
    /// downstream-owned file.
    #[must_use]
    pub const fn case_manifest_root(&self) -> ContentHash {
        self.case_manifest_root
    }

    /// Exact source-owned Deferred reason registry.
    #[must_use]
    pub const fn deferred_reason_registry_root(
        &self,
    ) -> BaseCoverageCloseDeferredReasonRegistryRootV1 {
        self.deferred_reason_registry_root
    }

    /// Sole exact Deferred reason.
    #[must_use]
    pub const fn deferred_reason(&self) -> DeferredReasonV1 {
        self.deferred_reason
    }

    /// Exact no-execution boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal observer-contract root.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseDeferredObservationContractRootV1 {
        self.root
    }
}

/// Additive V2 immutable downstream contribution.
///
/// The V1 type remains a compatibility declaration. V2 adds the explicit
/// observer/deferred contract, toolchain/target/feature identities, distinct
/// per-cell and per-shard budgets, expected partitions, and retained-relative
/// artifact policy. It intentionally does not contain its own Deferred
/// envelope, any actual observation, any result, or any execution claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseDownstreamContributionV2 {
    source_owner: Box<str>,
    payload_root: ContentHash,
    expected_partitions_root: ContentHash,
    per_cell_budgets: BaseCoverageCloseContributionBudgetsV1,
    shard_budgets: BaseCoverageCloseContributionBudgetsV1,
    schema_root: ContentHash,
    log_schema_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    toolchain_root: ToolchainIdentityRootV2,
    target_root: ContentHash,
    feature_set_root: ContentHash,
    retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
    observer_contract: BaseCoverageCloseDeferredObservationContractV1,
    no_claim: StableTokenV2,
    root: ContentHash,
}

impl BaseCoverageCloseDownstreamContributionV2 {
    #[allow(
        clippy::too_many_arguments,
        reason = "V2 intentionally binds every result-free contribution component without an open option bag"
    )]
    pub(crate) fn new(
        source_owner: impl Into<String>,
        payload_root: ContentHash,
        expected_partitions_root: ContentHash,
        per_cell_budgets: BaseCoverageCloseContributionBudgetsV1,
        shard_budgets: BaseCoverageCloseContributionBudgetsV1,
        schema_root: ContentHash,
        log_schema_root: ContentHash,
        source_root: SourceIdentityRootV2,
        build_root: BuildIdentityRootV2,
        toolchain_root: ToolchainIdentityRootV2,
        target_root: ContentHash,
        feature_set_root: ContentHash,
        retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
        observer_contract: BaseCoverageCloseDeferredObservationContractV1,
        no_claim: StableTokenV2,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_owner = source_owner.into();
        validate_untyped_case_id(&source_owner)?;
        if source_owner == observer_contract.execution_owner()
            || source_owner == observer_contract.semantic_consumer()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_v2.owner_roles",
                "an upstream source owner distinct from downstream semantic and execution ownership",
                source_owner,
            ));
        }
        for budgets in [per_cell_budgets, shard_budgets] {
            if budgets.resolved().profile()
                != BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution
            {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.contribution_v2.budget_profile",
                    "the exact downstream source-contribution profile for both budget scopes",
                    budgets.resolved().profile().code(),
                ));
            }
        }
        if per_cell_budgets.max_child_processes() > shard_budgets.max_child_processes()
            || per_cell_budgets.max_parallel_children() > shard_budgets.max_parallel_children()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_v2.budget_scope",
                "per-cell child-process limits no greater than per-shard limits",
                per_cell_budgets.max_child_processes(),
            ));
        }
        let roots = [
            payload_root,
            expected_partitions_root,
            schema_root,
            log_schema_root,
            target_root,
            feature_set_root,
        ];
        if roots.iter().copied().collect::<BTreeSet<_>>().len() != roots.len() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution_v2.semantic_roots",
                "distinct payload, partition, schema, log, target, and feature roots",
                roots.len(),
            ));
        }
        let root = close_downstream_contribution_root_v2(
            &source_owner,
            payload_root,
            expected_partitions_root,
            per_cell_budgets,
            shard_budgets,
            schema_root,
            log_schema_root,
            &source_root,
            &build_root,
            &toolchain_root,
            target_root,
            feature_set_root,
            retained_artifact_policy,
            &observer_contract,
            &no_claim,
        )?;
        Ok(Self {
            source_owner: source_owner.into_boxed_str(),
            payload_root,
            expected_partitions_root,
            per_cell_budgets,
            shard_budgets,
            schema_root,
            log_schema_root,
            source_root,
            build_root,
            toolchain_root,
            target_root,
            feature_set_root,
            retained_artifact_policy,
            observer_contract,
            no_claim,
            root,
        })
    }

    /// Upstream source owner of the result-free payload.
    #[must_use]
    pub fn source_owner(&self) -> &str {
        &self.source_owner
    }

    /// Exact immutable result-free payload root.
    #[must_use]
    pub const fn payload_root(&self) -> ContentHash {
        self.payload_root
    }

    /// Exact declaration-side expected-partition root.
    #[must_use]
    pub const fn expected_partitions_root(&self) -> ContentHash {
        self.expected_partitions_root
    }

    /// Exact per-cell contribution budget.
    #[must_use]
    pub const fn per_cell_budgets(&self) -> BaseCoverageCloseContributionBudgetsV1 {
        self.per_cell_budgets
    }

    /// Exact per-shard contribution budget.
    #[must_use]
    pub const fn shard_budgets(&self) -> BaseCoverageCloseContributionBudgetsV1 {
        self.shard_budgets
    }

    /// Exact schema-declaration root.
    #[must_use]
    pub const fn schema_root(&self) -> ContentHash {
        self.schema_root
    }

    /// Exact detailed-log schema declaration root.
    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    /// Exact compatible source identity.
    #[must_use]
    pub const fn source_root(&self) -> &SourceIdentityRootV2 {
        &self.source_root
    }

    /// Exact declared build identity.
    #[must_use]
    pub const fn build_root(&self) -> &BuildIdentityRootV2 {
        &self.build_root
    }

    /// Exact declared toolchain identity.
    #[must_use]
    pub const fn toolchain_root(&self) -> &ToolchainIdentityRootV2 {
        &self.toolchain_root
    }

    /// Exact target-contract root.
    #[must_use]
    pub const fn target_root(&self) -> ContentHash {
        self.target_root
    }

    /// Exact feature-contract root.
    #[must_use]
    pub const fn feature_set_root(&self) -> ContentHash {
        self.feature_set_root
    }

    /// Declaration-side retained-relative-artifact policy.
    #[must_use]
    pub const fn retained_artifact_policy(
        &self,
    ) -> BaseCoverageCloseRetainedRelativeArtifactPolicyV1 {
        self.retained_artifact_policy
    }

    /// Exact designated observer contract.
    #[must_use]
    pub const fn observer_contract(&self) -> &BaseCoverageCloseDeferredObservationContractV1 {
        &self.observer_contract
    }

    /// Exact no-execution/no-authority boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Domain-separated V2 declaration root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// A contribution is structurally incapable of carrying execution actuals.
    #[must_use]
    pub const fn execution_actual_field_count(&self) -> u16 {
        0
    }
}

/// Separate Deferred evidence envelope for one immutable V2 contribution.
///
/// Construction accepts only the typed contribution and its exact embedded
/// observer contract. There is no open disposition or optional-actual field
/// surface, so this type cannot fabricate Observed/NotObserved execution data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseDeferredEvidenceEnvelopeV1 {
    evidence_kind: BaseCoverageCloseEvidenceKindV1,
    disposition: RuntimeObservationDispositionV1,
    disposition_root: BaseCoverageCloseRuntimeObservationDispositionRootV1,
    contribution_root: ContentHash,
    observer_contract_root: BaseCoverageCloseDeferredObservationContractRootV1,
    deferred_reason_registry_root: BaseCoverageCloseDeferredReasonRegistryRootV1,
    deferred_reason: DeferredReasonV1,
    retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
    no_claim: StableTokenV2,
    root: BaseCoverageCloseEvidenceEnvelopeRootV1,
}

impl BaseCoverageCloseDeferredEvidenceEnvelopeV1 {
    /// Construct the only legal declaration-side Deferred envelope shape.
    pub(crate) fn new(
        contribution: &BaseCoverageCloseDownstreamContributionV2,
        no_claim: StableTokenV2,
    ) -> Result<Self, ConstructionErrorV2> {
        let observer = contribution.observer_contract();
        if observer.deferred_reason()
            != DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.deferred_envelope.reason",
                "immutable-contribution-awaits-designated-release-owner",
                observer.deferred_reason().code(),
            ));
        }
        let disposition = RuntimeObservationDispositionV1::Deferred;
        let disposition_root = disposition.root()?;
        let root = close_deferred_evidence_envelope_root_v1(
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution,
            disposition,
            disposition_root,
            contribution.root(),
            observer.root(),
            observer.deferred_reason_registry_root(),
            observer.deferred_reason(),
            contribution.retained_artifact_policy(),
            &no_claim,
        )?;
        Ok(Self {
            evidence_kind: BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution,
            disposition,
            disposition_root,
            contribution_root: contribution.root(),
            observer_contract_root: observer.root(),
            deferred_reason_registry_root: observer.deferred_reason_registry_root(),
            deferred_reason: observer.deferred_reason(),
            retained_artifact_policy: contribution.retained_artifact_policy(),
            no_claim,
            root,
        })
    }

    /// Exact immutable-contribution evidence kind.
    #[must_use]
    pub const fn evidence_kind(&self) -> BaseCoverageCloseEvidenceKindV1 {
        self.evidence_kind
    }

    /// Exact Deferred runtime-observation disposition.
    #[must_use]
    pub const fn disposition(&self) -> RuntimeObservationDispositionV1 {
        self.disposition
    }

    /// Nominal root of the exact Deferred disposition.
    #[must_use]
    pub const fn disposition_root(&self) -> BaseCoverageCloseRuntimeObservationDispositionRootV1 {
        self.disposition_root
    }

    /// Exact immutable contribution root.
    #[must_use]
    pub const fn contribution_root(&self) -> ContentHash {
        self.contribution_root
    }

    /// Exact designated observer-contract root.
    #[must_use]
    pub const fn observer_contract_root(
        &self,
    ) -> BaseCoverageCloseDeferredObservationContractRootV1 {
        self.observer_contract_root
    }

    /// Exact Deferred reason-registry root.
    #[must_use]
    pub const fn deferred_reason_registry_root(
        &self,
    ) -> BaseCoverageCloseDeferredReasonRegistryRootV1 {
        self.deferred_reason_registry_root
    }

    /// Sole exact Deferred reason.
    #[must_use]
    pub const fn deferred_reason(&self) -> DeferredReasonV1 {
        self.deferred_reason
    }

    /// Declaration-side retained-relative-artifact policy.
    #[must_use]
    pub const fn retained_artifact_policy(
        &self,
    ) -> BaseCoverageCloseRetainedRelativeArtifactPolicyV1 {
        self.retained_artifact_policy
    }

    /// Exact no-execution boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal Deferred evidence-envelope root.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseEvidenceEnvelopeRootV1 {
        self.root
    }

    /// A Deferred contribution envelope carries no execution actuals.
    #[must_use]
    pub const fn execution_actual_field_count(&self) -> u16 {
        0
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "the nominal observer root binds every distinct ownership and route role"
)]
fn close_deferred_observation_contract_root_v1(
    route_id: &StableTokenV2,
    semantic_consumer: &str,
    execution_owner: &str,
    driver_owner: &str,
    posix_wrapper_owner: &str,
    windows_wrapper_owner: &str,
    capability_profile: BaseCoverageCloseCapabilityProfileV1,
    driver: &StableTokenV2,
    posix_route: &str,
    windows_route: &str,
    case_manifest_path: &str,
    case_manifest_root: ContentHash,
    deferred_reason_registry_root: BaseCoverageCloseDeferredReasonRegistryRootV1,
    deferred_reason: DeferredReasonV1,
    no_claim: &StableTokenV2,
) -> Result<BaseCoverageCloseDeferredObservationContractRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEDEFERREDOBSERVER\x01", 8 * 1024)?;
    frame.push_u16("observer.api_generation", 2)?;
    frame.push_u16("observer.wire_version", 1)?;
    frame.push_str("observer.route_id", route_id.as_str())?;
    frame.push_str("observer.semantic_consumer", semantic_consumer)?;
    frame.push_str("observer.execution_owner", execution_owner)?;
    frame.push_str("observer.driver_owner", driver_owner)?;
    frame.push_str("observer.posix_wrapper_owner", posix_wrapper_owner)?;
    frame.push_str("observer.windows_wrapper_owner", windows_wrapper_owner)?;
    frame.push_u16("observer.capability_profile", capability_profile.code())?;
    frame.push_str(
        "observer.capability_profile_id",
        capability_profile.stable_id(),
    )?;
    frame.push_str("observer.driver", driver.as_str())?;
    frame.push_str("observer.posix_route", posix_route)?;
    frame.push_str("observer.windows_route", windows_route)?;
    frame.push_str("observer.case_manifest_path", case_manifest_path)?;
    frame.push_bytes("observer.case_manifest_root", case_manifest_root.as_bytes())?;
    frame.push_bytes(
        "observer.deferred_reason_registry_root",
        deferred_reason_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_u16("observer.deferred_reason", deferred_reason.code())?;
    frame.push_str(
        "observer.deferred_reason_name",
        deferred_reason.descriptor().name(),
    )?;
    frame.push_str("observer.no_claim", no_claim.as_str())?;
    Ok(
        BaseCoverageCloseDeferredObservationContractRootV1::from_content_hash(
            frame.root(BaseCoverageCloseDeferredObservationContractRootV1::DESCRIPTOR.domain()),
        ),
    )
}

#[allow(
    clippy::too_many_arguments,
    reason = "the V2 contribution root intentionally binds every declaration component and excludes its envelope"
)]
fn close_downstream_contribution_root_v2(
    source_owner: &str,
    payload_root: ContentHash,
    expected_partitions_root: ContentHash,
    per_cell_budgets: BaseCoverageCloseContributionBudgetsV1,
    shard_budgets: BaseCoverageCloseContributionBudgetsV1,
    schema_root: ContentHash,
    log_schema_root: ContentHash,
    source_root: &SourceIdentityRootV2,
    build_root: &BuildIdentityRootV2,
    toolchain_root: &ToolchainIdentityRootV2,
    target_root: ContentHash,
    feature_set_root: ContentHash,
    retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
    observer_contract: &BaseCoverageCloseDeferredObservationContractV1,
    no_claim: &StableTokenV2,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECONTRIBUTION\x02", 32 * 1024)?;
    frame.push_u16("contribution.api_generation", 2)?;
    frame.push_u16("contribution.wire_version", 1)?;
    frame.push_u16("contribution.wire_predecessor_policy", 1)?;
    frame.push_str("contribution.source_owner", source_owner)?;
    frame.push_str(
        "contribution.route_id",
        observer_contract.route_id().as_str(),
    )?;
    frame.push_bytes("contribution.payload_root", payload_root.as_bytes())?;
    frame.push_bytes(
        "contribution.expected_partitions_root",
        expected_partitions_root.as_bytes(),
    )?;
    push_contribution_budget_scope_v2(&mut frame, 1, per_cell_budgets)?;
    push_contribution_budget_scope_v2(&mut frame, 2, shard_budgets)?;
    frame.push_bytes("contribution.schema_root", schema_root.as_bytes())?;
    frame.push_bytes("contribution.log_schema_root", log_schema_root.as_bytes())?;
    frame.push_bytes("contribution.source_root", source_root.digest().bytes())?;
    frame.push_bytes("contribution.build_root", build_root.digest().bytes())?;
    frame.push_bytes(
        "contribution.toolchain_root",
        toolchain_root.digest().bytes(),
    )?;
    frame.push_bytes("contribution.target_root", target_root.as_bytes())?;
    frame.push_bytes("contribution.feature_set_root", feature_set_root.as_bytes())?;
    frame.push_u16(
        "contribution.retained_artifact_policy",
        retained_artifact_policy.code(),
    )?;
    frame.push_str(
        "contribution.retained_artifact_policy_name",
        retained_artifact_policy.stable_name(),
    )?;
    frame.push_bytes(
        "contribution.observer_contract_root",
        observer_contract.root().content_hash().as_bytes(),
    )?;
    frame.push_bytes(
        "contribution.deferred_reason_registry_root",
        observer_contract
            .deferred_reason_registry_root()
            .content_hash()
            .as_bytes(),
    )?;
    frame.push_u16(
        "contribution.deferred_reason",
        observer_contract.deferred_reason().code(),
    )?;
    frame.push_str("contribution.no_claim", no_claim.as_str())?;
    Ok(frame.root(BASE_COVERAGE_CLOSE_DOWNSTREAM_CONTRIBUTION_DOMAIN_V2))
}

fn push_contribution_budget_scope_v2(
    frame: &mut CanonicalFrameV1,
    scope_code: u16,
    budgets: BaseCoverageCloseContributionBudgetsV1,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16("contribution.budget.scope", scope_code)?;
    frame.push_u16(
        "contribution.budget.profile",
        budgets.resolved().profile().code(),
    )?;
    for row in budgets.resolved().rows() {
        frame.push_u16("contribution.budget.axis", row.axis().code())?;
        frame.push_u16("contribution.budget.width", row.hard().width().code())?;
        push_contribution_budget_value_v2(frame, row.hard())?;
        push_contribution_budget_value_v2(frame, row.soft())?;
        frame.push_u16("contribution.budget.unit", row.unit().unit().tag())?;
        frame.push_presence(
            "contribution.budget.unit_registry.present",
            row.unit().registry_identity().is_some(),
        )?;
        if let Some(registry_identity) = row.unit().registry_identity() {
            frame.push_bytes(
                "contribution.budget.unit_registry",
                registry_identity.as_bytes(),
            )?;
        }
    }
    frame.push_u32(
        "contribution.budget.max_child_processes",
        budgets.max_child_processes(),
    )?;
    frame.push_u32(
        "contribution.budget.max_parallel_children",
        budgets.max_parallel_children(),
    )?;
    Ok(())
}

fn push_contribution_budget_value_v2(
    frame: &mut CanonicalFrameV1,
    value: BaseCoverageCloseBudgetValueV1,
) -> Result<(), ConstructionErrorV2> {
    match value {
        BaseCoverageCloseBudgetValueV1::U32(value) => {
            frame.push_u32("contribution.budget.value_u32", value)
        }
        BaseCoverageCloseBudgetValueV1::U64(value) => {
            frame.push_u64("contribution.budget.value_u64", value)
        }
        BaseCoverageCloseBudgetValueV1::U128(value) => {
            frame.push_u128("contribution.budget.value_u128", value)
        }
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "the Deferred envelope root binds the exact closed contribution-only evidence shape"
)]
fn close_deferred_evidence_envelope_root_v1(
    evidence_kind: BaseCoverageCloseEvidenceKindV1,
    disposition: RuntimeObservationDispositionV1,
    disposition_root: BaseCoverageCloseRuntimeObservationDispositionRootV1,
    contribution_root: ContentHash,
    observer_contract_root: BaseCoverageCloseDeferredObservationContractRootV1,
    deferred_reason_registry_root: BaseCoverageCloseDeferredReasonRegistryRootV1,
    deferred_reason: DeferredReasonV1,
    retained_artifact_policy: BaseCoverageCloseRetainedRelativeArtifactPolicyV1,
    no_claim: &StableTokenV2,
) -> Result<BaseCoverageCloseEvidenceEnvelopeRootV1, ConstructionErrorV2> {
    if evidence_kind != BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
        || disposition != RuntimeObservationDispositionV1::Deferred
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.deferred_envelope.kind_disposition",
            "immutable downstream contribution plus Deferred",
            disposition.code(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEENVELOPE\x01", 8 * 1024)?;
    frame.push_u16("envelope.api_generation", 2)?;
    frame.push_u16("envelope.wire_version", 1)?;
    frame.push_u16("envelope.evidence_kind", evidence_kind.code())?;
    frame.push_bytes("envelope.payload_root", contribution_root.as_bytes())?;
    frame.push_u16("envelope.disposition", disposition.code())?;
    frame.push_bytes(
        "envelope.disposition_root",
        disposition_root.content_hash().as_bytes(),
    )?;
    frame.push_bytes("envelope.contribution_root", contribution_root.as_bytes())?;
    frame.push_bytes(
        "envelope.observer_contract_root",
        observer_contract_root.content_hash().as_bytes(),
    )?;
    frame.push_bytes(
        "envelope.deferred_reason_registry_root",
        deferred_reason_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_u16("envelope.deferred_reason", deferred_reason.code())?;
    frame.push_u16(
        "envelope.retained_artifact_policy",
        retained_artifact_policy.code(),
    )?;
    frame.push_str("envelope.no_claim", no_claim.as_str())?;
    Ok(BaseCoverageCloseEvidenceEnvelopeRootV1::from_content_hash(
        frame.root(BaseCoverageCloseEvidenceEnvelopeRootV1::DESCRIPTOR.domain()),
    ))
}

/// One result-free AC53 close declaration bound one-to-one to a source case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseCellDeclarationV1 {
    source_ordinal: u32,
    source_case_id: Box<str>,
    source_class: BaseCoverageManifestClassV1,
    source_path: Box<str>,
    group: BaseCoverageCloseGroupV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    downstream_contribution: Option<BaseCoverageCloseDownstreamContributionV1>,
    five_explicits: BaseCoverageCloseFiveExplicitsV1,
}

impl BaseCoverageCloseCellDeclarationV1 {
    /// Construct a bounded result-free declaration.
    #[allow(
        clippy::too_many_arguments,
        reason = "AC53 requires every classification axis explicitly"
    )]
    pub fn new(
        source_ordinal: u32,
        source_case_id: impl Into<String>,
        source_class: BaseCoverageManifestClassV1,
        source_path: impl Into<String>,
        group: BaseCoverageCloseGroupV1,
        facet: BaseCoverageCloseFacetV1,
        execution_scope: BaseCoverageCloseExecutionScopeV1,
        partition: BaseCoverageClosePartitionV1,
        expected_decision: BaseCoverageCloseDecisionV1,
        expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new_with_downstream_contribution(
            source_ordinal,
            source_case_id,
            source_class,
            source_path,
            group,
            facet,
            execution_scope,
            partition,
            expected_decision,
            expected_reason,
            None,
        )
    }

    /// Construct a declaration with the result-free downstream payload that is
    /// mandatory exactly for `ImmutableDownstreamContribution` scope.
    #[allow(
        clippy::too_many_arguments,
        reason = "AC53 requires every classification axis explicitly"
    )]
    pub fn new_with_downstream_contribution(
        source_ordinal: u32,
        source_case_id: impl Into<String>,
        source_class: BaseCoverageManifestClassV1,
        source_path: impl Into<String>,
        group: BaseCoverageCloseGroupV1,
        facet: BaseCoverageCloseFacetV1,
        execution_scope: BaseCoverageCloseExecutionScopeV1,
        partition: BaseCoverageClosePartitionV1,
        expected_decision: BaseCoverageCloseDecisionV1,
        expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
        downstream_contribution: Option<BaseCoverageCloseDownstreamContributionV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_case_id = source_case_id.into();
        let source_path = source_path.into();
        let five_explicits = frozen_close_five_explicits_v1(
            source_ordinal,
            &source_case_id,
            source_class,
            &source_path,
            facet,
            execution_scope,
            downstream_contribution.as_ref(),
        )?;
        Self::new_with_five_explicits(
            source_ordinal,
            source_case_id,
            source_class,
            source_path,
            group,
            facet,
            execution_scope,
            partition,
            expected_decision,
            expected_reason,
            downstream_contribution,
            five_explicits,
        )
    }

    /// Construct one declaration with a caller-presented Five Explicits
    /// record. Exact reconstruction uses this surface to prove that changing
    /// any explicit changes and then mismatches the stable cell.
    #[allow(
        clippy::too_many_arguments,
        reason = "the controlling close contract retains every classification and explicit axis"
    )]
    pub fn new_with_five_explicits(
        source_ordinal: u32,
        source_case_id: impl Into<String>,
        source_class: BaseCoverageManifestClassV1,
        source_path: impl Into<String>,
        group: BaseCoverageCloseGroupV1,
        facet: BaseCoverageCloseFacetV1,
        execution_scope: BaseCoverageCloseExecutionScopeV1,
        partition: BaseCoverageClosePartitionV1,
        expected_decision: BaseCoverageCloseDecisionV1,
        expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
        downstream_contribution: Option<BaseCoverageCloseDownstreamContributionV1>,
        five_explicits: BaseCoverageCloseFiveExplicitsV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_case_id = source_case_id.into();
        let source_path = source_path.into();
        if source_ordinal == 0 {
            return Err(refusal(
                ConstructionErrorKindV2::Zero,
                "coverage.close.source_ordinal",
                "a nonzero one-based source-manifest ordinal",
                source_ordinal,
            ));
        }
        validate_case_id(source_class, &source_case_id)?;
        validate_source_path(&source_path)?;
        if facet.group() != group {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.group",
                "the sole group owning the presented facet",
                group.stable_name(),
            ));
        }
        validate_close_partition_shape(partition, expected_decision, expected_reason)?;
        validate_close_reason_scope(expected_reason, execution_scope)?;
        let contribution_required =
            execution_scope == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution;
        if contribution_required != downstream_contribution.is_some() {
            return Err(refusal(
                if contribution_required {
                    ConstructionErrorKindV2::Missing
                } else {
                    ConstructionErrorKindV2::Unexpected
                },
                "coverage.close.downstream_contribution",
                "present exactly for immutable downstream contribution scope",
                downstream_contribution.is_some(),
            ));
        }
        validate_close_five_explicits_for_declaration(
            &source_case_id,
            execution_scope,
            downstream_contribution.as_ref(),
            &five_explicits,
        )?;
        Ok(Self {
            source_ordinal,
            source_case_id: source_case_id.into_boxed_str(),
            source_class,
            source_path: source_path.into_boxed_str(),
            group,
            facet,
            execution_scope,
            partition,
            expected_decision,
            expected_reason,
            downstream_contribution,
            five_explicits,
        })
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
    pub const fn source_class(&self) -> BaseCoverageManifestClassV1 {
        self.source_class
    }

    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
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
    pub const fn downstream_contribution(
        &self,
    ) -> Option<&BaseCoverageCloseDownstreamContributionV1> {
        self.downstream_contribution.as_ref()
    }

    #[must_use]
    pub const fn five_explicits(&self) -> &BaseCoverageCloseFiveExplicitsV1 {
        &self.five_explicits
    }
}

/// One admitted immutable AC53 close-manifest cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseManifestCellV1 {
    declaration: BaseCoverageCloseCellDeclarationV1,
    root: ContentHash,
}

impl BaseCoverageCloseManifestCellV1 {
    #[must_use]
    pub const fn declaration(&self) -> &BaseCoverageCloseCellDeclarationV1 {
        &self.declaration
    }

    #[must_use]
    pub const fn source_ordinal(&self) -> u32 {
        self.declaration.source_ordinal
    }

    #[must_use]
    pub fn source_case_id(&self) -> &str {
        self.declaration.source_case_id()
    }

    #[must_use]
    pub const fn source_class(&self) -> BaseCoverageManifestClassV1 {
        self.declaration.source_class
    }

    #[must_use]
    pub fn source_path(&self) -> &str {
        self.declaration.source_path()
    }

    #[must_use]
    pub const fn group(&self) -> BaseCoverageCloseGroupV1 {
        self.declaration.group
    }

    #[must_use]
    pub const fn facet(&self) -> BaseCoverageCloseFacetV1 {
        self.declaration.facet
    }

    #[must_use]
    pub const fn execution_scope(&self) -> BaseCoverageCloseExecutionScopeV1 {
        self.declaration.execution_scope
    }

    #[must_use]
    pub const fn partition(&self) -> BaseCoverageClosePartitionV1 {
        self.declaration.partition
    }

    #[must_use]
    pub const fn expected_decision(&self) -> BaseCoverageCloseDecisionV1 {
        self.declaration.expected_decision
    }

    #[must_use]
    pub const fn expected_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.declaration.expected_reason
    }

    #[must_use]
    pub const fn downstream_contribution(
        &self,
    ) -> Option<&BaseCoverageCloseDownstreamContributionV1> {
        self.declaration.downstream_contribution()
    }

    #[must_use]
    pub const fn five_explicits(&self) -> &BaseCoverageCloseFiveExplicitsV1 {
        self.declaration.five_explicits()
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Source-authoritative AC53 full-set close manifest.
///
/// This type deliberately has no subset-selection API. Construction requires
/// the exact frozen base and exact E2E/logging/source/downstream extension
/// sequence, then materializes one result-free close cell per source cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseManifestV1 {
    source_manifest_root: ContentHash,
    reason_registry_root: ContentHash,
    cells: Box<[BaseCoverageCloseManifestCellV1]>,
    root: ContentHash,
}

impl BaseCoverageCloseManifestV1 {
    /// Construct the independently frozen complete source and close manifests.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let source = frozen_full_source_manifest_v1()?;
        Self::reconstruct_full(&source)
    }

    /// Reconstruct close authority from the exact complete source manifest.
    ///
    /// A base-only manifest, an alternate extension set, or any reordered,
    /// reclassified, stale, or path-mutated source row refuses.
    pub fn reconstruct_full(source: &BaseCoverageManifestV1) -> Result<Self, ConstructionErrorV2> {
        let expected_source = frozen_full_source_manifest_v1()?;
        validate_exact_full_source_manifest(&expected_source, source)?;
        let declarations = source
            .cases()
            .iter()
            .map(classify_close_case_v1)
            .collect::<Result<Vec<_>, _>>()?;
        close_manifest_from_declarations(source, declarations)
    }

    /// Reconstruct from a caller-presented exact full declaration sequence.
    ///
    /// This checks every classification field against the independently
    /// materialized oracle. It never accepts an empty or partial declaration
    /// slice and never derives semantics from execution results.
    pub fn reconstruct_exact_full(
        source: &BaseCoverageManifestV1,
        presented: &[BaseCoverageCloseCellDeclarationV1],
    ) -> Result<Self, ConstructionErrorV2> {
        let expected = Self::reconstruct_full(source)?;
        validate_exact_close_declaration_sequence(&expected, presented)?;
        Ok(expected)
    }

    #[must_use]
    pub const fn source_manifest_root(&self) -> ContentHash {
        self.source_manifest_root
    }

    #[must_use]
    pub const fn reason_registry_root(&self) -> ContentHash {
        self.reason_registry_root
    }

    #[must_use]
    pub fn cells(&self) -> &[BaseCoverageCloseManifestCellV1] {
        &self.cells
    }

    #[must_use]
    pub fn cell(&self, source_case_id: &str) -> Option<&BaseCoverageCloseManifestCellV1> {
        self.cells
            .iter()
            .find(|cell| cell.source_case_id() == source_case_id)
    }

    #[must_use]
    pub fn group_count(&self, group: BaseCoverageCloseGroupV1) -> usize {
        self.cells
            .iter()
            .filter(|cell| cell.group() == group)
            .count()
    }

    #[must_use]
    pub fn facet_count(&self, facet: BaseCoverageCloseFacetV1) -> usize {
        self.cells
            .iter()
            .filter(|cell| cell.facet() == facet)
            .count()
    }

    #[must_use]
    pub fn applicable_facet_count(&self, facet: BaseCoverageCloseFacetV1) -> usize {
        self.cells
            .iter()
            .filter(|cell| {
                cell.facet() == facet
                    && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
            })
            .count()
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Typed evidence identity for one presented close result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseEvidenceKindV1 {
    OwnedHarnessExecution = 1,
    InProcessProjectionExecution = 2,
    ImmutableDownstreamContribution = 3,
    ApplicabilityDeclaration = 4,
}

impl BaseCoverageCloseEvidenceKindV1 {
    pub const ALL: [Self; 4] = [
        Self::OwnedHarnessExecution,
        Self::InProcessProjectionExecution,
        Self::ImmutableDownstreamContribution,
        Self::ApplicabilityDeclaration,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Content-bound, non-authoritative evidence attached to one close result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseResultEvidenceV1 {
    kind: BaseCoverageCloseEvidenceKindV1,
    root: ContentHash,
    retained_artifact: Option<Box<str>>,
}

impl BaseCoverageCloseResultEvidenceV1 {
    /// Construct one typed evidence identity with an optional safe relative
    /// retained-artifact reference.
    pub fn new(
        kind: BaseCoverageCloseEvidenceKindV1,
        root: ContentHash,
        retained_artifact: Option<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        if let Some(path) = retained_artifact.as_deref() {
            validate_source_path(path)?;
        }
        if matches!(
            kind,
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
                | BaseCoverageCloseEvidenceKindV1::ApplicabilityDeclaration
        ) && retained_artifact.is_some()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close_result.retained_artifact",
                "no fabricated retained artifact for a declaration-only evidence kind",
                kind.code(),
            ));
        }
        Ok(Self {
            kind,
            root,
            retained_artifact: retained_artifact.map(String::into_boxed_str),
        })
    }

    pub fn owned_harness_execution(
        root: ContentHash,
        retained_artifact: Option<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            BaseCoverageCloseEvidenceKindV1::OwnedHarnessExecution,
            root,
            retained_artifact,
        )
    }

    pub fn in_process_projection_execution(
        root: ContentHash,
        retained_artifact: Option<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            BaseCoverageCloseEvidenceKindV1::InProcessProjectionExecution,
            root,
            retained_artifact,
        )
    }

    pub fn immutable_downstream_contribution(
        contribution: &BaseCoverageCloseDownstreamContributionV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution,
            contribution.root(),
            None,
        )
    }

    pub fn applicability_declaration(
        manifest: &BaseCoverageCloseManifestV1,
        reason: BaseCoverageCloseReasonCodeV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            BaseCoverageCloseEvidenceKindV1::ApplicabilityDeclaration,
            close_applicability_evidence_root(manifest.reason_registry_root(), reason)?,
            None,
        )
    }

    #[must_use]
    pub const fn kind(&self) -> BaseCoverageCloseEvidenceKindV1 {
        self.kind
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    #[must_use]
    pub fn retained_artifact(&self) -> Option<&str> {
        self.retained_artifact.as_deref()
    }
}

/// Terminal accounting status for one caller-presented close result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseCoverageCloseResultStatusV1 {
    Matched = 1,
    UnexpectedMismatch = 2,
    ExecutionFailure = 3,
    UnexplainedSkip = 4,
}

impl BaseCoverageCloseResultStatusV1 {
    pub const ALL: [Self; 4] = [
        Self::Matched,
        Self::UnexpectedMismatch,
        Self::ExecutionFailure,
        Self::UnexplainedSkip,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// One full-set result carrying its immutable expected classification and
/// caller-observed terminal classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageClosePresentedResultV1 {
    close_manifest_root: ContentHash,
    cell_root: ContentHash,
    source_case_id: Box<str>,
    group: BaseCoverageCloseGroupV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    status: BaseCoverageCloseResultStatusV1,
    observed_decision: Option<BaseCoverageCloseDecisionV1>,
    observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
    evidence: BaseCoverageCloseResultEvidenceV1,
    root: ContentHash,
}

impl BaseCoverageClosePresentedResultV1 {
    /// Construct one explicitly classified result.
    ///
    /// A `Matched` row must carry exact expected decision and reason data.
    /// An `UnexpectedMismatch` row must actually differ. Execution failures
    /// and unexplained skips carry no fabricated observed classification.
    #[allow(
        clippy::too_many_arguments,
        reason = "AC53 result rows retain every classification axis"
    )]
    pub fn new(
        close_manifest_root: ContentHash,
        cell_root: ContentHash,
        source_case_id: impl Into<String>,
        group: BaseCoverageCloseGroupV1,
        facet: BaseCoverageCloseFacetV1,
        execution_scope: BaseCoverageCloseExecutionScopeV1,
        partition: BaseCoverageClosePartitionV1,
        expected_decision: BaseCoverageCloseDecisionV1,
        expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
        status: BaseCoverageCloseResultStatusV1,
        observed_decision: Option<BaseCoverageCloseDecisionV1>,
        observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
        evidence: BaseCoverageCloseResultEvidenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_case_id = source_case_id.into();
        validate_untyped_case_id(&source_case_id)?;
        if facet.group() != group {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close_result.group",
                "the sole group owning the presented facet",
                group.stable_name(),
            ));
        }
        validate_close_partition_shape(partition, expected_decision, expected_reason)?;
        validate_close_reason_scope(expected_reason, execution_scope)?;
        validate_close_evidence_scope_shape(execution_scope, &evidence)?;
        match status {
            BaseCoverageCloseResultStatusV1::Matched => {
                if observed_decision != Some(expected_decision)
                    || observed_reason != expected_reason
                {
                    return Err(refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "coverage.close_result.matched_observation",
                        "the exact expected decision and registered reason",
                        source_case_id.as_str(),
                    ));
                }
            }
            BaseCoverageCloseResultStatusV1::UnexpectedMismatch => {
                if observed_decision.is_none() {
                    return Err(refusal(
                        ConstructionErrorKindV2::Missing,
                        "coverage.close_result.observed_decision",
                        "one exact observed decision for an unexpected mismatch",
                        source_case_id.as_str(),
                    ));
                }
                if observed_decision == Some(expected_decision)
                    && observed_reason == expected_reason
                {
                    return Err(refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "coverage.close_result.unexpected_mismatch",
                        "an observed decision or reason different from the expected cell",
                        source_case_id.as_str(),
                    ));
                }
            }
            BaseCoverageCloseResultStatusV1::ExecutionFailure
            | BaseCoverageCloseResultStatusV1::UnexplainedSkip => {
                if observed_decision.is_some() || observed_reason.is_some() {
                    return Err(refusal(
                        ConstructionErrorKindV2::Unexpected,
                        "coverage.close_result.observed_classification",
                        "no fabricated decision or reason after failure or skip",
                        source_case_id.as_str(),
                    ));
                }
            }
        }
        let root = close_presented_result_root(
            close_manifest_root,
            cell_root,
            &source_case_id,
            group,
            facet,
            execution_scope,
            partition,
            expected_decision,
            expected_reason,
            status,
            observed_decision,
            observed_reason,
            &evidence,
        )?;
        Ok(Self {
            close_manifest_root,
            cell_root,
            source_case_id: source_case_id.into_boxed_str(),
            group,
            facet,
            execution_scope,
            partition,
            expected_decision,
            expected_reason,
            status,
            observed_decision,
            observed_reason,
            evidence,
            root,
        })
    }

    /// Exact matched result fixture for one immutable close cell.
    pub fn matched(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        evidence: BaseCoverageCloseResultEvidenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            manifest.root,
            cell.root,
            cell.source_case_id(),
            cell.group(),
            cell.facet(),
            cell.execution_scope(),
            cell.partition(),
            cell.expected_decision(),
            cell.expected_reason(),
            BaseCoverageCloseResultStatusV1::Matched,
            Some(cell.expected_decision()),
            cell.expected_reason(),
            evidence,
        )
    }

    /// Explicit unexpected semantic mismatch for one immutable close cell.
    pub fn unexpected_mismatch(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        observed_decision: BaseCoverageCloseDecisionV1,
        observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
        evidence: BaseCoverageCloseResultEvidenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            manifest.root,
            cell.root,
            cell.source_case_id(),
            cell.group(),
            cell.facet(),
            cell.execution_scope(),
            cell.partition(),
            cell.expected_decision(),
            cell.expected_reason(),
            BaseCoverageCloseResultStatusV1::UnexpectedMismatch,
            Some(observed_decision),
            observed_reason,
            evidence,
        )
    }

    /// Explicit execution failure without a fabricated decision.
    pub fn execution_failure(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        evidence: BaseCoverageCloseResultEvidenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            manifest.root,
            cell.root,
            cell.source_case_id(),
            cell.group(),
            cell.facet(),
            cell.execution_scope(),
            cell.partition(),
            cell.expected_decision(),
            cell.expected_reason(),
            BaseCoverageCloseResultStatusV1::ExecutionFailure,
            None,
            None,
            evidence,
        )
    }

    /// Explicit unexplained skip without a fabricated decision.
    pub fn unexplained_skip(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        evidence: BaseCoverageCloseResultEvidenceV1,
    ) -> Result<Self, ConstructionErrorV2> {
        Self::new(
            manifest.root,
            cell.root,
            cell.source_case_id(),
            cell.group(),
            cell.facet(),
            cell.execution_scope(),
            cell.partition(),
            cell.expected_decision(),
            cell.expected_reason(),
            BaseCoverageCloseResultStatusV1::UnexplainedSkip,
            None,
            None,
            evidence,
        )
    }

    #[must_use]
    pub const fn close_manifest_root(&self) -> ContentHash {
        self.close_manifest_root
    }

    #[must_use]
    pub const fn cell_root(&self) -> ContentHash {
        self.cell_root
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
    pub const fn status(&self) -> BaseCoverageCloseResultStatusV1 {
        self.status
    }

    #[must_use]
    pub const fn observed_decision(&self) -> Option<BaseCoverageCloseDecisionV1> {
        self.observed_decision
    }

    #[must_use]
    pub const fn observed_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.observed_reason
    }

    /// Opaque caller-presented evidence identity.
    ///
    /// The close layer binds but does not verify, retain, or promote this root.
    #[must_use]
    pub const fn evidence(&self) -> &BaseCoverageCloseResultEvidenceV1 {
        &self.evidence
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Exact stable-ID/reason pair retained by Unsupported and inapplicable sets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseReasonMatchV1 {
    source_case_id: Box<str>,
    reason: BaseCoverageCloseReasonCodeV1,
}

impl BaseCoverageCloseReasonMatchV1 {
    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    #[must_use]
    pub const fn reason(&self) -> BaseCoverageCloseReasonCodeV1 {
        self.reason
    }
}

/// Checked full-set AC53 accounting report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseCoverageCloseReportV1 {
    close_manifest_root: ContentHash,
    results: Box<[BaseCoverageClosePresentedResultV1]>,
    eligible_positive_ids: Box<[Box<str>]>,
    matched_positive_ids: Box<[Box<str>]>,
    eligible_expected_refusal_ids: Box<[Box<str>]>,
    matched_expected_refusal_ids: Box<[Box<str>]>,
    eligible_expected_failure_ids: Box<[Box<str>]>,
    matched_expected_failure_ids: Box<[Box<str>]>,
    eligible_mutation_ids: Box<[Box<str>]>,
    matched_mutation_ids: Box<[Box<str>]>,
    expected_unsupported: Box<[BaseCoverageCloseReasonMatchV1]>,
    matched_unsupported: Box<[BaseCoverageCloseReasonMatchV1]>,
    expected_inapplicable: Box<[BaseCoverageCloseReasonMatchV1]>,
    matched_inapplicable: Box<[BaseCoverageCloseReasonMatchV1]>,
    unexpected_mismatch_ids: Box<[Box<str>]>,
    execution_failure_ids: Box<[Box<str>]>,
    unexplained_skip_ids: Box<[Box<str>]>,
    adversarial_eligible: u32,
    adversarial_matched: u32,
    first_divergence_id: Option<Box<str>>,
    first_divergence_root: Option<ContentHash>,
    root: ContentHash,
}

impl BaseCoverageCloseReportV1 {
    /// Exact-join one result for every full close-manifest cell.
    ///
    /// There is intentionally no selected-subset argument. Empty, partial,
    /// extra, duplicate, reordered, stale, reclassified, wrong-scope, or
    /// wrong-expected-reason rows refuse before a report is returned.
    pub fn reconstruct_full(
        manifest: &BaseCoverageCloseManifestV1,
        presented: &[BaseCoverageClosePresentedResultV1],
    ) -> Result<Self, ConstructionErrorV2> {
        validate_close_results(manifest, presented)?;

        let mut eligible_positive_ids = Vec::new();
        let mut matched_positive_ids = Vec::new();
        let mut eligible_expected_refusal_ids = Vec::new();
        let mut matched_expected_refusal_ids = Vec::new();
        let mut eligible_expected_failure_ids = Vec::new();
        let mut matched_expected_failure_ids = Vec::new();
        let mut eligible_mutation_ids = Vec::new();
        let mut matched_mutation_ids = Vec::new();
        let mut expected_unsupported = Vec::new();
        let mut matched_unsupported = Vec::new();
        let mut expected_inapplicable = Vec::new();
        let mut matched_inapplicable = Vec::new();
        let mut unexpected_mismatch_ids = Vec::new();
        let mut execution_failure_ids = Vec::new();
        let mut unexplained_skip_ids = Vec::new();
        let mut first_divergence_id = None;
        let mut first_divergence_root = None;

        for result in presented {
            let id = || result.source_case_id.clone();
            let matched = result.status == BaseCoverageCloseResultStatusV1::Matched;
            match result.partition {
                BaseCoverageClosePartitionV1::Positive => {
                    eligible_positive_ids.push(id());
                    if matched {
                        matched_positive_ids.push(id());
                    }
                }
                BaseCoverageClosePartitionV1::ExpectedRefusal => {
                    eligible_expected_refusal_ids.push(id());
                    if matched {
                        matched_expected_refusal_ids.push(id());
                    }
                }
                BaseCoverageClosePartitionV1::ExpectedFailure => {
                    eligible_expected_failure_ids.push(id());
                    if matched {
                        matched_expected_failure_ids.push(id());
                    }
                }
                BaseCoverageClosePartitionV1::Mutation => {
                    eligible_mutation_ids.push(id());
                    if matched {
                        matched_mutation_ids.push(id());
                    }
                }
                BaseCoverageClosePartitionV1::Unsupported => {
                    let reason = result.expected_reason.ok_or_else(|| {
                        refusal(
                            ConstructionErrorKindV2::Missing,
                            "coverage.close_report.unsupported_reason",
                            "one exact registered Unsupported reason",
                            result.source_case_id(),
                        )
                    })?;
                    expected_unsupported.push(BaseCoverageCloseReasonMatchV1 {
                        source_case_id: id(),
                        reason,
                    });
                    if matched {
                        matched_unsupported.push(BaseCoverageCloseReasonMatchV1 {
                            source_case_id: id(),
                            reason,
                        });
                    }
                }
                BaseCoverageClosePartitionV1::Inapplicable => {
                    let reason = result.expected_reason.ok_or_else(|| {
                        refusal(
                            ConstructionErrorKindV2::Missing,
                            "coverage.close_report.inapplicable_reason",
                            "one exact registered inapplicability reason",
                            result.source_case_id(),
                        )
                    })?;
                    expected_inapplicable.push(BaseCoverageCloseReasonMatchV1 {
                        source_case_id: id(),
                        reason,
                    });
                    if matched {
                        matched_inapplicable.push(BaseCoverageCloseReasonMatchV1 {
                            source_case_id: id(),
                            reason,
                        });
                    }
                }
            }
            match result.status {
                BaseCoverageCloseResultStatusV1::Matched => {}
                BaseCoverageCloseResultStatusV1::UnexpectedMismatch => {
                    unexpected_mismatch_ids.push(id());
                }
                BaseCoverageCloseResultStatusV1::ExecutionFailure => {
                    execution_failure_ids.push(id());
                }
                BaseCoverageCloseResultStatusV1::UnexplainedSkip => {
                    unexplained_skip_ids.push(id());
                }
            }
            if first_divergence_id.is_none()
                && result.status != BaseCoverageCloseResultStatusV1::Matched
            {
                first_divergence_id = Some(id());
                first_divergence_root = Some(result.root());
            }
        }

        let adversarial_eligible = checked_adversarial_total(
            eligible_expected_refusal_ids.len(),
            eligible_expected_failure_ids.len(),
            eligible_mutation_ids.len(),
            "coverage.close_report.adversarial_eligible",
        )?;
        let adversarial_matched = checked_adversarial_total(
            matched_expected_refusal_ids.len(),
            matched_expected_failure_ids.len(),
            matched_mutation_ids.len(),
            "coverage.close_report.adversarial_matched",
        )?;
        let root = close_report_root(
            manifest.root,
            presented,
            adversarial_eligible,
            adversarial_matched,
            first_divergence_id.as_deref(),
            first_divergence_root,
        )?;
        Ok(Self {
            close_manifest_root: manifest.root,
            results: presented.to_vec().into_boxed_slice(),
            eligible_positive_ids: eligible_positive_ids.into_boxed_slice(),
            matched_positive_ids: matched_positive_ids.into_boxed_slice(),
            eligible_expected_refusal_ids: eligible_expected_refusal_ids.into_boxed_slice(),
            matched_expected_refusal_ids: matched_expected_refusal_ids.into_boxed_slice(),
            eligible_expected_failure_ids: eligible_expected_failure_ids.into_boxed_slice(),
            matched_expected_failure_ids: matched_expected_failure_ids.into_boxed_slice(),
            eligible_mutation_ids: eligible_mutation_ids.into_boxed_slice(),
            matched_mutation_ids: matched_mutation_ids.into_boxed_slice(),
            expected_unsupported: expected_unsupported.into_boxed_slice(),
            matched_unsupported: matched_unsupported.into_boxed_slice(),
            expected_inapplicable: expected_inapplicable.into_boxed_slice(),
            matched_inapplicable: matched_inapplicable.into_boxed_slice(),
            unexpected_mismatch_ids: unexpected_mismatch_ids.into_boxed_slice(),
            execution_failure_ids: execution_failure_ids.into_boxed_slice(),
            unexplained_skip_ids: unexplained_skip_ids.into_boxed_slice(),
            adversarial_eligible,
            adversarial_matched,
            first_divergence_id,
            first_divergence_root,
            root,
        })
    }

    #[must_use]
    pub const fn close_manifest_root(&self) -> ContentHash {
        self.close_manifest_root
    }

    #[must_use]
    pub fn results(&self) -> &[BaseCoverageClosePresentedResultV1] {
        &self.results
    }

    #[must_use]
    pub fn eligible_positive_ids(&self) -> &[Box<str>] {
        &self.eligible_positive_ids
    }

    #[must_use]
    pub fn matched_positive_ids(&self) -> &[Box<str>] {
        &self.matched_positive_ids
    }

    #[must_use]
    pub fn eligible_expected_refusal_ids(&self) -> &[Box<str>] {
        &self.eligible_expected_refusal_ids
    }

    #[must_use]
    pub fn matched_expected_refusal_ids(&self) -> &[Box<str>] {
        &self.matched_expected_refusal_ids
    }

    #[must_use]
    pub fn eligible_expected_failure_ids(&self) -> &[Box<str>] {
        &self.eligible_expected_failure_ids
    }

    #[must_use]
    pub fn matched_expected_failure_ids(&self) -> &[Box<str>] {
        &self.matched_expected_failure_ids
    }

    #[must_use]
    pub fn eligible_mutation_ids(&self) -> &[Box<str>] {
        &self.eligible_mutation_ids
    }

    #[must_use]
    pub fn matched_mutation_ids(&self) -> &[Box<str>] {
        &self.matched_mutation_ids
    }

    #[must_use]
    pub fn expected_unsupported(&self) -> &[BaseCoverageCloseReasonMatchV1] {
        &self.expected_unsupported
    }

    #[must_use]
    pub fn matched_unsupported(&self) -> &[BaseCoverageCloseReasonMatchV1] {
        &self.matched_unsupported
    }

    #[must_use]
    pub fn expected_inapplicable(&self) -> &[BaseCoverageCloseReasonMatchV1] {
        &self.expected_inapplicable
    }

    #[must_use]
    pub fn matched_inapplicable(&self) -> &[BaseCoverageCloseReasonMatchV1] {
        &self.matched_inapplicable
    }

    #[must_use]
    pub fn unexpected_mismatch_ids(&self) -> &[Box<str>] {
        &self.unexpected_mismatch_ids
    }

    #[must_use]
    pub fn execution_failure_ids(&self) -> &[Box<str>] {
        &self.execution_failure_ids
    }

    #[must_use]
    pub fn unexplained_skip_ids(&self) -> &[Box<str>] {
        &self.unexplained_skip_ids
    }

    /// Checked sum of expected-refusal, expected-failure, and mutation cells.
    #[must_use]
    pub const fn adversarial_eligible(&self) -> u32 {
        self.adversarial_eligible
    }

    /// Checked sum of the three independently matched adversarial partitions.
    #[must_use]
    pub const fn adversarial_matched(&self) -> u32 {
        self.adversarial_matched
    }

    #[must_use]
    pub fn first_divergence_id(&self) -> Option<&str> {
        self.first_divergence_id.as_deref()
    }

    /// Root of the first non-matched presented result in canonical order.
    #[must_use]
    pub const fn first_divergence_root(&self) -> Option<ContentHash> {
        self.first_divergence_root
    }

    /// The sole full-set green equation.
    ///
    /// Full cardinality/order is structural. Green additionally requires a
    /// nonempty positive corpus, a nonempty checked adversarial corpus, exact
    /// equality in every independent matched partition and registered-reason
    /// set, and no mismatch, execution-failure, or unexplained-skip IDs.
    #[must_use]
    pub fn is_green(&self) -> bool {
        !self.eligible_positive_ids.is_empty()
            && self.adversarial_eligible != 0
            && self.matched_positive_ids == self.eligible_positive_ids
            && self.matched_expected_refusal_ids == self.eligible_expected_refusal_ids
            && self.matched_expected_failure_ids == self.eligible_expected_failure_ids
            && self.matched_mutation_ids == self.eligible_mutation_ids
            && !self.expected_unsupported.is_empty()
            && self.matched_unsupported == self.expected_unsupported
            && !self.expected_inapplicable.is_empty()
            && self.matched_inapplicable == self.expected_inapplicable
            && self.adversarial_matched == self.adversarial_eligible
            && self.unexpected_mismatch_ids.is_empty()
            && self.execution_failure_ids.is_empty()
            && self.unexplained_skip_ids.is_empty()
            && self.first_divergence_id.is_none()
            && self.first_divergence_root.is_none()
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

#[derive(Debug, Clone, Copy)]
struct RustTestClassTemplateV1 {
    class: BaseCoverageManifestClassV1,
    test_names: &'static [&'static str],
}

#[derive(Debug, Clone, Copy)]
struct RustTestModuleTemplateV1 {
    module: &'static str,
    source_path: &'static str,
    class_templates: &'static [RustTestClassTemplateV1],
}

const fn classified_tests(
    class: BaseCoverageManifestClassV1,
    test_names: &'static [&'static str],
) -> RustTestClassTemplateV1 {
    RustTestClassTemplateV1 { class, test_names }
}

// This is an independent, handwritten classification table. Each source test
// name occurs literally under exactly one evidence class; no substring, suffix,
// source-module, or runtime-result heuristic participates in classification.
const RUST_TEST_MODULE_TEMPLATES_V1: &[RustTestModuleTemplateV1] = &[
    RustTestModuleTemplateV1 {
        module: "budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["precise_refusal_carries_unit_expectation_observation_and_owner"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "intrinsic_nonzero_relations_and_timeout_arithmetic_refuse_precisely",
                    "profile_boundaries_are_exact_and_one_over_refuses",
                    "disposition_zero_rules_and_publication_equation_are_exact",
                    "logical_work_is_exact_u128_and_registered_units_require_exact_registry_membership",
                    "publication_sum_overflow_is_typed_even_without_concrete_storage",
                    "publication_accounting_accepts_zero_one_exact_and_maximum_boundaries",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &["independent_literal_oracle_freezes_all_18_fields_and_widths"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &["every_one_field_mutation_moves_the_canonical_projection_and_root"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &["every_output_grant_is_checked_against_the_admitted_limits"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "canonical",
        source_path: "crates/fs-evidence-runner/src/canonical.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "exact_frame_bound_accepts_and_one_over_refuses_atomically",
                    "count_and_frame_length_overflow_helpers_refuse_precisely",
                    "preflight_refusal_happens_before_the_encoding_pass",
                    "preflighted_frame_refuses_growth_beyond_the_exact_capacity",
                    "magic_is_part_of_the_bound",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["count_only_preflight_matches_the_exact_encoded_length"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "integer_and_presence_fields_have_independent_known_big_endian_bytes",
                    "byte_and_string_fields_use_exact_u32_length_prefixes",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &["preflighted_frame_refuses_a_divergent_second_pass"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["policy_root_binds_each_semantic_field_and_ignores_opaque_observations"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "registry_is_bounded_sorted_duplicate_free_and_nonzero",
                    "registration_and_command_policy_sets_are_exact",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &["least_privilege_matrix_rejects_every_one_right_mutant"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "catalog",
        source_path: "crates/fs-evidence-runner/src/catalog.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "api_wire_and_predecessor_literal_oracle",
                    "terminal_command_and_diagnostic_literal_oracles",
                    "value_digest_and_role_literal_oracles",
                    "publication_capability_and_extension_literal_oracles",
                    "decision_detail_registry_reconstructs_the_exact_base_and_family_inventory",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "decision_detail_registry_refuses_unknown_collision_reorder_and_count_cap",
                    "every_decision_detail_registry_field_moves_identity_and_exactness_refuses",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &[
                    "list_has_no_selectors_budgets_disposition_or_publication",
                    "run_selection_has_no_default_mode_and_records_caller_provenance",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "exact_command_table_accepts_all_and_only_frozen_cells",
                    "profile_disposition_and_publication_presence_cannot_drift",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &["every_cross_command_provenance_cell_refuses"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "construction",
        source_path: "crates/fs-evidence-runner/src/construction.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &["observed_rendering_is_utf8_bounded_without_recursive_diagnostics"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &["sensitive_observed_fuzz_corpus_never_echoes_through_any_rendering"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "coverage",
        source_path: "crates/fs-evidence-runner/src/coverage.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "registered_extension_capability_registry_boundaries_and_diagnostics_are_exact",
                    "registered_extension_capability_sets_boundaries_and_diagnostics_are_exact",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "registered_extension_capability_ids_and_nominal_descriptor_roles_are_exact",
                    "registered_extension_capability_canonical_roots_and_magics_are_exact",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "dependency",
        source_path: "crates/fs-evidence-runner/src/dependency.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "current_and_eventual_literal_oracles_are_exact",
                    "every_pinned_source_identity_component_is_exact",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "every_absent_extra_or_order_mutant_refuses",
                    "every_forbidden_route_or_owner_feature_source_mutant_refuses",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &[
                    "compiled_phase_one_manifest_has_only_the_owned_normal_row",
                    "current_package_is_a_root_member_and_imports_only_its_owned_dependency",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &[
                    "repairs_are_bounded_structured_and_non_executable",
                    "typed_replacements_require_the_same_inline_or_retained_shape",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "every_repair_kind_rank_and_display_boundary_is_constructible_or_refused_exactly",
                    "diagnostic_requires_contiguous_repairs_and_joint_feasibility",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "diagnostic_counts_ranks_namespaces_and_prerequisites_are_exact",
                    "registered_detail_projection_is_bounded_sealed_and_non_authoritative",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "duplicate_repair_rank_refuses_while_count_remains_within_limit",
                    "every_registered_detail_identity_field_moves_the_projection_root",
                    "every_diagnostic_field_mutation_moves_the_root",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &["complete_frame_and_every_enclosing_grant_are_jointly_feasible"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["registered_descriptors_enforce_names_ids_and_canonical_allowed_units"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "logical_extent_base_axis_unit_table_and_u128_extrema_are_exact",
                    "registry_category_caps_are_independent_at_zero_one_64_and_65",
                    "u16_max_is_unknown_unless_registered_in_the_exact_typed_namespace",
                    "unit_conversion_checks_dimensions_normalization_extrema_and_overflow",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["registry_is_typed_permutation_invariant_and_exact_set_reconstructible"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "artifact_codec_catalog_is_exact_closed_and_registry_independent",
                    "logical_extent_schema_and_root_bind_axis_value_unit_in_exact_order",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "registry_refuses_duplicate_collision_unknown_and_over_cap_data",
                    "every_descriptor_field_and_allowed_unit_mutation_moves_the_registry_root",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &["duration_and_registered_axis_conversions_are_exact_and_bounded"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["wrapper_parser_checks_nominal_metadata_before_text"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &["digest_width_domain_and_all_zero_presence_are_checked"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "lowercase_hex_form_is_exact_and_round_trips_every_byte",
                    "constructor_owner_handoff_reconstruction_is_set_exact_and_order_stable",
                    "constructor_owner_handoff_ordered_closeout_rejects_reordering_without_weakening_exact_set",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "descriptor_inventory_is_complete_unique_and_generation_conformant",
                    "every_nominal_wrapper_has_a_checked_presented_parser",
                    "constructor_owner_handoff_exactly_covers_every_nominal_descriptor",
                    "root_free_guard_inventory_is_exact_ordered_and_root_stable",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "wrappers_reject_role_and_domain_substitution",
                    "each_digest_input_field_moves_presented_identity",
                    "owner_schema_domain_and_role_mutations_move_the_handoff_root_and_refuse",
                    "presented_constructor_owner_handoff_refuses_every_stale_or_unknown_identity_field",
                    "root_free_guard_rejects_missing_extra_duplicate_and_reordered_rows",
                    "root_free_guard_rejects_owner_edge_and_context_collapsing_mutants",
                    "root_free_guard_rejects_every_fabricated_identity_and_widening_surface",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "limits",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "profile_dependent_ceilings_are_exact",
                    "every_field_has_exact_width_zero_minimum_ceiling_and_maximum_boundaries",
                    "executable_structural_minima_and_declared_nested_minima_are_enforced",
                    "nested_and_per_case_capacities_remain_jointly_feasible",
                    "lifecycle_equation_is_checked_for_zero_256_and_overflow_cases",
                    "protocol_stored_relations_and_abstract_envelope_bound_are_exact",
                    "exact_256_artifact_envelope_accepts_and_257_refuses_precisely",
                    "checked_addition_refuses_overflow_before_accounting",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["whole_publication_algebra_recomputes_every_total"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "independent_literal_oracle_covers_all_71_fields",
                    "canonical_projection_encodes_all_71_width_tags_and_payloads_exactly",
                    "registry_limit_tail_is_distinct_profile_equal_and_tightenable",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "semantic_root_is_nominal_and_moves_with_an_admitted_limit_change",
                    "every_one_over_mutation_refuses_and_fixed_fields_cannot_move",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "logging",
        source_path: "crates/fs-evidence-runner/src/logging.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["target_root_must_match_the_exact_target_token"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "collection_bounds_and_reproduction_shape_fail_closed",
                    "complete_close_validation_refuses_zero_and_one_over_bounds_before_joining",
                    "bounded_writer_zero_one_exact_and_one_over_budget_boundaries_are_terminal",
                    "bounded_writer_checked_overflow_and_detail_count_bounds_refuse_first",
                    "schema_impact_logging_checks_arithmetic_and_never_accepts_raw_hostile_values",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "feature_and_target_roots_are_deterministic_canonical_and_sensitive",
                    "normalized_prefix_suffix_and_embedded_sensitive_aliases_refuse",
                    "canonical_event_and_log_roots_are_order_independent_but_mutation_sensitive",
                    "complete_close_log_replay_root_movement_and_red_divergence_are_exact",
                    "bounded_writer_replays_complete_log_deterministically_and_moves_with_input",
                    "budget_exceeded_terminal_binds_every_field_and_overflow_document_replays",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "closed_field_catalog_is_total_unique_and_round_trips",
                    "logging_schema_root_is_exact_deterministic_and_mutation_sensitive",
                    "case_outcome_matrix_and_first_divergence_are_exact",
                    "script_mapping_and_retained_artifact_are_distinct_typed_concepts",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "arbitrary_field_names_and_duplicate_codes_refuse_before_admission",
                    "exact_event_matrices_refuse_missing_extra_and_wrong_typed_fields",
                    "duplicate_case_and_unexpected_mismatch_cannot_form_a_green_log",
                    "ac35_row_evidence_fields_are_typed_required_and_mutation_checked",
                    "detail_manifest_fields_are_required_and_green_reconcile_exactly",
                    "caller_controlled_logging_rejections_never_echo_through_any_rendering",
                    "complete_close_effects_applicability_and_artifacts_refuse_substitution",
                    "bounded_writer_refuses_missing_extra_duplicate_reordered_and_post_terminal_details",
                    "bounded_writer_preserves_prefix_and_never_reports_overflow_as_success",
                    "repair_manifest_and_budget_terminal_are_deterministic_redacted_and_no_echo",
                    "schema_impact_report_retains_the_first_typed_and_rooted_divergence",
                    "schema_impact_log_refuses_missing_extra_duplicate_reordered_and_count_gaps",
                    "schema_impact_log_case_manifest_refuses_duplicate_cases_and_mixed_fragment_sources",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &[
                    "full_log_reconciles_sequences_journeys_rows_results_cells_and_counts",
                    "positive_and_expected_refusal_partitions_are_distinct_and_exact",
                    "mixed_semantic_rows_reconcile_exact_terminal_partitions",
                    "source_closure_green_counts_are_matches_not_expected_refusals",
                    "schema_impact_log_reconciles_every_partition_reason_and_source_fragment",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["windows_ascii_aliases_refuse_and_non_ascii_is_explicitly_unsupported"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "byte_and_segment_boundaries_are_exact",
                    "unsafe_or_ambiguous_single_path_forms_refuse",
                    "content_store_rejects_both_exact_reserved_first_segment_prefixes",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "exact_utf8_bytes_are_preserved_without_normalization",
                    "canonical_order_is_segment_sequence_byte_order",
                    "duplicate_and_strict_segment_prefix_are_distinct_and_deterministic",
                    "set_adjudication_is_invariant_to_input_permutation",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &["posix_and_content_store_cells_are_exact_bytewise"],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "projection",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["exact_source_reconstruction_is_deterministic"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &["schema_impact_v2_projection_has_exact_rows_paths_counts_and_canonical_roots"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "schema_impact_v2_projection_rejects_missing_extra_reordered_path_count_and_root_mutants",
                    "phase_one_contribution_roots_move_and_budget_owner_cycle_boundaries_fail_closed",
                    "one_field_projection_mutation_moves_the_root",
                    "every_harness_context_field_moves_exactly_one_context_root",
                    "source_reconstruction_rejects_missing_entry",
                    "source_reconstruction_rejects_extra_entry",
                    "source_reconstruction_rejects_duplicate_entry",
                    "source_reconstruction_rejects_reordered_entries",
                    "source_reconstruction_rejects_mutated_bytes",
                    "source_reconstruction_rejects_owner_route_identity_and_policy_mutations",
                    "source_reconstruction_rejects_length_content_and_snapshot_mutations",
                    "source_reconstruction_rejects_resealed_bytes_and_dependency_identity",
                    "source_entry_root_moves_for_each_admitted_metadata_field",
                    "journey_routes_reject_wrong_owner_driver_script_or_manifest",
                    "result_join_rejects_missing_first_middle_and_last_rows",
                    "result_join_rejects_extra_duplicate_and_reordered_rows",
                    "result_join_rejects_stale_unmapped_and_cross_journey_rows",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &[
                    "manifest_exactly_maps_five_scripts_and_all_base_rows",
                    "all_real_constructor_rows_agree_and_logs_are_deterministic",
                    "source_closure_membership_and_order_are_exact_and_content_bound",
                    "ac38_coverage_manifest_source_of_truth_and_checked_report_are_exact",
                    "exact_result_join_reconstructs_the_frozen_journey",
                    "phase_one_contribution_is_exact_result_free_deferred_and_no_execution",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "publication",
        source_path: "crates/fs-evidence-runner/src/publication.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "publication_storage_zero_max_role_envelope_and_overflow_boundaries_are_exact",
                    "atomic_result_projection_rejects_wrong_presence_and_second_copy_pressure",
                    "result_and_failure_caps_are_inclusive",
                    "every_command_result_shape_and_nested_boundary_is_exact",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["physical_observations_have_no_semantic_selection_field"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "profile_protocol_and_target_are_one_exact_cell",
                    "compatibility_matrix_accepts_exactly_three_cells_per_destination_mode",
                    "whole_publication_projection_has_one_domain_separated_canonical_root",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "each_semantic_publication_mutation_moves_the_projection_root",
                    "every_publication_storage_field_moves_identity_and_bad_accounting_refuses",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "schema_impact_accepts_every_exact_maximum_and_refuses_edge_overflow",
                    "schema_impact_text_ceilings_are_exact_and_one_over_refuses",
                    "schema_impact_row_collection_maxima_and_preflight_precedence_are_exact",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &["schema_impact_manifest_enforces_ownership_reciprocal_slots_and_acyclic_dag"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "component_descriptor_bytes_and_roots_match_independent_literal_oracles",
                    "registry_row_and_manifest_roots_match_independent_literal_oracles",
                    "production_meta_manifest_matches_independent_literal_oracle",
                    "nominal_root_registry_and_compatible_snapshot_are_exact",
                    "schema_impact_closed_catalogs_wrappers_frames_and_rows_are_exact",
                    "schema_impact_disposition_authority_policy_and_surface_matrix_is_exhaustive",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "schema_impact_rejects_every_shape_slot_and_snapshot_mutant",
                    "schema_impact_component_row_and_manifest_root_movement_is_complete",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &[
                    "schema_impact_manifest_reconstructs_the_source_frozen_meta_schema_without_authority",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &["active_stop_states_require_the_exact_nominal_basis"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &["not_run_basis_validates_manifest_boundaries_and_exact_diagnostic"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "not_run_remaining_suffix_arithmetic_is_exhaustive_and_allocation_free",
                    "exhaustive_cartesian_matrix_and_error_precedence",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "not_run_causes_have_exact_codes_names_and_nominal_accessors",
                    "exact_diagnostic_mapping_literal_oracle",
                    "role_state_sets_and_lifecycle_basis_are_exact",
                ],
            ),
        ],
    },
    RustTestModuleTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        class_templates: &[
            classified_tests(
                BaseCoverageManifestClassV1::Unit,
                &[
                    "semantic_seed_errors_have_exact_actionable_non_authority_metadata",
                    "semantic_seed_policy_payloads_bind_case_manifest_registry_domain_and_versions",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "rational_reduces_sign_zero_and_i128_min_without_overflow",
                    "decimal_has_one_representation_and_refuses_range_crossing",
                    "every_integer_width_preserves_both_extrema_exactly",
                    "token_boundaries_and_segment_grammar_are_exact",
                    "text_and_opaque_bytes_enforce_exact_byte_caps",
                    "semantic_seed_material_enforces_31_32_33_byte_boundaries_and_exact_bytes",
                    "semantic_seed_cli_is_exact_nonambient_and_rejects_every_noncanonical_form",
                    "stable_case_identity_accepts_real_coverage_ids_and_rejects_noncanonical_mutants",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "semantic_seed_derivation_registry_exact_reconstruction_rejects_missing_extra_mutation_and_stale_root",
                    "invocation_derived_material_and_roots_move_with_case_domain_registry_and_versions",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "ieee_wrappers_preserve_special_encodings_and_nan_payloads",
                    "numeric_tags_are_exact_and_nonrecursive",
                    "units_require_positive_canonical_scale_and_keep_exponent_order",
                    "typed_value_and_presence_tags_are_exact",
                    "semantic_seed_descriptors_freeze_api_wire_domains_and_no_claims",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::Mutation,
                &[
                    "semantic_seed_derivation_registry_rejects_zero_duplicate_reorder_collision_and_one_over",
                    "fixed_manifest_provenance_moves_with_case_manifest_material_and_versions",
                    "semantic_seed_resolution_rejects_unknown_and_cross_case_provenance",
                    "semantic_seed_canonical_provenance_and_debug_are_material_redacted",
                ],
            ),
        ],
    },
];

#[derive(Debug, Clone, Copy)]
struct CompileFailTemplateV1 {
    module: &'static str,
    source_path: &'static str,
    case_name: &'static str,
    expected_error_code: &'static str,
}

const COMPILE_FAIL_TEMPLATES_V1: &[CompileFailTemplateV1] = &[
    CompileFailTemplateV1 {
        module: "budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
        case_name: "no-postmutation-of-runner-budgets",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
        case_name: "no-postmutation-of-admitted-runner-budgets",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
        case_name: "no-physical-acquisition-material-as-semantic-right",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
        case_name: "immutable-policy-root",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
        case_name: "immutable-command-intent",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
        case_name: "no-command-intent-as-authority-scope",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-terminal-extension",
        expected_error_code: "E0609",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-refusal-extension",
        expected_error_code: "E0609",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-authority-mint",
        expected_error_code: "E0599",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-executable-repair",
        expected_error_code: "E0609",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "extension_capability_id_rejects_base_capability_id",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "base_capability_id_rejects_extension_capability_id",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "artifact-codec-has-no-executable-encoder",
        expected_error_code: "E0599",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "private-registered-artifact-role-fields",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "typed-registered-namespaces-cannot-cross-substitute",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "private-logical-extent-fields",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "sealed-digest-domain",
        expected_error_code: "E0423",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "no-generic-digest-as-nominal-root",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-source-identity-constructor",
        expected_error_code: "E0624",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-lifecycle-log-constructor",
        expected_error_code: "E0624",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-durable-publication-constructor",
        expected_error_code: "E0624",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-authority-scope-constructor",
        expected_error_code: "E0624",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "no-standalone-root-for-root-free-evaluator-members",
        expected_error_code: "E0432",
    },
    CompileFailTemplateV1 {
        module: "limits",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
        case_name: "no-postmutation-of-runner-limits",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-raw-string-as-logical-path",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "private-logical-path-constructor",
        expected_error_code: "E0423",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-postmutation-of-logical-path",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-bundle-path-as-content-key",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-postmutation-of-content-store-key",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "projection",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
        case_name: "immutable-source-closure",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "canonical-frame-version-is-not-runner-wire-version",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "compatible-source-snapshot-fields-are-private",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "leaf-extension-source-admission-is-private",
        expected_error_code: "E0624",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "leaf-extension-cannot-fill-frozen-base-position",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "raw-role-parts-cannot-mint-registry-membership",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "no-postvalidation-schema-impact-row-mutation",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "no-caller-forged-schema-impact-manifest-ordinal",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "schema-impact-row-root-is-not-manifest-root",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "schema_impact",
        source_path: "crates/fs-evidence-runner/src/schema_impact.rs",
        case_name: "generic-content-hash-cannot-mint-schema-impact-row-root",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "nominal-cancelled-root",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-generic-digest-cause",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "cause-payload-required",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-profile-filter-cause",
        expected_error_code: "E0599",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "private-not-run-basis-fields",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-unvalidated-state-candidate-as-terminal",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-rational-fields",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-rational",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-decimal-fields",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-decimal",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f32-bits-has-no-implicit-cmp",
        expected_error_code: "E0599",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f32-bits-has-no-relational-order",
        expected_error_code: "E0369",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f32-bits-has-no-btree-order",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f32-bits-has-no-default-sort",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f32-bits-has-no-key-sort",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f64-bits-has-no-implicit-cmp",
        expected_error_code: "E0599",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f64-bits-has-no-relational-order",
        expected_error_code: "E0369",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f64-bits-has-no-btree-order",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f64-bits-has-no-default-sort",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "f64-bits-has-no-key-sort",
        expected_error_code: "E0277",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-unit-fields",
        expected_error_code: "E0451",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-unit",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-string-as-stable-token",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-stable-token",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-string-as-bounded-text",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-bounded-text",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-vec-as-opaque-bytes",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-opaque-bytes",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-digest-as-typed-absence",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "raw_seed_material_cannot_bypass_validation",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "seed_material_fields_are_private",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "nominal_seed_roots_are_not_cross_substitutable",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "raw_invocation_seed_text_cannot_bypass_parser",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "stable_case_identity_fields_are_private",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "registry_root_rejects_cross_role_substitution",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "fixed_manifest_requires_nominal_case_manifest_root",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "invocation_derived_binding_fields_are_private",
        expected_error_code: "E0616",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "unregistered_domain_row_cannot_enter_policy",
        expected_error_code: "E0308",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "raw_seed_material_cannot_enter_policy",
        expected_error_code: "E0308",
    },
];

fn is_exact_compiler_error_code_v1(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() == 5 && bytes[0] == b'E' && bytes[1..].iter().all(u8::is_ascii_digit)
}

const MANIFEST_CONTRACT_TEST_NAMES_V1: &[&str] = &[
    "frozen_manifest_enumerates_exact_ratified_inventory_and_external_classes",
    "exact_base_reconstruction_rejects_missing_extra_duplicate_and_reordering",
    "declaration_grammar_and_semantic_mutations_refuse_or_move_identity",
    "extension_constructor_requires_external_class_order_and_global_uniqueness",
    "selection_accepts_exact_subsets_and_rejects_unknown_duplicate_and_reordered_ids",
    "checked_join_accepts_empty_and_mixed_exact_results",
    "checked_join_rejects_missing_and_selected_extra_results",
    "checked_join_rejects_stale_unmapped_and_multiply_reported_ids",
    "checked_join_rejects_reordered_results",
    "result_selection_manifest_and_report_roots_bind_every_semantic_field",
    "full_set_close_manifest_exactly_covers_nine_groups_and_twenty_two_facets",
    "full_set_close_manifest_exact_reconstruction_rejects_all_sequence_and_semantic_mutants",
    "full_set_close_report_reconstructs_exact_green_partitions_and_adversarial_sum",
    "full_set_close_report_refuses_wrong_manifest_cell_reason_scope_and_order",
    "full_set_close_report_is_fail_closed_for_mismatch_execution_failure_skip_and_partial_rows",
    "full_set_close_external_rows_are_inapplicable_immutable_contributions_never_passes",
    "full_set_close_first_divergence_and_roots_bind_every_semantic_partition",
    "full_set_close_empty_or_base_only_manifest_cannot_mint_full_set_authority",
    "five_explicits_numeric_domain_and_unit_references_are_exact",
    "five_explicits_numeric_surface_bounds_order_and_registered_identity_are_exact",
    "five_explicits_budget_axis_catalog_widths_units_and_order_are_exact",
    "five_explicits_budget_bounds_soft_relations_and_shape_refusals_are_exact",
    "five_explicits_downstream_budgets_keep_soft_rows_and_process_shape_independent",
    "five_explicits_profiles_and_one_field_mutations_move_roots",
    "five_explicits_maximum_numeric_frame_is_feasible",
    "race_facet_is_registered_inapplicable_for_pure_single_threaded_validator",
    "trait_facet_is_registered_inapplicable_without_public_trait_contract",
    "cancellation_facet_is_registered_inapplicable_for_pure_bounded_validator",
    "release_built_no_mock_e2e_facet_is_registered_as_downstream_owned",
];

const CLOSE_PUBLICATION_STATE_ROW_IDS_V1: &[&str] = &[
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

const CLOSE_PUBLICATION_V2_ROW_IDS_V1: &[&str] = &[
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

const CLOSE_VERIFIER_V2_ROW_IDS_V1: &[&str] = &[
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

const CLOSE_CANONICAL_RUNNER_V2_ROW_IDS_V1: &[&str] = &[
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

const CLOSE_RJOQ_HANDOFF_V1_ROW_IDS_V1: &[&str] = &[
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

#[derive(Debug, Clone, Copy)]
struct CloseJourneyRowsV1 {
    journey_key: &'static str,
    row_ids: &'static [&'static str],
}

const CLOSE_JOURNEY_ROWS_V1: &[CloseJourneyRowsV1] = &[
    CloseJourneyRowsV1 {
        journey_key: "canonical-runner-v2",
        row_ids: CLOSE_CANONICAL_RUNNER_V2_ROW_IDS_V1,
    },
    CloseJourneyRowsV1 {
        journey_key: "publication-state-v2",
        row_ids: CLOSE_PUBLICATION_STATE_ROW_IDS_V1,
    },
    CloseJourneyRowsV1 {
        journey_key: "publication-v2",
        row_ids: CLOSE_PUBLICATION_V2_ROW_IDS_V1,
    },
    CloseJourneyRowsV1 {
        journey_key: "rjoq-handoff-v1",
        row_ids: CLOSE_RJOQ_HANDOFF_V1_ROW_IDS_V1,
    },
    CloseJourneyRowsV1 {
        journey_key: "verifier-v2",
        row_ids: CLOSE_VERIFIER_V2_ROW_IDS_V1,
    },
];

const CLOSE_SOURCE_CLOSURE_IDS_V1: &[&str] = &[
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

const CLOSE_EXTENSION_CLASS_COUNTS_V1: [(BaseCoverageManifestClassV1, usize); 6] = [
    (BaseCoverageManifestClassV1::ProjectionE2e, 98),
    (BaseCoverageManifestClassV1::RuntimeLogging, 1),
    (BaseCoverageManifestClassV1::SourceClosure, 15),
    (BaseCoverageManifestClassV1::ExternalE2eScript, 5),
    (BaseCoverageManifestClassV1::ExternalMutation, 1),
    (BaseCoverageManifestClassV1::ExternalGovernance, 1),
];

fn frozen_full_source_manifest_v1() -> Result<BaseCoverageManifestV1, ConstructionErrorV2> {
    let mut extensions = Vec::with_capacity(121);
    for journey in CLOSE_JOURNEY_ROWS_V1 {
        for row_id in journey.row_ids {
            extensions.push(BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                format!("projection-e2e:{}:{row_id}", journey.journey_key),
                "crates/fs-evidence-runner/src/projection.rs",
            )?);
        }
    }
    extensions.push(BaseCoverageCaseDeclarationV1::new(
        BaseCoverageManifestClassV1::RuntimeLogging,
        "runtime-logging:aggregate-closed-log",
        "crates/fs-evidence-runner/src/logging.rs",
    )?);
    for id in CLOSE_SOURCE_CLOSURE_IDS_V1 {
        extensions.push(BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::SourceClosure,
            *id,
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
        (left.class().code(), left.id(), left.source_path()).cmp(&(
            right.class().code(),
            right.id(),
            right.source_path(),
        ))
    });
    if extensions.len() != 121 {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.full_extension_count",
            "exactly 121 source-authoritative extension cells",
            extensions.len(),
        ));
    }
    for (class, expected_count) in CLOSE_EXTENSION_CLASS_COUNTS_V1 {
        let observed_count = extensions
            .iter()
            .filter(|extension| extension.class() == class)
            .count();
        if observed_count != expected_count {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.full_extension_class_count",
                "the independently frozen count for this extension class",
                format_args!("{}:{observed_count}", class.code()),
            ));
        }
    }
    BaseCoverageManifestV1::with_exact_extensions(&extensions)
}

fn validate_close_partition_shape(
    partition: BaseCoverageClosePartitionV1,
    decision: BaseCoverageCloseDecisionV1,
    reason: Option<BaseCoverageCloseReasonCodeV1>,
) -> Result<(), ConstructionErrorV2> {
    let valid_decision = match partition {
        BaseCoverageClosePartitionV1::Positive => decision == BaseCoverageCloseDecisionV1::Accept,
        BaseCoverageClosePartitionV1::ExpectedRefusal => {
            decision == BaseCoverageCloseDecisionV1::Refuse
        }
        BaseCoverageClosePartitionV1::ExpectedFailure => {
            decision == BaseCoverageCloseDecisionV1::Fail
        }
        BaseCoverageClosePartitionV1::Mutation => matches!(
            decision,
            BaseCoverageCloseDecisionV1::Accept | BaseCoverageCloseDecisionV1::Refuse
        ),
        BaseCoverageClosePartitionV1::Unsupported => {
            decision == BaseCoverageCloseDecisionV1::Unsupported
        }
        BaseCoverageClosePartitionV1::Inapplicable => {
            decision == BaseCoverageCloseDecisionV1::Inapplicable
        }
    };
    if !valid_decision {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.expected_decision",
            "the exact decision admitted by the accounting partition",
            decision.stable_name(),
        ));
    }
    let valid_reason = match partition {
        BaseCoverageClosePartitionV1::Unsupported => {
            reason == Some(BaseCoverageCloseReasonCodeV1::WindowsNonasciiAliasLocallyUnadjudicable)
        }
        BaseCoverageClosePartitionV1::Inapplicable => reason.is_some_and(|reason| {
            reason != BaseCoverageCloseReasonCodeV1::WindowsNonasciiAliasLocallyUnadjudicable
        }),
        _ => reason.is_none(),
    };
    if !valid_reason {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.expected_reason",
            "an exact registered reason only for Unsupported or inapplicable cells",
            reason.map_or(0, BaseCoverageCloseReasonCodeV1::code),
        ));
    }
    Ok(())
}

fn validate_close_reason_scope(
    reason: Option<BaseCoverageCloseReasonCodeV1>,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
) -> Result<(), ConstructionErrorV2> {
    let Some(reason) = reason else {
        return Ok(());
    };
    let descriptor_scope = reason.descriptor().execution_scope();
    let scope_matches = execution_scope == descriptor_scope
        || (reason == BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned
            && execution_scope == BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration);
    if !scope_matches {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.reason_execution_scope",
            "the registered reason scope or its exact local facet-declaration scope",
            execution_scope.stable_name(),
        ));
    }
    Ok(())
}

fn validate_close_evidence_scope_shape(
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    evidence: &BaseCoverageCloseResultEvidenceV1,
) -> Result<(), ConstructionErrorV2> {
    let expected = match execution_scope {
        BaseCoverageCloseExecutionScopeV1::CrateTest
        | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => {
            BaseCoverageCloseEvidenceKindV1::OwnedHarnessExecution
        }
        BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
            BaseCoverageCloseEvidenceKindV1::InProcessProjectionExecution
        }
        BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
            BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
        }
        BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
            BaseCoverageCloseEvidenceKindV1::ApplicabilityDeclaration
        }
    };
    if evidence.kind() != expected {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close_result.evidence_kind",
            "the exact evidence kind admitted by the execution scope",
            evidence.kind().code(),
        ));
    }
    Ok(())
}

fn validate_close_evidence_for_cell(
    manifest: &BaseCoverageCloseManifestV1,
    cell: &BaseCoverageCloseManifestCellV1,
    evidence: &BaseCoverageCloseResultEvidenceV1,
) -> Result<(), ConstructionErrorV2> {
    validate_close_evidence_scope_shape(cell.execution_scope(), evidence)?;
    match cell.execution_scope() {
        BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
            let contribution = cell.downstream_contribution().ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::Missing,
                    "coverage.close_result.downstream_contribution",
                    "the immutable contribution bound by this cell",
                    cell.source_case_id(),
                )
            })?;
            if evidence.root() != contribution.root() {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close_result.evidence_root",
                    "the exact immutable downstream contribution root",
                    cell.source_case_id(),
                ));
            }
        }
        BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
            let reason = cell.expected_reason().ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::Missing,
                    "coverage.close_result.applicability_reason",
                    "the exact registered facet-applicability reason",
                    cell.source_case_id(),
                )
            })?;
            let exact = close_applicability_evidence_root(manifest.reason_registry_root(), reason)?;
            if evidence.root() != exact {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close_result.evidence_root",
                    "the applicability root binding the reason registry and exact reason",
                    cell.source_case_id(),
                ));
            }
        }
        BaseCoverageCloseExecutionScopeV1::CrateTest
        | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest
        | BaseCoverageCloseExecutionScopeV1::InProcessProjection => {}
    }
    Ok(())
}

fn close_declaration(
    case: &BaseCoverageManifestCaseV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
) -> Result<BaseCoverageCloseCellDeclarationV1, ConstructionErrorV2> {
    BaseCoverageCloseCellDeclarationV1::new(
        case.ordinal(),
        case.id(),
        case.class(),
        case.source_path(),
        facet.group(),
        facet,
        execution_scope,
        partition,
        expected_decision,
        expected_reason,
    )
}

fn close_contribution_reference_root(
    domain: &'static str,
    kind: &'static str,
    case: &BaseCoverageManifestCaseV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECONTRIBREF\x01", 1_024)?;
    frame.push_str("coverage.close.contribution.reference_kind", kind)?;
    frame.push_str("coverage.close.contribution.source_case_id", case.id())?;
    frame.push_str(
        "coverage.close.contribution.source_path",
        case.source_path(),
    )?;
    Ok(frame.root(domain))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CloseDownstreamRouteV1 {
    owner: &'static str,
    driver: &'static str,
    script: &'static str,
    manifest_path: &'static str,
}

fn frozen_close_downstream_route(
    source_case_id: &str,
) -> Result<CloseDownstreamRouteV1, ConstructionErrorV2> {
    let route = match source_case_id {
        "external-e2e:publication-state-v2" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.2.2.2",
            driver: "e2e-evidence-runner-publication-state-v2-driver",
            script: "scripts/ci/e2e_evidence_runner_publication_state_v2.sh",
            manifest_path: "scripts/ci/manifests/evidence_runner_publication_state_v2_cases.v1.json",
        },
        "external-e2e:publication-v2" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.2.2.3.3",
            driver: "e2e-evidence-runner-publication-v2-driver",
            script: "scripts/ci/e2e_evidence_runner_publication_v2.sh",
            manifest_path: "scripts/ci/manifests/evidence_runner_publication_v2_cases.v1.json",
        },
        "external-e2e:verifier-v2" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.3.3.3.3",
            driver: "e2e-evidence-verifier-v2-driver",
            script: "scripts/ci/e2e_evidence_verifier_v2.sh",
            manifest_path: "scripts/ci/manifests/evidence_verifier_v2_cases.v1.json",
        },
        "external-e2e:canonical-runner-v2"
        | "external-mutation:base-contract-exact-result-join" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.4.1.4",
            driver: "canonical-evidence-runner-v2-e2e-driver",
            script: "scripts/ci/canonical_evidence_runner_v2.sh",
            manifest_path: "scripts/ci/manifests/canonical_evidence_runner_v2_cases.v1.json",
        },
        "external-e2e:rjoq-handoff-v1" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.5.3.1",
            driver: "verify-runner-rjoq-handoff-v1-driver",
            script: "scripts/ci/verify_runner_rjoq_handoff_v1.sh",
            manifest_path: "scripts/ci/manifests/runner_rjoq_handoff_verifier_v1_cases.v1.json",
        },
        "external-governance:live-source-dependency-closure" => CloseDownstreamRouteV1 {
            owner: "frankensim-epic-foundations-huq.24.1.3.1",
            driver: "runner-v2-tool-governance-e2e-driver",
            script: "scripts/ci/e2e_runner_v2_tool_governance.sh",
            manifest_path: "scripts/ci/manifests/runner_v2_tool_governance_cases.v1.json",
        },
        _ => {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "coverage.close.downstream_route",
                "one exact owner, driver, script, and immutable manifest route",
                source_case_id,
            ));
        }
    };
    Ok(route)
}

fn close_downstream_manifest_reference_root(
    case: &BaseCoverageManifestCaseV1,
    route: CloseDownstreamRouteV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEMANIFESTREF\x01", 2 * 1024)?;
    frame.push_str(
        "coverage.close.contribution.manifest.source_case_id",
        case.id(),
    )?;
    frame.push_str("coverage.close.contribution.manifest.owner", route.owner)?;
    frame.push_str("coverage.close.contribution.manifest.driver", route.driver)?;
    frame.push_str("coverage.close.contribution.manifest.script", route.script)?;
    frame.push_str(
        "coverage.close.contribution.manifest.path",
        route.manifest_path,
    )?;
    Ok(frame.root("org.frankensim.fs-evidence-runner.close-contribution-manifest-reference.v1"))
}

fn frozen_close_downstream_contribution(
    case: &BaseCoverageManifestCaseV1,
) -> Result<BaseCoverageCloseDownstreamContributionV1, ConstructionErrorV2> {
    const LITERAL_DOMAIN: &str =
        "org.frankensim.fs-evidence-runner.close-contribution-literal-oracle.v1";
    const SEMANTIC_DOMAIN: &str =
        "org.frankensim.fs-evidence-runner.close-contribution-semantic-input.v1";
    const SCHEMA_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-contribution-schema.v1";
    const LOG_SCHEMA_DOMAIN: &str =
        "org.frankensim.fs-evidence-runner.close-contribution-log-schema.v1";
    const SOURCE_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-contribution-source.v1";
    const BUILD_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-contribution-build.v1";
    let route = frozen_close_downstream_route(case.id())?;
    let budgets =
        BaseCoverageCloseContributionBudgetsV1::new(frozen_downstream_close_budget_set()?, 4, 2)?;
    let source_reference = close_contribution_reference_root(SOURCE_DOMAIN, "source", case)?;
    let source_root = SourceIdentityRootV2::parse_presented(
        SourceIdentityRootV2::DESCRIPTOR.role(),
        SourceIdentityRootV2::DESCRIPTOR.domain(),
        &source_reference.to_hex(),
    )
    .map_err(|_| {
        refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.contribution.source_root",
            "an exact presented SourceIdentityRootV2",
            case.id(),
        )
    })?;
    let build_reference = close_contribution_reference_root(BUILD_DOMAIN, "build", case)?;
    let build_root = BuildIdentityRootV2::parse_presented(
        BuildIdentityRootV2::DESCRIPTOR.role(),
        BuildIdentityRootV2::DESCRIPTOR.domain(),
        &build_reference.to_hex(),
    )
    .map_err(|_| {
        refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.contribution.build_root",
            "an exact presented BuildIdentityRootV2",
            case.id(),
        )
    })?;
    BaseCoverageCloseDownstreamContributionV1::new(
        close_contribution_reference_root(LITERAL_DOMAIN, "literal-expectation-oracle", case)?,
        close_contribution_reference_root(SEMANTIC_DOMAIN, "semantic-input", case)?,
        budgets,
        close_contribution_reference_root(SCHEMA_DOMAIN, "schema", case)?,
        close_contribution_reference_root(LOG_SCHEMA_DOMAIN, "log-schema", case)?,
        source_root,
        build_root,
        route.owner,
        StableTokenV2::new(route.driver).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.contribution.downstream_driver",
                "one exact stable downstream driver ID",
                route.driver,
            )
        })?,
        route.script,
        route.manifest_path,
        close_downstream_manifest_reference_root(case, route)?,
        "downstream-contribution-is-not-execution-proof",
    )
}

fn close_downstream_declaration(
    case: &BaseCoverageManifestCaseV1,
    facet: BaseCoverageCloseFacetV1,
    partition: BaseCoverageClosePartitionV1,
    expected_reason: BaseCoverageCloseReasonCodeV1,
) -> Result<BaseCoverageCloseCellDeclarationV1, ConstructionErrorV2> {
    let contribution = frozen_close_downstream_contribution(case)?;
    BaseCoverageCloseCellDeclarationV1::new_with_downstream_contribution(
        case.ordinal(),
        case.id(),
        case.class(),
        case.source_path(),
        facet.group(),
        facet,
        BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution,
        partition,
        BaseCoverageCloseDecisionV1::Inapplicable,
        Some(expected_reason),
        Some(contribution),
    )
}

fn projection_row_id_for_exact_source_case(source_case_id: &str) -> Option<&'static str> {
    for journey in CLOSE_JOURNEY_ROWS_V1 {
        for row_id in journey.row_ids {
            let exact = format!("projection-e2e:{}:{row_id}", journey.journey_key);
            if source_case_id == exact {
                return Some(*row_id);
            }
        }
    }
    None
}

#[allow(
    clippy::too_many_lines,
    reason = "the controlling full-set close manifest uses explicit stable-ID overrides rather than execution-result inference"
)]
fn classify_close_case_v1(
    case: &BaseCoverageManifestCaseV1,
) -> Result<BaseCoverageCloseCellDeclarationV1, ConstructionErrorV2> {
    use BaseCoverageCloseDecisionV1::{Accept, Fail, Inapplicable, Refuse, Unsupported};
    use BaseCoverageCloseExecutionScopeV1::{
        CompileFailDoctest, CrateTest, FacetApplicabilityDeclaration, InProcessProjection,
    };
    use BaseCoverageCloseFacetV1::{
        Api, Boundary, Cancellation, CompileFail, DetailedDeterministicLogging, Fault, Fuzz,
        ImmutableE2eContribution, Literal, Metamorphic, Model, Mutation, Property, Race, Redaction,
        ReleaseBuiltNoMockE2e, Resource, SourceClosure, State, Trait, Typestate, Unit,
    };
    use BaseCoverageClosePartitionV1::{
        ExpectedFailure, ExpectedRefusal, Inapplicable as InapplicablePartition,
        Mutation as MutationPartition, Positive, Unsupported as UnsupportedPartition,
    };
    use BaseCoverageCloseReasonCodeV1::{
        CancellationNotApplicablePureBoundedValidator,
        RaceNotApplicablePureSingleThreadedValidator, ReleaseExecutionDownstreamOwned,
        TraitNotApplicableNoPublicTraitContract, WindowsNonasciiAliasLocallyUnadjudicable,
    };

    let id = case.id();
    match id {
        "manifest-contract:coverage:race_facet_is_registered_inapplicable_for_pure_single_threaded_validator" =>
        {
            return close_declaration(
                case,
                Race,
                FacetApplicabilityDeclaration,
                InapplicablePartition,
                Inapplicable,
                Some(RaceNotApplicablePureSingleThreadedValidator),
            );
        }
        "manifest-contract:coverage:trait_facet_is_registered_inapplicable_without_public_trait_contract" =>
        {
            return close_declaration(
                case,
                Trait,
                FacetApplicabilityDeclaration,
                InapplicablePartition,
                Inapplicable,
                Some(TraitNotApplicableNoPublicTraitContract),
            );
        }
        "manifest-contract:coverage:cancellation_facet_is_registered_inapplicable_for_pure_bounded_validator" =>
        {
            return close_declaration(
                case,
                Cancellation,
                FacetApplicabilityDeclaration,
                InapplicablePartition,
                Inapplicable,
                Some(CancellationNotApplicablePureBoundedValidator),
            );
        }
        "manifest-contract:coverage:release_built_no_mock_e2e_facet_is_registered_as_downstream_owned" =>
        {
            return close_declaration(
                case,
                ReleaseBuiltNoMockE2e,
                FacetApplicabilityDeclaration,
                InapplicablePartition,
                Inapplicable,
                Some(ReleaseExecutionDownstreamOwned),
            );
        }
        "manifest-contract:coverage:full_set_close_report_is_fail_closed_for_mismatch_execution_failure_skip_and_partial_rows" =>
        {
            return close_declaration(case, Fault, CrateTest, ExpectedFailure, Fail, None);
        }
        "manifest-contract:coverage:full_set_close_report_refuses_wrong_manifest_cell_reason_scope_and_order" =>
        {
            return close_declaration(case, Api, CrateTest, ExpectedRefusal, Refuse, None);
        }
        "boundary:budget:publication_sum_overflow_is_typed_even_without_concrete_storage" => {
            return close_declaration(case, Resource, CrateTest, ExpectedRefusal, Refuse, None);
        }
        "mutation:construction:sensitive_observed_fuzz_corpus_never_echoes_through_any_rendering" =>
        {
            return close_declaration(case, Fuzz, CrateTest, MutationPartition, Refuse, None);
        }
        "property-metamorphic:logging:normalized_prefix_suffix_and_embedded_sensitive_aliases_refuse" =>
        {
            return close_declaration(case, Redaction, CrateTest, ExpectedRefusal, Refuse, None);
        }
        "property-metamorphic:state:exhaustive_cartesian_matrix_and_error_precedence"
        | "property-metamorphic:state:not_run_remaining_suffix_arithmetic_is_exhaustive_and_allocation_free" =>
        {
            return close_declaration(case, Model, CrateTest, Positive, Accept, None);
        }
        "compile-fail:state:no-unvalidated-state-candidate-as-terminal"
        | "compile-fail:state:private-not-run-basis-fields" => {
            return close_declaration(
                case,
                Typestate,
                CompileFailDoctest,
                ExpectedFailure,
                Fail,
                None,
            );
        }
        "property-metamorphic:extension:registry_is_typed_permutation_invariant_and_exact_set_reconstructible"
        | "property-metamorphic:path:set_adjudication_is_invariant_to_input_permutation"
        | "property-metamorphic:logging:canonical_event_and_log_roots_are_order_independent_but_mutation_sensitive" =>
        {
            return close_declaration(case, Metamorphic, CrateTest, Positive, Accept, None);
        }
        "unit:identity:wrapper_parser_checks_nominal_metadata_before_text"
        | "unit:diagnostic:typed_replacements_require_the_same_inline_or_retained_shape"
        | "manifest-contract:coverage:full_set_close_manifest_exact_reconstruction_rejects_all_sequence_and_semantic_mutants" =>
        {
            return close_declaration(case, Api, CrateTest, Positive, Accept, None);
        }
        _ => {}
    }

    match case.class() {
        BaseCoverageManifestClassV1::Unit => {
            if case.source_path() == "crates/fs-evidence-runner/src/logging.rs" {
                close_declaration(
                    case,
                    DetailedDeterministicLogging,
                    CrateTest,
                    Positive,
                    Accept,
                    None,
                )
            } else if case.source_path() == "crates/fs-evidence-runner/src/state.rs" {
                close_declaration(case, State, CrateTest, Positive, Accept, None)
            } else {
                close_declaration(case, Unit, CrateTest, Positive, Accept, None)
            }
        }
        BaseCoverageManifestClassV1::Boundary => {
            if case.source_path() == "crates/fs-evidence-runner/src/logging.rs" {
                close_declaration(
                    case,
                    DetailedDeterministicLogging,
                    CrateTest,
                    Positive,
                    Accept,
                    None,
                )
            } else if case.source_path() == "crates/fs-evidence-runner/src/state.rs" {
                close_declaration(case, State, CrateTest, Positive, Accept, None)
            } else {
                close_declaration(case, Boundary, CrateTest, Positive, Accept, None)
            }
        }
        BaseCoverageManifestClassV1::PropertyMetamorphic => {
            if case.source_path() == "crates/fs-evidence-runner/src/logging.rs" {
                close_declaration(
                    case,
                    DetailedDeterministicLogging,
                    CrateTest,
                    Positive,
                    Accept,
                    None,
                )
            } else if case.source_path() == "crates/fs-evidence-runner/src/state.rs" {
                close_declaration(case, Model, CrateTest, Positive, Accept, None)
            } else {
                close_declaration(case, Property, CrateTest, Positive, Accept, None)
            }
        }
        BaseCoverageManifestClassV1::SchemaDescriptor => {
            if case.source_path() == "crates/fs-evidence-runner/src/logging.rs" {
                close_declaration(
                    case,
                    DetailedDeterministicLogging,
                    CrateTest,
                    Positive,
                    Accept,
                    None,
                )
            } else if case.source_path() == "crates/fs-evidence-runner/src/state.rs" {
                close_declaration(case, State, CrateTest, Positive, Accept, None)
            } else {
                close_declaration(case, Literal, CrateTest, Positive, Accept, None)
            }
        }
        BaseCoverageManifestClassV1::Mutation => {
            if case.source_path() == "crates/fs-evidence-runner/src/logging.rs" {
                close_declaration(
                    case,
                    DetailedDeterministicLogging,
                    CrateTest,
                    MutationPartition,
                    Accept,
                    None,
                )
            } else {
                close_declaration(case, Mutation, CrateTest, MutationPartition, Accept, None)
            }
        }
        BaseCoverageManifestClassV1::NoMockIntegration => {
            close_declaration(case, Api, CrateTest, Positive, Accept, None)
        }
        BaseCoverageManifestClassV1::CompileFailDoctest => close_declaration(
            case,
            CompileFail,
            CompileFailDoctest,
            ExpectedFailure,
            Fail,
            None,
        ),
        BaseCoverageManifestClassV1::ManifestContract => {
            close_declaration(case, Unit, CrateTest, Positive, Accept, None)
        }
        BaseCoverageManifestClassV1::ProjectionE2e => {
            let row_id = projection_row_id_for_exact_source_case(id).ok_or_else(|| {
                refusal(
                    ConstructionErrorKindV2::UnknownCode,
                    "coverage.close.projection_source_case_id",
                    "one exact independently frozen projection journey row",
                    id,
                )
            })?;
            match row_id {
                "windows-unicode-alias" => close_declaration(
                    case,
                    ImmutableE2eContribution,
                    InProcessProjection,
                    UnsupportedPartition,
                    Unsupported,
                    Some(WindowsNonasciiAliasLocallyUnadjudicable),
                ),
                "identity-mutation" => close_declaration(
                    case,
                    Mutation,
                    InProcessProjection,
                    MutationPartition,
                    Accept,
                    None,
                ),
                "unknown-catalog-code"
                | "overlong-stable-token"
                | "reserved-content-store-prefix"
                | "budget-child-relation"
                | "publication-cross-cell"
                | "capability-extra-right"
                | "state-usage-in-lifecycle"
                | "diagnostic-rank-gap"
                | "atomic-result-presence" => close_declaration(
                    case,
                    ImmutableE2eContribution,
                    InProcessProjection,
                    ExpectedRefusal,
                    Refuse,
                    None,
                ),
                "catalog-literals"
                | "canonical-rational"
                | "logical-path"
                | "limit-catalog"
                | "budget-admission"
                | "publication-selection"
                | "capability-least-privilege"
                | "state-pass"
                | "diagnostic"
                | "no-claim-nominality"
                | "atomic-result"
                | "publication-storage"
                | "command-list" => close_declaration(
                    case,
                    ImmutableE2eContribution,
                    InProcessProjection,
                    Positive,
                    Accept,
                    None,
                ),
                _ => Err(refusal(
                    ConstructionErrorKindV2::UnknownCode,
                    "coverage.close.projection_row_id",
                    "one exact classified projection row ID",
                    row_id,
                )),
            }
        }
        BaseCoverageManifestClassV1::RuntimeLogging => close_declaration(
            case,
            DetailedDeterministicLogging,
            InProcessProjection,
            Positive,
            Accept,
            None,
        ),
        BaseCoverageManifestClassV1::SourceClosure => {
            if id == "source-closure:exact-positive" {
                close_declaration(
                    case,
                    SourceClosure,
                    InProcessProjection,
                    Positive,
                    Accept,
                    None,
                )
            } else {
                close_declaration(
                    case,
                    SourceClosure,
                    InProcessProjection,
                    MutationPartition,
                    Refuse,
                    None,
                )
            }
        }
        BaseCoverageManifestClassV1::ExternalE2eScript => close_downstream_declaration(
            case,
            ImmutableE2eContribution,
            InapplicablePartition,
            ReleaseExecutionDownstreamOwned,
        ),
        BaseCoverageManifestClassV1::ExternalMutation => close_downstream_declaration(
            case,
            Mutation,
            InapplicablePartition,
            ReleaseExecutionDownstreamOwned,
        ),
        BaseCoverageManifestClassV1::ExternalGovernance => close_downstream_declaration(
            case,
            SourceClosure,
            InapplicablePartition,
            ReleaseExecutionDownstreamOwned,
        ),
    }
}

fn validate_exact_full_source_manifest(
    expected: &BaseCoverageManifestV1,
    presented: &BaseCoverageManifestV1,
) -> Result<(), ConstructionErrorV2> {
    if presented.cases.len() < expected.cases.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.source_manifest",
            "the complete exact full source manifest",
            presented.cases.len(),
        ));
    }
    if presented.cases.len() > expected.cases.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close.source_manifest",
            "no source case beyond the complete exact full manifest",
            presented.cases.len(),
        ));
    }
    let expected_ids = expected
        .cases
        .iter()
        .map(BaseCoverageManifestCaseV1::id)
        .collect::<BTreeSet<_>>();
    for (index, (expected_case, presented_case)) in expected
        .cases
        .iter()
        .zip(presented.cases.iter())
        .enumerate()
    {
        if expected_case == presented_case {
            continue;
        }
        let kind = if expected_ids.contains(presented_case.id()) {
            ConstructionErrorKindV2::OutOfOrder
        } else {
            ConstructionErrorKindV2::Incompatible
        };
        return Err(refusal(
            kind,
            "coverage.close.source_manifest_case",
            "the exact source case at this full-manifest ordinal",
            format_args!("{index}:{}", presented_case.id()),
        ));
    }
    if presented.root != expected.root {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.source_manifest_root",
            "the root of the exact complete source sequence",
            presented.root.to_hex(),
        ));
    }
    Ok(())
}

fn validate_exact_close_declaration_sequence(
    expected: &BaseCoverageCloseManifestV1,
    presented: &[BaseCoverageCloseCellDeclarationV1],
) -> Result<(), ConstructionErrorV2> {
    if presented.len() < expected.cells.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.declarations",
            "one exact declaration per full source-manifest cell",
            presented.len(),
        ));
    }
    if presented.len() > expected.cells.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close.declarations",
            "no declaration beyond the full source manifest",
            presented.len(),
        ));
    }
    let mut ids = BTreeSet::new();
    for declaration in presented {
        if !ids.insert(declaration.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.close.source_case_id",
                "one declaration for every full source case",
                declaration.source_case_id(),
            ));
        }
    }
    for (index, (cell, declaration)) in expected.cells.iter().zip(presented).enumerate() {
        if cell.declaration == *declaration {
            continue;
        }
        let kind = if expected
            .cells
            .iter()
            .any(|expected_cell| expected_cell.source_case_id() == declaration.source_case_id())
        {
            if cell.source_case_id() != declaration.source_case_id() {
                ConstructionErrorKindV2::OutOfOrder
            } else {
                ConstructionErrorKindV2::Incompatible
            }
        } else {
            ConstructionErrorKindV2::UnknownCode
        };
        return Err(refusal(
            kind,
            "coverage.close.declaration",
            "the exact result-free declaration at this ordinal",
            format_args!("{index}:{}", declaration.source_case_id()),
        ));
    }
    Ok(())
}

fn close_manifest_from_declarations(
    source: &BaseCoverageManifestV1,
    declarations: Vec<BaseCoverageCloseCellDeclarationV1>,
) -> Result<BaseCoverageCloseManifestV1, ConstructionErrorV2> {
    if declarations.len() != source.cases.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.cell_count",
            "exactly one close declaration per full source case",
            declarations.len(),
        ));
    }
    let reason_registry_root = close_reason_registry_root()?;
    let mut cells = Vec::with_capacity(declarations.len());
    for (source_case, declaration) in source.cases.iter().zip(declarations) {
        if declaration.source_ordinal != source_case.ordinal()
            || declaration.source_case_id() != source_case.id()
            || declaration.source_class != source_case.class()
            || declaration.source_path() != source_case.source_path()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.source_binding",
                "one exact close declaration for the corresponding source cell",
                declaration.source_case_id(),
            ));
        }
        let root = close_cell_root(source.root, reason_registry_root, &declaration)?;
        cells.push(BaseCoverageCloseManifestCellV1 { declaration, root });
    }
    validate_close_corpus(&cells)?;
    let root = close_manifest_root(source.root, reason_registry_root, &cells)?;
    Ok(BaseCoverageCloseManifestV1 {
        source_manifest_root: source.root,
        reason_registry_root,
        cells: cells.into_boxed_slice(),
        root,
    })
}

fn validate_close_corpus(
    cells: &[BaseCoverageCloseManifestCellV1],
) -> Result<(), ConstructionErrorV2> {
    if cells.is_empty() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.cells",
            "the nonempty exact full-set close corpus",
            0,
        ));
    }
    let mut ids = BTreeSet::new();
    for (index, cell) in cells.iter().enumerate() {
        if cell.source_ordinal()
            != u32::try_from(index + 1).map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "coverage.close.source_ordinal",
                    "a one-based u32 ordinal",
                    index + 1,
                )
            })?
        {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.close.source_ordinal",
                "contiguous full source-manifest order",
                cell.source_ordinal(),
            ));
        }
        if !ids.insert(cell.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.close.source_case_id",
                "one close cell per full source case",
                cell.source_case_id(),
            ));
        }
        validate_close_partition_shape(
            cell.partition(),
            cell.expected_decision(),
            cell.expected_reason(),
        )?;
        validate_close_reason_scope(cell.expected_reason(), cell.execution_scope())?;
        let contribution_required = cell.execution_scope()
            == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution;
        if contribution_required != cell.downstream_contribution().is_some() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.downstream_contribution",
                "present exactly for immutable downstream contribution scope",
                cell.source_case_id(),
            ));
        }
    }
    for group in BaseCoverageCloseGroupV1::ALL {
        if !cells.iter().any(|cell| cell.group() == group) {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.group",
                "a nonzero corpus for every controlling group",
                group.stable_name(),
            ));
        }
    }
    for facet in BaseCoverageCloseFacetV1::ALL {
        let facet_cells = cells
            .iter()
            .filter(|cell| cell.facet() == facet)
            .collect::<Vec<_>>();
        if facet_cells.is_empty() {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.facet",
                "a nonzero applicable corpus or registered inapplicability declaration",
                facet.stable_name(),
            ));
        }
        let applicable_count = facet_cells
            .iter()
            .filter(|cell| cell.partition() != BaseCoverageClosePartitionV1::Inapplicable)
            .count();
        let exact_inapplicable_reason = match facet {
            BaseCoverageCloseFacetV1::Race => {
                Some(BaseCoverageCloseReasonCodeV1::RaceNotApplicablePureSingleThreadedValidator)
            }
            BaseCoverageCloseFacetV1::Trait => {
                Some(BaseCoverageCloseReasonCodeV1::TraitNotApplicableNoPublicTraitContract)
            }
            BaseCoverageCloseFacetV1::Cancellation => {
                Some(BaseCoverageCloseReasonCodeV1::CancellationNotApplicablePureBoundedValidator)
            }
            BaseCoverageCloseFacetV1::ReleaseBuiltNoMockE2e => {
                Some(BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned)
            }
            _ => None,
        };
        if let Some(reason) = exact_inapplicable_reason {
            if applicable_count != 0
                || facet_cells.len() != 1
                || facet_cells[0].expected_reason() != Some(reason)
            {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.close.inapplicable_facet",
                    "one exact registered declaration and no hidden eligible cell",
                    facet.stable_name(),
                ));
            }
        } else if applicable_count == 0 {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.applicable_facet",
                "a separately nonzero applicable corpus",
                facet.stable_name(),
            ));
        }
    }
    for partition in [
        BaseCoverageClosePartitionV1::Positive,
        BaseCoverageClosePartitionV1::ExpectedRefusal,
        BaseCoverageClosePartitionV1::ExpectedFailure,
        BaseCoverageClosePartitionV1::Mutation,
        BaseCoverageClosePartitionV1::Unsupported,
        BaseCoverageClosePartitionV1::Inapplicable,
    ] {
        if !cells.iter().any(|cell| cell.partition() == partition) {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.partition",
                "a separately nonzero exact corpus for every partition",
                partition.stable_name(),
            ));
        }
    }
    for reason in BaseCoverageCloseReasonCodeV1::ALL {
        if !cells
            .iter()
            .any(|cell| cell.expected_reason() == Some(reason))
        {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.close.reason",
                "at least one exact source cell for every registered reason",
                reason.code(),
            ));
        }
    }
    Ok(())
}

fn validate_close_results(
    manifest: &BaseCoverageCloseManifestV1,
    presented: &[BaseCoverageClosePresentedResultV1],
) -> Result<(), ConstructionErrorV2> {
    if presented.len() < manifest.cells.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close_results",
            "one result for every full close-manifest cell",
            presented.len(),
        ));
    }
    if presented.len() > manifest.cells.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close_results",
            "no result beyond the full close-manifest cells",
            presented.len(),
        ));
    }
    let mut ids = BTreeSet::new();
    for result in presented {
        if !ids.insert(result.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.close_result.source_case_id",
                "one result per full close-manifest cell",
                result.source_case_id(),
            ));
        }
    }
    for (index, (cell, result)) in manifest.cells.iter().zip(presented).enumerate() {
        if result.close_manifest_root != manifest.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close_result.close_manifest_root",
                "the exact current full close-manifest root",
                result.close_manifest_root.to_hex(),
            ));
        }
        if result.source_case_id() != cell.source_case_id() {
            let mapped = manifest.cell(result.source_case_id()).is_some();
            return Err(refusal(
                if mapped {
                    ConstructionErrorKindV2::OutOfOrder
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "coverage.close_result.source_case_id",
                "the exact full-manifest source ID at this ordinal",
                format_args!("{index}:{}", result.source_case_id()),
            ));
        }
        if result.cell_root != cell.root
            || result.group != cell.group()
            || result.facet != cell.facet()
            || result.execution_scope != cell.execution_scope()
            || result.partition != cell.partition()
            || result.expected_decision != cell.expected_decision()
            || result.expected_reason != cell.expected_reason()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close_result.expected_classification",
                "the exact cell root, group, facet, scope, partition, decision, and registered reason",
                result.source_case_id(),
            ));
        }
        validate_close_partition_shape(
            result.partition,
            result.expected_decision,
            result.expected_reason,
        )?;
        validate_close_reason_scope(result.expected_reason, result.execution_scope)?;
        validate_close_evidence_for_cell(manifest, cell, &result.evidence)?;
        match result.status {
            BaseCoverageCloseResultStatusV1::Matched => {
                if result.observed_decision != Some(result.expected_decision)
                    || result.observed_reason != result.expected_reason
                {
                    return Err(refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "coverage.close_result.matched_observation",
                        "the exact expected decision and registered reason",
                        result.source_case_id(),
                    ));
                }
            }
            BaseCoverageCloseResultStatusV1::UnexpectedMismatch => {
                if result.observed_decision.is_none()
                    || (result.observed_decision == Some(result.expected_decision)
                        && result.observed_reason == result.expected_reason)
                {
                    return Err(refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "coverage.close_result.unexpected_mismatch",
                        "one actually divergent observed classification",
                        result.source_case_id(),
                    ));
                }
            }
            BaseCoverageCloseResultStatusV1::ExecutionFailure
            | BaseCoverageCloseResultStatusV1::UnexplainedSkip => {
                if result.observed_decision.is_some() || result.observed_reason.is_some() {
                    return Err(refusal(
                        ConstructionErrorKindV2::Unexpected,
                        "coverage.close_result.observed_classification",
                        "no fabricated classification after failure or skip",
                        result.source_case_id(),
                    ));
                }
            }
        }
        let reconstructed_root = close_presented_result_root(
            result.close_manifest_root,
            result.cell_root,
            result.source_case_id(),
            result.group,
            result.facet,
            result.execution_scope,
            result.partition,
            result.expected_decision,
            result.expected_reason,
            result.status,
            result.observed_decision,
            result.observed_reason,
            &result.evidence,
        )?;
        if reconstructed_root != result.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close_result.root",
                "the root reconstructed from every presented result field",
                result.source_case_id(),
            ));
        }
    }
    Ok(())
}

fn checked_adversarial_total(
    expected_refusals: usize,
    expected_failures: usize,
    mutations: usize,
    field: &'static str,
) -> Result<u32, ConstructionErrorV2> {
    let refusals = u32::try_from(expected_refusals).map_err(|_| {
        refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "three independently checked u32 adversarial counts",
            expected_refusals,
        )
    })?;
    let failures = u32::try_from(expected_failures).map_err(|_| {
        refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "three independently checked u32 adversarial counts",
            expected_failures,
        )
    })?;
    let mutations = u32::try_from(mutations).map_err(|_| {
        refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "three independently checked u32 adversarial counts",
            mutations,
        )
    })?;
    refusals
        .checked_add(failures)
        .and_then(|sum| sum.checked_add(mutations))
        .ok_or_else(|| {
            refusal(
                ConstructionErrorKindV2::ArithmeticOverflow,
                field,
                "the checked sum of refusal, failure, and mutation counts",
                u64::from(refusals) + u64::from(failures) + u64::from(mutations),
            )
        })
}

fn validate_close_capability_set(
    field: &'static str,
    values: &[StableTokenV2],
) -> Result<(), ConstructionErrorV2> {
    if values.len() > 64 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            "coverage.close.five_explicits.capability_set",
            "at most 64 exact capability IDs per set",
            values.len(),
        ));
    }
    for pair in values.windows(2) {
        if pair[0] == pair[1] {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                field,
                "unique capability IDs in canonical order",
                pair[1].as_str(),
            ));
        }
        if pair[0] > pair[1] {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                field,
                "capability IDs in canonical lexical order",
                pair[1].as_str(),
            ));
        }
    }
    Ok(())
}

fn exact_close_capability_definition(
    code: u16,
) -> Option<(
    &'static str,
    BaseCoverageCloseCapabilityPolicyV1,
    &'static str,
)> {
    let index = usize::from(code.checked_sub(1)?);
    let definition = EXACT_CLOSE_CAPABILITY_DEFINITIONS_V1.get(index)?;
    Some((definition.stable_id, definition.policy, definition.no_claim))
}

fn validate_close_capability_registry_rows(
    rows: &[BaseCoverageCloseCapabilityDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    if rows.len() < BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.capability_registry.rows",
            "the exact five source capability descriptors",
            rows.len(),
        ));
    }
    if rows.len() > BASE_COVERAGE_CLOSE_CAPABILITY_DESCRIPTOR_COUNT_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close.capability_registry.rows",
            "no capability descriptor beyond the exact source registry",
            rows.len(),
        ));
    }
    for (index, row) in rows.iter().enumerate() {
        let Some((stable_id, policy, no_claim)) = exact_close_capability_definition(row.id.code())
        else {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "coverage.close.capability_registry.id",
                "one of the exact five registered capability IDs",
                row.id.code(),
            ));
        };
        if row.stable_id.as_str() != stable_id
            || row.owner.as_str() != BASE_COVERAGE_CLOSE_CAPABILITY_OWNER_V1
            || row.policy != policy
            || row.no_claim.as_str() != no_claim
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_registry.row",
                "the exact stable ID, owner, policy, and no-claim for this registry ID",
                row.id.code(),
            ));
        }
        let expected_root = close_capability_descriptor_root(
            row.id,
            &row.stable_id,
            &row.owner,
            row.policy,
            &row.no_claim,
        )?;
        if row.root != expected_root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_registry.descriptor_root",
                "the nominal root of the complete exact descriptor row",
                row.id.code(),
            ));
        }
        if let Some(previous) = index.checked_sub(1).map(|previous| &rows[previous]) {
            if previous.id == row.id {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.capability_registry.id",
                    "unique capability IDs in exact registry order",
                    row.id.code(),
                ));
            }
            if previous.id > row.id {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.close.capability_registry.id",
                    "strict ascending capability registry IDs",
                    row.id.code(),
                ));
            }
        }
        for previous in &rows[..index] {
            if previous.stable_id == row.stable_id {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.capability_registry.stable_id",
                    "unique stable capability IDs",
                    row.id.code(),
                ));
            }
            if previous.root == row.root {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "coverage.close.capability_registry.descriptor_root",
                    "unique nominal descriptor roots",
                    row.id.code(),
                ));
            }
        }
    }
    Ok(())
}

fn validate_close_capability_profile_rows(
    capability_registry: &BaseCoverageCloseCapabilityRegistryV1,
    rows: &[BaseCoverageCloseCapabilityProfileDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    if rows.len() < BaseCoverageCloseCapabilityProfileV1::ALL.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.capability_profile_registry.rows",
            "the exact five source capability profiles",
            rows.len(),
        ));
    }
    if rows.len() > BaseCoverageCloseCapabilityProfileV1::ALL.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close.capability_profile_registry.rows",
            "no profile beyond the exact source capability-profile registry",
            rows.len(),
        ));
    }
    for (index, row) in rows.iter().enumerate() {
        let expected_profile = BaseCoverageCloseCapabilityProfileV1::ALL[index];
        if row.profile != expected_profile {
            return Err(refusal(
                if rows[..index]
                    .iter()
                    .any(|previous| previous.profile == row.profile)
                {
                    ConstructionErrorKindV2::Duplicate
                } else {
                    ConstructionErrorKindV2::OutOfOrder
                },
                "coverage.close.capability_profile_registry.profile",
                "each exact capability profile once in registry order",
                row.profile.code(),
            ));
        }
        if row.stable_id.as_str() != expected_profile.stable_id()
            || row.no_claim.as_str() != BASE_COVERAGE_CLOSE_CAPABILITY_CONTRACT_NO_CLAIM_V1
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_profile_registry.row",
                "the exact stable profile ID and no-claim",
                row.profile.code(),
            ));
        }
        validate_close_capability_id_set(
            "coverage.close.capability_profile_registry.required",
            capability_registry,
            &row.required,
        )?;
        validate_close_capability_id_set(
            "coverage.close.capability_profile_registry.permitted",
            capability_registry,
            &row.permitted,
        )?;
        let exact = expected_profile.required_codes();
        if row.required.len() != exact.len()
            || row.permitted.len() != exact.len()
            || row.required.as_ref() != row.permitted.as_ref()
            || !row
                .required
                .iter()
                .zip(exact)
                .all(|(observed, expected)| observed.code() == *expected)
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.capability_profile_registry.sets",
                "required equals permitted and matches the exact profile table",
                row.profile.code(),
            ));
        }
    }
    Ok(())
}

fn validate_close_capability_id_set(
    field: &'static str,
    registry: &BaseCoverageCloseCapabilityRegistryV1,
    values: &[BaseCoverageCloseCapabilityIdV1],
) -> Result<(), ConstructionErrorV2> {
    if values.len() > BASE_COVERAGE_CLOSE_CAPABILITY_SET_MAX_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "at most the exact five base semantic capability IDs",
            values.len(),
        ));
    }
    for (index, value) in values.iter().enumerate() {
        if registry.descriptor(*value).is_none() {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                field,
                "an ID registered by the exact close capability registry",
                value.code(),
            ));
        }
        if let Some(previous) = index.checked_sub(1).map(|previous| values[previous]) {
            if previous == *value {
                return Err(refusal(
                    ConstructionErrorKindV2::Duplicate,
                    field,
                    "unique capability IDs in canonical registry order",
                    value.code(),
                ));
            }
            if previous > *value {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    field,
                    "capability IDs in canonical registry order",
                    value.code(),
                ));
            }
        }
    }
    Ok(())
}

fn is_base_close_capability_stable_id(value: &str) -> bool {
    EXACT_CLOSE_CAPABILITY_DEFINITIONS_V1
        .iter()
        .any(|definition| definition.stable_id == value)
}

fn close_capability_set_is_subset(
    subset: &[BaseCoverageCloseCapabilityIdV1],
    superset: &[BaseCoverageCloseCapabilityIdV1],
) -> bool {
    subset.iter().all(|id| superset.binary_search(id).is_ok())
}

fn validate_close_numeric_explicit_sequence(
    field: &'static str,
    values: &[BaseCoverageCloseNumericExplicitV1],
) -> Result<(), ConstructionErrorV2> {
    if values.len() > BASE_COVERAGE_CLOSE_NUMERIC_EXPLICIT_MAX_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "an explicit exact-empty through 64-row typed numeric profile",
            values.len(),
        ));
    }
    for pair in values.windows(2) {
        if pair[0].name == pair[1].name {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                field,
                "unique numeric explicit names in canonical order",
                pair[1].name.as_str(),
            ));
        }
        if pair[0].name > pair[1].name {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                field,
                "numeric explicit names in canonical lexical order",
                pair[1].name.as_str(),
            ));
        }
    }
    Ok(())
}

// The current base leaf's property, metamorphic, mutation, and fuzz-labelled
// tests are deterministic exhaustive corpora; none consumes an RNG. This
// exact-set table is intentionally source-authoritative rather than inferred
// from a coverage facet. A future genuinely randomized cell must add its
// stable ID here and bind explicit material plus generator/minimizer versions.
const CLOSE_SEMANTICALLY_SEEDED_CASE_IDS_V1: &[&str] = &[];

fn close_case_uses_semantic_seed(source_case_id: &str) -> bool {
    CLOSE_SEMANTICALLY_SEEDED_CASE_IDS_V1
        .binary_search(&source_case_id)
        .is_ok()
}

fn expected_close_target_and_profile(
    execution_scope: BaseCoverageCloseExecutionScopeV1,
) -> (BaseCoverageCloseTargetV1, BaseCoverageCloseProfileV1) {
    match execution_scope {
        BaseCoverageCloseExecutionScopeV1::CrateTest => (
            BaseCoverageCloseTargetV1::TargetIndependentPureValidation,
            BaseCoverageCloseProfileV1::CrateTest,
        ),
        BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => (
            BaseCoverageCloseTargetV1::DeclaredRustTarget,
            BaseCoverageCloseProfileV1::CompileFailDoctest,
        ),
        BaseCoverageCloseExecutionScopeV1::InProcessProjection => (
            BaseCoverageCloseTargetV1::DeclaredHostTarget,
            BaseCoverageCloseProfileV1::InProcessProjection,
        ),
        BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => (
            BaseCoverageCloseTargetV1::DownstreamPlatformMatrix,
            BaseCoverageCloseProfileV1::DownstreamRelease,
        ),
        BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => (
            BaseCoverageCloseTargetV1::TargetIndependentPureValidation,
            BaseCoverageCloseProfileV1::ApplicabilityDeclaration,
        ),
    }
}

fn validate_close_five_explicits_for_declaration(
    source_case_id: &str,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    downstream_contribution: Option<&BaseCoverageCloseDownstreamContributionV1>,
    five: &BaseCoverageCloseFiveExplicitsV1,
) -> Result<(), ConstructionErrorV2> {
    if five.versions.api_generation != RUNNER_SPEC_V2_API_GENERATION
        || five.versions.wire_version != RUNNER_V2_WIRE_VERSION
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.api_wire",
            "RunnerSpecV2 API generation 2 and wire version 1",
            five.versions.api_generation.code(),
        ));
    }
    let (target, profile) = expected_close_target_and_profile(execution_scope);
    if five.versions.target != target || five.versions.profile != profile {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.target_profile",
            "the exact target/profile pair for the execution scope",
            five.versions.profile.code(),
        ));
    }
    if !five.numeric_inputs.is_empty()
        || !five.numeric_grants.is_empty()
        || !five.numeric_observations.is_empty()
    {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.close.five_explicits.numeric_profiles",
            "the current source-authoritative exact-empty semantic input, grant, and observation profiles",
            five.numeric_inputs.len() + five.numeric_grants.len() + five.numeric_observations.len(),
        ));
    }
    if close_case_uses_semantic_seed(source_case_id) != five.seed.material().is_some() {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.seed_applicability",
            "material plus generator/minimizer exactly for a source-registered randomized cell",
            five.seed.material().is_some(),
        ));
    }
    if let BaseCoverageCloseSeedExplicitV1::Applicable {
        generator_version,
        minimizer_version,
        ..
    } = &five.seed
        && (generator_version.as_str() == minimizer_version.as_str()
            || generator_version.as_str().is_empty()
            || minimizer_version.as_str().is_empty())
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.seed_versions",
            "distinct exact generator and minimizer versions",
            generator_version.as_str(),
        ));
    }
    if matches!(
        &five.seed,
        BaseCoverageCloseSeedExplicitV1::Inapplicable {
            reason: SeedInapplicableCodeV1::NoRandomnessByContract
        }
    ) && close_case_uses_semantic_seed(source_case_id)
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.seed_inapplicable",
            "inapplicable only for a cell that consumes no randomness",
            source_case_id,
        ));
    }

    let capability_sets = &five.capabilities;
    if let Some(contribution) = downstream_contribution {
        if five.budgets.profile != BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution
            || five.budgets != contribution.budgets.resolved
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.downstream_budget_profile",
                "the exact source-named downstream profile and contribution-resolved rows",
                five.budgets.profile.name(),
            ));
        }
        if capability_sets.required.len() != 1
            || capability_sets.required[0].as_str() != contribution.downstream_owner()
            || !capability_sets.granted.is_empty()
            || !capability_sets.observed.is_empty()
            || !capability_sets.returned.is_empty()
            || !capability_sets.revoked.is_empty()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.downstream_capabilities",
                "one exact downstream requirement and empty unobserved effect sets",
                capability_sets.required.len(),
            ));
        }
        if five.versions.source_root != contribution.source_root
            || five.versions.build_root != contribution.build_root
            || five.versions.schema_root != contribution.schema_root
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.downstream_versions",
                "the exact contribution source, build, and schema roots",
                contribution.downstream_owner(),
            ));
        }
    } else {
        if five.budgets.profile != BaseCoverageCloseBudgetProfileV1::LocalSourceValidation
            || five.budgets != frozen_local_close_budget_set()?
        {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.local_budget_profile",
                "the exact source-named local profile and resolved rows",
                five.budgets.profile.name(),
            ));
        }
        if !capability_sets.required.is_empty()
            || !capability_sets.granted.is_empty()
            || !capability_sets.observed.is_empty()
            || !capability_sets.returned.is_empty()
            || !capability_sets.revoked.is_empty()
        {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.close.five_explicits.local_capabilities",
                "explicit empty capability sets for pure local validation",
                capability_sets.required.len(),
            ));
        }
    }
    Ok(())
}

fn close_five_reference_root(
    domain: &'static str,
    source_ordinal: u32,
    source_case_id: &str,
    source_class: BaseCoverageManifestClassV1,
    source_path: &str,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEFIVEREF\x01", 2 * 1024)?;
    frame.push_u32(
        "coverage.close.five_reference.source_ordinal",
        source_ordinal,
    )?;
    frame.push_str(
        "coverage.close.five_reference.source_case_id",
        source_case_id,
    )?;
    frame.push_u16(
        "coverage.close.five_reference.source_class",
        source_class.code(),
    )?;
    frame.push_str("coverage.close.five_reference.source_path", source_path)?;
    frame.push_u16("coverage.close.five_reference.facet", facet.code())?;
    frame.push_u16(
        "coverage.close.five_reference.execution_scope",
        execution_scope.code(),
    )?;
    Ok(frame.root(domain))
}

fn fixed_close_budget_unit(
    unit: LogicalUnitV2,
) -> Result<BaseCoverageCloseLogicalUnitReferenceV1, ConstructionErrorV2> {
    BaseCoverageCloseLogicalUnitReferenceV1::fixed(unit)
}

fn frozen_local_close_budget_set() -> Result<BaseCoverageCloseBudgetSetV1, ConstructionErrorV2> {
    BaseCoverageCloseBudgetSetV1::new(
        BaseCoverageCloseBudgetProfileV1::LocalSourceValidation,
        vec![
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                fixed_close_budget_unit(LogicalUnitV2::Nanoseconds)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Memory,
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                fixed_close_budget_unit(LogicalUnitV2::LogicalBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::LogicalWork,
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                fixed_close_budget_unit(LogicalUnitV2::Operations)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Processes,
                BaseCoverageCloseBudgetValueV1::U32(1),
                BaseCoverageCloseBudgetValueV1::U32(0),
                fixed_close_budget_unit(LogicalUnitV2::Count)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Artifacts,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Output,
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Logs,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
        ],
    )
}

fn frozen_downstream_close_budget_set() -> Result<BaseCoverageCloseBudgetSetV1, ConstructionErrorV2>
{
    BaseCoverageCloseBudgetSetV1::new(
        BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
        vec![
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                fixed_close_budget_unit(LogicalUnitV2::Nanoseconds)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Memory,
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                fixed_close_budget_unit(LogicalUnitV2::LogicalBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::LogicalWork,
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                fixed_close_budget_unit(LogicalUnitV2::Operations)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Processes,
                BaseCoverageCloseBudgetValueV1::U32(8),
                BaseCoverageCloseBudgetValueV1::U32(6),
                fixed_close_budget_unit(LogicalUnitV2::Count)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Artifacts,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Output,
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
            BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::Logs,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                fixed_close_budget_unit(LogicalUnitV2::EncodedBytes)?,
            )?,
        ],
    )
}

/// Reconstruct one exact downstream contribution budget with an independently
/// declared child-process shape.
///
/// This crate-private surface lets source-frozen V2 projection declarations
/// reuse the unchanged V1 seven-axis budget vocabulary without exposing a
/// public constructor that could pretend to be a source-owned route.
pub(crate) fn frozen_downstream_close_contribution_budgets_v1(
    max_child_processes: u32,
    max_parallel_children: u32,
) -> Result<BaseCoverageCloseContributionBudgetsV1, ConstructionErrorV2> {
    BaseCoverageCloseContributionBudgetsV1::new(
        frozen_downstream_close_budget_set()?,
        max_child_processes,
        max_parallel_children,
    )
}

#[allow(
    clippy::too_many_arguments,
    reason = "the frozen Five Explicits derive from every stable cell identity axis"
)]
fn frozen_close_five_explicits_v1(
    source_ordinal: u32,
    source_case_id: &str,
    source_class: BaseCoverageManifestClassV1,
    source_path: &str,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    downstream_contribution: Option<&BaseCoverageCloseDownstreamContributionV1>,
) -> Result<BaseCoverageCloseFiveExplicitsV1, ConstructionErrorV2> {
    const SCHEMA_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-five-schema-version.v1";
    const SOURCE_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-five-source-version.v1";
    const BUILD_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-five-build-version.v1";
    const TOOLCHAIN_DOMAIN: &str =
        "org.frankensim.fs-evidence-runner.close-five-toolchain-version.v1";
    const FEATURE_DOMAIN: &str = "org.frankensim.fs-evidence-runner.close-five-feature-set.v1";

    let reference = |domain| {
        close_five_reference_root(
            domain,
            source_ordinal,
            source_case_id,
            source_class,
            source_path,
            facet,
            execution_scope,
        )
    };
    let seed = if close_case_uses_semantic_seed(source_case_id) {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.close.five_explicits.seed_material",
            "source-authoritative explicit seed material for every registered randomized cell",
            source_case_id,
        ));
    } else {
        BaseCoverageCloseSeedExplicitV1::Inapplicable {
            reason: SeedInapplicableCodeV1::NoRandomnessByContract,
        }
    };

    let budgets = downstream_contribution
        .map_or_else(frozen_local_close_budget_set, |contribution| {
            Ok(contribution.budgets.resolved)
        })?;

    let source_root = if let Some(contribution) = downstream_contribution {
        contribution.source_root.clone()
    } else {
        let root = reference(SOURCE_DOMAIN)?;
        SourceIdentityRootV2::parse_presented(
            SourceIdentityRootV2::DESCRIPTOR.role(),
            SourceIdentityRootV2::DESCRIPTOR.domain(),
            &root.to_hex(),
        )
        .map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.source_root",
                "one exact presented source identity",
                source_case_id,
            )
        })?
    };
    let build_root = if let Some(contribution) = downstream_contribution {
        contribution.build_root.clone()
    } else {
        let root = reference(BUILD_DOMAIN)?;
        BuildIdentityRootV2::parse_presented(
            BuildIdentityRootV2::DESCRIPTOR.role(),
            BuildIdentityRootV2::DESCRIPTOR.domain(),
            &root.to_hex(),
        )
        .map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.build_root",
                "one exact presented build identity",
                source_case_id,
            )
        })?
    };
    let toolchain_reference = reference(TOOLCHAIN_DOMAIN)?;
    let toolchain_root = ToolchainIdentityRootV2::parse_presented(
        ToolchainIdentityRootV2::DESCRIPTOR.role(),
        ToolchainIdentityRootV2::DESCRIPTOR.domain(),
        &toolchain_reference.to_hex(),
    )
    .map_err(|_| {
        refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close.five_explicits.toolchain_root",
            "one exact presented toolchain identity",
            source_case_id,
        )
    })?;
    let (target, profile) = expected_close_target_and_profile(execution_scope);
    let versions = BaseCoverageCloseVersionSetV1::new(
        RUNNER_SPEC_V2_API_GENERATION,
        RUNNER_V2_WIRE_VERSION,
        downstream_contribution.map_or(reference(SCHEMA_DOMAIN)?, |value| value.schema_root),
        source_root,
        build_root,
        toolchain_root,
        target,
        profile,
        reference(FEATURE_DOMAIN)?,
    );

    let required = downstream_contribution
        .map(|contribution| {
            StableTokenV2::new(contribution.downstream_owner()).map(|token| vec![token])
        })
        .transpose()
        .map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.required_capability",
                "one exact stable downstream owner capability ID",
                source_case_id,
            )
        })?
        .unwrap_or_default();
    let capabilities =
        BaseCoverageCloseCapabilitySetsV1::new(required, vec![], vec![], vec![], vec![])?;
    // The current source oracle declares no semantic numeric case inputs,
    // numeric grants, or expected numeric observations. Classification
    // ordinals and result-row cardinality are bound elsewhere and must not
    // masquerade as semantics.
    let numeric_inputs = vec![];
    let numeric_grants = vec![];
    let numeric_observations = vec![];
    let no_claim = StableTokenV2::new("five-explicits-prove-no-execution-science-or-admission")
        .map_err(|_| {
            refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.close.five_explicits.no_claim",
                "one exact stable no-claim token",
                source_ordinal,
            )
        })?;
    BaseCoverageCloseFiveExplicitsV1::new(
        numeric_inputs,
        numeric_grants,
        numeric_observations,
        seed,
        budgets,
        versions,
        capabilities,
        no_claim,
    )
}

fn push_close_logical_unit(
    frame: &mut CanonicalFrameV1,
    tag_field: &'static str,
    presence_field: &'static str,
    id_field: &'static str,
    unit: LogicalUnitV2,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(tag_field, unit.tag())?;
    frame.push_presence(presence_field, unit.registered_id().is_some())?;
    if let Some(id) = unit.registered_id() {
        frame.push_u16(id_field, id)?;
    }
    Ok(())
}

fn push_close_logical_unit_reference(
    frame: &mut CanonicalFrameV1,
    reference: BaseCoverageCloseLogicalUnitReferenceV1,
) -> Result<(), ConstructionErrorV2> {
    push_close_logical_unit(
        frame,
        "coverage.close.five.logical_unit_tag",
        "coverage.close.five.logical_unit_id_present",
        "coverage.close.five.logical_unit_id",
        reference.unit,
    )?;
    frame.push_presence(
        "coverage.close.five.logical_unit_registry_identity_present",
        reference.registry_identity.is_some(),
    )?;
    if let Some(identity) = reference.registry_identity {
        frame.push_bytes(
            "coverage.close.five.logical_unit_registry_identity",
            identity.as_bytes(),
        )?;
    }
    Ok(())
}

fn push_close_numeric_value(
    frame: &mut CanonicalFrameV1,
    value: &NumericValueV2,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16("coverage.close.five.numeric_value_tag", value.wire_tag())?;
    match value {
        NumericValueV2::I8(value) => frame.push_i8("coverage.close.five.numeric_i8", *value)?,
        NumericValueV2::I16(value) => frame.push_i16("coverage.close.five.numeric_i16", *value)?,
        NumericValueV2::I32(value) => frame.push_i32("coverage.close.five.numeric_i32", *value)?,
        NumericValueV2::I64(value) => frame.push_i64("coverage.close.five.numeric_i64", *value)?,
        NumericValueV2::I128(value) => {
            frame.push_i128("coverage.close.five.numeric_i128", *value)?;
        }
        NumericValueV2::U8(value) => frame.push_u8("coverage.close.five.numeric_u8", *value)?,
        NumericValueV2::U16(value) => frame.push_u16("coverage.close.five.numeric_u16", *value)?,
        NumericValueV2::U32(value) => frame.push_u32("coverage.close.five.numeric_u32", *value)?,
        NumericValueV2::U64(value) => frame.push_u64("coverage.close.five.numeric_u64", *value)?,
        NumericValueV2::U128(value) => {
            frame.push_u128("coverage.close.five.numeric_u128", *value)?;
        }
        NumericValueV2::Rational(value) => {
            frame.push_i128(
                "coverage.close.five.numeric_rational_numerator",
                value.numerator(),
            )?;
            frame.push_u128(
                "coverage.close.five.numeric_rational_denominator",
                value.denominator(),
            )?;
        }
        NumericValueV2::Decimal(value) => {
            frame.push_i128(
                "coverage.close.five.numeric_decimal_coefficient",
                value.coefficient(),
            )?;
            frame.push_i32("coverage.close.five.numeric_decimal_scale", value.scale())?;
        }
        NumericValueV2::F32Bits(value) => {
            frame.push_u32("coverage.close.five.numeric_f32_bits", value.bits())?;
        }
        NumericValueV2::F64Bits(value) => {
            frame.push_u64("coverage.close.five.numeric_f64_bits", value.bits())?;
        }
    }
    Ok(())
}

fn push_close_numeric_unit(
    frame: &mut CanonicalFrameV1,
    unit: BaseCoverageCloseNumericUnitV1,
) -> Result<(), ConstructionErrorV2> {
    match unit {
        BaseCoverageCloseNumericUnitV1::Physical(unit) => {
            frame.push_u16("coverage.close.five.numeric_unit_domain", 1)?;
            let scale = unit.scale();
            frame.push_i128(
                "coverage.close.five.numeric_physical_scale_numerator",
                scale.numerator(),
            )?;
            frame.push_u128(
                "coverage.close.five.numeric_physical_scale_denominator",
                scale.denominator(),
            )?;
            for exponent in unit.exponents().as_array() {
                frame.push_i16(
                    "coverage.close.five.numeric_physical_dimension_exponent",
                    *exponent,
                )?;
            }
        }
        BaseCoverageCloseNumericUnitV1::Logical(reference) => {
            frame.push_u16("coverage.close.five.numeric_unit_domain", 2)?;
            push_close_logical_unit_reference(frame, reference)?;
        }
    }
    Ok(())
}

fn push_close_capability_ids(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    values: &[BaseCoverageCloseCapabilityIdV1],
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(
        field,
        u16::try_from(values.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                field,
                "a u16 bounded capability-ID count",
                values.len(),
            )
        })?,
    )?;
    for value in values {
        frame.push_u16(field, value.code())?;
    }
    Ok(())
}

fn close_capability_descriptor_root(
    id: BaseCoverageCloseCapabilityIdV1,
    stable_id: &StableTokenV2,
    owner: &StableTokenV2,
    policy: BaseCoverageCloseCapabilityPolicyV1,
    no_claim: &StableTokenV2,
) -> Result<BaseCoverageCloseCapabilityDescriptorRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECAPDESC\x01", 1024)?;
    frame.push_u16(
        "coverage.close.capability_descriptor.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.capability_descriptor.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_descriptor.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_u16("coverage.close.capability_descriptor.id", id.code())?;
    frame.push_str(
        "coverage.close.capability_descriptor.stable_id",
        stable_id.as_str(),
    )?;
    frame.push_str("coverage.close.capability_descriptor.owner", owner.as_str())?;
    frame.push_u16(
        "coverage.close.capability_descriptor.policy_code",
        policy.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_descriptor.policy_name",
        policy.stable_name(),
    )?;
    frame.push_str(
        "coverage.close.capability_descriptor.no_claim",
        no_claim.as_str(),
    )?;
    Ok(
        BaseCoverageCloseCapabilityDescriptorRootV1::from_content_hash(
            frame.root(BaseCoverageCloseCapabilityDescriptorRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_registered_extension_capability_descriptor_root(
    id: BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
    stable_id: &BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1,
    owner: &BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1,
    scope: &BaseCoverageCloseRegisteredExtensionCapabilityScopeV1,
    no_claim: &BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1,
) -> Result<BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEEXTCAPDESC\x01", 2 * 1024)?;
    frame.push_u16(
        "coverage.close.extension_capability_descriptor.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_descriptor.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_descriptor.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_descriptor.id",
        id.code(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_descriptor.stable_id",
        stable_id.as_str(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_descriptor.owner",
        owner.as_str(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_descriptor.scope",
        scope.as_str(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_descriptor.no_claim",
        no_claim.as_str(),
    )?;
    Ok(
        BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1::from_content_hash(
            frame.root(
                BaseCoverageCloseRegisteredExtensionCapabilityDescriptorRootV1::DESCRIPTOR.domain(),
            ),
        ),
    )
}

fn close_registered_extension_capability_registry_root(
    rows: &[BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1],
) -> Result<BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEEXTCAPREG\x01", 8 * 1024)?;
    frame.push_u16(
        "coverage.close.extension_capability_registry.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_registry.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_registry.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_registry.count",
        u16::try_from(rows.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.extension_capability_registry.count",
                "a u16 bounded extension capability count",
                rows.len(),
            )
        })?,
    )?;
    for row in rows {
        frame.push_bytes(
            "coverage.close.extension_capability_registry.descriptor_root",
            row.root().content_hash().as_bytes(),
        )?;
    }
    frame.push_str(
        "coverage.close.extension_capability_registry.no_claim",
        BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1::DESCRIPTOR.no_claim(),
    )?;
    Ok(
        BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1::from_content_hash(
            frame.root(
                BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1::DESCRIPTOR.domain(),
            ),
        ),
    )
}

fn close_registered_extension_capability_set_root(
    registry_root: BaseCoverageCloseRegisteredExtensionCapabilityRegistryRootV1,
    values: &[BaseCoverageCloseRegisteredExtensionCapabilityIdV1],
) -> Result<BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEEXTCAPSET\x01", 2 * 1024)?;
    frame.push_u16(
        "coverage.close.extension_capability_set.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_set.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.extension_capability_set.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_bytes(
        "coverage.close.extension_capability_set.registry_root",
        registry_root.content_hash().as_bytes(),
    )?;
    frame.push_u16(
        "coverage.close.extension_capability_set.count",
        u16::try_from(values.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.extension_capability_set.count",
                "a u16 bounded extension capability count",
                values.len(),
            )
        })?,
    )?;
    for value in values {
        frame.push_u16(
            "coverage.close.extension_capability_set.value",
            value.code(),
        )?;
    }
    frame.push_str(
        "coverage.close.extension_capability_set.no_claim",
        BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1::DESCRIPTOR.no_claim(),
    )?;
    Ok(
        BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1::from_content_hash(
            frame
                .root(BaseCoverageCloseRegisteredExtensionCapabilitySetRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_runtime_observation_disposition_root(
    disposition: RuntimeObservationDispositionV1,
) -> Result<BaseCoverageCloseRuntimeObservationDispositionRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEOBSDISPOSITION\x01", 128)?;
    frame.push_u16(
        "coverage.close.runtime_observation.disposition",
        disposition.code(),
    )?;
    Ok(
        BaseCoverageCloseRuntimeObservationDispositionRootV1::from_content_hash(
            frame.root(BaseCoverageCloseRuntimeObservationDispositionRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_not_observed_reason_registry_root(
    descriptors: &[NotObservedReasonDescriptorV1],
) -> Result<BaseCoverageCloseNotObservedReasonRegistryRootV1, ConstructionErrorV2> {
    if descriptors.len() != NotObservedReasonV1::ALL.len() {
        return Err(refusal(
            if descriptors.len() < NotObservedReasonV1::ALL.len() {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "coverage.close.not_observed_reason_registry.rows",
            "the exact three source-owned NotObserved reason rows",
            descriptors.len(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSCLOSENOTOBSREASONS\x01", 4 * 1024)?;
    frame.push_u32(
        "coverage.close.not_observed_reason_registry.count",
        u32::try_from(descriptors.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.not_observed_reason_registry.count",
                "a u32 exact reason count",
                descriptors.len(),
            )
        })?,
    )?;
    for (expected, descriptor) in NotObservedReasonV1::ALL.iter().zip(descriptors) {
        if descriptor.reason() != *expected {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.close.not_observed_reason_registry.reason",
                "the exact NotObserved reason rows in wire-code order",
                descriptor.reason().code(),
            ));
        }
        frame.push_u16(
            "coverage.close.not_observed_reason_registry.reason",
            descriptor.reason().code(),
        )?;
        frame.push_str(
            "coverage.close.not_observed_reason_registry.name",
            descriptor.name(),
        )?;
        frame.push_str(
            "coverage.close.not_observed_reason_registry.owner",
            descriptor.owner(),
        )?;
        frame.push_str(
            "coverage.close.not_observed_reason_registry.scope",
            descriptor.scope(),
        )?;
        frame.push_str(
            "coverage.close.not_observed_reason_registry.prerequisite",
            descriptor.prerequisite(),
        )?;
        frame.push_u16(
            "coverage.close.not_observed_reason_registry.diagnostic",
            descriptor.diagnostic().code(),
        )?;
        frame.push_str(
            "coverage.close.not_observed_reason_registry.no_claim",
            descriptor.no_claim(),
        )?;
    }
    Ok(
        BaseCoverageCloseNotObservedReasonRegistryRootV1::from_content_hash(
            frame.root(BaseCoverageCloseNotObservedReasonRegistryRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_deferred_reason_registry_root(
    descriptors: &[DeferredReasonDescriptorV1],
) -> Result<BaseCoverageCloseDeferredReasonRegistryRootV1, ConstructionErrorV2> {
    if descriptors.len() != DeferredReasonV1::ALL.len() {
        return Err(refusal(
            if descriptors.len() < DeferredReasonV1::ALL.len() {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "coverage.close.deferred_reason_registry.rows",
            "the exact one source-owned Deferred reason row",
            descriptors.len(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEDEFERREDREASONS\x01", 2 * 1024)?;
    frame.push_u32(
        "coverage.close.deferred_reason_registry.count",
        u32::try_from(descriptors.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.deferred_reason_registry.count",
                "a u32 exact reason count",
                descriptors.len(),
            )
        })?,
    )?;
    for (expected, descriptor) in DeferredReasonV1::ALL.iter().zip(descriptors) {
        if descriptor.reason() != *expected {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.close.deferred_reason_registry.reason",
                "the exact Deferred reason rows in wire-code order",
                descriptor.reason().code(),
            ));
        }
        frame.push_u16(
            "coverage.close.deferred_reason_registry.reason",
            descriptor.reason().code(),
        )?;
        frame.push_str(
            "coverage.close.deferred_reason_registry.name",
            descriptor.name(),
        )?;
        frame.push_str(
            "coverage.close.deferred_reason_registry.owner",
            descriptor.owner(),
        )?;
        frame.push_str(
            "coverage.close.deferred_reason_registry.scope",
            descriptor.scope(),
        )?;
        frame.push_str(
            "coverage.close.deferred_reason_registry.prerequisite",
            descriptor.prerequisite(),
        )?;
        frame.push_u16(
            "coverage.close.deferred_reason_registry.diagnostic",
            descriptor.diagnostic().code(),
        )?;
        frame.push_str(
            "coverage.close.deferred_reason_registry.no_claim",
            descriptor.no_claim(),
        )?;
    }
    Ok(
        BaseCoverageCloseDeferredReasonRegistryRootV1::from_content_hash(
            frame.root(BaseCoverageCloseDeferredReasonRegistryRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_capability_registry_root(
    rows: &[BaseCoverageCloseCapabilityDescriptorV1],
) -> Result<BaseCoverageCloseCapabilityRegistryRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECAPREG\x01", 8 * 1024)?;
    frame.push_u16(
        "coverage.close.capability_registry.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.capability_registry.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_registry.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_u16(
        "coverage.close.capability_registry.row_count",
        u16::try_from(rows.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.capability_registry.row_count",
                "a u16 bounded exact registry row count",
                rows.len(),
            )
        })?,
    )?;
    for row in rows {
        frame.push_u16("coverage.close.capability_registry.id", row.id.code())?;
        frame.push_bytes(
            "coverage.close.capability_registry.descriptor_root",
            row.root.content_hash().as_bytes(),
        )?;
    }
    Ok(
        BaseCoverageCloseCapabilityRegistryRootV1::from_content_hash(
            frame.root(BaseCoverageCloseCapabilityRegistryRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_capability_profile_registry_root(
    capability_registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    rows: &[BaseCoverageCloseCapabilityProfileDescriptorV1],
) -> Result<BaseCoverageCloseCapabilityProfileRegistryRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECAPPROFILES\x01", 8 * 1024)?;
    frame.push_u16(
        "coverage.close.capability_profile_registry.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.capability_profile_registry.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_profile_registry.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_bytes(
        "coverage.close.capability_profile_registry.capability_registry_root",
        capability_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_u16(
        "coverage.close.capability_profile_registry.row_count",
        u16::try_from(rows.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.capability_profile_registry.row_count",
                "a u16 bounded exact profile row count",
                rows.len(),
            )
        })?,
    )?;
    for row in rows {
        frame.push_u16(
            "coverage.close.capability_profile_registry.profile_code",
            row.profile.code(),
        )?;
        frame.push_str(
            "coverage.close.capability_profile_registry.profile_id",
            row.stable_id.as_str(),
        )?;
        push_close_capability_ids(
            &mut frame,
            "coverage.close.capability_profile_registry.required",
            &row.required,
        )?;
        push_close_capability_ids(
            &mut frame,
            "coverage.close.capability_profile_registry.permitted",
            &row.permitted,
        )?;
        frame.push_str(
            "coverage.close.capability_profile_registry.no_claim",
            row.no_claim.as_str(),
        )?;
    }
    Ok(
        BaseCoverageCloseCapabilityProfileRegistryRootV1::from_content_hash(
            frame.root(BaseCoverageCloseCapabilityProfileRegistryRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_capability_contract_root(
    capability_registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    profile_registry_root: BaseCoverageCloseCapabilityProfileRegistryRootV1,
    profile: BaseCoverageCloseCapabilityProfileV1,
    required: &[BaseCoverageCloseCapabilityIdV1],
    permitted: &[BaseCoverageCloseCapabilityIdV1],
    no_claim: &StableTokenV2,
) -> Result<BaseCoverageCloseCapabilityContractRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECAPCONTRACT\x01", 4 * 1024)?;
    frame.push_u16(
        "coverage.close.capability_contract.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.capability_contract.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_contract.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_bytes(
        "coverage.close.capability_contract.capability_registry_root",
        capability_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.capability_contract.profile_registry_root",
        profile_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_u16(
        "coverage.close.capability_contract.profile_code",
        profile.code(),
    )?;
    frame.push_str(
        "coverage.close.capability_contract.profile_id",
        profile.stable_id(),
    )?;
    push_close_capability_ids(
        &mut frame,
        "coverage.close.capability_contract.required",
        required,
    )?;
    push_close_capability_ids(
        &mut frame,
        "coverage.close.capability_contract.permitted",
        permitted,
    )?;
    frame.push_str(
        "coverage.close.capability_contract.no_claim",
        no_claim.as_str(),
    )?;
    Ok(
        BaseCoverageCloseCapabilityContractRootV1::from_content_hash(
            frame.root(BaseCoverageCloseCapabilityContractRootV1::DESCRIPTOR.domain()),
        ),
    )
}

#[allow(
    clippy::too_many_arguments,
    reason = "the five independently inspectable semantic sets are the AC56 contract"
)]
fn close_observed_capability_sets_root(
    capability_registry_root: BaseCoverageCloseCapabilityRegistryRootV1,
    capability_contract_root: BaseCoverageCloseCapabilityContractRootV1,
    required: &[BaseCoverageCloseCapabilityIdV1],
    granted: &[BaseCoverageCloseCapabilityIdV1],
    observed: &[BaseCoverageCloseCapabilityIdV1],
    returned: &[BaseCoverageCloseCapabilityIdV1],
    revoked: &[BaseCoverageCloseCapabilityIdV1],
) -> Result<BaseCoverageCloseObservedCapabilitySetsRootV1, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECAPOBSERVED\x01", 4 * 1024)?;
    frame.push_u16(
        "coverage.close.observed_capabilities.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "coverage.close.observed_capabilities.wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "coverage.close.observed_capabilities.predecessor_policy",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_bytes(
        "coverage.close.observed_capabilities.capability_registry_root",
        capability_registry_root.content_hash().as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.observed_capabilities.capability_contract_root",
        capability_contract_root.content_hash().as_bytes(),
    )?;
    for (partition, values) in [
        ("required", required),
        ("granted", granted),
        ("observed", observed),
        ("returned", returned),
        ("revoked", revoked),
    ] {
        frame.push_str("coverage.close.observed_capabilities.partition", partition)?;
        push_close_capability_ids(
            &mut frame,
            "coverage.close.observed_capabilities.ids",
            values,
        )?;
    }
    frame.push_str(
        "coverage.close.observed_capabilities.no_claim",
        BaseCoverageCloseObservedCapabilitySetsRootV1::DESCRIPTOR.no_claim(),
    )?;
    Ok(
        BaseCoverageCloseObservedCapabilitySetsRootV1::from_content_hash(
            frame.root(BaseCoverageCloseObservedCapabilitySetsRootV1::DESCRIPTOR.domain()),
        ),
    )
}

fn close_numeric_profile_root(
    partition: BaseCoverageCloseNumericPartitionV1,
    values: &[BaseCoverageCloseNumericExplicitV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(
        b"FSCLOSENUMERIC\x01",
        BASE_COVERAGE_CLOSE_NUMERIC_PROFILE_FRAME_MAX_BYTES_V1,
    )?;
    frame.push_u16(
        "coverage.close.numeric_profile.partition_code",
        partition.code(),
    )?;
    frame.push_str(
        "coverage.close.numeric_profile.partition_name",
        partition.name(),
    )?;
    frame.push_u32(
        "coverage.close.numeric_profile.count",
        u32::try_from(values.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.numeric_profile.count",
                "a u32 numeric explicit count",
                values.len(),
            )
        })?,
    )?;
    for value in values {
        frame.push_str(
            "coverage.close.numeric_profile.numeric_name",
            value.name.as_str(),
        )?;
        push_close_numeric_value(&mut frame, &value.value)?;
        push_close_numeric_unit(&mut frame, value.unit)?;
    }
    Ok(frame.root(BASE_COVERAGE_CLOSE_NUMERIC_PROFILE_DOMAIN_V1))
}

fn push_close_budget_value(
    frame: &mut CanonicalFrameV1,
    value: BaseCoverageCloseBudgetValueV1,
) -> Result<(), ConstructionErrorV2> {
    match value {
        BaseCoverageCloseBudgetValueV1::U32(value) => {
            frame.push_u32("coverage.close.five.budget_value_u32", value)?;
        }
        BaseCoverageCloseBudgetValueV1::U64(value) => {
            frame.push_u64("coverage.close.five.budget_value_u64", value)?;
        }
        BaseCoverageCloseBudgetValueV1::U128(value) => {
            frame.push_u128("coverage.close.five.budget_value_u128", value)?;
        }
    }
    Ok(())
}

fn push_close_budget_set(
    frame: &mut CanonicalFrameV1,
    budgets: BaseCoverageCloseBudgetSetV1,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(
        "coverage.close.five.budget_profile_code",
        budgets.profile.code(),
    )?;
    frame.push_str(
        "coverage.close.five.budget_profile_name",
        budgets.profile.name(),
    )?;
    frame.push_u16(
        "coverage.close.five.budget_axis_count",
        BASE_COVERAGE_CLOSE_BUDGET_AXIS_COUNT_V1 as u16,
    )?;
    for budget in budgets.rows {
        frame.push_u16("coverage.close.five.budget_axis_code", budget.axis.code())?;
        frame.push_str("coverage.close.five.budget_axis_name", budget.axis.name())?;
        frame.push_u16(
            "coverage.close.five.budget_width_tag",
            budget.axis.width().code(),
        )?;
        push_close_budget_value(frame, budget.hard)?;
        push_close_budget_value(frame, budget.soft)?;
        push_close_logical_unit_reference(frame, budget.unit)?;
    }
    Ok(())
}

#[allow(
    clippy::too_many_arguments,
    reason = "the root binds every Five Explicits collection separately"
)]
fn close_five_explicits_root(
    numeric_inputs_root: ContentHash,
    numeric_grants_root: ContentHash,
    numeric_observations_root: ContentHash,
    seed: &BaseCoverageCloseSeedExplicitV1,
    budgets: BaseCoverageCloseBudgetSetV1,
    versions: &BaseCoverageCloseVersionSetV1,
    capabilities: &BaseCoverageCloseCapabilitySetsV1,
    no_claim: &StableTokenV2,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(
        b"FSCLOSEFIVE\x01",
        BASE_COVERAGE_CLOSE_FIVE_EXPLICITS_FRAME_MAX_BYTES_V1,
    )?;
    frame.push_bytes(
        "coverage.close.five.numeric_inputs_root",
        numeric_inputs_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.five.numeric_grants_root",
        numeric_grants_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.five.numeric_observations_root",
        numeric_observations_root.as_bytes(),
    )?;

    match seed {
        BaseCoverageCloseSeedExplicitV1::Applicable {
            material,
            generator_version,
            minimizer_version,
        } => {
            frame.push_u16("coverage.close.five.seed_tag", 1)?;
            frame.push_bytes(
                "coverage.close.five.seed_material_root",
                material.root().content_hash().as_bytes(),
            )?;
            frame.push_str(
                "coverage.close.five.generator_version",
                generator_version.as_str(),
            )?;
            frame.push_str(
                "coverage.close.five.minimizer_version",
                minimizer_version.as_str(),
            )?;
        }
        BaseCoverageCloseSeedExplicitV1::Inapplicable { reason } => {
            frame.push_u16("coverage.close.five.seed_tag", 0)?;
            frame.push_u16("coverage.close.five.seed_reason_code", reason.code())?;
            frame.push_str("coverage.close.five.seed_reason_name", reason.name())?;
            frame.push_str("coverage.close.five.seed_reason_owner", reason.owner())?;
            frame.push_str("coverage.close.five.seed_reason_scope", reason.scope())?;
            frame.push_str(
                "coverage.close.five.seed_reason_prerequisite",
                reason.prerequisite(),
            )?;
            frame.push_str(
                "coverage.close.five.seed_reason_no_claim",
                reason.no_claim(),
            )?;
        }
    }
    push_close_budget_set(&mut frame, budgets)?;
    frame.push_u16(
        "coverage.close.five.api_generation",
        versions.api_generation.code(),
    )?;
    frame.push_u16(
        "coverage.close.five.wire_version",
        versions.wire_version.code(),
    )?;
    frame.push_bytes(
        "coverage.close.five.schema_root",
        versions.schema_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.five.source_root",
        versions.source_root.bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.five.build_root",
        versions.build_root.bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.five.toolchain_root",
        versions.toolchain_root.bytes(),
    )?;
    frame.push_u16("coverage.close.five.target", versions.target.code())?;
    frame.push_u16("coverage.close.five.profile", versions.profile.code())?;
    frame.push_bytes(
        "coverage.close.five.feature_set_root",
        versions.feature_set_root.as_bytes(),
    )?;
    for (partition, values) in [
        ("required", capabilities.required()),
        ("granted", capabilities.granted()),
        ("observed", capabilities.observed()),
        ("returned", capabilities.returned()),
        ("revoked", capabilities.revoked()),
    ] {
        frame.push_str("coverage.close.five.capability_partition", partition)?;
        frame.push_u32(
            "coverage.close.five.capability_count",
            u32::try_from(values.len()).map_err(|_| {
                refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "coverage.close.five.capability_count",
                    "a u32 capability count",
                    values.len(),
                )
            })?,
        )?;
        for value in values {
            frame.push_str("coverage.close.five.capability_id", value.as_str())?;
        }
    }
    frame.push_str("coverage.close.five.no_claim", no_claim.as_str())?;
    Ok(frame.root(BASE_COVERAGE_CLOSE_FIVE_EXPLICITS_DOMAIN_V1))
}

#[allow(
    clippy::too_many_arguments,
    reason = "the contribution root binds every AC53 downstream field"
)]
fn close_downstream_contribution_root(
    literal_expectation_oracle_root: ContentHash,
    semantic_input_root: ContentHash,
    budgets: BaseCoverageCloseContributionBudgetsV1,
    schema_root: ContentHash,
    log_schema_root: ContentHash,
    source_root: &SourceIdentityRootV2,
    build_root: &BuildIdentityRootV2,
    downstream_owner: &str,
    downstream_driver: &StableTokenV2,
    downstream_script: &str,
    downstream_manifest_path: &str,
    downstream_manifest_root: ContentHash,
    no_claim: &str,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECONTRIBUTION\x01", 4 * 1024)?;
    frame.push_bytes(
        "coverage.close.contribution.literal_expectation_oracle_root",
        literal_expectation_oracle_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.contribution.semantic_input_root",
        semantic_input_root.as_bytes(),
    )?;
    push_close_budget_set(&mut frame, budgets.resolved)?;
    frame.push_u32(
        "coverage.close.contribution.max_child_processes",
        budgets.max_child_processes,
    )?;
    frame.push_u32(
        "coverage.close.contribution.max_parallel_children",
        budgets.max_parallel_children,
    )?;
    frame.push_bytes(
        "coverage.close.contribution.schema_root",
        schema_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.contribution.log_schema_root",
        log_schema_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.contribution.source_root",
        source_root.bytes(),
    )?;
    frame.push_bytes("coverage.close.contribution.build_root", build_root.bytes())?;
    frame.push_str(
        "coverage.close.contribution.downstream_owner",
        downstream_owner,
    )?;
    frame.push_str(
        "coverage.close.contribution.downstream_driver",
        downstream_driver.as_str(),
    )?;
    frame.push_str(
        "coverage.close.contribution.downstream_script",
        downstream_script,
    )?;
    frame.push_str(
        "coverage.close.contribution.downstream_manifest_path",
        downstream_manifest_path,
    )?;
    frame.push_bytes(
        "coverage.close.contribution.downstream_manifest_root",
        downstream_manifest_root.as_bytes(),
    )?;
    frame.push_str("coverage.close.contribution.no_claim", no_claim)?;
    Ok(frame.root(BASE_COVERAGE_CLOSE_DOWNSTREAM_CONTRIBUTION_DOMAIN_V1))
}

fn close_reason_registry_root() -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEREASONS\x01", 16 * 1024)?;
    frame.push_u32(
        "coverage.close.reason_count",
        u32::try_from(BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.reason_count",
                "a u32 reason count",
                BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1.len(),
            )
        })?,
    )?;
    for descriptor in BASE_COVERAGE_CLOSE_REASON_DESCRIPTORS_V1 {
        frame.push_u16("coverage.close.reason_code", descriptor.code.code())?;
        frame.push_str("coverage.close.reason_name", descriptor.name)?;
        frame.push_str("coverage.close.reason_owner", descriptor.owner)?;
        frame.push_u16(
            "coverage.close.reason_execution_scope",
            descriptor.execution_scope.code(),
        )?;
        frame.push_str(
            "coverage.close.reason_prerequisite",
            descriptor.prerequisite,
        )?;
        frame.push_str("coverage.close.reason_no_claim", descriptor.no_claim)?;
    }
    Ok(frame.root(BASE_COVERAGE_CLOSE_REASON_REGISTRY_DOMAIN_V1))
}

fn close_applicability_evidence_root(
    reason_registry_root: ContentHash,
    reason: BaseCoverageCloseReasonCodeV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEAPPLICABILITY\x01", 512)?;
    frame.push_bytes(
        "coverage.close.applicability.reason_registry_root",
        reason_registry_root.as_bytes(),
    )?;
    frame.push_u16("coverage.close.applicability.reason", reason.code())?;
    Ok(frame.root(BASE_COVERAGE_CLOSE_APPLICABILITY_EVIDENCE_DOMAIN_V1))
}

fn close_cell_root(
    source_manifest_root: ContentHash,
    reason_registry_root: ContentHash,
    declaration: &BaseCoverageCloseCellDeclarationV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSECELL\x01", 4 * 1024)?;
    frame.push_bytes(
        "coverage.close.cell.source_manifest_root",
        source_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.cell.reason_registry_root",
        reason_registry_root.as_bytes(),
    )?;
    frame.push_u32(
        "coverage.close.cell.source_ordinal",
        declaration.source_ordinal,
    )?;
    frame.push_str(
        "coverage.close.cell.source_case_id",
        declaration.source_case_id(),
    )?;
    frame.push_u16(
        "coverage.close.cell.source_class",
        declaration.source_class.code(),
    )?;
    frame.push_str("coverage.close.cell.source_path", declaration.source_path())?;
    frame.push_u16("coverage.close.cell.group", declaration.group.code())?;
    frame.push_u16("coverage.close.cell.facet", declaration.facet.code())?;
    frame.push_u16(
        "coverage.close.cell.execution_scope",
        declaration.execution_scope.code(),
    )?;
    frame.push_u16(
        "coverage.close.cell.partition",
        declaration.partition.code(),
    )?;
    frame.push_u16(
        "coverage.close.cell.expected_decision",
        declaration.expected_decision.code(),
    )?;
    push_optional_close_reason(
        &mut frame,
        "coverage.close.cell.expected_reason_presence",
        "coverage.close.cell.expected_reason",
        declaration.expected_reason,
    )?;
    frame.push_u16(
        "coverage.close.cell.downstream_contribution_presence",
        u16::from(declaration.downstream_contribution.is_some()),
    )?;
    if let Some(contribution) = declaration.downstream_contribution.as_ref() {
        frame.push_bytes(
            "coverage.close.cell.downstream_contribution_root",
            contribution.root().as_bytes(),
        )?;
    }
    frame.push_bytes(
        "coverage.close.cell.five_explicits_root",
        declaration.five_explicits.root().as_bytes(),
    )?;
    Ok(frame.root(BASE_COVERAGE_CLOSE_CELL_DOMAIN_V1))
}

fn close_manifest_root(
    source_manifest_root: ContentHash,
    reason_registry_root: ContentHash,
    cells: &[BaseCoverageCloseManifestCellV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEMANIFEST\x01", 128 * 1024)?;
    frame.push_bytes(
        "coverage.close.manifest.source_manifest_root",
        source_manifest_root.as_bytes(),
    )?;
    frame.push_bytes(
        "coverage.close.manifest.reason_registry_root",
        reason_registry_root.as_bytes(),
    )?;
    frame.push_u32(
        "coverage.close.manifest.cell_count",
        u32::try_from(cells.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close.manifest.cell_count",
                "a u32 full-set cell count",
                cells.len(),
            )
        })?,
    )?;
    for cell in cells {
        frame.push_bytes("coverage.close.manifest.cell_root", cell.root().as_bytes())?;
    }
    Ok(frame.root(BASE_COVERAGE_CLOSE_MANIFEST_DOMAIN_V1))
}

fn push_optional_close_reason(
    frame: &mut CanonicalFrameV1,
    presence_field: &'static str,
    value_field: &'static str,
    reason: Option<BaseCoverageCloseReasonCodeV1>,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(presence_field, u16::from(reason.is_some()))?;
    if let Some(reason) = reason {
        frame.push_u16(value_field, reason.code())?;
    }
    Ok(())
}

fn push_optional_close_decision(
    frame: &mut CanonicalFrameV1,
    presence_field: &'static str,
    value_field: &'static str,
    decision: Option<BaseCoverageCloseDecisionV1>,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(presence_field, u16::from(decision.is_some()))?;
    if let Some(decision) = decision {
        frame.push_u16(value_field, decision.code())?;
    }
    Ok(())
}

#[allow(
    clippy::too_many_arguments,
    reason = "the presented-result root binds every AC53 result field"
)]
fn close_presented_result_root(
    close_manifest_root: ContentHash,
    cell_root: ContentHash,
    source_case_id: &str,
    group: BaseCoverageCloseGroupV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    status: BaseCoverageCloseResultStatusV1,
    observed_decision: Option<BaseCoverageCloseDecisionV1>,
    observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
    evidence: &BaseCoverageCloseResultEvidenceV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSCLOSERESULT\x01", 4 * 1024)?;
    frame.push_bytes(
        "coverage.close_result.close_manifest_root",
        close_manifest_root.as_bytes(),
    )?;
    frame.push_bytes("coverage.close_result.cell_root", cell_root.as_bytes())?;
    frame.push_str("coverage.close_result.source_case_id", source_case_id)?;
    frame.push_u16("coverage.close_result.group", group.code())?;
    frame.push_u16("coverage.close_result.facet", facet.code())?;
    frame.push_u16(
        "coverage.close_result.execution_scope",
        execution_scope.code(),
    )?;
    frame.push_u16("coverage.close_result.partition", partition.code())?;
    frame.push_u16(
        "coverage.close_result.expected_decision",
        expected_decision.code(),
    )?;
    push_optional_close_reason(
        &mut frame,
        "coverage.close_result.expected_reason_presence",
        "coverage.close_result.expected_reason",
        expected_reason,
    )?;
    frame.push_u16("coverage.close_result.status", status.code())?;
    push_optional_close_decision(
        &mut frame,
        "coverage.close_result.observed_decision_presence",
        "coverage.close_result.observed_decision",
        observed_decision,
    )?;
    push_optional_close_reason(
        &mut frame,
        "coverage.close_result.observed_reason_presence",
        "coverage.close_result.observed_reason",
        observed_reason,
    )?;
    frame.push_bytes(
        "coverage.close_result.evidence_root",
        evidence.root().as_bytes(),
    )?;
    frame.push_u16(
        "coverage.close_result.evidence_kind",
        evidence.kind().code(),
    )?;
    frame.push_u16(
        "coverage.close_result.retained_artifact_presence",
        u16::from(evidence.retained_artifact().is_some()),
    )?;
    if let Some(path) = evidence.retained_artifact() {
        frame.push_str("coverage.close_result.retained_artifact", path)?;
    }
    Ok(frame.root(BASE_COVERAGE_CLOSE_PRESENTED_RESULT_DOMAIN_V1))
}

fn close_report_root(
    close_manifest_root: ContentHash,
    results: &[BaseCoverageClosePresentedResultV1],
    adversarial_eligible: u32,
    adversarial_matched: u32,
    first_divergence_id: Option<&str>,
    first_divergence_root: Option<ContentHash>,
) -> Result<ContentHash, ConstructionErrorV2> {
    if first_divergence_id.is_some() != first_divergence_root.is_some() {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.close_report.first_divergence",
            "an ID and presented-result root that are both present or both absent",
            first_divergence_id.is_some(),
        ));
    }
    let mut frame = CanonicalFrameV1::new(b"FSCLOSEREPORT\x01", 128 * 1024)?;
    frame.push_bytes(
        "coverage.close_report.close_manifest_root",
        close_manifest_root.as_bytes(),
    )?;
    frame.push_u32(
        "coverage.close_report.result_count",
        u32::try_from(results.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.close_report.result_count",
                "a u32 full-set result count",
                results.len(),
            )
        })?,
    )?;
    for result in results {
        frame.push_bytes(
            "coverage.close_report.result_root",
            result.root().as_bytes(),
        )?;
    }
    frame.push_u32(
        "coverage.close_report.adversarial_eligible",
        adversarial_eligible,
    )?;
    frame.push_u32(
        "coverage.close_report.adversarial_matched",
        adversarial_matched,
    )?;
    frame.push_u16(
        "coverage.close_report.first_divergence_presence",
        u16::from(first_divergence_id.is_some()),
    )?;
    if let (Some(id), Some(root)) = (first_divergence_id, first_divergence_root) {
        frame.push_str("coverage.close_report.first_divergence_id", id)?;
        frame.push_bytes(
            "coverage.close_report.first_divergence_root",
            root.as_bytes(),
        )?;
    }
    Ok(frame.root(BASE_COVERAGE_CLOSE_REPORT_DOMAIN_V1))
}

fn frozen_base_declarations() -> Result<Vec<BaseCoverageCaseDeclarationV1>, ConstructionErrorV2> {
    let capacity = BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1
        + BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1
        + BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1;
    let mut declarations = Vec::with_capacity(capacity);
    declarations.extend(frozen_rust_test_declarations()?);
    for case in COMPILE_FAIL_TEMPLATES_V1 {
        if !is_exact_compiler_error_code_v1(case.expected_error_code) {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.compile_fail.expected_error_code",
                "one exact E followed by four ASCII digits",
                case.expected_error_code,
            ));
        }
        let is_root_free_missing_type_contract = case.module == "identity"
            && case.case_name == "no-standalone-root-for-root-free-evaluator-members";
        if (case.expected_error_code == "E0432") != is_root_free_missing_type_contract {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.compile_fail.unresolved_import_owner",
                "E0432 owned only by the root-free evaluator-member absence contract",
                case.case_name,
            ));
        }
        declarations.push(BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::CompileFailDoctest,
            format!("compile-fail:{}:{}", case.module, case.case_name),
            case.source_path,
        )?);
    }
    for test_name in MANIFEST_CONTRACT_TEST_NAMES_V1 {
        declarations.push(BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::ManifestContract,
            format!("manifest-contract:coverage:{test_name}"),
            "crates/fs-evidence-runner/src/coverage.rs",
        )?);
    }
    if declarations.len() != capacity {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.base.case_count",
            "the exact 217 classified Rust tests, 78 compile-fail cases, and 29 manifest-contract cases",
            declarations.len(),
        ));
    }
    validate_unique_declarations(&declarations)?;
    Ok(declarations)
}

fn frozen_rust_test_declarations() -> Result<Vec<BaseCoverageCaseDeclarationV1>, ConstructionErrorV2>
{
    let mut declarations = Vec::with_capacity(BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1);
    let mut class_counts = BTreeMap::<BaseCoverageManifestClassV1, usize>::new();
    let mut source_tests = BTreeSet::new();
    let mut previous_module = None;

    for module in RUST_TEST_MODULE_TEMPLATES_V1 {
        if previous_module.is_some_and(|previous| module.module <= previous) {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.rust_test.module",
                "strict source-module order",
                module.module,
            ));
        }
        previous_module = Some(module.module);
        if module.class_templates.is_empty() {
            return Err(refusal(
                ConstructionErrorKindV2::Missing,
                "coverage.rust_test.class_templates",
                "at least one explicit evidence-class group per source module",
                module.module,
            ));
        }

        let mut previous_class = None;
        for class_template in module.class_templates {
            if !BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES
                .contains(&class_template.class)
            {
                return Err(refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "coverage.rust_test.class",
                    "one of the six frozen Rust-test evidence classes",
                    class_template.class.code(),
                ));
            }
            if previous_class.is_some_and(|previous| class_template.class.code() <= previous) {
                return Err(refusal(
                    ConstructionErrorKindV2::OutOfOrder,
                    "coverage.rust_test.class",
                    "strict nonrepeating evidence-class order within each source module",
                    class_template.class.code(),
                ));
            }
            previous_class = Some(class_template.class.code());
            if class_template.test_names.is_empty() {
                return Err(refusal(
                    ConstructionErrorKindV2::Missing,
                    "coverage.rust_test.test_names",
                    "at least one explicitly classified test name",
                    module.module,
                ));
            }

            for test_name in class_template.test_names {
                if !source_tests.insert((module.module, *test_name)) {
                    return Err(refusal(
                        ConstructionErrorKindV2::Duplicate,
                        "coverage.rust_test.source_identity",
                        "one explicit evidence-class assignment per source test",
                        format_args!("{}:{test_name}", module.module),
                    ));
                }
                *class_counts.entry(class_template.class).or_default() += 1;
                declarations.push(BaseCoverageCaseDeclarationV1::new(
                    class_template.class,
                    format!(
                        "{}{}:{test_name}",
                        class_template.class.stable_prefix(),
                        module.module
                    ),
                    module.source_path,
                )?);
            }
        }
    }

    for (class, expected_count) in BASE_COVERAGE_RUST_TEST_CLASS_COUNTS_V1 {
        let observed_count = class_counts.get(&class).copied().unwrap_or_default();
        if observed_count != expected_count {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.rust_test.class_count",
                "the exact authoritative count for this Rust-test evidence class",
                format_args!("{}:{observed_count}", class.code()),
            ));
        }
    }
    if class_counts.len() != BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES.len()
        || declarations.len() != BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.rust_test.total",
            "exactly 217 tests partitioned once across all six required evidence classes",
            declarations.len(),
        ));
    }
    Ok(declarations)
}

fn manifest_from_declarations(
    declarations: Vec<BaseCoverageCaseDeclarationV1>,
) -> Result<BaseCoverageManifestV1, ConstructionErrorV2> {
    validate_unique_declarations(&declarations)?;
    let mut cases = Vec::with_capacity(declarations.len());
    for (index, declaration) in declarations.into_iter().enumerate() {
        let ordinal = u32::try_from(index + 1).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.manifest.ordinal",
                "a one-based u32 ordinal",
                index + 1,
            )
        })?;
        cases.push(BaseCoverageManifestCaseV1 {
            ordinal,
            declaration,
        });
    }
    let root = manifest_root(&cases)?;
    Ok(BaseCoverageManifestV1 {
        cases: cases.into_boxed_slice(),
        root,
    })
}

fn validate_exact_declaration_sequence(
    expected: &[BaseCoverageCaseDeclarationV1],
    presented: &[BaseCoverageCaseDeclarationV1],
) -> Result<(), ConstructionErrorV2> {
    validate_unique_declarations(presented)?;
    if presented.len() < expected.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.base.declarations",
            "the complete exact frozen declaration sequence",
            presented.len(),
        ));
    }
    if presented.len() > expected.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.base.declarations",
            "no declaration beyond the exact frozen sequence",
            presented.len(),
        ));
    }
    for (index, (observed, exact)) in presented.iter().zip(expected).enumerate() {
        if observed == exact {
            continue;
        }
        let kind = if expected.contains(observed) {
            ConstructionErrorKindV2::OutOfOrder
        } else {
            ConstructionErrorKindV2::Incompatible
        };
        return Err(refusal(
            kind,
            "coverage.base.declaration",
            "the exact declaration at this ordinal",
            format_args!("{index}:{}", observed.id()),
        ));
    }
    Ok(())
}

fn validate_extensions(
    extensions: &[BaseCoverageCaseDeclarationV1],
) -> Result<(), ConstructionErrorV2> {
    if extensions.len() > BASE_COVERAGE_EXTENSION_CASES_MAX_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            "coverage.extensions",
            "at most 1024 exact extension declarations",
            extensions.len(),
        ));
    }
    let base_ids = frozen_base_declarations()?
        .into_iter()
        .map(|case| case.id)
        .collect::<BTreeSet<_>>();
    let mut seen = BTreeSet::new();
    let mut previous = None;
    for extension in extensions {
        if !extension.class.is_extension() {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.extension.class",
                "one of the six exact extension classes",
                extension.class.code(),
            ));
        }
        if base_ids.contains(extension.id()) || !seen.insert(extension.id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.extension.source_case_id",
                "a globally unique extension source-case ID",
                extension.id(),
            ));
        }
        let key = (
            extension.class.code(),
            extension.id(),
            extension.source_path(),
        );
        if previous.is_some_and(|prior| key <= prior) {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.extensions",
                "strict (class, ID, source path) order",
                extension.id(),
            ));
        }
        previous = Some(key);
    }
    Ok(())
}

fn validate_unique_declarations(
    declarations: &[BaseCoverageCaseDeclarationV1],
) -> Result<(), ConstructionErrorV2> {
    let mut ids = BTreeSet::new();
    for declaration in declarations {
        if !ids.insert(declaration.id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.manifest.source_case_id",
                "one globally unique declaration per source-case ID",
                declaration.id(),
            ));
        }
    }
    Ok(())
}

fn validate_selection_against_manifest(
    manifest: &BaseCoverageManifestV1,
    selection: &BaseCoverageExecutableSubsetV1,
) -> Result<(), ConstructionErrorV2> {
    if selection.manifest_root != manifest.root {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.selection.manifest_root",
            "the current exact manifest root",
            selection.manifest_root.to_hex(),
        ));
    }
    let selected = selection
        .source_case_ids
        .iter()
        .map(Box::as_ref)
        .collect::<Vec<_>>();
    let reconstructed = manifest.select_exact(&selected)?;
    if reconstructed.root != selection.root {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.selection.root",
            "the root reconstructed from the exact manifest and selected IDs",
            selection.root.to_hex(),
        ));
    }
    Ok(())
}

fn validate_presented_results(
    manifest: &BaseCoverageManifestV1,
    selection: &BaseCoverageExecutableSubsetV1,
    presented: &[BaseCoveragePresentedResultV1],
) -> Result<(), ConstructionErrorV2> {
    let manifest_ids = manifest
        .cases
        .iter()
        .map(BaseCoverageManifestCaseV1::id)
        .collect::<BTreeSet<_>>();
    let selected_ids = selection
        .source_case_ids
        .iter()
        .map(Box::as_ref)
        .collect::<BTreeSet<_>>();
    let mut reported = BTreeSet::new();
    for result in presented {
        if result.manifest_root != manifest.root {
            return Err(refusal(
                ConstructionErrorKindV2::Incompatible,
                "coverage.result.manifest_root",
                "the current exact manifest root",
                result.manifest_root.to_hex(),
            ));
        }
        if !manifest_ids.contains(result.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::UnknownCode,
                "coverage.result.source_case_id",
                "an ID mapped by the exact manifest",
                result.source_case_id(),
            ));
        }
        if !reported.insert(result.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Duplicate,
                "coverage.result.source_case_id",
                "exactly one presented result per selected source-case ID",
                result.source_case_id(),
            ));
        }
        if !selected_ids.contains(result.source_case_id()) {
            return Err(refusal(
                ConstructionErrorKindV2::Unexpected,
                "coverage.result.source_case_id",
                "an ID in the caller-selected executable subset",
                result.source_case_id(),
            ));
        }
    }
    if presented.len() < selection.source_case_ids.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.results",
            "one result for every selected source-case ID",
            presented.len(),
        ));
    }
    if presented.len() > selection.source_case_ids.len() {
        return Err(refusal(
            ConstructionErrorKindV2::Unexpected,
            "coverage.results",
            "no result beyond the selected source-case IDs",
            presented.len(),
        ));
    }
    for (index, (result, selected_id)) in presented
        .iter()
        .zip(selection.source_case_ids.iter())
        .enumerate()
    {
        if result.source_case_id() != selected_id.as_ref() {
            return Err(refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "coverage.results",
                "exact selected manifest order",
                format_args!("{index}:{}", result.source_case_id()),
            ));
        }
    }
    Ok(())
}

fn validate_case_id(
    class: BaseCoverageManifestClassV1,
    id: &str,
) -> Result<(), ConstructionErrorV2> {
    validate_untyped_case_id(id)?;
    if !id.starts_with(class.stable_prefix()) {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.source_case_id",
            "the stable prefix owned by its closed coverage class",
            id,
        ));
    }
    Ok(())
}

fn validate_untyped_case_id(id: &str) -> Result<(), ConstructionErrorV2> {
    if id.is_empty() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.source_case_id",
            "a nonempty stable source-case ID",
            0,
        ));
    }
    if id.len() > BASE_COVERAGE_CASE_ID_MAX_BYTES_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            "coverage.source_case_id",
            "at most 160 UTF-8 bytes",
            id.len(),
        ));
    }
    if !id.bytes().all(|byte| {
        byte.is_ascii_lowercase()
            || byte.is_ascii_digit()
            || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
    }) || id.starts_with(':')
        || id.ends_with(':')
        || id.contains("::")
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.source_case_id",
            "lowercase ASCII stable-ID grammar without empty colon segments",
            id,
        ));
    }
    Ok(())
}

fn validate_source_path(path: &str) -> Result<(), ConstructionErrorV2> {
    if path.is_empty() {
        return Err(refusal(
            ConstructionErrorKindV2::Missing,
            "coverage.source_path",
            "a nonempty workspace-relative source path",
            0,
        ));
    }
    if path.len() > BASE_COVERAGE_SOURCE_PATH_MAX_BYTES_V1 {
        return Err(refusal(
            ConstructionErrorKindV2::TooLarge,
            "coverage.source_path",
            "at most 240 UTF-8 bytes",
            path.len(),
        ));
    }
    if path.starts_with('/')
        || path.ends_with('/')
        || path.contains('\\')
        || path.chars().any(char::is_control)
        || path
            .split('/')
            .any(|segment| segment.is_empty() || segment == "." || segment == "..")
    {
        return Err(refusal(
            ConstructionErrorKindV2::Incompatible,
            "coverage.source_path",
            "an exact clean workspace-relative UTF-8 path",
            path,
        ));
    }
    Ok(())
}

fn manifest_root(cases: &[BaseCoverageManifestCaseV1]) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASECOVERAGEMANIFEST\x01", 512 * 1024)?;
    frame.push_u32(
        "coverage.manifest.case_count",
        u32::try_from(cases.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.manifest.case_count",
                "a u32 case count",
                cases.len(),
            )
        })?,
    )?;
    for case in cases {
        frame.push_u32("coverage.manifest.ordinal", case.ordinal)?;
        frame.push_u16("coverage.manifest.class", case.class().code())?;
        frame.push_str("coverage.manifest.source_case_id", case.id())?;
        frame.push_str("coverage.manifest.source_path", case.source_path())?;
    }
    Ok(frame.root(BASE_COVERAGE_MANIFEST_DOMAIN_V1))
}

fn selection_root(
    manifest_root: ContentHash,
    source_case_ids: &[Box<str>],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASECOVERAGESELECTION\x01", 192 * 1024)?;
    frame.push_bytes("coverage.selection.manifest_root", manifest_root.as_bytes())?;
    frame.push_u32(
        "coverage.selection.case_count",
        u32::try_from(source_case_ids.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.selection.case_count",
                "a u32 selected case count",
                source_case_ids.len(),
            )
        })?,
    )?;
    for id in source_case_ids {
        frame.push_str("coverage.selection.source_case_id", id)?;
    }
    Ok(frame.root(BASE_COVERAGE_SELECTION_DOMAIN_V1))
}

fn presented_result_root(
    manifest_root: ContentHash,
    source_case_id: &str,
    outcome: BaseCoveragePresentedOutcomeV1,
    evidence_root: ContentHash,
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASECOVERAGERESULT\x01", 512)?;
    frame.push_bytes("coverage.result.manifest_root", manifest_root.as_bytes())?;
    frame.push_str("coverage.result.source_case_id", source_case_id)?;
    frame.push_u16("coverage.result.outcome", outcome.code())?;
    frame.push_bytes("coverage.result.evidence_root", evidence_root.as_bytes())?;
    Ok(frame.root(BASE_COVERAGE_PRESENTED_RESULT_DOMAIN_V1))
}

fn checked_report_root(
    selection_root: ContentHash,
    results: &[BaseCoveragePresentedResultV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASECOVERAGEREPORT\x01", 128 * 1024)?;
    frame.push_bytes("coverage.report.selection_root", selection_root.as_bytes())?;
    frame.push_u32(
        "coverage.report.result_count",
        u32::try_from(results.len()).map_err(|_| {
            refusal(
                ConstructionErrorKindV2::TooLarge,
                "coverage.report.result_count",
                "a u32 result count",
                results.len(),
            )
        })?,
    )?;
    for result in results {
        frame.push_bytes("coverage.report.result_root", result.root.as_bytes())?;
    }
    Ok(frame.root(BASE_COVERAGE_CHECKED_REPORT_DOMAIN_V1))
}

fn refusal(
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
    observed: impl std::fmt::Display,
) -> ConstructionErrorV2 {
    // This helper is used at structural boundaries where many observed values
    // originate in caller-presented IDs, paths, or reconstructed rows.  The
    // non-wire construction refusal therefore records only provenance; exact
    // typed expected/observed data belongs in the actionable diagnostic and
    // close-log projections.
    let _ = observed;
    ConstructionErrorV2::new_redacted(
        kind,
        field,
        expected,
        ConstructionObservedDataClassV2::CallerControlledText,
    )
}

#[cfg(test)]
mod tests {
    #[allow(
        deprecated,
        reason = "one compatibility assertion freezes the misleading historical aggregate alias"
    )]
    use super::{
        BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
        BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
        BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1,
        BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1, BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1,
        BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1,
        BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1, BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1,
        BASE_COVERAGE_RUST_TEST_CLASS_COUNTS_V1, BASE_COVERAGE_UNIT_CASE_COUNT_V1,
        BaseCoverageCaseDeclarationV1, BaseCoverageCheckedReportV1, BaseCoverageCloseBudgetAxisV1,
        BaseCoverageCloseBudgetProfileV1, BaseCoverageCloseBudgetSetV1,
        BaseCoverageCloseBudgetValueV1, BaseCoverageCloseBudgetWidthV1,
        BaseCoverageCloseCapabilityContractV1, BaseCoverageCloseCapabilityDescriptorV1,
        BaseCoverageCloseCapabilityIdV1, BaseCoverageCloseCapabilityPolicyV1,
        BaseCoverageCloseCapabilityProfileRegistryV1, BaseCoverageCloseCapabilityProfileV1,
        BaseCoverageCloseCapabilityRegistryV1, BaseCoverageCloseCellDeclarationV1,
        BaseCoverageCloseContributionBudgetsV1, BaseCoverageCloseDecisionV1,
        BaseCoverageCloseDownstreamContributionV1, BaseCoverageCloseEvidenceKindV1,
        BaseCoverageCloseExecutionScopeV1, BaseCoverageCloseFacetV1,
        BaseCoverageCloseFiveExplicitsV1, BaseCoverageCloseGroupV1,
        BaseCoverageCloseLogicalUnitReferenceV1, BaseCoverageCloseManifestCellV1,
        BaseCoverageCloseManifestV1, BaseCoverageCloseNumericExplicitV1,
        BaseCoverageCloseNumericPartitionV1, BaseCoverageCloseNumericUnitV1,
        BaseCoverageCloseObservedCapabilitySetsV1, BaseCoverageClosePartitionV1,
        BaseCoverageClosePresentedResultV1, BaseCoverageCloseReasonCodeV1,
        BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1,
        BaseCoverageCloseRegisteredExtensionCapabilityIdV1,
        BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1,
        BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1,
        BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1,
        BaseCoverageCloseRegisteredExtensionCapabilityScopeV1,
        BaseCoverageCloseRegisteredExtensionCapabilitySetV1,
        BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1, BaseCoverageCloseReportV1,
        BaseCoverageCloseResultEvidenceV1, BaseCoverageCloseTypedBudgetV1,
        BaseCoverageManifestClassV1, BaseCoverageManifestV1, BaseCoveragePresentedOutcomeV1,
        BaseCoveragePresentedResultV1, COMPILE_FAIL_TEMPLATES_V1,
        CanonicalSchemaImpactDispositionV1, CanonicalSchemaMigrationPolicyV1,
        DeferredReasonRegistryV1, DeferredReasonV1, NotObservedReasonRegistryV1,
        NotObservedReasonV1, RuntimeObservationDispositionV1,
        base_coverage_close_capability_profile_for_source_case_v1,
        base_coverage_close_nominal_root_descriptors_v1, is_exact_compiler_error_code_v1,
    };
    use crate::{
        ConstructionErrorKindV2,
        catalog::{DiagnosticCodeV2, LogicalUnitV2, WirePredecessorPolicyV1},
        value::{
            DecimalV2, F32BitsV2, F64BitsV2, NumericValueV2, RationalV2, StableTokenV2, UnitV2,
        },
    };
    use fs_blake3::{ContentHash, hash_domain};
    use std::collections::{BTreeMap, BTreeSet};

    const COMPILE_FAIL_RUSTDOC_SOURCES_V1: [(&str, &str); 12] = [
        (
            "crates/fs-evidence-runner/src/budget.rs",
            include_str!("budget.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/capability.rs",
            include_str!("capability.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/command.rs",
            include_str!("command.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/diagnostic.rs",
            include_str!("diagnostic.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/extension.rs",
            include_str!("extension.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/identity.rs",
            include_str!("identity.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/limits.rs",
            include_str!("limits.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/path.rs",
            include_str!("path.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/projection.rs",
            include_str!("projection.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/schema_impact.rs",
            include_str!("schema_impact.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/state.rs",
            include_str!("state.rs"),
        ),
        (
            "crates/fs-evidence-runner/src/value.rs",
            include_str!("value.rs"),
        ),
    ];

    const COMPILE_FAIL_ERROR_CODE_DISTRIBUTION_V1: [(&str, usize); 10] = [
        ("E0277", 9),
        ("E0308", 24),
        ("E0369", 2),
        ("E0423", 2),
        ("E0432", 1),
        ("E0451", 7),
        ("E0599", 5),
        ("E0609", 3),
        ("E0616", 20),
        ("E0624", 5),
    ];

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct CompileFailFenceV1 {
        source_path: &'static str,
        line: usize,
        expected_error_code: Box<str>,
    }

    fn compile_fail_fences_v1(
        sources: &[(&'static str, &'static str)],
        expected_count: usize,
    ) -> Result<Vec<CompileFailFenceV1>, String> {
        let mut fences = Vec::new();
        for &(source_path, source) in sources {
            for (index, line) in source.lines().enumerate() {
                let Some(doc_line) = line.trim_start().strip_prefix("///") else {
                    continue;
                };
                let header = doc_line.trim_start();
                if !header.starts_with("```compile_fail") {
                    continue;
                }
                if header == "```compile_fail" {
                    return Err(format!(
                        "bare compile_fail fence at {source_path}:{}",
                        index + 1
                    ));
                }
                let Some(error_code) = header.strip_prefix("```compile_fail,") else {
                    return Err(format!(
                        "malformed compile_fail fence at {source_path}:{}: {header}",
                        index + 1
                    ));
                };
                if !is_exact_compiler_error_code_v1(error_code) {
                    return Err(format!(
                        "malformed or multi-code compile_fail fence at {source_path}:{}: {header}",
                        index + 1
                    ));
                }
                fences.push(CompileFailFenceV1 {
                    source_path,
                    line: index + 1,
                    expected_error_code: error_code.into(),
                });
            }
        }
        if fences.len() != expected_count {
            return Err(format!(
                "compile_fail fence count mismatch: expected {expected_count}, observed {}",
                fences.len()
            ));
        }
        Ok(fences)
    }

    fn compile_fail_distribution_v1<'a>(
        codes: impl IntoIterator<Item = &'a str>,
    ) -> Vec<(&'a str, usize)> {
        let mut distribution = BTreeMap::new();
        for code in codes {
            *distribution.entry(code).or_insert(0) += 1;
        }
        distribution.into_iter().collect()
    }

    fn root(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.fs-evidence-runner.coverage-test.v1",
            label.as_bytes(),
        )
    }

    fn fixed_logical(unit: LogicalUnitV2) -> BaseCoverageCloseLogicalUnitReferenceV1 {
        BaseCoverageCloseLogicalUnitReferenceV1::fixed(unit).expect("fixed logical unit")
    }

    fn numeric(
        name: impl Into<String>,
        value: NumericValueV2,
        unit: BaseCoverageCloseNumericUnitV1,
    ) -> BaseCoverageCloseNumericExplicitV1 {
        BaseCoverageCloseNumericExplicitV1::new(
            StableTokenV2::new(name).expect("stable numeric name"),
            value,
            unit,
        )
    }

    fn extension_stable_id(
        value: impl Into<String>,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1 {
        BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1::new(
            StableTokenV2::new(value).expect("stable extension capability ID"),
        )
        .expect("namespaced non-base extension capability ID")
    }

    fn extension_owner(
        value: impl Into<String>,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1 {
        BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1::new(
            StableTokenV2::new(value).expect("stable extension owner"),
        )
        .expect("namespaced extension owner")
    }

    fn extension_scope(
        value: impl Into<String>,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityScopeV1 {
        BaseCoverageCloseRegisteredExtensionCapabilityScopeV1::new(
            StableTokenV2::new(value).expect("stable extension scope"),
        )
        .expect("namespaced extension scope")
    }

    fn extension_no_claim(
        value: impl Into<String>,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1 {
        BaseCoverageCloseRegisteredExtensionCapabilityNoClaimV1::new(
            StableTokenV2::new(value).expect("stable extension no-claim"),
        )
    }

    fn extension_descriptor(
        code: u16,
    ) -> BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1 {
        BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(code)
                .expect("bounded extension capability ID"),
            extension_stable_id(format!(
                "org.example.fs-evidence-runner.extension.capability-{code}"
            )),
            extension_owner("org.example.fs-evidence-runner.extension-owner"),
            extension_scope("org.example.fs-evidence-runner.extension-scope"),
            extension_no_claim("extension-contract-proves-no-acquisition-effect-or-authority"),
        )
        .expect("valid extension descriptor")
    }

    fn five_from_template(
        numeric_inputs: Vec<BaseCoverageCloseNumericExplicitV1>,
        numeric_observations: Vec<BaseCoverageCloseNumericExplicitV1>,
        budgets: BaseCoverageCloseBudgetSetV1,
    ) -> Result<BaseCoverageCloseFiveExplicitsV1, crate::ConstructionErrorV2> {
        five_from_template_with_grants(numeric_inputs, vec![], numeric_observations, budgets)
    }

    fn five_from_template_with_grants(
        numeric_inputs: Vec<BaseCoverageCloseNumericExplicitV1>,
        numeric_grants: Vec<BaseCoverageCloseNumericExplicitV1>,
        numeric_observations: Vec<BaseCoverageCloseNumericExplicitV1>,
        budgets: BaseCoverageCloseBudgetSetV1,
    ) -> Result<BaseCoverageCloseFiveExplicitsV1, crate::ConstructionErrorV2> {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("template manifest");
        let template = manifest.cells()[0].five_explicits();
        BaseCoverageCloseFiveExplicitsV1::new(
            numeric_inputs,
            numeric_grants,
            numeric_observations,
            template.seed().clone(),
            budgets,
            template.versions().clone(),
            template.capabilities().clone(),
            template.no_claim().clone(),
        )
    }

    fn budget_row(
        axis: BaseCoverageCloseBudgetAxisV1,
        hard: BaseCoverageCloseBudgetValueV1,
        soft: BaseCoverageCloseBudgetValueV1,
        unit: LogicalUnitV2,
    ) -> Result<BaseCoverageCloseTypedBudgetV1, crate::ConstructionErrorV2> {
        BaseCoverageCloseTypedBudgetV1::new(axis, hard, soft, fixed_logical(unit))
    }

    fn one_over_budget_value(
        value: BaseCoverageCloseBudgetValueV1,
    ) -> Option<BaseCoverageCloseBudgetValueV1> {
        match value {
            BaseCoverageCloseBudgetValueV1::U32(value) => value
                .checked_add(1)
                .map(BaseCoverageCloseBudgetValueV1::U32),
            BaseCoverageCloseBudgetValueV1::U64(value) => value
                .checked_add(1)
                .map(BaseCoverageCloseBudgetValueV1::U64),
            BaseCoverageCloseBudgetValueV1::U128(value) => value
                .checked_add(1)
                .map(BaseCoverageCloseBudgetValueV1::U128),
        }
    }

    fn extension(
        class: BaseCoverageManifestClassV1,
        id: &str,
        path: &str,
    ) -> BaseCoverageCaseDeclarationV1 {
        BaseCoverageCaseDeclarationV1::new(class, id, path).expect("valid extension fixture")
    }

    fn result(
        manifest: &BaseCoverageManifestV1,
        id: &str,
        outcome: BaseCoveragePresentedOutcomeV1,
        evidence: &str,
    ) -> BaseCoveragePresentedResultV1 {
        BaseCoveragePresentedResultV1::new(manifest.root(), id, outcome, root(evidence))
            .expect("valid result fixture")
    }

    fn close_evidence(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        label: &str,
    ) -> BaseCoverageCloseResultEvidenceV1 {
        match cell.execution_scope() {
            BaseCoverageCloseExecutionScopeV1::CrateTest
            | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => {
                BaseCoverageCloseResultEvidenceV1::owned_harness_execution(root(label), None)
                    .expect("owned evidence")
            }
            BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
                BaseCoverageCloseResultEvidenceV1::in_process_projection_execution(
                    root(label),
                    None,
                )
                .expect("projection evidence")
            }
            BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
                BaseCoverageCloseResultEvidenceV1::immutable_downstream_contribution(
                    cell.downstream_contribution().expect("contribution"),
                )
                .expect("immutable contribution")
            }
            BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
                BaseCoverageCloseResultEvidenceV1::applicability_declaration(
                    manifest,
                    cell.expected_reason().expect("applicability reason"),
                )
                .expect("applicability evidence")
            }
        }
    }

    fn matched_close_results(
        manifest: &BaseCoverageCloseManifestV1,
    ) -> Vec<BaseCoverageClosePresentedResultV1> {
        manifest
            .cells()
            .iter()
            .enumerate()
            .map(|(index, cell)| {
                BaseCoverageClosePresentedResultV1::matched(
                    manifest,
                    cell,
                    close_evidence(manifest, cell, &format!("matched-{index}")),
                )
                .expect("matched result")
            })
            .collect()
    }

    const fn divergent_decision(
        expected: BaseCoverageCloseDecisionV1,
    ) -> BaseCoverageCloseDecisionV1 {
        match expected {
            BaseCoverageCloseDecisionV1::Accept => BaseCoverageCloseDecisionV1::Refuse,
            BaseCoverageCloseDecisionV1::Refuse => BaseCoverageCloseDecisionV1::Accept,
            BaseCoverageCloseDecisionV1::Fail => BaseCoverageCloseDecisionV1::Accept,
            BaseCoverageCloseDecisionV1::Unsupported => BaseCoverageCloseDecisionV1::Refuse,
            BaseCoverageCloseDecisionV1::Inapplicable => BaseCoverageCloseDecisionV1::Refuse,
        }
    }

    #[test]
    #[allow(
        deprecated,
        reason = "this test intentionally verifies the compatibility-only aggregate alias"
    )]
    fn frozen_manifest_enumerates_exact_ratified_inventory_and_external_classes() {
        let first = BaseCoverageManifestV1::frozen().expect("frozen manifest");
        let second = BaseCoverageManifestV1::frozen().expect("deterministic manifest");
        assert_eq!(first, second);
        assert_eq!(BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1, 130);
        assert_eq!(BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1, 87);
        assert_eq!(BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1, 217);
        assert_eq!(
            BASE_COVERAGE_UNIT_CASE_COUNT_V1,
            BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1
        );
        let mut classified_total = 0;
        for (class, exact_count) in BASE_COVERAGE_RUST_TEST_CLASS_COUNTS_V1 {
            assert!(BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES.contains(&class));
            assert_ne!(
                exact_count, 0,
                "{class:?} must remain independently nonzero"
            );
            assert_eq!(first.case_count(class), exact_count, "{class:?}");
            classified_total += exact_count;
        }
        assert_eq!(classified_total, BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1);
        assert_eq!(
            first.case_count(BaseCoverageManifestClassV1::CompileFailDoctest),
            BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1
        );
        assert_eq!(BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1, 78);
        let compile_fail_fences = compile_fail_fences_v1(
            &COMPILE_FAIL_RUSTDOC_SOURCES_V1,
            BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1,
        )
        .expect("all 78 Rustdoc compile-fail fences are singly cause-gated");
        assert_eq!(
            compile_fail_distribution_v1(
                compile_fail_fences
                    .iter()
                    .map(|fence| fence.expected_error_code.as_ref()),
            ),
            COMPILE_FAIL_ERROR_CODE_DISTRIBUTION_V1.to_vec(),
            "live Rustdoc error-code distribution"
        );
        assert_eq!(
            compile_fail_distribution_v1(
                COMPILE_FAIL_TEMPLATES_V1
                    .iter()
                    .map(|case| case.expected_error_code),
            ),
            COMPILE_FAIL_ERROR_CODE_DISTRIBUTION_V1.to_vec(),
            "source-authoritative compile-fail oracle distribution"
        );
        let mut live_source_counts = BTreeMap::new();
        for fence in &compile_fail_fences {
            *live_source_counts.entry(fence.source_path).or_insert(0) += 1;
        }
        let mut oracle_source_counts = BTreeMap::new();
        for case in COMPILE_FAIL_TEMPLATES_V1 {
            assert!(
                is_exact_compiler_error_code_v1(case.expected_error_code),
                "{}:{} has malformed expected error code {}",
                case.module,
                case.case_name,
                case.expected_error_code
            );
            *oracle_source_counts.entry(case.source_path).or_insert(0) += 1;
        }
        assert_eq!(
            live_source_counts, oracle_source_counts,
            "every source-owned Rustdoc fence has exactly one oracle row"
        );
        let live_unresolved_imports = compile_fail_fences
            .iter()
            .filter(|fence| fence.expected_error_code.as_ref() == "E0432")
            .collect::<Vec<_>>();
        assert_eq!(live_unresolved_imports.len(), 1);
        assert_eq!(
            live_unresolved_imports[0].source_path,
            "crates/fs-evidence-runner/src/identity.rs"
        );
        let oracle_unresolved_imports = COMPILE_FAIL_TEMPLATES_V1
            .iter()
            .filter(|case| case.expected_error_code == "E0432")
            .collect::<Vec<_>>();
        assert_eq!(oracle_unresolved_imports.len(), 1);
        assert_eq!(oracle_unresolved_imports[0].module, "identity");
        assert_eq!(
            oracle_unresolved_imports[0].case_name,
            "no-standalone-root-for-root-free-evaluator-members"
        );
        assert_eq!(
            oracle_unresolved_imports[0].source_path,
            "crates/fs-evidence-runner/src/identity.rs"
        );

        let bare = compile_fail_fences_v1(
            &[(
                "synthetic.rs",
                "/// ```compile_fail\n/// let _ = 1;\n/// ```\n",
            )],
            1,
        )
        .unwrap_err();
        assert!(bare.contains("bare compile_fail fence at synthetic.rs:1"));
        let malformed = compile_fail_fences_v1(
            &[(
                "synthetic.rs",
                "/// ```compile_fail,E616\n/// let _ = 1;\n/// ```\n",
            )],
            1,
        )
        .unwrap_err();
        assert!(malformed.contains("malformed or multi-code compile_fail fence"));
        let multi_code = compile_fail_fences_v1(
            &[(
                "synthetic.rs",
                "/// ```compile_fail,E0616,E0308\n/// let _ = 1;\n/// ```\n",
            )],
            1,
        )
        .unwrap_err();
        assert!(multi_code.contains("malformed or multi-code compile_fail fence"));
        assert_eq!(
            first.case_count(BaseCoverageManifestClassV1::ManifestContract),
            BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1
        );
        assert_eq!(first.cases().len(), 324);
        assert_eq!(
            BaseCoverageManifestClassV1::ALL,
            [
                BaseCoverageManifestClassV1::Unit,
                BaseCoverageManifestClassV1::CompileFailDoctest,
                BaseCoverageManifestClassV1::ManifestContract,
                BaseCoverageManifestClassV1::ProjectionE2e,
                BaseCoverageManifestClassV1::RuntimeLogging,
                BaseCoverageManifestClassV1::SourceClosure,
                BaseCoverageManifestClassV1::ExternalE2eScript,
                BaseCoverageManifestClassV1::ExternalMutation,
                BaseCoverageManifestClassV1::ExternalGovernance,
                BaseCoverageManifestClassV1::Boundary,
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                BaseCoverageManifestClassV1::SchemaDescriptor,
                BaseCoverageManifestClassV1::Mutation,
                BaseCoverageManifestClassV1::NoMockIntegration,
            ]
        );
        assert_eq!(
            BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES,
            [
                BaseCoverageManifestClassV1::Unit,
                BaseCoverageManifestClassV1::Boundary,
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                BaseCoverageManifestClassV1::SchemaDescriptor,
                BaseCoverageManifestClassV1::Mutation,
                BaseCoverageManifestClassV1::NoMockIntegration,
            ]
        );
        assert_eq!(
            BaseCoverageManifestClassV1::EXTERNALLY_OWNED,
            [
                BaseCoverageManifestClassV1::ExternalE2eScript,
                BaseCoverageManifestClassV1::ExternalMutation,
                BaseCoverageManifestClassV1::ExternalGovernance,
            ]
        );
        for (index, case) in first.cases().iter().enumerate() {
            assert_eq!(case.ordinal(), u32::try_from(index + 1).unwrap());
            assert_eq!(first.case(case.id()), Some(case));
            assert!(
                case.id().starts_with(case.class().stable_prefix()),
                "{}:{:?}",
                case.id(),
                case.class()
            );
        }
    }

    #[test]
    fn exact_base_reconstruction_rejects_missing_extra_duplicate_and_reordering() {
        let frozen = BaseCoverageManifestV1::frozen().expect("frozen");
        let exact = frozen
            .cases()
            .iter()
            .map(|case| case.declaration().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&exact).unwrap(),
            frozen
        );

        let mut missing = exact.clone();
        missing.pop();
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&missing)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut extra = exact.clone();
        extra.push(extension(
            BaseCoverageManifestClassV1::Unit,
            "unit:coverage:extra",
            "crates/fs-evidence-runner/src/coverage.rs",
        ));
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&extra)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let mut duplicate = exact.clone();
        duplicate[1] = duplicate[0].clone();
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&duplicate)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reordered = exact;
        reordered.swap(0, 1);
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&reordered)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let exact = frozen
            .cases()
            .iter()
            .map(|case| case.declaration().clone())
            .collect::<Vec<_>>();
        let boundary_index = exact
            .iter()
            .position(|case| case.class() == BaseCoverageManifestClassV1::Boundary)
            .expect("nonzero boundary class");
        let boundary = &exact[boundary_index];
        let mut reclassified = exact.clone();
        reclassified[boundary_index] = BaseCoverageCaseDeclarationV1::new(
            BaseCoverageManifestClassV1::Unit,
            format!(
                "unit:{}",
                boundary
                    .id()
                    .strip_prefix("boundary:")
                    .expect("exact boundary prefix")
            ),
            boundary.source_path(),
        )
        .expect("valid but semantically wrong reclassification");
        assert_eq!(
            BaseCoverageManifestV1::reconstruct_exact_base(&reclassified)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn declaration_grammar_and_semantic_mutations_refuse_or_move_identity() {
        for id in [
            "",
            "projection-e2e:",
            "Projection-e2e:case",
            "projection-e2e::case",
            "projection e2e:case",
        ] {
            assert!(
                BaseCoverageCaseDeclarationV1::new(
                    BaseCoverageManifestClassV1::ProjectionE2e,
                    id,
                    "tests/e2e.rs",
                )
                .is_err(),
                "{id:?}"
            );
        }
        let exact_id = format!(
            "projection-e2e:{}",
            "a".repeat(160 - "projection-e2e:".len())
        );
        assert_eq!(exact_id.len(), 160);
        assert!(
            BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                exact_id,
                "tests/e2e.rs",
            )
            .is_ok()
        );
        assert!(
            BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                format!(
                    "projection-e2e:{}",
                    "a".repeat(161 - "projection-e2e:".len())
                ),
                "tests/e2e.rs",
            )
            .is_err()
        );
        assert!(
            BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                "projection-e2e:exact-source-path-bound",
                "a".repeat(240),
            )
            .is_ok()
        );
        assert!(
            BaseCoverageCaseDeclarationV1::new(
                BaseCoverageManifestClassV1::ProjectionE2e,
                "projection-e2e:one-over-source-path-bound",
                "a".repeat(241),
            )
            .is_err()
        );
        for path in ["", "/absolute", "../escape", "a/../b", "a//b", "a\\b"] {
            assert!(
                BaseCoverageCaseDeclarationV1::new(
                    BaseCoverageManifestClassV1::ProjectionE2e,
                    "projection-e2e:case",
                    path,
                )
                .is_err(),
                "{path:?}"
            );
        }

        let first = extension(
            BaseCoverageManifestClassV1::ProjectionE2e,
            "projection-e2e:a",
            "tests/a.rs",
        );
        let changed_path = extension(
            BaseCoverageManifestClassV1::ProjectionE2e,
            "projection-e2e:a",
            "tests/b.rs",
        );
        assert_ne!(
            BaseCoverageManifestV1::with_exact_extensions(&[first])
                .unwrap()
                .root(),
            BaseCoverageManifestV1::with_exact_extensions(&[changed_path])
                .unwrap()
                .root()
        );
    }

    #[test]
    fn extension_constructor_requires_external_class_order_and_global_uniqueness() {
        for class in BaseCoverageManifestClassV1::RUST_TEST_EVIDENCE_CLASSES {
            let crate_owned = extension(
                class,
                &format!("{}coverage:forged-extension", class.stable_prefix()),
                "tests/e2e.rs",
            );
            assert_eq!(
                BaseCoverageManifestV1::with_exact_extensions(&[crate_owned])
                    .unwrap_err()
                    .kind(),
                ConstructionErrorKindV2::Incompatible,
                "{class:?}"
            );
        }

        let a = extension(
            BaseCoverageManifestClassV1::ProjectionE2e,
            "projection-e2e:a",
            "tests/e2e.rs",
        );
        let b = extension(
            BaseCoverageManifestClassV1::RuntimeLogging,
            "runtime-logging:b",
            "tests/e2e.rs",
        );
        assert!(BaseCoverageManifestV1::with_exact_extensions(&[a.clone(), b.clone()]).is_ok());
        assert_eq!(
            BaseCoverageManifestV1::with_exact_extensions(&[b, a.clone()])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(
            BaseCoverageManifestV1::with_exact_extensions(&[a.clone(), a])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let maximum = (0..1_024)
            .map(|index| {
                extension(
                    BaseCoverageManifestClassV1::ExternalGovernance,
                    &format!("external-governance:case-{index:04}"),
                    "scripts/check-governance.sh",
                )
            })
            .collect::<Vec<_>>();
        assert!(BaseCoverageManifestV1::with_exact_extensions(&maximum).is_ok());
        let mut one_over = maximum;
        one_over.push(extension(
            BaseCoverageManifestClassV1::ExternalGovernance,
            "external-governance:case-1024",
            "scripts/check-governance.sh",
        ));
        assert_eq!(
            BaseCoverageManifestV1::with_exact_extensions(&one_over)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
    }

    #[test]
    fn selection_accepts_exact_subsets_and_rejects_unknown_duplicate_and_reordered_ids() {
        let manifest = BaseCoverageManifestV1::frozen().expect("frozen");
        let first = manifest.cases()[0].id();
        let last = manifest.cases().last().unwrap().id();
        let selected = manifest
            .select_exact(&[first, last])
            .expect("ordered subset");
        assert_eq!(
            selected.source_case_ids(),
            &[Box::<str>::from(first), Box::<str>::from(last)]
        );
        assert!(manifest.select_exact(&[]).is_ok());
        assert_eq!(
            manifest
                .select_exact(&["unit:unknown:case"])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        assert_eq!(
            manifest.select_exact(&[first, first]).unwrap_err().kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            manifest.select_exact(&[last, first]).unwrap_err().kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
    }

    #[test]
    fn checked_join_accepts_empty_and_mixed_exact_results() {
        let manifest = BaseCoverageManifestV1::frozen().expect("frozen");
        let empty = manifest.select_exact(&[]).expect("empty selection");
        let empty_report =
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &empty, &[]).unwrap();
        assert!(empty_report.is_green());
        assert!(empty_report.results().is_empty());

        let ids = [manifest.cases()[0].id(), manifest.cases()[1].id()];
        let selected = manifest.select_exact(&ids).expect("selection");
        let results = [
            result(
                &manifest,
                ids[0],
                BaseCoveragePresentedOutcomeV1::PositiveMatched,
                "positive",
            ),
            result(
                &manifest,
                ids[1],
                BaseCoveragePresentedOutcomeV1::ExpectedRefusalMatched,
                "refusal",
            ),
        ];
        let report =
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &results).unwrap();
        assert!(report.is_green());
        assert_eq!(report.positive_matched(), 1);
        assert_eq!(report.expected_refusals_matched(), 1);
        assert_eq!(report.expected_unsupported_matched(), 0);
        assert_eq!(report.unexpected_mismatches(), 0);
    }

    #[test]
    fn checked_join_rejects_missing_and_selected_extra_results() {
        let manifest = BaseCoverageManifestV1::frozen().expect("frozen");
        let ids = [manifest.cases()[0].id(), manifest.cases()[1].id()];
        let selected = manifest.select_exact(&ids[..1]).expect("one selected");
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &[])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let results = [
            result(
                &manifest,
                ids[0],
                BaseCoveragePresentedOutcomeV1::PositiveMatched,
                "a",
            ),
            result(
                &manifest,
                ids[1],
                BaseCoveragePresentedOutcomeV1::PositiveMatched,
                "b",
            ),
        ];
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &results)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
    }

    #[test]
    fn checked_join_rejects_stale_unmapped_and_multiply_reported_ids() {
        let manifest = BaseCoverageManifestV1::frozen().expect("frozen");
        let id = manifest.cases()[0].id();
        let selected = manifest.select_exact(&[id]).expect("selection");
        let stale = BaseCoveragePresentedResultV1::new(
            root("stale-manifest"),
            id,
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            root("evidence"),
        )
        .unwrap();
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &[stale])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let unmapped = BaseCoveragePresentedResultV1::new(
            manifest.root(),
            "unit:unknown:case",
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            root("evidence"),
        )
        .unwrap();
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &[unmapped])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        let repeated = result(
            &manifest,
            id,
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            "evidence",
        );
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(
                &manifest,
                &selected,
                &[repeated.clone(), repeated],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
    }

    #[test]
    fn checked_join_rejects_reordered_results() {
        let manifest = BaseCoverageManifestV1::frozen().expect("frozen");
        let ids = [manifest.cases()[0].id(), manifest.cases()[1].id()];
        let selected = manifest.select_exact(&ids).expect("selection");
        let first = result(
            &manifest,
            ids[0],
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            "first",
        );
        let second = result(
            &manifest,
            ids[1],
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            "second",
        );
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest, &selected, &[second, first],)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
    }

    #[test]
    fn result_selection_manifest_and_report_roots_bind_every_semantic_field() {
        let extension_a = extension(
            BaseCoverageManifestClassV1::ProjectionE2e,
            "projection-e2e:a",
            "tests/a.rs",
        );
        let extension_b = extension(
            BaseCoverageManifestClassV1::ProjectionE2e,
            "projection-e2e:b",
            "tests/a.rs",
        );
        let manifest_a =
            BaseCoverageManifestV1::with_exact_extensions(&[extension_a]).expect("manifest a");
        let manifest_b =
            BaseCoverageManifestV1::with_exact_extensions(&[extension_b]).expect("manifest b");
        assert_ne!(manifest_a.root(), manifest_b.root());

        let id = manifest_a.cases().last().unwrap().id();
        let selection_a = manifest_a.select_exact(&[id]).expect("selection a");
        let selection_empty = manifest_a.select_exact(&[]).expect("empty");
        assert_ne!(selection_a.root(), selection_empty.root());

        let positive = result(
            &manifest_a,
            id,
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            "evidence-a",
        );
        let mismatch = result(
            &manifest_a,
            id,
            BaseCoveragePresentedOutcomeV1::UnexpectedMismatch,
            "evidence-a",
        );
        let other_evidence = result(
            &manifest_a,
            id,
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            "evidence-b",
        );
        let stale_manifest = BaseCoveragePresentedResultV1::new(
            manifest_b.root(),
            id,
            BaseCoveragePresentedOutcomeV1::PositiveMatched,
            root("evidence-a"),
        )
        .unwrap();
        assert_ne!(positive.root(), mismatch.root());
        assert_ne!(positive.root(), other_evidence.root());
        assert_ne!(positive.root(), stale_manifest.root());

        let green =
            BaseCoverageCheckedReportV1::reconstruct(&manifest_a, &selection_a, &[positive])
                .unwrap();
        let red = BaseCoverageCheckedReportV1::reconstruct(&manifest_a, &selection_a, &[mismatch])
            .unwrap();
        assert!(green.is_green());
        assert!(!red.is_green());
        assert_ne!(green.root(), red.root());
        assert_eq!(
            BaseCoverageCheckedReportV1::reconstruct(&manifest_b, &selection_a, &[])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn full_set_close_manifest_exactly_covers_nine_groups_and_twenty_two_facets() {
        let source = super::frozen_full_source_manifest_v1().expect("full source manifest");
        assert_eq!(source.cases().len(), 445);
        assert_eq!(
            [
                source.case_count(BaseCoverageManifestClassV1::ProjectionE2e),
                source.case_count(BaseCoverageManifestClassV1::RuntimeLogging),
                source.case_count(BaseCoverageManifestClassV1::SourceClosure),
                source.case_count(BaseCoverageManifestClassV1::ExternalE2eScript),
                source.case_count(BaseCoverageManifestClassV1::ExternalMutation),
                source.case_count(BaseCoverageManifestClassV1::ExternalGovernance),
            ],
            [98, 1, 15, 5, 1, 1]
        );
        assert_eq!(
            BaseCoverageManifestClassV1::ALL
                .iter()
                .copied()
                .filter(|class| class.is_extension())
                .map(|class| source.case_count(class))
                .sum::<usize>(),
            121
        );
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        assert_eq!(manifest.cells().len(), 445);
        assert_eq!(manifest.cells().len(), source.cases().len());
        assert_eq!(manifest.source_manifest_root(), source.root());
        for group in BaseCoverageCloseGroupV1::ALL {
            assert_ne!(manifest.group_count(group), 0, "{group:?}");
        }
        for facet in BaseCoverageCloseFacetV1::ALL {
            assert_ne!(manifest.facet_count(facet), 0, "{facet:?}");
            let inapplicable = matches!(
                facet,
                BaseCoverageCloseFacetV1::Race
                    | BaseCoverageCloseFacetV1::Trait
                    | BaseCoverageCloseFacetV1::Cancellation
                    | BaseCoverageCloseFacetV1::ReleaseBuiltNoMockE2e
            );
            assert_eq!(
                manifest.applicable_facet_count(facet) == 0,
                inapplicable,
                "{facet:?}"
            );
        }
        assert_eq!(
            super::base_coverage_close_reason_descriptors_v1()
                .iter()
                .map(|reason| (reason.code().code(), reason.name()))
                .collect::<Vec<_>>(),
            vec![
                (1, "race-not-applicable-pure-single-threaded-validator"),
                (2, "trait-not-applicable-no-public-trait-contract"),
                (3, "cancellation-not-applicable-pure-bounded-validator"),
                (4, "release-execution-downstream-owned"),
                (5, "windows-nonascii-alias-locally-unadjudicable"),
            ]
        );
        let downstream = manifest
            .cells()
            .iter()
            .filter(|cell| {
                cell.execution_scope()
                    == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution
            })
            .collect::<Vec<_>>();
        assert_eq!(downstream.len(), 7);
        for cell in downstream {
            let contribution = cell.downstream_contribution().expect("typed contribution");
            assert_eq!(
                contribution
                    .budgets()
                    .resolved()
                    .logical_work()
                    .unit()
                    .unit(),
                LogicalUnitV2::Operations
            );
            assert!(
                contribution.budgets().max_parallel_children()
                    <= contribution.budgets().max_child_processes()
            );
            let expected_script = match cell.source_class() {
                BaseCoverageManifestClassV1::ExternalE2eScript => cell.source_path(),
                BaseCoverageManifestClassV1::ExternalMutation => {
                    assert_eq!(
                        cell.source_path(),
                        "crates/fs-evidence-runner/src/projection.rs"
                    );
                    "scripts/ci/canonical_evidence_runner_v2.sh"
                }
                BaseCoverageManifestClassV1::ExternalGovernance => {
                    assert_eq!(
                        cell.source_path(),
                        "crates/fs-evidence-runner/src/dependency.rs"
                    );
                    "scripts/ci/e2e_runner_v2_tool_governance.sh"
                }
                class => {
                    unreachable!("downstream contribution has non-external class {class:?}")
                }
            };
            assert_eq!(contribution.downstream_script(), expected_script);
        }
        assert!(manifest.cells().iter().all(|cell| {
            (cell.execution_scope()
                == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution)
                == cell.downstream_contribution().is_some()
        }));
    }

    #[test]
    fn full_set_close_manifest_exact_reconstruction_rejects_all_sequence_and_semantic_mutants() {
        let source = super::frozen_full_source_manifest_v1().expect("full source");
        let manifest =
            BaseCoverageCloseManifestV1::reconstruct_full(&source).expect("full close manifest");
        let exact = manifest
            .cells()
            .iter()
            .map(|cell| cell.declaration().clone())
            .collect::<Vec<_>>();
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &exact).unwrap(),
            manifest
        );

        let mut missing = exact.clone();
        missing.pop();
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &missing)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut extra = exact.clone();
        extra.push(exact[0].clone());
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &extra)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let mut duplicate = exact.clone();
        duplicate[1] = duplicate[0].clone();
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &duplicate)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reordered = exact.clone();
        reordered.swap(0, 1);
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &reordered)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let mut reclassified = exact.clone();
        reclassified[0].facet = BaseCoverageCloseFacetV1::Boundary;
        reclassified[0].group = BaseCoverageCloseFacetV1::Boundary.group();
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &reclassified)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let race_index = exact
            .iter()
            .position(|cell| cell.facet() == BaseCoverageCloseFacetV1::Race)
            .unwrap();
        let mut wrong_scope = exact.clone();
        wrong_scope[race_index].execution_scope = BaseCoverageCloseExecutionScopeV1::CrateTest;
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &wrong_scope)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let downstream_index = exact
            .iter()
            .position(|cell| cell.downstream_contribution().is_some())
            .unwrap();
        let mut missing_contribution = exact;
        missing_contribution[downstream_index].downstream_contribution = None;
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_exact_full(&source, &missing_contribution)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn full_set_close_report_reconstructs_exact_green_partitions_and_adversarial_sum() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let results = matched_close_results(&manifest);
        let report =
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &results).expect("green report");
        assert!(report.is_green());
        assert_eq!(report.results().len(), manifest.cells().len());
        assert!(!report.eligible_positive_ids().is_empty());
        assert!(!report.eligible_expected_refusal_ids().is_empty());
        assert!(!report.eligible_expected_failure_ids().is_empty());
        assert!(!report.eligible_mutation_ids().is_empty());
        assert_eq!(
            report.matched_positive_ids(),
            report.eligible_positive_ids()
        );
        assert_eq!(
            report.matched_expected_refusal_ids(),
            report.eligible_expected_refusal_ids()
        );
        assert_eq!(
            report.matched_expected_failure_ids(),
            report.eligible_expected_failure_ids()
        );
        assert_eq!(
            report.matched_mutation_ids(),
            report.eligible_mutation_ids()
        );
        assert_eq!(
            usize::try_from(report.adversarial_eligible()).unwrap(),
            report.eligible_expected_refusal_ids().len()
                + report.eligible_expected_failure_ids().len()
                + report.eligible_mutation_ids().len()
        );
        assert_eq!(report.adversarial_matched(), report.adversarial_eligible());
        assert_eq!(report.matched_unsupported(), report.expected_unsupported());
        assert_eq!(
            report.matched_inapplicable(),
            report.expected_inapplicable()
        );
        assert!(report.first_divergence_id().is_none());
        assert!(report.first_divergence_root().is_none());
    }

    #[test]
    fn full_set_close_report_refuses_wrong_manifest_cell_reason_scope_and_order() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let exact = matched_close_results(&manifest);

        let mut stale = exact.clone();
        stale[0].close_manifest_root = root("stale-close-manifest");
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &stale)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_cell = exact.clone();
        wrong_cell[0].cell_root = root("wrong-cell");
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &wrong_cell)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let reason_index = exact
            .iter()
            .position(|result| {
                result.expected_reason()
                    == Some(BaseCoverageCloseReasonCodeV1::WindowsNonasciiAliasLocallyUnadjudicable)
            })
            .unwrap();
        let mut wrong_reason = exact.clone();
        wrong_reason[reason_index].expected_reason =
            Some(BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned);
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &wrong_reason)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut wrong_scope = exact.clone();
        wrong_scope[reason_index].execution_scope = BaseCoverageCloseExecutionScopeV1::CrateTest;
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &wrong_scope)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut reordered = exact.clone();
        reordered.swap(0, 1);
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &reordered)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let downstream_index = exact
            .iter()
            .position(|result| {
                result.evidence().kind()
                    == BaseCoverageCloseEvidenceKindV1::ImmutableDownstreamContribution
            })
            .unwrap();
        let mut wrong_contribution = exact;
        wrong_contribution[downstream_index].evidence.root = root("wrong-contribution");
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &wrong_contribution)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn full_set_close_report_is_fail_closed_for_mismatch_execution_failure_skip_and_partial_rows() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let mut results = matched_close_results(&manifest);
        let mismatch_index = manifest
            .cells()
            .iter()
            .position(|cell| cell.partition() == BaseCoverageClosePartitionV1::Positive)
            .unwrap();
        let failure_index = mismatch_index + 1;
        let skip_index = mismatch_index + 2;
        let mismatch_cell = &manifest.cells()[mismatch_index];
        results[mismatch_index] = BaseCoverageClosePresentedResultV1::unexpected_mismatch(
            &manifest,
            mismatch_cell,
            divergent_decision(mismatch_cell.expected_decision()),
            None,
            close_evidence(&manifest, mismatch_cell, "mismatch"),
        )
        .unwrap();
        let failure_cell = &manifest.cells()[failure_index];
        results[failure_index] = BaseCoverageClosePresentedResultV1::execution_failure(
            &manifest,
            failure_cell,
            close_evidence(&manifest, failure_cell, "execution-failure"),
        )
        .unwrap();
        let skip_cell = &manifest.cells()[skip_index];
        results[skip_index] = BaseCoverageClosePresentedResultV1::unexplained_skip(
            &manifest,
            skip_cell,
            close_evidence(&manifest, skip_cell, "skip"),
        )
        .unwrap();
        let first_root = results[mismatch_index].root();
        let report =
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &results).expect("red report");
        assert!(!report.is_green());
        assert_eq!(
            report.unexpected_mismatch_ids(),
            &[Box::<str>::from(mismatch_cell.source_case_id())]
        );
        assert_eq!(
            report.execution_failure_ids(),
            &[Box::<str>::from(failure_cell.source_case_id())]
        );
        assert_eq!(
            report.unexplained_skip_ids(),
            &[Box::<str>::from(skip_cell.source_case_id())]
        );
        assert_eq!(
            report.first_divergence_id(),
            Some(mismatch_cell.source_case_id())
        );
        assert_eq!(report.first_divergence_root(), Some(first_root));
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &[])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        results.pop();
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &results)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
    }

    #[test]
    fn full_set_close_external_rows_are_inapplicable_immutable_contributions_never_passes() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let external = manifest
            .cells()
            .iter()
            .filter(|cell| cell.source_class() == BaseCoverageManifestClassV1::ExternalE2eScript)
            .collect::<Vec<_>>();
        assert_eq!(external.len(), 5);
        for cell in external {
            assert_eq!(
                cell.facet(),
                BaseCoverageCloseFacetV1::ImmutableE2eContribution
            );
            assert_eq!(
                cell.execution_scope(),
                BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution
            );
            assert_eq!(cell.partition(), BaseCoverageClosePartitionV1::Inapplicable);
            assert_eq!(
                cell.expected_decision(),
                BaseCoverageCloseDecisionV1::Inapplicable
            );
            assert_eq!(
                cell.expected_reason(),
                Some(BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned)
            );
            let contribution = cell.downstream_contribution().unwrap();
            let expected_owner = match cell.source_case_id() {
                "external-e2e:publication-state-v2" => "frankensim-epic-foundations-huq.24.2.2.2",
                "external-e2e:publication-v2" => "frankensim-epic-foundations-huq.24.2.2.3.3",
                "external-e2e:verifier-v2" => "frankensim-epic-foundations-huq.24.3.3.3.3",
                "external-e2e:canonical-runner-v2" => "frankensim-epic-foundations-huq.24.4.1.4",
                "external-e2e:rjoq-handoff-v1" => "frankensim-epic-foundations-huq.24.5.3.1",
                other => panic!("unexpected external ID {other}"),
            };
            assert_eq!(contribution.downstream_owner(), expected_owner);
            assert_eq!(
                contribution.no_claim(),
                "downstream-contribution-is-not-execution-proof"
            );
        }
    }

    #[test]
    fn full_set_close_first_divergence_and_roots_bind_every_semantic_partition() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let first = manifest
            .cells()
            .iter()
            .position(|cell| cell.execution_scope() == BaseCoverageCloseExecutionScopeV1::CrateTest)
            .unwrap();
        let mut a = matched_close_results(&manifest);
        let report_a = BaseCoverageCloseReportV1::reconstruct_full(&manifest, &a).unwrap();
        a[first] = BaseCoverageClosePresentedResultV1::matched(
            &manifest,
            &manifest.cells()[first],
            BaseCoverageCloseResultEvidenceV1::owned_harness_execution(
                root("changed-evidence"),
                Some("artifacts/close/owned.log".to_owned()),
            )
            .unwrap(),
        )
        .unwrap();
        let report_b = BaseCoverageCloseReportV1::reconstruct_full(&manifest, &a).unwrap();
        assert!(report_a.is_green());
        assert!(report_b.is_green());
        assert_ne!(
            report_a.results()[first].root(),
            report_b.results()[first].root()
        );
        assert_ne!(report_a.root(), report_b.root());
        assert_eq!(
            report_b.results()[first].evidence().retained_artifact(),
            Some("artifacts/close/owned.log")
        );
    }

    #[test]
    fn full_set_close_empty_or_base_only_manifest_cannot_mint_full_set_authority() {
        let base = BaseCoverageManifestV1::frozen().expect("base source manifest");
        assert_eq!(
            BaseCoverageCloseManifestV1::reconstruct_full(&base)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let full = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        assert_eq!(
            BaseCoverageCloseReportV1::reconstruct_full(&full, &[])
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
    }

    #[test]
    fn five_explicits_numeric_domain_and_unit_references_are_exact() {
        let physical = BaseCoverageCloseNumericUnitV1::physical(
            UnitV2::from_parts(1, 1_000, [1, 0, -1, 0, 0, 0, 0]).expect("canonical physical unit"),
        );
        let logical = BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Operations)
            .expect("fixed logical unit");
        let values = vec![
            NumericValueV2::I8(i8::MIN),
            NumericValueV2::I16(i16::MIN),
            NumericValueV2::I32(i32::MIN),
            NumericValueV2::I64(i64::MIN),
            NumericValueV2::I128(i128::MIN),
            NumericValueV2::U8(u8::MAX),
            NumericValueV2::U16(u16::MAX),
            NumericValueV2::U32(u32::MAX),
            NumericValueV2::U64(u64::MAX),
            NumericValueV2::U128(u128::MAX),
            NumericValueV2::Rational(RationalV2::new(-7, 13).expect("rational")),
            NumericValueV2::Decimal(DecimalV2::new(-12345, 7).expect("decimal")),
            NumericValueV2::F32Bits(F32BitsV2::from_bits(0x7fc0_0042)),
            NumericValueV2::F64Bits(F64BitsV2::from_bits(0x7ff8_0000_0000_0042)),
        ];
        let inputs = values
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                numeric(
                    format!("numeric-{index:02}"),
                    value,
                    if index < 7 { physical } else { logical },
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            inputs
                .iter()
                .map(|value| value.value().wire_tag())
                .collect::<Vec<_>>(),
            (1_u16..=14).collect::<Vec<_>>()
        );
        let five = five_from_template(
            inputs.clone(),
            vec![],
            super::frozen_local_close_budget_set().expect("local budget profile"),
        )
        .expect("all numeric variants fit the Five Explicits frame");
        assert_eq!(five.numeric_inputs().len(), 14);
        assert!(five.numeric_grants().is_empty());
        assert!(five.numeric_observations().is_empty());
        let variant_roots = inputs
            .into_iter()
            .map(|value| {
                five_from_template(
                    vec![value],
                    vec![],
                    super::frozen_local_close_budget_set().expect("local budget profile"),
                )
                .expect("one exact numeric variant")
                .numeric_inputs_root()
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(
            variant_roots.len(),
            14,
            "every NumericValueV2 variant tag and payload enters the profile root"
        );
        let payload_mutations = vec![
            (NumericValueV2::I8(-2), NumericValueV2::I8(-1)),
            (NumericValueV2::I16(-2), NumericValueV2::I16(-1)),
            (NumericValueV2::I32(-2), NumericValueV2::I32(-1)),
            (NumericValueV2::I64(-2), NumericValueV2::I64(-1)),
            (NumericValueV2::I128(-2), NumericValueV2::I128(-1)),
            (NumericValueV2::U8(1), NumericValueV2::U8(2)),
            (NumericValueV2::U16(1), NumericValueV2::U16(2)),
            (NumericValueV2::U32(1), NumericValueV2::U32(2)),
            (NumericValueV2::U64(1), NumericValueV2::U64(2)),
            (NumericValueV2::U128(1), NumericValueV2::U128(2)),
            (
                NumericValueV2::Rational(RationalV2::new(1, 3).expect("rational a")),
                NumericValueV2::Rational(RationalV2::new(2, 3).expect("rational b")),
            ),
            (
                NumericValueV2::Decimal(DecimalV2::new(1, 2).expect("decimal a")),
                NumericValueV2::Decimal(DecimalV2::new(1, 3).expect("decimal b")),
            ),
            (
                NumericValueV2::F32Bits(F32BitsV2::from_bits(1)),
                NumericValueV2::F32Bits(F32BitsV2::from_bits(2)),
            ),
            (
                NumericValueV2::F64Bits(F64BitsV2::from_bits(1)),
                NumericValueV2::F64Bits(F64BitsV2::from_bits(2)),
            ),
        ];
        for (index, (first, second)) in payload_mutations.into_iter().enumerate() {
            let first = five_from_template(
                vec![numeric(format!("payload-{index:02}"), first, logical)],
                vec![],
                super::frozen_local_close_budget_set().expect("local budget profile"),
            )
            .expect("first payload");
            let second = five_from_template(
                vec![numeric(format!("payload-{index:02}"), second, logical)],
                vec![],
                super::frozen_local_close_budget_set().expect("local budget profile"),
            )
            .expect("second payload");
            assert_ne!(
                first.numeric_inputs_root(),
                second.numeric_inputs_root(),
                "every NumericValueV2 payload field enters the profile root"
            );
        }

        let registered = LogicalUnitV2::from_tag(16, Some(7)).expect("registered unit syntax");
        assert_eq!(
            BaseCoverageCloseLogicalUnitReferenceV1::new(registered, None)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseCoverageCloseLogicalUnitReferenceV1::new(
                LogicalUnitV2::Operations,
                Some(root("unexpected-registry")),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let registered_unit =
            BaseCoverageCloseNumericUnitV1::registered_logical(registered, root("registry-a"))
                .expect("registry-bound logical unit");
        assert!(matches!(
            registered_unit,
            BaseCoverageCloseNumericUnitV1::Logical(reference)
                if reference.unit() == registered
                    && reference.registry_identity() == Some(root("registry-a"))
        ));
    }

    #[test]
    fn five_explicits_numeric_surface_bounds_order_and_registered_identity_are_exact() {
        let budgets = super::frozen_local_close_budget_set().expect("local budget profile");
        let empty = five_from_template(vec![], vec![], budgets).expect("exact-empty profiles");
        assert_ne!(empty.numeric_inputs_root(), empty.numeric_grants_root());
        assert_ne!(
            empty.numeric_inputs_root(),
            empty.numeric_observations_root()
        );
        assert_ne!(
            empty.numeric_grants_root(),
            empty.numeric_observations_root()
        );
        assert!(
            five_from_template(
                vec![numeric(
                    "one",
                    NumericValueV2::U8(1),
                    BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count)
                        .expect("logical unit"),
                )],
                vec![],
                budgets,
            )
            .is_ok()
        );
        let one = numeric(
            "one",
            NumericValueV2::U8(1),
            BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count).expect("logical unit"),
        );
        assert!(five_from_template_with_grants(vec![], vec![one.clone()], vec![], budgets).is_ok());
        assert!(five_from_template_with_grants(vec![], vec![], vec![one], budgets).is_ok());

        let exact = (0_u128..64)
            .map(|index| {
                numeric(
                    format!("n{index:02}"),
                    NumericValueV2::U128(u128::from(index)),
                    BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count)
                        .expect("logical unit"),
                )
            })
            .collect::<Vec<_>>();
        assert!(
            five_from_template_with_grants(exact.clone(), exact.clone(), exact.clone(), budgets)
                .is_ok()
        );
        let mut one_over = exact.clone();
        one_over.push(numeric(
            "n64",
            NumericValueV2::U128(64),
            BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count).expect("logical unit"),
        ));
        assert_eq!(
            five_from_template(one_over, vec![], budgets)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let mut one_over_grants = exact.clone();
        one_over_grants.push(numeric(
            "n64",
            NumericValueV2::U128(64),
            BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count).expect("logical unit"),
        ));
        assert_eq!(
            five_from_template_with_grants(vec![], one_over_grants, vec![], budgets)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let mut one_over_observations = exact.clone();
        one_over_observations.push(numeric(
            "n64",
            NumericValueV2::U128(64),
            BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count).expect("logical unit"),
        ));
        assert_eq!(
            five_from_template_with_grants(vec![], vec![], one_over_observations, budgets)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let duplicate = vec![exact[0].clone(), exact[0].clone()];
        assert_eq!(
            five_from_template_with_grants(vec![], duplicate, vec![], budgets)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let reordered = vec![exact[1].clone(), exact[0].clone()];
        assert_eq!(
            five_from_template_with_grants(vec![], vec![], reordered, budgets)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let registered = LogicalUnitV2::from_tag(16, Some(9)).expect("registered unit syntax");
        let with_registry = |identity| {
            vec![numeric(
                "registered-value",
                NumericValueV2::U64(17),
                BaseCoverageCloseNumericUnitV1::registered_logical(registered, identity)
                    .expect("registry-bound unit"),
            )]
        };
        let a = five_from_template(with_registry(root("registry-a")), vec![], budgets)
            .expect("registry a");
        let b = five_from_template(with_registry(root("registry-b")), vec![], budgets)
            .expect("registry b");
        assert_ne!(a.root(), b.root());

        let close_manifest =
            BaseCoverageCloseManifestV1::frozen().expect("source-authoritative close manifest");
        let local_cell = close_manifest
            .cells()
            .iter()
            .find(|cell| cell.downstream_contribution().is_none())
            .expect("one locally owned cell");
        let declared = local_cell.five_explicits();
        let unexpected_profile = |partition| {
            let value = numeric(
                "unexpected-value",
                NumericValueV2::U64(1),
                BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count)
                    .expect("logical unit"),
            );
            let (inputs, grants, observations) = match partition {
                BaseCoverageCloseNumericPartitionV1::Inputs => (vec![value], vec![], vec![]),
                BaseCoverageCloseNumericPartitionV1::Grants => (vec![], vec![value], vec![]),
                BaseCoverageCloseNumericPartitionV1::Observations => (vec![], vec![], vec![value]),
            };
            BaseCoverageCloseFiveExplicitsV1::new(
                inputs,
                grants,
                observations,
                declared.seed().clone(),
                declared.budgets(),
                declared.versions().clone(),
                declared.capabilities().clone(),
                declared.no_claim().clone(),
            )
            .expect("internally valid nonempty numeric profile")
        };
        for partition in [
            BaseCoverageCloseNumericPartitionV1::Inputs,
            BaseCoverageCloseNumericPartitionV1::Grants,
            BaseCoverageCloseNumericPartitionV1::Observations,
        ] {
            assert_eq!(
                BaseCoverageCloseCellDeclarationV1::new_with_five_explicits(
                    local_cell.source_ordinal(),
                    local_cell.source_case_id(),
                    local_cell.source_class(),
                    local_cell.source_path(),
                    local_cell.group(),
                    local_cell.facet(),
                    local_cell.execution_scope(),
                    local_cell.partition(),
                    local_cell.expected_decision(),
                    local_cell.expected_reason(),
                    None,
                    unexpected_profile(partition),
                )
                .unwrap_err()
                .kind(),
                ConstructionErrorKindV2::Unexpected,
                "exact reconstruction refuses every numeric profile absent from the source oracle"
            );
        }
    }

    #[test]
    fn five_explicits_budget_axis_catalog_widths_units_and_order_are_exact() {
        assert_eq!(
            BaseCoverageCloseBudgetAxisV1::ALL.map(BaseCoverageCloseBudgetAxisV1::code),
            [1, 2, 3, 4, 5, 6, 7]
        );
        assert_eq!(
            BaseCoverageCloseBudgetAxisV1::ALL.map(BaseCoverageCloseBudgetAxisV1::width),
            [
                BaseCoverageCloseBudgetWidthV1::U64,
                BaseCoverageCloseBudgetWidthV1::U64,
                BaseCoverageCloseBudgetWidthV1::U128,
                BaseCoverageCloseBudgetWidthV1::U32,
                BaseCoverageCloseBudgetWidthV1::U64,
                BaseCoverageCloseBudgetWidthV1::U64,
                BaseCoverageCloseBudgetWidthV1::U64,
            ]
        );
        assert_eq!(
            [
                BaseCoverageCloseBudgetWidthV1::U32.code(),
                BaseCoverageCloseBudgetWidthV1::U64.code(),
                BaseCoverageCloseBudgetWidthV1::U128.code(),
            ],
            [1, 2, 3]
        );
        assert_eq!(
            [
                BaseCoverageCloseBudgetProfileV1::LocalSourceValidation.code(),
                BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution.code(),
            ],
            [1, 2]
        );
        assert_eq!(
            [
                BaseCoverageCloseBudgetProfileV1::LocalSourceValidation.name(),
                BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution.name(),
            ],
            [
                "base-close-local-source-validation-v1",
                "base-close-downstream-source-contribution-v1",
            ]
        );
        assert_eq!(
            BaseCoverageCloseBudgetAxisV1::ALL.map(BaseCoverageCloseBudgetAxisV1::fixed_unit),
            [
                Some(LogicalUnitV2::Nanoseconds),
                Some(LogicalUnitV2::LogicalBytes),
                None,
                Some(LogicalUnitV2::Count),
                Some(LogicalUnitV2::EncodedBytes),
                Some(LogicalUnitV2::EncodedBytes),
                Some(LogicalUnitV2::EncodedBytes),
            ]
        );
        let profile = super::frozen_local_close_budget_set().expect("local budget profile");
        assert_eq!(
            profile.profile(),
            BaseCoverageCloseBudgetProfileV1::LocalSourceValidation
        );
        assert_eq!(
            profile.rows().map(BaseCoverageCloseTypedBudgetV1::hard),
            [
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U32(1),
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
            ]
        );
        assert_eq!(
            profile.rows().map(BaseCoverageCloseTypedBudgetV1::soft),
            [
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                BaseCoverageCloseBudgetValueV1::U32(0),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
            ]
        );
        for (expected, row) in BaseCoverageCloseBudgetAxisV1::ALL
            .iter()
            .zip(profile.rows())
        {
            assert_eq!(row.axis(), *expected);
            assert_eq!(row.hard().width(), expected.width());
            assert_eq!(row.soft().width(), expected.width());
            if let Some(expected_unit) = expected.fixed_unit() {
                assert_eq!(row.unit(), fixed_logical(expected_unit));
            }
        }
        let downstream =
            super::frozen_downstream_close_budget_set().expect("downstream budget profile");
        assert_eq!(
            downstream.rows().map(BaseCoverageCloseTypedBudgetV1::hard),
            [
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U32(8),
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
            ]
        );
        assert_eq!(
            downstream.rows().map(BaseCoverageCloseTypedBudgetV1::soft),
            [
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                BaseCoverageCloseBudgetValueV1::U32(6),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
            ]
        );
        for (expected, row) in BaseCoverageCloseBudgetAxisV1::ALL
            .iter()
            .zip(downstream.rows())
        {
            assert_eq!(row.axis(), *expected);
            assert_eq!(row.hard().width(), expected.width());
            assert_eq!(row.soft().width(), expected.width());
            if let Some(expected_unit) = expected.fixed_unit() {
                assert_eq!(row.unit(), fixed_logical(expected_unit));
            }
        }
    }

    #[test]
    fn five_explicits_budget_bounds_soft_relations_and_shape_refusals_are_exact() {
        assert_eq!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(0),
                BaseCoverageCloseBudgetValueV1::U64(0),
                LogicalUnitV2::Nanoseconds,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(1),
                BaseCoverageCloseBudgetValueV1::U64(0),
                LogicalUnitV2::Nanoseconds,
            )
            .is_ok()
        );
        assert!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(1),
                BaseCoverageCloseBudgetValueV1::U64(1),
                LogicalUnitV2::Nanoseconds,
            )
            .is_ok()
        );
        assert_eq!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(1),
                BaseCoverageCloseBudgetValueV1::U64(2),
                LogicalUnitV2::Nanoseconds,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U32(1),
                BaseCoverageCloseBudgetValueV1::U32(1),
                LogicalUnitV2::Nanoseconds,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            budget_row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(1),
                BaseCoverageCloseBudgetValueV1::U64(1),
                LogicalUnitV2::Seconds,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        for profile in [
            BaseCoverageCloseBudgetProfileV1::LocalSourceValidation,
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
        ] {
            let source = match profile {
                BaseCoverageCloseBudgetProfileV1::LocalSourceValidation => {
                    super::frozen_local_close_budget_set().expect("local profile")
                }
                BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution => {
                    super::frozen_downstream_close_budget_set().expect("downstream profile")
                }
            };
            for axis in BaseCoverageCloseBudgetAxisV1::ALL {
                let index = usize::from(axis.code() - 1);
                let ceiling = profile.hard_ceiling(axis);
                let mut exact = source.rows().to_vec();
                exact[index] = BaseCoverageCloseTypedBudgetV1::new(
                    axis,
                    ceiling,
                    ceiling,
                    exact[index].unit(),
                )
                .expect("exact governing ceiling row");
                assert!(BaseCoverageCloseBudgetSetV1::new(profile, exact).is_ok());
                if let Some(one_over) = one_over_budget_value(ceiling) {
                    let mut over = source.rows().to_vec();
                    over[index] = BaseCoverageCloseTypedBudgetV1::new(
                        axis,
                        one_over,
                        ceiling,
                        over[index].unit(),
                    )
                    .expect("one-over retains width, unit, and soft relation");
                    assert_eq!(
                        BaseCoverageCloseBudgetSetV1::new(profile, over)
                            .unwrap_err()
                            .kind(),
                        ConstructionErrorKindV2::TooLarge
                    );
                }
            }
        }
        assert_eq!(
            BaseCoverageCloseBudgetProfileV1::LocalSourceValidation
                .hard_ceiling(BaseCoverageCloseBudgetAxisV1::LogicalWork),
            BaseCoverageCloseBudgetValueV1::U128(u128::MAX)
        );
        assert_eq!(
            one_over_budget_value(BaseCoverageCloseBudgetValueV1::U128(u128::MAX)),
            None,
            "the exact u128 maximum is admitted and has no representable one-over value"
        );

        let source = super::frozen_local_close_budget_set().expect("local profile");
        let mut missing = source.rows().to_vec();
        missing.pop();
        assert_eq!(
            BaseCoverageCloseBudgetSetV1::new(source.profile(), missing)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut extra = source.rows().to_vec();
        extra.push(source.logs());
        assert_eq!(
            BaseCoverageCloseBudgetSetV1::new(source.profile(), extra)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let mut duplicate = source.rows().to_vec();
        duplicate[1] = duplicate[0];
        assert_eq!(
            BaseCoverageCloseBudgetSetV1::new(source.profile(), duplicate)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reordered = source.rows().to_vec();
        reordered.swap(0, 1);
        assert_eq!(
            BaseCoverageCloseBudgetSetV1::new(source.profile(), reordered)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
    }

    #[test]
    fn five_explicits_downstream_budgets_keep_soft_rows_and_process_shape_independent() {
        let resolved =
            super::frozen_downstream_close_budget_set().expect("downstream budget profile");
        assert_eq!(
            resolved.profile(),
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution
        );
        assert_eq!(
            resolved.processes().hard(),
            BaseCoverageCloseBudgetValueV1::U32(8)
        );
        assert_eq!(
            resolved.processes().soft(),
            BaseCoverageCloseBudgetValueV1::U32(6)
        );
        let contribution =
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 4, 2).expect("shape");
        assert_eq!(contribution.resolved(), resolved);
        assert_eq!(contribution.max_child_processes(), 4);
        assert_eq!(contribution.max_parallel_children(), 2);
        assert_ne!(
            resolved.processes().hard(),
            BaseCoverageCloseBudgetValueV1::U32(contribution.max_child_processes())
        );
        assert_ne!(
            resolved.processes().soft(),
            BaseCoverageCloseBudgetValueV1::U32(contribution.max_parallel_children())
        );

        let mut changed_soft_rows = resolved.rows().to_vec();
        changed_soft_rows[3] = budget_row(
            BaseCoverageCloseBudgetAxisV1::Processes,
            BaseCoverageCloseBudgetValueV1::U32(8),
            BaseCoverageCloseBudgetValueV1::U32(5),
            LogicalUnitV2::Count,
        )
        .expect("independent soft row");
        let changed_soft = BaseCoverageCloseBudgetSetV1::new(
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
            changed_soft_rows,
        )
        .expect("changed explicit soft profile");
        assert!(
            BaseCoverageCloseContributionBudgetsV1::new(changed_soft, 4, 2).is_ok(),
            "process shape does not synthesize or overwrite the independent soft row"
        );
        let resolved_five =
            five_from_template_with_grants(vec![], vec![], vec![], resolved).expect("resolved");
        let changed_soft_five =
            five_from_template_with_grants(vec![], vec![], vec![], changed_soft)
                .expect("changed explicit soft");
        assert_ne!(
            resolved_five.root(),
            changed_soft_five.root(),
            "an independently declared soft value enters the Five Explicits root"
        );
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("close manifest");
        let template = manifest
            .cells()
            .iter()
            .find_map(|cell| cell.downstream_contribution())
            .expect("downstream contribution template");
        let rebuild_contribution = |budgets| {
            BaseCoverageCloseDownstreamContributionV1::new(
                template.literal_expectation_oracle_root(),
                template.semantic_input_root(),
                budgets,
                template.schema_root(),
                template.log_schema_root(),
                template.source_root().clone(),
                template.build_root().clone(),
                template.downstream_owner(),
                template.downstream_driver().clone(),
                template.downstream_script(),
                template.downstream_manifest_path(),
                template.downstream_manifest_root(),
                template.no_claim(),
            )
            .expect("reconstructed downstream contribution")
        };
        let baseline_contribution = rebuild_contribution(contribution);
        let soft_contribution = rebuild_contribution(
            BaseCoverageCloseContributionBudgetsV1::new(changed_soft, 4, 2)
                .expect("changed soft contribution budgets"),
        );
        let total_shape_contribution = rebuild_contribution(
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 5, 2)
                .expect("changed total-child shape"),
        );
        let parallel_shape_contribution = rebuild_contribution(
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 4, 3)
                .expect("changed parallel-child shape"),
        );
        let mut changed_hard_rows = resolved.rows().to_vec();
        changed_hard_rows[3] = budget_row(
            BaseCoverageCloseBudgetAxisV1::Processes,
            BaseCoverageCloseBudgetValueV1::U32(7),
            BaseCoverageCloseBudgetValueV1::U32(6),
            LogicalUnitV2::Count,
        )
        .expect("independent hard row");
        let changed_hard = BaseCoverageCloseBudgetSetV1::new(
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
            changed_hard_rows,
        )
        .expect("changed explicit hard profile");
        let hard_contribution = rebuild_contribution(
            BaseCoverageCloseContributionBudgetsV1::new(changed_hard, 4, 2)
                .expect("changed hard contribution budgets"),
        );
        for mutation in [
            &soft_contribution,
            &hard_contribution,
            &total_shape_contribution,
            &parallel_shape_contribution,
        ] {
            assert_ne!(
                baseline_contribution.root(),
                mutation.root(),
                "downstream hard, soft, total-child, and parallel-child fields enter the contribution root independently"
            );
        }

        assert_eq!(
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 4, 5)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 257, 1)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            BaseCoverageCloseContributionBudgetsV1::new(resolved, 65, 65)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut too_small_process_rows = resolved.rows().to_vec();
        too_small_process_rows[3] = budget_row(
            BaseCoverageCloseBudgetAxisV1::Processes,
            BaseCoverageCloseBudgetValueV1::U32(3),
            BaseCoverageCloseBudgetValueV1::U32(2),
            LogicalUnitV2::Count,
        )
        .expect("bounded process row");
        let too_small_process = BaseCoverageCloseBudgetSetV1::new(
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
            too_small_process_rows,
        )
        .expect("bounded process profile");
        assert_eq!(
            BaseCoverageCloseContributionBudgetsV1::new(too_small_process, 4, 2)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn five_explicits_profiles_and_one_field_mutations_move_roots() {
        let descriptor_oracle = [
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
        let descriptors = base_coverage_close_nominal_root_descriptors_v1();
        assert_eq!(
            descriptors.len(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
        );
        assert_eq!(descriptors.len(), descriptor_oracle.len());
        let mut schema_names = BTreeSet::new();
        let mut domains = BTreeSet::new();
        let mut no_claims = BTreeSet::new();
        for (descriptor, (schema_name, domain, no_claim)) in
            descriptors.iter().zip(descriptor_oracle)
        {
            assert_eq!(descriptor.schema_name(), schema_name);
            assert_eq!(descriptor.domain(), domain);
            assert_eq!(descriptor.no_claim(), no_claim);
            assert_eq!(descriptor.api_generation().code(), 2);
            assert_eq!(descriptor.wire_version().code(), 1);
            assert_eq!(
                descriptor.predecessor_policy(),
                WirePredecessorPolicyV1::NoPredecessor
            );
            assert_eq!(descriptor.predecessor_policy().predecessor(), None);
            assert!(
                schema_names.insert(descriptor.schema_name()),
                "duplicate nominal-root schema name {}",
                descriptor.schema_name()
            );
            assert!(
                domains.insert(descriptor.domain()),
                "duplicate nominal-root domain {}",
                descriptor.domain()
            );
            assert!(
                no_claims.insert(descriptor.no_claim()),
                "duplicate nominal-root no-claim {}",
                descriptor.no_claim()
            );
        }
        assert_eq!(
            schema_names.len(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
        );
        assert_eq!(
            domains.len(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
        );
        assert_eq!(
            no_claims.len(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
        );

        let disposition_oracle = [
            (RuntimeObservationDispositionV1::Observed, 1, "observed"),
            (
                RuntimeObservationDispositionV1::NotObserved,
                2,
                "not-observed",
            ),
            (RuntimeObservationDispositionV1::Deferred, 3, "deferred"),
            (
                RuntimeObservationDispositionV1::Inapplicable,
                4,
                "inapplicable",
            ),
        ];
        let mut disposition_roots = BTreeSet::new();
        for (disposition, code, name) in disposition_oracle {
            assert_eq!(disposition.code(), code);
            assert_eq!(disposition.stable_name(), name);
            assert_eq!(
                RuntimeObservationDispositionV1::try_from_code(code)
                    .expect("exact disposition code"),
                disposition
            );
            assert!(
                disposition_roots.insert(
                    disposition
                        .root()
                        .expect("nominal disposition root")
                        .content_hash()
                ),
                "every disposition has a distinct nominal identity"
            );
        }
        assert_eq!(
            RuntimeObservationDispositionV1::try_from_code(0)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            RuntimeObservationDispositionV1::try_from_code(5)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let not_observed_oracle = [
            (
                NotObservedReasonV1::ExecutionFailedBeforeCompleteness,
                1,
                "execution-failed-before-completeness",
                "fs-evidence-runner.coverage",
                "execution-owned-attempt",
                "a-complete-terminal-attempt-and-runtime-observation",
                DiagnosticCodeV2::RunnerInternalError,
                "not-observed-execution-failure-proves-no-complete-observation-or-success",
            ),
            (
                NotObservedReasonV1::ObservationChannelFailedBeforeCompleteness,
                2,
                "observation-channel-failed-before-completeness",
                "fs-evidence-runner.coverage",
                "execution-owned-observation-channel",
                "a-complete-redacted-runtime-observation-channel",
                DiagnosticCodeV2::RunnerNoData,
                "not-observed-channel-failure-proves-no-complete-observation-or-success",
            ),
            (
                NotObservedReasonV1::UnexpectedUnstartedOrSkipped,
                3,
                "unexpected-unstarted-or-skipped",
                "fs-evidence-runner.coverage",
                "execution-owned-dispatch",
                "a-complete-execution-owned-attempt",
                DiagnosticCodeV2::RunnerNotRun,
                "not-observed-unstarted-or-skipped-proves-no-execution-observation-or-success",
            ),
        ];
        let not_observed_registry =
            NotObservedReasonRegistryV1::frozen().expect("NotObserved registry");
        assert_eq!(
            not_observed_registry.descriptors().len(),
            not_observed_oracle.len()
        );
        for (descriptor, (reason, code, name, owner, scope, prerequisite, diagnostic, no_claim)) in
            not_observed_registry
                .descriptors()
                .iter()
                .zip(not_observed_oracle)
        {
            assert_eq!(descriptor.reason(), reason);
            assert_eq!(reason.code(), code);
            assert_eq!(
                NotObservedReasonV1::try_from_code(code).expect("exact NotObserved code"),
                reason
            );
            assert_eq!(descriptor.name(), name);
            assert_eq!(descriptor.owner(), owner);
            assert_eq!(descriptor.scope(), scope);
            assert_eq!(descriptor.prerequisite(), prerequisite);
            assert_eq!(descriptor.diagnostic(), diagnostic);
            assert_eq!(descriptor.no_claim(), no_claim);
        }
        assert_eq!(
            NotObservedReasonV1::try_from_code(0).unwrap_err().kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            NotObservedReasonV1::try_from_code(4).unwrap_err().kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        assert_eq!(
            NotObservedReasonRegistryV1::frozen()
                .expect("deterministic NotObserved registry")
                .root(),
            not_observed_registry.root()
        );

        let deferred_registry = DeferredReasonRegistryV1::frozen().expect("Deferred registry");
        let deferred_descriptor = deferred_registry
            .descriptor(DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner);
        assert_eq!(deferred_registry.descriptors().len(), 1);
        assert_eq!(
            deferred_descriptor.reason(),
            DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        );
        assert_eq!(deferred_descriptor.reason().code(), 1);
        assert_eq!(
            DeferredReasonV1::try_from_code(1).expect("exact Deferred code"),
            DeferredReasonV1::ImmutableContributionAwaitsDesignatedReleaseOwner
        );
        assert_eq!(
            deferred_descriptor.name(),
            "immutable-contribution-awaits-designated-release-owner"
        );
        assert_eq!(deferred_descriptor.owner(), "fs-evidence-runner.coverage");
        assert_eq!(
            deferred_descriptor.scope(),
            "immutable-downstream-contribution"
        );
        assert_eq!(
            deferred_descriptor.prerequisite(),
            "release-built-no-mock-designated-owner-execution"
        );
        assert_eq!(
            deferred_descriptor.diagnostic(),
            DiagnosticCodeV2::RunnerNotRun
        );
        assert_eq!(
            deferred_descriptor.no_claim(),
            "deferred-contribution-proves-no-designated-owner-execution-or-success"
        );
        assert_eq!(
            DeferredReasonV1::try_from_code(0).unwrap_err().kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            DeferredReasonV1::try_from_code(2).unwrap_err().kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        assert_ne!(
            deferred_registry.root().content_hash(),
            not_observed_registry.root().content_hash(),
            "reason registries use distinct nominal domains and frames"
        );

        let schema_impact_oracle = [
            (
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                1,
                "new-v1-no-predecessor",
            ),
            (
                CanonicalSchemaImpactDispositionV1::UnchangedV1,
                2,
                "unchanged-v1",
            ),
            (
                CanonicalSchemaImpactDispositionV1::MigratedV1ToV2,
                3,
                "migrated-v1-to-v2",
            ),
            (
                CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1,
                4,
                "decode-only-legacy-v1",
            ),
            (
                CanonicalSchemaImpactDispositionV1::RetiredV1,
                5,
                "retired-v1",
            ),
            (
                CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame,
                6,
                "inapplicable-no-canonical-frame",
            ),
        ];
        for (disposition, code, name) in schema_impact_oracle {
            assert_eq!(disposition.code(), code);
            assert_eq!(disposition.stable_name(), name);
            assert_eq!(
                CanonicalSchemaImpactDispositionV1::try_from_code(code)
                    .expect("exact schema-impact code"),
                disposition
            );
        }
        assert_eq!(
            CanonicalSchemaImpactDispositionV1::try_from_code(0)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            CanonicalSchemaImpactDispositionV1::try_from_code(7)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let migration_policy_oracle = [
            (
                CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor,
                1,
                "no-schema-predecessor",
            ),
            (
                CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence,
                2,
                "v1-decode-only-compatibility-evidence",
            ),
            (CanonicalSchemaMigrationPolicyV1::V1Retired, 3, "v1-retired"),
        ];
        for (policy, code, name) in migration_policy_oracle {
            assert_eq!(policy.code(), code);
            assert_eq!(policy.stable_name(), name);
            assert_eq!(
                CanonicalSchemaMigrationPolicyV1::try_from_code(code)
                    .expect("exact migration-policy code"),
                policy
            );
        }
        assert_eq!(
            CanonicalSchemaMigrationPolicyV1::try_from_code(0)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            CanonicalSchemaMigrationPolicyV1::try_from_code(4)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let capability_registry =
            BaseCoverageCloseCapabilityRegistryV1::frozen().expect("capability registry");
        let capability_oracle = [
            (
                1,
                "fs-evidence-runner.close.control-input.read",
                BaseCoverageCloseCapabilityPolicyV1::DeclaredControlInputRead,
                "declared-control-input-read-proves-no-content-authenticity-or-execution",
            ),
            (
                2,
                "fs-evidence-runner.close.release-process.control",
                BaseCoverageCloseCapabilityPolicyV1::VersionBoundReleaseProcessControl,
                "declared-process-control-proves-no-launch-drain-success-or-version-match",
            ),
            (
                3,
                "fs-evidence-runner.close.retained-evidence.write",
                BaseCoverageCloseCapabilityPolicyV1::AttemptConfinedRetainedEvidenceWrite,
                "declared-retention-write-proves-no-completeness-durability-or-validity",
            ),
            (
                4,
                "fs-evidence-runner.close.evidence-input.read",
                BaseCoverageCloseCapabilityPolicyV1::DeclaredEvidenceInputRead,
                "declared-evidence-input-read-proves-no-evidence-validity-or-verification",
            ),
            (
                5,
                "fs-evidence-runner.close.publication-output.commit",
                BaseCoverageCloseCapabilityPolicyV1::SelectionBoundPublicationOutputCommit,
                "declared-publication-output-proves-no-write-durability-receipt-or-authority",
            ),
        ];
        assert_eq!(capability_registry.rows().len(), capability_oracle.len());
        let mut capability_roots = BTreeSet::new();
        for (row, (code, stable_id, policy, no_claim)) in
            capability_registry.rows().iter().zip(capability_oracle)
        {
            assert_eq!(row.id().code(), code);
            assert_eq!(row.stable_id().as_str(), stable_id);
            assert_eq!(row.owner().as_str(), "fs-evidence-runner.coverage");
            assert_eq!(row.policy(), policy);
            assert_eq!(row.policy().code(), code);
            assert_eq!(row.no_claim().as_str(), no_claim);
            assert!(capability_roots.insert(row.root()));
        }
        assert_eq!(capability_roots.len(), capability_oracle.len());
        assert_eq!(
            capability_registry
                .descriptor_by_stable_id("fs-evidence-runner.close.evidence-input.read")
                .expect("stable-ID lookup")
                .id()
                .code(),
            4
        );
        assert_eq!(
            capability_registry
                .descriptor_by_policy(
                    BaseCoverageCloseCapabilityPolicyV1::SelectionBoundPublicationOutputCommit,
                )
                .expect("policy lookup")
                .id()
                .code(),
            5
        );
        assert!(
            capability_registry
                .descriptor_by_stable_id("frankensim-epic-foundations-huq.24.4.1.4")
                .is_none(),
            "a Bead owner/route identifier cannot become a capability"
        );
        assert_eq!(
            capability_registry
                .reconstruct_exact(capability_registry.rows(), capability_registry.root())
                .expect("exact capability registry"),
            capability_registry
        );

        let mut duplicate_capability_rows = capability_registry.rows().to_vec();
        duplicate_capability_rows[1] = duplicate_capability_rows[0].clone();
        assert_eq!(
            capability_registry
                .reconstruct_exact(&duplicate_capability_rows, capability_registry.root())
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reordered_capability_rows = capability_registry.rows().to_vec();
        reordered_capability_rows.swap(0, 1);
        assert_eq!(
            capability_registry
                .reconstruct_exact(&reordered_capability_rows, capability_registry.root())
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(
            capability_registry
                .reconstruct_exact(
                    &capability_registry.rows()[..capability_registry.rows().len() - 1],
                    capability_registry.root(),
                )
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Missing
        );
        let mut unknown_capability_rows = capability_registry.rows().to_vec();
        unknown_capability_rows[0].id =
            BaseCoverageCloseCapabilityIdV1::new(6).expect("nonzero unknown ID");
        assert_eq!(
            capability_registry
                .reconstruct_exact(&unknown_capability_rows, capability_registry.root())
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        let moved_descriptor = BaseCoverageCloseCapabilityDescriptorV1::frozen(
            1,
            "fs-evidence-runner.close.control-input.read",
            BaseCoverageCloseCapabilityPolicyV1::DeclaredControlInputRead,
            "declared-control-input-read-proves-no-content-authenticity-or-authority",
        )
        .expect("one-field descriptor mutation");
        assert_ne!(
            moved_descriptor.root(),
            capability_registry.rows()[0].root()
        );

        let capability_profiles =
            BaseCoverageCloseCapabilityProfileRegistryV1::frozen(&capability_registry)
                .expect("capability profile registry");
        let profile_oracle = [
            (
                BaseCoverageCloseCapabilityProfileV1::None,
                "fs-evidence-runner.close-capability.none.v1",
                &[][..],
            ),
            (
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
                "fs-evidence-runner.close-capability.release-control.v1",
                &[1, 2, 3][..],
            ),
            (
                BaseCoverageCloseCapabilityProfileV1::ReleasePublication,
                "fs-evidence-runner.close-capability.release-publication.v1",
                &[1, 2, 3, 5][..],
            ),
            (
                BaseCoverageCloseCapabilityProfileV1::ReleaseVerification,
                "fs-evidence-runner.close-capability.release-verification.v1",
                &[1, 2, 3, 4][..],
            ),
            (
                BaseCoverageCloseCapabilityProfileV1::ReleaseCanonical,
                "fs-evidence-runner.close-capability.release-canonical.v1",
                &[1, 2, 3, 4, 5][..],
            ),
        ];
        assert_eq!(capability_profiles.rows().len(), profile_oracle.len());
        for (row, (profile, stable_id, codes)) in
            capability_profiles.rows().iter().zip(profile_oracle)
        {
            assert_eq!(row.profile(), profile);
            assert_eq!(row.stable_id().as_str(), stable_id);
            let observed_codes = row
                .required()
                .iter()
                .map(|id| id.code())
                .collect::<Vec<_>>();
            assert_eq!(observed_codes.as_slice(), codes);
            assert_eq!(row.required(), row.permitted());
            assert_eq!(
                row.no_claim().as_str(),
                "capability-contract-proves-no-acquisition-effect-success-or-authority"
            );
        }
        assert_eq!(
            capability_profiles
                .descriptor_by_stable_id(
                    "fs-evidence-runner.close-capability.release-verification.v1",
                )
                .expect("stable profile lookup")
                .profile(),
            BaseCoverageCloseCapabilityProfileV1::ReleaseVerification
        );
        assert_eq!(
            capability_profiles.profile_by_stable_id(
                "fs-evidence-runner.close-capability.release-publication.v1",
            ),
            Some(BaseCoverageCloseCapabilityProfileV1::ReleasePublication)
        );
        assert!(
            capability_profiles
                .profile_by_stable_id("frankensim-epic-foundations-huq.24.4.1.4")
                .is_none(),
            "a Bead owner/route identifier cannot become a capability profile"
        );
        assert_ne!(
            capability_registry.root().content_hash(),
            capability_profiles.root().content_hash()
        );
        assert_eq!(
            capability_profiles
                .reconstruct_exact(
                    &capability_registry,
                    capability_profiles.rows(),
                    capability_profiles.root(),
                )
                .expect("exact capability profile registry"),
            capability_profiles
        );
        let mut duplicate_profile_rows = capability_profiles.rows().to_vec();
        duplicate_profile_rows[1] = duplicate_profile_rows[0].clone();
        assert_eq!(
            capability_profiles
                .reconstruct_exact(
                    &capability_registry,
                    &duplicate_profile_rows,
                    capability_profiles.root(),
                )
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reordered_profile_rows = capability_profiles.rows().to_vec();
        reordered_profile_rows.swap(0, 1);
        assert_eq!(
            capability_profiles
                .reconstruct_exact(
                    &capability_registry,
                    &reordered_profile_rows,
                    capability_profiles.root(),
                )
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let mapping_oracle = [
            (
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:publication-state-v2",
                BaseCoverageCloseCapabilityProfileV1::ReleasePublication,
            ),
            (
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:publication-v2",
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
            ),
            (
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:verifier-v2",
                BaseCoverageCloseCapabilityProfileV1::ReleaseVerification,
            ),
            (
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:rjoq-handoff-v1",
                BaseCoverageCloseCapabilityProfileV1::ReleaseVerification,
            ),
            (
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:canonical-runner-v2",
                BaseCoverageCloseCapabilityProfileV1::ReleaseCanonical,
            ),
            (
                BaseCoverageManifestClassV1::ExternalMutation,
                "external-mutation:base-contract-exact-result-join",
                BaseCoverageCloseCapabilityProfileV1::ReleaseCanonical,
            ),
            (
                BaseCoverageManifestClassV1::ExternalGovernance,
                "external-governance:live-source-dependency-closure",
                BaseCoverageCloseCapabilityProfileV1::ReleaseControl,
            ),
        ];
        let mut mapped_profiles = BTreeSet::from([BaseCoverageCloseCapabilityProfileV1::None]);
        for (source_class, source_case_id, expected_profile) in mapping_oracle {
            let observed = base_coverage_close_capability_profile_for_source_case_v1(
                source_class,
                source_case_id,
                BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution,
            )
            .expect("literal downstream capability mapping");
            assert_eq!(observed, expected_profile);
            mapped_profiles.insert(observed);
        }
        for (source_class, source_case_id, scope) in [
            (
                BaseCoverageManifestClassV1::Unit,
                "unit:coverage:local-capability",
                BaseCoverageCloseExecutionScopeV1::CrateTest,
            ),
            (
                BaseCoverageManifestClassV1::CompileFailDoctest,
                "compile-fail:coverage:local-capability",
                BaseCoverageCloseExecutionScopeV1::CompileFailDoctest,
            ),
            (
                BaseCoverageManifestClassV1::ProjectionE2e,
                "projection-e2e:coverage:local-capability",
                BaseCoverageCloseExecutionScopeV1::InProcessProjection,
            ),
            (
                BaseCoverageManifestClassV1::ManifestContract,
                "manifest-contract:coverage:local-capability",
                BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration,
            ),
        ] {
            assert_eq!(
                base_coverage_close_capability_profile_for_source_case_v1(
                    source_class,
                    source_case_id,
                    scope,
                )
                .expect("literal local capability mapping"),
                BaseCoverageCloseCapabilityProfileV1::None
            );
        }
        assert_eq!(
            mapped_profiles,
            BTreeSet::from(BaseCoverageCloseCapabilityProfileV1::ALL)
        );
        assert_eq!(
            base_coverage_close_capability_profile_for_source_case_v1(
                BaseCoverageManifestClassV1::ExternalE2eScript,
                "external-e2e:unknown-capability-route",
                BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let mut contracts = Vec::new();
        for profile in BaseCoverageCloseCapabilityProfileV1::ALL {
            let contract = BaseCoverageCloseCapabilityContractV1::for_profile(
                &capability_registry,
                &capability_profiles,
                profile,
            )
            .expect("exact capability contract");
            assert_eq!(contract.profile(), profile);
            assert_eq!(contract.required(), contract.permitted());
            assert_eq!(contract.registry_root(), capability_registry.root());
            assert_eq!(contract.profile_registry_root(), capability_profiles.root());
            assert_eq!(
                contract.no_claim().as_str(),
                "capability-contract-proves-no-acquisition-effect-success-or-authority"
            );
            contracts.push(contract);
        }
        assert!(contracts[0].required().is_empty());
        assert_eq!(
            contracts
                .iter()
                .map(|contract| contract.root())
                .collect::<BTreeSet<_>>()
                .len(),
            BaseCoverageCloseCapabilityProfileV1::ALL.len()
        );
        let release_control = &contracts[1];
        let reconstructed_release_control = release_control
            .reconstruct_exact(
                &capability_registry,
                &capability_profiles,
                release_control.profile(),
                release_control.required(),
                release_control.permitted(),
                release_control.root(),
            )
            .expect("exact release-control contract");
        assert_eq!(&reconstructed_release_control, release_control);
        let mut reversed_required = release_control.required().to_vec();
        reversed_required.swap(0, 1);
        assert_eq!(
            release_control
                .reconstruct_exact(
                    &capability_registry,
                    &capability_profiles,
                    release_control.profile(),
                    &reversed_required,
                    release_control.permitted(),
                    release_control.root(),
                )
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let unknown_id = BaseCoverageCloseCapabilityIdV1::new(6).expect("nonzero unknown ID");
        assert_eq!(
            release_control
                .reconstruct_exact(
                    &capability_registry,
                    &capability_profiles,
                    release_control.profile(),
                    &[unknown_id],
                    &[unknown_id],
                    release_control.root(),
                )
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let ids = |codes: &[u16]| {
            codes
                .iter()
                .copied()
                .map(BaseCoverageCloseCapabilityIdV1::new)
                .collect::<Result<Vec<_>, _>>()
                .expect("nonzero capability IDs")
        };
        for (label, values) in [
            ("empty", ids(&[])),
            ("one", ids(&[1])),
            ("exact-five", ids(&[1, 2, 3, 4, 5])),
        ] {
            super::validate_close_capability_id_set(
                "coverage.close.capability.boundary_fixture",
                &capability_registry,
                &values,
            )
            .unwrap_or_else(|error| panic!("{label} base capability set must be legal: {error}"));
        }
        assert_eq!(
            super::validate_close_capability_id_set(
                "coverage.close.capability.boundary_fixture",
                &capability_registry,
                &vec![BaseCoverageCloseCapabilityIdV1::new(1).expect("id"); 6],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        assert_eq!(
            BASE_COVERAGE_CLOSE_REGISTERED_EXTENSION_CAPABILITY_MAX_V1,
            64
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(0)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(65)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1::new(
                StableTokenV2::new("fs-evidence-runner.close.control-input.read")
                    .expect("base namespace")
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let empty_extension_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![])
                .expect("explicit empty extension registry");
        let empty_extension_set = BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
            &empty_extension_registry,
            vec![],
        )
        .expect("explicit empty extension set");
        assert!(empty_extension_registry.rows().is_empty());
        assert!(empty_extension_set.values().is_empty());
        assert_eq!(
            empty_extension_set.registry_root(),
            empty_extension_registry.root()
        );

        let one_extension_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(1),
            ])
            .expect("one-row extension registry");
        let one_extension_set = BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
            &one_extension_registry,
            vec![BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(1).expect("extension ID")],
        )
        .expect("one-member extension set");
        assert_eq!(one_extension_registry.rows().len(), 1);
        assert_eq!(one_extension_set.values().len(), 1);

        let extension_rows = (1..=64).map(extension_descriptor).collect::<Vec<_>>();
        let maximum_extension_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(extension_rows.clone())
                .expect("64-row extension registry");
        let maximum_extension_ids = (1..=64)
            .map(|code| {
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(code)
                    .expect("bounded extension ID")
            })
            .collect::<Vec<_>>();
        let maximum_extension_set = BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
            &maximum_extension_registry,
            maximum_extension_ids.clone(),
        )
        .expect("64-member extension set");
        assert_eq!(maximum_extension_registry.rows().len(), 64);
        assert_eq!(maximum_extension_set.values().len(), 64);
        assert_eq!(
            maximum_extension_registry
                .descriptor(
                    BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(64)
                        .expect("maximum extension ID"),
                )
                .expect("maximum registered descriptor")
                .stable_id()
                .as_str(),
            "org.example.fs-evidence-runner.extension.capability-64"
        );

        let mut sixty_five_extension_rows = extension_rows.clone();
        sixty_five_extension_rows.push(extension_rows[63].clone());
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(
                sixty_five_extension_rows,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let mut sixty_five_extension_ids = maximum_extension_ids.clone();
        sixty_five_extension_ids.push(maximum_extension_ids[63]);
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &maximum_extension_registry,
                sixty_five_extension_ids,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut reordered_extension_rows = extension_rows.clone();
        reordered_extension_rows.swap(0, 1);
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(
                reordered_extension_rows,
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let duplicate_namespace_row =
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(2).expect("extension ID"),
                extension_rows[0].stable_id().clone(),
                extension_owner("org.example.fs-evidence-runner.extension-owner"),
                extension_scope("org.example.fs-evidence-runner.extension-scope"),
                extension_no_claim("extension-contract-proves-no-acquisition-effect-or-authority"),
            )
            .expect("descriptor-level namespace is valid before registry uniqueness");
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_rows[0].clone(),
                duplicate_namespace_row,
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let extension_id_one =
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(1).expect("extension ID");
        let extension_id_two =
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(2).expect("extension ID");
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &maximum_extension_registry,
                vec![extension_id_one, extension_id_one],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &maximum_extension_registry,
                vec![extension_id_two, extension_id_one],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &one_extension_registry,
                vec![extension_id_two],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let moved_extension_descriptor =
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                extension_id_one,
                extension_rows[0].stable_id().clone(),
                extension_rows[0].owner().clone(),
                extension_scope("org.example.fs-evidence-runner.changed-extension-scope"),
                extension_rows[0].no_claim().clone(),
            )
            .expect("one-field extension descriptor mutation");
        assert_ne!(
            moved_extension_descriptor.root(),
            extension_rows[0].root(),
            "every extension descriptor field enters its nominal root"
        );
        assert_ne!(
            empty_extension_registry.root(),
            one_extension_registry.root(),
            "extension registry membership enters its nominal root"
        );
        assert_ne!(
            empty_extension_set.root(),
            one_extension_set.root(),
            "extension set membership enters its nominal root"
        );
        assert_ne!(
            maximum_extension_registry.root().content_hash(),
            capability_registry.root().content_hash(),
            "extension and base capability registries remain nominally disjoint"
        );

        let observed = BaseCoverageCloseObservedCapabilitySetsV1::new(
            &capability_registry,
            release_control,
            ids(&[1, 2, 3]),
            ids(&[1, 2, 3]),
            ids(&[1]),
            ids(&[1, 2]),
            ids(&[3]),
        )
        .expect("reconciled observed capability sets");
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::reconstruct_exact(
                &capability_registry,
                release_control,
                observed.required(),
                observed.granted(),
                observed.observed(),
                observed.returned(),
                observed.revoked(),
                observed.root(),
            )
            .expect("exact observed capability reconstruction"),
            observed
        );
        let moved_observed = BaseCoverageCloseObservedCapabilitySetsV1::new(
            &capability_registry,
            release_control,
            ids(&[1, 2, 3]),
            ids(&[1, 2, 3]),
            ids(&[1, 2]),
            ids(&[1, 2]),
            ids(&[3]),
        )
        .expect("one-field observed-set mutation");
        assert_ne!(observed.root(), moved_observed.root());
        let empty_observed = BaseCoverageCloseObservedCapabilitySetsV1::new(
            &capability_registry,
            &contracts[0],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
        .expect("explicit empty local observation");
        assert!(empty_observed.required().is_empty());
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2]),
                ids(&[1, 2, 3]),
                ids(&[1]),
                ids(&[1, 2]),
                ids(&[3]),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2, 3]),
                ids(&[1, 2]),
                ids(&[1]),
                ids(&[1, 2]),
                vec![],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2, 3]),
                ids(&[1, 2, 3]),
                ids(&[1, 4]),
                ids(&[1, 2]),
                ids(&[3]),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2, 3]),
                ids(&[1, 2, 3]),
                ids(&[1]),
                ids(&[1, 2]),
                ids(&[2, 3]),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2, 3]),
                ids(&[1, 2, 3]),
                ids(&[1]),
                ids(&[1]),
                ids(&[3]),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseCoverageCloseObservedCapabilitySetsV1::new(
                &capability_registry,
                release_control,
                ids(&[1, 2, 3]),
                vec![BaseCoverageCloseCapabilityIdV1::new(1).expect("id"); 6],
                vec![],
                vec![],
                vec![],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        assert_eq!(
            [
                BaseCoverageCloseNumericPartitionV1::Inputs.code(),
                BaseCoverageCloseNumericPartitionV1::Grants.code(),
                BaseCoverageCloseNumericPartitionV1::Observations.code(),
            ],
            [1, 2, 3]
        );
        assert_eq!(
            [
                BaseCoverageCloseNumericPartitionV1::Inputs.name(),
                BaseCoverageCloseNumericPartitionV1::Grants.name(),
                BaseCoverageCloseNumericPartitionV1::Observations.name(),
            ],
            [
                "semantic-numeric-inputs",
                "semantic-numeric-grants",
                "expected-numeric-observations",
            ]
        );
        let local = super::frozen_local_close_budget_set().expect("local budget profile");
        let baseline = five_from_template(vec![], vec![], local).expect("baseline");
        let physical = BaseCoverageCloseNumericUnitV1::physical(
            UnitV2::from_parts(1, 1, [0; 7]).expect("dimensionless physical unit"),
        );
        let logical = BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Dimensionless)
            .expect("dimensionless logical unit");
        let value_a = five_from_template(
            vec![numeric("semantic-value", NumericValueV2::U64(1), physical)],
            vec![],
            local,
        )
        .expect("value a");
        let value_b = five_from_template(
            vec![numeric("semantic-value", NumericValueV2::U64(2), physical)],
            vec![],
            local,
        )
        .expect("value b");
        let logical_unit = five_from_template(
            vec![numeric("semantic-value", NumericValueV2::U64(1), logical)],
            vec![],
            local,
        )
        .expect("logical unit");
        let scaled_physical_unit = five_from_template(
            vec![numeric(
                "semantic-value",
                NumericValueV2::U64(1),
                BaseCoverageCloseNumericUnitV1::physical(
                    UnitV2::from_parts(2, 1, [0; 7]).expect("scaled physical unit"),
                ),
            )],
            vec![],
            local,
        )
        .expect("scaled physical unit");
        let dimensioned_physical_unit = five_from_template(
            vec![numeric(
                "semantic-value",
                NumericValueV2::U64(1),
                BaseCoverageCloseNumericUnitV1::physical(
                    UnitV2::from_parts(1, 1, [1, 0, 0, 0, 0, 0, 0])
                        .expect("dimensioned physical unit"),
                ),
            )],
            vec![],
            local,
        )
        .expect("dimensioned physical unit");
        let other_logical_unit = five_from_template(
            vec![numeric(
                "semantic-value",
                NumericValueV2::U64(1),
                BaseCoverageCloseNumericUnitV1::logical(LogicalUnitV2::Count)
                    .expect("count logical unit"),
            )],
            vec![],
            local,
        )
        .expect("other logical unit");
        let grant_value = five_from_template_with_grants(
            vec![],
            vec![numeric("semantic-value", NumericValueV2::U64(1), physical)],
            vec![],
            local,
        )
        .expect("grant-only mutation");
        let observation_value = five_from_template_with_grants(
            vec![],
            vec![],
            vec![numeric("semantic-value", NumericValueV2::U64(1), physical)],
            local,
        )
        .expect("observation-only mutation");
        assert_ne!(baseline.root(), value_a.root());
        assert_ne!(value_a.root(), value_b.root());
        assert_ne!(value_a.root(), logical_unit.root());
        assert_ne!(value_a.root(), scaled_physical_unit.root());
        assert_ne!(value_a.root(), dimensioned_physical_unit.root());
        assert_ne!(logical_unit.root(), other_logical_unit.root());
        assert_ne!(baseline.root(), grant_value.root());
        assert_ne!(baseline.root(), observation_value.root());
        assert_ne!(value_a.root(), grant_value.root());
        assert_ne!(grant_value.root(), observation_value.root());
        assert_ne!(
            value_a.numeric_inputs_root(),
            baseline.numeric_inputs_root()
        );
        assert_eq!(
            value_a.numeric_grants_root(),
            baseline.numeric_grants_root()
        );
        assert_eq!(
            value_a.numeric_observations_root(),
            baseline.numeric_observations_root()
        );
        assert_eq!(
            grant_value.numeric_inputs_root(),
            baseline.numeric_inputs_root()
        );
        assert_ne!(
            grant_value.numeric_grants_root(),
            baseline.numeric_grants_root()
        );
        assert_eq!(
            grant_value.numeric_observations_root(),
            baseline.numeric_observations_root()
        );
        assert_eq!(
            observation_value.numeric_inputs_root(),
            baseline.numeric_inputs_root()
        );
        assert_eq!(
            observation_value.numeric_grants_root(),
            baseline.numeric_grants_root()
        );
        assert_ne!(
            observation_value.numeric_observations_root(),
            baseline.numeric_observations_root()
        );

        let mut budget_rows = local.rows().to_vec();
        budget_rows[0] = budget_row(
            BaseCoverageCloseBudgetAxisV1::Time,
            BaseCoverageCloseBudgetValueV1::U64(59_000_000_000),
            BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
            LogicalUnitV2::Nanoseconds,
        )
        .expect("one-field budget mutation");
        let changed_budget = BaseCoverageCloseBudgetSetV1::new(local.profile(), budget_rows)
            .expect("changed budget");
        let budget_mutant =
            five_from_template(vec![], vec![], changed_budget).expect("budget mutant");
        assert_ne!(baseline.root(), budget_mutant.root());

        let downstream_name = BaseCoverageCloseBudgetSetV1::new(
            BaseCoverageCloseBudgetProfileV1::DownstreamSourceContribution,
            local.rows().to_vec(),
        )
        .expect("same resolved rows under a distinct source profile");
        let profile_mutant =
            five_from_template(vec![], vec![], downstream_name).expect("profile mutant");
        assert_ne!(baseline.root(), profile_mutant.root());

        let registered = LogicalUnitV2::from_tag(16, Some(33)).expect("registered unit syntax");
        let registered_budget = |registry_identity| {
            let mut rows = local.rows().to_vec();
            rows[2] = BaseCoverageCloseTypedBudgetV1::new(
                BaseCoverageCloseBudgetAxisV1::LogicalWork,
                rows[2].hard(),
                rows[2].soft(),
                BaseCoverageCloseLogicalUnitReferenceV1::new(registered, Some(registry_identity))
                    .expect("registry-bound work unit"),
            )
            .expect("registered logical-work row");
            BaseCoverageCloseBudgetSetV1::new(local.profile(), rows)
                .expect("registered logical-work profile")
        };
        let registry_a =
            five_from_template(vec![], vec![], registered_budget(root("budget-registry-a")))
                .expect("budget registry a");
        let registry_b =
            five_from_template(vec![], vec![], registered_budget(root("budget-registry-b")))
                .expect("budget registry b");
        assert_ne!(
            registry_a.root(),
            registry_b.root(),
            "logical-work registry identity enters the Five Explicits root"
        );
    }

    #[test]
    fn registered_extension_capability_ids_and_nominal_descriptor_roles_are_exact() {
        assert_eq!(
            BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            44
        );
        assert_eq!(BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1, 47);
        assert_eq!(
            base_coverage_close_nominal_root_descriptors_v1()
                [BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1..]
                .iter()
                .map(|descriptor| descriptor.schema_name())
                .collect::<Vec<_>>(),
            [
                "registered-extension-capability-descriptor",
                "registered-extension-capability-registry",
                "registered-extension-capability-set",
            ]
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(0)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(65)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let base_registry =
            BaseCoverageCloseCapabilityRegistryV1::frozen().expect("exact base registry");
        for row in base_registry.rows() {
            assert_eq!(
                BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1::new(
                    row.stable_id().clone()
                )
                .unwrap_err()
                .kind(),
                ConstructionErrorKindV2::Incompatible,
                "{} must remain reserved to the base registry",
                row.stable_id().as_str()
            );
        }
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityStableIdV1::new(
                StableTokenV2::new("unnamespaced").expect("stable token")
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityOwnerV1::new(
                StableTokenV2::new("unnamespaced").expect("stable token")
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityScopeV1::new(
                StableTokenV2::new("unnamespaced").expect("stable token")
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let duplicated_role = "org.example.extension.duplicated-role";
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(1).expect("extension ID"),
                extension_stable_id(duplicated_role),
                extension_owner(duplicated_role),
                extension_scope("org.example.extension.scope"),
                extension_no_claim("extension-proves-no-authority"),
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let maximum_stable_id = format!("a.{}", "b".repeat(126));
        let maximum_owner = format!("c.{}", "d".repeat(126));
        let maximum_scope = format!("e.{}", "f".repeat(126));
        let maximum_no_claim = "g".repeat(128);
        for value in [
            maximum_stable_id.as_str(),
            maximum_owner.as_str(),
            maximum_scope.as_str(),
            maximum_no_claim.as_str(),
        ] {
            assert_eq!(value.len(), 128);
        }
        let maximum = BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(1).expect("extension ID"),
            extension_stable_id(maximum_stable_id),
            extension_owner(maximum_owner),
            extension_scope(maximum_scope),
            extension_no_claim(maximum_no_claim),
        )
        .expect("maximum-width descriptor");
        assert_eq!(maximum.stable_id().as_str().len(), 128);
        assert_eq!(maximum.owner().as_str().len(), 128);
        assert_eq!(maximum.scope().as_str().len(), 128);
        assert_eq!(maximum.no_claim().as_str().len(), 128);

        let baseline = extension_descriptor(1);
        let mutations = [
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(2).expect("extension ID"),
                baseline.stable_id().clone(),
                baseline.owner().clone(),
                baseline.scope().clone(),
                baseline.no_claim().clone(),
            )
            .expect("ID mutation"),
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                baseline.id(),
                extension_stable_id("org.example.extension.changed-capability"),
                baseline.owner().clone(),
                baseline.scope().clone(),
                baseline.no_claim().clone(),
            )
            .expect("stable-ID mutation"),
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                baseline.id(),
                baseline.stable_id().clone(),
                extension_owner("org.example.extension.changed-owner"),
                baseline.scope().clone(),
                baseline.no_claim().clone(),
            )
            .expect("owner mutation"),
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                baseline.id(),
                baseline.stable_id().clone(),
                baseline.owner().clone(),
                extension_scope("org.example.extension.changed-scope"),
                baseline.no_claim().clone(),
            )
            .expect("scope mutation"),
            BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
                baseline.id(),
                baseline.stable_id().clone(),
                baseline.owner().clone(),
                baseline.scope().clone(),
                extension_no_claim("changed-extension-proves-no-authority"),
            )
            .expect("no-claim mutation"),
        ];
        for mutation in mutations {
            assert_ne!(
                mutation.root(),
                baseline.root(),
                "each independently nominal descriptor field must move the root"
            );
        }
    }

    #[test]
    fn registered_extension_capability_registry_boundaries_and_diagnostics_are_exact() {
        let empty = BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![])
            .expect("exact-empty registry");
        let one = BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
            extension_descriptor(1),
        ])
        .expect("one-row registry");
        let rows = (1..=64).map(extension_descriptor).collect::<Vec<_>>();
        let maximum = BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(rows.clone())
            .expect("64-row registry");
        assert!(empty.rows().is_empty());
        assert_eq!(one.rows().len(), 1);
        assert_eq!(maximum.rows().len(), 64);
        assert_eq!(
            maximum
                .descriptor(
                    BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(64)
                        .expect("maximum extension ID")
                )
                .expect("maximum descriptor")
                .id()
                .code(),
            64
        );

        let mut one_over = rows.clone();
        one_over.push(rows[63].clone());
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(one_over)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(2)
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(1),
                extension_descriptor(3),
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(2),
                extension_descriptor(1),
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(1),
                extension_descriptor(1),
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                extension_descriptor(2),
                extension_descriptor(2),
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate,
            "duplicate diagnosis precedes simultaneous missing/order faults"
        );

        let first = extension_descriptor(1);
        let duplicate_namespace = BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
            BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(2).expect("extension ID"),
            first.stable_id().clone(),
            extension_owner("org.example.extension.second-owner"),
            extension_scope("org.example.extension.second-scope"),
            extension_no_claim("second-extension-proves-no-authority"),
        )
        .expect("descriptor-local namespace remains valid");
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![
                first.clone(),
                duplicate_namespace,
            ])
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let changed_scope = BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
            first.id(),
            first.stable_id().clone(),
            first.owner().clone(),
            extension_scope("org.example.extension.changed-registry-scope"),
            first.no_claim().clone(),
        )
        .expect("same-cardinality row mutation");
        let original_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![first])
                .expect("original one-row registry");
        let mutated_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![changed_scope])
                .expect("mutated one-row registry");
        assert_ne!(original_registry.root(), mutated_registry.root());
        assert_ne!(
            maximum.root().content_hash(),
            BaseCoverageCloseCapabilityRegistryV1::frozen()
                .expect("base registry")
                .root()
                .content_hash()
        );
    }

    #[test]
    fn registered_extension_capability_sets_boundaries_and_diagnostics_are_exact() {
        let rows = (1..=64).map(extension_descriptor).collect::<Vec<_>>();
        let registry = BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(rows.clone())
            .expect("64-row registry");
        let ids = (1..=64)
            .map(|code| {
                BaseCoverageCloseRegisteredExtensionCapabilityIdV1::new(code)
                    .expect("bounded extension ID")
            })
            .collect::<Vec<_>>();
        let empty = BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(&registry, vec![])
            .expect("exact-empty set");
        let one = BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(&registry, vec![ids[0]])
            .expect("one-member set");
        let maximum =
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(&registry, ids.clone())
                .expect("64-member set");
        assert!(empty.values().is_empty());
        assert_eq!(one.values(), &[ids[0]]);
        assert_eq!(maximum.values().len(), 64);
        assert_ne!(
            empty.root(),
            one.root(),
            "membership must move the root under one unchanged registry"
        );

        let mut one_over = ids.clone();
        one_over.push(ids[63]);
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(&registry, one_over)
                .unwrap_err()
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &registry,
                vec![ids[0], ids[0]],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &registry,
                vec![ids[0], ids[1], ids[0]],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::Duplicate,
            "a non-adjacent duplicate must not degrade into OutOfOrder"
        );
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &registry,
                vec![ids[1], ids[0]],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let one_row_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![rows[0].clone()])
                .expect("one-row registry");
        assert_eq!(
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &one_row_registry,
                vec![ids[1]],
            )
            .unwrap_err()
            .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let changed_scope = BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1::new(
            rows[0].id(),
            rows[0].stable_id().clone(),
            rows[0].owner().clone(),
            extension_scope("org.example.extension.changed-set-registry-scope"),
            rows[0].no_claim().clone(),
        )
        .expect("registry mutation");
        let changed_registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![changed_scope])
                .expect("changed registry");
        let same_member_changed_registry =
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(
                &changed_registry,
                vec![ids[0]],
            )
            .expect("same member under changed registry");
        assert_ne!(one.root(), same_member_changed_registry.root());
    }

    #[test]
    fn registered_extension_capability_canonical_roots_and_magics_are_exact() {
        fn descriptor_oracle(
            magic: &'static [u8],
            row: &BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1,
        ) -> ContentHash {
            let mut frame = super::CanonicalFrameV1::new(magic, 2 * 1024)
                .expect("independent descriptor frame");
            frame.push_u16("api", 2).expect("api");
            frame.push_u16("wire", 1).expect("wire");
            frame
                .push_str("predecessor", "no-predecessor")
                .expect("predecessor");
            frame.push_u16("id", row.id().code()).expect("id");
            frame
                .push_str("stable_id", row.stable_id().as_str())
                .expect("stable ID");
            frame
                .push_str("owner", row.owner().as_str())
                .expect("owner");
            frame
                .push_str("scope", row.scope().as_str())
                .expect("scope");
            frame
                .push_str("no_claim", row.no_claim().as_str())
                .expect("no-claim");
            frame.root(
                "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-descriptor.v1",
            )
        }

        fn registry_oracle(
            magic: &'static [u8],
            rows: &[BaseCoverageCloseRegisteredExtensionCapabilityDescriptorV1],
        ) -> ContentHash {
            let mut frame =
                super::CanonicalFrameV1::new(magic, 8 * 1024).expect("independent registry frame");
            frame.push_u16("api", 2).expect("api");
            frame.push_u16("wire", 1).expect("wire");
            frame
                .push_str("predecessor", "no-predecessor")
                .expect("predecessor");
            frame
                .push_u16("count", u16::try_from(rows.len()).expect("bounded count"))
                .expect("count");
            for row in rows {
                frame
                    .push_bytes("descriptor_root", row.root().content_hash().as_bytes())
                    .expect("descriptor root");
            }
            frame
                .push_str(
                    "no_claim",
                    "registered-extension-capability-registry-root-proves-declared-extension-membership-not-base-membership-acquisition-or-authority",
                )
                .expect("no-claim");
            frame.root(
                "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-registry.v1",
            )
        }

        fn set_oracle(
            magic: &'static [u8],
            registry: &BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1,
            values: &[BaseCoverageCloseRegisteredExtensionCapabilityIdV1],
        ) -> ContentHash {
            let mut frame =
                super::CanonicalFrameV1::new(magic, 2 * 1024).expect("independent set frame");
            frame.push_u16("api", 2).expect("api");
            frame.push_u16("wire", 1).expect("wire");
            frame
                .push_str("predecessor", "no-predecessor")
                .expect("predecessor");
            frame
                .push_bytes("registry_root", registry.root().content_hash().as_bytes())
                .expect("registry root");
            frame
                .push_u16("count", u16::try_from(values.len()).expect("bounded count"))
                .expect("count");
            for value in values {
                frame.push_u16("value", value.code()).expect("value");
            }
            frame
                .push_str(
                    "no_claim",
                    "registered-extension-capability-set-root-proves-structural-extension-membership-not-base-membership-acquisition-or-authority",
                )
                .expect("no-claim");
            frame.root(
                "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-set.v1",
            )
        }

        const DESCRIPTOR_MAGIC: &[u8] = b"FSCLOSEEXTCAPDESC\x01";
        const REGISTRY_MAGIC: &[u8] = b"FSCLOSEEXTCAPREG\x01";
        const SET_MAGIC: &[u8] = b"FSCLOSEEXTCAPSET\x01";
        assert_eq!(DESCRIPTOR_MAGIC.last(), Some(&1));
        assert_eq!(REGISTRY_MAGIC.last(), Some(&1));
        assert_eq!(SET_MAGIC.last(), Some(&1));

        let row = extension_descriptor(1);
        let registry =
            BaseCoverageCloseRegisteredExtensionCapabilityRegistryV1::new(vec![row.clone()])
                .expect("one-row registry");
        let values = [row.id()];
        let set =
            BaseCoverageCloseRegisteredExtensionCapabilitySetV1::new(&registry, values.to_vec())
                .expect("one-member set");

        assert_eq!(
            row.root().content_hash(),
            descriptor_oracle(DESCRIPTOR_MAGIC, &row)
        );
        assert_eq!(
            registry.root().content_hash(),
            registry_oracle(REGISTRY_MAGIC, registry.rows())
        );
        assert_eq!(
            set.root().content_hash(),
            set_oracle(SET_MAGIC, &registry, &values)
        );

        for mutation in [
            b"FSCLOSEEXTCAPDESC\x02".as_slice(),
            b"FSCLOSEEXTCAPDESC".as_slice(),
            b"FSCLOSEEXTCAPDESC\x01x".as_slice(),
        ] {
            assert_ne!(row.root().content_hash(), descriptor_oracle(mutation, &row));
        }
        for mutation in [
            b"FSCLOSEEXTCAPREG\x02".as_slice(),
            b"FSCLOSEEXTCAPREG".as_slice(),
            b"FSCLOSEEXTCAPREG\x01x".as_slice(),
        ] {
            assert_ne!(
                registry.root().content_hash(),
                registry_oracle(mutation, registry.rows())
            );
        }
        for mutation in [
            b"FSCLOSEEXTCAPSET\x02".as_slice(),
            b"FSCLOSEEXTCAPSET".as_slice(),
            b"FSCLOSEEXTCAPSET\x01x".as_slice(),
        ] {
            assert_ne!(
                set.root().content_hash(),
                set_oracle(mutation, &registry, &values)
            );
        }
    }

    #[test]
    fn five_explicits_maximum_numeric_frame_is_feasible() {
        let physical = BaseCoverageCloseNumericUnitV1::physical(
            UnitV2::from_parts(i128::MAX, u128::MAX, [i16::MIN, i16::MAX, -1, 0, 1, 2, 3])
                .expect("maximum-shape physical unit"),
        );
        let numeric_value =
            NumericValueV2::Rational(RationalV2::new(i128::MAX, u128::MAX).expect("rational"));
        let maximum_partition = (0..64)
            .map(|index| {
                let name = format!("n{index:02}{}", "a".repeat(125));
                assert_eq!(name.len(), 128);
                numeric(name, numeric_value.clone(), physical)
            })
            .collect::<Vec<_>>();
        let budgets = super::frozen_local_close_budget_set().expect("local budget profile");
        let first = five_from_template_with_grants(
            maximum_partition.clone(),
            maximum_partition.clone(),
            maximum_partition.clone(),
            budgets,
        )
        .expect("maximum canonical three-profile Five Explicits frame");
        let second = five_from_template_with_grants(
            maximum_partition.clone(),
            maximum_partition.clone(),
            maximum_partition,
            budgets,
        )
        .expect("deterministic maximum three-profile frame replay");
        assert_eq!(first.numeric_inputs().len(), 64);
        assert_eq!(first.numeric_grants().len(), 64);
        assert_eq!(first.numeric_observations().len(), 64);
        assert_ne!(first.numeric_inputs_root(), first.numeric_grants_root());
        assert_ne!(
            first.numeric_inputs_root(),
            first.numeric_observations_root()
        );
        assert_ne!(
            first.numeric_grants_root(),
            first.numeric_observations_root()
        );
        assert_eq!(first.root(), second.root());
    }

    #[test]
    fn race_facet_is_registered_inapplicable_for_pure_single_threaded_validator() {
        assert_registered_inapplicable_facet(
            BaseCoverageCloseFacetV1::Race,
            BaseCoverageCloseReasonCodeV1::RaceNotApplicablePureSingleThreadedValidator,
            "race-not-applicable-pure-single-threaded-validator",
        );
    }

    #[test]
    fn trait_facet_is_registered_inapplicable_without_public_trait_contract() {
        assert_registered_inapplicable_facet(
            BaseCoverageCloseFacetV1::Trait,
            BaseCoverageCloseReasonCodeV1::TraitNotApplicableNoPublicTraitContract,
            "trait-not-applicable-no-public-trait-contract",
        );
    }

    #[test]
    fn cancellation_facet_is_registered_inapplicable_for_pure_bounded_validator() {
        assert_registered_inapplicable_facet(
            BaseCoverageCloseFacetV1::Cancellation,
            BaseCoverageCloseReasonCodeV1::CancellationNotApplicablePureBoundedValidator,
            "cancellation-not-applicable-pure-bounded-validator",
        );
    }

    #[test]
    fn release_built_no_mock_e2e_facet_is_registered_as_downstream_owned() {
        assert_registered_inapplicable_facet(
            BaseCoverageCloseFacetV1::ReleaseBuiltNoMockE2e,
            BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned,
            "release-execution-downstream-owned",
        );
    }

    fn assert_registered_inapplicable_facet(
        facet: BaseCoverageCloseFacetV1,
        reason: BaseCoverageCloseReasonCodeV1,
        exact_name: &str,
    ) {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("full close manifest");
        let cells = manifest
            .cells()
            .iter()
            .filter(|cell| cell.facet() == facet)
            .collect::<Vec<_>>();
        assert_eq!(cells.len(), 1);
        assert_eq!(manifest.applicable_facet_count(facet), 0);
        assert_eq!(
            cells[0].partition(),
            BaseCoverageClosePartitionV1::Inapplicable
        );
        assert_eq!(cells[0].expected_reason(), Some(reason));
        assert_eq!(
            cells[0].execution_scope(),
            BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration
        );
        assert_eq!(reason.descriptor().name(), exact_name);
        assert!(!reason.descriptor().owner().is_empty());
        assert!(!reason.descriptor().prerequisite().is_empty());
        assert!(!reason.descriptor().no_claim().is_empty());
    }
}
