//! Stage-A declaration and pure evaluator for Runner V2 base values.
//!
//! The declaration is canonical source-authoritative contract data. The
//! evaluator is deliberately different: it produces one fresh, complete,
//! bounded, rootless handoff and never receives an expected result, subset,
//! callback, capability, route, or attempt identity from its caller.

use core::{cmp::Ordering, fmt};

use super::super::handoff::{
    RunnerV2LocalWorkPackageHandoffV1, RunnerV2RawCellObservationV1, RunnerV2RawDiagnosticV1,
    RunnerV2RawOutcomeKindV1, RunnerV2RawReasonV1, RunnerV2RawRepairV1,
    RunnerV2SafeNumericObservationV1,
};
use crate::canonical::{CanonicalFrameSinkV1, CanonicalFrameV1};
use crate::catalog::{
    DiagnosticCodeV2, DigestRoleV2, LogicalUnitV2, RUNNER_SPEC_V2_API_GENERATION,
    RUNNER_V2_PREDECESSOR_POLICY, RUNNER_V2_WIRE_VERSION, RepairActionKindV2, RetryabilityV2,
    RunProfileV2, RunnerApiGeneration, RunnerWireVersion, WirePredecessorPolicyV1,
};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::coverage::{
    BaseCoverageCloseBudgetAxisV1, BaseCoverageCloseBudgetProfileV1, BaseCoverageCloseBudgetSetV1,
    BaseCoverageCloseBudgetValueV1, BaseCoverageCloseCapabilityContractV1,
    BaseCoverageCloseCapabilityProfileRegistryV1, BaseCoverageCloseCapabilityProfileV1,
    BaseCoverageCloseCapabilityRegistryV1, BaseCoverageCloseLogicalUnitReferenceV1,
    BaseCoverageCloseProfileV1, BaseCoverageCloseSeedExplicitV1, BaseCoverageCloseTargetV1,
    BaseCoverageCloseTypedBudgetV1, CanonicalSchemaImpactDispositionV1,
    CanonicalSchemaMigrationPolicyV1,
};
use crate::identity::{
    BuildIdentityRootV2, DigestValueV2, SourceIdentityRootV2, ToolchainIdentityRootV2,
};
use crate::limits::{
    RUNNER_LIMIT_FIELD_COUNT_V2, RunnerFamilyLimitRequirementsV2, RunnerLimitExpectationV2,
    RunnerLimitFieldV2, RunnerLimitMinimumRuleV2, RunnerLimitTightenabilityV2, RunnerLimitUnitV2,
    RunnerLimitValueV2, RunnerLimitWidthV2, RunnerLimitsCandidateV2, RunnerLimitsV2,
    RunnerLimitsViolationKindV2, RunnerLimitsViolationV2, checked_lifecycle_record_requirement,
};
use crate::path::LogicalBundlePathV1;
use crate::value::{F32BitsV2, F64BitsV2, SeedInapplicableCodeV1, StableTokenV2, TypedOptionV1};
use fs_blake3::{ContentHash, hash_domain};

/// Stable package identity for the first foundational Runner V2 work package.
pub const RUNNER_V2_BASE_VALUES_PACKAGE_ID_V1: &str = "runner-v2.work-package.24-1-1-1-1.v1";
/// Frozen future public entry point; work package `.7` implements the wrapper.
pub const RUNNER_V2_BASE_VALUES_PUBLIC_ENTRY_POINT_V1: &str =
    "fs_evidence_runner::runner_v2::work_packages::run_24_1_1_1_1_v1";
/// Exact local route declaration owned by this work package.
pub const RUNNER_V2_BASE_VALUES_LOCAL_ROUTE_ID_V1: &str =
    "runner-v2.route.24-1-1-1-1.local.work-package.v1";
/// Exact LocalInProcess route count for this work package.
pub const RUNNER_V2_BASE_VALUES_LOCAL_IN_PROCESS_ROUTE_COUNT_V1: usize = 1;
/// Exact ExecutionOwned route count for this Stage-A child.
pub const RUNNER_V2_BASE_VALUES_EXECUTION_OWNED_ROUTE_COUNT_V1: usize = 0;
/// Exact ContributionOnly route count for this Stage-A child.
pub const RUNNER_V2_BASE_VALUES_CONTRIBUTION_ONLY_ROUTE_COUNT_V1: usize = 0;
/// Exact Inapplicable route count for this Stage-A child.
pub const RUNNER_V2_BASE_VALUES_INAPPLICABLE_ROUTE_COUNT_V1: usize = 0;
/// Source-owned no-claim attached to every Stage-A component.
pub const RUNNER_V2_BASE_VALUES_NO_CLAIM_V1: &str =
    "stage-a-proves-no-runtime-observation-retention-authority-or-release-e2e";
/// Exact number of boundary kinds for every limit field.
pub const RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1: usize = 12;
/// Exact number of executable cases in the Stage-A admission fixture.
pub const RUNNER_V2_LIMIT_FIXTURE_CASE_COUNT_V1: usize = 1;
/// Checked lifecycle minimum implied by the frozen one-case, zero-row fixture.
pub const RUNNER_V2_LIMIT_FIXTURE_LIFECYCLE_MINIMUM_V1: u32 = 5;
/// Exact limit-cell cardinality: 71 fields by 12 boundary kinds.
pub const RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1: usize =
    RUNNER_LIMIT_FIELD_COUNT_V2 * RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1;
/// Exact number of non-limit declaration/evaluator cells.
pub const RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1: usize = 15;
/// Exact complete raw-cell cardinality for this package.
pub const RUNNER_V2_BASE_VALUES_CELL_COUNT_V1: usize =
    RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1 + RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1;
/// Exact number of retained pre-Runner-V2 domain obligations.
pub const RUNNER_V2_RETAINED_DOMAIN_OBLIGATION_COUNT_V1: usize = 50;
/// Exact number of opposite-width auxiliary limit mutations.
pub const RUNNER_V2_LIMIT_MUTATION_OBLIGATION_COUNT_V1: usize = 71;
/// Exact Refused-compatible raw-reason contract cardinality.
pub const RUNNER_V2_BASE_VALUES_REFUSAL_REASON_COUNT_V1: usize = 12;
/// Exact number of deferred common-contract requirements.
pub const RUNNER_V2_COMMON_REQUIREMENT_COUNT_V1: usize = 31;
/// Existing broad source inventory cardinality before Runner V2 work packages.
pub const RUNNER_V2_EXISTING_SOURCE_COUNT_V1: usize = 27;
/// Exact number of frozen future source additions.
pub const RUNNER_V2_FUTURE_SOURCE_COUNT_V1: usize = 13;
/// Exact eventual broad source inventory cardinality.
pub const RUNNER_V2_FINAL_SOURCE_COUNT_V1: usize =
    RUNNER_V2_EXISTING_SOURCE_COUNT_V1 + RUNNER_V2_FUTURE_SOURCE_COUNT_V1;
/// Exact number of source members realized by this child.
pub const RUNNER_V2_BASE_VALUES_OWNER_SOURCE_COUNT_V1: usize = 2;
/// Exact current source dependency closure for Stage-A semantics.
pub const RUNNER_V2_BASE_VALUES_DEPENDENCY_SOURCE_COUNT_V1: usize = 16;
/// Exact number of canonical Stage-A schemas deferred to work package `.3`.
pub const RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1: usize = 42;
/// Exact number of separately classified rootless Stage-A semantic types.
pub const RUNNER_V2_BASE_VALUES_ROOTLESS_SCHEMA_COUNT_V1: usize = 1;
/// Exact owned schema inventory: 42 canonical schemas plus one rootless handoff.
pub const RUNNER_V2_BASE_VALUES_OWNED_SCHEMA_COUNT_V1: usize =
    RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1
        + RUNNER_V2_BASE_VALUES_ROOTLESS_SCHEMA_COUNT_V1;

/// One source-declared canonical schema name deferred to the `.3` registry owner.
///
/// This nominal wrapper is deliberately distinct from
/// [`RunnerV2RootlessHandoffSchemaNameV1`]. A rootless semantic type cannot be
/// passed to an API that accepts a canonical schema name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2CanonicalSchemaNameV1(StableTokenV2);

impl RunnerV2CanonicalSchemaNameV1 {
    fn from_source_literal(
        field: &'static str,
        name: &'static str,
    ) -> Result<Self, ConstructionErrorV2> {
        Ok(Self(stage_a_token(field, name)?))
    }

    /// Exact stable canonical schema name.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}

/// The sole source-declared rootless Stage-A semantic type.
///
/// This wrapper has no conversion to [`RunnerV2CanonicalSchemaNameV1`], a
/// canonical frame, or a nominal root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RootlessHandoffSchemaNameV1(StableTokenV2);

impl RunnerV2RootlessHandoffSchemaNameV1 {
    fn from_source_literal(
        field: &'static str,
        name: &'static str,
    ) -> Result<Self, ConstructionErrorV2> {
        Ok(Self(stage_a_token(field, name)?))
    }

    /// Exact stable rootless semantic-type name.
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.0.as_str()
    }
}

const _: () = {
    #[allow(
        dead_code,
        reason = "distinct impls are an always-compiled nominal-type separation witness"
    )]
    trait RunnerV2StageASchemaNameSeparationWitnessV1 {}

    impl RunnerV2StageASchemaNameSeparationWitnessV1 for RunnerV2CanonicalSchemaNameV1 {}
    impl RunnerV2StageASchemaNameSeparationWitnessV1 for RunnerV2RootlessHandoffSchemaNameV1 {}
};

/// One ranked, non-executable repair for a Stage-A inventory mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunnerV2StageAInventoryRepairV1 {
    rank: u8,
    kind: RepairActionKindV2,
    target: &'static str,
}

impl RunnerV2StageAInventoryRepairV1 {
    /// One-based contiguous rank.
    #[must_use]
    pub const fn rank(self) -> u8 {
        self.rank
    }

    /// Closed non-executable repair class.
    #[must_use]
    pub const fn kind(self) -> RepairActionKindV2 {
        self.kind
    }

    /// Stable semantic repair target.
    #[must_use]
    pub const fn target(self) -> &'static str {
        self.target
    }
}

/// Actionable first-divergence report for one exact Stage-A inventory.
///
/// This source-validation diagnostic is separate from the rootless evaluator
/// handoff. It carries no runtime actual, attempt, authority, executable
/// command, physical locator, or retained-log claim. Unregistered observed
/// text is replaced with a fixed placeholder before it reaches this value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAInventoryMismatchV1 {
    kind: ConstructionErrorKindV2,
    inventory: &'static str,
    first_mismatch_index0: usize,
    expected_ordinal1: usize,
    expected_identity: String,
    observed_safe_identity: String,
    observed_identity_redacted: bool,
    component: &'static str,
    expected_semantic_value: String,
    observed_safe_value: String,
    observed_value_redacted: bool,
    semantic_owner: &'static str,
    expected_count: usize,
    observed_count: usize,
    repairs: [RunnerV2StageAInventoryRepairV1; 1],
}

impl RunnerV2StageAInventoryMismatchV1 {
    /// Closed mismatch class.
    #[must_use]
    pub const fn kind(&self) -> ConstructionErrorKindV2 {
        self.kind
    }

    /// Stable inventory identity.
    #[must_use]
    pub const fn inventory(&self) -> &'static str {
        self.inventory
    }

    /// Exact zero-based machine join key for the first divergence.
    #[must_use]
    pub const fn first_mismatch_index0(&self) -> usize {
        self.first_mismatch_index0
    }

    /// Exact one-based ordinal intended for human-facing diagnostics.
    #[must_use]
    pub const fn expected_ordinal1(&self) -> usize {
        self.expected_ordinal1
    }

    /// Exact expected stable identity or the explicit end sentinel.
    #[must_use]
    pub fn expected_identity(&self) -> &str {
        &self.expected_identity
    }

    /// Safe observed stable identity, missing sentinel, or redaction marker.
    #[must_use]
    pub fn observed_safe_identity(&self) -> &str {
        &self.observed_safe_identity
    }

    /// Whether unregistered observed input was redacted.
    #[must_use]
    pub const fn observed_identity_redacted(&self) -> bool {
        self.observed_identity_redacted
    }

    /// First row component whose semantic value diverged.
    #[must_use]
    pub const fn component(&self) -> &'static str {
        self.component
    }

    /// Exact expected semantic value for the first divergent component.
    #[must_use]
    pub fn expected_semantic_value(&self) -> &str {
        &self.expected_semantic_value
    }

    /// Safe observed semantic value or a fixed redaction marker.
    #[must_use]
    pub fn observed_safe_value(&self) -> &str {
        &self.observed_safe_value
    }

    /// Whether the observed semantic value was redacted.
    #[must_use]
    pub const fn observed_value_redacted(&self) -> bool {
        self.observed_value_redacted
    }

    /// Source-authoritative semantic owner.
    #[must_use]
    pub const fn semantic_owner(&self) -> &'static str {
        self.semantic_owner
    }

    /// Exact expected inventory cardinality.
    #[must_use]
    pub const fn expected_count(&self) -> usize {
        self.expected_count
    }

    /// Exact observed inventory cardinality.
    #[must_use]
    pub const fn observed_count(&self) -> usize {
        self.observed_count
    }

    /// Complete ranked non-executable remediation list.
    #[must_use]
    pub const fn repairs(&self) -> &[RunnerV2StageAInventoryRepairV1] {
        &self.repairs
    }

    fn as_construction_error(&self) -> ConstructionErrorV2 {
        stage_a_error(
            self.kind,
            self.inventory,
            "the exact Stage-A inventory at the reported zero-based first divergence",
            self.first_mismatch_index0,
        )
    }
}

impl fmt::Display for RunnerV2StageAInventoryMismatchV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{} mismatch at index {} (ordinal {}), component {}: identity expected {}; identity observed {}; semantic expected {}; semantic observed {}; owner {}; count {}/{}; repair 1 {}",
            self.inventory,
            self.first_mismatch_index0,
            self.expected_ordinal1,
            self.component,
            self.expected_identity,
            self.observed_safe_identity,
            self.expected_semantic_value,
            self.observed_safe_value,
            self.semantic_owner,
            self.observed_count,
            self.expected_count,
            self.repairs[0].target,
        )
    }
}

impl std::error::Error for RunnerV2StageAInventoryMismatchV1 {}

const STAGE_A_FRAME_MAX_BYTES_V1: usize = 1024 * 1024;
const KIB: u64 = 1024;
const MIB: u64 = 1024 * KIB;
const STAGE_A_CANONICAL_SCHEMA_NAMES_V1: [&str; RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1] = [
    "runner-v2-stage-a-declaration-root-v1",
    "runner-v2-stage-a-oracle-root-v1",
    "runner-v2-stage-a-case-manifest-root-v1",
    "runner-v2-stage-a-schema-inventory-root-v1",
    "runner-v2-stage-a-feature-declaration-root-v1",
    "runner-v2-stage-a-five-explicits-root-v1",
    "runner-v2-stage-a-source-member-root-v1",
    "runner-v2-limit-boundary-kind-v1",
    "runner-v2-limit-fixture-declaration-v1",
    "runner-v2-limit-companion-normalization-v1",
    "runner-v2-stage-a-expected-partition-v1",
    "runner-v2-stage-a-cell-group-v1",
    "runner-v2-retained-domain-facet-v1",
    "runner-v2-retained-domain-obligation-v1",
    "runner-v2-limit-mutation-obligation-v1",
    "runner-v2-stage-a-version-requirements-v1",
    "runner-v2-stage-a-five-explicits-v1",
    "runner-v2-contract-plane-set-v1",
    "runner-v2-common-fulfillment-stage-v1",
    "runner-v2-unavailable-common-root-v1",
    "runner-v2-common-contract-requirement-v1",
    "runner-v2-future-source-requirement-v1",
    "runner-v2-owner-source-member-v1",
    "runner-v2-dependency-source-member-v1",
    "runner-v2-schema-impact-deferral-v1",
    "runner-v2-rootless-ac58-fragment-v1",
    "runner-v2-local-route-class-v1",
    "runner-v2-local-route-declaration-v1",
    "runner-v2-stage-a-inapplicability-declaration-v1",
    "runner-v2-stage-a-oracle-row-v1",
    "runner-v2-stage-a-projection-row-v1",
    "runner-v2-stage-a-meta-operation-v1",
    "runner-v2-stage-a-cell-operation-v1",
    "runner-v2-stage-a-cell-declaration-v1",
    "runner-v2-base-values-stage-a-declaration-v1",
    "runner-v2-raw-outcome-reason-contract-v1",
    "runner-v2-safe-numeric-value-v1",
    "runner-v2-safe-numeric-unit-v1",
    "runner-v2-safe-numeric-observation-v1",
    "runner-v2-raw-repair-v1",
    "runner-v2-raw-diagnostic-v1",
    "runner-v2-raw-cell-observation-v1",
];
const STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1: &str = "runner-v2-local-work-package-handoff-v1";

const _: [(); 42] = [(); STAGE_A_CANONICAL_SCHEMA_NAMES_V1.len()];
const _: [(); 1] = [(); RUNNER_V2_BASE_VALUES_ROOTLESS_SCHEMA_COUNT_V1];
const _: [(); 43] = [(); RUNNER_V2_BASE_VALUES_OWNED_SCHEMA_COUNT_V1];
const _: [(); 12] = [(); RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1];
const _: [(); 71] = [(); RUNNER_LIMIT_FIELD_COUNT_V2];
const _: [(); 852] = [(); RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1];
const _: [(); 15] = [(); RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1];
const _: [(); 867] = [(); RUNNER_V2_BASE_VALUES_CELL_COUNT_V1];
const _: [(); 50] = [(); RUNNER_V2_RETAINED_DOMAIN_OBLIGATION_COUNT_V1];
const _: [(); 71] = [(); RUNNER_V2_LIMIT_MUTATION_OBLIGATION_COUNT_V1];
const _: [(); 12] = [(); RUNNER_V2_BASE_VALUES_REFUSAL_REASON_COUNT_V1];
const _: [(); 2] = [(); RUNNER_V2_BASE_VALUES_OWNER_SOURCE_COUNT_V1];
const _: [(); 16] = [(); RUNNER_V2_BASE_VALUES_DEPENDENCY_SOURCE_COUNT_V1];

/// Exact boundary grammar applied to every one of the 71 limit fields.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RunnerV2LimitBoundaryKindV1 {
    /// Exact zero value in the field's declared primitive width.
    Zero = 1,
    /// Exact one value in the field's declared primitive width.
    One = 2,
    /// Field-specific structural minimum.
    StructuralMinimum = 3,
    /// One step below the field-specific structural minimum.
    OneBelowStructuralMinimum = 4,
    /// Exact Smoke-profile ceiling.
    SmokeCeiling = 5,
    /// Canonical jointly feasible tightening of the Smoke ceiling.
    SmokeTightened = 6,
    /// One step above the Smoke-profile ceiling.
    SmokeOneOver = 7,
    /// Exact Full-profile ceiling.
    FullCeiling = 8,
    /// Canonical jointly feasible tightening of the Full ceiling.
    FullTightened = 9,
    /// One step above the Full-profile ceiling.
    FullOneOver = 10,
    /// Maximum value representable by the field's declared width.
    RepresentationalMaximum = 11,
    /// Checked refusal when stepping beyond the representational maximum.
    CheckedRepresentationalOverflowRefusal = 12,
}

impl RunnerV2LimitBoundaryKindV1 {
    /// Every boundary kind in frozen source order.
    pub const ALL: [Self; RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1] = [
        Self::Zero,
        Self::One,
        Self::StructuralMinimum,
        Self::OneBelowStructuralMinimum,
        Self::SmokeCeiling,
        Self::SmokeTightened,
        Self::SmokeOneOver,
        Self::FullCeiling,
        Self::FullTightened,
        Self::FullOneOver,
        Self::RepresentationalMaximum,
        Self::CheckedRepresentationalOverflowRefusal,
    ];

    /// Frozen one-based code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable source name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Zero => "zero",
            Self::One => "one",
            Self::StructuralMinimum => "structural-minimum",
            Self::OneBelowStructuralMinimum => "one-below-structural-minimum",
            Self::SmokeCeiling => "smoke-ceiling",
            Self::SmokeTightened => "smoke-tightened",
            Self::SmokeOneOver => "smoke-one-over",
            Self::FullCeiling => "full-ceiling",
            Self::FullTightened => "full-tightened",
            Self::FullOneOver => "full-one-over",
            Self::RepresentationalMaximum => "representational-maximum",
            Self::CheckedRepresentationalOverflowRefusal => {
                "checked-representational-overflow-refusal"
            }
        }
    }

    /// Profile against which this exact boundary is admitted.
    #[must_use]
    pub const fn profile(self) -> RunProfileV2 {
        match self {
            Self::SmokeCeiling | Self::SmokeTightened | Self::SmokeOneOver => RunProfileV2::Smoke,
            Self::Zero
            | Self::One
            | Self::StructuralMinimum
            | Self::OneBelowStructuralMinimum
            | Self::FullCeiling
            | Self::FullTightened
            | Self::FullOneOver
            | Self::RepresentationalMaximum
            | Self::CheckedRepresentationalOverflowRefusal => RunProfileV2::Full,
        }
    }

    /// Whether this boundary uses the canonical one-step tightening fixture.
    #[must_use]
    pub const fn is_tightened(self) -> bool {
        matches!(self, Self::SmokeTightened | Self::FullTightened)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RunnerV2LimitBoundaryDefinitionV1 {
    ordinal: u16,
    kind: RunnerV2LimitBoundaryKindV1,
    stable_name: &'static str,
}

const STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1: [RunnerV2LimitBoundaryDefinitionV1;
    RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1] = [
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 1,
        kind: RunnerV2LimitBoundaryKindV1::Zero,
        stable_name: "zero",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 2,
        kind: RunnerV2LimitBoundaryKindV1::One,
        stable_name: "one",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 3,
        kind: RunnerV2LimitBoundaryKindV1::StructuralMinimum,
        stable_name: "structural-minimum",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 4,
        kind: RunnerV2LimitBoundaryKindV1::OneBelowStructuralMinimum,
        stable_name: "one-below-structural-minimum",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 5,
        kind: RunnerV2LimitBoundaryKindV1::SmokeCeiling,
        stable_name: "smoke-ceiling",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 6,
        kind: RunnerV2LimitBoundaryKindV1::SmokeTightened,
        stable_name: "smoke-tightened",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 7,
        kind: RunnerV2LimitBoundaryKindV1::SmokeOneOver,
        stable_name: "smoke-one-over",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 8,
        kind: RunnerV2LimitBoundaryKindV1::FullCeiling,
        stable_name: "full-ceiling",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 9,
        kind: RunnerV2LimitBoundaryKindV1::FullTightened,
        stable_name: "full-tightened",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 10,
        kind: RunnerV2LimitBoundaryKindV1::FullOneOver,
        stable_name: "full-one-over",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 11,
        kind: RunnerV2LimitBoundaryKindV1::RepresentationalMaximum,
        stable_name: "representational-maximum",
    },
    RunnerV2LimitBoundaryDefinitionV1 {
        ordinal: 12,
        kind: RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal,
        stable_name: "checked-representational-overflow-refusal",
    },
];

const _: [(); 12] = [(); STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1.len()];

/// Canonical executable-family fixture used by every Stage-A limit cell.
///
/// Its one case has zero family rows, so the checked lifecycle requirement is
/// `3 + (2 + 0) = 5`. This is immutable declaration input, never an evaluator
/// result or a runtime-actual claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2LimitFixtureDeclarationV1 {
    executable: bool,
    family_rows_by_case: Box<[u32]>,
    declared_minimums_present_empty: bool,
    lifecycle_document_structural_minimum: u32,
    no_claim: StableTokenV2,
}

impl RunnerV2LimitFixtureDeclarationV1 {
    /// Whether the canonical family participates in executable admission.
    #[must_use]
    pub const fn executable(&self) -> bool {
        self.executable
    }

    /// Exact source-ordered family-row count for each canonical case.
    #[must_use]
    pub fn family_rows_by_case(&self) -> &[u32] {
        &self.family_rows_by_case
    }

    /// Whether the separate declared-minimum list is present and exactly empty.
    #[must_use]
    pub const fn declared_minimums_present_empty(&self) -> bool {
        self.declared_minimums_present_empty
    }

    /// Checked lifecycle-document minimum for this exact fixture.
    #[must_use]
    pub const fn lifecycle_document_structural_minimum(&self) -> u32 {
        self.lifecycle_document_structural_minimum
    }

    /// Exact fixture no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// One exact companion value applied only to keep a tightened limit candidate
/// jointly feasible while the target field is evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunnerV2LimitCompanionNormalizationV1 {
    field: RunnerLimitFieldV2,
    value: RunnerLimitValueV2,
}

impl RunnerV2LimitCompanionNormalizationV1 {
    /// Companion field changed by the canonical tightened fixture.
    #[must_use]
    pub const fn field(self) -> RunnerLimitFieldV2 {
        self.field
    }

    /// Exact same-width companion value.
    #[must_use]
    pub const fn value(self) -> RunnerLimitValueV2 {
        self.value
    }
}

/// Exact declaration-side expected partition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RunnerV2StageAExpectedPartitionV1 {
    /// Applicable operation expected to accept.
    EligiblePositive = 1,
    /// Applicable operation expected to return a typed refusal.
    ExpectedRefusal = 2,
    /// Applicable operation expected to return a modeled failure.
    ExpectedFailure = 3,
    /// Deliberate mutation expected to be detected.
    Mutation = 4,
    /// Registered operation that the implementation does not support.
    Unsupported = 5,
    /// Registered facet with no runtime operation.
    Inapplicable = 6,
}

impl RunnerV2StageAExpectedPartitionV1 {
    /// Frozen code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Exact verification facet assigned to one stable Stage-A cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RunnerV2StageACellGroupV1 {
    /// Literal, unit, and exact boundary verification.
    LiteralUnitBoundary = 1,
    /// Property-law and metamorphic verification.
    PropertyMetamorphic = 2,
    /// State, typestate, and model verification.
    StateModel = 3,
    /// Mutation and bounded fuzz-style refusal verification.
    MutationFuzz = 4,
    /// Public-API and compile-fail verification.
    ApiCompileFail = 5,
    /// Fault, resource, and cancellation-model verification.
    FaultResourceCancellation = 6,
    /// Real in-process integration without a mock evaluator.
    NoMockLocalIntegration = 7,
    /// Deterministic diagnostic logging and redaction verification.
    DetailedLoggingRedaction = 8,
    /// Structured reproduction declaration verification.
    ReproductionDeclaration = 9,
    /// Exact source-membership and source-identity verification.
    SourceClosure = 10,
}

impl RunnerV2StageACellGroupV1 {
    /// Frozen code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Facet of one retained pre-Runner-V2 base-value obligation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum RunnerV2RetainedDomainFacetV1 {
    /// Integer, rational, decimal, IEEE, and typed-literal behavior.
    NumericLiteral = 1,
    /// Physical and logical unit behavior.
    Unit = 2,
    /// Token, text, opaque-byte, and logical-path behavior.
    TokenTextPath = 3,
    /// Closed catalogs and nominal role/domain behavior.
    CatalogAndNominalIdentity = 4,
    /// Determinism, round-trip, and metamorphic behavior.
    PropertyAndMetamorphic = 5,
    /// Malformed, mutation, drift, and substitution refusals.
    MutationAndRefusal = 6,
    /// API privacy, compile-fail, and no-authority boundaries.
    ApiAndCompileFail = 7,
    /// Checked overflow, bounded resources, redaction, and integration.
    FaultResourceAndIntegration = 8,
}

impl RunnerV2RetainedDomainFacetV1 {
    /// Frozen facet code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// One exact retained domain obligation outside the 867 new Stage-A cells.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RetainedDomainObligationV1 {
    ordinal: u16,
    stable_id: StableTokenV2,
    facet: RunnerV2RetainedDomainFacetV1,
}

impl RunnerV2RetainedDomainObligationV1 {
    /// Exact one-based order.
    #[must_use]
    pub const fn ordinal(&self) -> u16 {
        self.ordinal
    }

    /// Stable test-obligation identity.
    #[must_use]
    pub const fn stable_id(&self) -> &StableTokenV2 {
        &self.stable_id
    }

    /// Exact retained coverage facet.
    #[must_use]
    pub const fn facet(&self) -> RunnerV2RetainedDomainFacetV1 {
        self.facet
    }
}

/// One explicit wrong-width mutation obligation for one limit field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2LimitMutationObligationV1 {
    ordinal: u16,
    stable_id: StableTokenV2,
    field: RunnerLimitFieldV2,
    field_name: StableTokenV2,
    declared_width: RunnerLimitWidthV2,
    opposite_width_zero: RunnerLimitValueV2,
    unit: RunnerLimitUnitV2,
    expected_reason: RunnerV2RawReasonV1,
    diagnostic_owner: StableTokenV2,
    repair_rank: u8,
    repair_kind: RepairActionKindV2,
    repair_target: StableTokenV2,
}

impl RunnerV2LimitMutationObligationV1 {
    /// Exact one-based mutation order.
    #[must_use]
    pub const fn ordinal(&self) -> u16 {
        self.ordinal
    }

    /// Stable mutation identity.
    #[must_use]
    pub const fn stable_id(&self) -> &StableTokenV2 {
        &self.stable_id
    }

    /// Field whose opposite primitive width must refuse.
    #[must_use]
    pub const fn field(&self) -> RunnerLimitFieldV2 {
        self.field
    }

    /// Independently declared stable field name.
    #[must_use]
    pub const fn field_name(&self) -> &StableTokenV2 {
        &self.field_name
    }

    /// Primitive width admitted by the field.
    #[must_use]
    pub const fn declared_width(&self) -> RunnerLimitWidthV2 {
        self.declared_width
    }

    /// Exact zero value in the opposite primitive width.
    #[must_use]
    pub const fn opposite_width_zero(&self) -> RunnerLimitValueV2 {
        self.opposite_width_zero
    }

    /// Exact semantic unit retained by the refusal.
    #[must_use]
    pub const fn unit(&self) -> RunnerLimitUnitV2 {
        self.unit
    }

    /// Exact raw reason produced by the wrong-width refusal.
    #[must_use]
    pub const fn expected_reason(&self) -> RunnerV2RawReasonV1 {
        self.expected_reason
    }

    /// Stable semantic owner of the refusal diagnostic.
    #[must_use]
    pub const fn diagnostic_owner(&self) -> &StableTokenV2 {
        &self.diagnostic_owner
    }

    /// One-based rank of the primary non-executable repair.
    #[must_use]
    pub const fn repair_rank(&self) -> u8 {
        self.repair_rank
    }

    /// Closed repair class for a wrong primitive width.
    #[must_use]
    pub const fn repair_kind(&self) -> RepairActionKindV2 {
        self.repair_kind
    }

    /// Stable semantic target of the repair.
    #[must_use]
    pub const fn repair_target(&self) -> &StableTokenV2 {
        &self.repair_target
    }
}

macro_rules! nominal_content_root {
    ($name:ident, $doc:literal) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name(ContentHash);

        impl $name {
            fn from_content_hash(root: ContentHash) -> Self {
                Self(root)
            }

            /// Exact 32-byte nominal identity.
            #[must_use]
            pub fn bytes(&self) -> &[u8; 32] {
                self.0.as_bytes()
            }
        }
    };
}

nominal_content_root!(
    RunnerV2StageADeclarationRootV1,
    "Nominal identity of the complete Stage-A base-values declaration."
);
nominal_content_root!(
    RunnerV2StageAOracleRootV1,
    "Nominal identity of one independent declaration-side oracle row."
);
nominal_content_root!(
    RunnerV2StageACaseManifestRootV1,
    "Nominal identity of one result-free operation and fixture manifest."
);
nominal_content_root!(
    RunnerV2StageASchemaInventoryRootV1,
    "Nominal identity of the direct Stage-A schema inventory."
);
nominal_content_root!(
    RunnerV2StageAFeatureDeclarationRootV1,
    "Nominal identity of the exact Stage-A feature declaration."
);
nominal_content_root!(
    RunnerV2StageAFiveExplicitsRootV1,
    "Nominal identity of declaration-side Stage-A Five Explicits."
);
nominal_content_root!(
    RunnerV2StageASourceMemberRootV1,
    "Nominal content identity of one child-owned source member."
);

/// Direct version/provenance requirements for Stage A.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAVersionRequirementsV1 {
    api_generation: RunnerApiGeneration,
    wire_version: RunnerWireVersion,
    predecessor_policy: WirePredecessorPolicyV1,
    source_identity: SourceIdentityRootV2,
    build_identity: BuildIdentityRootV2,
    toolchain_identity: ToolchainIdentityRootV2,
    schema_inventory_root: RunnerV2StageASchemaInventoryRootV1,
    feature_declaration_root: RunnerV2StageAFeatureDeclarationRootV1,
    target: BaseCoverageCloseTargetV1,
    profile: BaseCoverageCloseProfileV1,
}

impl RunnerV2StageAVersionRequirementsV1 {
    /// Declared Runner API generation.
    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        self.api_generation
    }

    /// Declared Runner wire version.
    #[must_use]
    pub const fn wire_version(&self) -> RunnerWireVersion {
        self.wire_version
    }

    /// Exact wire predecessor policy.
    #[must_use]
    pub const fn predecessor_policy(&self) -> WirePredecessorPolicyV1 {
        self.predecessor_policy
    }

    /// Identity of the complete current dependency-source closure.
    #[must_use]
    pub const fn source_identity(&self) -> &SourceIdentityRootV2 {
        &self.source_identity
    }

    /// Identity of the exact build-input declaration.
    #[must_use]
    pub const fn build_identity(&self) -> &BuildIdentityRootV2 {
        &self.build_identity
    }

    /// Identity of the exact declared Rust toolchain input.
    #[must_use]
    pub const fn toolchain_identity(&self) -> &ToolchainIdentityRootV2 {
        &self.toolchain_identity
    }

    /// Nominal identity of the direct Stage-A schema inventory.
    #[must_use]
    pub const fn schema_inventory_root(&self) -> RunnerV2StageASchemaInventoryRootV1 {
        self.schema_inventory_root
    }

    /// Nominal identity of the exact feature declaration.
    #[must_use]
    pub const fn feature_declaration_root(&self) -> RunnerV2StageAFeatureDeclarationRootV1 {
        self.feature_declaration_root
    }

    /// Declared target binding for pure Stage-A validation.
    #[must_use]
    pub const fn target(&self) -> BaseCoverageCloseTargetV1 {
        self.target
    }

    /// Declared execution-profile binding for Stage-A tests.
    #[must_use]
    pub const fn profile(&self) -> BaseCoverageCloseProfileV1 {
        self.profile
    }
}

/// Declaration-only Five Explicits for the pure Stage-A package.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAFiveExplicitsV1 {
    numeric_inputs_present_empty: bool,
    numeric_grants_present_empty: bool,
    expected_numeric_observations_present_empty: bool,
    seed: BaseCoverageCloseSeedExplicitV1,
    budgets: BaseCoverageCloseBudgetSetV1,
    versions: RunnerV2StageAVersionRequirementsV1,
    capability_registry: BaseCoverageCloseCapabilityRegistryV1,
    capability_profile_registry: BaseCoverageCloseCapabilityProfileRegistryV1,
    capability_contract: BaseCoverageCloseCapabilityContractV1,
    no_claim: StableTokenV2,
    root: RunnerV2StageAFiveExplicitsRootV1,
}

impl RunnerV2StageAFiveExplicitsV1 {
    /// Whether semantic numeric inputs are explicitly present and empty.
    #[must_use]
    pub const fn numeric_inputs_present_empty(&self) -> bool {
        self.numeric_inputs_present_empty
    }

    /// Whether semantic numeric grants are explicitly present and empty.
    #[must_use]
    pub const fn numeric_grants_present_empty(&self) -> bool {
        self.numeric_grants_present_empty
    }

    /// Whether expected numeric observations are explicitly present and empty.
    #[must_use]
    pub const fn expected_numeric_observations_present_empty(&self) -> bool {
        self.expected_numeric_observations_present_empty
    }

    /// Exact source-declared seed or seed-inapplicability policy.
    #[must_use]
    pub const fn seed(&self) -> &BaseCoverageCloseSeedExplicitV1 {
        &self.seed
    }

    /// Exact seven-axis declaration-side budget set.
    #[must_use]
    pub const fn budgets(&self) -> &BaseCoverageCloseBudgetSetV1 {
        &self.budgets
    }

    /// Direct version and provenance requirements.
    #[must_use]
    pub const fn versions(&self) -> &RunnerV2StageAVersionRequirementsV1 {
        &self.versions
    }

    /// Frozen semantic capability registry.
    #[must_use]
    pub const fn capability_registry(&self) -> &BaseCoverageCloseCapabilityRegistryV1 {
        &self.capability_registry
    }

    /// Frozen declaration-side capability-profile registry.
    #[must_use]
    pub const fn capability_profile_registry(
        &self,
    ) -> &BaseCoverageCloseCapabilityProfileRegistryV1 {
        &self.capability_profile_registry
    }

    /// Exact declaration-side `None` capability contract.
    #[must_use]
    pub const fn capability_contract(&self) -> &BaseCoverageCloseCapabilityContractV1 {
        &self.capability_contract
    }

    /// Exact Stage-A no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal identity of the complete declaration-side Five Explicits.
    #[must_use]
    pub const fn root(&self) -> RunnerV2StageAFiveExplicitsRootV1 {
        self.root
    }
}

/// Three declaration/execution/retention planes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum RunnerV2ContractPlaneV1 {
    /// Attempt-invariant canonical declaration or observation plane.
    Canonical = 1,
    /// Attempt-specific execution plane.
    Execution = 2,
    /// Attempt-specific retained-evidence plane.
    Retention = 4,
}

/// Exact closed set of included planes for one deferred common requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RunnerV2ContractPlaneSetV1(u8);

impl RunnerV2ContractPlaneSetV1 {
    const fn from_mask(mask: u8) -> Self {
        Self(mask)
    }

    /// Whether this exact set contains a plane.
    #[must_use]
    pub const fn contains(self, plane: RunnerV2ContractPlaneV1) -> bool {
        self.0 & plane as u8 != 0
    }

    /// Exact bit mask over the three closed planes.
    #[must_use]
    pub const fn mask(self) -> u8 {
        self.0
    }
}

/// Later foundational owner that realizes a deferred common contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum RunnerV2CommonFulfillmentStageV1 {
    /// Runtime evidence, actual-explicits, and reconciliation owner `.4`.
    RuntimeEvidence = 4,
    /// Route, owner, dispatcher, and shard/resume owner `.5`.
    RoutesAndDispatch = 5,
    /// Logging, reproduction, telemetry, and operator-view owner `.6`.
    LoggingAndReproduction = 6,
}

impl RunnerV2CommonFulfillmentStageV1 {
    /// Frozen owner-stage code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Uninhabited future root payload.
///
/// A requirement can carry only `TypedOptionV1::Absent` until its later owner
/// defines the nominal type and `.7` resolves it. There is no zero-digest or
/// generic-hash escape hatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RunnerV2UnavailableCommonRootV1 {}

/// Typed deferral of canonical schema-impact rows to the dedicated `.3` owner.
///
/// The rootless handoff is classified locally by [`RunnerV2RootlessAc58FragmentV1`];
/// every other owned nominal schema is enumerated here without fabricating the
/// future `.3` manifest root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2SchemaImpactDeferralV1 {
    resolution_owner: StableTokenV2,
    canonical_schema_names: Box<[RunnerV2CanonicalSchemaNameV1]>,
    future_manifest_root: TypedOptionV1<RunnerV2UnavailableCommonRootV1>,
    no_claim: StableTokenV2,
}

impl RunnerV2SchemaImpactDeferralV1 {
    /// Dedicated schema-registry work package that resolves the rows.
    #[must_use]
    pub const fn resolution_owner(&self) -> &StableTokenV2 {
        &self.resolution_owner
    }

    /// Exact schemas requiring one canonical AC58 row from `.3`.
    #[must_use]
    pub fn canonical_schema_names(&self) -> &[RunnerV2CanonicalSchemaNameV1] {
        &self.canonical_schema_names
    }

    /// Typed-absent future `.3` manifest identity.
    #[must_use]
    pub const fn future_manifest_root(&self) -> &TypedOptionV1<RunnerV2UnavailableCommonRootV1> {
        &self.future_manifest_root
    }

    /// Exact no-claim for the deferral.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// One exact source-authoritative requirement for a later common contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2CommonContractRequirementV1 {
    ordinal: u16,
    slot_id: StableTokenV2,
    api_generation: RunnerApiGeneration,
    wire_version: RunnerWireVersion,
    predecessor_policy: WirePredecessorPolicyV1,
    semantic_owner: StableTokenV2,
    realization_owner: StableTokenV2,
    future_nominal_role: StableTokenV2,
    future_domain: StableTokenV2,
    included_planes: RunnerV2ContractPlaneSetV1,
    fulfillment_stage: RunnerV2CommonFulfillmentStageV1,
    resolution_owner: StableTokenV2,
    future_root: TypedOptionV1<RunnerV2UnavailableCommonRootV1>,
    no_claim: StableTokenV2,
}

impl RunnerV2CommonContractRequirementV1 {
    /// Exact one-based requirement order.
    #[must_use]
    pub const fn ordinal(&self) -> u16 {
        self.ordinal
    }

    /// Stable identity of the deferred common-contract slot.
    #[must_use]
    pub const fn slot_id(&self) -> &StableTokenV2 {
        &self.slot_id
    }

    /// Runner API generation required by the future contract.
    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        self.api_generation
    }

    /// Runner wire version required by the future contract.
    #[must_use]
    pub const fn wire_version(&self) -> RunnerWireVersion {
        self.wire_version
    }

    /// Exact predecessor policy required by the future contract.
    #[must_use]
    pub const fn predecessor_policy(&self) -> WirePredecessorPolicyV1 {
        self.predecessor_policy
    }

    /// Parent that owns the requirement's semantics.
    #[must_use]
    pub const fn semantic_owner(&self) -> &StableTokenV2 {
        &self.semantic_owner
    }

    /// Later foundational work package that defines the contract.
    #[must_use]
    pub const fn realization_owner(&self) -> &StableTokenV2 {
        &self.realization_owner
    }

    /// Future nominal root role, without a fabricated root value.
    #[must_use]
    pub const fn future_nominal_role(&self) -> &StableTokenV2 {
        &self.future_nominal_role
    }

    /// Future canonical domain, without a fabricated frame.
    #[must_use]
    pub const fn future_domain(&self) -> &StableTokenV2 {
        &self.future_domain
    }

    /// Exact nonempty set of canonical, execution, and retention planes.
    #[must_use]
    pub const fn included_planes(&self) -> RunnerV2ContractPlaneSetV1 {
        self.included_planes
    }

    /// Foundational stage that must fulfill the requirement.
    #[must_use]
    pub const fn fulfillment_stage(&self) -> RunnerV2CommonFulfillmentStageV1 {
        self.fulfillment_stage
    }

    /// Final integration owner that resolves the future root.
    #[must_use]
    pub const fn resolution_owner(&self) -> &StableTokenV2 {
        &self.resolution_owner
    }

    /// Typed-absent future nominal root.
    #[must_use]
    pub const fn future_root(&self) -> &TypedOptionV1<RunnerV2UnavailableCommonRootV1> {
        &self.future_root
    }

    /// Exact requirement no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// One exact future broad-source member requirement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2FutureSourceRequirementV1 {
    final_ordinal: u16,
    path: LogicalBundlePathV1,
    future_content_root: TypedOptionV1<RunnerV2UnavailableCommonRootV1>,
}

impl RunnerV2FutureSourceRequirementV1 {
    /// Exact one-based ordinal in the eventual broad source inventory.
    #[must_use]
    pub const fn final_ordinal(&self) -> u16 {
        self.final_ordinal
    }

    /// Exact repository-relative future source path.
    #[must_use]
    pub const fn path(&self) -> &LogicalBundlePathV1 {
        &self.path
    }

    /// Typed-absent future content identity.
    #[must_use]
    pub const fn future_content_root(&self) -> &TypedOptionV1<RunnerV2UnavailableCommonRootV1> {
        &self.future_content_root
    }
}

/// One realized child-owned source member.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2OwnerSourceMemberV1 {
    path: LogicalBundlePathV1,
    content_root: RunnerV2StageASourceMemberRootV1,
}

/// One current source dependency whose content can change Stage-A semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2DependencySourceMemberV1 {
    path: LogicalBundlePathV1,
    content_root: RunnerV2StageASourceMemberRootV1,
}

impl RunnerV2DependencySourceMemberV1 {
    /// Exact repository-relative dependency path.
    #[must_use]
    pub const fn path(&self) -> &LogicalBundlePathV1 {
        &self.path
    }

    /// Domain-separated identity of the exact dependency bytes.
    #[must_use]
    pub const fn content_root(&self) -> RunnerV2StageASourceMemberRootV1 {
        self.content_root
    }
}

impl RunnerV2OwnerSourceMemberV1 {
    /// Exact repository-relative child-owned path.
    #[must_use]
    pub const fn path(&self) -> &LogicalBundlePathV1 {
        &self.path
    }

    /// Domain-separated identity of the exact child-owned bytes.
    #[must_use]
    pub const fn content_root(&self) -> RunnerV2StageASourceMemberRootV1 {
        self.content_root
    }
}

/// Lightweight Stage-A AC58 declaration for the rootless handoff.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2RootlessAc58FragmentV1 {
    semantic_type: RunnerV2RootlessHandoffSchemaNameV1,
    disposition: CanonicalSchemaImpactDispositionV1,
    migration_policy: CanonicalSchemaMigrationPolicyV1,
    authority_surfaces_present_empty: bool,
    no_claim: StableTokenV2,
}

impl RunnerV2RootlessAc58FragmentV1 {
    /// Exact semantic type classified by this lightweight fragment.
    #[must_use]
    pub const fn semantic_type(&self) -> &RunnerV2RootlessHandoffSchemaNameV1 {
        &self.semantic_type
    }

    /// AC58 disposition proving that no canonical frame applies.
    #[must_use]
    pub const fn disposition(&self) -> CanonicalSchemaImpactDispositionV1 {
        self.disposition
    }

    /// Explicit no-predecessor migration policy.
    #[must_use]
    pub const fn migration_policy(&self) -> CanonicalSchemaMigrationPolicyV1 {
        self.migration_policy
    }

    /// Whether the authority-surface set is explicitly present and empty.
    #[must_use]
    pub const fn authority_surfaces_present_empty(&self) -> bool {
        self.authority_surfaces_present_empty
    }

    /// Exact rootless-fragment no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// Exact single local route declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum RunnerV2LocalRouteClassV1 {
    /// Pure local in-process evaluator route; execution evidence is later-owned.
    LocalOnly = 1,
}

impl RunnerV2LocalRouteClassV1 {
    /// Frozen route-class code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// Exact single local route declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2LocalRouteDeclarationV1 {
    route_id: StableTokenV2,
    class: RunnerV2LocalRouteClassV1,
    public_entry_point: &'static str,
    execution_owner: StableTokenV2,
    capability_profile: BaseCoverageCloseCapabilityProfileV1,
    external_driver: TypedOptionV1<RunnerV2UnavailableCommonRootV1>,
    no_claim: StableTokenV2,
}

impl RunnerV2LocalRouteDeclarationV1 {
    /// Stable identity of the sole local route.
    #[must_use]
    pub const fn route_id(&self) -> &StableTokenV2 {
        &self.route_id
    }

    /// Exact `LocalOnly` route class.
    #[must_use]
    pub const fn class(&self) -> RunnerV2LocalRouteClassV1 {
        self.class
    }

    /// Frozen future public wrapper entry point.
    #[must_use]
    pub const fn public_entry_point(&self) -> &'static str {
        self.public_entry_point
    }

    /// Later work package that owns actual route execution.
    #[must_use]
    pub const fn execution_owner(&self) -> &StableTokenV2 {
        &self.execution_owner
    }

    /// Exact declaration-side capability profile.
    #[must_use]
    pub const fn capability_profile(&self) -> BaseCoverageCloseCapabilityProfileV1 {
        self.capability_profile
    }

    /// Typed-absent external driver for this local-only route.
    #[must_use]
    pub const fn external_driver(&self) -> &TypedOptionV1<RunnerV2UnavailableCommonRootV1> {
        &self.external_driver
    }

    /// Exact route no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// Explicit Stage-A shard or resume applicability.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAInapplicabilityDeclarationV1 {
    stable_id: StableTokenV2,
    reason: RunnerV2RawReasonV1,
    owner: StableTokenV2,
    prerequisite: StableTokenV2,
    no_claim: StableTokenV2,
}

impl RunnerV2StageAInapplicabilityDeclarationV1 {
    /// Stable identity of the inapplicability declaration.
    #[must_use]
    pub const fn stable_id(&self) -> &StableTokenV2 {
        &self.stable_id
    }

    /// Closed reason why the facet has no Stage-A runtime operation.
    #[must_use]
    pub const fn reason(&self) -> RunnerV2RawReasonV1 {
        self.reason
    }

    /// Source owner of the inapplicability decision.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Exact prerequisite that makes the declaration applicable.
    #[must_use]
    pub const fn prerequisite(&self) -> &StableTokenV2 {
        &self.prerequisite
    }

    /// Exact inapplicability no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

/// One declaration-side independent expected-oracle row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunnerV2StageAOracleNumericValueV1 {
    /// One exact heterogeneous Runner-limit value.
    Limit(RunnerLimitValueV2),
    /// One exact bounded count.
    Count(u64),
}

/// Unit of one independently declared oracle numeric value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunnerV2StageAOracleNumericUnitV1 {
    /// One field-specific Runner-limit unit.
    Limit(RunnerLimitUnitV2),
    /// The closed logical count unit.
    LogicalCount,
}

/// One independently declared expected numeric observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAOracleNumericV1 {
    name: StableTokenV2,
    value: RunnerV2StageAOracleNumericValueV1,
    unit: RunnerV2StageAOracleNumericUnitV1,
}

impl RunnerV2StageAOracleNumericV1 {
    /// Stable expected observation name.
    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Independently declared exact expected value.
    #[must_use]
    pub const fn value(&self) -> RunnerV2StageAOracleNumericValueV1 {
        self.value
    }

    /// Independently declared exact expected unit.
    #[must_use]
    pub const fn unit(&self) -> RunnerV2StageAOracleNumericUnitV1 {
        self.unit
    }
}

/// One independently declared expected repair descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAOracleRepairV1 {
    rank: u8,
    kind: RepairActionKindV2,
    target: StableTokenV2,
}

impl RunnerV2StageAOracleRepairV1 {
    /// One-based contiguous expected rank.
    #[must_use]
    pub const fn rank(&self) -> u8 {
        self.rank
    }

    /// Closed expected repair class.
    #[must_use]
    pub const fn kind(&self) -> RepairActionKindV2 {
        self.kind
    }

    /// Stable expected semantic target.
    #[must_use]
    pub const fn target(&self) -> &StableTokenV2 {
        &self.target
    }
}

/// Complete independently declared expected diagnostic projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAOracleDiagnosticV1 {
    code: DiagnosticCodeV2,
    owner: StableTokenV2,
    retryability: RetryabilityV2,
    prerequisites: Box<[StableTokenV2]>,
    repairs: Box<[RunnerV2StageAOracleRepairV1]>,
}

impl RunnerV2StageAOracleDiagnosticV1 {
    /// Closed expected diagnostic code.
    #[must_use]
    pub const fn code(&self) -> DiagnosticCodeV2 {
        self.code
    }

    /// Stable expected diagnostic owner.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Closed expected retryability.
    #[must_use]
    pub const fn retryability(&self) -> RetryabilityV2 {
        self.retryability
    }

    /// Exact source-ordered expected prerequisites.
    #[must_use]
    pub const fn prerequisites(&self) -> &[StableTokenV2] {
        &self.prerequisites
    }

    /// Exact ranked expected repair descriptors.
    #[must_use]
    pub const fn repairs(&self) -> &[RunnerV2StageAOracleRepairV1] {
        &self.repairs
    }
}

/// One declaration-side independent expected-oracle row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAOracleRowV1 {
    cell_id: StableTokenV2,
    expected_outcome: RunnerV2RawOutcomeKindV1,
    expected_reason: RunnerV2RawReasonV1,
    expected_partition: RunnerV2StageAExpectedPartitionV1,
    expected_numeric: Box<[RunnerV2StageAOracleNumericV1]>,
    expected_diagnostic: Option<RunnerV2StageAOracleDiagnosticV1>,
    root: RunnerV2StageAOracleRootV1,
}

/// One result-free parent-projection declaration for one stable cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageAProjectionRowV1 {
    ordinal: u16,
    cell_id: StableTokenV2,
    consumer_route: StableTokenV2,
    consumer_owner: StableTokenV2,
    dispatcher: StableTokenV2,
    posix_script: LogicalBundlePathV1,
    windows_script: LogicalBundlePathV1,
    expected_partition: RunnerV2StageAExpectedPartitionV1,
    case_manifest_root: RunnerV2StageACaseManifestRootV1,
    no_claim: StableTokenV2,
}

impl RunnerV2StageAProjectionRowV1 {
    /// Exact one-based projection-row order.
    #[must_use]
    pub const fn ordinal(&self) -> u16 {
        self.ordinal
    }

    /// Stable identity of the projected Stage-A cell.
    #[must_use]
    pub const fn cell_id(&self) -> &StableTokenV2 {
        &self.cell_id
    }

    /// Future route that consumes the declaration.
    #[must_use]
    pub const fn consumer_route(&self) -> &StableTokenV2 {
        &self.consumer_route
    }

    /// Final integration owner that consumes the declaration.
    #[must_use]
    pub const fn consumer_owner(&self) -> &StableTokenV2 {
        &self.consumer_owner
    }

    /// Future closed dispatcher entry point.
    #[must_use]
    pub const fn dispatcher(&self) -> &StableTokenV2 {
        &self.dispatcher
    }

    /// Future POSIX E2E wrapper path.
    #[must_use]
    pub const fn posix_script(&self) -> &LogicalBundlePathV1 {
        &self.posix_script
    }

    /// Future native-Windows E2E wrapper path.
    #[must_use]
    pub const fn windows_script(&self) -> &LogicalBundlePathV1 {
        &self.windows_script
    }

    /// Exact declaration-side expected partition.
    #[must_use]
    pub const fn expected_partition(&self) -> RunnerV2StageAExpectedPartitionV1 {
        self.expected_partition
    }

    /// Result-free identity of the exact cell operation and fixture.
    #[must_use]
    pub const fn case_manifest_root(&self) -> RunnerV2StageACaseManifestRootV1 {
        self.case_manifest_root
    }

    /// Exact projection no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }
}

impl RunnerV2StageAOracleRowV1 {
    /// Stable cell identity joined to this oracle.
    #[must_use]
    pub const fn cell_id(&self) -> &StableTokenV2 {
        &self.cell_id
    }

    /// Independently declared expected raw outcome.
    #[must_use]
    pub const fn expected_outcome(&self) -> RunnerV2RawOutcomeKindV1 {
        self.expected_outcome
    }

    /// Independently declared expected raw reason.
    #[must_use]
    pub const fn expected_reason(&self) -> RunnerV2RawReasonV1 {
        self.expected_reason
    }

    /// Independently declared expected verification partition.
    #[must_use]
    pub const fn expected_partition(&self) -> RunnerV2StageAExpectedPartitionV1 {
        self.expected_partition
    }

    /// Complete independently declared expected numeric projection.
    #[must_use]
    pub const fn expected_numeric(&self) -> &[RunnerV2StageAOracleNumericV1] {
        &self.expected_numeric
    }

    /// Complete independently declared expected diagnostic projection.
    #[must_use]
    pub const fn expected_diagnostic(&self) -> Option<&RunnerV2StageAOracleDiagnosticV1> {
        self.expected_diagnostic.as_ref()
    }

    /// Nominal identity of the complete independent oracle row.
    #[must_use]
    pub const fn root(&self) -> RunnerV2StageAOracleRootV1 {
        self.root
    }
}

/// One bounded non-limit operation in the Stage-A evaluator corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunnerV2StageAMetaOperationV1 {
    /// Typed absence is distinct from a present zero digest.
    TypedAbsenceDistinctFromZero,
    /// Binary32 ordering is available only through its named policy.
    F32NamedTotalOrder,
    /// Binary64 ordering is available only through its named policy.
    F64NamedTotalOrder,
    /// The capability contract is exactly the None profile.
    CapabilityNoneContract,
    /// Deferred common requirements have exact membership and order.
    CommonRequirementsExact,
    /// Reordering deferred requirements refuses.
    CommonRequirementReorderedRefusal,
    /// Future source paths have exact membership and order.
    FutureSourcesExact,
    /// The handoff has an explicit rootless AC58 disposition.
    RootlessAc58,
    /// The child source fragment has exactly two content-rooted members.
    OwnerSourceFragment,
    /// Exactly one local route is declared.
    LocalRoute,
    /// Sensitive observations are redacted without input echo.
    DiagnosticRedaction,
    /// Reproduction is declared for its later common owner.
    ReproductionDeclaration,
    /// Prohibited implicit float ordering remains compile-fail-only.
    CompileFailOrderingSurface,
    /// Sharding is inapplicable to this complete local evaluator.
    ShardInapplicable,
    /// Resume is inapplicable because each invocation recomputes all cells.
    ResumeInapplicable,
}

impl RunnerV2StageAMetaOperationV1 {
    /// Frozen one-based operation code.
    #[must_use]
    pub const fn code(self) -> u16 {
        meta_operation_code_v1(self)
    }
}

/// Exact result-free operation declared for one Stage-A cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunnerV2StageACellOperationV1 {
    /// One heterogeneous Runner-limit boundary admission.
    Limit {
        /// Exact Runner-limit field.
        field: RunnerLimitFieldV2,
        /// Exact boundary kind.
        boundary: RunnerV2LimitBoundaryKindV1,
        /// Typed input, absent only for inapplicable or overflow rows.
        value: TypedOptionV1<RunnerLimitValueV2>,
    },
    /// One bounded non-limit semantic operation.
    Meta(RunnerV2StageAMetaOperationV1),
}

/// One exact stable Stage-A cell declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2StageACellDeclarationV1 {
    ordinal: u16,
    cell_id: StableTokenV2,
    group: RunnerV2StageACellGroupV1,
    operation: RunnerV2StageACellOperationV1,
    companion_normalization: Box<[RunnerV2LimitCompanionNormalizationV1]>,
    oracle_root: RunnerV2StageAOracleRootV1,
    case_manifest_root: RunnerV2StageACaseManifestRootV1,
}

impl RunnerV2StageACellDeclarationV1 {
    /// Exact one-based cell order.
    #[must_use]
    pub const fn ordinal(&self) -> u16 {
        self.ordinal
    }

    /// Stable source-owned cell identity.
    #[must_use]
    pub const fn cell_id(&self) -> &StableTokenV2 {
        &self.cell_id
    }

    /// Exact verification group assigned to the cell.
    #[must_use]
    pub const fn group(&self) -> RunnerV2StageACellGroupV1 {
        self.group
    }

    /// Exact result-free operation declared for this cell.
    #[must_use]
    pub const fn operation(&self) -> RunnerV2StageACellOperationV1 {
        self.operation
    }

    /// Exact source-ordered companion normalization for this cell.
    #[must_use]
    pub fn companion_normalization(&self) -> &[RunnerV2LimitCompanionNormalizationV1] {
        &self.companion_normalization
    }

    /// Nominal identity of the independent expected-oracle row.
    #[must_use]
    pub const fn oracle_root(&self) -> RunnerV2StageAOracleRootV1 {
        self.oracle_root
    }

    /// Result-free nominal identity of the exact operation and shared fixture.
    #[must_use]
    pub const fn case_manifest_root(&self) -> RunnerV2StageACaseManifestRootV1 {
        self.case_manifest_root
    }
}

/// Complete source-authoritative Stage-A declaration for base values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseValuesStageADeclarationV1 {
    package_id: StableTokenV2,
    cells: Box<[RunnerV2StageACellDeclarationV1]>,
    oracles: Box<[RunnerV2StageAOracleRowV1]>,
    projections: Box<[RunnerV2StageAProjectionRowV1]>,
    limit_fixture: RunnerV2LimitFixtureDeclarationV1,
    retained_domain_obligations: Box<[RunnerV2RetainedDomainObligationV1]>,
    limit_mutation_obligations: Box<[RunnerV2LimitMutationObligationV1]>,
    five_explicits: RunnerV2StageAFiveExplicitsV1,
    route: RunnerV2LocalRouteDeclarationV1,
    common_requirements: Box<[RunnerV2CommonContractRequirementV1]>,
    future_sources: Box<[RunnerV2FutureSourceRequirementV1]>,
    owner_source_fragment: Box<[RunnerV2OwnerSourceMemberV1]>,
    dependency_source_closure: Box<[RunnerV2DependencySourceMemberV1]>,
    schema_impact_deferral: RunnerV2SchemaImpactDeferralV1,
    ac58: RunnerV2RootlessAc58FragmentV1,
    shard: RunnerV2StageAInapplicabilityDeclarationV1,
    resume: RunnerV2StageAInapplicabilityDeclarationV1,
    no_claim: StableTokenV2,
    root: RunnerV2StageADeclarationRootV1,
}

impl RunnerV2BaseValuesStageADeclarationV1 {
    /// Stable identity of this foundational work package.
    #[must_use]
    pub const fn package_id(&self) -> &StableTokenV2 {
        &self.package_id
    }

    /// Complete source-ordered cell declaration set.
    #[must_use]
    pub fn cells(&self) -> &[RunnerV2StageACellDeclarationV1] {
        &self.cells
    }

    /// Complete independent expected-oracle set in cell order.
    #[must_use]
    pub fn oracles(&self) -> &[RunnerV2StageAOracleRowV1] {
        &self.oracles
    }

    /// Complete result-free parent-projection set in cell order.
    #[must_use]
    pub fn projections(&self) -> &[RunnerV2StageAProjectionRowV1] {
        &self.projections
    }

    /// Canonical executable-family input shared by every limit cell.
    #[must_use]
    pub const fn limit_fixture(&self) -> &RunnerV2LimitFixtureDeclarationV1 {
        &self.limit_fixture
    }

    /// Exact legacy-domain coverage obligations preserved outside the 867 cells.
    #[must_use]
    pub fn retained_domain_obligations(&self) -> &[RunnerV2RetainedDomainObligationV1] {
        &self.retained_domain_obligations
    }

    /// One exact opposite-width mutation obligation for every limit field.
    #[must_use]
    pub fn limit_mutation_obligations(&self) -> &[RunnerV2LimitMutationObligationV1] {
        &self.limit_mutation_obligations
    }

    /// Complete declaration-side Five Explicits.
    #[must_use]
    pub const fn five_explicits(&self) -> &RunnerV2StageAFiveExplicitsV1 {
        &self.five_explicits
    }

    /// Sole local-only route declaration.
    #[must_use]
    pub const fn route(&self) -> &RunnerV2LocalRouteDeclarationV1 {
        &self.route
    }

    /// Exact ordered requirements for common contracts owned by `.4`–`.6`.
    #[must_use]
    pub fn common_requirements(&self) -> &[RunnerV2CommonContractRequirementV1] {
        &self.common_requirements
    }

    /// Exact typed-absent additions to the eventual broad source inventory.
    #[must_use]
    pub fn future_sources(&self) -> &[RunnerV2FutureSourceRequirementV1] {
        &self.future_sources
    }

    /// Exact two-member content-rooted child ownership fragment.
    #[must_use]
    pub fn owner_source_fragment(&self) -> &[RunnerV2OwnerSourceMemberV1] {
        &self.owner_source_fragment
    }

    /// Exact content-rooted source dependencies that can change Stage A.
    #[must_use]
    pub fn dependency_source_closure(&self) -> &[RunnerV2DependencySourceMemberV1] {
        &self.dependency_source_closure
    }

    /// Exact canonical-schema rows deferred to `.3` without a fabricated root.
    #[must_use]
    pub const fn schema_impact_deferral(&self) -> &RunnerV2SchemaImpactDeferralV1 {
        &self.schema_impact_deferral
    }

    /// Lightweight local AC58 classification of the rootless handoff.
    #[must_use]
    pub const fn ac58(&self) -> &RunnerV2RootlessAc58FragmentV1 {
        &self.ac58
    }

    /// Explicit declaration that Stage-A sharding is inapplicable.
    #[must_use]
    pub const fn shard(&self) -> &RunnerV2StageAInapplicabilityDeclarationV1 {
        &self.shard
    }

    /// Explicit declaration that Stage-A resume is inapplicable.
    #[must_use]
    pub const fn resume(&self) -> &RunnerV2StageAInapplicabilityDeclarationV1 {
        &self.resume
    }

    /// Exact declaration-wide no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal identity of the complete source-authoritative declaration.
    #[must_use]
    pub const fn root(&self) -> RunnerV2StageADeclarationRootV1 {
        self.root
    }
}

#[derive(Debug, Clone, Copy)]
struct RunnerV2LimitLiteralV1 {
    field: RunnerLimitFieldV2,
    width: RunnerLimitWidthV2,
    unit: RunnerLimitUnitV2,
    tightenability: RunnerLimitTightenabilityV2,
    minimum_rule: RunnerLimitMinimumRuleV2,
    smoke: RunnerLimitValueV2,
    full: RunnerLimitValueV2,
}

const fn opposite_width_zero_v1(width: RunnerLimitWidthV2) -> RunnerLimitValueV2 {
    match width {
        RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U64(0),
        RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U32(0),
    }
}

macro_rules! limit_literal_u32 {
    ($field:ident, $unit:ident, $tight:ident, $minimum:ident, $smoke:expr, $full:expr) => {
        RunnerV2LimitLiteralV1 {
            field: RunnerLimitFieldV2::$field,
            width: RunnerLimitWidthV2::U32,
            unit: RunnerLimitUnitV2::$unit,
            tightenability: RunnerLimitTightenabilityV2::$tight,
            minimum_rule: RunnerLimitMinimumRuleV2::$minimum,
            smoke: RunnerLimitValueV2::U32($smoke),
            full: RunnerLimitValueV2::U32($full),
        }
    };
}

macro_rules! limit_literal_u64 {
    ($field:ident, $unit:ident, $tight:ident, $minimum:ident, $smoke:expr, $full:expr) => {
        RunnerV2LimitLiteralV1 {
            field: RunnerLimitFieldV2::$field,
            width: RunnerLimitWidthV2::U64,
            unit: RunnerLimitUnitV2::$unit,
            tightenability: RunnerLimitTightenabilityV2::$tight,
            minimum_rule: RunnerLimitMinimumRuleV2::$minimum,
            smoke: RunnerLimitValueV2::U64($smoke),
            full: RunnerLimitValueV2::U64($full),
        }
    };
}

// This literal table is deliberately independent of
// `limits::RUNNER_LIMIT_DESCRIPTORS_V2` and `RunnerLimitsCandidateV2::base`.
// Conformance tests compare both directions; the oracle never reads either
// production table while deciding an expected result.
const STAGE_A_LIMIT_LITERALS_V1: [RunnerV2LimitLiteralV1; RUNNER_LIMIT_FIELD_COUNT_V2] = [
    limit_literal_u32!(ArgvTokens, Count, Tightenable, AtLeastOne, 64, 64),
    limit_literal_u64!(
        ArgvTokenBytes,
        LogicalBytes,
        Tightenable,
        AtLeastOne,
        8 * KIB,
        8 * KIB
    ),
    limit_literal_u64!(
        ArgvAggregateBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        64 * KIB,
        64 * KIB
    ),
    limit_literal_u64!(
        LifecycleRecordEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        16 * KIB,
        16 * KIB
    ),
    limit_literal_u32!(
        CaseLifecycleRecords,
        Records,
        Tightenable,
        ExecutableCaseAtLeastTwoRecords,
        256,
        256
    ),
    limit_literal_u64!(
        CaseLifecycleEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        256 * KIB,
        256 * KIB
    ),
    limit_literal_u32!(FamilyRowsPerCase, Rows, Tightenable, ZeroAllowed, 254, 254),
    limit_literal_u32!(
        InvocationCases,
        Count,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        256,
        256
    ),
    limit_literal_u32!(
        LifecycleDocumentRecords,
        Records,
        Tightenable,
        CheckedLifecycleEquation,
        4096,
        4096
    ),
    limit_literal_u64!(
        LifecycleDocumentEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        4 * MIB,
        4 * MIB
    ),
    limit_literal_u64!(
        CommandResultStdoutBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        5 * MIB,
        5 * MIB
    ),
    limit_literal_u64!(
        ChildStdoutBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        4 * MIB,
        4 * MIB
    ),
    limit_literal_u64!(
        CombinedChildStdoutBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        16 * MIB,
        128 * MIB
    ),
    limit_literal_u64!(
        ChildStderrBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        64 * KIB,
        64 * KIB
    ),
    limit_literal_u64!(
        CombinedChildStderrBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        256 * KIB,
        256 * KIB
    ),
    limit_literal_u64!(
        ManifestEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        MIB,
        MIB
    ),
    limit_literal_u32!(NestingDepth, Depth, Tightenable, AtLeastOne, 32, 32),
    limit_literal_u32!(
        ComparisonNodes,
        Nodes,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        256,
        256
    ),
    limit_literal_u32!(
        EffectNodes,
        Nodes,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        256,
        256
    ),
    limit_literal_u64!(
        TextBytes,
        LogicalBytes,
        Tightenable,
        AtLeastOne,
        8 * KIB,
        8 * KIB
    ),
    limit_literal_u64!(
        StableTokenBytes,
        LogicalBytes,
        Tightenable,
        AtLeastOne,
        128,
        128
    ),
    limit_literal_u64!(
        BundleRelativePathBytes,
        LogicalBytes,
        Tightenable,
        AtLeastOne,
        240,
        240
    ),
    limit_literal_u32!(
        DiagnosticsPerCase,
        Diagnostics,
        Tightenable,
        AtLeastOne,
        32,
        32
    ),
    limit_literal_u32!(
        DiagnosticsPerRun,
        Diagnostics,
        Tightenable,
        AtLeastOne,
        256,
        256
    ),
    limit_literal_u32!(
        PrerequisitesPerDiagnostic,
        Prerequisites,
        Tightenable,
        ZeroAllowed,
        16,
        16
    ),
    limit_literal_u32!(
        RepairsPerDiagnostic,
        Repairs,
        Tightenable,
        AtLeastOne,
        16,
        16
    ),
    limit_literal_u32!(Artifacts, Artifacts, Tightenable, ZeroAllowed, 256, 256),
    limit_literal_u64!(
        ArtifactEncodedBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        64 * MIB,
        64 * MIB
    ),
    limit_literal_u64!(
        ArtifactExpandedBytes,
        ExpandedBytes,
        Tightenable,
        ZeroAllowed,
        64 * MIB,
        64 * MIB
    ),
    limit_literal_u64!(
        ArtifactStoredBytes,
        StoredBytes,
        Tightenable,
        ZeroAllowed,
        64 * MIB + 4 * KIB,
        64 * MIB + 4 * KIB
    ),
    limit_literal_u64!(
        BundleEncodedBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        64 * MIB,
        512 * MIB
    ),
    limit_literal_u64!(
        BundleExpandedBytes,
        ExpandedBytes,
        Tightenable,
        ZeroAllowed,
        64 * MIB,
        512 * MIB
    ),
    limit_literal_u64!(
        ArtifactStoredAggregateBytes,
        StoredBytes,
        Tightenable,
        ZeroAllowed,
        65 * MIB,
        513 * MIB
    ),
    limit_literal_u64!(
        SystemPublicationStoredBytes,
        StoredBytes,
        Tightenable,
        ZeroAllowed,
        8 * MIB,
        8 * MIB
    ),
    limit_literal_u64!(
        PublicationStoredBytes,
        StoredBytes,
        Tightenable,
        ZeroAllowed,
        73 * MIB,
        521 * MIB
    ),
    limit_literal_u64!(
        ChildStreamDiscardBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        MIB,
        MIB
    ),
    limit_literal_u32!(
        ModesPerFamily,
        Count,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        64,
        64
    ),
    limit_literal_u32!(
        ExtensionDiagnosticsPerFamily,
        Diagnostics,
        Tightenable,
        ZeroAllowed,
        256,
        256
    ),
    limit_literal_u32!(
        ArtifactRolesPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RootPoliciesPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RegisteredUnitsPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        DigestDomainsPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        ExtensionSchemasPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        ExecutableDescriptorsPerFamily,
        Count,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        64,
        64
    ),
    limit_literal_u32!(MapEntries, Count, Tightenable, ZeroAllowed, 256, 256),
    limit_literal_u32!(
        GenericArrayItems,
        Count,
        Tightenable,
        ZeroAllowed,
        4096,
        4096
    ),
    limit_literal_u32!(PathSegments, Segments, Tightenable, AtLeastOne, 32, 32),
    limit_literal_u32!(IntegerDigits, Digits, Fixed, Fixed, 39, 39),
    limit_literal_u64!(RationalComponentBytes, EncodedBytes, Fixed, Fixed, 16, 16),
    limit_literal_u64!(DecimalCoefficientBytes, EncodedBytes, Fixed, Fixed, 16, 16),
    limit_literal_u32!(DecimalAbsoluteScale, DecimalScale, Fixed, Fixed, 6144, 6144),
    limit_literal_u32!(
        LogicalExtentsPerArtifact,
        Count,
        Tightenable,
        ZeroAllowed,
        16,
        16
    ),
    limit_literal_u32!(
        ObservationKeysPerCase,
        Count,
        Tightenable,
        ZeroAllowed,
        256,
        256
    ),
    limit_literal_u32!(
        DecisionDetailNamespaces,
        Namespaces,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(OutputClasses, Classes, Tightenable, ZeroAllowed, 64, 64),
    limit_literal_u64!(
        OpaqueValueBytes,
        LogicalBytes,
        Tightenable,
        AtLeastOne,
        8192,
        8192
    ),
    limit_literal_u64!(
        RetainedUnknownExtensionBytes,
        EncodedBytes,
        Tightenable,
        ZeroAllowed,
        65_536,
        65_536
    ),
    limit_literal_u32!(ExpressionEdges, Count, Tightenable, ZeroAllowed, 512, 512),
    limit_literal_u32!(
        MemoizedEvaluationVisits,
        Visits,
        Tightenable,
        ExecutableFamilyAtLeastOne,
        4096,
        4096
    ),
    limit_literal_u64!(
        RepairActionEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        1024,
        1024
    ),
    limit_literal_u64!(
        ActionableDiagnosticEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        8192,
        8192
    ),
    limit_literal_u64!(
        FailureStderrEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        16_384,
        16_384
    ),
    limit_literal_u64!(
        RunnerCatalogEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        MIB,
        MIB
    ),
    limit_literal_u64!(
        PublishedBundleReceiptEncodedBytes,
        EncodedBytes,
        Tightenable,
        AtLeastOne,
        MIB,
        MIB
    ),
    limit_literal_u64!(
        ContentStoreEnvelopeNonPayloadBytes,
        StoredBytes,
        Tightenable,
        ZeroAllowed,
        4096,
        4096
    ),
    limit_literal_u32!(
        RegisteredExtentAxesPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RegisteredObservationKeysPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        4096,
        4096
    ),
    limit_literal_u32!(
        RegisteredAuthorityScopesPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RegisteredExternalRootClassesPerFamily,
        Classes,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RegisteredEvaluationUnitsPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        64,
        64
    ),
    limit_literal_u32!(
        RegisteredResourceIdentitiesPerFamily,
        Count,
        Tightenable,
        ZeroAllowed,
        256,
        256
    ),
];
const _: [(); 71] = [(); STAGE_A_LIMIT_LITERALS_V1.len()];

// This independent literal name table deliberately does not read
// `RunnerLimitFieldV2::descriptor`. Its conformance test catches drift in
// either direction, while expected diagnostic repair targets remain
// independently specified.
const STAGE_A_INDEPENDENT_LIMIT_NAMES_V1: [&str; RUNNER_LIMIT_FIELD_COUNT_V2] = [
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

#[derive(Debug, Clone, Copy)]
struct CommonRequirementDefinitionV1 {
    slot_id: &'static str,
    realization_owner: &'static str,
    role: &'static str,
    domain: &'static str,
    planes: u8,
    stage: RunnerV2CommonFulfillmentStageV1,
}

const C: u8 = RunnerV2ContractPlaneV1::Canonical as u8;
const E: u8 = RunnerV2ContractPlaneV1::Execution as u8;
const R: u8 = RunnerV2ContractPlaneV1::Retention as u8;

const COMMON_REQUIREMENT_DEFINITIONS_V1: [CommonRequirementDefinitionV1;
    RUNNER_V2_COMMON_REQUIREMENT_COUNT_V1] = [
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.attempt-identity-contract.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "attempt-identity-root-v1",
        domain: "runner-v2-attempt-identity-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.canonical-runtime-observation-projection.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "canonical-runtime-observation-projection-root-v1",
        domain: "runner-v2-canonical-runtime-observation-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.runtime-evidence-envelope.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "runtime-evidence-envelope-root-v1",
        domain: "runner-v2-runtime-evidence-envelope-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.actual-five-explicits.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "actual-five-explicits-root-v1",
        domain: "runner-v2-actual-five-explicits-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.completeness-disposition.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "completeness-disposition-root-v1",
        domain: "runner-v2-completeness-disposition-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.safe-partial-evidence.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "safe-partial-evidence-root-v1",
        domain: "runner-v2-safe-partial-evidence-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.capability-reconciliation.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "capability-reconciliation-root-v1",
        domain: "runner-v2-capability-reconciliation-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.resource-reconciliation.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "resource-reconciliation-root-v1",
        domain: "runner-v2-resource-reconciliation-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.attempt-retention-receipt.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "attempt-retention-receipt-root-v1",
        domain: "runner-v2-attempt-retention-receipt-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.atomic-retention-finalization.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.4",
        role: "atomic-retention-finalization-root-v1",
        domain: "runner-v2-atomic-retention-finalization-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::RuntimeEvidence,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.route-schema.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "route-schema-root-v1",
        domain: "runner-v2-route-schema-v1",
        planes: C | E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.owner-matrix.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "owner-matrix-root-v1",
        domain: "runner-v2-owner-matrix-v1",
        planes: C | E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.route-registry.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "route-registry-root-v1",
        domain: "runner-v2-route-registry-v1",
        planes: C | E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.deferred-route-registry.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "deferred-route-registry-root-v1",
        domain: "runner-v2-deferred-route-registry-v1",
        planes: C | E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.dispatcher.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "dispatcher-root-v1",
        domain: "runner-v2-dispatcher-v1",
        planes: E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.native-bootstrap.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "native-bootstrap-root-v1",
        domain: "runner-v2-native-bootstrap-v1",
        planes: E,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.execution-source-binding.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "execution-source-binding-root-v1",
        domain: "runner-v2-execution-source-binding-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.retention-scope.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "retention-scope-root-v1",
        domain: "runner-v2-retention-scope-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.finalization-protocol.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "finalization-protocol-root-v1",
        domain: "runner-v2-finalization-protocol-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.shard-contract.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "shard-contract-root-v1",
        domain: "runner-v2-shard-contract-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.resume-contract.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.5",
        role: "resume-contract-root-v1",
        domain: "runner-v2-resume-contract-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.command-schema.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "command-schema-root-v1",
        domain: "runner-v2-command-schema-v1",
        planes: C | E,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.jsonl-event-schema.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "jsonl-event-schema-root-v1",
        domain: "runner-v2-jsonl-event-schema-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.terminal-reservation.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "terminal-reservation-root-v1",
        domain: "runner-v2-terminal-reservation-v1",
        planes: E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.redaction-policy.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "redaction-policy-root-v1",
        domain: "runner-v2-redaction-policy-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.first-divergence.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "first-divergence-root-v1",
        domain: "runner-v2-first-divergence-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.reproduction-schema.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "reproduction-schema-root-v1",
        domain: "runner-v2-reproduction-schema-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.relative-artifact-schema.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "relative-artifact-schema-root-v1",
        domain: "runner-v2-relative-artifact-schema-v1",
        planes: C | E | R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.raw-audit-binding.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "raw-audit-binding-root-v1",
        domain: "runner-v2-raw-audit-binding-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.stage-telemetry.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "stage-telemetry-root-v1",
        domain: "runner-v2-stage-telemetry-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
    CommonRequirementDefinitionV1 {
        slot_id: "runner-v2.common.operator-view.v1",
        realization_owner: "frankensim-epic-foundations-huq.24.1.1.1.6",
        role: "operator-view-root-v1",
        domain: "runner-v2-operator-view-v1",
        planes: R,
        stage: RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction,
    },
];

const FUTURE_SOURCE_PATHS_V1: [&str; RUNNER_V2_FUTURE_SOURCE_COUNT_V1] = [
    "crates/fs-evidence-runner/src/runner_v2.rs",
    "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/diagnostics.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/schema_registry.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/runtime_evidence.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/routes.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/detailed_logging.rs",
    "crates/fs-evidence-runner/src/runner_v2/work_packages/execution.rs",
    "crates/fs-evidence-runner/tests/runner_v2_base_work_packages.rs",
    "scripts/ci/runner_v2_base_work_packages_e2e.sh",
    "scripts/ci/runner_v2_base_work_packages_e2e.ps1",
];

#[derive(Clone, Copy)]
struct StageASourceDeclarationV1 {
    path: &'static str,
    bytes: &'static [u8],
}

macro_rules! stage_a_workspace_source_v1 {
    ($path:literal) => {
        StageASourceDeclarationV1 {
            path: $path,
            bytes: include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../", $path)),
        }
    };
}

macro_rules! stage_a_source_declarations_v1 {
    ($declarations:ident, $paths:ident, $count:expr, [$($path:literal),+ $(,)?]) => {
        const $declarations: [StageASourceDeclarationV1; $count] = [
            $(stage_a_workspace_source_v1!($path)),+
        ];
        const $paths: [&str; $count] = [$($path),+];
    };
}

stage_a_source_declarations_v1!(
    OWNER_SOURCE_DECLARATIONS_V1,
    OWNER_SOURCE_PATHS_V1,
    RUNNER_V2_BASE_VALUES_OWNER_SOURCE_COUNT_V1,
    [
        "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
    ]
);

stage_a_source_declarations_v1!(
    DEPENDENCY_SOURCE_DECLARATIONS_V1,
    DEPENDENCY_SOURCE_PATHS_V1,
    RUNNER_V2_BASE_VALUES_DEPENDENCY_SOURCE_COUNT_V1,
    [
        "crates/fs-blake3/src/lib.rs",
        "crates/fs-evidence-runner/src/lib.rs",
        "crates/fs-evidence-runner/src/canonical.rs",
        "crates/fs-evidence-runner/src/catalog.rs",
        "crates/fs-evidence-runner/src/construction.rs",
        "crates/fs-evidence-runner/src/coverage.rs",
        "crates/fs-evidence-runner/src/identity.rs",
        "crates/fs-evidence-runner/src/limits.rs",
        "crates/fs-evidence-runner/src/path.rs",
        "crates/fs-evidence-runner/src/projection.rs",
        "crates/fs-evidence-runner/src/schema_impact.rs",
        "crates/fs-evidence-runner/src/value.rs",
        "crates/fs-evidence-runner/src/runner_v2.rs",
        "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
    ]
);

const RETAINED_DOMAIN_OBLIGATION_DEFINITIONS_V1: [(&str, RunnerV2RetainedDomainFacetV1);
    RUNNER_V2_RETAINED_DOMAIN_OBLIGATION_COUNT_V1] = [
    (
        "signed-integer-width-and-extremes",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "unsigned-integer-width-and-extremes",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "integer-two-to-the-fifty-three-boundaries",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "rational-canonical-normalization",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "rational-zero-denominator-refusal",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "decimal-canonical-normalization",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "decimal-scale-minimum-maximum-and-one-over",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "binary32-exact-bit-identity",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "binary64-exact-bit-identity",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "binary32-named-ieee-total-order",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "binary64-named-ieee-total-order",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "ieee-signed-zero-nan-infinity-and-subnormal",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "duration-byte-and-rate-literal-catalogs",
        RunnerV2RetainedDomainFacetV1::NumericLiteral,
    ),
    (
        "physical-unit-positive-canonical-scale",
        RunnerV2RetainedDomainFacetV1::Unit,
    ),
    (
        "unit-dimension-and-scale-equivalence",
        RunnerV2RetainedDomainFacetV1::Unit,
    ),
    (
        "logical-unit-catalog-exact-order",
        RunnerV2RetainedDomainFacetV1::Unit,
    ),
    (
        "physical-logical-unit-substitution-refusal",
        RunnerV2RetainedDomainFacetV1::Unit,
    ),
    (
        "stable-token-empty-minimum-maximum-one-over",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "stable-token-separator-canonicality",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "text-empty-minimum-maximum-one-over",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "opaque-bytes-empty-maximum-one-over",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "logical-path-segment-empty-maximum-one-over",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "path-alias-and-normalization-refusal",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "reserved-prefix-refusal",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "non-ascii-and-platform-alias-boundaries",
        RunnerV2RetainedDomainFacetV1::TokenTextPath,
    ),
    (
        "closed-literal-catalog-code-name-order",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "unknown-catalog-code-refusal",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "nominal-root-role-and-domain-binding",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "raw-hash-as-nominal-root-refusal",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "typed-absence-versus-present-zero-digest",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "source-build-toolchain-role-substitution-refusal",
        RunnerV2RetainedDomainFacetV1::CatalogAndNominalIdentity,
    ),
    (
        "deterministic-construction-repeatability",
        RunnerV2RetainedDomainFacetV1::PropertyAndMetamorphic,
    ),
    (
        "canonical-round-trip-identity",
        RunnerV2RetainedDomainFacetV1::PropertyAndMetamorphic,
    ),
    (
        "unit-rescaling-invariance",
        RunnerV2RetainedDomainFacetV1::PropertyAndMetamorphic,
    ),
    (
        "named-ieee-order-repeatability",
        RunnerV2RetainedDomainFacetV1::PropertyAndMetamorphic,
    ),
    (
        "public-constructor-malformed-input-refusal",
        RunnerV2RetainedDomainFacetV1::MutationAndRefusal,
    ),
    (
        "missing-extra-duplicate-and-reordered-refusal",
        RunnerV2RetainedDomainFacetV1::MutationAndRefusal,
    ),
    (
        "cross-role-and-cross-domain-substitution-refusal",
        RunnerV2RetainedDomainFacetV1::MutationAndRefusal,
    ),
    (
        "profile-version-source-and-feature-drift",
        RunnerV2RetainedDomainFacetV1::MutationAndRefusal,
    ),
    (
        "ambient-unpinned-and-mixed-snapshot-refusal",
        RunnerV2RetainedDomainFacetV1::MutationAndRefusal,
    ),
    (
        "validated-wrapper-private-field-boundaries",
        RunnerV2RetainedDomainFacetV1::ApiAndCompileFail,
    ),
    (
        "binary32-no-ord-or-partial-ord",
        RunnerV2RetainedDomainFacetV1::ApiAndCompileFail,
    ),
    (
        "binary64-no-ord-or-partial-ord",
        RunnerV2RetainedDomainFacetV1::ApiAndCompileFail,
    ),
    (
        "ieee-wrapper-no-ordered-collection-or-sort",
        RunnerV2RetainedDomainFacetV1::ApiAndCompileFail,
    ),
    (
        "rootless-handoff-no-canonical-or-authority-surface",
        RunnerV2RetainedDomainFacetV1::ApiAndCompileFail,
    ),
    (
        "checked-arithmetic-overflow-refusal",
        RunnerV2RetainedDomainFacetV1::FaultResourceAndIntegration,
    ),
    (
        "bounded-allocation-and-one-over-refusal",
        RunnerV2RetainedDomainFacetV1::FaultResourceAndIntegration,
    ),
    (
        "real-no-mock-local-evaluator-integration",
        RunnerV2RetainedDomainFacetV1::FaultResourceAndIntegration,
    ),
    (
        "diagnostic-redaction-and-forbidden-value-no-echo",
        RunnerV2RetainedDomainFacetV1::FaultResourceAndIntegration,
    ),
    (
        "reproduction-and-compatible-source-closure-declaration",
        RunnerV2RetainedDomainFacetV1::FaultResourceAndIntegration,
    ),
];
const _: [(); 50] = [(); RETAINED_DOMAIN_OBLIGATION_DEFINITIONS_V1.len()];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MetaCellDefinitionV1 {
    ordinal: u16,
    id_suffix: &'static str,
    group: RunnerV2StageACellGroupV1,
    operation: RunnerV2StageAMetaOperationV1,
    expected_outcome: RunnerV2RawOutcomeKindV1,
    expected_reason: RunnerV2RawReasonV1,
    expected_partition: RunnerV2StageAExpectedPartitionV1,
}

const META_CELL_DEFINITIONS_V1: [MetaCellDefinitionV1; RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1] = [
    MetaCellDefinitionV1 {
        ordinal: 1,
        id_suffix: "typed-absence-distinct-from-zero",
        group: RunnerV2StageACellGroupV1::LiteralUnitBoundary,
        operation: RunnerV2StageAMetaOperationV1::TypedAbsenceDistinctFromZero,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 2,
        id_suffix: "f32-named-total-order",
        group: RunnerV2StageACellGroupV1::PropertyMetamorphic,
        operation: RunnerV2StageAMetaOperationV1::F32NamedTotalOrder,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 3,
        id_suffix: "f64-named-total-order",
        group: RunnerV2StageACellGroupV1::PropertyMetamorphic,
        operation: RunnerV2StageAMetaOperationV1::F64NamedTotalOrder,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 4,
        id_suffix: "capability-none-contract",
        group: RunnerV2StageACellGroupV1::StateModel,
        operation: RunnerV2StageAMetaOperationV1::CapabilityNoneContract,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 5,
        id_suffix: "common-requirements-exact",
        group: RunnerV2StageACellGroupV1::StateModel,
        operation: RunnerV2StageAMetaOperationV1::CommonRequirementsExact,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 6,
        id_suffix: "common-requirement-reordered-refusal",
        group: RunnerV2StageACellGroupV1::MutationFuzz,
        operation: RunnerV2StageAMetaOperationV1::CommonRequirementReorderedRefusal,
        expected_outcome: RunnerV2RawOutcomeKindV1::Refused,
        expected_reason: RunnerV2RawReasonV1::ExactMembershipMismatch,
        expected_partition: RunnerV2StageAExpectedPartitionV1::Mutation,
    },
    MetaCellDefinitionV1 {
        ordinal: 7,
        id_suffix: "future-sources-exact",
        group: RunnerV2StageACellGroupV1::SourceClosure,
        operation: RunnerV2StageAMetaOperationV1::FutureSourcesExact,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 8,
        id_suffix: "rootless-ac58",
        group: RunnerV2StageACellGroupV1::ApiCompileFail,
        operation: RunnerV2StageAMetaOperationV1::RootlessAc58,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 9,
        id_suffix: "owner-source-fragment",
        group: RunnerV2StageACellGroupV1::SourceClosure,
        operation: RunnerV2StageAMetaOperationV1::OwnerSourceFragment,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 10,
        id_suffix: "local-route",
        group: RunnerV2StageACellGroupV1::NoMockLocalIntegration,
        operation: RunnerV2StageAMetaOperationV1::LocalRoute,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 11,
        id_suffix: "diagnostic-redaction",
        group: RunnerV2StageACellGroupV1::DetailedLoggingRedaction,
        operation: RunnerV2StageAMetaOperationV1::DiagnosticRedaction,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 12,
        id_suffix: "reproduction-declaration",
        group: RunnerV2StageACellGroupV1::ReproductionDeclaration,
        operation: RunnerV2StageAMetaOperationV1::ReproductionDeclaration,
        expected_outcome: RunnerV2RawOutcomeKindV1::Accepted,
        expected_reason: RunnerV2RawReasonV1::ExactCheckedValue,
        expected_partition: RunnerV2StageAExpectedPartitionV1::EligiblePositive,
    },
    MetaCellDefinitionV1 {
        ordinal: 13,
        id_suffix: "compile-fail-ordering-surface",
        group: RunnerV2StageACellGroupV1::ApiCompileFail,
        operation: RunnerV2StageAMetaOperationV1::CompileFailOrderingSurface,
        expected_outcome: RunnerV2RawOutcomeKindV1::Inapplicable,
        expected_reason: RunnerV2RawReasonV1::PureDeclarationFacet,
        expected_partition: RunnerV2StageAExpectedPartitionV1::Inapplicable,
    },
    MetaCellDefinitionV1 {
        ordinal: 14,
        id_suffix: "shard-inapplicable",
        group: RunnerV2StageACellGroupV1::FaultResourceCancellation,
        operation: RunnerV2StageAMetaOperationV1::ShardInapplicable,
        expected_outcome: RunnerV2RawOutcomeKindV1::Inapplicable,
        expected_reason: RunnerV2RawReasonV1::ShardInapplicable,
        expected_partition: RunnerV2StageAExpectedPartitionV1::Inapplicable,
    },
    MetaCellDefinitionV1 {
        ordinal: 15,
        id_suffix: "resume-inapplicable",
        group: RunnerV2StageACellGroupV1::FaultResourceCancellation,
        operation: RunnerV2StageAMetaOperationV1::ResumeInapplicable,
        expected_outcome: RunnerV2RawOutcomeKindV1::Inapplicable,
        expected_reason: RunnerV2RawReasonV1::ResumeInapplicable,
        expected_partition: RunnerV2StageAExpectedPartitionV1::Inapplicable,
    },
];

#[derive(Debug, Clone)]
struct SourceOperationRowV1 {
    cell_id: StableTokenV2,
    group: RunnerV2StageACellGroupV1,
    operation: RunnerV2StageACellOperationV1,
}

#[derive(Debug, Clone, Copy)]
struct IndependentLimitViolationV1 {
    field: RunnerLimitFieldV2,
    expected: RunnerLimitExpectationV2,
    observed: RunnerLimitValueV2,
    unit: RunnerLimitUnitV2,
    reason: RunnerV2RawReasonV1,
}

#[derive(Clone)]
struct IndependentOracleProjectionV1 {
    outcome: RunnerV2RawOutcomeKindV1,
    reason: RunnerV2RawReasonV1,
    numeric: Vec<RunnerV2StageAOracleNumericV1>,
    diagnostic: Option<RunnerV2StageAOracleDiagnosticV1>,
}

/// Construct the complete immutable Stage-A declaration.
///
/// This function never invokes the evaluator and never binds an invocation
/// result, attempt, AC57 value, actual Five Explicits, retained artifact,
/// terminal log, reproduction instance, receipt, telemetry, or authority.
pub fn declare_24_1_1_1_1_v1() -> Result<RunnerV2BaseValuesStageADeclarationV1, ConstructionErrorV2>
{
    validate_limit_boundary_definitions_exact_v1(&STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1)?;
    validate_limit_literal_definitions_exact_v1(&STAGE_A_LIMIT_LITERALS_V1)?;
    validate_meta_definitions_exact_v1(&META_CELL_DEFINITIONS_V1)?;
    let package_id = stage_a_token(
        "runner_v2.base_values.package_id",
        RUNNER_V2_BASE_VALUES_PACKAGE_ID_V1,
    )?;
    let operations = source_operation_rows_v1()?;
    let oracles = independent_oracle_rows_v1()?;
    validate_operation_oracle_join_v1(&operations, &oracles)?;
    let limit_fixture = build_limit_fixture_v1()?;
    let dependency_source_closure = build_dependency_source_closure_v1()?;
    validate_dependency_source_closure_exact_v1(&dependency_source_closure)?;
    let five_explicits = build_five_explicits_v1(&dependency_source_closure)?;

    let cells = operations
        .iter()
        .zip(&oracles)
        .enumerate()
        .map(|(index, (operation, oracle))| {
            let companion_normalization = declared_companion_normalization_v1(operation.operation)?;
            let case_manifest_root = case_manifest_root_v1(
                operation,
                &companion_normalization,
                &limit_fixture,
                &five_explicits,
            )?;
            Ok(RunnerV2StageACellDeclarationV1 {
                ordinal: u16::try_from(index + 1).map_err(|_| {
                    stage_a_error(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "runner_v2.base_values.cell.ordinal",
                        "one-based cell ordinal representable as u16",
                        index,
                    )
                })?,
                cell_id: operation.cell_id.clone(),
                group: operation.group,
                operation: operation.operation,
                companion_normalization: companion_normalization.into_boxed_slice(),
                oracle_root: oracle.root,
                case_manifest_root,
            })
        })
        .collect::<Result<Vec<_>, ConstructionErrorV2>>()?;

    let projections = build_projection_rows_v1(&cells, &oracles)?;
    let retained_domain_obligations = build_retained_domain_obligations_v1()?;
    validate_retained_domain_obligations_exact_v1(&retained_domain_obligations)?;
    let limit_mutation_obligations = build_limit_mutation_obligations_v1()?;
    validate_limit_mutation_obligations_exact_v1(&limit_mutation_obligations)?;
    let route = build_route_v1()?;
    let common_requirements = build_common_requirements_v1()?;
    validate_common_requirements_exact_v1(&common_requirements)?;
    let future_sources = build_future_sources_v1()?;
    validate_future_sources_exact_v1(&future_sources)?;
    let owner_source_fragment = build_owner_source_fragment_v1()?;
    validate_owner_source_fragment_exact_v1(&owner_source_fragment)?;
    let schema_impact_deferral = build_schema_impact_deferral_v1()?;
    let ac58 = build_ac58_v1()?;
    let shard = build_inapplicability_v1(
        "runner-v2.base-values.shard-inapplicable.v1",
        RunnerV2RawReasonV1::ShardInapplicable,
        "complete-single-pass-local-evaluator",
    )?;
    let resume = build_inapplicability_v1(
        "runner-v2.base-values.resume-inapplicable.v1",
        RunnerV2RawReasonV1::ResumeInapplicable,
        "fresh-invocation-recomputes-complete-package",
    )?;
    let no_claim = stage_a_token(
        "runner_v2.base_values.no_claim",
        RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
    )?;

    let root = stage_a_declaration_root_v1(
        &package_id,
        &cells,
        &oracles,
        &projections,
        &limit_fixture,
        &retained_domain_obligations,
        &limit_mutation_obligations,
        &five_explicits,
        &route,
        &common_requirements,
        &future_sources,
        &owner_source_fragment,
        &dependency_source_closure,
        &schema_impact_deferral,
        &ac58,
        &shard,
        &resume,
        &no_claim,
    )?;

    Ok(RunnerV2BaseValuesStageADeclarationV1 {
        package_id,
        cells: cells.into_boxed_slice(),
        oracles: oracles.into_boxed_slice(),
        projections: projections.into_boxed_slice(),
        limit_fixture,
        retained_domain_obligations: retained_domain_obligations.into_boxed_slice(),
        limit_mutation_obligations: limit_mutation_obligations.into_boxed_slice(),
        five_explicits,
        route,
        common_requirements: common_requirements.into_boxed_slice(),
        future_sources: future_sources.into_boxed_slice(),
        owner_source_fragment: owner_source_fragment.into_boxed_slice(),
        dependency_source_closure: dependency_source_closure.into_boxed_slice(),
        schema_impact_deferral,
        ac58,
        shard,
        resume,
        no_claim,
        root,
    })
}

/// Freshly evaluate every source-declared base-values cell exactly once.
///
/// `Err` means the complete bounded report itself could not be constructed.
/// Ordinary domain refusals, modeled failures, and inapplicable facets are
/// retained as raw rows in the successful rootless handoff.
pub(crate) fn evaluate_24_1_1_1_1_cell_v1()
-> Result<RunnerV2LocalWorkPackageHandoffV1, ConstructionErrorV2> {
    let operations = source_operation_rows_v1()?;
    let limit_fixture = build_limit_fixture_v1()?;
    let declared_ids = operations
        .iter()
        .map(|row| row.cell_id.clone())
        .collect::<Vec<_>>();
    let rows = operations
        .iter()
        .map(|source| evaluate_source_operation_v1(source, &limit_fixture))
        .collect::<Result<Vec<_>, _>>()?;
    RunnerV2LocalWorkPackageHandoffV1::new(
        stage_a_token(
            "runner_v2.base_values.package_id",
            RUNNER_V2_BASE_VALUES_PACKAGE_ID_V1,
        )?,
        &declared_ids,
        rows,
    )
}

fn source_operation_rows_v1() -> Result<Vec<SourceOperationRowV1>, ConstructionErrorV2> {
    let mut rows = Vec::with_capacity(RUNNER_V2_BASE_VALUES_CELL_COUNT_V1);
    for literal in STAGE_A_LIMIT_LITERALS_V1 {
        for boundary in RunnerV2LimitBoundaryKindV1::ALL {
            rows.push(SourceOperationRowV1 {
                cell_id: limit_cell_id_v1(literal.field, boundary)?,
                group: RunnerV2StageACellGroupV1::LiteralUnitBoundary,
                operation: RunnerV2StageACellOperationV1::Limit {
                    field: literal.field,
                    boundary,
                    value: boundary_value_v1(literal, boundary),
                },
            });
        }
    }
    for (index, definition) in META_CELL_DEFINITIONS_V1.iter().enumerate() {
        rows.push(SourceOperationRowV1 {
            cell_id: meta_cell_id_v1(index, definition.id_suffix)?,
            group: definition.group,
            operation: RunnerV2StageACellOperationV1::Meta(definition.operation),
        });
    }
    if rows.len() != RUNNER_V2_BASE_VALUES_CELL_COUNT_V1 {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.operations",
            "the exact 852 limit plus 15 non-limit source operations",
            rows.len(),
        ));
    }
    let observed_ids = rows
        .iter()
        .map(|row| row.cell_id.as_str().to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_limit_cell_inventory_v1(
        &observed_ids[..RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1],
    )
    .map_err(|mismatch| mismatch.as_construction_error())?;
    validate_stage_a_complete_cell_inventory_v1(&observed_ids)
        .map_err(|mismatch| mismatch.as_construction_error())?;
    validate_unique_cell_ids_v1(rows.iter().map(|row| &row.cell_id))?;
    Ok(rows)
}

fn independent_oracle_rows_v1() -> Result<Vec<RunnerV2StageAOracleRowV1>, ConstructionErrorV2> {
    let mut rows = Vec::with_capacity(RUNNER_V2_BASE_VALUES_CELL_COUNT_V1);
    for literal in STAGE_A_LIMIT_LITERALS_V1 {
        for boundary in RunnerV2LimitBoundaryKindV1::ALL {
            let cell_id = limit_cell_id_v1(literal.field, boundary)?;
            let projection = independent_limit_oracle_projection_v1(literal, boundary)?;
            let expected_partition = match projection.outcome {
                RunnerV2RawOutcomeKindV1::Accepted => {
                    RunnerV2StageAExpectedPartitionV1::EligiblePositive
                }
                RunnerV2RawOutcomeKindV1::Refused => {
                    RunnerV2StageAExpectedPartitionV1::ExpectedRefusal
                }
                RunnerV2RawOutcomeKindV1::Failed => {
                    RunnerV2StageAExpectedPartitionV1::ExpectedFailure
                }
                RunnerV2RawOutcomeKindV1::Unsupported => {
                    RunnerV2StageAExpectedPartitionV1::Unsupported
                }
                RunnerV2RawOutcomeKindV1::Inapplicable => {
                    RunnerV2StageAExpectedPartitionV1::Inapplicable
                }
            };
            rows.push(build_oracle_row_v1(
                cell_id,
                projection.outcome,
                projection.reason,
                expected_partition,
                projection.numeric,
                projection.diagnostic,
            )?);
        }
    }
    for (index, definition) in META_CELL_DEFINITIONS_V1.iter().enumerate() {
        let projection = independent_meta_oracle_projection_v1(definition.operation)?;
        if projection.outcome != definition.expected_outcome
            || projection.reason != definition.expected_reason
        {
            return Err(stage_a_error(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.meta_oracle",
                "independent outcome and reason matching the frozen meta declaration",
                index,
            ));
        }
        rows.push(build_oracle_row_v1(
            meta_cell_id_v1(index, definition.id_suffix)?,
            projection.outcome,
            projection.reason,
            definition.expected_partition,
            projection.numeric,
            projection.diagnostic,
        )?);
    }
    if rows.len() != RUNNER_V2_BASE_VALUES_CELL_COUNT_V1 {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.oracles",
            "the exact independent oracle row for every source operation",
            rows.len(),
        ));
    }
    validate_unique_cell_ids_v1(rows.iter().map(|row| &row.cell_id))?;
    Ok(rows)
}

fn build_oracle_row_v1(
    cell_id: StableTokenV2,
    expected_outcome: RunnerV2RawOutcomeKindV1,
    expected_reason: RunnerV2RawReasonV1,
    expected_partition: RunnerV2StageAExpectedPartitionV1,
    expected_numeric: Vec<RunnerV2StageAOracleNumericV1>,
    expected_diagnostic: Option<RunnerV2StageAOracleDiagnosticV1>,
) -> Result<RunnerV2StageAOracleRowV1, ConstructionErrorV2> {
    let frame = CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-ORACLE\x01", 16 * 1024, |sink| {
        sink.push_str("runner_v2.base_values.oracle.cell_id", cell_id.as_str())?;
        sink.push_u16(
            "runner_v2.base_values.oracle.outcome",
            expected_outcome.code(),
        )?;
        sink.push_u16(
            "runner_v2.base_values.oracle.reason",
            expected_reason.code(),
        )?;
        sink.push_u16(
            "runner_v2.base_values.oracle.partition",
            expected_partition.code(),
        )?;
        sink.push_u32(
            "runner_v2.base_values.oracle.numeric_count",
            checked_u32_v1(
                "runner_v2.base_values.oracle.numeric_count",
                expected_numeric.len(),
            )?,
        )?;
        for numeric in &expected_numeric {
            push_oracle_numeric_v1(sink, numeric)?;
        }
        sink.push_presence(
            "runner_v2.base_values.oracle.diagnostic",
            expected_diagnostic.is_some(),
        )?;
        if let Some(diagnostic) = &expected_diagnostic {
            push_oracle_diagnostic_v1(sink, diagnostic)?;
        }
        Ok(())
    })?;
    Ok(RunnerV2StageAOracleRowV1 {
        cell_id,
        expected_outcome,
        expected_reason,
        expected_partition,
        expected_numeric: expected_numeric.into_boxed_slice(),
        expected_diagnostic,
        root: RunnerV2StageAOracleRootV1::from_content_hash(
            frame.root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.oracle.v1"),
        ),
    })
}

fn independent_limit_oracle_projection_v1(
    literal: RunnerV2LimitLiteralV1,
    boundary: RunnerV2LimitBoundaryKindV1,
) -> Result<IndependentOracleProjectionV1, ConstructionErrorV2> {
    let mut numeric = vec![oracle_count_v1(
        "field-ordinal",
        u64::from(literal.field.ordinal()),
    )?];
    let (outcome, reason, diagnostic) = match independent_boundary_value_v1(literal, boundary) {
        TypedOptionV1::Absent
            if matches!(
                boundary,
                RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal
            ) =>
        {
            numeric.push(oracle_limit_v1(
                "representational-maximum",
                representational_maximum_v1(literal.width),
                literal.unit,
            )?);
            let outcome = RunnerV2RawOutcomeKindV1::Refused;
            let reason = RunnerV2RawReasonV1::CheckedRepresentationalOverflow;
            (
                outcome,
                reason,
                Some(independent_oracle_diagnostic_v1(
                    outcome,
                    reason,
                    "fs-evidence-runner.runner-limits",
                    Some((
                        RepairActionKindV2::ChangeArguments,
                        independent_limit_name_v1(literal.field),
                    )),
                )?),
            )
        }
        TypedOptionV1::Absent => {
            let outcome = RunnerV2RawOutcomeKindV1::Inapplicable;
            let reason = RunnerV2RawReasonV1::PureDeclarationFacet;
            (
                outcome,
                reason,
                Some(independent_oracle_diagnostic_v1(
                    outcome,
                    reason,
                    "fs-evidence-runner.runner-v2.base-values",
                    None,
                )?),
            )
        }
        TypedOptionV1::Present(value) => {
            numeric.push(oracle_limit_v1("observed-value", value, literal.unit)?);
            let mut candidate = independent_base_vector_v1(boundary.profile());
            set_independent_value_v1(&mut candidate, literal.field, value);
            if boundary.is_tightened() {
                normalize_tightened_independent_v1(&mut candidate, literal.field);
            }
            match independent_candidate_violation_v1(boundary.profile(), &candidate) {
                None => (
                    RunnerV2RawOutcomeKindV1::Accepted,
                    RunnerV2RawReasonV1::ExactCheckedValue,
                    None,
                ),
                Some(violation) => {
                    append_independent_violation_numeric_v1(&mut numeric, violation)?;
                    let repair = if matches!(
                        violation.reason,
                        RunnerV2RawReasonV1::CheckedRepresentationalOverflow
                    ) {
                        None
                    } else {
                        Some((
                            match violation.expected {
                                RunnerLimitExpectationV2::AtMost(_) => {
                                    RepairActionKindV2::ReduceResourceDemand
                                }
                                RunnerLimitExpectationV2::Width(_)
                                | RunnerLimitExpectationV2::AtLeast(_)
                                | RunnerLimitExpectationV2::Exactly(_)
                                | RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
                                    RepairActionKindV2::UpdatePolicyOrCapability
                                }
                            },
                            independent_limit_name_v1(violation.field),
                        ))
                    };
                    (
                        RunnerV2RawOutcomeKindV1::Refused,
                        violation.reason,
                        Some(independent_oracle_diagnostic_v1(
                            RunnerV2RawOutcomeKindV1::Refused,
                            violation.reason,
                            "fs-evidence-runner.runner-limits",
                            repair,
                        )?),
                    )
                }
            }
        }
    };
    numeric.sort_by(|left, right| left.name.cmp(&right.name));
    Ok(IndependentOracleProjectionV1 {
        outcome,
        reason,
        numeric,
        diagnostic,
    })
}

fn independent_meta_oracle_projection_v1(
    operation: RunnerV2StageAMetaOperationV1,
) -> Result<IndependentOracleProjectionV1, ConstructionErrorV2> {
    let (outcome, reason, observed_count) = match operation {
        RunnerV2StageAMetaOperationV1::TypedAbsenceDistinctFromZero
        | RunnerV2StageAMetaOperationV1::F32NamedTotalOrder
        | RunnerV2StageAMetaOperationV1::F64NamedTotalOrder
        | RunnerV2StageAMetaOperationV1::CapabilityNoneContract
        | RunnerV2StageAMetaOperationV1::RootlessAc58
        | RunnerV2StageAMetaOperationV1::LocalRoute
        | RunnerV2StageAMetaOperationV1::DiagnosticRedaction
        | RunnerV2StageAMetaOperationV1::ReproductionDeclaration => (
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            1,
        ),
        RunnerV2StageAMetaOperationV1::CommonRequirementsExact => (
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            31,
        ),
        RunnerV2StageAMetaOperationV1::CommonRequirementReorderedRefusal => (
            RunnerV2RawOutcomeKindV1::Refused,
            RunnerV2RawReasonV1::ExactMembershipMismatch,
            31,
        ),
        RunnerV2StageAMetaOperationV1::FutureSourcesExact => (
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            13,
        ),
        RunnerV2StageAMetaOperationV1::OwnerSourceFragment => (
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            2,
        ),
        RunnerV2StageAMetaOperationV1::CompileFailOrderingSurface => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::PureDeclarationFacet,
            1,
        ),
        RunnerV2StageAMetaOperationV1::ShardInapplicable => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::ShardInapplicable,
            1,
        ),
        RunnerV2StageAMetaOperationV1::ResumeInapplicable => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::ResumeInapplicable,
            1,
        ),
    };
    let diagnostic = if matches!(outcome, RunnerV2RawOutcomeKindV1::Accepted) {
        None
    } else {
        Some(independent_oracle_diagnostic_v1(
            outcome,
            reason,
            "fs-evidence-runner.runner-v2.base-values",
            None,
        )?)
    };
    Ok(IndependentOracleProjectionV1 {
        outcome,
        reason,
        numeric: vec![oracle_count_v1("observed-count", observed_count)?],
        diagnostic,
    })
}

fn oracle_count_v1(
    name: &'static str,
    value: u64,
) -> Result<RunnerV2StageAOracleNumericV1, ConstructionErrorV2> {
    Ok(RunnerV2StageAOracleNumericV1 {
        name: stage_a_token("runner_v2.base_values.oracle.numeric.name", name)?,
        value: RunnerV2StageAOracleNumericValueV1::Count(value),
        unit: RunnerV2StageAOracleNumericUnitV1::LogicalCount,
    })
}

fn oracle_limit_v1(
    name: &'static str,
    value: RunnerLimitValueV2,
    unit: RunnerLimitUnitV2,
) -> Result<RunnerV2StageAOracleNumericV1, ConstructionErrorV2> {
    Ok(RunnerV2StageAOracleNumericV1 {
        name: stage_a_token("runner_v2.base_values.oracle.numeric.name", name)?,
        value: RunnerV2StageAOracleNumericValueV1::Limit(value),
        unit: RunnerV2StageAOracleNumericUnitV1::Limit(unit),
    })
}

fn append_independent_violation_numeric_v1(
    numeric: &mut Vec<RunnerV2StageAOracleNumericV1>,
    violation: IndependentLimitViolationV1,
) -> Result<(), ConstructionErrorV2> {
    numeric.push(oracle_count_v1(
        "violation-field-ordinal",
        u64::from(violation.field.ordinal()),
    )?);
    numeric.push(oracle_limit_v1(
        "violation-observed",
        violation.observed,
        violation.unit,
    )?);
    match violation.expected {
        RunnerLimitExpectationV2::AtMost(value)
        | RunnerLimitExpectationV2::AtLeast(value)
        | RunnerLimitExpectationV2::Exactly(value) => {
            numeric.push(oracle_limit_v1("expected-bound", value, violation.unit)?);
        }
        RunnerLimitExpectationV2::Width(width) => {
            numeric.push(oracle_count_v1(
                "expected-width-bits",
                match width {
                    RunnerLimitWidthV2::U32 => 32,
                    RunnerLimitWidthV2::U64 => 64,
                },
            )?);
        }
        RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
            numeric.push(oracle_count_v1("expected-minimum-ordinal-step", 1)?);
        }
    }
    Ok(())
}

fn independent_oracle_diagnostic_v1(
    outcome: RunnerV2RawOutcomeKindV1,
    reason: RunnerV2RawReasonV1,
    owner: &'static str,
    repair: Option<(RepairActionKindV2, &'static str)>,
) -> Result<RunnerV2StageAOracleDiagnosticV1, ConstructionErrorV2> {
    let (code, retryability) = match outcome {
        RunnerV2RawOutcomeKindV1::Accepted => {
            return Err(stage_a_error(
                ConstructionErrorKindV2::Unexpected,
                "runner_v2.base_values.oracle.diagnostic",
                "no diagnostic for an accepted oracle projection",
                outcome.code(),
            ));
        }
        RunnerV2RawOutcomeKindV1::Refused => (
            DiagnosticCodeV2::RunnerRefused,
            RetryabilityV2::AfterInputChange,
        ),
        RunnerV2RawOutcomeKindV1::Failed => {
            (DiagnosticCodeV2::RunnerInternalError, RetryabilityV2::Never)
        }
        RunnerV2RawOutcomeKindV1::Unsupported => (
            DiagnosticCodeV2::RunnerUnsupported,
            RetryabilityV2::AfterPrerequisiteChange,
        ),
        RunnerV2RawOutcomeKindV1::Inapplicable => {
            (DiagnosticCodeV2::RunnerNotRun, RetryabilityV2::Never)
        }
    };
    let prerequisites = if matches!(
        outcome,
        RunnerV2RawOutcomeKindV1::Inapplicable | RunnerV2RawOutcomeKindV1::Unsupported
    ) {
        vec![stage_a_token(
            "runner_v2.base_values.oracle.diagnostic.prerequisite",
            match (outcome, reason) {
                (
                    RunnerV2RawOutcomeKindV1::Inapplicable,
                    RunnerV2RawReasonV1::ShardInapplicable,
                ) => "registered-shard-policy",
                (
                    RunnerV2RawOutcomeKindV1::Inapplicable,
                    RunnerV2RawReasonV1::ResumeInapplicable,
                ) => "registered-resume-policy",
                (RunnerV2RawOutcomeKindV1::Unsupported, _) => "registered-evaluator-implementation",
                _ => "runtime-operation-applicability",
            },
        )?]
    } else {
        Vec::new()
    };
    let repairs = repair
        .map(|(kind, target)| {
            Ok(RunnerV2StageAOracleRepairV1 {
                rank: 1,
                kind,
                target: stage_a_token(
                    "runner_v2.base_values.oracle.diagnostic.repair_target",
                    target,
                )?,
            })
        })
        .into_iter()
        .collect::<Result<Vec<_>, ConstructionErrorV2>>()?;
    Ok(RunnerV2StageAOracleDiagnosticV1 {
        code,
        owner: stage_a_token("runner_v2.base_values.oracle.diagnostic.owner", owner)?,
        retryability,
        prerequisites: prerequisites.into_boxed_slice(),
        repairs: repairs.into_boxed_slice(),
    })
}

fn push_oracle_numeric_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    numeric: &RunnerV2StageAOracleNumericV1,
) -> Result<(), ConstructionErrorV2> {
    sink.push_str(
        "runner_v2.base_values.oracle.numeric.name",
        numeric.name.as_str(),
    )?;
    match numeric.value {
        RunnerV2StageAOracleNumericValueV1::Limit(value) => {
            sink.push_u16("runner_v2.base_values.oracle.numeric.value_kind", 1)?;
            push_limit_value_v1(sink, value)?;
        }
        RunnerV2StageAOracleNumericValueV1::Count(value) => {
            sink.push_u16("runner_v2.base_values.oracle.numeric.value_kind", 2)?;
            sink.push_bytes(
                "runner_v2.base_values.oracle.numeric.count",
                &value.to_be_bytes(),
            )?;
        }
    }
    match numeric.unit {
        RunnerV2StageAOracleNumericUnitV1::Limit(unit) => {
            sink.push_u16("runner_v2.base_values.oracle.numeric.unit_kind", 1)?;
            sink.push_str(
                "runner_v2.base_values.oracle.numeric.limit_unit",
                independent_limit_unit_name_v1(unit),
            )
        }
        RunnerV2StageAOracleNumericUnitV1::LogicalCount => {
            sink.push_u16("runner_v2.base_values.oracle.numeric.unit_kind", 2)?;
            sink.push_str("runner_v2.base_values.oracle.numeric.logical_unit", "count")
        }
    }
}

fn push_oracle_diagnostic_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    diagnostic: &RunnerV2StageAOracleDiagnosticV1,
) -> Result<(), ConstructionErrorV2> {
    sink.push_u16(
        "runner_v2.base_values.oracle.diagnostic.code",
        diagnostic.code.code(),
    )?;
    sink.push_str(
        "runner_v2.base_values.oracle.diagnostic.owner",
        diagnostic.owner.as_str(),
    )?;
    sink.push_u16(
        "runner_v2.base_values.oracle.diagnostic.retryability",
        diagnostic.retryability.code(),
    )?;
    sink.push_u32(
        "runner_v2.base_values.oracle.diagnostic.prerequisite_count",
        checked_u32_v1(
            "runner_v2.base_values.oracle.diagnostic.prerequisite_count",
            diagnostic.prerequisites.len(),
        )?,
    )?;
    for prerequisite in &diagnostic.prerequisites {
        sink.push_str(
            "runner_v2.base_values.oracle.diagnostic.prerequisite",
            prerequisite.as_str(),
        )?;
    }
    sink.push_u32(
        "runner_v2.base_values.oracle.diagnostic.repair_count",
        checked_u32_v1(
            "runner_v2.base_values.oracle.diagnostic.repair_count",
            diagnostic.repairs.len(),
        )?,
    )?;
    for repair in &diagnostic.repairs {
        sink.push_u8(
            "runner_v2.base_values.oracle.diagnostic.repair_rank",
            repair.rank,
        )?;
        sink.push_u16(
            "runner_v2.base_values.oracle.diagnostic.repair_kind",
            repair.kind.code(),
        )?;
        sink.push_str(
            "runner_v2.base_values.oracle.diagnostic.repair_target",
            repair.target.as_str(),
        )?;
    }
    Ok(())
}

const fn independent_limit_name_v1(field: RunnerLimitFieldV2) -> &'static str {
    STAGE_A_INDEPENDENT_LIMIT_NAMES_V1[field.ordinal() as usize - 1]
}

const fn independent_limit_unit_name_v1(unit: RunnerLimitUnitV2) -> &'static str {
    match unit {
        RunnerLimitUnitV2::Count => "count",
        RunnerLimitUnitV2::Records => "records",
        RunnerLimitUnitV2::Rows => "rows",
        RunnerLimitUnitV2::EncodedBytes => "encoded-bytes",
        RunnerLimitUnitV2::ExpandedBytes => "expanded-bytes",
        RunnerLimitUnitV2::StoredBytes => "stored-bytes",
        RunnerLimitUnitV2::LogicalBytes => "logical-bytes",
        RunnerLimitUnitV2::Depth => "depth",
        RunnerLimitUnitV2::Nodes => "nodes",
        RunnerLimitUnitV2::Digits => "digits",
        RunnerLimitUnitV2::Segments => "segments",
        RunnerLimitUnitV2::Diagnostics => "diagnostics",
        RunnerLimitUnitV2::Prerequisites => "prerequisites",
        RunnerLimitUnitV2::Repairs => "repairs",
        RunnerLimitUnitV2::Artifacts => "artifacts",
        RunnerLimitUnitV2::Namespaces => "namespaces",
        RunnerLimitUnitV2::Classes => "classes",
        RunnerLimitUnitV2::Visits => "visits",
        RunnerLimitUnitV2::DecimalScale => "decimal-scale",
    }
}

const fn runner_v2_raw_reason_name_v1(reason: RunnerV2RawReasonV1) -> &'static str {
    match reason {
        RunnerV2RawReasonV1::ExactCheckedValue => "exact-checked-value",
        RunnerV2RawReasonV1::BelowStructuralMinimum => "below-structural-minimum",
        RunnerV2RawReasonV1::AboveProfileCeiling => "above-profile-ceiling",
        RunnerV2RawReasonV1::FixedRepresentationChanged => "fixed-representation-changed",
        RunnerV2RawReasonV1::WrongPrimitiveWidth => "wrong-primitive-width",
        RunnerV2RawReasonV1::CheckedRepresentationalOverflow => "checked-representational-overflow",
        RunnerV2RawReasonV1::JointFeasibilityViolation => "joint-feasibility-violation",
        RunnerV2RawReasonV1::UnknownClosedValue => "unknown-closed-value",
        RunnerV2RawReasonV1::MalformedOrNoncanonical => "malformed-or-noncanonical",
        RunnerV2RawReasonV1::RequiredValueAbsent => "required-value-absent",
        RunnerV2RawReasonV1::UnexpectedValuePresent => "unexpected-value-present",
        RunnerV2RawReasonV1::ExactMembershipMismatch => "exact-membership-mismatch",
        RunnerV2RawReasonV1::PureDeclarationFacet => "pure-declaration-facet",
        RunnerV2RawReasonV1::CancellationInapplicable => "cancellation-inapplicable",
        RunnerV2RawReasonV1::ShardInapplicable => "shard-inapplicable",
        RunnerV2RawReasonV1::ResumeInapplicable => "resume-inapplicable",
        RunnerV2RawReasonV1::SourceDeclarationMismatch => "source-declaration-mismatch",
        RunnerV2RawReasonV1::InternalInvariantFailure => "internal-invariant-failure",
        RunnerV2RawReasonV1::UnsupportedClosedValue => "unsupported-closed-value",
    }
}

const fn limit_width_name_v1(width: RunnerLimitWidthV2) -> &'static str {
    match width {
        RunnerLimitWidthV2::U32 => "u32",
        RunnerLimitWidthV2::U64 => "u64",
    }
}

fn limit_value_safe_name_v1(value: RunnerLimitValueV2) -> String {
    match value {
        RunnerLimitValueV2::U32(value) => format!("u32:{value}"),
        RunnerLimitValueV2::U64(value) => format!("u64:{value}"),
    }
}

fn validate_operation_oracle_join_v1(
    operations: &[SourceOperationRowV1],
    oracles: &[RunnerV2StageAOracleRowV1],
) -> Result<(), ConstructionErrorV2> {
    if operations.len() < oracles.len() {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.base_values.operation_oracle_join",
            "one oracle for every operation",
            oracles.len(),
        ));
    }
    if operations.len() > oracles.len() {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Missing,
            "runner_v2.base_values.operation_oracle_join",
            "one oracle for every operation",
            oracles.len(),
        ));
    }
    for (index, (operation, oracle)) in operations.iter().zip(oracles).enumerate() {
        if operation.cell_id != oracle.cell_id {
            return Err(stage_a_error(
                ConstructionErrorKindV2::OutOfOrder,
                "runner_v2.base_values.operation_oracle_join",
                "the exact source-ordered operation/oracle cell identity",
                index,
            ));
        }
    }
    Ok(())
}

fn build_projection_rows_v1(
    cells: &[RunnerV2StageACellDeclarationV1],
    oracles: &[RunnerV2StageAOracleRowV1],
) -> Result<Vec<RunnerV2StageAProjectionRowV1>, ConstructionErrorV2> {
    if cells.len() != oracles.len() {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.projection_join",
            "one cell manifest for every independent oracle",
            cells.len(),
        ));
    }
    let consumer_route = stage_a_token(
        "runner_v2.base_values.projection.consumer_route",
        "runner-v2.route.24-1-1-1-7.execution-owner.v1",
    )?;
    let consumer_owner = stage_a_token(
        "runner_v2.base_values.projection.consumer_owner",
        "frankensim-epic-foundations-huq.24.1.1.1.7",
    )?;
    let dispatcher = stage_a_token(
        "runner_v2.base_values.projection.dispatcher",
        "runner-v2-work-package-dispatcher-v1",
    )?;
    let posix_script = stage_a_path(
        "runner_v2.base_values.projection.posix_script",
        "scripts/ci/runner_v2_base_work_packages_e2e.sh",
    )?;
    let windows_script = stage_a_path(
        "runner_v2.base_values.projection.windows_script",
        "scripts/ci/runner_v2_base_work_packages_e2e.ps1",
    )?;
    let no_claim = stage_a_token(
        "runner_v2.base_values.projection.no_claim",
        RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
    )?;

    oracles
        .iter()
        .zip(cells)
        .enumerate()
        .map(|(index, (oracle, cell))| {
            if oracle.cell_id != cell.cell_id {
                return Err(stage_a_error(
                    ConstructionErrorKindV2::OutOfOrder,
                    "runner_v2.base_values.projection_join",
                    "the exact source-ordered cell and oracle identity",
                    index,
                ));
            }
            Ok(RunnerV2StageAProjectionRowV1 {
                ordinal: u16::try_from(index + 1).map_err(|_| {
                    stage_a_error(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "runner_v2.base_values.projection.ordinal",
                        "one-based projection ordinal representable as u16",
                        index,
                    )
                })?,
                cell_id: oracle.cell_id.clone(),
                consumer_route: consumer_route.clone(),
                consumer_owner: consumer_owner.clone(),
                dispatcher: dispatcher.clone(),
                posix_script: posix_script.clone(),
                windows_script: windows_script.clone(),
                expected_partition: oracle.expected_partition,
                case_manifest_root: cell.case_manifest_root,
                no_claim: no_claim.clone(),
            })
        })
        .collect()
}

fn declared_companion_normalization_v1(
    operation: RunnerV2StageACellOperationV1,
) -> Result<Vec<RunnerV2LimitCompanionNormalizationV1>, ConstructionErrorV2> {
    let RunnerV2StageACellOperationV1::Limit {
        field,
        boundary,
        value,
    } = operation
    else {
        return Ok(Vec::new());
    };
    if !boundary.is_tightened() {
        return Ok(Vec::new());
    }
    let TypedOptionV1::Present(target_value) = value else {
        return Ok(Vec::new());
    };
    let mut rows = Vec::new();
    let mut push = |field, value| {
        rows.push(RunnerV2LimitCompanionNormalizationV1 { field, value });
    };
    match field {
        RunnerLimitFieldV2::ArgvAggregateBytes => {
            push(
                RunnerLimitFieldV2::ArgvTokenBytes,
                RunnerLimitValueV2::U64(1),
            );
        }
        RunnerLimitFieldV2::CaseLifecycleRecords => {
            let RunnerLimitValueV2::U32(records) = target_value else {
                return Err(stage_a_error(
                    ConstructionErrorKindV2::Incompatible,
                    "runner_v2.base_values.companion_normalization",
                    "u32 case-lifecycle target",
                    field.ordinal(),
                ));
            };
            push(
                RunnerLimitFieldV2::FamilyRowsPerCase,
                RunnerLimitValueV2::U32(records.saturating_sub(2)),
            );
        }
        RunnerLimitFieldV2::CaseLifecycleEncodedBytes => {
            push_declared_diagnostic_chain_v1(&mut push, true);
        }
        RunnerLimitFieldV2::LifecycleDocumentEncodedBytes => {
            push_declared_diagnostic_chain_v1(&mut push, true);
            push(
                RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
                RunnerLimitValueV2::U64(1),
            );
        }
        RunnerLimitFieldV2::CommandResultStdoutBytes => {
            push_declared_diagnostic_chain_v1(&mut push, true);
            for companion in [
                RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
                RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
                RunnerLimitFieldV2::RunnerCatalogEncodedBytes,
                RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes,
            ] {
                push(companion, RunnerLimitValueV2::U64(1));
            }
        }
        RunnerLimitFieldV2::CombinedChildStdoutBytes => push(
            RunnerLimitFieldV2::ChildStdoutBytes,
            RunnerLimitValueV2::U64(0),
        ),
        RunnerLimitFieldV2::CombinedChildStderrBytes => push(
            RunnerLimitFieldV2::ChildStderrBytes,
            RunnerLimitValueV2::U64(0),
        ),
        RunnerLimitFieldV2::BundleEncodedBytes
        | RunnerLimitFieldV2::BundleExpandedBytes
        | RunnerLimitFieldV2::ArtifactStoredAggregateBytes => {
            push(RunnerLimitFieldV2::Artifacts, RunnerLimitValueV2::U32(0))
        }
        RunnerLimitFieldV2::PublicationStoredBytes => push(
            RunnerLimitFieldV2::SystemPublicationStoredBytes,
            RunnerLimitValueV2::U64(0),
        ),
        RunnerLimitFieldV2::LifecycleRecordEncodedBytes
        | RunnerLimitFieldV2::FailureStderrEncodedBytes => {
            push_declared_diagnostic_chain_v1(&mut push, false);
        }
        RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes => push(
            RunnerLimitFieldV2::RepairActionEncodedBytes,
            RunnerLimitValueV2::U64(1),
        ),
        _ => {}
    }
    drop(push);
    rows.sort_by_key(|row| row.field.ordinal());
    if rows.windows(2).any(|pair| pair[0].field == pair[1].field) {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Duplicate,
            "runner_v2.base_values.companion_normalization",
            "at most one exact companion per field",
            field.ordinal(),
        ));
    }
    Ok(rows)
}

fn push_declared_diagnostic_chain_v1(
    push: &mut impl FnMut(RunnerLimitFieldV2, RunnerLimitValueV2),
    include_lifecycle_record: bool,
) {
    push(
        RunnerLimitFieldV2::RepairActionEncodedBytes,
        RunnerLimitValueV2::U64(1),
    );
    push(
        RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
        RunnerLimitValueV2::U64(1),
    );
    if include_lifecycle_record {
        push(
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            RunnerLimitValueV2::U64(1),
        );
    }
}

fn case_manifest_root_v1(
    source: &SourceOperationRowV1,
    companion_normalization: &[RunnerV2LimitCompanionNormalizationV1],
    fixture: &RunnerV2LimitFixtureDeclarationV1,
    five_explicits: &RunnerV2StageAFiveExplicitsV1,
) -> Result<RunnerV2StageACaseManifestRootV1, ConstructionErrorV2> {
    let frame =
        CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-CASE-MANIFEST\x01", 4096, |sink| {
            sink.push_str(
                "runner_v2.base_values.case_manifest.cell_id",
                source.cell_id.as_str(),
            )?;
            sink.push_u16(
                "runner_v2.base_values.case_manifest.group",
                source.group.code(),
            )?;
            push_operation_v1(sink, source.operation)?;
            sink.push_u32(
                "runner_v2.base_values.case_manifest.companion_count",
                checked_u32_v1(
                    "runner_v2.base_values.case_manifest.companion_count",
                    companion_normalization.len(),
                )?,
            )?;
            for companion in companion_normalization {
                sink.push_u16(
                    "runner_v2.base_values.case_manifest.companion_field",
                    companion.field.ordinal(),
                )?;
                push_limit_value_v1(sink, companion.value)?;
            }
            sink.push_presence(
                "runner_v2.base_values.case_manifest.fixture_executable",
                fixture.executable,
            )?;
            sink.push_u32(
                "runner_v2.base_values.case_manifest.fixture_case_count",
                checked_u32_v1(
                    "runner_v2.base_values.case_manifest.fixture_case_count",
                    fixture.family_rows_by_case.len(),
                )?,
            )?;
            for rows in &fixture.family_rows_by_case {
                sink.push_u32(
                    "runner_v2.base_values.case_manifest.fixture_family_rows",
                    *rows,
                )?;
            }
            sink.push_presence(
                "runner_v2.base_values.case_manifest.declared_minimums_present_empty",
                fixture.declared_minimums_present_empty,
            )?;
            sink.push_u32(
                "runner_v2.base_values.case_manifest.lifecycle_minimum",
                fixture.lifecycle_document_structural_minimum,
            )?;
            sink.push_str(
                "runner_v2.base_values.case_manifest.no_claim",
                fixture.no_claim.as_str(),
            )?;
            sink.push_fixed_bytes_32(
                "runner_v2.base_values.case_manifest.five_explicits",
                five_explicits.root.bytes(),
            )
        })?;
    Ok(RunnerV2StageACaseManifestRootV1::from_content_hash(
        frame.root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.case-manifest.v1"),
    ))
}

fn boundary_value_v1(
    literal: RunnerV2LimitLiteralV1,
    boundary: RunnerV2LimitBoundaryKindV1,
) -> TypedOptionV1<RunnerLimitValueV2> {
    let value = match boundary {
        RunnerV2LimitBoundaryKindV1::Zero => Some(same_width_value_v1(literal.width, 0)),
        RunnerV2LimitBoundaryKindV1::One => Some(same_width_value_v1(literal.width, 1)),
        RunnerV2LimitBoundaryKindV1::StructuralMinimum => structural_minimum_v1(literal),
        RunnerV2LimitBoundaryKindV1::OneBelowStructuralMinimum => {
            structural_minimum_v1(literal).and_then(checked_sub_one_v1)
        }
        RunnerV2LimitBoundaryKindV1::SmokeCeiling => Some(literal.smoke),
        RunnerV2LimitBoundaryKindV1::SmokeTightened => {
            canonical_tightened_value_v1(literal, literal.smoke)
        }
        RunnerV2LimitBoundaryKindV1::SmokeOneOver => checked_add_one_v1(literal.smoke),
        RunnerV2LimitBoundaryKindV1::FullCeiling => Some(literal.full),
        RunnerV2LimitBoundaryKindV1::FullTightened => {
            canonical_tightened_value_v1(literal, literal.full)
        }
        RunnerV2LimitBoundaryKindV1::FullOneOver => checked_add_one_v1(literal.full),
        RunnerV2LimitBoundaryKindV1::RepresentationalMaximum => {
            Some(representational_maximum_v1(literal.width))
        }
        RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal => None,
    };
    value.map_or(TypedOptionV1::Absent, TypedOptionV1::Present)
}

const fn structural_minimum_v1(literal: RunnerV2LimitLiteralV1) -> Option<RunnerLimitValueV2> {
    match literal.minimum_rule {
        RunnerLimitMinimumRuleV2::ZeroAllowed => Some(same_width_value_v1(literal.width, 0)),
        RunnerLimitMinimumRuleV2::AtLeastOne
        | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne => {
            Some(same_width_value_v1(literal.width, 1))
        }
        RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords => {
            Some(RunnerLimitValueV2::U32(2))
        }
        // The frozen Stage-A fixture is one executable case with zero family
        // rows, so the exact lifecycle equation is 3 + (2 + 0) = 5.
        RunnerLimitMinimumRuleV2::CheckedLifecycleEquation => Some(RunnerLimitValueV2::U32(5)),
        RunnerLimitMinimumRuleV2::Fixed => None,
    }
}

const fn canonical_tightened_value_v1(
    literal: RunnerV2LimitLiteralV1,
    ceiling: RunnerLimitValueV2,
) -> Option<RunnerLimitValueV2> {
    if matches!(literal.tightenability, RunnerLimitTightenabilityV2::Fixed) {
        return None;
    }
    let minimum = match structural_minimum_v1(literal) {
        Some(value) => value,
        None => return None,
    };
    if ceiling.as_u128() <= minimum.as_u128() {
        None
    } else {
        checked_sub_one_v1(ceiling)
    }
}

const fn same_width_value_v1(width: RunnerLimitWidthV2, value: u32) -> RunnerLimitValueV2 {
    match width {
        RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(value),
        RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(value as u64),
    }
}

const fn checked_sub_one_v1(value: RunnerLimitValueV2) -> Option<RunnerLimitValueV2> {
    match value {
        RunnerLimitValueV2::U32(value) => match value.checked_sub(1) {
            Some(value) => Some(RunnerLimitValueV2::U32(value)),
            None => None,
        },
        RunnerLimitValueV2::U64(value) => match value.checked_sub(1) {
            Some(value) => Some(RunnerLimitValueV2::U64(value)),
            None => None,
        },
    }
}

const fn checked_add_one_v1(value: RunnerLimitValueV2) -> Option<RunnerLimitValueV2> {
    match value {
        RunnerLimitValueV2::U32(value) => match value.checked_add(1) {
            Some(value) => Some(RunnerLimitValueV2::U32(value)),
            None => None,
        },
        RunnerLimitValueV2::U64(value) => match value.checked_add(1) {
            Some(value) => Some(RunnerLimitValueV2::U64(value)),
            None => None,
        },
    }
}

const fn representational_maximum_v1(width: RunnerLimitWidthV2) -> RunnerLimitValueV2 {
    match width {
        RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(u32::MAX),
        RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(u64::MAX),
    }
}

// This path deliberately does not call `boundary_value_v1`,
// `structural_minimum_v1`, or `canonical_tightened_value_v1`. Declaration
// construction and expected-result construction must not share the formula
// that selects the tested input.
fn independent_boundary_value_v1(
    literal: RunnerV2LimitLiteralV1,
    boundary: RunnerV2LimitBoundaryKindV1,
) -> TypedOptionV1<RunnerLimitValueV2> {
    let present = match boundary {
        RunnerV2LimitBoundaryKindV1::Zero => Some(match literal.width {
            RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(0),
            RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(0),
        }),
        RunnerV2LimitBoundaryKindV1::One => Some(match literal.width {
            RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(1),
            RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(1),
        }),
        RunnerV2LimitBoundaryKindV1::StructuralMinimum => {
            independent_structural_minimum_v1(literal)
        }
        RunnerV2LimitBoundaryKindV1::OneBelowStructuralMinimum => {
            match independent_structural_minimum_v1(literal) {
                Some(RunnerLimitValueV2::U32(value)) => {
                    value.checked_sub(1).map(RunnerLimitValueV2::U32)
                }
                Some(RunnerLimitValueV2::U64(value)) => {
                    value.checked_sub(1).map(RunnerLimitValueV2::U64)
                }
                None => None,
            }
        }
        RunnerV2LimitBoundaryKindV1::SmokeCeiling => Some(literal.smoke),
        RunnerV2LimitBoundaryKindV1::SmokeTightened => {
            independent_tightened_value_v1(literal, literal.smoke)
        }
        RunnerV2LimitBoundaryKindV1::SmokeOneOver => match literal.smoke {
            RunnerLimitValueV2::U32(value) => value.checked_add(1).map(RunnerLimitValueV2::U32),
            RunnerLimitValueV2::U64(value) => value.checked_add(1).map(RunnerLimitValueV2::U64),
        },
        RunnerV2LimitBoundaryKindV1::FullCeiling => Some(literal.full),
        RunnerV2LimitBoundaryKindV1::FullTightened => {
            independent_tightened_value_v1(literal, literal.full)
        }
        RunnerV2LimitBoundaryKindV1::FullOneOver => match literal.full {
            RunnerLimitValueV2::U32(value) => value.checked_add(1).map(RunnerLimitValueV2::U32),
            RunnerLimitValueV2::U64(value) => value.checked_add(1).map(RunnerLimitValueV2::U64),
        },
        RunnerV2LimitBoundaryKindV1::RepresentationalMaximum => Some(match literal.width {
            RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(u32::MAX),
            RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(u64::MAX),
        }),
        RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal => None,
    };
    present.map_or(TypedOptionV1::Absent, TypedOptionV1::Present)
}

fn independent_structural_minimum_v1(
    literal: RunnerV2LimitLiteralV1,
) -> Option<RunnerLimitValueV2> {
    match (literal.minimum_rule, literal.width) {
        (RunnerLimitMinimumRuleV2::ZeroAllowed, RunnerLimitWidthV2::U32) => {
            Some(RunnerLimitValueV2::U32(0))
        }
        (RunnerLimitMinimumRuleV2::ZeroAllowed, RunnerLimitWidthV2::U64) => {
            Some(RunnerLimitValueV2::U64(0))
        }
        (
            RunnerLimitMinimumRuleV2::AtLeastOne
            | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
            RunnerLimitWidthV2::U32,
        ) => Some(RunnerLimitValueV2::U32(1)),
        (
            RunnerLimitMinimumRuleV2::AtLeastOne
            | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
            RunnerLimitWidthV2::U64,
        ) => Some(RunnerLimitValueV2::U64(1)),
        (RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords, _) => {
            Some(RunnerLimitValueV2::U32(2))
        }
        (RunnerLimitMinimumRuleV2::CheckedLifecycleEquation, _) => Some(RunnerLimitValueV2::U32(
            RUNNER_V2_LIMIT_FIXTURE_LIFECYCLE_MINIMUM_V1,
        )),
        (RunnerLimitMinimumRuleV2::Fixed, _) => None,
    }
}

fn independent_tightened_value_v1(
    literal: RunnerV2LimitLiteralV1,
    ceiling: RunnerLimitValueV2,
) -> Option<RunnerLimitValueV2> {
    if literal.tightenability == RunnerLimitTightenabilityV2::Fixed {
        return None;
    }
    let minimum = independent_structural_minimum_v1(literal)?;
    if ceiling.as_u128() <= minimum.as_u128() {
        return None;
    }
    match ceiling {
        RunnerLimitValueV2::U32(value) => value.checked_sub(1).map(RunnerLimitValueV2::U32),
        RunnerLimitValueV2::U64(value) => value.checked_sub(1).map(RunnerLimitValueV2::U64),
    }
}

fn independent_base_vector_v1(profile: RunProfileV2) -> [RunnerLimitValueV2; 71] {
    let mut values = [RunnerLimitValueV2::U32(0); RUNNER_LIMIT_FIELD_COUNT_V2];
    for literal in STAGE_A_LIMIT_LITERALS_V1 {
        values[usize::from(literal.field.ordinal() - 1)] = match profile {
            RunProfileV2::Smoke => literal.smoke,
            RunProfileV2::Full => literal.full,
        };
    }
    values
}

fn set_independent_value_v1(
    candidate: &mut [RunnerLimitValueV2; 71],
    field: RunnerLimitFieldV2,
    value: RunnerLimitValueV2,
) {
    candidate[usize::from(field.ordinal() - 1)] = value;
}

fn independent_value_v1(
    candidate: &[RunnerLimitValueV2; 71],
    field: RunnerLimitFieldV2,
) -> RunnerLimitValueV2 {
    candidate[usize::from(field.ordinal() - 1)]
}

fn normalize_tightened_independent_v1(
    candidate: &mut [RunnerLimitValueV2; 71],
    field: RunnerLimitFieldV2,
) {
    let set_u32 =
        |candidate: &mut [RunnerLimitValueV2; 71], field: RunnerLimitFieldV2, value: u32| {
            set_independent_value_v1(candidate, field, RunnerLimitValueV2::U32(value));
        };
    let set_u64 =
        |candidate: &mut [RunnerLimitValueV2; 71], field: RunnerLimitFieldV2, value: u64| {
            set_independent_value_v1(candidate, field, RunnerLimitValueV2::U64(value));
        };
    match field {
        RunnerLimitFieldV2::ArgvAggregateBytes => {
            set_u64(candidate, RunnerLimitFieldV2::ArgvTokenBytes, 1);
        }
        RunnerLimitFieldV2::CaseLifecycleRecords => {
            if let RunnerLimitValueV2::U32(records) = independent_value_v1(candidate, field) {
                set_u32(
                    candidate,
                    RunnerLimitFieldV2::FamilyRowsPerCase,
                    records.saturating_sub(2),
                );
            }
        }
        RunnerLimitFieldV2::CaseLifecycleEncodedBytes => {
            lower_diagnostic_chain_v1(candidate, true, false);
        }
        RunnerLimitFieldV2::LifecycleDocumentEncodedBytes => {
            lower_diagnostic_chain_v1(candidate, true, false);
            set_u64(candidate, RunnerLimitFieldV2::CaseLifecycleEncodedBytes, 1);
        }
        RunnerLimitFieldV2::CommandResultStdoutBytes => {
            lower_diagnostic_chain_v1(candidate, true, false);
            set_u64(candidate, RunnerLimitFieldV2::CaseLifecycleEncodedBytes, 1);
            set_u64(
                candidate,
                RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
                1,
            );
            set_u64(candidate, RunnerLimitFieldV2::RunnerCatalogEncodedBytes, 1);
            set_u64(
                candidate,
                RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes,
                1,
            );
        }
        RunnerLimitFieldV2::CombinedChildStdoutBytes => {
            set_u64(candidate, RunnerLimitFieldV2::ChildStdoutBytes, 0);
        }
        RunnerLimitFieldV2::CombinedChildStderrBytes => {
            set_u64(candidate, RunnerLimitFieldV2::ChildStderrBytes, 0);
        }
        RunnerLimitFieldV2::BundleEncodedBytes
        | RunnerLimitFieldV2::BundleExpandedBytes
        | RunnerLimitFieldV2::ArtifactStoredAggregateBytes => {
            set_u32(candidate, RunnerLimitFieldV2::Artifacts, 0);
        }
        RunnerLimitFieldV2::PublicationStoredBytes => {
            set_u64(
                candidate,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
                0,
            );
        }
        RunnerLimitFieldV2::LifecycleRecordEncodedBytes
        | RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes
        | RunnerLimitFieldV2::FailureStderrEncodedBytes => {
            lower_diagnostic_chain_v1(candidate, false, true);
        }
        _ => {}
    }
}

fn lower_diagnostic_chain_v1(
    candidate: &mut [RunnerLimitValueV2; 71],
    lower_lifecycle_record: bool,
    keep_target: bool,
) {
    set_independent_value_v1(
        candidate,
        RunnerLimitFieldV2::RepairActionEncodedBytes,
        RunnerLimitValueV2::U64(1),
    );
    if !keep_target
        || !matches!(
            independent_value_v1(
                candidate,
                RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes
            ),
            RunnerLimitValueV2::U64(_)
        )
    {
        set_independent_value_v1(
            candidate,
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            RunnerLimitValueV2::U64(1),
        );
    } else if lower_lifecycle_record {
        set_independent_value_v1(
            candidate,
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            RunnerLimitValueV2::U64(1),
        );
    }
    if lower_lifecycle_record {
        set_independent_value_v1(
            candidate,
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            RunnerLimitValueV2::U64(1),
        );
    }
}

fn independent_candidate_violation_v1(
    profile: RunProfileV2,
    candidate: &[RunnerLimitValueV2; 71],
) -> Option<IndependentLimitViolationV1> {
    let base = independent_base_vector_v1(profile);
    for literal in STAGE_A_LIMIT_LITERALS_V1 {
        let observed = independent_value_v1(candidate, literal.field);
        let ceiling = independent_value_v1(&base, literal.field);
        if observed.width() != literal.width {
            return Some(independent_violation_v1(
                literal.field,
                RunnerLimitExpectationV2::Width(literal.width),
                observed,
                RunnerV2RawReasonV1::WrongPrimitiveWidth,
            ));
        }
        match literal.tightenability {
            RunnerLimitTightenabilityV2::Fixed if observed != ceiling => {
                return Some(independent_violation_v1(
                    literal.field,
                    RunnerLimitExpectationV2::Exactly(ceiling),
                    observed,
                    RunnerV2RawReasonV1::FixedRepresentationChanged,
                ));
            }
            RunnerLimitTightenabilityV2::Tightenable if observed.as_u128() > ceiling.as_u128() => {
                return Some(independent_violation_v1(
                    literal.field,
                    RunnerLimitExpectationV2::AtMost(ceiling),
                    observed,
                    RunnerV2RawReasonV1::AboveProfileCeiling,
                ));
            }
            RunnerLimitTightenabilityV2::Fixed | RunnerLimitTightenabilityV2::Tightenable => {}
        }
        let minimum = match literal.minimum_rule {
            RunnerLimitMinimumRuleV2::AtLeastOne
            | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne => {
                Some(same_width_value_v1(literal.width, 1))
            }
            RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords => {
                Some(RunnerLimitValueV2::U32(2))
            }
            RunnerLimitMinimumRuleV2::ZeroAllowed
            | RunnerLimitMinimumRuleV2::CheckedLifecycleEquation
            | RunnerLimitMinimumRuleV2::Fixed => None,
        };
        if let Some(minimum) = minimum
            && observed.as_u128() < minimum.as_u128()
        {
            return Some(independent_violation_v1(
                literal.field,
                RunnerLimitExpectationV2::AtLeast(minimum),
                observed,
                RunnerV2RawReasonV1::BelowStructuralMinimum,
            ));
        }
    }

    for (outer, inner) in [
        (
            RunnerLimitFieldV2::ArgvAggregateBytes,
            RunnerLimitFieldV2::ArgvTokenBytes,
        ),
        (
            RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
            RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::CombinedChildStdoutBytes,
            RunnerLimitFieldV2::ChildStdoutBytes,
        ),
        (
            RunnerLimitFieldV2::CombinedChildStderrBytes,
            RunnerLimitFieldV2::ChildStderrBytes,
        ),
        (
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            RunnerLimitFieldV2::RepairActionEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::FailureStderrEncodedBytes,
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            RunnerLimitFieldV2::RunnerCatalogEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes,
        ),
        (
            RunnerLimitFieldV2::PublicationStoredBytes,
            RunnerLimitFieldV2::SystemPublicationStoredBytes,
        ),
    ] {
        let outer_value = independent_value_v1(candidate, outer);
        let inner_value = independent_value_v1(candidate, inner);
        if outer_value.as_u128() < inner_value.as_u128() {
            return Some(independent_violation_v1(
                outer,
                RunnerLimitExpectationV2::AtLeast(inner_value),
                outer_value,
                RunnerV2RawReasonV1::JointFeasibilityViolation,
            ));
        }
    }

    if independent_value_v1(candidate, RunnerLimitFieldV2::Artifacts).as_u128() > 0 {
        for (outer, inner) in [
            (
                RunnerLimitFieldV2::BundleEncodedBytes,
                RunnerLimitFieldV2::ArtifactEncodedBytes,
            ),
            (
                RunnerLimitFieldV2::BundleExpandedBytes,
                RunnerLimitFieldV2::ArtifactExpandedBytes,
            ),
            (
                RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
                RunnerLimitFieldV2::ArtifactStoredBytes,
            ),
        ] {
            let outer_value = independent_value_v1(candidate, outer);
            let inner_value = independent_value_v1(candidate, inner);
            if outer_value.as_u128() < inner_value.as_u128() {
                return Some(independent_violation_v1(
                    outer,
                    RunnerLimitExpectationV2::AtLeast(inner_value),
                    outer_value,
                    RunnerV2RawReasonV1::JointFeasibilityViolation,
                ));
            }
        }
    }

    let invocation_value = independent_value_v1(candidate, RunnerLimitFieldV2::InvocationCases);
    let family_rows_value = independent_value_v1(candidate, RunnerLimitFieldV2::FamilyRowsPerCase);
    let case_records_value =
        independent_value_v1(candidate, RunnerLimitFieldV2::CaseLifecycleRecords);
    let document_records_value =
        independent_value_v1(candidate, RunnerLimitFieldV2::LifecycleDocumentRecords);
    let RunnerLimitValueV2::U32(invocation_cases) = invocation_value else {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitExpectationV2::Width(RunnerLimitWidthV2::U32),
            invocation_value,
            RunnerV2RawReasonV1::WrongPrimitiveWidth,
        ));
    };
    let RunnerLimitValueV2::U32(family_rows) = family_rows_value else {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::FamilyRowsPerCase,
            RunnerLimitExpectationV2::Width(RunnerLimitWidthV2::U32),
            family_rows_value,
            RunnerV2RawReasonV1::WrongPrimitiveWidth,
        ));
    };
    let RunnerLimitValueV2::U32(case_records) = case_records_value else {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::CaseLifecycleRecords,
            RunnerLimitExpectationV2::Width(RunnerLimitWidthV2::U32),
            case_records_value,
            RunnerV2RawReasonV1::WrongPrimitiveWidth,
        ));
    };
    let RunnerLimitValueV2::U32(document_records) = document_records_value else {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            RunnerLimitExpectationV2::Width(RunnerLimitWidthV2::U32),
            document_records_value,
            RunnerV2RawReasonV1::WrongPrimitiveWidth,
        ));
    };

    if invocation_cases > 0 {
        let Some(required_case_records) = family_rows.checked_add(2) else {
            return Some(independent_violation_v1(
                RunnerLimitFieldV2::CaseLifecycleRecords,
                RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(u32::MAX)),
                RunnerLimitValueV2::U32(family_rows),
                RunnerV2RawReasonV1::CheckedRepresentationalOverflow,
            ));
        };
        if case_records < required_case_records {
            return Some(independent_violation_v1(
                RunnerLimitFieldV2::CaseLifecycleRecords,
                RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(required_case_records)),
                RunnerLimitValueV2::U32(case_records),
                RunnerV2RawReasonV1::JointFeasibilityViolation,
            ));
        }
        if document_records < case_records {
            return Some(independent_violation_v1(
                RunnerLimitFieldV2::LifecycleDocumentRecords,
                RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(case_records)),
                RunnerLimitValueV2::U32(document_records),
                RunnerV2RawReasonV1::JointFeasibilityViolation,
            ));
        }
    }

    if invocation_cases < RUNNER_V2_LIMIT_FIXTURE_CASE_COUNT_V1 as u32 {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(invocation_cases)),
            RunnerLimitValueV2::U32(RUNNER_V2_LIMIT_FIXTURE_CASE_COUNT_V1 as u32),
            RunnerV2RawReasonV1::BelowStructuralMinimum,
        ));
    }
    if document_records < RUNNER_V2_LIMIT_FIXTURE_LIFECYCLE_MINIMUM_V1 {
        return Some(independent_violation_v1(
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(
                RUNNER_V2_LIMIT_FIXTURE_LIFECYCLE_MINIMUM_V1,
            )),
            RunnerLimitValueV2::U32(document_records),
            RunnerV2RawReasonV1::BelowStructuralMinimum,
        ));
    }
    None
}

fn independent_violation_v1(
    field: RunnerLimitFieldV2,
    expected: RunnerLimitExpectationV2,
    observed: RunnerLimitValueV2,
    reason: RunnerV2RawReasonV1,
) -> IndependentLimitViolationV1 {
    IndependentLimitViolationV1 {
        field,
        expected,
        observed,
        unit: STAGE_A_LIMIT_LITERALS_V1[usize::from(field.ordinal() - 1)].unit,
        reason,
    }
}

fn evaluate_source_operation_v1(
    source: &SourceOperationRowV1,
    limit_fixture: &RunnerV2LimitFixtureDeclarationV1,
) -> Result<RunnerV2RawCellObservationV1, ConstructionErrorV2> {
    match source.operation {
        RunnerV2StageACellOperationV1::Limit {
            field,
            boundary,
            value,
        } => evaluate_limit_operation_v1(
            source.cell_id.clone(),
            field,
            boundary,
            value,
            limit_fixture,
        ),
        RunnerV2StageACellOperationV1::Meta(operation) => {
            evaluate_meta_operation_v1(source.cell_id.clone(), operation)
        }
    }
}

fn evaluate_limit_operation_v1(
    cell_id: StableTokenV2,
    field: RunnerLimitFieldV2,
    boundary: RunnerV2LimitBoundaryKindV1,
    value: TypedOptionV1<RunnerLimitValueV2>,
    limit_fixture: &RunnerV2LimitFixtureDeclarationV1,
) -> Result<RunnerV2RawCellObservationV1, ConstructionErrorV2> {
    let literal = STAGE_A_LIMIT_LITERALS_V1[usize::from(field.ordinal() - 1)];
    let mut numeric = vec![RunnerV2SafeNumericObservationV1::count(
        stage_a_token("runner_v2.base_values.observation.name", "field-ordinal")?,
        u64::from(field.ordinal()),
    )];

    let (outcome, reason, diagnostic) = match value {
        TypedOptionV1::Absent
            if matches!(
                boundary,
                RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal
            ) =>
        {
            let maximum = representational_maximum_v1(literal.width);
            numeric.push(RunnerV2SafeNumericObservationV1::limit(
                stage_a_token(
                    "runner_v2.base_values.observation.name",
                    "representational-maximum",
                )?,
                maximum,
                literal.unit,
            ));
            if checked_add_one_v1(maximum).is_none() {
                let reason = RunnerV2RawReasonV1::CheckedRepresentationalOverflow;
                (
                    RunnerV2RawOutcomeKindV1::Refused,
                    reason,
                    Some(raw_diagnostic_v1(
                        RunnerV2RawOutcomeKindV1::Refused,
                        reason,
                        "fs-evidence-runner.runner-limits",
                        Some(field),
                    )?),
                )
            } else {
                let reason = RunnerV2RawReasonV1::InternalInvariantFailure;
                (
                    RunnerV2RawOutcomeKindV1::Failed,
                    reason,
                    Some(raw_diagnostic_v1(
                        RunnerV2RawOutcomeKindV1::Failed,
                        reason,
                        "fs-evidence-runner.runner-v2.base-values",
                        Some(field),
                    )?),
                )
            }
        }
        TypedOptionV1::Absent => {
            let reason = RunnerV2RawReasonV1::PureDeclarationFacet;
            (
                RunnerV2RawOutcomeKindV1::Inapplicable,
                reason,
                Some(raw_diagnostic_v1(
                    RunnerV2RawOutcomeKindV1::Inapplicable,
                    reason,
                    "fs-evidence-runner.runner-v2.base-values",
                    Some(field),
                )?),
            )
        }
        TypedOptionV1::Present(value) => {
            numeric.push(RunnerV2SafeNumericObservationV1::limit(
                stage_a_token("runner_v2.base_values.observation.name", "observed-value")?,
                value,
                literal.unit,
            ));
            let mut candidate = RunnerLimitsV2::base(boundary.profile()).to_candidate();
            if let Err(violation) = candidate.set_value(field, value) {
                let reason = limit_violation_reason_v1(&violation);
                append_limit_violation_observations_v1(&mut numeric, &violation)?;
                (
                    RunnerV2RawOutcomeKindV1::Refused,
                    reason,
                    Some(raw_limit_diagnostic_v1(&violation, reason)?),
                )
            } else {
                if boundary.is_tightened() {
                    normalize_tightened_candidate_v1(&mut candidate, field)?;
                }
                let requirements = RunnerFamilyLimitRequirementsV2 {
                    executable: limit_fixture.executable,
                    family_rows_by_case: &limit_fixture.family_rows_by_case,
                    declared_minimums: &[],
                };
                match RunnerLimitsV2::admit_family(boundary.profile(), candidate, requirements) {
                    Ok(_) => (
                        RunnerV2RawOutcomeKindV1::Accepted,
                        RunnerV2RawReasonV1::ExactCheckedValue,
                        None,
                    ),
                    Err(violation) => {
                        let reason = limit_violation_reason_v1(&violation);
                        append_limit_violation_observations_v1(&mut numeric, &violation)?;
                        (
                            RunnerV2RawOutcomeKindV1::Refused,
                            reason,
                            Some(raw_limit_diagnostic_v1(&violation, reason)?),
                        )
                    }
                }
            }
        }
    };

    numeric.sort_by(|left, right| left.name().cmp(right.name()));
    RunnerV2RawCellObservationV1::new(cell_id, outcome, reason, numeric, diagnostic)
}

fn normalize_tightened_candidate_v1(
    candidate: &mut RunnerLimitsCandidateV2,
    field: RunnerLimitFieldV2,
) -> Result<(), ConstructionErrorV2> {
    let case_lifecycle_records = if field == RunnerLimitFieldV2::CaseLifecycleRecords {
        let RunnerLimitValueV2::U32(records) = candidate.value(field) else {
            return Err(stage_a_error(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.tightened_fixture",
                "u32 case-lifecycle record value",
                field.ordinal(),
            ));
        };
        Some(records)
    } else {
        None
    };
    let mut set = |field, value| {
        candidate.set_value(field, value).map_err(|violation| {
            stage_a_error(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.tightened_fixture",
                "the exact source-frozen same-width companion value",
                violation.field().ordinal(),
            )
        })
    };
    match field {
        RunnerLimitFieldV2::ArgvAggregateBytes => {
            set(
                RunnerLimitFieldV2::ArgvTokenBytes,
                RunnerLimitValueV2::U64(1),
            )?;
        }
        RunnerLimitFieldV2::CaseLifecycleRecords => {
            let records = case_lifecycle_records.ok_or_else(|| {
                stage_a_error(
                    ConstructionErrorKindV2::Incompatible,
                    "runner_v2.base_values.tightened_fixture",
                    "captured u32 case-lifecycle record value",
                    field.ordinal(),
                )
            })?;
            set(
                RunnerLimitFieldV2::FamilyRowsPerCase,
                RunnerLimitValueV2::U32(records.saturating_sub(2)),
            )?;
        }
        RunnerLimitFieldV2::CaseLifecycleEncodedBytes => {
            set_production_diagnostic_chain_v1(&mut set, true, true)?;
        }
        RunnerLimitFieldV2::LifecycleDocumentEncodedBytes => {
            set_production_diagnostic_chain_v1(&mut set, true, true)?;
            set(
                RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
                RunnerLimitValueV2::U64(1),
            )?;
        }
        RunnerLimitFieldV2::CommandResultStdoutBytes => {
            set_production_diagnostic_chain_v1(&mut set, true, true)?;
            for companion in [
                RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
                RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
                RunnerLimitFieldV2::RunnerCatalogEncodedBytes,
                RunnerLimitFieldV2::PublishedBundleReceiptEncodedBytes,
            ] {
                set(companion, RunnerLimitValueV2::U64(1))?;
            }
        }
        RunnerLimitFieldV2::CombinedChildStdoutBytes => {
            set(
                RunnerLimitFieldV2::ChildStdoutBytes,
                RunnerLimitValueV2::U64(0),
            )?;
        }
        RunnerLimitFieldV2::CombinedChildStderrBytes => {
            set(
                RunnerLimitFieldV2::ChildStderrBytes,
                RunnerLimitValueV2::U64(0),
            )?;
        }
        RunnerLimitFieldV2::BundleEncodedBytes
        | RunnerLimitFieldV2::BundleExpandedBytes
        | RunnerLimitFieldV2::ArtifactStoredAggregateBytes => {
            set(RunnerLimitFieldV2::Artifacts, RunnerLimitValueV2::U32(0))?;
        }
        RunnerLimitFieldV2::PublicationStoredBytes => {
            set(
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
                RunnerLimitValueV2::U64(0),
            )?;
        }
        RunnerLimitFieldV2::LifecycleRecordEncodedBytes
        | RunnerLimitFieldV2::FailureStderrEncodedBytes => {
            set(
                RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
                RunnerLimitValueV2::U64(1),
            )?;
            set(
                RunnerLimitFieldV2::RepairActionEncodedBytes,
                RunnerLimitValueV2::U64(1),
            )?;
        }
        RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes => {
            set(
                RunnerLimitFieldV2::RepairActionEncodedBytes,
                RunnerLimitValueV2::U64(1),
            )?;
        }
        _ => {}
    }
    Ok(())
}

fn set_production_diagnostic_chain_v1(
    set: &mut impl FnMut(RunnerLimitFieldV2, RunnerLimitValueV2) -> Result<(), ConstructionErrorV2>,
    lifecycle_record: bool,
    actionable: bool,
) -> Result<(), ConstructionErrorV2> {
    set(
        RunnerLimitFieldV2::RepairActionEncodedBytes,
        RunnerLimitValueV2::U64(1),
    )?;
    if actionable {
        set(
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            RunnerLimitValueV2::U64(1),
        )?;
    }
    if lifecycle_record {
        set(
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            RunnerLimitValueV2::U64(1),
        )?;
    }
    Ok(())
}

fn limit_violation_reason_v1(violation: &RunnerLimitsViolationV2) -> RunnerV2RawReasonV1 {
    match violation.kind() {
        RunnerLimitsViolationKindV2::WrongWidth => RunnerV2RawReasonV1::WrongPrimitiveWidth,
        RunnerLimitsViolationKindV2::ExceedsBaseCeiling => RunnerV2RawReasonV1::AboveProfileCeiling,
        RunnerLimitsViolationKindV2::FixedFieldChanged => {
            RunnerV2RawReasonV1::FixedRepresentationChanged
        }
        RunnerLimitsViolationKindV2::BelowStructuralMinimum
        | RunnerLimitsViolationKindV2::DeclaredMinimumUnmet
        | RunnerLimitsViolationKindV2::ExecutableCaseSetEmpty
        | RunnerLimitsViolationKindV2::CaseCountExceeded
        | RunnerLimitsViolationKindV2::FamilyRowsExceeded
        | RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient => {
            RunnerV2RawReasonV1::BelowStructuralMinimum
        }
        RunnerLimitsViolationKindV2::ArithmeticOverflow => {
            RunnerV2RawReasonV1::CheckedRepresentationalOverflow
        }
        RunnerLimitsViolationKindV2::JointFeasibilityViolation
        | RunnerLimitsViolationKindV2::ProtocolStoredLengthMismatch
        | RunnerLimitsViolationKindV2::EnvelopeOverheadExceeded
        | RunnerLimitsViolationKindV2::ArtifactCountExceeded
        | RunnerLimitsViolationKindV2::SystemObjectSetMismatch
        | RunnerLimitsViolationKindV2::AggregateMismatch => {
            RunnerV2RawReasonV1::JointFeasibilityViolation
        }
        RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder
        | RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum
        | RunnerLimitsViolationKindV2::NonExecutableCaseSetPresent => {
            RunnerV2RawReasonV1::MalformedOrNoncanonical
        }
    }
}

fn append_limit_violation_observations_v1(
    numeric: &mut Vec<RunnerV2SafeNumericObservationV1>,
    violation: &RunnerLimitsViolationV2,
) -> Result<(), ConstructionErrorV2> {
    numeric.push(RunnerV2SafeNumericObservationV1::count(
        stage_a_token(
            "runner_v2.base_values.observation.name",
            "violation-field-ordinal",
        )?,
        u64::from(violation.field().ordinal()),
    ));
    numeric.push(RunnerV2SafeNumericObservationV1::limit(
        stage_a_token(
            "runner_v2.base_values.observation.name",
            "violation-observed",
        )?,
        violation.observed(),
        violation.unit(),
    ));
    match violation.expected() {
        RunnerLimitExpectationV2::AtMost(value)
        | RunnerLimitExpectationV2::AtLeast(value)
        | RunnerLimitExpectationV2::Exactly(value) => {
            numeric.push(RunnerV2SafeNumericObservationV1::limit(
                stage_a_token("runner_v2.base_values.observation.name", "expected-bound")?,
                value,
                violation.unit(),
            ));
        }
        RunnerLimitExpectationV2::Width(width) => {
            numeric.push(RunnerV2SafeNumericObservationV1::count(
                stage_a_token(
                    "runner_v2.base_values.observation.name",
                    "expected-width-bits",
                )?,
                match width {
                    RunnerLimitWidthV2::U32 => 32,
                    RunnerLimitWidthV2::U64 => 64,
                },
            ));
        }
        RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
            numeric.push(RunnerV2SafeNumericObservationV1::count(
                stage_a_token(
                    "runner_v2.base_values.observation.name",
                    "expected-minimum-ordinal-step",
                )?,
                1,
            ));
        }
    }
    Ok(())
}

fn raw_limit_diagnostic_v1(
    violation: &RunnerLimitsViolationV2,
    reason: RunnerV2RawReasonV1,
) -> Result<RunnerV2RawDiagnosticV1, ConstructionErrorV2> {
    let repair = RunnerV2RawRepairV1::new(
        violation.repair_rank(),
        violation.repair_kind(),
        stage_a_token(
            "runner_v2.base_values.diagnostic.repair_target",
            violation.repair_target(),
        )?,
    );
    RunnerV2RawDiagnosticV1::new(
        DiagnosticCodeV2::RunnerRefused,
        stage_a_token("runner_v2.base_values.diagnostic.owner", violation.owner())?,
        RetryabilityV2::AfterInputChange,
        Vec::new(),
        if matches!(reason, RunnerV2RawReasonV1::CheckedRepresentationalOverflow) {
            Vec::new()
        } else {
            vec![repair]
        },
    )
}

fn raw_diagnostic_v1(
    outcome: RunnerV2RawOutcomeKindV1,
    reason: RunnerV2RawReasonV1,
    owner: &'static str,
    field: Option<RunnerLimitFieldV2>,
) -> Result<RunnerV2RawDiagnosticV1, ConstructionErrorV2> {
    let (code, retryability) = match outcome {
        RunnerV2RawOutcomeKindV1::Accepted => {
            return Err(stage_a_error(
                ConstructionErrorKindV2::Unexpected,
                "runner_v2.base_values.diagnostic",
                "no diagnostic for an accepted raw cell",
                outcome.code(),
            ));
        }
        RunnerV2RawOutcomeKindV1::Refused => (
            DiagnosticCodeV2::RunnerRefused,
            RetryabilityV2::AfterInputChange,
        ),
        RunnerV2RawOutcomeKindV1::Failed => {
            (DiagnosticCodeV2::RunnerInternalError, RetryabilityV2::Never)
        }
        RunnerV2RawOutcomeKindV1::Unsupported => (
            DiagnosticCodeV2::RunnerUnsupported,
            RetryabilityV2::AfterPrerequisiteChange,
        ),
        RunnerV2RawOutcomeKindV1::Inapplicable => {
            (DiagnosticCodeV2::RunnerNotRun, RetryabilityV2::Never)
        }
    };
    let repairs = if let Some(field) =
        field.filter(|_| matches!(outcome, RunnerV2RawOutcomeKindV1::Refused))
    {
        vec![RunnerV2RawRepairV1::new(
            1,
            RepairActionKindV2::ChangeArguments,
            stage_a_token(
                "runner_v2.base_values.diagnostic.repair_target",
                field.descriptor().name,
            )?,
        )]
    } else {
        Vec::new()
    };
    let prerequisites = if matches!(
        outcome,
        RunnerV2RawOutcomeKindV1::Inapplicable | RunnerV2RawOutcomeKindV1::Unsupported
    ) {
        vec![stage_a_token(
            "runner_v2.base_values.diagnostic.prerequisite",
            match (outcome, reason) {
                (
                    RunnerV2RawOutcomeKindV1::Inapplicable,
                    RunnerV2RawReasonV1::ShardInapplicable,
                ) => "registered-shard-policy",
                (
                    RunnerV2RawOutcomeKindV1::Inapplicable,
                    RunnerV2RawReasonV1::ResumeInapplicable,
                ) => "registered-resume-policy",
                (RunnerV2RawOutcomeKindV1::Unsupported, _) => "registered-evaluator-implementation",
                _ => "runtime-operation-applicability",
            },
        )?]
    } else {
        Vec::new()
    };
    RunnerV2RawDiagnosticV1::new(
        code,
        stage_a_token("runner_v2.base_values.diagnostic.owner", owner)?,
        retryability,
        prerequisites,
        repairs,
    )
}

fn evaluate_meta_operation_v1(
    cell_id: StableTokenV2,
    operation: RunnerV2StageAMetaOperationV1,
) -> Result<RunnerV2RawCellObservationV1, ConstructionErrorV2> {
    let (outcome, reason, observed_count) = match operation {
        RunnerV2StageAMetaOperationV1::TypedAbsenceDistinctFromZero => {
            let absent: TypedOptionV1<DigestValueV2> = TypedOptionV1::Absent;
            let zero = DigestValueV2::from_array(
                DigestRoleV2::Source,
                SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
                [0_u8; 32],
            );
            let present = TypedOptionV1::Present(zero);
            bool_meta_result_v1(absent.wire_tag() == 0 && present.wire_tag() == 1)
        }
        RunnerV2StageAMetaOperationV1::F32NamedTotalOrder => bool_meta_result_v1(
            F32BitsV2::from_bits(0x8000_0000).ieee_total_cmp_v1(F32BitsV2::from_bits(0))
                == Ordering::Less
                && F32BitsV2::from_bits(0x7fc0_0001)
                    .ieee_total_cmp_v1(F32BitsV2::from_bits(0x7fc0_0002))
                    != Ordering::Equal,
        ),
        RunnerV2StageAMetaOperationV1::F64NamedTotalOrder => bool_meta_result_v1(
            F64BitsV2::from_bits(0x8000_0000_0000_0000).ieee_total_cmp_v1(F64BitsV2::from_bits(0))
                == Ordering::Less
                && F64BitsV2::from_bits(0x7ff8_0000_0000_0001)
                    .ieee_total_cmp_v1(F64BitsV2::from_bits(0x7ff8_0000_0000_0002))
                    != Ordering::Equal,
        ),
        RunnerV2StageAMetaOperationV1::CapabilityNoneContract => {
            let (registry, profile_registry, contract) = build_capability_none_v1()?;
            bool_meta_result_v1(
                registry.rows().len() == 5
                    && profile_registry.rows().len() == 5
                    && contract.profile() == BaseCoverageCloseCapabilityProfileV1::None
                    && contract.required().is_empty()
                    && contract.permitted().is_empty(),
            )
        }
        RunnerV2StageAMetaOperationV1::CommonRequirementsExact => {
            let rows = build_common_requirements_v1()?;
            let accepted = validate_common_requirements_exact_v1(&rows).is_ok();
            (
                if accepted {
                    RunnerV2RawOutcomeKindV1::Accepted
                } else {
                    RunnerV2RawOutcomeKindV1::Failed
                },
                if accepted {
                    RunnerV2RawReasonV1::ExactCheckedValue
                } else {
                    RunnerV2RawReasonV1::InternalInvariantFailure
                },
                rows.len(),
            )
        }
        RunnerV2StageAMetaOperationV1::CommonRequirementReorderedRefusal => {
            let mut rows = build_common_requirements_v1()?;
            rows.swap(0, 1);
            let refused = validate_common_requirements_exact_v1(&rows).is_err();
            (
                if refused {
                    RunnerV2RawOutcomeKindV1::Refused
                } else {
                    RunnerV2RawOutcomeKindV1::Failed
                },
                if refused {
                    RunnerV2RawReasonV1::ExactMembershipMismatch
                } else {
                    RunnerV2RawReasonV1::InternalInvariantFailure
                },
                rows.len(),
            )
        }
        RunnerV2StageAMetaOperationV1::FutureSourcesExact => {
            let rows = build_future_sources_v1()?;
            let accepted = validate_future_sources_exact_v1(&rows).is_ok();
            (
                if accepted {
                    RunnerV2RawOutcomeKindV1::Accepted
                } else {
                    RunnerV2RawOutcomeKindV1::Failed
                },
                if accepted {
                    RunnerV2RawReasonV1::ExactCheckedValue
                } else {
                    RunnerV2RawReasonV1::InternalInvariantFailure
                },
                rows.len(),
            )
        }
        RunnerV2StageAMetaOperationV1::RootlessAc58 => {
            let row = build_ac58_v1()?;
            bool_meta_result_v1(
                row.disposition == CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame
                    && row.authority_surfaces_present_empty,
            )
        }
        RunnerV2StageAMetaOperationV1::OwnerSourceFragment => {
            let rows = build_owner_source_fragment_v1()?;
            let accepted = validate_owner_source_fragment_exact_v1(&rows).is_ok();
            (
                if accepted {
                    RunnerV2RawOutcomeKindV1::Accepted
                } else {
                    RunnerV2RawOutcomeKindV1::Failed
                },
                if accepted {
                    RunnerV2RawReasonV1::ExactCheckedValue
                } else {
                    RunnerV2RawReasonV1::InternalInvariantFailure
                },
                rows.len(),
            )
        }
        RunnerV2StageAMetaOperationV1::LocalRoute => {
            let route = build_route_v1()?;
            bool_meta_result_v1(
                route.route_id.as_str() == RUNNER_V2_BASE_VALUES_LOCAL_ROUTE_ID_V1
                    && route.class == RunnerV2LocalRouteClassV1::LocalOnly
                    && route.execution_owner.as_str()
                        == "frankensim-epic-foundations-huq.24.1.1.1.7"
                    && route.capability_profile == BaseCoverageCloseCapabilityProfileV1::None
                    && matches!(route.external_driver, TypedOptionV1::Absent),
            )
        }
        RunnerV2StageAMetaOperationV1::DiagnosticRedaction => {
            let forbidden_value =
                String::from("runner-v2-sensitive-canary-value-must-never-be-retained");
            let error = consume_and_redact_stage_a_value_v1(forbidden_value);
            let display = format!("{error}");
            let debug = format!("{error:?}");
            bool_meta_result_v1(
                error.observed() == "<redacted:sensitive-or-ambient>"
                    && !error.observed().contains("runner-v2-sensitive-canary")
                    && !display.contains("runner-v2-sensitive-canary")
                    && !debug.contains("runner-v2-sensitive-canary"),
            )
        }
        RunnerV2StageAMetaOperationV1::ReproductionDeclaration => {
            let requirements = build_common_requirements_v1()?;
            bool_meta_result_v1(requirements.iter().any(|row| {
                row.slot_id.as_str() == "runner-v2.common.reproduction-schema.v1"
                    && matches!(row.future_root, TypedOptionV1::Absent)
            }))
        }
        RunnerV2StageAMetaOperationV1::CompileFailOrderingSurface => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::PureDeclarationFacet,
            1,
        ),
        RunnerV2StageAMetaOperationV1::ShardInapplicable => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::ShardInapplicable,
            1,
        ),
        RunnerV2StageAMetaOperationV1::ResumeInapplicable => (
            RunnerV2RawOutcomeKindV1::Inapplicable,
            RunnerV2RawReasonV1::ResumeInapplicable,
            1,
        ),
    };
    let diagnostic = if matches!(outcome, RunnerV2RawOutcomeKindV1::Accepted) {
        None
    } else {
        Some(raw_diagnostic_v1(
            outcome,
            reason,
            "fs-evidence-runner.runner-v2.base-values",
            None,
        )?)
    };
    RunnerV2RawCellObservationV1::new(
        cell_id,
        outcome,
        reason,
        vec![RunnerV2SafeNumericObservationV1::count(
            stage_a_token("runner_v2.base_values.observation.name", "observed-count")?,
            u64::try_from(observed_count).map_err(|_| {
                stage_a_error(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.base_values.observation.count",
                    "observed count representable as u64",
                    observed_count,
                )
            })?,
        )],
        diagnostic,
    )
}

fn consume_and_redact_stage_a_value_v1(forbidden_value: String) -> ConstructionErrorV2 {
    drop(forbidden_value);
    ConstructionErrorV2::new_redacted(
        ConstructionErrorKindV2::Incompatible,
        "runner_v2.base_values.redaction",
        "one redacted semantic class",
        crate::construction::ConstructionObservedDataClassV2::SensitiveOrAmbient,
    )
}

const fn bool_meta_result_v1(
    accepted: bool,
) -> (RunnerV2RawOutcomeKindV1, RunnerV2RawReasonV1, usize) {
    if accepted {
        (
            RunnerV2RawOutcomeKindV1::Accepted,
            RunnerV2RawReasonV1::ExactCheckedValue,
            1,
        )
    } else {
        (
            RunnerV2RawOutcomeKindV1::Failed,
            RunnerV2RawReasonV1::InternalInvariantFailure,
            0,
        )
    }
}

fn build_retained_domain_obligations_v1()
-> Result<Vec<RunnerV2RetainedDomainObligationV1>, ConstructionErrorV2> {
    RETAINED_DOMAIN_OBLIGATION_DEFINITIONS_V1
        .iter()
        .enumerate()
        .map(|(index, (suffix, facet))| {
            let ordinal = u16::try_from(index + 1).map_err(|_| {
                stage_a_error(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.base_values.retained_domain.ordinal",
                    "one-based retained-domain ordinal representable as u16",
                    index,
                )
            })?;
            Ok(RunnerV2RetainedDomainObligationV1 {
                ordinal,
                stable_id: stage_a_token_owned(
                    "runner_v2.base_values.retained_domain.stable_id",
                    format!(
                        "runner-v2.base-values.retained-{:03}-{}.v1",
                        ordinal, suffix
                    ),
                )?,
                facet: *facet,
            })
        })
        .collect()
}

fn retained_domain_inventory_ids_v1() -> Vec<String> {
    RETAINED_DOMAIN_OBLIGATION_DEFINITIONS_V1
        .iter()
        .enumerate()
        .map(|(index, (suffix, _))| {
            format!(
                "runner-v2.base-values.retained-{:03}-{}.v1",
                index + 1,
                suffix
            )
        })
        .collect()
}

fn validate_stage_a_retained_domain_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let expected = retained_domain_inventory_ids_v1();
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.retained_domain",
        &expected,
        rows,
        &[],
        "restore-exact-retained-domain-obligation-catalog",
    )
}

fn validate_retained_domain_obligations_diagnostic_v1(
    rows: &[RunnerV2RetainedDomainObligationV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| row.stable_id.as_str().to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_retained_domain_inventory_v1(&observed_ids)?;

    for (index, (row, (suffix, facet))) in rows
        .iter()
        .zip(RETAINED_DOMAIN_OBLIGATION_DEFINITIONS_V1)
        .enumerate()
    {
        let expected_ordinal = u16::try_from(index + 1).map_err(|_| {
            stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "runner_v2.base_values.retained_domain.ordinal",
                index,
                format!("retained-domain-row-{}", index + 1),
                row.stable_id.as_str().to_owned(),
                false,
                "ordinal",
                "one-based ordinal representable as u16".to_owned(),
                (index + 1).to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-retained-domain-obligation-catalog",
            )
        })?;
        let expected_id = format!(
            "runner-v2.base-values.retained-{:03}-{}.v1",
            expected_ordinal, suffix
        );
        if row.ordinal != expected_ordinal {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.retained_domain",
                index,
                expected_id.clone(),
                row.stable_id.as_str().to_owned(),
                false,
                "ordinal",
                expected_ordinal.to_string(),
                row.ordinal.to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-retained-domain-obligation-catalog",
            ));
        }
        if row.facet != facet {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.retained_domain",
                index,
                expected_id,
                row.stable_id.as_str().to_owned(),
                false,
                "facet",
                facet.code().to_string(),
                row.facet.code().to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-retained-domain-obligation-catalog",
            ));
        }
    }
    Ok(())
}

fn validate_retained_domain_obligations_exact_v1(
    rows: &[RunnerV2RetainedDomainObligationV1],
) -> Result<(), ConstructionErrorV2> {
    validate_retained_domain_obligations_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn build_limit_mutation_obligations_v1()
-> Result<Vec<RunnerV2LimitMutationObligationV1>, ConstructionErrorV2> {
    let diagnostic_owner = stage_a_token(
        "runner_v2.base_values.limit_mutation.diagnostic_owner",
        "fs-evidence-runner.runner-limits",
    )?;
    STAGE_A_LIMIT_LITERALS_V1
        .into_iter()
        .map(|literal| {
            let field_name = stage_a_token(
                "runner_v2.base_values.limit_mutation.field_name",
                independent_limit_name_v1(literal.field),
            )?;
            Ok(RunnerV2LimitMutationObligationV1 {
                ordinal: literal.field.ordinal(),
                stable_id: stage_a_token_owned(
                    "runner_v2.base_values.limit_mutation.stable_id",
                    format!(
                        "runner-v2.base-values.limit-{:03}.wrong-width-mutation.v1",
                        literal.field.ordinal()
                    ),
                )?,
                field: literal.field,
                field_name: field_name.clone(),
                declared_width: literal.width,
                opposite_width_zero: opposite_width_zero_v1(literal.width),
                unit: literal.unit,
                expected_reason: RunnerV2RawReasonV1::WrongPrimitiveWidth,
                diagnostic_owner: diagnostic_owner.clone(),
                repair_rank: 1,
                repair_kind: RepairActionKindV2::UpdatePolicyOrCapability,
                repair_target: field_name,
            })
        })
        .collect()
}

fn limit_mutation_inventory_ids_v1() -> Vec<String> {
    STAGE_A_LIMIT_LITERALS_V1
        .iter()
        .map(|literal| {
            format!(
                "runner-v2.base-values.limit-{:03}.wrong-width-mutation.v1",
                literal.field.ordinal()
            )
        })
        .collect()
}

fn validate_stage_a_limit_mutation_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let expected = limit_mutation_inventory_ids_v1();
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.limit_mutations",
        &expected,
        rows,
        &[],
        "restore-exact-limit-wrong-width-mutation-catalog",
    )
}

fn validate_limit_mutation_obligations_diagnostic_v1(
    rows: &[RunnerV2LimitMutationObligationV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| row.stable_id.as_str().to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_limit_mutation_inventory_v1(&observed_ids)?;

    for (index, (row, literal)) in rows.iter().zip(STAGE_A_LIMIT_LITERALS_V1).enumerate() {
        let expected_id = format!(
            "runner-v2.base-values.limit-{:03}.wrong-width-mutation.v1",
            literal.field.ordinal()
        );
        let expected_name = independent_limit_name_v1(literal.field);
        let mismatch = |component: &'static str,
                        expected_semantic_value: String,
                        observed_safe_value: String| {
            stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.limit_mutations",
                index,
                expected_id.clone(),
                row.stable_id.as_str().to_owned(),
                false,
                component,
                expected_semantic_value,
                observed_safe_value,
                false,
                rows.len(),
                rows.len(),
                "restore-exact-limit-wrong-width-mutation-catalog",
            )
        };
        if row.ordinal != literal.field.ordinal() {
            return Err(mismatch(
                "ordinal",
                literal.field.ordinal().to_string(),
                row.ordinal.to_string(),
            ));
        }
        if row.field != literal.field {
            return Err(mismatch(
                "field",
                format!("{}:{expected_name}", literal.field.ordinal()),
                format!("{}:{}", row.field.ordinal(), row.field.descriptor().name),
            ));
        }
        if row.field_name.as_str() != expected_name {
            return Err(mismatch(
                "field-name",
                expected_name.to_owned(),
                row.field_name.as_str().to_owned(),
            ));
        }
        if row.declared_width != literal.width {
            return Err(mismatch(
                "declared-width",
                limit_width_name_v1(literal.width).to_owned(),
                limit_width_name_v1(row.declared_width).to_owned(),
            ));
        }
        let expected_opposite_zero = opposite_width_zero_v1(literal.width);
        if row.opposite_width_zero != expected_opposite_zero {
            return Err(mismatch(
                "opposite-width-zero",
                limit_value_safe_name_v1(expected_opposite_zero),
                limit_value_safe_name_v1(row.opposite_width_zero),
            ));
        }
        if row.unit != literal.unit {
            return Err(mismatch(
                "unit",
                independent_limit_unit_name_v1(literal.unit).to_owned(),
                independent_limit_unit_name_v1(row.unit).to_owned(),
            ));
        }
        if row.expected_reason != RunnerV2RawReasonV1::WrongPrimitiveWidth {
            return Err(mismatch(
                "expected-reason",
                format!(
                    "{:02}:{}",
                    RunnerV2RawReasonV1::WrongPrimitiveWidth.code(),
                    runner_v2_raw_reason_name_v1(RunnerV2RawReasonV1::WrongPrimitiveWidth)
                ),
                format!(
                    "{:02}:{}",
                    row.expected_reason.code(),
                    runner_v2_raw_reason_name_v1(row.expected_reason)
                ),
            ));
        }
        if row.diagnostic_owner.as_str() != "fs-evidence-runner.runner-limits" {
            return Err(mismatch(
                "diagnostic-owner",
                "fs-evidence-runner.runner-limits".to_owned(),
                row.diagnostic_owner.as_str().to_owned(),
            ));
        }
        if row.repair_rank != 1 {
            return Err(mismatch(
                "repair-rank",
                "1".to_owned(),
                row.repair_rank.to_string(),
            ));
        }
        if row.repair_kind != RepairActionKindV2::UpdatePolicyOrCapability {
            return Err(mismatch(
                "repair-kind",
                RepairActionKindV2::UpdatePolicyOrCapability
                    .code()
                    .to_string(),
                row.repair_kind.code().to_string(),
            ));
        }
        if row.repair_target.as_str() != expected_name {
            return Err(mismatch(
                "repair-target",
                expected_name.to_owned(),
                row.repair_target.as_str().to_owned(),
            ));
        }
    }
    Ok(())
}

fn validate_limit_mutation_obligations_exact_v1(
    rows: &[RunnerV2LimitMutationObligationV1],
) -> Result<(), ConstructionErrorV2> {
    validate_limit_mutation_obligations_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn build_common_requirements_v1()
-> Result<Vec<RunnerV2CommonContractRequirementV1>, ConstructionErrorV2> {
    let semantic_owner = stage_a_token(
        "runner_v2.base_values.requirement.semantic_owner",
        "frankensim-epic-foundations-huq.24.1.1.1",
    )?;
    let resolution_owner = stage_a_token(
        "runner_v2.base_values.requirement.resolution_owner",
        "frankensim-epic-foundations-huq.24.1.1.1.7",
    )?;
    let no_claim = stage_a_token(
        "runner_v2.base_values.requirement.no_claim",
        RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
    )?;
    COMMON_REQUIREMENT_DEFINITIONS_V1
        .iter()
        .enumerate()
        .map(|(index, definition)| {
            Ok(RunnerV2CommonContractRequirementV1 {
                ordinal: u16::try_from(index + 1).map_err(|_| {
                    stage_a_error(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "runner_v2.base_values.requirement.ordinal",
                        "one-based requirement ordinal representable as u16",
                        index,
                    )
                })?,
                slot_id: stage_a_token(
                    "runner_v2.base_values.requirement.slot_id",
                    definition.slot_id,
                )?,
                api_generation: RUNNER_SPEC_V2_API_GENERATION,
                wire_version: RUNNER_V2_WIRE_VERSION,
                predecessor_policy: RUNNER_V2_PREDECESSOR_POLICY,
                semantic_owner: semantic_owner.clone(),
                realization_owner: stage_a_token(
                    "runner_v2.base_values.requirement.realization_owner",
                    definition.realization_owner,
                )?,
                future_nominal_role: stage_a_token(
                    "runner_v2.base_values.requirement.future_nominal_role",
                    definition.role,
                )?,
                future_domain: stage_a_token(
                    "runner_v2.base_values.requirement.future_domain",
                    definition.domain,
                )?,
                included_planes: RunnerV2ContractPlaneSetV1::from_mask(definition.planes),
                fulfillment_stage: definition.stage,
                resolution_owner: resolution_owner.clone(),
                future_root: TypedOptionV1::Absent,
                no_claim: no_claim.clone(),
            })
        })
        .collect()
}

fn common_requirement_inventory_ids_v1() -> Vec<String> {
    COMMON_REQUIREMENT_DEFINITIONS_V1
        .iter()
        .map(|definition| definition.slot_id.to_owned())
        .collect()
}

fn validate_stage_a_common_requirement_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.requirements",
        &common_requirement_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-common-requirement-catalog",
    )
}

fn common_requirement_semantic_mismatch_v1(
    rows: &[RunnerV2CommonContractRequirementV1],
    index: usize,
    component: &'static str,
    expected: String,
    observed: String,
) -> RunnerV2StageAInventoryMismatchV1 {
    let identity = COMMON_REQUIREMENT_DEFINITIONS_V1[index].slot_id.to_owned();
    stage_a_inventory_mismatch_v1(
        ConstructionErrorKindV2::Incompatible,
        "runner_v2.base_values.requirements",
        index,
        identity.clone(),
        identity,
        false,
        component,
        expected,
        observed,
        false,
        rows.len(),
        rows.len(),
        "restore-exact-common-requirement-catalog",
    )
}

fn validate_common_requirements_diagnostic_v1(
    rows: &[RunnerV2CommonContractRequirementV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| row.slot_id.as_str().to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_common_requirement_inventory_v1(&observed_ids)?;

    for (index, (row, definition)) in rows
        .iter()
        .zip(COMMON_REQUIREMENT_DEFINITIONS_V1)
        .enumerate()
    {
        let expected_ordinal = u16::try_from(index + 1).map_err(|_| {
            stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "runner_v2.base_values.requirements",
                index,
                definition.slot_id.to_owned(),
                row.slot_id.as_str().to_owned(),
                false,
                "ordinal",
                "one-based requirement ordinal representable as u16".to_owned(),
                (index + 1).to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-common-requirement-catalog",
            )
        })?;
        if row.ordinal != expected_ordinal {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "ordinal",
                expected_ordinal.to_string(),
                row.ordinal.to_string(),
            ));
        }
        if row.api_generation != RUNNER_SPEC_V2_API_GENERATION {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "api-generation",
                format!("{:?}", RUNNER_SPEC_V2_API_GENERATION),
                format!("{:?}", row.api_generation),
            ));
        }
        if row.wire_version != RUNNER_V2_WIRE_VERSION {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "wire-version",
                format!("{:?}", RUNNER_V2_WIRE_VERSION),
                format!("{:?}", row.wire_version),
            ));
        }
        if row.predecessor_policy != RUNNER_V2_PREDECESSOR_POLICY {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "predecessor-policy",
                format!("{:?}", RUNNER_V2_PREDECESSOR_POLICY),
                format!("{:?}", row.predecessor_policy),
            ));
        }
        if row.semantic_owner.as_str() != STAGE_A_INVENTORY_SEMANTIC_OWNER_V1 {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "semantic-owner",
                STAGE_A_INVENTORY_SEMANTIC_OWNER_V1.to_owned(),
                row.semantic_owner.as_str().to_owned(),
            ));
        }
        if row.realization_owner.as_str() != definition.realization_owner {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "realization-owner",
                definition.realization_owner.to_owned(),
                row.realization_owner.as_str().to_owned(),
            ));
        }
        if row.future_nominal_role.as_str() != definition.role {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "future-nominal-role",
                definition.role.to_owned(),
                row.future_nominal_role.as_str().to_owned(),
            ));
        }
        if row.future_domain.as_str() != definition.domain {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "future-domain",
                definition.domain.to_owned(),
                row.future_domain.as_str().to_owned(),
            ));
        }
        if row.included_planes.mask() != definition.planes
            || definition.planes == 0
            || definition.planes & !(C | E | R) != 0
        {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "included-planes",
                format!("0b{:03b}", definition.planes),
                format!("0b{:03b}", row.included_planes.mask()),
            ));
        }
        if row.fulfillment_stage != definition.stage {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "fulfillment-stage",
                definition.stage.code().to_string(),
                row.fulfillment_stage.code().to_string(),
            ));
        }
        if row.resolution_owner.as_str() != "frankensim-epic-foundations-huq.24.1.1.1.7" {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "resolution-owner",
                "frankensim-epic-foundations-huq.24.1.1.1.7".to_owned(),
                row.resolution_owner.as_str().to_owned(),
            ));
        }
        if !matches!(row.future_root, TypedOptionV1::Absent) {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "future-root-presence",
                "absent".to_owned(),
                "present".to_owned(),
            ));
        }
        if row.no_claim.as_str() != RUNNER_V2_BASE_VALUES_NO_CLAIM_V1 {
            return Err(common_requirement_semantic_mismatch_v1(
                rows,
                index,
                "no-claim",
                RUNNER_V2_BASE_VALUES_NO_CLAIM_V1.to_owned(),
                row.no_claim.as_str().to_owned(),
            ));
        }
    }
    Ok(())
}

fn validate_common_requirements_exact_v1(
    rows: &[RunnerV2CommonContractRequirementV1],
) -> Result<(), ConstructionErrorV2> {
    validate_common_requirements_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn build_future_sources_v1() -> Result<Vec<RunnerV2FutureSourceRequirementV1>, ConstructionErrorV2>
{
    FUTURE_SOURCE_PATHS_V1
        .iter()
        .enumerate()
        .map(|(index, path)| {
            Ok(RunnerV2FutureSourceRequirementV1 {
                final_ordinal: u16::try_from(RUNNER_V2_EXISTING_SOURCE_COUNT_V1 + index + 1)
                    .map_err(|_| {
                        stage_a_error(
                            ConstructionErrorKindV2::ArithmeticOverflow,
                            "runner_v2.base_values.future_source.ordinal",
                            "one final source ordinal representable as u16",
                            index,
                        )
                    })?,
                path: stage_a_path("runner_v2.base_values.future_source.path", path)?,
                future_content_root: TypedOptionV1::Absent,
            })
        })
        .collect()
}

fn future_source_inventory_ids_v1() -> Vec<String> {
    FUTURE_SOURCE_PATHS_V1
        .iter()
        .map(|path| (*path).to_owned())
        .collect()
}

fn validate_stage_a_future_source_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.future_sources",
        &future_source_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-future-source-catalog",
    )
}

fn future_source_semantic_mismatch_v1(
    rows: &[RunnerV2FutureSourceRequirementV1],
    index: usize,
    component: &'static str,
    expected: String,
    observed: String,
) -> RunnerV2StageAInventoryMismatchV1 {
    let identity = FUTURE_SOURCE_PATHS_V1[index].to_owned();
    stage_a_inventory_mismatch_v1(
        ConstructionErrorKindV2::Incompatible,
        "runner_v2.base_values.future_sources",
        index,
        identity.clone(),
        identity,
        false,
        component,
        expected,
        observed,
        false,
        rows.len(),
        rows.len(),
        "restore-exact-future-source-catalog",
    )
}

fn validate_future_sources_diagnostic_v1(
    rows: &[RunnerV2FutureSourceRequirementV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| row.path.as_str().to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_future_source_inventory_v1(&observed_ids)?;

    for (index, (row, expected_path)) in rows.iter().zip(FUTURE_SOURCE_PATHS_V1).enumerate() {
        let expected_ordinal = u16::try_from(RUNNER_V2_EXISTING_SOURCE_COUNT_V1 + index + 1)
            .map_err(|_| {
                stage_a_inventory_mismatch_v1(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.base_values.future_sources",
                    index,
                    expected_path.to_owned(),
                    row.path.as_str().to_owned(),
                    false,
                    "final-ordinal",
                    "one final source ordinal representable as u16".to_owned(),
                    (RUNNER_V2_EXISTING_SOURCE_COUNT_V1 + index + 1).to_string(),
                    false,
                    rows.len(),
                    rows.len(),
                    "restore-exact-future-source-catalog",
                )
            })?;
        if row.final_ordinal != expected_ordinal {
            return Err(future_source_semantic_mismatch_v1(
                rows,
                index,
                "final-ordinal",
                expected_ordinal.to_string(),
                row.final_ordinal.to_string(),
            ));
        }
        if !matches!(row.future_content_root, TypedOptionV1::Absent) {
            return Err(future_source_semantic_mismatch_v1(
                rows,
                index,
                "future-content-root-presence",
                "absent".to_owned(),
                "present".to_owned(),
            ));
        }
    }
    Ok(())
}

fn validate_future_sources_exact_v1(
    rows: &[RunnerV2FutureSourceRequirementV1],
) -> Result<(), ConstructionErrorV2> {
    validate_future_sources_diagnostic_v1(rows).map_err(|mismatch| mismatch.as_construction_error())
}

fn build_owner_source_fragment_v1() -> Result<Vec<RunnerV2OwnerSourceMemberV1>, ConstructionErrorV2>
{
    OWNER_SOURCE_DECLARATIONS_V1
        .into_iter()
        .map(|declaration| {
            Ok(RunnerV2OwnerSourceMemberV1 {
                path: stage_a_path("runner_v2.base_values.owner_source.path", declaration.path)?,
                content_root: RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                    "org.frankensim.fs-evidence-runner.runner-v2.stage-a.source-member.v1",
                    declaration.bytes,
                )),
            })
        })
        .collect()
}

fn independently_expected_owner_source_root_v1(
    index: usize,
) -> Option<RunnerV2StageASourceMemberRootV1> {
    let bytes: &[u8] = match index {
        0 => include_bytes!("../handoff.rs"),
        1 => include_bytes!("base_values.rs"),
        _ => return None,
    };
    Some(RunnerV2StageASourceMemberRootV1::from_content_hash(
        hash_domain(
            "org.frankensim.fs-evidence-runner.runner-v2.stage-a.source-member.v1",
            bytes,
        ),
    ))
}

fn source_root_safe_name_v1(root: RunnerV2StageASourceMemberRootV1) -> String {
    format!("{:02x?}", root.bytes())
}

fn validate_source_members_diagnostic_v1(
    inventory: &'static str,
    repair_target: &'static str,
    expected_paths: &[&str],
    observed_paths: &[String],
    observed_roots: &[RunnerV2StageASourceMemberRootV1],
    expected_root: impl Fn(usize) -> Option<RunnerV2StageASourceMemberRootV1>,
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let expected_ids = expected_paths
        .iter()
        .map(|path| (*path).to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_exact_inventory_v1(
        inventory,
        &expected_ids,
        observed_paths,
        &[],
        repair_target,
    )?;
    for (index, ((expected_path, observed_path), observed_root)) in expected_paths
        .iter()
        .zip(observed_paths)
        .zip(observed_roots)
        .enumerate()
    {
        let Some(expected_root) = expected_root(index) else {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                inventory,
                index,
                (*expected_path).to_owned(),
                observed_path.clone(),
                false,
                "content-root-oracle",
                "independent-content-root-present".to_owned(),
                "independent-content-root-missing".to_owned(),
                false,
                expected_paths.len(),
                observed_paths.len(),
                repair_target,
            ));
        };
        if *observed_root != expected_root {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                inventory,
                index,
                (*expected_path).to_owned(),
                observed_path.clone(),
                false,
                "content-root",
                source_root_safe_name_v1(expected_root),
                source_root_safe_name_v1(*observed_root),
                false,
                expected_paths.len(),
                observed_paths.len(),
                repair_target,
            ));
        }
    }
    Ok(())
}

fn validate_owner_source_fragment_diagnostic_v1(
    rows: &[RunnerV2OwnerSourceMemberV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let paths = rows
        .iter()
        .map(|row| row.path.as_str().to_owned())
        .collect::<Vec<_>>();
    let roots = rows.iter().map(|row| row.content_root).collect::<Vec<_>>();
    validate_source_members_diagnostic_v1(
        "runner_v2.base_values.owner_sources",
        "restore-exact-owner-source-fragment",
        &OWNER_SOURCE_PATHS_V1,
        &paths,
        &roots,
        independently_expected_owner_source_root_v1,
    )
}

fn validate_owner_source_fragment_exact_v1(
    rows: &[RunnerV2OwnerSourceMemberV1],
) -> Result<(), ConstructionErrorV2> {
    validate_owner_source_fragment_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn build_dependency_source_closure_v1()
-> Result<Vec<RunnerV2DependencySourceMemberV1>, ConstructionErrorV2> {
    DEPENDENCY_SOURCE_DECLARATIONS_V1
        .into_iter()
        .map(|declaration| {
            Ok(RunnerV2DependencySourceMemberV1 {
                path: stage_a_path(
                    "runner_v2.base_values.dependency_source.path",
                    declaration.path,
                )?,
                content_root: RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                    "org.frankensim.fs-evidence-runner.runner-v2.stage-a.dependency-source-member.v1",
                    declaration.bytes,
                )),
            })
        })
        .collect()
}

fn independently_expected_dependency_source_root_v1(
    index: usize,
) -> Option<RunnerV2StageASourceMemberRootV1> {
    let bytes: &[u8] = match index {
        0 => include_bytes!("../../../../fs-blake3/src/lib.rs"),
        1 => include_bytes!("../../lib.rs"),
        2 => include_bytes!("../../canonical.rs"),
        3 => include_bytes!("../../catalog.rs"),
        4 => include_bytes!("../../construction.rs"),
        5 => include_bytes!("../../coverage.rs"),
        6 => include_bytes!("../../identity.rs"),
        7 => include_bytes!("../../limits.rs"),
        8 => include_bytes!("../../path.rs"),
        9 => include_bytes!("../../projection.rs"),
        10 => include_bytes!("../../schema_impact.rs"),
        11 => include_bytes!("../../value.rs"),
        12 => include_bytes!("../../runner_v2.rs"),
        13 => include_bytes!("../handoff.rs"),
        14 => include_bytes!("../work_packages.rs"),
        15 => include_bytes!("base_values.rs"),
        _ => return None,
    };
    Some(RunnerV2StageASourceMemberRootV1::from_content_hash(
        hash_domain(
            "org.frankensim.fs-evidence-runner.runner-v2.stage-a.dependency-source-member.v1",
            bytes,
        ),
    ))
}

fn validate_dependency_source_closure_diagnostic_v1(
    rows: &[RunnerV2DependencySourceMemberV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let paths = rows
        .iter()
        .map(|row| row.path.as_str().to_owned())
        .collect::<Vec<_>>();
    let roots = rows.iter().map(|row| row.content_root).collect::<Vec<_>>();
    validate_source_members_diagnostic_v1(
        "runner_v2.base_values.dependency_sources",
        "restore-exact-dependency-source-closure",
        &DEPENDENCY_SOURCE_PATHS_V1,
        &paths,
        &roots,
        independently_expected_dependency_source_root_v1,
    )
}

fn validate_dependency_source_closure_exact_v1(
    rows: &[RunnerV2DependencySourceMemberV1],
) -> Result<(), ConstructionErrorV2> {
    validate_dependency_source_closure_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

const STAGE_A_INVENTORY_MISSING_SENTINEL_V1: &str = "<missing>";
const STAGE_A_INVENTORY_END_SENTINEL_V1: &str = "<end-of-inventory>";
const STAGE_A_INVENTORY_REDACTED_SENTINEL_V1: &str = "<redacted:unregistered-inventory-member>";
const STAGE_A_INVENTORY_SEMANTIC_OWNER_V1: &str = "frankensim-epic-foundations-huq.24.1.1.1";

#[allow(
    clippy::too_many_arguments,
    reason = "the constructor makes every actionable diagnostic component explicit"
)]
fn stage_a_inventory_mismatch_v1(
    kind: ConstructionErrorKindV2,
    inventory: &'static str,
    first_mismatch_index0: usize,
    expected_identity: String,
    observed_safe_identity: String,
    observed_identity_redacted: bool,
    component: &'static str,
    expected_semantic_value: String,
    observed_safe_value: String,
    observed_value_redacted: bool,
    expected_count: usize,
    observed_count: usize,
    repair_target: &'static str,
) -> RunnerV2StageAInventoryMismatchV1 {
    RunnerV2StageAInventoryMismatchV1 {
        kind,
        inventory,
        first_mismatch_index0,
        expected_ordinal1: first_mismatch_index0 + 1,
        expected_identity,
        observed_safe_identity,
        observed_identity_redacted,
        component,
        expected_semantic_value,
        observed_safe_value,
        observed_value_redacted,
        semantic_owner: STAGE_A_INVENTORY_SEMANTIC_OWNER_V1,
        expected_count,
        observed_count,
        repairs: [RunnerV2StageAInventoryRepairV1 {
            rank: 1,
            kind: RepairActionKindV2::ChangeArguments,
            target: repair_target,
        }],
    }
}

fn validate_stage_a_exact_inventory_v1(
    inventory: &'static str,
    expected: &[String],
    observed: &[String],
    additionally_safe_observed: &[String],
    repair_target: &'static str,
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let first_mismatch_index0 = expected
        .iter()
        .zip(observed)
        .position(|(expected_row, observed_row)| expected_row != observed_row)
        .or_else(|| (expected.len() != observed.len()).then(|| expected.len().min(observed.len())));
    let Some(first_mismatch_index0) = first_mismatch_index0 else {
        return Ok(());
    };

    let expected_row = expected.get(first_mismatch_index0);
    let observed_row = observed.get(first_mismatch_index0);
    let kind = match expected.len().cmp(&observed.len()) {
        Ordering::Greater => ConstructionErrorKindV2::Missing,
        Ordering::Less => observed_row.map_or(ConstructionErrorKindV2::Unexpected, |row| {
            let expected_occurrences = expected
                .iter()
                .filter(|expected_row| *expected_row == row)
                .count();
            let observed_occurrences = observed
                .iter()
                .filter(|observed_row| *observed_row == row)
                .count();
            if observed_occurrences > expected_occurrences {
                ConstructionErrorKindV2::Duplicate
            } else {
                ConstructionErrorKindV2::Unexpected
            }
        }),
        Ordering::Equal => observed_row.map_or(ConstructionErrorKindV2::Missing, |row| {
            let expected_occurrences = expected
                .iter()
                .filter(|expected_row| *expected_row == row)
                .count();
            let observed_occurrences = observed
                .iter()
                .filter(|observed_row| *observed_row == row)
                .count();
            if observed_occurrences > expected_occurrences {
                ConstructionErrorKindV2::Duplicate
            } else if expected.contains(row) {
                ConstructionErrorKindV2::OutOfOrder
            } else {
                ConstructionErrorKindV2::Incompatible
            }
        }),
    };

    let expected_identity = expected_row
        .cloned()
        .unwrap_or_else(|| STAGE_A_INVENTORY_END_SENTINEL_V1.to_owned());
    let (observed_safe_identity, observed_identity_redacted) = match observed_row {
        None => (STAGE_A_INVENTORY_MISSING_SENTINEL_V1.to_owned(), false),
        Some(row)
            if expected.contains(row)
                || additionally_safe_observed.iter().any(|safe| safe == row) =>
        {
            (row.clone(), false)
        }
        Some(_) => (STAGE_A_INVENTORY_REDACTED_SENTINEL_V1.to_owned(), true),
    };

    Err(stage_a_inventory_mismatch_v1(
        kind,
        inventory,
        first_mismatch_index0,
        expected_identity.clone(),
        observed_safe_identity.clone(),
        observed_identity_redacted,
        "identity",
        expected_identity,
        observed_safe_identity,
        observed_identity_redacted,
        expected.len(),
        observed.len(),
        repair_target,
    ))
}

fn limit_boundary_inventory_ids_v1() -> Vec<String> {
    RunnerV2LimitBoundaryKindV1::ALL
        .iter()
        .map(|kind| kind.stable_name().to_owned())
        .collect()
}

fn validate_stage_a_limit_boundary_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.limit_boundaries",
        &limit_boundary_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-limit-boundary-catalog",
    )
}

fn validate_limit_boundary_definitions_diagnostic_v1(
    rows: &[RunnerV2LimitBoundaryDefinitionV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| row.stable_name.to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_limit_boundary_inventory_v1(&observed_ids)?;

    for (index, (row, expected_kind)) in rows
        .iter()
        .zip(RunnerV2LimitBoundaryKindV1::ALL)
        .enumerate()
    {
        let expected_identity = expected_kind.stable_name().to_owned();
        let expected_ordinal = u16::try_from(index + 1).map_err(|_| {
            stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "runner_v2.base_values.limit_boundaries",
                index,
                expected_identity.clone(),
                row.stable_name.to_owned(),
                false,
                "ordinal",
                "one-based ordinal representable as u16".to_owned(),
                (index + 1).to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-limit-boundary-catalog",
            )
        })?;
        if row.ordinal != expected_ordinal {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.limit_boundaries",
                index,
                expected_identity.clone(),
                row.stable_name.to_owned(),
                false,
                "ordinal",
                expected_ordinal.to_string(),
                row.ordinal.to_string(),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-limit-boundary-catalog",
            ));
        }
        if row.kind != expected_kind {
            return Err(stage_a_inventory_mismatch_v1(
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.limit_boundaries",
                index,
                expected_identity,
                row.stable_name.to_owned(),
                false,
                "kind",
                format!(
                    "{:02}:{}",
                    expected_kind.code(),
                    expected_kind.stable_name()
                ),
                format!("{:02}:{}", row.kind.code(), row.kind.stable_name()),
                false,
                rows.len(),
                rows.len(),
                "restore-exact-limit-boundary-catalog",
            ));
        }
    }
    Ok(())
}

fn validate_limit_boundary_definitions_exact_v1(
    rows: &[RunnerV2LimitBoundaryDefinitionV1],
) -> Result<(), ConstructionErrorV2> {
    validate_limit_boundary_definitions_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn limit_field_inventory_ids_v1() -> Vec<String> {
    STAGE_A_INDEPENDENT_LIMIT_NAMES_V1
        .iter()
        .map(|name| (*name).to_owned())
        .collect()
}

fn validate_stage_a_limit_field_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.limit_fields",
        &limit_field_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-limit-field-and-width-catalog",
    )
}

fn limit_literal_semantic_mismatch_v1(
    rows: &[RunnerV2LimitLiteralV1],
    index: usize,
    component: &'static str,
    expected: String,
    observed: String,
) -> RunnerV2StageAInventoryMismatchV1 {
    let identity = STAGE_A_INDEPENDENT_LIMIT_NAMES_V1[index].to_owned();
    stage_a_inventory_mismatch_v1(
        ConstructionErrorKindV2::Incompatible,
        "runner_v2.base_values.limit_fields",
        index,
        identity.clone(),
        identity,
        false,
        component,
        expected,
        observed,
        false,
        rows.len(),
        rows.len(),
        "restore-exact-limit-field-and-width-catalog",
    )
}

fn validate_limit_literal_definitions_diagnostic_v1(
    rows: &[RunnerV2LimitLiteralV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| independent_limit_name_v1(row.field).to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_limit_field_inventory_v1(&observed_ids)?;

    for (index, (row, expected)) in rows.iter().zip(STAGE_A_LIMIT_LITERALS_V1).enumerate() {
        if row.field != expected.field {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "field",
                format!(
                    "{}:{}",
                    expected.field.ordinal(),
                    independent_limit_name_v1(expected.field)
                ),
                format!(
                    "{}:{}",
                    row.field.ordinal(),
                    independent_limit_name_v1(row.field)
                ),
            ));
        }
        if row.width != expected.width {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "width",
                limit_width_name_v1(expected.width).to_owned(),
                limit_width_name_v1(row.width).to_owned(),
            ));
        }
        if row.unit != expected.unit {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "unit",
                independent_limit_unit_name_v1(expected.unit).to_owned(),
                independent_limit_unit_name_v1(row.unit).to_owned(),
            ));
        }
        if row.tightenability != expected.tightenability {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "tightenability",
                format!("{:?}", expected.tightenability),
                format!("{:?}", row.tightenability),
            ));
        }
        if row.minimum_rule != expected.minimum_rule {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "minimum-rule",
                format!("{:?}", expected.minimum_rule),
                format!("{:?}", row.minimum_rule),
            ));
        }
        if row.smoke != expected.smoke {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "smoke-ceiling",
                limit_value_safe_name_v1(expected.smoke),
                limit_value_safe_name_v1(row.smoke),
            ));
        }
        if row.full != expected.full {
            return Err(limit_literal_semantic_mismatch_v1(
                rows,
                index,
                "full-ceiling",
                limit_value_safe_name_v1(expected.full),
                limit_value_safe_name_v1(row.full),
            ));
        }
    }
    Ok(())
}

fn validate_limit_literal_definitions_exact_v1(
    rows: &[RunnerV2LimitLiteralV1],
) -> Result<(), ConstructionErrorV2> {
    validate_limit_literal_definitions_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn meta_inventory_ids_v1() -> Vec<String> {
    META_CELL_DEFINITIONS_V1
        .iter()
        .map(|row| {
            format!(
                "runner-v2.base-values.meta-{:03}-{}.v1",
                row.ordinal, row.id_suffix
            )
        })
        .collect()
}

fn validate_stage_a_meta_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.meta_cells",
        &meta_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-meta-cell-catalog",
    )
}

fn meta_definition_semantic_mismatch_v1(
    rows: &[MetaCellDefinitionV1],
    index: usize,
    component: &'static str,
    expected: String,
    observed: String,
) -> RunnerV2StageAInventoryMismatchV1 {
    let identity = meta_inventory_ids_v1()[index].clone();
    stage_a_inventory_mismatch_v1(
        ConstructionErrorKindV2::Incompatible,
        "runner_v2.base_values.meta_cells",
        index,
        identity.clone(),
        identity,
        false,
        component,
        expected,
        observed,
        false,
        rows.len(),
        rows.len(),
        "restore-exact-meta-cell-catalog",
    )
}

fn validate_meta_definitions_diagnostic_v1(
    rows: &[MetaCellDefinitionV1],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let observed_ids = rows
        .iter()
        .map(|row| {
            format!(
                "runner-v2.base-values.meta-{:03}-{}.v1",
                row.ordinal, row.id_suffix
            )
        })
        .collect::<Vec<_>>();
    validate_stage_a_meta_inventory_v1(&observed_ids)?;

    for (index, (row, expected)) in rows.iter().zip(META_CELL_DEFINITIONS_V1).enumerate() {
        if row.ordinal != expected.ordinal {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "ordinal",
                expected.ordinal.to_string(),
                row.ordinal.to_string(),
            ));
        }
        if row.group != expected.group {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "group",
                expected.group.code().to_string(),
                row.group.code().to_string(),
            ));
        }
        if row.operation != expected.operation {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "operation",
                format!("{:?}", expected.operation),
                format!("{:?}", row.operation),
            ));
        }
        if row.expected_outcome != expected.expected_outcome {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "expected-outcome",
                format!("{:?}", expected.expected_outcome),
                format!("{:?}", row.expected_outcome),
            ));
        }
        if row.expected_reason != expected.expected_reason {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "expected-reason",
                format!(
                    "{:02}:{}",
                    expected.expected_reason.code(),
                    runner_v2_raw_reason_name_v1(expected.expected_reason)
                ),
                format!(
                    "{:02}:{}",
                    row.expected_reason.code(),
                    runner_v2_raw_reason_name_v1(row.expected_reason)
                ),
            ));
        }
        if row.expected_partition != expected.expected_partition {
            return Err(meta_definition_semantic_mismatch_v1(
                rows,
                index,
                "expected-partition",
                expected.expected_partition.code().to_string(),
                row.expected_partition.code().to_string(),
            ));
        }
    }
    Ok(())
}

fn validate_meta_definitions_exact_v1(
    rows: &[MetaCellDefinitionV1],
) -> Result<(), ConstructionErrorV2> {
    validate_meta_definitions_diagnostic_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn limit_cell_inventory_ids_v1() -> Vec<String> {
    let mut rows = Vec::with_capacity(RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1);
    for literal in STAGE_A_LIMIT_LITERALS_V1 {
        for boundary in RunnerV2LimitBoundaryKindV1::ALL {
            rows.push(format!(
                "runner-v2.base-values.limit-{:03}.boundary-{:02}-{}.v1",
                literal.field.ordinal(),
                boundary.code(),
                boundary.stable_name(),
            ));
        }
    }
    rows
}

fn validate_stage_a_limit_cell_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.limit_cells",
        &limit_cell_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-71-by-12-limit-cell-inventory",
    )
}

fn complete_cell_inventory_ids_v1() -> Vec<String> {
    let mut rows = limit_cell_inventory_ids_v1();
    rows.extend(meta_inventory_ids_v1());
    rows
}

fn validate_stage_a_complete_cell_inventory_v1(
    rows: &[String],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.complete_cells",
        &complete_cell_inventory_ids_v1(),
        rows,
        &[],
        "restore-exact-852-plus-15-cell-inventory",
    )
}

fn validate_stage_a_canonical_schema_inventory_v1(
    rows: &[&str],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let expected = STAGE_A_CANONICAL_SCHEMA_NAMES_V1
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    let observed = rows
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    let additionally_safe_observed = vec![STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1.to_owned()];
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.schema_inventory.canonical",
        &expected,
        &observed,
        &additionally_safe_observed,
        "restore-exact-canonical-schema-inventory",
    )
}

fn validate_stage_a_canonical_schema_names_exact_v1(
    rows: &[&str],
) -> Result<(), ConstructionErrorV2> {
    validate_stage_a_canonical_schema_inventory_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn validate_stage_a_rootless_schema_inventory_v1(
    rows: &[&str],
) -> Result<(), RunnerV2StageAInventoryMismatchV1> {
    let expected = vec![STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1.to_owned()];
    let observed = rows
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    let additionally_safe_observed = STAGE_A_CANONICAL_SCHEMA_NAMES_V1
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    validate_stage_a_exact_inventory_v1(
        "runner_v2.base_values.schema_inventory.rootless_handoff",
        &expected,
        &observed,
        &additionally_safe_observed,
        "restore-sole-rootless-handoff-schema-inventory",
    )
}

fn validate_stage_a_rootless_schema_names_exact_v1(
    rows: &[&str],
) -> Result<(), ConstructionErrorV2> {
    validate_stage_a_rootless_schema_inventory_v1(rows)
        .map_err(|mismatch| mismatch.as_construction_error())
}

fn validate_stage_a_schema_partition_exact_v1(
    canonical_rows: &[&str],
    rootless_rows: &[&str],
) -> Result<(), ConstructionErrorV2> {
    validate_stage_a_canonical_schema_names_exact_v1(canonical_rows)?;
    validate_stage_a_rootless_schema_names_exact_v1(rootless_rows)?;
    let owned_count = canonical_rows
        .len()
        .checked_add(rootless_rows.len())
        .ok_or_else(|| {
            stage_a_error(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "runner_v2.base_values.schema_inventory.owned_count",
                "a checked 42 canonical plus one rootless inventory",
                usize::MAX,
            )
        })?;
    if owned_count != RUNNER_V2_BASE_VALUES_OWNED_SCHEMA_COUNT_V1 {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.schema_inventory.owned_count",
            "exactly 43 separately classified owned schemas",
            owned_count,
        ));
    }
    Ok(())
}

fn build_schema_impact_deferral_v1() -> Result<RunnerV2SchemaImpactDeferralV1, ConstructionErrorV2>
{
    validate_stage_a_schema_partition_exact_v1(
        &STAGE_A_CANONICAL_SCHEMA_NAMES_V1,
        &[STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1],
    )?;
    let canonical_schema_names = STAGE_A_CANONICAL_SCHEMA_NAMES_V1
        .iter()
        .copied()
        .map(|name| {
            RunnerV2CanonicalSchemaNameV1::from_source_literal(
                "runner_v2.base_values.schema_deferral.schema",
                name,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    if canonical_schema_names.len() != RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1 {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.schema_deferral",
            "exactly 42 canonical schemas deferred once to work package .3",
            canonical_schema_names.len(),
        ));
    }
    Ok(RunnerV2SchemaImpactDeferralV1 {
        resolution_owner: stage_a_token(
            "runner_v2.base_values.schema_deferral.resolution_owner",
            "frankensim-epic-foundations-huq.24.1.1.1.3",
        )?,
        canonical_schema_names: canonical_schema_names.into_boxed_slice(),
        future_manifest_root: TypedOptionV1::Absent,
        no_claim: stage_a_token(
            "runner_v2.base_values.schema_deferral.no_claim",
            RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
        )?,
    })
}

fn build_ac58_v1() -> Result<RunnerV2RootlessAc58FragmentV1, ConstructionErrorV2> {
    validate_stage_a_rootless_schema_names_exact_v1(&[STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1])?;
    Ok(RunnerV2RootlessAc58FragmentV1 {
        semantic_type: RunnerV2RootlessHandoffSchemaNameV1::from_source_literal(
            "runner_v2.base_values.ac58.semantic_type",
            STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1,
        )?,
        disposition: CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame,
        migration_policy: CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor,
        authority_surfaces_present_empty: true,
        no_claim: stage_a_token(
            "runner_v2.base_values.ac58.no_claim",
            RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
        )?,
    })
}

fn build_limit_fixture_v1() -> Result<RunnerV2LimitFixtureDeclarationV1, ConstructionErrorV2> {
    let family_rows_by_case = vec![0_u32].into_boxed_slice();
    let lifecycle_document_structural_minimum =
        checked_lifecycle_record_requirement(&family_rows_by_case).map_err(|violation| {
            stage_a_error(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "runner_v2.base_values.limit_fixture.lifecycle_minimum",
                "the checked one-case lifecycle requirement",
                violation.observed().as_u128(),
            )
        })?;
    if family_rows_by_case.len() != RUNNER_V2_LIMIT_FIXTURE_CASE_COUNT_V1
        || lifecycle_document_structural_minimum != RUNNER_V2_LIMIT_FIXTURE_LIFECYCLE_MINIMUM_V1
    {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.limit_fixture",
            "one executable zero-row case with lifecycle minimum five",
            lifecycle_document_structural_minimum,
        ));
    }
    Ok(RunnerV2LimitFixtureDeclarationV1 {
        executable: true,
        family_rows_by_case,
        declared_minimums_present_empty: true,
        lifecycle_document_structural_minimum,
        no_claim: stage_a_token(
            "runner_v2.base_values.limit_fixture.no_claim",
            RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
        )?,
    })
}

fn build_route_v1() -> Result<RunnerV2LocalRouteDeclarationV1, ConstructionErrorV2> {
    Ok(RunnerV2LocalRouteDeclarationV1 {
        route_id: stage_a_token(
            "runner_v2.base_values.route.id",
            RUNNER_V2_BASE_VALUES_LOCAL_ROUTE_ID_V1,
        )?,
        class: RunnerV2LocalRouteClassV1::LocalOnly,
        public_entry_point: RUNNER_V2_BASE_VALUES_PUBLIC_ENTRY_POINT_V1,
        execution_owner: stage_a_token(
            "runner_v2.base_values.route.execution_owner",
            "frankensim-epic-foundations-huq.24.1.1.1.7",
        )?,
        capability_profile: BaseCoverageCloseCapabilityProfileV1::None,
        external_driver: TypedOptionV1::Absent,
        no_claim: stage_a_token(
            "runner_v2.base_values.route.no_claim",
            RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
        )?,
    })
}

fn build_inapplicability_v1(
    stable_id: &'static str,
    reason: RunnerV2RawReasonV1,
    prerequisite: &'static str,
) -> Result<RunnerV2StageAInapplicabilityDeclarationV1, ConstructionErrorV2> {
    Ok(RunnerV2StageAInapplicabilityDeclarationV1 {
        stable_id: stage_a_token("runner_v2.base_values.inapplicable.id", stable_id)?,
        reason,
        owner: stage_a_token(
            "runner_v2.base_values.inapplicable.owner",
            "fs-evidence-runner.runner-v2.base-values",
        )?,
        prerequisite: stage_a_token(
            "runner_v2.base_values.inapplicable.prerequisite",
            prerequisite,
        )?,
        no_claim: stage_a_token(
            "runner_v2.base_values.inapplicable.no_claim",
            RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
        )?,
    })
}

fn build_capability_none_v1() -> Result<
    (
        BaseCoverageCloseCapabilityRegistryV1,
        BaseCoverageCloseCapabilityProfileRegistryV1,
        BaseCoverageCloseCapabilityContractV1,
    ),
    ConstructionErrorV2,
> {
    let registry = BaseCoverageCloseCapabilityRegistryV1::frozen()?;
    let profile_registry = BaseCoverageCloseCapabilityProfileRegistryV1::frozen(&registry)?;
    let contract = BaseCoverageCloseCapabilityContractV1::for_profile(
        &registry,
        &profile_registry,
        BaseCoverageCloseCapabilityProfileV1::None,
    )?;
    if !contract.required().is_empty() || !contract.permitted().is_empty() {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.capability_none",
            "exact-empty required and permitted declaration sets",
            contract.required().len(),
        ));
    }
    Ok((registry, profile_registry, contract))
}

fn build_local_budget_set_v1() -> Result<BaseCoverageCloseBudgetSetV1, ConstructionErrorV2> {
    let row = |axis, hard, soft, unit| {
        Ok(BaseCoverageCloseTypedBudgetV1::new(
            axis,
            hard,
            soft,
            BaseCoverageCloseLogicalUnitReferenceV1::fixed(unit)?,
        )?)
    };
    BaseCoverageCloseBudgetSetV1::new(
        BaseCoverageCloseBudgetProfileV1::LocalSourceValidation,
        vec![
            row(
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                LogicalUnitV2::Nanoseconds,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::Memory,
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                LogicalUnitV2::LogicalBytes,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::LogicalWork,
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                LogicalUnitV2::Operations,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::Processes,
                BaseCoverageCloseBudgetValueV1::U32(1),
                BaseCoverageCloseBudgetValueV1::U32(0),
                LogicalUnitV2::Count,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::Artifacts,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                LogicalUnitV2::EncodedBytes,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::Output,
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                LogicalUnitV2::EncodedBytes,
            )?,
            row(
                BaseCoverageCloseBudgetAxisV1::Logs,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                LogicalUnitV2::EncodedBytes,
            )?,
        ],
    )
}

fn cargo_manifest_declares_feature_v1(manifest: &str) -> bool {
    let mut dependency_scope = false;
    for raw_line in manifest.lines() {
        let syntax = unquoted_toml_syntax_v1(raw_line);
        let line = syntax.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(table) = line
            .strip_prefix('[')
            .and_then(|line| line.strip_suffix(']'))
        {
            let table = table.trim();
            if table == "features" {
                return true;
            }
            dependency_scope = cargo_dependency_table_v1(table);
            continue;
        }
        if dependency_scope {
            let compact = line
                .chars()
                .filter(|character| !character.is_ascii_whitespace())
                .collect::<String>();
            if compact
                .split([',', '{', '}'])
                .any(|clause| clause == "optional=true")
            {
                return true;
            }
        }
    }
    false
}

fn cargo_dependency_table_v1(table: &str) -> bool {
    table == "dependencies"
        || table == "build-dependencies"
        || table.starts_with("dependencies.")
        || table.starts_with("build-dependencies.")
        || table.ends_with(".dependencies")
        || table.ends_with(".build-dependencies")
        || table.contains(".dependencies.")
        || table.contains(".build-dependencies.")
}

fn unquoted_toml_syntax_v1(line: &str) -> String {
    let mut syntax = String::with_capacity(line.len());
    let mut single_quoted = false;
    let mut double_quoted = false;
    let mut escaped = false;
    for character in line.chars() {
        if double_quoted {
            if escaped {
                escaped = false;
            } else if character == '\\' {
                escaped = true;
            } else if character == '"' {
                double_quoted = false;
            }
            syntax.push(' ');
            continue;
        }
        if single_quoted {
            if character == '\'' {
                single_quoted = false;
            }
            syntax.push(' ');
            continue;
        }
        match character {
            '#' => break,
            '"' => {
                double_quoted = true;
                syntax.push(' ');
            }
            '\'' => {
                single_quoted = true;
                syntax.push(' ');
            }
            _ => syntax.push(character),
        }
    }
    syntax
}

fn build_version_requirements_v1(
    dependency_sources: &[RunnerV2DependencySourceMemberV1],
) -> Result<RunnerV2StageAVersionRequirementsV1, ConstructionErrorV2> {
    validate_dependency_source_closure_exact_v1(dependency_sources)?;
    let source_frame =
        CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-SOURCE\x01", 16 * 1024, |sink| {
            for member in dependency_sources {
                sink.push_str(
                    "runner_v2.base_values.version.source_path",
                    member.path.as_str(),
                )?;
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.version.source_member",
                    member.content_root.bytes(),
                )?;
            }
            Ok(())
        })?;
    let source_identity = source_nominal_root_v1(
        source_frame.root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.source-closure.v1"),
    )?;

    let build_frame =
        CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-BUILD\x01", 1024 * 1024, |sink| {
            sink.push_bytes(
                "runner_v2.base_values.version.crate_manifest",
                include_bytes!("../../../Cargo.toml"),
            )?;
            sink.push_bytes(
                "runner_v2.base_values.version.hash_crate_manifest",
                include_bytes!("../../../../fs-blake3/Cargo.toml"),
            )?;
            sink.push_bytes(
                "runner_v2.base_values.version.workspace_manifest",
                include_bytes!("../../../../../Cargo.toml"),
            )?;
            sink.push_bytes(
                "runner_v2.base_values.version.workspace_lock",
                include_bytes!("../../../../../Cargo.lock"),
            )?;
            sink.push_bytes(
                "runner_v2.base_values.version.constellation_lock",
                include_bytes!("../../../../../constellation.lock"),
            )?;
            sink.push_bytes(
                "runner_v2.base_values.version.contract",
                include_bytes!("../../../CONTRACT.md"),
            )
        })?;
    let build_identity = build_nominal_root_v1(
        build_frame.root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.build-inputs.v1"),
    )?;

    let toolchain_identity = toolchain_nominal_root_v1(hash_domain(
        "org.frankensim.fs-evidence-runner.runner-v2.stage-a.toolchain.v1",
        include_bytes!("../../../../../rust-toolchain.toml"),
    ))?;

    let schema_frame =
        CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-SCHEMAS\x01", 16 * 1024, |sink| {
            sink.push_u32(
                "runner_v2.base_values.version.canonical_schema_count",
                checked_u32_v1(
                    "runner_v2.base_values.version.canonical_schema_count",
                    STAGE_A_CANONICAL_SCHEMA_NAMES_V1.len(),
                )?,
            )?;
            for schema in STAGE_A_CANONICAL_SCHEMA_NAMES_V1 {
                sink.push_str("runner_v2.base_values.version.canonical_schema", schema)?;
            }
            sink.push_u32(
                "runner_v2.base_values.version.rootless_schema_count",
                checked_u32_v1(
                    "runner_v2.base_values.version.rootless_schema_count",
                    RUNNER_V2_BASE_VALUES_ROOTLESS_SCHEMA_COUNT_V1,
                )?,
            )?;
            sink.push_str(
                "runner_v2.base_values.version.rootless_schema",
                STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1,
            )?;
            Ok(())
        })?;
    let schema_inventory_root = RunnerV2StageASchemaInventoryRootV1::from_content_hash(
        schema_frame
            .root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.schema-inventory.v1"),
    );
    let crate_manifest = include_str!("../../../Cargo.toml");
    if cargo_manifest_declares_feature_v1(crate_manifest) {
        return Err(stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.version.feature_declaration",
            "no explicit feature table or implicit optional-dependency feature for the exact-empty Stage-A feature set",
            1_u16,
        ));
    }
    let feature_declaration_root =
        RunnerV2StageAFeatureDeclarationRootV1::from_content_hash(hash_domain(
            "org.frankensim.fs-evidence-runner.runner-v2.stage-a.features.v1",
            b"exact-empty-feature-set;validated-no-features-table-or-optional-dependencies-in-crate-manifest",
        ));

    Ok(RunnerV2StageAVersionRequirementsV1 {
        api_generation: RUNNER_SPEC_V2_API_GENERATION,
        wire_version: RUNNER_V2_WIRE_VERSION,
        predecessor_policy: RUNNER_V2_PREDECESSOR_POLICY,
        source_identity,
        build_identity,
        toolchain_identity,
        schema_inventory_root,
        feature_declaration_root,
        target: BaseCoverageCloseTargetV1::TargetIndependentPureValidation,
        profile: BaseCoverageCloseProfileV1::CrateTest,
    })
}

fn build_five_explicits_v1(
    dependency_sources: &[RunnerV2DependencySourceMemberV1],
) -> Result<RunnerV2StageAFiveExplicitsV1, ConstructionErrorV2> {
    let budgets = build_local_budget_set_v1()?;
    let versions = build_version_requirements_v1(dependency_sources)?;
    let (capability_registry, capability_profile_registry, capability_contract) =
        build_capability_none_v1()?;
    let no_claim = stage_a_token(
        "runner_v2.base_values.five.no_claim",
        RUNNER_V2_BASE_VALUES_NO_CLAIM_V1,
    )?;
    let seed = BaseCoverageCloseSeedExplicitV1::Inapplicable {
        reason: SeedInapplicableCodeV1::NoRandomnessByContract,
    };

    let frame = CanonicalFrameV1::preflighted(b"FSRUNNER-STAGE-A-FIVE\x01", 16 * 1024, |sink| {
        sink.push_presence("runner_v2.base_values.five.inputs_present", true)?;
        sink.push_u32("runner_v2.base_values.five.input_count", 0)?;
        sink.push_presence("runner_v2.base_values.five.grants_present", true)?;
        sink.push_u32("runner_v2.base_values.five.grant_count", 0)?;
        sink.push_presence(
            "runner_v2.base_values.five.expected_observations_present",
            true,
        )?;
        sink.push_u32("runner_v2.base_values.five.expected_observation_count", 0)?;
        sink.push_u16(
            "runner_v2.base_values.five.seed_inapplicable",
            SeedInapplicableCodeV1::NoRandomnessByContract.code(),
        )?;
        push_budget_set_v1(sink, &budgets)?;
        sink.push_u16(
            "runner_v2.base_values.five.api_generation",
            versions.api_generation.code(),
        )?;
        sink.push_u16(
            "runner_v2.base_values.five.wire_version",
            versions.wire_version.code(),
        )?;
        sink.push_str(
            "runner_v2.base_values.five.predecessor",
            versions.predecessor_policy.name(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.source",
            versions.source_identity.bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.build",
            versions.build_identity.bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.toolchain",
            versions.toolchain_identity.bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.schema_inventory",
            versions.schema_inventory_root.bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.features",
            versions.feature_declaration_root.bytes(),
        )?;
        sink.push_u16("runner_v2.base_values.five.target", versions.target.code())?;
        sink.push_u16(
            "runner_v2.base_values.five.profile",
            versions.profile.code(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.capability_registry",
            capability_registry.root().content_hash().as_bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.capability_profile_registry",
            capability_profile_registry.root().content_hash().as_bytes(),
        )?;
        sink.push_fixed_bytes_32(
            "runner_v2.base_values.five.capability_contract",
            capability_contract.root().content_hash().as_bytes(),
        )?;
        sink.push_u32(
            "runner_v2.base_values.five.required_capability_count",
            u32::try_from(capability_contract.required().len()).map_err(|_| {
                stage_a_error(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.base_values.five.required_capability_count",
                    "capability count representable as u32",
                    capability_contract.required().len(),
                )
            })?,
        )?;
        sink.push_u32(
            "runner_v2.base_values.five.permitted_capability_count",
            u32::try_from(capability_contract.permitted().len()).map_err(|_| {
                stage_a_error(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "runner_v2.base_values.five.permitted_capability_count",
                    "capability count representable as u32",
                    capability_contract.permitted().len(),
                )
            })?,
        )?;
        sink.push_str("runner_v2.base_values.five.no_claim", no_claim.as_str())
    })?;
    let root = RunnerV2StageAFiveExplicitsRootV1::from_content_hash(
        frame.root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.five-explicits.v1"),
    );

    Ok(RunnerV2StageAFiveExplicitsV1 {
        numeric_inputs_present_empty: true,
        numeric_grants_present_empty: true,
        expected_numeric_observations_present_empty: true,
        seed,
        budgets,
        versions,
        capability_registry,
        capability_profile_registry,
        capability_contract,
        no_claim,
        root,
    })
}

fn push_budget_set_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    budgets: &BaseCoverageCloseBudgetSetV1,
) -> Result<(), ConstructionErrorV2> {
    sink.push_u16(
        "runner_v2.base_values.five.budget_profile",
        budgets.profile().code(),
    )?;
    for budget in budgets.rows() {
        sink.push_u16(
            "runner_v2.base_values.five.budget_axis",
            budget.axis().code(),
        )?;
        push_budget_value_v1(sink, budget.hard())?;
        push_budget_value_v1(sink, budget.soft())?;
        sink.push_u16(
            "runner_v2.base_values.five.budget_unit",
            budget.unit().unit().tag(),
        )?;
    }
    Ok(())
}

fn push_budget_value_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    value: BaseCoverageCloseBudgetValueV1,
) -> Result<(), ConstructionErrorV2> {
    match value {
        BaseCoverageCloseBudgetValueV1::U32(value) => {
            sink.push_u16("runner_v2.base_values.five.budget_width", 1)?;
            sink.push_bytes(
                "runner_v2.base_values.five.budget_value",
                &value.to_be_bytes(),
            )
        }
        BaseCoverageCloseBudgetValueV1::U64(value) => {
            sink.push_u16("runner_v2.base_values.five.budget_width", 2)?;
            sink.push_bytes(
                "runner_v2.base_values.five.budget_value",
                &value.to_be_bytes(),
            )
        }
        BaseCoverageCloseBudgetValueV1::U128(value) => {
            sink.push_u16("runner_v2.base_values.five.budget_width", 3)?;
            sink.push_bytes(
                "runner_v2.base_values.five.budget_value",
                &value.to_be_bytes(),
            )
        }
    }
}

#[allow(
    clippy::too_many_arguments,
    reason = "the Stage-A declaration root exact-binds every independently owned component"
)]
fn stage_a_declaration_root_v1(
    package_id: &StableTokenV2,
    cells: &[RunnerV2StageACellDeclarationV1],
    oracles: &[RunnerV2StageAOracleRowV1],
    projections: &[RunnerV2StageAProjectionRowV1],
    limit_fixture: &RunnerV2LimitFixtureDeclarationV1,
    retained_domain_obligations: &[RunnerV2RetainedDomainObligationV1],
    limit_mutation_obligations: &[RunnerV2LimitMutationObligationV1],
    five_explicits: &RunnerV2StageAFiveExplicitsV1,
    route: &RunnerV2LocalRouteDeclarationV1,
    common_requirements: &[RunnerV2CommonContractRequirementV1],
    future_sources: &[RunnerV2FutureSourceRequirementV1],
    owner_source_fragment: &[RunnerV2OwnerSourceMemberV1],
    dependency_source_closure: &[RunnerV2DependencySourceMemberV1],
    schema_impact_deferral: &RunnerV2SchemaImpactDeferralV1,
    ac58: &RunnerV2RootlessAc58FragmentV1,
    shard: &RunnerV2StageAInapplicabilityDeclarationV1,
    resume: &RunnerV2StageAInapplicabilityDeclarationV1,
    no_claim: &StableTokenV2,
) -> Result<RunnerV2StageADeclarationRootV1, ConstructionErrorV2> {
    let frame = CanonicalFrameV1::preflighted(
        b"FSRUNNER-STAGE-A-BASE-VALUES\x01",
        STAGE_A_FRAME_MAX_BYTES_V1,
        |sink| {
            sink.push_str("runner_v2.base_values.package_id", package_id.as_str())?;
            sink.push_u32(
                "runner_v2.base_values.cell_count",
                checked_u32_v1("runner_v2.base_values.cell_count", cells.len())?,
            )?;
            for cell in cells {
                sink.push_u16("runner_v2.base_values.cell.ordinal", cell.ordinal)?;
                sink.push_str("runner_v2.base_values.cell.id", cell.cell_id.as_str())?;
                sink.push_u16("runner_v2.base_values.cell.group", cell.group.code())?;
                push_operation_v1(sink, cell.operation)?;
                sink.push_u32(
                    "runner_v2.base_values.cell.companion_count",
                    checked_u32_v1(
                        "runner_v2.base_values.cell.companion_count",
                        cell.companion_normalization.len(),
                    )?,
                )?;
                for companion in &cell.companion_normalization {
                    sink.push_u16(
                        "runner_v2.base_values.cell.companion_field",
                        companion.field.ordinal(),
                    )?;
                    push_limit_value_v1(sink, companion.value)?;
                }
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.cell.oracle_root",
                    cell.oracle_root.bytes(),
                )?;
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.cell.case_manifest_root",
                    cell.case_manifest_root.bytes(),
                )?;
            }

            sink.push_u32(
                "runner_v2.base_values.oracle_count",
                checked_u32_v1("runner_v2.base_values.oracle_count", oracles.len())?,
            )?;
            for oracle in oracles {
                sink.push_str(
                    "runner_v2.base_values.oracle.cell_id",
                    oracle.cell_id.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.oracle.outcome",
                    oracle.expected_outcome.code(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.oracle.reason",
                    oracle.expected_reason.code(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.oracle.partition",
                    oracle.expected_partition.code(),
                )?;
                sink.push_fixed_bytes_32("runner_v2.base_values.oracle.root", oracle.root.bytes())?;
            }

            sink.push_u32(
                "runner_v2.base_values.projection_count",
                checked_u32_v1("runner_v2.base_values.projection_count", projections.len())?,
            )?;
            for projection in projections {
                sink.push_u16(
                    "runner_v2.base_values.projection.ordinal",
                    projection.ordinal,
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.cell_id",
                    projection.cell_id.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.consumer_route",
                    projection.consumer_route.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.consumer_owner",
                    projection.consumer_owner.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.dispatcher",
                    projection.dispatcher.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.posix_script",
                    projection.posix_script.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.windows_script",
                    projection.windows_script.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.projection.partition",
                    projection.expected_partition.code(),
                )?;
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.projection.case_manifest",
                    projection.case_manifest_root.bytes(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.projection.no_claim",
                    projection.no_claim.as_str(),
                )?;
            }

            sink.push_presence(
                "runner_v2.base_values.limit_fixture.executable",
                limit_fixture.executable,
            )?;
            sink.push_u32(
                "runner_v2.base_values.limit_fixture.case_count",
                checked_u32_v1(
                    "runner_v2.base_values.limit_fixture.case_count",
                    limit_fixture.family_rows_by_case.len(),
                )?,
            )?;
            for rows in &limit_fixture.family_rows_by_case {
                sink.push_u32("runner_v2.base_values.limit_fixture.family_rows", *rows)?;
            }
            sink.push_presence(
                "runner_v2.base_values.limit_fixture.declared_minimums_present_empty",
                limit_fixture.declared_minimums_present_empty,
            )?;
            sink.push_u32(
                "runner_v2.base_values.limit_fixture.lifecycle_minimum",
                limit_fixture.lifecycle_document_structural_minimum,
            )?;
            sink.push_str(
                "runner_v2.base_values.limit_fixture.no_claim",
                limit_fixture.no_claim.as_str(),
            )?;

            sink.push_u32(
                "runner_v2.base_values.retained_domain_count",
                checked_u32_v1(
                    "runner_v2.base_values.retained_domain_count",
                    retained_domain_obligations.len(),
                )?,
            )?;
            for obligation in retained_domain_obligations {
                sink.push_u16(
                    "runner_v2.base_values.retained_domain.ordinal",
                    obligation.ordinal,
                )?;
                sink.push_str(
                    "runner_v2.base_values.retained_domain.stable_id",
                    obligation.stable_id.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.retained_domain.facet",
                    obligation.facet.code(),
                )?;
            }
            sink.push_u32(
                "runner_v2.base_values.limit_mutation_count",
                checked_u32_v1(
                    "runner_v2.base_values.limit_mutation_count",
                    limit_mutation_obligations.len(),
                )?,
            )?;
            for obligation in limit_mutation_obligations {
                sink.push_u16(
                    "runner_v2.base_values.limit_mutation.ordinal",
                    obligation.ordinal,
                )?;
                sink.push_str(
                    "runner_v2.base_values.limit_mutation.stable_id",
                    obligation.stable_id.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.limit_mutation.field",
                    obligation.field.ordinal(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.limit_mutation.field_name",
                    obligation.field_name.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.limit_mutation.declared_width",
                    match obligation.declared_width {
                        RunnerLimitWidthV2::U32 => 32,
                        RunnerLimitWidthV2::U64 => 64,
                    },
                )?;
                push_limit_value_v1(sink, obligation.opposite_width_zero)?;
                sink.push_str(
                    "runner_v2.base_values.limit_mutation.unit",
                    independent_limit_unit_name_v1(obligation.unit),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.limit_mutation.expected_reason",
                    obligation.expected_reason.code(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.limit_mutation.diagnostic_owner",
                    obligation.diagnostic_owner.as_str(),
                )?;
                sink.push_u8(
                    "runner_v2.base_values.limit_mutation.repair_rank",
                    obligation.repair_rank,
                )?;
                sink.push_u16(
                    "runner_v2.base_values.limit_mutation.repair_kind",
                    obligation.repair_kind.code(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.limit_mutation.repair_target",
                    obligation.repair_target.as_str(),
                )?;
            }

            sink.push_fixed_bytes_32(
                "runner_v2.base_values.five_root",
                five_explicits.root.bytes(),
            )?;
            sink.push_str("runner_v2.base_values.route_id", route.route_id.as_str())?;
            sink.push_u16("runner_v2.base_values.route_class", route.class.code())?;
            sink.push_str(
                "runner_v2.base_values.entry_point",
                route.public_entry_point,
            )?;
            sink.push_str(
                "runner_v2.base_values.route_execution_owner",
                route.execution_owner.as_str(),
            )?;
            sink.push_u32(
                "runner_v2.base_values.local_in_process_route_count",
                checked_u32_v1(
                    "runner_v2.base_values.local_in_process_route_count",
                    RUNNER_V2_BASE_VALUES_LOCAL_IN_PROCESS_ROUTE_COUNT_V1,
                )?,
            )?;
            sink.push_u32(
                "runner_v2.base_values.execution_owned_route_count",
                checked_u32_v1(
                    "runner_v2.base_values.execution_owned_route_count",
                    RUNNER_V2_BASE_VALUES_EXECUTION_OWNED_ROUTE_COUNT_V1,
                )?,
            )?;
            sink.push_u32(
                "runner_v2.base_values.contribution_only_route_count",
                checked_u32_v1(
                    "runner_v2.base_values.contribution_only_route_count",
                    RUNNER_V2_BASE_VALUES_CONTRIBUTION_ONLY_ROUTE_COUNT_V1,
                )?,
            )?;
            sink.push_u32(
                "runner_v2.base_values.inapplicable_route_count",
                checked_u32_v1(
                    "runner_v2.base_values.inapplicable_route_count",
                    RUNNER_V2_BASE_VALUES_INAPPLICABLE_ROUTE_COUNT_V1,
                )?,
            )?;
            sink.push_u16(
                "runner_v2.base_values.capability_profile",
                route.capability_profile.code(),
            )?;
            sink.push_presence(
                "runner_v2.base_values.external_driver",
                matches!(route.external_driver, TypedOptionV1::Present(_)),
            )?;

            sink.push_u32(
                "runner_v2.base_values.requirement_count",
                checked_u32_v1(
                    "runner_v2.base_values.requirement_count",
                    common_requirements.len(),
                )?,
            )?;
            for requirement in common_requirements {
                sink.push_u16(
                    "runner_v2.base_values.requirement.ordinal",
                    requirement.ordinal,
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.slot_id",
                    requirement.slot_id.as_str(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.requirement.api_generation",
                    requirement.api_generation.code(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.requirement.wire_version",
                    requirement.wire_version.code(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.predecessor",
                    requirement.predecessor_policy.name(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.semantic_owner",
                    requirement.semantic_owner.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.realization_owner",
                    requirement.realization_owner.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.role",
                    requirement.future_nominal_role.as_str(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.domain",
                    requirement.future_domain.as_str(),
                )?;
                sink.push_u8(
                    "runner_v2.base_values.requirement.planes",
                    requirement.included_planes.mask(),
                )?;
                sink.push_u16(
                    "runner_v2.base_values.requirement.stage",
                    requirement.fulfillment_stage.code(),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.resolution_owner",
                    requirement.resolution_owner.as_str(),
                )?;
                sink.push_presence(
                    "runner_v2.base_values.requirement.future_root",
                    matches!(requirement.future_root, TypedOptionV1::Present(_)),
                )?;
                sink.push_str(
                    "runner_v2.base_values.requirement.no_claim",
                    requirement.no_claim.as_str(),
                )?;
            }

            sink.push_u32(
                "runner_v2.base_values.future_source_count",
                checked_u32_v1(
                    "runner_v2.base_values.future_source_count",
                    future_sources.len(),
                )?,
            )?;
            sink.push_u32(
                "runner_v2.base_values.final_source_count",
                checked_u32_v1(
                    "runner_v2.base_values.final_source_count",
                    RUNNER_V2_FINAL_SOURCE_COUNT_V1,
                )?,
            )?;
            for source in future_sources {
                sink.push_u16(
                    "runner_v2.base_values.future_source.ordinal",
                    source.final_ordinal,
                )?;
                sink.push_str(
                    "runner_v2.base_values.future_source.path",
                    source.path.as_str(),
                )?;
                sink.push_presence(
                    "runner_v2.base_values.future_source.root",
                    matches!(source.future_content_root, TypedOptionV1::Present(_)),
                )?;
            }

            sink.push_u32(
                "runner_v2.base_values.owner_source_count",
                checked_u32_v1(
                    "runner_v2.base_values.owner_source_count",
                    owner_source_fragment.len(),
                )?,
            )?;
            for source in owner_source_fragment {
                sink.push_str(
                    "runner_v2.base_values.owner_source.path",
                    source.path.as_str(),
                )?;
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.owner_source.root",
                    source.content_root.bytes(),
                )?;
            }

            sink.push_u32(
                "runner_v2.base_values.dependency_source_count",
                checked_u32_v1(
                    "runner_v2.base_values.dependency_source_count",
                    dependency_source_closure.len(),
                )?,
            )?;
            for source in dependency_source_closure {
                sink.push_str(
                    "runner_v2.base_values.dependency_source.path",
                    source.path.as_str(),
                )?;
                sink.push_fixed_bytes_32(
                    "runner_v2.base_values.dependency_source.root",
                    source.content_root.bytes(),
                )?;
            }

            sink.push_str(
                "runner_v2.base_values.schema_deferral.resolution_owner",
                schema_impact_deferral.resolution_owner.as_str(),
            )?;
            sink.push_u32(
                "runner_v2.base_values.schema_deferral.schema_count",
                checked_u32_v1(
                    "runner_v2.base_values.schema_deferral.schema_count",
                    schema_impact_deferral.canonical_schema_names.len(),
                )?,
            )?;
            for schema in &schema_impact_deferral.canonical_schema_names {
                sink.push_str(
                    "runner_v2.base_values.schema_deferral.schema",
                    schema.as_str(),
                )?;
            }
            sink.push_presence(
                "runner_v2.base_values.schema_deferral.future_manifest_root",
                matches!(
                    schema_impact_deferral.future_manifest_root,
                    TypedOptionV1::Present(_)
                ),
            )?;
            sink.push_str(
                "runner_v2.base_values.schema_deferral.no_claim",
                schema_impact_deferral.no_claim.as_str(),
            )?;

            sink.push_str(
                "runner_v2.base_values.ac58.semantic_type",
                ac58.semantic_type.as_str(),
            )?;
            sink.push_u16(
                "runner_v2.base_values.ac58.disposition",
                ac58.disposition.code(),
            )?;
            sink.push_u16(
                "runner_v2.base_values.ac58.migration",
                ac58.migration_policy.code(),
            )?;
            sink.push_presence(
                "runner_v2.base_values.ac58.authority_surfaces_present_empty",
                ac58.authority_surfaces_present_empty,
            )?;

            push_inapplicability_v1(sink, shard)?;
            push_inapplicability_v1(sink, resume)?;
            sink.push_str("runner_v2.base_values.no_claim", no_claim.as_str())
        },
    )?;
    Ok(RunnerV2StageADeclarationRootV1::from_content_hash(
        frame
            .root("org.frankensim.fs-evidence-runner.runner-v2.stage-a.base-values-declaration.v1"),
    ))
}

fn push_operation_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    operation: RunnerV2StageACellOperationV1,
) -> Result<(), ConstructionErrorV2> {
    match operation {
        RunnerV2StageACellOperationV1::Limit {
            field,
            boundary,
            value,
        } => {
            sink.push_u16("runner_v2.base_values.operation.kind", 1)?;
            sink.push_u16(
                "runner_v2.base_values.operation.limit_field",
                field.ordinal(),
            )?;
            sink.push_u16("runner_v2.base_values.operation.boundary", boundary.code())?;
            sink.push_presence(
                "runner_v2.base_values.operation.value",
                matches!(value, TypedOptionV1::Present(_)),
            )?;
            if let TypedOptionV1::Present(value) = value {
                push_limit_value_v1(sink, value)?;
            }
            Ok(())
        }
        RunnerV2StageACellOperationV1::Meta(operation) => {
            sink.push_u16("runner_v2.base_values.operation.kind", 2)?;
            sink.push_u16(
                "runner_v2.base_values.operation.meta",
                meta_operation_code_v1(operation),
            )
        }
    }
}

fn push_limit_value_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    value: RunnerLimitValueV2,
) -> Result<(), ConstructionErrorV2> {
    match value {
        RunnerLimitValueV2::U32(value) => {
            sink.push_u16("runner_v2.base_values.limit_value.width", 1)?;
            sink.push_bytes(
                "runner_v2.base_values.limit_value.value",
                &value.to_be_bytes(),
            )
        }
        RunnerLimitValueV2::U64(value) => {
            sink.push_u16("runner_v2.base_values.limit_value.width", 2)?;
            sink.push_bytes(
                "runner_v2.base_values.limit_value.value",
                &value.to_be_bytes(),
            )
        }
    }
}

fn push_inapplicability_v1(
    sink: &mut dyn CanonicalFrameSinkV1,
    declaration: &RunnerV2StageAInapplicabilityDeclarationV1,
) -> Result<(), ConstructionErrorV2> {
    sink.push_str(
        "runner_v2.base_values.inapplicability.id",
        declaration.stable_id.as_str(),
    )?;
    sink.push_u16(
        "runner_v2.base_values.inapplicability.reason",
        declaration.reason.code(),
    )?;
    sink.push_str(
        "runner_v2.base_values.inapplicability.owner",
        declaration.owner.as_str(),
    )?;
    sink.push_str(
        "runner_v2.base_values.inapplicability.prerequisite",
        declaration.prerequisite.as_str(),
    )?;
    sink.push_str(
        "runner_v2.base_values.inapplicability.no_claim",
        declaration.no_claim.as_str(),
    )
}

const fn meta_operation_code_v1(operation: RunnerV2StageAMetaOperationV1) -> u16 {
    match operation {
        RunnerV2StageAMetaOperationV1::TypedAbsenceDistinctFromZero => 1,
        RunnerV2StageAMetaOperationV1::F32NamedTotalOrder => 2,
        RunnerV2StageAMetaOperationV1::F64NamedTotalOrder => 3,
        RunnerV2StageAMetaOperationV1::CapabilityNoneContract => 4,
        RunnerV2StageAMetaOperationV1::CommonRequirementsExact => 5,
        RunnerV2StageAMetaOperationV1::CommonRequirementReorderedRefusal => 6,
        RunnerV2StageAMetaOperationV1::FutureSourcesExact => 7,
        RunnerV2StageAMetaOperationV1::RootlessAc58 => 8,
        RunnerV2StageAMetaOperationV1::OwnerSourceFragment => 9,
        RunnerV2StageAMetaOperationV1::LocalRoute => 10,
        RunnerV2StageAMetaOperationV1::DiagnosticRedaction => 11,
        RunnerV2StageAMetaOperationV1::ReproductionDeclaration => 12,
        RunnerV2StageAMetaOperationV1::CompileFailOrderingSurface => 13,
        RunnerV2StageAMetaOperationV1::ShardInapplicable => 14,
        RunnerV2StageAMetaOperationV1::ResumeInapplicable => 15,
    }
}

fn source_nominal_root_v1(
    content: ContentHash,
) -> Result<SourceIdentityRootV2, ConstructionErrorV2> {
    SourceIdentityRootV2::from_digest(DigestValueV2::from_array(
        DigestRoleV2::Source,
        SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
        *content.as_bytes(),
    ))
    .map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.version.source_identity",
            "the exact SourceIdentityRootV2 role and domain",
            0_u16,
        )
    })
}

fn build_nominal_root_v1(content: ContentHash) -> Result<BuildIdentityRootV2, ConstructionErrorV2> {
    BuildIdentityRootV2::from_digest(DigestValueV2::from_array(
        DigestRoleV2::Build,
        BuildIdentityRootV2::DESCRIPTOR.domain_witness(),
        *content.as_bytes(),
    ))
    .map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.version.build_identity",
            "the exact BuildIdentityRootV2 role and domain",
            0_u16,
        )
    })
}

fn toolchain_nominal_root_v1(
    content: ContentHash,
) -> Result<ToolchainIdentityRootV2, ConstructionErrorV2> {
    ToolchainIdentityRootV2::from_digest(DigestValueV2::from_array(
        DigestRoleV2::Toolchain,
        ToolchainIdentityRootV2::DESCRIPTOR.domain_witness(),
        *content.as_bytes(),
    ))
    .map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.version.toolchain_identity",
            "the exact ToolchainIdentityRootV2 role and domain",
            0_u16,
        )
    })
}

fn limit_cell_id_v1(
    field: RunnerLimitFieldV2,
    boundary: RunnerV2LimitBoundaryKindV1,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    stage_a_token_owned(
        "runner_v2.base_values.limit_cell_id",
        format!(
            "runner-v2.base-values.limit-{:03}.boundary-{:02}-{}.v1",
            field.ordinal(),
            boundary.code(),
            boundary.stable_name(),
        ),
    )
}

fn meta_cell_id_v1(
    index: usize,
    suffix: &'static str,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    stage_a_token_owned(
        "runner_v2.base_values.meta_cell_id",
        format!("runner-v2.base-values.meta-{:03}-{}.v1", index + 1, suffix),
    )
}

fn validate_unique_cell_ids_v1<'a>(
    values: impl IntoIterator<Item = &'a StableTokenV2>,
) -> Result<(), ConstructionErrorV2> {
    let mut seen = std::collections::BTreeSet::new();
    for (index, value) in values.into_iter().enumerate() {
        if !seen.insert(value.as_str()) {
            return Err(stage_a_error(
                ConstructionErrorKindV2::Duplicate,
                "runner_v2.base_values.cell_id",
                "one occurrence of every stable cell identity",
                index,
            ));
        }
    }
    Ok(())
}

fn checked_u32_v1(field: &'static str, value: usize) -> Result<u32, ConstructionErrorV2> {
    u32::try_from(value).map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::ArithmeticOverflow,
            field,
            "a source-bounded count representable as u32",
            value,
        )
    })
}

fn stage_a_token(
    field: &'static str,
    value: &'static str,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    stage_a_token_owned(field, value.to_owned())
}

fn stage_a_token_owned(
    field: &'static str,
    value: String,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    let length = value.len();
    StableTokenV2::new(value).map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            field,
            "one exact bounded stable token",
            length,
        )
    })
}

fn stage_a_path(
    field: &'static str,
    value: &str,
) -> Result<LogicalBundlePathV1, ConstructionErrorV2> {
    LogicalBundlePathV1::new(value).map_err(|_| {
        stage_a_error(
            ConstructionErrorKindV2::Incompatible,
            field,
            "one exact validated logical relative path",
            value.len(),
        )
    })
}

fn stage_a_error(
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
    observed: impl Into<crate::construction::ConstructionObservedV2>,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(kind, field, expected, observed)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::runner_v2::handoff::{RunnerV2SafeNumericUnitV1, RunnerV2SafeNumericValueV1};

    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestBoundaryExpectedV1 {
        ordinal: u16,
        boundary: RunnerV2LimitBoundaryKindV1,
        code: u16,
        stable_name: &'static str,
    }

    macro_rules! test_boundary_v1 {
        ($ordinal:literal, $boundary:ident, $code:literal, $name:literal) => {
            TestBoundaryExpectedV1 {
                ordinal: $ordinal,
                boundary: RunnerV2LimitBoundaryKindV1::$boundary,
                code: $code,
                stable_name: $name,
            }
        };
    }

    const TEST_BOUNDARIES_EXPECTED_V1: [TestBoundaryExpectedV1; 12] = [
        test_boundary_v1!(1, Zero, 1, "zero"),
        test_boundary_v1!(2, One, 2, "one"),
        test_boundary_v1!(3, StructuralMinimum, 3, "structural-minimum"),
        test_boundary_v1!(
            4,
            OneBelowStructuralMinimum,
            4,
            "one-below-structural-minimum"
        ),
        test_boundary_v1!(5, SmokeCeiling, 5, "smoke-ceiling"),
        test_boundary_v1!(6, SmokeTightened, 6, "smoke-tightened"),
        test_boundary_v1!(7, SmokeOneOver, 7, "smoke-one-over"),
        test_boundary_v1!(8, FullCeiling, 8, "full-ceiling"),
        test_boundary_v1!(9, FullTightened, 9, "full-tightened"),
        test_boundary_v1!(10, FullOneOver, 10, "full-one-over"),
        test_boundary_v1!(11, RepresentationalMaximum, 11, "representational-maximum"),
        test_boundary_v1!(
            12,
            CheckedRepresentationalOverflowRefusal,
            12,
            "checked-representational-overflow-refusal"
        ),
    ];

    const _: [(); 12] = [(); TEST_BOUNDARIES_EXPECTED_V1.len()];

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestMetaExpectedV1 {
        ordinal: u16,
        stable_id: &'static str,
        id_suffix: &'static str,
        group: RunnerV2StageACellGroupV1,
        operation: RunnerV2StageAMetaOperationV1,
        outcome: RunnerV2RawOutcomeKindV1,
        reason: RunnerV2RawReasonV1,
        partition: RunnerV2StageAExpectedPartitionV1,
    }

    macro_rules! test_meta_v1 {
        (
            $ordinal:literal,
            $id:literal,
            $suffix:literal,
            $group:ident,
            $operation:ident,
            $outcome:ident,
            $reason:ident,
            $partition:ident
        ) => {
            TestMetaExpectedV1 {
                ordinal: $ordinal,
                stable_id: $id,
                id_suffix: $suffix,
                group: RunnerV2StageACellGroupV1::$group,
                operation: RunnerV2StageAMetaOperationV1::$operation,
                outcome: RunnerV2RawOutcomeKindV1::$outcome,
                reason: RunnerV2RawReasonV1::$reason,
                partition: RunnerV2StageAExpectedPartitionV1::$partition,
            }
        };
    }

    const TEST_META_EXPECTED_V1: [TestMetaExpectedV1; 15] = [
        test_meta_v1!(
            1,
            "runner-v2.base-values.meta-001-typed-absence-distinct-from-zero.v1",
            "typed-absence-distinct-from-zero",
            LiteralUnitBoundary,
            TypedAbsenceDistinctFromZero,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            2,
            "runner-v2.base-values.meta-002-f32-named-total-order.v1",
            "f32-named-total-order",
            PropertyMetamorphic,
            F32NamedTotalOrder,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            3,
            "runner-v2.base-values.meta-003-f64-named-total-order.v1",
            "f64-named-total-order",
            PropertyMetamorphic,
            F64NamedTotalOrder,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            4,
            "runner-v2.base-values.meta-004-capability-none-contract.v1",
            "capability-none-contract",
            StateModel,
            CapabilityNoneContract,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            5,
            "runner-v2.base-values.meta-005-common-requirements-exact.v1",
            "common-requirements-exact",
            StateModel,
            CommonRequirementsExact,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            6,
            "runner-v2.base-values.meta-006-common-requirement-reordered-refusal.v1",
            "common-requirement-reordered-refusal",
            MutationFuzz,
            CommonRequirementReorderedRefusal,
            Refused,
            ExactMembershipMismatch,
            Mutation
        ),
        test_meta_v1!(
            7,
            "runner-v2.base-values.meta-007-future-sources-exact.v1",
            "future-sources-exact",
            SourceClosure,
            FutureSourcesExact,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            8,
            "runner-v2.base-values.meta-008-rootless-ac58.v1",
            "rootless-ac58",
            ApiCompileFail,
            RootlessAc58,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            9,
            "runner-v2.base-values.meta-009-owner-source-fragment.v1",
            "owner-source-fragment",
            SourceClosure,
            OwnerSourceFragment,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            10,
            "runner-v2.base-values.meta-010-local-route.v1",
            "local-route",
            NoMockLocalIntegration,
            LocalRoute,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            11,
            "runner-v2.base-values.meta-011-diagnostic-redaction.v1",
            "diagnostic-redaction",
            DetailedLoggingRedaction,
            DiagnosticRedaction,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            12,
            "runner-v2.base-values.meta-012-reproduction-declaration.v1",
            "reproduction-declaration",
            ReproductionDeclaration,
            ReproductionDeclaration,
            Accepted,
            ExactCheckedValue,
            EligiblePositive
        ),
        test_meta_v1!(
            13,
            "runner-v2.base-values.meta-013-compile-fail-ordering-surface.v1",
            "compile-fail-ordering-surface",
            ApiCompileFail,
            CompileFailOrderingSurface,
            Inapplicable,
            PureDeclarationFacet,
            Inapplicable
        ),
        test_meta_v1!(
            14,
            "runner-v2.base-values.meta-014-shard-inapplicable.v1",
            "shard-inapplicable",
            FaultResourceCancellation,
            ShardInapplicable,
            Inapplicable,
            ShardInapplicable,
            Inapplicable
        ),
        test_meta_v1!(
            15,
            "runner-v2.base-values.meta-015-resume-inapplicable.v1",
            "resume-inapplicable",
            FaultResourceCancellation,
            ResumeInapplicable,
            Inapplicable,
            ResumeInapplicable,
            Inapplicable
        ),
    ];

    const _: [(); 15] = [(); TEST_META_EXPECTED_V1.len()];

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestCommonExpectedV1 {
        ordinal: u16,
        slot_id: &'static str,
        realization_owner: &'static str,
        role: &'static str,
        domain: &'static str,
        plane_mask: u8,
        stage: RunnerV2CommonFulfillmentStageV1,
    }

    macro_rules! test_common_v1 {
        (
            $ordinal:literal,
            $slot:literal,
            $owner:literal,
            $role:literal,
            $domain:literal,
            $mask:literal,
            $stage:ident
        ) => {
            TestCommonExpectedV1 {
                ordinal: $ordinal,
                slot_id: $slot,
                realization_owner: $owner,
                role: $role,
                domain: $domain,
                plane_mask: $mask,
                stage: RunnerV2CommonFulfillmentStageV1::$stage,
            }
        };
    }

    const TEST_COMMON_EXPECTED_V1: [TestCommonExpectedV1; 31] = [
        test_common_v1!(
            1,
            "runner-v2.common.attempt-identity-contract.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "attempt-identity-root-v1",
            "runner-v2-attempt-identity-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            2,
            "runner-v2.common.canonical-runtime-observation-projection.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "canonical-runtime-observation-projection-root-v1",
            "runner-v2-canonical-runtime-observation-v1",
            0b111,
            RuntimeEvidence
        ),
        test_common_v1!(
            3,
            "runner-v2.common.runtime-evidence-envelope.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "runtime-evidence-envelope-root-v1",
            "runner-v2-runtime-evidence-envelope-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            4,
            "runner-v2.common.actual-five-explicits.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "actual-five-explicits-root-v1",
            "runner-v2-actual-five-explicits-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            5,
            "runner-v2.common.completeness-disposition.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "completeness-disposition-root-v1",
            "runner-v2-completeness-disposition-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            6,
            "runner-v2.common.safe-partial-evidence.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "safe-partial-evidence-root-v1",
            "runner-v2-safe-partial-evidence-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            7,
            "runner-v2.common.capability-reconciliation.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "capability-reconciliation-root-v1",
            "runner-v2-capability-reconciliation-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            8,
            "runner-v2.common.resource-reconciliation.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "resource-reconciliation-root-v1",
            "runner-v2-resource-reconciliation-v1",
            0b110,
            RuntimeEvidence
        ),
        test_common_v1!(
            9,
            "runner-v2.common.attempt-retention-receipt.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "attempt-retention-receipt-root-v1",
            "runner-v2-attempt-retention-receipt-v1",
            0b100,
            RuntimeEvidence
        ),
        test_common_v1!(
            10,
            "runner-v2.common.atomic-retention-finalization.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.4",
            "atomic-retention-finalization-root-v1",
            "runner-v2-atomic-retention-finalization-v1",
            0b100,
            RuntimeEvidence
        ),
        test_common_v1!(
            11,
            "runner-v2.common.route-schema.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "route-schema-root-v1",
            "runner-v2-route-schema-v1",
            0b011,
            RoutesAndDispatch
        ),
        test_common_v1!(
            12,
            "runner-v2.common.owner-matrix.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "owner-matrix-root-v1",
            "runner-v2-owner-matrix-v1",
            0b011,
            RoutesAndDispatch
        ),
        test_common_v1!(
            13,
            "runner-v2.common.route-registry.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "route-registry-root-v1",
            "runner-v2-route-registry-v1",
            0b011,
            RoutesAndDispatch
        ),
        test_common_v1!(
            14,
            "runner-v2.common.deferred-route-registry.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "deferred-route-registry-root-v1",
            "runner-v2-deferred-route-registry-v1",
            0b011,
            RoutesAndDispatch
        ),
        test_common_v1!(
            15,
            "runner-v2.common.dispatcher.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "dispatcher-root-v1",
            "runner-v2-dispatcher-v1",
            0b010,
            RoutesAndDispatch
        ),
        test_common_v1!(
            16,
            "runner-v2.common.native-bootstrap.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "native-bootstrap-root-v1",
            "runner-v2-native-bootstrap-v1",
            0b010,
            RoutesAndDispatch
        ),
        test_common_v1!(
            17,
            "runner-v2.common.execution-source-binding.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "execution-source-binding-root-v1",
            "runner-v2-execution-source-binding-v1",
            0b110,
            RoutesAndDispatch
        ),
        test_common_v1!(
            18,
            "runner-v2.common.retention-scope.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "retention-scope-root-v1",
            "runner-v2-retention-scope-v1",
            0b100,
            RoutesAndDispatch
        ),
        test_common_v1!(
            19,
            "runner-v2.common.finalization-protocol.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "finalization-protocol-root-v1",
            "runner-v2-finalization-protocol-v1",
            0b100,
            RoutesAndDispatch
        ),
        test_common_v1!(
            20,
            "runner-v2.common.shard-contract.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "shard-contract-root-v1",
            "runner-v2-shard-contract-v1",
            0b111,
            RoutesAndDispatch
        ),
        test_common_v1!(
            21,
            "runner-v2.common.resume-contract.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.5",
            "resume-contract-root-v1",
            "runner-v2-resume-contract-v1",
            0b111,
            RoutesAndDispatch
        ),
        test_common_v1!(
            22,
            "runner-v2.common.command-schema.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "command-schema-root-v1",
            "runner-v2-command-schema-v1",
            0b011,
            LoggingAndReproduction
        ),
        test_common_v1!(
            23,
            "runner-v2.common.jsonl-event-schema.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "jsonl-event-schema-root-v1",
            "runner-v2-jsonl-event-schema-v1",
            0b111,
            LoggingAndReproduction
        ),
        test_common_v1!(
            24,
            "runner-v2.common.terminal-reservation.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "terminal-reservation-root-v1",
            "runner-v2-terminal-reservation-v1",
            0b110,
            LoggingAndReproduction
        ),
        test_common_v1!(
            25,
            "runner-v2.common.redaction-policy.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "redaction-policy-root-v1",
            "runner-v2-redaction-policy-v1",
            0b111,
            LoggingAndReproduction
        ),
        test_common_v1!(
            26,
            "runner-v2.common.first-divergence.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "first-divergence-root-v1",
            "runner-v2-first-divergence-v1",
            0b111,
            LoggingAndReproduction
        ),
        test_common_v1!(
            27,
            "runner-v2.common.reproduction-schema.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "reproduction-schema-root-v1",
            "runner-v2-reproduction-schema-v1",
            0b111,
            LoggingAndReproduction
        ),
        test_common_v1!(
            28,
            "runner-v2.common.relative-artifact-schema.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "relative-artifact-schema-root-v1",
            "runner-v2-relative-artifact-schema-v1",
            0b111,
            LoggingAndReproduction
        ),
        test_common_v1!(
            29,
            "runner-v2.common.raw-audit-binding.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "raw-audit-binding-root-v1",
            "runner-v2-raw-audit-binding-v1",
            0b100,
            LoggingAndReproduction
        ),
        test_common_v1!(
            30,
            "runner-v2.common.stage-telemetry.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "stage-telemetry-root-v1",
            "runner-v2-stage-telemetry-v1",
            0b100,
            LoggingAndReproduction
        ),
        test_common_v1!(
            31,
            "runner-v2.common.operator-view.v1",
            "frankensim-epic-foundations-huq.24.1.1.1.6",
            "operator-view-root-v1",
            "runner-v2-operator-view-v1",
            0b100,
            LoggingAndReproduction
        ),
    ];

    const _: [(); 31] = [(); TEST_COMMON_EXPECTED_V1.len()];

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestFutureSourceExpectedV1 {
        final_ordinal: u16,
        path: &'static str,
    }

    const TEST_FUTURE_SOURCES_EXPECTED_V1: [TestFutureSourceExpectedV1; 13] = [
        TestFutureSourceExpectedV1 {
            final_ordinal: 28,
            path: "crates/fs-evidence-runner/src/runner_v2.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 29,
            path: "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 30,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 31,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 32,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/diagnostics.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 33,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/schema_registry.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 34,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/runtime_evidence.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 35,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/routes.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 36,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/detailed_logging.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 37,
            path: "crates/fs-evidence-runner/src/runner_v2/work_packages/execution.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 38,
            path: "crates/fs-evidence-runner/tests/runner_v2_base_work_packages.rs",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 39,
            path: "scripts/ci/runner_v2_base_work_packages_e2e.sh",
        },
        TestFutureSourceExpectedV1 {
            final_ordinal: 40,
            path: "scripts/ci/runner_v2_base_work_packages_e2e.ps1",
        },
    ];

    const TEST_OWNER_SOURCE_PATHS_EXPECTED_V1: [&str; 2] = [
        "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
    ];

    const TEST_DEPENDENCY_SOURCE_PATHS_EXPECTED_V1: [&str; 16] = [
        "crates/fs-blake3/src/lib.rs",
        "crates/fs-evidence-runner/src/lib.rs",
        "crates/fs-evidence-runner/src/canonical.rs",
        "crates/fs-evidence-runner/src/catalog.rs",
        "crates/fs-evidence-runner/src/construction.rs",
        "crates/fs-evidence-runner/src/coverage.rs",
        "crates/fs-evidence-runner/src/identity.rs",
        "crates/fs-evidence-runner/src/limits.rs",
        "crates/fs-evidence-runner/src/path.rs",
        "crates/fs-evidence-runner/src/projection.rs",
        "crates/fs-evidence-runner/src/schema_impact.rs",
        "crates/fs-evidence-runner/src/value.rs",
        "crates/fs-evidence-runner/src/runner_v2.rs",
        "crates/fs-evidence-runner/src/runner_v2/handoff.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages.rs",
        "crates/fs-evidence-runner/src/runner_v2/work_packages/base_values.rs",
    ];

    const _: [(); 13] = [(); TEST_FUTURE_SOURCES_EXPECTED_V1.len()];
    const _: [(); 2] = [(); TEST_OWNER_SOURCE_PATHS_EXPECTED_V1.len()];
    const _: [(); 16] = [(); TEST_DEPENDENCY_SOURCE_PATHS_EXPECTED_V1.len()];

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestLimitMutationExpectedV1 {
        ordinal: u16,
        stable_id: &'static str,
        field: RunnerLimitFieldV2,
        field_name: &'static str,
        declared_width: RunnerLimitWidthV2,
        opposite_zero: RunnerLimitValueV2,
        unit: RunnerLimitUnitV2,
    }

    macro_rules! test_limit_mutation_v1 {
        (
            $ordinal:literal,
            $id:literal,
            $field:ident,
            $name:literal,
            $width:ident,
            $opposite:ident,
            $unit:ident
        ) => {
            TestLimitMutationExpectedV1 {
                ordinal: $ordinal,
                stable_id: $id,
                field: RunnerLimitFieldV2::$field,
                field_name: $name,
                declared_width: RunnerLimitWidthV2::$width,
                opposite_zero: RunnerLimitValueV2::$opposite(0),
                unit: RunnerLimitUnitV2::$unit,
            }
        };
    }

    const TEST_LIMIT_MUTATIONS_EXPECTED_V1: [TestLimitMutationExpectedV1; 71] = [
        test_limit_mutation_v1!(
            1,
            "runner-v2.base-values.limit-001.wrong-width-mutation.v1",
            ArgvTokens,
            "argv_tokens",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            2,
            "runner-v2.base-values.limit-002.wrong-width-mutation.v1",
            ArgvTokenBytes,
            "argv_token_bytes",
            U64,
            U32,
            LogicalBytes
        ),
        test_limit_mutation_v1!(
            3,
            "runner-v2.base-values.limit-003.wrong-width-mutation.v1",
            ArgvAggregateBytes,
            "argv_aggregate_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            4,
            "runner-v2.base-values.limit-004.wrong-width-mutation.v1",
            LifecycleRecordEncodedBytes,
            "lifecycle_record_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            5,
            "runner-v2.base-values.limit-005.wrong-width-mutation.v1",
            CaseLifecycleRecords,
            "case_lifecycle_records",
            U32,
            U64,
            Records
        ),
        test_limit_mutation_v1!(
            6,
            "runner-v2.base-values.limit-006.wrong-width-mutation.v1",
            CaseLifecycleEncodedBytes,
            "case_lifecycle_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            7,
            "runner-v2.base-values.limit-007.wrong-width-mutation.v1",
            FamilyRowsPerCase,
            "family_rows_per_case",
            U32,
            U64,
            Rows
        ),
        test_limit_mutation_v1!(
            8,
            "runner-v2.base-values.limit-008.wrong-width-mutation.v1",
            InvocationCases,
            "invocation_cases",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            9,
            "runner-v2.base-values.limit-009.wrong-width-mutation.v1",
            LifecycleDocumentRecords,
            "lifecycle_document_records",
            U32,
            U64,
            Records
        ),
        test_limit_mutation_v1!(
            10,
            "runner-v2.base-values.limit-010.wrong-width-mutation.v1",
            LifecycleDocumentEncodedBytes,
            "lifecycle_document_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            11,
            "runner-v2.base-values.limit-011.wrong-width-mutation.v1",
            CommandResultStdoutBytes,
            "command_result_stdout_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            12,
            "runner-v2.base-values.limit-012.wrong-width-mutation.v1",
            ChildStdoutBytes,
            "child_stdout_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            13,
            "runner-v2.base-values.limit-013.wrong-width-mutation.v1",
            CombinedChildStdoutBytes,
            "combined_child_stdout_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            14,
            "runner-v2.base-values.limit-014.wrong-width-mutation.v1",
            ChildStderrBytes,
            "child_stderr_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            15,
            "runner-v2.base-values.limit-015.wrong-width-mutation.v1",
            CombinedChildStderrBytes,
            "combined_child_stderr_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            16,
            "runner-v2.base-values.limit-016.wrong-width-mutation.v1",
            ManifestEncodedBytes,
            "manifest_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            17,
            "runner-v2.base-values.limit-017.wrong-width-mutation.v1",
            NestingDepth,
            "nesting_depth",
            U32,
            U64,
            Depth
        ),
        test_limit_mutation_v1!(
            18,
            "runner-v2.base-values.limit-018.wrong-width-mutation.v1",
            ComparisonNodes,
            "comparison_nodes",
            U32,
            U64,
            Nodes
        ),
        test_limit_mutation_v1!(
            19,
            "runner-v2.base-values.limit-019.wrong-width-mutation.v1",
            EffectNodes,
            "effect_nodes",
            U32,
            U64,
            Nodes
        ),
        test_limit_mutation_v1!(
            20,
            "runner-v2.base-values.limit-020.wrong-width-mutation.v1",
            TextBytes,
            "text_bytes",
            U64,
            U32,
            LogicalBytes
        ),
        test_limit_mutation_v1!(
            21,
            "runner-v2.base-values.limit-021.wrong-width-mutation.v1",
            StableTokenBytes,
            "stable_token_bytes",
            U64,
            U32,
            LogicalBytes
        ),
        test_limit_mutation_v1!(
            22,
            "runner-v2.base-values.limit-022.wrong-width-mutation.v1",
            BundleRelativePathBytes,
            "bundle_relative_path_bytes",
            U64,
            U32,
            LogicalBytes
        ),
        test_limit_mutation_v1!(
            23,
            "runner-v2.base-values.limit-023.wrong-width-mutation.v1",
            DiagnosticsPerCase,
            "diagnostics_per_case",
            U32,
            U64,
            Diagnostics
        ),
        test_limit_mutation_v1!(
            24,
            "runner-v2.base-values.limit-024.wrong-width-mutation.v1",
            DiagnosticsPerRun,
            "diagnostics_per_run",
            U32,
            U64,
            Diagnostics
        ),
        test_limit_mutation_v1!(
            25,
            "runner-v2.base-values.limit-025.wrong-width-mutation.v1",
            PrerequisitesPerDiagnostic,
            "prerequisites_per_diagnostic",
            U32,
            U64,
            Prerequisites
        ),
        test_limit_mutation_v1!(
            26,
            "runner-v2.base-values.limit-026.wrong-width-mutation.v1",
            RepairsPerDiagnostic,
            "repairs_per_diagnostic",
            U32,
            U64,
            Repairs
        ),
        test_limit_mutation_v1!(
            27,
            "runner-v2.base-values.limit-027.wrong-width-mutation.v1",
            Artifacts,
            "artifacts",
            U32,
            U64,
            Artifacts
        ),
        test_limit_mutation_v1!(
            28,
            "runner-v2.base-values.limit-028.wrong-width-mutation.v1",
            ArtifactEncodedBytes,
            "artifact_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            29,
            "runner-v2.base-values.limit-029.wrong-width-mutation.v1",
            ArtifactExpandedBytes,
            "artifact_expanded_bytes",
            U64,
            U32,
            ExpandedBytes
        ),
        test_limit_mutation_v1!(
            30,
            "runner-v2.base-values.limit-030.wrong-width-mutation.v1",
            ArtifactStoredBytes,
            "artifact_stored_bytes",
            U64,
            U32,
            StoredBytes
        ),
        test_limit_mutation_v1!(
            31,
            "runner-v2.base-values.limit-031.wrong-width-mutation.v1",
            BundleEncodedBytes,
            "bundle_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            32,
            "runner-v2.base-values.limit-032.wrong-width-mutation.v1",
            BundleExpandedBytes,
            "bundle_expanded_bytes",
            U64,
            U32,
            ExpandedBytes
        ),
        test_limit_mutation_v1!(
            33,
            "runner-v2.base-values.limit-033.wrong-width-mutation.v1",
            ArtifactStoredAggregateBytes,
            "artifact_stored_aggregate_bytes",
            U64,
            U32,
            StoredBytes
        ),
        test_limit_mutation_v1!(
            34,
            "runner-v2.base-values.limit-034.wrong-width-mutation.v1",
            SystemPublicationStoredBytes,
            "system_publication_stored_bytes",
            U64,
            U32,
            StoredBytes
        ),
        test_limit_mutation_v1!(
            35,
            "runner-v2.base-values.limit-035.wrong-width-mutation.v1",
            PublicationStoredBytes,
            "publication_stored_bytes",
            U64,
            U32,
            StoredBytes
        ),
        test_limit_mutation_v1!(
            36,
            "runner-v2.base-values.limit-036.wrong-width-mutation.v1",
            ChildStreamDiscardBytes,
            "child_stream_discard_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            37,
            "runner-v2.base-values.limit-037.wrong-width-mutation.v1",
            ModesPerFamily,
            "modes_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            38,
            "runner-v2.base-values.limit-038.wrong-width-mutation.v1",
            ExtensionDiagnosticsPerFamily,
            "extension_diagnostics_per_family",
            U32,
            U64,
            Diagnostics
        ),
        test_limit_mutation_v1!(
            39,
            "runner-v2.base-values.limit-039.wrong-width-mutation.v1",
            ArtifactRolesPerFamily,
            "artifact_roles_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            40,
            "runner-v2.base-values.limit-040.wrong-width-mutation.v1",
            RootPoliciesPerFamily,
            "root_policies_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            41,
            "runner-v2.base-values.limit-041.wrong-width-mutation.v1",
            RegisteredUnitsPerFamily,
            "registered_units_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            42,
            "runner-v2.base-values.limit-042.wrong-width-mutation.v1",
            DigestDomainsPerFamily,
            "digest_domains_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            43,
            "runner-v2.base-values.limit-043.wrong-width-mutation.v1",
            ExtensionSchemasPerFamily,
            "extension_schemas_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            44,
            "runner-v2.base-values.limit-044.wrong-width-mutation.v1",
            ExecutableDescriptorsPerFamily,
            "executable_descriptors_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            45,
            "runner-v2.base-values.limit-045.wrong-width-mutation.v1",
            MapEntries,
            "map_entries",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            46,
            "runner-v2.base-values.limit-046.wrong-width-mutation.v1",
            GenericArrayItems,
            "generic_array_items",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            47,
            "runner-v2.base-values.limit-047.wrong-width-mutation.v1",
            PathSegments,
            "path_segments",
            U32,
            U64,
            Segments
        ),
        test_limit_mutation_v1!(
            48,
            "runner-v2.base-values.limit-048.wrong-width-mutation.v1",
            IntegerDigits,
            "integer_digits",
            U32,
            U64,
            Digits
        ),
        test_limit_mutation_v1!(
            49,
            "runner-v2.base-values.limit-049.wrong-width-mutation.v1",
            RationalComponentBytes,
            "rational_component_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            50,
            "runner-v2.base-values.limit-050.wrong-width-mutation.v1",
            DecimalCoefficientBytes,
            "decimal_coefficient_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            51,
            "runner-v2.base-values.limit-051.wrong-width-mutation.v1",
            DecimalAbsoluteScale,
            "decimal_absolute_scale",
            U32,
            U64,
            DecimalScale
        ),
        test_limit_mutation_v1!(
            52,
            "runner-v2.base-values.limit-052.wrong-width-mutation.v1",
            LogicalExtentsPerArtifact,
            "logical_extents_per_artifact",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            53,
            "runner-v2.base-values.limit-053.wrong-width-mutation.v1",
            ObservationKeysPerCase,
            "observation_keys_per_case",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            54,
            "runner-v2.base-values.limit-054.wrong-width-mutation.v1",
            DecisionDetailNamespaces,
            "decision_detail_namespaces",
            U32,
            U64,
            Namespaces
        ),
        test_limit_mutation_v1!(
            55,
            "runner-v2.base-values.limit-055.wrong-width-mutation.v1",
            OutputClasses,
            "output_classes",
            U32,
            U64,
            Classes
        ),
        test_limit_mutation_v1!(
            56,
            "runner-v2.base-values.limit-056.wrong-width-mutation.v1",
            OpaqueValueBytes,
            "opaque_value_bytes",
            U64,
            U32,
            LogicalBytes
        ),
        test_limit_mutation_v1!(
            57,
            "runner-v2.base-values.limit-057.wrong-width-mutation.v1",
            RetainedUnknownExtensionBytes,
            "retained_unknown_extension_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            58,
            "runner-v2.base-values.limit-058.wrong-width-mutation.v1",
            ExpressionEdges,
            "expression_edges",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            59,
            "runner-v2.base-values.limit-059.wrong-width-mutation.v1",
            MemoizedEvaluationVisits,
            "memoized_evaluation_visits",
            U32,
            U64,
            Visits
        ),
        test_limit_mutation_v1!(
            60,
            "runner-v2.base-values.limit-060.wrong-width-mutation.v1",
            RepairActionEncodedBytes,
            "repair_action_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            61,
            "runner-v2.base-values.limit-061.wrong-width-mutation.v1",
            ActionableDiagnosticEncodedBytes,
            "actionable_diagnostic_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            62,
            "runner-v2.base-values.limit-062.wrong-width-mutation.v1",
            FailureStderrEncodedBytes,
            "failure_stderr_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            63,
            "runner-v2.base-values.limit-063.wrong-width-mutation.v1",
            RunnerCatalogEncodedBytes,
            "runner_catalog_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            64,
            "runner-v2.base-values.limit-064.wrong-width-mutation.v1",
            PublishedBundleReceiptEncodedBytes,
            "published_bundle_receipt_encoded_bytes",
            U64,
            U32,
            EncodedBytes
        ),
        test_limit_mutation_v1!(
            65,
            "runner-v2.base-values.limit-065.wrong-width-mutation.v1",
            ContentStoreEnvelopeNonPayloadBytes,
            "content_store_envelope_non_payload_bytes",
            U64,
            U32,
            StoredBytes
        ),
        test_limit_mutation_v1!(
            66,
            "runner-v2.base-values.limit-066.wrong-width-mutation.v1",
            RegisteredExtentAxesPerFamily,
            "registered_extent_axes_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            67,
            "runner-v2.base-values.limit-067.wrong-width-mutation.v1",
            RegisteredObservationKeysPerFamily,
            "registered_observation_keys_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            68,
            "runner-v2.base-values.limit-068.wrong-width-mutation.v1",
            RegisteredAuthorityScopesPerFamily,
            "registered_authority_scopes_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            69,
            "runner-v2.base-values.limit-069.wrong-width-mutation.v1",
            RegisteredExternalRootClassesPerFamily,
            "registered_external_root_classes_per_family",
            U32,
            U64,
            Classes
        ),
        test_limit_mutation_v1!(
            70,
            "runner-v2.base-values.limit-070.wrong-width-mutation.v1",
            RegisteredEvaluationUnitsPerFamily,
            "registered_evaluation_units_per_family",
            U32,
            U64,
            Count
        ),
        test_limit_mutation_v1!(
            71,
            "runner-v2.base-values.limit-071.wrong-width-mutation.v1",
            RegisteredResourceIdentitiesPerFamily,
            "registered_resource_identities_per_family",
            U32,
            U64,
            Count
        ),
    ];

    const _: [(); 71] = [(); TEST_LIMIT_MUTATIONS_EXPECTED_V1.len()];

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct TestRetainedExpectedV1 {
        ordinal: u16,
        stable_id: &'static str,
        facet: RunnerV2RetainedDomainFacetV1,
    }

    macro_rules! test_retained_v1 {
        ($ordinal:literal, $id:literal, $facet:ident) => {
            TestRetainedExpectedV1 {
                ordinal: $ordinal,
                stable_id: $id,
                facet: RunnerV2RetainedDomainFacetV1::$facet,
            }
        };
    }

    const TEST_RETAINED_EXPECTED_V1: [TestRetainedExpectedV1; 50] = [
        test_retained_v1!(
            1,
            "runner-v2.base-values.retained-001-signed-integer-width-and-extremes.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            2,
            "runner-v2.base-values.retained-002-unsigned-integer-width-and-extremes.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            3,
            "runner-v2.base-values.retained-003-integer-two-to-the-fifty-three-boundaries.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            4,
            "runner-v2.base-values.retained-004-rational-canonical-normalization.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            5,
            "runner-v2.base-values.retained-005-rational-zero-denominator-refusal.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            6,
            "runner-v2.base-values.retained-006-decimal-canonical-normalization.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            7,
            "runner-v2.base-values.retained-007-decimal-scale-minimum-maximum-and-one-over.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            8,
            "runner-v2.base-values.retained-008-binary32-exact-bit-identity.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            9,
            "runner-v2.base-values.retained-009-binary64-exact-bit-identity.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            10,
            "runner-v2.base-values.retained-010-binary32-named-ieee-total-order.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            11,
            "runner-v2.base-values.retained-011-binary64-named-ieee-total-order.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            12,
            "runner-v2.base-values.retained-012-ieee-signed-zero-nan-infinity-and-subnormal.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            13,
            "runner-v2.base-values.retained-013-duration-byte-and-rate-literal-catalogs.v1",
            NumericLiteral
        ),
        test_retained_v1!(
            14,
            "runner-v2.base-values.retained-014-physical-unit-positive-canonical-scale.v1",
            Unit
        ),
        test_retained_v1!(
            15,
            "runner-v2.base-values.retained-015-unit-dimension-and-scale-equivalence.v1",
            Unit
        ),
        test_retained_v1!(
            16,
            "runner-v2.base-values.retained-016-logical-unit-catalog-exact-order.v1",
            Unit
        ),
        test_retained_v1!(
            17,
            "runner-v2.base-values.retained-017-physical-logical-unit-substitution-refusal.v1",
            Unit
        ),
        test_retained_v1!(
            18,
            "runner-v2.base-values.retained-018-stable-token-empty-minimum-maximum-one-over.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            19,
            "runner-v2.base-values.retained-019-stable-token-separator-canonicality.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            20,
            "runner-v2.base-values.retained-020-text-empty-minimum-maximum-one-over.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            21,
            "runner-v2.base-values.retained-021-opaque-bytes-empty-maximum-one-over.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            22,
            "runner-v2.base-values.retained-022-logical-path-segment-empty-maximum-one-over.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            23,
            "runner-v2.base-values.retained-023-path-alias-and-normalization-refusal.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            24,
            "runner-v2.base-values.retained-024-reserved-prefix-refusal.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            25,
            "runner-v2.base-values.retained-025-non-ascii-and-platform-alias-boundaries.v1",
            TokenTextPath
        ),
        test_retained_v1!(
            26,
            "runner-v2.base-values.retained-026-closed-literal-catalog-code-name-order.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            27,
            "runner-v2.base-values.retained-027-unknown-catalog-code-refusal.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            28,
            "runner-v2.base-values.retained-028-nominal-root-role-and-domain-binding.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            29,
            "runner-v2.base-values.retained-029-raw-hash-as-nominal-root-refusal.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            30,
            "runner-v2.base-values.retained-030-typed-absence-versus-present-zero-digest.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            31,
            "runner-v2.base-values.retained-031-source-build-toolchain-role-substitution-refusal.v1",
            CatalogAndNominalIdentity
        ),
        test_retained_v1!(
            32,
            "runner-v2.base-values.retained-032-deterministic-construction-repeatability.v1",
            PropertyAndMetamorphic
        ),
        test_retained_v1!(
            33,
            "runner-v2.base-values.retained-033-canonical-round-trip-identity.v1",
            PropertyAndMetamorphic
        ),
        test_retained_v1!(
            34,
            "runner-v2.base-values.retained-034-unit-rescaling-invariance.v1",
            PropertyAndMetamorphic
        ),
        test_retained_v1!(
            35,
            "runner-v2.base-values.retained-035-named-ieee-order-repeatability.v1",
            PropertyAndMetamorphic
        ),
        test_retained_v1!(
            36,
            "runner-v2.base-values.retained-036-public-constructor-malformed-input-refusal.v1",
            MutationAndRefusal
        ),
        test_retained_v1!(
            37,
            "runner-v2.base-values.retained-037-missing-extra-duplicate-and-reordered-refusal.v1",
            MutationAndRefusal
        ),
        test_retained_v1!(
            38,
            "runner-v2.base-values.retained-038-cross-role-and-cross-domain-substitution-refusal.v1",
            MutationAndRefusal
        ),
        test_retained_v1!(
            39,
            "runner-v2.base-values.retained-039-profile-version-source-and-feature-drift.v1",
            MutationAndRefusal
        ),
        test_retained_v1!(
            40,
            "runner-v2.base-values.retained-040-ambient-unpinned-and-mixed-snapshot-refusal.v1",
            MutationAndRefusal
        ),
        test_retained_v1!(
            41,
            "runner-v2.base-values.retained-041-validated-wrapper-private-field-boundaries.v1",
            ApiAndCompileFail
        ),
        test_retained_v1!(
            42,
            "runner-v2.base-values.retained-042-binary32-no-ord-or-partial-ord.v1",
            ApiAndCompileFail
        ),
        test_retained_v1!(
            43,
            "runner-v2.base-values.retained-043-binary64-no-ord-or-partial-ord.v1",
            ApiAndCompileFail
        ),
        test_retained_v1!(
            44,
            "runner-v2.base-values.retained-044-ieee-wrapper-no-ordered-collection-or-sort.v1",
            ApiAndCompileFail
        ),
        test_retained_v1!(
            45,
            "runner-v2.base-values.retained-045-rootless-handoff-no-canonical-or-authority-surface.v1",
            ApiAndCompileFail
        ),
        test_retained_v1!(
            46,
            "runner-v2.base-values.retained-046-checked-arithmetic-overflow-refusal.v1",
            FaultResourceAndIntegration
        ),
        test_retained_v1!(
            47,
            "runner-v2.base-values.retained-047-bounded-allocation-and-one-over-refusal.v1",
            FaultResourceAndIntegration
        ),
        test_retained_v1!(
            48,
            "runner-v2.base-values.retained-048-real-no-mock-local-evaluator-integration.v1",
            FaultResourceAndIntegration
        ),
        test_retained_v1!(
            49,
            "runner-v2.base-values.retained-049-diagnostic-redaction-and-forbidden-value-no-echo.v1",
            FaultResourceAndIntegration
        ),
        test_retained_v1!(
            50,
            "runner-v2.base-values.retained-050-reproduction-and-compatible-source-closure-declaration.v1",
            FaultResourceAndIntegration
        ),
    ];

    const _: [(); 50] = [(); TEST_RETAINED_EXPECTED_V1.len()];
    const _: [(); 50] = [(); 13 + 4 + 8 + 6 + 4 + 5 + 5 + 5];

    fn assert_refuses<T>(result: Result<T, ConstructionErrorV2>) {
        assert!(result.is_err(), "mutated exact-set input must refuse");
    }

    fn assert_inventory_mismatch_v1(
        result: Result<(), RunnerV2StageAInventoryMismatchV1>,
        kind: ConstructionErrorKindV2,
        inventory: &'static str,
        repair_target: &'static str,
        index0: usize,
        expected_count: usize,
        observed_count: usize,
        expected_identity: &str,
        observed_safe_identity: &str,
        observed_redacted: bool,
    ) {
        let mismatch = result.expect_err("mutated exact inventory must refuse");
        assert_eq!(mismatch.kind(), kind);
        assert_eq!(mismatch.inventory(), inventory);
        assert_eq!(mismatch.first_mismatch_index0(), index0);
        assert_eq!(mismatch.expected_ordinal1(), index0 + 1);
        assert_eq!(mismatch.expected_identity(), expected_identity);
        assert_eq!(mismatch.observed_safe_identity(), observed_safe_identity);
        assert_eq!(mismatch.observed_identity_redacted(), observed_redacted);
        assert_eq!(mismatch.component(), "identity");
        assert_eq!(mismatch.expected_semantic_value(), expected_identity);
        assert_eq!(mismatch.observed_safe_value(), observed_safe_identity);
        assert_eq!(mismatch.observed_value_redacted(), observed_redacted);
        assert_eq!(
            mismatch.semantic_owner(),
            STAGE_A_INVENTORY_SEMANTIC_OWNER_V1
        );
        assert_eq!(mismatch.expected_count(), expected_count);
        assert_eq!(mismatch.observed_count(), observed_count);
        assert_eq!(mismatch.repairs().len(), 1);
        assert_eq!(mismatch.repairs()[0].rank(), 1);
        assert_eq!(
            mismatch.repairs()[0].kind(),
            RepairActionKindV2::ChangeArguments
        );
        assert_eq!(mismatch.repairs()[0].target(), repair_target);
        assert!(
            mismatch
                .to_string()
                .contains(&format!("ordinal {}", index0 + 1))
        );
    }

    fn assert_every_position_inventory_mutations_v1(
        expected: &[String],
        validate: impl Fn(&[String]) -> Result<(), RunnerV2StageAInventoryMismatchV1>,
        inventory: &'static str,
        repair_target: &'static str,
    ) {
        assert!(
            !expected.is_empty(),
            "an exact Stage-A inventory must be nonempty"
        );
        validate(expected).expect("independent literal inventory must be accepted");

        for index in 0..expected.len() {
            let mut missing = expected.to_vec();
            missing.remove(index);
            assert_inventory_mismatch_v1(
                validate(&missing),
                ConstructionErrorKindV2::Missing,
                inventory,
                repair_target,
                index,
                expected.len(),
                expected.len() - 1,
                &expected[index],
                missing
                    .get(index)
                    .map_or(STAGE_A_INVENTORY_MISSING_SENTINEL_V1, String::as_str),
                false,
            );

            let mut substituted = expected.to_vec();
            substituted[index] =
                "runner-v2-stage-a-unregistered-substitution-must-be-redacted".to_owned();
            assert_inventory_mismatch_v1(
                validate(&substituted),
                ConstructionErrorKindV2::Incompatible,
                inventory,
                repair_target,
                index,
                expected.len(),
                expected.len(),
                &expected[index],
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );
        }

        for insertion_index in 0..=expected.len() {
            let mut extra = expected.to_vec();
            extra.insert(
                insertion_index,
                "runner-v2-stage-a-unregistered-extra-must-be-redacted".to_owned(),
            );
            assert_inventory_mismatch_v1(
                validate(&extra),
                ConstructionErrorKindV2::Unexpected,
                inventory,
                repair_target,
                insertion_index,
                expected.len(),
                expected.len() + 1,
                expected
                    .get(insertion_index)
                    .map_or(STAGE_A_INVENTORY_END_SENTINEL_V1, String::as_str),
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );
        }

        for insertion_index in 0..=expected.len() {
            let mut duplicate = expected.to_vec();
            let duplicate_identity = if expected.len() == 1 {
                expected[0].clone()
            } else if insertion_index == expected.len() {
                expected[expected.len() - 1].clone()
            } else {
                expected[(insertion_index + 1) % expected.len()].clone()
            };
            duplicate.insert(insertion_index, duplicate_identity);
            let first_mismatch_index0 = if expected.len() == 1 {
                1
            } else {
                insertion_index
            };
            assert_inventory_mismatch_v1(
                validate(&duplicate),
                ConstructionErrorKindV2::Duplicate,
                inventory,
                repair_target,
                first_mismatch_index0,
                expected.len(),
                expected.len() + 1,
                expected
                    .get(first_mismatch_index0)
                    .map_or(STAGE_A_INVENTORY_END_SENTINEL_V1, String::as_str),
                &duplicate[first_mismatch_index0],
                false,
            );
        }

        if expected.len() == 1 {
            assert!(
                expected.windows(2).next().is_none(),
                "a one-row inventory has no nonidentity reordering"
            );
        } else {
            for index in 0..expected.len() - 1 {
                let mut reordered = expected.to_vec();
                reordered.swap(index, index + 1);
                assert_inventory_mismatch_v1(
                    validate(&reordered),
                    ConstructionErrorKindV2::OutOfOrder,
                    inventory,
                    repair_target,
                    index,
                    expected.len(),
                    expected.len(),
                    &expected[index],
                    &reordered[index],
                    false,
                );
            }
        }
    }

    fn assert_every_position_source_path_mutations_v1<T: Clone>(
        exact: &[T],
        expected_paths: &[&str],
        inventory: &'static str,
        repair_target: &'static str,
        validate: impl Fn(&[T]) -> Result<(), RunnerV2StageAInventoryMismatchV1>,
        path_of: impl Fn(&T) -> &str,
        set_path: impl Fn(&mut T, LogicalBundlePathV1),
    ) {
        assert_eq!(exact.len(), expected_paths.len());
        assert!(exact.len() >= 2);
        validate(exact).expect("exact typed source inventory");

        for index in 0..exact.len() {
            let mut missing = exact.to_vec();
            missing.remove(index);
            assert_inventory_mismatch_v1(
                validate(&missing),
                ConstructionErrorKindV2::Missing,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len() - 1,
                expected_paths[index],
                missing
                    .get(index)
                    .map_or(STAGE_A_INVENTORY_MISSING_SENTINEL_V1, |row| path_of(row)),
                false,
            );

            let mut substituted = exact.to_vec();
            set_path(
                &mut substituted[index],
                stage_a_path(
                    "test.source.path",
                    "crates/fs-evidence-runner/src/runner_v2/unregistered_source_member.rs",
                )
                .expect("valid but unregistered source path"),
            );
            assert_inventory_mismatch_v1(
                validate(&substituted),
                ConstructionErrorKindV2::Incompatible,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len(),
                expected_paths[index],
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );
        }

        for insertion_index in 0..=exact.len() {
            let mut extra_row = exact[0].clone();
            set_path(
                &mut extra_row,
                stage_a_path(
                    "test.source.path",
                    "crates/fs-evidence-runner/src/runner_v2/unregistered_extra_source_member.rs",
                )
                .expect("valid but unregistered extra source path"),
            );
            let mut extra = exact.to_vec();
            extra.insert(insertion_index, extra_row);
            assert_inventory_mismatch_v1(
                validate(&extra),
                ConstructionErrorKindV2::Unexpected,
                inventory,
                repair_target,
                insertion_index,
                exact.len(),
                exact.len() + 1,
                expected_paths
                    .get(insertion_index)
                    .copied()
                    .unwrap_or(STAGE_A_INVENTORY_END_SENTINEL_V1),
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );

            let duplicate_source_index = if insertion_index == exact.len() {
                exact.len() - 1
            } else {
                (insertion_index + 1) % exact.len()
            };
            let mut duplicate = exact.to_vec();
            duplicate.insert(insertion_index, exact[duplicate_source_index].clone());
            assert_inventory_mismatch_v1(
                validate(&duplicate),
                ConstructionErrorKindV2::Duplicate,
                inventory,
                repair_target,
                insertion_index,
                exact.len(),
                exact.len() + 1,
                expected_paths
                    .get(insertion_index)
                    .copied()
                    .unwrap_or(STAGE_A_INVENTORY_END_SENTINEL_V1),
                path_of(&duplicate[insertion_index]),
                false,
            );
        }

        for index in 0..exact.len() - 1 {
            let mut reordered = exact.to_vec();
            reordered.swap(index, index + 1);
            assert_inventory_mismatch_v1(
                validate(&reordered),
                ConstructionErrorKindV2::OutOfOrder,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len(),
                expected_paths[index],
                path_of(&reordered[index]),
                false,
            );
        }
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the helper keeps the exact typed membership-mutation contract explicit"
    )]
    fn assert_every_position_typed_identity_mutations_v1<T: Clone>(
        exact: &[T],
        expected_ids: &[&str],
        inventory: &'static str,
        repair_target: &'static str,
        substitution_identity: &'static str,
        extra_identity: &'static str,
        validate: impl Fn(&[T]) -> Result<(), RunnerV2StageAInventoryMismatchV1>,
        identity_of: impl Fn(&T) -> &str,
        set_identity: impl Fn(&mut T, &'static str),
    ) {
        assert_eq!(exact.len(), expected_ids.len());
        assert!(exact.len() >= 2);
        validate(exact).expect("exact typed inventory");

        for index in 0..exact.len() {
            let mut missing = exact.to_vec();
            missing.remove(index);
            assert_inventory_mismatch_v1(
                validate(&missing),
                ConstructionErrorKindV2::Missing,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len() - 1,
                expected_ids[index],
                missing
                    .get(index)
                    .map_or(STAGE_A_INVENTORY_MISSING_SENTINEL_V1, |row| {
                        identity_of(row)
                    }),
                false,
            );

            let mut substituted = exact.to_vec();
            set_identity(&mut substituted[index], substitution_identity);
            assert_inventory_mismatch_v1(
                validate(&substituted),
                ConstructionErrorKindV2::Incompatible,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len(),
                expected_ids[index],
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );
        }

        for insertion_index in 0..=exact.len() {
            let mut extra_row = exact[0].clone();
            set_identity(&mut extra_row, extra_identity);
            let mut extra = exact.to_vec();
            extra.insert(insertion_index, extra_row);
            assert_inventory_mismatch_v1(
                validate(&extra),
                ConstructionErrorKindV2::Unexpected,
                inventory,
                repair_target,
                insertion_index,
                exact.len(),
                exact.len() + 1,
                expected_ids
                    .get(insertion_index)
                    .copied()
                    .unwrap_or(STAGE_A_INVENTORY_END_SENTINEL_V1),
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );

            let duplicate_source_index = if insertion_index == exact.len() {
                exact.len() - 1
            } else {
                (insertion_index + 1) % exact.len()
            };
            let mut duplicate = exact.to_vec();
            duplicate.insert(insertion_index, exact[duplicate_source_index].clone());
            assert_inventory_mismatch_v1(
                validate(&duplicate),
                ConstructionErrorKindV2::Duplicate,
                inventory,
                repair_target,
                insertion_index,
                exact.len(),
                exact.len() + 1,
                expected_ids
                    .get(insertion_index)
                    .copied()
                    .unwrap_or(STAGE_A_INVENTORY_END_SENTINEL_V1),
                identity_of(&duplicate[insertion_index]),
                false,
            );
        }

        for index in 0..exact.len() - 1 {
            let mut reordered = exact.to_vec();
            reordered.swap(index, index + 1);
            assert_inventory_mismatch_v1(
                validate(&reordered),
                ConstructionErrorKindV2::OutOfOrder,
                inventory,
                repair_target,
                index,
                exact.len(),
                exact.len(),
                expected_ids[index],
                identity_of(&reordered[index]),
                false,
            );
        }
    }

    fn assert_semantic_inventory_mismatch_v1(
        result: Result<(), RunnerV2StageAInventoryMismatchV1>,
        inventory: &'static str,
        repair_target: &'static str,
        index0: usize,
        count: usize,
        identity: &str,
        component: &'static str,
        expected_semantic_value: &str,
        observed_safe_value: &str,
    ) {
        let mismatch = result.expect_err("mutated inventory semantics must refuse");
        assert_eq!(mismatch.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(mismatch.inventory(), inventory);
        assert_eq!(mismatch.first_mismatch_index0(), index0);
        assert_eq!(mismatch.expected_ordinal1(), index0 + 1);
        assert_eq!(mismatch.expected_identity(), identity);
        assert_eq!(mismatch.observed_safe_identity(), identity);
        assert!(!mismatch.observed_identity_redacted());
        assert_eq!(mismatch.component(), component);
        assert_eq!(mismatch.expected_semantic_value(), expected_semantic_value);
        assert_eq!(mismatch.observed_safe_value(), observed_safe_value);
        assert!(!mismatch.observed_value_redacted());
        assert_eq!(
            mismatch.semantic_owner(),
            STAGE_A_INVENTORY_SEMANTIC_OWNER_V1
        );
        assert_eq!(mismatch.expected_count(), count);
        assert_eq!(mismatch.observed_count(), count);
        assert_eq!(mismatch.repairs().len(), 1);
        assert_eq!(mismatch.repairs()[0].rank(), 1);
        assert_eq!(
            mismatch.repairs()[0].kind(),
            RepairActionKindV2::ChangeArguments
        );
        assert_eq!(mismatch.repairs()[0].target(), repair_target);
    }

    fn assert_schema_error_v1(
        result: Result<(), ConstructionErrorV2>,
        kind: ConstructionErrorKindV2,
        field: &'static str,
        first_mismatch: usize,
    ) {
        let error = result.expect_err("mutated schema inventory must refuse");
        assert_eq!(error.kind(), kind);
        assert_eq!(error.field(), field);
        assert_eq!(error.observed(), first_mismatch.to_string());
    }

    fn assert_all_canonical_schema_position_mutations_v1(expected: &[&'static str]) {
        let expected = expected
            .iter()
            .map(|name| (*name).to_owned())
            .collect::<Vec<_>>();
        assert_eq!(expected.len(), 42);
        assert_every_position_inventory_mutations_v1(
            &expected,
            |rows| {
                let borrowed = rows.iter().map(String::as_str).collect::<Vec<_>>();
                validate_stage_a_canonical_schema_inventory_v1(&borrowed)
            },
            "runner_v2.base_values.schema_inventory.canonical",
            "restore-exact-canonical-schema-inventory",
        );
    }

    fn oracle_root_for_test_v1(
        cell_id: StableTokenV2,
        expected_outcome: RunnerV2RawOutcomeKindV1,
        expected_reason: RunnerV2RawReasonV1,
        expected_partition: RunnerV2StageAExpectedPartitionV1,
        projection: IndependentOracleProjectionV1,
    ) -> RunnerV2StageAOracleRootV1 {
        build_oracle_row_v1(
            cell_id,
            expected_outcome,
            expected_reason,
            expected_partition,
            projection.numeric,
            projection.diagnostic,
        )
        .expect("bounded test oracle row")
        .root()
    }

    #[test]
    fn boundary_field_meta_and_complete_cell_inventories_are_independent_and_mutation_complete() {
        assert_eq!(
            STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1.len(),
            TEST_BOUNDARIES_EXPECTED_V1.len()
        );
        for (row, expected) in STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1
            .iter()
            .zip(TEST_BOUNDARIES_EXPECTED_V1)
        {
            assert_eq!(row.ordinal, expected.ordinal);
            assert_eq!(row.kind, expected.boundary);
            assert_eq!(row.kind.code(), expected.code);
            assert_eq!(row.stable_name, expected.stable_name);
            assert_eq!(row.kind.stable_name(), expected.stable_name);
        }
        validate_limit_boundary_definitions_diagnostic_v1(&STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1)
            .expect("exact boundary definitions");
        let independent_boundary_ids = TEST_BOUNDARIES_EXPECTED_V1
            .iter()
            .map(|row| row.stable_name.to_owned())
            .collect::<Vec<_>>();
        assert_every_position_inventory_mutations_v1(
            &independent_boundary_ids,
            validate_stage_a_limit_boundary_inventory_v1,
            "runner_v2.base_values.limit_boundaries",
            "restore-exact-limit-boundary-catalog",
        );
        for (index, expected) in TEST_BOUNDARIES_EXPECTED_V1.iter().enumerate() {
            let mut wrong_ordinal = STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1;
            wrong_ordinal[index].ordinal = if expected.ordinal == 12 {
                1
            } else {
                expected.ordinal + 1
            };
            assert_semantic_inventory_mismatch_v1(
                validate_limit_boundary_definitions_diagnostic_v1(&wrong_ordinal),
                "runner_v2.base_values.limit_boundaries",
                "restore-exact-limit-boundary-catalog",
                index,
                12,
                expected.stable_name,
                "ordinal",
                &expected.ordinal.to_string(),
                &wrong_ordinal[index].ordinal.to_string(),
            );

            let mut wrong_kind = STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1;
            wrong_kind[index].kind = TEST_BOUNDARIES_EXPECTED_V1[(index + 1) % 12].boundary;
            assert_semantic_inventory_mismatch_v1(
                validate_limit_boundary_definitions_diagnostic_v1(&wrong_kind),
                "runner_v2.base_values.limit_boundaries",
                "restore-exact-limit-boundary-catalog",
                index,
                12,
                expected.stable_name,
                "kind",
                &format!("{:02}:{}", expected.code, expected.stable_name),
                &format!(
                    "{:02}:{}",
                    wrong_kind[index].kind.code(),
                    wrong_kind[index].kind.stable_name()
                ),
            );

            let mut wrong_name = STAGE_A_LIMIT_BOUNDARY_DEFINITIONS_V1;
            wrong_name[index].stable_name = "runner-v2-unregistered-boundary-name-must-be-redacted";
            assert_inventory_mismatch_v1(
                validate_limit_boundary_definitions_diagnostic_v1(&wrong_name),
                ConstructionErrorKindV2::Incompatible,
                "runner_v2.base_values.limit_boundaries",
                "restore-exact-limit-boundary-catalog",
                index,
                12,
                12,
                expected.stable_name,
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1,
                true,
            );
        }

        assert_eq!(
            STAGE_A_LIMIT_LITERALS_V1.len(),
            TEST_LIMIT_MUTATIONS_EXPECTED_V1.len()
        );
        for (row, expected) in STAGE_A_LIMIT_LITERALS_V1
            .iter()
            .zip(TEST_LIMIT_MUTATIONS_EXPECTED_V1)
        {
            assert_eq!(row.field, expected.field);
            assert_eq!(independent_limit_name_v1(row.field), expected.field_name);
            assert_eq!(row.width, expected.declared_width);
            assert_eq!(row.unit, expected.unit);
        }
        validate_limit_literal_definitions_diagnostic_v1(&STAGE_A_LIMIT_LITERALS_V1)
            .expect("exact limit literal definitions");
        let independent_field_ids = TEST_LIMIT_MUTATIONS_EXPECTED_V1
            .iter()
            .map(|row| row.field_name.to_owned())
            .collect::<Vec<_>>();
        assert_every_position_inventory_mutations_v1(
            &independent_field_ids,
            validate_stage_a_limit_field_inventory_v1,
            "runner_v2.base_values.limit_fields",
            "restore-exact-limit-field-and-width-catalog",
        );
        let changed_same_width = |value| match value {
            RunnerLimitValueV2::U32(0) => RunnerLimitValueV2::U32(1),
            RunnerLimitValueV2::U32(_) => RunnerLimitValueV2::U32(0),
            RunnerLimitValueV2::U64(0) => RunnerLimitValueV2::U64(1),
            RunnerLimitValueV2::U64(_) => RunnerLimitValueV2::U64(0),
        };
        for (index, expected) in TEST_LIMIT_MUTATIONS_EXPECTED_V1.iter().enumerate() {
            let exact = STAGE_A_LIMIT_LITERALS_V1;

            let mut wrong_field = exact;
            wrong_field[index].field = TEST_LIMIT_MUTATIONS_EXPECTED_V1[(index + 1) % 71].field;
            let mismatch = validate_limit_literal_definitions_diagnostic_v1(&wrong_field)
                .expect_err("field substitution must refuse");
            assert_eq!(mismatch.first_mismatch_index0(), index);
            assert_eq!(mismatch.expected_identity(), expected.field_name);
            assert_eq!(
                mismatch.observed_safe_identity(),
                TEST_LIMIT_MUTATIONS_EXPECTED_V1[(index + 1) % 71].field_name
            );
            assert_eq!(mismatch.component(), "identity");
            assert!(!mismatch.observed_identity_redacted());

            let mut wrong_width = exact;
            wrong_width[index].width = match expected.declared_width {
                RunnerLimitWidthV2::U32 => RunnerLimitWidthV2::U64,
                RunnerLimitWidthV2::U64 => RunnerLimitWidthV2::U32,
            };
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_width),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "width",
                limit_width_name_v1(expected.declared_width),
                limit_width_name_v1(wrong_width[index].width),
            );

            let mut wrong_unit = exact;
            wrong_unit[index].unit = if expected.unit == RunnerLimitUnitV2::Count {
                RunnerLimitUnitV2::Records
            } else {
                RunnerLimitUnitV2::Count
            };
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_unit),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "unit",
                independent_limit_unit_name_v1(expected.unit),
                independent_limit_unit_name_v1(wrong_unit[index].unit),
            );

            let mut wrong_tightenability = exact;
            wrong_tightenability[index].tightenability = match exact[index].tightenability {
                RunnerLimitTightenabilityV2::Tightenable => RunnerLimitTightenabilityV2::Fixed,
                RunnerLimitTightenabilityV2::Fixed => RunnerLimitTightenabilityV2::Tightenable,
            };
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_tightenability),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "tightenability",
                &format!("{:?}", exact[index].tightenability),
                &format!("{:?}", wrong_tightenability[index].tightenability),
            );

            let mut wrong_minimum = exact;
            wrong_minimum[index].minimum_rule =
                if matches!(exact[index].minimum_rule, RunnerLimitMinimumRuleV2::Fixed) {
                    RunnerLimitMinimumRuleV2::ZeroAllowed
                } else {
                    RunnerLimitMinimumRuleV2::Fixed
                };
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_minimum),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "minimum-rule",
                &format!("{:?}", exact[index].minimum_rule),
                &format!("{:?}", wrong_minimum[index].minimum_rule),
            );

            let mut wrong_smoke = exact;
            wrong_smoke[index].smoke = changed_same_width(exact[index].smoke);
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_smoke),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "smoke-ceiling",
                &limit_value_safe_name_v1(exact[index].smoke),
                &limit_value_safe_name_v1(wrong_smoke[index].smoke),
            );

            let mut wrong_full = exact;
            wrong_full[index].full = changed_same_width(exact[index].full);
            assert_semantic_inventory_mismatch_v1(
                validate_limit_literal_definitions_diagnostic_v1(&wrong_full),
                "runner_v2.base_values.limit_fields",
                "restore-exact-limit-field-and-width-catalog",
                index,
                71,
                expected.field_name,
                "full-ceiling",
                &limit_value_safe_name_v1(exact[index].full),
                &limit_value_safe_name_v1(wrong_full[index].full),
            );
        }

        assert_eq!(META_CELL_DEFINITIONS_V1.len(), TEST_META_EXPECTED_V1.len());
        for (row, expected) in META_CELL_DEFINITIONS_V1.iter().zip(TEST_META_EXPECTED_V1) {
            assert_eq!(row.ordinal, expected.ordinal);
            assert_eq!(row.id_suffix, expected.id_suffix);
            assert_eq!(row.group, expected.group);
            assert_eq!(row.operation, expected.operation);
            assert_eq!(row.expected_outcome, expected.outcome);
            assert_eq!(row.expected_reason, expected.reason);
            assert_eq!(row.expected_partition, expected.partition);
        }
        validate_meta_definitions_diagnostic_v1(&META_CELL_DEFINITIONS_V1)
            .expect("exact meta definitions");
        let independent_meta_ids = TEST_META_EXPECTED_V1
            .iter()
            .map(|row| row.stable_id.to_owned())
            .collect::<Vec<_>>();
        assert_every_position_inventory_mutations_v1(
            &independent_meta_ids,
            validate_stage_a_meta_inventory_v1,
            "runner_v2.base_values.meta_cells",
            "restore-exact-meta-cell-catalog",
        );
        for (index, expected) in TEST_META_EXPECTED_V1.iter().enumerate() {
            let exact = META_CELL_DEFINITIONS_V1;

            let mut wrong_ordinal = exact;
            wrong_ordinal[index].ordinal = if expected.ordinal == 15 {
                1
            } else {
                expected.ordinal + 1
            };
            let mismatch = validate_meta_definitions_diagnostic_v1(&wrong_ordinal)
                .expect_err("meta ordinal substitution must refuse");
            assert_eq!(mismatch.first_mismatch_index0(), index);
            assert_eq!(mismatch.expected_identity(), expected.stable_id);
            assert_eq!(mismatch.component(), "identity");

            let mut wrong_suffix = exact;
            wrong_suffix[index].id_suffix = "runner-v2-unregistered-meta-suffix-must-be-redacted";
            let mismatch = validate_meta_definitions_diagnostic_v1(&wrong_suffix)
                .expect_err("meta suffix substitution must refuse");
            assert_eq!(mismatch.first_mismatch_index0(), index);
            assert_eq!(mismatch.expected_identity(), expected.stable_id);
            assert_eq!(
                mismatch.observed_safe_identity(),
                STAGE_A_INVENTORY_REDACTED_SENTINEL_V1
            );
            assert!(mismatch.observed_identity_redacted());

            let mut wrong_group = exact;
            wrong_group[index].group =
                if expected.group == RunnerV2StageACellGroupV1::LiteralUnitBoundary {
                    RunnerV2StageACellGroupV1::StateModel
                } else {
                    RunnerV2StageACellGroupV1::LiteralUnitBoundary
                };
            assert_semantic_inventory_mismatch_v1(
                validate_meta_definitions_diagnostic_v1(&wrong_group),
                "runner_v2.base_values.meta_cells",
                "restore-exact-meta-cell-catalog",
                index,
                15,
                expected.stable_id,
                "group",
                &expected.group.code().to_string(),
                &wrong_group[index].group.code().to_string(),
            );

            let mut wrong_operation = exact;
            wrong_operation[index].operation = TEST_META_EXPECTED_V1[(index + 1) % 15].operation;
            assert_semantic_inventory_mismatch_v1(
                validate_meta_definitions_diagnostic_v1(&wrong_operation),
                "runner_v2.base_values.meta_cells",
                "restore-exact-meta-cell-catalog",
                index,
                15,
                expected.stable_id,
                "operation",
                &format!("{:?}", expected.operation),
                &format!("{:?}", wrong_operation[index].operation),
            );

            let mut wrong_outcome = exact;
            wrong_outcome[index].expected_outcome =
                if expected.outcome == RunnerV2RawOutcomeKindV1::Accepted {
                    RunnerV2RawOutcomeKindV1::Refused
                } else {
                    RunnerV2RawOutcomeKindV1::Accepted
                };
            assert_semantic_inventory_mismatch_v1(
                validate_meta_definitions_diagnostic_v1(&wrong_outcome),
                "runner_v2.base_values.meta_cells",
                "restore-exact-meta-cell-catalog",
                index,
                15,
                expected.stable_id,
                "expected-outcome",
                &format!("{:?}", expected.outcome),
                &format!("{:?}", wrong_outcome[index].expected_outcome),
            );

            let mut wrong_reason = exact;
            wrong_reason[index].expected_reason = RunnerV2RawReasonV1::UnknownClosedValue;
            assert_semantic_inventory_mismatch_v1(
                validate_meta_definitions_diagnostic_v1(&wrong_reason),
                "runner_v2.base_values.meta_cells",
                "restore-exact-meta-cell-catalog",
                index,
                15,
                expected.stable_id,
                "expected-reason",
                &format!(
                    "{:02}:{}",
                    expected.reason.code(),
                    runner_v2_raw_reason_name_v1(expected.reason)
                ),
                "08:unknown-closed-value",
            );

            let mut wrong_partition = exact;
            wrong_partition[index].expected_partition =
                RunnerV2StageAExpectedPartitionV1::Unsupported;
            assert_semantic_inventory_mismatch_v1(
                validate_meta_definitions_diagnostic_v1(&wrong_partition),
                "runner_v2.base_values.meta_cells",
                "restore-exact-meta-cell-catalog",
                index,
                15,
                expected.stable_id,
                "expected-partition",
                &expected.partition.code().to_string(),
                &RunnerV2StageAExpectedPartitionV1::Unsupported
                    .code()
                    .to_string(),
            );
        }

        let independent_limit_cell_ids = TEST_LIMIT_MUTATIONS_EXPECTED_V1
            .iter()
            .flat_map(|field| {
                TEST_BOUNDARIES_EXPECTED_V1.iter().map(move |boundary| {
                    format!(
                        "runner-v2.base-values.limit-{:03}.boundary-{:02}-{}.v1",
                        field.ordinal, boundary.code, boundary.stable_name
                    )
                })
            })
            .collect::<Vec<_>>();
        assert_eq!(independent_limit_cell_ids.len(), 852);
        assert_every_position_inventory_mutations_v1(
            &independent_limit_cell_ids,
            validate_stage_a_limit_cell_inventory_v1,
            "runner_v2.base_values.limit_cells",
            "restore-exact-71-by-12-limit-cell-inventory",
        );

        let mut independent_complete_cell_ids = independent_limit_cell_ids;
        independent_complete_cell_ids.extend(
            TEST_META_EXPECTED_V1
                .iter()
                .map(|row| row.stable_id.to_owned()),
        );
        assert_eq!(independent_complete_cell_ids.len(), 867);
        assert_every_position_inventory_mutations_v1(
            &independent_complete_cell_ids,
            validate_stage_a_complete_cell_inventory_v1,
            "runner_v2.base_values.complete_cells",
            "restore-exact-852-plus-15-cell-inventory",
        );
    }

    #[test]
    fn limit_literal_table_exact_matches_all_descriptors_and_both_profiles() {
        assert_eq!(STAGE_A_LIMIT_LITERALS_V1.len(), RUNNER_LIMIT_FIELD_COUNT_V2);
        let mut smoke_full_differences = Vec::new();
        let mut tightenable = 0;
        let mut fixed = 0;

        for (index, (field, literal)) in RunnerLimitFieldV2::ALL
            .iter()
            .copied()
            .zip(STAGE_A_LIMIT_LITERALS_V1)
            .enumerate()
        {
            let descriptor = field.descriptor();
            assert_eq!(literal.field, field);
            assert_eq!(descriptor.ordinal, u16::try_from(index + 1).unwrap());
            assert_eq!(descriptor.field, field);
            assert_eq!(STAGE_A_INDEPENDENT_LIMIT_NAMES_V1[index], descriptor.name);
            assert_eq!(literal.width, descriptor.width);
            assert_eq!(literal.unit, descriptor.unit);
            assert_eq!(literal.tightenability, descriptor.tightenability);
            assert_eq!(literal.minimum_rule, descriptor.minimum_rule);
            assert_eq!(
                literal.smoke,
                RunnerLimitsV2::base(RunProfileV2::Smoke).value(field)
            );
            assert_eq!(
                literal.full,
                RunnerLimitsV2::base(RunProfileV2::Full).value(field)
            );
            if literal.smoke != literal.full {
                smoke_full_differences.push(field.ordinal());
            }
            match literal.tightenability {
                RunnerLimitTightenabilityV2::Tightenable => tightenable += 1,
                RunnerLimitTightenabilityV2::Fixed => fixed += 1,
            }
        }

        assert_eq!(smoke_full_differences, vec![13, 31, 32, 33, 35]);
        assert_eq!(tightenable, 67);
        assert_eq!(fixed, 4);
    }

    #[test]
    fn declaration_is_exact_867_cell_cartesian_product_with_distinct_manifests() {
        let declaration = declare_24_1_1_1_1_v1().expect("complete declaration");
        let repeated = declare_24_1_1_1_1_v1().expect("repeat declaration");
        assert_eq!(declaration, repeated);
        assert_eq!(
            declaration.package_id().as_str(),
            RUNNER_V2_BASE_VALUES_PACKAGE_ID_V1
        );
        assert_eq!(declaration.cells().len(), 867);
        assert_eq!(declaration.oracles().len(), 867);
        assert_eq!(declaration.projections().len(), 867);

        let mut cell_ids = BTreeSet::new();
        let mut oracle_roots = BTreeSet::new();
        let mut case_manifest_roots = BTreeSet::new();
        for (index, ((cell, oracle), projection)) in declaration
            .cells()
            .iter()
            .zip(declaration.oracles())
            .zip(declaration.projections())
            .enumerate()
        {
            assert_eq!(cell.ordinal(), u16::try_from(index + 1).unwrap());
            assert_eq!(projection.ordinal(), cell.ordinal());
            assert_eq!(cell.cell_id(), oracle.cell_id());
            assert_eq!(cell.cell_id(), projection.cell_id());
            assert_eq!(cell.oracle_root(), oracle.root());
            assert_eq!(cell.case_manifest_root(), projection.case_manifest_root());
            assert_ne!(
                cell.oracle_root().bytes(),
                cell.case_manifest_root().bytes()
            );
            assert!(cell_ids.insert(cell.cell_id().as_str()));
            assert!(oracle_roots.insert(*cell.oracle_root().bytes()));
            assert!(case_manifest_roots.insert(*cell.case_manifest_root().bytes()));
            assert_eq!(
                projection.posix_script().as_str(),
                "scripts/ci/runner_v2_base_work_packages_e2e.sh"
            );
            assert_eq!(
                projection.windows_script().as_str(),
                "scripts/ci/runner_v2_base_work_packages_e2e.ps1"
            );
            assert_eq!(projection.expected_partition(), oracle.expected_partition());
        }
        assert_eq!(cell_ids.len(), 867);
        assert_eq!(oracle_roots.len(), 867);
        assert_eq!(case_manifest_roots.len(), 867);

        for (field_index, field) in RunnerLimitFieldV2::ALL.iter().copied().enumerate() {
            let literal = STAGE_A_LIMIT_LITERALS_V1[field_index];
            for (boundary_index, boundary) in
                RunnerV2LimitBoundaryKindV1::ALL.iter().copied().enumerate()
            {
                let index = field_index * RUNNER_V2_LIMIT_BOUNDARY_KIND_COUNT_V1 + boundary_index;
                let cell = &declaration.cells()[index];
                assert_eq!(
                    cell.cell_id(),
                    &limit_cell_id_v1(field, boundary).expect("stable limit cell ID")
                );
                assert_eq!(cell.group(), RunnerV2StageACellGroupV1::LiteralUnitBoundary);
                match cell.operation() {
                    RunnerV2StageACellOperationV1::Limit {
                        field: actual_field,
                        boundary: actual_boundary,
                        value,
                    } => {
                        assert_eq!(actual_field, field);
                        assert_eq!(actual_boundary, boundary);
                        assert_eq!(value, independent_boundary_value_v1(literal, boundary));
                    }
                    RunnerV2StageACellOperationV1::Meta(_) => {
                        panic!("limit Cartesian row became a meta operation")
                    }
                }
            }
        }

        assert_eq!(
            declaration.cells()[RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1].ordinal(),
            853
        );
        assert!(
            declaration.cells()[RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1..]
                .iter()
                .all(|cell| matches!(cell.operation(), RunnerV2StageACellOperationV1::Meta(_)))
        );
    }

    #[test]
    fn oracle_root_binds_every_declared_header_numeric_and_diagnostic_component() {
        let literal = STAGE_A_LIMIT_LITERALS_V1[0];
        let boundary = RunnerV2LimitBoundaryKindV1::CheckedRepresentationalOverflowRefusal;
        let cell_id = limit_cell_id_v1(literal.field, boundary).expect("stable cell ID");
        let projection =
            independent_limit_oracle_projection_v1(literal, boundary).expect("refusal oracle");
        assert!(projection.numeric.len() >= 2);
        assert_eq!(
            projection
                .diagnostic
                .as_ref()
                .expect("refusal diagnostic")
                .repairs
                .len(),
            1
        );

        let base_root = oracle_root_for_test_v1(
            cell_id.clone(),
            projection.outcome,
            projection.reason,
            RunnerV2StageAExpectedPartitionV1::ExpectedRefusal,
            projection.clone(),
        );
        let root_for_projection = |mutated: IndependentOracleProjectionV1| {
            oracle_root_for_test_v1(
                cell_id.clone(),
                projection.outcome,
                projection.reason,
                RunnerV2StageAExpectedPartitionV1::ExpectedRefusal,
                mutated,
            )
        };

        assert_ne!(
            base_root,
            oracle_root_for_test_v1(
                stage_a_token(
                    "runner_v2.base_values.test.oracle.cell_id",
                    "runner-v2-oracle-root-sensitivity-other-cell.v1",
                )
                .expect("alternate cell ID"),
                projection.outcome,
                projection.reason,
                RunnerV2StageAExpectedPartitionV1::ExpectedRefusal,
                projection.clone(),
            ),
            "cell identity must affect the oracle root"
        );
        assert_ne!(
            base_root,
            oracle_root_for_test_v1(
                cell_id.clone(),
                RunnerV2RawOutcomeKindV1::Failed,
                projection.reason,
                RunnerV2StageAExpectedPartitionV1::ExpectedRefusal,
                projection.clone(),
            ),
            "outcome must affect the oracle root"
        );
        assert_ne!(
            base_root,
            oracle_root_for_test_v1(
                cell_id.clone(),
                projection.outcome,
                RunnerV2RawReasonV1::InternalInvariantFailure,
                RunnerV2StageAExpectedPartitionV1::ExpectedRefusal,
                projection.clone(),
            ),
            "reason must affect the oracle root"
        );
        assert_ne!(
            base_root,
            oracle_root_for_test_v1(
                cell_id.clone(),
                projection.outcome,
                projection.reason,
                RunnerV2StageAExpectedPartitionV1::ExpectedFailure,
                projection.clone(),
            ),
            "partition must affect the oracle root"
        );

        let mut numeric_count = projection.clone();
        numeric_count.numeric.pop();
        assert_ne!(
            base_root,
            root_for_projection(numeric_count),
            "numeric cardinality must affect the oracle root"
        );

        let mut numeric_order = projection.clone();
        numeric_order.numeric.swap(0, 1);
        assert_ne!(
            base_root,
            root_for_projection(numeric_order),
            "numeric source order must affect the oracle root"
        );

        let mut numeric_name = projection.clone();
        numeric_name.numeric[0].name = stage_a_token(
            "runner_v2.base_values.test.oracle.numeric.name",
            "mutated-field-ordinal",
        )
        .expect("alternate numeric name");
        assert_ne!(
            base_root,
            root_for_projection(numeric_name),
            "numeric name must affect the oracle root"
        );

        let mut numeric_value = projection.clone();
        numeric_value.numeric[0].value = RunnerV2StageAOracleNumericValueV1::Count(9_999);
        assert_ne!(
            base_root,
            root_for_projection(numeric_value),
            "numeric value must affect the oracle root"
        );

        let mut numeric_unit = projection.clone();
        numeric_unit.numeric[0].unit =
            RunnerV2StageAOracleNumericUnitV1::Limit(RunnerLimitUnitV2::Records);
        assert_ne!(
            base_root,
            root_for_projection(numeric_unit),
            "numeric unit must affect the oracle root"
        );

        let mut diagnostic_presence = projection.clone();
        diagnostic_presence.diagnostic = None;
        assert_ne!(
            base_root,
            root_for_projection(diagnostic_presence),
            "diagnostic presence must affect the oracle root"
        );

        let mut diagnostic_code = projection.clone();
        diagnostic_code
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .code = DiagnosticCodeV2::RunnerBlocked;
        assert_ne!(
            base_root,
            root_for_projection(diagnostic_code),
            "diagnostic code must affect the oracle root"
        );

        let mut diagnostic_owner = projection.clone();
        diagnostic_owner
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .owner = stage_a_token(
            "runner_v2.base_values.test.oracle.diagnostic.owner",
            "fs-evidence-runner.runner-v2.root-sensitivity",
        )
        .expect("alternate diagnostic owner");
        assert_ne!(
            base_root,
            root_for_projection(diagnostic_owner),
            "diagnostic owner must affect the oracle root"
        );

        let mut diagnostic_retryability = projection.clone();
        diagnostic_retryability
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .retryability = RetryabilityV2::AfterEnvironmentChange;
        assert_ne!(
            base_root,
            root_for_projection(diagnostic_retryability),
            "diagnostic retryability must affect the oracle root"
        );

        let mut prerequisite_first = projection.clone();
        prerequisite_first
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .prerequisites = vec![
            stage_a_token(
                "runner_v2.base_values.test.oracle.diagnostic.prerequisite",
                "first-root-sensitivity-prerequisite",
            )
            .expect("first prerequisite"),
        ]
        .into_boxed_slice();
        let prerequisite_first_root = root_for_projection(prerequisite_first.clone());
        assert_ne!(
            base_root, prerequisite_first_root,
            "diagnostic prerequisite cardinality must affect the oracle root"
        );
        prerequisite_first
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .prerequisites[0] = stage_a_token(
            "runner_v2.base_values.test.oracle.diagnostic.prerequisite",
            "second-root-sensitivity-prerequisite",
        )
        .expect("second prerequisite");
        assert_ne!(
            prerequisite_first_root,
            root_for_projection(prerequisite_first),
            "diagnostic prerequisite identity must affect the oracle root"
        );

        let mut repair_count = projection.clone();
        repair_count
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .repairs = Vec::new().into_boxed_slice();
        assert_ne!(
            base_root,
            root_for_projection(repair_count),
            "repair cardinality must affect the oracle root"
        );

        let mut repair_rank = projection.clone();
        repair_rank.diagnostic.as_mut().expect("diagnostic").repairs[0].rank = 2;
        assert_ne!(
            base_root,
            root_for_projection(repair_rank),
            "repair rank must affect the oracle root"
        );

        let mut repair_kind = projection.clone();
        repair_kind.diagnostic.as_mut().expect("diagnostic").repairs[0].kind =
            RepairActionKindV2::SupplyEvidence;
        assert_ne!(
            base_root,
            root_for_projection(repair_kind),
            "repair kind must affect the oracle root"
        );

        let mut repair_target = projection.clone();
        repair_target
            .diagnostic
            .as_mut()
            .expect("diagnostic")
            .repairs[0]
            .target = stage_a_token(
            "runner_v2.base_values.test.oracle.diagnostic.repair_target",
            "mutated-root-sensitivity-target",
        )
        .expect("alternate repair target");
        assert_ne!(
            base_root,
            root_for_projection(repair_target),
            "repair target must affect the oracle root"
        );
    }

    #[test]
    fn meta_oracle_fixture_is_independent_and_literal_exact() {
        let expected = [
            (
                "typed-absence-distinct-from-zero",
                RunnerV2StageACellGroupV1::LiteralUnitBoundary,
                RunnerV2StageAMetaOperationV1::TypedAbsenceDistinctFromZero,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "f32-named-total-order",
                RunnerV2StageACellGroupV1::PropertyMetamorphic,
                RunnerV2StageAMetaOperationV1::F32NamedTotalOrder,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "f64-named-total-order",
                RunnerV2StageACellGroupV1::PropertyMetamorphic,
                RunnerV2StageAMetaOperationV1::F64NamedTotalOrder,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "capability-none-contract",
                RunnerV2StageACellGroupV1::StateModel,
                RunnerV2StageAMetaOperationV1::CapabilityNoneContract,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "common-requirements-exact",
                RunnerV2StageACellGroupV1::StateModel,
                RunnerV2StageAMetaOperationV1::CommonRequirementsExact,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "common-requirement-reordered-refusal",
                RunnerV2StageACellGroupV1::MutationFuzz,
                RunnerV2StageAMetaOperationV1::CommonRequirementReorderedRefusal,
                RunnerV2RawOutcomeKindV1::Refused,
                RunnerV2RawReasonV1::ExactMembershipMismatch,
                RunnerV2StageAExpectedPartitionV1::Mutation,
            ),
            (
                "future-sources-exact",
                RunnerV2StageACellGroupV1::SourceClosure,
                RunnerV2StageAMetaOperationV1::FutureSourcesExact,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "rootless-ac58",
                RunnerV2StageACellGroupV1::ApiCompileFail,
                RunnerV2StageAMetaOperationV1::RootlessAc58,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "owner-source-fragment",
                RunnerV2StageACellGroupV1::SourceClosure,
                RunnerV2StageAMetaOperationV1::OwnerSourceFragment,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "local-route",
                RunnerV2StageACellGroupV1::NoMockLocalIntegration,
                RunnerV2StageAMetaOperationV1::LocalRoute,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "diagnostic-redaction",
                RunnerV2StageACellGroupV1::DetailedLoggingRedaction,
                RunnerV2StageAMetaOperationV1::DiagnosticRedaction,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "reproduction-declaration",
                RunnerV2StageACellGroupV1::ReproductionDeclaration,
                RunnerV2StageAMetaOperationV1::ReproductionDeclaration,
                RunnerV2RawOutcomeKindV1::Accepted,
                RunnerV2RawReasonV1::ExactCheckedValue,
                RunnerV2StageAExpectedPartitionV1::EligiblePositive,
            ),
            (
                "compile-fail-ordering-surface",
                RunnerV2StageACellGroupV1::ApiCompileFail,
                RunnerV2StageAMetaOperationV1::CompileFailOrderingSurface,
                RunnerV2RawOutcomeKindV1::Inapplicable,
                RunnerV2RawReasonV1::PureDeclarationFacet,
                RunnerV2StageAExpectedPartitionV1::Inapplicable,
            ),
            (
                "shard-inapplicable",
                RunnerV2StageACellGroupV1::FaultResourceCancellation,
                RunnerV2StageAMetaOperationV1::ShardInapplicable,
                RunnerV2RawOutcomeKindV1::Inapplicable,
                RunnerV2RawReasonV1::ShardInapplicable,
                RunnerV2StageAExpectedPartitionV1::Inapplicable,
            ),
            (
                "resume-inapplicable",
                RunnerV2StageACellGroupV1::FaultResourceCancellation,
                RunnerV2StageAMetaOperationV1::ResumeInapplicable,
                RunnerV2RawOutcomeKindV1::Inapplicable,
                RunnerV2RawReasonV1::ResumeInapplicable,
                RunnerV2StageAExpectedPartitionV1::Inapplicable,
            ),
        ];
        let expected_ids = [
            "runner-v2.base-values.meta-001-typed-absence-distinct-from-zero.v1",
            "runner-v2.base-values.meta-002-f32-named-total-order.v1",
            "runner-v2.base-values.meta-003-f64-named-total-order.v1",
            "runner-v2.base-values.meta-004-capability-none-contract.v1",
            "runner-v2.base-values.meta-005-common-requirements-exact.v1",
            "runner-v2.base-values.meta-006-common-requirement-reordered-refusal.v1",
            "runner-v2.base-values.meta-007-future-sources-exact.v1",
            "runner-v2.base-values.meta-008-rootless-ac58.v1",
            "runner-v2.base-values.meta-009-owner-source-fragment.v1",
            "runner-v2.base-values.meta-010-local-route.v1",
            "runner-v2.base-values.meta-011-diagnostic-redaction.v1",
            "runner-v2.base-values.meta-012-reproduction-declaration.v1",
            "runner-v2.base-values.meta-013-compile-fail-ordering-surface.v1",
            "runner-v2.base-values.meta-014-shard-inapplicable.v1",
            "runner-v2.base-values.meta-015-resume-inapplicable.v1",
        ];

        assert_eq!(expected.len(), RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1);
        for (index, (definition, expected)) in
            META_CELL_DEFINITIONS_V1.iter().zip(expected).enumerate()
        {
            assert_eq!(definition.id_suffix, expected.0);
            assert_eq!(definition.group, expected.1);
            assert_eq!(definition.operation, expected.2);
            assert_eq!(definition.expected_outcome, expected.3);
            assert_eq!(definition.expected_reason, expected.4);
            assert_eq!(definition.expected_partition, expected.5);
            assert_eq!(
                meta_cell_id_v1(index, definition.id_suffix)
                    .unwrap()
                    .as_str(),
                expected_ids[index]
            );
        }

        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        let handoff = evaluate_24_1_1_1_1_cell_v1().expect("fresh complete handoff");
        let meta_start = RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1;
        let declaration_ids = declaration.cells()[meta_start..]
            .iter()
            .map(|cell| cell.cell_id().as_str())
            .collect::<Vec<_>>();
        let oracle_ids = declaration.oracles()[meta_start..]
            .iter()
            .map(|oracle| oracle.cell_id().as_str())
            .collect::<Vec<_>>();
        let projection_ids = declaration.projections()[meta_start..]
            .iter()
            .map(|projection| projection.cell_id().as_str())
            .collect::<Vec<_>>();
        let handoff_ids = handoff.cells()[meta_start..]
            .iter()
            .map(|cell| cell.cell_id().as_str())
            .collect::<Vec<_>>();
        assert_eq!(declaration_ids.as_slice(), expected_ids.as_slice());
        assert_eq!(oracle_ids.as_slice(), expected_ids.as_slice());
        assert_eq!(projection_ids.as_slice(), expected_ids.as_slice());
        assert_eq!(handoff_ids.as_slice(), expected_ids.as_slice());
        assert_eq!(
            expected_ids.iter().copied().collect::<BTreeSet<_>>().len(),
            RUNNER_V2_BASE_VALUES_META_CELL_COUNT_V1
        );
    }

    #[test]
    fn fresh_evaluator_exact_joins_all_cells_and_matches_independent_oracles() {
        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        let first = evaluate_24_1_1_1_1_cell_v1().expect("first full evaluation");
        let second = evaluate_24_1_1_1_1_cell_v1().expect("second full evaluation");
        assert_eq!(first, second);
        assert_eq!(first.package_id(), declaration.package_id());
        assert_eq!(first.cells().len(), RUNNER_V2_BASE_VALUES_CELL_COUNT_V1);
        first
            .validate_exact_cell_order(
                &declaration
                    .cells()
                    .iter()
                    .map(|cell| cell.cell_id().clone())
                    .collect::<Vec<_>>(),
            )
            .expect("complete exact join");

        let mut accepted = 0;
        let mut refused = 0;
        let mut inapplicable = 0;
        let mut numeric_total = 0;
        let mut diagnostic_total = 0;
        let mut repair_total = 0;
        let mut prerequisite_total = 0;
        let mut limit_accepted = 0;
        let mut limit_refused = 0;
        let mut limit_inapplicable = 0;
        let mut limit_numeric_total = 0;
        let mut limit_diagnostic_total = 0;
        let mut limit_repair_total = 0;
        let mut limit_prerequisite_total = 0;
        let mut above_ceiling = 0;
        let mut fixed_changed = 0;
        let mut below_minimum = 0;
        let mut joint_feasibility = 0;
        let mut checked_overflow = 0;
        for (cell_index, ((actual, oracle), declaration_cell)) in first
            .cells()
            .iter()
            .zip(declaration.oracles())
            .zip(declaration.cells())
            .enumerate()
        {
            assert_eq!(actual.cell_id(), declaration_cell.cell_id());
            assert_eq!(actual.outcome(), oracle.expected_outcome());
            assert_eq!(actual.reason(), oracle.expected_reason());
            assert_eq!(
                actual.numeric().len(),
                oracle.expected_numeric().len(),
                "numeric count mismatch for {}",
                actual.cell_id().as_str()
            );
            for (numeric_index, (actual_numeric, expected_numeric)) in actual
                .numeric()
                .iter()
                .zip(oracle.expected_numeric())
                .enumerate()
            {
                assert_eq!(
                    actual_numeric.name(),
                    expected_numeric.name(),
                    "numeric name mismatch for {} at index {numeric_index}",
                    actual.cell_id().as_str()
                );
                match expected_numeric.value() {
                    RunnerV2StageAOracleNumericValueV1::Limit(expected) => assert_eq!(
                        actual_numeric.value(),
                        &RunnerV2SafeNumericValueV1::Limit(expected),
                        "numeric value mismatch for {} at index {numeric_index}",
                        actual.cell_id().as_str()
                    ),
                    RunnerV2StageAOracleNumericValueV1::Count(expected) => assert_eq!(
                        actual_numeric.value(),
                        &RunnerV2SafeNumericValueV1::Count(expected),
                        "numeric count mismatch for {} at index {numeric_index}",
                        actual.cell_id().as_str()
                    ),
                }
                match expected_numeric.unit() {
                    RunnerV2StageAOracleNumericUnitV1::Limit(expected) => assert_eq!(
                        actual_numeric.unit(),
                        RunnerV2SafeNumericUnitV1::Limit(expected),
                        "numeric unit mismatch for {} at index {numeric_index}",
                        actual.cell_id().as_str()
                    ),
                    RunnerV2StageAOracleNumericUnitV1::LogicalCount => assert_eq!(
                        actual_numeric.unit(),
                        RunnerV2SafeNumericUnitV1::Logical(LogicalUnitV2::Count),
                        "numeric count unit mismatch for {} at index {numeric_index}",
                        actual.cell_id().as_str()
                    ),
                }
            }
            assert!(
                actual
                    .numeric()
                    .windows(2)
                    .all(|pair| pair[0].name() < pair[1].name())
            );
            match (actual.diagnostic(), oracle.expected_diagnostic()) {
                (None, None) => {}
                (Some(actual_diagnostic), Some(expected_diagnostic)) => {
                    assert_eq!(
                        actual_diagnostic.code(),
                        expected_diagnostic.code(),
                        "diagnostic code mismatch for {}",
                        actual.cell_id().as_str()
                    );
                    assert_eq!(
                        actual_diagnostic.owner(),
                        expected_diagnostic.owner(),
                        "diagnostic owner mismatch for {}",
                        actual.cell_id().as_str()
                    );
                    assert_eq!(
                        actual_diagnostic.retryability(),
                        expected_diagnostic.retryability(),
                        "diagnostic retryability mismatch for {}",
                        actual.cell_id().as_str()
                    );
                    assert_eq!(
                        actual_diagnostic.prerequisites(),
                        expected_diagnostic.prerequisites(),
                        "diagnostic prerequisite mismatch for {}",
                        actual.cell_id().as_str()
                    );
                    assert_eq!(
                        actual_diagnostic.repairs().len(),
                        expected_diagnostic.repairs().len(),
                        "diagnostic repair count mismatch for {}",
                        actual.cell_id().as_str()
                    );
                    for (repair_index, (actual_repair, expected_repair)) in actual_diagnostic
                        .repairs()
                        .iter()
                        .zip(expected_diagnostic.repairs())
                        .enumerate()
                    {
                        assert_eq!(
                            actual_repair.rank(),
                            expected_repair.rank(),
                            "repair rank mismatch for {} at index {repair_index}",
                            actual.cell_id().as_str()
                        );
                        assert_eq!(
                            actual_repair.kind(),
                            expected_repair.kind(),
                            "repair kind mismatch for {} at index {repair_index}",
                            actual.cell_id().as_str()
                        );
                        assert_eq!(
                            actual_repair.target(),
                            expected_repair.target(),
                            "repair target mismatch for {} at index {repair_index}",
                            actual.cell_id().as_str()
                        );
                    }
                }
                _ => panic!(
                    "diagnostic presence mismatch for {}",
                    actual.cell_id().as_str()
                ),
            }
            numeric_total += actual.numeric().len();
            if let Some(diagnostic) = actual.diagnostic() {
                diagnostic_total += 1;
                repair_total += diagnostic.repairs().len();
                prerequisite_total += diagnostic.prerequisites().len();
            }
            if cell_index < RUNNER_V2_BASE_VALUES_LIMIT_CELL_COUNT_V1 {
                limit_numeric_total += actual.numeric().len();
                if let Some(diagnostic) = actual.diagnostic() {
                    limit_diagnostic_total += 1;
                    limit_repair_total += diagnostic.repairs().len();
                    limit_prerequisite_total += diagnostic.prerequisites().len();
                }
                match actual.outcome() {
                    RunnerV2RawOutcomeKindV1::Accepted => limit_accepted += 1,
                    RunnerV2RawOutcomeKindV1::Refused => {
                        limit_refused += 1;
                        match actual.reason() {
                            RunnerV2RawReasonV1::AboveProfileCeiling => above_ceiling += 1,
                            RunnerV2RawReasonV1::FixedRepresentationChanged => fixed_changed += 1,
                            RunnerV2RawReasonV1::BelowStructuralMinimum => below_minimum += 1,
                            RunnerV2RawReasonV1::JointFeasibilityViolation => {
                                joint_feasibility += 1
                            }
                            RunnerV2RawReasonV1::CheckedRepresentationalOverflow => {
                                checked_overflow += 1
                            }
                            other => panic!(
                                "unexpected limit refusal reason {other:?} for {}",
                                actual.cell_id().as_str()
                            ),
                        }
                    }
                    RunnerV2RawOutcomeKindV1::Inapplicable => limit_inapplicable += 1,
                    RunnerV2RawOutcomeKindV1::Failed | RunnerV2RawOutcomeKindV1::Unsupported => {
                        panic!("unexpected limit outcome for {}", actual.cell_id().as_str())
                    }
                }
            }
            match actual.outcome() {
                RunnerV2RawOutcomeKindV1::Accepted => {
                    accepted += 1;
                    assert!(actual.diagnostic().is_none());
                }
                RunnerV2RawOutcomeKindV1::Refused => {
                    refused += 1;
                    let diagnostic = actual.diagnostic().expect("refusal diagnostic");
                    assert_eq!(diagnostic.code(), DiagnosticCodeV2::RunnerRefused);
                    assert_eq!(diagnostic.retryability(), RetryabilityV2::AfterInputChange);
                }
                RunnerV2RawOutcomeKindV1::Inapplicable => {
                    inapplicable += 1;
                    let diagnostic = actual.diagnostic().expect("inapplicable diagnostic");
                    assert_eq!(diagnostic.code(), DiagnosticCodeV2::RunnerNotRun);
                    assert_eq!(diagnostic.retryability(), RetryabilityV2::Never);
                    assert!(!diagnostic.prerequisites().is_empty());
                    assert!(diagnostic.repairs().is_empty());
                }
                RunnerV2RawOutcomeKindV1::Failed | RunnerV2RawOutcomeKindV1::Unsupported => {
                    panic!("source-frozen Stage-A implementation produced an unexpected outcome")
                }
            }
        }
        assert_eq!(accepted, 422);
        assert_eq!(refused, 389);
        assert_eq!(inapplicable, 56);
        assert_eq!(accepted + refused + inapplicable, 867);
        assert_eq!(numeric_total, 2_617);
        assert_eq!(diagnostic_total, 445);
        assert_eq!(repair_total, 388);
        assert_eq!(prerequisite_total, 56);
        assert_eq!(limit_accepted, 411);
        assert_eq!(limit_refused, 388);
        assert_eq!(limit_inapplicable, 53);
        assert_eq!(limit_numeric_total, 2_602);
        assert_eq!(limit_diagnostic_total, 441);
        assert_eq!(limit_repair_total, 388);
        assert_eq!(limit_prerequisite_total, 53);
        assert_eq!(above_ceiling, 201);
        assert_eq!(fixed_changed, 20);
        assert_eq!(below_minimum, 59);
        assert_eq!(joint_feasibility, 37);
        assert_eq!(checked_overflow, 71);
    }

    #[test]
    fn tightened_companion_declarations_match_production_and_remain_jointly_feasible() {
        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        for cell in declaration.cells() {
            let RunnerV2StageACellOperationV1::Limit {
                field,
                boundary,
                value,
            } = cell.operation()
            else {
                assert!(cell.companion_normalization().is_empty());
                continue;
            };
            assert!(
                cell.companion_normalization()
                    .windows(2)
                    .all(|pair| pair[0].field().ordinal() < pair[1].field().ordinal())
            );
            if !boundary.is_tightened() {
                assert!(cell.companion_normalization().is_empty());
                continue;
            }
            let TypedOptionV1::Present(value) = value else {
                assert!(cell.companion_normalization().is_empty());
                continue;
            };

            let mut candidate = RunnerLimitsV2::base(boundary.profile()).to_candidate();
            candidate
                .set_value(field, value)
                .expect("declared input has the exact field width");
            normalize_tightened_candidate_v1(&mut candidate, field)
                .expect("production companion normalization");
            for companion in cell.companion_normalization() {
                assert_eq!(candidate.value(companion.field()), companion.value());
                assert_eq!(
                    companion.value().width(),
                    companion.field().descriptor().width
                );
            }
            let admitted = RunnerLimitsV2::admit_family(
                boundary.profile(),
                candidate,
                RunnerFamilyLimitRequirementsV2 {
                    executable: declaration.limit_fixture().executable(),
                    family_rows_by_case: declaration.limit_fixture().family_rows_by_case(),
                    declared_minimums: &[],
                },
            );
            assert!(
                admitted.is_ok(),
                "tightened cell {} must remain feasible",
                cell.cell_id().as_str()
            );
        }
    }

    #[test]
    fn all_71_wrong_width_mutations_refuse_with_exact_field_width_and_unit() {
        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        assert_eq!(declaration.limit_mutation_obligations().len(), 71);
        for (index, (obligation, expected)) in declaration
            .limit_mutation_obligations()
            .iter()
            .zip(TEST_LIMIT_MUTATIONS_EXPECTED_V1)
            .enumerate()
        {
            assert_eq!(expected.ordinal, u16::try_from(index + 1).unwrap());
            assert_eq!(obligation.ordinal(), expected.ordinal);
            assert_eq!(obligation.stable_id().as_str(), expected.stable_id);
            assert_eq!(obligation.field(), expected.field);
            assert_eq!(obligation.field_name().as_str(), expected.field_name);
            assert_eq!(obligation.declared_width(), expected.declared_width);
            assert_eq!(obligation.opposite_width_zero(), expected.opposite_zero);
            assert_eq!(obligation.unit(), expected.unit);
            assert_eq!(
                obligation.expected_reason(),
                RunnerV2RawReasonV1::WrongPrimitiveWidth
            );
            assert_eq!(
                obligation.diagnostic_owner().as_str(),
                "fs-evidence-runner.runner-limits"
            );
            assert_eq!(obligation.repair_rank(), 1);
            assert_eq!(
                obligation.repair_kind(),
                RepairActionKindV2::UpdatePolicyOrCapability
            );
            assert_eq!(obligation.repair_target().as_str(), expected.field_name);

            let mut candidate = RunnerLimitsV2::base(RunProfileV2::Full).to_candidate();
            let violation = candidate
                .set_value(expected.field, expected.opposite_zero)
                .expect_err("opposite primitive width must refuse");
            assert_eq!(violation.kind(), RunnerLimitsViolationKindV2::WrongWidth);
            assert_eq!(violation.field(), expected.field);
            assert_eq!(violation.observed(), expected.opposite_zero);
            assert_eq!(violation.unit(), expected.unit);
            assert_eq!(
                violation.expected(),
                RunnerLimitExpectationV2::Width(expected.declared_width)
            );
            assert_eq!(violation.owner(), "fs-evidence-runner.runner-limits");
            assert_eq!(violation.repair_rank(), 1);
            assert_eq!(
                violation.repair_kind(),
                RepairActionKindV2::UpdatePolicyOrCapability
            );
            assert_eq!(violation.repair_target(), expected.field_name);
            let raw_reason = limit_violation_reason_v1(&violation);
            assert_eq!(raw_reason, RunnerV2RawReasonV1::WrongPrimitiveWidth);
            let diagnostic =
                raw_limit_diagnostic_v1(&violation, raw_reason).expect("bounded raw diagnostic");
            assert_eq!(diagnostic.code(), DiagnosticCodeV2::RunnerRefused);
            assert_eq!(
                diagnostic.owner().as_str(),
                "fs-evidence-runner.runner-limits"
            );
            assert_eq!(diagnostic.retryability(), RetryabilityV2::AfterInputChange);
            assert!(diagnostic.prerequisites().is_empty());
            assert_eq!(diagnostic.repairs().len(), 1);
            assert_eq!(diagnostic.repairs()[0].rank(), 1);
            assert_eq!(
                diagnostic.repairs()[0].kind(),
                RepairActionKindV2::UpdatePolicyOrCapability
            );
            assert_eq!(
                diagnostic.repairs()[0].target().as_str(),
                expected.field_name
            );
        }

        let exact = build_limit_mutation_obligations_v1().unwrap();
        validate_limit_mutation_obligations_exact_v1(&exact).unwrap();
        let independent_ids = TEST_LIMIT_MUTATIONS_EXPECTED_V1
            .iter()
            .map(|row| row.stable_id.to_owned())
            .collect::<Vec<_>>();
        assert_every_position_inventory_mutations_v1(
            &independent_ids,
            validate_stage_a_limit_mutation_inventory_v1,
            "runner_v2.base_values.limit_mutations",
            "restore-exact-limit-wrong-width-mutation-catalog",
        );

        for (index, expected) in TEST_LIMIT_MUTATIONS_EXPECTED_V1.iter().enumerate() {
            let assert_component = |rows: &[RunnerV2LimitMutationObligationV1],
                                    component: &'static str,
                                    expected_value: String,
                                    observed_value: String| {
                assert_semantic_inventory_mismatch_v1(
                    validate_limit_mutation_obligations_diagnostic_v1(rows),
                    "runner_v2.base_values.limit_mutations",
                    "restore-exact-limit-wrong-width-mutation-catalog",
                    index,
                    71,
                    expected.stable_id,
                    component,
                    &expected_value,
                    &observed_value,
                );
            };

            let mut wrong_ordinal = exact.clone();
            wrong_ordinal[index].ordinal = if expected.ordinal == 71 {
                1
            } else {
                expected.ordinal + 1
            };
            assert_component(
                &wrong_ordinal,
                "ordinal",
                expected.ordinal.to_string(),
                wrong_ordinal[index].ordinal.to_string(),
            );

            let mut wrong_field = exact.clone();
            wrong_field[index].field = TEST_LIMIT_MUTATIONS_EXPECTED_V1[(index + 1) % 71].field;
            assert_component(
                &wrong_field,
                "field",
                format!("{}:{}", expected.ordinal, expected.field_name),
                format!(
                    "{}:{}",
                    wrong_field[index].field.ordinal(),
                    wrong_field[index].field.descriptor().name
                ),
            );

            let mut wrong_field_name = exact.clone();
            wrong_field_name[index].field_name =
                stage_a_token("test.field_name", "substituted-limit-field").unwrap();
            assert_component(
                &wrong_field_name,
                "field-name",
                expected.field_name.to_owned(),
                "substituted-limit-field".to_owned(),
            );

            let mut wrong_width = exact.clone();
            wrong_width[index].declared_width = match expected.declared_width {
                RunnerLimitWidthV2::U32 => RunnerLimitWidthV2::U64,
                RunnerLimitWidthV2::U64 => RunnerLimitWidthV2::U32,
            };
            assert_component(
                &wrong_width,
                "declared-width",
                limit_width_name_v1(expected.declared_width).to_owned(),
                limit_width_name_v1(wrong_width[index].declared_width).to_owned(),
            );

            let mut wrong_opposite_zero = exact.clone();
            wrong_opposite_zero[index].opposite_width_zero = match expected.opposite_zero {
                RunnerLimitValueV2::U32(_) => RunnerLimitValueV2::U32(1),
                RunnerLimitValueV2::U64(_) => RunnerLimitValueV2::U64(1),
            };
            assert_component(
                &wrong_opposite_zero,
                "opposite-width-zero",
                limit_value_safe_name_v1(expected.opposite_zero),
                limit_value_safe_name_v1(wrong_opposite_zero[index].opposite_width_zero),
            );

            let mut wrong_unit = exact.clone();
            wrong_unit[index].unit = if expected.unit == RunnerLimitUnitV2::Count {
                RunnerLimitUnitV2::Records
            } else {
                RunnerLimitUnitV2::Count
            };
            assert_component(
                &wrong_unit,
                "unit",
                independent_limit_unit_name_v1(expected.unit).to_owned(),
                independent_limit_unit_name_v1(wrong_unit[index].unit).to_owned(),
            );

            let mut wrong_reason = exact.clone();
            wrong_reason[index].expected_reason = RunnerV2RawReasonV1::UnknownClosedValue;
            assert_component(
                &wrong_reason,
                "expected-reason",
                "05:wrong-primitive-width".to_owned(),
                "08:unknown-closed-value".to_owned(),
            );

            let mut wrong_owner = exact.clone();
            wrong_owner[index].diagnostic_owner =
                stage_a_token("test.owner", "fs-evidence-runner.wrong-owner").unwrap();
            assert_component(
                &wrong_owner,
                "diagnostic-owner",
                "fs-evidence-runner.runner-limits".to_owned(),
                "fs-evidence-runner.wrong-owner".to_owned(),
            );

            let mut wrong_rank = exact.clone();
            wrong_rank[index].repair_rank = 2;
            assert_component(&wrong_rank, "repair-rank", "1".to_owned(), "2".to_owned());

            let mut wrong_kind = exact.clone();
            wrong_kind[index].repair_kind = RepairActionKindV2::ChangeArguments;
            assert_component(
                &wrong_kind,
                "repair-kind",
                RepairActionKindV2::UpdatePolicyOrCapability
                    .code()
                    .to_string(),
                RepairActionKindV2::ChangeArguments.code().to_string(),
            );

            let mut wrong_target = exact.clone();
            wrong_target[index].repair_target =
                stage_a_token("test.target", "substituted-repair-target").unwrap();
            assert_component(
                &wrong_target,
                "repair-target",
                expected.field_name.to_owned(),
                "substituted-repair-target".to_owned(),
            );
        }
    }

    #[test]
    fn retained_domain_manifest_is_exact_50_with_frozen_facet_counts_and_mutations() {
        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        let rows = declaration.retained_domain_obligations();
        assert_eq!(rows.len(), 50);
        let mut counts = [0_usize; 8];
        let mut ids = BTreeSet::new();
        for (index, (row, expected)) in rows.iter().zip(TEST_RETAINED_EXPECTED_V1).enumerate() {
            assert_eq!(expected.ordinal, u16::try_from(index + 1).unwrap());
            assert_eq!(row.ordinal(), expected.ordinal);
            assert_eq!(row.stable_id().as_str(), expected.stable_id);
            assert_eq!(row.facet(), expected.facet);
            assert!(ids.insert(row.stable_id().as_str()));
            counts[usize::from(row.facet().code() - 1)] += 1;
        }
        assert_eq!(counts, [13, 4, 8, 6, 4, 5, 5, 5]);

        let exact = build_retained_domain_obligations_v1().unwrap();
        validate_retained_domain_obligations_exact_v1(&exact).unwrap();
        let independent_ids = TEST_RETAINED_EXPECTED_V1
            .iter()
            .map(|row| row.stable_id.to_owned())
            .collect::<Vec<_>>();
        assert_every_position_inventory_mutations_v1(
            &independent_ids,
            validate_stage_a_retained_domain_inventory_v1,
            "runner_v2.base_values.retained_domain",
            "restore-exact-retained-domain-obligation-catalog",
        );

        for (index, expected) in TEST_RETAINED_EXPECTED_V1.iter().enumerate() {
            let mut wrong_ordinal = exact.clone();
            wrong_ordinal[index].ordinal = if expected.ordinal == 50 {
                1
            } else {
                expected.ordinal + 1
            };
            assert_semantic_inventory_mismatch_v1(
                validate_retained_domain_obligations_diagnostic_v1(&wrong_ordinal),
                "runner_v2.base_values.retained_domain",
                "restore-exact-retained-domain-obligation-catalog",
                index,
                50,
                expected.stable_id,
                "ordinal",
                &expected.ordinal.to_string(),
                &wrong_ordinal[index].ordinal.to_string(),
            );

            let mut wrong_facet = exact.clone();
            wrong_facet[index].facet =
                if expected.facet == RunnerV2RetainedDomainFacetV1::NumericLiteral {
                    RunnerV2RetainedDomainFacetV1::Unit
                } else {
                    RunnerV2RetainedDomainFacetV1::NumericLiteral
                };
            assert_semantic_inventory_mismatch_v1(
                validate_retained_domain_obligations_diagnostic_v1(&wrong_facet),
                "runner_v2.base_values.retained_domain",
                "restore-exact-retained-domain-obligation-catalog",
                index,
                50,
                expected.stable_id,
                "facet",
                &expected.facet.code().to_string(),
                &wrong_facet[index].facet.code().to_string(),
            );
        }
    }

    #[test]
    fn common_future_owner_and_dependency_catalogs_fail_closed_under_mutation() {
        let common = build_common_requirements_v1().expect("common requirements");
        assert_eq!(common.len(), 31);
        validate_common_requirements_exact_v1(&common).unwrap();
        for (row, expected) in common.iter().zip(TEST_COMMON_EXPECTED_V1) {
            assert_eq!(row.ordinal(), expected.ordinal);
            assert_eq!(row.slot_id().as_str(), expected.slot_id);
            assert_eq!(row.api_generation(), RUNNER_SPEC_V2_API_GENERATION);
            assert_eq!(row.wire_version(), RUNNER_V2_WIRE_VERSION);
            assert_eq!(row.predecessor_policy(), RUNNER_V2_PREDECESSOR_POLICY);
            assert_eq!(
                row.semantic_owner().as_str(),
                "frankensim-epic-foundations-huq.24.1.1.1"
            );
            assert_eq!(row.realization_owner().as_str(), expected.realization_owner);
            assert_eq!(row.future_nominal_role().as_str(), expected.role);
            assert_eq!(row.future_domain().as_str(), expected.domain);
            assert_eq!(row.included_planes().mask(), expected.plane_mask);
            assert_eq!(row.fulfillment_stage(), expected.stage);
            assert_eq!(
                row.resolution_owner().as_str(),
                "frankensim-epic-foundations-huq.24.1.1.1.7"
            );
            assert!(matches!(row.future_root(), TypedOptionV1::Absent));
            assert_eq!(row.no_claim().as_str(), RUNNER_V2_BASE_VALUES_NO_CLAIM_V1);
        }
        let independent_common_ids = TEST_COMMON_EXPECTED_V1
            .iter()
            .map(|row| row.slot_id)
            .collect::<Vec<_>>();
        assert_every_position_typed_identity_mutations_v1(
            &common,
            &independent_common_ids,
            "runner_v2.base_values.requirements",
            "restore-exact-common-requirement-catalog",
            "runner-v2.unregistered-common-requirement-substitution.v1",
            "runner-v2.unregistered-common-requirement-extra.v1",
            validate_common_requirements_diagnostic_v1,
            |row| row.slot_id.as_str(),
            |row, identity| {
                row.slot_id = stage_a_token("test.common.slot_id", identity).unwrap();
            },
        );
        for (index, expected) in TEST_COMMON_EXPECTED_V1.iter().enumerate() {
            let assert_component = |rows: &[RunnerV2CommonContractRequirementV1],
                                    component: &'static str,
                                    expected_value: String,
                                    observed_value: String| {
                assert_semantic_inventory_mismatch_v1(
                    validate_common_requirements_diagnostic_v1(rows),
                    "runner_v2.base_values.requirements",
                    "restore-exact-common-requirement-catalog",
                    index,
                    31,
                    expected.slot_id,
                    component,
                    &expected_value,
                    &observed_value,
                );
            };

            let mut wrong_ordinal = common.clone();
            wrong_ordinal[index].ordinal = if expected.ordinal == 31 {
                1
            } else {
                expected.ordinal + 1
            };
            assert_component(
                &wrong_ordinal,
                "ordinal",
                expected.ordinal.to_string(),
                wrong_ordinal[index].ordinal.to_string(),
            );

            let mut wrong_semantic_owner = common.clone();
            wrong_semantic_owner[index].semantic_owner =
                stage_a_token("test.owner", "frankensim-wrong-semantic-owner").unwrap();
            assert_component(
                &wrong_semantic_owner,
                "semantic-owner",
                STAGE_A_INVENTORY_SEMANTIC_OWNER_V1.to_owned(),
                "frankensim-wrong-semantic-owner".to_owned(),
            );

            let mut wrong_realization_owner = common.clone();
            wrong_realization_owner[index].realization_owner =
                stage_a_token("test.owner", "frankensim-wrong-realization-owner").unwrap();
            assert_component(
                &wrong_realization_owner,
                "realization-owner",
                expected.realization_owner.to_owned(),
                "frankensim-wrong-realization-owner".to_owned(),
            );

            let mut wrong_role = common.clone();
            wrong_role[index].future_nominal_role =
                stage_a_token("test.role", "runner-v2-wrong-future-role-v1").unwrap();
            assert_component(
                &wrong_role,
                "future-nominal-role",
                expected.role.to_owned(),
                "runner-v2-wrong-future-role-v1".to_owned(),
            );

            let mut wrong_domain = common.clone();
            wrong_domain[index].future_domain =
                stage_a_token("test.domain", "runner-v2-wrong-future-domain-v1").unwrap();
            assert_component(
                &wrong_domain,
                "future-domain",
                expected.domain.to_owned(),
                "runner-v2-wrong-future-domain-v1".to_owned(),
            );

            let mut wrong_planes = common.clone();
            wrong_planes[index].included_planes =
                RunnerV2ContractPlaneSetV1::from_mask(expected.plane_mask ^ 0b001);
            assert_component(
                &wrong_planes,
                "included-planes",
                format!("0b{:03b}", expected.plane_mask),
                format!("0b{:03b}", wrong_planes[index].included_planes.mask()),
            );

            let mut wrong_stage = common.clone();
            wrong_stage[index].fulfillment_stage = match expected.stage {
                RunnerV2CommonFulfillmentStageV1::RuntimeEvidence => {
                    RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch
                }
                RunnerV2CommonFulfillmentStageV1::RoutesAndDispatch => {
                    RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction
                }
                RunnerV2CommonFulfillmentStageV1::LoggingAndReproduction => {
                    RunnerV2CommonFulfillmentStageV1::RuntimeEvidence
                }
            };
            assert_component(
                &wrong_stage,
                "fulfillment-stage",
                expected.stage.code().to_string(),
                wrong_stage[index].fulfillment_stage.code().to_string(),
            );

            let mut wrong_resolution_owner = common.clone();
            wrong_resolution_owner[index].resolution_owner =
                stage_a_token("test.owner", "frankensim-wrong-resolution-owner").unwrap();
            assert_component(
                &wrong_resolution_owner,
                "resolution-owner",
                "frankensim-epic-foundations-huq.24.1.1.1.7".to_owned(),
                "frankensim-wrong-resolution-owner".to_owned(),
            );

            let mut wrong_no_claim = common.clone();
            wrong_no_claim[index].no_claim =
                stage_a_token("test.no_claim", "runner-v2-wrong-no-claim").unwrap();
            assert_component(
                &wrong_no_claim,
                "no-claim",
                RUNNER_V2_BASE_VALUES_NO_CLAIM_V1.to_owned(),
                "runner-v2-wrong-no-claim".to_owned(),
            );
        }
        assert!(common.iter().all(|row| {
            row.api_generation() == RUNNER_SPEC_V2_API_GENERATION
                && row.wire_version() == RUNNER_V2_WIRE_VERSION
                && row.predecessor_policy() == RUNNER_V2_PREDECESSOR_POLICY
                && row.semantic_owner().as_str() == "frankensim-epic-foundations-huq.24.1.1.1"
                && row.resolution_owner().as_str() == "frankensim-epic-foundations-huq.24.1.1.1.7"
                && row.included_planes().mask() != 0
                && matches!(row.future_root(), TypedOptionV1::Absent)
        }));
        let mut missing = common.clone();
        missing.pop();
        assert_refuses(validate_common_requirements_exact_v1(&missing));
        let mut extra = common.clone();
        extra.push(common[0].clone());
        assert_refuses(validate_common_requirements_exact_v1(&extra));
        let mut reordered = common.clone();
        reordered.swap(0, 1);
        assert_refuses(validate_common_requirements_exact_v1(&reordered));
        let mut wrong_owner = common.clone();
        wrong_owner[0].realization_owner =
            stage_a_token("test.owner", "frankensim-epic-foundations-huq.24.1.1.1.6").unwrap();
        assert_refuses(validate_common_requirements_exact_v1(&wrong_owner));
        let mut wrong_role = common.clone();
        wrong_role[0].future_nominal_role =
            stage_a_token("test.role", "wrong-nominal-role-v1").unwrap();
        assert_refuses(validate_common_requirements_exact_v1(&wrong_role));
        let mut wrong_plane = common;
        wrong_plane[0].included_planes = RunnerV2ContractPlaneSetV1::from_mask(7);
        assert_refuses(validate_common_requirements_exact_v1(&wrong_plane));

        let future = build_future_sources_v1().expect("future sources");
        assert_eq!(future.len(), 13);
        validate_future_sources_exact_v1(&future).unwrap();
        for (row, expected) in future.iter().zip(TEST_FUTURE_SOURCES_EXPECTED_V1) {
            assert_eq!(row.final_ordinal(), expected.final_ordinal);
            assert_eq!(row.path().as_str(), expected.path);
            assert!(matches!(row.future_content_root(), TypedOptionV1::Absent));
        }
        let independent_future_ids = TEST_FUTURE_SOURCES_EXPECTED_V1
            .iter()
            .map(|row| row.path)
            .collect::<Vec<_>>();
        assert_every_position_typed_identity_mutations_v1(
            &future,
            &independent_future_ids,
            "runner_v2.base_values.future_sources",
            "restore-exact-future-source-catalog",
            "crates/fs-evidence-runner/src/runner_v2/unregistered_substitution.rs",
            "crates/fs-evidence-runner/src/runner_v2/unregistered_extra.rs",
            validate_future_sources_diagnostic_v1,
            |row| row.path.as_str(),
            |row, identity| {
                row.path = stage_a_path("test.future.path", identity).unwrap();
            },
        );
        for (index, expected) in TEST_FUTURE_SOURCES_EXPECTED_V1.iter().enumerate() {
            let mut wrong_ordinal = future.clone();
            wrong_ordinal[index].final_ordinal = if expected.final_ordinal == 40 {
                28
            } else {
                expected.final_ordinal + 1
            };
            assert_semantic_inventory_mismatch_v1(
                validate_future_sources_diagnostic_v1(&wrong_ordinal),
                "runner_v2.base_values.future_sources",
                "restore-exact-future-source-catalog",
                index,
                13,
                expected.path,
                "final-ordinal",
                &expected.final_ordinal.to_string(),
                &wrong_ordinal[index].final_ordinal.to_string(),
            );
        }
        let mut missing = future.clone();
        missing.pop();
        assert_refuses(validate_future_sources_exact_v1(&missing));
        let mut extra = future.clone();
        extra.push(future[0].clone());
        assert_refuses(validate_future_sources_exact_v1(&extra));
        let mut reordered = future.clone();
        reordered.swap(0, 1);
        assert_refuses(validate_future_sources_exact_v1(&reordered));
        let mut wrong_path = future;
        wrong_path[0].path = stage_a_path("test.path", "wrong/future.rs").unwrap();
        assert_refuses(validate_future_sources_exact_v1(&wrong_path));

        let owner = build_owner_source_fragment_v1().expect("owner source");
        validate_owner_source_fragment_exact_v1(&owner).unwrap();
        for (index, (row, expected_path)) in owner
            .iter()
            .zip(TEST_OWNER_SOURCE_PATHS_EXPECTED_V1)
            .enumerate()
        {
            assert_eq!(row.path().as_str(), expected_path);
            assert_eq!(
                row.content_root(),
                independently_expected_owner_source_root_v1(index).unwrap()
            );
        }
        assert_every_position_source_path_mutations_v1(
            &owner,
            &TEST_OWNER_SOURCE_PATHS_EXPECTED_V1,
            "runner_v2.base_values.owner_sources",
            "restore-exact-owner-source-fragment",
            validate_owner_source_fragment_diagnostic_v1,
            |row| row.path.as_str(),
            |row, path| row.path = path,
        );
        for (index, expected_path) in TEST_OWNER_SOURCE_PATHS_EXPECTED_V1.iter().enumerate() {
            let expected_root = independently_expected_owner_source_root_v1(index).unwrap();
            let mut wrong_root = owner.clone();
            wrong_root[index].content_root = owner[(index + 1) % 2].content_root;
            assert_semantic_inventory_mismatch_v1(
                validate_owner_source_fragment_diagnostic_v1(&wrong_root),
                "runner_v2.base_values.owner_sources",
                "restore-exact-owner-source-fragment",
                index,
                2,
                expected_path,
                "content-root",
                &source_root_safe_name_v1(expected_root),
                &source_root_safe_name_v1(wrong_root[index].content_root),
            );
        }
        let mut missing = owner.clone();
        missing.pop();
        assert_refuses(validate_owner_source_fragment_exact_v1(&missing));
        let mut extra = owner.clone();
        extra.push(owner[0].clone());
        assert_refuses(validate_owner_source_fragment_exact_v1(&extra));
        let mut reordered = owner.clone();
        reordered.swap(0, 1);
        assert_refuses(validate_owner_source_fragment_exact_v1(&reordered));
        let mut wrong_path = owner;
        wrong_path[0].path = stage_a_path("test.path", "wrong/owner.rs").unwrap();
        assert_refuses(validate_owner_source_fragment_exact_v1(&wrong_path));
        let mut wrong_root = build_owner_source_fragment_v1().unwrap();
        wrong_root[0].content_root = RunnerV2StageASourceMemberRootV1::from_content_hash(
            hash_domain("org.frankensim.test.owner-source.v1", b"wrong owner source"),
        );
        assert_refuses(validate_owner_source_fragment_exact_v1(&wrong_root));
        let mut swapped_roots = build_owner_source_fragment_v1().unwrap();
        let first_root = swapped_roots[0].content_root;
        swapped_roots[0].content_root = swapped_roots[1].content_root;
        swapped_roots[1].content_root = first_root;
        assert_refuses(validate_owner_source_fragment_exact_v1(&swapped_roots));
        let mut resealed_owner = build_owner_source_fragment_v1().unwrap();
        resealed_owner[0].content_root =
            RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                "org.frankensim.fs-evidence-runner.runner-v2.stage-a.source-member.v1",
                b"mutated owner bytes with a correctly resealed candidate root",
            ));
        assert_refuses(validate_owner_source_fragment_exact_v1(&resealed_owner));

        let dependency = build_dependency_source_closure_v1().expect("dependency closure");
        validate_dependency_source_closure_exact_v1(&dependency).unwrap();
        assert_eq!(dependency.len(), 16);
        for (index, (row, expected_path)) in dependency
            .iter()
            .zip(TEST_DEPENDENCY_SOURCE_PATHS_EXPECTED_V1)
            .enumerate()
        {
            assert_eq!(row.path().as_str(), expected_path);
            assert_eq!(
                row.content_root(),
                independently_expected_dependency_source_root_v1(index).unwrap()
            );
        }
        assert_every_position_source_path_mutations_v1(
            &dependency,
            &TEST_DEPENDENCY_SOURCE_PATHS_EXPECTED_V1,
            "runner_v2.base_values.dependency_sources",
            "restore-exact-dependency-source-closure",
            validate_dependency_source_closure_diagnostic_v1,
            |row| row.path.as_str(),
            |row, path| row.path = path,
        );
        for (index, expected_path) in TEST_DEPENDENCY_SOURCE_PATHS_EXPECTED_V1.iter().enumerate() {
            let expected_root = independently_expected_dependency_source_root_v1(index).unwrap();
            let mut wrong_root = dependency.clone();
            wrong_root[index].content_root = dependency[(index + 1) % 16].content_root;
            assert_semantic_inventory_mismatch_v1(
                validate_dependency_source_closure_diagnostic_v1(&wrong_root),
                "runner_v2.base_values.dependency_sources",
                "restore-exact-dependency-source-closure",
                index,
                16,
                expected_path,
                "content-root",
                &source_root_safe_name_v1(expected_root),
                &source_root_safe_name_v1(wrong_root[index].content_root),
            );
        }
        let _exact_versions = build_version_requirements_v1(&dependency).unwrap();
        let mut changed_content = dependency.clone();
        changed_content[0].content_root =
            RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                "org.frankensim.test.source-movement.v1",
                b"one changed dependency",
            ));
        assert_refuses(validate_dependency_source_closure_exact_v1(
            &changed_content,
        ));
        assert_refuses(build_version_requirements_v1(&changed_content));
        let mut swapped_roots = dependency.clone();
        let first_root = swapped_roots[0].content_root;
        swapped_roots[0].content_root = swapped_roots[1].content_root;
        swapped_roots[1].content_root = first_root;
        assert_refuses(validate_dependency_source_closure_exact_v1(&swapped_roots));
        let mut resealed_dependency = dependency.clone();
        resealed_dependency[0].content_root =
            RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                "org.frankensim.fs-evidence-runner.runner-v2.stage-a.dependency-source-member.v1",
                b"mutated dependency bytes with a correctly resealed candidate root",
            ));
        assert_refuses(validate_dependency_source_closure_exact_v1(
            &resealed_dependency,
        ));
        let mut wrong_domain = dependency.clone();
        wrong_domain[0].content_root =
            RunnerV2StageASourceMemberRootV1::from_content_hash(hash_domain(
                "org.frankensim.fs-evidence-runner.runner-v2.stage-a.source-member.v1",
                DEPENDENCY_SOURCE_DECLARATIONS_V1[0].bytes,
            ));
        assert_refuses(validate_dependency_source_closure_exact_v1(&wrong_domain));
        let mut missing = dependency.clone();
        missing.pop();
        assert_refuses(validate_dependency_source_closure_exact_v1(&missing));
        let mut extra = dependency.clone();
        extra.push(dependency[0].clone());
        assert_refuses(validate_dependency_source_closure_exact_v1(&extra));
        let mut reordered = dependency.clone();
        reordered.swap(0, 1);
        assert_refuses(validate_dependency_source_closure_exact_v1(&reordered));
        let mut wrong_path = dependency;
        wrong_path[0].path = stage_a_path("test.path", "wrong/dependency.rs").unwrap();
        assert_refuses(validate_dependency_source_closure_exact_v1(&wrong_path));
    }

    #[test]
    fn exact_empty_feature_guard_rejects_explicit_and_implicit_cargo_features() {
        assert!(!cargo_manifest_declares_feature_v1(include_str!(
            "../../../Cargo.toml"
        )));
        assert!(cargo_manifest_declares_feature_v1(
            "[features]\ndefault = []\n"
        ));
        assert!(cargo_manifest_declares_feature_v1(
            "[dependencies]\nexample = { version = \"1\", optional = true }\n"
        ));
        assert!(cargo_manifest_declares_feature_v1(
            "[dependencies.example]\noptional = true\n"
        ));
        assert!(cargo_manifest_declares_feature_v1(
            "[build-dependencies]\nexample = { optional=true }\n"
        ));
        assert!(cargo_manifest_declares_feature_v1(
            "[target.'cfg(unix)'.dependencies]\nexample = { optional = true }\n"
        ));
        assert!(cargo_manifest_declares_feature_v1(
            "[target.'cfg(windows)'.build-dependencies.example]\noptional=true\n"
        ));
        assert!(!cargo_manifest_declares_feature_v1(
            "[dependencies]\n# example = { optional = true }\n"
        ));
        assert!(!cargo_manifest_declares_feature_v1(
            "[dependencies]\nexample = { note = \"optional=true\" }\n"
        ));
        assert!(!cargo_manifest_declares_feature_v1(
            "[dependencies]\nexample = { optional = false }\n"
        ));
        assert!(!cargo_manifest_declares_feature_v1(
            "[package.metadata]\noptional = true\n"
        ));
    }

    #[test]
    fn diagnostic_redaction_consumes_a_real_canary_without_echo() {
        let canary = "runner-v2-sensitive-canary-direct-test-must-never-be-retained";
        let error = consume_and_redact_stage_a_value_v1(canary.to_owned());
        assert_eq!(error.observed(), "<redacted:sensitive-or-ambient>");
        assert!(!error.observed().contains(canary));
        assert!(!format!("{error}").contains(canary));
        assert!(!format!("{error:?}").contains(canary));

        let cell = evaluate_meta_operation_v1(
            stage_a_token(
                "test.cell_id",
                "runner-v2.base-values.meta-redaction-direct-test.v1",
            )
            .unwrap(),
            RunnerV2StageAMetaOperationV1::DiagnosticRedaction,
        )
        .unwrap();
        assert_eq!(cell.outcome(), RunnerV2RawOutcomeKindV1::Accepted);
        assert!(!format!("{cell:?}").contains("runner-v2-sensitive-canary"));
    }

    #[test]
    fn five_explicits_route_fixture_and_schema_deferral_are_exact() {
        let declaration = declare_24_1_1_1_1_v1().expect("declaration");
        let fixture = declaration.limit_fixture();
        assert!(fixture.executable());
        assert_eq!(fixture.family_rows_by_case(), &[0]);
        assert!(fixture.declared_minimums_present_empty());
        assert_eq!(fixture.lifecycle_document_structural_minimum(), 5);

        let five = declaration.five_explicits();
        assert!(five.numeric_inputs_present_empty());
        assert!(five.numeric_grants_present_empty());
        assert!(five.expected_numeric_observations_present_empty());
        assert_eq!(
            five.seed().inapplicable_reason(),
            Some(SeedInapplicableCodeV1::NoRandomnessByContract)
        );
        assert_eq!(five.budgets().rows().len(), 7);
        let expected_budgets = [
            (
                BaseCoverageCloseBudgetAxisV1::Time,
                BaseCoverageCloseBudgetValueV1::U64(60_000_000_000),
                BaseCoverageCloseBudgetValueV1::U64(45_000_000_000),
                LogicalUnitV2::Nanoseconds,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::Memory,
                BaseCoverageCloseBudgetValueV1::U64(536_870_912),
                BaseCoverageCloseBudgetValueV1::U64(402_653_184),
                LogicalUnitV2::LogicalBytes,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::LogicalWork,
                BaseCoverageCloseBudgetValueV1::U128(1_000_000),
                BaseCoverageCloseBudgetValueV1::U128(750_000),
                LogicalUnitV2::Operations,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::Processes,
                BaseCoverageCloseBudgetValueV1::U32(1),
                BaseCoverageCloseBudgetValueV1::U32(0),
                LogicalUnitV2::Count,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::Artifacts,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                LogicalUnitV2::EncodedBytes,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::Output,
                BaseCoverageCloseBudgetValueV1::U64(5_242_880),
                BaseCoverageCloseBudgetValueV1::U64(4_194_304),
                LogicalUnitV2::EncodedBytes,
            ),
            (
                BaseCoverageCloseBudgetAxisV1::Logs,
                BaseCoverageCloseBudgetValueV1::U64(67_108_864),
                BaseCoverageCloseBudgetValueV1::U64(50_331_648),
                LogicalUnitV2::EncodedBytes,
            ),
        ];
        for (row, expected) in five.budgets().rows().iter().zip(expected_budgets) {
            assert_eq!(row.axis(), expected.0);
            assert_eq!(row.hard(), expected.1);
            assert_eq!(row.soft(), expected.2);
            assert_eq!(row.unit().unit(), expected.3);
        }
        assert_eq!(five.capability_registry().rows().len(), 5);
        assert_eq!(five.capability_profile_registry().rows().len(), 5);
        assert_eq!(
            five.capability_contract().profile(),
            BaseCoverageCloseCapabilityProfileV1::None
        );
        assert!(five.capability_contract().required().is_empty());
        assert!(five.capability_contract().permitted().is_empty());
        assert!(!cargo_manifest_declares_feature_v1(include_str!(
            "../../../Cargo.toml"
        )));

        let route = declaration.route();
        assert_eq!(
            route.route_id().as_str(),
            RUNNER_V2_BASE_VALUES_LOCAL_ROUTE_ID_V1
        );
        assert_eq!(route.class(), RunnerV2LocalRouteClassV1::LocalOnly);
        assert_eq!(
            route.public_entry_point(),
            RUNNER_V2_BASE_VALUES_PUBLIC_ENTRY_POINT_V1
        );
        assert_eq!(
            route.execution_owner().as_str(),
            "frankensim-epic-foundations-huq.24.1.1.1.7"
        );
        assert!(matches!(route.external_driver(), TypedOptionV1::Absent));
        assert_eq!(RUNNER_V2_BASE_VALUES_LOCAL_IN_PROCESS_ROUTE_COUNT_V1, 1);
        assert_eq!(RUNNER_V2_BASE_VALUES_EXECUTION_OWNED_ROUTE_COUNT_V1, 0);
        assert_eq!(RUNNER_V2_BASE_VALUES_CONTRIBUTION_ONLY_ROUTE_COUNT_V1, 0);
        assert_eq!(RUNNER_V2_BASE_VALUES_INAPPLICABLE_ROUTE_COUNT_V1, 0);

        let expected_canonical_schema_names = [
            "runner-v2-stage-a-declaration-root-v1",
            "runner-v2-stage-a-oracle-root-v1",
            "runner-v2-stage-a-case-manifest-root-v1",
            "runner-v2-stage-a-schema-inventory-root-v1",
            "runner-v2-stage-a-feature-declaration-root-v1",
            "runner-v2-stage-a-five-explicits-root-v1",
            "runner-v2-stage-a-source-member-root-v1",
            "runner-v2-limit-boundary-kind-v1",
            "runner-v2-limit-fixture-declaration-v1",
            "runner-v2-limit-companion-normalization-v1",
            "runner-v2-stage-a-expected-partition-v1",
            "runner-v2-stage-a-cell-group-v1",
            "runner-v2-retained-domain-facet-v1",
            "runner-v2-retained-domain-obligation-v1",
            "runner-v2-limit-mutation-obligation-v1",
            "runner-v2-stage-a-version-requirements-v1",
            "runner-v2-stage-a-five-explicits-v1",
            "runner-v2-contract-plane-set-v1",
            "runner-v2-common-fulfillment-stage-v1",
            "runner-v2-unavailable-common-root-v1",
            "runner-v2-common-contract-requirement-v1",
            "runner-v2-future-source-requirement-v1",
            "runner-v2-owner-source-member-v1",
            "runner-v2-dependency-source-member-v1",
            "runner-v2-schema-impact-deferral-v1",
            "runner-v2-rootless-ac58-fragment-v1",
            "runner-v2-local-route-class-v1",
            "runner-v2-local-route-declaration-v1",
            "runner-v2-stage-a-inapplicability-declaration-v1",
            "runner-v2-stage-a-oracle-row-v1",
            "runner-v2-stage-a-projection-row-v1",
            "runner-v2-stage-a-meta-operation-v1",
            "runner-v2-stage-a-cell-operation-v1",
            "runner-v2-stage-a-cell-declaration-v1",
            "runner-v2-base-values-stage-a-declaration-v1",
            "runner-v2-raw-outcome-reason-contract-v1",
            "runner-v2-safe-numeric-value-v1",
            "runner-v2-safe-numeric-unit-v1",
            "runner-v2-safe-numeric-observation-v1",
            "runner-v2-raw-repair-v1",
            "runner-v2-raw-diagnostic-v1",
            "runner-v2-raw-cell-observation-v1",
        ];
        let expected_rootless_schema_names = ["runner-v2-local-work-package-handoff-v1"];
        assert_eq!(
            STAGE_A_CANONICAL_SCHEMA_NAMES_V1,
            expected_canonical_schema_names
        );
        assert_eq!(
            STAGE_A_ROOTLESS_HANDOFF_SCHEMA_NAME_V1,
            expected_rootless_schema_names[0]
        );
        validate_stage_a_schema_partition_exact_v1(
            &expected_canonical_schema_names,
            &expected_rootless_schema_names,
        )
        .expect("separately typed literal schema inventories");
        let canonical_schema_names = STAGE_A_CANONICAL_SCHEMA_NAMES_V1
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        assert_eq!(
            STAGE_A_CANONICAL_SCHEMA_NAMES_V1.len(),
            RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1
        );
        assert_eq!(
            canonical_schema_names.len(),
            RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1
        );
        assert_eq!(expected_rootless_schema_names.len(), 1);
        assert_eq!(
            STAGE_A_CANONICAL_SCHEMA_NAMES_V1.len() + expected_rootless_schema_names.len(),
            RUNNER_V2_BASE_VALUES_OWNED_SCHEMA_COUNT_V1
        );
        let deferral = declaration.schema_impact_deferral();
        assert_eq!(
            deferral.resolution_owner().as_str(),
            "frankensim-epic-foundations-huq.24.1.1.1.3"
        );
        assert_eq!(
            deferral.canonical_schema_names().len(),
            RUNNER_V2_BASE_VALUES_CANONICAL_SCHEMA_COUNT_V1
        );
        assert!(
            deferral
                .canonical_schema_names()
                .iter()
                .any(|name| name.as_str() == "runner-v2-raw-outcome-reason-contract-v1")
        );
        assert_eq!(
            deferral
                .canonical_schema_names()
                .iter()
                .filter(|name| name.as_str() == "runner-v2-schema-impact-deferral-v1")
                .count(),
            1
        );
        assert_eq!(
            deferral
                .canonical_schema_names()
                .iter()
                .filter(|name| name.as_str() == "runner-v2-raw-outcome-reason-contract-v1")
                .count(),
            1
        );
        assert!(matches!(
            deferral.future_manifest_root(),
            TypedOptionV1::Absent
        ));
        assert!(
            !deferral
                .canonical_schema_names()
                .iter()
                .any(|name| name.as_str() == "runner-v2-local-work-package-handoff-v1")
        );
        assert_eq!(
            deferral
                .canonical_schema_names()
                .iter()
                .map(RunnerV2CanonicalSchemaNameV1::as_str)
                .collect::<Vec<_>>(),
            expected_canonical_schema_names
        );
        assert_all_canonical_schema_position_mutations_v1(&expected_canonical_schema_names);
        let independent_rootless_schema_names =
            vec!["runner-v2-local-work-package-handoff-v1".to_owned()];
        assert_every_position_inventory_mutations_v1(
            &independent_rootless_schema_names,
            |rows| {
                let borrowed = rows.iter().map(String::as_str).collect::<Vec<_>>();
                validate_stage_a_rootless_schema_inventory_v1(&borrowed)
            },
            "runner_v2.base_values.schema_inventory.rootless_handoff",
            "restore-sole-rootless-handoff-schema-inventory",
        );

        let mut rootless_as_canonical = expected_canonical_schema_names;
        rootless_as_canonical[0] = expected_rootless_schema_names[0];
        assert_schema_error_v1(
            validate_stage_a_canonical_schema_names_exact_v1(&rootless_as_canonical),
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.schema_inventory.canonical",
            0,
        );
        assert_schema_error_v1(
            validate_stage_a_canonical_schema_names_exact_v1(
                &expected_canonical_schema_names[..41],
            ),
            ConstructionErrorKindV2::Missing,
            "runner_v2.base_values.schema_inventory.canonical",
            41,
        );
        assert_schema_error_v1(
            validate_stage_a_rootless_schema_names_exact_v1(&[]),
            ConstructionErrorKindV2::Missing,
            "runner_v2.base_values.schema_inventory.rootless_handoff",
            0,
        );
        assert_schema_error_v1(
            validate_stage_a_rootless_schema_names_exact_v1(&[expected_canonical_schema_names[0]]),
            ConstructionErrorKindV2::Incompatible,
            "runner_v2.base_values.schema_inventory.rootless_handoff",
            0,
        );
        assert_schema_error_v1(
            validate_stage_a_rootless_schema_names_exact_v1(&[
                expected_rootless_schema_names[0],
                expected_rootless_schema_names[0],
            ]),
            ConstructionErrorKindV2::Duplicate,
            "runner_v2.base_values.schema_inventory.rootless_handoff",
            1,
        );
        assert_schema_error_v1(
            validate_stage_a_rootless_schema_names_exact_v1(&[
                expected_rootless_schema_names[0],
                "runner-v2-unregistered-extra-rootless-v1",
            ]),
            ConstructionErrorKindV2::Unexpected,
            "runner_v2.base_values.schema_inventory.rootless_handoff",
            1,
        );
        assert!(
            expected_rootless_schema_names.windows(2).next().is_none(),
            "a one-row rootless inventory has no nonidentity reordering; cross-partition movement is tested separately"
        );
        assert_eq!(
            declaration.ac58().semantic_type().as_str(),
            expected_rootless_schema_names[0]
        );
        assert_eq!(
            declaration.ac58().disposition(),
            CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame
        );
        assert_eq!(
            declaration.ac58().migration_policy(),
            CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor
        );
        assert!(declaration.ac58().authority_surfaces_present_empty());
    }

    #[test]
    fn float_named_total_order_covers_zero_infinities_subnormals_and_nan_payloads() {
        assert_eq!(
            F32BitsV2::from_bits(0x8000_0000).ieee_total_cmp_v1(F32BitsV2::from_bits(0)),
            Ordering::Less
        );
        assert_eq!(
            F32BitsV2::from_bits(1).ieee_total_cmp_v1(F32BitsV2::from_bits(0)),
            Ordering::Greater
        );
        assert_eq!(
            F32BitsV2::from_bits(f32::NEG_INFINITY.to_bits())
                .ieee_total_cmp_v1(F32BitsV2::from_bits(f32::INFINITY.to_bits())),
            Ordering::Less
        );
        assert_ne!(
            F32BitsV2::from_bits(0x7fc0_0001).ieee_total_cmp_v1(F32BitsV2::from_bits(0x7fc0_0002)),
            Ordering::Equal
        );
        assert_eq!(
            F64BitsV2::from_bits(0x8000_0000_0000_0000).ieee_total_cmp_v1(F64BitsV2::from_bits(0)),
            Ordering::Less
        );
        assert_eq!(
            F64BitsV2::from_bits(1).ieee_total_cmp_v1(F64BitsV2::from_bits(0)),
            Ordering::Greater
        );
        assert_eq!(
            F64BitsV2::from_bits(f64::NEG_INFINITY.to_bits())
                .ieee_total_cmp_v1(F64BitsV2::from_bits(f64::INFINITY.to_bits())),
            Ordering::Less
        );
        assert_ne!(
            F64BitsV2::from_bits(0x7ff8_0000_0000_0001)
                .ieee_total_cmp_v1(F64BitsV2::from_bits(0x7ff8_0000_0000_0002)),
            Ordering::Equal
        );
    }
}
