//! Five source-closed, non-wire base E2E projections.
//!
//! These projections exercise real public constructors and validators in
//! process. They do not create or execute the downstream-owned shell scripts,
//! emit lifecycle records, publish bundles, or mint authority.

use crate::budget::{
    RunnerBudgetFieldV2, RunnerBudgetsCandidateV2, RunnerBudgetsV2,
};
use crate::canonical::CanonicalFrameV1;
use crate::capability::{
    NarrowedPolicyViewV2, OverlapPolicyRegistrationV2, RootCapabilityPolicyV2,
    RootPolicyRegistryProjectionV2, expected_rights, validate_policy_against_selection_v2,
};
use crate::catalog::{
    ArtifactDispositionV2, ArtifactRoleV2, DestinationAdmissionModeV2, DiagnosticCodeV2,
    DigestRoleV2, LifecycleRecordKindV2, LogicalExtentAxisV2, LogicalUnitV2,
    NotRunCauseCodeV2, OverlapPolicyRelationV2, PlatformPathProfileV2, ProofExitV2,
    PublicationProtocolV2, RefusedReasonV2, RepairActionKindV2, RetryabilityV2,
    RootCapabilityAccessV2, RootCapabilityRightV2, RootClassV2, RunProfileV2,
    RunnerApiGeneration, RunnerCommandV2, RunnerWireVersion, StateBearingRecordRoleV2,
    TypedOptionTagV1, TypedValueTagV2, WirePredecessorPolicyV1,
};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::diagnostic::{
    ActionableDiagnosticV2, DiagnosticCodeRefV2, DiagnosticEnvelopeGrantsV2, DiagnosticValueV2,
    RepairActionV2,
};
use crate::identity::{
    BuildIdentityRootV2, CancelledStopRootV2, DrainedInternalErrorRootV2, NoClaimScopeRootV1,
    SourceIdentityRootV2, TimedOutStopRootV2, ToolchainIdentityRootV2,
};
use crate::limits::{
    ArtifactStorageProjectionV2, PublicationStorageProjectionV2, RUNNER_LIMIT_DESCRIPTORS_V2,
    RunnerFamilyLimitRequirementsV2, RunnerLimitFieldV2, RunnerLimitTightenabilityV2,
    RunnerLimitValueV2, RunnerLimitsV2, RunnerLimitsViolationKindV2,
    SystemObjectStorageProjectionV2, SystemPublicationObjectRoleV2,
};
use crate::logging::{
    BaseE2eLogEventV1, BaseE2eLogFieldV1, BaseE2eLogKindV1, BaseE2eLogV1, BaseE2eOutcomeV1,
    SymbolicReproductionArgV1,
};
use crate::path::{
    ContentStoreObjectKeyV1, LogicalBundlePathV1, PathSetAdjudicationV1,
    adjudicate_logical_bundle_path_set,
};
use crate::publication::{
    PublicationSelectionV2, PublicationTargetV2, SymbolicCommandResultPlanV2,
};
use crate::state::{
    NotRunBasisV2, NotRunCauseV2, PresentedDrainRootKindV2, StateValidationInputV2,
    validate_state_v2,
};
use crate::value::{RationalV2, StableTokenV2, TypedValueV2};
use fs_blake3::{ContentHash, hash_domain};

/// Overall non-wire projection root domain.
pub const BASE_E2E_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-projection.v1";
/// Per-journey non-wire projection root domain.
pub const BASE_E2E_JOURNEY_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-journey-projection.v1";
/// Domain for one exact embedded source file's raw bytes.
pub const BASE_SOURCE_FILE_CONTENT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-file-content.v1";
/// Domain for one path-, length-, and content-bound source entry.
pub const BASE_SOURCE_FILE_ENTRY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-file-entry.v1";
/// Domain for the exact ordered base-schema source closure.
pub const BASE_SOURCE_CLOSURE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-source-closure.v1";
/// Domain for the immutable, result-free coverage-source inventory.
pub const BASE_COVERAGE_INVENTORY_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-coverage-inventory.v1";

#[derive(Clone, Copy)]
struct EmbeddedSourceFileV1 {
    path: &'static str,
    bytes: &'static [u8],
}

// This is the exact bytewise-lexicographic source set owned or consumed by the
// base leaf. `include_bytes!` makes the compiled projection move whenever any
// source, contract, manifest, or lock input changes.
const EMBEDDED_SOURCE_FILES_V1: [EmbeddedSourceFileV1; 21] = [
    EmbeddedSourceFileV1 {
        path: "Cargo.lock",
        bytes: include_bytes!("../../../Cargo.lock"),
    },
    EmbeddedSourceFileV1 {
        path: "Cargo.toml",
        bytes: include_bytes!("../../../Cargo.toml"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/CONTRACT.md",
        bytes: include_bytes!("../CONTRACT.md"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/Cargo.toml",
        bytes: include_bytes!("../Cargo.toml"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/budget.rs",
        bytes: include_bytes!("budget.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/canonical.rs",
        bytes: include_bytes!("canonical.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/capability.rs",
        bytes: include_bytes!("capability.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/catalog.rs",
        bytes: include_bytes!("catalog.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/command.rs",
        bytes: include_bytes!("command.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/construction.rs",
        bytes: include_bytes!("construction.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/dependency.rs",
        bytes: include_bytes!("dependency.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/diagnostic.rs",
        bytes: include_bytes!("diagnostic.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/identity.rs",
        bytes: include_bytes!("identity.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/lib.rs",
        bytes: include_bytes!("lib.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/limits.rs",
        bytes: include_bytes!("limits.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/logging.rs",
        bytes: include_bytes!("logging.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/path.rs",
        bytes: include_bytes!("path.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/projection.rs",
        bytes: include_bytes!("projection.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/publication.rs",
        bytes: include_bytes!("publication.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/state.rs",
        bytes: include_bytes!("state.rs"),
    },
    EmbeddedSourceFileV1 {
        path: "crates/fs-evidence-runner/src/value.rs",
        bytes: include_bytes!("value.rs"),
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
    bytes: Vec<u8>,
}

impl BaseSourceClosureInputV1 {
    /// Constructs one raw reconstruction input without validating it.
    #[must_use]
    pub fn presented(path: impl Into<String>, bytes: impl Into<Vec<u8>>) -> Self {
        Self {
            path: path.into(),
            bytes: bytes.into(),
        }
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
}

/// One exact path-, length-, and content-bound source-closure entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseSourceClosureEntryV1 {
    path: &'static str,
    encoded_bytes: u64,
    content_root: ContentHash,
    entry_root: ContentHash,
}

impl BaseSourceClosureEntryV1 {
    /// Exact workspace-relative source path.
    #[must_use]
    pub const fn path(self) -> &'static str {
        self.path
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

    /// Domain-separated root binding path, length, and content root.
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
/// ```compile_fail
/// use fs_evidence_runner::projection::RunnerV2BaseSourceClosureV1;
///
/// let mut closure = RunnerV2BaseSourceClosureV1::frozen().unwrap();
/// closure.entries.clear();
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseSourceClosureV1 {
    entries: Box<[BaseSourceClosureEntryV1]>,
    root: ContentHash,
}

impl RunnerV2BaseSourceClosureV1 {
    /// Reconstructs the closure from the exact compile-time embedded inputs.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let inputs = EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| BaseSourceClosureInputV1::presented(file.path, file.bytes.to_vec()))
            .collect::<Vec<_>>();
        Self::reconstruct(&inputs)
    }

    /// Checks and reconstructs the one exact ordered source closure.
    ///
    /// Duplicate, missing, extra, reordered, or byte-mutated inputs refuse
    /// before a closure root is returned.
    pub fn reconstruct(inputs: &[BaseSourceClosureInputV1]) -> Result<Self, ConstructionErrorV2> {
        let mut paths = std::collections::BTreeSet::new();
        for input in inputs {
            if !paths.insert(input.path()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_source_closure.path",
                    "one unique path per exact source entry",
                    input.path(),
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

        let expected_paths = EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| file.path)
            .collect::<std::collections::BTreeSet<_>>();
        let mut entries = Vec::with_capacity(EMBEDDED_SOURCE_FILES_V1.len());
        for (ordinal, (input, expected)) in inputs.iter().zip(EMBEDDED_SOURCE_FILES_V1).enumerate()
        {
            if input.path() != expected.path {
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
                return Err(ConstructionErrorV2::new(
                    kind,
                    "base_source_closure.path",
                    expectation,
                    format_args!("{ordinal}:{}", input.path()),
                ));
            }
            if input.bytes() != expected.bytes {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_source_closure.bytes",
                    "the exact compile-time included source bytes",
                    expected.path,
                ));
            }
            entries.push(source_closure_entry(expected.path, expected.bytes)?);
        }
        let root = source_closure_root(&entries)?;
        Ok(Self {
            entries: entries.into_boxed_slice(),
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
}

/// Result-free source-coverage classes retained by the immutable manifest.
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

/// One immutable source case. It records no execution or pass/fail result.
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

/// Immutable source-case inventory with no execution-result fields.
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
    VerifierV1 = 3,
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
        Self::VerifierV1,
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
            Self::VerifierV1 => "verifier-v1",
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
            Self::VerifierV1 => "scripts/ci/e2e_evidence_verifier_v1.sh",
            Self::CanonicalRunnerV2 => "scripts/ci/canonical_evidence_runner_v2.sh",
            Self::RjoqHandoffV1 => "scripts/ci/verify_runner_rjoq_handoff_v1.sh",
        }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
enum BaseE2eCaseKindV1 {
    CatalogLiterals = 1,
    UnknownCatalogCode = 2,
    CanonicalRational = 3,
    OverlongStableToken = 4,
    LogicalPath = 5,
    ReservedContentStorePrefix = 6,
    WindowsUnicodeAlias = 7,
    LimitCatalog = 8,
    BudgetAdmission = 9,
    BudgetChildRelation = 10,
    PublicationSelection = 11,
    PublicationCrossCell = 12,
    CapabilityLeastPrivilege = 13,
    CapabilityExtraRight = 14,
    StatePass = 15,
    StateUsageInLifecycle = 16,
    Diagnostic = 17,
    DiagnosticRankGap = 18,
    IdentityMutation = 19,
    NoClaimNominality = 20,
    AtomicResult = 21,
    AtomicResultPresence = 22,
    PublicationStorage = 23,
    CommandList = 24,
}

impl BaseE2eCaseKindV1 {
    const fn code(self) -> u16 {
        self as u16
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

/// One immutable source-closed projection row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eProjectionRowV1 {
    id: StableTokenV2,
    kind: BaseE2eCaseKindV1,
    expected: BaseE2eExpectedDecisionV1,
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
}

/// One journey-keyed immutable projection and root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eJourneyProjectionV1 {
    journey: BaseE2eJourneyV1,
    script_path: LogicalBundlePathV1,
    rows: Box<[BaseE2eProjectionRowV1]>,
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

    /// Exact row set.
    #[must_use]
    pub fn rows(&self) -> &[BaseE2eProjectionRowV1] {
        &self.rows
    }

    /// Immutable journey projection root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Complete five-journey, non-wire projection manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerV2BaseE2eProjectionV1 {
    journeys: Box<[BaseE2eJourneyProjectionV1]>,
    source_closure: RunnerV2BaseSourceClosureV1,
    coverage_inventory: BaseCoverageInventoryV1,
    root: ContentHash,
}

impl RunnerV2BaseE2eProjectionV1 {
    /// Construct the exact five journey projections and roots.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        let mut journeys = Vec::with_capacity(BaseE2eJourneyV1::ALL.len());
        for journey in BaseE2eJourneyV1::ALL {
            let script_path = LogicalBundlePathV1::new(journey.script_path()).map_err(|error| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_projection.script_path",
                    "the frozen logical relative script path",
                    format_args!("{error:?}"),
                )
            })?;
            let rows = BASE_CASE_TEMPLATES_V1
                .iter()
                .map(|template| {
                    StableTokenV2::new(template.id)
                        .map(|id| BaseE2eProjectionRowV1 {
                            id,
                            kind: template.kind,
                            expected: template.expected,
                        })
                        .map_err(|error| {
                            ConstructionErrorV2::new(
                                ConstructionErrorKindV2::Incompatible,
                                "base_e2e_projection.row_id",
                                "a frozen stable token",
                                format_args!("{error:?}"),
                            )
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let root = journey_root(journey, &script_path, &rows)?;
            journeys.push(BaseE2eJourneyProjectionV1 {
                journey,
                script_path,
                rows: rows.into_boxed_slice(),
                root,
            });
        }
        let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
        let coverage_inventory = coverage_inventory(&journeys)?;
        let root = projection_root(&journeys, source_closure.root(), coverage_inventory.root())?;
        Ok(Self {
            journeys: journeys.into_boxed_slice(),
            source_closure,
            coverage_inventory,
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
    #[must_use]
    pub const fn coverage_inventory(&self) -> &BaseCoverageInventoryV1 {
        &self.coverage_inventory
    }

    /// Complete immutable projection root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
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
    features: Box<[StableTokenV2]>,
    no_claim_scope: NoClaimScopeRootV1,
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
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_harness.features",
                    "unique stable feature tokens",
                    feature.as_str(),
                ));
            }
        }
        features.sort();
        Ok(Self {
            source,
            build,
            toolchain,
            target,
            features: features.into_boxed_slice(),
            no_claim_scope,
        })
    }
}

/// Exact aggregate projection execution counts and typed log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eProjectionReportV1 {
    eligible: u32,
    passed: u32,
    failed: u32,
    unsupported: u32,
    projection_e2e_checked: u32,
    logging_events_checked: u32,
    source_closure_eligible: u32,
    source_closure_passed: u32,
    source_closure_failed: u32,
    projection_root: ContentHash,
    source_closure_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    source_closure_paths: Box<[LogicalBundlePathV1]>,
    log: BaseE2eLogV1,
}

impl BaseE2eProjectionReportV1 {
    /// Eligible locally adjudicable rows.
    #[must_use]
    pub const fn eligible(&self) -> u32 {
        self.eligible
    }

    /// Eligible rows whose actual decision matched.
    #[must_use]
    pub const fn passed(&self) -> u32 {
        self.passed
    }

    /// Rows whose actual decision disagreed.
    #[must_use]
    pub const fn failed(&self) -> u32 {
        self.failed
    }

    /// Explicitly unsupported platform-owned cells.
    #[must_use]
    pub const fn unsupported(&self) -> u32 {
        self.unsupported
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

    /// Eligible source-closure checks executed by this report.
    #[must_use]
    pub const fn source_closure_eligible(&self) -> u32 {
        self.source_closure_eligible
    }

    /// Source-closure checks whose observed decision matched.
    #[must_use]
    pub const fn source_closure_passed(&self) -> u32 {
        self.source_closure_passed
    }

    /// Source-closure checks whose observed decision disagreed.
    #[must_use]
    pub const fn source_closure_failed(&self) -> u32 {
        self.source_closure_failed
    }

    /// Projection root executed.
    #[must_use]
    pub const fn projection_root(&self) -> ContentHash {
        self.projection_root
    }

    /// Exact source-closure root reconstructed before row execution.
    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
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
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BaseE2eCaseExecutionV1 {
    decision: BaseE2eExpectedDecisionV1,
    checked_cells: u32,
    first_failed_cell: Option<String>,
}

impl BaseE2eCaseExecutionV1 {
    fn accepted(checked_cells: u32) -> Self {
        Self {
            decision: BaseE2eExpectedDecisionV1::Accept,
            checked_cells,
            first_failed_cell: None,
        }
    }

    fn refused(checked_cells: u32) -> Self {
        Self {
            decision: BaseE2eExpectedDecisionV1::Refuse,
            checked_cells,
            first_failed_cell: None,
        }
    }

    fn unsupported(checked_cells: u32) -> Self {
        Self {
            decision: BaseE2eExpectedDecisionV1::Unsupported,
            checked_cells,
            first_failed_cell: None,
        }
    }

    fn with_failure(
        decision: BaseE2eExpectedDecisionV1,
        checked_cells: u32,
        first_failed_cell: impl Into<String>,
    ) -> Self {
        Self {
            decision,
            checked_cells,
            first_failed_cell: Some(first_failed_cell.into()),
        }
    }
}

/// Run every frozen row through real in-process public constructors and
/// validators with deterministic detailed logging.
pub fn run_base_e2e_projection_v1(
    projection: &RunnerV2BaseE2eProjectionV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<BaseE2eProjectionReportV1, ConstructionErrorV2> {
    let expected = RunnerV2BaseE2eProjectionV1::frozen()?;
    if projection != &expected {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_projection.root",
            "the exact immutable frozen projection",
            projection.root.to_hex(),
        ));
    }

    let (source_closure_eligible, source_closure_passed, source_closure_failed) =
        run_source_closure_checks(&projection.source_closure);
    let mut sequence = 0_u32;
    let mut eligible = 0_u32;
    let mut passed = 0_u32;
    let mut failed = 0_u32;
    let mut unsupported = 0_u32;
    let mut projection_e2e_checked = 0_u32;
    let mut events = Vec::new();
    for journey in &projection.journeys {
        events.push(log_event(
            sequence,
            journey,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            harness,
            Vec::new(),
        )?);
        sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
        let mut journey_eligible = 0_u32;
        let mut journey_passed = 0_u32;
        let mut journey_failed = 0_u32;
        let mut journey_unsupported = 0_u32;
        for row in &journey.rows {
            let execution = execute_case(row.kind, harness);
            let actual = execution.decision;
            projection_e2e_checked = projection_e2e_checked
                .checked_add(execution.checked_cells)
                .ok_or_else(sequence_overflow)?;
            let agrees = actual == row.expected;
            let outcome = if row.expected == BaseE2eExpectedDecisionV1::Unsupported && agrees {
                unsupported = unsupported.checked_add(1).ok_or_else(sequence_overflow)?;
                journey_unsupported = journey_unsupported
                    .checked_add(1)
                    .ok_or_else(sequence_overflow)?;
                BaseE2eOutcomeV1::Unsupported
            } else {
                eligible = eligible.checked_add(1).ok_or_else(sequence_overflow)?;
                journey_eligible = journey_eligible
                    .checked_add(1)
                    .ok_or_else(sequence_overflow)?;
                if agrees {
                    passed = passed.checked_add(1).ok_or_else(sequence_overflow)?;
                    journey_passed = journey_passed
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                    BaseE2eOutcomeV1::Passed
                } else {
                    failed = failed.checked_add(1).ok_or_else(sequence_overflow)?;
                    journey_failed = journey_failed
                        .checked_add(1)
                        .ok_or_else(sequence_overflow)?;
                    BaseE2eOutcomeV1::Failed
                }
            };
            let mut case_fields = vec![
                field("checked-cells", TypedValueV2::U32(execution.checked_cells))?,
                field("expected", TypedValueV2::Token(token(row.expected.name())?))?,
                field("observed", TypedValueV2::Token(token(actual.name())?))?,
            ];
            if let Some(first_failed_cell) = execution.first_failed_cell {
                case_fields.push(field(
                    "first-failed-cell",
                    TypedValueV2::Token(token(&first_failed_cell)?),
                )?);
            }
            case_fields.extend(case_detail_fields(row.kind, harness)?);
            events.push(log_event(
                sequence,
                journey,
                Some(row),
                BaseE2eLogKindV1::CaseTerminal,
                outcome,
                harness,
                case_fields,
            )?);
            sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
        }
        events.push(log_event(
            sequence,
            journey,
            None,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eOutcomeV1::NotApplicable,
            harness,
            count_fields(
                journey_eligible,
                journey_passed,
                journey_failed,
                journey_unsupported,
            )?,
        )?);
        sequence = sequence.checked_add(1).ok_or_else(sequence_overflow)?;
    }
    let summary_journey = &projection.journeys[0];
    let logging_events_checked = u32::try_from(events.len() + 1).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_projection.logging_events_checked",
            "a u32 deterministic event count",
            events.len() + 1,
        )
    })?;
    let mut summary_fields = count_fields(eligible, passed, failed, unsupported)?;
    summary_fields.extend([
        field(
            "coverage-source-cases",
            TypedValueV2::U32(
                u32::try_from(projection.coverage_inventory.cases.len())
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
            TypedValueV2::U32(source_closure_eligible),
        )?,
        field(
            "source-closure-failed",
            TypedValueV2::U32(source_closure_failed),
        )?,
        field(
            "source-closure-passed",
            TypedValueV2::U32(source_closure_passed),
        )?,
        field(
            "source-closure-root",
            opaque_root(projection.source_closure.root())?,
        )?,
    ]);
    events.push(log_event(
        sequence,
        summary_journey,
        None,
        BaseE2eLogKindV1::ProjectionSummary,
        BaseE2eOutcomeV1::NotApplicable,
        harness,
        summary_fields,
    )?);
    let log = BaseE2eLogV1::new(events)?;
    let source_closure_paths = projection
        .source_closure
        .entries()
        .iter()
        .map(|entry| {
            LogicalBundlePathV1::new(entry.path()).map_err(|error| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_projection.source_closure_path",
                    "a source-closure-relative path",
                    format_args!("{error:?}"),
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_boxed_slice();
    Ok(BaseE2eProjectionReportV1 {
        eligible,
        passed,
        failed,
        unsupported,
        projection_e2e_checked,
        logging_events_checked,
        source_closure_eligible,
        source_closure_passed,
        source_closure_failed,
        projection_root: projection.root,
        source_closure_root: projection.source_closure.root(),
        source_root: harness.source.clone(),
        build_root: harness.build.clone(),
        source_closure_paths,
        log,
    })
}

fn execute_case(
    kind: BaseE2eCaseKindV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> BaseE2eCaseExecutionV1 {
    match kind {
        BaseE2eCaseKindV1::CatalogLiterals => aggregate_accept(catalog_literal_matrix()),
        BaseE2eCaseKindV1::UnknownCatalogCode => refuse_if(
            ProofExitV2::from_code(1).is_err(),
            "catalog.unknown-code",
        ),
        BaseE2eCaseKindV1::CanonicalRational => accept_if(
            RationalV2::new(6, 8).ok() == RationalV2::new(3, 4).ok(),
            "value.rational-equivalence",
        ),
        BaseE2eCaseKindV1::OverlongStableToken => refuse_if(
            StableTokenV2::new("a".repeat(129)).is_err(),
            "value.overlong-token",
        ),
        BaseE2eCaseKindV1::LogicalPath => accept_if(
            LogicalBundlePathV1::new("runner/seal").is_ok(),
            "path.logical",
        ),
        BaseE2eCaseKindV1::ReservedContentStorePrefix => refuse_if(
            ContentStoreObjectKeyV1::new("__runner_private/object").is_err(),
            "path.reserved-prefix",
        ),
        BaseE2eCaseKindV1::WindowsUnicodeAlias => {
            let paths = [LogicalBundlePathV1::new("résumé/a").expect("valid UTF-8 path")];
            match adjudicate_logical_bundle_path_set(
                PlatformPathProfileV2::WindowsHandleRelativeV1,
                &paths,
            ) {
                PathSetAdjudicationV1::UnsupportedWindowsNonAsciiAlias { .. } => {
                    BaseE2eCaseExecutionV1::unsupported(1)
                }
                _ => BaseE2eCaseExecutionV1::with_failure(
                    BaseE2eExpectedDecisionV1::Refuse,
                    1,
                    "path.windows-unicode-alias",
                ),
            }
        }
        BaseE2eCaseKindV1::LimitCatalog => aggregate_accept(limit_matrix()),
        BaseE2eCaseKindV1::BudgetAdmission => aggregate_accept(budget_matrix()),
        BaseE2eCaseKindV1::BudgetChildRelation => {
            let mut candidate = durable_budget_candidate();
            candidate.max_parallel_children = candidate.max_child_processes + 1;
            refuse_if(
                RunnerBudgetsV2::try_new(candidate).is_err(),
                "budget.parallel-children",
            )
        }
        BaseE2eCaseKindV1::PublicationSelection => {
            aggregate_accept(publication_selection_matrix())
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
                .is_err(),
                "publication.cross-cell",
            )
        }
        BaseE2eCaseKindV1::CapabilityLeastPrivilege => {
            aggregate_accept(capability_valid_matrix(harness.no_claim_scope.clone()))
        }
        BaseE2eCaseKindV1::CapabilityExtraRight => {
            aggregate_refusal(capability_invalid_matrix(harness.no_claim_scope.clone()))
        }
        BaseE2eCaseKindV1::StatePass => aggregate_accept(state_and_not_run_matrix()),
        BaseE2eCaseKindV1::StateUsageInLifecycle => refuse_if(
            validate_state_v2(StateValidationInputV2::new(
                StateBearingRecordRoleV2::ExecutedCaseTerminal,
                ProofExitV2::Usage,
                None,
                Some(DiagnosticCodeV2::RunnerUsage),
                None,
            ))
            .is_err(),
            "state.usage-in-lifecycle",
        ),
        BaseE2eCaseKindV1::Diagnostic => {
            aggregate_accept(diagnostic_matrix(harness.no_claim_scope.clone()))
        }
        BaseE2eCaseKindV1::DiagnosticRankGap => refuse_if(
            diagnostic(harness.no_claim_scope.clone(), 2).is_err(),
            "diagnostic.rank-gap",
        ),
        BaseE2eCaseKindV1::IdentityMutation => {
            aggregate_accept(identity_mutation_matrix(harness))
        }
        BaseE2eCaseKindV1::NoClaimNominality => {
            aggregate_accept(no_claim_matrix(&harness.no_claim_scope))
        }
        BaseE2eCaseKindV1::AtomicResult => accept_if(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::List,
                32,
                0,
                128,
                0,
                1024,
                1024,
            )
            .is_ok(),
            "result.atomic",
        ),
        BaseE2eCaseKindV1::AtomicResultPresence => refuse_if(
            SymbolicCommandResultPlanV2::new(
                RunnerCommandV2::Run,
                32,
                128,
                1,
                128,
                1024,
                1024,
            )
            .is_err(),
            "result.atomic-presence",
        ),
        BaseE2eCaseKindV1::PublicationStorage => {
            accept_if(publication_storage().is_ok(), "publication.storage")
        }
        BaseE2eCaseKindV1::CommandList => aggregate_accept(command_matrix()),
    }
}

fn accept_if(condition: bool, failed_cell: &'static str) -> BaseE2eCaseExecutionV1 {
    if condition {
        BaseE2eCaseExecutionV1::accepted(1)
    } else {
        BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Refuse,
            1,
            failed_cell,
        )
    }
}

fn refuse_if(condition: bool, failed_cell: &'static str) -> BaseE2eCaseExecutionV1 {
    if condition {
        BaseE2eCaseExecutionV1::refused(1)
    } else {
        BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Accept,
            1,
            failed_cell,
        )
    }
}

fn aggregate_accept(result: Result<u32, (u32, String)>) -> BaseE2eCaseExecutionV1 {
    match result {
        Ok(checked_cells) => BaseE2eCaseExecutionV1::accepted(checked_cells),
        Err((checked_cells, failed_cell)) => BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Refuse,
            checked_cells,
            failed_cell,
        ),
    }
}

fn aggregate_refusal(result: Result<u32, (u32, String)>) -> BaseE2eCaseExecutionV1 {
    match result {
        Ok(checked_cells) => BaseE2eCaseExecutionV1::refused(checked_cells),
        Err((checked_cells, failed_cell)) => BaseE2eCaseExecutionV1::with_failure(
            BaseE2eExpectedDecisionV1::Accept,
            checked_cells,
            failed_cell,
        ),
    }
}

fn case_detail_fields(
    kind: BaseE2eCaseKindV1,
    harness: &BaseE2eHarnessIdentityV1,
) -> Result<Vec<BaseE2eLogFieldV1>, ConstructionErrorV2> {
    let fields = match kind {
        BaseE2eCaseKindV1::CatalogLiterals => vec![field(
            "catalog-literal-cells",
            TypedValueV2::U32(184),
        )?],
        BaseE2eCaseKindV1::LimitCatalog => vec![
            field("limit-field-count", TypedValueV2::U32(65))?,
            field("limit-profile-cells", TypedValueV2::U32(130))?,
        ],
        BaseE2eCaseKindV1::BudgetAdmission => vec![
            field("budget-field-count", TypedValueV2::U32(18))?,
            field("logical-unit-count", TypedValueV2::U32(16))?,
        ],
        BaseE2eCaseKindV1::CapabilityLeastPrivilege
        | BaseE2eCaseKindV1::CapabilityExtraRight => vec![
            field("capability-valid-cells", TypedValueV2::U32(12))?,
            field("capability-mutant-cells", TypedValueV2::U32(390))?,
            field("capability-right-count", TypedValueV2::U32(10))?,
        ],
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
            field(
                "diagnostic-expected",
                TypedValueV2::U64(4),
            )?,
            field(
                "diagnostic-observed",
                TypedValueV2::U64(5),
            )?,
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
        BaseE2eCaseKindV1::NoClaimNominality => vec![field(
            "no-claim-scope",
            TypedValueV2::Digest(harness.no_claim_scope.digest().clone()),
        )?],
        _ => Vec::new(),
    };
    Ok(fields)
}

fn catalog_literal_matrix() -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;

    macro_rules! check_closed_catalog {
        ($catalog:ty, $label:literal) => {
            for value in <$catalog>::ALL {
                checked += 1;
                if <$catalog>::from_code(value.code()) != Ok(value) {
                    return Err((checked, format!("catalog.{}.{}", $label, value.code())));
                }
            }
        };
    }

    macro_rules! check_registered_catalog {
        ($catalog:ty, $label:literal) => {
            for descriptor in <$catalog>::ALL {
                checked += 1;
                let registered_id = descriptor.requires_registered_id().then_some(7);
                let value = match <$catalog>::from_tag(descriptor.tag(), registered_id) {
                    Ok(value) => value,
                    Err(_) => {
                        return Err((
                            checked,
                            format!("catalog.{}.{}", $label, descriptor.tag()),
                        ));
                    }
                };
                if value.tag() != descriptor.tag()
                    || value.name() != descriptor.name()
                    || value.registered_id() != registered_id
                {
                    return Err((
                        checked,
                        format!("catalog.{}.{}", $label, descriptor.tag()),
                    ));
                }
            }
        };
    }

    check_closed_catalog!(RunnerApiGeneration, "api-generation");
    check_closed_catalog!(RunnerWireVersion, "wire-version");
    checked += 1;
    if WirePredecessorPolicyV1::ALL != [WirePredecessorPolicyV1::NoPredecessor]
        || WirePredecessorPolicyV1::NoPredecessor.predecessor().is_some()
    {
        return Err((checked, "catalog.wire-predecessor.1".to_owned()));
    }
    check_closed_catalog!(ProofExitV2, "proof-exit");
    check_closed_catalog!(RefusedReasonV2, "refused-reason");
    check_closed_catalog!(RunnerCommandV2, "runner-command");
    check_closed_catalog!(RunProfileV2, "run-profile");
    check_closed_catalog!(ArtifactDispositionV2, "artifact-disposition");
    check_closed_catalog!(PlatformPathProfileV2, "path-profile");
    check_closed_catalog!(LifecycleRecordKindV2, "record-kind");
    check_closed_catalog!(StateBearingRecordRoleV2, "record-role");
    check_closed_catalog!(DiagnosticCodeV2, "diagnostic-code");
    check_closed_catalog!(RetryabilityV2, "retryability");
    check_closed_catalog!(RepairActionKindV2, "repair-kind");
    check_closed_catalog!(NotRunCauseCodeV2, "not-run-cause");
    check_closed_catalog!(TypedValueTagV2, "typed-value");
    check_closed_catalog!(TypedOptionTagV1, "typed-option");
    check_closed_catalog!(DigestRoleV2, "digest-role");
    check_closed_catalog!(PublicationProtocolV2, "publication-protocol");
    check_closed_catalog!(DestinationAdmissionModeV2, "destination-mode");
    check_closed_catalog!(RootCapabilityAccessV2, "capability-access");
    check_closed_catalog!(RootCapabilityRightV2, "capability-right");
    check_closed_catalog!(OverlapPolicyRelationV2, "overlap-relation");
    check_registered_catalog!(RootClassV2, "root-class");
    check_registered_catalog!(LogicalUnitV2, "logical-unit");
    check_registered_catalog!(ArtifactRoleV2, "artifact-role");
    check_registered_catalog!(LogicalExtentAxisV2, "logical-axis");

    if checked != 184 {
        return Err((checked, "catalog.total-count".to_owned()));
    }
    Ok(checked)
}

fn limit_matrix() -> Result<u32, (u32, String)> {
    if RUNNER_LIMIT_DESCRIPTORS_V2.len() != 65 {
        return Err((0, "limit.descriptor-count".to_owned()));
    }
    let mut checked = 0_u32;
    for profile in RunProfileV2::ALL {
        let admitted = RunnerLimitsV2::base(profile);
        for field in RunnerLimitFieldV2::ALL {
            checked += 1;
            let descriptor = field.descriptor();
            let failure = || format!("limit.{}.{}", profile.name(), descriptor.name);
            if descriptor.field != field
                || descriptor.ordinal != field.ordinal()
                || RunnerLimitFieldV2::from_ordinal(descriptor.ordinal) != Some(field)
                || admitted.value(field).width() != descriptor.width
            {
                return Err((checked, failure()));
            }

            let mut one_over = admitted.to_candidate();
            let one_over_value = match one_over.value(field) {
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
            if one_over.set_value(field, one_over_value).is_err() {
                return Err((checked, failure()));
            }
            let error = match RunnerLimitsV2::admit_family(
                profile,
                one_over,
                RunnerFamilyLimitRequirementsV2::NONE,
            ) {
                Ok(_) => return Err((checked, failure())),
                Err(error) => error,
            };
            let expected_kind = match descriptor.tightenability {
                RunnerLimitTightenabilityV2::Fixed => {
                    RunnerLimitsViolationKindV2::FixedFieldChanged
                }
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling
                }
            };
            if error.kind() != expected_kind || error.field() != field {
                return Err((checked, failure()));
            }
        }
    }
    if checked != 130 {
        return Err((checked, "limit.total-count".to_owned()));
    }
    Ok(checked)
}

fn budget_matrix() -> Result<u32, (u32, String)> {
    let base = RunnerBudgetsV2::try_new(durable_budget_candidate())
        .map_err(|_| (0, "budget.base-construction".to_owned()))?;
    let base_root = base.semantic_root();
    let mut checked = 0_u32;

    for field in RunnerBudgetFieldV2::ALL {
        checked += 1;
        let descriptor = field.descriptor();
        let failure = || format!("budget.field.{}", descriptor.name);
        if descriptor.field != field
            || descriptor.ordinal != field.ordinal()
            || RunnerBudgetFieldV2::from_ordinal(descriptor.ordinal) != Some(field)
        {
            return Err((checked, failure()));
        }
        let mutated = RunnerBudgetsV2::try_new(mutated_budget_candidate(field))
            .map_err(|_| (checked, failure()))?;
        if mutated.value(field) == base.value(field)
            || mutated.semantic_root().bytes() == base_root.bytes()
        {
            return Err((checked, failure()));
        }
    }

    for descriptor in LogicalUnitV2::ALL {
        checked += 1;
        let unit = LogicalUnitV2::from_tag(
            descriptor.tag(),
            descriptor.requires_registered_id().then_some(7),
        )
        .map_err(|_| {
            (
                checked,
                format!("budget.logical-unit.{}", descriptor.tag()),
            )
        })?;
        let mut candidate = durable_budget_candidate();
        candidate.logical_work_unit = unit;
        let value = RunnerBudgetsV2::try_new(candidate).map_err(|_| {
            (
                checked,
                format!("budget.logical-unit.{}", descriptor.tag()),
            )
        })?;
        if value.logical_work_unit() != unit {
            return Err((
                checked,
                format!("budget.logical-unit.{}", descriptor.tag()),
            ));
        }
    }

    for profile in RunProfileV2::ALL {
        checked += 1;
        if base
            .admit(
                profile,
                ArtifactDispositionV2::DurableBundleRequired,
                &RunnerLimitsV2::base(profile),
            )
            .is_err()
        {
            return Err((checked, format!("budget.admission.{}", profile.name())));
        }
    }

    if checked != 36 {
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

fn admitted_durable_budget() -> Result<(), ()> {
    let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
    RunnerBudgetsV2::try_new(durable_budget_candidate())
        .map_err(|_| ())?
        .admit(
            RunProfileV2::Smoke,
            ArtifactDispositionV2::DurableBundleRequired,
            &limits,
        )
        .map(|_| ())
        .map_err(|_| ())
}

fn publication_selection_matrix() -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    let mut roots = Vec::new();
    for profile in PlatformPathProfileV2::ALL {
        for mode in DestinationAdmissionModeV2::ALL {
            checked += 1;
            let selection = selection_for_profile(profile, mode)
                .map_err(|_| (checked, format!("publication.{}.{}", profile.name(), mode.name())))?;
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
    let logical = LogicalBundlePathV1::new("runner/seal").map_err(|error| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.publication_path",
            "a valid logical path",
            format_args!("{error:?}"),
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
                ContentStoreObjectKeyV1::new("objects/seal").map_err(|error| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "projection.content_store_key",
                        "a valid logical object key",
                        format_args!("{error:?}"),
                    )
                })?,
            ),
        ),
    };
    PublicationSelectionV2::new(profile, protocol, mode, target)
}

fn publication_selection(path: &str) -> Result<PublicationSelectionV2, ConstructionErrorV2> {
    let path = LogicalBundlePathV1::new(path).map_err(|error| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.publication_path",
            "a valid logical path",
            format_args!("{error:?}"),
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

fn capability_valid_matrix(
    no_claim_scope: NoClaimScopeRootV1,
) -> Result<u32, (u32, String)> {
    let registry =
        capability_registry().map_err(|_| (0, "capability.registry".to_owned()))?;
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
                let policy = RootCapabilityPolicyV2::new(
                    root_class_for_access(access),
                    profile,
                    access,
                    expected_rights(profile, access, mode),
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

fn capability_invalid_matrix(
    no_claim_scope: NoClaimScopeRootV1,
) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    for profile in PlatformPathProfileV2::ALL {
        for access in RootCapabilityAccessV2::ALL {
            for mode in DestinationAdmissionModeV2::ALL {
                let exact = expected_rights(profile, access, mode);
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
        Err(_) => return true,
    };
    match access {
        RootCapabilityAccessV2::ReadOnlyInput => {
            NarrowedPolicyViewV2::for_read_only(&policy).is_err()
        }
        RootCapabilityAccessV2::DurableOutput => selection_for_profile(profile, mode)
            .and_then(|selection| NarrowedPolicyViewV2::for_publication(&policy, &selection))
            .is_err(),
    }
}

const fn root_class_for_access(access: RootCapabilityAccessV2) -> RootClassV2 {
    match access {
        RootCapabilityAccessV2::ReadOnlyInput => RootClassV2::InputArtifactRoot,
        RootCapabilityAccessV2::DurableOutput => RootClassV2::OutputArtifactRoot,
    }
}

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
                        .is_ok();
                        let expected =
                            expected_state_cell(role, state, *reason, *diagnostic, drain);
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

    let (cancelled, timed_out, internal_error) = presented_stop_fixture_roots()
        .map_err(|_| (checked, "not-run.fixture-roots".to_owned()))?;
    let causes = [
        NotRunCauseV2::PriorCancelled(cancelled),
        NotRunCauseV2::PriorTimedOut(timed_out),
        NotRunCauseV2::PriorControlledInternalError(internal_error),
    ];
    for cause in causes {
        let code = cause.code();
        checked += 1;
        let first = NotRunBasisV2::new(cause.clone(), 0, 1)
            .map_err(|_| (checked, format!("not-run.{}.first", code)))?;
        if first.remaining_case_count(1) != Ok(1)
            || first.diagnostic() != DiagnosticCodeV2::RunnerNotRun
            || first.state() != ProofExitV2::NotRun
        {
            return Err((checked, format!("not-run.{}.first", code)));
        }

        checked += 1;
        let last = NotRunBasisV2::new(cause.clone(), 255, 256)
            .map_err(|_| (checked, format!("not-run.{}.last", code)))?;
        if last.remaining_case_count(256) != Ok(1) {
            return Err((checked, format!("not-run.{}.last", code)));
        }

        checked += 1;
        if NotRunBasisV2::new(cause.clone(), 256, 256).is_ok() {
            return Err((checked, format!("not-run.{}.one-over", code)));
        }

        checked += 1;
        if NotRunBasisV2::new(cause, 0, 0).is_ok() {
            return Err((checked, format!("not-run.{}.empty", code)));
        }
    }

    if checked != 32_460 {
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
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            field,
            "a nominal presented stop fixture root",
            "fixture construction failed",
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

fn diagnostic_matrix(
    no_claim_scope: NoClaimScopeRootV1,
) -> Result<u32, (u32, String)> {
    let mut checked = 0_u32;
    for code in DiagnosticCodeV2::ALL {
        checked += 1;
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(code),
            RetryabilityV2::AfterInputChange,
            RepairActionKindV2::ChangeArguments,
            1,
        )
        .map_err(|_| (checked, format!("diagnostic.code.{}", code.code())))?;
        if value.code().code() != code.code() {
            return Err((checked, format!("diagnostic.code.{}", code.code())));
        }
    }

    checked += 1;
    let registered = DiagnosticCodeRefV2::registered(7, 9)
        .map_err(|_| (checked, "diagnostic.registered".to_owned()))?;
    let registered_value = diagnostic_fixture(
        no_claim_scope.clone(),
        registered,
        RetryabilityV2::AfterPrerequisiteChange,
        RepairActionKindV2::ContactOwner,
        1,
    )
    .map_err(|_| (checked, "diagnostic.registered".to_owned()))?;
    if registered_value.code().registered_namespace() != Some(7)
        || registered_value.code().code() != 9
    {
        return Err((checked, "diagnostic.registered".to_owned()));
    }

    for retryability in RetryabilityV2::ALL {
        checked += 1;
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            retryability,
            RepairActionKindV2::ChangeArguments,
            1,
        )
        .map_err(|_| {
            (
                checked,
                format!("diagnostic.retryability.{}", retryability.code()),
            )
        })?;
        if value.retryability() != retryability {
            return Err((
                checked,
                format!("diagnostic.retryability.{}", retryability.code()),
            ));
        }
    }

    for kind in RepairActionKindV2::ALL {
        checked += 1;
        let value = diagnostic_fixture(
            no_claim_scope.clone(),
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            RetryabilityV2::AfterInputChange,
            kind,
            1,
        )
        .map_err(|_| (checked, format!("diagnostic.repair-kind.{}", kind.code())))?;
        if value.repairs()[0].kind() != kind {
            return Err((checked, format!("diagnostic.repair-kind.{}", kind.code())));
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

fn publication_storage() -> Result<(), ()> {
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
    limits
        .validate_publication_storage(PublicationStorageProjectionV2 {
            artifacts: &artifacts,
            system_objects: &system_objects,
            artifact_encoded_bytes: 1,
            artifact_stored_bytes: 1,
            system_publication_stored_bytes: 6,
            publication_stored_bytes: 7,
        })
        .map_err(|_| ())
}

fn source_closure_entry(
    path: &'static str,
    bytes: &[u8],
) -> Result<BaseSourceClosureEntryV1, ConstructionErrorV2> {
    let encoded_bytes = u64::try_from(bytes.len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_source_closure.encoded_bytes",
            "a u64 source byte length",
            bytes.len(),
        )
    })?;
    let content_root = hash_domain(BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, bytes);
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCEENTRY\x01", 2048)?;
    frame.push_str("source.path", path)?;
    frame.push_u64("source.encoded_bytes", encoded_bytes)?;
    frame.push_bytes("source.content_root", content_root.as_bytes())?;
    let entry_root = frame.root(BASE_SOURCE_FILE_ENTRY_DOMAIN_V1);
    Ok(BaseSourceClosureEntryV1 {
        path,
        encoded_bytes,
        content_root,
        entry_root,
    })
}

fn source_closure_root(
    entries: &[BaseSourceClosureEntryV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASESOURCECLOSURE\x01", 16 * 1024)?;
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
        frame.push_u64("source.encoded_bytes", entry.encoded_bytes())?;
        frame.push_bytes("source.content_root", entry.content_root().as_bytes())?;
        frame.push_bytes("source.entry_root", entry.entry_root().as_bytes())?;
    }
    Ok(frame.root(BASE_SOURCE_CLOSURE_DOMAIN_V1))
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
                class.code(),
            ));
        }
    }

    let expected_e2e_rows = BaseE2eJourneyV1::ALL
        .len()
        .checked_mul(BASE_CASE_TEMPLATES_V1.len())
        .ok_or_else(sequence_overflow)?;
    let actual_e2e_rows = cases
        .iter()
        .filter(|source_case| source_case.class == BaseCoverageClassV1::ProjectionE2e)
        .count();
    if actual_e2e_rows != expected_e2e_rows || actual_e2e_rows != 120 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "coverage.projection_e2e_source_cases",
            "exactly five journeys times twenty-four rows, or 120 source cases",
            actual_e2e_rows,
        ));
    }

    let mut seen = std::collections::BTreeSet::new();
    for source_case in &cases {
        if !seen.insert((source_case.class.code(), source_case.id.as_str())) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                "coverage.source_case_id",
                "unique IDs within every coverage class",
                source_case.id.as_str(),
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
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "coverage.source_path",
            "a member of the exact embedded source closure",
            source_path,
        ));
    }
    let id = token(id)?;
    let source_path = LogicalBundlePathV1::new(source_path).map_err(|error| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "coverage.source_path",
            "a validated relative source path",
            format_args!("{error:?}"),
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
    script_path: &LogicalBundlePathV1,
    rows: &[BaseE2eProjectionRowV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSBASEJOURNEY\x01", 64 * 1024)?;
    frame.push_u16("projection.journey", journey.code())?;
    frame.push_str("projection.script_path", script_path.as_str())?;
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
        frame.push_str("projection.row_id", row.id.as_str())?;
        frame.push_u16("projection.case_kind", row.kind.code())?;
        frame.push_u16("projection.expected", row.expected.code())?;
    }
    Ok(frame.root(BASE_E2E_JOURNEY_PROJECTION_DOMAIN_V1))
}

fn projection_root(
    journeys: &[BaseE2eJourneyProjectionV1],
    source_closure_root: ContentHash,
    coverage_inventory_root: ContentHash,
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
    Ok(frame.root(BASE_E2E_PROJECTION_DOMAIN_V1))
}

fn log_event(
    sequence: u32,
    journey: &BaseE2eJourneyProjectionV1,
    row: Option<&BaseE2eProjectionRowV1>,
    kind: BaseE2eLogKindV1,
    outcome: BaseE2eOutcomeV1,
    harness: &BaseE2eHarnessIdentityV1,
    mut fields: Vec<BaseE2eLogFieldV1>,
) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
    fields.extend([
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
        field(
            "projection-root",
            TypedValueV2::OpaqueBytes(
                crate::value::OpaqueBytesV2::new(journey.root.as_bytes().to_vec())
                    .expect("32-byte root fits opaque value"),
            ),
        )?,
    ]);
    BaseE2eLogEventV1::new(
        sequence,
        token(journey.journey.key())?,
        row.map(|row| row.id.clone()),
        kind,
        outcome,
        fields,
        Some(journey.script_path.clone()),
        vec![
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(token(journey.journey.key())?),
        ],
    )
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

fn field(
    name: &'static str,
    value: TypedValueV2,
) -> Result<BaseE2eLogFieldV1, ConstructionErrorV2> {
    Ok(BaseE2eLogFieldV1::new(token(name)?, value))
}

fn token(value: &str) -> Result<StableTokenV2, ConstructionErrorV2> {
    StableTokenV2::new(value).map_err(|error| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "projection.stable_token",
            "a valid stable token",
            format_args!("{error:?}"),
        )
    })
}

fn sequence_overflow() -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::ArithmeticOverflow,
        "projection.logical_sequence",
        "a checked u32 sequence",
        "overflow",
    )
}

#[cfg(test)]
mod tests {
    use super::{
        BASE_CASE_TEMPLATES_V1, BASE_SOURCE_FILE_CONTENT_DOMAIN_V1, BaseCoverageClassV1,
        BaseE2eHarnessIdentityV1, BaseE2eJourneyV1, BaseSourceClosureInputV1,
        EMBEDDED_SOURCE_FILES_V1, RunnerV2BaseE2eProjectionV1, RunnerV2BaseSourceClosureV1,
        run_base_e2e_projection_v1,
    };
    use crate::catalog::DigestRoleV2;
    use crate::construction::ConstructionErrorKindV2;
    use crate::identity::{
        BuildIdentityRootV2, NoClaimScopeRootV1, SourceIdentityRootV2, ToolchainIdentityRootV2,
    };
    use crate::value::StableTokenV2;
    use fs_blake3::hash_domain;

    const EXPECTED_SOURCE_PATHS_V1: [&str; 21] = [
        "Cargo.lock",
        "Cargo.toml",
        "crates/fs-evidence-runner/CONTRACT.md",
        "crates/fs-evidence-runner/Cargo.toml",
        "crates/fs-evidence-runner/src/budget.rs",
        "crates/fs-evidence-runner/src/canonical.rs",
        "crates/fs-evidence-runner/src/capability.rs",
        "crates/fs-evidence-runner/src/catalog.rs",
        "crates/fs-evidence-runner/src/command.rs",
        "crates/fs-evidence-runner/src/construction.rs",
        "crates/fs-evidence-runner/src/dependency.rs",
        "crates/fs-evidence-runner/src/diagnostic.rs",
        "crates/fs-evidence-runner/src/identity.rs",
        "crates/fs-evidence-runner/src/lib.rs",
        "crates/fs-evidence-runner/src/limits.rs",
        "crates/fs-evidence-runner/src/logging.rs",
        "crates/fs-evidence-runner/src/path.rs",
        "crates/fs-evidence-runner/src/projection.rs",
        "crates/fs-evidence-runner/src/publication.rs",
        "crates/fs-evidence-runner/src/state.rs",
        "crates/fs-evidence-runner/src/value.rs",
    ];

    fn presented<T>(
        role: DigestRoleV2,
        domain: &str,
        parser: impl FnOnce(DigestRoleV2, &str, &str) -> Result<T, crate::identity::IdentityError>,
    ) -> T {
        parser(role, domain, &"00".repeat(32)).expect("presented fixture")
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

    fn frozen_source_inputs() -> Vec<BaseSourceClosureInputV1> {
        EMBEDDED_SOURCE_FILES_V1
            .iter()
            .map(|file| BaseSourceClosureInputV1::presented(file.path, file.bytes.to_vec()))
            .collect()
    }

    #[test]
    fn manifest_exactly_maps_five_scripts_and_all_base_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        assert_eq!(projection.journeys().len(), BaseE2eJourneyV1::ALL.len());
        let mut row_count = 0_usize;
        for (index, journey) in projection.journeys().iter().enumerate() {
            assert_eq!(journey.journey(), BaseE2eJourneyV1::ALL[index]);
            assert_eq!(
                journey.script_path().as_str(),
                BaseE2eJourneyV1::ALL[index].script_path()
            );
            assert_eq!(journey.rows().len(), BASE_CASE_TEMPLATES_V1.len());
            row_count += journey.rows().len();
        }
        assert_eq!(row_count, 120);
    }

    #[test]
    fn all_real_constructor_rows_agree_and_logs_are_deterministic() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let first = run_base_e2e_projection_v1(&projection, &harness()).expect("projection run");
        let second = run_base_e2e_projection_v1(&projection, &harness()).expect("projection run");
        assert_eq!(first, second);
        assert_eq!(first.failed(), 0);
        assert!(first.eligible() > 0);
        assert!(first.passed() > 0);
        assert!(first.unsupported() > 0);
        assert!(!first.log().events().is_empty());
    }

    #[test]
    fn one_field_projection_mutation_moves_the_root() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let first = projection.journeys()[0].root();
        let second = projection.journeys()[1].root();
        assert_ne!(first, second);
    }

    #[test]
    fn source_closure_membership_and_order_are_exact_and_content_bound() {
        let closure = RunnerV2BaseSourceClosureV1::frozen().expect("frozen source closure");
        assert_eq!(closure.entries().len(), EXPECTED_SOURCE_PATHS_V1.len());
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
        }
        assert!(
            EXPECTED_SOURCE_PATHS_V1
                .windows(2)
                .all(|pair| pair[0].as_bytes() < pair[1].as_bytes())
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

    #[test]
    fn coverage_inventory_is_result_free_nonzero_and_exactly_maps_120_e2e_rows() {
        let projection = RunnerV2BaseE2eProjectionV1::frozen().expect("frozen projection");
        let inventory = projection.coverage_inventory();
        let expected_counts = [
            (BaseCoverageClassV1::Unit, 14),
            (BaseCoverageClassV1::Boundary, 5),
            (BaseCoverageClassV1::PropertyMetamorphic, 4),
            (BaseCoverageClassV1::CompileFailDoctest, 2),
            (BaseCoverageClassV1::SchemaDescriptor, 4),
            (BaseCoverageClassV1::Mutation, 6),
            (BaseCoverageClassV1::Integration, 4),
            (BaseCoverageClassV1::ProjectionE2e, 120),
            (BaseCoverageClassV1::Logging, 3),
            (BaseCoverageClassV1::SourceClosure, 2),
        ];
        assert_eq!(inventory.cases().len(), 164);
        for (class, expected_count) in expected_counts {
            assert_eq!(inventory.source_case_count(class), expected_count);
            assert!(inventory.source_case_count(class) > 0);
            let ordinals = inventory
                .cases()
                .iter()
                .filter(|source_case| source_case.class() == class)
                .map(|source_case| source_case.ordinal())
                .collect::<Vec<_>>();
            assert_eq!(
                ordinals,
                (1..=u32::try_from(expected_count).expect("fixture count fits u32"))
                    .collect::<Vec<_>>()
            );
        }
        let closure_paths = projection
            .source_closure()
            .entries()
            .iter()
            .map(|entry| entry.path())
            .collect::<std::collections::BTreeSet<_>>();
        assert!(
            inventory
                .cases()
                .iter()
                .all(|source_case| { closure_paths.contains(source_case.source_path().as_str()) })
        );
    }
}
