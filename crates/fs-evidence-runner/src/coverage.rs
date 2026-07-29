//! Exact, result-free coverage declarations and checked result accounting.
//!
//! This module deliberately separates three things that are easy to conflate:
//!
//! 1. the immutable source-case manifest;
//! 2. a caller-selected, manifest-ordered executable subset; and
//! 3. presented result records joined back to that exact subset.
//!
//! The frozen base enumerates the current 160 non-manifest Rust tests and all
//! 47
//! `compile_fail` contracts. The historical aggregate is recorded as 130
//! ratified cases plus a thirty-case delta, but this source does not pretend
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
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
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

/// Exact number of Rust tests in the ratified pre-manifest base inventory.
pub const BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1: usize = 130;

/// Exact Rust-test delta added by the same implementation train.
pub const BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1: usize = 30;

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
pub const BASE_COVERAGE_UNIT_CLASS_CASE_COUNT_V1: usize = 10;

/// Exact number of boundary and checked-arithmetic cases.
pub const BASE_COVERAGE_BOUNDARY_CASE_COUNT_V1: usize = 39;

/// Exact number of property and metamorphic cases.
pub const BASE_COVERAGE_PROPERTY_METAMORPHIC_CASE_COUNT_V1: usize = 17;

/// Exact number of schema and descriptor cases.
pub const BASE_COVERAGE_SCHEMA_DESCRIPTOR_CASE_COUNT_V1: usize = 39;

/// Exact number of mutation and malformed-presentation cases.
pub const BASE_COVERAGE_MUTATION_CASE_COUNT_V1: usize = 41;

/// Exact number of no-mock, in-process public-API integration cases.
pub const BASE_COVERAGE_NO_MOCK_INTEGRATION_CASE_COUNT_V1: usize = 14;

/// Exact number of compile-fail contracts in the ratified base inventory.
pub const BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1: usize = 47;

/// Exact number of unit tests that protect this manifest contract itself.
pub const BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1: usize = 10;

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
    /// The 47 source-frozen Rustdoc compile-fail contracts.
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
    /// Construct the frozen 160-Rust-test, 47-compile-fail, and
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
                    "magic_is_part_of_the_bound",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "integer_and_presence_fields_have_independent_known_big_endian_bytes",
                    "byte_and_string_fields_use_exact_u32_length_prefixes",
                ],
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
        class_templates: &[classified_tests(
            BaseCoverageManifestClassV1::Boundary,
            &["observed_rendering_is_utf8_bounded_without_recursive_diagnostics"],
        )],
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
                &["collection_bounds_and_reproduction_shape_fail_closed"],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::PropertyMetamorphic,
                &[
                    "feature_and_target_roots_are_deterministic_canonical_and_sensitive",
                    "normalized_prefix_suffix_and_embedded_sensitive_aliases_refuse",
                    "canonical_event_and_log_roots_are_order_independent_but_mutation_sensitive",
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
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::NoMockIntegration,
                &[
                    "full_log_reconciles_sequences_journeys_rows_results_cells_and_counts",
                    "positive_and_expected_refusal_partitions_are_distinct_and_exact",
                    "mixed_semantic_rows_reconcile_exact_terminal_partitions",
                    "source_closure_green_counts_are_matches_not_expected_refusals",
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
                BaseCoverageManifestClassV1::Mutation,
                &[
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
                BaseCoverageManifestClassV1::Boundary,
                &[
                    "rational_reduces_sign_zero_and_i128_min_without_overflow",
                    "decimal_has_one_representation_and_refuses_range_crossing",
                    "every_integer_width_preserves_both_extrema_exactly",
                    "token_boundaries_and_segment_grammar_are_exact",
                    "text_and_opaque_bytes_enforce_exact_byte_caps",
                ],
            ),
            classified_tests(
                BaseCoverageManifestClassV1::SchemaDescriptor,
                &[
                    "ieee_wrappers_preserve_special_encodings_and_nan_payloads",
                    "numeric_tags_are_exact_and_nonrecursive",
                    "units_require_positive_canonical_scale_and_keep_exponent_order",
                    "typed_value_and_presence_tags_are_exact",
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
}

const COMPILE_FAIL_TEMPLATES_V1: &[CompileFailTemplateV1] = &[
    CompileFailTemplateV1 {
        module: "budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
        case_name: "no-postmutation-of-runner-budgets",
    },
    CompileFailTemplateV1 {
        module: "budget",
        source_path: "crates/fs-evidence-runner/src/budget.rs",
        case_name: "no-postmutation-of-admitted-runner-budgets",
    },
    CompileFailTemplateV1 {
        module: "capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
        case_name: "no-physical-acquisition-material-as-semantic-right",
    },
    CompileFailTemplateV1 {
        module: "capability",
        source_path: "crates/fs-evidence-runner/src/capability.rs",
        case_name: "immutable-policy-root",
    },
    CompileFailTemplateV1 {
        module: "command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
        case_name: "immutable-command-intent",
    },
    CompileFailTemplateV1 {
        module: "command",
        source_path: "crates/fs-evidence-runner/src/command.rs",
        case_name: "no-command-intent-as-authority-scope",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-terminal-extension",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-refusal-extension",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-authority-mint",
    },
    CompileFailTemplateV1 {
        module: "diagnostic",
        source_path: "crates/fs-evidence-runner/src/diagnostic.rs",
        case_name: "no-executable-repair",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "artifact-codec-has-no-executable-encoder",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "private-registered-artifact-role-fields",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "typed-registered-namespaces-cannot-cross-substitute",
    },
    CompileFailTemplateV1 {
        module: "extension",
        source_path: "crates/fs-evidence-runner/src/extension.rs",
        case_name: "private-logical-extent-fields",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "sealed-digest-domain",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "no-generic-digest-as-nominal-root",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-source-identity-constructor",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-lifecycle-log-constructor",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-durable-publication-constructor",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "private-authority-scope-constructor",
    },
    CompileFailTemplateV1 {
        module: "identity",
        source_path: "crates/fs-evidence-runner/src/identity.rs",
        case_name: "no-standalone-root-for-root-free-evaluator-members",
    },
    CompileFailTemplateV1 {
        module: "limits",
        source_path: "crates/fs-evidence-runner/src/limits.rs",
        case_name: "no-postmutation-of-runner-limits",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-raw-string-as-logical-path",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "private-logical-path-constructor",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-postmutation-of-logical-path",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-bundle-path-as-content-key",
    },
    CompileFailTemplateV1 {
        module: "path",
        source_path: "crates/fs-evidence-runner/src/path.rs",
        case_name: "no-postmutation-of-content-store-key",
    },
    CompileFailTemplateV1 {
        module: "projection",
        source_path: "crates/fs-evidence-runner/src/projection.rs",
        case_name: "immutable-source-closure",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "nominal-cancelled-root",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-generic-digest-cause",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "cause-payload-required",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-profile-filter-cause",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "private-not-run-basis-fields",
    },
    CompileFailTemplateV1 {
        module: "state",
        source_path: "crates/fs-evidence-runner/src/state.rs",
        case_name: "no-unvalidated-state-candidate-as-terminal",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-rational-fields",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-rational",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-decimal-fields",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-decimal",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "private-unit-fields",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-unit",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-string-as-stable-token",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-stable-token",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-string-as-bounded-text",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-bounded-text",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-vec-as-opaque-bytes",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-postmutation-of-opaque-bytes",
    },
    CompileFailTemplateV1 {
        module: "value",
        source_path: "crates/fs-evidence-runner/src/value.rs",
        case_name: "no-digest-as-typed-absence",
    },
];

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
];

fn frozen_base_declarations() -> Result<Vec<BaseCoverageCaseDeclarationV1>, ConstructionErrorV2> {
    let capacity = BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1
        + BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1
        + BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1;
    let mut declarations = Vec::with_capacity(capacity);
    declarations.extend(frozen_rust_test_declarations()?);
    for case in COMPILE_FAIL_TEMPLATES_V1 {
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
            "the exact 160 classified Rust tests, 47 compile-fail cases, and 10 manifest-contract cases",
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
            "exactly 160 tests partitioned once across all six required evidence classes",
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
    ConstructionErrorV2::new(kind, field, expected, observed)
}

#[cfg(test)]
mod tests {
    #[allow(
        deprecated,
        reason = "one compatibility assertion freezes the misleading historical aggregate alias"
    )]
    use super::{
        BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1, BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1,
        BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1,
        BASE_COVERAGE_PREEXISTING_UNIT_CASE_COUNT_V1, BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1,
        BASE_COVERAGE_RUST_TEST_CLASS_COUNTS_V1, BASE_COVERAGE_UNIT_CASE_COUNT_V1,
        BaseCoverageCaseDeclarationV1, BaseCoverageCheckedReportV1, BaseCoverageManifestClassV1,
        BaseCoverageManifestV1, BaseCoveragePresentedOutcomeV1, BaseCoveragePresentedResultV1,
    };
    use crate::ConstructionErrorKindV2;
    use fs_blake3::{ContentHash, hash_domain};

    fn root(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.fs-evidence-runner.coverage-test.v1",
            label.as_bytes(),
        )
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
        assert_eq!(BASE_COVERAGE_POST_RATIFICATION_UNIT_CASE_DELTA_V1, 30);
        assert_eq!(BASE_COVERAGE_RUST_TEST_CASE_COUNT_V1, 160);
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
        assert_eq!(BASE_COVERAGE_COMPILE_FAIL_CASE_COUNT_V1, 47);
        assert_eq!(
            first.case_count(BaseCoverageManifestClassV1::ManifestContract),
            BASE_COVERAGE_MANIFEST_CONTRACT_CASE_COUNT_V1
        );
        assert_eq!(first.cases().len(), 217);
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
}
