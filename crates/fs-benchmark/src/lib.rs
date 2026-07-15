//! fs-benchmark — the wedge-vertical benchmark & trace corpus (plan addendum,
//! Proposal 7). Layer: UTIL (versioned data + measurement helpers).
//!
//! Governance Rule 2 (doctrine): "a proposal whose kill measurement was never
//! instrumented counts as killed — unmeasured survival is not survival." Many
//! addendum kill criteria measure against "the wedge vertical's benchmark set /
//! recorded traces / merge trials", yet no other bead owns that artifact. THIS
//! is the single, shared, versioned, DETERMINISTIC corpus and the fail-closed
//! evaluator-schema boundary. A declared schema does not count as
//! instrumentation until all of its typed retained evidence reconstructs.
//!
//! It bundles, for the conjugate-heat-transfer (electronics cooling) wedge:
//! a [`QueryCase`] set (QoIs, units, tolerances, and evidence-backed reference
//! answers), an optimizer [`DesignTask`] set, recorded [`EditTrace`]s with
//! retained numerator/denominator evidence, an [`MmsCase`] elliptic battery,
//! and swarm [`MergeTrial`]s. [`reconstruct_metric`] resolves typed metric
//! inputs rather than trusting detached numbers. The schema-versioned,
//! length-framed [`corpus_digest`] binds every semantic field.

pub use fs_blake3::ContentHash;
pub use fs_evidence::{
    AdmissionDecision, AdmissionReceipt, AdmissionRejection, AdmissionVerifier, AdmittedColor,
    COLOR_ALGEBRA_VERSION, Color, ColorRank, NoAdmissionVerifier,
};

use fs_blake3::hash_domain;
use fs_evidence::validate_color_payload;
use std::sync::LazyLock;

/// The corpus version (measurements are only comparable within a version).
pub const BENCHMARK_VERSION: u32 = 3;
/// Canonical length-framed corpus-identity schema.
pub const CORPUS_IDENTITY_SCHEMA_VERSION: u32 = 1;
/// Canonical retained-evidence identity schema.
pub const EVIDENCE_IDENTITY_SCHEMA_VERSION: u32 = 1;

const CORPUS_IDENTITY_DOMAIN: &str = "org.frankensim.fs-benchmark.corpus.v1";
const EVIDENCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-benchmark.evidence.v1";
const METRIC_IDENTITY_DOMAIN: &str = "org.frankensim.fs-benchmark.metric.v1";
const QUERY_ADMISSION_NODE_DOMAIN: &str = "org.frankensim.fs-benchmark.query-admission-node.v1";
const METRIC_IDENTITY_SCHEMA_VERSION: u32 = 1;
const MAX_EXACT_F64_COUNT: u64 = 1_u64 << 53;

/// How an evidence reference must be resolved.
///
/// There is deliberately no permissive or "best effort" policy. A row resolves
/// only when exactly one retained record has the requested id and its content
/// digest matches the consumer's independently retained identity. Resolution
/// establishes integrity, not scientific admission; positive color authority
/// still requires [`AdmittedColor`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceResolutionPolicy {
    /// Exact id and exact semantic content; consumers separately bind context.
    ExactContentV1,
}

impl EvidenceResolutionPolicy {
    const fn tag(self) -> u8 {
        match self {
            Self::ExactContentV1 => 1,
        }
    }
}

/// A reference to retained evidence. Resolution is fail-closed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceRef {
    /// Stable retained-evidence id.
    pub id: &'static str,
    /// Exact semantic digest expected by the consuming corpus row.
    pub expected_digest: ContentHash,
    /// Required resolution policy.
    pub policy: EvidenceResolutionPolicy,
}

impl EvidenceRef {
    /// Construct an exact retained-evidence reference.
    #[must_use]
    pub const fn exact(id: &'static str, expected_digest: ContentHash) -> Self {
        Self {
            id,
            expected_digest,
            policy: EvidenceResolutionPolicy::ExactContentV1,
        }
    }
}

/// Semantic role authenticated as part of a retained evidence record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EvidenceRole {
    /// Query answer, tolerance, cost, and epistemic declaration.
    QueryReference,
    /// Generic speedup baseline cost.
    BaselineCost,
    /// Generic speedup candidate cost.
    CandidateCost,
    /// Generic rate numerator.
    RateNumerator,
    /// Generic rate denominator.
    RateDenominator,
    /// Proposal 8 fixed-planner baseline cost.
    PlannerBaselineCost,
    /// Proposal 8 candidate-planner cost.
    PlannerCandidateCost,
    /// Proposal 1 adjoint task wins.
    AdjointWins,
    /// Proposal 1 matched-budget task comparisons.
    AdjointComparisons,
    /// Proposal 2 hash-memoization baseline cost.
    SkipBaselineCost,
    /// Proposal 2 certified-skip candidate cost.
    SkipCandidateCost,
    /// Proposal D guard-endpoint catches.
    GuardEndpointCatches,
    /// Proposal D guard-endpoint trials.
    GuardEndpointTrials,
    /// Proposal D random-design catches.
    GuardRandomCatches,
    /// Proposal D random-design trials.
    GuardRandomTrials,
    /// Proposal F robust results not dominated on realized cost.
    RobustNotDominated,
    /// Proposal F paired realized-cost comparisons.
    RobustComparisons,
    /// Proposal 9 accepted speculative solves.
    SpeculationAccepts,
    /// Proposal 9 attempted speculative solves.
    SpeculationAttempts,
    /// Proposal 9 cold-start cost.
    ColdStartCost,
    /// Proposal 9 warm-start cost.
    WarmStartCost,
    /// Proposal 10 unresolved candidate conflicts.
    MergeConflicts,
    /// Proposal 10 escalations.
    MergeEscalations,
    /// Proposal 10 refusals.
    MergeRefusals,
    /// Proposal 10 type conflicts.
    MergeTypeConflicts,
    /// Proposal 10 retained realistic merge attempts.
    MergeAttempts,
    /// Proposal A RB-certified query volume.
    CoveredQueryVolume,
    /// Proposal A total query volume.
    TotalQueryVolume,
    /// Non-promoting edit skip-count diagnostic.
    EditDiagnosticSkips,
    /// Non-promoting edit total-op diagnostic.
    EditDiagnosticTotalOps,
    /// Non-promoting synthetic merge-conflict diagnostic.
    MergeDiagnosticConflicts,
    /// Non-promoting synthetic merge-total diagnostic.
    MergeDiagnosticTotal,
}

impl EvidenceRole {
    const fn tag(self) -> u8 {
        match self {
            Self::QueryReference => 1,
            Self::BaselineCost => 2,
            Self::CandidateCost => 3,
            Self::RateNumerator => 4,
            Self::RateDenominator => 5,
            Self::PlannerBaselineCost => 6,
            Self::PlannerCandidateCost => 7,
            Self::AdjointWins => 8,
            Self::AdjointComparisons => 9,
            Self::SkipBaselineCost => 10,
            Self::SkipCandidateCost => 11,
            Self::GuardEndpointCatches => 12,
            Self::GuardEndpointTrials => 13,
            Self::GuardRandomCatches => 14,
            Self::GuardRandomTrials => 15,
            Self::RobustNotDominated => 16,
            Self::RobustComparisons => 17,
            Self::SpeculationAccepts => 18,
            Self::SpeculationAttempts => 19,
            Self::ColdStartCost => 20,
            Self::WarmStartCost => 21,
            Self::MergeConflicts => 22,
            Self::MergeEscalations => 23,
            Self::MergeRefusals => 24,
            Self::MergeTypeConflicts => 25,
            Self::MergeAttempts => 26,
            Self::CoveredQueryVolume => 27,
            Self::TotalQueryVolume => 28,
            Self::EditDiagnosticSkips => 29,
            Self::EditDiagnosticTotalOps => 30,
            Self::MergeDiagnosticConflicts => 31,
            Self::MergeDiagnosticTotal => 32,
        }
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::QueryReference => "query-reference",
            Self::BaselineCost => "baseline-cost",
            Self::CandidateCost => "candidate-cost",
            Self::RateNumerator => "rate-numerator",
            Self::RateDenominator => "rate-denominator",
            Self::PlannerBaselineCost => "planner-baseline-cost",
            Self::PlannerCandidateCost => "planner-candidate-cost",
            Self::AdjointWins => "adjoint-wins",
            Self::AdjointComparisons => "adjoint-comparisons",
            Self::SkipBaselineCost => "skip-baseline-cost",
            Self::SkipCandidateCost => "skip-candidate-cost",
            Self::GuardEndpointCatches => "guard-endpoint-catches",
            Self::GuardEndpointTrials => "guard-endpoint-trials",
            Self::GuardRandomCatches => "guard-random-catches",
            Self::GuardRandomTrials => "guard-random-trials",
            Self::RobustNotDominated => "robust-not-dominated",
            Self::RobustComparisons => "robust-comparisons",
            Self::SpeculationAccepts => "speculation-accepts",
            Self::SpeculationAttempts => "speculation-attempts",
            Self::ColdStartCost => "cold-start-cost",
            Self::WarmStartCost => "warm-start-cost",
            Self::MergeConflicts => "merge-conflicts",
            Self::MergeEscalations => "merge-escalations",
            Self::MergeRefusals => "merge-refusals",
            Self::MergeTypeConflicts => "merge-type-conflicts",
            Self::MergeAttempts => "merge-attempts",
            Self::CoveredQueryVolume => "covered-query-volume",
            Self::TotalQueryVolume => "total-query-volume",
            Self::EditDiagnosticSkips => "edit-diagnostic-skips",
            Self::EditDiagnosticTotalOps => "edit-diagnostic-total-ops",
            Self::MergeDiagnosticConflicts => "merge-diagnostic-conflicts",
            Self::MergeDiagnosticTotal => "merge-diagnostic-total",
        }
    }

    const fn datum_tag(self) -> u8 {
        match self {
            Self::QueryReference => 3,
            Self::BaselineCost
            | Self::CandidateCost
            | Self::PlannerBaselineCost
            | Self::PlannerCandidateCost
            | Self::SkipBaselineCost
            | Self::SkipCandidateCost
            | Self::ColdStartCost
            | Self::WarmStartCost
            | Self::CoveredQueryVolume
            | Self::TotalQueryVolume => 1,
            Self::RateNumerator
            | Self::RateDenominator
            | Self::AdjointWins
            | Self::AdjointComparisons
            | Self::GuardEndpointCatches
            | Self::GuardEndpointTrials
            | Self::GuardRandomCatches
            | Self::GuardRandomTrials
            | Self::RobustNotDominated
            | Self::RobustComparisons
            | Self::SpeculationAccepts
            | Self::SpeculationAttempts
            | Self::MergeConflicts
            | Self::MergeEscalations
            | Self::MergeRefusals
            | Self::MergeTypeConflicts
            | Self::MergeAttempts
            | Self::EditDiagnosticSkips
            | Self::EditDiagnosticTotalOps
            | Self::MergeDiagnosticConflicts
            | Self::MergeDiagnosticTotal => 2,
        }
    }

    const fn proposal_units(self) -> Option<&'static str> {
        match self {
            Self::PlannerBaselineCost
            | Self::PlannerCandidateCost
            | Self::SkipBaselineCost
            | Self::SkipCandidateCost
            | Self::ColdStartCost
            | Self::WarmStartCost => Some("work-unit"),
            Self::AdjointWins
            | Self::AdjointComparisons
            | Self::GuardEndpointCatches
            | Self::GuardEndpointTrials
            | Self::GuardRandomCatches
            | Self::GuardRandomTrials
            | Self::RobustNotDominated
            | Self::RobustComparisons
            | Self::SpeculationAccepts
            | Self::SpeculationAttempts
            | Self::MergeConflicts
            | Self::MergeEscalations
            | Self::MergeRefusals
            | Self::MergeTypeConflicts
            | Self::MergeAttempts => Some("count"),
            Self::CoveredQueryVolume | Self::TotalQueryVolume => Some("query-volume"),
            Self::QueryReference
            | Self::BaselineCost
            | Self::CandidateCost
            | Self::RateNumerator
            | Self::RateDenominator
            | Self::EditDiagnosticSkips
            | Self::EditDiagnosticTotalOps
            | Self::MergeDiagnosticConflicts
            | Self::MergeDiagnosticTotal => None,
        }
    }
}

/// Typed payload stored in a retained evidence record.
#[derive(Debug, Clone, PartialEq)]
pub enum EvidenceDatum {
    /// A finite scalar with no epistemic-rank claim.
    Scalar(f64),
    /// A non-negative integer count.
    Count(u64),
    /// A complete query reference declaration.
    ///
    /// Positive colors are candidates only until `admission_receipt` is
    /// authenticated by an injected [`AdmissionVerifier`]. Estimated colors
    /// carry no receipt and remain explicitly non-authoritative.
    QueryReference {
        /// Retained answer.
        answer: f64,
        /// Retained requested tolerance, in the query QoI units.
        tolerance: f64,
        /// Retained reference compute cost.
        reference_cost: f64,
        /// Explicit cost units.
        reference_cost_units: &'static str,
        /// Full epistemic candidate/declaration payload.
        color: Color,
        /// Admission receipt required for positive colors.
        admission_receipt: Option<AdmissionReceipt>,
    },
}

/// One retained evidence record.
///
/// The stable id is a locator. A consuming row or proposal manifest retains a
/// separate [`EvidenceRef`] for the remaining semantic content, so a record
/// never signs its own fields.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceRecord {
    /// Stable evidence id.
    pub id: &'static str,
    /// Dataset row this evidence supports.
    pub subject_id: &'static str,
    /// Quantity represented by `datum`.
    pub quantity: &'static str,
    /// Explicit units or the literal `count`.
    pub units: &'static str,
    /// Exact evaluator semantics that produced the datum.
    pub evaluator_semantics: &'static str,
    /// Exact proposal/governance semantics under which it is interpreted.
    pub proposal_semantics: &'static str,
    /// Authenticated semantic role; metric positions cannot relabel this field.
    pub role: EvidenceRole,
    /// Typed retained value.
    pub datum: EvidenceDatum,
}

impl EvidenceRecord {
    /// Construct a retained finite scalar.
    #[must_use]
    pub fn scalar(
        id: &'static str,
        subject_id: &'static str,
        quantity: &'static str,
        units: &'static str,
        evaluator_semantics: &'static str,
        proposal_semantics: &'static str,
        role: EvidenceRole,
        value: f64,
    ) -> Self {
        Self {
            id,
            subject_id,
            quantity,
            units,
            evaluator_semantics,
            proposal_semantics,
            role,
            datum: EvidenceDatum::Scalar(value),
        }
    }

    /// Construct a retained count.
    #[must_use]
    pub fn count(
        id: &'static str,
        subject_id: &'static str,
        quantity: &'static str,
        evaluator_semantics: &'static str,
        proposal_semantics: &'static str,
        role: EvidenceRole,
        value: u64,
    ) -> Self {
        Self {
            id,
            subject_id,
            quantity,
            units: "count",
            evaluator_semantics,
            proposal_semantics,
            role,
            datum: EvidenceDatum::Count(value),
        }
    }

    /// Construct a complete retained query-reference declaration/candidate.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn query_reference(
        id: &'static str,
        subject_id: &'static str,
        quantity: &'static str,
        units: &'static str,
        evaluator_semantics: &'static str,
        proposal_semantics: &'static str,
        answer: f64,
        tolerance: f64,
        reference_cost: f64,
        reference_cost_units: &'static str,
        color: Color,
        admission_receipt: Option<AdmissionReceipt>,
    ) -> Self {
        Self {
            id,
            subject_id,
            quantity,
            units,
            evaluator_semantics,
            proposal_semantics,
            role: EvidenceRole::QueryReference,
            datum: EvidenceDatum::QueryReference {
                answer,
                tolerance,
                reference_cost,
                reference_cost_units,
                color,
                admission_receipt,
            },
        }
    }

    /// Reconstruct this record's semantic digest.
    #[must_use]
    pub fn semantic_digest(&self) -> ContentHash {
        evidence_semantic_digest(self)
    }

    /// Context identity a positive query receipt's node hash must authenticate.
    /// It binds the locator, answer, tolerance, cost, units, semantics, role,
    /// and full candidate color, while excluding the receipt itself.
    pub fn query_admission_node_hash(&self) -> Result<ContentHash, EvidenceError> {
        query_admission_node_hash(self)
    }
}

/// One conjugate-heat-transfer design query with a retained reference answer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QueryCase {
    /// Stable id.
    pub id: &'static str,
    /// The quantity of interest.
    pub qoi: &'static str,
    /// Explicit QoI units.
    pub units: &'static str,
    /// The requested tolerance.
    pub tolerance: f64,
    /// The retained reference answer.
    pub reference_answer: f64,
    /// The reference compute cost to reach it.
    pub reference_cost: f64,
    /// Explicit reference-cost units.
    pub reference_cost_units: &'static str,
    /// Retained evidence that supplies the complete color candidate. Positive
    /// rank still requires an injected admission verifier at resolution time.
    pub reference_evidence: EvidenceRef,
    /// Exact evaluator semantics required from the retained record.
    pub reference_evaluator_semantics: &'static str,
    /// Exact proposal semantics required from the retained record.
    pub reference_proposal_semantics: &'static str,
}

/// One optimizer benchmark design task with a known optimum (for adjoint vs
/// derivative-free comparison).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignTask {
    /// Stable id.
    pub id: &'static str,
    /// Design dimension.
    pub dimension: usize,
    /// The known optimal objective value.
    pub optimum: f64,
}

/// A recorded agent design-iteration edit trace with a known-correct skip set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EditTrace {
    /// Stable id.
    pub id: &'static str,
    /// Total ops in the DAG for this edit.
    pub total_ops: usize,
    /// Ops that are CERTIFIABLY skippable (the known-correct skip set size).
    pub correct_skips: usize,
    /// Evidence for `correct_skips`.
    pub correct_skips_evidence: EvidenceRef,
    /// Evidence for `total_ops`, the metric denominator.
    pub total_ops_evidence: EvidenceRef,
    /// Exact evaluator semantics for the known-correct skip set.
    pub evaluator_semantics: &'static str,
    /// Exact proposal semantics for interpreting the skip fraction.
    pub proposal_semantics: &'static str,
}

/// An elliptic manufactured-solution case with a high-precision reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MmsCase {
    /// Stable id.
    pub id: &'static str,
    /// The exact solution sampled at the domain center (reference value).
    pub exact_center: f64,
}

/// A synthetic swarm-concurrency fixture with its candidate-remainder count.
///
/// These small rows exercise the corpus shape and measurement API. They are not
/// retained realistic merge traces and do not discharge Proposal 10's broader
/// unresolved-merge gate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeTrial {
    /// Stable id.
    pub id: &'static str,
    /// Total merges attempted.
    pub total_merges: usize,
    /// Synthetic merges assigned a fixed-iteration candidate remainder.
    pub candidate_remainder_conflicts: usize,
    /// Evidence for `candidate_remainder_conflicts`.
    pub conflict_count_evidence: EvidenceRef,
    /// Evidence for `total_merges`, the metric denominator.
    pub total_merges_evidence: EvidenceRef,
    /// Exact evaluator semantics for the synthetic diagnostic.
    pub evaluator_semantics: &'static str,
    /// Exact proposal semantics for interpreting the diagnostic rate.
    pub proposal_semantics: &'static str,
}

const QUERY_REFERENCE_EVALUATOR: &str = "cht-reference-solve-v1";
const QUERY_REFERENCE_PROPOSAL: &str = "benchmark-query-reference-v1";
const EDIT_SKIP_EVALUATOR: &str = "edit-dag-known-skip-set-v1";
const EDIT_SKIP_PROPOSAL: &str = "addendum-proposal-2:certified-skip-yield-v1";
const MERGE_EVALUATOR: &str = "fixed-iteration-candidate-remainder-count-v1";
const MERGE_PROPOSAL: &str = "addendum-proposal-10:synthetic-diagnostic-v1";

fn legacy_estimated_reference(
    id: &'static str,
    subject_id: &'static str,
    quantity: &'static str,
    units: &'static str,
    answer: f64,
    tolerance: f64,
    reference_cost: f64,
) -> EvidenceRecord {
    EvidenceRecord::query_reference(
        id,
        subject_id,
        quantity,
        units,
        QUERY_REFERENCE_EVALUATOR,
        QUERY_REFERENCE_PROPOSAL,
        answer,
        tolerance,
        reference_cost,
        "work-unit",
        Color::Estimated {
            estimator: format!("fs-benchmark:legacy-unadmitted:{subject_id}:v1"),
            dispersion: f64::INFINITY,
        },
        None,
    )
}

static RETAINED_EVIDENCE: LazyLock<Vec<EvidenceRecord>> = LazyLock::new(|| {
    vec![
        legacy_estimated_reference(
            "evidence:query:cht-q1:reference:v1",
            "cht-q1",
            "max junction temperature",
            "K",
            358.2,
            0.5,
            1000.0,
        ),
        legacy_estimated_reference(
            "evidence:query:cht-q2:reference:v1",
            "cht-q2",
            "board-to-ambient pressure drop",
            "Pa",
            42.0,
            5.0,
            600.0,
        ),
        legacy_estimated_reference(
            "evidence:query:cht-q3:reference:v1",
            "cht-q3",
            "hotspot thermal margin",
            "K",
            7.4,
            1.0,
            250.0,
        ),
        EvidenceRecord::count(
            "evidence:edit:raise-fin-count:correct-skips:v1",
            "raise-fin-count",
            "correct-skips",
            EDIT_SKIP_EVALUATOR,
            EDIT_SKIP_PROPOSAL,
            EvidenceRole::EditDiagnosticSkips,
            96,
        ),
        EvidenceRecord::count(
            "evidence:edit:raise-fin-count:total-ops:v1",
            "raise-fin-count",
            "total-ops",
            EDIT_SKIP_EVALUATOR,
            EDIT_SKIP_PROPOSAL,
            EvidenceRole::EditDiagnosticTotalOps,
            120,
        ),
        EvidenceRecord::count(
            "evidence:edit:move-inlet:correct-skips:v1",
            "move-inlet",
            "correct-skips",
            EDIT_SKIP_EVALUATOR,
            EDIT_SKIP_PROPOSAL,
            EvidenceRole::EditDiagnosticSkips,
            40,
        ),
        EvidenceRecord::count(
            "evidence:edit:move-inlet:total-ops:v1",
            "move-inlet",
            "total-ops",
            EDIT_SKIP_EVALUATOR,
            EDIT_SKIP_PROPOSAL,
            EvidenceRole::EditDiagnosticTotalOps,
            80,
        ),
        EvidenceRecord::count(
            "evidence:merge:two-agent-fin-vs-duct:conflicts:v1",
            "two-agent-fin-vs-duct",
            "candidate-remainder-conflicts",
            MERGE_EVALUATOR,
            MERGE_PROPOSAL,
            EvidenceRole::MergeDiagnosticConflicts,
            6,
        ),
        EvidenceRecord::count(
            "evidence:merge:two-agent-fin-vs-duct:total:v1",
            "two-agent-fin-vs-duct",
            "total-merges",
            MERGE_EVALUATOR,
            MERGE_PROPOSAL,
            EvidenceRole::MergeDiagnosticTotal,
            40,
        ),
        EvidenceRecord::count(
            "evidence:merge:three-agent-layout:conflicts:v1",
            "three-agent-layout",
            "candidate-remainder-conflicts",
            MERGE_EVALUATOR,
            MERGE_PROPOSAL,
            EvidenceRole::MergeDiagnosticConflicts,
            13,
        ),
        EvidenceRecord::count(
            "evidence:merge:three-agent-layout:total:v1",
            "three-agent-layout",
            "total-merges",
            MERGE_EVALUATOR,
            MERGE_PROPOSAL,
            EvidenceRole::MergeDiagnosticTotal,
            60,
        ),
    ]
});

fn built_in_evidence_ref(id: &'static str) -> EvidenceRef {
    let record = RETAINED_EVIDENCE
        .iter()
        .find(|record| record.id == id)
        .expect("built-in evidence id must exist");
    EvidenceRef::exact(record.id, record.semantic_digest())
}

static QUERY_SET: LazyLock<Vec<QueryCase>> = LazyLock::new(|| {
    vec![
        QueryCase {
            id: "cht-q1",
            qoi: "max junction temperature",
            units: "K",
            tolerance: 0.5,
            reference_answer: 358.2,
            reference_cost: 1000.0,
            reference_cost_units: "work-unit",
            reference_evidence: built_in_evidence_ref("evidence:query:cht-q1:reference:v1"),
            reference_evaluator_semantics: QUERY_REFERENCE_EVALUATOR,
            reference_proposal_semantics: QUERY_REFERENCE_PROPOSAL,
        },
        QueryCase {
            id: "cht-q2",
            qoi: "board-to-ambient pressure drop",
            units: "Pa",
            tolerance: 5.0,
            reference_answer: 42.0,
            reference_cost: 600.0,
            reference_cost_units: "work-unit",
            reference_evidence: built_in_evidence_ref("evidence:query:cht-q2:reference:v1"),
            reference_evaluator_semantics: QUERY_REFERENCE_EVALUATOR,
            reference_proposal_semantics: QUERY_REFERENCE_PROPOSAL,
        },
        QueryCase {
            id: "cht-q3",
            qoi: "hotspot thermal margin",
            units: "K",
            tolerance: 1.0,
            reference_answer: 7.4,
            reference_cost: 250.0,
            reference_cost_units: "work-unit",
            reference_evidence: built_in_evidence_ref("evidence:query:cht-q3:reference:v1"),
            reference_evaluator_semantics: QUERY_REFERENCE_EVALUATOR,
            reference_proposal_semantics: QUERY_REFERENCE_PROPOSAL,
        },
    ]
});

const DESIGN_TASKS: [DesignTask; 3] = [
    DesignTask {
        id: "heatsink-fin-pitch",
        dimension: 4,
        optimum: 351.0,
    },
    DesignTask {
        id: "duct-routing",
        dimension: 8,
        optimum: 38.5,
    },
    DesignTask {
        id: "component-placement",
        dimension: 16,
        optimum: 9.1,
    },
];

static EDIT_TRACES: LazyLock<Vec<EditTrace>> = LazyLock::new(|| {
    vec![
        EditTrace {
            id: "raise-fin-count",
            total_ops: 120,
            correct_skips: 96,
            correct_skips_evidence: built_in_evidence_ref(
                "evidence:edit:raise-fin-count:correct-skips:v1",
            ),
            total_ops_evidence: built_in_evidence_ref("evidence:edit:raise-fin-count:total-ops:v1"),
            evaluator_semantics: EDIT_SKIP_EVALUATOR,
            proposal_semantics: EDIT_SKIP_PROPOSAL,
        },
        EditTrace {
            id: "move-inlet",
            total_ops: 80,
            correct_skips: 40,
            correct_skips_evidence: built_in_evidence_ref(
                "evidence:edit:move-inlet:correct-skips:v1",
            ),
            total_ops_evidence: built_in_evidence_ref("evidence:edit:move-inlet:total-ops:v1"),
            evaluator_semantics: EDIT_SKIP_EVALUATOR,
            proposal_semantics: EDIT_SKIP_PROPOSAL,
        },
    ]
});

const MMS_BATTERY: [MmsCase; 2] = [
    MmsCase {
        id: "poisson-sin",
        exact_center: 1.0,
    },
    MmsCase {
        id: "aniso-diffusion",
        exact_center: 0.25,
    },
];

static MERGE_TRIALS: LazyLock<Vec<MergeTrial>> = LazyLock::new(|| {
    vec![
        MergeTrial {
            id: "two-agent-fin-vs-duct",
            total_merges: 40,
            candidate_remainder_conflicts: 6,
            conflict_count_evidence: built_in_evidence_ref(
                "evidence:merge:two-agent-fin-vs-duct:conflicts:v1",
            ),
            total_merges_evidence: built_in_evidence_ref(
                "evidence:merge:two-agent-fin-vs-duct:total:v1",
            ),
            evaluator_semantics: MERGE_EVALUATOR,
            proposal_semantics: MERGE_PROPOSAL,
        },
        MergeTrial {
            id: "three-agent-layout",
            total_merges: 60,
            candidate_remainder_conflicts: 13,
            conflict_count_evidence: built_in_evidence_ref(
                "evidence:merge:three-agent-layout:conflicts:v1",
            ),
            total_merges_evidence: built_in_evidence_ref(
                "evidence:merge:three-agent-layout:total:v1",
            ),
            evaluator_semantics: MERGE_EVALUATOR,
            proposal_semantics: MERGE_PROPOSAL,
        },
    ]
});

/// The conjugate-heat-transfer query set.
#[must_use]
pub fn query_set() -> &'static [QueryCase] {
    QUERY_SET.as_slice()
}

/// The optimizer benchmark design tasks.
#[must_use]
pub fn design_tasks() -> &'static [DesignTask] {
    &DESIGN_TASKS
}

/// The recorded edit traces.
#[must_use]
pub fn edit_traces() -> &'static [EditTrace] {
    EDIT_TRACES.as_slice()
}

/// The elliptic MMS battery.
#[must_use]
pub fn mms_battery() -> &'static [MmsCase] {
    &MMS_BATTERY
}

/// The swarm merge trials.
#[must_use]
pub fn merge_trials() -> &'static [MergeTrial] {
    MERGE_TRIALS.as_slice()
}

/// The retained evidence catalog backing declarations and typed diagnostics.
#[must_use]
pub fn retained_evidence() -> &'static [EvidenceRecord] {
    RETAINED_EVIDENCE.as_slice()
}

// -- Measurement helpers (each addendum kill measurement) ------------------

/// An evidence-resolution failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceError {
    /// No retained record has the requested id.
    Missing { id: &'static str },
    /// More than one retained record has the requested id.
    Duplicate { id: &'static str },
    /// Two different ids claim the same authoritative semantic content.
    DuplicateDigest {
        id: &'static str,
        other_id: &'static str,
    },
    /// A record's semantic content does not match the consumer's identity.
    Tampered { id: &'static str },
    /// A required string field is empty.
    EmptyField {
        id: &'static str,
        field: &'static str,
    },
    /// A retained scalar is non-finite.
    NonFinite { id: &'static str },
    /// A query-reference color payload is structurally invalid.
    MalformedColor { id: &'static str },
    /// A positive color candidate has no admission receipt.
    MissingAdmissionReceipt { id: &'static str },
    /// An Estimated declaration incorrectly carries an admission receipt.
    UnexpectedAdmissionReceipt { id: &'static str },
    /// The injected admission authority refused a positive color candidate.
    AdmissionRejected {
        id: &'static str,
        rejection: AdmissionRejection,
    },
    /// A positive receipt authenticates a different query context node.
    AdmissionContextMismatch { id: &'static str },
    /// A resolved record does not match its consumer's exact context.
    ContextMismatch {
        id: &'static str,
        field: &'static str,
    },
}

impl core::fmt::Display for EvidenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Missing { id } => write!(f, "retained evidence '{id}' is missing"),
            Self::Duplicate { id } => write!(f, "retained evidence id '{id}' is ambiguous"),
            Self::DuplicateDigest { id, other_id } => write!(
                f,
                "retained evidence '{id}' duplicates the content identity of '{other_id}'"
            ),
            Self::Tampered { id } => {
                write!(
                    f,
                    "retained evidence '{id}' mismatches its expected content identity"
                )
            }
            Self::EmptyField { id, field } => {
                write!(f, "retained evidence '{id}' has empty field '{field}'")
            }
            Self::NonFinite { id } => {
                write!(f, "retained evidence '{id}' has a non-finite scalar")
            }
            Self::MalformedColor { id } => {
                write!(f, "retained evidence '{id}' has a malformed color payload")
            }
            Self::MissingAdmissionReceipt { id } => {
                write!(
                    f,
                    "retained evidence '{id}' has a positive color without a receipt"
                )
            }
            Self::UnexpectedAdmissionReceipt { id } => write!(
                f,
                "retained evidence '{id}' attaches admission authority to an Estimated declaration"
            ),
            Self::AdmissionRejected { id, rejection } => {
                write!(f, "retained evidence '{id}' was not admitted: {rejection}")
            }
            Self::AdmissionContextMismatch { id } => write!(
                f,
                "retained evidence '{id}' receipt does not bind its complete query context"
            ),
            Self::ContextMismatch { id, field } => {
                write!(
                    f,
                    "retained evidence '{id}' mismatches context field '{field}'"
                )
            }
        }
    }
}

/// A metric reconstruction failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MetricError {
    /// Retained evidence could not be resolved exactly.
    Evidence(EvidenceError),
    /// A raw scalar input was non-finite.
    NonFiniteInput { role: &'static str },
    /// Finite inputs produced a non-finite quotient.
    NonFiniteResult,
    /// A scalar domain requires a non-negative value.
    NegativeInput { role: &'static str },
    /// The denominator is zero and therefore cannot support a metric.
    ZeroDenominator,
    /// A count numerator exceeds its denominator.
    NumeratorExceedsDenominator,
    /// A count cannot be represented exactly by the returned `f64` metric.
    InexactCount { role: &'static str },
    /// Both metric roles resolve the same retained evidence.
    SameEvidence,
    /// Numerator and denominator evidence have different units.
    UnitsMismatch,
    /// The typed payload is incompatible with the requested metric.
    WrongDatum {
        id: &'static str,
        expected: &'static str,
    },
    /// Evaluator/proposal semantics are absent or disagree with the evidence.
    SemanticsMismatch { field: &'static str },
    /// The retained record's authenticated role is wrong for its metric slot.
    RoleMismatch {
        expected: EvidenceRole,
        actual: EvidenceRole,
    },
    /// A proposal evaluator's predeclared evidence role is absent.
    MissingProposalRole { role: EvidenceRole },
    /// The metric definition kind is incompatible with its fixed schema.
    SchemaKindMismatch,
}

impl From<EvidenceError> for MetricError {
    fn from(value: EvidenceError) -> Self {
        Self::Evidence(value)
    }
}

impl core::fmt::Display for MetricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Evidence(error) => core::fmt::Display::fmt(error, f),
            Self::NonFiniteInput { role } => write!(f, "metric {role} is non-finite"),
            Self::NonFiniteResult => f.write_str("metric quotient is non-finite"),
            Self::NegativeInput { role } => write!(f, "metric {role} is negative"),
            Self::ZeroDenominator => f.write_str("metric denominator is zero"),
            Self::NumeratorExceedsDenominator => {
                f.write_str("metric numerator exceeds its denominator")
            }
            Self::InexactCount { role } => {
                write!(f, "metric {role} exceeds the exact f64 integer range")
            }
            Self::SameEvidence => {
                f.write_str("metric numerator and denominator use the same evidence")
            }
            Self::UnitsMismatch => f.write_str("metric evidence units do not match"),
            Self::WrongDatum { id, expected } => {
                write!(f, "metric evidence '{id}' is not a {expected}")
            }
            Self::SemanticsMismatch { field } => {
                write!(f, "metric evidence mismatches '{field}' semantics")
            }
            Self::RoleMismatch { expected, actual } => write!(
                f,
                "metric role '{}' cannot occupy '{}'",
                actual.name(),
                expected.name()
            ),
            Self::MissingProposalRole { role } => {
                write!(f, "proposal evidence role '{}' is missing", role.name())
            }
            Self::SchemaKindMismatch => {
                f.write_str("metric definition kind does not match its fixed schema")
            }
        }
    }
}

/// One metric role with an explicit evidence locator and expected context.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MetricEvidenceRole {
    /// Exact retained-evidence reference.
    pub reference: EvidenceRef,
    /// Dataset row the evidence must describe.
    pub subject_id: &'static str,
    /// Quantity this metric role requires.
    pub quantity: &'static str,
    /// Units this metric role requires.
    pub units: &'static str,
}

impl MetricEvidenceRole {
    /// Construct one exact metric evidence role.
    #[must_use]
    pub const fn exact(
        reference: EvidenceRef,
        subject_id: &'static str,
        quantity: &'static str,
        units: &'static str,
    ) -> Self {
        Self {
            reference,
            subject_id,
            quantity,
            units,
        }
    }
}

/// Typed reconstruction requested from retained evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricDefinition {
    /// `baseline_cost / candidate_cost`, reconstructed from finite scalars.
    Speedup {
        baseline: MetricEvidenceRole,
        candidate: MetricEvidenceRole,
    },
    /// `numerator / denominator`, reconstructed from integer counts.
    Rate {
        numerator: MetricEvidenceRole,
        denominator: MetricEvidenceRole,
    },
    /// `numerator / denominator`, reconstructed from non-negative scalars.
    Fraction {
        numerator: MetricEvidenceRole,
        denominator: MetricEvidenceRole,
    },
}

/// Fixed evaluator schema that owns the semantic roles of both metric slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricSchema {
    /// Generic scalar speedup.
    GenericSpeedup,
    /// Generic integer rate.
    GenericRate,
    /// Non-promoting edit skip fraction.
    EditSkipDiagnostic,
    /// Non-promoting synthetic merge-conflict fraction.
    MergeConflictDiagnostic,
    /// Proposal 8 planner speedup.
    PlannerSpeedup,
    /// Proposal 1 adjoint task win rate.
    AdjointWinRate,
    /// Proposal 2 certified-skip speedup.
    CertifiedSkipSpeedup,
    /// Proposal D guard endpoint catch rate.
    GuardEndpointRate,
    /// Proposal D random-design catch rate.
    GuardRandomRate,
    /// Proposal F not-dominated rate.
    RobustNotDominatedRate,
    /// Proposal 9 speculative accept rate.
    SpeculationAcceptRate,
    /// Proposal 9 warm-start speedup.
    WarmStartSpeedup,
    /// Proposal A covered query-volume fraction.
    QueryVolumeFraction,
}

impl MetricSchema {
    const fn tag(self) -> u8 {
        match self {
            Self::GenericSpeedup => 1,
            Self::GenericRate => 2,
            Self::EditSkipDiagnostic => 3,
            Self::MergeConflictDiagnostic => 4,
            Self::PlannerSpeedup => 5,
            Self::AdjointWinRate => 6,
            Self::CertifiedSkipSpeedup => 7,
            Self::GuardEndpointRate => 8,
            Self::GuardRandomRate => 9,
            Self::RobustNotDominatedRate => 10,
            Self::SpeculationAcceptRate => 11,
            Self::WarmStartSpeedup => 12,
            Self::QueryVolumeFraction => 13,
        }
    }

    const fn kind(self) -> MetricKind {
        match self {
            Self::GenericSpeedup
            | Self::PlannerSpeedup
            | Self::CertifiedSkipSpeedup
            | Self::WarmStartSpeedup => MetricKind::Speedup,
            Self::GenericRate
            | Self::EditSkipDiagnostic
            | Self::MergeConflictDiagnostic
            | Self::AdjointWinRate
            | Self::GuardEndpointRate
            | Self::GuardRandomRate
            | Self::RobustNotDominatedRate
            | Self::SpeculationAcceptRate => MetricKind::Rate,
            Self::QueryVolumeFraction => MetricKind::Fraction,
        }
    }

    const fn expected_roles(self) -> (EvidenceRole, EvidenceRole) {
        match self {
            Self::GenericSpeedup => (EvidenceRole::BaselineCost, EvidenceRole::CandidateCost),
            Self::GenericRate => (EvidenceRole::RateNumerator, EvidenceRole::RateDenominator),
            Self::EditSkipDiagnostic => (
                EvidenceRole::EditDiagnosticSkips,
                EvidenceRole::EditDiagnosticTotalOps,
            ),
            Self::MergeConflictDiagnostic => (
                EvidenceRole::MergeDiagnosticConflicts,
                EvidenceRole::MergeDiagnosticTotal,
            ),
            Self::PlannerSpeedup => (
                EvidenceRole::PlannerBaselineCost,
                EvidenceRole::PlannerCandidateCost,
            ),
            Self::AdjointWinRate => (EvidenceRole::AdjointWins, EvidenceRole::AdjointComparisons),
            Self::CertifiedSkipSpeedup => (
                EvidenceRole::SkipBaselineCost,
                EvidenceRole::SkipCandidateCost,
            ),
            Self::GuardEndpointRate => (
                EvidenceRole::GuardEndpointCatches,
                EvidenceRole::GuardEndpointTrials,
            ),
            Self::GuardRandomRate => (
                EvidenceRole::GuardRandomCatches,
                EvidenceRole::GuardRandomTrials,
            ),
            Self::RobustNotDominatedRate => (
                EvidenceRole::RobustNotDominated,
                EvidenceRole::RobustComparisons,
            ),
            Self::SpeculationAcceptRate => (
                EvidenceRole::SpeculationAccepts,
                EvidenceRole::SpeculationAttempts,
            ),
            Self::WarmStartSpeedup => (EvidenceRole::ColdStartCost, EvidenceRole::WarmStartCost),
            Self::QueryVolumeFraction => (
                EvidenceRole::CoveredQueryVolume,
                EvidenceRole::TotalQueryVolume,
            ),
        }
    }
}

/// Exact evaluator and proposal semantics for a typed metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MetricRequest {
    /// Fixed metric schema; callers cannot redefine positional roles.
    pub schema: MetricSchema,
    /// Exact evaluator algorithm/version.
    pub evaluator_semantics: &'static str,
    /// Exact proposal/governance interpretation.
    pub proposal_semantics: &'static str,
    /// Typed numerator/denominator role.
    pub definition: MetricDefinition,
}

/// The kind of reconstructed metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetricKind {
    /// Baseline/candidate scalar speedup.
    Speedup,
    /// Count/count rate.
    Rate,
    /// Non-negative scalar fraction.
    Fraction,
}

/// A reconstructed metric tied to both retained evidence records.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReconstructedMetric {
    /// Metric kind.
    pub kind: MetricKind,
    /// Reconstructed finite value.
    pub value: f64,
    /// Numerator/baseline evidence.
    pub numerator_evidence: EvidenceRef,
    /// Denominator/candidate evidence.
    pub denominator_evidence: EvidenceRef,
    /// Length-framed identity over the complete typed request.
    pub identity: ContentHash,
}

/// A query reference after exact retained-evidence resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedReference {
    /// Retained reference answer.
    answer: f64,
    /// Admitted positive evidence or an explicit non-authoritative declaration.
    authority: ReferenceAuthority,
    /// Digest that was reconstructed and checked.
    evidence_digest: ContentHash,
    /// Complete query context bound into any admitted authority.
    query_context_hash: ContentHash,
}

impl ResolvedReference {
    /// Retained reference answer.
    #[must_use]
    pub fn answer(&self) -> f64 {
        self.answer
    }

    /// Epistemic standing attached to this exact resolved context.
    #[must_use]
    pub fn authority(&self) -> &ReferenceAuthority {
        &self.authority
    }

    /// Digest reconstructed from the retained evidence record.
    #[must_use]
    pub fn evidence_digest(&self) -> ContentHash {
        self.evidence_digest
    }

    /// Complete query-context identity carried by this resolution.
    #[must_use]
    pub fn query_context_hash(&self) -> ContentHash {
        self.query_context_hash
    }
}

/// Opaque, structurally valid Estimated declaration. Its private field keeps a
/// caller from placing a positive `Color` behind the non-authoritative variant.
#[derive(Debug, Clone, PartialEq)]
pub struct EstimatedReference {
    color: Color,
}

/// Opaque admitted color whose provenance remains attached to the complete
/// query context authenticated by its receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedReference {
    color: AdmittedColor,
    query_context_hash: ContentHash,
}

impl AdmittedReference {
    /// The admitted rank.
    #[must_use]
    pub fn rank(&self) -> ColorRank {
        self.color.rank()
    }

    /// Complete query context authenticated by the admission receipt.
    #[must_use]
    pub fn query_context_hash(&self) -> ContentHash {
        self.query_context_hash
    }
}

impl EstimatedReference {
    /// The structurally validated Estimated declaration.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// Estimated declarations never report a positive rank.
    #[must_use]
    pub const fn rank(&self) -> ColorRank {
        ColorRank::Estimated
    }
}

/// Epistemic standing returned by query-reference resolution.
#[derive(Debug, Clone, PartialEq)]
pub enum ReferenceAuthority {
    /// Positive evidence authenticated by an injected admission capability.
    Admitted(AdmittedReference),
    /// Structurally valid Estimated declaration with no promotion authority.
    EstimatedDeclaration(EstimatedReference),
}

impl ReferenceAuthority {
    /// The represented rank. Positive ranks can only come from `AdmittedColor`.
    #[must_use]
    pub fn rank(&self) -> ColorRank {
        match self {
            Self::Admitted(color) => color.rank(),
            Self::EstimatedDeclaration(declaration) => declaration.rank(),
        }
    }

    /// Whether an admission authority authenticated this reference.
    #[must_use]
    pub const fn is_admitted(&self) -> bool {
        matches!(self, Self::Admitted(_))
    }
}

/// A borrowed corpus view, used for deterministic replay and adversarial
/// validation without any global mutation.
#[derive(Debug, Clone, Copy)]
pub struct BenchmarkCorpus<'a> {
    /// Corpus semantic version.
    pub version: u32,
    /// Canonical identity schema version.
    pub identity_schema_version: u32,
    /// Query rows.
    pub query_set: &'a [QueryCase],
    /// Design-task rows.
    pub design_tasks: &'a [DesignTask],
    /// Edit-trace rows.
    pub edit_traces: &'a [EditTrace],
    /// MMS rows.
    pub mms_battery: &'a [MmsCase],
    /// Merge-trial rows.
    pub merge_trials: &'a [MergeTrial],
    /// Retained evidence records.
    pub retained_evidence: &'a [EvidenceRecord],
    /// Independently retained proposal-role references. Records without a
    /// matching reference remain unavailable even when present in the catalog.
    pub proposal_evidence: &'a [ProposalEvidenceRef],
    /// Instrumented proposal rows.
    pub instrumented_proposals: &'a [InstrumentedProposal],
}

/// The built-in immutable benchmark corpus.
#[must_use]
pub fn benchmark_corpus() -> BenchmarkCorpus<'static> {
    BenchmarkCorpus {
        version: BENCHMARK_VERSION,
        identity_schema_version: CORPUS_IDENTITY_SCHEMA_VERSION,
        query_set: query_set(),
        design_tasks: &DESIGN_TASKS,
        edit_traces: edit_traces(),
        mms_battery: &MMS_BATTERY,
        merge_trials: merge_trials(),
        retained_evidence: retained_evidence(),
        proposal_evidence: &[],
        instrumented_proposals: &INSTRUMENTED,
    }
}

fn validate_evidence_record(record: &EvidenceRecord) -> Result<(), EvidenceError> {
    for (field, value) in [
        ("id", record.id),
        ("subject-id", record.subject_id),
        ("quantity", record.quantity),
        ("units", record.units),
        ("evaluator-semantics", record.evaluator_semantics),
        ("proposal-semantics", record.proposal_semantics),
    ] {
        if value.trim().is_empty() {
            return Err(EvidenceError::EmptyField {
                id: record.id,
                field,
            });
        }
    }
    let datum_tag = match &record.datum {
        EvidenceDatum::Scalar(_) => 1,
        EvidenceDatum::Count(_) => 2,
        EvidenceDatum::QueryReference { .. } => 3,
    };
    if datum_tag != record.role.datum_tag() {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "role-datum",
        });
    }
    if matches!(&record.datum, EvidenceDatum::Count(_)) && record.units != "count" {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "count-units",
        });
    }
    match &record.datum {
        EvidenceDatum::Scalar(value) if !value.is_finite() => {
            Err(EvidenceError::NonFinite { id: record.id })
        }
        EvidenceDatum::QueryReference {
            answer,
            tolerance,
            reference_cost,
            reference_cost_units,
            color,
            admission_receipt,
        } => {
            if !answer.is_finite() || !tolerance.is_finite() || !reference_cost.is_finite() {
                return Err(EvidenceError::NonFinite { id: record.id });
            }
            if *tolerance <= 0.0 {
                return Err(EvidenceError::ContextMismatch {
                    id: record.id,
                    field: "tolerance",
                });
            }
            if *reference_cost <= 0.0 {
                return Err(EvidenceError::ContextMismatch {
                    id: record.id,
                    field: "reference-cost",
                });
            }
            if reference_cost_units.trim().is_empty() {
                return Err(EvidenceError::EmptyField {
                    id: record.id,
                    field: "reference-cost-units",
                });
            }
            validate_color_payload(color)
                .map_err(|_| EvidenceError::MalformedColor { id: record.id })?;
            if let Color::Verified { lo, hi } = color {
                if *answer < *lo || *answer > *hi {
                    return Err(EvidenceError::ContextMismatch {
                        id: record.id,
                        field: "verified-answer",
                    });
                }
            }
            match (color.rank(), admission_receipt) {
                (ColorRank::Estimated, None) => Ok(()),
                (ColorRank::Estimated, Some(_)) => {
                    Err(EvidenceError::UnexpectedAdmissionReceipt { id: record.id })
                }
                (_, None) => Err(EvidenceError::MissingAdmissionReceipt { id: record.id }),
                (_, Some(receipt)) => {
                    if receipt.node_hash() == record.query_admission_node_hash()? {
                        Ok(())
                    } else {
                        Err(EvidenceError::AdmissionContextMismatch { id: record.id })
                    }
                }
            }
        }
        _ => Ok(()),
    }
}

fn validate_evidence_registry(corpus: &BenchmarkCorpus<'_>) -> Result<(), EvidenceError> {
    for (index, record) in corpus.retained_evidence.iter().enumerate() {
        validate_evidence_record(record)?;
        let digest = record.semantic_digest();
        for other in &corpus.retained_evidence[..index] {
            if other.id == record.id {
                return Err(EvidenceError::Duplicate { id: record.id });
            }
            if other.semantic_digest() == digest {
                return Err(EvidenceError::DuplicateDigest {
                    id: record.id,
                    other_id: other.id,
                });
            }
        }
    }
    Ok(())
}

/// Resolve exactly one retained record and verify its expected content digest.
/// Any duplicate id or content identity anywhere in the borrowed registry
/// refuses the resolution, so callers cannot consume a partially ambiguous
/// corpus.
pub fn resolve_evidence(
    corpus: &BenchmarkCorpus<'_>,
    reference: EvidenceRef,
) -> Result<EvidenceRecord, EvidenceError> {
    match reference.policy {
        EvidenceResolutionPolicy::ExactContentV1 => {
            validate_evidence_registry(corpus)?;
        }
    }
    let mut matches = corpus
        .retained_evidence
        .iter()
        .cloned()
        .filter(|record| record.id == reference.id);
    let record = matches
        .next()
        .ok_or(EvidenceError::Missing { id: reference.id })?;
    if matches.next().is_some() {
        return Err(EvidenceError::Duplicate { id: reference.id });
    }
    let actual_digest = record.semantic_digest();
    if actual_digest != reference.expected_digest {
        return Err(EvidenceError::Tampered { id: reference.id });
    }
    Ok(record)
}

fn require_context(
    record: &EvidenceRecord,
    subject_id: &'static str,
    quantity: &'static str,
    units: &'static str,
    evaluator_semantics: &'static str,
    proposal_semantics: &'static str,
) -> Result<(), EvidenceError> {
    for (field, matches) in [
        ("subject-id", record.subject_id == subject_id),
        ("quantity", record.quantity == quantity),
        ("units", record.units == units),
        (
            "evaluator-semantics",
            record.evaluator_semantics == evaluator_semantics,
        ),
        (
            "proposal-semantics",
            record.proposal_semantics == proposal_semantics,
        ),
    ] {
        if !matches {
            return Err(EvidenceError::ContextMismatch {
                id: record.id,
                field,
            });
        }
    }
    Ok(())
}

/// Resolve a query under the default deny-all admission policy.
///
/// Estimated declarations may resolve, but positive evidence is refused at the
/// UTIL layer. Call [`resolve_query_reference_with_verifier`] only with the
/// actual admission authority.
pub fn resolve_query_reference(
    corpus: &BenchmarkCorpus<'_>,
    query: &QueryCase,
) -> Result<ResolvedReference, EvidenceError> {
    resolve_query_reference_with_verifier(corpus, query, &NoAdmissionVerifier)
}

/// Resolve a query reference and authenticate any positive color candidate.
pub fn resolve_query_reference_with_verifier(
    corpus: &BenchmarkCorpus<'_>,
    query: &QueryCase,
    verifier: &dyn AdmissionVerifier,
) -> Result<ResolvedReference, EvidenceError> {
    let record = resolve_evidence(corpus, query.reference_evidence)?;
    require_context(
        &record,
        query.id,
        query.qoi,
        query.units,
        query.reference_evaluator_semantics,
        query.reference_proposal_semantics,
    )?;
    let EvidenceDatum::QueryReference {
        answer,
        tolerance,
        reference_cost,
        reference_cost_units,
        color,
        admission_receipt,
    } = &record.datum
    else {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "query-reference-datum",
        });
    };
    let query_context_hash = record.query_admission_node_hash()?;
    for (field, matches) in [
        (
            "reference-answer",
            answer.to_bits() == query.reference_answer.to_bits(),
        ),
        (
            "tolerance",
            tolerance.to_bits() == query.tolerance.to_bits(),
        ),
        (
            "reference-cost",
            reference_cost.to_bits() == query.reference_cost.to_bits(),
        ),
        (
            "reference-cost-units",
            *reference_cost_units == query.reference_cost_units,
        ),
    ] {
        if !matches {
            return Err(EvidenceError::ContextMismatch {
                id: record.id,
                field,
            });
        }
    }
    let authority = match color.rank() {
        ColorRank::Estimated => ReferenceAuthority::EstimatedDeclaration(EstimatedReference {
            color: color.clone(),
        }),
        ColorRank::Validated | ColorRank::Verified => {
            let receipt = (*admission_receipt)
                .ok_or(EvidenceError::MissingAdmissionReceipt { id: record.id })?;
            let admitted = AdmittedColor::from_receipt(color.clone(), receipt, verifier).map_err(
                |rejection| EvidenceError::AdmissionRejected {
                    id: record.id,
                    rejection,
                },
            )?;
            ReferenceAuthority::Admitted(AdmittedReference {
                color: admitted,
                query_context_hash,
            })
        }
    };
    Ok(ResolvedReference {
        answer: *answer,
        authority,
        evidence_digest: record.semantic_digest(),
        query_context_hash,
    })
}

fn metric_context(
    request: &MetricRequest,
    numerator_role: MetricEvidenceRole,
    denominator_role: MetricEvidenceRole,
    numerator: &EvidenceRecord,
    denominator: &EvidenceRecord,
) -> Result<(), MetricError> {
    if request.evaluator_semantics.trim().is_empty() {
        return Err(MetricError::SemanticsMismatch {
            field: "evaluator-semantics",
        });
    }
    if request.proposal_semantics.trim().is_empty() {
        return Err(MetricError::SemanticsMismatch {
            field: "proposal-semantics",
        });
    }
    if numerator_role.reference.id == denominator_role.reference.id {
        return Err(MetricError::SameEvidence);
    }
    fn require_metric_role(
        request: &MetricRequest,
        role: MetricEvidenceRole,
        record: &EvidenceRecord,
        expected_role: EvidenceRole,
    ) -> Result<(), MetricError> {
        for (field, value) in [
            ("subject-id", role.subject_id),
            ("quantity", role.quantity),
            ("units", role.units),
        ] {
            if value.trim().is_empty() {
                return Err(MetricError::SemanticsMismatch { field });
            }
        }
        require_context(
            record,
            role.subject_id,
            role.quantity,
            role.units,
            request.evaluator_semantics,
            request.proposal_semantics,
        )?;
        if record.role != expected_role {
            return Err(MetricError::RoleMismatch {
                expected: expected_role,
                actual: record.role,
            });
        }
        Ok(())
    }
    let (expected_numerator, expected_denominator) = request.schema.expected_roles();
    require_metric_role(request, numerator_role, numerator, expected_numerator)?;
    require_metric_role(request, denominator_role, denominator, expected_denominator)?;
    if numerator_role.subject_id != denominator_role.subject_id {
        return Err(MetricError::SemanticsMismatch {
            field: "subject-id",
        });
    }
    if numerator_role.units != denominator_role.units {
        return Err(MetricError::UnitsMismatch);
    }
    Ok(())
}

/// Reconstruct a metric from typed numerator and denominator evidence.
pub fn reconstruct_metric(
    corpus: &BenchmarkCorpus<'_>,
    request: MetricRequest,
) -> Result<ReconstructedMetric, MetricError> {
    let (kind, numerator_role, denominator_role) = match request.definition {
        MetricDefinition::Speedup {
            baseline,
            candidate,
        } => (MetricKind::Speedup, baseline, candidate),
        MetricDefinition::Rate {
            numerator,
            denominator,
        } => (MetricKind::Rate, numerator, denominator),
        MetricDefinition::Fraction {
            numerator,
            denominator,
        } => (MetricKind::Fraction, numerator, denominator),
    };
    if kind != request.schema.kind() {
        return Err(MetricError::SchemaKindMismatch);
    }
    if numerator_role.reference.id == denominator_role.reference.id {
        return Err(MetricError::SameEvidence);
    }
    let numerator = resolve_evidence(corpus, numerator_role.reference)?;
    let denominator = resolve_evidence(corpus, denominator_role.reference)?;
    metric_context(
        &request,
        numerator_role,
        denominator_role,
        &numerator,
        &denominator,
    )?;

    let value = match kind {
        MetricKind::Speedup => {
            let EvidenceDatum::Scalar(baseline) = &numerator.datum else {
                return Err(MetricError::WrongDatum {
                    id: numerator.id,
                    expected: "scalar",
                });
            };
            let EvidenceDatum::Scalar(candidate) = &denominator.datum else {
                return Err(MetricError::WrongDatum {
                    id: denominator.id,
                    expected: "scalar",
                });
            };
            speedup(*baseline, *candidate)?
        }
        MetricKind::Rate => {
            let EvidenceDatum::Count(numerator_count) = &numerator.datum else {
                return Err(MetricError::WrongDatum {
                    id: numerator.id,
                    expected: "count",
                });
            };
            let EvidenceDatum::Count(denominator_count) = &denominator.datum else {
                return Err(MetricError::WrongDatum {
                    id: denominator.id,
                    expected: "count",
                });
            };
            exact_count_rate(*numerator_count, *denominator_count)?
        }
        MetricKind::Fraction => {
            let EvidenceDatum::Scalar(numerator_value) = &numerator.datum else {
                return Err(MetricError::WrongDatum {
                    id: numerator.id,
                    expected: "scalar",
                });
            };
            let EvidenceDatum::Scalar(denominator_value) = &denominator.datum else {
                return Err(MetricError::WrongDatum {
                    id: denominator.id,
                    expected: "scalar",
                });
            };
            scalar_fraction(*numerator_value, *denominator_value)?
        }
    };

    Ok(ReconstructedMetric {
        kind,
        value,
        numerator_evidence: numerator_role.reference,
        denominator_evidence: denominator_role.reference,
        identity: metric_identity(&request),
    })
}

/// Typed skip-fraction request for an edit trace.
#[must_use]
pub const fn edit_skip_metric(trace: &EditTrace) -> MetricRequest {
    MetricRequest {
        schema: MetricSchema::EditSkipDiagnostic,
        evaluator_semantics: trace.evaluator_semantics,
        proposal_semantics: trace.proposal_semantics,
        definition: MetricDefinition::Rate {
            numerator: MetricEvidenceRole::exact(
                trace.correct_skips_evidence,
                trace.id,
                "correct-skips",
                "count",
            ),
            denominator: MetricEvidenceRole::exact(
                trace.total_ops_evidence,
                trace.id,
                "total-ops",
                "count",
            ),
        },
    }
}

/// Typed candidate-remainder-rate request for a merge trial.
#[must_use]
pub const fn merge_conflict_metric(trial: &MergeTrial) -> MetricRequest {
    MetricRequest {
        schema: MetricSchema::MergeConflictDiagnostic,
        evaluator_semantics: trial.evaluator_semantics,
        proposal_semantics: trial.proposal_semantics,
        definition: MetricDefinition::Rate {
            numerator: MetricEvidenceRole::exact(
                trial.conflict_count_evidence,
                trial.id,
                "candidate-remainder-conflicts",
                "count",
            ),
            denominator: MetricEvidenceRole::exact(
                trial.total_merges_evidence,
                trial.id,
                "total-merges",
                "count",
            ),
        },
    }
}

/// Speedup `baseline / candidate` (Proposal 8 planner, Proposal 2 skip-yield).
/// Non-finite, negative, and zero-denominator inputs are refused.
pub fn speedup(baseline_cost: f64, candidate_cost: f64) -> Result<f64, MetricError> {
    if !baseline_cost.is_finite() {
        return Err(MetricError::NonFiniteInput { role: "baseline" });
    }
    if !candidate_cost.is_finite() {
        return Err(MetricError::NonFiniteInput { role: "candidate" });
    }
    if baseline_cost < 0.0 {
        return Err(MetricError::NegativeInput { role: "baseline" });
    }
    if candidate_cost < 0.0 {
        return Err(MetricError::NegativeInput { role: "candidate" });
    }
    if candidate_cost == 0.0 {
        return Err(MetricError::ZeroDenominator);
    }
    let quotient = baseline_cost / candidate_cost;
    if quotient.is_finite() {
        Ok(quotient)
    } else {
        Err(MetricError::NonFiniteResult)
    }
}

fn scalar_fraction(numerator: f64, denominator: f64) -> Result<f64, MetricError> {
    if !numerator.is_finite() {
        return Err(MetricError::NonFiniteInput { role: "numerator" });
    }
    if !denominator.is_finite() {
        return Err(MetricError::NonFiniteInput {
            role: "denominator",
        });
    }
    if numerator < 0.0 {
        return Err(MetricError::NegativeInput { role: "numerator" });
    }
    if denominator < 0.0 {
        return Err(MetricError::NegativeInput {
            role: "denominator",
        });
    }
    if denominator == 0.0 {
        return Err(MetricError::ZeroDenominator);
    }
    if numerator > denominator {
        return Err(MetricError::NumeratorExceedsDenominator);
    }
    let quotient = numerator / denominator;
    if quotient.is_finite() {
        Ok(quotient)
    } else {
        Err(MetricError::NonFiniteResult)
    }
}

/// Win rate: the fraction of true outcomes. Empty input is refused because it
/// provides no denominator evidence.
pub fn win_rate(outcomes: &[bool]) -> Result<f64, MetricError> {
    rate(
        outcomes.iter().filter(|&&outcome| outcome).count(),
        outcomes.len(),
    )
}

/// A generic count rate. Zero totals and counts above the total are refused.
pub fn rate(count: usize, total: usize) -> Result<f64, MetricError> {
    let count =
        u64::try_from(count).map_err(|_| MetricError::InexactCount { role: "numerator" })?;
    let total = u64::try_from(total).map_err(|_| MetricError::InexactCount {
        role: "denominator",
    })?;
    exact_count_rate(count, total)
}

#[allow(clippy::cast_precision_loss)]
fn exact_count_rate(count: u64, total: u64) -> Result<f64, MetricError> {
    if total == 0 {
        return Err(MetricError::ZeroDenominator);
    }
    if count > total {
        return Err(MetricError::NumeratorExceedsDenominator);
    }
    if count > MAX_EXACT_F64_COUNT {
        return Err(MetricError::InexactCount { role: "numerator" });
    }
    if total > MAX_EXACT_F64_COUNT {
        return Err(MetricError::InexactCount {
            role: "denominator",
        });
    }
    Ok(count as f64 / total as f64)
}

/// The candidate-remainder diagnostic rate of a synthetic merge fixture.
pub fn conflict_rate(trial: &MergeTrial) -> Result<f64, MetricError> {
    rate(trial.candidate_remainder_conflicts, trial.total_merges)
}

/// The accept-rate `accepts / attempts` (Proposal 9 certified speculation).
pub fn accept_rate(accepts: usize, attempts: usize) -> Result<f64, MetricError> {
    rate(accepts, attempts)
}

// -- Governance Rule 2 evaluator schemas -----------------------------------

/// Corpus dataset selected by an evaluator schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatasetKind {
    /// Conjugate-heat-transfer queries.
    QuerySet,
    /// Optimization design tasks.
    DesignTasks,
    /// Recorded edit traces.
    EditTraces,
    /// Manufactured-solution battery.
    MmsBattery,
    /// Retained merge trials.
    MergeTrials,
}

impl DatasetKind {
    const fn tag(self) -> u8 {
        match self {
            Self::QuerySet => 1,
            Self::DesignTasks => 2,
            Self::EditTraces => 3,
            Self::MmsBattery => 4,
            Self::MergeTrials => 5,
        }
    }

    /// Stable schema name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::QuerySet => "query_set",
            Self::DesignTasks => "design_tasks",
            Self::EditTraces => "edit_traces",
            Self::MmsBattery => "mms_battery",
            Self::MergeTrials => "merge_trials",
        }
    }
}

/// Formula owned by a proposal evaluator schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvaluatorFormula {
    /// One fixed two-role metric schema.
    Metric(MetricSchema),
    /// `(conflicts + escalations + refusals + type conflicts) / attempts`.
    MergeOutcomeRate,
}

impl EvaluatorFormula {
    const fn tag(self) -> u8 {
        match self {
            Self::Metric(schema) => schema.tag(),
            Self::MergeOutcomeRate => 255,
        }
    }
}

/// Closed set of addendum evaluator schemas. Each variant fixes its dataset,
/// formula, and authenticated evidence roles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProposalEvaluator {
    /// Proposal 8 planner comparison.
    Planner8,
    /// Proposal 1 adjoint comparison.
    Adjoint1,
    /// Proposal 2 certified skip comparison.
    CertifiedSkip2,
    /// Proposal D endpoint-guard comparison.
    GuardD,
    /// Proposal F robust-design comparison.
    RobustF,
    /// Proposal 9 speculation comparison.
    Speculation9,
    /// Proposal 10 merge outcome comparison.
    Merge10,
    /// Proposal A query-volume coverage.
    CoverageA,
}

impl ProposalEvaluator {
    const fn tag(self) -> u8 {
        match self {
            Self::Planner8 => 1,
            Self::Adjoint1 => 2,
            Self::CertifiedSkip2 => 3,
            Self::GuardD => 4,
            Self::RobustF => 5,
            Self::Speculation9 => 6,
            Self::Merge10 => 7,
            Self::CoverageA => 8,
        }
    }

    /// Exact formulas required by this evaluator.
    #[must_use]
    pub const fn formulas(self) -> &'static [EvaluatorFormula] {
        match self {
            Self::Planner8 => &[EvaluatorFormula::Metric(MetricSchema::PlannerSpeedup)],
            Self::Adjoint1 => &[EvaluatorFormula::Metric(MetricSchema::AdjointWinRate)],
            Self::CertifiedSkip2 => &[EvaluatorFormula::Metric(MetricSchema::CertifiedSkipSpeedup)],
            Self::GuardD => &[
                EvaluatorFormula::Metric(MetricSchema::GuardEndpointRate),
                EvaluatorFormula::Metric(MetricSchema::GuardRandomRate),
            ],
            Self::RobustF => &[EvaluatorFormula::Metric(
                MetricSchema::RobustNotDominatedRate,
            )],
            Self::Speculation9 => &[
                EvaluatorFormula::Metric(MetricSchema::SpeculationAcceptRate),
                EvaluatorFormula::Metric(MetricSchema::WarmStartSpeedup),
            ],
            Self::Merge10 => &[EvaluatorFormula::MergeOutcomeRate],
            Self::CoverageA => &[EvaluatorFormula::Metric(MetricSchema::QueryVolumeFraction)],
        }
    }

    /// Exact authenticated roles required for every selected dataset row.
    #[must_use]
    pub const fn required_roles(self) -> &'static [EvidenceRole] {
        match self {
            Self::Planner8 => &[
                EvidenceRole::PlannerBaselineCost,
                EvidenceRole::PlannerCandidateCost,
            ],
            Self::Adjoint1 => &[EvidenceRole::AdjointWins, EvidenceRole::AdjointComparisons],
            Self::CertifiedSkip2 => &[
                EvidenceRole::SkipBaselineCost,
                EvidenceRole::SkipCandidateCost,
            ],
            Self::GuardD => &[
                EvidenceRole::GuardEndpointCatches,
                EvidenceRole::GuardEndpointTrials,
                EvidenceRole::GuardRandomCatches,
                EvidenceRole::GuardRandomTrials,
            ],
            Self::RobustF => &[
                EvidenceRole::RobustNotDominated,
                EvidenceRole::RobustComparisons,
            ],
            Self::Speculation9 => &[
                EvidenceRole::SpeculationAccepts,
                EvidenceRole::SpeculationAttempts,
                EvidenceRole::ColdStartCost,
                EvidenceRole::WarmStartCost,
            ],
            Self::Merge10 => &[
                EvidenceRole::MergeConflicts,
                EvidenceRole::MergeEscalations,
                EvidenceRole::MergeRefusals,
                EvidenceRole::MergeTypeConflicts,
                EvidenceRole::MergeAttempts,
            ],
            Self::CoverageA => &[
                EvidenceRole::CoveredQueryVolume,
                EvidenceRole::TotalQueryVolume,
            ],
        }
    }
}

/// One proposal evaluator schema. A row is only actually instrumented when
/// [`evaluate_proposal`] returns [`ProposalEvaluation::Available`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InstrumentedProposal {
    /// The proposal id.
    pub proposal: &'static str,
    /// Typed dataset selected by the evaluator.
    pub dataset: DatasetKind,
    /// Typed evaluator, formula, and required-role schema.
    pub evaluator: ProposalEvaluator,
    /// The kill metric.
    pub kill_metric: &'static str,
    /// Exact evaluator algorithm/version for reconstructing the metric.
    pub evaluator_semantics: &'static str,
    /// Exact proposal/governance interpretation of the result.
    pub proposal_semantics: &'static str,
}

/// Independent locator for one proposal evaluator input. Quantity, units,
/// datum kind, evaluator semantics, and proposal semantics come from the
/// closed [`ProposalEvaluator`] / [`EvidenceRole`] schema; this manifest only
/// selects the exact retained identity allowed to occupy that role.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProposalEvidenceRef {
    /// Canonical proposal id.
    pub proposal: &'static str,
    /// Selected dataset row.
    pub subject_id: &'static str,
    /// Authenticated evaluator role.
    pub role: EvidenceRole,
    /// Independently retained exact evidence identity.
    pub reference: EvidenceRef,
}

/// The canonical proposal evaluator schemas. These are measurement contracts,
/// not proof that their required evidence exists.
const INSTRUMENTED: [InstrumentedProposal; 8] = [
    InstrumentedProposal {
        proposal: "8",
        dataset: DatasetKind::QuerySet,
        evaluator: ProposalEvaluator::Planner8,
        kill_metric: "planner >=2x vs baseline at equal certified accuracy",
        evaluator_semantics: "equal-accuracy-cost-ratio-v1",
        proposal_semantics: "addendum-proposal-8:planner-kill-v1",
    },
    InstrumentedProposal {
        proposal: "1",
        dataset: DatasetKind::DesignTasks,
        evaluator: ProposalEvaluator::Adjoint1,
        kill_metric: "adjoint beats derivative-free on >=70% of tasks",
        evaluator_semantics: "matched-budget-task-win-rate-v1",
        proposal_semantics: "addendum-proposal-1:adjoint-kill-v1",
    },
    InstrumentedProposal {
        proposal: "2",
        dataset: DatasetKind::EditTraces,
        evaluator: ProposalEvaluator::CertifiedSkip2,
        kill_metric: "certified skip-yield >=2x vs hash memoization",
        evaluator_semantics: "equal-edit-cost-ratio-v1",
        proposal_semantics: "addendum-proposal-2:skip-yield-kill-v1",
    },
    InstrumentedProposal {
        proposal: "D",
        dataset: DatasetKind::DesignTasks,
        evaluator: ProposalEvaluator::GuardD,
        kill_metric: "guard endpoint catch-rate > random-design catch-rate",
        evaluator_semantics: "matched-budget-endpoint-catch-rate-v1",
        proposal_semantics: "addendum-proposal-d:guard-kill-v1",
    },
    InstrumentedProposal {
        proposal: "F",
        dataset: DatasetKind::DesignTasks,
        evaluator: ProposalEvaluator::RobustF,
        kill_metric: "robust optima not dominated by nominal+safety on realized cost",
        evaluator_semantics: "retained-realized-cost-dominance-v1",
        proposal_semantics: "addendum-proposal-f:robustness-kill-v1",
    },
    InstrumentedProposal {
        proposal: "9",
        dataset: DatasetKind::MmsBattery,
        evaluator: ProposalEvaluator::Speculation9,
        kill_metric: "accept-rate >30% AND warm-start >=1.5x on the elliptic/FEEC battery",
        evaluator_semantics: "elliptic-accept-and-warm-start-v1",
        proposal_semantics: "addendum-proposal-9:speculation-kill-v1",
    },
    InstrumentedProposal {
        proposal: "10",
        dataset: DatasetKind::MergeTrials,
        evaluator: ProposalEvaluator::Merge10,
        kill_metric: "unresolved merge outcomes <25% on retained realistic traces",
        evaluator_semantics: "retained-merge-outcome-rate-v1",
        proposal_semantics: "addendum-proposal-10:merge-kill-v1",
    },
    InstrumentedProposal {
        proposal: "A",
        dataset: DatasetKind::QuerySet,
        evaluator: ProposalEvaluator::CoverageA,
        kill_metric: "RB-certified regions cover >=20% of wedge query volume",
        evaluator_semantics: "typed-query-volume-coverage-v1",
        proposal_semantics: "addendum-proposal-a:ladder-kill-v1",
    },
];

/// The declared evaluator schemas. Declaration alone is not instrumentation.
#[must_use]
pub fn instrumented_proposals() -> &'static [InstrumentedProposal] {
    &INSTRUMENTED
}

/// One exact evidence role absent from one selected dataset row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MissingEvidenceRole {
    /// Proposal whose evaluator is unavailable.
    pub proposal: &'static str,
    /// Selected dataset row.
    pub subject_id: &'static str,
    /// Missing authenticated role.
    pub role: EvidenceRole,
}

/// A fail-closed proposal-evaluation diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProposalRefusal {
    /// A mutable proposal row disagrees with its closed typed schema.
    SchemaMismatch {
        proposal: &'static str,
        field: &'static str,
    },
    /// The selected dataset repeats a subject id and would reuse evidence.
    DuplicateSubject {
        proposal: &'static str,
        subject_id: &'static str,
    },
    /// More than one record claims the same required role for a dataset row.
    AmbiguousRole {
        proposal: &'static str,
        subject_id: &'static str,
        role: EvidenceRole,
    },
    /// A selected retained record is structurally invalid.
    Evidence {
        proposal: &'static str,
        subject_id: &'static str,
        error: EvidenceError,
    },
    /// Evidence existed but could not reconstruct the fixed formula.
    Metric {
        proposal: &'static str,
        subject_id: &'static str,
        error: MetricError,
    },
}

impl core::fmt::Display for ProposalRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SchemaMismatch { proposal, field } => {
                write!(
                    f,
                    "proposal {proposal} mismatches evaluator schema field '{field}'"
                )
            }
            Self::DuplicateSubject {
                proposal,
                subject_id,
            } => write!(
                f,
                "proposal {proposal} dataset repeats subject '{subject_id}'"
            ),
            Self::AmbiguousRole {
                proposal,
                subject_id,
                role,
            } => write!(
                f,
                "proposal {proposal} row '{subject_id}' has ambiguous role '{}'",
                role.name()
            ),
            Self::Evidence {
                proposal,
                subject_id,
                error,
            } => write!(
                f,
                "proposal {proposal} row '{subject_id}' evidence refused: {error}"
            ),
            Self::Metric {
                proposal,
                subject_id,
                error,
            } => write!(
                f,
                "proposal {proposal} row '{subject_id}' metric refused: {error}"
            ),
        }
    }
}

/// One reconstructed proposal metric.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProposalMetric {
    /// Dataset row measured.
    pub subject_id: &'static str,
    /// Fixed evaluator formula.
    pub formula: EvaluatorFormula,
    /// Reconstructed finite result.
    pub value: f64,
    /// Canonical identity of every formula input.
    pub identity: ContentHash,
}

/// Availability of one proposal's real kill measurement.
#[derive(Debug, Clone, PartialEq)]
pub enum ProposalEvaluation {
    /// Every selected row had all required roles and every formula reconstructed.
    Available { metrics: Vec<ProposalMetric> },
    /// Required retained measurements do not yet exist.
    Unavailable { missing: Vec<MissingEvidenceRole> },
    /// Evidence or schema was present but ambiguous or invalid.
    Refused { diagnostics: Vec<ProposalRefusal> },
}

fn proposal_schema(proposal: ProposalEvaluator) -> &'static InstrumentedProposal {
    match proposal {
        ProposalEvaluator::Planner8 => &INSTRUMENTED[0],
        ProposalEvaluator::Adjoint1 => &INSTRUMENTED[1],
        ProposalEvaluator::CertifiedSkip2 => &INSTRUMENTED[2],
        ProposalEvaluator::GuardD => &INSTRUMENTED[3],
        ProposalEvaluator::RobustF => &INSTRUMENTED[4],
        ProposalEvaluator::Speculation9 => &INSTRUMENTED[5],
        ProposalEvaluator::Merge10 => &INSTRUMENTED[6],
        ProposalEvaluator::CoverageA => &INSTRUMENTED[7],
    }
}

fn proposal_schema_refusals(proposal: &InstrumentedProposal) -> Vec<ProposalRefusal> {
    let expected = proposal_schema(proposal.evaluator);
    let mut diagnostics = Vec::new();
    for (field, matches) in [
        ("proposal", proposal.proposal == expected.proposal),
        ("dataset", proposal.dataset == expected.dataset),
        ("kill-metric", proposal.kill_metric == expected.kill_metric),
        (
            "evaluator-semantics",
            proposal.evaluator_semantics == expected.evaluator_semantics,
        ),
        (
            "proposal-semantics",
            proposal.proposal_semantics == expected.proposal_semantics,
        ),
    ] {
        if !matches {
            diagnostics.push(ProposalRefusal::SchemaMismatch {
                proposal: proposal.proposal,
                field,
            });
        }
    }
    diagnostics
}

fn dataset_subjects(corpus: &BenchmarkCorpus<'_>, dataset: DatasetKind) -> Vec<&'static str> {
    match dataset {
        DatasetKind::QuerySet => corpus.query_set.iter().map(|row| row.id).collect(),
        DatasetKind::DesignTasks => corpus.design_tasks.iter().map(|row| row.id).collect(),
        DatasetKind::EditTraces => corpus.edit_traces.iter().map(|row| row.id).collect(),
        DatasetKind::MmsBattery => corpus.mms_battery.iter().map(|row| row.id).collect(),
        DatasetKind::MergeTrials => corpus.merge_trials.iter().map(|row| row.id).collect(),
    }
}

struct ResolvedProposalEvidence {
    record: EvidenceRecord,
    reference: EvidenceRef,
}

fn proposal_evidence(
    corpus: &BenchmarkCorpus<'_>,
    proposal: &InstrumentedProposal,
    subject_id: &'static str,
    role: EvidenceRole,
) -> Result<Option<ResolvedProposalEvidence>, ProposalRefusal> {
    let mut matches = corpus.proposal_evidence.iter().filter(|entry| {
        entry.proposal == proposal.proposal && entry.subject_id == subject_id && entry.role == role
    });
    let Some(entry) = matches.next() else {
        return Ok(None);
    };
    if matches.next().is_some() {
        return Err(ProposalRefusal::AmbiguousRole {
            proposal: proposal.proposal,
            subject_id,
            role,
        });
    }
    let record =
        resolve_evidence(corpus, entry.reference).map_err(|error| ProposalRefusal::Evidence {
            proposal: proposal.proposal,
            subject_id,
            error,
        })?;
    let units = role
        .proposal_units()
        .ok_or(ProposalRefusal::SchemaMismatch {
            proposal: proposal.proposal,
            field: "non-proposal-role",
        })?;
    require_context(
        &record,
        subject_id,
        role.name(),
        units,
        proposal.evaluator_semantics,
        proposal.proposal_semantics,
    )
    .map_err(|error| ProposalRefusal::Evidence {
        proposal: proposal.proposal,
        subject_id,
        error,
    })?;
    if record.role != role {
        return Err(ProposalRefusal::Evidence {
            proposal: proposal.proposal,
            subject_id,
            error: EvidenceError::ContextMismatch {
                id: record.id,
                field: "evidence-role",
            },
        });
    }
    Ok(Some(ResolvedProposalEvidence {
        record,
        reference: entry.reference,
    }))
}

fn proposal_metric_role(evidence: &ResolvedProposalEvidence) -> MetricEvidenceRole {
    MetricEvidenceRole::exact(
        evidence.reference,
        evidence.record.subject_id,
        evidence.record.quantity,
        evidence.record.units,
    )
}

fn evaluate_pair_formula(
    corpus: &BenchmarkCorpus<'_>,
    proposal: &InstrumentedProposal,
    subject_id: &'static str,
    schema: MetricSchema,
) -> Result<ProposalMetric, ProposalRefusal> {
    let (left_role, right_role) = schema.expected_roles();
    let left = proposal_evidence(corpus, proposal, subject_id, left_role)?.ok_or_else(|| {
        ProposalRefusal::Metric {
            proposal: proposal.proposal,
            subject_id,
            error: MetricError::MissingProposalRole { role: left_role },
        }
    })?;
    let right = proposal_evidence(corpus, proposal, subject_id, right_role)?.ok_or_else(|| {
        ProposalRefusal::Metric {
            proposal: proposal.proposal,
            subject_id,
            error: MetricError::MissingProposalRole { role: right_role },
        }
    })?;
    let left = proposal_metric_role(&left);
    let right = proposal_metric_role(&right);
    let definition = match schema.kind() {
        MetricKind::Speedup => MetricDefinition::Speedup {
            baseline: left,
            candidate: right,
        },
        MetricKind::Rate => MetricDefinition::Rate {
            numerator: left,
            denominator: right,
        },
        MetricKind::Fraction => MetricDefinition::Fraction {
            numerator: left,
            denominator: right,
        },
    };
    let metric = reconstruct_metric(
        corpus,
        MetricRequest {
            schema,
            evaluator_semantics: proposal.evaluator_semantics,
            proposal_semantics: proposal.proposal_semantics,
            definition,
        },
    )
    .map_err(|error| ProposalRefusal::Metric {
        proposal: proposal.proposal,
        subject_id,
        error,
    })?;
    Ok(ProposalMetric {
        subject_id,
        formula: EvaluatorFormula::Metric(schema),
        value: metric.value,
        identity: metric.identity,
    })
}

fn evaluate_merge_formula(
    corpus: &BenchmarkCorpus<'_>,
    proposal: &InstrumentedProposal,
    subject_id: &'static str,
) -> Result<ProposalMetric, ProposalRefusal> {
    let mut adverse = 0_u64;
    let mut inputs = Vec::new();
    for role in [
        EvidenceRole::MergeConflicts,
        EvidenceRole::MergeEscalations,
        EvidenceRole::MergeRefusals,
        EvidenceRole::MergeTypeConflicts,
    ] {
        let record = proposal_evidence(corpus, proposal, subject_id, role)?.ok_or_else(|| {
            ProposalRefusal::Metric {
                proposal: proposal.proposal,
                subject_id,
                error: MetricError::MissingProposalRole { role },
            }
        })?;
        let EvidenceDatum::Count(value) = &record.record.datum else {
            return Err(ProposalRefusal::Metric {
                proposal: proposal.proposal,
                subject_id,
                error: MetricError::WrongDatum {
                    id: record.record.id,
                    expected: "count",
                },
            });
        };
        adverse = adverse
            .checked_add(*value)
            .ok_or_else(|| ProposalRefusal::Metric {
                proposal: proposal.proposal,
                subject_id,
                error: MetricError::InexactCount { role: "numerator" },
            })?;
        inputs.push(record);
    }
    let attempts = proposal_evidence(corpus, proposal, subject_id, EvidenceRole::MergeAttempts)?
        .ok_or_else(|| ProposalRefusal::Metric {
            proposal: proposal.proposal,
            subject_id,
            error: MetricError::MissingProposalRole {
                role: EvidenceRole::MergeAttempts,
            },
        })?;
    let EvidenceDatum::Count(total) = &attempts.record.datum else {
        return Err(ProposalRefusal::Metric {
            proposal: proposal.proposal,
            subject_id,
            error: MetricError::WrongDatum {
                id: attempts.record.id,
                expected: "count",
            },
        });
    };
    inputs.push(attempts);
    let value = exact_count_rate(adverse, *total).map_err(|error| ProposalRefusal::Metric {
        proposal: proposal.proposal,
        subject_id,
        error,
    })?;
    Ok(ProposalMetric {
        subject_id,
        formula: EvaluatorFormula::MergeOutcomeRate,
        value,
        identity: merge_metric_identity(proposal, subject_id, &inputs),
    })
}

/// Reconstruct one proposal's fixed kill measurement. Missing evidence is
/// reported exactly and never upgraded to instrumentation.
#[must_use]
pub fn evaluate_proposal(
    corpus: &BenchmarkCorpus<'_>,
    proposal: &InstrumentedProposal,
) -> ProposalEvaluation {
    let diagnostics = proposal_schema_refusals(proposal);
    if !diagnostics.is_empty() {
        return ProposalEvaluation::Refused { diagnostics };
    }
    if let Err(error) = validate_evidence_registry(corpus) {
        return ProposalEvaluation::Refused {
            diagnostics: vec![ProposalRefusal::Evidence {
                proposal: proposal.proposal,
                subject_id: "<registry>",
                error,
            }],
        };
    }
    let subjects = dataset_subjects(corpus, proposal.dataset);
    if subjects.is_empty() {
        return ProposalEvaluation::Refused {
            diagnostics: vec![ProposalRefusal::SchemaMismatch {
                proposal: proposal.proposal,
                field: "dataset-empty",
            }],
        };
    }
    let mut seen_subjects = Vec::new();
    for subject_id in &subjects {
        if seen_subjects.contains(subject_id) {
            return ProposalEvaluation::Refused {
                diagnostics: vec![ProposalRefusal::DuplicateSubject {
                    proposal: proposal.proposal,
                    subject_id: *subject_id,
                }],
            };
        }
        seen_subjects.push(*subject_id);
    }
    let mut missing = Vec::new();
    for subject_id in &subjects {
        let subject_id = *subject_id;
        for role in proposal.evaluator.required_roles() {
            match proposal_evidence(corpus, proposal, subject_id, *role) {
                Ok(Some(_)) => {}
                Ok(None) => missing.push(MissingEvidenceRole {
                    proposal: proposal.proposal,
                    subject_id,
                    role: *role,
                }),
                Err(diagnostic) => {
                    return ProposalEvaluation::Refused {
                        diagnostics: vec![diagnostic],
                    };
                }
            }
        }
    }
    if !missing.is_empty() {
        return ProposalEvaluation::Unavailable { missing };
    }

    let mut metrics = Vec::new();
    for subject_id in subjects {
        for formula in proposal.evaluator.formulas() {
            let result = match formula {
                EvaluatorFormula::Metric(schema) => {
                    evaluate_pair_formula(corpus, proposal, subject_id, *schema)
                }
                EvaluatorFormula::MergeOutcomeRate => {
                    evaluate_merge_formula(corpus, proposal, subject_id)
                }
            };
            match result {
                Ok(metric) => metrics.push(metric),
                Err(diagnostic) => {
                    return ProposalEvaluation::Refused {
                        diagnostics: vec![diagnostic],
                    };
                }
            }
        }
    }
    ProposalEvaluation::Available { metrics }
}

// -- Determinism -----------------------------------------------------------

fn query_admission_node_hash(record: &EvidenceRecord) -> Result<ContentHash, EvidenceError> {
    if record.role != EvidenceRole::QueryReference {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "evidence-role",
        });
    }
    let EvidenceDatum::QueryReference {
        answer,
        tolerance,
        reference_cost,
        reference_cost_units,
        color,
        ..
    } = &record.datum
    else {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "query-reference-datum",
        });
    };
    let mut encoder = CanonicalEncoder::new();
    encoder.u64("schema-version", 1);
    encoder.str("evidence-id", record.id);
    encoder.str("subject-id", record.subject_id);
    encoder.str("quantity", record.quantity);
    encoder.str("units", record.units);
    encoder.str("evaluator-semantics", record.evaluator_semantics);
    encoder.str("proposal-semantics", record.proposal_semantics);
    encoder.u8("evidence-role", record.role.tag());
    encoder.u64("answer-bits", answer.to_bits());
    encoder.u64("tolerance-bits", tolerance.to_bits());
    encoder.u64("reference-cost-bits", reference_cost.to_bits());
    encoder.str("reference-cost-units", reference_cost_units);
    encoder.color(color);
    Ok(encoder.finish(QUERY_ADMISSION_NODE_DOMAIN))
}

fn evidence_semantic_digest(record: &EvidenceRecord) -> ContentHash {
    let mut encoder = CanonicalEncoder::new();
    encoder.u64(
        "schema-version",
        u64::from(EVIDENCE_IDENTITY_SCHEMA_VERSION),
    );
    encoder.str("subject-id", record.subject_id);
    encoder.str("quantity", record.quantity);
    encoder.str("units", record.units);
    encoder.str("evaluator-semantics", record.evaluator_semantics);
    encoder.str("proposal-semantics", record.proposal_semantics);
    encoder.u8("evidence-role", record.role.tag());
    match &record.datum {
        EvidenceDatum::Scalar(value) => {
            encoder.u8("datum-tag", 1);
            encoder.u64("scalar-bits", value.to_bits());
        }
        EvidenceDatum::Count(value) => {
            encoder.u8("datum-tag", 2);
            encoder.u64("count", *value);
        }
        EvidenceDatum::QueryReference {
            answer,
            tolerance,
            reference_cost,
            reference_cost_units,
            color,
            admission_receipt,
        } => {
            encoder.u8("datum-tag", 3);
            encoder.u64("answer-bits", answer.to_bits());
            encoder.u64("tolerance-bits", tolerance.to_bits());
            encoder.u64("reference-cost-bits", reference_cost.to_bits());
            encoder.str("reference-cost-units", reference_cost_units);
            encoder.color(color);
            encoder.optional_admission_receipt(admission_receipt.as_ref());
        }
    }
    encoder.finish(EVIDENCE_IDENTITY_DOMAIN)
}

struct CanonicalEncoder {
    bytes: Vec<u8>,
}

impl CanonicalEncoder {
    fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    fn frame(&mut self, bytes: &[u8]) {
        self.bytes
            .extend_from_slice(&(bytes.len() as u64).to_le_bytes());
        self.bytes.extend_from_slice(bytes);
    }

    fn marker(&mut self, marker: &str) {
        self.frame(marker.as_bytes());
    }

    fn str(&mut self, name: &str, value: &str) {
        self.frame(name.as_bytes());
        self.frame(value.as_bytes());
    }

    fn u64(&mut self, name: &str, value: u64) {
        self.frame(name.as_bytes());
        self.frame(&value.to_le_bytes());
    }

    fn u8(&mut self, name: &str, value: u8) {
        self.frame(name.as_bytes());
        self.frame(&[value]);
    }

    fn content_hash(&mut self, name: &str, value: ContentHash) {
        self.frame(name.as_bytes());
        self.frame(value.as_bytes());
    }

    fn color(&mut self, color: &Color) {
        self.marker("color");
        match color {
            Color::Verified { lo, hi } => {
                self.u8("color-tag", 1);
                self.u64("verified-lo-bits", lo.to_bits());
                self.u64("verified-hi-bits", hi.to_bits());
            }
            Color::Validated { regime, dataset } => {
                self.u8("color-tag", 2);
                self.str("validated-dataset", dataset);
                self.u64("validated-axis-count", regime.bounds().len() as u64);
                for (axis, (lo, hi)) in regime.bounds() {
                    self.marker("validated-axis");
                    self.str("axis", axis);
                    self.u64("lo-bits", lo.to_bits());
                    self.u64("hi-bits", hi.to_bits());
                }
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                self.u8("color-tag", 3);
                self.str("estimated-estimator", estimator);
                self.u64("estimated-dispersion-bits", dispersion.to_bits());
            }
        }
    }

    fn optional_admission_receipt(&mut self, receipt: Option<&AdmissionReceipt>) {
        match receipt {
            None => self.u8("admission-receipt-tag", 0),
            Some(receipt) => {
                self.u8("admission-receipt-tag", 1);
                self.content_hash("admission-node-hash", receipt.node_hash());
                self.u64(
                    "admission-row-schema-version",
                    u64::from(receipt.row_schema_version()),
                );
                self.u64(
                    "admission-color-algebra-version",
                    u64::from(receipt.color_algebra_version()),
                );
                self.content_hash("admission-policy-fingerprint", receipt.policy_fingerprint());
            }
        }
    }

    fn evidence_ref(&mut self, role: &str, reference: EvidenceRef) {
        self.marker(role);
        self.str("evidence-id", reference.id);
        self.content_hash("expected-digest", reference.expected_digest);
        self.u8("resolution-policy", reference.policy.tag());
    }

    fn evidence_role(&mut self, role: &str, evidence: MetricEvidenceRole) {
        self.marker(role);
        self.evidence_ref("reference", evidence.reference);
        self.str("subject-id", evidence.subject_id);
        self.str("quantity", evidence.quantity);
        self.str("units", evidence.units);
    }

    fn finish(self, domain: &str) -> ContentHash {
        hash_domain(domain, &self.bytes)
    }
}

fn metric_identity(request: &MetricRequest) -> ContentHash {
    let mut hasher = CanonicalEncoder::new();
    hasher.u64("schema-version", u64::from(METRIC_IDENTITY_SCHEMA_VERSION));
    hasher.str("evaluator-semantics", request.evaluator_semantics);
    hasher.str("proposal-semantics", request.proposal_semantics);
    hasher.u8("metric-schema", request.schema.tag());
    match request.definition {
        MetricDefinition::Speedup {
            baseline,
            candidate,
        } => {
            hasher.u8("metric-kind", 1);
            hasher.evidence_role("baseline", baseline);
            hasher.evidence_role("candidate", candidate);
        }
        MetricDefinition::Rate {
            numerator,
            denominator,
        } => {
            hasher.u8("metric-kind", 2);
            hasher.evidence_role("numerator", numerator);
            hasher.evidence_role("denominator", denominator);
        }
        MetricDefinition::Fraction {
            numerator,
            denominator,
        } => {
            hasher.u8("metric-kind", 3);
            hasher.evidence_role("numerator", numerator);
            hasher.evidence_role("denominator", denominator);
        }
    }
    hasher.finish(METRIC_IDENTITY_DOMAIN)
}

fn merge_metric_identity(
    proposal: &InstrumentedProposal,
    subject_id: &'static str,
    inputs: &[ResolvedProposalEvidence],
) -> ContentHash {
    let mut hasher = CanonicalEncoder::new();
    hasher.u64("schema-version", u64::from(METRIC_IDENTITY_SCHEMA_VERSION));
    hasher.str("proposal", proposal.proposal);
    hasher.u8("dataset", proposal.dataset.tag());
    hasher.u8("evaluator", proposal.evaluator.tag());
    hasher.u8("formula", EvaluatorFormula::MergeOutcomeRate.tag());
    hasher.str("subject-id", subject_id);
    hasher.str("evaluator-semantics", proposal.evaluator_semantics);
    hasher.str("proposal-semantics", proposal.proposal_semantics);
    hasher.u64("input-count", inputs.len() as u64);
    for evidence in inputs {
        hasher.marker("input");
        hasher.u8("evidence-role", evidence.record.role.tag());
        hasher.evidence_ref("evidence", evidence.reference);
    }
    hasher.finish(METRIC_IDENTITY_DOMAIN)
}

/// Schema-versioned, length-framed identity over every corpus semantic field.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn corpus_digest_for(corpus: &BenchmarkCorpus<'_>) -> ContentHash {
    let mut hasher = CanonicalEncoder::new();
    hasher.u64(
        "identity-schema-version",
        corpus.identity_schema_version as u64,
    );
    hasher.u64("corpus-version", corpus.version as u64);

    hasher.u64("query-count", corpus.query_set.len() as u64);
    for query in corpus.query_set {
        hasher.marker("query-row");
        hasher.str("id", query.id);
        hasher.str("qoi", query.qoi);
        hasher.str("units", query.units);
        hasher.u64("tolerance-bits", query.tolerance.to_bits());
        hasher.u64("reference-answer-bits", query.reference_answer.to_bits());
        hasher.u64("reference-cost-bits", query.reference_cost.to_bits());
        hasher.str("reference-cost-units", query.reference_cost_units);
        hasher.evidence_ref("reference-evidence", query.reference_evidence);
        hasher.str(
            "reference-evaluator-semantics",
            query.reference_evaluator_semantics,
        );
        hasher.str(
            "reference-proposal-semantics",
            query.reference_proposal_semantics,
        );
    }

    hasher.u64("design-task-count", corpus.design_tasks.len() as u64);
    for task in corpus.design_tasks {
        hasher.marker("design-task-row");
        hasher.str("id", task.id);
        hasher.u64("dimension", task.dimension as u64);
        hasher.u64("optimum-bits", task.optimum.to_bits());
    }

    hasher.u64("edit-trace-count", corpus.edit_traces.len() as u64);
    for trace in corpus.edit_traces {
        hasher.marker("edit-trace-row");
        hasher.str("id", trace.id);
        hasher.u64("total-ops", trace.total_ops as u64);
        hasher.u64("correct-skips", trace.correct_skips as u64);
        hasher.evidence_ref("correct-skips-evidence", trace.correct_skips_evidence);
        hasher.evidence_ref("total-ops-evidence", trace.total_ops_evidence);
        hasher.str("evaluator-semantics", trace.evaluator_semantics);
        hasher.str("proposal-semantics", trace.proposal_semantics);
    }

    hasher.u64("mms-count", corpus.mms_battery.len() as u64);
    for case in corpus.mms_battery {
        hasher.marker("mms-row");
        hasher.str("id", case.id);
        hasher.u64("exact-center-bits", case.exact_center.to_bits());
    }

    hasher.u64("merge-trial-count", corpus.merge_trials.len() as u64);
    for trial in corpus.merge_trials {
        hasher.marker("merge-trial-row");
        hasher.str("id", trial.id);
        hasher.u64("total-merges", trial.total_merges as u64);
        hasher.u64(
            "candidate-remainder-conflicts",
            trial.candidate_remainder_conflicts as u64,
        );
        hasher.evidence_ref("conflict-count-evidence", trial.conflict_count_evidence);
        hasher.evidence_ref("total-merges-evidence", trial.total_merges_evidence);
        hasher.str("evaluator-semantics", trial.evaluator_semantics);
        hasher.str("proposal-semantics", trial.proposal_semantics);
    }

    hasher.u64(
        "retained-evidence-count",
        corpus.retained_evidence.len() as u64,
    );
    for record in corpus.retained_evidence {
        hasher.marker("retained-evidence-row");
        hasher.str("id", record.id);
        hasher.str("subject-id", record.subject_id);
        hasher.str("quantity", record.quantity);
        hasher.str("units", record.units);
        hasher.str("evaluator-semantics", record.evaluator_semantics);
        hasher.str("proposal-semantics", record.proposal_semantics);
        hasher.u8("evidence-role", record.role.tag());
        match &record.datum {
            EvidenceDatum::Scalar(value) => {
                hasher.u8("datum-tag", 1);
                hasher.u64("scalar-bits", value.to_bits());
            }
            EvidenceDatum::Count(value) => {
                hasher.u8("datum-tag", 2);
                hasher.u64("count", *value);
            }
            EvidenceDatum::QueryReference {
                answer,
                tolerance,
                reference_cost,
                reference_cost_units,
                color,
                admission_receipt,
            } => {
                hasher.u8("datum-tag", 3);
                hasher.u64("answer-bits", answer.to_bits());
                hasher.u64("tolerance-bits", tolerance.to_bits());
                hasher.u64("reference-cost-bits", reference_cost.to_bits());
                hasher.str("reference-cost-units", reference_cost_units);
                hasher.color(color);
                hasher.optional_admission_receipt(admission_receipt.as_ref());
            }
        }
    }

    hasher.u64(
        "proposal-evidence-count",
        corpus.proposal_evidence.len() as u64,
    );
    for entry in corpus.proposal_evidence {
        hasher.marker("proposal-evidence-row");
        hasher.str("proposal", entry.proposal);
        hasher.str("subject-id", entry.subject_id);
        hasher.u8("evidence-role", entry.role.tag());
        hasher.evidence_ref("evidence", entry.reference);
    }

    hasher.u64(
        "instrumented-proposal-count",
        corpus.instrumented_proposals.len() as u64,
    );
    for proposal in corpus.instrumented_proposals {
        hasher.marker("instrumented-proposal-row");
        hasher.str("proposal", proposal.proposal);
        hasher.u8("dataset", proposal.dataset.tag());
        hasher.u8("evaluator", proposal.evaluator.tag());
        hasher.u64("formula-count", proposal.evaluator.formulas().len() as u64);
        for formula in proposal.evaluator.formulas() {
            hasher.marker("formula");
            hasher.u8("formula-tag", formula.tag());
        }
        hasher.u64(
            "required-role-count",
            proposal.evaluator.required_roles().len() as u64,
        );
        for role in proposal.evaluator.required_roles() {
            hasher.marker("required-role");
            hasher.u8("evidence-role", role.tag());
            hasher.str("canonical-quantity", role.name());
            hasher.u8("canonical-datum-tag", role.datum_tag());
            match role.proposal_units() {
                Some(units) => {
                    hasher.u8("canonical-units-tag", 1);
                    hasher.str("canonical-units", units);
                }
                None => hasher.u8("canonical-units-tag", 0),
            }
        }
        hasher.str("kill-metric", proposal.kill_metric);
        hasher.str("evaluator-semantics", proposal.evaluator_semantics);
        hasher.str("proposal-semantics", proposal.proposal_semantics);
    }
    hasher.finish(CORPUS_IDENTITY_DOMAIN)
}

/// Deterministic identity of the built-in corpus.
#[must_use]
pub fn corpus_digest() -> ContentHash {
    corpus_digest_for(&benchmark_corpus())
}

// -- Audit -----------------------------------------------------------------

/// The corpus completeness audit.
#[derive(Debug, Clone, PartialEq)]
pub struct CorpusAudit {
    /// The corpus version.
    pub version: u32,
    /// Number of proposal schemas with fully available retained measurements.
    pub instrumented: usize,
    /// Gaps (human-readable) — empty when the corpus is complete.
    pub gaps: Vec<String>,
}

impl CorpusAudit {
    /// Is the corpus complete (no gaps)?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }
}

fn dataset_len(corpus: &BenchmarkCorpus<'_>, dataset: DatasetKind) -> usize {
    match dataset {
        DatasetKind::QuerySet => corpus.query_set.len(),
        DatasetKind::DesignTasks => corpus.design_tasks.len(),
        DatasetKind::EditTraces => corpus.edit_traces.len(),
        DatasetKind::MmsBattery => corpus.mms_battery.len(),
        DatasetKind::MergeTrials => corpus.merge_trials.len(),
    }
}

fn note_row_id(gaps: &mut Vec<String>, seen: &mut Vec<&'static str>, kind: &str, id: &'static str) {
    if id.trim().is_empty() {
        gaps.push(format!("{kind} id is empty"));
    } else if seen.contains(&id) {
        gaps.push(format!("duplicate corpus row id '{id}'"));
    } else {
        seen.push(id);
    }
}

fn count_binding(
    corpus: &BenchmarkCorpus<'_>,
    reference: EvidenceRef,
    subject_id: &'static str,
    quantity: &'static str,
    evaluator_semantics: &'static str,
    proposal_semantics: &'static str,
    role: EvidenceRole,
    expected: usize,
) -> Result<(), EvidenceError> {
    let record = resolve_evidence(corpus, reference)?;
    require_context(
        &record,
        subject_id,
        quantity,
        "count",
        evaluator_semantics,
        proposal_semantics,
    )?;
    if record.role != role {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "evidence-role",
        });
    }
    let EvidenceDatum::Count(value) = &record.datum else {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "count-datum",
        });
    };
    if u64::try_from(expected).ok() != Some(*value) {
        return Err(EvidenceError::ContextMismatch {
            id: record.id,
            field: "count-value",
        });
    }
    Ok(())
}

/// Audit an arbitrary borrowed corpus view under the safe deny-all admission
/// policy.
///
/// Structurally valid Estimated query declarations can pass this audit, but
/// positive query evidence is refused because no admission authority is
/// available. Composition roots that hold the real authority must call
/// [`audit_corpus_with_verifier`] instead.
#[must_use]
pub fn audit_corpus(corpus: &BenchmarkCorpus<'_>) -> CorpusAudit {
    audit_corpus_with_verifier(corpus, &NoAdmissionVerifier)
}

/// Audit an arbitrary borrowed corpus view with an injected query-evidence
/// admission authority.
///
/// The verifier is consulted only for positive query colors after the retained
/// record, receipt context, and complete query binding have passed their local
/// fail-closed checks.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn audit_corpus_with_verifier(
    corpus: &BenchmarkCorpus<'_>,
    verifier: &dyn AdmissionVerifier,
) -> CorpusAudit {
    let mut gaps = Vec::new();
    if corpus.version != BENCHMARK_VERSION {
        gaps.push(format!(
            "unsupported corpus version {} (expected {BENCHMARK_VERSION})",
            corpus.version
        ));
    }
    if corpus.identity_schema_version != CORPUS_IDENTITY_SCHEMA_VERSION {
        gaps.push(format!(
            "unsupported identity schema {} (expected {CORPUS_IDENTITY_SCHEMA_VERSION})",
            corpus.identity_schema_version
        ));
    }
    if corpus.query_set.is_empty() {
        gaps.push("query set is empty".to_string());
    }
    if corpus.design_tasks.is_empty() {
        gaps.push("design tasks are empty".to_string());
    }
    if corpus.edit_traces.is_empty() {
        gaps.push("edit traces are empty".to_string());
    }
    if corpus.mms_battery.is_empty() {
        gaps.push("MMS battery is empty".to_string());
    }
    if corpus.merge_trials.is_empty() {
        gaps.push("merge trials are empty".to_string());
    }

    let mut row_ids = Vec::new();
    for query in corpus.query_set {
        note_row_id(&mut gaps, &mut row_ids, "query", query.id);
        for (field, value) in [
            ("qoi", query.qoi),
            ("units", query.units),
            ("reference cost units", query.reference_cost_units),
            ("reference evaluator", query.reference_evaluator_semantics),
            ("reference proposal", query.reference_proposal_semantics),
        ] {
            if value.trim().is_empty() {
                gaps.push(format!("query '{}' has empty {field}", query.id));
            }
        }
        if !query.tolerance.is_finite() || query.tolerance <= 0.0 {
            gaps.push(format!("query '{}' has invalid tolerance", query.id));
        }
        if !query.reference_answer.is_finite() {
            gaps.push(format!(
                "query '{}' has non-finite reference answer",
                query.id
            ));
        }
        if !query.reference_cost.is_finite() || query.reference_cost <= 0.0 {
            gaps.push(format!("query '{}' has invalid reference cost", query.id));
        }
        if let Err(error) = resolve_query_reference_with_verifier(corpus, query, verifier) {
            gaps.push(format!("query '{}': {error}", query.id));
        }
    }

    for task in corpus.design_tasks {
        note_row_id(&mut gaps, &mut row_ids, "design task", task.id);
        if task.dimension == 0 {
            gaps.push(format!("design task '{}' has zero dimension", task.id));
        }
        if !task.optimum.is_finite() {
            gaps.push(format!("design task '{}' has non-finite optimum", task.id));
        }
    }

    for trace in corpus.edit_traces {
        note_row_id(&mut gaps, &mut row_ids, "edit trace", trace.id);
        if trace.total_ops == 0 {
            gaps.push(format!("edit trace '{}' has zero total ops", trace.id));
        }
        if trace.correct_skips > trace.total_ops {
            gaps.push(format!("edit trace '{}' skips exceed total ops", trace.id));
        }
        for result in [
            count_binding(
                corpus,
                trace.correct_skips_evidence,
                trace.id,
                "correct-skips",
                trace.evaluator_semantics,
                trace.proposal_semantics,
                EvidenceRole::EditDiagnosticSkips,
                trace.correct_skips,
            ),
            count_binding(
                corpus,
                trace.total_ops_evidence,
                trace.id,
                "total-ops",
                trace.evaluator_semantics,
                trace.proposal_semantics,
                EvidenceRole::EditDiagnosticTotalOps,
                trace.total_ops,
            ),
        ] {
            if let Err(error) = result {
                gaps.push(format!("edit trace '{}': {error}", trace.id));
            }
        }
        if let Err(error) = reconstruct_metric(corpus, edit_skip_metric(trace)) {
            gaps.push(format!("edit trace '{}' metric: {error}", trace.id));
        }
    }

    for case in corpus.mms_battery {
        note_row_id(&mut gaps, &mut row_ids, "MMS case", case.id);
        if !case.exact_center.is_finite() {
            gaps.push(format!(
                "MMS case '{}' has non-finite exact center",
                case.id
            ));
        }
    }

    for trial in corpus.merge_trials {
        note_row_id(&mut gaps, &mut row_ids, "merge trial", trial.id);
        if trial.total_merges == 0 {
            gaps.push(format!("merge trial '{}' has zero total merges", trial.id));
        }
        if trial.candidate_remainder_conflicts > trial.total_merges {
            gaps.push(format!(
                "merge trial '{}' conflicts exceed total merges",
                trial.id
            ));
        }
        for result in [
            count_binding(
                corpus,
                trial.conflict_count_evidence,
                trial.id,
                "candidate-remainder-conflicts",
                trial.evaluator_semantics,
                trial.proposal_semantics,
                EvidenceRole::MergeDiagnosticConflicts,
                trial.candidate_remainder_conflicts,
            ),
            count_binding(
                corpus,
                trial.total_merges_evidence,
                trial.id,
                "total-merges",
                trial.evaluator_semantics,
                trial.proposal_semantics,
                EvidenceRole::MergeDiagnosticTotal,
                trial.total_merges,
            ),
        ] {
            if let Err(error) = result {
                gaps.push(format!("merge trial '{}': {error}", trial.id));
            }
        }
        if let Err(error) = reconstruct_metric(corpus, merge_conflict_metric(trial)) {
            gaps.push(format!("merge trial '{}' metric: {error}", trial.id));
        }
    }

    let mut evidence_ids = Vec::new();
    let mut evidence_digests = Vec::new();
    for record in corpus.retained_evidence {
        if evidence_ids.contains(&record.id) {
            gaps.push(format!("duplicate retained evidence id '{}'", record.id));
        } else {
            evidence_ids.push(record.id);
        }
        let digest = record.semantic_digest();
        if let Some((_, other_id)) = evidence_digests
            .iter()
            .find(|(seen_digest, _)| *seen_digest == digest)
        {
            gaps.push(format!(
                "duplicate retained evidence digest for '{}' and '{}'",
                other_id, record.id
            ));
        } else {
            evidence_digests.push((digest, record.id));
        }
        if let Err(error) = validate_evidence_record(record) {
            gaps.push(error.to_string());
        }
    }

    let mut proposal_manifest_keys = Vec::new();
    for entry in corpus.proposal_evidence {
        let key = (entry.proposal, entry.subject_id, entry.role);
        if proposal_manifest_keys.contains(&key) {
            gaps.push(format!(
                "duplicate proposal evidence role '{}' for {} row '{}'",
                entry.role.name(),
                entry.proposal,
                entry.subject_id
            ));
        } else {
            proposal_manifest_keys.push(key);
        }
        let Some(schema) = INSTRUMENTED
            .iter()
            .find(|proposal| proposal.proposal == entry.proposal)
        else {
            gaps.push(format!(
                "proposal evidence references unknown proposal '{}'",
                entry.proposal
            ));
            continue;
        };
        if !schema.evaluator.required_roles().contains(&entry.role) {
            gaps.push(format!(
                "proposal {} does not accept role '{}'",
                entry.proposal,
                entry.role.name()
            ));
        }
        if !dataset_subjects(corpus, schema.dataset).contains(&entry.subject_id) {
            gaps.push(format!(
                "proposal {} evidence references unknown dataset row '{}'",
                entry.proposal, entry.subject_id
            ));
        }
    }

    let mut proposal_ids = Vec::new();
    let mut instrumented = 0;
    if corpus.instrumented_proposals.len() != INSTRUMENTED.len() {
        gaps.push(format!(
            "proposal evaluator schema count is {} (expected {})",
            corpus.instrumented_proposals.len(),
            INSTRUMENTED.len()
        ));
    }
    for expected in &INSTRUMENTED {
        if !corpus
            .instrumented_proposals
            .iter()
            .any(|proposal| proposal.evaluator == expected.evaluator)
        {
            gaps.push(format!(
                "proposal {} evaluator schema is missing",
                expected.proposal
            ));
        }
    }
    for proposal in corpus.instrumented_proposals {
        if proposal_ids.contains(&proposal.proposal) {
            gaps.push(format!(
                "duplicate instrumented proposal '{}'",
                proposal.proposal
            ));
        } else {
            proposal_ids.push(proposal.proposal);
        }
        for (field, value) in [
            ("proposal", proposal.proposal),
            ("kill metric", proposal.kill_metric),
            ("evaluator semantics", proposal.evaluator_semantics),
            ("proposal semantics", proposal.proposal_semantics),
        ] {
            if value.trim().is_empty() {
                gaps.push(format!(
                    "instrumented proposal '{}' has empty {field}",
                    proposal.proposal
                ));
            }
        }
        if dataset_len(corpus, proposal.dataset) == 0 {
            gaps.push(format!(
                "proposal {} references empty dataset '{}'",
                proposal.proposal,
                proposal.dataset.name()
            ));
        }
        match evaluate_proposal(corpus, proposal) {
            ProposalEvaluation::Available { .. } => instrumented += 1,
            ProposalEvaluation::Unavailable { missing } => {
                for missing in missing {
                    gaps.push(format!(
                        "proposal {} unavailable for '{}': missing role '{}'",
                        missing.proposal,
                        missing.subject_id,
                        missing.role.name()
                    ));
                }
            }
            ProposalEvaluation::Refused { diagnostics } => {
                gaps.extend(
                    diagnostics
                        .into_iter()
                        .map(|diagnostic| diagnostic.to_string()),
                );
            }
        }
    }

    CorpusAudit {
        version: corpus.version,
        instrumented,
        gaps,
    }
}

/// Audit the built-in corpus.
#[must_use]
pub fn audit() -> CorpusAudit {
    audit_corpus(&benchmark_corpus())
}
