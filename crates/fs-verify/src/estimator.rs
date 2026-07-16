//! The VERIFIER: equilibrated-flux a-posteriori bounds (Prager–Synge,
//! 1D elliptic class), interval-evaluated to VERIFIED color.
//!
//! The rigor structure: ANY σ with `σ′ = −f` yields the guaranteed
//! bound `‖(u − u_h)′‖ ≤ ‖σ − u_h′‖` — so the free constant in
//! `σ = c − F` is optimized in plain f64 for TIGHTNESS while the bound
//! itself is evaluated with outward-rounded intervals over exact Gauss
//! quadrature (polynomial data ⇒ the quadrature identity is exact;
//! only rounding needs enclosing). Malformed inputs and unusable
//! enclosures FAIL CLOSED as structured refusals: no color, ever.

use crate::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS, MMS_PROBLEM_IDENTITY_VERSION,
    MmsProblem, gauss5, require_converged, true_energy_error, try_zeroed, validate_candidate,
    validate_problem,
};
use crate::interval::Iv;
use fs_blake3::{Blake3, ContentHash, hash_domain};
use fs_evidence::Color;
use std::fmt::Write as _;

/// Largest mesh admitted by the synchronous v0 verifier.
pub const MAX_VERIFIER_MESH_NODES: usize = MAX_FEM1D_MESH_NODES;
/// Exactness envelope for the manufactured solution: degree at most five.
pub const MAX_VERIFIER_POLY_COEFFICIENTS: usize = MAX_FEM1D_POLY_COEFFICIENTS;
/// Semantic version of the verifier's complete bounded-work accounting.
pub const VERIFIER_WORK_PLAN_VERSION: u32 = 2;
/// Semantic version of the verifier's callback/checkpoint schedule.
pub const VERIFIER_POLL_POLICY_VERSION: u32 = 1;
/// Maximum completed logical work between verifier work-boundary callbacks.
pub const VERIFIER_POLL_STRIDE_WORK_UNITS: u128 = 256;
/// Semantic version of the candidate-bound reconstructed-flux identity.
pub const VERIFIER_FLUX_IDENTITY_VERSION: u32 = 2;
/// Semantic schema for production verifier receipts.
pub const VERIFIER_RECEIPT_SCHEMA_VERSION: u32 = 1;
/// Exact theorem implemented by the production verifier receipt.
pub const VERIFIER_RECEIPT_THEOREM: &str = "prager-synge/equilibrated-flux/elliptic-1d/v1";
/// Exact outward-rounded arithmetic policy implemented by the verifier.
pub const VERIFIER_RECEIPT_ARITHMETIC: &str = "outward-rounded-f64-interval/gauss5-enclosed/v1";
/// Exact operator named by the verifier receipt.
pub const VERIFIER_RECEIPT_OPERATOR: &str = "poisson-1d/homogeneous-dirichlet/p1/v1";
/// Exact quantity certified by this verifier.
pub const VERIFIER_RECEIPT_QOI: &str = "fem1d-energy-error";
/// Semantic units of [`VERIFIER_RECEIPT_QOI`].
pub const VERIFIER_RECEIPT_UNITS: &str = "energy-norm";

const VERIFIER_RECEIPT_HASH_DOMAIN: &str = "fs-verify:verifier-receipt:v1";
const VERIFIER_RECEIPT_MAGIC: &[u8; 8] = b"FSVRCP01";
const MAX_VERIFIER_RECEIPT_CANONICAL_BYTES: usize = 16 * 1024;
const MAX_VERIFIER_RECEIPT_STRING_BYTES: usize = 1024;
const MAX_VERIFIER_RECEIPT_HYPOTHESES: usize = 16;
const VERIFIER_RECEIPT_HYPOTHESES: [&str; 4] = [
    "canonical degree-at-most-five manufactured class",
    "homogeneous Dirichlet endpoints",
    "finite strictly increasing mesh over [0,1]",
    "equilibrated flux has exact derivative minus forcing",
];

/// Owner-local verifier-receipt identity declaration consumed by
/// `xtask check-identities`.
#[allow(dead_code)]
pub const VERIFIER_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-verify:verifier-receipt",
    "version_const=VERIFIER_RECEIPT_SCHEMA_VERSION",
    "version=1",
    "domain=fs-verify:verifier-receipt:v1",
    "domain_const=VERIFIER_RECEIPT_HASH_DOMAIN",
    "encoder=VerifierReceipt::calculated_artifact_root",
    "encoder_helpers=VerifierReceipt::canonical_bytes_inner,receipt_push,receipt_push_string,receipt_push_hash,receipt_phase_tag,receipt_checkpoint_tag",
    "schema_constants=VERIFIER_RECEIPT_SCHEMA_VERSION,VERIFIER_RECEIPT_HASH_DOMAIN,VERIFIER_RECEIPT_MAGIC,MAX_VERIFIER_RECEIPT_CANONICAL_BYTES,MAX_VERIFIER_RECEIPT_STRING_BYTES,MAX_VERIFIER_RECEIPT_HYPOTHESES",
    "schema_functions=VerifierReceipt::from_successful_report,VerifierProducerSourceIdentity::current,current_verifier_feature_set,source_set_root,framed_parts_root,f64_sequence_root,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-verify:fem1d-mms-problem",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=VerifierReceipt,VerifierProducerSourceIdentity",
    "source_fields=VerifierReceipt.schema_version:semantic,VerifierReceipt.theorem:semantic,VerifierReceipt.producer:derived:hashed-exclusively-through-its-owned-VerifierProducerSourceIdentity-sub-fields,VerifierReceipt.problem_identity_version:semantic,VerifierReceipt.problem_root:semantic,VerifierReceipt.candidate_root:semantic,VerifierReceipt.mesh_root:semantic,VerifierReceipt.operator_root:semantic,VerifierReceipt.coefficient_root:semantic,VerifierReceipt.query_root:semantic,VerifierReceipt.qoi:semantic,VerifierReceipt.units:semantic,VerifierReceipt.flux_hash:semantic,VerifierReceipt.verifier_family:semantic,VerifierReceipt.arithmetic:semantic,VerifierReceipt.hypotheses:semantic,VerifierReceipt.bound_lo_bits:semantic,VerifierReceipt.bound_hi_bits:semantic,VerifierReceipt.tolerance_bits:semantic,VerifierReceipt.accepted:semantic,VerifierReceipt.work_plan:semantic,VerifierReceipt.observed_completed_work:semantic,VerifierReceipt.observed_planned_work:semantic,VerifierReceipt.final_phase:semantic,VerifierReceipt.final_checkpoint:semantic,VerifierReceipt.publication_observed:semantic,VerifierReceipt.artifact_root:derived:recomputed-from-canonical-fields,VerifierProducerSourceIdentity.crate_name:semantic,VerifierProducerSourceIdentity.crate_version:semantic,VerifierProducerSourceIdentity.features:semantic,VerifierProducerSourceIdentity.producer_source_root:semantic,VerifierProducerSourceIdentity.dependency_source_root:semantic,VerifierProducerSourceIdentity.workspace_manifest_root:semantic,VerifierProducerSourceIdentity.workspace_lock_root:semantic,VerifierProducerSourceIdentity.toolchain_root:semantic",
    "source_bindings=VerifierReceipt.schema_version>schema-version,VerifierReceipt.theorem>theorem,VerifierReceipt.problem_identity_version>problem-identity-version,VerifierReceipt.problem_root>problem-root,VerifierReceipt.candidate_root>candidate-root,VerifierReceipt.mesh_root>mesh-root,VerifierReceipt.operator_root>operator-root,VerifierReceipt.coefficient_root>coefficient-root,VerifierReceipt.query_root>query-root,VerifierReceipt.qoi>qoi,VerifierReceipt.units>units,VerifierReceipt.flux_hash>flux-hash,VerifierReceipt.verifier_family>verifier-family,VerifierReceipt.arithmetic>arithmetic,VerifierReceipt.hypotheses>ordered-hypotheses,VerifierReceipt.bound_lo_bits>bound-lo-bits,VerifierReceipt.bound_hi_bits>bound-hi-bits,VerifierReceipt.tolerance_bits>tolerance-bits,VerifierReceipt.accepted>accepted,VerifierReceipt.work_plan>work-plan,VerifierReceipt.observed_completed_work>observed-completed-work,VerifierReceipt.observed_planned_work>observed-planned-work,VerifierReceipt.final_phase>final-phase,VerifierReceipt.final_checkpoint>final-checkpoint,VerifierReceipt.publication_observed>publication-observed,VerifierProducerSourceIdentity.crate_name>producer-crate,VerifierProducerSourceIdentity.crate_version>producer-version,VerifierProducerSourceIdentity.features>producer-features,VerifierProducerSourceIdentity.producer_source_root>producer-source-root,VerifierProducerSourceIdentity.dependency_source_root>dependency-source-root,VerifierProducerSourceIdentity.workspace_manifest_root>workspace-manifest-root,VerifierProducerSourceIdentity.workspace_lock_root>workspace-lock-root,VerifierProducerSourceIdentity.toolchain_root>toolchain-root",
    "external_semantic_fields=canonical-magic,digest-domain",
    "semantic_fields=canonical-magic,digest-domain,schema-version,theorem,producer-crate,producer-version,producer-features,producer-source-root,dependency-source-root,workspace-manifest-root,workspace-lock-root,toolchain-root,problem-identity-version,problem-root,candidate-root,mesh-root,operator-root,coefficient-root,query-root,qoi,units,flux-hash,verifier-family,arithmetic,ordered-hypotheses,bound-lo-bits,bound-hi-bits,tolerance-bits,accepted,work-plan,observed-completed-work,observed-planned-work,final-phase,final-checkpoint,publication-observed",
    "excluded_fields=artifact-root:derived-recomputed-not-canonical-input,statement:derived-display-only,producer-label:derived-display-only,binary-artifact-identity:no-claim",
    "consumers=VerifierReceipt::artifact_root,PresentedVerifierReceipt::from_retained_bytes,admit_verifier_receipt,fs-ir::planner::VerifierCertificate,fs-flywheel-e2e",
    "mutations=canonical-magic:crates/fs-verify/src/estimator.rs#production_receipt_retention_requires_independent_root_and_replay,digest-domain:crates/fs-verify/src/estimator.rs#production_receipt_retention_requires_independent_root_and_replay,schema-version:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,theorem:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,producer-crate:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,producer-version:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,producer-features:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,producer-source-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,dependency-source-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,workspace-manifest-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,workspace-lock-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,toolchain-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,problem-identity-version:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,problem-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,candidate-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,mesh-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,operator-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,coefficient-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,query-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,qoi:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,units:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,flux-hash:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,verifier-family:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,arithmetic:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,ordered-hypotheses:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,bound-lo-bits:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,bound-hi-bits:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,tolerance-bits:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,accepted:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,work-plan:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,observed-completed-work:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,observed-planned-work:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,final-phase:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,final-checkpoint:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,publication-observed:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay",
    "nonsemantic_mutations=artifact-root:crates/fs-verify/src/estimator.rs#every_receipt_semantic_field_moves_root_and_fails_exact_replay,statement:crates/fs-verify/src/estimator.rs#production_receipt_retention_requires_independent_root_and_replay,producer-label:crates/fs-verify/src/estimator.rs#producer_identity_reports_the_exact_compiled_feature_set,binary-artifact-identity:crates/fs-verify/src/estimator.rs#producer_identity_reports_the_exact_compiled_feature_set",
    "field_guard=classify_verifier_receipt_identity_fields",
    "transport_guard=PresentedVerifierReceipt::from_retained_bytes",
    "version_guard=crates/fs-verify/src/estimator.rs#production_receipt_retention_requires_independent_root_and_replay",
    "coupling_surface=fs-verify:verifier-receipt",
];

/// One phase of the bounded equilibrated-flux verification workflow.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VerifierPhase {
    /// Input and derived-polynomial validation.
    Validation,
    /// Rounded optimizer for the free equilibrated-flux constant.
    Tightness,
    /// Outward-rounded equilibrated-bound construction.
    Equilibrated,
    /// Deterministic reconstructed-flux identity.
    Hash,
    /// Final report construction and publication.
    Finalization,
}

/// Why a verifier progress callback is being invoked.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VerifierCheckpointKind {
    /// Mandatory callback before a phase begins.
    PhaseEntry,
    /// Invocation-global multiple of [`VERIFIER_POLL_STRIDE_WORK_UNITS`].
    WorkBoundary,
    /// Mandatory callback after a structured refusal has inspected all work it reports.
    RefusalFlush,
    /// Mandatory final callback after the complete report is ready to publish.
    Publication,
}

/// Immutable invocation-global verifier progress passed to a callback by value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VerifierProgress {
    /// Callback reason.
    pub kind: VerifierCheckpointKind,
    /// Phase active at this callback.
    pub phase: VerifierPhase,
    /// Exact logical work completed across all phases.
    pub completed_work_units: u128,
    /// Complete constant-time preflighted work for this invocation.
    pub planned_work_units: u128,
}

/// Complete checked logical-work shape for one verifier invocation.
///
/// Counts are exact credited logical progress, not instruction counters. The
/// fixed-cap polynomial helpers process at most six coefficients atomically and
/// are credited only after success; this can conservatively lag physical work
/// by one such micro-tile but never overstates completed work.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VerifierWorkPlan {
    validation: u128,
    tightness: u128,
    equilibrated: u128,
    hash: u128,
    finalization: u128,
    total: u128,
}

impl VerifierWorkPlan {
    /// Preflight the complete work shape without running a callback.
    ///
    /// Shape refusals therefore have exact zero-work semantics. Content
    /// validation happens later through [`verify_with_checkpoint`].
    ///
    /// # Errors
    /// Returns a structured shape refusal for an inadmissible mesh, candidate,
    /// polynomial envelope, or checked work-plan overflow.
    pub fn for_inputs(problem: &MmsProblem, candidate: &[f64]) -> Result<Self, VerifierRefusal> {
        let mesh_nodes = problem.mesh().len();
        if !(2..=MAX_VERIFIER_MESH_NODES).contains(&mesh_nodes) {
            return Err(VerifierRefusal::MeshNodeCount);
        }
        if candidate.len() != mesh_nodes {
            return Err(VerifierRefusal::CandidateLength);
        }
        for (role, polynomial) in [
            (VerifierPolynomial::ExactSolution, problem.exact_solution()),
            (VerifierPolynomial::Forcing, problem.forcing()),
            (
                VerifierPolynomial::ForcingAntiderivative,
                problem.rounded_forcing_antiderivative(),
            ),
        ] {
            if !(1..=MAX_VERIFIER_POLY_COEFFICIENTS).contains(&polynomial.coefficients().len()) {
                return Err(VerifierRefusal::PolynomialCoefficientCount { polynomial: role });
            }
        }

        let mesh_nodes =
            u128::try_from(mesh_nodes).map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let cells = mesh_nodes
            .checked_sub(1)
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let exact_coefficients = u128::try_from(problem.exact_solution().coefficients().len())
            .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let forcing_coefficients = u128::try_from(problem.forcing().coefficients().len())
            .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let antiderivative_coefficients = u128::try_from(
            problem
                .rounded_forcing_antiderivative()
                .coefficients()
                .len(),
        )
        .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;

        let validation_work_units = 3_u128
            .checked_add(
                mesh_nodes
                    .checked_mul(2)
                    .ok_or(VerifierRefusal::WorkPlanOverflow)?,
            )
            .and_then(|work| work.checked_add(cells))
            .and_then(|work| work.checked_add(exact_coefficients.checked_mul(3)?))
            .and_then(|work| work.checked_add(forcing_coefficients.checked_mul(3)?))
            .and_then(|work| work.checked_add(antiderivative_coefficients.checked_mul(2)?))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let tightness_work_units = cells;
        let equilibrated_work_units = cells;
        let hash_work_units = 7_u128
            .checked_add(
                mesh_nodes
                    .checked_mul(2)
                    .ok_or(VerifierRefusal::WorkPlanOverflow)?,
            )
            .and_then(|work| work.checked_add(forcing_coefficients))
            .and_then(|work| work.checked_add(antiderivative_coefficients))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let finalization_work_units = 1;
        let planned_work_units = validation_work_units
            .checked_add(tightness_work_units)
            .and_then(|work| work.checked_add(equilibrated_work_units))
            .and_then(|work| work.checked_add(hash_work_units))
            .and_then(|work| work.checked_add(finalization_work_units))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        Ok(Self {
            validation: validation_work_units,
            tightness: tightness_work_units,
            equilibrated: equilibrated_work_units,
            hash: hash_work_units,
            finalization: finalization_work_units,
            total: planned_work_units,
        })
    }

    /// Validation work in the plan.
    #[must_use]
    pub const fn validation_work_units(self) -> u128 {
        self.validation
    }

    /// Tightness-optimizer work in the plan.
    #[must_use]
    pub const fn tightness_work_units(self) -> u128 {
        self.tightness
    }

    /// Equilibrated-bound work in the plan.
    #[must_use]
    pub const fn equilibrated_work_units(self) -> u128 {
        self.equilibrated
    }

    /// Flux-identity work in the plan.
    #[must_use]
    pub const fn hash_work_units(self) -> u128 {
        self.hash
    }

    /// Final report/publication work in the plan.
    #[must_use]
    pub const fn finalization_work_units(self) -> u128 {
        self.finalization
    }

    /// Total work in the plan.
    #[must_use]
    pub const fn planned_work_units(self) -> u128 {
        self.total
    }

    /// Stable phase counts used by downstream evidence identities.
    #[must_use]
    pub const fn identity_fields(self) -> [u128; 6] {
        [
            self.validation,
            self.tightness,
            self.equilibrated,
            self.hash,
            self.finalization,
            self.total,
        ]
    }
}

/// Estimator families (Proposal D's independence escalation needs at
/// least two registered per class).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimatorFamily {
    /// Equilibrated flux (guaranteed, constant-free — the verifier).
    EquilibratedFlux,
    /// Hierarchical (refined-mesh comparison — independent, NOT
    /// guaranteed; the falsifier's cross-check).
    Hierarchical,
}

/// Polynomial role carried by a structured verifier refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierPolynomial {
    /// Manufactured exact-solution metadata (`u`).
    ExactSolution,
    /// Canonical forcing (`f = -u''`).
    Forcing,
    /// Canonical zero-constant antiderivative of the forcing (`big_f`).
    ForcingAntiderivative,
}

/// Stable, structured reason why no verifier authority was issued.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierRefusal {
    /// Mesh length is outside the bounded `2..=MAX_VERIFIER_MESH_NODES` class.
    MeshNodeCount,
    /// One polynomial is empty or exceeds the degree-five exactness class.
    PolynomialCoefficientCount {
        /// Polynomial whose resource envelope was violated.
        polynomial: VerifierPolynomial,
    },
    /// Candidate and mesh lengths differ.
    CandidateLength,
    /// Tolerance is non-finite or non-positive.
    InvalidTolerance,
    /// Mesh endpoints are not canonical `+0.0` and `1.0`.
    MeshDomain,
    /// A mesh coordinate is non-finite or the mesh is not strictly increasing.
    MeshCoordinates,
    /// A candidate value is non-finite.
    CandidateNonFinite,
    /// Candidate endpoints are not canonical homogeneous `+0.0` values.
    CandidateBoundary,
    /// One polynomial contains a non-finite coefficient.
    PolynomialNonFinite {
        /// Polynomial containing the non-finite coefficient.
        polynomial: VerifierPolynomial,
    },
    /// The exact-solution polynomial does not vanish canonically at both ends.
    ExactSolutionBoundary,
    /// A public derived polynomial differs from the canonical value recomputed from `u`.
    DerivedPolynomialMismatch {
        /// Public derived polynomial that did not match its canonical value.
        polynomial: VerifierPolynomial,
    },
    /// The optional tightness constant could not be computed finitely.
    NonFiniteTightness,
    /// Interval construction produced a non-finite, reversed, or unusable enclosure.
    InvalidEnclosure,
    /// Complete verifier work-shape arithmetic overflowed.
    WorkPlanOverflow,
    /// Executed logical work did not match the complete preflighted plan.
    WorkPlanMismatch,
}

impl VerifierRefusal {
    /// Stable identifier for diagnostics and ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            Self::MeshNodeCount => "mesh-node-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-coefficient-count",
            Self::CandidateLength => "candidate-length",
            Self::InvalidTolerance => "invalid-tolerance",
            Self::MeshDomain => "mesh-domain",
            Self::MeshCoordinates => "mesh-coordinates",
            Self::CandidateNonFinite => "candidate-non-finite",
            Self::CandidateBoundary => "candidate-boundary",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-non-finite",
            Self::ExactSolutionBoundary => "exact-solution-boundary",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::Forcing,
            } => "derived-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "derived-big-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "derived-u-mismatch",
            Self::NonFiniteTightness => "non-finite-tightness",
            Self::InvalidEnclosure => "invalid-enclosure",
            Self::WorkPlanOverflow => "work-plan-overflow",
            Self::WorkPlanMismatch => "work-plan-mismatch",
        }
    }
}

impl EstimatorFamily {
    /// Stable id for ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            EstimatorFamily::EquilibratedFlux => "equilibrated-flux-1d",
            EstimatorFamily::Hierarchical => "hierarchical-h2",
        }
    }
}

/// The verifier's verdict on one candidate.
#[derive(Debug, Clone)]
pub struct VerifierReport {
    /// The certified error-bound enclosure (energy norm).
    pub bound: Iv,
    /// Accept ⟺ `bound.hi ≤ tolerance` for an admitted finite report.
    pub accept: bool,
    /// The verified color carried by an ACCEPT (`None` on reject or refusal).
    pub color: Option<Color>,
    /// The tolerance tested against (feeds the planner).
    pub tolerance: f64,
    /// Estimator family id.
    pub family: &'static str,
    /// FNV hash of the reconstructed flux (ledger identity).
    pub flux_hash: u64,
    /// Structured refusal (`None` only when a finite bound was produced).
    pub refusal: Option<VerifierRefusal>,
}

impl VerifierReport {
    /// The review-round-3 ledger row (structured, never stdout).
    #[must_use]
    pub fn to_row(&self, problem: &str, oracle_error: f64) -> String {
        let problem = json_escape(problem);
        let family = json_escape(self.family);
        let bound_lo = finite_scientific(self.bound.lo);
        let bound_hi = finite_scientific(self.bound.hi);
        let oracle = finite_scientific(oracle_error);
        let tolerance = finite_scientific(self.tolerance);
        let effectivity = if self.refusal.is_none()
            && oracle_error.is_finite()
            && oracle_error > 0.0
            && self.bound.hi.is_finite()
        {
            finite_fixed(self.bound.hi / oracle_error)
        } else if self.refusal.is_none() && oracle_error == 0.0 {
            "1.0000".to_string()
        } else {
            "null".to_string()
        };
        let refusal = self.refusal.map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", reason.id()),
        );
        let verdict = if self.refusal.is_some() {
            "refused"
        } else if self.accept {
            "accept"
        } else {
            "reject"
        };
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"problem\":\"{problem}\",\"estimator_family_id\":\"{}\",\
             \"flux_hash\":\"{:016X}\",\"bound_lo\":{bound_lo},\"bound_hi\":{bound_hi},\
             \"oracle_true_error\":{oracle},\"effectivity\":{effectivity},\
             \"verdict\":\"{verdict}\",\"tolerance\":{tolerance},\"refusal\":{refusal}}}",
            family, self.flux_hash,
        );
        s
    }
}

/// Collision-resistant address of one canonical production verifier receipt.
///
/// This is an fs-verify-owned nominal wrapper. A root authenticates bytes, not
/// scientific authority; callers must still pass the decoded receipt through
/// [`admit_verifier_receipt`] for an exact independent replay.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VerifierArtifactRoot(ContentHash);

impl VerifierArtifactRoot {
    /// Construct a presented root from exact retained bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(ContentHash(bytes))
    }

    /// Raw digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal transport.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }

    /// Workspace lower-layer digest value.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.0
    }
}

impl core::fmt::Display for VerifierArtifactRoot {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Display::fmt(&self.0, f)
    }
}

impl core::fmt::Debug for VerifierArtifactRoot {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "VerifierArtifactRoot({self})")
    }
}

/// Source-cone identity declared by the verifier producer.
///
/// This deliberately does not claim to be a binary attestation. It binds the
/// exact in-tree manifests and source bytes for fs-verify plus its production
/// fs-evidence/fs-obs/fs-blake3 dependency cone, along with workspace, lock,
/// feature, and toolchain inputs. External artifact attestation remains a
/// separate authority boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifierProducerSourceIdentity {
    crate_name: String,
    crate_version: String,
    features: String,
    producer_source_root: ContentHash,
    dependency_source_root: ContentHash,
    workspace_manifest_root: ContentHash,
    workspace_lock_root: ContentHash,
    toolchain_root: ContentHash,
}

impl VerifierProducerSourceIdentity {
    fn current() -> Result<Self, VerifierReceiptError> {
        Ok(Self {
            crate_name: try_owned("fs-verify")?,
            crate_version: try_owned(crate::VERSION)?,
            features: try_owned(current_verifier_feature_set())?,
            producer_source_root: source_set_root(
                "fs-verify/producer-source-cone/v1",
                &[
                    include_bytes!("../Cargo.toml"),
                    include_bytes!("lib.rs"),
                    include_bytes!("economics.rs"),
                    include_bytes!("estimator.rs"),
                    include_bytes!("fem1d.rs"),
                    include_bytes!("interval.rs"),
                    include_bytes!("zoo.rs"),
                ],
            ),
            dependency_source_root: source_set_root(
                "fs-verify/production-dependency-source-cone/v1",
                &[
                    include_bytes!("../../fs-evidence/Cargo.toml"),
                    include_bytes!("../../fs-evidence/src/lib.rs"),
                    include_bytes!("../../fs-evidence/src/admitted.rs"),
                    include_bytes!("../../fs-evidence/src/cards.rs"),
                    include_bytes!("../../fs-evidence/src/color.rs"),
                    include_bytes!("../../fs-evidence/src/discrepancy.rs"),
                    include_bytes!("../../fs-evidence/src/falsify.rs"),
                    include_bytes!("../../fs-evidence/src/vv.rs"),
                    include_bytes!("../../fs-evidence/src/vv/codec.rs"),
                    include_bytes!("../../fs-evidence/src/vv/model.rs"),
                    include_bytes!("../../fs-obs/Cargo.toml"),
                    include_bytes!("../../fs-obs/src/lib.rs"),
                    include_bytes!("../../fs-obs/src/ident.rs"),
                    include_bytes!("../../fs-blake3/Cargo.toml"),
                    include_bytes!("../../fs-blake3/src/lib.rs"),
                    include_bytes!("../../fs-blake3/src/identity.rs"),
                ],
            ),
            workspace_manifest_root: hash_domain(
                "fs-verify:workspace-manifest:v1",
                include_bytes!("../../../Cargo.toml"),
            ),
            workspace_lock_root: hash_domain(
                "fs-verify:workspace-lock:v1",
                include_bytes!("../../../Cargo.lock"),
            ),
            toolchain_root: hash_domain(
                "fs-verify:toolchain:v1",
                include_bytes!("../../../rust-toolchain.toml"),
            ),
        })
    }

    /// Producer crate name.
    #[must_use]
    pub fn crate_name(&self) -> &str {
        &self.crate_name
    }

    /// Producer crate version.
    #[must_use]
    pub fn crate_version(&self) -> &str {
        &self.crate_version
    }

    /// Exact enabled production feature set represented by the receipt API.
    #[must_use]
    pub fn features(&self) -> &str {
        &self.features
    }

    /// Root of the fs-verify manifest and complete source tree.
    #[must_use]
    pub const fn producer_source_root(&self) -> ContentHash {
        self.producer_source_root
    }

    /// Root of the complete production dependency source cone.
    #[must_use]
    pub const fn dependency_source_root(&self) -> ContentHash {
        self.dependency_source_root
    }

    /// Root of the workspace manifest governing the producer build.
    #[must_use]
    pub const fn workspace_manifest_root(&self) -> ContentHash {
        self.workspace_manifest_root
    }

    /// Root of the workspace dependency lock input.
    #[must_use]
    pub const fn workspace_lock_root(&self) -> ContentHash {
        self.workspace_lock_root
    }

    /// Root of the pinned toolchain input.
    #[must_use]
    pub const fn toolchain_root(&self) -> ContentHash {
        self.toolchain_root
    }

    /// Honest source-identity label; this is not a binary-attestation label.
    #[must_use]
    pub fn label(&self) -> String {
        format!("{}-source@{}", self.crate_name, self.crate_version)
    }
}

const fn current_verifier_feature_set() -> &'static str {
    if cfg!(feature = "certified-speculation") {
        "certified-speculation"
    } else {
        "none"
    }
}

/// Actual callback progress observed before a verifier attempt stopped.
///
/// Cancellation and partial progress are telemetry only. This type cannot be
/// converted into a verifier receipt or admitted scientific authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VerifierAttemptTelemetry {
    work_plan: [u128; 6],
    last_progress: Option<VerifierProgress>,
    publication_observed: bool,
}

impl VerifierAttemptTelemetry {
    /// Exact preflighted work shape.
    #[must_use]
    pub const fn work_plan(&self) -> [u128; 6] {
        self.work_plan
    }

    /// Last real verifier callback observed before cancellation.
    #[must_use]
    pub const fn last_progress(&self) -> Option<VerifierProgress> {
        self.last_progress
    }

    /// Whether the publication callback completed. Cancelled attempts are
    /// always false and cannot mint receipts.
    #[must_use]
    pub const fn publication_observed(&self) -> bool {
        self.publication_observed
    }
}

/// Fail-closed receipt production, transport, or replay error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerifierReceiptError {
    /// The scientific verifier refused the supplied inputs or enclosure.
    VerifierRefused(VerifierRefusal),
    /// A real cancellation callback stopped the attempt before publication.
    Cancelled(VerifierAttemptTelemetry),
    /// The verifier report and observed callback transcript disagreed.
    Protocol(&'static str),
    /// Receipt storage could not be reserved.
    AllocationFailed,
    /// Retained bytes exceeded the fixed transport cap.
    ReceiptTooLarge {
        /// Supplied byte count.
        requested: usize,
        /// Fixed receipt cap.
        cap: usize,
    },
    /// Retained bytes were malformed or non-canonical.
    MalformedRetained(&'static str),
    /// The independently supplied BLAKE3 root did not authenticate the bytes.
    ArtifactRootMismatch,
    /// The presented receipt differed from a fresh exact verifier replay.
    ReplayMismatch,
}

impl core::fmt::Display for VerifierReceiptError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::VerifierRefused(reason) => write!(f, "verifier refused: {}", reason.id()),
            Self::Cancelled(telemetry) => write!(
                f,
                "verifier cancelled before publication after {:?}",
                telemetry.last_progress
            ),
            Self::Protocol(stage) => write!(f, "verifier receipt protocol mismatch: {stage}"),
            Self::AllocationFailed => f.write_str("verifier receipt allocation failed"),
            Self::ReceiptTooLarge { requested, cap } => {
                write!(f, "verifier receipt has {requested} bytes above cap {cap}")
            }
            Self::MalformedRetained(stage) => {
                write!(f, "malformed retained verifier receipt: {stage}")
            }
            Self::ArtifactRootMismatch => f.write_str("retained verifier receipt root mismatch"),
            Self::ReplayMismatch => f.write_str("verifier receipt differs from exact replay"),
        }
    }
}

impl std::error::Error for VerifierReceiptError {}

/// Immutable production proof record emitted only by the real verifier.
///
/// Fields are private. This type is returned only by a successful in-process
/// verifier publication. Retained bytes decode into the distinct
/// [`PresentedVerifierReceipt`] transport type instead.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifierReceipt {
    schema_version: u32,
    theorem: String,
    producer: VerifierProducerSourceIdentity,
    problem_identity_version: u32,
    problem_root: ContentHash,
    candidate_root: ContentHash,
    mesh_root: ContentHash,
    operator_root: ContentHash,
    coefficient_root: ContentHash,
    query_root: ContentHash,
    qoi: String,
    units: String,
    flux_hash: u64,
    verifier_family: String,
    arithmetic: String,
    hypotheses: Vec<String>,
    bound_lo_bits: u64,
    bound_hi_bits: u64,
    tolerance_bits: u64,
    accepted: bool,
    work_plan: [u128; 6],
    observed_completed_work: u128,
    observed_planned_work: u128,
    final_phase: VerifierPhase,
    final_checkpoint: VerifierCheckpointKind,
    publication_observed: bool,
    artifact_root: VerifierArtifactRoot,
}

/// Root-authenticated retained verifier data without scientific authority.
///
/// The decoded production record is deliberately private. Presented transport
/// exposes identity fields needed to route and replay it, but not its claimed
/// bounds, acceptance bit, statement, color, or the underlying
/// [`VerifierReceipt`]. Exact independent replay through
/// [`admit_verifier_receipt`] is the only public route to those result fields.
///
/// ```compile_fail
/// use fs_verify::estimator::PresentedVerifierReceipt;
///
/// fn cannot_read_unadmitted_bounds(receipt: &PresentedVerifierReceipt) {
///     let _ = receipt.bound_lo();
///     let _ = receipt.bound_hi();
///     let _ = receipt.accepted();
/// }
/// ```
///
/// ```compile_fail
/// use fs_verify::estimator::PresentedVerifierReceipt;
///
/// fn cannot_mint_unadmitted_color(receipt: &PresentedVerifierReceipt) {
///     let _ = receipt.color();
/// }
/// ```
#[derive(Clone, PartialEq, Eq)]
pub struct PresentedVerifierReceipt {
    receipt: VerifierReceipt,
}

impl core::fmt::Debug for PresentedVerifierReceipt {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter
            .debug_tuple("PresentedVerifierReceipt")
            .field(&self.artifact_root())
            .finish()
    }
}

/// Opaque capability proving one exact receipt survived independent replay.
pub struct AdmittedVerifierReceipt<'a> {
    presented: &'a PresentedVerifierReceipt,
}

#[allow(dead_code)]
fn classify_verifier_receipt_identity_fields(receipt: &VerifierReceipt) {
    let VerifierReceipt {
        schema_version,
        theorem,
        producer,
        problem_identity_version,
        problem_root,
        candidate_root,
        mesh_root,
        operator_root,
        coefficient_root,
        query_root,
        qoi,
        units,
        flux_hash,
        verifier_family,
        arithmetic,
        hypotheses,
        bound_lo_bits,
        bound_hi_bits,
        tolerance_bits,
        accepted,
        work_plan,
        observed_completed_work,
        observed_planned_work,
        final_phase,
        final_checkpoint,
        publication_observed,
        artifact_root,
    } = receipt;
    let VerifierProducerSourceIdentity {
        crate_name,
        crate_version,
        features,
        producer_source_root,
        dependency_source_root,
        workspace_manifest_root,
        workspace_lock_root,
        toolchain_root,
    } = producer;
    let _ = (
        schema_version,
        theorem,
        crate_name,
        crate_version,
        features,
        producer_source_root,
        dependency_source_root,
        workspace_manifest_root,
        workspace_lock_root,
        toolchain_root,
        problem_identity_version,
        problem_root,
        candidate_root,
        mesh_root,
        operator_root,
        coefficient_root,
        query_root,
        qoi,
        units,
        flux_hash,
        verifier_family,
        arithmetic,
        hypotheses,
        bound_lo_bits,
        bound_hi_bits,
        tolerance_bits,
        accepted,
        work_plan,
        observed_completed_work,
        observed_planned_work,
        final_phase,
        final_checkpoint,
        publication_observed,
        artifact_root,
    );
}

fn try_owned(value: &str) -> Result<String, VerifierReceiptError> {
    let mut owned = String::new();
    owned
        .try_reserve_exact(value.len())
        .map_err(|_| VerifierReceiptError::AllocationFailed)?;
    owned.push_str(value);
    Ok(owned)
}

/// Hash an ordered, length-framed byte-part sequence without retaining the
/// potentially large source/candidate preimage. The frame prefix and exact
/// domain bytes make this use of plain streaming BLAKE3 unambiguous.
fn framed_parts_root(domain: &str, parts: &[&[u8]]) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(b"fs-verify:length-framed-parts:v1");
    hasher.update(&(domain.len() as u128).to_le_bytes());
    hasher.update(domain.as_bytes());
    hasher.update(&(parts.len() as u128).to_le_bytes());
    for part in parts {
        hasher.update(&(part.len() as u128).to_le_bytes());
        hasher.update(part);
    }
    hasher.finalize()
}

fn source_set_root(domain: &str, sources: &[&[u8]]) -> ContentHash {
    framed_parts_root(domain, sources)
}

fn f64_sequence_root(domain: &str, values: &[f64]) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(b"fs-verify:f64-sequence:v1");
    hasher.update(&(domain.len() as u128).to_le_bytes());
    hasher.update(domain.as_bytes());
    hasher.update(&(values.len() as u128).to_le_bytes());
    for value in values {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

fn receipt_phase_tag(phase: VerifierPhase) -> u8 {
    match phase {
        VerifierPhase::Validation => 0,
        VerifierPhase::Tightness => 1,
        VerifierPhase::Equilibrated => 2,
        VerifierPhase::Hash => 3,
        VerifierPhase::Finalization => 4,
    }
}

fn receipt_checkpoint_tag(checkpoint: VerifierCheckpointKind) -> u8 {
    match checkpoint {
        VerifierCheckpointKind::PhaseEntry => 0,
        VerifierCheckpointKind::WorkBoundary => 1,
        VerifierCheckpointKind::RefusalFlush => 2,
        VerifierCheckpointKind::Publication => 3,
    }
}

fn receipt_push(bytes: &mut Vec<u8>, value: &[u8]) -> Result<(), VerifierReceiptError> {
    let requested =
        bytes
            .len()
            .checked_add(value.len())
            .ok_or(VerifierReceiptError::ReceiptTooLarge {
                requested: usize::MAX,
                cap: MAX_VERIFIER_RECEIPT_CANONICAL_BYTES,
            })?;
    if requested > MAX_VERIFIER_RECEIPT_CANONICAL_BYTES {
        return Err(VerifierReceiptError::ReceiptTooLarge {
            requested,
            cap: MAX_VERIFIER_RECEIPT_CANONICAL_BYTES,
        });
    }
    bytes
        .try_reserve_exact(value.len())
        .map_err(|_| VerifierReceiptError::AllocationFailed)?;
    bytes.extend_from_slice(value);
    Ok(())
}

fn receipt_push_string(bytes: &mut Vec<u8>, value: &str) -> Result<(), VerifierReceiptError> {
    if value.len() > MAX_VERIFIER_RECEIPT_STRING_BYTES {
        return Err(VerifierReceiptError::ReceiptTooLarge {
            requested: value.len(),
            cap: MAX_VERIFIER_RECEIPT_STRING_BYTES,
        });
    }
    let len = u64::try_from(value.len())
        .map_err(|_| VerifierReceiptError::MalformedRetained("string length"))?;
    receipt_push(bytes, &len.to_le_bytes())?;
    receipt_push(bytes, value.as_bytes())
}

fn receipt_push_hash(bytes: &mut Vec<u8>, root: ContentHash) -> Result<(), VerifierReceiptError> {
    receipt_push(bytes, root.as_bytes())
}

struct ReceiptCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ReceiptCursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], VerifierReceiptError> {
        let end = self
            .offset
            .checked_add(count)
            .ok_or(VerifierReceiptError::MalformedRetained(
                "field length overflow",
            ))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(VerifierReceiptError::MalformedRetained("truncated field"))?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, VerifierReceiptError> {
        Ok(self.take(1)?[0])
    }

    fn bool(&mut self) -> Result<bool, VerifierReceiptError> {
        match self.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(VerifierReceiptError::MalformedRetained("boolean tag")),
        }
    }

    fn u32(&mut self) -> Result<u32, VerifierReceiptError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| VerifierReceiptError::MalformedRetained("u32"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, VerifierReceiptError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| VerifierReceiptError::MalformedRetained("u64"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn u128(&mut self) -> Result<u128, VerifierReceiptError> {
        let bytes: [u8; 16] = self
            .take(16)?
            .try_into()
            .map_err(|_| VerifierReceiptError::MalformedRetained("u128"))?;
        Ok(u128::from_le_bytes(bytes))
    }

    fn hash(&mut self) -> Result<ContentHash, VerifierReceiptError> {
        let bytes: [u8; 32] = self
            .take(32)?
            .try_into()
            .map_err(|_| VerifierReceiptError::MalformedRetained("content hash"))?;
        Ok(ContentHash(bytes))
    }

    fn string(&mut self) -> Result<String, VerifierReceiptError> {
        let len = usize::try_from(self.u64()?)
            .map_err(|_| VerifierReceiptError::MalformedRetained("string length"))?;
        if len > MAX_VERIFIER_RECEIPT_STRING_BYTES {
            return Err(VerifierReceiptError::MalformedRetained("oversized string"));
        }
        let value = core::str::from_utf8(self.take(len)?)
            .map_err(|_| VerifierReceiptError::MalformedRetained("string utf8"))?;
        try_owned(value)
    }

    fn phase(&mut self) -> Result<VerifierPhase, VerifierReceiptError> {
        match self.u8()? {
            0 => Ok(VerifierPhase::Validation),
            1 => Ok(VerifierPhase::Tightness),
            2 => Ok(VerifierPhase::Equilibrated),
            3 => Ok(VerifierPhase::Hash),
            4 => Ok(VerifierPhase::Finalization),
            _ => Err(VerifierReceiptError::MalformedRetained("phase tag")),
        }
    }

    fn checkpoint(&mut self) -> Result<VerifierCheckpointKind, VerifierReceiptError> {
        match self.u8()? {
            0 => Ok(VerifierCheckpointKind::PhaseEntry),
            1 => Ok(VerifierCheckpointKind::WorkBoundary),
            2 => Ok(VerifierCheckpointKind::RefusalFlush),
            3 => Ok(VerifierCheckpointKind::Publication),
            _ => Err(VerifierReceiptError::MalformedRetained("checkpoint tag")),
        }
    }

    fn finish(self) -> Result<(), VerifierReceiptError> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err(VerifierReceiptError::MalformedRetained("trailing bytes"))
        }
    }
}

impl VerifierReceipt {
    fn canonical_bytes_inner(&self) -> Result<Vec<u8>, VerifierReceiptError> {
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(2_048)
            .map_err(|_| VerifierReceiptError::AllocationFailed)?;
        receipt_push(&mut bytes, VERIFIER_RECEIPT_MAGIC)?;
        receipt_push(&mut bytes, &self.schema_version.to_le_bytes())?;
        receipt_push_string(&mut bytes, &self.theorem)?;
        receipt_push_string(&mut bytes, &self.producer.crate_name)?;
        receipt_push_string(&mut bytes, &self.producer.crate_version)?;
        receipt_push_string(&mut bytes, &self.producer.features)?;
        for root in [
            self.producer.producer_source_root,
            self.producer.dependency_source_root,
            self.producer.workspace_manifest_root,
            self.producer.workspace_lock_root,
            self.producer.toolchain_root,
        ] {
            receipt_push_hash(&mut bytes, root)?;
        }
        receipt_push(&mut bytes, &self.problem_identity_version.to_le_bytes())?;
        for root in [
            self.problem_root,
            self.candidate_root,
            self.mesh_root,
            self.operator_root,
            self.coefficient_root,
            self.query_root,
        ] {
            receipt_push_hash(&mut bytes, root)?;
        }
        receipt_push_string(&mut bytes, &self.qoi)?;
        receipt_push_string(&mut bytes, &self.units)?;
        receipt_push(&mut bytes, &self.flux_hash.to_le_bytes())?;
        receipt_push_string(&mut bytes, &self.verifier_family)?;
        receipt_push_string(&mut bytes, &self.arithmetic)?;
        let hypothesis_count = u32::try_from(self.hypotheses.len())
            .map_err(|_| VerifierReceiptError::MalformedRetained("hypothesis count"))?;
        if self.hypotheses.len() > MAX_VERIFIER_RECEIPT_HYPOTHESES {
            return Err(VerifierReceiptError::MalformedRetained("hypothesis count"));
        }
        receipt_push(&mut bytes, &hypothesis_count.to_le_bytes())?;
        for hypothesis in &self.hypotheses {
            receipt_push_string(&mut bytes, hypothesis)?;
        }
        for bits in [self.bound_lo_bits, self.bound_hi_bits, self.tolerance_bits] {
            receipt_push(&mut bytes, &bits.to_le_bytes())?;
        }
        receipt_push(&mut bytes, &[u8::from(self.accepted)])?;
        for work in self.work_plan {
            receipt_push(&mut bytes, &work.to_le_bytes())?;
        }
        receipt_push(&mut bytes, &self.observed_completed_work.to_le_bytes())?;
        receipt_push(&mut bytes, &self.observed_planned_work.to_le_bytes())?;
        receipt_push(&mut bytes, &[receipt_phase_tag(self.final_phase)])?;
        receipt_push(&mut bytes, &[receipt_checkpoint_tag(self.final_checkpoint)])?;
        receipt_push(&mut bytes, &[u8::from(self.publication_observed)])?;
        Ok(bytes)
    }

    fn calculated_artifact_root(&self) -> Result<VerifierArtifactRoot, VerifierReceiptError> {
        Ok(VerifierArtifactRoot(hash_domain(
            VERIFIER_RECEIPT_HASH_DOMAIN,
            &self.canonical_bytes_inner()?,
        )))
    }

    fn from_successful_report(
        problem: &MmsProblem,
        candidate: &[f64],
        tolerance: f64,
        plan: VerifierWorkPlan,
        progress: VerifierProgress,
        report: &VerifierReport,
    ) -> Result<Self, VerifierReceiptError> {
        if report.refusal.is_some() {
            return Err(VerifierReceiptError::Protocol("refused report"));
        }
        if !report.bound.lo.is_finite()
            || !report.bound.hi.is_finite()
            || report.bound.lo < 0.0
            || report.bound.lo > report.bound.hi
            || report.tolerance.to_bits() != tolerance.to_bits()
            || report.family != EstimatorFamily::EquilibratedFlux.id()
            || report.accept != (report.bound.hi <= tolerance)
        {
            return Err(VerifierReceiptError::Protocol("report fields"));
        }
        match (&report.color, report.accept) {
            // declared-color-ok: guarded multi-line pattern read validates report shape and constructs no positive evidence (6pf9)
            (Some(Color::Verified { lo, hi }), true)
                if lo.to_bits() == 0.0_f64.to_bits()
                    && hi.to_bits() == report.bound.hi.to_bits() => {}
            (None, false) => {}
            _ => return Err(VerifierReceiptError::Protocol("report color")),
        }
        if progress.kind != VerifierCheckpointKind::Publication
            || progress.phase != VerifierPhase::Finalization
            || progress.completed_work_units != plan.planned_work_units()
            || progress.planned_work_units != plan.planned_work_units()
        {
            return Err(VerifierReceiptError::Protocol("publication transcript"));
        }

        let producer = VerifierProducerSourceIdentity::current()?;
        let problem_root = hash_domain(
            "fs-verify:mms-problem-strong-root:v1",
            problem.canonical_bytes(),
        );
        let candidate_root = f64_sequence_root("fs-verify:candidate-nodal:v1", candidate);
        let mesh_root = f64_sequence_root("fs-verify:mesh:v1", problem.mesh());
        let operator_root = hash_domain(
            "fs-verify:operator:v1",
            VERIFIER_RECEIPT_OPERATOR.as_bytes(),
        );
        let coefficient_root = hash_domain(
            "fs-verify:mms-class-strong-root:v1",
            problem.class().canonical_bytes(),
        );
        let tolerance_bytes = tolerance.to_bits().to_le_bytes();
        let query_root = framed_parts_root(
            "fs-verify:verification-query:v1",
            &[
                problem_root.as_bytes(),
                candidate_root.as_bytes(),
                mesh_root.as_bytes(),
                operator_root.as_bytes(),
                coefficient_root.as_bytes(),
                &tolerance_bytes,
                VERIFIER_RECEIPT_QOI.as_bytes(),
                VERIFIER_RECEIPT_UNITS.as_bytes(),
            ],
        );
        let mut hypotheses = Vec::new();
        hypotheses
            .try_reserve_exact(VERIFIER_RECEIPT_HYPOTHESES.len())
            .map_err(|_| VerifierReceiptError::AllocationFailed)?;
        for hypothesis in VERIFIER_RECEIPT_HYPOTHESES {
            hypotheses.push(try_owned(hypothesis)?);
        }
        let mut receipt = Self {
            schema_version: VERIFIER_RECEIPT_SCHEMA_VERSION,
            theorem: try_owned(VERIFIER_RECEIPT_THEOREM)?,
            producer,
            problem_identity_version: MMS_PROBLEM_IDENTITY_VERSION,
            problem_root,
            candidate_root,
            mesh_root,
            operator_root,
            coefficient_root,
            query_root,
            qoi: try_owned(VERIFIER_RECEIPT_QOI)?,
            units: try_owned(VERIFIER_RECEIPT_UNITS)?,
            flux_hash: report.flux_hash,
            verifier_family: try_owned(EstimatorFamily::EquilibratedFlux.id())?,
            arithmetic: try_owned(VERIFIER_RECEIPT_ARITHMETIC)?,
            hypotheses,
            // This is the theorem's true-error interval. The interval
            // evaluator's positive lower enclosure is not a lower bound on
            // the unknown true error; zero is the only proved lower endpoint.
            bound_lo_bits: 0.0_f64.to_bits(),
            bound_hi_bits: report.bound.hi.to_bits(),
            tolerance_bits: tolerance.to_bits(),
            accepted: report.accept,
            work_plan: plan.identity_fields(),
            observed_completed_work: progress.completed_work_units,
            observed_planned_work: progress.planned_work_units,
            final_phase: progress.phase,
            final_checkpoint: progress.kind,
            publication_observed: true,
            artifact_root: VerifierArtifactRoot(ContentHash([0; 32])),
        };
        receipt.artifact_root = receipt.calculated_artifact_root()?;
        Ok(receipt)
    }

    /// Exact retained byte representation. The derived artifact root is not
    /// recursively encoded; it authenticates these bytes from outside.
    ///
    /// # Errors
    /// Returns a bounded allocation/size error. Receipt production has already
    /// exercised this path once, so a later failure grants no authority.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, VerifierReceiptError> {
        self.canonical_bytes_inner()
    }

    fn decode_retained_bytes(
        bytes: &[u8],
        expected_root: VerifierArtifactRoot,
    ) -> Result<Self, VerifierReceiptError> {
        if bytes.len() > MAX_VERIFIER_RECEIPT_CANONICAL_BYTES {
            return Err(VerifierReceiptError::ReceiptTooLarge {
                requested: bytes.len(),
                cap: MAX_VERIFIER_RECEIPT_CANONICAL_BYTES,
            });
        }
        let actual_root = VerifierArtifactRoot(hash_domain(VERIFIER_RECEIPT_HASH_DOMAIN, bytes));
        if actual_root != expected_root {
            return Err(VerifierReceiptError::ArtifactRootMismatch);
        }
        let mut cursor = ReceiptCursor::new(bytes);
        if cursor.take(VERIFIER_RECEIPT_MAGIC.len())? != VERIFIER_RECEIPT_MAGIC {
            return Err(VerifierReceiptError::MalformedRetained("magic"));
        }
        let schema_version = cursor.u32()?;
        if schema_version != VERIFIER_RECEIPT_SCHEMA_VERSION {
            return Err(VerifierReceiptError::MalformedRetained("schema version"));
        }
        let theorem = cursor.string()?;
        let producer = VerifierProducerSourceIdentity {
            crate_name: cursor.string()?,
            crate_version: cursor.string()?,
            features: cursor.string()?,
            producer_source_root: cursor.hash()?,
            dependency_source_root: cursor.hash()?,
            workspace_manifest_root: cursor.hash()?,
            workspace_lock_root: cursor.hash()?,
            toolchain_root: cursor.hash()?,
        };
        let problem_identity_version = cursor.u32()?;
        let problem_root = cursor.hash()?;
        let candidate_root = cursor.hash()?;
        let mesh_root = cursor.hash()?;
        let operator_root = cursor.hash()?;
        let coefficient_root = cursor.hash()?;
        let query_root = cursor.hash()?;
        let qoi = cursor.string()?;
        let units = cursor.string()?;
        let flux_hash = cursor.u64()?;
        let verifier_family = cursor.string()?;
        let arithmetic = cursor.string()?;
        let hypothesis_count = usize::try_from(cursor.u32()?)
            .map_err(|_| VerifierReceiptError::MalformedRetained("hypothesis count"))?;
        if hypothesis_count > MAX_VERIFIER_RECEIPT_HYPOTHESES {
            return Err(VerifierReceiptError::MalformedRetained("hypothesis count"));
        }
        let mut hypotheses = Vec::new();
        hypotheses
            .try_reserve_exact(hypothesis_count)
            .map_err(|_| VerifierReceiptError::AllocationFailed)?;
        for _ in 0..hypothesis_count {
            hypotheses.push(cursor.string()?);
        }
        let bound_lo_bits = cursor.u64()?;
        let bound_hi_bits = cursor.u64()?;
        let tolerance_bits = cursor.u64()?;
        let accepted = cursor.bool()?;
        let mut work_plan = [0_u128; 6];
        for work in &mut work_plan {
            *work = cursor.u128()?;
        }
        let observed_completed_work = cursor.u128()?;
        let observed_planned_work = cursor.u128()?;
        let final_phase = cursor.phase()?;
        let final_checkpoint = cursor.checkpoint()?;
        let publication_observed = cursor.bool()?;
        cursor.finish()?;
        let receipt = Self {
            schema_version,
            theorem,
            producer,
            problem_identity_version,
            problem_root,
            candidate_root,
            mesh_root,
            operator_root,
            coefficient_root,
            query_root,
            qoi,
            units,
            flux_hash,
            verifier_family,
            arithmetic,
            hypotheses,
            bound_lo_bits,
            bound_hi_bits,
            tolerance_bits,
            accepted,
            work_plan,
            observed_completed_work,
            observed_planned_work,
            final_phase,
            final_checkpoint,
            publication_observed,
            artifact_root: expected_root,
        };
        if receipt.canonical_bytes_inner()? != bytes {
            return Err(VerifierReceiptError::MalformedRetained(
                "non-canonical encoding",
            ));
        }
        Ok(receipt)
    }

    /// Receipt schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact theorem identifier.
    #[must_use]
    pub fn theorem(&self) -> &str {
        &self.theorem
    }

    /// Honest producer source-cone identity (not binary attestation).
    #[must_use]
    pub const fn producer(&self) -> &VerifierProducerSourceIdentity {
        &self.producer
    }

    /// Retained lower-layer problem-identity schema version.
    #[must_use]
    pub const fn problem_identity_version(&self) -> u32 {
        self.problem_identity_version
    }

    /// Strong root of the exact canonical manufactured problem.
    #[must_use]
    pub const fn problem_root(&self) -> ContentHash {
        self.problem_root
    }

    /// Strong root of exact candidate values.
    #[must_use]
    pub const fn candidate_root(&self) -> ContentHash {
        self.candidate_root
    }

    /// Strong root of exact mesh values.
    #[must_use]
    pub const fn mesh_root(&self) -> ContentHash {
        self.mesh_root
    }

    /// Exact operator identity.
    #[must_use]
    pub const fn operator_root(&self) -> ContentHash {
        self.operator_root
    }

    /// Strong root of the exact manufactured coefficient class.
    #[must_use]
    pub const fn coefficient_root(&self) -> ContentHash {
        self.coefficient_root
    }

    /// Exact verification query/QoI/tolerance identity.
    #[must_use]
    pub const fn query_root(&self) -> ContentHash {
        self.query_root
    }

    /// Certified quantity identifier.
    #[must_use]
    pub fn qoi(&self) -> &str {
        &self.qoi
    }

    /// Certified quantity units.
    #[must_use]
    pub fn units(&self) -> &str {
        &self.units
    }

    /// Candidate-bound reconstructed-flux identity.
    #[must_use]
    pub const fn flux_hash(&self) -> u64 {
        self.flux_hash
    }

    /// Production verifier family.
    #[must_use]
    pub fn verifier_family(&self) -> &str {
        &self.verifier_family
    }

    /// Arithmetic policy encoded in the receipt.
    #[must_use]
    pub fn arithmetic(&self) -> &str {
        &self.arithmetic
    }

    /// Exact theorem hypotheses.
    #[must_use]
    pub fn hypotheses(&self) -> &[String] {
        &self.hypotheses
    }

    /// Proved lower endpoint of the true-error interval.
    #[must_use]
    pub fn bound_lo(&self) -> f64 {
        f64::from_bits(self.bound_lo_bits)
    }

    /// Proved upper endpoint of the true-error interval.
    #[must_use]
    pub fn bound_hi(&self) -> f64 {
        f64::from_bits(self.bound_hi_bits)
    }

    /// Exact tested tolerance.
    #[must_use]
    pub fn tolerance(&self) -> f64 {
        f64::from_bits(self.tolerance_bits)
    }

    /// Whether the proved upper endpoint met the tested tolerance.
    #[must_use]
    pub const fn accepted(&self) -> bool {
        self.accepted
    }

    /// Exact preflighted logical work shape.
    #[must_use]
    pub const fn work_plan(&self) -> [u128; 6] {
        self.work_plan
    }

    /// Actual completed work at successful publication.
    #[must_use]
    pub const fn observed_completed_work(&self) -> u128 {
        self.observed_completed_work
    }

    /// Planned work reported by the successful publication callback.
    #[must_use]
    pub const fn observed_planned_work(&self) -> u128 {
        self.observed_planned_work
    }

    /// True only for a completed production publication transcript.
    #[must_use]
    pub const fn publication_observed(&self) -> bool {
        self.publication_observed
    }

    /// Collision-resistant address of the exact retained bytes.
    #[must_use]
    pub const fn artifact_root(&self) -> VerifierArtifactRoot {
        self.artifact_root
    }

    /// Lower-owned deterministic statement used by evidence packages.
    #[must_use]
    pub fn statement(&self) -> String {
        format!(
            "{} in {} is certified within [{:.17e}, {:.17e}] by {} at tolerance {:.17e}",
            self.qoi(),
            self.units(),
            self.bound_lo(),
            self.bound_hi(),
            self.theorem(),
            self.tolerance(),
        )
    }
}

impl PresentedVerifierReceipt {
    /// Parse exact retained bytes authenticated by an independently supplied
    /// collision-resistant root. This yields presented transport data only;
    /// callers must use [`admit_verifier_receipt`] before reading any claimed
    /// scientific result.
    ///
    /// # Errors
    /// Fails before decoding on a root mismatch, and otherwise rejects size,
    /// schema, UTF-8, enum-tag, boolean, truncation, and trailing-byte defects.
    pub fn from_retained_bytes(
        bytes: &[u8],
        expected_root: VerifierArtifactRoot,
    ) -> Result<Self, VerifierReceiptError> {
        Ok(Self {
            receipt: VerifierReceipt::decode_retained_bytes(bytes, expected_root)?,
        })
    }

    /// Collision-resistant identity of the exact retained bytes.
    #[must_use]
    pub const fn artifact_root(&self) -> VerifierArtifactRoot {
        self.receipt.artifact_root
    }

    /// Production verifier-family identity used for replay routing.
    #[must_use]
    pub fn verifier_family(&self) -> &str {
        &self.receipt.verifier_family
    }

    /// Honest producer source-cone identity (not binary attestation).
    #[must_use]
    pub const fn producer(&self) -> &VerifierProducerSourceIdentity {
        &self.receipt.producer
    }

    /// Strong identity of the exact canonical manufactured problem.
    #[must_use]
    pub const fn problem_root(&self) -> ContentHash {
        self.receipt.problem_root
    }

    /// Strong identity of exact candidate values.
    #[must_use]
    pub const fn candidate_root(&self) -> ContentHash {
        self.receipt.candidate_root
    }

    /// Strong identity of exact mesh values.
    #[must_use]
    pub const fn mesh_root(&self) -> ContentHash {
        self.receipt.mesh_root
    }

    /// Exact operator identity.
    #[must_use]
    pub const fn operator_root(&self) -> ContentHash {
        self.receipt.operator_root
    }

    /// Strong identity of the manufactured coefficient class.
    #[must_use]
    pub const fn coefficient_root(&self) -> ContentHash {
        self.receipt.coefficient_root
    }

    /// Exact verification-query identity.
    #[must_use]
    pub const fn query_root(&self) -> ContentHash {
        self.receipt.query_root
    }

    /// Exact preflighted logical work-shape identity.
    #[must_use]
    pub const fn work_plan(&self) -> [u128; 6] {
        self.receipt.work_plan
    }
}

impl<'a> AdmittedVerifierReceipt<'a> {
    fn authenticated(&self) -> &VerifierReceipt {
        &self.presented.receipt
    }

    /// Whether the replay proved that the upper endpoint met the tested
    /// tolerance.
    #[must_use]
    pub fn accepted(&self) -> bool {
        self.authenticated().accepted
    }

    /// Replay-authenticated lower endpoint of the true-error interval.
    #[must_use]
    pub fn bound_lo(&self) -> f64 {
        self.authenticated().bound_lo()
    }

    /// Replay-authenticated upper endpoint of the true-error interval.
    #[must_use]
    pub fn bound_hi(&self) -> f64 {
        self.authenticated().bound_hi()
    }

    /// Replay-authenticated quantity identifier.
    #[must_use]
    pub fn qoi(&self) -> &str {
        self.authenticated().qoi()
    }

    /// Replay-authenticated producer source-cone identity.
    #[must_use]
    pub fn producer(&self) -> &VerifierProducerSourceIdentity {
        self.authenticated().producer()
    }

    /// Replay-authenticated reconstructed-flux identity.
    #[must_use]
    pub fn flux_hash(&self) -> u64 {
        self.authenticated().flux_hash()
    }

    /// Replay-authenticated verifier-family identifier.
    #[must_use]
    pub fn verifier_family(&self) -> &str {
        self.authenticated().verifier_family()
    }

    /// Address of the exact retained bytes that survived replay.
    #[must_use]
    pub fn artifact_root(&self) -> VerifierArtifactRoot {
        self.authenticated().artifact_root()
    }

    /// Lower-owned deterministic statement, available only after replay.
    #[must_use]
    pub fn statement(&self) -> String {
        self.authenticated().statement()
    }

    /// Verified color only when the replay-authenticated receipt discharged
    /// its exact tested tolerance.
    #[must_use]
    pub fn color(&self) -> Option<Color> {
        // declared-color-ok: exact independent replay matched the presented receipt before this opaque admitted wrapper reconstructs its certified bound (6pf9)
        self.accepted().then(|| Color::Verified {
            lo: self.bound_lo(),
            hi: self.bound_hi(),
        })
    }
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

fn finite_scientific(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

fn finite_fixed(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.4}")
    } else {
        "null".to_string()
    }
}

fn fnv_extend(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

enum VerifierRunError<E> {
    Callback(E),
    Refusal(VerifierRefusal),
}

struct VerifierDriver<F> {
    callback: F,
    plan: VerifierWorkPlan,
    completed_work_units: u128,
    phase: VerifierPhase,
}

impl<F> VerifierDriver<F> {
    fn new(plan: VerifierWorkPlan, callback: F) -> Self {
        Self {
            callback,
            plan,
            completed_work_units: 0,
            phase: VerifierPhase::Validation,
        }
    }

    fn emit<E>(&mut self, kind: VerifierCheckpointKind) -> Result<(), E>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        (self.callback)(VerifierProgress {
            kind,
            phase: self.phase,
            completed_work_units: self.completed_work_units,
            planned_work_units: self.plan.total,
        })
    }

    fn enter<E>(&mut self, phase: VerifierPhase) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.phase = phase;
        self.emit(VerifierCheckpointKind::PhaseEntry)
            .map_err(VerifierRunError::Callback)
    }

    fn complete_one<E>(&mut self) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(1)
            .filter(|completed| *completed <= self.plan.total)
            .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
        if self
            .completed_work_units
            .is_multiple_of(VERIFIER_POLL_STRIDE_WORK_UNITS)
        {
            self.emit(VerifierCheckpointKind::WorkBoundary)
                .map_err(VerifierRunError::Callback)?;
        }
        Ok(())
    }

    fn refusal_flush<E>(&mut self) -> Result<(), E>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.emit(VerifierCheckpointKind::RefusalFlush)
    }

    fn require_completed<E>(&self, expected: u128) -> Result<(), VerifierRunError<E>> {
        if self.completed_work_units == expected {
            Ok(())
        } else {
            Err(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))
        }
    }

    fn publication<E>(&mut self) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.emit(VerifierCheckpointKind::Publication)
            .map_err(VerifierRunError::Callback)
    }
}

#[allow(clippy::too_many_lines)] // One auditable refusal order with exact partial-work semantics.
fn validate_inputs_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<(crate::fem1d::Poly, crate::fem1d::Poly), VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        let valid_count =
            (1..=MAX_VERIFIER_POLY_COEFFICIENTS).contains(&polynomial.coefficients().len());
        driver.complete_one()?;
        if !valid_count {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::PolynomialCoefficientCount { polynomial: role },
            ));
        }
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(VerifierRunError::Refusal(VerifierRefusal::InvalidTolerance));
    }
    if problem.mesh().first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || problem.mesh().last().map(|value| value.to_bits()) != Some(1.0_f64.to_bits())
    {
        return Err(VerifierRunError::Refusal(VerifierRefusal::MeshDomain));
    }
    for value in problem.mesh() {
        let finite = value.is_finite();
        driver.complete_one()?;
        if !finite {
            return Err(VerifierRunError::Refusal(VerifierRefusal::MeshCoordinates));
        }
    }
    for pair in problem.mesh().windows(2) {
        let increasing = pair[0] < pair[1];
        driver.complete_one()?;
        if !increasing {
            return Err(VerifierRunError::Refusal(VerifierRefusal::MeshCoordinates));
        }
    }
    for value in candidate {
        let finite = value.is_finite();
        driver.complete_one()?;
        if !finite {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::CandidateNonFinite,
            ));
        }
    }
    if candidate.first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || candidate.last().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
    {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::CandidateBoundary,
        ));
    }
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        for value in polynomial.coefficients() {
            let finite = value.is_finite();
            driver.complete_one()?;
            if !finite {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::PolynomialNonFinite { polynomial: role },
                ));
            }
        }
    }
    // The fixed 34-limb boundary certifier is an atomic U<=6 micro-tile. Credit
    // its coefficient units only after the result exists, so cancellation can
    // conservatively under-report by at most one tile but never over-report.
    let first_is_zero = problem
        .exact_solution()
        .coefficients()
        .first()
        .map(|value| value.to_bits())
        == Some(0.0_f64.to_bits());
    let zero_at_one = problem.exact_solution().is_exactly_zero_at_one();
    let exact_solution_has_boundary = first_is_zero && zero_at_one;
    for _ in problem.exact_solution().coefficients() {
        driver.complete_one()?;
    }
    if !exact_solution_has_boundary {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::ExactSolutionBoundary,
        ));
    }

    let expected_f = problem
        .exact_solution()
        .derive()
        .and_then(|derivative| derivative.derive())
        .map(crate::fem1d::Poly::neg)
        .map_err(|_| {
            VerifierRunError::Refusal(VerifierRefusal::PolynomialNonFinite {
                polynomial: VerifierPolynomial::Forcing,
            })
        })?;
    for _ in problem.exact_solution().coefficients() {
        driver.complete_one()?;
    }
    if problem.forcing().coefficients().len() != expected_f.coefficients().len() {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::Forcing,
            },
        ));
    }
    for (declared, expected) in problem
        .forcing()
        .coefficients()
        .iter()
        .zip(expected_f.coefficients())
    {
        let equal = declared.to_bits() == expected.to_bits();
        driver.complete_one()?;
        if !equal {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::DerivedPolynomialMismatch {
                    polynomial: VerifierPolynomial::Forcing,
                },
            ));
        }
    }
    let expected_big_f = expected_f.antiderive().map_err(|_| {
        VerifierRunError::Refusal(VerifierRefusal::PolynomialNonFinite {
            polynomial: VerifierPolynomial::ForcingAntiderivative,
        })
    })?;
    for _ in problem.forcing().coefficients() {
        driver.complete_one()?;
    }
    if problem
        .rounded_forcing_antiderivative()
        .coefficients()
        .len()
        != expected_big_f.coefficients().len()
    {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            },
        ));
    }
    for (declared, expected) in problem
        .rounded_forcing_antiderivative()
        .coefficients()
        .iter()
        .zip(expected_big_f.coefficients())
    {
        let equal = declared.to_bits() == expected.to_bits();
        driver.complete_one()?;
        if !equal {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::DerivedPolynomialMismatch {
                    polynomial: VerifierPolynomial::ForcingAntiderivative,
                },
            ));
        }
    }
    Ok((expected_f, expected_big_f))
}

fn tightness_constant_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    big_f: &crate::fem1d::Poly,
    driver: &mut VerifierDriver<F>,
) -> Result<f64, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut mean = 0.0;
    for element in 0..problem.mesh().len() - 1 {
        let (x0, x1) = (problem.mesh()[element], problem.mesh()[element + 1]);
        let h = x1 - x0;
        let slope = (candidate[element + 1] - candidate[element]) / h;
        if !h.is_finite() || h <= 0.0 || !slope.is_finite() {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::NonFiniteTightness,
            ));
        }
        for (point, weight) in gauss5(x0, x1) {
            let value = big_f.eval(point) + slope;
            let contribution = weight * value;
            if !point.is_finite()
                || !weight.is_finite()
                || !value.is_finite()
                || !contribution.is_finite()
            {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::NonFiniteTightness,
                ));
            }
            mean += contribution;
            if !mean.is_finite() {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::NonFiniteTightness,
                ));
            }
        }
        driver.complete_one()?;
    }
    Ok(mean)
}

fn finite_interval(interval: Iv) -> Result<Iv, VerifierRefusal> {
    if interval.lo.is_finite() && interval.hi.is_finite() && interval.lo <= interval.hi {
        Ok(interval)
    } else {
        Err(VerifierRefusal::InvalidEnclosure)
    }
}

fn interval_element_geometry(x0: f64, x1: f64) -> Result<(Iv, Iv, Iv), VerifierRefusal> {
    let x0 = Iv::point(x0);
    let x1 = Iv::point(x1);
    let h = finite_interval(x1.sub(x0))?;
    if h.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    let midpoint = finite_interval(x0.add(x1).scale_pos(0.5))?;
    let half = finite_interval(h.scale_pos(0.5))?;
    if half.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((h, midpoint, half))
}

fn interval_candidate_slope(first: f64, second: f64, h: Iv) -> Result<Iv, VerifierRefusal> {
    let difference = finite_interval(Iv::point(second).sub(Iv::point(first)))?;
    finite_interval(difference.div_pos(h))
}

fn interval_quadrature_geometry(
    midpoint: Iv,
    half: Iv,
    node_constant: f64,
    weight_constant: f64,
) -> Result<(Iv, Iv), VerifierRefusal> {
    let node = finite_interval(midpoint.add(half.mul(iv_c(node_constant))))?;
    let weight = finite_interval(half.mul(iv_c(weight_constant)))?;
    if weight.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((node, weight))
}

fn interval_antiderivative_coefficient(
    coefficient: f64,
    exponent: usize,
) -> Result<Iv, VerifierRefusal> {
    finite_interval(Iv::point(coefficient).div_pos(Iv::point(exponent as f64)))
}

fn interval_forcing_antiderivative(
    forcing: &crate::fem1d::Poly,
    x: Iv,
) -> Result<Iv, VerifierRefusal> {
    // F(x) = x * Horner(f_k / (k + 1)). Coefficient division is itself
    // intervalized: the rounded coefficients in `big_f` are replay metadata,
    // not point enclosures of the exact antiderivative of the authoritative f.
    let mut accumulated = Iv::zero();
    for (degree, coefficient) in forcing.coefficients().iter().copied().enumerate().rev() {
        let antiderivative_coefficient =
            interval_antiderivative_coefficient(coefficient, degree + 1)?;
        accumulated = finite_interval(accumulated.mul(x).add(antiderivative_coefficient))?;
    }
    finite_interval(x.mul(accumulated))
}

fn equilibrated_bound_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    forcing: &crate::fem1d::Poly,
    c_star: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<Iv, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut eta_sq = Iv::zero();
    for element in 0..problem.mesh().len() - 1 {
        let (h, midpoint, half) =
            interval_element_geometry(problem.mesh()[element], problem.mesh()[element + 1])
                .map_err(VerifierRunError::Refusal)?;
        let slope = interval_candidate_slope(candidate[element], candidate[element + 1], h)
            .map_err(VerifierRunError::Refusal)?;
        for (node_constant, weight_constant) in GAUSS5_REF {
            let (node, weight) =
                interval_quadrature_geometry(midpoint, half, node_constant, weight_constant)
                    .map_err(VerifierRunError::Refusal)?;
            let antiderivative = interval_forcing_antiderivative(forcing, node)
                .map_err(VerifierRunError::Refusal)?;
            let residual = finite_interval(Iv::point(c_star).sub(antiderivative).sub(slope))
                .map_err(VerifierRunError::Refusal)?;
            let contribution =
                finite_interval(weight.mul(residual.sq())).map_err(VerifierRunError::Refusal)?;
            eta_sq =
                finite_interval(eta_sq.add(contribution)).map_err(VerifierRunError::Refusal)?;
        }
        driver.complete_one()?;
    }
    let bound = finite_interval(eta_sq.sqrt()).map_err(VerifierRunError::Refusal)?;
    if bound.lo < 0.0 {
        Err(VerifierRunError::Refusal(VerifierRefusal::InvalidEnclosure))
    } else {
        Ok(bound)
    }
}

fn flux_hash_with_checkpoint<F, E>(
    c_star: f64,
    mesh: &[f64],
    candidate: &[f64],
    forcing: &crate::fem1d::Poly,
    antiderivative: &crate::fem1d::Poly,
    driver: &mut VerifierDriver<F>,
) -> Result<u64, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut hash = fnv_extend(
        0xcbf2_9ce4_8422_2325,
        b"fs-verify/equilibrated-flux-reconstruction/v2",
    );
    driver.complete_one()?;
    hash = fnv_extend(hash, &VERIFIER_FLUX_IDENTITY_VERSION.to_le_bytes());
    driver.complete_one()?;
    hash = fnv_extend(hash, &c_star.to_bits().to_le_bytes());
    driver.complete_one()?;
    for values in [mesh, candidate] {
        let length = u64::try_from(values.len())
            .map_err(|_| VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
        hash = fnv_extend(hash, &length.to_le_bytes());
        driver.complete_one()?;
        for value in values {
            hash = fnv_extend(hash, &value.to_bits().to_le_bytes());
            driver.complete_one()?;
        }
    }
    for polynomial in [forcing, antiderivative] {
        let length = u64::try_from(polynomial.coefficients().len())
            .map_err(|_| VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
        hash = fnv_extend(hash, &length.to_le_bytes());
        driver.complete_one()?;
        for coefficient in polynomial.coefficients() {
            hash = fnv_extend(hash, &coefficient.to_bits().to_le_bytes());
            driver.complete_one()?;
        }
    }
    Ok(hash)
}

fn refused(tolerance: f64, reason: VerifierRefusal) -> VerifierReport {
    VerifierReport {
        bound: Iv {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        },
        accept: false,
        color: None,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: 0,
        refusal: Some(reason),
    }
}

fn run_verifier<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<VerifierReport, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    driver.enter(VerifierPhase::Validation)?;
    let (canonical_f, canonical_big_f) =
        validate_inputs_with_checkpoint(problem, candidate, tolerance, driver)?;
    driver.require_completed(driver.plan.validation)?;

    driver.enter(VerifierPhase::Tightness)?;
    // Any finite c is sound. This rounded optimizer affects tightness only.
    let c_star = tightness_constant_with_checkpoint(problem, candidate, &canonical_big_f, driver)?;
    let after_tightness = driver
        .plan
        .validation
        .checked_add(driver.plan.tightness)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_tightness)?;

    driver.enter(VerifierPhase::Equilibrated)?;
    let bound =
        equilibrated_bound_with_checkpoint(problem, candidate, &canonical_f, c_star, driver)?;
    let after_equilibrated = after_tightness
        .checked_add(driver.plan.equilibrated)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_equilibrated)?;

    driver.enter(VerifierPhase::Hash)?;
    let flux_hash = flux_hash_with_checkpoint(
        c_star,
        problem.mesh(),
        candidate,
        &canonical_f,
        &canonical_big_f,
        driver,
    )?;
    let after_hash = after_equilibrated
        .checked_add(driver.plan.hash)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_hash)?;

    let accept = bound.hi <= tolerance;
    let color = if accept {
        // declared-color-ok: the rigorous equilibrated-flux accept declares a fresh report; retained authority requires exact receipt replay (6pf9)
        Some(Color::Verified {
            lo: 0.0,
            hi: bound.hi,
        })
    } else {
        None
    };
    driver.enter(VerifierPhase::Finalization)?;
    let report = VerifierReport {
        bound,
        accept,
        color,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash,
        refusal: None,
    };
    driver.complete_one()?;
    driver.require_completed(driver.plan.total)?;
    driver.publication()?;
    Ok(report)
}

/// The equilibrated-flux VERIFIER with an explicit sparse progress callback.
///
/// The callback runs at every phase entry, each invocation-global multiple of
/// [`VERIFIER_POLL_STRIDE_WORK_UNITS`], every structured-refusal flush, and the
/// final publication gate. Callback failure wins over any pending scientific
/// refusal and no report is returned. Shape refusals happen during constant-time
/// preflight and invoke no callback.
///
/// # Errors
/// Returns the callback's error unchanged. Scientific input or arithmetic
/// failures remain fail-closed [`VerifierReport`] values with a structured
/// [`VerifierReport::refusal`].
pub fn verify_with_checkpoint<E, F>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    callback: F,
) -> Result<VerifierReport, E>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let plan = match VerifierWorkPlan::for_inputs(problem, candidate) {
        Ok(plan) => plan,
        Err(reason) => return Ok(refused(tolerance, reason)),
    };
    let mut driver = VerifierDriver::new(plan, callback);
    match run_verifier(problem, candidate, tolerance, &mut driver) {
        Ok(report) => Ok(report),
        Err(VerifierRunError::Callback(error)) => Err(error),
        Err(VerifierRunError::Refusal(reason)) => {
            driver.refusal_flush()?;
            Ok(refused(tolerance, reason))
        }
    }
}

/// The equilibrated-flux VERIFIER: certify (or reject) a candidate's
/// nodal values against `tolerance`. The returned bound is a TRUE
/// upper bound on `‖(u − u_h)′‖` whenever the candidate satisfies the
/// boundary conditions; the enclosure is rigorous by outward rounding.
///
/// This convenience wrapper is bitwise equivalent to
/// [`verify_with_checkpoint`] with an infallible no-op callback.
#[must_use]
pub fn verify(problem: &MmsProblem, candidate: &[f64], tolerance: f64) -> VerifierReport {
    match verify_with_checkpoint(problem, candidate, tolerance, |_| {
        Ok::<(), core::convert::Infallible>(())
    }) {
        Ok(report) => report,
        Err(never) => match never {},
    }
}

/// Run the real verifier and emit an immutable production receipt only after
/// its successful publication callback. The cancellation predicate observes
/// actual verifier checkpoints; returning true stops the callback immediately
/// and yields telemetry rather than a promotable partial receipt.
///
/// # Errors
/// Returns a structured scientific refusal, real cancellation telemetry, or a
/// fail-closed receipt construction/protocol error. No receipt exists on any
/// error path.
pub fn verify_with_receipt_cancellable<F>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    mut should_cancel: F,
) -> Result<VerifierReceipt, VerifierReceiptError>
where
    F: FnMut(VerifierProgress) -> bool,
{
    let plan = VerifierWorkPlan::for_inputs(problem, candidate)
        .map_err(VerifierReceiptError::VerifierRefused)?;
    let mut last_progress = None;
    let mut publication_observed = false;
    let result = verify_with_checkpoint(problem, candidate, tolerance, |progress| {
        last_progress = Some(progress);
        if should_cancel(progress) {
            Err(())
        } else {
            if progress.kind == VerifierCheckpointKind::Publication {
                publication_observed = true;
            }
            Ok(())
        }
    });
    let report = match result {
        Ok(report) => report,
        Err(()) => {
            return Err(VerifierReceiptError::Cancelled(VerifierAttemptTelemetry {
                work_plan: plan.identity_fields(),
                last_progress,
                publication_observed: false,
            }));
        }
    };
    if let Some(reason) = report.refusal {
        return Err(VerifierReceiptError::VerifierRefused(reason));
    }
    if !publication_observed {
        return Err(VerifierReceiptError::Protocol(
            "missing successful publication callback",
        ));
    }
    let progress = last_progress.ok_or(VerifierReceiptError::Protocol(
        "missing verifier callback transcript",
    ))?;
    VerifierReceipt::from_successful_report(problem, candidate, tolerance, plan, progress, &report)
}

/// Run the real verifier to completion and return its exact production-owned
/// receipt. This is the receipt-producing counterpart to [`verify`].
///
/// # Errors
/// Returns a structured verifier refusal or fail-closed receipt construction
/// error. Above-tolerance finite reports still produce an unaccepted receipt;
/// after exact replay they yield `None` from
/// [`AdmittedVerifierReceipt::color`].
pub fn verify_with_receipt(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
) -> Result<VerifierReceipt, VerifierReceiptError> {
    verify_with_receipt_cancellable(problem, candidate, tolerance, |_| false)
}

/// Admit a presented verifier receipt only after authenticating its retained
/// bytes and independently replaying the exact production verifier over the
/// supplied problem, candidate, and tolerance.
///
/// # Errors
/// Returns [`VerifierReceiptError::ArtifactRootMismatch`] when the retained
/// fields do not match the stored root, or a replay/refusal error when current
/// production output differs from the presented receipt.
pub fn admit_verifier_receipt<'a>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    receipt: &'a PresentedVerifierReceipt,
) -> Result<AdmittedVerifierReceipt<'a>, VerifierReceiptError> {
    if receipt.receipt.calculated_artifact_root()? != receipt.receipt.artifact_root {
        return Err(VerifierReceiptError::ArtifactRootMismatch);
    }
    let replay = verify_with_receipt(problem, candidate, tolerance)?;
    if replay != receipt.receipt {
        return Err(VerifierReceiptError::ReplayMismatch);
    }
    Ok(AdmittedVerifierReceipt { presented: receipt })
}

const GAUSS5_REF: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_08),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.0, 0.568_888_888_888_888_9),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.906_179_845_938_664, 0.236_926_885_056_189_08),
];

/// One-ulp-widened constant (the tabulated Gauss data carries ~1 ulp
/// of transcription error; widening keeps enclosures honest).
fn iv_c(v: f64) -> Iv {
    Iv {
        lo: crate::interval::down(v),
        hi: crate::interval::up(v),
    }
}

/// The INDEPENDENT second family: hierarchical estimate from a
/// uniformly refined solve (`h/2`). Not guaranteed — the falsifier's
/// cross-check, never a color source.
///
/// # Errors
/// Returns [`Fem1dError`] for malformed inputs, refinement overflow/resource
/// excess, allocation failure, or a non-finite estimate.
pub fn hierarchical_estimate(problem: &MmsProblem, candidate: &[f64]) -> Result<f64, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let fine_nodes = problem
        .mesh()
        .len()
        .checked_mul(2)
        .and_then(|nodes| nodes.checked_sub(1))
        .ok_or(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: usize::MAX,
            limit: MAX_FEM1D_MESH_NODES,
        })?;
    if fine_nodes > MAX_FEM1D_MESH_NODES {
        return Err(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: fine_nodes,
            limit: MAX_FEM1D_MESH_NODES,
        });
    }
    let mut fine_mesh = Vec::new();
    fine_mesh
        .try_reserve_exact(fine_nodes)
        .map_err(|_| Fem1dError::AllocationFailed {
            stage: "hierarchical mesh",
            requested: fine_nodes,
        })?;
    for w in problem.mesh().windows(2) {
        fine_mesh.push(w[0]);
        fine_mesh.push(f64::midpoint(w[0], w[1]));
    }
    fine_mesh.push(problem.mesh()[problem.mesh().len() - 1]);
    let fine = problem.with_mesh(fine_mesh)?;
    let fine_u = crate::fem1d::solve_p1(&fine)?;
    // ‖u_{h/2}′ − u_h′‖ over the fine mesh.
    let mut acc = 0.0;
    for e in 0..fine.mesh().len() - 1 {
        let (x0, x1) = (fine.mesh()[e], fine.mesh()[e + 1]);
        let h = x1 - x0;
        let fine_slope = (fine_u[e + 1] - fine_u[e]) / h;
        // The coarse element containing this fine element.
        let coarse_e = e / 2;
        let ch = problem.mesh()[coarse_e + 1] - problem.mesh()[coarse_e];
        let coarse_slope = (candidate[coarse_e + 1] - candidate[coarse_e]) / ch;
        let d = fine_slope - coarse_slope;
        let updated = (h * d).mul_add(d, acc);
        if !fine_slope.is_finite()
            || !coarse_slope.is_finite()
            || !d.is_finite()
            || !updated.is_finite()
        {
            return Err(Fem1dError::NonFiniteIntermediate {
                stage: "hierarchical estimate",
                index: Some(e),
            });
        }
        acc = updated;
    }
    let estimate = acc.sqrt();
    if estimate.is_finite() {
        Ok(estimate)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "hierarchical estimate",
            index: None,
        })
    }
}

/// The nonlinear WARM-START fallback: the candidate is accepted only
/// as a starting point; the measured value is iteration savings and
/// the color is ESTIMATED, never verified (the honest R1 boundary).
#[derive(Debug, Clone)]
pub struct WarmStartReport {
    /// Newton iterations from a cold start (zero).
    pub cold_iterations: u32,
    /// Newton iterations from the candidate.
    pub warm_iterations: u32,
    /// The color of the claim (always `Estimated`).
    pub color: Color,
}

/// Measure warm-start savings on the toy nonlinear class.
///
/// # Errors
/// Returns [`Fem1dError`] when either run is malformed, unusable, or does not
/// converge within the admitted budget. Nonconvergence never becomes savings.
pub fn warm_start(
    problem: &MmsProblem,
    candidate: &[f64],
    max_iter: u32,
) -> Result<WarmStartReport, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let zero = try_zeroed("cold nonlinear start", problem.mesh().len())?;
    let cold = crate::fem1d::solve_nonlinear(problem, &zero, max_iter)?;
    require_converged(&cold, "cold nonlinear solve")?;
    let warm = crate::fem1d::solve_nonlinear(problem, candidate, max_iter)?;
    require_converged(&warm, "warm nonlinear solve")?;
    Ok(WarmStartReport {
        cold_iterations: cold.iterations,
        warm_iterations: warm.iterations,
        color: Color::Estimated {
            estimator: "warm-start-iteration-savings".to_string(),
            dispersion: f64::INFINITY,
        },
    })
}

/// Convenience for the batteries: effectivity of a report against the
/// oracle.
///
/// # Errors
/// Returns [`Fem1dError`] when the independent oracle or report bound is not a
/// usable finite value. Oracle failure is never mapped to effectivity `1.0`.
pub fn effectivity(
    problem: &MmsProblem,
    candidate: &[f64],
    report: &VerifierReport,
) -> Result<f64, Fem1dError> {
    if report.refusal.is_some() {
        return Err(Fem1dError::InvalidScalar {
            field: "verifier report",
            reason: "refused reports have no defined effectivity",
        });
    }
    let truth = true_energy_error(problem, candidate)?;
    if !report.bound.hi.is_finite() || report.bound.hi < 0.0 {
        return Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity report bound",
            index: None,
        });
    }
    if truth == 0.0 {
        return Err(Fem1dError::InvalidScalar {
            field: "oracle true error",
            reason: "effectivity is undefined for a zero denominator",
        });
    }
    let value = report.bound.hi / truth;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity ratio",
            index: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fem1d::{Poly, solve_p1};
    use fs_math::dd::Dd;

    fn receipt_fixture() -> (MmsProblem, Vec<f64>, f64, VerifierReceipt) {
        let problem = MmsProblem::new(
            "verifier-receipt-fixture",
            Poly::new(vec![0.0, 1.0, -1.0]).expect("admitted polynomial"),
            vec![0.0, 0.125, 0.375, 0.75, 1.0],
        )
        .expect("admitted manufactured problem");
        let candidate = solve_p1(&problem).expect("finite P1 fixture");
        let tolerance = 1.0;
        let receipt = verify_with_receipt(&problem, &candidate, tolerance)
            .expect("production verifier receipt");
        assert!(receipt.accepted());
        (problem, candidate, tolerance, receipt)
    }

    fn flipped(root: ContentHash) -> ContentHash {
        let mut bytes = *root.as_bytes();
        bytes[0] ^= 1;
        ContentHash(bytes)
    }

    #[test]
    fn producer_identity_reports_the_exact_compiled_feature_set() {
        let (_, _, _, receipt) = receipt_fixture();
        assert_eq!(
            receipt.producer().features(),
            current_verifier_feature_set()
        );
        assert_eq!(
            receipt.producer().features(),
            if cfg!(feature = "certified-speculation") {
                "certified-speculation"
            } else {
                "none"
            }
        );
    }

    #[test]
    fn production_receipt_retention_requires_independent_root_and_replay() {
        let (problem, candidate, tolerance, receipt) = receipt_fixture();
        let bytes = receipt.canonical_bytes().expect("canonical receipt bytes");
        let presented =
            PresentedVerifierReceipt::from_retained_bytes(&bytes, receipt.artifact_root())
                .expect("root-authenticated canonical receipt");
        assert_eq!(presented.artifact_root(), receipt.artifact_root());
        let presented_debug = format!("{presented:?}");
        assert_eq!(
            presented_debug,
            format!("PresentedVerifierReceipt({:?})", receipt.artifact_root()),
            "presented Debug is an artifact-identity allowlist"
        );
        for forbidden_field in [
            "receipt:",
            "theorem:",
            "qoi:",
            "units:",
            "hypotheses:",
            "bound_lo_bits:",
            "bound_hi_bits:",
            "tolerance_bits:",
            "accepted:",
        ] {
            assert!(
                !presented_debug.contains(forbidden_field),
                "presented Debug must redact scientific field {forbidden_field}: {presented_debug}"
            );
        }
        let admitted = admit_verifier_receipt(&problem, &candidate, tolerance, &presented)
            .expect("exact independent replay admits the receipt");
        assert_eq!(admitted.artifact_root(), receipt.artifact_root());
        assert_eq!(admitted.qoi(), receipt.qoi());
        assert_eq!(admitted.statement(), receipt.statement());
        assert_eq!(
            admitted.color(),
            Some(Color::Verified {
                lo: receipt.bound_lo(),
                hi: receipt.bound_hi(),
            }),
            "only the replay-admitted retained receipt exposes color authority"
        );

        let rejected_tolerance = receipt.bound_hi() / 2.0;
        assert!(rejected_tolerance.is_finite() && rejected_tolerance > 0.0);
        let rejected = verify_with_receipt(&problem, &candidate, rejected_tolerance)
            .expect("above-tolerance verification still publishes a receipt");
        assert!(!rejected.accepted());
        let rejected_bytes = rejected
            .canonical_bytes()
            .expect("canonical rejected receipt bytes");
        let rejected_presented = PresentedVerifierReceipt::from_retained_bytes(
            &rejected_bytes,
            rejected.artifact_root(),
        )
        .expect("root-authenticated rejected transport");
        let rejected_admitted = admit_verifier_receipt(
            &problem,
            &candidate,
            rejected_tolerance,
            &rejected_presented,
        )
        .expect("exact independent replay admits the rejected result data");
        assert!(!rejected_admitted.accepted());
        assert_eq!(
            rejected_admitted.color(),
            None,
            "replay admission cannot color an above-tolerance result"
        );

        let wrong_root = VerifierArtifactRoot(flipped(receipt.artifact_root().content_hash()));
        assert_eq!(
            PresentedVerifierReceipt::from_retained_bytes(&bytes, wrong_root),
            Err(VerifierReceiptError::ArtifactRootMismatch)
        );
        let mut corrupted = bytes.clone();
        corrupted[0] ^= 1;
        assert_eq!(
            PresentedVerifierReceipt::from_retained_bytes(&corrupted, receipt.artifact_root()),
            Err(VerifierReceiptError::ArtifactRootMismatch)
        );
        let corrupted_root =
            VerifierArtifactRoot(hash_domain(VERIFIER_RECEIPT_HASH_DOMAIN, &corrupted));
        assert_eq!(
            PresentedVerifierReceipt::from_retained_bytes(&corrupted, corrupted_root),
            Err(VerifierReceiptError::MalformedRetained("magic"))
        );

        let mut future = receipt.clone();
        future.schema_version += 1;
        future.artifact_root = future
            .calculated_artifact_root()
            .expect("future-version root");
        let future_bytes = future.canonical_bytes().expect("future-version bytes");
        assert_eq!(
            PresentedVerifierReceipt::from_retained_bytes(&future_bytes, future.artifact_root()),
            Err(VerifierReceiptError::MalformedRetained("schema version"))
        );

        // A foreign receipt can be internally self-consistent and root-valid.
        // It is still only presented bytes until the independent production
        // replay rejects its changed scientific semantics.
        let mut foreign = receipt.clone();
        foreign.qoi.push_str("-foreign");
        foreign.artifact_root = foreign
            .calculated_artifact_root()
            .expect("foreign root calculation");
        let foreign_bytes = foreign
            .canonical_bytes()
            .expect("foreign canonical receipt bytes");
        let foreign_presented =
            PresentedVerifierReceipt::from_retained_bytes(&foreign_bytes, foreign.artifact_root())
                .expect("self-consistent foreign bytes are presentable");
        assert!(matches!(
            admit_verifier_receipt(&problem, &candidate, tolerance, &foreign_presented),
            Err(VerifierReceiptError::ReplayMismatch)
        ));

        let mut trailing = bytes.clone();
        trailing.push(0);
        let trailing_root =
            VerifierArtifactRoot(hash_domain(VERIFIER_RECEIPT_HASH_DOMAIN, &trailing));
        assert_eq!(
            PresentedVerifierReceipt::from_retained_bytes(&trailing, trailing_root),
            Err(VerifierReceiptError::MalformedRetained("trailing bytes"))
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn every_receipt_semantic_field_moves_root_and_fails_exact_replay() {
        let (problem, candidate, tolerance, receipt) = receipt_fixture();
        let original_root = receipt.artifact_root();

        macro_rules! semantic_mutation {
            ($field:literal, $mutation:expr) => {{
                let mut changed = receipt.clone();
                $mutation(&mut changed);
                let moved_root = changed
                    .calculated_artifact_root()
                    .expect("mutated canonical root");
                assert_ne!(moved_root, original_root, "{} must move root", $field);
                changed.artifact_root = moved_root;
                let presented = PresentedVerifierReceipt { receipt: changed };
                assert!(
                    matches!(
                        admit_verifier_receipt(&problem, &candidate, tolerance, &presented),
                        Err(VerifierReceiptError::ReplayMismatch)
                    ),
                    "{} must fail independent replay",
                    $field
                );
            }};
        }

        semantic_mutation!("schema version", |r: &mut VerifierReceipt| r
            .schema_version +=
            1);
        semantic_mutation!("theorem", |r: &mut VerifierReceipt| r.theorem.push('x'));
        semantic_mutation!("producer crate", |r: &mut VerifierReceipt| r
            .producer
            .crate_name
            .push('x'));
        semantic_mutation!("producer version", |r: &mut VerifierReceipt| r
            .producer
            .crate_version
            .push('x'));
        semantic_mutation!("producer features", |r: &mut VerifierReceipt| r
            .producer
            .features
            .push('x'));
        semantic_mutation!("producer source", |r: &mut VerifierReceipt| r
            .producer
            .producer_source_root =
            flipped(r.producer.producer_source_root));
        semantic_mutation!("dependency source", |r: &mut VerifierReceipt| r
            .producer
            .dependency_source_root =
            flipped(r.producer.dependency_source_root));
        semantic_mutation!("workspace manifest", |r: &mut VerifierReceipt| r
            .producer
            .workspace_manifest_root =
            flipped(r.producer.workspace_manifest_root));
        semantic_mutation!("workspace lock", |r: &mut VerifierReceipt| r
            .producer
            .workspace_lock_root =
            flipped(r.producer.workspace_lock_root));
        semantic_mutation!("toolchain", |r: &mut VerifierReceipt| r
            .producer
            .toolchain_root =
            flipped(r.producer.toolchain_root));
        semantic_mutation!("problem identity version", |r: &mut VerifierReceipt| r
            .problem_identity_version +=
            1);
        semantic_mutation!("problem root", |r: &mut VerifierReceipt| r.problem_root =
            flipped(r.problem_root));
        semantic_mutation!("candidate root", |r: &mut VerifierReceipt| r
            .candidate_root =
            flipped(r.candidate_root));
        semantic_mutation!("mesh root", |r: &mut VerifierReceipt| r.mesh_root =
            flipped(r.mesh_root));
        semantic_mutation!("operator root", |r: &mut VerifierReceipt| r.operator_root =
            flipped(r.operator_root));
        semantic_mutation!("coefficient root", |r: &mut VerifierReceipt| r
            .coefficient_root =
            flipped(r.coefficient_root));
        semantic_mutation!("query root", |r: &mut VerifierReceipt| r.query_root =
            flipped(r.query_root));
        semantic_mutation!("QoI", |r: &mut VerifierReceipt| r.qoi.push('x'));
        semantic_mutation!("units", |r: &mut VerifierReceipt| r.units.push('x'));
        semantic_mutation!("flux hash", |r: &mut VerifierReceipt| r.flux_hash ^= 1);
        semantic_mutation!("verifier family", |r: &mut VerifierReceipt| r
            .verifier_family
            .push('x'));
        semantic_mutation!("arithmetic", |r: &mut VerifierReceipt| r
            .arithmetic
            .push('x'));
        semantic_mutation!("hypotheses", |r: &mut VerifierReceipt| {
            r.hypotheses.pop();
        });
        semantic_mutation!("lower endpoint", |r: &mut VerifierReceipt| r
            .bound_lo_bits ^=
            1);
        semantic_mutation!("upper endpoint", |r: &mut VerifierReceipt| r
            .bound_hi_bits ^=
            1);
        semantic_mutation!("tolerance", |r: &mut VerifierReceipt| r.tolerance_bits ^= 1);
        semantic_mutation!("acceptance", |r: &mut VerifierReceipt| r.accepted =
            !r.accepted);
        for index in 0..6 {
            semantic_mutation!("work-plan component", |r: &mut VerifierReceipt| r
                .work_plan[index] +=
                1);
        }
        semantic_mutation!("observed completed work", |r: &mut VerifierReceipt| r
            .observed_completed_work +=
            1);
        semantic_mutation!("observed planned work", |r: &mut VerifierReceipt| r
            .observed_planned_work +=
            1);
        semantic_mutation!("final phase", |r: &mut VerifierReceipt| r.final_phase =
            VerifierPhase::Hash);
        semantic_mutation!("final checkpoint", |r: &mut VerifierReceipt| r
            .final_checkpoint =
            VerifierCheckpointKind::PhaseEntry);
        semantic_mutation!("publication", |r: &mut VerifierReceipt| r
            .publication_observed =
            false);

        let mut bad_derived_root = receipt.clone();
        bad_derived_root.artifact_root =
            VerifierArtifactRoot(flipped(original_root.content_hash()));
        assert_eq!(
            bad_derived_root
                .calculated_artifact_root()
                .expect("derived root is excluded from canonical bytes"),
            original_root
        );
        let bad_derived_root = PresentedVerifierReceipt {
            receipt: bad_derived_root,
        };
        assert!(matches!(
            admit_verifier_receipt(&problem, &candidate, tolerance, &bad_derived_root),
            Err(VerifierReceiptError::ArtifactRootMismatch)
        ));
    }

    #[test]
    fn cancellation_is_real_checkpoint_telemetry_and_never_a_receipt() {
        let (problem, candidate, tolerance, receipt) = receipt_fixture();
        let cancelled =
            verify_with_receipt_cancellable(&problem, &candidate, tolerance, |progress| {
                progress.kind == VerifierCheckpointKind::PhaseEntry
                    && progress.phase == VerifierPhase::Equilibrated
            });
        let Err(VerifierReceiptError::Cancelled(telemetry)) = cancelled else {
            panic!("a real phase-entry cancellation must return telemetry only");
        };
        let progress = telemetry
            .last_progress()
            .expect("real callback was observed");
        assert_eq!(progress.kind, VerifierCheckpointKind::PhaseEntry);
        assert_eq!(progress.phase, VerifierPhase::Equilibrated);
        assert!(!telemetry.publication_observed());
        assert_eq!(telemetry.work_plan(), receipt.work_plan());

        // Even cancelling at the fully computed publication gate mints no
        // receipt: publication authority exists only after the callback wins.
        let at_publication =
            verify_with_receipt_cancellable(&problem, &candidate, tolerance, |progress| {
                progress.kind == VerifierCheckpointKind::Publication
            });
        let Err(VerifierReceiptError::Cancelled(telemetry)) = at_publication else {
            panic!("publication-gate cancellation must return telemetry only");
        };
        assert_eq!(
            telemetry.last_progress().map(|progress| progress.kind),
            Some(VerifierCheckpointKind::Publication)
        );
        assert!(!telemetry.publication_observed());
    }

    #[test]
    fn candidate_identity_moves_both_strong_root_and_flux_identity() {
        let (problem, candidate, tolerance, receipt) = receipt_fixture();
        let mut changed_candidate = candidate.clone();
        changed_candidate[1] = f64::from_bits(changed_candidate[1].to_bits() + 1);
        let changed = verify_with_receipt(&problem, &changed_candidate, tolerance)
            .expect("changed finite candidate still produces a receipt");
        assert_ne!(changed.candidate_root(), receipt.candidate_root());
        assert_ne!(changed.flux_hash(), receipt.flux_hash());
        assert_ne!(changed.artifact_root(), receipt.artifact_root());
    }

    fn contains_dd(interval: Iv, exact: Dd) -> bool {
        !exact.lt(Dd::from_f64(interval.lo)) && !Dd::from_f64(interval.hi).lt(exact)
    }

    #[test]
    fn intervalized_element_inputs_cover_independent_dd_oracle() {
        // These decimal-looking f64 inputs make every legacy point computation
        // below round away a nonzero residual. The double-double oracle therefore
        // detects removal of intervalized mesh, slope, node, or weight arithmetic.
        let (x0, x1) = (0.1, 0.4);
        let (candidate0, candidate1) = (0.1, 0.2);
        let (node_constant, weight_constant) = GAUSS5_REF[0];
        let (dx0, dx1) = (Dd::from_f64(x0), Dd::from_f64(x1));
        let half_constant = Dd::from_f64(0.5);

        let (h, midpoint, half) = interval_element_geometry(x0, x1).unwrap();
        let exact_h = dx1 - dx0;
        let exact_midpoint = (dx0 + dx1) * half_constant;
        let exact_half = exact_h * half_constant;
        assert_ne!(exact_h, Dd::from_f64(x1 - x0));
        assert_ne!(exact_midpoint, Dd::from_f64(f64::midpoint(x0, x1)));
        assert_ne!(exact_half, Dd::from_f64((x1 - x0) * 0.5));
        assert!(contains_dd(h, exact_h));
        assert!(contains_dd(midpoint, exact_midpoint));
        assert!(contains_dd(half, exact_half));

        let slope = interval_candidate_slope(candidate0, candidate1, h).unwrap();
        let exact_difference = Dd::from_f64(candidate1) - Dd::from_f64(candidate0);
        let rounded_slope = (candidate1 - candidate0) / (x1 - x0);
        assert_ne!(Dd::from_f64(rounded_slope) * exact_h, exact_difference);
        assert!(!(exact_difference).lt(Dd::from_f64(slope.lo) * exact_h));
        assert!(!(Dd::from_f64(slope.hi) * exact_h).lt(exact_difference));

        let (node, weight) =
            interval_quadrature_geometry(midpoint, half, node_constant, weight_constant).unwrap();
        let exact_node = exact_midpoint + exact_half * Dd::from_f64(node_constant);
        let exact_weight = exact_half * Dd::from_f64(weight_constant);
        let rounded_node = f64::midpoint(x0, x1) + (x1 - x0) * 0.5 * node_constant;
        let rounded_weight = (x1 - x0) * 0.5 * weight_constant;
        assert_ne!(exact_node, Dd::from_f64(rounded_node));
        assert_ne!(exact_weight, Dd::from_f64(rounded_weight));
        assert!(contains_dd(node, exact_node));
        assert!(contains_dd(weight, exact_weight));

        // `1/3` is not representable. The coefficient interval must reach the
        // side of the rounded quotient selected by the exact FMA residual;
        // treating the rounded antiderivative coefficient as a point fails it.
        let coefficient = interval_antiderivative_coefficient(1.0, 3).unwrap();
        let rounded = 1.0_f64 / 3.0;
        let residual = rounded.mul_add(3.0, -1.0);
        if residual > 0.0 {
            assert!(coefficient.lo <= crate::interval::down(rounded));
        } else if residual < 0.0 {
            assert!(coefficient.hi >= crate::interval::up(rounded));
        } else {
            assert!(coefficient.lo <= rounded && rounded <= coefficient.hi);
        }
    }

    #[test]
    fn gauss_constants_enclose_independent_truth_brackets() {
        // Each bit pair is the adjacent-f64 bracket around the corresponding
        // high-precision Gauss-Legendre constant, derived independently from
        // the decimal reference values. Fifteen-digit literals miss some
        // weights by up to eight ulps, so this locks the certified quadrature
        // inputs rather than merely checking that `iv_c` widens its input.
        let positive_constants = [
            (
                GAUSS5_REF[4].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[3].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
            (
                GAUSS5_REF[0].1,
                0x3fce_539e_c36e_038c,
                0x3fce_539e_c36e_038d,
            ),
            (
                GAUSS5_REF[1].1,
                0x3fde_a1da_25ae_415a,
                0x3fde_a1da_25ae_415b,
            ),
            (
                GAUSS5_REF[2].1,
                0x3fe2_3456_789a_bcdf,
                0x3fe2_3456_789a_bce0,
            ),
        ];
        for (constant, lower_bits, upper_bits) in positive_constants {
            let interval = iv_c(constant);
            assert!(interval.lo <= f64::from_bits(lower_bits));
            assert!(interval.hi >= f64::from_bits(upper_bits));
        }

        for (constant, positive_lower_bits, positive_upper_bits) in [
            (
                GAUSS5_REF[0].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[1].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
        ] {
            let interval = iv_c(constant);
            assert!(interval.lo <= -f64::from_bits(positive_upper_bits));
            assert!(interval.hi >= -f64::from_bits(positive_lower_bits));
        }
    }
}
