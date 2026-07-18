//! CMA-ES as natural-gradient IGO (plan §9.3, Bet 6): weighted
//! recombination with log-rank weights, rank-µ + rank-1 covariance
//! updates, cumulative step-size adaptation — the standard Hansen
//! couplings, which ARE the natural-gradient couplings on the Gaussian
//! statistical manifold. Rank-based selection makes the evolution
//! invariant to monotone transforms of the objective BY CONSTRUCTION
//! (property-tested bitwise, not cited).
//!
//! Determinism: sampling from a keyed Philox stream, `total_cmp` ranking
//! with lowest-index tie-breaks, fixed eigendecomposition cadence via
//! the landed cyclic Jacobi — the trajectory is a pure function of the
//! seed.

use fs_blake3::Blake3;
use fs_la::eigen::jacobi_eigh;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_rand::StreamKey;

/// Kernel id for CMA sampling streams (stable registry).
const K_CMA: u32 = 0xD1F0;

/// Domain stride used by the versioned checked BIPOP seed derivation.
/// The checked entry point refuses before any callback when the complete
/// conservative restart range could wrap this coordinate.
const BIPOP_RESTART_SEED_STRIDE: u64 = 0x9E37_79B9;

/// Largest zero-based population-doubling rung launched by BIPOP.
const BIPOP_LARGE_RUN_CAP: u32 = 8;

/// Local generation envelope assigned to one restart.
const BIPOP_GENERATIONS_PER_RESTART: usize = 250;

/// Sum of the large-population scales that can finance a later small restart.
/// Production terminates immediately after rung eight, so only rungs `0..=7`
/// contribute: `2^0 + ... + 2^7 = 2^8 - 1`.
const BIPOP_PRE_TERMINAL_LARGE_RUN_SCALE_SUM: usize = (1usize << BIPOP_LARGE_RUN_CAP as usize) - 1;

/// A full `250*lambda` local cap admits the initial callback and at most 249
/// complete populations: `1 + 249*lambda` callbacks in total.
const BIPOP_FULL_ENVELOPE_GENERATIONS: usize = BIPOP_GENERATIONS_PER_RESTART - 1;

/// One `(seed, kernel, tile)` Philox stream has exactly `2^64` distinct
/// counter positions. Consuming all of them once is valid; requesting one more
/// would reuse position zero because `Stream` advances its `u64` index modulo
/// this cardinality.
const FS_RAND_STREAM_COUNTER_CARDINALITY: u128 = 1u128 << 64;

/// Schema version for [`BipopAdmission`].
pub const BIPOP_ADMISSION_SCHEMA_VERSION: u32 = 3;

/// Schema version for [`BipopRestartRecord`].
pub const BIPOP_RESTART_SCHEMA_VERSION: u32 = 1;

/// Schema version for the root, callback trace, restart ledger, and full-study
/// identity retained by [`BipopReport`].
///
/// Version 3 makes the earliest IEEE-754 total-order minimum authoritative
/// inside each CMA restart, separates that representative from numeric target
/// witnesses, and composes the length-framed trace-identity v2.
pub const BIPOP_REPORT_SCHEMA_VERSION: u32 = 3;

/// Schema version for each borrowed [`BipopEvaluationRecord`].
pub const BIPOP_EVALUATION_SCHEMA_VERSION: u32 = 1;

/// Canonical fs-obs identity kind for one exact BIPOP root-input preimage.
pub const BIPOP_ROOT_IDENTITY_KIND: &str = "fs-dfo-bipop-root-v1";

/// Schema version for the canonical BIPOP root-input identity.
pub const BIPOP_ROOT_IDENTITY_SCHEMA_VERSION: u32 = 1;

/// Schema version for the allocation-free streaming trace identity.
///
/// Version 2 length-frames the domain before the fixed-width trace header, so
/// cross-domain separation follows from the byte grammar rather than from a
/// registry convention about domain suffixes.
pub const BIPOP_TRACE_IDENTITY_SCHEMA_VERSION: u32 = 2;

/// Domain prefix for the production BIPOP callback-trace BLAKE3.
pub const BIPOP_TRACE_IDENTITY_DOMAIN: &str = "frankensim.fs-dfo.bipop-callback-trace.v2";

/// Schema version for the complete production BIPOP study identity.
pub const BIPOP_STUDY_IDENTITY_SCHEMA_VERSION: u32 = 1;

/// Domain prefix for the complete root, callback, and restart-ledger identity.
pub const BIPOP_STUDY_IDENTITY_DOMAIN: &str = "frankensim.fs-dfo.bipop-full-study.v1";

/// Owner-local declaration for the canonical BIPOP root-input identity.
pub const BIPOP_ROOT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-dfo:bipop-root-inputs",
    "version_const=BIPOP_ROOT_IDENTITY_SCHEMA_VERSION",
    "version=1",
    "domain=fs-dfo-bipop-root-v1",
    "domain_const=BIPOP_ROOT_IDENTITY_KIND",
    "encoder=bipop_root_identity",
    "encoder_helpers=none",
    "schema_constants=BIPOP_ROOT_IDENTITY_SCHEMA_VERSION,BIPOP_ROOT_IDENTITY_KIND,BIPOP_REPORT_SCHEMA_VERSION,BIPOP_ADMISSION_SCHEMA_VERSION",
    "schema_functions=build_bipop_root_inputs",
    "schema_dependencies=fs-obs:replay-identity-frame",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=BipopRootIdentitySource",
    "source_fields=BipopRootIdentitySource.start:semantic,BipopRootIdentitySource.sigma:semantic,BipopRootIdentitySource.total_budget:semantic,BipopRootIdentitySource.target:semantic,BipopRootIdentitySource.seed:semantic",
    "source_bindings=BipopRootIdentitySource.start>dimension+start-coordinate-bits,BipopRootIdentitySource.sigma>initial-sigma-bits,BipopRootIdentitySource.total_budget>total-budget,BipopRootIdentitySource.target>target-presence+target-bits,BipopRootIdentitySource.seed>root-seed",
    "external_semantic_fields=root-domain,root-schema-version,canonical-field-order,report-schema-version,admission-schema-version",
    "semantic_fields=root-domain,root-schema-version,canonical-field-order,report-schema-version,admission-schema-version,dimension,start-coordinate-bits,initial-sigma-bits,total-budget,target-presence,target-bits,root-seed",
    "excluded_fields=none",
    "consumers=build_bipop_root_inputs,BipopReport::root_inputs,BipopReport::validate_ledger,BipopReport::admit_study_identity_with_replay",
    "mutations=root-domain:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,root-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,canonical-field-order:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,report-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,admission-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,dimension:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,start-coordinate-bits:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,initial-sigma-bits:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,total-budget:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,target-presence:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,target-bits:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently,root-seed:crates/fs-dfo/src/cma.rs#bipop_root_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_bipop_root_identity_fields",
    "transport_guard=BipopReport::validate_ledger",
    "version_guard=crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-dfo:bipop-root-inputs",
];

/// Owner-local declaration for the root-bound production callback trace.
pub const BIPOP_TRACE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-dfo:bipop-callback-trace",
    "version_const=BIPOP_TRACE_IDENTITY_SCHEMA_VERSION",
    "version=2",
    "domain=frankensim.fs-dfo.bipop-callback-trace.v2",
    "domain_const=BIPOP_TRACE_IDENTITY_DOMAIN",
    "encoder=bipop_trace_identity",
    "encoder_helpers=none",
    "schema_constants=BIPOP_TRACE_IDENTITY_SCHEMA_VERSION,BIPOP_TRACE_IDENTITY_DOMAIN,BIPOP_EVALUATION_SCHEMA_VERSION",
    "schema_functions=build_bipop_trace_identity,crates/fs-blake3/src/lib.rs#Blake3::new,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes",
    "schema_dependencies=fs-dfo:bipop-root-inputs",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=BipopTraceIdentitySource",
    "source_fields=BipopTraceIdentitySource.root_canonical_bytes:semantic,BipopTraceIdentitySource.dimension:semantic,BipopTraceIdentitySource.rows:semantic,BipopTraceIdentitySource.points:semantic",
    "source_bindings=BipopTraceIdentitySource.root_canonical_bytes>root-canonical-bytes,BipopTraceIdentitySource.dimension>dimension,BipopTraceIdentitySource.rows>row-count+row-schema-version+restart-ordinal+local-offset+objective-bits+row-order,BipopTraceIdentitySource.points>point-coordinate-bits+point-order",
    "external_semantic_fields=trace-domain,trace-schema-version,canonical-field-order",
    "semantic_fields=trace-domain,trace-schema-version,canonical-field-order,root-canonical-bytes,dimension,row-count,row-schema-version,restart-ordinal,local-offset,objective-bits,row-order,point-coordinate-bits,point-order",
    "excluded_fields=none",
    "consumers=build_bipop_trace_identity,BipopReport::trace_identity,BipopReport::validate_ledger,BipopReport::validate_study_identity",
    "mutations=trace-domain:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,trace-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,canonical-field-order:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,root-canonical-bytes:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,dimension:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,row-count:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,row-schema-version:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,restart-ordinal:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,local-offset:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,objective-bits:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,row-order:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,point-coordinate-bits:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently,point-order:crates/fs-dfo/src/cma.rs#bipop_trace_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_bipop_trace_identity_fields",
    "transport_guard=BipopReport::validate_ledger",
    "version_guard=crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-dfo:bipop-callback-trace",
];

/// Owner-local declaration for the complete BIPOP study identity.
pub const BIPOP_STUDY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-dfo:bipop-full-study",
    "version_const=BIPOP_STUDY_IDENTITY_SCHEMA_VERSION",
    "version=1",
    "domain=frankensim.fs-dfo.bipop-full-study.v1",
    "domain_const=BIPOP_STUDY_IDENTITY_DOMAIN",
    "encoder=bipop_study_identity",
    "encoder_helpers=study_hash_field,study_hash_u64,study_hash_usize,study_hash_u32,study_hash_f64,study_hash_flag,bipop_lane_tag,cma_stop_reason_tag",
    "schema_constants=BIPOP_STUDY_IDENTITY_SCHEMA_VERSION,BIPOP_STUDY_IDENTITY_DOMAIN,BIPOP_REPORT_SCHEMA_VERSION,BIPOP_ADMISSION_SCHEMA_VERSION,BIPOP_RESTART_SCHEMA_VERSION,BIPOP_EVALUATION_SCHEMA_VERSION,BIPOP_TRACE_IDENTITY_SCHEMA_VERSION,crates/fs-rand/src/lib.rs#STREAM_SEMANTICS_VERSION,crates/fs-la/src/eigen.rs#JACOBI_EIGH_ADMISSION_SCHEMA_VERSION,K_CMA,BIPOP_RESTART_SEED_STRIDE,BIPOP_LARGE_RUN_CAP,BIPOP_GENERATIONS_PER_RESTART",
    "schema_functions=build_bipop_study_identity,crates/fs-blake3/src/lib.rs#Blake3::new,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes",
    "schema_dependencies=fs-dfo:bipop-root-inputs,fs-dfo:bipop-callback-trace",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=BipopStudyIdentitySource",
    "source_fields=BipopStudyIdentitySource.report_schema_version:semantic,BipopStudyIdentitySource.root:semantic,BipopStudyIdentitySource.root_content_identity:semantic,BipopStudyIdentitySource.schedule:semantic,BipopStudyIdentitySource.total_evals:semantic,BipopStudyIdentitySource.records:semantic,BipopStudyIdentitySource.best_restart:semantic,BipopStudyIdentitySource.best:semantic,BipopStudyIdentitySource.retained_trace_identity:semantic,BipopStudyIdentitySource.callback_content_identity:semantic",
    "source_bindings=BipopStudyIdentitySource.report_schema_version>report-schema-version,BipopStudyIdentitySource.root>retained-root-canonical-bytes,BipopStudyIdentitySource.root_content_identity>recomputed-root-canonical-bytes,BipopStudyIdentitySource.schedule>schedule,BipopStudyIdentitySource.total_evals>total-evals,BipopStudyIdentitySource.records>restart-records,BipopStudyIdentitySource.best_restart>best-restart,BipopStudyIdentitySource.best>best-projection,BipopStudyIdentitySource.retained_trace_identity>retained-trace-identity,BipopStudyIdentitySource.callback_content_identity>callback-content-identity",
    "external_semantic_fields=study-domain,study-schema-version,canonical-field-order,admission-schema-version,restart-schema-version,evaluation-schema-version,trace-schema-version,rand-stream-semantics-version,jacobi-admission-schema-version,cma-stream-kernel,restart-seed-stride,large-run-cap,generations-per-restart",
    "semantic_fields=study-domain,study-schema-version,canonical-field-order,report-schema-version,admission-schema-version,restart-schema-version,evaluation-schema-version,trace-schema-version,rand-stream-semantics-version,jacobi-admission-schema-version,cma-stream-kernel,restart-seed-stride,large-run-cap,generations-per-restart,retained-root-canonical-bytes,recomputed-root-canonical-bytes,schedule,total-evals,restart-records,best-restart,best-projection,retained-trace-identity,callback-content-identity",
    "excluded_fields=none",
    "consumers=build_bipop_study_identity,BipopReport::study_identity,BipopReport::validate_study_identity,BipopReport::admit_study_identity,BipopReport::admit_study_identity_with_replay",
    "mutations=study-domain:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,study-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,canonical-field-order:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,report-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,admission-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,restart-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,evaluation-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,trace-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,rand-stream-semantics-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,jacobi-admission-schema-version:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,cma-stream-kernel:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,restart-seed-stride:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,large-run-cap:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,generations-per-restart:crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed,retained-root-canonical-bytes:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,recomputed-root-canonical-bytes:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,schedule:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,total-evals:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,restart-records:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,best-restart:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,best-projection:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,retained-trace-identity:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently,callback-content-identity:crates/fs-dfo/src/cma.rs#bipop_study_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_bipop_study_identity_fields",
    "transport_guard=BipopReport::admit_study_identity",
    "version_guard=crates/fs-dfo/src/cma.rs#bipop_identity_versions_and_domains_fail_closed",
    "coupling_surface=fs-dfo:bipop-full-study",
];

/// Sealed result of callback-free BIPOP input and arithmetic admission.
///
/// The restart bound combines two independent proofs. One objective evaluation
/// is the minimum spend of a launched restart, so `total_budget - 1` is a
/// budget bound. Independently, production admits nine large rungs (`0..=8`)
/// and terminates immediately after publishing rung eight. Only rungs `0..=7`
/// can therefore finance later small restarts. A full local envelope on rung
/// `r` spends at most `1 + 249*base_lambda*2^r` callbacks. Before the `j`th
/// small restart can launch, its predecessor small spend still trails the
/// cumulative large spend. While more than `base_lambda` aggregate callbacks
/// remain, every continuing small restart completes at least one population
/// and therefore spends at least `base_lambda + 1` callbacks. At most
/// `ceil((8 + 249*base_lambda*(2^8-1))/(base_lambda+1))` such restarts can
/// launch. Once no more than `base_lambda` callbacks remain, at most that many
/// further launches are possible: every continuing launch spends at least one
/// callback, while a zero-callback generated-start refusal terminates the run.
/// Adding the eight preceding large ordinals gives the scheduler bound below.
/// Admission uses the minimum of this theorem and the budget bound.
///
/// The receipt also proves that neither the shared restart-perturbation stream
/// nor any per-restart CMA stream can reuse a Philox counter coordinate. Holding
/// it proves only that these formulas were representable; it is not an
/// authenticated identity for the start, target, sigma, seed, or callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BipopAdmission {
    schema_version: u32,
    stream_semantics_version: u32,
    jacobi_admission: Option<fs_la::eigen::JacobiEighAdmission>,
    dimension: usize,
    total_budget: usize,
    base_lambda: usize,
    max_large_lambda: usize,
    max_local_budget: usize,
    max_restart_ordinal: u64,
    max_matrix_entries: usize,
    max_population_entries: usize,
    max_restart_stream_blocks: u128,
    max_cma_stream_blocks: u128,
}

impl BipopAdmission {
    /// Admission schema used for the receipt.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// fs-rand stream semantics whose exact two-block normal consumption is
    /// assumed by both counter-range proofs.
    #[must_use]
    pub fn stream_semantics_version(&self) -> u32 {
        self.stream_semantics_version
    }

    /// Exact fs-la dense-Jacobi capability receipt when the aggregate budget
    /// can reach at least one complete CMA generation. Initial-evaluation-only
    /// schedules return `None` because they never invoke that dependency.
    #[must_use]
    pub fn jacobi_admission(&self) -> Option<fs_la::eigen::JacobiEighAdmission> {
        self.jacobi_admission
    }

    /// Decision-vector dimension admitted for every restart.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Hard aggregate callback budget.
    #[must_use]
    pub fn total_budget(&self) -> usize {
        self.total_budget
    }

    /// Standard population at the first large and every small restart.
    #[must_use]
    pub fn base_lambda(&self) -> usize {
        self.base_lambda
    }

    /// Conservative population representability bound for the admitted ladder.
    #[must_use]
    pub fn max_large_lambda(&self) -> usize {
        self.max_large_lambda
    }

    /// Largest pre-minimum local budget formula (`lambda * 250`).
    #[must_use]
    pub fn max_local_budget(&self) -> usize {
        self.max_local_budget
    }

    /// Conservative largest restart ordinal under both budget and ladder caps.
    #[must_use]
    pub fn max_restart_ordinal(&self) -> u64 {
        self.max_restart_ordinal
    }

    /// Largest dense square-matrix element count needed by one CMA run.
    #[must_use]
    pub fn max_matrix_entries(&self) -> usize {
        self.max_matrix_entries
    }

    /// Largest population-coordinate element count needed by one CMA run.
    #[must_use]
    pub fn max_population_entries(&self) -> usize {
        self.max_population_entries
    }

    /// Most Philox counter blocks consumed by the shared restart stream.
    #[must_use]
    pub fn max_restart_stream_blocks(&self) -> u128 {
        self.max_restart_stream_blocks
    }

    /// Most Philox counter blocks consumed by any one CMA restart stream.
    #[must_use]
    pub fn max_cma_stream_blocks(&self) -> u128 {
        self.max_cma_stream_blocks
    }
}

/// Structured refusal from [`admit_bipop`] or [`try_bipop_cmaes`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BipopError {
    /// A zero-dimensional decision cannot enter CMA.
    EmptyStart,
    /// One supplied coordinate is NaN or infinite.
    NonFiniteStart {
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// Sigma must be finite and strictly positive (including refusing ±0).
    InvalidSigma {
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// `Some(target)` must be finite; `None` is the typed no-target spelling.
    NonFiniteTarget {
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// No restart can be launched under a zero callback budget.
    ZeroBudget,
    /// A dimension-derived storage formula is not representable.
    DimensionOverflow {
        /// Stable formula label.
        what: &'static str,
    },
    /// A native field cannot be represented by the canonical identity
    /// schema's unsigned 64-bit integer encoding.
    IdentityFieldOverflow {
        /// Stable field label.
        what: &'static str,
    },
    /// Capacity for the complete retained callback trace could not be
    /// reserved before the first objective invocation.
    TraceAllocationFailed {
        /// Evaluation rows requested by the hard budget.
        evaluations: usize,
        /// Decision-coordinate entries requested by dimension times budget.
        point_entries: usize,
    },
    /// The budget/ladder restart envelope cannot be represented as a u64 ordinal.
    RestartOrdinalOverflow {
        /// Requested aggregate budget.
        total_budget: usize,
    },
    /// The conservative derived-seed range would wrap u64.
    SeedRangeOverflow {
        /// Root seed.
        seed: u64,
        /// Largest conservatively reachable restart ordinal.
        max_restart_ordinal: u64,
    },
    /// The large-population ladder cannot be represented.
    PopulationOverflow {
        /// Zero-based large-run rung.
        large_run: u32,
    },
    /// `lambda * 250` cannot be represented.
    LocalBudgetOverflow {
        /// Population whose local budget overflowed.
        lambda: usize,
    },
    /// A Philox stream would need more than its `2^64` distinct positions.
    RandomCounterRangeOverflow {
        /// Stable stream-domain label.
        stream: &'static str,
        /// Exact requested block count, or `None` if even the `u128` product
        /// was not representable.
        required_blocks: Option<u128>,
    },
    /// The shared fs-la dense-Jacobi authority refused this decision dimension.
    EigensolverAdmission {
        /// Exact dependency refusal; no duplicated capability formula is used.
        error: fs_la::eigen::JacobiEighAdmissionError,
    },
    /// Checked scheduler accounting failed at a named restart boundary.
    ArithmeticOverflow {
        /// Restart ordinal being prepared or finalized.
        restart: u64,
        /// Stable formula label.
        what: &'static str,
    },
    /// A finite admitted root plus a finite perturbation produced a non-finite
    /// restart coordinate; the affected restart was not invoked.
    GeneratedStartNonFinite {
        /// Restart ordinal.
        restart: u64,
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// CMA generated a non-finite decision coordinate. Earlier callbacks or
    /// restarts may already have completed, but the affected candidate was not
    /// passed to the objective.
    GeneratedCandidateNonFinite {
        /// Restart ordinal.
        restart: u64,
        /// One-based CMA generation.
        generation: usize,
        /// Zero-based candidate within the population.
        candidate: usize,
        /// Coordinate index.
        component: usize,
        /// Exact IEEE-754 payload.
        bits: u64,
    },
    /// A CMA run violated the local hard budget assumed by the scheduler.
    InternalBudgetViolation {
        /// Restart ordinal.
        restart: u64,
        /// Callbacks reported by CMA.
        spent: usize,
        /// Local cap supplied to CMA.
        allocated: usize,
    },
    /// Aggregate trace accounting exceeded the admitted hard budget.
    InternalAggregateBudgetViolation {
        /// Restart ordinal whose report crossed the boundary.
        restart: u64,
        /// Aggregate trace end after the restart.
        total_spent: usize,
        /// Admitted aggregate cap.
        total_budget: usize,
    },
    /// An admitted scheduler reached a state ruled out by its preflight and
    /// loop invariants.
    InternalInvariant {
        /// Stable invariant label.
        what: &'static str,
    },
    /// Generated evidence failed its own structural validator.
    GeneratedLedgerInvalid {
        /// Restart index when the invariant is local.
        restart: Option<usize>,
        /// Stable validator invariant.
        invariant: &'static str,
    },
}

impl core::fmt::Display for BipopError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match *self {
            Self::EmptyStart => {
                formatter.write_str("BIPOP start must contain at least one coordinate")
            }
            Self::NonFiniteStart { component, bits } => write!(
                formatter,
                "BIPOP start component {component} is non-finite (bits 0x{bits:016x})"
            ),
            Self::InvalidSigma { bits } => write!(
                formatter,
                "BIPOP sigma must be finite and strictly positive (bits 0x{bits:016x})"
            ),
            Self::NonFiniteTarget { bits } => write!(
                formatter,
                "BIPOP target must be finite when present (bits 0x{bits:016x}); use None for no target"
            ),
            Self::ZeroBudget => formatter.write_str("BIPOP total callback budget must be positive"),
            Self::DimensionOverflow { what } => {
                write!(formatter, "BIPOP dimension formula `{what}` overflowed")
            }
            Self::IdentityFieldOverflow { what } => write!(
                formatter,
                "BIPOP identity field `{what}` does not fit its u64 encoding"
            ),
            Self::TraceAllocationFailed {
                evaluations,
                point_entries,
            } => write!(
                formatter,
                "BIPOP could not reserve retained trace capacity for {evaluations} evaluations and {point_entries} decision coordinates"
            ),
            Self::RestartOrdinalOverflow { total_budget } => write!(
                formatter,
                "BIPOP restart envelope for budget {total_budget} cannot be represented by the u64 ordinal"
            ),
            Self::SeedRangeOverflow {
                seed,
                max_restart_ordinal,
            } => write!(
                formatter,
                "BIPOP seed range from root {seed} through restart {max_restart_ordinal} would wrap"
            ),
            Self::PopulationOverflow { large_run } => write!(
                formatter,
                "BIPOP population ladder overflowed at large run {large_run}"
            ),
            Self::LocalBudgetOverflow { lambda } => write!(
                formatter,
                "BIPOP local budget formula overflowed for population {lambda}"
            ),
            Self::RandomCounterRangeOverflow {
                stream,
                required_blocks,
            } => match required_blocks {
                Some(required_blocks) => write!(
                    formatter,
                    "BIPOP {stream} stream requires {required_blocks} Philox blocks, exceeding the 2^64-position counter domain"
                ),
                None => write!(
                    formatter,
                    "BIPOP {stream} stream block count overflowed the u128 admission accumulator"
                ),
            },
            Self::EigensolverAdmission { error } => {
                write!(formatter, "BIPOP eigensolver admission refused: {error}")
            }
            Self::ArithmeticOverflow { restart, what } => write!(
                formatter,
                "BIPOP scheduler formula `{what}` overflowed at restart {restart}"
            ),
            Self::GeneratedStartNonFinite {
                restart,
                component,
                bits,
            } => write!(
                formatter,
                "BIPOP restart {restart} start component {component} became non-finite (bits 0x{bits:016x})"
            ),
            Self::GeneratedCandidateNonFinite {
                restart,
                generation,
                candidate,
                component,
                bits,
            } => write!(
                formatter,
                "BIPOP restart {restart} generation {generation} candidate {candidate} component {component} became non-finite (bits 0x{bits:016x})"
            ),
            Self::InternalBudgetViolation {
                restart,
                spent,
                allocated,
            } => write!(
                formatter,
                "BIPOP restart {restart} spent {spent} callbacks under local cap {allocated}"
            ),
            Self::InternalAggregateBudgetViolation {
                restart,
                total_spent,
                total_budget,
            } => write!(
                formatter,
                "BIPOP restart {restart} advanced aggregate callbacks to {total_spent} under hard cap {total_budget}"
            ),
            Self::InternalInvariant { what } => {
                write!(formatter, "BIPOP internal invariant failed: {what}")
            }
            Self::GeneratedLedgerInvalid { restart, invariant } => match restart {
                Some(restart) => write!(
                    formatter,
                    "generated BIPOP restart {restart} violates {invariant}"
                ),
                None => write!(formatter, "generated BIPOP ledger violates {invariant}"),
            },
        }
    }
}

impl std::error::Error for BipopError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::EigensolverAdmission { error } => Some(error),
            _ => None,
        }
    }
}

fn checked_random_counter_blocks(
    stream: &'static str,
    factors: &[u128],
) -> Result<u128, BipopError> {
    let required_blocks = factors
        .iter()
        .try_fold(1u128, |product, factor| product.checked_mul(*factor));
    let Some(required_blocks) = required_blocks else {
        return Err(BipopError::RandomCounterRangeOverflow {
            stream,
            required_blocks: None,
        });
    };
    if required_blocks > FS_RAND_STREAM_COUNTER_CARDINALITY {
        return Err(BipopError::RandomCounterRangeOverflow {
            stream,
            required_blocks: Some(required_blocks),
        });
    }
    Ok(required_blocks)
}

fn restart_stream_block_bound(
    dimension: usize,
    max_restart_ordinal: u64,
) -> Result<u128, BipopError> {
    checked_random_counter_blocks(
        "restart-perturbation",
        &[2, dimension as u128, u128::from(max_restart_ordinal)],
    )
}

fn cma_stream_block_bound(
    dimension: usize,
    lambda: usize,
    max_generations: usize,
) -> Result<u128, BipopError> {
    checked_random_counter_blocks(
        "CMA-candidate",
        &[
            2,
            dimension as u128,
            lambda as u128,
            max_generations as u128,
        ],
    )
}

fn checked_dense_matrix_allocation(max_matrix_entries: usize) -> Result<usize, BipopError> {
    let bytes = max_matrix_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "dense covariance bytes",
        })?;
    if bytes > isize::MAX as usize {
        return Err(BipopError::DimensionOverflow {
            what: "dense covariance address space",
        });
    }
    Ok(bytes)
}

fn scheduler_max_restart_ordinal(
    base_lambda: usize,
    total_budget: usize,
) -> Result<u64, BipopError> {
    let preterminal_large_spend = (base_lambda as u128)
        .checked_mul(BIPOP_FULL_ENVELOPE_GENERATIONS as u128)
        .and_then(|spend| spend.checked_mul(BIPOP_PRE_TERMINAL_LARGE_RUN_SCALE_SUM as u128))
        .and_then(|spend| spend.checked_add(u128::from(BIPOP_LARGE_RUN_CAP)))
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let continuing_small_spend = (base_lambda as u128)
        .checked_add(1)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let max_full_generation_small_count = preterminal_large_spend
        .checked_add(continuing_small_spend - 1)
        .map(|numerator| numerator / continuing_small_spend)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let max_small_restart_count = max_full_generation_small_count
        .checked_add(base_lambda as u128)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let scheduler_max_restart_ordinal = u128::from(BIPOP_LARGE_RUN_CAP)
        .checked_add(max_small_restart_count)
        .ok_or(BipopError::RestartOrdinalOverflow { total_budget })?;
    let budget_max_restart_ordinal = (total_budget - 1) as u128;
    u64::try_from(budget_max_restart_ordinal.min(scheduler_max_restart_ordinal))
        .map_err(|_| BipopError::RestartOrdinalOverflow { total_budget })
}

/// Validate BIPOP inputs and every conservative arithmetic/storage envelope
/// before allocating scheduler state or invoking the objective.
///
/// `target = None` disables target stopping. This avoids using a non-finite
/// floating-point sentinel on the checked API.
///
/// # Errors
/// Returns [`BipopError`] for malformed input or an unrepresentable envelope.
pub fn admit_bipop(
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
) -> Result<BipopAdmission, BipopError> {
    if x0.is_empty() {
        return Err(BipopError::EmptyStart);
    }
    for (component, value) in x0.iter().enumerate() {
        if !value.is_finite() {
            return Err(BipopError::NonFiniteStart {
                component,
                bits: value.to_bits(),
            });
        }
    }
    if !sigma0.is_finite() || sigma0 <= 0.0 {
        return Err(BipopError::InvalidSigma {
            bits: sigma0.to_bits(),
        });
    }
    if let Some(target) = target {
        if !target.is_finite() {
            return Err(BipopError::NonFiniteTarget {
                bits: target.to_bits(),
            });
        }
    }
    if total_budget == 0 {
        return Err(BipopError::ZeroBudget);
    }

    let dimension = x0.len();
    let lambda_offset = (3.0 * fs_math::det::ln(dimension as f64)).floor();
    if !lambda_offset.is_finite() || lambda_offset < 0.0 || lambda_offset > usize::MAX as f64 {
        return Err(BipopError::DimensionOverflow {
            what: "base population",
        });
    }
    let base_lambda =
        4usize
            .checked_add(lambda_offset as usize)
            .ok_or(BipopError::DimensionOverflow {
                what: "base population",
            })?;
    let jacobi_admission = if total_budget > base_lambda {
        Some(
            fs_la::eigen::admit_jacobi_eigh(dimension)
                .map_err(|error| BipopError::EigensolverAdmission { error })?,
        )
    } else {
        None
    };

    let max_restart_ordinal = scheduler_max_restart_ordinal(base_lambda, total_budget)?;
    let seed_delta = max_restart_ordinal
        .checked_mul(BIPOP_RESTART_SEED_STRIDE)
        .ok_or(BipopError::SeedRangeOverflow {
            seed,
            max_restart_ordinal,
        })?;
    seed.checked_add(seed_delta)
        .ok_or(BipopError::SeedRangeOverflow {
            seed,
            max_restart_ordinal,
        })?;

    let budget_ladder_cap = u32::try_from(max_restart_ordinal)
        .unwrap_or(u32::MAX)
        .min(BIPOP_LARGE_RUN_CAP);
    let scale = 1usize
        .checked_shl(budget_ladder_cap)
        .ok_or(BipopError::PopulationOverflow {
            large_run: budget_ladder_cap,
        })?;
    let max_large_lambda =
        base_lambda
            .checked_mul(scale)
            .ok_or(BipopError::PopulationOverflow {
                large_run: budget_ladder_cap,
            })?;
    let max_local_budget = max_large_lambda
        .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
        .ok_or(BipopError::LocalBudgetOverflow {
            lambda: max_large_lambda,
        })?;
    let max_matrix_entries =
        dimension
            .checked_mul(dimension)
            .ok_or(BipopError::DimensionOverflow {
                what: "dense covariance entries",
            })?;
    let max_population_entries =
        dimension
            .checked_mul(max_large_lambda)
            .ok_or(BipopError::DimensionOverflow {
                what: "population coordinate entries",
            })?;
    checked_dense_matrix_allocation(max_matrix_entries)?;
    max_population_entries
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "population coordinate bytes",
        })?;

    let max_restart_stream_blocks = restart_stream_block_bound(dimension, max_restart_ordinal)?;
    let mut max_cma_stream_blocks = 0u128;
    for large_run in 0..=budget_ladder_cap {
        let scale = 1usize
            .checked_shl(large_run)
            .ok_or(BipopError::PopulationOverflow { large_run })?;
        let lambda = base_lambda
            .checked_mul(scale)
            .ok_or(BipopError::PopulationOverflow { large_run })?;
        let local_envelope = lambda
            .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
            .ok_or(BipopError::LocalBudgetOverflow { lambda })?;
        let allocated_budget = local_envelope.min(total_budget);
        // CMA spends one callback at the start, then admits only complete
        // populations: `1 + generations * lambda <= allocated_budget`.
        let max_generations = (allocated_budget - 1) / lambda;
        let blocks = cma_stream_block_bound(dimension, lambda, max_generations)?;
        max_cma_stream_blocks = max_cma_stream_blocks.max(blocks);
    }

    Ok(BipopAdmission {
        schema_version: BIPOP_ADMISSION_SCHEMA_VERSION,
        stream_semantics_version: fs_rand::STREAM_SEMANTICS_VERSION,
        jacobi_admission,
        dimension,
        total_budget,
        base_lambda,
        max_large_lambda,
        max_local_budget,
        max_restart_ordinal,
        max_matrix_entries,
        max_population_entries,
        max_restart_stream_blocks,
        max_cma_stream_blocks,
    })
}

/// Tunables (defaults follow Hansen's standard settings).
#[derive(Debug, Clone)]
pub struct CmaParams {
    /// Population size λ (default 4 + ⌊3·ln n⌋).
    pub lambda: usize,
    /// Initial step size σ₀.
    pub sigma0: f64,
    /// Evaluation budget.
    pub max_evals: usize,
    /// Target objective value (stop when reached).
    pub f_target: f64,
    /// Generations between eigendecompositions (SPD refresh cadence).
    pub eigen_interval: usize,
}

impl CmaParams {
    /// Standard defaults for dimension `n`.
    #[must_use]
    pub fn standard(n: usize, sigma0: f64, max_evals: usize, f_target: f64) -> CmaParams {
        let lambda = 4 + (3.0 * fs_math::det::ln(n as f64)).floor() as usize;
        CmaParams {
            lambda,
            sigma0,
            max_evals,
            f_target,
            eigen_interval: 1,
        }
    }
}

/// Why one CMA run stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmaStopReason {
    /// The requested objective target was reached.
    TargetReached,
    /// The local budget could not admit another complete population.
    BudgetExhausted,
    /// TolX/TolFun stopped a run.
    Stagnated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CmaGeneratedPointError {
    generation: usize,
    candidate: usize,
    component: usize,
    bits: u64,
}

/// Run evidence.
#[derive(Debug, Clone)]
pub struct CmaReport {
    /// Best point found.
    pub x_best: Vec<f64>,
    /// Earliest best objective under exact [`f64::total_cmp`] ordering.
    pub f_best: f64,
    /// Objective evaluations consumed.
    pub evals: usize,
    /// Generations run.
    pub generations: usize,
    /// Whether any evaluated objective numerically reached the finite target.
    ///
    /// This is deliberately independent of [`Self::f_best`]: IEEE-754 total
    /// order can select a negative NaN ahead of a finite target witness.
    pub converged: bool,
    /// Final step size (diagnostic).
    pub sigma: f64,
}

/// Full-covariance CMA-ES from `x0`. Deterministic per `seed`.
///
/// # Panics
/// Panics before the affected callback if the initial or a generated decision
/// point contains a non-finite coordinate. The fallible BIPOP surface projects
/// the same guard into [`BipopError::GeneratedCandidateNonFinite`].
#[must_use]
pub fn cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
) -> CmaReport {
    cmaes_with_stop(f, x0, params, seed).0
}

#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
fn cmaes_with_stop<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
) -> (CmaReport, CmaStopReason) {
    cmaes_with_stop_target(f, x0, params, seed, Some(params.f_target)).unwrap_or_else(|error| {
        panic!(
            "CMA generated a non-finite query at generation {} candidate {} component {} (bits 0x{:016x})",
            error.generation, error.candidate, error.component, error.bits
        )
    })
}

#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
fn cmaes_with_stop_target<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
    target: Option<f64>,
) -> Result<(CmaReport, CmaStopReason), CmaGeneratedPointError> {
    let n = x0.len();
    assert!(n >= 1, "dimension must be positive");
    if let Some((component, value)) = x0.iter().enumerate().find(|(_, value)| !value.is_finite()) {
        return Err(CmaGeneratedPointError {
            generation: 0,
            candidate: 0,
            component,
            bits: value.to_bits(),
        });
    }
    let lambda = params.lambda.max(4);
    let mu = lambda / 2;
    // Log-rank recombination weights (Hansen standard).
    let raw: Vec<f64> = (0..mu)
        .map(|i| {
            fs_math::det::ln(f64::midpoint(lambda as f64, 1.0)) - fs_math::det::ln(i as f64 + 1.0)
        })
        .collect();
    let wsum: f64 = raw.iter().sum();
    let weights: Vec<f64> = raw.iter().map(|w| w / wsum).collect();
    let mu_eff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let nf = n as f64;
    // Standard strategy parameters (the IGO/natural-gradient couplings).
    let cc = (4.0 + mu_eff / nf) / (nf + 4.0 + 2.0 * mu_eff / nf);
    let cs = (mu_eff + 2.0) / (nf + mu_eff + 5.0);
    let c1 = 2.0 / ((nf + 1.3) * (nf + 1.3) + mu_eff);
    let cmu =
        (1.0 - c1).min(2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((nf + 2.0) * (nf + 2.0) + mu_eff));
    let damps = 1.0 + 2.0 * (fs_math::det::sqrt((mu_eff - 1.0) / (nf + 1.0)) - 1.0).max(0.0) + cs;
    // E‖N(0,I)‖ (Hansen's approximation).
    let chi_n = fs_math::det::sqrt(nf) * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));

    let mut mean = x0.to_vec();
    let mut sigma = params.sigma0;
    let mut cov = vec![0.0f64; n * n];
    for i in 0..n {
        cov[i * n + i] = 1.0;
    }
    let mut p_c = vec![0.0f64; n];
    let mut p_s = vec![0.0f64; n];
    // Eigen state: C = B·diag(d²)·Bᵀ; sqrt factors refreshed on cadence.
    let mut b_mat = cov.clone();
    let mut d_sqrt = vec![1.0f64; n];
    let mut stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 0,
    }
    .stream();

    let mut x_best = mean.clone();
    let mut f_best = f(&mean);
    let mut evals = 1usize;
    let mut generations = 0usize;
    if target.is_some_and(|target| f_best <= target) {
        return Ok((
            CmaReport {
                x_best,
                f_best,
                evals,
                generations,
                converged: true,
                sigma,
            },
            CmaStopReason::TargetReached,
        ));
    }
    let mut stop_reason = CmaStopReason::BudgetExhausted;
    // TolFun stagnation: generations since a meaningful f_best improvement.
    let mut gens_since_improve = 0usize;

    let mut zs: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut ys: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut fitness: Vec<f64> = vec![0.0; lambda];

    while evals
        .checked_add(lambda)
        .is_some_and(|next_generation| next_generation <= params.max_evals)
    {
        generations += 1;
        let mut generation_reached_target = false;
        // Refresh eigendecomposition on cadence (SPD maintenance).
        if generations % params.eigen_interval.max(1) == 1 || params.eigen_interval <= 1 {
            // Symmetrize (roundoff hygiene) then eigh; floor eigenvalues.
            for i in 0..n {
                for j in i + 1..n {
                    let avg = f64::midpoint(cov[i * n + j], cov[j * n + i]);
                    cov[i * n + j] = avg;
                    cov[j * n + i] = avg;
                }
            }
            let (vals, vecs) = jacobi_eigh(&cov, n);
            let vmax = vals.last().copied().unwrap_or(1.0).max(f64::MIN_POSITIVE);
            for (k, &v) in vals.iter().enumerate() {
                d_sqrt[k] = fs_math::det::sqrt(v.max(1e-14 * vmax));
            }
            b_mat.copy_from_slice(&vecs);
        }
        // Sample λ candidates: x = m + σ·B·diag(d)·z.
        for (k, z) in zs.iter_mut().enumerate() {
            for zi in z.iter_mut() {
                *zi = stream.next_normal();
            }
            let y = &mut ys[k];
            for i in 0..n {
                let mut acc = 0.0f64;
                for j in 0..n {
                    acc = (b_mat[i * n + j] * d_sqrt[j]).mul_add(z[j], acc);
                }
                y[i] = acc;
            }
            let x: Vec<f64> = mean
                .iter()
                .zip(y.iter())
                .map(|(m, yi)| sigma.mul_add(*yi, *m))
                .collect();
            if let Some((component, value)) =
                x.iter().enumerate().find(|(_, value)| !value.is_finite())
            {
                return Err(CmaGeneratedPointError {
                    generation: generations,
                    candidate: k,
                    component,
                    bits: value.to_bits(),
                });
            }
            fitness[k] = f(&x);
            evals += 1;
            generation_reached_target |= target.is_some_and(|target| fitness[k] <= target);
            if fitness[k].total_cmp(&f_best).is_lt() {
                if f_best - fitness[k] > 1e-12 * (1.0 + f_best.abs()) {
                    gens_since_improve = 0;
                }
                f_best = fitness[k];
                x_best = x;
            }
        }
        gens_since_improve += 1;
        if generation_reached_target {
            return Ok((
                CmaReport {
                    x_best,
                    f_best,
                    evals,
                    generations,
                    converged: true,
                    sigma,
                },
                CmaStopReason::TargetReached,
            ));
        }
        // Rank (total_cmp, lowest index on ties — P2).
        let mut order: Vec<usize> = (0..lambda).collect();
        order.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]).then(a.cmp(&b)));
        // Weighted recombination in y-space.
        let mut y_w = vec![0.0f64; n];
        for (w, &idx) in weights.iter().zip(&order) {
            for i in 0..n {
                y_w[i] = w.mul_add(ys[idx][i], y_w[i]);
            }
        }
        // Mean update.
        for i in 0..n {
            mean[i] = sigma.mul_add(y_w[i], mean[i]);
        }
        // CSA path: p_s ← (1−cs)p_s + √(cs(2−cs)µeff)·C^{−1/2}·y_w,
        // with C^{−1/2} = B·diag(1/d)·Bᵀ.
        let mut c_inv_half_yw = vec![0.0f64; n];
        for i in 0..n {
            // t = Bᵀ y_w
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[j * n + i].mul_add(y_w[j], acc);
            }
            c_inv_half_yw[i] = acc / d_sqrt[i];
        }
        let mut tmp = vec![0.0f64; n];
        for i in 0..n {
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[i * n + j].mul_add(c_inv_half_yw[j], acc);
            }
            tmp[i] = acc;
        }
        let csn = fs_math::det::sqrt(cs * (2.0 - cs) * mu_eff);
        for i in 0..n {
            p_s[i] = (1.0 - cs).mul_add(p_s[i], csn * tmp[i]);
        }
        let ps_norm = fs_math::det::sqrt(p_s.iter().map(|t| t * t).sum::<f64>());
        // Step-size update (the natural-gradient-consistent coupling).
        sigma *= fs_math::det::exp((cs / damps) * (ps_norm / chi_n - 1.0));
        // STAGNATION STOP: once the search distribution has collapsed
        // (σ·√λmax(C) negligible vs σ₀) the run is dead — keep sampling
        // and it just burns budget polishing whatever basin it's in.
        // BIPOP's restart ladder DEPENDS on dead runs terminating
        // (measured during bring-up: without this, failed runs consumed
        // their entire 120k budget at f ≈ 1 on rastrigin).
        let spread = sigma * d_sqrt.iter().fold(0.0f64, |m, &d| m.max(d));
        if spread < 1e-12 * params.sigma0 || gens_since_improve > 120 {
            // TolX OR TolFun: σ-collapse alone fires too slowly inside a
            // per-run budget (measured: a λ=150 local-basin run burned
            // 120k evals with f stalled for hundreds of generations) —
            // the f-stall criterion is what actually frees the budget.
            stop_reason = CmaStopReason::Stagnated;
            break;
        }
        // Rank-1 path with stall indicator h_σ.
        let h_sig = ps_norm
            / fs_math::det::sqrt(
                1.0 - fs_math::det::powi(
                    1.0 - cs,
                    2 * i32::try_from(generations.min(100_000)).expect("generation count"),
                ),
            )
            < (1.4 + 2.0 / (nf + 1.0)) * chi_n;
        let ccn = fs_math::det::sqrt(cc * (2.0 - cc) * mu_eff);
        for i in 0..n {
            let h = if h_sig { ccn * y_w[i] } else { 0.0 };
            p_c[i] = (1.0 - cc).mul_add(p_c[i], h);
        }
        // Covariance update: rank-1 + rank-µ.
        let delta_h = if h_sig { 0.0 } else { cc * (2.0 - cc) };
        for i in 0..n {
            for j in 0..n {
                let mut rank_mu = 0.0f64;
                for (w, &idx) in weights.iter().zip(&order) {
                    rank_mu = (w * ys[idx][i]).mul_add(ys[idx][j], rank_mu);
                }
                let rank1 = p_c[i] * p_c[j];
                cov[i * n + j] = (1.0 - c1 - cmu).mul_add(
                    cov[i * n + j],
                    c1.mul_add(rank1 + delta_h * cov[i * n + j], cmu * rank_mu),
                );
            }
        }
    }
    Ok((
        CmaReport {
            x_best,
            f_best,
            evals,
            generations,
            converged: false,
            sigma,
        },
        stop_reason,
    ))
}

/// Which BIPOP budget lane launched a restart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BipopLane {
    /// The doubling population ladder.
    Large,
    /// The base-population interleave.
    Small,
}

#[derive(Debug, Clone, Copy)]
struct BipopRootIdentitySource<'a> {
    start: &'a [f64],
    sigma: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
}

#[allow(dead_code)]
fn classify_bipop_root_identity_fields(source: &BipopRootIdentitySource<'_>) {
    let BipopRootIdentitySource {
        start: _,
        sigma: _,
        total_budget: _,
        target: _,
        seed: _,
    } = source;
}

fn bipop_root_identity(source: &BipopRootIdentitySource<'_>) -> Result<ReplayIdentity, BipopError> {
    let dimension = u64::try_from(source.start.len())
        .map_err(|_| BipopError::IdentityFieldOverflow { what: "dimension" })?;
    let identity_budget =
        u64::try_from(source.total_budget).map_err(|_| BipopError::IdentityFieldOverflow {
            what: "total budget",
        })?;
    let mut builder = IdentityBuilder::new(BIPOP_ROOT_IDENTITY_KIND)
        .u64(
            "report-schema-version",
            u64::from(BIPOP_REPORT_SCHEMA_VERSION),
        )
        .u64(
            "admission-schema-version",
            u64::from(BIPOP_ADMISSION_SCHEMA_VERSION),
        )
        .u64("dimension", dimension)
        .u64("total-budget", identity_budget)
        .u64("root-seed", source.seed)
        .f64_bits("initial-sigma", source.sigma)
        .flag("target-present", source.target.is_some());
    if let Some(target) = source.target {
        builder = builder.f64_bits("target", target);
    }
    for coordinate in source.start {
        builder = builder.f64_bits("start-coordinate", *coordinate);
    }
    Ok(builder.finish())
}

/// Exact causal inputs for one BIPOP study.
///
/// The nested [`ReplayIdentity`] is a canonical, typed fs-obs preimage over
/// every field. It authenticates the retained fields against accidental or
/// stale mutation; callers that require an external study identity must still
/// compare it with their separately retained expected root.
#[derive(Debug, Clone)]
pub struct BipopRootInputs {
    start: Vec<f64>,
    sigma: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
    identity: ReplayIdentity,
}

impl BipopRootInputs {
    /// Exact initial point supplied to BIPOP.
    #[must_use]
    pub fn start(&self) -> &[f64] {
        &self.start
    }

    /// Initial CMA step size, retained by exact IEEE-754 bits.
    #[must_use]
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Hard aggregate callback budget.
    #[must_use]
    pub fn total_budget(&self) -> usize {
        self.total_budget
    }

    /// Typed finite target, or `None` when target stopping is disabled.
    #[must_use]
    pub fn target(&self) -> Option<f64> {
        self.target
    }

    /// Root seed from which restart and CMA streams are derived.
    #[must_use]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Canonical identity binding every retained root-input bit.
    #[must_use]
    pub fn identity(&self) -> &ReplayIdentity {
        &self.identity
    }
}

/// Strong streaming identity over the exact root-bound callback trace.
///
/// The BLAKE3 preimage starts with the domain byte length plus exact domain,
/// then includes the trace schema, canonical root-input identity preimage,
/// dimension, row count, every row's restart ownership/local offset/objective
/// bits, and every decision bit in global callback order. The preimage is
/// streamed and is not duplicated in memory beside the retained trace.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BipopTraceIdentity {
    schema_version: u32,
    rows: usize,
    dimension: usize,
    digest: [u8; 32],
}

impl BipopTraceIdentity {
    /// Trace-identity schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Number of callbacks committed by the identity.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Decision dimension committed for every row.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Raw 256-bit BLAKE3 digest.
    #[must_use]
    pub fn digest(&self) -> &[u8; 32] {
        &self.digest
    }
}

fn build_bipop_root_inputs(
    start: &[f64],
    sigma: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
) -> Result<BipopRootInputs, BipopError> {
    let source = BipopRootIdentitySource {
        start,
        sigma,
        total_budget,
        target,
        seed,
    };
    Ok(BipopRootInputs {
        start: start.to_vec(),
        sigma,
        total_budget,
        target,
        seed,
        identity: bipop_root_identity(&source)?,
    })
}

#[derive(Debug, Clone)]
struct BipopEvaluationRow {
    schema_version: u32,
    restart: u64,
    local_offset: usize,
    objective: f64,
}

#[derive(Debug, Clone, Copy)]
struct BipopTraceIdentitySource<'a> {
    root_canonical_bytes: &'a [u8],
    dimension: usize,
    rows: &'a [BipopEvaluationRow],
    points: &'a [f64],
}

#[allow(dead_code)]
fn classify_bipop_trace_identity_fields(source: &BipopTraceIdentitySource<'_>) {
    let BipopTraceIdentitySource {
        root_canonical_bytes: _,
        dimension: _,
        rows: _,
        points: _,
    } = source;
}

fn build_bipop_trace_identity(
    root: &BipopRootInputs,
    rows: &[BipopEvaluationRow],
    points: &[f64],
) -> Result<BipopTraceIdentity, BipopError> {
    let source = BipopTraceIdentitySource {
        root_canonical_bytes: root.identity.canonical_bytes(),
        dimension: root.start.len(),
        rows,
        points,
    };
    bipop_trace_identity(&source)
}

fn bipop_trace_identity(
    source: &BipopTraceIdentitySource<'_>,
) -> Result<BipopTraceIdentity, BipopError> {
    if source.dimension == 0 {
        return Err(BipopError::InternalInvariant {
            what: "callback trace dimension must be positive",
        });
    }
    let encoded_dimension =
        u64::try_from(source.dimension).map_err(|_| BipopError::IdentityFieldOverflow {
            what: "trace dimension",
        })?;
    let encoded_rows =
        u64::try_from(source.rows.len()).map_err(|_| BipopError::IdentityFieldOverflow {
            what: "trace row count",
        })?;
    let encoded_domain_bytes = u64::try_from(BIPOP_TRACE_IDENTITY_DOMAIN.len()).map_err(|_| {
        BipopError::IdentityFieldOverflow {
            what: "trace identity domain byte length",
        }
    })?;
    let encoded_root_bytes = u64::try_from(source.root_canonical_bytes.len()).map_err(|_| {
        BipopError::IdentityFieldOverflow {
            what: "root identity byte length",
        }
    })?;
    let expected_points =
        source
            .rows
            .len()
            .checked_mul(source.dimension)
            .ok_or(BipopError::DimensionOverflow {
                what: "callback trace coordinate entries",
            })?;
    if source.points.len() != expected_points {
        return Err(BipopError::InternalInvariant {
            what: "callback trace points must form a dense row-major matrix",
        });
    }

    let mut hasher = Blake3::new();
    hasher.update(&encoded_domain_bytes.to_le_bytes());
    hasher.update(BIPOP_TRACE_IDENTITY_DOMAIN.as_bytes());
    hasher.update(&BIPOP_TRACE_IDENTITY_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&encoded_dimension.to_le_bytes());
    hasher.update(&encoded_rows.to_le_bytes());
    hasher.update(&encoded_root_bytes.to_le_bytes());
    hasher.update(source.root_canonical_bytes);
    for (row, point) in source
        .rows
        .iter()
        .zip(source.points.chunks_exact(source.dimension))
    {
        let encoded_local =
            u64::try_from(row.local_offset).map_err(|_| BipopError::IdentityFieldOverflow {
                what: "trace local offset",
            })?;
        hasher.update(&row.schema_version.to_le_bytes());
        hasher.update(&row.restart.to_le_bytes());
        hasher.update(&encoded_local.to_le_bytes());
        hasher.update(&row.objective.to_bits().to_le_bytes());
        for coordinate in point {
            hasher.update(&coordinate.to_bits().to_le_bytes());
        }
    }
    Ok(BipopTraceIdentity {
        schema_version: BIPOP_TRACE_IDENTITY_SCHEMA_VERSION,
        rows: source.rows.len(),
        dimension: source.dimension,
        digest: *hasher.finalize().as_bytes(),
    })
}

/// Borrowed view of one exact production objective invocation.
///
/// Global ordering is the row order in [`BipopReport::evaluations`]. The
/// decision slice borrows the report's flat coordinate store, avoiding one
/// allocation per callback while retaining every queried bit.
#[derive(Debug, Clone, Copy)]
pub struct BipopEvaluationRecord<'a> {
    schema_version: u32,
    global_offset: usize,
    restart: u64,
    local_offset: usize,
    point: &'a [f64],
    objective: f64,
}

impl BipopEvaluationRecord<'_> {
    /// Evaluation-record schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Zero-based position in the complete study trace.
    #[must_use]
    pub fn global_offset(&self) -> usize {
        self.global_offset
    }

    /// Restart ordinal that owns this callback.
    #[must_use]
    pub fn restart(&self) -> u64 {
        self.restart
    }

    /// Zero-based callback offset inside the owning restart.
    #[must_use]
    pub fn local_offset(&self) -> usize {
        self.local_offset
    }

    /// Exact finite decision point supplied to the objective.
    #[must_use]
    pub fn point(&self) -> &[f64] {
        self.point
    }

    /// Exact objective-result bits returned by user code.
    ///
    /// Non-finite objective outputs remain data in the current CMA contract;
    /// they are retained without normalization or a false finiteness claim.
    #[must_use]
    pub fn objective(&self) -> f64 {
        self.objective
    }
}

/// One immutable, versioned BIPOP restart receipt.
///
/// Point and objective values retain their exact `f64` bits. The aggregate
/// trace interval is half-open and indexes the exact retained production
/// callback records exposed by [`BipopReport::evaluations`].
#[derive(Debug, Clone)]
pub struct BipopRestartRecord {
    schema_version: u32,
    ordinal: u64,
    lane: BipopLane,
    lambda: usize,
    allocated_budget: usize,
    seed: u64,
    start: Vec<f64>,
    trace_start: usize,
    trace_end: usize,
    stop_reason: CmaStopReason,
    report: CmaReport,
}

impl BipopRestartRecord {
    /// Restart-record schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Zero-based restart ordinal.
    #[must_use]
    pub fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Large or small BIPOP budget lane.
    #[must_use]
    pub fn lane(&self) -> BipopLane {
        self.lane
    }

    /// Population size used by this restart.
    #[must_use]
    pub fn lambda(&self) -> usize {
        self.lambda
    }

    /// Local evaluation cap assigned to this restart.
    #[must_use]
    pub fn allocated_budget(&self) -> usize {
        self.allocated_budget
    }

    /// CMA stream seed derived for this restart.
    #[must_use]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Exact start point supplied to this restart.
    #[must_use]
    pub fn start(&self) -> &[f64] {
        &self.start
    }

    /// Start of this restart's half-open aggregate evaluation interval.
    #[must_use]
    pub fn trace_start(&self) -> usize {
        self.trace_start
    }

    /// End of this restart's half-open aggregate evaluation interval.
    #[must_use]
    pub fn trace_end(&self) -> usize {
        self.trace_end
    }

    /// Causal terminal classification retained from the CMA run.
    #[must_use]
    pub fn stop_reason(&self) -> CmaStopReason {
        self.stop_reason
    }

    /// Complete CMA result for this restart.
    #[must_use]
    pub fn report(&self) -> &CmaReport {
        &self.report
    }
}

/// Structured refusal from [`BipopReport::validate_ledger`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BipopLedgerError {
    restart: Option<usize>,
    invariant: &'static str,
}

impl BipopLedgerError {
    fn global(invariant: &'static str) -> Self {
        Self {
            restart: None,
            invariant,
        }
    }

    fn at(restart: usize, invariant: &'static str) -> Self {
        Self {
            restart: Some(restart),
            invariant,
        }
    }

    /// Restart index associated with the refusal, if it is local.
    #[must_use]
    pub fn restart(&self) -> Option<usize> {
        self.restart
    }

    /// Stable invariant name suitable for structured diagnostics.
    #[must_use]
    pub fn invariant(&self) -> &'static str {
        self.invariant
    }
}

impl core::fmt::Display for BipopLedgerError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self.restart {
            Some(restart) => write!(
                formatter,
                "BIPOP restart {restart} violates {}",
                self.invariant
            ),
            None => write!(formatter, "BIPOP ledger violates {}", self.invariant),
        }
    }
}

impl std::error::Error for BipopLedgerError {}

fn cma_reports_match_bits(left: &CmaReport, right: &CmaReport) -> bool {
    left.f_best.to_bits() == right.f_best.to_bits()
        && left.evals == right.evals
        && left.generations == right.generations
        && left.converged == right.converged
        && left.sigma.to_bits() == right.sigma.to_bits()
        && left.x_best.len() == right.x_best.len()
        && left
            .x_best
            .iter()
            .zip(&right.x_best)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

fn f64_slices_match_bits(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

/// Strong identity over one complete BIPOP result payload.
///
/// The digest composes the exact root-input preimage, the retained and
/// independently recomputed callback-trace receipts, every compatibility
/// projection, and every ordered restart-record bit. The two trace receipts
/// are deliberately distinct fields: stale trace content and a stale retained
/// trace receipt must both perturb the full-study identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BipopStudyIdentity {
    schema_version: u32,
    restarts: usize,
    evaluations: usize,
    digest: [u8; 32],
}

impl BipopStudyIdentity {
    /// Full-study identity schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Number of ordered restart records bound by the digest.
    #[must_use]
    pub fn restarts(&self) -> usize {
        self.restarts
    }

    /// Number of ordered objective callbacks bound by the digest.
    #[must_use]
    pub fn evaluations(&self) -> usize {
        self.evaluations
    }

    /// Raw 256-bit BLAKE3 digest.
    #[must_use]
    pub fn digest(&self) -> &[u8; 32] {
        &self.digest
    }
}

impl core::fmt::Display for BipopStudyIdentity {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "v{}:r{}:e{}:",
            self.schema_version, self.restarts, self.evaluations
        )?;
        for byte in self.digest {
            write!(formatter, "{byte:02x}")?;
        }
        Ok(())
    }
}

/// Typed refusal from production full-study identity validation/admission.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BipopStudyAdmissionError {
    /// The retained identity does not match the report payload currently held.
    PayloadIdentityMismatch {
        /// Identity retained in the report.
        declared: BipopStudyIdentity,
        /// Identity recomputed from the current report payload.
        computed: BipopStudyIdentity,
    },
    /// An internally valid report names a different study than the caller's
    /// separately retained reference.
    ReferenceIdentityMismatch {
        /// External study identity required by the caller.
        expected: BipopStudyIdentity,
        /// Internally valid identity carried by this report.
        found: BipopStudyIdentity,
    },
    /// A native cardinality could not enter the fixed-width identity frame.
    IdentityEncoding {
        /// Stable field label.
        field: &'static str,
    },
    /// The retained payload cannot form the canonical identity preimage.
    PayloadInvalid {
        /// Stable failed preimage invariant.
        invariant: &'static str,
    },
    /// Structural or semantic ledger admission failed before reference matching.
    Ledger {
        /// First deterministic ledger refusal.
        error: BipopLedgerError,
    },
}

impl core::fmt::Display for BipopStudyAdmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PayloadIdentityMismatch { declared, computed } => write!(
                formatter,
                "BIPOP PayloadIdentityMismatch: declared {declared}, computed {computed}"
            ),
            Self::ReferenceIdentityMismatch { expected, found } => write!(
                formatter,
                "BIPOP ReferenceIdentityMismatch: expected {expected}, found {found}"
            ),
            Self::IdentityEncoding { field } => write!(
                formatter,
                "BIPOP study identity field `{field}` does not fit its canonical encoding"
            ),
            Self::PayloadInvalid { invariant } => write!(
                formatter,
                "BIPOP study identity payload violates `{invariant}`"
            ),
            Self::Ledger { error } => write!(formatter, "BIPOP study admission refused: {error}"),
        }
    }
}

impl std::error::Error for BipopStudyAdmissionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Ledger { error } => Some(error),
            _ => None,
        }
    }
}

/// Typed refusal from replay-backed BIPOP semantic admission.
///
/// This is intentionally distinct from [`BipopStudyAdmissionError`]. The cheap
/// identity/ledger gate authenticates retained bytes against an external
/// reference without invoking user code. Replay-backed admission performs that
/// gate first, then executes a caller-supplied objective oracle from the exact
/// retained root and requires the independently produced complete study to have
/// the same identity.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BipopReplayAdmissionError {
    /// The callback-free payload, ledger, or external-reference gate refused.
    EvidenceAdmission {
        /// Exact cheap-admission refusal. The replay oracle was not invoked.
        error: BipopStudyAdmissionError,
    },
    /// Re-executing the retained root against the supplied oracle failed.
    ReplayExecution {
        /// Exact production execution refusal from the replay attempt.
        error: BipopError,
    },
    /// Replay completed and validated, but produced a different complete study.
    SemanticMismatch {
        /// Identity of the admitted retained report.
        retained: BipopStudyIdentity,
        /// Identity independently produced by replaying the objective oracle.
        replayed: BipopStudyIdentity,
    },
}

impl core::fmt::Display for BipopReplayAdmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EvidenceAdmission { error } => {
                write!(
                    formatter,
                    "BIPOP replay admission refused retained evidence: {error}"
                )
            }
            Self::ReplayExecution { error } => {
                write!(
                    formatter,
                    "BIPOP semantic replay failed to execute: {error}"
                )
            }
            Self::SemanticMismatch { retained, replayed } => write!(
                formatter,
                "BIPOP semantic replay mismatch: retained {retained}, replayed {replayed}"
            ),
        }
    }
}

impl std::error::Error for BipopReplayAdmissionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::EvidenceAdmission { error } => Some(error),
            Self::ReplayExecution { error } => Some(error),
            Self::SemanticMismatch { .. } => None,
        }
    }
}

fn study_identity_usize(value: usize, what: &'static str) -> Result<u64, BipopError> {
    u64::try_from(value).map_err(|_| BipopError::IdentityFieldOverflow { what })
}

fn study_hash_field(
    hasher: &mut Blake3,
    label: &'static str,
    value: &[u8],
) -> Result<(), BipopError> {
    let label_len = study_identity_usize(label.len(), "study field label length")?;
    let value_len = study_identity_usize(value.len(), "study field value length")?;
    hasher.update(&label_len.to_le_bytes());
    hasher.update(label.as_bytes());
    hasher.update(&value_len.to_le_bytes());
    hasher.update(value);
    Ok(())
}

fn study_hash_u64(hasher: &mut Blake3, label: &'static str, value: u64) -> Result<(), BipopError> {
    study_hash_field(hasher, label, &value.to_le_bytes())
}

fn study_hash_usize(
    hasher: &mut Blake3,
    label: &'static str,
    value: usize,
) -> Result<(), BipopError> {
    study_hash_u64(hasher, label, study_identity_usize(value, label)?)
}

fn study_hash_u32(hasher: &mut Blake3, label: &'static str, value: u32) -> Result<(), BipopError> {
    study_hash_field(hasher, label, &value.to_le_bytes())
}

fn study_hash_f64(hasher: &mut Blake3, label: &'static str, value: f64) -> Result<(), BipopError> {
    study_hash_u64(hasher, label, value.to_bits())
}

fn study_hash_flag(
    hasher: &mut Blake3,
    label: &'static str,
    value: bool,
) -> Result<(), BipopError> {
    study_hash_field(hasher, label, &[u8::from(value)])
}

fn bipop_lane_tag(lane: BipopLane) -> u8 {
    match lane {
        BipopLane::Large => 0,
        BipopLane::Small => 1,
    }
}

fn cma_stop_reason_tag(reason: CmaStopReason) -> u8 {
    match reason {
        CmaStopReason::TargetReached => 0,
        CmaStopReason::BudgetExhausted => 1,
        CmaStopReason::Stagnated => 2,
    }
}

#[derive(Debug, Clone, Copy)]
struct BipopStudyIdentitySource<'a> {
    report_schema_version: u32,
    root: &'a BipopRootInputs,
    root_content_identity: &'a ReplayIdentity,
    schedule: &'a [usize],
    total_evals: usize,
    records: &'a [BipopRestartRecord],
    best_restart: usize,
    best: &'a CmaReport,
    retained_trace_identity: BipopTraceIdentity,
    callback_content_identity: BipopTraceIdentity,
}

#[allow(dead_code)]
fn classify_bipop_study_identity_fields(source: &BipopStudyIdentitySource<'_>) {
    let BipopStudyIdentitySource {
        report_schema_version: _,
        root: _,
        root_content_identity: _,
        schedule: _,
        total_evals: _,
        records: _,
        best_restart: _,
        best: _,
        retained_trace_identity: _,
        callback_content_identity: _,
    } = source;
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn build_bipop_study_identity(
    report_schema_version: u32,
    root: &BipopRootInputs,
    root_content_identity: &ReplayIdentity,
    schedule: &[usize],
    total_evals: usize,
    records: &[BipopRestartRecord],
    best_restart: usize,
    best: &CmaReport,
    retained_trace_identity: BipopTraceIdentity,
    callback_content_identity: BipopTraceIdentity,
) -> Result<BipopStudyIdentity, BipopError> {
    let source = BipopStudyIdentitySource {
        report_schema_version,
        root,
        root_content_identity,
        schedule,
        total_evals,
        records,
        best_restart,
        best,
        retained_trace_identity,
        callback_content_identity,
    };
    bipop_study_identity(&source)
}

#[allow(clippy::too_many_lines)] // One ordered encoder keeps every retained projection visibly bound.
fn bipop_study_identity(
    source: &BipopStudyIdentitySource<'_>,
) -> Result<BipopStudyIdentity, BipopError> {
    let BipopStudyIdentitySource {
        report_schema_version,
        root,
        root_content_identity,
        schedule,
        total_evals,
        records,
        best_restart,
        best,
        retained_trace_identity,
        callback_content_identity,
    } = *source;
    let mut hasher = Blake3::new();
    study_hash_field(
        &mut hasher,
        "domain",
        BIPOP_STUDY_IDENTITY_DOMAIN.as_bytes(),
    )?;
    study_hash_u32(
        &mut hasher,
        "study-schema-version",
        BIPOP_STUDY_IDENTITY_SCHEMA_VERSION,
    )?;
    study_hash_u32(&mut hasher, "report-schema-version", report_schema_version)?;
    study_hash_u32(
        &mut hasher,
        "admission-schema-version",
        BIPOP_ADMISSION_SCHEMA_VERSION,
    )?;
    study_hash_u32(
        &mut hasher,
        "restart-schema-version",
        BIPOP_RESTART_SCHEMA_VERSION,
    )?;
    study_hash_u32(
        &mut hasher,
        "evaluation-schema-version",
        BIPOP_EVALUATION_SCHEMA_VERSION,
    )?;
    study_hash_u32(
        &mut hasher,
        "trace-schema-version",
        BIPOP_TRACE_IDENTITY_SCHEMA_VERSION,
    )?;
    study_hash_u32(
        &mut hasher,
        "rand-stream-semantics-version",
        fs_rand::STREAM_SEMANTICS_VERSION,
    )?;
    study_hash_u32(
        &mut hasher,
        "jacobi-admission-schema-version",
        fs_la::eigen::JACOBI_EIGH_ADMISSION_SCHEMA_VERSION,
    )?;
    study_hash_u32(&mut hasher, "cma-stream-kernel", K_CMA)?;
    study_hash_u64(
        &mut hasher,
        "restart-seed-stride",
        BIPOP_RESTART_SEED_STRIDE,
    )?;
    study_hash_u32(&mut hasher, "large-run-cap", BIPOP_LARGE_RUN_CAP)?;
    study_hash_usize(
        &mut hasher,
        "generations-per-restart",
        BIPOP_GENERATIONS_PER_RESTART,
    )?;
    study_hash_field(
        &mut hasher,
        "retained-root-canonical-bytes",
        root.identity.canonical_bytes(),
    )?;
    study_hash_field(
        &mut hasher,
        "root-content-canonical-bytes",
        root_content_identity.canonical_bytes(),
    )?;
    study_hash_u32(
        &mut hasher,
        "retained-trace-schema-version",
        retained_trace_identity.schema_version,
    )?;
    study_hash_usize(
        &mut hasher,
        "retained-trace-rows",
        retained_trace_identity.rows,
    )?;
    study_hash_usize(
        &mut hasher,
        "retained-trace-dimension",
        retained_trace_identity.dimension,
    )?;
    study_hash_field(
        &mut hasher,
        "retained-trace-digest",
        &retained_trace_identity.digest,
    )?;
    study_hash_u32(
        &mut hasher,
        "callback-content-schema-version",
        callback_content_identity.schema_version,
    )?;
    study_hash_usize(
        &mut hasher,
        "callback-content-rows",
        callback_content_identity.rows,
    )?;
    study_hash_usize(
        &mut hasher,
        "callback-content-dimension",
        callback_content_identity.dimension,
    )?;
    study_hash_field(
        &mut hasher,
        "callback-content-digest",
        &callback_content_identity.digest,
    )?;
    study_hash_usize(&mut hasher, "total-evals", total_evals)?;
    study_hash_usize(&mut hasher, "schedule-length", schedule.len())?;
    for (index, lambda) in schedule.iter().copied().enumerate() {
        study_hash_usize(&mut hasher, "schedule-index", index)?;
        study_hash_usize(&mut hasher, "schedule-lambda", lambda)?;
    }
    study_hash_usize(&mut hasher, "restart-record-count", records.len())?;
    study_hash_usize(&mut hasher, "best-restart", best_restart)?;
    study_hash_usize(&mut hasher, "best-x-length", best.x_best.len())?;
    for (coordinate, value) in best.x_best.iter().copied().enumerate() {
        study_hash_usize(&mut hasher, "best-coordinate-index", coordinate)?;
        study_hash_f64(&mut hasher, "best-coordinate", value)?;
    }
    study_hash_f64(&mut hasher, "best-objective", best.f_best)?;
    study_hash_usize(&mut hasher, "best-evaluations", best.evals)?;
    study_hash_usize(&mut hasher, "best-generations", best.generations)?;
    study_hash_flag(&mut hasher, "best-converged", best.converged)?;
    study_hash_f64(&mut hasher, "best-sigma", best.sigma)?;

    for (record_index, record) in records.iter().enumerate() {
        study_hash_usize(&mut hasher, "record-index", record_index)?;
        study_hash_u32(&mut hasher, "record-schema-version", record.schema_version)?;
        study_hash_u64(&mut hasher, "record-ordinal", record.ordinal)?;
        study_hash_field(&mut hasher, "record-lane", &[bipop_lane_tag(record.lane)])?;
        study_hash_usize(&mut hasher, "record-lambda", record.lambda)?;
        study_hash_usize(
            &mut hasher,
            "record-allocated-budget",
            record.allocated_budget,
        )?;
        study_hash_u64(&mut hasher, "record-seed", record.seed)?;
        study_hash_usize(&mut hasher, "record-start-length", record.start.len())?;
        for (coordinate, value) in record.start.iter().copied().enumerate() {
            study_hash_usize(&mut hasher, "record-start-coordinate-index", coordinate)?;
            study_hash_f64(&mut hasher, "record-start-coordinate", value)?;
        }
        study_hash_usize(&mut hasher, "record-trace-start", record.trace_start)?;
        study_hash_usize(&mut hasher, "record-trace-end", record.trace_end)?;
        study_hash_field(
            &mut hasher,
            "record-stop-reason",
            &[cma_stop_reason_tag(record.stop_reason)],
        )?;
        study_hash_usize(
            &mut hasher,
            "record-report-x-length",
            record.report.x_best.len(),
        )?;
        for (coordinate, value) in record.report.x_best.iter().copied().enumerate() {
            study_hash_usize(&mut hasher, "record-report-coordinate-index", coordinate)?;
            study_hash_f64(&mut hasher, "record-report-coordinate", value)?;
        }
        study_hash_f64(&mut hasher, "record-report-objective", record.report.f_best)?;
        study_hash_usize(
            &mut hasher,
            "record-report-evaluations",
            record.report.evals,
        )?;
        study_hash_usize(
            &mut hasher,
            "record-report-generations",
            record.report.generations,
        )?;
        study_hash_flag(
            &mut hasher,
            "record-report-converged",
            record.report.converged,
        )?;
        study_hash_f64(&mut hasher, "record-report-sigma", record.report.sigma)?;
    }

    Ok(BipopStudyIdentity {
        schema_version: BIPOP_STUDY_IDENTITY_SCHEMA_VERSION,
        restarts: records.len(),
        evaluations: total_evals,
        digest: *hasher.finalize().as_bytes(),
    })
}

/// BIPOP restart evidence.
#[derive(Debug, Clone)]
pub struct BipopReport {
    /// Compatibility projection of [`Self::best_record`].
    pub best: CmaReport,
    /// Compatibility projection of every restart's population size.
    pub schedule: Vec<usize>,
    /// Compatibility projection of the terminal aggregate trace offset.
    pub total_evals: usize,
    schema_version: u32,
    root: BipopRootInputs,
    records: Vec<BipopRestartRecord>,
    best_restart: usize,
    trace_rows: Vec<BipopEvaluationRow>,
    trace_points: Vec<f64>,
    trace_identity: BipopTraceIdentity,
    study_identity: BipopStudyIdentity,
}

impl BipopReport {
    /// Production report schema governing the root, restart ledger, and trace.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact, canonically identified inputs that caused this study.
    #[must_use]
    pub fn root_inputs(&self) -> &BipopRootInputs {
        &self.root
    }

    /// Ordered immutable restart ledger.
    #[must_use]
    pub fn records(&self) -> &[BipopRestartRecord] {
        &self.records
    }

    /// Admitted hard aggregate callback budget that governs every record's
    /// local allocation and the ledger's terminal completeness.
    #[must_use]
    pub fn total_budget(&self) -> usize {
        self.root.total_budget
    }

    /// Index of the earliest restart attaining the best objective under
    /// `f64::total_cmp`.
    #[must_use]
    pub fn best_restart(&self) -> usize {
        self.best_restart
    }

    /// Named record from which [`Self::best`] is projected.
    #[must_use]
    pub fn best_record(&self) -> Option<&BipopRestartRecord> {
        self.records.get(self.best_restart)
    }

    /// One retained objective invocation by global trace offset.
    #[must_use]
    pub fn evaluation(&self, global_offset: usize) -> Option<BipopEvaluationRecord<'_>> {
        let row = self.trace_rows.get(global_offset)?;
        let dimension = self.root.start.len();
        let point_start = global_offset.checked_mul(dimension)?;
        let point_end = point_start.checked_add(dimension)?;
        let point = self.trace_points.get(point_start..point_end)?;
        Some(BipopEvaluationRecord {
            schema_version: row.schema_version,
            global_offset,
            restart: row.restart,
            local_offset: row.local_offset,
            point,
            objective: row.objective,
        })
    }

    /// Complete ordered production callback trace.
    #[must_use]
    pub fn evaluations(&self) -> impl ExactSizeIterator<Item = BipopEvaluationRecord<'_>> + '_ {
        (0..self.trace_rows.len()).map(|index| {
            self.evaluation(index)
                .expect("private BIPOP trace layout must remain contiguous")
        })
    }

    /// Strong identity over the exact root-bound callback trace.
    #[must_use]
    pub fn trace_identity(&self) -> BipopTraceIdentity {
        self.trace_identity
    }

    /// Strong identity over the complete root, callback trace, and restart ledger.
    #[must_use]
    pub fn study_identity(&self) -> BipopStudyIdentity {
        self.study_identity
    }

    fn computed_study_identity(&self) -> Result<BipopStudyIdentity, BipopError> {
        let root_content = build_bipop_root_inputs(
            &self.root.start,
            self.root.sigma,
            self.root.total_budget,
            self.root.target,
            self.root.seed,
        )?;
        let callback_content_identity =
            build_bipop_trace_identity(&self.root, &self.trace_rows, &self.trace_points)?;
        build_bipop_study_identity(
            self.schema_version,
            &self.root,
            &root_content.identity,
            &self.schedule,
            self.total_evals,
            &self.records,
            self.best_restart,
            &self.best,
            self.trace_identity,
            callback_content_identity,
        )
    }

    /// Recompute and validate the retained full-study identity.
    ///
    /// This check distinguishes a stale payload from comparison against a
    /// different external study. It does not replace [`Self::validate_ledger`],
    /// which supplies the independent structural and semantic checks.
    ///
    /// # Errors
    /// Returns [`BipopStudyAdmissionError::PayloadIdentityMismatch`] when the
    /// current report bytes do not match their retained identity, or a payload
    /// layout/encoding refusal when no canonical preimage can be formed.
    pub fn validate_study_identity(&self) -> Result<(), BipopStudyAdmissionError> {
        let computed = self
            .computed_study_identity()
            .map_err(|error| match error {
                BipopError::IdentityFieldOverflow { what }
                | BipopError::DimensionOverflow { what } => {
                    BipopStudyAdmissionError::IdentityEncoding { field: what }
                }
                BipopError::InternalInvariant { what } => {
                    BipopStudyAdmissionError::PayloadInvalid { invariant: what }
                }
                _ => BipopStudyAdmissionError::PayloadInvalid {
                    invariant: "study identity reconstruction",
                },
            })?;
        if computed == self.study_identity {
            Ok(())
        } else {
            Err(BipopStudyAdmissionError::PayloadIdentityMismatch {
                declared: self.study_identity,
                computed,
            })
        }
    }

    /// Admit this report against one separately retained full-study identity.
    ///
    /// Validation order is intentional: stale payloads are classified before
    /// semantic ledger failures, and an internally valid but different study is
    /// classified only after both checks succeed.
    ///
    /// This is the cheap, callback-free evidence gate. It proves that all
    /// retained fields match their strong identity and the structural ledger,
    /// but it does not reconstruct candidate generation, final CMA state, or
    /// objective semantics. Call [`Self::admit_study_identity_with_replay`] when
    /// those causal semantics must be checked against an executable oracle.
    ///
    /// # Errors
    /// Returns a typed payload, ledger, encoding, or reference mismatch.
    pub fn admit_study_identity(
        &self,
        expected: BipopStudyIdentity,
    ) -> Result<(), BipopStudyAdmissionError> {
        self.validate_study_identity()?;
        self.validate_ledger()
            .map_err(|error| BipopStudyAdmissionError::Ledger { error })?;
        if self.study_identity == expected {
            Ok(())
        } else {
            Err(BipopStudyAdmissionError::ReferenceIdentityMismatch {
                expected,
                found: self.study_identity,
            })
        }
    }

    /// Admit retained evidence and then replay its exact root against an oracle.
    ///
    /// The callback-free [`Self::admit_study_identity`] gate always runs first.
    /// Consequently, stale payloads, invalid ledgers, and wrong external
    /// references invoke `objective_oracle` zero times. After that gate succeeds,
    /// the method re-executes BIPOP with the retained start, sigma, hard budget,
    /// target, and seed. The production executor validates the replayed ledger;
    /// equality of the complete study identities then binds every replayed root,
    /// callback, restart, report, and compatibility-projection bit.
    ///
    /// The supplied oracle is executable authority: it must implement the same
    /// deterministic objective semantics claimed by the retained study. Panics
    /// are not caught, and replay consumes the full retained callback budget (or
    /// terminates under the same recorded rule). This method is therefore an
    /// expensive semantic gate rather than a replacement for cheap admission.
    ///
    /// # Errors
    /// Returns [`BipopReplayAdmissionError::EvidenceAdmission`] before invoking
    /// the oracle when cheap admission fails,
    /// [`BipopReplayAdmissionError::ReplayExecution`] when replay cannot
    /// complete, or [`BipopReplayAdmissionError::SemanticMismatch`] when a valid
    /// replay produces different complete evidence.
    pub fn admit_study_identity_with_replay<F>(
        &self,
        expected: BipopStudyIdentity,
        objective_oracle: &mut F,
    ) -> Result<(), BipopReplayAdmissionError>
    where
        F: FnMut(&[f64]) -> f64,
    {
        self.admit_study_identity(expected)
            .map_err(|error| BipopReplayAdmissionError::EvidenceAdmission { error })?;
        let replay = try_bipop_cmaes(
            objective_oracle,
            &self.root.start,
            self.root.sigma,
            self.root.total_budget,
            self.root.target,
            self.root.seed,
        )
        .map_err(|error| BipopReplayAdmissionError::ReplayExecution { error })?;
        let replayed = replay.study_identity;
        if replayed == self.study_identity {
            Ok(())
        } else {
            Err(BipopReplayAdmissionError::SemanticMismatch {
                retained: self.study_identity,
                replayed,
            })
        }
    }

    /// Recheck the ordered ledger and every compatibility projection.
    ///
    /// This is a structural validator over retained evidence. It recomputes the
    /// canonical root identity, restart derivations, complete callback-trace
    /// projection, strong trace identity, and complete study identity. A
    /// consumer asserting one specific study should use
    /// [`Self::admit_study_identity`] with its separately retained reference.
    ///
    /// # Errors
    /// Returns a [`BipopLedgerError`] naming the first deterministic invariant
    /// violation.
    #[allow(clippy::too_many_lines)] // one ordered pass mirrors the versioned record schema
    pub fn validate_ledger(&self) -> Result<(), BipopLedgerError> {
        if self.schema_version != BIPOP_REPORT_SCHEMA_VERSION {
            return Err(BipopLedgerError::global("report-schema-version"));
        }
        let first = self
            .records
            .first()
            .ok_or_else(|| BipopLedgerError::global("nonempty"))?;
        if self.schedule.len() != self.records.len() {
            return Err(BipopLedgerError::global("schedule-length"));
        }
        if self.root.total_budget == 0 {
            return Err(BipopLedgerError::global("positive-total-budget"));
        }
        let base_lambda = first.lambda;
        let base_seed = self.root.seed;
        let point_dim = self.root.start.len();
        if point_dim == 0 {
            return Err(BipopLedgerError::at(0, "positive-point-dimension"));
        }
        if self.root.start.iter().any(|value| !value.is_finite()) {
            return Err(BipopLedgerError::global("finite-root-start"));
        }
        if !self.root.sigma.is_finite() || self.root.sigma <= 0.0 {
            return Err(BipopLedgerError::global("positive-finite-root-sigma"));
        }
        if self.root.target.is_some_and(|target| !target.is_finite()) {
            return Err(BipopLedgerError::global("finite-root-target"));
        }
        let expected_root = build_bipop_root_inputs(
            &self.root.start,
            self.root.sigma,
            self.root.total_budget,
            self.root.target,
            self.root.seed,
        )
        .map_err(|_| BipopLedgerError::global("root-identity-encoding"))?;
        if self.root.identity != expected_root.identity {
            return Err(BipopLedgerError::global("root-identity"));
        }
        // Classify a non-finite first retained decision point under the same
        // stable invariant as every later restart before reporting the broader
        // root-projection mismatch that it necessarily also causes.
        if first.start.iter().any(|value| !value.is_finite()) {
            return Err(BipopLedgerError::at(0, "finite-start"));
        }
        if !f64_slices_match_bits(&first.start, &self.root.start) {
            return Err(BipopLedgerError::at(0, "root-start-projection"));
        }
        let matrix_entries = point_dim
            .checked_mul(point_dim)
            .ok_or_else(|| BipopLedgerError::at(0, "dense-matrix-admission"))?;
        checked_dense_matrix_allocation(matrix_entries)
            .map_err(|_| BipopLedgerError::at(0, "dense-matrix-admission"))?;
        checked_random_counter_blocks(
            "restart-perturbation",
            &[2, point_dim as u128, (self.records.len() - 1) as u128],
        )
        .map_err(|_| BipopLedgerError::global("restart-counter-range"))?;
        let expected_base_lambda = 4 + (3.0 * fs_math::det::ln(point_dim as f64)).floor() as usize;
        if base_lambda != expected_base_lambda {
            return Err(BipopLedgerError::at(0, "base-population"));
        }
        // Mirror callback-free admission, not merely the realized trace. An
        // early target hit cannot retroactively authenticate a problem whose
        // hard budget could have entered the refused Jacobi dependency.
        if self.root.total_budget > expected_base_lambda {
            fs_la::eigen::admit_jacobi_eigh(point_dim)
                .map_err(|_| BipopLedgerError::at(0, "eigensolver-admission"))?;
        }
        let admitted_max_ordinal =
            scheduler_max_restart_ordinal(base_lambda, self.root.total_budget)
                .map_err(|_| BipopLedgerError::global("restart-envelope"))?;
        let last_ordinal = u64::try_from(self.records.len() - 1)
            .map_err(|_| BipopLedgerError::global("restart-envelope"))?;
        if last_ordinal > admitted_max_ordinal {
            return Err(BipopLedgerError::global("restart-envelope"));
        }
        let admitted_seed_delta = admitted_max_ordinal
            .checked_mul(BIPOP_RESTART_SEED_STRIDE)
            .ok_or_else(|| BipopLedgerError::global("admission-seed-range"))?;
        base_seed
            .checked_add(admitted_seed_delta)
            .ok_or_else(|| BipopLedgerError::global("admission-seed-range"))?;
        let mut cursor = 0usize;
        let mut large_budget_used = 0usize;
        let mut small_budget_used = 0usize;
        let mut large_runs = 0u32;
        let mut expected_restart_stream = StreamKey {
            seed: self.root.seed,
            kernel: K_CMA,
            tile: 1,
        }
        .stream();

        for (index, record) in self.records.iter().enumerate() {
            // Production stops immediately after publishing rung eight, so no
            // tenth-large or trailing-small record can be authentic even if
            // its local lane arithmetic is otherwise self-consistent.
            if large_runs > BIPOP_LARGE_RUN_CAP {
                return Err(BipopLedgerError::at(index, "large-run-cap"));
            }
            if record.schema_version != BIPOP_RESTART_SCHEMA_VERSION {
                return Err(BipopLedgerError::at(index, "schema-version"));
            }
            let expected_ordinal = u64::try_from(index)
                .map_err(|_| BipopLedgerError::at(index, "ordinal-overflow"))?;
            if record.ordinal != expected_ordinal {
                return Err(BipopLedgerError::at(index, "ordinal"));
            }
            let expected_seed = expected_ordinal
                .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                .and_then(|delta| base_seed.checked_add(delta))
                .ok_or_else(|| BipopLedgerError::at(index, "derived-seed-overflow"))?;
            if record.seed != expected_seed {
                return Err(BipopLedgerError::at(index, "derived-seed"));
            }
            if record.start.len() != point_dim || record.report.x_best.len() != point_dim {
                return Err(BipopLedgerError::at(index, "point-dimension"));
            }
            if record.start.iter().any(|value| !value.is_finite()) {
                return Err(BipopLedgerError::at(index, "finite-start"));
            }
            if index > 0 {
                for (component, actual) in record.start.iter().enumerate() {
                    let expected = self.root.sigma.mul_add(
                        expected_restart_stream.next_normal(),
                        self.root.start[component],
                    );
                    if actual.to_bits() != expected.to_bits() {
                        return Err(BipopLedgerError::at(index, "restart-start"));
                    }
                }
            }
            if record.report.x_best.iter().any(|value| !value.is_finite()) {
                return Err(BipopLedgerError::at(index, "finite-best-point"));
            }
            if !record.report.sigma.is_finite() || record.report.sigma.is_sign_negative() {
                return Err(BipopLedgerError::at(index, "finite-nonnegative-sigma"));
            }

            let expected_lane = if large_budget_used <= small_budget_used {
                BipopLane::Large
            } else {
                BipopLane::Small
            };
            if record.lane != expected_lane {
                return Err(BipopLedgerError::at(index, "lane-selection"));
            }
            let expected_lambda = match expected_lane {
                BipopLane::Large => 1usize
                    .checked_shl(large_runs)
                    .and_then(|scale| base_lambda.checked_mul(scale))
                    .ok_or_else(|| BipopLedgerError::at(index, "population-overflow"))?,
                BipopLane::Small => base_lambda,
            };
            if record.lambda != expected_lambda || self.schedule[index] != record.lambda {
                return Err(BipopLedgerError::at(index, "population-schedule"));
            }
            if record.trace_start != cursor {
                return Err(BipopLedgerError::at(index, "trace-start"));
            }
            let expected_end = cursor
                .checked_add(record.report.evals)
                .ok_or_else(|| BipopLedgerError::at(index, "trace-overflow"))?;
            if record.trace_end != expected_end {
                return Err(BipopLedgerError::at(index, "trace-end"));
            }
            if record.report.evals > record.allocated_budget {
                return Err(BipopLedgerError::at(index, "local-budget"));
            }
            if record.allocated_budget == 0 {
                return Err(BipopLedgerError::at(index, "positive-local-budget"));
            }
            let accounted_evals = record
                .report
                .generations
                .checked_mul(record.lambda)
                .and_then(|samples| samples.checked_add(1))
                .ok_or_else(|| BipopLedgerError::at(index, "evaluation-overflow"))?;
            if record.report.evals != accounted_evals {
                return Err(BipopLedgerError::at(index, "generation-accounting"));
            }
            match record.stop_reason {
                CmaStopReason::TargetReached => {
                    if !record.report.converged || self.root.target.is_none() {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
                CmaStopReason::BudgetExhausted => {
                    let next_generation = record
                        .report
                        .evals
                        .checked_add(record.lambda)
                        .ok_or_else(|| BipopLedgerError::at(index, "evaluation-overflow"))?;
                    if record.report.converged || next_generation <= record.allocated_budget {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
                CmaStopReason::Stagnated => {
                    if record.report.converged || record.report.generations == 0 {
                        return Err(BipopLedgerError::at(index, "terminal-reason"));
                    }
                }
            }
            if index + 1 < self.records.len() && record.report.converged {
                return Err(BipopLedgerError::at(index, "continued-after-convergence"));
            }
            // Validate the report's intrinsic generation/terminal arithmetic
            // first so an overflowing forged report cannot hide behind a second
            // scheduler violation. Then bind the local cap to the production
            // `lambda * 250` envelope used by the admission theorem.
            let local_envelope = record
                .lambda
                .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
                .ok_or_else(|| BipopLedgerError::at(index, "local-budget-overflow"))?;
            let remaining = self
                .root
                .total_budget
                .checked_sub(cursor)
                .ok_or_else(|| BipopLedgerError::at(index, "aggregate-budget"))?;
            let expected_allocated_budget = local_envelope.min(remaining);
            if record.allocated_budget != expected_allocated_budget {
                return Err(BipopLedgerError::at(index, "allocated-budget"));
            }
            cma_stream_block_bound(point_dim, record.lambda, record.report.generations)
                .map_err(|_| BipopLedgerError::at(index, "candidate-counter-range"))?;

            cursor = expected_end;
            match record.lane {
                BipopLane::Large => {
                    large_budget_used = large_budget_used
                        .checked_add(record.report.evals)
                        .ok_or_else(|| BipopLedgerError::at(index, "lane-budget-overflow"))?;
                    large_runs = large_runs
                        .checked_add(1)
                        .ok_or_else(|| BipopLedgerError::at(index, "large-run-overflow"))?;
                }
                BipopLane::Small => {
                    small_budget_used = small_budget_used
                        .checked_add(record.report.evals)
                        .ok_or_else(|| BipopLedgerError::at(index, "lane-budget-overflow"))?;
                }
            }
        }

        if cursor != self.total_evals {
            return Err(BipopLedgerError::global("total-evaluations"));
        }
        let last = self
            .records
            .last()
            .ok_or_else(|| BipopLedgerError::global("nonempty"))?;
        if cursor < self.root.total_budget
            && !last.report.converged
            && large_runs <= BIPOP_LARGE_RUN_CAP
        {
            return Err(BipopLedgerError::global("nonterminal-prefix"));
        }
        let mut expected_best = 0usize;
        for index in 1..self.records.len() {
            if self.records[index]
                .report
                .f_best
                .total_cmp(&self.records[expected_best].report.f_best)
                .is_lt()
            {
                expected_best = index;
            }
        }
        if self.best_restart != expected_best {
            return Err(BipopLedgerError::global("best-restart"));
        }
        if !cma_reports_match_bits(&self.best, &self.records[expected_best].report) {
            return Err(BipopLedgerError::global("best-projection"));
        }
        if self.trace_rows.len() != self.total_evals {
            return Err(BipopLedgerError::global("trace-length"));
        }
        let expected_point_entries = self
            .total_evals
            .checked_mul(point_dim)
            .ok_or_else(|| BipopLedgerError::global("trace-point-overflow"))?;
        if self.trace_points.len() != expected_point_entries {
            return Err(BipopLedgerError::global("trace-point-layout"));
        }
        for (restart_index, record) in self.records.iter().enumerate() {
            let rows = self
                .trace_rows
                .get(record.trace_start..record.trace_end)
                .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-row-range"))?;
            let target_witnessed = self
                .root
                .target
                .is_some_and(|target| rows.iter().any(|row| row.objective <= target));
            let target_terminal = record.stop_reason == CmaStopReason::TargetReached;
            if target_witnessed != target_terminal || record.report.converged != target_terminal {
                return Err(BipopLedgerError::at(restart_index, "terminal-reason"));
            }
            let mut local_best = 0usize;
            for (local_offset, row) in rows.iter().enumerate() {
                let global_offset = record
                    .trace_start
                    .checked_add(local_offset)
                    .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-overflow"))?;
                if row.schema_version != BIPOP_EVALUATION_SCHEMA_VERSION {
                    return Err(BipopLedgerError::at(restart_index, "trace-schema-version"));
                }
                if row.restart != record.ordinal {
                    return Err(BipopLedgerError::at(restart_index, "trace-restart"));
                }
                if row.local_offset != local_offset {
                    return Err(BipopLedgerError::at(restart_index, "trace-local-offset"));
                }
                let point_start = global_offset
                    .checked_mul(point_dim)
                    .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-overflow"))?;
                let point_end = point_start
                    .checked_add(point_dim)
                    .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-overflow"))?;
                let point = self
                    .trace_points
                    .get(point_start..point_end)
                    .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-range"))?;
                if point.iter().any(|coordinate| !coordinate.is_finite()) {
                    return Err(BipopLedgerError::at(restart_index, "finite-trace-point"));
                }
                if local_offset == 0 && !f64_slices_match_bits(point, &record.start) {
                    return Err(BipopLedgerError::at(
                        restart_index,
                        "trace-start-projection",
                    ));
                }
                if row.objective.total_cmp(&rows[local_best].objective).is_lt() {
                    local_best = local_offset;
                }
            }
            let best_global = record
                .trace_start
                .checked_add(local_best)
                .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-overflow"))?;
            let best_point_start = best_global
                .checked_mul(point_dim)
                .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-overflow"))?;
            let best_point_end = best_point_start
                .checked_add(point_dim)
                .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-overflow"))?;
            let best_point = self
                .trace_points
                .get(best_point_start..best_point_end)
                .ok_or_else(|| BipopLedgerError::at(restart_index, "trace-point-range"))?;
            if rows[local_best].objective.to_bits() != record.report.f_best.to_bits()
                || !f64_slices_match_bits(best_point, &record.report.x_best)
            {
                return Err(BipopLedgerError::at(restart_index, "trace-best-projection"));
            }
        }
        let expected_trace_identity =
            build_bipop_trace_identity(&self.root, &self.trace_rows, &self.trace_points)
                .map_err(|_| BipopLedgerError::global("trace-identity-encoding"))?;
        if self.trace_identity != expected_trace_identity {
            return Err(BipopLedgerError::global("trace-identity"));
        }
        let expected_study_identity = build_bipop_study_identity(
            self.schema_version,
            &self.root,
            &expected_root.identity,
            &self.schedule,
            self.total_evals,
            &self.records,
            self.best_restart,
            &self.best,
            self.trace_identity,
            expected_trace_identity,
        )
        .map_err(|_| BipopLedgerError::global("study-identity-encoding"))?;
        if self.study_identity != expected_study_identity {
            return Err(BipopLedgerError::global("study-identity"));
        }
        Ok(())
    }
}

/// Fallible BIPOP-CMA-ES under one admitted hard callback budget.
///
/// `target = None` disables target stopping. Every raw-input, conservative
/// seed-range, dimension, population, and local-budget formula is checked
/// before the first callback. Per-restart start generation and accounting are
/// checked again before that restart becomes visible to the callback.
/// A finite target is reached when any callback is numerically at or below it;
/// that witness is independent of the exact total-order representative retained
/// as [`CmaReport::f_best`].
///
/// # Errors
/// Returns [`BipopError`] without invoking `f` for raw-input or preflight
/// refusal. A later generated-start or generated-candidate refusal can follow
/// completed callbacks or earlier restarts, but the affected decision point is
/// never passed to the objective or published.
#[allow(clippy::too_many_lines)] // scheduler and record publication are one atomic state machine
pub fn try_bipop_cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    target: Option<f64>,
    seed: u64,
) -> Result<BipopReport, BipopError> {
    let admission = admit_bipop(x0, sigma0, total_budget, target, seed)?;
    run_admitted_bipop(f, x0, sigma0, target, seed, admission)
}

/// Legacy BIPOP spelling retained as a checked compatibility projection.
///
/// Finite targets map to `Some(target)`. Historical `-∞` means no target and
/// maps to `None`; all other malformed inputs now refuse before callbacks and
/// panic at this legacy boundary. New callers should use [`try_bipop_cmaes`]
/// for typed refusal.
#[must_use]
pub fn bipop_cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    f_target: f64,
    seed: u64,
) -> BipopReport {
    let target = if f_target.to_bits() == f64::NEG_INFINITY.to_bits() {
        None
    } else {
        Some(f_target)
    };
    try_bipop_cmaes(f, x0, sigma0, total_budget, target, seed)
        .unwrap_or_else(|error| panic!("BIPOP input or scheduler refused: {error}"))
}

#[allow(clippy::too_many_lines)] // scheduler and record publication are one atomic state machine
fn run_admitted_bipop<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    target: Option<f64>,
    seed: u64,
    admission: BipopAdmission,
) -> Result<BipopReport, BipopError> {
    let base_lambda = admission.base_lambda;
    let total_budget = admission.total_budget;
    let dimension = admission.dimension;
    let root = build_bipop_root_inputs(x0, sigma0, total_budget, target, seed)?;
    let trace_point_capacity =
        total_budget
            .checked_mul(dimension)
            .ok_or(BipopError::DimensionOverflow {
                what: "callback trace coordinate entries",
            })?;
    let trace_point_bytes = trace_point_capacity
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(BipopError::DimensionOverflow {
            what: "callback trace coordinate bytes",
        })?;
    let trace_row_bytes = total_budget
        .checked_mul(core::mem::size_of::<BipopEvaluationRow>())
        .ok_or(BipopError::DimensionOverflow {
            what: "callback trace row bytes",
        })?;
    if trace_point_bytes > isize::MAX as usize || trace_row_bytes > isize::MAX as usize {
        return Err(BipopError::DimensionOverflow {
            what: "callback trace address space",
        });
    }
    let mut trace_rows: Vec<BipopEvaluationRow> = Vec::new();
    let mut trace_points: Vec<f64> = Vec::new();
    trace_rows
        .try_reserve_exact(total_budget)
        .map_err(|_| BipopError::TraceAllocationFailed {
            evaluations: total_budget,
            point_entries: trace_point_capacity,
        })?;
    trace_points
        .try_reserve_exact(trace_point_capacity)
        .map_err(|_| BipopError::TraceAllocationFailed {
            evaluations: total_budget,
            point_entries: trace_point_capacity,
        })?;
    let mut records: Vec<BipopRestartRecord> = Vec::new();
    let mut total_evals = 0usize;
    let mut best_restart: Option<usize> = None;
    let mut large_runs = 0u32;
    let mut restart = 0u64;
    let mut small_budget_used = 0usize;
    let mut large_budget_used = 0usize;
    let mut restart_stream_blocks_used = 0u128;
    // Deterministic restart-start perturbations (a restart from the SAME
    // point with a tiny sigma is just a polish run and cannot escape a
    // local basin — measured during bring-up on rastrigin).
    let mut restart_stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 1,
    }
    .stream();
    while total_evals < total_budget {
        if restart > admission.max_restart_ordinal {
            return Err(BipopError::InternalInvariant {
                what: "restart ordinal exceeded its admitted scheduler bound",
            });
        }
        // BIPOP rule: run LARGE next if its cumulative budget lags.
        let run_large = large_budget_used <= small_budget_used;
        let lambda = if run_large {
            let scale = 1usize
                .checked_shl(large_runs)
                .ok_or(BipopError::PopulationOverflow {
                    large_run: large_runs,
                })?;
            base_lambda
                .checked_mul(scale)
                .ok_or(BipopError::PopulationOverflow {
                    large_run: large_runs,
                })?
        } else {
            base_lambda
        };
        // Per-run budget scales with the population (≈250 generations):
        // handing a small-λ run half the TOTAL budget just polishes one
        // local minimum expensively — the doubling ladder must be reached
        // (measured during bring-up on rastrigin).
        let local_envelope = lambda
            .checked_mul(BIPOP_GENERATIONS_PER_RESTART)
            .ok_or(BipopError::LocalBudgetOverflow { lambda })?;
        if lambda > admission.max_large_lambda || local_envelope > admission.max_local_budget {
            return Err(BipopError::InternalInvariant {
                what: "runtime population or local budget exceeded admission",
            });
        }
        let remaining =
            total_budget
                .checked_sub(total_evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "remaining aggregate budget",
                })?;
        let budget = local_envelope.min(remaining);
        if budget == 0 {
            return Err(BipopError::InternalInvariant {
                what: "a launched restart must receive at least one callback",
            });
        }
        let max_generations = (budget - 1) / lambda;
        let planned_cma_stream_blocks = cma_stream_block_bound(dimension, lambda, max_generations)?;
        if planned_cma_stream_blocks > admission.max_cma_stream_blocks {
            return Err(BipopError::InternalInvariant {
                what: "runtime CMA stream exceeded its admitted counter bound",
            });
        }
        // Reserve every scheduler addition against the full local cap before
        // this restart becomes visible to the callback. The post-run checks
        // below then narrow these envelopes to the actual callback count.
        let trace_start = total_evals;
        let trace_cap = trace_start
            .checked_add(budget)
            .ok_or(BipopError::ArithmeticOverflow {
                restart,
                what: "aggregate trace envelope",
            })?;
        if trace_cap > total_budget {
            return Err(BipopError::InternalAggregateBudgetViolation {
                restart,
                total_spent: trace_cap,
                total_budget,
            });
        }
        let (large_budget_cap, small_budget_cap) = if run_large {
            (
                large_budget_used
                    .checked_add(budget)
                    .ok_or(BipopError::ArithmeticOverflow {
                        restart,
                        what: "large-lane budget envelope",
                    })?,
                small_budget_used,
            )
        } else {
            (
                large_budget_used,
                small_budget_used
                    .checked_add(budget)
                    .ok_or(BipopError::ArithmeticOverflow {
                        restart,
                        what: "small-lane budget envelope",
                    })?,
            )
        };
        let next_restart = restart
            .checked_add(1)
            .ok_or(BipopError::ArithmeticOverflow {
                restart,
                what: "restart ordinal",
            })?;
        let large_runs_after = if run_large {
            large_runs
                .checked_add(1)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "large-run count",
                })?
        } else {
            large_runs
        };
        let params = CmaParams {
            lambda,
            sigma0,
            max_evals: budget,
            f_target: target.unwrap_or(f64::NEG_INFINITY),
            eigen_interval: 1,
        };
        // Restarts after the first launch from a perturbed start point.
        let start: Vec<f64> = if restart == 0 {
            x0.to_vec()
        } else {
            let blocks_for_start = restart_stream_block_bound(dimension, 1)?;
            let next_restart_stream_blocks = restart_stream_blocks_used
                .checked_add(blocks_for_start)
                .ok_or(BipopError::RandomCounterRangeOverflow {
                    stream: "restart-perturbation",
                    required_blocks: None,
                })?;
            if next_restart_stream_blocks > admission.max_restart_stream_blocks
                || next_restart_stream_blocks > FS_RAND_STREAM_COUNTER_CARDINALITY
            {
                return Err(BipopError::InternalInvariant {
                    what: "runtime restart stream exceeded its admitted counter bound",
                });
            }
            let start = x0
                .iter()
                .map(|&v| sigma0.mul_add(restart_stream.next_normal(), v))
                .collect();
            restart_stream_blocks_used = next_restart_stream_blocks;
            start
        };
        for (component, value) in start.iter().enumerate() {
            if !value.is_finite() {
                return Err(BipopError::GeneratedStartNonFinite {
                    restart,
                    component,
                    bits: value.to_bits(),
                });
            }
        }
        let derived_seed = restart
            .checked_mul(BIPOP_RESTART_SEED_STRIDE)
            .and_then(|delta| seed.checked_add(delta))
            .ok_or(BipopError::SeedRangeOverflow {
                seed,
                max_restart_ordinal: restart,
            })?;
        if trace_rows.len() != trace_start
            || trace_points.len()
                != trace_start
                    .checked_mul(dimension)
                    .ok_or(BipopError::ArithmeticOverflow {
                        restart,
                        what: "callback trace point offset",
                    })?
        {
            return Err(BipopError::InternalInvariant {
                what: "retained callback trace must precede its restart record",
            });
        }
        let trace_row_start = trace_rows.len();
        let trace_point_start = trace_points.len();
        let result = {
            let mut traced_objective = |point: &[f64]| {
                let local_offset = trace_rows.len() - trace_row_start;
                let objective = f(point);
                trace_points.extend_from_slice(point);
                trace_rows.push(BipopEvaluationRow {
                    schema_version: BIPOP_EVALUATION_SCHEMA_VERSION,
                    restart,
                    local_offset,
                    objective,
                });
                objective
            };
            cmaes_with_stop_target(&mut traced_objective, &start, &params, derived_seed, target)
        };
        let (rep, stop_reason) =
            result.map_err(|error| BipopError::GeneratedCandidateNonFinite {
                restart,
                generation: error.generation,
                candidate: error.candidate,
                component: error.component,
                bits: error.bits,
            })?;
        let retained_rows = trace_rows.len() - trace_row_start;
        let retained_points = trace_points.len() - trace_point_start;
        let expected_retained_points =
            rep.evals
                .checked_mul(dimension)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "callback trace point count",
                })?;
        if retained_rows != rep.evals || retained_points != expected_retained_points {
            return Err(BipopError::InternalInvariant {
                what: "CMA report must project the retained callback trace",
            });
        }
        let actual_cma_stream_blocks = cma_stream_block_bound(dimension, lambda, rep.generations)?;
        if actual_cma_stream_blocks > planned_cma_stream_blocks {
            return Err(BipopError::InternalInvariant {
                what: "actual CMA stream consumption exceeded its preflight",
            });
        }
        if rep.evals > budget {
            return Err(BipopError::InternalBudgetViolation {
                restart,
                spent: rep.evals,
                allocated: budget,
            });
        }
        let trace_end =
            trace_start
                .checked_add(rep.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "aggregate trace end",
                })?;
        if trace_end > total_budget {
            return Err(BipopError::InternalAggregateBudgetViolation {
                restart,
                total_spent: trace_end,
                total_budget,
            });
        }
        if trace_end > trace_cap {
            return Err(BipopError::InternalInvariant {
                what: "actual trace exceeded its admitted envelope",
            });
        }
        let record_index = records.len();
        let is_better = best_restart.is_none_or(|best_index| {
            rep.f_best
                .total_cmp(&records[best_index].report.f_best)
                .is_lt()
        });
        records.push(BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: restart,
            lane: if run_large {
                BipopLane::Large
            } else {
                BipopLane::Small
            },
            lambda,
            allocated_budget: budget,
            seed: derived_seed,
            start,
            trace_start,
            trace_end,
            stop_reason,
            report: rep,
        });
        if is_better {
            best_restart = Some(record_index);
        }
        total_evals = trace_end;
        if run_large {
            large_budget_used = large_budget_used
                .checked_add(records[record_index].report.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "large-lane evaluation total",
                })?;
            if large_budget_used > large_budget_cap {
                return Err(BipopError::InternalInvariant {
                    what: "actual large-lane total exceeded its admitted envelope",
                });
            }
        } else {
            small_budget_used = small_budget_used
                .checked_add(records[record_index].report.evals)
                .ok_or(BipopError::ArithmeticOverflow {
                    restart,
                    what: "small-lane evaluation total",
                })?;
            if small_budget_used > small_budget_cap {
                return Err(BipopError::InternalInvariant {
                    what: "actual small-lane total exceeded its admitted envelope",
                });
            }
        }
        large_runs = large_runs_after;
        if records[record_index].report.converged {
            break;
        }
        restart = next_restart;
        if large_runs > BIPOP_LARGE_RUN_CAP {
            // Cap the LADDER, not total restarts: small runs are cheap
            // and interleave freely; counting them against the cap
            // stalled the ladder at λ ≈ 64 (measured during bring-up).
            break;
        }
    }
    let best_restart = best_restart.ok_or(BipopError::InternalInvariant {
        what: "positive admitted budget must launch one restart",
    })?;
    let schedule: Vec<usize> = records.iter().map(BipopRestartRecord::lambda).collect();
    let best = records[best_restart].report.clone();
    let trace_identity = build_bipop_trace_identity(&root, &trace_rows, &trace_points)?;
    let study_identity = build_bipop_study_identity(
        BIPOP_REPORT_SCHEMA_VERSION,
        &root,
        &root.identity,
        &schedule,
        total_evals,
        &records,
        best_restart,
        &best,
        trace_identity,
        trace_identity,
    )?;
    let report = BipopReport {
        best,
        schedule,
        total_evals,
        schema_version: BIPOP_REPORT_SCHEMA_VERSION,
        root,
        records,
        best_restart,
        trace_rows,
        trace_points,
        trace_identity,
        study_identity,
    };
    report
        .validate_ledger()
        .map_err(|error| BipopError::GeneratedLedgerInvalid {
            restart: error.restart,
            invariant: error.invariant,
        })?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_root(
        start: &[f64],
        sigma: f64,
        total_budget: usize,
        target: Option<f64>,
        seed: u64,
    ) -> BipopRootInputs {
        build_bipop_root_inputs(start, sigma, total_budget, target, seed)
            .expect("test root must be identity-representable")
    }

    fn report_without_trace(
        best: CmaReport,
        schedule: Vec<usize>,
        total_evals: usize,
        root: BipopRootInputs,
        records: Vec<BipopRestartRecord>,
        best_restart: usize,
    ) -> BipopReport {
        let trace_rows = Vec::new();
        let trace_points = Vec::new();
        let trace_identity = build_bipop_trace_identity(&root, &trace_rows, &trace_points)
            .expect("empty test trace identity must be representable");
        let study_identity = build_bipop_study_identity(
            BIPOP_REPORT_SCHEMA_VERSION,
            &root,
            &root.identity,
            &schedule,
            total_evals,
            &records,
            best_restart,
            &best,
            trace_identity,
            trace_identity,
        )
        .expect("test study identity must be representable");
        BipopReport {
            best,
            schedule,
            total_evals,
            schema_version: BIPOP_REPORT_SCHEMA_VERSION,
            root,
            records,
            best_restart,
            trace_rows,
            trace_points,
            trace_identity,
            study_identity,
        }
    }

    fn reseal_private_identities(report: &mut BipopReport) {
        report.root.identity = build_bipop_root_inputs(
            &report.root.start,
            report.root.sigma,
            report.root.total_budget,
            report.root.target,
            report.root.seed,
        )
        .expect("mutated test root remains representable")
        .identity;
        report.trace_identity =
            build_bipop_trace_identity(&report.root, &report.trace_rows, &report.trace_points)
                .expect("mutated test trace remains representable");
        report.study_identity = report
            .computed_study_identity()
            .expect("mutated test study remains representable");
    }

    #[test]
    fn bipop_root_identity_fields_move_independently() {
        let start = [1.0, -2.0];
        let source = BipopRootIdentitySource {
            start: &start,
            sigma: 0.5,
            total_budget: 20,
            target: Some(-1.0),
            seed: 7,
        };
        let base = bipop_root_identity(&source).expect("canonical root identity");

        let changed_start = [1.0, -2.5];
        let mut changed = source;
        changed.start = &changed_start;
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        let changed_dimension = [1.0, -2.0, 3.0];
        changed = source;
        changed.start = &changed_dimension;
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        changed = source;
        changed.sigma = f64::from_bits(source.sigma.to_bits() ^ 1);
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        changed = source;
        changed.total_budget += 1;
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        changed = source;
        changed.target = None;
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        changed = source;
        changed.target = Some(f64::from_bits((-1.0_f64).to_bits() ^ 1));
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);
        changed = source;
        changed.seed ^= 1;
        assert_ne!(bipop_root_identity(&changed).unwrap(), base);

        let reordered = IdentityBuilder::new(BIPOP_ROOT_IDENTITY_KIND)
            .u64(
                "admission-schema-version",
                u64::from(BIPOP_ADMISSION_SCHEMA_VERSION),
            )
            .u64(
                "report-schema-version",
                u64::from(BIPOP_REPORT_SCHEMA_VERSION),
            )
            .u64("dimension", 2)
            .u64("total-budget", 20)
            .u64("root-seed", 7)
            .f64_bits("initial-sigma", 0.5)
            .flag("target-present", true)
            .f64_bits("target", -1.0)
            .f64_bits("start-coordinate", 1.0)
            .f64_bits("start-coordinate", -2.0)
            .finish();
        assert_ne!(reordered, base, "canonical field order is semantic");
    }

    #[test]
    fn bipop_trace_identity_fields_move_independently() {
        let root = test_root(&[1.0, -2.0], 0.5, 20, None, 7);
        let rows = [
            BipopEvaluationRow {
                schema_version: BIPOP_EVALUATION_SCHEMA_VERSION,
                restart: 0,
                local_offset: 0,
                objective: 5.0,
            },
            BipopEvaluationRow {
                schema_version: BIPOP_EVALUATION_SCHEMA_VERSION,
                restart: 0,
                local_offset: 1,
                objective: 4.0,
            },
        ];
        let points = [1.0, -2.0, 0.5, -1.5];
        let source = BipopTraceIdentitySource {
            root_canonical_bytes: root.identity.canonical_bytes(),
            dimension: 2,
            rows: &rows,
            points: &points,
        };
        let base = bipop_trace_identity(&source).expect("canonical trace identity");

        let mut changed = source;
        changed.root_canonical_bytes = b"different-root";
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        let one_dimensional_points = [1.0, 0.5];
        changed = source;
        changed.dimension = 1;
        changed.points = &one_dimensional_points;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed = source;
        changed.rows = &rows[..1];
        changed.points = &points[..2];
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        let mut changed_rows = rows.clone();
        changed_rows[0].schema_version ^= 1;
        changed = source;
        changed.rows = &changed_rows;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed_rows = rows.clone();
        changed_rows[0].restart = 1;
        changed.rows = &changed_rows;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed_rows = rows.clone();
        changed_rows[1].local_offset = 2;
        changed.rows = &changed_rows;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed_rows = rows.clone();
        changed_rows[0].objective = f64::from_bits(rows[0].objective.to_bits() ^ 1);
        changed.rows = &changed_rows;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed_rows = rows.clone();
        changed_rows.swap(0, 1);
        changed.rows = &changed_rows;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        let mut changed_points = points;
        changed_points[2] = f64::from_bits(changed_points[2].to_bits() ^ 1);
        changed = source;
        changed.points = &changed_points;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
        changed_points = points;
        changed_points.swap(0, 2);
        changed = source;
        changed.points = &changed_points;
        assert_ne!(bipop_trace_identity(&changed).unwrap(), base);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One source-level test moves every complete-study projection class.
    fn bipop_study_identity_fields_move_independently() {
        let mut objective = |point: &[f64]| point.iter().map(|value| value * value).sum::<f64>();
        let report = try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, 11)
            .expect("study owner fixture admits");
        assert!(
            report.records.len() > 1,
            "fixture needs two restart records"
        );
        let root_content = build_bipop_root_inputs(
            &report.root.start,
            report.root.sigma,
            report.root.total_budget,
            report.root.target,
            report.root.seed,
        )
        .expect("recompute root identity");
        let callback_content =
            build_bipop_trace_identity(&report.root, &report.trace_rows, &report.trace_points)
                .expect("recompute callback identity");
        let source = BipopStudyIdentitySource {
            report_schema_version: report.schema_version,
            root: &report.root,
            root_content_identity: &root_content.identity,
            schedule: &report.schedule,
            total_evals: report.total_evals,
            records: &report.records,
            best_restart: report.best_restart,
            best: &report.best,
            retained_trace_identity: report.trace_identity,
            callback_content_identity: callback_content,
        };
        let base = bipop_study_identity(&source).expect("canonical study identity");
        assert_eq!(base, report.study_identity);

        let mut changed = source;
        changed.report_schema_version ^= 1;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let alternative_root = test_root(&[2.0, -1.0], 0.75, 20, None, 12);
        changed = source;
        changed.root = &alternative_root;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        changed = source;
        changed.root_content_identity = &alternative_root.identity;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let mut schedule = report.schedule.clone();
        schedule[0] += 1;
        changed = source;
        changed.schedule = &schedule;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        changed = source;
        changed.total_evals += 1;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let mut records = report.records.clone();
        records[0].seed ^= 1;
        changed = source;
        changed.records = &records;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        changed = source;
        changed.best_restart = (report.best_restart + 1) % report.records.len();
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let mut best = report.best.clone();
        best.sigma = f64::from_bits(best.sigma.to_bits() ^ 1);
        changed = source;
        changed.best = &best;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let mut retained_trace = report.trace_identity;
        retained_trace.digest[0] ^= 1;
        changed = source;
        changed.retained_trace_identity = retained_trace;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
        let mut recomputed_trace = callback_content;
        recomputed_trace.digest[0] ^= 1;
        changed = source;
        changed.callback_content_identity = recomputed_trace;
        assert_ne!(bipop_study_identity(&changed).unwrap(), base);
    }

    #[test]
    fn bipop_identity_versions_and_domains_fail_closed() {
        assert_eq!(BIPOP_ROOT_IDENTITY_SCHEMA_VERSION, 1);
        assert!(BIPOP_ROOT_IDENTITY_KIND.ends_with("-v1"));
        assert!(BIPOP_TRACE_IDENTITY_DOMAIN.ends_with(".v2"));
        assert!(BIPOP_STUDY_IDENTITY_DOMAIN.ends_with(".v1"));

        let make_report = || {
            let mut objective =
                |point: &[f64]| point.iter().map(|value| value * value).sum::<f64>();
            try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, 11)
                .expect("version refusal fixture admits")
        };
        let mut report = make_report();
        report.schema_version ^= 1;
        assert_eq!(
            report.validate_ledger().unwrap_err().invariant(),
            "report-schema-version"
        );

        report = make_report();
        report.root.identity = IdentityBuilder::new("fs-dfo-bipop-root-v2").finish();
        assert_eq!(
            report.validate_ledger().unwrap_err().invariant(),
            "root-identity"
        );

        report = make_report();
        report.trace_rows[0].schema_version ^= 1;
        assert_eq!(
            report.validate_ledger().unwrap_err().invariant(),
            "trace-schema-version"
        );

        report = make_report();
        report.trace_identity.schema_version ^= 1;
        assert_eq!(
            report.validate_ledger().unwrap_err().invariant(),
            "trace-identity"
        );

        report = make_report();
        report.study_identity.schema_version ^= 1;
        assert!(matches!(
            report.validate_study_identity(),
            Err(BipopStudyAdmissionError::PayloadIdentityMismatch { .. })
        ));
    }

    /// G0: both random domains count every semantic axis, accept exactly the
    /// full `2^64` Philox coordinate set once, and refuse the first reuse. The
    /// small exact KATs are mutation guards for the Box-Muller factor, decision
    /// dimension, population/generation product, and restart ordinal.
    #[test]
    fn random_counter_bounds_count_every_axis_and_refuse_reuse() {
        assert_eq!(
            restart_stream_block_bound(3, 5).expect("small restart bound"),
            2 * 3 * 5
        );
        assert_eq!(
            cma_stream_block_bound(3, 7, 11).expect("small CMA bound"),
            2 * 3 * 7 * 11
        );
        assert_eq!(
            checked_random_counter_blocks(
                "exact-cardinality",
                &[2, FS_RAND_STREAM_COUNTER_CARDINALITY / 2],
            )
            .expect("every counter coordinate may be consumed once"),
            FS_RAND_STREAM_COUNTER_CARDINALITY
        );

        let required_blocks = FS_RAND_STREAM_COUNTER_CARDINALITY + 2;
        assert_eq!(
            checked_random_counter_blocks(
                "first-reuse",
                &[2, FS_RAND_STREAM_COUNTER_CARDINALITY / 2 + 1],
            )
            .expect_err("the first coordinate reuse must refuse"),
            BipopError::RandomCounterRangeOverflow {
                stream: "first-reuse",
                required_blocks: Some(required_blocks),
            }
        );
        assert_eq!(
            checked_random_counter_blocks("u128-overflow", &[u128::MAX, 2])
                .expect_err("the proof accumulator must also fail closed"),
            BipopError::RandomCounterRangeOverflow {
                stream: "u128-overflow",
                required_blocks: None,
            }
        );
    }

    /// G0: checked multiplication alone is not sufficient for a Rust `Vec`.
    /// One dense covariance allocation must also fit the target's `isize`
    /// address-difference domain before any objective callback is possible.
    #[test]
    fn dense_covariance_allocation_refuses_first_unaddressable_size() {
        let element_bytes = core::mem::size_of::<f64>();
        let max_entries = (isize::MAX as usize) / element_bytes;
        let exact_bytes = checked_dense_matrix_allocation(max_entries)
            .expect("the last whole f64 allocation below isize::MAX admits");
        assert!(exact_bytes <= isize::MAX as usize);

        let first_unaddressable = max_entries
            .checked_add(1)
            .expect("one entry beyond the f64 address-space boundary fits usize");
        assert_eq!(
            checked_dense_matrix_allocation(first_unaddressable),
            Err(BipopError::DimensionOverflow {
                what: "dense covariance address space",
            })
        );
        assert_eq!(
            checked_dense_matrix_allocation(usize::MAX),
            Err(BipopError::DimensionOverflow {
                what: "dense covariance bytes",
            }),
            "the checked byte-product branch is independently retained"
        );
    }

    /// G0/G4: the infallible direct-CMA compatibility surface remains a panic
    /// projection, but it shares the finite-query guard with fallible BIPOP and
    /// never exposes the overflowing candidate to user code. This retained seed
    /// has a positive first normal, so `MAX + MAX*z` deterministically overflows
    /// at generation one, candidate zero, component zero.
    #[test]
    fn legacy_cma_panics_before_nonfinite_candidate_callback() {
        let calls = std::cell::Cell::new(0usize);
        let mut objective = |point: &[f64]| {
            calls.set(calls.get() + 1);
            assert!(point.iter().all(|value| value.is_finite()));
            0.0
        };
        let params = CmaParams::standard(1, f64::MAX, 5, f64::NEG_INFINITY);
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            cmaes(&mut objective, &[f64::MAX], &params, 0xB1_90_00_02)
        }));

        assert!(
            panic.is_err(),
            "legacy CMA must project typed refusal as panic"
        );
        assert_eq!(
            calls.get(),
            1,
            "only the finite initial point may reach the callback"
        );
    }

    /// G0/G3: production stops immediately after its ninth large run (rungs
    /// zero through eight). A locally self-consistent trailing small record is
    /// therefore a forged history and must fail the independent validator.
    #[test]
    fn ledger_refuses_any_record_after_ninth_large_run() {
        let base_lambda = 4usize;
        let total_budget = 258_056usize;
        let mut records = Vec::new();
        let mut schedule = Vec::new();
        let mut cursor = 0usize;
        let mut restart_stream = StreamKey {
            seed: 0,
            kernel: K_CMA,
            tile: 1,
        }
        .stream();
        let mut push_record =
            |lane: BipopLane, lambda: usize, allocated_budget: usize, generations: usize| {
                let index = records.len();
                let ordinal = u64::try_from(index).expect("fixture ordinal fits");
                let evals = generations
                    .checked_mul(lambda)
                    .and_then(|samples| samples.checked_add(1))
                    .expect("fixture evaluation count fits");
                let trace_end = cursor
                    .checked_add(evals)
                    .expect("fixture trace offset fits");
                let seed = ordinal
                    .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                    .expect("fixture seed fits");
                let report = CmaReport {
                    x_best: vec![0.0],
                    f_best: 0.0,
                    evals,
                    generations,
                    converged: false,
                    sigma: 1.0,
                };
                records.push(BipopRestartRecord {
                    schema_version: BIPOP_RESTART_SCHEMA_VERSION,
                    ordinal,
                    lane,
                    lambda,
                    allocated_budget,
                    seed,
                    start: if index == 0 {
                        vec![0.0]
                    } else {
                        vec![restart_stream.next_normal()]
                    },
                    trace_start: cursor,
                    trace_end,
                    stop_reason: CmaStopReason::Stagnated,
                    report,
                });
                schedule.push(lambda);
                cursor = trace_end;
            };

        // Each first-eight large rung is paired with enough small-lane
        // generations to make the cumulative spends exactly equal. Every
        // allocated cap is the production local envelope under this large
        // retained aggregate budget.
        for large_run in 0..BIPOP_LARGE_RUN_CAP {
            let scale = 1usize << large_run;
            let large_lambda = base_lambda * scale;
            push_record(
                BipopLane::Large,
                large_lambda,
                large_lambda * BIPOP_GENERATIONS_PER_RESTART,
                1,
            );
            push_record(
                BipopLane::Small,
                base_lambda,
                base_lambda * BIPOP_GENERATIONS_PER_RESTART,
                scale,
            );
        }

        // Rung eight consumes exactly the remaining aggregate allocation and
        // is the ninth/terminal large run. The following locally plausible
        // small record is forbidden before any of its fields can self-justify.
        let final_large_lambda = base_lambda * (1usize << BIPOP_LARGE_RUN_CAP);
        push_record(
            BipopLane::Large,
            final_large_lambda,
            final_large_lambda * BIPOP_GENERATIONS_PER_RESTART,
            1,
        );
        push_record(
            BipopLane::Small,
            base_lambda,
            base_lambda * BIPOP_GENERATIONS_PER_RESTART,
            1,
        );
        assert_eq!(records.len(), 18);
        assert_eq!(cursor, 3_086);
        let report = report_without_trace(
            records[0].report.clone(),
            schedule,
            cursor,
            test_root(&[0.0], 1.0, total_budget, None, 0),
            records,
            0,
        );

        let error = report
            .validate_ledger()
            .expect_err("post-cap history must fail closed");
        assert_eq!(error.restart(), Some(17));
        assert_eq!(error.invariant(), "large-run-cap");
    }

    /// G0/G3: a forged record cannot mint more large-lane spend than the
    /// production `lambda * 250` envelope. Without this independent check, an
    /// oversized first large record could justify arbitrarily many small
    /// records and invalidate the scheduler theorem used by admission.
    #[test]
    fn ledger_refuses_allocated_budget_above_local_envelope() {
        let lambda = 4usize;
        let allocated_budget = lambda * BIPOP_GENERATIONS_PER_RESTART + 1;
        let generations = (allocated_budget - 1) / lambda;
        let run = CmaReport {
            x_best: vec![0.0],
            f_best: 0.0,
            evals: allocated_budget,
            generations,
            converged: false,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget,
            seed: 0,
            start: vec![0.0],
            trace_start: 0,
            trace_end: allocated_budget,
            stop_reason: CmaStopReason::BudgetExhausted,
            report: run.clone(),
        };
        let report = report_without_trace(
            run.clone(),
            vec![lambda],
            allocated_budget,
            test_root(&[0.0], 1.0, allocated_budget, None, 0),
            vec![record],
            0,
        );

        let error = report
            .validate_ledger()
            .expect_err("an oversized local cap is not production-authentic");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "allocated-budget");
    }

    /// G0/G3: a prefix can be locally self-consistent yet impossible as a
    /// terminal report under its retained aggregate budget. The validator must
    /// not accept truncation merely because every surviving interval closes.
    #[test]
    fn ledger_refuses_nonterminal_prefix_under_retained_budget() {
        let mut objective = |_point: &[f64]| 0.0;
        let mut report = try_bipop_cmaes(&mut objective, &[0.0], 0.5, 20, None, 0)
            .expect("production fixture admits");
        assert!(report.records.len() > 1, "fixture must contain a suffix");
        report.records.pop();
        report.schedule.pop();
        report.total_evals = report
            .records
            .last()
            .expect("the truncated prefix remains nonempty")
            .trace_end;

        let error = report
            .validate_ledger()
            .expect_err("a nonterminal prefix is not a complete run");
        assert_eq!(error.restart(), None);
        assert_eq!(error.invariant(), "nonterminal-prefix");
    }

    /// G0/G3: target convergence is an immediate scheduler terminal. A later
    /// record cannot be authenticated by giving that record another locally
    /// plausible converged report and relying on the final-record check.
    #[test]
    fn ledger_refuses_any_record_after_target_convergence() {
        let lambda = 4usize;
        let total_budget = 1_001usize;
        let run = CmaReport {
            x_best: vec![0.0],
            f_best: 0.0,
            evals: 1,
            generations: 0,
            converged: true,
            sigma: 1.0,
        };
        let make_record = |index: usize, lane: BipopLane| BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: u64::try_from(index).expect("fixture ordinal fits"),
            lane,
            lambda,
            allocated_budget: 1_000,
            seed: u64::try_from(index)
                .expect("fixture ordinal fits")
                .checked_mul(BIPOP_RESTART_SEED_STRIDE)
                .expect("fixture seed fits"),
            start: vec![0.0],
            trace_start: index,
            trace_end: index + 1,
            stop_reason: CmaStopReason::TargetReached,
            report: run.clone(),
        };
        let report = report_without_trace(
            run.clone(),
            vec![lambda, lambda],
            2,
            test_root(&[0.0], 1.0, total_budget, Some(0.0), 0),
            vec![
                make_record(0, BipopLane::Large),
                make_record(1, BipopLane::Small),
            ],
            0,
        );

        let error = report
            .validate_ledger()
            .expect_err("production stops immediately at the first target hit");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "continued-after-convergence");
    }

    /// G0/G3: a lucky target hit on the initial callback cannot bypass the
    /// dependency authority that preflight required because the admitted hard
    /// budget could have entered a complete generation.
    #[test]
    fn ledger_early_target_cannot_bypass_reachable_jacobi_authority() {
        let dimension = 4_096usize;
        let lambda = 28usize;
        let total_budget = 29usize;
        let run = CmaReport {
            x_best: vec![0.0; dimension],
            f_best: 0.0,
            evals: 1,
            generations: 0,
            converged: true,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget: total_budget,
            seed: 0,
            start: vec![0.0; dimension],
            trace_start: 0,
            trace_end: 1,
            stop_reason: CmaStopReason::TargetReached,
            report: run.clone(),
        };
        let report = report_without_trace(
            run,
            vec![lambda],
            1,
            test_root(&vec![0.0; dimension], 1.0, total_budget, Some(0.0), 0),
            vec![record],
            0,
        );

        let error = report
            .validate_ledger()
            .expect_err("early convergence cannot weaken callback-free admission");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "eigensolver-admission");
    }

    /// G0/G3: retained decision points must preserve the finite-query guards
    /// enforced by production. Objective outputs remain intentionally outside
    /// this check.
    #[test]
    fn ledger_refuses_nonfinite_start_and_best_points() {
        let make_report = || {
            let mut objective = |_point: &[f64]| 0.0;
            try_bipop_cmaes(&mut objective, &[0.0], 0.5, 1, None, 0)
                .expect("one-callback fixture admits")
        };

        let mut bad_start = make_report();
        bad_start.records[0].start[0] = f64::NAN;
        let error = bad_start
            .validate_ledger()
            .expect_err("non-finite retained start must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-start");

        let mut bad_best = make_report();
        bad_best.records[0].report.x_best[0] = f64::INFINITY;
        let error = bad_best
            .validate_ledger()
            .expect_err("non-finite retained best point must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-best-point");

        let mut bad_sigma = make_report();
        bad_sigma.records[0].report.sigma = f64::INFINITY;
        let error = bad_sigma
            .validate_ledger()
            .expect_err("non-finite derived step size must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-nonnegative-sigma");

        let mut negative_sigma = make_report();
        negative_sigma.records[0].report.sigma = -0.0;
        let error = negative_sigma
            .validate_ledger()
            .expect_err("a forged negative-sign step size must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "finite-nonnegative-sigma");
    }

    /// G3: stale and structurally reassigned callback evidence cannot survive
    /// the production validator, even when every value remains finite.
    #[test]
    fn ledger_refuses_trace_truncation_reordering_reassignment_and_bit_mutation() {
        let make_report = || {
            let mut objective =
                |point: &[f64]| point.iter().map(|value| value * value).sum::<f64>();
            try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, 7)
                .expect("trace mutation fixture admits")
        };

        let mut truncated = make_report();
        truncated.trace_rows.pop();
        let error = truncated
            .validate_ledger()
            .expect_err("truncated trace must refuse");
        assert_eq!(error.invariant(), "trace-length");

        let mut reordered = make_report();
        reordered.trace_rows.swap(0, 1);
        let error = reordered
            .validate_ledger()
            .expect_err("reordered trace rows must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "trace-local-offset");

        let mut duplicated = make_report();
        duplicated.trace_rows[1] = duplicated.trace_rows[0].clone();
        let error = duplicated
            .validate_ledger()
            .expect_err("duplicated trace metadata must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "trace-local-offset");

        let mut reassigned = make_report();
        reassigned.trace_rows[0].restart = 1;
        let error = reassigned
            .validate_ledger()
            .expect_err("reassigned trace row must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "trace-restart");

        let mut changed_point = make_report();
        let first_record = &changed_point.records[0];
        let first_noninitial = first_record
            .trace_start
            .checked_add(1)
            .expect("fixture trace offset fits");
        let nonbest_global = (first_noninitial..first_record.trace_end)
            .find(|global| {
                changed_point.trace_rows[*global].objective.to_bits()
                    != first_record.report.f_best.to_bits()
            })
            .expect("fixture retains a non-best callback");
        let coordinate = nonbest_global
            .checked_mul(changed_point.root.start.len())
            .expect("fixture point offset fits");
        let original = changed_point.trace_points[coordinate];
        changed_point.trace_points[coordinate] = f64::from_bits(original.to_bits() ^ 1);
        let error = changed_point
            .validate_ledger()
            .expect_err("stale trace identity must catch a non-best point mutation");
        assert_eq!(error.restart(), None);
        assert_eq!(error.invariant(), "trace-identity");

        let mut changed_objective = make_report();
        let first_record = &changed_objective.records[0];
        let nonbest_global = (first_record.trace_start..first_record.trace_end)
            .find(|global| {
                changed_objective.trace_rows[*global].objective > first_record.report.f_best
            })
            .expect("fixture retains a strictly non-best callback");
        let original = changed_objective.trace_rows[nonbest_global].objective;
        changed_objective.trace_rows[nonbest_global].objective = f64::from_bits(
            original
                .to_bits()
                .checked_add(1)
                .expect("finite fixture ULP"),
        );
        let error = changed_objective
            .validate_ledger()
            .expect_err("stale trace identity must catch a non-best objective mutation");
        assert_eq!(error.restart(), None);
        assert_eq!(error.invariant(), "trace-identity");
    }

    /// G3: recomputing the nested identities does not make semantically
    /// incompatible root fields authentic for the retained run.
    #[test]
    fn ledger_refuses_resealed_root_start_sigma_budget_target_and_seed_mutations() {
        let make_report = || {
            let mut objective =
                |point: &[f64]| point.iter().map(|value| value * value).sum::<f64>();
            try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, Some(-1.0), 7)
                .expect("root mutation fixture admits")
        };

        let mut changed_start = make_report();
        changed_start.root.start[0] = f64::from_bits(changed_start.root.start[0].to_bits() ^ 1);
        reseal_private_identities(&mut changed_start);
        let error = changed_start
            .validate_ledger()
            .expect_err("resealed root start mutation must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "root-start-projection");

        let mut changed_sigma = make_report();
        changed_sigma.root.sigma *= 2.0;
        reseal_private_identities(&mut changed_sigma);
        let error = changed_sigma
            .validate_ledger()
            .expect_err("resealed root sigma mutation must refuse");
        assert_eq!(error.restart(), Some(1));
        assert_eq!(error.invariant(), "restart-start");

        let mut changed_budget = make_report();
        changed_budget.root.total_budget += 1;
        reseal_private_identities(&mut changed_budget);
        let error = changed_budget
            .validate_ledger()
            .expect_err("resealed root budget mutation must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "allocated-budget");

        let mut changed_target = make_report();
        changed_target.root.target = Some(f64::MAX);
        reseal_private_identities(&mut changed_target);
        let error = changed_target
            .validate_ledger()
            .expect_err("resealed root target mutation must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "terminal-reason");

        let mut changed_seed = make_report();
        changed_seed.root.seed ^= 1;
        reseal_private_identities(&mut changed_seed);
        let error = changed_seed
            .validate_ledger()
            .expect_err("resealed root seed mutation must refuse");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "derived-seed");
    }

    /// G3/G5: a stale but structurally plausible ledger field is a payload
    /// mismatch; resealing makes it a different reference identity; resealing a
    /// semantically impossible field still reaches the independent ledger gate.
    #[test]
    fn study_identity_separates_stale_reference_and_semantic_mutations() {
        let mut objective = |point: &[f64]| point.iter().map(|value| value * value).sum::<f64>();
        let canonical = try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, 11)
            .expect("study identity mutation fixture admits");
        assert!(
            canonical.records.len() > 1,
            "fixture needs a non-best restart"
        );
        let reference = canonical.study_identity;
        let nonbest = if canonical.best_restart == 0 { 1 } else { 0 };
        let assert_payload_mismatch = |report: &BipopReport| match report.validate_study_identity()
        {
            Err(BipopStudyAdmissionError::PayloadIdentityMismatch { declared, computed }) => {
                assert_eq!(declared, reference);
                assert_ne!(computed, reference);
            }
            other => panic!("expected payload identity mismatch, found {other:?}"),
        };

        let mut stale_root = canonical.clone();
        stale_root.root.seed ^= 1;
        assert_payload_mismatch(&stale_root);

        let mut stale_trace_content = canonical.clone();
        let noninitial_coordinate = stale_trace_content
            .root
            .start
            .len()
            .checked_add(1)
            .expect("fixture coordinate offset fits");
        stale_trace_content.trace_points[noninitial_coordinate] =
            f64::from_bits(stale_trace_content.trace_points[noninitial_coordinate].to_bits() ^ 1);
        assert_payload_mismatch(&stale_trace_content);

        let mut stale_trace_receipt = canonical.clone();
        stale_trace_receipt.trace_identity.digest[0] ^= 1;
        assert_payload_mismatch(&stale_trace_receipt);

        let mut stale_schedule = canonical.clone();
        stale_schedule.schedule[0] = stale_schedule.schedule[0]
            .checked_add(1)
            .expect("fixture population fits");
        assert_payload_mismatch(&stale_schedule);

        let mut stale_best_projection = canonical.clone();
        stale_best_projection.best.sigma = f64::from_bits(
            stale_best_projection
                .best
                .sigma
                .to_bits()
                .checked_add(1)
                .expect("finite fixture ULP"),
        );
        assert_payload_mismatch(&stale_best_projection);

        let mut stale = canonical.clone();
        let original_sigma = stale.records[nonbest].report.sigma;
        stale.records[nonbest].report.sigma = f64::from_bits(
            original_sigma
                .to_bits()
                .checked_add(1)
                .expect("finite fixture ULP"),
        );
        let computed = stale
            .computed_study_identity()
            .expect("mutated payload remains identity-representable");
        assert_eq!(
            stale.validate_study_identity(),
            Err(BipopStudyAdmissionError::PayloadIdentityMismatch {
                declared: reference,
                computed,
            })
        );
        let error = stale
            .validate_ledger()
            .expect_err("stale full-study identity must fail the ledger gate");
        assert_eq!(error.invariant(), "study-identity");

        stale.study_identity = computed;
        stale
            .validate_ledger()
            .expect("finite non-best diagnostic sigma has no independent replay oracle");
        assert_eq!(
            stale.admit_study_identity(reference),
            Err(BipopStudyAdmissionError::ReferenceIdentityMismatch {
                expected: reference,
                found: computed,
            })
        );

        let mut semantic_forgery = canonical;
        semantic_forgery.records[0].seed ^= 1;
        semantic_forgery.study_identity = semantic_forgery
            .computed_study_identity()
            .expect("semantic forgery remains identity-representable");
        let error = semantic_forgery
            .admit_study_identity(reference)
            .expect_err("resealing cannot bypass derived-seed admission");
        match error {
            BipopStudyAdmissionError::Ledger { error } => {
                assert_eq!(error.restart(), Some(0));
                assert_eq!(error.invariant(), "derived-seed");
            }
            other => panic!("expected semantic ledger refusal, found {other:?}"),
        }
    }

    /// G3: a correctly resealed `TargetReached` receipt still needs a numeric
    /// witness in its retained callback trace; the total-order best alone is
    /// neither necessary nor sufficient evidence of that terminal condition.
    #[test]
    fn ledger_refuses_resealed_target_terminal_without_callback_witness() {
        let negative_nan = f64::from_bits(0xfff8_0000_0000_0024);
        let values = [negative_nan, 4.0, 0.25, 3.0, 2.0];
        let mut calls = 0usize;
        let mut objective = |_point: &[f64]| {
            let value = values[calls];
            calls += 1;
            value
        };
        let mut report = try_bipop_cmaes(&mut objective, &[1.0], 0.5, 5, Some(0.5), 11)
            .expect("target-witness mutation fixture admits");
        assert_eq!(report.best.f_best.to_bits(), negative_nan.to_bits());
        assert_eq!(report.records[0].stop_reason, CmaStopReason::TargetReached);

        let witness = report
            .trace_rows
            .iter_mut()
            .find(|row| row.objective.to_bits() == 0.25_f64.to_bits())
            .expect("fixture retains its finite target witness");
        witness.objective = 0.75;
        reseal_private_identities(&mut report);

        let error = report
            .validate_ledger()
            .expect_err("resealing cannot invent a missing numeric target witness");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "terminal-reason");
    }

    /// G3/G5: cheap admission intentionally authenticates retained bytes and
    /// structural ledger semantics without replaying objective/CMA state. The
    /// replay-backed gate must reject correctly resealed mutations from every
    /// class that requires causal re-execution to distinguish.
    #[test]
    fn replay_admission_refuses_resealed_sigma_point_objective_and_stop_reason() {
        let make_report = || {
            let mut objective = |point: &[f64]| {
                point
                    .iter()
                    .map(|coordinate| coordinate * coordinate)
                    .sum::<f64>()
            };
            try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, 11)
                .expect("semantic replay mutation fixture admits")
        };
        let assert_replay_refuses =
            |mut forged: BipopReport, canonical_identity: BipopStudyIdentity| {
                reseal_private_identities(&mut forged);
                let forged_identity = forged.study_identity;
                forged
                    .admit_study_identity(forged_identity)
                    .expect("resealed fixture remains cheap-admission consistent");
                let mut calls = 0usize;
                let mut oracle = |point: &[f64]| {
                    calls += 1;
                    point
                        .iter()
                        .map(|coordinate| coordinate * coordinate)
                        .sum::<f64>()
                };
                assert_eq!(
                    forged.admit_study_identity_with_replay(forged_identity, &mut oracle),
                    Err(BipopReplayAdmissionError::SemanticMismatch {
                        retained: forged_identity,
                        replayed: canonical_identity,
                    })
                );
                assert_eq!(calls, forged.total_evals);
            };

        let canonical = make_report();
        let canonical_identity = canonical.study_identity;
        let nonbest_restart = if canonical.best_restart == 0 { 1 } else { 0 };

        let mut changed_sigma = canonical.clone();
        changed_sigma.records[nonbest_restart].report.sigma = f64::from_bits(
            changed_sigma.records[nonbest_restart]
                .report
                .sigma
                .to_bits()
                .checked_add(1)
                .expect("finite sigma fixture admits one ULP"),
        );
        assert_replay_refuses(changed_sigma, canonical_identity);

        let first = &canonical.records[0];
        let nonwinning_global = (first.trace_start + 1..first.trace_end)
            .find(|global| {
                canonical.trace_rows[*global]
                    .objective
                    .total_cmp(&first.report.f_best)
                    .is_gt()
            })
            .expect("fixture retains one noninitial nonwinning callback");

        let mut changed_point = canonical.clone();
        let point_coordinate = nonwinning_global
            .checked_mul(changed_point.root.start.len())
            .expect("fixture coordinate offset fits");
        changed_point.trace_points[point_coordinate] =
            f64::from_bits(changed_point.trace_points[point_coordinate].to_bits() ^ 1);
        assert_replay_refuses(changed_point, canonical_identity);

        let mut changed_objective = canonical.clone();
        let objective_bits = changed_objective.trace_rows[nonwinning_global]
            .objective
            .to_bits();
        changed_objective.trace_rows[nonwinning_global].objective = f64::from_bits(
            objective_bits
                .checked_add(1)
                .expect("finite objective ULP fits"),
        );
        assert_replay_refuses(changed_objective, canonical_identity);

        let mut changed_stop_reason = canonical;
        let plausible_stagnation = changed_stop_reason
            .records
            .iter_mut()
            .find(|record| {
                record.stop_reason == CmaStopReason::BudgetExhausted
                    && record.report.generations > 0
            })
            .expect("fixture contains a generated budget terminal");
        plausible_stagnation.stop_reason = CmaStopReason::Stagnated;
        assert_replay_refuses(changed_stop_reason, canonical_identity);
    }

    /// G0: an overflowing hypothetical next generation is not evidence that a
    /// run exhausted its local budget. The validator must refuse the arithmetic
    /// boundary instead of treating `checked_add(None)` as "no generation fits."
    #[test]
    fn ledger_refuses_wrapped_next_generation_accounting() {
        let lambda = 6usize;
        let evals = usize::MAX - 2;
        let generations = (evals - 1) / lambda;
        assert_eq!(generations * lambda + 1, evals);

        let run = CmaReport {
            x_best: vec![0.0, 0.0],
            f_best: 0.0,
            evals,
            generations,
            converged: false,
            sigma: 1.0,
        };
        let record = BipopRestartRecord {
            schema_version: BIPOP_RESTART_SCHEMA_VERSION,
            ordinal: 0,
            lane: BipopLane::Large,
            lambda,
            allocated_budget: usize::MAX,
            seed: 7,
            start: vec![0.0, 0.0],
            trace_start: 0,
            trace_end: evals,
            stop_reason: CmaStopReason::BudgetExhausted,
            report: run.clone(),
        };
        let report = report_without_trace(
            run,
            vec![lambda],
            evals,
            test_root(&[0.0, 0.0], 1.0, usize::MAX, None, 7),
            vec![record],
            0,
        );

        let error = report
            .validate_ledger()
            .expect_err("overflowing next-generation accounting must fail closed");
        assert_eq!(error.restart(), Some(0));
        assert_eq!(error.invariant(), "evaluation-overflow");
    }
}
