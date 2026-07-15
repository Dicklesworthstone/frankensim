//! GREEDY FIDELITY-LADDER PLANNER (addendum Proposal 8, bead lmp4.16;
//! [F] — behind the `ladder-planner` feature): a LADDER WALK, not a
//! general planner. Governance Rule 1 forbids opening general planning
//! as a research program; the search space is deliberately collapsed to
//! the fidelity-refinement lattice and solved greedily over the operator
//! menu `{cache, speculate, solve-rung, DWR-refine, climb}` with costs
//! LEARNED from telemetry (cold estimates fall back to a conservative
//! default). All the intelligence is inherited from the flywheel
//! underneath: certified verification (Proposal 9's verifier), the
//! content-addressed cache (Proposal 2), the fidelity-ladder registry
//! (Proposal 3), and colors on the returned interval.
//!
//! Determinism (G5): fixed operator order, deterministic tie-breaks
//! (refine preferred over climb on equal predicted cost) — a replayed
//! query reproduces the same operator sequence and interval trajectory.
//! The cannot-discharge boundary hands off cleanly to refusal semantics
//! with the best achieved certified interval — never a false in-budget
//! answer.

use core::fmt;
use std::collections::BTreeMap;

use fs_evidence::{Color, NumericalCertificate, color_leaf_identity_reason, verified_from};
use fs_verify::estimator::{VerifierReceipt, verify_with_receipt};
use fs_verify::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MmsClass, MmsProblem, Poly, gauss5, solve_p1,
};

const MAX_EXACT_CELLS: u128 = 1_u128 << 53;
const MAX_EXACT_CELLS_F64: f64 = 9_007_199_254_740_992.0;
/// Maximum cells a single v0 planner mesh may own.
///
/// Exact floating-point accounting is a weaker constraint than bounded memory:
/// this ceiling is checked before every uniform or adaptive mesh allocation.
/// One mesh with `n` cells owns `n + 1` nodes, so this is derived from the
/// lower-layer MMS node envelope rather than maintained independently.
pub const MAX_PLANNER_CELLS: usize = MAX_FEM1D_MESH_NODES - 1;
/// Maximum polynomial coefficients admitted by the v0 family boundary.
/// This is the equilibrated verifier's exact five-point-Gauss envelope:
/// `(c - F - slope)^2` has degree `2 * (degree(u) - 1) <= 9`.
pub const MAX_FAMILY_COEFFICIENTS: usize = fs_verify::fem1d::MAX_FEM1D_POLY_COEFFICIENTS;
/// Maximum entries in one fidelity ladder.
pub const MAX_LADDER_RUNGS: usize = 4_096;
/// Maximum coefficient-by-cell-by-quadrature-point work admitted to one
/// synchronous v0 solve or verification. This couples the otherwise
/// independent family and mesh caps.
pub const MAX_POLYNOMIAL_CELL_WORK: usize = 16_000_000;
const VERIFIER_GAUSS_POINTS: usize = 5;

/// Semantic version of the retained planner-cache key transport.
///
/// Version 3 is the already-shipped `fs-ir-ladder:v3:` grammar. This constant
/// makes that existing version explicit; it does not rotate or re-key it.
pub const PLANNER_CACHE_KEY_VERSION: u32 = 3;
/// Exact domain and prefix of every canonical planner-cache key.
pub const PLANNER_CACHE_KEY_DOMAIN: &str = "fs-ir-ladder:v3:";
const PLANNER_CACHE_KEY_PREFIX_STEM: &str = "fs-ir-ladder:v";

/// Owner-local declaration consumed by `xtask check-identities`.
pub const PLANNER_CACHE_KEY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ir:planner-cache-key",
    "version_const=PLANNER_CACHE_KEY_VERSION",
    "version=3",
    "domain=fs-ir-ladder:v3:",
    "domain_const=PLANNER_CACHE_KEY_DOMAIN",
    "encoder=cache_key",
    "encoder_helpers=cache_key_with_prefix,ProblemFamily::scaled_class",
    "schema_constants=PLANNER_CACHE_KEY_VERSION,PLANNER_CACHE_KEY_DOMAIN,PLANNER_CACHE_KEY_PREFIX_STEM",
    "schema_functions=ProblemFamily::base,ProblemFamily::kernel,validate_finite,canonicalize_zero,canonical_f64_bits",
    "schema_dependencies=fs-verify:fem1d-mms-class",
    "digest=none-exact-canonical-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=ProblemFamily",
    "source_fields=ProblemFamily.base_class:semantic",
    "source_bindings=ProblemFamily.base_class>base-class-canonical-identity",
    "external_semantic_fields=domain-prefix,identity-version,theta-scaled-class-exact-bits",
    "semantic_fields=domain-prefix,identity-version,base-class-canonical-identity,theta-scaled-class-exact-bits",
    "excluded_fields=base-class-signed-zero-and-trailing-zero-spelling:normalized-before-ProblemFamily-admission,theta-signed-zero-sign:normalized-by-scaled-MmsClass-admission,tolerance:lookup-filter-and-independent-reverification-only,budget-and-ladder:execution-policy-not-answer-identity",
    "consumers=plan_observed,AnswerCache,MemCache,retained-planner-cache",
    "mutations=domain-prefix:crates/fs-ir/src/planner.rs#planner_cache_schema_and_domain_move_identity,identity-version:crates/fs-ir/src/planner.rs#planner_cache_schema_and_domain_move_identity,base-class-canonical-identity:crates/fs-ir/src/planner.rs#planner_cache_base_class_identity_moves_key,theta-scaled-class-exact-bits:crates/fs-ir/src/planner.rs#planner_cache_theta_exact_bits_move_key",
    "nonsemantic_mutations=base-class-signed-zero-and-trailing-zero-spelling:crates/fs-ir/src/planner.rs#planner_cache_intentional_normalizations_do_not_move_identity,theta-signed-zero-sign:crates/fs-ir/src/planner.rs#planner_cache_intentional_normalizations_do_not_move_identity,tolerance:crates/fs-ir/src/planner.rs#planner_cache_execution_policy_does_not_move_answer_identity,budget-and-ladder:crates/fs-ir/src/planner.rs#planner_cache_execution_policy_does_not_move_answer_identity",
    "field_guard=classify_planner_cache_identity_fields",
    "transport_guard=admit_planner_cache_key",
    "version_guard=crates/fs-ir/src/planner.rs#planner_cache_key_versions_fail_closed",
    "coupling_surface=fs-ir:planner-cache-key",
];

/// One planner operator (the menu).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PlanOp {
    /// Proposal 2: the content-addressed answer cache.
    CacheLookup,
    /// Proposal 9: verify a prolongated coarse answer without solving.
    Speculate,
    /// Solve at the current rung (uniform mesh).
    SolveRung,
    /// Refine ONLY where the residual indicators concentrate.
    DwrRefine,
    /// Move to the next rung (uniform, finer everywhere).
    Climb,
}

impl PlanOp {
    /// Stable name for logs and the cost table.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            PlanOp::CacheLookup => "cache",
            PlanOp::Speculate => "speculate",
            PlanOp::SolveRung => "solve-rung",
            PlanOp::DwrRefine => "dwr-refine",
            PlanOp::Climb => "climb",
        }
    }
}

/// A structured refusal at the planner's public validity boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanError {
    /// A scalar planner input violates its finite/range contract.
    InvalidScalar {
        /// Input name.
        field: &'static str,
        /// Stable rejection reason.
        reason: &'static str,
    },
    /// A required ladder or coefficient sequence is empty.
    EmptySequence {
        /// Sequence name.
        field: &'static str,
    },
    /// One sequence entry violates ordering, range, or finiteness.
    InvalidSequenceEntry {
        /// Sequence name.
        field: &'static str,
        /// Offending entry.
        index: usize,
        /// Stable rejection reason.
        reason: &'static str,
    },
    /// The problem-family definition is not usable by this kernel.
    InvalidFamily {
        /// Family field (`kernel`, `base`, or `boundary`).
        field: &'static str,
        /// Offending coefficient, when applicable.
        index: Option<usize>,
        /// Stable rejection reason.
        reason: &'static str,
    },
    /// A mesh violates the 1-D verifier's topology/domain contract.
    InvalidMesh {
        /// Offending node, when applicable.
        index: Option<usize>,
        /// Stable rejection reason.
        reason: &'static str,
    },
    /// A nodal candidate cannot be evaluated on its declared mesh.
    InvalidCandidate {
        /// Offending node, when applicable.
        index: Option<usize>,
        /// Stable rejection reason.
        reason: &'static str,
    },
    /// Cost telemetry overflowed while updating an operator mean.
    CostOverflow {
        /// Operator whose telemetry overflowed.
        op: PlanOp,
    },
    /// A resource-driving public value exceeds the v0 planner envelope.
    ResourceLimit {
        /// Resource or public field name.
        field: &'static str,
        /// Requested count.
        requested: usize,
        /// Maximum admitted count.
        limit: usize,
    },
    /// A bounded vector allocation failed before numerical work began.
    AllocationFailed {
        /// Stable allocation stage.
        stage: &'static str,
        /// Elements requested from the allocator.
        requested: usize,
    },
    /// A checked numerical stage produced an unusable result.
    NumericalFailure {
        /// Stable computation stage.
        stage: &'static str,
    },
    /// The bounded FEM solver refused or failed with structured context.
    Fem1dFailure {
        /// Planner stage invoking the solver.
        stage: &'static str,
        /// Original FEM boundary failure.
        source: Fem1dError,
    },
}

impl fmt::Display for PlanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidScalar { field, reason } => {
                write!(f, "invalid planner scalar `{field}`: {reason}")
            }
            Self::EmptySequence { field } => {
                write!(f, "planner sequence `{field}` must not be empty")
            }
            Self::InvalidSequenceEntry {
                field,
                index,
                reason,
            } => write!(f, "invalid `{field}` entry {index}: {reason}"),
            Self::InvalidFamily {
                field,
                index: Some(index),
                reason,
            } => write!(f, "invalid family `{field}` entry {index}: {reason}"),
            Self::InvalidFamily {
                field,
                index: None,
                reason,
            } => write!(f, "invalid family field `{field}`: {reason}"),
            Self::InvalidMesh {
                index: Some(index),
                reason,
            } => write!(f, "invalid mesh node {index}: {reason}"),
            Self::InvalidMesh {
                index: None,
                reason,
            } => write!(f, "invalid mesh: {reason}"),
            Self::InvalidCandidate {
                index: Some(index),
                reason,
            } => write!(f, "invalid candidate node {index}: {reason}"),
            Self::InvalidCandidate {
                index: None,
                reason,
            } => write!(f, "invalid candidate: {reason}"),
            Self::CostOverflow { op } => {
                write!(f, "cost telemetry overflowed for `{}`", op.name())
            }
            Self::ResourceLimit {
                field,
                requested,
                limit,
            } => write!(
                f,
                "planner resource `{field}` requested {requested}, limit is {limit}"
            ),
            Self::AllocationFailed { stage, requested } => {
                write!(
                    f,
                    "planner allocation failed during {stage} for {requested} elements"
                )
            }
            Self::NumericalFailure { stage } => {
                write!(f, "planner numerical failure during {stage}")
            }
            Self::Fem1dFailure { stage, source } => {
                write!(f, "planner FEM failure during {stage}: {source}")
            }
        }
    }
}

impl std::error::Error for PlanError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fem1dFailure { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Fail-closed reason for refusing a retained planner-cache key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlannerCacheKeyAdmissionError {
    /// The key does not use the planner-cache namespace/version grammar.
    MalformedPrefix,
    /// The key declares a well-formed but unsupported schema version.
    UnsupportedVersion {
        /// Version declared by retained state.
        declared: u32,
        /// Exact version supported by this build.
        supported: u32,
    },
    /// The version parses as current but is not in its one canonical spelling.
    NonCanonicalPrefix,
    /// The lower-layer canonical bytes are not canonical lowercase hex.
    MalformedPayload,
}

impl fmt::Display for PlannerCacheKeyAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MalformedPrefix => f.write_str("malformed planner-cache key prefix"),
            Self::UnsupportedVersion {
                declared,
                supported,
            } => write!(
                f,
                "planner-cache key v{declared} is unsupported; this build accepts exactly v{supported}"
            ),
            Self::NonCanonicalPrefix => {
                f.write_str("planner-cache key version prefix is not canonically encoded")
            }
            Self::MalformedPayload => {
                f.write_str("planner-cache key payload is not non-empty lowercase hexadecimal")
            }
        }
    }
}

impl std::error::Error for PlannerCacheKeyAdmissionError {}

/// Admit a retained cache key under the exact current prefix and return its
/// canonical lower-layer payload.
///
/// Unknown versions, aliases such as `v03`, changed domains, uppercase hex,
/// and malformed payloads are refused rather than guessed or normalized.
///
/// # Errors
/// Returns [`PlannerCacheKeyAdmissionError`] for any non-current or
/// non-canonical retained key.
pub fn admit_planner_cache_key(key: &str) -> Result<&str, PlannerCacheKeyAdmissionError> {
    let versioned = key
        .strip_prefix(PLANNER_CACHE_KEY_PREFIX_STEM)
        .ok_or(PlannerCacheKeyAdmissionError::MalformedPrefix)?;
    let (version_text, payload) = versioned
        .split_once(':')
        .ok_or(PlannerCacheKeyAdmissionError::MalformedPrefix)?;
    if version_text.is_empty() || !version_text.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(PlannerCacheKeyAdmissionError::MalformedPrefix);
    }
    let declared = version_text
        .parse::<u32>()
        .map_err(|_| PlannerCacheKeyAdmissionError::MalformedPrefix)?;
    if declared != PLANNER_CACHE_KEY_VERSION {
        return Err(PlannerCacheKeyAdmissionError::UnsupportedVersion {
            declared,
            supported: PLANNER_CACHE_KEY_VERSION,
        });
    }
    if !key.starts_with(PLANNER_CACHE_KEY_DOMAIN) {
        return Err(PlannerCacheKeyAdmissionError::NonCanonicalPrefix);
    }
    if payload.is_empty()
        || !payload.len().is_multiple_of(2)
        || !payload
            .bytes()
            .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f'))
    {
        return Err(PlannerCacheKeyAdmissionError::MalformedPayload);
    }
    Ok(payload)
}

/// A cache-declared answer. The planner independently re-verifies it before any
/// certificate or discharge decision is emitted.
#[derive(Debug, Clone, PartialEq)]
pub struct CachedAnswer {
    nodal: Vec<f64>,
    bound: f64,
    mesh: Vec<f64>,
}

impl CachedAnswer {
    /// Construct a structurally checked cache record.
    ///
    /// This does not trust the claimed certificate: [`plan`] independently
    /// replays the verifier before using any returned cache entry.
    ///
    /// # Errors
    /// Returns [`PlanError`] for an invalid bound, mesh, or nodal vector.
    pub fn new(mut nodal: Vec<f64>, bound: f64, mut mesh: Vec<f64>) -> Result<Self, PlanError> {
        validate_bound(bound, "cached_bound")?;
        validate_mesh(&mesh)?;
        validate_candidate(&mesh, &nodal)?;
        canonicalize_mesh_zeros(&mut mesh);
        canonicalize_mesh_zeros(&mut nodal);
        Ok(Self { nodal, bound, mesh })
    }

    /// Cached nodal values.
    #[must_use]
    pub fn nodal(&self) -> &[f64] {
        &self.nodal
    }

    /// Cache-declared bound. The planner never trusts it without replay.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Cached mesh.
    #[must_use]
    pub fn mesh(&self) -> &[f64] {
        &self.mesh
    }
}

/// The Proposal-2 cache seam (implemented over the content-addressed
/// store; the planner only needs lookup/insert semantics).
pub trait AnswerCache {
    /// A candidate for `key` whose declared bound is ≤ `tol`, if any. The
    /// implementation is not a trust root; [`plan`] re-verifies the candidate.
    /// Implementations backed by retained state must fail closed through
    /// [`admit_planner_cache_key`] before interpreting a key.
    fn lookup(&self, key: &str, tol: f64) -> Option<CachedAnswer>;
    /// Record an answer that the planner has just independently verified. A
    /// retained implementation must refuse non-current/non-canonical keys.
    fn insert(&mut self, key: &str, answer: CachedAnswer);
}

/// The trivial in-memory cache.
#[derive(Debug, Default)]
pub struct MemCache {
    items: BTreeMap<String, CachedAnswer>,
}

impl AnswerCache for MemCache {
    fn lookup(&self, key: &str, tol: f64) -> Option<CachedAnswer> {
        admit_planner_cache_key(key).ok()?;
        self.items.get(key).filter(|a| a.bound <= tol).cloned()
    }

    fn insert(&mut self, key: &str, answer: CachedAnswer) {
        if admit_planner_cache_key(key).is_ok() {
            self.items.insert(key.to_string(), answer);
        }
    }
}

/// LEARNED cost table: mean observed cost (cells solved) per operator.
/// Cold entries fall back to the conservative default.
#[derive(Debug)]
pub struct CostTable {
    seen: BTreeMap<&'static str, (f64, u64)>,
    cold_default: f64,
}

impl CostTable {
    /// Construct a table with a finite, strictly positive cold fallback.
    ///
    /// # Errors
    /// Returns [`PlanError`] when the fallback could poison operator choice.
    pub fn new(cold_default: f64) -> Result<CostTable, PlanError> {
        validate_positive_finite(cold_default, "cold_default")?;
        Ok(CostTable {
            seen: BTreeMap::new(),
            cold_default,
        })
    }

    /// Record one finite, strictly positive observed cost.
    ///
    /// # Errors
    /// Returns [`PlanError`] without mutating the table when the sample or
    /// accumulator is unusable.
    pub fn record(&mut self, op: PlanOp, cost: f64) -> Result<(), PlanError> {
        validate_positive_finite(cost, "cost_sample")?;
        let e = self.seen.entry(op.name()).or_insert((0.0, 0));
        let sum = e.0 + cost;
        let count = e.1.checked_add(1).ok_or(PlanError::CostOverflow { op })?;
        if !sum.is_finite() {
            return Err(PlanError::CostOverflow { op });
        }
        *e = (sum, count);
        Ok(())
    }

    /// Predict an operator's cost (learned mean, else the default).
    #[must_use]
    pub fn predict(&self, op: PlanOp) -> f64 {
        self.seen
            .get(op.name())
            .filter(|(_, n)| *n > 0)
            .map_or(self.cold_default, |(sum, n)| {
                #[allow(clippy::cast_precision_loss)]
                {
                    sum / *n as f64
                }
            })
    }

    /// The checked cold fallback used for unseen operators.
    #[must_use]
    pub fn cold_default(&self) -> f64 {
        self.cold_default
    }
}

/// One executed step in the plan (the audit trail).
#[derive(Debug, Clone, PartialEq)]
pub struct OpLog {
    /// Which operator ran.
    pub op: PlanOp,
    /// Cells involved (the cost unit).
    pub cost: f64,
    /// The certified bound after the step (∞ before any verify).
    pub bound_after: f64,
    /// Full verifier authority after the step, when one was produced.
    pub certificate_after: Option<VerifierCertificate>,
}

/// An equilibrated enclosure together with the verifier identity that minted it.
///
/// Fields are private so callers cannot detach a `Verified` color from the
/// verifier family and reconstructed-flux identity checked by the planner.
#[derive(Debug, Clone, PartialEq)]
pub struct VerifierCertificate {
    receipt: VerifierReceipt,
    color: Color,
}

impl VerifierCertificate {
    /// Certified energy-error half width.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.receipt.bound_hi()
    }

    /// Evidence color minted through `fs-evidence`'s guarded enclosure door.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// Stable verifier-family identity.
    #[must_use]
    pub fn verifier_family(&self) -> &str {
        self.receipt.verifier_family()
    }

    /// Identity of the verifier's reconstructed flux.
    #[must_use]
    pub fn flux_hash(&self) -> u64 {
        self.receipt.flux_hash()
    }

    /// Immutable lower-layer receipt emitted by the exact verification run
    /// that minted this certificate.
    #[must_use]
    pub fn receipt(&self) -> &VerifierReceipt {
        &self.receipt
    }
}

/// The planner verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanOutcome {
    /// The query is discharged: the certified bound meets tolerance.
    Discharged {
        /// The nodal answer on its mesh.
        nodal: Vec<f64>,
        /// The final mesh.
        mesh: Vec<f64>,
        /// The certified energy bound (VERIFIED color: an equilibrated
        /// enclosure, never a DWR guess).
        bound: f64,
        /// Full verifier authority for `bound`.
        certificate: VerifierCertificate,
        /// The executed operator sequence.
        ops: Vec<OpLog>,
        /// Total cost spent (cells).
        cost: f64,
    },
    /// The budget could not discharge the query: hand off to refusal
    /// semantics with the BEST ACHIEVED certified interval — never a
    /// false in-budget answer.
    RefusedWithBest {
        /// The best certified bound achieved.
        best_bound: f64,
        /// Full verifier authority for `best_bound`.
        best_certificate: VerifierCertificate,
        /// The nodal answer that achieved it.
        best_nodal: Vec<f64>,
        /// Its mesh.
        best_mesh: Vec<f64>,
        /// The executed operator sequence.
        ops: Vec<OpLog>,
        /// Total cost spent (never exceeds the admitted budget).
        cost: f64,
        /// What to tell the caller.
        reason: String,
    },
    /// The request was valid, but its budget could not fund even one solve and
    /// no independently verified cache answer existed. No interval or color is
    /// fabricated.
    RefusedWithoutAnswer {
        /// Operators attempted before the refusal (normally cache lookup only).
        ops: Vec<OpLog>,
        /// Total solved-cell cost spent.
        cost: f64,
        /// Teaching refusal.
        reason: String,
    },
}

/// The 1-D elliptic problem family the v0 planner discharges (the
/// verifier's kernel class): exact solution `theta`-scaled.
#[derive(Debug, Clone, PartialEq)]
pub struct ProblemFamily {
    base_class: MmsClass,
}

#[allow(dead_code)]
fn classify_planner_cache_identity_fields(family: &ProblemFamily) {
    let ProblemFamily { base_class: _ } = family;
}

impl ProblemFamily {
    /// Construct a checked 1-D manufactured-problem family.
    ///
    /// # Errors
    /// Returns [`PlanError`] for a malformed kernel identity, empty/non-finite
    /// polynomial, or non-homogeneous boundary values.
    pub fn new(base: Poly, kernel: impl AsRef<str>) -> Result<Self, PlanError> {
        let kernel = kernel.as_ref();
        if let Some(reason) = color_leaf_identity_reason(kernel) {
            return Err(PlanError::InvalidFamily {
                field: "kernel",
                index: None,
                reason,
            });
        }
        let base_class = MmsClass::new(kernel, base).map_err(|source| PlanError::Fem1dFailure {
            stage: "problem-family admission",
            source,
        })?;
        Ok(Self { base_class })
    }

    /// Exact-solution polynomial at `theta = 1`.
    #[must_use]
    pub fn base(&self) -> &Poly {
        self.base_class.exact_solution()
    }

    /// Stable kernel identity.
    #[must_use]
    pub fn kernel(&self) -> &str {
        self.base_class.name()
    }

    /// Canonical lower-layer bytes of the unscaled family class.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        self.base_class.canonical_bytes()
    }

    fn scaled_class(&self, theta: f64) -> Result<MmsClass, PlanError> {
        validate_finite(theta, "theta")?;
        let coefficients = self.base().coefficients();
        let mut scaled = Vec::new();
        scaled
            .try_reserve_exact(coefficients.len())
            .map_err(|_| PlanError::AllocationFailed {
                stage: "problem-family scaling",
                requested: coefficients.len(),
            })?;
        for (index, coefficient) in coefficients.iter().enumerate() {
            let value = coefficient * theta;
            if !value.is_finite() {
                return Err(PlanError::InvalidFamily {
                    field: "scaled_base",
                    index: Some(index),
                    reason: "coefficient scaling is non-finite",
                });
            }
            scaled.push(canonicalize_zero(value));
        }
        let scaled = Poly::new(scaled).map_err(|source| PlanError::Fem1dFailure {
            stage: "scaled problem-family polynomial",
            source,
        })?;
        MmsClass::new(self.kernel(), scaled).map_err(|source| PlanError::Fem1dFailure {
            stage: "scaled problem-family admission",
            source,
        })
    }

    /// Instantiate the checked problem at `theta` on an arbitrary valid mesh.
    ///
    /// # Errors
    /// Returns [`PlanError`] for invalid `theta`/mesh values or coefficient
    /// overflow while scaling/differentiating the family.
    pub fn at(&self, theta: f64, mesh: Vec<f64>) -> Result<MmsProblem, PlanError> {
        validate_mesh(&mesh)?;
        validate_family_cell_work(self, mesh.len() - 1)?;
        let class = self.scaled_class(theta)?;
        MmsProblem::from_class(class, mesh).map_err(|source| PlanError::Fem1dFailure {
            stage: "problem-family mesh admission",
            source,
        })
    }
}

pub(crate) fn validate_finite(value: f64, field: &'static str) -> Result<(), PlanError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PlanError::InvalidScalar {
            field,
            reason: "must be finite",
        })
    }
}

pub(crate) fn validate_positive_finite(value: f64, field: &'static str) -> Result<(), PlanError> {
    if !value.is_finite() {
        Err(PlanError::InvalidScalar {
            field,
            reason: "must be finite",
        })
    } else if value <= 0.0 {
        Err(PlanError::InvalidScalar {
            field,
            reason: "must be strictly positive",
        })
    } else {
        Ok(())
    }
}

pub(crate) fn validate_bound(bound: f64, field: &'static str) -> Result<(), PlanError> {
    if !bound.is_finite() {
        Err(PlanError::InvalidScalar {
            field,
            reason: "must be finite",
        })
    } else if bound < 0.0 {
        Err(PlanError::InvalidScalar {
            field,
            reason: "must be non-negative",
        })
    } else {
        Ok(())
    }
}

pub(crate) fn validate_budget(value: f64, field: &'static str) -> Result<(), PlanError> {
    validate_positive_finite(value, field)?;
    if value > MAX_EXACT_CELLS_F64 {
        Err(PlanError::InvalidScalar {
            field,
            reason: "exceeds exact f64 cell accounting",
        })
    } else {
        Ok(())
    }
}

fn validate_family_cell_work(family: &ProblemFamily, cells: usize) -> Result<(), PlanError> {
    let requested = family
        .base()
        .coefficients()
        .len()
        .checked_mul(cells)
        .and_then(|work| work.checked_mul(VERIFIER_GAUSS_POINTS))
        .ok_or(PlanError::ResourceLimit {
            field: "polynomial_cell_work",
            requested: usize::MAX,
            limit: MAX_POLYNOMIAL_CELL_WORK,
        })?;
    if requested > MAX_POLYNOMIAL_CELL_WORK {
        Err(PlanError::ResourceLimit {
            field: "polynomial_cell_work",
            requested,
            limit: MAX_POLYNOMIAL_CELL_WORK,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn validate_rung_cells(rung_cells: &[usize]) -> Result<(), PlanError> {
    if rung_cells.is_empty() {
        return Err(PlanError::EmptySequence {
            field: "rung_cells",
        });
    }
    if rung_cells.len() > MAX_LADDER_RUNGS {
        return Err(PlanError::ResourceLimit {
            field: "rung_cells",
            requested: rung_cells.len(),
            limit: MAX_LADDER_RUNGS,
        });
    }
    for (index, &cells) in rung_cells.iter().enumerate() {
        if cells == 0 {
            return Err(PlanError::InvalidSequenceEntry {
                field: "rung_cells",
                index,
                reason: "cell count must be non-zero",
            });
        }
        if widened_usize(cells) > MAX_EXACT_CELLS {
            return Err(PlanError::InvalidSequenceEntry {
                field: "rung_cells",
                index,
                reason: "cell count exceeds exact f64 budget accounting",
            });
        }
        if cells > MAX_PLANNER_CELLS {
            return Err(PlanError::ResourceLimit {
                field: "rung_cells",
                requested: cells,
                limit: MAX_PLANNER_CELLS,
            });
        }
        if index > 0 && cells <= rung_cells[index - 1] {
            return Err(PlanError::InvalidSequenceEntry {
                field: "rung_cells",
                index,
                reason: "cell counts must be strictly increasing",
            });
        }
    }
    Ok(())
}

fn validate_mesh(mesh: &[f64]) -> Result<(), PlanError> {
    if mesh.len() < 2 {
        return Err(PlanError::InvalidMesh {
            index: None,
            reason: "at least two boundary nodes are required",
        });
    }
    let cells = mesh.len() - 1;
    if widened_usize(cells) > MAX_EXACT_CELLS {
        return Err(PlanError::InvalidMesh {
            index: None,
            reason: "cell count exceeds exact f64 budget accounting",
        });
    }
    if cells > MAX_PLANNER_CELLS {
        return Err(PlanError::ResourceLimit {
            field: "mesh_cells",
            requested: cells,
            limit: MAX_PLANNER_CELLS,
        });
    }
    for (index, node) in mesh.iter().enumerate() {
        if !node.is_finite() {
            return Err(PlanError::InvalidMesh {
                index: Some(index),
                reason: "node must be finite",
            });
        }
    }
    if canonical_f64_bits(mesh[0]) != 0 {
        return Err(PlanError::InvalidMesh {
            index: Some(0),
            reason: "first node must equal 0",
        });
    }
    if canonical_f64_bits(*mesh.last().ok_or(PlanError::InvalidMesh {
        index: None,
        reason: "at least two boundary nodes are required",
    })?) != 1.0_f64.to_bits()
    {
        return Err(PlanError::InvalidMesh {
            index: Some(mesh.len() - 1),
            reason: "last node must equal 1",
        });
    }
    for (index, pair) in mesh.windows(2).enumerate() {
        if pair[0].partial_cmp(&pair[1]) != Some(core::cmp::Ordering::Less) {
            return Err(PlanError::InvalidMesh {
                index: Some(index + 1),
                reason: "nodes must be strictly increasing",
            });
        }
        // Align with fs-verify's `admit_mesh` (fem1d `NonFiniteReciprocal`): a
        // cell so narrow that `1/(b-a)` overflows to ±inf blows up the P1
        // stiffness assembly. fs-verify rejects it, so fs-ir must too — else a
        // CachedAnswer with such a mesh is admitted here but rejected at replay
        // (a divergence against this hardening's own "aligned with admit_mesh"
        // goal), and it can never be discharged into a certificate.
        if !(1.0 / (pair[1] - pair[0])).is_finite() {
            return Err(PlanError::InvalidMesh {
                index: Some(index + 1),
                reason: "cell is too narrow: 1/(b-a) is non-finite (fs-verify admit_mesh rejects it)",
            });
        }
    }
    Ok(())
}

fn validate_candidate(mesh: &[f64], nodal: &[f64]) -> Result<(), PlanError> {
    if nodal.len() != mesh.len() {
        return Err(PlanError::InvalidCandidate {
            index: None,
            reason: "nodal length must equal mesh-node length",
        });
    }
    for (index, value) in nodal.iter().enumerate() {
        if !value.is_finite() {
            return Err(PlanError::InvalidCandidate {
                index: Some(index),
                reason: "nodal value must be finite",
            });
        }
    }
    if canonical_f64_bits(nodal[0]) != 0
        || canonical_f64_bits(*nodal.last().ok_or(PlanError::InvalidCandidate {
            index: None,
            reason: "nodal values must include both boundaries",
        })?) != 0
    {
        return Err(PlanError::InvalidCandidate {
            index: None,
            reason: "homogeneous boundary values must equal zero",
        });
    }
    Ok(())
}

fn canonical_f64_bits(value: f64) -> u64 {
    const SIGN_BIT: u64 = 1_u64 << 63;
    match value.to_bits() {
        SIGN_BIT => 0,
        bits => bits,
    }
}

#[allow(clippy::cast_lossless)]
const fn widened_usize(value: usize) -> u128 {
    value as u128
}

fn canonicalize_zero(value: f64) -> f64 {
    f64::from_bits(canonical_f64_bits(value))
}

fn canonicalize_mesh_zeros(mesh: &mut [f64]) {
    for node in mesh {
        *node = canonicalize_zero(*node);
    }
}

fn meshes_have_same_nodes(left: &[f64], right: &[f64]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| canonical_f64_bits(*left) == canonical_f64_bits(*right))
}

fn uniform_mesh(cells: usize) -> Result<Vec<f64>, PlanError> {
    validate_rung_cells(&[cells])?;
    let nodes = cells.checked_add(1).ok_or(PlanError::ResourceLimit {
        field: "mesh_nodes",
        requested: cells,
        limit: MAX_PLANNER_CELLS + 1,
    })?;
    let mut mesh = Vec::new();
    mesh.try_reserve_exact(nodes)
        .map_err(|_| PlanError::AllocationFailed {
            stage: "uniform mesh",
            requested: nodes,
        })?;
    #[allow(clippy::cast_precision_loss)]
    for k in 0..=cells {
        mesh.push(k as f64 / cells as f64);
    }
    validate_mesh(&mesh)?;
    Ok(mesh)
}

/// Per-element energy-residual indicators (the same integrand the
/// equilibrated verifier bounds, localized): `∫_K (c* − F − u′)²`.
fn element_indicators(problem: &MmsProblem, nodal: &[f64]) -> Result<Vec<f64>, PlanError> {
    let m = problem.mesh();
    validate_mesh(m)?;
    validate_candidate(m, nodal)?;
    // The verifier's optimal constant.
    let mut c_star = 0.0f64;
    for (element, nodes) in m.windows(2).enumerate() {
        let h = nodes[1] - nodes[0];
        let slope = (nodal[element + 1] - nodal[element]) / h;
        for (gx, gw) in gauss5(nodes[0], nodes[1]) {
            c_star += gw * (problem.rounded_forcing_antiderivative().eval(gx) + slope);
        }
    }
    if !c_star.is_finite() {
        return Err(PlanError::NumericalFailure {
            stage: "residual indicator flux constant",
        });
    }
    let indicator_count = m.len() - 1;
    let mut indicators = Vec::new();
    indicators
        .try_reserve_exact(indicator_count)
        .map_err(|_| PlanError::AllocationFailed {
            stage: "residual indicators",
            requested: indicator_count,
        })?;
    for (element, nodes) in m.windows(2).enumerate() {
        let h = nodes[1] - nodes[0];
        let slope = (nodal[element + 1] - nodal[element]) / h;
        let mut acc = 0.0f64;
        for (gx, gw) in gauss5(nodes[0], nodes[1]) {
            let r = c_star - problem.rounded_forcing_antiderivative().eval(gx) - slope;
            acc += gw * r * r;
        }
        if !acc.is_finite() || acc < 0.0 {
            return Err(PlanError::NumericalFailure {
                stage: "residual indicator integration",
            });
        }
        indicators.push(acc);
    }
    Ok(indicators)
}

/// EQUIDISTRIBUTION refinement (the textbook optimal-mesh criterion):
/// split every element whose squared-residual contribution exceeds the
/// per-element target `tol²/n`, with a PER-ELEMENT depth from its own
/// gap (splitting an element into 2^d pieces cuts its contribution by
/// ~4^d in this residual model). Deterministic; converges in a couple
/// of solve rounds instead of crawling at the tail.
fn refine_to_target(mesh: &[f64], indicators: &[f64], tol: f64) -> Result<Vec<f64>, PlanError> {
    validate_mesh(mesh)?;
    validate_positive_finite(tol, "tolerance")?;
    if indicators.len() + 1 != mesh.len() {
        return Err(PlanError::InvalidSequenceEntry {
            field: "indicators",
            index: indicators.len(),
            reason: "indicator count must equal mesh cell count",
        });
    }
    if let Some(index) = indicators
        .iter()
        .position(|indicator| !indicator.is_finite() || *indicator < 0.0)
    {
        return Err(PlanError::InvalidSequenceEntry {
            field: "indicators",
            index,
            reason: "indicator must be finite and non-negative",
        });
    }
    let target = refinement_target(indicators.len(), tol);
    let refined_cells = refined_cell_count(indicators, target)?;
    let refined_nodes = refined_cells
        .checked_add(1)
        .ok_or(PlanError::ResourceLimit {
            field: "adaptive_mesh_nodes",
            requested: refined_cells,
            limit: MAX_PLANNER_CELLS + 1,
        })?;
    let mut out = Vec::new();
    out.try_reserve_exact(refined_nodes)
        .map_err(|_| PlanError::AllocationFailed {
            stage: "adaptive mesh",
            requested: refined_nodes,
        })?;
    for (element, nodes) in mesh.windows(2).enumerate() {
        out.push(nodes[0]);
        let pieces = refinement_pieces(indicators[element], target);
        if pieces > 1 {
            let (a, b) = (nodes[0], nodes[1]);
            for k in 1..pieces {
                #[allow(clippy::cast_precision_loss)]
                out.push(a + (b - a) * k as f64 / pieces as f64);
            }
        }
    }
    out.push(1.0);
    canonicalize_mesh_zeros(&mut out);
    validate_mesh(&out)?;
    Ok(out)
}

fn refinement_pieces(indicator: f64, target: f64) -> usize {
    if indicator <= target {
        return 1;
    }
    let gap = indicator / target;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let depth = if gap.is_finite() {
        ((gap.log2() / 2.0).ceil() as u32).clamp(1, 5)
    } else {
        5
    };
    1usize << depth
}

fn refinement_target(indicator_count: usize, tol: f64) -> f64 {
    let count = indicator_count.max(1);
    #[allow(clippy::cast_precision_loss)]
    let target = tol * tol / count as f64;
    target.max(f64::MIN_POSITIVE)
}

fn refined_cell_count(indicators: &[f64], target: f64) -> Result<usize, PlanError> {
    let mut cells = 0usize;
    for indicator in indicators {
        cells = cells
            .checked_add(refinement_pieces(*indicator, target))
            .ok_or(PlanError::ResourceLimit {
                field: "adaptive_mesh_cells",
                requested: usize::MAX,
                limit: MAX_PLANNER_CELLS,
            })?;
        if cells > MAX_PLANNER_CELLS {
            return Err(PlanError::ResourceLimit {
                field: "adaptive_mesh_cells",
                requested: cells,
                limit: MAX_PLANNER_CELLS,
            });
        }
    }
    Ok(cells)
}

/// Deterministic P1 interpolation from an arbitrary valid coarse mesh onto an
/// arbitrary valid target mesh over the same `[0,1]` domain.
fn prolong_linear(
    coarse_mesh: &[f64],
    coarse_nodal: &[f64],
    target_mesh: &[f64],
) -> Result<Vec<f64>, PlanError> {
    validate_mesh(coarse_mesh)?;
    validate_candidate(coarse_mesh, coarse_nodal)?;
    validate_mesh(target_mesh)?;

    let mut output = Vec::new();
    output
        .try_reserve_exact(target_mesh.len())
        .map_err(|_| PlanError::AllocationFailed {
            stage: "linear prolongation",
            requested: target_mesh.len(),
        })?;
    let mut coarse_element = 0usize;
    for &x in target_mesh {
        while coarse_element + 1 < coarse_mesh.len() - 1 && x > coarse_mesh[coarse_element + 1] {
            coarse_element += 1;
        }
        let a = coarse_mesh[coarse_element];
        let b = coarse_mesh[coarse_element + 1];
        let value = if canonical_f64_bits(x) == canonical_f64_bits(a) {
            coarse_nodal[coarse_element]
        } else if canonical_f64_bits(x) == canonical_f64_bits(b) {
            coarse_nodal[coarse_element + 1]
        } else {
            let fraction = (x - a) / (b - a);
            coarse_nodal[coarse_element]
                + fraction * (coarse_nodal[coarse_element + 1] - coarse_nodal[coarse_element])
        };
        if !value.is_finite() {
            return Err(PlanError::NumericalFailure {
                stage: "linear prolongation",
            });
        }
        output.push(canonicalize_zero(value));
    }
    validate_candidate(target_mesh, &output)?;
    Ok(output)
}

#[derive(Debug, Clone)]
struct CheckedVerification {
    certificate: VerifierCertificate,
    accept: bool,
}

fn checked_verify(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
) -> Result<CheckedVerification, PlanError> {
    validate_candidate(problem.mesh(), candidate)?;
    let receipt = verify_with_receipt(problem, candidate, tolerance).map_err(|_| {
        PlanError::NumericalFailure {
            stage: "verifier receipt production",
        }
    })?;
    let lo = receipt.bound_lo();
    let hi = receipt.bound_hi();
    if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || lo > hi {
        return Err(PlanError::NumericalFailure {
            stage: "verifier enclosure",
        });
    }
    if canonical_f64_bits(receipt.tolerance()) != canonical_f64_bits(tolerance) {
        return Err(PlanError::NumericalFailure {
            stage: "verifier tolerance consistency",
        });
    }
    let expected_accept = hi <= tolerance;
    if receipt.accepted() != expected_accept {
        return Err(PlanError::NumericalFailure {
            stage: "verifier acceptance consistency",
        });
    }
    let guarded_color = verified_from(&NumericalCertificate::enclosure(0.0, hi)).map_err(|_| {
        PlanError::NumericalFailure {
            stage: "verifier evidence-color admission",
        }
    })?;
    let receipt_color = receipt.color();
    match (receipt_color.as_ref(), expected_accept) {
        (
            Some(Color::Verified {
                lo: color_lo,
                hi: color_hi,
            }),
            true,
        ) if canonical_f64_bits(*color_lo) == 0
            && canonical_f64_bits(*color_hi) == canonical_f64_bits(hi) => {}
        (None, false) => {}
        _ => {
            return Err(PlanError::NumericalFailure {
                stage: "verifier color consistency",
            });
        }
    }
    Ok(CheckedVerification {
        certificate: VerifierCertificate {
            receipt,
            color: receipt_color.unwrap_or(guarded_color),
        },
        accept: expected_accept,
    })
}

fn cache_key(family: &ProblemFamily, theta: f64) -> Result<String, PlanError> {
    cache_key_with_prefix(family, theta, PLANNER_CACHE_KEY_DOMAIN)
}

fn cache_key_with_prefix(
    family: &ProblemFamily,
    theta: f64,
    prefix: &str,
) -> Result<String, PlanError> {
    let class = family.scaled_class(theta)?;
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let bytes = class.canonical_bytes();
    let requested = bytes
        .len()
        .checked_mul(2)
        .and_then(|encoded| encoded.checked_add(prefix.len()))
        .ok_or(PlanError::AllocationFailed {
            stage: "planner canonical cache key",
            requested: usize::MAX,
        })?;
    let mut key = String::new();
    key.try_reserve_exact(requested)
        .map_err(|_| PlanError::AllocationFailed {
            stage: "planner canonical cache key",
            requested,
        })?;
    key.push_str(prefix);
    for byte in bytes {
        key.push(char::from(HEX[usize::from(byte >> 4)]));
        key.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    Ok(key)
}

fn solved_cell_cost(mesh: &[f64]) -> Result<f64, PlanError> {
    validate_mesh(mesh)?;
    cell_count_cost(mesh.len() - 1)
}

fn cell_count_cost(cells: usize) -> Result<f64, PlanError> {
    if cells == 0 || cells > MAX_PLANNER_CELLS {
        return Err(PlanError::ResourceLimit {
            field: "operator_cells",
            requested: cells,
            limit: MAX_PLANNER_CELLS,
        });
    }
    #[allow(clippy::cast_precision_loss)]
    let cost = cells as f64;
    validate_positive_finite(cost, "operator_cost")?;
    Ok(cost)
}

fn can_afford(spent: f64, next_cost: f64, budget: f64) -> bool {
    next_cost <= budget - spent
}

#[derive(Debug, Clone, Copy)]
struct PendingTransition {
    op: PlanOp,
    observed_cost: f64,
}

impl PendingTransition {
    fn new(op: PlanOp) -> Self {
        Self {
            op,
            observed_cost: 0.0,
        }
    }

    fn observe(&mut self, cost: f64) -> Result<(), PlanError> {
        validate_positive_finite(cost, "transition_observed_cost")?;
        let observed_cost = self.observed_cost + cost;
        if !observed_cost.is_finite() {
            return Err(PlanError::CostOverflow { op: self.op });
        }
        self.observed_cost = observed_cost;
        Ok(())
    }
}

fn record_observed_transition(
    costs: &mut CostTable,
    pending: &mut Option<PendingTransition>,
) -> Result<(), PlanError> {
    if pending
        .as_ref()
        .is_some_and(|transition| transition.observed_cost > 0.0)
    {
        let transition = pending.take().ok_or(PlanError::NumericalFailure {
            stage: "pending transition accounting",
        })?;
        costs.record(transition.op, transition.observed_cost)?;
    }
    Ok(())
}

/// Whether an operational planner observer permits later work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PlanControl {
    Continue,
    Stop,
}

/// Operational checkpoints exposed to the anytime adapter.
pub(crate) trait PlanObserver {
    /// Called before any paid operation or its resource allocation begins.
    fn before_work(
        &mut self,
        spent: f64,
        next_cost: f64,
        best: Option<(&VerifierCertificate, &[f64])>,
        costs: &CostTable,
    ) -> Result<PlanControl, PlanError>;

    /// Called immediately after a verifier certificate and its current-work
    /// telemetry are committed, before any later transition, cache insertion,
    /// telemetry, allocation, or numerical work.
    fn certified(
        &mut self,
        spent: f64,
        best: (&VerifierCertificate, &[f64]),
        costs: &CostTable,
    ) -> Result<PlanControl, PlanError>;
}

struct ContinueObserver;

impl PlanObserver for ContinueObserver {
    fn before_work(
        &mut self,
        _spent: f64,
        _next_cost: f64,
        _best: Option<(&VerifierCertificate, &[f64])>,
        _costs: &CostTable,
    ) -> Result<PlanControl, PlanError> {
        Ok(PlanControl::Continue)
    }

    fn certified(
        &mut self,
        _spent: f64,
        _best: (&VerifierCertificate, &[f64]),
        _costs: &CostTable,
    ) -> Result<PlanControl, PlanError> {
        Ok(PlanControl::Continue)
    }
}

pub(crate) struct ObservedPlan {
    pub(crate) outcome: PlanOutcome,
    pub(crate) stopped: bool,
}

fn finished_plan(outcome: PlanOutcome) -> ObservedPlan {
    ObservedPlan {
        outcome,
        stopped: false,
    }
}

fn refuse(
    best: Option<(VerifierCertificate, Vec<f64>, Vec<f64>)>,
    ops: Vec<OpLog>,
    cost: f64,
    reason: String,
) -> PlanOutcome {
    match best {
        Some((best_certificate, best_nodal, best_mesh)) => PlanOutcome::RefusedWithBest {
            best_bound: best_certificate.bound(),
            best_certificate,
            best_nodal,
            best_mesh,
            ops,
            cost,
            reason,
        },
        None => PlanOutcome::RefusedWithoutAnswer { ops, cost, reason },
    }
}

fn best_ref(
    best: Option<&(VerifierCertificate, Vec<f64>, Vec<f64>)>,
) -> Option<(&VerifierCertificate, &[f64])> {
    best.map(|(certificate, _, mesh)| (certificate, mesh.as_slice()))
}

fn stopped_plan(
    best: Option<(VerifierCertificate, Vec<f64>, Vec<f64>)>,
    ops: Vec<OpLog>,
    cost: f64,
) -> ObservedPlan {
    ObservedPlan {
        outcome: refuse(
            best,
            ops,
            cost,
            "execution stopped by the anytime observer before later planner work".to_string(),
        ),
        stopped: true,
    }
}

/// The greedy ladder walk. `rung_cells` is the fidelity lattice
/// (coarsest first); `budget_cells` is the cost budget in solved cells.
///
/// # Errors
/// Returns [`PlanError`] before issuing a certificate when any public input,
/// cache replay, telemetry update, or numerical intermediate is unusable.
#[allow(clippy::too_many_lines)]
pub fn plan(
    family: &ProblemFamily,
    theta: f64,
    tol: f64,
    budget_cells: f64,
    rung_cells: &[usize],
    cache: &mut dyn AnswerCache,
    costs: &mut CostTable,
) -> Result<PlanOutcome, PlanError> {
    let observed = plan_observed(
        family,
        theta,
        tol,
        budget_cells,
        rung_cells,
        cache,
        costs,
        &mut ContinueObserver,
    )?;
    if observed.stopped {
        return Err(PlanError::NumericalFailure {
            stage: "non-stoppable planner observer",
        });
    }
    Ok(observed.outcome)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn plan_observed(
    family: &ProblemFamily,
    theta: f64,
    tol: f64,
    budget_cells: f64,
    rung_cells: &[usize],
    cache: &mut dyn AnswerCache,
    costs: &mut CostTable,
    observer: &mut dyn PlanObserver,
) -> Result<ObservedPlan, PlanError> {
    validate_finite(theta, "theta")?;
    validate_positive_finite(tol, "tolerance")?;
    validate_budget(budget_cells, "budget_cells")?;
    validate_rung_cells(rung_cells)?;
    for &cells in rung_cells {
        validate_family_cell_work(family, cells)?;
    }
    // Re-instantiation below checks theta-scaled derived coefficients too.
    let key = cache_key(family, theta)?;
    let mut ops: Vec<OpLog> = Vec::new();
    let mut spent = 0.0f64;
    let mut best: Option<(VerifierCertificate, Vec<f64>, Vec<f64>)> = None;
    // ---- Operator 1: the cache (zero solved-cell cost on a hit). Cache
    // metadata is never evidence: independently replay the verifier.
    if let Some(hit) = cache.lookup(&key, tol) {
        let replay = family
            .at(theta, hit.mesh.clone())
            .and_then(|problem| checked_verify(&problem, &hit.nodal, tol));
        if let Ok(checked) = replay
            && checked.accept
        {
            let certificate = checked.certificate;
            ops.push(OpLog {
                op: PlanOp::CacheLookup,
                cost: 0.0,
                bound_after: certificate.bound(),
                certificate_after: Some(certificate.clone()),
            });
            best = Some((certificate.clone(), hit.nodal.clone(), hit.mesh.clone()));
            let best_certificate = best_ref(best.as_ref()).ok_or(PlanError::NumericalFailure {
                stage: "missing best certificate after cache replay",
            })?;
            if observer.certified(0.0, best_certificate, costs)? == PlanControl::Stop {
                return Ok(stopped_plan(best, ops, 0.0));
            }
            return Ok(finished_plan(PlanOutcome::Discharged {
                nodal: hit.nodal,
                mesh: hit.mesh,
                bound: certificate.bound(),
                certificate,
                ops,
                cost: 0.0,
            }));
        }
    }
    ops.push(OpLog {
        op: PlanOp::CacheLookup,
        cost: 0.0,
        bound_after: f64::INFINITY,
        certificate_after: None,
    });
    let initial_cost = cell_count_cost(rung_cells[0])?;
    validate_family_cell_work(family, rung_cells[0])?;
    if observer.before_work(spent, initial_cost, best_ref(best.as_ref()), costs)?
        == PlanControl::Stop
    {
        return Ok(stopped_plan(best, ops, spent));
    }
    if !can_afford(spent, initial_cost, budget_cells) {
        return Ok(finished_plan(refuse(
            best,
            ops,
            spent,
            format!(
                "budget {budget_cells} cells cannot fund the initial {initial_cost:.0}-cell solve; \
                 refusal occurs before allocation, so no mesh is allocated and no uncertified \
                 answer is returned"
            ),
        )));
    }
    let mut rung = 0usize;
    let mut mesh = uniform_mesh(rung_cells[0])?;
    let mut carried: Option<Vec<f64>> = None; // prolongated candidate
    let mut pending_transition: Option<PendingTransition> = None;
    loop {
        // ---- Operator 2: speculate — verify a carried candidate
        // WITHOUT solving (prolongation from the previous rung).
        if let Some(cand) = carried.take() {
            let vcost = 0.2 * solved_cell_cost(&mesh)?;
            if observer.before_work(spent, vcost, best_ref(best.as_ref()), costs)?
                == PlanControl::Stop
            {
                return Ok(stopped_plan(best, ops, spent));
            }
            if !can_afford(spent, vcost, budget_cells) {
                return Ok(finished_plan(refuse(
                    best,
                    ops,
                    spent,
                    format!(
                        "budget {budget_cells} cells cannot fund speculative verification \
                         costing {vcost:.1}; no unverified candidate is returned — hand off to \
                         refusal/anytime semantics"
                    ),
                )));
            }
            let problem = family.at(theta, mesh.clone())?;
            spent += vcost;
            let checked = checked_verify(&problem, &cand, tol)?;
            let transition = pending_transition
                .as_mut()
                .ok_or(PlanError::NumericalFailure {
                    stage: "speculation without a pending climb",
                })?;
            if transition.op != PlanOp::Climb {
                return Err(PlanError::NumericalFailure {
                    stage: "speculation attached to a non-climb transition",
                });
            }
            transition.observe(vcost)?;
            record_observed_transition(costs, &mut pending_transition)?;
            costs.record(PlanOp::Speculate, vcost)?;
            let certificate = checked.certificate;
            ops.push(OpLog {
                op: PlanOp::Speculate,
                cost: vcost,
                bound_after: certificate.bound(),
                certificate_after: Some(certificate.clone()),
            });
            let better = best
                .as_ref()
                .is_none_or(|(best, _, _)| certificate.bound() < best.bound());
            if better {
                best = Some((certificate.clone(), cand.clone(), mesh.clone()));
            }
            let best_certificate = best_ref(best.as_ref()).ok_or(PlanError::NumericalFailure {
                stage: "missing best certificate after speculation",
            })?;
            if observer.certified(spent, best_certificate, costs)? == PlanControl::Stop {
                return Ok(stopped_plan(best, ops, spent));
            }
            if checked.accept {
                cache.insert(
                    &key,
                    CachedAnswer::new(cand.clone(), certificate.bound(), mesh.clone())?,
                );
                return Ok(finished_plan(PlanOutcome::Discharged {
                    nodal: cand,
                    mesh,
                    bound: certificate.bound(),
                    certificate,
                    ops,
                    cost: spent,
                }));
            }
        }
        // ---- Operator 3: solve at the current rung.
        let scost = solved_cell_cost(&mesh)?;
        if observer.before_work(spent, scost, best_ref(best.as_ref()), costs)? == PlanControl::Stop
        {
            return Ok(stopped_plan(best, ops, spent));
        }
        if !can_afford(spent, scost, budget_cells) {
            // A climb may already have paid for speculative verification; that
            // observed work is real. A zero-cost refine/climb whose first
            // downstream compute never ran remains unobserved and is dropped.
            record_observed_transition(costs, &mut pending_transition)?;
            return Ok(finished_plan(refuse(
                best,
                ops,
                spent,
                format!(
                    "budget {budget_cells} cells cannot fund the next {scost:.0}-cell solve; \
                     no uncertified answer is returned — hand off to refusal/anytime semantics"
                ),
            )));
        }
        let problem = family.at(theta, mesh.clone())?;
        spent += scost;
        let nodal = solve_p1(&problem).map_err(|source| PlanError::Fem1dFailure {
            stage: "solve rung",
            source,
        })?;
        validate_candidate(&mesh, &nodal)?;
        let checked = checked_verify(&problem, &nodal, tol)?;
        if let Some(transition) = pending_transition.as_mut() {
            transition.observe(scost)?;
        }
        record_observed_transition(costs, &mut pending_transition)?;
        costs.record(PlanOp::SolveRung, scost)?;
        let certificate = checked.certificate;
        ops.push(OpLog {
            op: PlanOp::SolveRung,
            cost: scost,
            bound_after: certificate.bound(),
            certificate_after: Some(certificate.clone()),
        });
        let better = best
            .as_ref()
            .is_none_or(|(best, _, _)| certificate.bound() < best.bound());
        if better {
            best = Some((certificate.clone(), nodal.clone(), mesh.clone()));
        }
        let best_certificate = best_ref(best.as_ref()).ok_or(PlanError::NumericalFailure {
            stage: "missing best certificate after solve",
        })?;
        if observer.certified(spent, best_certificate, costs)? == PlanControl::Stop {
            return Ok(stopped_plan(best, ops, spent));
        }
        if checked.accept {
            cache.insert(
                &key,
                CachedAnswer::new(nodal.clone(), certificate.bound(), mesh.clone())?,
            );
            return Ok(finished_plan(PlanOutcome::Discharged {
                nodal,
                mesh,
                bound: certificate.bound(),
                certificate,
                ops,
                cost: spent,
            }));
        }
        // Predictions select a deterministic zero-cost transition. They are
        // not an admission oracle: the exact verification/solve cost is
        // checked at the top of the next iteration before work begins.
        let next_refine = costs.predict(PlanOp::DwrRefine);
        let next_climb = costs.predict(PlanOp::Climb);
        let climb_available = rung + 1 < rung_cells.len();
        let choose_refine = !climb_available || next_refine <= next_climb;
        // ---- Greedy choice: refine-where-indicated vs climb, by
        // learned predicted cost; deterministic tie-break prefers
        // DwrRefine (the cheaper-in-principle local move).
        if choose_refine {
            let indicators = element_indicators(&problem, &nodal)?;
            let target = refinement_target(indicators.len(), tol);
            let refined_cells = refined_cell_count(&indicators, target)?;
            let refined_cost = cell_count_cost(refined_cells)?;
            validate_family_cell_work(family, refined_cells)?;
            if observer.before_work(spent, refined_cost, best_ref(best.as_ref()), costs)?
                == PlanControl::Stop
            {
                return Ok(stopped_plan(best, ops, spent));
            }
            if !can_afford(spent, refined_cost, budget_cells) {
                return Ok(finished_plan(refuse(
                    best,
                    ops,
                    spent,
                    format!(
                        "budget {budget_cells} cells cannot fund the next {refined_cost:.0}-cell \
                         adaptive solve; refusal occurs before any refined mesh is allocated"
                    ),
                )));
            }
            let refined = refine_to_target(&mesh, &indicators, tol)?;
            if meshes_have_same_nodes(&mesh, &refined) {
                return Ok(finished_plan(refuse(
                    best,
                    ops,
                    spent,
                    "refinement produced no new representable mesh node; refusing instead of \
                     repeating the same solve"
                        .to_string(),
                )));
            }
            ops.push(OpLog {
                op: PlanOp::DwrRefine,
                cost: 0.0, // the refine itself is bookkeeping; the solve pays
                bound_after: f64::INFINITY,
                certificate_after: None,
            });
            mesh = refined;
            carried = None;
            pending_transition = Some(PendingTransition::new(PlanOp::DwrRefine));
        } else {
            rung += 1;
            let speculative_cost = 0.2 * cell_count_cost(rung_cells[rung])?;
            validate_family_cell_work(family, rung_cells[rung])?;
            if observer.before_work(spent, speculative_cost, best_ref(best.as_ref()), costs)?
                == PlanControl::Stop
            {
                return Ok(stopped_plan(best, ops, spent));
            }
            if !can_afford(spent, speculative_cost, budget_cells) {
                return Ok(finished_plan(refuse(
                    best,
                    ops,
                    spent,
                    format!(
                        "budget {budget_cells} cells cannot fund {speculative_cost:.1}-cell \
                         verification on rung {}; refusal occurs before any fine mesh is allocated",
                        rung_cells[rung]
                    ),
                )));
            }
            let fine_mesh = uniform_mesh(rung_cells[rung])?;
            // Interpolate over the actual current mesh. It may be non-dyadic,
            // nonuniform, or locally finer than this uniform target after DWR.
            let cand = prolong_linear(&mesh, &nodal, &fine_mesh)?;
            ops.push(OpLog {
                op: PlanOp::Climb,
                cost: 0.0,
                bound_after: f64::INFINITY,
                certificate_after: None,
            });
            mesh = fine_mesh;
            carried = Some(cand);
            pending_transition = Some(PendingTransition::new(PlanOp::Climb));
        }
    }
}

/// The FIXED BASELINE the kill criterion measures against: a single
/// mid-rung solve, then UNIFORM refinement until the tolerance is met
/// (no cache, no speculation, no locality).
///
/// # Errors
/// Returns [`PlanError`] for invalid inputs, cell-count overflow, or an
/// unusable solver/verifier result.
pub fn baseline_uniform(
    family: &ProblemFamily,
    theta: f64,
    tol: f64,
    mid_rung_cells: usize,
    max_doublings: usize,
) -> Result<(f64, f64), PlanError> {
    validate_finite(theta, "theta")?;
    validate_positive_finite(tol, "tolerance")?;
    validate_rung_cells(&[mid_rung_cells])?;
    let mut cells = mid_rung_cells;
    let mut spent = 0.0f64;
    for _ in 0..=max_doublings {
        validate_family_cell_work(family, cells)?;
        let mesh = uniform_mesh(cells)?;
        let problem = family.at(theta, mesh)?;
        spent += solved_cell_cost(problem.mesh())?;
        if !spent.is_finite() {
            return Err(PlanError::NumericalFailure {
                stage: "baseline cost accumulation",
            });
        }
        let nodal = solve_p1(&problem).map_err(|source| PlanError::Fem1dFailure {
            stage: "uniform baseline solve",
            source,
        })?;
        validate_candidate(problem.mesh(), &nodal)?;
        let checked = checked_verify(&problem, &nodal, tol)?;
        if checked.accept {
            return Ok((spent, checked.certificate.bound()));
        }
        cells = cells.checked_mul(2).ok_or(PlanError::NumericalFailure {
            stage: "baseline rung doubling",
        })?;
        if widened_usize(cells) > MAX_EXACT_CELLS {
            return Err(PlanError::InvalidSequenceEntry {
                field: "baseline_cells",
                index: 0,
                reason: "cell count exceeds exact f64 budget accounting",
            });
        }
    }
    Ok((spent, f64::INFINITY))
}

#[cfg(test)]
mod tests {
    use super::{
        CostTable, MemCache, PLANNER_CACHE_KEY_DOMAIN, PLANNER_CACHE_KEY_VERSION, PlanOutcome,
        PlannerCacheKeyAdmissionError, ProblemFamily, admit_planner_cache_key, cache_key,
        cache_key_with_prefix, plan, prolong_linear,
    };
    use fs_verify::fem1d::Poly;

    fn family(coefficients: Vec<f64>, kernel: &str) -> ProblemFamily {
        ProblemFamily::new(
            Poly::new(coefficients).expect("valid planner test polynomial"),
            kernel,
        )
        .expect("valid planner test family")
    }

    #[test]
    fn adaptive_mesh_prolongates_to_unrelated_uniform_nodes() {
        let coarse_mesh = [0.0, 0.1, 0.4, 0.9, 1.0];
        let coarse_nodal = [0.0, 0.2, 0.6, 0.1, 0.0];
        let target_mesh = [0.0, 0.25, 0.5, 0.75, 1.0];
        let prolonged = prolong_linear(&coarse_mesh, &coarse_nodal, &target_mesh).unwrap();
        let expected = [0.0, 0.4, 0.5, 0.25, 0.0];
        for (actual, expected) in prolonged.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-15);
        }
    }

    #[test]
    fn non_dyadic_prolongation_reproduces_coarse_nodes_exactly() {
        let coarse_mesh = [0.0, 0.3, 0.55, 1.0];
        let coarse_nodal = [0.0, 0.7, -0.2, 0.0];
        let prolonged = prolong_linear(&coarse_mesh, &coarse_nodal, &coarse_mesh).unwrap();
        assert_eq!(prolonged, coarse_nodal);
    }

    #[test]
    fn plan_retains_verifier_receipt_through_audit_and_outcome() {
        let receipt_family = family(vec![0.0, 0.2, -0.2, 0.0, 1.0, -1.0], "receipt-propagation");
        let mut costs = CostTable::new(200.0).expect("valid cold cost");
        let outcome = plan(
            &receipt_family,
            1.0,
            0.05,
            2_000.0,
            &[12, 24, 48, 96],
            &mut MemCache::default(),
            &mut costs,
        )
        .expect("receipt-propagation plan");
        let (certificate, ops) = match &outcome {
            PlanOutcome::Discharged {
                certificate, ops, ..
            } => (certificate, ops),
            PlanOutcome::RefusedWithBest {
                best_certificate,
                ops,
                ..
            } => (best_certificate, ops),
            PlanOutcome::RefusedWithoutAnswer { reason, .. } => {
                panic!("the propagation fixture must produce a receipt: {reason}")
            }
        };
        let audited = ops
            .iter()
            .rev()
            .find_map(|entry| entry.certificate_after.as_ref())
            .expect("audit trail retains the verifier receipt");

        assert_eq!(audited.receipt(), certificate.receipt());
        assert_eq!(
            audited.receipt().artifact_root(),
            certificate.receipt().artifact_root()
        );
        assert_eq!(
            audited.receipt().candidate_root(),
            certificate.receipt().candidate_root()
        );
    }

    #[test]
    fn cache_key_remains_v3_prefix_plus_lower_class_hex() {
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let key = cache_key(&ordinary, 1.0).unwrap();
        let scaled = ordinary.scaled_class(1.0).unwrap();
        let encoded = key
            .strip_prefix(PLANNER_CACHE_KEY_DOMAIN)
            .expect("cache-key schema prefix");
        assert_eq!(encoded.len(), scaled.canonical_bytes().len() * 2);
        for (pair, expected) in encoded
            .as_bytes()
            .chunks_exact(2)
            .zip(scaled.canonical_bytes())
        {
            let pair = core::str::from_utf8(pair).expect("ASCII hexadecimal cache key");
            assert_eq!(u8::from_str_radix(pair, 16).unwrap(), *expected);
        }
    }

    #[test]
    fn planner_cache_intentional_normalizations_do_not_move_identity() {
        let normalized = family(vec![-0.0, 1.0, -1.0, -0.0], "elliptic");
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        assert_eq!(normalized.canonical_bytes(), ordinary.canonical_bytes());
        let normalized_key = cache_key(&normalized, 1.0).unwrap();
        let ordinary_key = cache_key(&ordinary, 1.0).unwrap();
        assert_eq!(normalized_key, ordinary_key);
        assert_eq!(
            cache_key(&ordinary, -0.0).unwrap(),
            cache_key(&ordinary, 0.0).unwrap()
        );
    }

    #[test]
    fn planner_cache_execution_policy_does_not_move_answer_identity() {
        fn key_under_policy(
            family: &ProblemFamily,
            theta: f64,
            _tolerance: f64,
            _budget_cells: usize,
            _rung_cells: &[usize],
        ) -> String {
            cache_key(family, theta).expect("valid cache identity")
        }

        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let strict = key_under_policy(&ordinary, 1.0, 1.0e-12, 1_024, &[8, 16, 32]);
        let permissive = key_under_policy(&ordinary, 1.0, 1.0e-3, 65_536, &[64, 256]);
        assert_eq!(strict, permissive);
    }

    #[test]
    fn planner_cache_base_class_identity_moves_key() {
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let renamed = family(vec![0.0, 1.0, -1.0], "elliptic-renamed");
        let rescaled = family(vec![0.0, 2.0, -2.0], "elliptic");
        assert_ne!(ordinary.canonical_bytes(), renamed.canonical_bytes());
        assert_ne!(ordinary.canonical_bytes(), rescaled.canonical_bytes());
        assert_ne!(
            cache_key(&ordinary, 1.0).unwrap(),
            cache_key(&renamed, 1.0).unwrap()
        );
        assert_ne!(
            cache_key(&ordinary, 1.0).unwrap(),
            cache_key(&rescaled, 1.0).unwrap()
        );
    }

    #[test]
    fn planner_cache_theta_exact_bits_move_key() {
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let theta = 1.0_f64;
        let next_theta = f64::from_bits(theta.to_bits() + 1);
        assert_ne!(theta.to_bits(), next_theta.to_bits());
        assert_eq!(format!("{theta:.6e}"), format!("{next_theta:.6e}"));
        assert_ne!(
            cache_key(&ordinary, theta).unwrap(),
            cache_key(&ordinary, next_theta).unwrap()
        );
    }

    #[test]
    fn planner_cache_schema_and_domain_move_identity() {
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let current = cache_key(&ordinary, 1.0).unwrap();
        let next_schema = cache_key_with_prefix(&ordinary, 1.0, "fs-ir-ladder:v4:").unwrap();
        let next_domain = cache_key_with_prefix(&ordinary, 1.0, "fs-ir-shadow-ladder:v3:").unwrap();
        assert_ne!(current, next_schema);
        assert_ne!(current, next_domain);
        assert_ne!(next_schema, next_domain);
    }

    #[test]
    fn planner_cache_key_replay_is_deterministic() {
        let first = family(vec![0.0, 1.0, -1.0], "elliptic");
        let replay = family(vec![-0.0, 1.0, -1.0, 0.0], "elliptic");
        assert_eq!(
            cache_key(&first, 1.25).unwrap(),
            cache_key(&replay, 1.25).unwrap()
        );
    }

    #[test]
    fn planner_cache_key_versions_fail_closed() {
        let ordinary = family(vec![0.0, 1.0, -1.0], "elliptic");
        let current = cache_key(&ordinary, 1.0).unwrap();
        assert_eq!(PLANNER_CACHE_KEY_VERSION, 3);
        assert_eq!(
            admit_planner_cache_key(&current).unwrap(),
            current.strip_prefix(PLANNER_CACHE_KEY_DOMAIN).unwrap()
        );

        let stale = cache_key_with_prefix(&ordinary, 1.0, "fs-ir-ladder:v2:").unwrap();
        assert_eq!(
            admit_planner_cache_key(&stale),
            Err(PlannerCacheKeyAdmissionError::UnsupportedVersion {
                declared: 2,
                supported: PLANNER_CACHE_KEY_VERSION,
            })
        );
        let aliased = cache_key_with_prefix(&ordinary, 1.0, "fs-ir-ladder:v03:").unwrap();
        assert_eq!(
            admit_planner_cache_key(&aliased),
            Err(PlannerCacheKeyAdmissionError::NonCanonicalPrefix)
        );
        let wrong_domain =
            cache_key_with_prefix(&ordinary, 1.0, "fs-ir-shadow-ladder:v3:").unwrap();
        assert_eq!(
            admit_planner_cache_key(&wrong_domain),
            Err(PlannerCacheKeyAdmissionError::MalformedPrefix)
        );
        assert_eq!(
            admit_planner_cache_key("fs-ir-ladder:v3:ABC0"),
            Err(PlannerCacheKeyAdmissionError::MalformedPayload)
        );
    }
}
