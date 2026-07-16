//! SHEAF-ADJUDICATED THREE-WAY MERGE (addendum Proposal 10, bead
//! lmp4.12 — THE CROWN JEWEL; [M], behind the `sheaf-merge` feature
//! until its Gauntlet tier + kill metric are green): the sheaf
//! machinery built for watertightness supplies useful operators, but the
//! fixed-iteration decomposition is a guarded merge heuristic rather than an
//! H¹ certifier:
//!
//! - the union cochain receives a deterministic fixed-iteration gauge
//!   reconciliation → applied automatically only when the nominal post-gauge
//!   residual independently passes the supplied tolerance;
//! - a dominant decomposition remainder is retained only as a CANDIDATE
//!   merge conflict, localized to its supra-tolerance interface cells with
//!   both caller-supplied parent labels. It proves neither H¹, non-exactness,
//!   authenticated provenance, nor that no local repair exists;
//! - NON-GEOMETRIC assignment payloads (load cases, materials) currently
//!   refuse before decomposition. Without the base assignment map, even two
//!   different branch values do not prove a three-way conflict: either value
//!   may be unchanged from the base. No assignment is silently dropped;
//! - trust is CONDITIONED on the complex's spectral gap (Proposal 5,
//!   risk R5): merges in degraded-gap regions are flagged
//!   low-confidence.

use crate::sheaf_repair::{
    AdmittedSheafSkeleton, RepairAccountant, RepairAdmission, SheafRepairBudget, SheafRepairError,
    SheafRepairUsage, SheafSkeleton, SheafSkeletonError, admit_repair_budget, apply_gauge,
    bounded_d0, hodge_decompose, hodge_decompose_accounted, planned_cost, validate_bounded_cochain,
    zeroed_output_bounded,
};
use fs_exec::Cx;
use std::collections::{BTreeMap, BTreeSet};

/// One branch's edits: a mismatch cochain plus non-geometric keyed
/// assignments (load cases, materials — the typed layer's inputs).
#[derive(Debug, Clone)]
pub struct BranchState {
    /// Caller-supplied parent label (for example, a commit or agent label).
    /// This string is not an authenticated provenance identity.
    pub provenance: String,
    /// The branch's interface mismatch cochain.
    pub mismatch: Vec<f64>,
    /// Non-geometric keyed assignments.
    pub assignments: BTreeMap<String, String>,
}

/// A dominant fixed-iteration decomposition remainder that prevents automatic
/// merge at the requested tolerance.
///
/// This candidate proves neither non-exactness nor an H¹/topology claim.
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateRemainderConflict {
    /// Supra-tolerance interface cells (patch pairs) with magnitudes,
    /// strongest first. This is not the complete mathematical support.
    pub cells: Vec<((usize, usize), f64)>,
    /// Caller-supplied parent labels (X, Y); not authenticated identities.
    pub parents: (String, String),
}

/// A pairwise assignment-difference candidate: the two branch payloads carry
/// different values for the same key. Without a base assignment map this does
/// not prove that both branches edited or that a three-way conflict exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TypeConflict {
    /// The colliding key.
    pub key: String,
    /// X's value.
    pub x_value: String,
    /// Y's value.
    pub y_value: String,
}

/// Heuristic merge-confidence diagnostic conditioned on the current spectral
/// gap proxy (Proposal 5 / R5). This is not certification authority.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Confidence {
    /// The measured proxy meets the supplied threshold. The fixed-iteration
    /// eigensolve does not certify convergence or separation.
    Normal {
        /// The algebraic-connectivity gap.
        gap: f64,
    },
    /// Degraded measured proxy: the candidate-remainder/gauge-fit split is
    /// especially fragile here — treat the merge as provisional.
    LowGap {
        /// The measured gap.
        gap: f64,
        /// The threshold it fell below.
        threshold: f64,
    },
}

/// Local nominal post-gauge threshold observation attached to auto-resolved
/// merges. The v1 fields are not a portable replay/content-bound receipt.
///
/// This is a deterministic `f64` threshold check, not an interval
/// watertightness certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeResidualReceipt {
    /// Worst post-reconciliation interface mismatch.
    pub post_norm: f64,
    /// The tolerance it passed.
    pub tol: f64,
}

/// The merge verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum MergeOutcome {
    /// Inputs could not support an honest merge decision.
    Refused {
        /// Stable, actionable refusal reason.
        reason: &'static str,
    },
    /// A boundary fast path fired (X == Y, or one branch unchanged):
    /// no decomposition needed. This proves edit equality only and carries no
    /// post-state residual receipt.
    Trivial {
        /// Which fast path.
        reason: &'static str,
        /// The merged cochain.
        merged: Vec<f64>,
    },
    /// The deterministic gauge reconciliation reached the supplied nominal
    /// residual tolerance.
    Resolved {
        /// The reconciled cochain.
        merged: Vec<f64>,
        /// The gauge offsets applied.
        gauge: Vec<f64>,
        /// The checked local nominal residual observation.
        residual_receipt: MergeResidualReceipt,
        /// Gap-conditioned trust.
        confidence: Confidence,
    },
    /// Operational numerical merge conflicts. Candidate remainders make no
    /// topology claim. The assignment field is retained for the future
    /// base-aware API; the current entry point refuses all assignment payloads
    /// before producing this variant.
    Conflicted {
        /// Dominant fixed-iteration remainders.
        candidate_remainders: Vec<CandidateRemainderConflict>,
        /// Base-aware keyed-assignment conflicts (currently always empty).
        type_conflicts: Vec<TypeConflict>,
        /// Gap-conditioned trust.
        confidence: Confidence,
    },
    /// The Sev-0 guard: reconciliation could not reach the nominal residual
    /// tolerance (for example, a coexact residue) — escalated unresolved
    /// rather than reported as resolved.
    EscalatedUnresolved {
        /// The post-reconciliation norm that failed.
        post_norm: f64,
        /// The tolerance it failed against.
        tol: f64,
        /// Diagnostic squared-norm ratios from the fixed-iteration split; not
        /// a certified orthogonal energy partition.
        fractions: (f64, f64, f64),
        /// Gap-conditioned trust, retained on every numerical adjudication.
        confidence: Confidence,
    },
}

/// Complete resource envelope for one admitted numerical merge.
///
/// `repair` governs the shared Hodge operator schedule, scalar work, ambient
/// cost, and cancellation polling. The remaining fields bound the dense
/// spectral proxy and every retained merge payload before publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafMergeBudget {
    /// Hodge sweeps, operator applications, total work, and poll cadence.
    pub repair: SheafRepairBudget,
    /// Maximum deterministic Jacobi sweeps for the spectral-gap proxy.
    pub spectral_sweeps: usize,
    /// Maximum conservative simultaneously live scalar-slot envelope for the
    /// whole merge, including the dense flat Laplacian.
    pub max_scalar_slots: usize,
    /// Maximum requested logical bytes across the published outcome's vector
    /// and provenance payloads. Allocator-internal spare capacity is outside
    /// this portable counter.
    pub max_output_bytes: usize,
    /// Maximum localized interface cells retained in a candidate conflict.
    pub max_conflict_cells: usize,
    /// Maximum UTF-8 bytes retained across both caller-supplied parent labels.
    pub max_provenance_bytes: usize,
}

/// Enforced and measured consumption for one bounded merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SheafMergeUsage {
    /// Shared Hodge/operator/work/ambient accounting. Its admitted scalar and
    /// work fields cover the whole merge rather than only decomposition.
    pub execution: SheafRepairUsage,
    /// Jacobi sweeps completed before the spectral proxy converged or hit its
    /// deterministic sweep ceiling.
    pub spectral_sweeps_completed: usize,
    /// Conservative logical output-byte envelope admitted before work began.
    pub admitted_output_bytes: usize,
    /// Requested logical payload bytes charged for the published outcome;
    /// allocator-internal spare capacity is not measured.
    pub reserved_output_bytes: usize,
    /// Localized conflict cells retained in the published outcome.
    pub conflict_cells: usize,
    /// Parent-label UTF-8 bytes retained in the published outcome.
    pub provenance_bytes: usize,
}

/// A merge verdict published only after all bounded work and the final
/// cancellation checkpoint complete.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedMergeOutcome {
    /// The legacy-compatible numerical verdict payload.
    pub outcome: MergeOutcome,
    /// Exact caller-admitted envelope.
    pub budget: SheafMergeBudget,
    /// Enforced and measured resource consumption.
    pub usage: SheafMergeUsage,
}

/// Aggregate resource envelope for one seeded candidate-remainder diagnostic.
///
/// `merge` bounds each individual adjudication. The remaining fields bound the
/// complete trial set under one [`Cx`] and one accountant; in particular,
/// discarded per-trial verdict payloads are charged cumulatively rather than
/// being treated as free temporary allocations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateRateBudget {
    /// Per-trial numerical, spectral, and retained-verdict envelope.
    pub merge: SheafMergeBudget,
    /// Maximum number of seeded trials in this diagnostic transaction.
    pub max_trials: usize,
    /// Maximum edit-incidence and Hodge operator evaluations across every
    /// trial.
    pub max_operator_evaluations: usize,
    /// Maximum scalar work items across edit generation and every merge.
    pub max_work_items: usize,
    /// Maximum conservative simultaneously live scalar-slot envelope,
    /// including the base, both generated branches, and one merge.
    pub max_scalar_slots: usize,
    /// Maximum cumulative logical bytes requested by discarded per-trial
    /// verdict payloads. Allocator-internal spare capacity is excluded.
    pub max_total_output_bytes: usize,
}

/// Enforced and measured consumption for one bounded seeded diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateRateUsage {
    /// Aggregate Hodge/operator/work/ambient accounting for the full trial set.
    pub execution: SheafRepairUsage,
    /// Exact aggregate operator schedule admitted before edit generation.
    pub admitted_operator_evaluations: usize,
    /// Total Jacobi sweeps completed across all per-trial spectral proxies.
    pub spectral_sweeps_completed: usize,
    /// Trials completely adjudicated before final aggregate publication.
    pub trials_completed: usize,
    /// Completed trials classified with a non-empty candidate remainder.
    pub conflicts: usize,
    /// Conservative cumulative logical verdict-byte envelope admitted before
    /// edit allocation or numerical work began.
    pub admitted_output_bytes: usize,
    /// Actual cumulative logical verdict bytes requested across discarded
    /// per-trial outcomes.
    pub requested_output_bytes: usize,
}

/// A seeded conflict-rate observation published only after every trial and the
/// final cancellation checkpoint complete.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedCandidateRate {
    /// Candidate-remainder conflicts divided by completed trials.
    pub rate: f64,
    /// Seed that deterministically generated both edit streams.
    pub seed: u64,
    /// Exact caller-admitted aggregate envelope.
    pub budget: CandidateRateBudget,
    /// Enforced and measured aggregate consumption.
    pub usage: CandidateRateUsage,
}

/// Structured refusal from the admitted, cancellation-aware merge path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SheafMergeError {
    /// Shared incidence, Hodge, work, ambient, or cancellation refusal.
    Repair(SheafRepairError),
    /// A merge-specific budget field that must be positive was zero.
    InvalidBudget {
        /// Stable field name.
        field: &'static str,
    },
    /// One numerical threshold violated its finite/range contract.
    InvalidThreshold {
        /// Stable field name.
        field: &'static str,
    },
    /// Weight cardinality differs from retained edge cardinality.
    WeightLength {
        /// Required number of weights.
        expected: usize,
        /// Supplied number of weights.
        actual: usize,
    },
    /// One edge weight is non-finite or negative.
    InvalidWeight {
        /// First invalid caller-order weight.
        index: usize,
    },
    /// The v1 entry point cannot adjudicate keyed assignments without a base
    /// assignment map.
    AssignmentsUnsupported,
    /// The conservative whole-merge scalar envelope exceeds its cap.
    MemoryBudgetExceeded {
        /// Required scalar slots.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// A retained merge payload exceeds its cardinality or byte cap.
    OutputBudgetExceeded {
        /// Stable resource name.
        resource: &'static str,
        /// Required cardinality or bytes.
        required: u128,
        /// Caller-admitted ceiling.
        cap: usize,
    },
    /// Checked merge-admission arithmetic exceeded `u128` or `usize`.
    BudgetArithmeticOverflow {
        /// Stable preflight stage.
        stage: &'static str,
    },
    /// Finite merge or spectral arithmetic overflowed.
    NumericalOverflow {
        /// Stable arithmetic stage.
        stage: &'static str,
    },
}

impl From<SheafRepairError> for SheafMergeError {
    fn from(source: SheafRepairError) -> Self {
        Self::Repair(source)
    }
}

impl From<SheafSkeletonError> for SheafMergeError {
    fn from(source: SheafSkeletonError) -> Self {
        Self::Repair(source.into())
    }
}

impl core::fmt::Display for SheafMergeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Repair(source) => write!(f, "{source}"),
            Self::InvalidBudget { field } => {
                write!(f, "sheaf merge budget field {field} must be positive")
            }
            Self::InvalidThreshold { field } => {
                write!(f, "sheaf merge threshold {field} is invalid")
            }
            Self::WeightLength { expected, actual } => write!(
                f,
                "sheaf merge requires {expected} edge weights, got {actual}"
            ),
            Self::InvalidWeight { index } => {
                write!(
                    f,
                    "sheaf merge edge weight {index} must be finite and non-negative"
                )
            }
            Self::AssignmentsUnsupported => write!(
                f,
                "sheaf merge requires a base assignment map for keyed payloads"
            ),
            Self::MemoryBudgetExceeded { required, cap } => write!(
                f,
                "sheaf merge scalar envelope requires {required} slots above cap {cap}"
            ),
            Self::OutputBudgetExceeded {
                resource,
                required,
                cap,
            } => write!(
                f,
                "sheaf merge output {resource} requires {required} above cap {cap}"
            ),
            Self::BudgetArithmeticOverflow { stage } => {
                write!(f, "sheaf merge budget arithmetic overflowed during {stage}")
            }
            Self::NumericalOverflow { stage } => {
                write!(f, "sheaf merge arithmetic overflowed during {stage}")
            }
        }
    }
}

impl std::error::Error for SheafMergeError {}

/// Invalid or refused input to the seeded candidate-remainder diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateRateError {
    /// A rate requires a non-empty trial set.
    ZeroTrials,
    /// Edit magnitudes must be finite and non-negative.
    InvalidEditScale,
    /// The supplied skeleton cannot support the incidence operators.
    MalformedSkeleton,
    /// The requested trial count exceeds the aggregate admitted ceiling.
    TrialBudgetExceeded {
        /// Requested seeded trials.
        requested: usize,
        /// Admitted trial ceiling.
        cap: usize,
    },
    /// A generated merge refused rather than producing an adjudication.
    TrialRefused {
        /// Stable refusal reason from [`three_way_merge`].
        reason: &'static str,
    },
    /// Structured refusal from the bounded aggregate path.
    Bounded(SheafMergeError),
}

impl From<SheafMergeError> for CandidateRateError {
    fn from(source: SheafMergeError) -> Self {
        Self::Bounded(source)
    }
}

impl From<SheafRepairError> for CandidateRateError {
    fn from(source: SheafRepairError) -> Self {
        Self::Bounded(source.into())
    }
}

fn norm_inf(v: &[f64]) -> f64 {
    if v.iter().any(|value| !value.is_finite()) {
        f64::INFINITY
    } else {
        v.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
    }
}

fn skeleton_is_well_formed(skeleton: &SheafSkeleton) -> bool {
    if skeleton.n_patches == 0 {
        return false;
    }
    let mut edges = BTreeSet::new();
    for &(u, v) in &skeleton.edges {
        if u >= v || v >= skeleton.n_patches || !edges.insert((u, v)) {
            return false;
        }
    }
    let mut triangles = BTreeSet::new();
    for &(a, b, c) in &skeleton.triangles {
        if a >= b
            || b >= c
            || c >= skeleton.n_patches
            || !triangles.insert((a, b, c))
            || !edges.contains(&(a, b))
            || !edges.contains(&(b, c))
            || !edges.contains(&(a, c))
        {
            return false;
        }
    }
    true
}

/// Dense symmetric eigenvalues by cyclic Jacobi (small matrices —
/// patch counts, not DOF counts).
#[allow(clippy::needless_range_loop)] // rotations touch (k,p),(k,q) pairs
fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    for _sweep in 0..64 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p][q] * a[p][q];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut ev: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    ev.sort_by(f64::total_cmp);
    ev
}

/// The spectral gap (algebraic connectivity λ₂) of the weighted
/// patch-adjacency Laplacian — the Proposal-5 trust signal. Weights
/// default to 1 (e.g. interface sample counts belong here).
#[must_use]
pub fn spectral_gap(skeleton: &SheafSkeleton, weights: Option<&[f64]>) -> f64 {
    if !skeleton_is_well_formed(skeleton)
        || weights.is_some_and(|values| {
            values.len() != skeleton.edges.len()
                || values
                    .iter()
                    .any(|weight| !weight.is_finite() || *weight < 0.0)
        })
    {
        return f64::NAN;
    }
    let n = skeleton.n_patches;
    let mut lap = vec![vec![0.0f64; n]; n];
    for (k, &(u, v)) in skeleton.edges.iter().enumerate() {
        let w = weights.map_or(1.0, |ws| ws[k]);
        lap[u][u] += w;
        lap[v][v] += w;
        lap[u][v] -= w;
        lap[v][u] -= w;
    }
    if lap.iter().flatten().any(|entry| !entry.is_finite()) {
        return f64::NAN;
    }
    let ev = jacobi_eigenvalues(lap);
    // Algebraic connectivity is the SECOND sorted eigenvalue. Searching for
    // the first positive eigenvalue would falsely call a disconnected graph
    // well connected by skipping its additional zero modes.
    let Some(&lambda_2) = ev.get(1) else {
        return 0.0;
    };
    if !lambda_2.is_finite() {
        lambda_2
    } else if lambda_2 <= 1e-9 {
        0.0
    } else {
        lambda_2
    }
}

fn merge_add(left: u128, right: u128, stage: &'static str) -> Result<u128, SheafMergeError> {
    left.checked_add(right)
        .ok_or(SheafMergeError::BudgetArithmeticOverflow { stage })
}

fn merge_mul(left: u128, right: u128, stage: &'static str) -> Result<u128, SheafMergeError> {
    left.checked_mul(right)
        .ok_or(SheafMergeError::BudgetArithmeticOverflow { stage })
}

#[derive(Clone, Copy)]
struct MergeAdmission {
    execution: RepairAdmission,
    output_bytes: usize,
}

fn checked_merge_scalar_envelope(
    skeleton: &AdmittedSheafSkeleton,
    hodge_slots: usize,
) -> Result<u128, SheafMergeError> {
    let patches = skeleton.n_patches() as u128;
    let edges = skeleton.edges().len() as u128;
    let spectral = merge_mul(patches, patches, "spectral-scalar-envelope")?;
    let numerical = merge_add(
        hodge_slots as u128,
        merge_mul(edges, 2, "merge-live-edge-scalars")?,
        "merge-numerical-scalar-envelope",
    )?;
    Ok(spectral.max(numerical).max(edges))
}

fn checked_merge_output_envelope(
    skeleton: &AdmittedSheafSkeleton,
    budget: SheafMergeBudget,
) -> Result<u128, SheafMergeError> {
    let scalar_bytes = core::mem::size_of::<f64>() as u128;
    let edges = skeleton.edges().len() as u128;
    let patches = skeleton.n_patches() as u128;
    let trivial = merge_mul(edges, scalar_bytes, "trivial-output-bytes")?;
    let resolved = merge_mul(
        merge_add(edges, patches, "resolved-output-scalars")?,
        scalar_bytes,
        "resolved-output-bytes",
    )?;
    let conflict_cells = merge_mul(
        budget.max_conflict_cells as u128,
        core::mem::size_of::<((usize, usize), f64)>() as u128,
        "conflict-cell-output-bytes",
    )?;
    let conflict = merge_add(
        merge_add(
            conflict_cells,
            budget.max_provenance_bytes as u128,
            "conflict-provenance-output-bytes",
        )?,
        core::mem::size_of::<CandidateRemainderConflict>() as u128,
        "conflict-container-output-bytes",
    )?;
    Ok(trivial.max(resolved).max(conflict))
}

fn checked_merge_work_envelope(
    skeleton: &AdmittedSheafSkeleton,
    hodge_work: usize,
    budget: SheafMergeBudget,
    has_weights: bool,
) -> Result<u128, SheafMergeError> {
    let patches = skeleton.n_patches() as u128;
    let edges = skeleton.edges().len() as u128;
    let weight_pass = if has_weights { 1 } else { 0 };
    let validation = merge_mul(edges, 3 + weight_pass, "merge-validation-work")?;
    let equality = merge_mul(edges, 3, "merge-equality-work")?;
    let matrix = merge_mul(patches, patches, "spectral-matrix-work")?;
    let pairs = merge_mul(patches, patches.saturating_sub(1), "spectral-pair-work")? / 2;
    let per_pair = merge_add(
        2,
        merge_mul(patches, 2, "spectral-rotation-span")?,
        "spectral-pair-span",
    )?;
    let spectral_sweeps = merge_mul(
        merge_mul(pairs, per_pair, "spectral-sweep-work")?,
        budget.spectral_sweeps as u128,
        "spectral-total-work",
    )?;
    let spectral = merge_add(
        merge_add(matrix, edges, "spectral-assembly-work")?,
        merge_add(spectral_sweeps, patches, "spectral-selection-work")?,
        "spectral-work",
    )?;
    let union_and_gauge = merge_mul(edges, 4, "merge-union-gauge-work")?;
    let residual_norms = merge_mul(edges, 3, "merge-residual-norm-work")?;
    let conflict_scan = edges;
    let retained_cells = (budget.max_conflict_cells.min(skeleton.edges().len())) as u128;
    let conflict_sort = merge_mul(
        retained_cells,
        retained_cells.saturating_sub(1),
        "merge-conflict-sort-work",
    )? / 2;
    let conflict_output = merge_add(
        merge_add(conflict_scan, edges, "merge-conflict-retain-work")?,
        merge_add(
            conflict_sort,
            budget.max_provenance_bytes as u128,
            "merge-conflict-provenance-work",
        )?,
        "merge-conflict-output-work",
    )?;
    let nontrivial = merge_add(
        merge_add(
            merge_add(spectral, union_and_gauge, "merge-numerical-work")?,
            hodge_work as u128,
            "merge-hodge-work",
        )?,
        merge_add(residual_norms, conflict_output, "merge-post-hodge-work")?,
        "merge-nontrivial-work",
    )?;
    let prefix = merge_add(validation, equality, "merge-prefix-work")?;
    let trivial_copy = edges;
    merge_add(prefix, nontrivial.max(trivial_copy), "merge-work-envelope")
}

fn admit_merge_budget(
    skeleton: &AdmittedSheafSkeleton,
    budget: SheafMergeBudget,
    has_weights: bool,
) -> Result<MergeAdmission, SheafMergeError> {
    if budget.spectral_sweeps == 0 {
        return Err(SheafMergeError::InvalidBudget {
            field: "spectral_sweeps",
        });
    }
    let hodge = admit_repair_budget(skeleton, budget.repair)?;
    let scalar_slots = checked_merge_scalar_envelope(skeleton, hodge.scalar_slots)?;
    if scalar_slots > budget.max_scalar_slots as u128 {
        return Err(SheafMergeError::MemoryBudgetExceeded {
            required: scalar_slots,
            cap: budget.max_scalar_slots,
        });
    }
    let work_items = checked_merge_work_envelope(skeleton, hodge.work_items, budget, has_weights)?;
    if work_items > budget.repair.max_work_items as u128 {
        return Err(SheafRepairError::WorkItemBudgetExceeded {
            stage: "merge-work-preflight",
            required: work_items,
            cap: budget.repair.max_work_items,
        }
        .into());
    }
    let output_bytes = checked_merge_output_envelope(skeleton, budget)?;
    if output_bytes > budget.max_output_bytes as u128 {
        return Err(SheafMergeError::OutputBudgetExceeded {
            resource: "output-bytes",
            required: output_bytes,
            cap: budget.max_output_bytes,
        });
    }
    Ok(MergeAdmission {
        execution: RepairAdmission {
            scalar_slots: usize::try_from(scalar_slots).map_err(|_| {
                SheafMergeError::BudgetArithmeticOverflow {
                    stage: "merge-scalar-publication",
                }
            })?,
            operator_evaluations: hodge.operator_evaluations,
            work_items: usize::try_from(work_items).map_err(|_| {
                SheafMergeError::BudgetArithmeticOverflow {
                    stage: "merge-work-publication",
                }
            })?,
        },
        output_bytes: usize::try_from(output_bytes).map_err(|_| {
            SheafMergeError::BudgetArithmeticOverflow {
                stage: "merge-output-publication",
            }
        })?,
    })
}

#[allow(clippy::too_many_lines)] // mirrors the legacy cyclic-Jacobi loop under one accountant
fn bounded_spectral_gap(
    skeleton: &AdmittedSheafSkeleton,
    weights: Option<&[f64]>,
    sweeps: usize,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(f64, usize), SheafMergeError> {
    let n = skeleton.n_patches();
    let matrix_len = n
        .checked_mul(n)
        .ok_or(SheafMergeError::BudgetArithmeticOverflow {
            stage: "spectral-matrix-length",
        })?;
    let mut lap = zeroed_output_bounded(matrix_len, "spectral-matrix-allocation", accountant)?;
    for (edge, &(u, v)) in skeleton.edges().iter().enumerate() {
        accountant.consume_item("spectral-assembly")?;
        let weight = weights.map_or(1.0, |values| values[edge]);
        let uu = u * n + u;
        let vv = v * n + v;
        let uv = u * n + v;
        let vu = v * n + u;
        lap[uu] += weight;
        lap[vv] += weight;
        lap[uv] -= weight;
        lap[vu] -= weight;
        if [lap[uu], lap[vv], lap[uv], lap[vu]]
            .into_iter()
            .any(|entry| !entry.is_finite())
        {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "spectral-assembly",
            });
        }
    }

    let mut completed_sweeps = 0usize;
    for _ in 0..sweeps {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                accountant.consume_item("spectral-off-diagonal")?;
                let value = lap[p * n + q];
                // Match legacy Jacobi: a finite entry may square to +∞. The
                // sum is only a convergence sentinel, so +∞ means "continue
                // rotating" and is not itself a failed numerical result.
                off += value * value;
            }
        }
        completed_sweeps =
            completed_sweeps
                .checked_add(1)
                .ok_or(SheafMergeError::BudgetArithmeticOverflow {
                    stage: "spectral-sweeps-completed",
                })?;
        if off < 1e-24 {
            accountant.checkpoint("spectral-sweep")?;
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                accountant.consume_item("spectral-rotation")?;
                let pq = p * n + q;
                if lap[pq].abs() < 1e-300 {
                    continue;
                }
                let theta = (lap[q * n + q] - lap[p * n + p]) / (2.0 * lap[pq]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // Preserve the legacy finite-result semantics: `theta²` may
                // benignly overflow to +∞, yielding t=0, c=1, s=0. Reject
                // only when the applied rotation itself is non-finite.
                if [t, c, s].into_iter().any(|value| !value.is_finite()) {
                    return Err(SheafMergeError::NumericalOverflow {
                        stage: "spectral-rotation",
                    });
                }
                for k in 0..n {
                    accountant.consume_item("spectral-column-rotation")?;
                    let kp = k * n + p;
                    let kq = k * n + q;
                    let (akp, akq) = (lap[kp], lap[kq]);
                    lap[kp] = c * akp - s * akq;
                    lap[kq] = s * akp + c * akq;
                    if !(lap[kp].is_finite() && lap[kq].is_finite()) {
                        return Err(SheafMergeError::NumericalOverflow {
                            stage: "spectral-column-rotation",
                        });
                    }
                }
                for k in 0..n {
                    accountant.consume_item("spectral-row-rotation")?;
                    let pk = p * n + k;
                    let qk = q * n + k;
                    let (apk, aqk) = (lap[pk], lap[qk]);
                    lap[pk] = c * apk - s * aqk;
                    lap[qk] = s * apk + c * aqk;
                    if !(lap[pk].is_finite() && lap[qk].is_finite()) {
                        return Err(SheafMergeError::NumericalOverflow {
                            stage: "spectral-row-rotation",
                        });
                    }
                }
            }
        }
        accountant.checkpoint("spectral-sweep")?;
    }

    let mut smallest = None::<f64>;
    let mut second = None::<f64>;
    for index in 0..n {
        accountant.consume_item("spectral-eigenvalue-selection")?;
        let value = lap[index * n + index];
        if !value.is_finite() {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "spectral-eigenvalue-selection",
            });
        }
        if smallest.is_none_or(|first| value.total_cmp(&first).is_lt()) {
            second = smallest;
            smallest = Some(value);
        } else if second.is_none_or(|current| value.total_cmp(&current).is_lt()) {
            second = Some(value);
        }
    }
    let lambda_2 = second.unwrap_or(0.0);
    let gap = if lambda_2 <= 1e-9 { 0.0 } else { lambda_2 };
    Ok((gap, completed_sweeps))
}

/// Detect pairwise keyed-assignment differences. Without a base assignment
/// map these are only conflict candidates, not authoritative three-way
/// conflicts (either branch may equal the base). Coupling-graph legality of a
/// future merged assignment set is fs-iface's contract at its own layer.
#[must_use]
pub fn type_conflicts(x: &BranchState, y: &BranchState) -> Vec<TypeConflict> {
    let mut out = Vec::new();
    for (k, xv) in &x.assignments {
        if let Some(yv) = y.assignments.get(k)
            && xv != yv
        {
            out.push(TypeConflict {
                key: k.clone(),
                x_value: xv.clone(),
                y_value: yv.clone(),
            });
        }
    }
    out
}

/// The three-way merge. `base` is the common ancestor's mismatch
/// cochain; `tol` is the nominal residual tolerance the reconciled state
/// must satisfy to be reported resolved; `gap_threshold` conditions trust.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn three_way_merge(
    skeleton: &SheafSkeleton,
    base: &[f64],
    x: &BranchState,
    y: &BranchState,
    weights: Option<&[f64]>,
    tol: f64,
    gap_threshold: f64,
) -> MergeOutcome {
    if !skeleton_is_well_formed(skeleton) {
        return MergeOutcome::Refused {
            reason: "malformed sheaf skeleton",
        };
    }
    if base.len() != skeleton.edges.len()
        || x.mismatch.len() != base.len()
        || y.mismatch.len() != base.len()
    {
        return MergeOutcome::Refused {
            reason: "cochain length mismatch",
        };
    }
    if !tol.is_finite() || tol < 0.0 || !gap_threshold.is_finite() || gap_threshold <= 0.0 {
        return MergeOutcome::Refused {
            reason: "residual tolerance must be finite non-negative and gap threshold finite positive",
        };
    }
    if base
        .iter()
        .chain(&x.mismatch)
        .chain(&y.mismatch)
        .any(|value| !value.is_finite())
    {
        return MergeOutcome::Refused {
            reason: "cochains must be finite",
        };
    }
    if weights.is_some_and(|values| {
        values.len() != skeleton.edges.len()
            || values
                .iter()
                .any(|weight| !weight.is_finite() || *weight < 0.0)
    }) {
        return MergeOutcome::Refused {
            reason: "weights must match edges and be finite non-negative values",
        };
    }
    // This v1 API carries only the base mismatch cochain, not the base
    // assignment map, and none of its success variants can return merged
    // assignments. Even pairwise-different branch values do not prove a
    // three-way conflict because either may equal the unknown base. Refuse all
    // assignment payloads rather than inventing conflict authority or dropping
    // state.
    if !x.assignments.is_empty() || !y.assignments.is_empty() {
        return MergeOutcome::Refused {
            reason: "base-aware assignment merge is not represented",
        };
    }
    // Boundary fast paths: no decomposition, no false ceremony.
    // Compare in place so an O(1) result does not allocate whole bit-vectors.
    let same_bits = |left: &[f64], right: &[f64]| {
        left.len() == right.len()
            && left
                .iter()
                .zip(right)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    };
    if same_bits(&x.mismatch, &y.mismatch) {
        return MergeOutcome::Trivial {
            reason: "branches identical",
            merged: x.mismatch.clone(),
        };
    }
    if same_bits(&x.mismatch, base) {
        return MergeOutcome::Trivial {
            reason: "X unchanged from base",
            merged: y.mismatch.clone(),
        };
    }
    if same_bits(&y.mismatch, base) {
        return MergeOutcome::Trivial {
            reason: "Y unchanged from base",
            merged: x.mismatch.clone(),
        };
    }
    // Nontrivial numerical adjudication is the first path that consumes gap
    // confidence; trivial equality must not pay for the dense Jacobi proxy.
    let gap = spectral_gap(skeleton, weights);
    if !gap.is_finite() {
        return MergeOutcome::Refused {
            reason: "spectral-gap arithmetic is non-finite",
        };
    }
    let confidence = if gap < gap_threshold {
        Confidence::LowGap {
            gap,
            threshold: gap_threshold,
        }
    } else {
        Confidence::Normal { gap }
    };
    // The naive union of edits at the cochain level: X + Y − B.
    let union: Vec<f64> = x
        .mismatch
        .iter()
        .zip(&y.mismatch)
        .zip(base)
        .map(|((a, b), c)| a + b - c)
        .collect();
    if union.iter().any(|value| !value.is_finite()) {
        return MergeOutcome::Refused {
            reason: "merged cochain arithmetic is non-finite",
        };
    }
    let split = match hodge_decompose(skeleton, &union) {
        Ok(split) => split,
        Err(
            SheafRepairError::NumericalOverflow { .. }
            | SheafRepairError::Skeleton(SheafSkeletonError::NumericalOverflow { .. }),
        ) => {
            return MergeOutcome::Refused {
                reason: "decomposition arithmetic is non-finite",
            };
        }
        Err(_) => {
            return MergeOutcome::Refused {
                reason: "decomposition refused malformed, non-finite, or exhausted input",
            };
        }
    };
    if split
        .exact
        .iter()
        .chain(&split.potential)
        .chain(&split.coexact)
        .chain(&split.harmonic)
        .any(|value| !value.is_finite())
        || [split.fractions.0, split.fractions.1, split.fractions.2]
            .into_iter()
            .any(|value| !value.is_finite())
    {
        return MergeOutcome::Refused {
            reason: "decomposition arithmetic is non-finite",
        };
    }
    // Deterministic fixed-iteration gauge reconciliation. Resolution is
    // authorized only by the nominal residual check below, not by a
    // cohomology classification.
    let merged = match apply_gauge(skeleton, &union, &split.potential) {
        Ok(merged) => merged,
        Err(
            SheafRepairError::NumericalOverflow { .. }
            | SheafRepairError::Skeleton(SheafSkeletonError::NumericalOverflow { .. }),
        ) => {
            return MergeOutcome::Refused {
                reason: "post-gauge arithmetic is non-finite",
            };
        }
        Err(_) => {
            return MergeOutcome::Refused {
                reason: "post-gauge construction refused malformed, non-finite, or exhausted input",
            };
        }
    };
    if merged.iter().any(|value| !value.is_finite()) {
        return MergeOutcome::Refused {
            reason: "post-gauge arithmetic is non-finite",
        };
    }
    // Sev-0 post-state check: nominal values only. This is not the interval
    // watertightness certificate from the base sheaf path.
    let post_norm = norm_inf(&merged);
    if post_norm <= tol {
        return MergeOutcome::Resolved {
            merged,
            gauge: split.potential,
            residual_receipt: MergeResidualReceipt { post_norm, tol },
            confidence,
        };
    }
    // The nominal residual check failed. A dominant fixed-iteration remainder
    // is retained as an operational merge-conflict candidate only; anything
    // else (for example a coexact circulation) escalates unresolved.
    let harmonic_norm = norm_inf(&split.harmonic);
    let coexact_norm = norm_inf(&split.coexact);
    if harmonic_norm > tol && harmonic_norm >= coexact_norm {
        let mut cells: Vec<((usize, usize), f64)> = skeleton
            .edges
            .iter()
            .zip(&split.harmonic)
            .filter(|(_, h)| h.abs() > tol)
            .map(|(&e, &h)| (e, h.abs()))
            .collect();
        cells.sort_by(|a, b| b.1.total_cmp(&a.1));
        return MergeOutcome::Conflicted {
            candidate_remainders: vec![CandidateRemainderConflict {
                cells,
                parents: (x.provenance.clone(), y.provenance.clone()),
            }],
            type_conflicts: Vec::new(),
            confidence,
        };
    }
    MergeOutcome::EscalatedUnresolved {
        post_norm,
        tol,
        fractions: split.fractions,
        confidence,
    }
}

fn retained_bytes<T>(len: usize, stage: &'static str) -> Result<usize, SheafMergeError> {
    len.checked_mul(core::mem::size_of::<T>())
        .ok_or(SheafMergeError::BudgetArithmeticOverflow { stage })
}

fn output_vec_with_capacity<T>(
    len: usize,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<T>, SheafMergeError> {
    let bytes = retained_bytes::<T>(len, stage)?;
    accountant.reserve_plan_bytes(stage, bytes)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(len)
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    accountant.checkpoint(stage)?;
    Ok(output)
}

fn retain_existing_output<T>(
    len: usize,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<(), SheafMergeError> {
    accountant.reserve_plan_bytes(stage, retained_bytes::<T>(len, stage)?)?;
    Ok(())
}

fn accounted_same_bits(
    left: &[f64],
    right: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<bool, SheafMergeError> {
    if left.len() != right.len() {
        return Ok(false);
    }
    for (a, b) in left.iter().zip(right) {
        accountant.consume_item("merge-equality")?;
        if a.to_bits() != b.to_bits() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn accounted_output_copy(
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafMergeError> {
    let mut output = output_vec_with_capacity(values.len(), stage, accountant)?;
    for value in values {
        accountant.consume_item(stage)?;
        output.push(*value);
    }
    Ok(output)
}

fn bounded_merge_norm_inf(
    values: &[f64],
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<f64, SheafMergeError> {
    let mut largest = 0.0f64;
    for value in values {
        accountant.consume_item(stage)?;
        if !value.is_finite() {
            return Err(SheafMergeError::NumericalOverflow { stage });
        }
        largest = largest.max(value.abs());
    }
    Ok(largest)
}

fn bounded_merge_union(
    base: &[f64],
    x: &[f64],
    y: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafMergeError> {
    let mut union = zeroed_output_bounded(x.len(), "merge-union-allocation", accountant)?;
    for (((value, x_value), y_value), base_value) in union.iter_mut().zip(x).zip(y).zip(base) {
        accountant.consume_item("merge-union")?;
        *value = x_value + y_value - base_value;
        if !value.is_finite() {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "merge-union",
            });
        }
    }
    Ok(union)
}

fn bounded_merge_apply_gauge(
    skeleton: &AdmittedSheafSkeleton,
    mismatch: &[f64],
    gauge: &[f64],
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafMergeError> {
    if gauge.len() != skeleton.n_patches() {
        return Err(SheafSkeletonError::CochainLength {
            role: "gauge",
            expected: skeleton.n_patches(),
            actual: gauge.len(),
        }
        .into());
    }
    let mut merged =
        zeroed_output_bounded(skeleton.edges().len(), "merge-gauge-allocation", accountant)?;
    for (edge, &(u, v)) in skeleton.edges().iter().enumerate() {
        accountant.consume_item("merge-gauge-application")?;
        let correction = gauge[v] - gauge[u];
        merged[edge] = mismatch[edge] - correction;
        if !(correction.is_finite() && merged[edge].is_finite()) {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "merge-gauge-application",
            });
        }
    }
    Ok(merged)
}

fn accounted_parent(
    value: &str,
    stage: &'static str,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<String, SheafMergeError> {
    accountant.reserve_plan_bytes(stage, value.len())?;
    let mut output = String::new();
    output
        .try_reserve_exact(value.len())
        .map_err(|_| SheafSkeletonError::ResourceExhausted { stage })?;
    accountant.checkpoint(stage)?;
    for character in value.chars() {
        for _ in 0..character.len_utf8() {
            accountant.consume_item(stage)?;
        }
        output.push(character);
    }
    accountant.checkpoint(stage)?;
    Ok(output)
}

struct AccountedMergeOutcome {
    outcome: MergeOutcome,
    spectral_sweeps_completed: usize,
    conflict_cells: usize,
    provenance_bytes: usize,
}

fn finish_accounted_merge(
    outcome: MergeOutcome,
    spectral_sweeps_completed: usize,
    conflict_cells: usize,
    provenance_bytes: usize,
) -> AccountedMergeOutcome {
    AccountedMergeOutcome {
        outcome,
        spectral_sweeps_completed,
        conflict_cells,
        provenance_bytes,
    }
}

/// Run a three-way numerical merge over sealed incidence under one explicit
/// spectral, Hodge, work, memory, output, deadline, and cancellation envelope.
///
/// The returned [`MergeOutcome`] retains the legacy heuristic/no-claim
/// semantics. A refusal returns no partial verdict or retained output.
///
/// # Errors
/// Returns [`SheafMergeError`] for invalid thresholds, cochains, weights, or
/// unsupported assignments; insufficient work, memory, output, or ambient
/// budget; cancellation; allocation refusal; or finite arithmetic overflow.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn three_way_merge_bounded(
    skeleton: &AdmittedSheafSkeleton,
    base: &[f64],
    x: &BranchState,
    y: &BranchState,
    weights: Option<&[f64]>,
    tol: f64,
    gap_threshold: f64,
    budget: SheafMergeBudget,
    cx: &Cx<'_>,
) -> Result<BoundedMergeOutcome, SheafMergeError> {
    let edge_count = skeleton.edges().len();
    for (role, values) in [
        ("base", base),
        ("x-mismatch", x.mismatch.as_slice()),
        ("y-mismatch", y.mismatch.as_slice()),
    ] {
        if values.len() != edge_count {
            return Err(SheafSkeletonError::CochainLength {
                role,
                expected: edge_count,
                actual: values.len(),
            }
            .into());
        }
    }
    if !tol.is_finite() || tol < 0.0 {
        return Err(SheafMergeError::InvalidThreshold { field: "tol" });
    }
    if !gap_threshold.is_finite() || gap_threshold <= 0.0 {
        return Err(SheafMergeError::InvalidThreshold {
            field: "gap_threshold",
        });
    }
    let admission = admit_merge_budget(skeleton, budget, weights.is_some())?;
    let mut accountant = RepairAccountant::new(
        cx,
        budget.repair,
        planned_cost(admission.execution)?,
        budget.max_output_bytes,
        0,
    )?;
    let accounted = three_way_merge_accounted(
        skeleton,
        base,
        x,
        y,
        weights,
        tol,
        gap_threshold,
        budget,
        &mut accountant,
    )?;
    accountant.checkpoint("merge-publication")?;
    let usage = SheafMergeUsage {
        execution: accountant.usage(
            admission.execution.scalar_slots,
            admission.execution.work_items,
        ),
        spectral_sweeps_completed: accounted.spectral_sweeps_completed,
        admitted_output_bytes: admission.output_bytes,
        reserved_output_bytes: accountant.reserved_plan_bytes(),
        conflict_cells: accounted.conflict_cells,
        provenance_bytes: accounted.provenance_bytes,
    };
    Ok(BoundedMergeOutcome {
        outcome: accounted.outcome,
        budget,
        usage,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn three_way_merge_accounted(
    skeleton: &AdmittedSheafSkeleton,
    base: &[f64],
    x: &BranchState,
    y: &BranchState,
    weights: Option<&[f64]>,
    tol: f64,
    gap_threshold: f64,
    budget: SheafMergeBudget,
    mut accountant: &mut RepairAccountant<'_, '_>,
) -> Result<AccountedMergeOutcome, SheafMergeError> {
    let edge_count = skeleton.edges().len();
    accountant.checkpoint("merge-admission")?;
    validate_bounded_cochain(
        base,
        edge_count,
        "base",
        "merge-base-validation",
        &mut accountant,
    )?;
    validate_bounded_cochain(
        &x.mismatch,
        edge_count,
        "x-mismatch",
        "merge-x-validation",
        &mut accountant,
    )?;
    validate_bounded_cochain(
        &y.mismatch,
        edge_count,
        "y-mismatch",
        "merge-y-validation",
        &mut accountant,
    )?;
    if let Some(values) = weights {
        if values.len() != edge_count {
            return Err(SheafMergeError::WeightLength {
                expected: edge_count,
                actual: values.len(),
            });
        }
        for (index, weight) in values.iter().enumerate() {
            accountant.consume_item("merge-weight-validation")?;
            if !weight.is_finite() || *weight < 0.0 {
                return Err(SheafMergeError::InvalidWeight { index });
            }
        }
    }
    if !x.assignments.is_empty() || !y.assignments.is_empty() {
        return Err(SheafMergeError::AssignmentsUnsupported);
    }

    if accounted_same_bits(&x.mismatch, &y.mismatch, &mut accountant)? {
        let merged = accounted_output_copy(
            &x.mismatch,
            "merge-trivial-identical-output",
            &mut accountant,
        )?;
        return Ok(finish_accounted_merge(
            MergeOutcome::Trivial {
                reason: "branches identical",
                merged,
            },
            0,
            0,
            0,
        ));
    }
    if accounted_same_bits(&x.mismatch, base, &mut accountant)? {
        let merged = accounted_output_copy(&y.mismatch, "merge-trivial-x-output", &mut accountant)?;
        return Ok(finish_accounted_merge(
            MergeOutcome::Trivial {
                reason: "X unchanged from base",
                merged,
            },
            0,
            0,
            0,
        ));
    }
    if accounted_same_bits(&y.mismatch, base, &mut accountant)? {
        let merged = accounted_output_copy(&x.mismatch, "merge-trivial-y-output", &mut accountant)?;
        return Ok(finish_accounted_merge(
            MergeOutcome::Trivial {
                reason: "Y unchanged from base",
                merged,
            },
            0,
            0,
            0,
        ));
    }

    let (gap, spectral_sweeps_completed) =
        bounded_spectral_gap(skeleton, weights, budget.spectral_sweeps, &mut accountant)?;
    let confidence = if gap < gap_threshold {
        Confidence::LowGap {
            gap,
            threshold: gap_threshold,
        }
    } else {
        Confidence::Normal { gap }
    };
    let union = bounded_merge_union(base, &x.mismatch, &y.mismatch, &mut accountant)?;
    let split = hodge_decompose_accounted(skeleton, &union, &mut accountant)?;
    let merged = bounded_merge_apply_gauge(skeleton, &union, &split.potential, &mut accountant)?;
    let post_norm = bounded_merge_norm_inf(&merged, "merge-post-norm", &mut accountant)?;
    if post_norm <= tol {
        retain_existing_output::<f64>(
            merged.len(),
            "merge-resolved-cochain-output",
            &mut accountant,
        )?;
        retain_existing_output::<f64>(
            split.potential.len(),
            "merge-resolved-gauge-output",
            &mut accountant,
        )?;
        return Ok(finish_accounted_merge(
            MergeOutcome::Resolved {
                merged,
                gauge: split.potential,
                residual_receipt: MergeResidualReceipt { post_norm, tol },
                confidence,
            },
            spectral_sweeps_completed,
            0,
            0,
        ));
    }

    let harmonic_norm =
        bounded_merge_norm_inf(&split.harmonic, "merge-harmonic-norm", &mut accountant)?;
    let coexact_norm =
        bounded_merge_norm_inf(&split.coexact, "merge-coexact-norm", &mut accountant)?;
    if harmonic_norm > tol && harmonic_norm >= coexact_norm {
        let mut required_cells = 0usize;
        for value in &split.harmonic {
            accountant.consume_item("merge-conflict-count")?;
            if value.abs() > tol {
                required_cells = required_cells.checked_add(1).ok_or(
                    SheafMergeError::BudgetArithmeticOverflow {
                        stage: "merge-conflict-count",
                    },
                )?;
            }
        }
        if required_cells > budget.max_conflict_cells {
            return Err(SheafMergeError::OutputBudgetExceeded {
                resource: "conflict-cells",
                required: required_cells as u128,
                cap: budget.max_conflict_cells,
            });
        }
        let provenance_bytes = x.provenance.len().checked_add(y.provenance.len()).ok_or(
            SheafMergeError::BudgetArithmeticOverflow {
                stage: "merge-provenance-bytes",
            },
        )?;
        if provenance_bytes > budget.max_provenance_bytes {
            return Err(SheafMergeError::OutputBudgetExceeded {
                resource: "provenance-bytes",
                required: provenance_bytes as u128,
                cap: budget.max_provenance_bytes,
            });
        }
        let mut cells = output_vec_with_capacity::<((usize, usize), f64)>(
            required_cells,
            "merge-conflict-cell-output",
            &mut accountant,
        )?;
        for (&edge, &value) in skeleton.edges().iter().zip(&split.harmonic) {
            accountant.consume_item("merge-conflict-retain")?;
            if value.abs() > tol {
                cells.push((edge, value.abs()));
            }
        }
        for index in 1..cells.len() {
            let mut cursor = index;
            while cursor > 0 {
                accountant.consume_item("merge-conflict-sort")?;
                if cells[cursor - 1].1.total_cmp(&cells[cursor].1).is_ge() {
                    break;
                }
                cells.swap(cursor - 1, cursor);
                cursor -= 1;
            }
        }
        let parent_x = accounted_parent(&x.provenance, "merge-parent-x-output", &mut accountant)?;
        let parent_y = accounted_parent(&y.provenance, "merge-parent-y-output", &mut accountant)?;
        let mut candidate_remainders = output_vec_with_capacity::<CandidateRemainderConflict>(
            1,
            "merge-conflict-container-output",
            &mut accountant,
        )?;
        candidate_remainders.push(CandidateRemainderConflict {
            cells,
            parents: (parent_x, parent_y),
        });
        return Ok(finish_accounted_merge(
            MergeOutcome::Conflicted {
                candidate_remainders,
                type_conflicts: Vec::new(),
                confidence,
            },
            spectral_sweeps_completed,
            required_cells,
            provenance_bytes,
        ));
    }

    Ok(finish_accounted_merge(
        MergeOutcome::EscalatedUnresolved {
            post_norm,
            tol,
            fractions: split.fractions,
            confidence,
        },
        spectral_sweeps_completed,
        0,
        0,
    ))
}

/// Run seeded random three-way merges and measure the candidate-remainder
/// conflict rate.
///
/// This is one diagnostic input to Proposal 10's broader kill criterion. It
/// deliberately counts only [`CandidateRemainderConflict`] outcomes; callers
/// must separately retain escalations, refusals, and type conflicts before
/// claiming the full fraction of merges that could not auto-reconcile. It does
/// not count certified H¹ classes.
pub fn candidate_remainder_conflict_rate(
    skeleton: &SheafSkeleton,
    trials: usize,
    edit_scale: f64,
    seed: u64,
) -> Result<f64, CandidateRateError> {
    if trials == 0 {
        return Err(CandidateRateError::ZeroTrials);
    }
    if !edit_scale.is_finite() || edit_scale < 0.0 {
        return Err(CandidateRateError::InvalidEditScale);
    }
    if !skeleton_is_well_formed(skeleton) {
        return Err(CandidateRateError::MalformedSkeleton);
    }
    let mut state = seed;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let n_edges = skeleton.edges.len();
    let n_patches = skeleton.n_patches;
    let mut conflicts = 0usize;
    for _ in 0..trials {
        // Seeded synthetic edits: each branch re-gauges random patch values
        // (coboundary-style work), then receives small interface noise.
        let mut edit = |scale: f64| -> Result<Vec<f64>, CandidateRateError> {
            let offsets: Vec<f64> = (0..n_patches).map(|_| scale * lcg()).collect();
            let mut m = skeleton
                .d0(&offsets)
                .map_err(|_| CandidateRateError::TrialRefused {
                    reason: "seeded edit incidence refused",
                })?;
            // A small amount of independent interface noise.
            for v in &mut m {
                *v += 0.01 * scale * lcg();
                if !v.is_finite() {
                    return Err(CandidateRateError::TrialRefused {
                        reason: "seeded edit arithmetic is non-finite",
                    });
                }
            }
            Ok(m)
        };
        let base = vec![0.0f64; n_edges];
        let x = BranchState {
            provenance: "trial-x".to_string(),
            mismatch: edit(edit_scale)?,
            assignments: BTreeMap::new(),
        };
        let y = BranchState {
            provenance: "trial-y".to_string(),
            mismatch: edit(edit_scale)?,
            assignments: BTreeMap::new(),
        };
        let out = three_way_merge(skeleton, &base, &x, &y, None, 0.05 * edit_scale, 1e-6);
        match out {
            MergeOutcome::Conflicted {
                candidate_remainders,
                ..
            } if !candidate_remainders.is_empty() => conflicts += 1,
            MergeOutcome::Refused { reason } => {
                return Err(CandidateRateError::TrialRefused { reason });
            }
            _ => {}
        }
    }
    #[allow(clippy::cast_precision_loss)]
    {
        Ok(conflicts as f64 / trials as f64)
    }
}

#[derive(Clone, Copy)]
struct CandidateRateAdmission {
    execution: RepairAdmission,
    output_bytes: usize,
}

fn candidate_rate_usize(value: u128, stage: &'static str) -> Result<usize, SheafMergeError> {
    usize::try_from(value).map_err(|_| SheafMergeError::BudgetArithmeticOverflow { stage })
}

fn admit_candidate_rate_budget(
    skeleton: &AdmittedSheafSkeleton,
    trials: usize,
    budget: CandidateRateBudget,
) -> Result<CandidateRateAdmission, CandidateRateError> {
    if trials > budget.max_trials {
        return Err(CandidateRateError::TrialBudgetExceeded {
            requested: trials,
            cap: budget.max_trials,
        });
    }
    let per_trial = admit_merge_budget(skeleton, budget.merge, false)?;
    let trial_count = trials as u128;
    let patches = skeleton.n_patches() as u128;
    let edges = skeleton.edges().len() as u128;

    let per_trial_operators = merge_add(
        per_trial.execution.operator_evaluations as u128,
        2,
        "candidate-rate-edit-operators",
    )?;
    let operator_evaluations = merge_mul(
        per_trial_operators,
        trial_count,
        "candidate-rate-operator-envelope",
    )?;
    if operator_evaluations > budget.max_operator_evaluations as u128 {
        return Err(SheafRepairError::WorkBudgetExceeded {
            required: operator_evaluations,
            cap: budget.max_operator_evaluations,
        }
        .into());
    }

    // Each generated branch zero-initializes and fills P offsets (2P), applies
    // one bounded d0 (2E), then adds E noise values. Two branches therefore
    // consume 4P+6E work items and two real operator evaluations. One
    // additional item records the completed trial classification; the base is
    // initialized once for the whole transaction.
    let generated_offsets = merge_mul(patches, 4, "candidate-rate-offset-work")?;
    let generated_edges = merge_mul(edges, 6, "candidate-rate-edge-work")?;
    let generated_edits = merge_add(
        generated_offsets,
        generated_edges,
        "candidate-rate-edit-work",
    )?;
    let per_trial_work = merge_add(
        merge_add(
            per_trial.execution.work_items as u128,
            generated_edits,
            "candidate-rate-trial-work",
        )?,
        1,
        "candidate-rate-classification-work",
    )?;
    let work_items = merge_add(
        edges,
        merge_mul(
            per_trial_work,
            trial_count,
            "candidate-rate-total-trial-work",
        )?,
        "candidate-rate-work-envelope",
    )?;
    if work_items > budget.max_work_items as u128 {
        return Err(SheafRepairError::WorkItemBudgetExceeded {
            stage: "candidate-rate-work-preflight",
            required: work_items,
            cap: budget.max_work_items,
        }
        .into());
    }

    // The base and both completed branches remain live while one merge runs.
    // During Y generation its P offsets are additionally live, so take the
    // larger of that span and the admitted per-merge scalar envelope.
    let scalar_slots = merge_add(
        merge_mul(edges, 3, "candidate-rate-live-branches")?,
        patches.max(per_trial.execution.scalar_slots as u128),
        "candidate-rate-scalar-envelope",
    )?;
    if scalar_slots > budget.max_scalar_slots as u128 {
        return Err(SheafMergeError::MemoryBudgetExceeded {
            required: scalar_slots,
            cap: budget.max_scalar_slots,
        }
        .into());
    }

    // Trial verdicts are discarded after classification, but their requested
    // logical bytes remain cumulatively charged to bound allocation volume.
    let output_bytes = merge_mul(
        per_trial.output_bytes as u128,
        trial_count,
        "candidate-rate-output-envelope",
    )?;
    if output_bytes > budget.max_total_output_bytes as u128 {
        return Err(SheafMergeError::OutputBudgetExceeded {
            resource: "candidate-rate-output-bytes",
            required: output_bytes,
            cap: budget.max_total_output_bytes,
        }
        .into());
    }

    Ok(CandidateRateAdmission {
        execution: RepairAdmission {
            scalar_slots: candidate_rate_usize(scalar_slots, "candidate-rate-scalar-publication")?,
            operator_evaluations: candidate_rate_usize(
                operator_evaluations,
                "candidate-rate-operator-publication",
            )?,
            work_items: candidate_rate_usize(work_items, "candidate-rate-work-publication")?,
        },
        output_bytes: candidate_rate_usize(output_bytes, "candidate-rate-output-publication")?,
    })
}

fn candidate_rate_lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

fn bounded_candidate_edit(
    skeleton: &AdmittedSheafSkeleton,
    scale: f64,
    state: &mut u64,
    accountant: &mut RepairAccountant<'_, '_>,
) -> Result<Vec<f64>, SheafMergeError> {
    let mut offsets = zeroed_output_bounded(
        skeleton.n_patches(),
        "candidate-offset-allocation",
        accountant,
    )?;
    for offset in &mut offsets {
        accountant.consume_item("candidate-offset-generation")?;
        *offset = scale * candidate_rate_lcg(state);
        if !offset.is_finite() {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "candidate-offset-generation",
            });
        }
    }
    accountant.checkpoint("candidate-offset-generation")?;

    let mut mismatch = bounded_d0(skeleton, &offsets, "candidate-edit-d0", accountant)?;
    for value in &mut mismatch {
        accountant.consume_item("candidate-mismatch-generation")?;
        *value += 0.01 * scale * candidate_rate_lcg(state);
        if !value.is_finite() {
            return Err(SheafMergeError::NumericalOverflow {
                stage: "candidate-mismatch-generation",
            });
        }
    }
    accountant.checkpoint("candidate-mismatch-generation")?;
    Ok(mismatch)
}

/// Run the seeded candidate-remainder diagnostic over sealed incidence under
/// one aggregate work, memory, output, deadline, and cancellation transaction.
///
/// This preserves [`candidate_remainder_conflict_rate`]'s deterministic edit
/// stream and diagnostic-only numerator. Every trial shares one accountant;
/// no independently admitted per-trial `Cx` can reset aggregate authority.
/// Discarded verdict payloads remain cumulatively charged, and no rate or
/// partial histogram is published before the final cancellation checkpoint.
///
/// # Errors
/// Returns [`CandidateRateError`] for an empty or oversized trial set, invalid
/// edit scale, insufficient aggregate or per-merge resources, cancellation,
/// deadline/ambient refusal, allocation refusal, or finite arithmetic overflow.
#[allow(clippy::too_many_lines)]
pub fn candidate_remainder_conflict_rate_bounded(
    skeleton: &AdmittedSheafSkeleton,
    trials: usize,
    edit_scale: f64,
    seed: u64,
    budget: CandidateRateBudget,
    cx: &Cx<'_>,
) -> Result<BoundedCandidateRate, CandidateRateError> {
    if trials == 0 {
        return Err(CandidateRateError::ZeroTrials);
    }
    if !edit_scale.is_finite() || edit_scale < 0.0 {
        return Err(CandidateRateError::InvalidEditScale);
    }
    let admission = admit_candidate_rate_budget(skeleton, trials, budget)?;
    let execution_budget = SheafRepairBudget {
        sweeps: budget.merge.repair.sweeps,
        max_operator_evaluations: budget.max_operator_evaluations,
        max_work_items: budget.max_work_items,
        max_scalar_slots: budget.max_scalar_slots,
        poll_stride: budget.merge.repair.poll_stride,
    };
    let mut accountant = RepairAccountant::new(
        cx,
        execution_budget,
        planned_cost(admission.execution)?,
        budget.merge.max_output_bytes,
        0,
    )?;
    accountant.checkpoint("candidate-rate-admission")?;

    let base = zeroed_output_bounded(
        skeleton.edges().len(),
        "candidate-base-allocation",
        &mut accountant,
    )?;
    let mut state = seed;
    let mut conflicts = 0usize;
    let mut trials_completed = 0usize;
    let mut spectral_sweeps_completed = 0usize;
    let mut requested_output_bytes = 0usize;
    for _ in 0..trials {
        let output_bytes_before = accountant.reserved_plan_bytes();
        let x = BranchState {
            // Parent labels cannot affect the numerical classification and no
            // per-trial verdict escapes this aggregate diagnostic.
            provenance: String::new(),
            mismatch: bounded_candidate_edit(skeleton, edit_scale, &mut state, &mut accountant)?,
            assignments: BTreeMap::new(),
        };
        let y = BranchState {
            provenance: String::new(),
            mismatch: bounded_candidate_edit(skeleton, edit_scale, &mut state, &mut accountant)?,
            assignments: BTreeMap::new(),
        };
        let trial = three_way_merge_accounted(
            skeleton,
            &base,
            &x,
            &y,
            None,
            0.05 * edit_scale,
            1e-6,
            budget.merge,
            &mut accountant,
        )?;
        spectral_sweeps_completed = spectral_sweeps_completed
            .checked_add(trial.spectral_sweeps_completed)
            .ok_or(SheafMergeError::BudgetArithmeticOverflow {
                stage: "candidate-rate-spectral-sweeps",
            })?;
        accountant.consume_item("candidate-rate-classification")?;
        let is_conflict = matches!(
            &trial.outcome,
            MergeOutcome::Conflicted {
                candidate_remainders,
                ..
            } if !candidate_remainders.is_empty()
        );
        let output_bytes_after = accountant.reserved_plan_bytes();
        let trial_output_bytes = output_bytes_after.checked_sub(output_bytes_before).ok_or(
            SheafMergeError::BudgetArithmeticOverflow {
                stage: "candidate-rate-trial-output-delta",
            },
        )?;
        requested_output_bytes = requested_output_bytes
            .checked_add(trial_output_bytes)
            .ok_or(SheafMergeError::BudgetArithmeticOverflow {
                stage: "candidate-rate-requested-output",
            })?;
        if requested_output_bytes > budget.max_total_output_bytes {
            return Err(SheafMergeError::OutputBudgetExceeded {
                resource: "candidate-rate-output-bytes",
                required: requested_output_bytes as u128,
                cap: budget.max_total_output_bytes,
            }
            .into());
        }
        drop(trial);
        accountant.release_plan_bytes("candidate-trial-output-release", trial_output_bytes)?;
        if accountant.reserved_plan_bytes() != output_bytes_before {
            return Err(SheafMergeError::BudgetArithmeticOverflow {
                stage: "candidate-rate-output-release",
            }
            .into());
        }
        if is_conflict {
            conflicts =
                conflicts
                    .checked_add(1)
                    .ok_or(SheafMergeError::BudgetArithmeticOverflow {
                        stage: "candidate-rate-conflicts",
                    })?;
        }
        trials_completed =
            trials_completed
                .checked_add(1)
                .ok_or(SheafMergeError::BudgetArithmeticOverflow {
                    stage: "candidate-rate-trials-completed",
                })?;
        accountant.checkpoint("candidate-trial-complete")?;
    }
    accountant.checkpoint("candidate-rate-publication")?;

    #[allow(clippy::cast_precision_loss)]
    let rate = conflicts as f64 / trials as f64;
    let usage = CandidateRateUsage {
        execution: accountant.usage(
            admission.execution.scalar_slots,
            admission.execution.work_items,
        ),
        admitted_operator_evaluations: admission.execution.operator_evaluations,
        spectral_sweeps_completed,
        trials_completed,
        conflicts,
        admitted_output_bytes: admission.output_bytes,
        requested_output_bytes,
    };
    Ok(BoundedCandidateRate {
        rate,
        seed,
        budget,
        usage,
    })
}
