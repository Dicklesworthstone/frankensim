//! fs-truss-e2e — TrussPath: a deterministic truss iterate with an advisory,
//! endpoint-checked load path. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! A structural optimizer returns member sizes. This returns a deterministic
//! truss iterate plus explicit numerical diagnostics and an advisory account of
//! how the load travels through it, composing crates never designed to meet:
//!
//! - **Ground-structure optimization** ([`fs_truss`]): a Michell ground
//!   structure (all admissible candidate bars) is iterated toward minimum
//!   volume and equilibrium by a first-order PDHG solver. A separate outward
//!   certificate repairs the primal to exact feasibility and independently
//!   scales and checks the dual before publishing optimum bounds.
//! - **The critical load path** ([`fs_tropical`]): the active bars form a
//!   directed acyclic graph oriented by distance-to-support; a MAX-PLUS
//!   (tropical) critical-path computation finds a connected chain of active
//!   bars from the load node to a support, and names a bottleneck only when the
//!   rounded task graph has a unique heaviest chain.
//! - **Honest colors** ([`fs_evidence`]): optimality becomes `Verified` only
//!   from the retained outward primal/dual certificate. A load path becomes
//!   `Verified` only when the same solver receipt supplies separated member-
//!   force intervals, every active-set decision and support-ward orientation
//!   separates, and interval path/bottleneck comparisons are unique. Otherwise
//!   the rounded advisory path remains `Estimated`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::{Color, ProvenanceHash};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_tropical::{MAX_TASK_DAG_EDGES, MAX_TASK_DAG_NODES, TaskDag, TropicalError};
use fs_truss::{
    GroundLimits, GroundRules, GroundStructure, LayoutCase, LayoutCertificateError,
    LayoutCertificateIdentity, LayoutCertificateLimits, LayoutCertificateProblem,
    LayoutCertificateStatus, LayoutLimits, LayoutLp, PdhgError, PdhgSettings,
    TrussConstructionError,
};
use std::collections::BTreeSet;

/// Maximum grid nodes admitted to the cubic ground-structure constructor.
pub const MAX_TRUSS_CAMPAIGN_NODES: usize = 256;
/// Maximum cubic node-triplet checks admitted before ground construction.
pub const MAX_TRUSS_GROUND_CHECKS: usize = 262_144;
/// Maximum candidate members retained for one campaign solve.
pub const MAX_TRUSS_CANDIDATE_MEMBERS: usize = 512;
/// Maximum conservative scalar operations admitted to the fixed PDHG solve.
pub const MAX_TRUSS_PDHG_SCALAR_STEPS: usize = 1 << 27;
/// Version of the exact-input load-path certificate receipt.
pub const LOAD_PATH_CERTIFICATE_VERSION: u32 = 1;
/// Relative member-force threshold used by the certified active-set policy.
pub const LOAD_PATH_ACTIVE_RELATIVE_THRESHOLD: f64 = 1e-3;
/// Positive scale floor used before applying the relative threshold.
pub const LOAD_PATH_ACTIVE_FORCE_FLOOR: f64 = 1e-12;

const TRUSS_PDHG_MAX_ITERS: usize = 60_000;
const TRUSS_PDHG_CHECK_EVERY: usize = 500;

/// Structured TrussPath refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrussError {
    /// One public campaign field is outside its bounded numerical domain.
    InvalidInput {
        /// Stable field name.
        field: &'static str,
        /// Stable domain requirement.
        requirement: &'static str,
    },
    /// The bounded grid produced no admissible candidate member.
    NoCandidateMembers,
    /// The thresholded active graph has no connected multi-bar load-to-support path.
    NoCompleteLoadPath,
    /// A deterministic construction or solver work budget was exceeded.
    WorkBudget {
        /// Bounded resource.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Observed request, saturated on arithmetic overflow.
        observed: usize,
    },
    /// Solver-derived path data violated its checked domain.
    InvalidLoadPath {
        /// Stable diagnosis.
        reason: &'static str,
    },
    /// The checked PDHG solver refused its controls or warm-start state.
    Solver(PdhgError),
    /// The optimum-certificate attempt encountered malformed state, allocation
    /// failure, or cancellation. A sound numerical refusal is not an error and
    /// remains `Estimated`.
    Certificate(LayoutCertificateError),
    /// Ground-structure or LP construction refused before publishing output.
    Construction(TrussConstructionError),
    /// Tropical analysis refused solver-derived task data.
    Tropical(TropicalError),
}

impl core::fmt::Display for TrussError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { field, requirement } => {
                write!(formatter, "truss campaign {field} {requirement}")
            }
            Self::NoCandidateMembers => {
                formatter.write_str("truss campaign has no candidate members")
            }
            Self::NoCompleteLoadPath => formatter.write_str(
                "truss campaign has no connected multi-bar active path from load to support",
            ),
            Self::WorkBudget {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "truss campaign {resource} work {observed} exceeds limit {limit}"
            ),
            Self::InvalidLoadPath { reason } => {
                write!(formatter, "truss load-path input {reason}")
            }
            Self::Solver(error) => write!(formatter, "truss solver refused: {error}"),
            Self::Certificate(error) => {
                write!(formatter, "truss optimum certificate refused: {error}")
            }
            Self::Construction(error) => {
                write!(formatter, "truss construction refused: {error}")
            }
            Self::Tropical(error) => write!(formatter, "truss load-path analysis refused: {error}"),
        }
    }
}

impl std::error::Error for TrussError {}

impl From<TropicalError> for TrussError {
    fn from(value: TropicalError) -> Self {
        Self::Tropical(value)
    }
}

impl From<PdhgError> for TrussError {
    fn from(value: PdhgError) -> Self {
        Self::Solver(value)
    }
}

impl From<LayoutCertificateError> for TrussError {
    fn from(value: LayoutCertificateError) -> Self {
        Self::Certificate(value)
    }
}

impl From<TrussConstructionError> for TrussError {
    fn from(value: TrussConstructionError) -> Self {
        match value {
            TrussConstructionError::WorkBudget {
                resource,
                limit,
                observed,
            } => Self::WorkBudget {
                resource,
                limit,
                observed,
            },
            other => Self::Construction(other),
        }
    }
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct TrussReport {
    /// Candidate bars in the ground structure.
    pub num_members: usize,
    /// Bars in the separated certificate active set, or the rounded advisory
    /// active set when certification is unavailable.
    pub num_active: usize,
    /// Approximate primal volume of the returned PDHG iterate.
    pub total_volume: f64,
    /// Relative primal/dual objective separation diagnostic from PDHG.
    pub gap: f64,
    /// The equilibrium residual `‖Ax−b‖/‖b‖`.
    pub eq_residual: f64,
    /// PDHG iterations run.
    pub iters: usize,
    /// Did the iterative solver meet its gap and equilibrium-residual targets?
    pub solver_converged: bool,
    /// The advisory load path as original bar indices (load → support).
    pub critical_path: Vec<usize>,
    /// The volume carried by the critical path (tropical makespan).
    pub critical_path_volume: f64,
    /// The uniquely heaviest bar on a unique advisory path (original index).
    pub bottleneck_member: Option<usize>,
    /// Certified optimum bounds, or an honest diagnostic estimate when the
    /// outward proof is unavailable.
    pub optimality_color: Color,
    /// Load-path evidence. `Verified` is reachable only through
    /// [`LoadPathCertificateStatus::Certified`].
    pub load_path_color: Color,
    /// Retained proof or an exact reason the advisory path was not promoted.
    pub load_path_status: LoadPathCertificateStatus,
}

/// Checked advisory load-path analysis shared by native and browser campaigns.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadPathAnalysis {
    /// Original member indices, ordered from the load node to a support.
    pub members: Vec<usize>,
    /// Rounded sum of the selected member weights.
    pub weight: f64,
    /// A uniquely heaviest member when both path and weight ranking are unique.
    pub bottleneck_member: Option<usize>,
    /// Whether directed rounding separates the selected path from all rivals.
    pub path_is_unique: bool,
}

/// Sound reason an otherwise well-formed advisory path could not be promoted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LoadPathCertificateRefusal {
    /// The underlying optimum certificate was numerically unavailable.
    OptimumCertificateUnavailable,
    /// The optimum receipt did not bind the supplied LP, iterates, and settings.
    SolverIdentityMismatch,
    /// One member-force box straddled the interval active-set threshold.
    ActiveSetUnseparated {
        /// Physical member index.
        member: usize,
    },
    /// An active member's endpoint distances could not be strictly ordered.
    OrientationUnseparated {
        /// Physical member index.
        member: usize,
    },
    /// No complete multi-bar load-to-support path survived certified filtering.
    NoCompleteLoadPath,
    /// The selected path's lower weight did not exceed every rival upper weight.
    CriticalPathUnseparated,
    /// The selected bottleneck's lower weight did not exceed every peer upper weight.
    BottleneckUnseparated,
    /// Outward arithmetic became non-finite.
    NonFiniteArithmetic {
        /// Stable proof stage.
        stage: &'static str,
    },
}

impl LoadPathCertificateRefusal {
    const fn estimator(&self) -> &'static str {
        match self {
            Self::OptimumCertificateUnavailable => {
                "interval-load-path-optimum-certificate-unavailable-v1"
            }
            Self::SolverIdentityMismatch => "interval-load-path-solver-identity-mismatch-v1",
            Self::ActiveSetUnseparated { .. } => "interval-load-path-active-set-unseparated-v1",
            Self::OrientationUnseparated { .. } => "interval-load-path-orientation-unseparated-v1",
            Self::NoCompleteLoadPath => "interval-load-path-incomplete-v1",
            Self::CriticalPathUnseparated => "interval-load-path-critical-tie-v1",
            Self::BottleneckUnseparated => "interval-load-path-bottleneck-tie-v1",
            Self::NonFiniteArithmetic { .. } => "interval-load-path-non-finite-v1",
        }
    }
}

/// Exact native/browser input identity retained by a load-path proof.
///
/// The solver identities are collision-resistant BLAKE3 receipts from
/// `fs-truss`. Geometry, endpoint, and threshold identity is retained in full
/// canonical form, avoiding a second hash trust boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
struct LoadPathInputIdentity {
    version: u32,
    problem: LayoutCertificateIdentity,
    solver_input: LayoutCertificateIdentity,
    solver_certificate: LayoutCertificateIdentity,
    nodes: Vec<[u64; 2]>,
    members: Vec<(usize, usize)>,
    load_node: usize,
    support_nodes: Vec<usize>,
    relative_threshold: u64,
    force_floor: u64,
}

/// Private-by-construction interval certificate for one load-to-support path.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadPathCertificate {
    identity: LoadPathInputIdentity,
    analysis: LoadPathAnalysis,
    path_weight: Interval,
    active_threshold: Interval,
    active_members: Vec<usize>,
    member_weights: Vec<Interval>,
    replay_golden: u64,
}

impl LoadPathCertificate {
    /// Certified advisory path and its deterministic representative weight.
    #[must_use]
    pub const fn analysis(&self) -> &LoadPathAnalysis {
        &self.analysis
    }

    /// Outward interval enclosing the complete selected path weight.
    #[must_use]
    pub const fn path_weight_bounds(&self) -> Interval {
        self.path_weight
    }

    /// Outward active-force threshold used for every membership comparison.
    #[must_use]
    pub const fn active_threshold(&self) -> Interval {
        self.active_threshold
    }

    /// Exactly separated active member identities.
    #[must_use]
    pub fn active_members(&self) -> &[usize] {
        &self.active_members
    }

    /// Outward material-volume product for every physical member.
    #[must_use]
    pub fn member_weight_bounds(&self) -> &[Interval] {
        &self.member_weights
    }

    /// Collision-resistant identity of the complete lower-level solver proof.
    #[must_use]
    pub const fn solver_certificate_identity(&self) -> LayoutCertificateIdentity {
        self.identity.solver_certificate
    }

    /// Deterministic 64-bit replay sentinel over the exact receipt.
    ///
    /// This legacy FNV golden detects native/browser drift; it is not authority
    /// and cannot replace the retained exact fields or the BLAKE3 solver receipt.
    #[must_use]
    pub const fn replay_golden(&self) -> u64 {
        self.replay_golden
    }

    /// Re-run the bounded proof and require exact receipt equality.
    ///
    /// # Errors
    /// Returns a structured malformed-input, cancellation, or allocation
    /// error. A clean `false` means the supplied geometry, endpoints, solver
    /// state, threshold implementation, or proof output differs.
    #[allow(clippy::too_many_arguments)] // Exact proof identity has no implicit inputs.
    pub fn verifies_for(
        &self,
        problem: &LayoutCertificateProblem<'_>,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        optimum_status: &LayoutCertificateStatus,
        nodes: &[[f64; 2]],
        members: &[(usize, usize)],
        load_node: usize,
        support_nodes: &[usize],
        cx: &Cx<'_>,
    ) -> Result<bool, TrussError> {
        let replayed = certify_load_path(
            problem,
            x,
            y,
            settings,
            optimum_status,
            nodes,
            members,
            load_node,
            support_nodes,
            cx,
        )?;
        Ok(matches!(
            &replayed,
            LoadPathCertificateStatus::Certified(candidate) if candidate == self
        ))
    }
}

/// Result of a well-formed interval load-path certificate attempt.
#[derive(Debug, Clone, PartialEq)]
// Boxing the successful receipt would add an unreported allocation to a path
// whose other retained-state allocations are explicit.
#[allow(clippy::large_enum_variant)]
pub enum LoadPathCertificateStatus {
    /// Every promotion obligation separated and the exact receipt was retained.
    Certified(LoadPathCertificate),
    /// The rounded advisory path remains useful but cannot be called verified.
    Unavailable(LoadPathCertificateRefusal),
}

/// Convert a load-path certificate status into the sole evidence-promotion gate.
#[must_use]
pub fn load_path_color_from_certificate(status: &LoadPathCertificateStatus) -> Color {
    match status {
        LoadPathCertificateStatus::Certified(certificate) => {
            let bounds = certificate.path_weight_bounds();
            // declared-color-ok: candidate from the local path-weight interval certificate; admitted only at a consumer's authority boundary (6pf9)
            Color::Verified {
                lo: bounds.lo(),
                hi: bounds.hi(),
            }
        }
        LoadPathCertificateStatus::Unavailable(reason) => Color::Estimated {
            estimator: reason.estimator().to_string(),
            dispersion: f64::INFINITY,
        },
    }
}

/// Convert one outward certificate result into the sole optimality-promotion
/// gate shared by native and browser TrussPath consumers.
///
/// A finite PDHG gap or equilibrium residual is never sufficient for
/// `Verified`. An unavailable proof, or a proof whose private receipt no longer
/// validates, falls back to an explicitly diagnostic `Estimated` color.
///
/// # Errors
/// Returns a structured cancellation or allocation error if the bounded,
/// context-binding receipt preflight cannot finish.
#[allow(clippy::too_many_arguments)] // Complete problem, solver state, diagnostics, and Cx binding.
pub fn optimality_color_from_certificate(
    problem: &LayoutCertificateProblem<'_>,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    status: &LayoutCertificateStatus,
    gap: f64,
    eq_residual: f64,
    cx: &Cx<'_>,
) -> Result<Color, LayoutCertificateError> {
    if let LayoutCertificateStatus::Certified(certificate) = status
        && certificate.verifies_for_problem(problem, x, y, settings, cx)?
    {
        let bounds = certificate.bounds();
        // declared-color-ok: candidate from the local compliance interval certificate; admitted only at a consumer's authority boundary (6pf9)
        return Ok(Color::Verified {
            lo: bounds.lower(),
            hi: bounds.upper(),
        });
    }

    Ok(estimated_optimality_color(gap, eq_residual))
}

/// Preserve finite PDHG diagnostics without implying a proved optimum bound.
///
/// Browser consumers use this same fallback when malformed state, allocation
/// failure, or cancellation prevents a complete certificate attempt.
#[must_use]
pub fn estimated_optimality_color(gap: f64, eq_residual: f64) -> Color {
    Color::Estimated {
        estimator: "pdhg-diagnostics-with-unavailable-optimum-certificate-v1".to_string(),
        dispersion: if gap.is_finite()
            && eq_residual.is_finite()
            && gap >= 0.0
            && eq_residual >= 0.0
        {
            gap.max(eq_residual)
        } else {
            f64::INFINITY
        },
    }
}

/// Divide a verified optimum interval by a finite positive physical scale
/// using outward arithmetic.
///
/// This preserves an existing proof; it never promotes weaker evidence. An
/// invalid scale or malformed interval is demoted to an infinite-dispersion
/// diagnostic estimate.
#[must_use]
pub fn rescale_optimality_color(color: &Color, positive_divisor: f64) -> Color {
    match color {
        // declared-color-ok: guarded match-arm pattern READ split across lines; destructures rank, constructs nothing (6pf9)
        Color::Verified { lo, hi }
            if positive_divisor.is_finite()
                && positive_divisor > 0.0
                && lo.is_finite()
                && hi.is_finite()
                && lo <= hi =>
        {
            let scaled = Interval::new(*lo, *hi) / Interval::point(positive_divisor);
            if scaled.lo().is_finite() && scaled.hi().is_finite() {
                // declared-color-ok: rescaled candidate keeps its source certificate's declared rank; admitted only at a consumer's authority boundary (6pf9)
                Color::Verified {
                    lo: scaled.lo(),
                    hi: scaled.hi(),
                }
            } else {
                estimated_optimality_color(f64::INFINITY, f64::INFINITY)
            }
        }
        Color::Verified { .. } => estimated_optimality_color(f64::INFINITY, f64::INFINITY),
        weaker => weaker.clone(),
    }
}

#[derive(Debug, Clone, Copy)]
struct OrientedMember {
    original: usize,
    from: usize,
    to: usize,
    weight: f64,
}

/// Extract a connected, strictly support-ward path from thresholded members.
///
/// Every retained task is reachable from `load_node` and can reach one of the
/// indexed supports. This prevents a heavy disconnected component or an
/// interior suffix from being mislabeled as a load-to-support chain.
///
/// # Errors
/// Refuses malformed indices, non-finite geometry/weights, duplicate active or
/// support identities, bounded-resource excess, and graphs without a connected
/// path of at least two bars.
#[allow(clippy::too_many_lines)] // One bounded load-to-support graph witness and verifier.
pub fn analyze_load_path(
    nodes: &[[f64; 2]],
    members: &[(usize, usize)],
    active: &[usize],
    weights: &[f64],
    load_node: usize,
    support_nodes: &[usize],
) -> Result<LoadPathAnalysis, TrussError> {
    if nodes.is_empty() || nodes.len() > MAX_TRUSS_CAMPAIGN_NODES {
        return Err(TrussError::InvalidLoadPath {
            reason: "node count must be within the campaign bound",
        });
    }
    if members.len() != weights.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "member and weight counts must match",
        });
    }
    if load_node >= nodes.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "load node is out of range",
        });
    }
    if support_nodes.is_empty() || support_nodes.len() > nodes.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "support count must be within 1..=node count",
        });
    }
    if active.len() > MAX_TASK_DAG_NODES {
        return Err(TrussError::WorkBudget {
            resource: "active member count",
            limit: MAX_TASK_DAG_NODES,
            observed: active.len(),
        });
    }
    if nodes
        .iter()
        .flatten()
        .any(|coordinate| !coordinate.is_finite())
    {
        return Err(TrussError::InvalidLoadPath {
            reason: "node coordinates must be finite",
        });
    }

    let supports: BTreeSet<usize> = support_nodes.iter().copied().collect();
    if supports.len() != support_nodes.len()
        || supports.iter().any(|&node| node >= nodes.len())
        || supports.contains(&load_node)
    {
        return Err(TrussError::InvalidLoadPath {
            reason: "supports must be unique, in range, and exclude the load node",
        });
    }
    let active_set: BTreeSet<usize> = active.iter().copied().collect();
    if active_set.len() != active.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "active member identities must be unique",
        });
    }

    let support_points: Vec<[f64; 2]> = supports.iter().map(|&index| nodes[index]).collect();
    let distance: Vec<f64> = nodes
        .iter()
        .map(|point| {
            support_points
                .iter()
                .map(|support| (point[0] - support[0]).hypot(point[1] - support[1]))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    if distance.iter().any(|value| !value.is_finite()) {
        return Err(TrussError::InvalidLoadPath {
            reason: "distance-to-support must be finite",
        });
    }

    let mut oriented = Vec::with_capacity(active.len());
    for &member in active {
        let Some(&(a, b)) = members.get(member) else {
            return Err(TrussError::InvalidLoadPath {
                reason: "active member is out of range",
            });
        };
        if a >= nodes.len() || b >= nodes.len() || a == b {
            return Err(TrussError::InvalidLoadPath {
                reason: "active member endpoints must be distinct and in range",
            });
        }
        let weight = weights[member];
        if !weight.is_finite() || weight <= 0.0 {
            return Err(TrussError::InvalidLoadPath {
                reason: "active member weights must be finite and positive",
            });
        }
        let (from, to) = if distance[a] > distance[b] {
            (a, b)
        } else if distance[b] > distance[a] {
            (b, a)
        } else {
            // Equal-distance members do not make strictly support-ward progress.
            continue;
        };
        oriented.push(OrientedMember {
            original: member,
            from,
            to,
            weight,
        });
    }

    let mut reachable = vec![false; nodes.len()];
    reachable[load_node] = true;
    loop {
        let mut changed = false;
        for member in &oriented {
            if reachable[member.from] && !reachable[member.to] {
                reachable[member.to] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let mut reaches_support = vec![false; nodes.len()];
    for &support in &supports {
        reaches_support[support] = true;
    }
    loop {
        let mut changed = false;
        for member in oriented.iter().rev() {
            if reaches_support[member.to] && !reaches_support[member.from] {
                reaches_support[member.from] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    oriented.retain(|member| reachable[member.from] && reaches_support[member.to]);
    oriented.sort_by(|a, b| {
        distance[b.from]
            .total_cmp(&distance[a.from])
            .then(a.original.cmp(&b.original))
    });
    if oriented.len() < 2 {
        return Err(TrussError::NoCompleteLoadPath);
    }

    let mut starts_at = vec![Vec::new(); nodes.len()];
    for (index, member) in oriented.iter().enumerate() {
        starts_at[member.from].push(index);
    }
    let mut dag = TaskDag::new(oriented.iter().map(|member| member.weight).collect());
    let mut edge_count = 0usize;
    for (index, member) in oriented.iter().enumerate() {
        for &successor in &starts_at[member.to] {
            edge_count = edge_count.checked_add(1).ok_or(TrussError::WorkBudget {
                resource: "load-path edge count",
                limit: MAX_TASK_DAG_EDGES,
                observed: usize::MAX,
            })?;
            if edge_count > MAX_TASK_DAG_EDGES {
                return Err(TrussError::WorkBudget {
                    resource: "load-path edge count",
                    limit: MAX_TASK_DAG_EDGES,
                    observed: edge_count,
                });
            }
            dag = dag.with_edge(index, successor);
        }
    }

    let critical = dag.critical_path()?;
    let path: Vec<OrientedMember> = critical.path.iter().map(|&index| oriented[index]).collect();
    if path.len() < 2
        || path.first().is_none_or(|member| member.from != load_node)
        || path
            .last()
            .is_none_or(|member| !supports.contains(&member.to))
        || path.windows(2).any(|pair| pair[0].to != pair[1].from)
    {
        return Err(TrussError::NoCompleteLoadPath);
    }
    let bottleneck_member = dag.bottleneck()?.map(|index| oriented[index].original);
    Ok(LoadPathAnalysis {
        members: path.iter().map(|member| member.original).collect(),
        weight: critical.makespan,
        bottleneck_member,
        path_is_unique: critical.path_is_unique,
    })
}

#[derive(Debug, Clone, Copy)]
struct CertifiedOrientedMember {
    original: usize,
    from: usize,
    to: usize,
    weight: Interval,
}

fn path_poll(cx: &Cx<'_>, stage: &'static str) -> Result<(), TrussError> {
    cx.checkpoint()
        .map_err(|_| TrussError::Certificate(LayoutCertificateError::Cancelled { stage }))
}

fn finite_interval(interval: Interval) -> bool {
    interval.lo().is_finite() && interval.hi().is_finite() && interval.lo() <= interval.hi()
}

fn distance_interval(a: [f64; 2], b: [f64; 2]) -> Interval {
    if a[0].to_bits() == b[0].to_bits() && a[1].to_bits() == b[1].to_bits() {
        return Interval::point(0.0);
    }
    let dx = if a[0].to_bits() == b[0].to_bits() {
        Interval::point(0.0)
    } else {
        Interval::point(a[0]) - Interval::point(b[0])
    };
    let dy = if a[1].to_bits() == b[1].to_bits() {
        Interval::point(0.0)
    } else {
        Interval::point(a[1]) - Interval::point(b[1])
    };
    let abs_x = dx.abs();
    let abs_y = dy.abs();
    (abs_x * abs_x + abs_y * abs_y).sqrt()
}

fn nearest_support_distances(
    nodes: &[[f64; 2]],
    supports: &BTreeSet<usize>,
    cx: &Cx<'_>,
) -> Result<Vec<Interval>, TrussError> {
    let mut distances = Vec::new();
    distances.try_reserve_exact(nodes.len()).map_err(|_| {
        TrussError::Certificate(LayoutCertificateError::AllocationFailed {
            resource: "load-path support distances",
            requested: nodes.len(),
        })
    })?;
    for (node_index, &node) in nodes.iter().enumerate() {
        if node_index.is_multiple_of(64) {
            path_poll(cx, "load-path support-distance enclosure")?;
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::INFINITY;
        for &support in supports {
            let distance = distance_interval(node, nodes[support]);
            lo = lo.min(distance.lo());
            hi = hi.min(distance.hi());
        }
        if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || lo > hi {
            return Err(TrussError::InvalidLoadPath {
                reason: "support-distance enclosure must be finite and nonnegative",
            });
        }
        distances.push(Interval::new(lo, hi));
    }
    Ok(distances)
}

fn checked_upper_sum(left: f64, right: f64) -> Option<f64> {
    let sum = (Interval::point(left) + Interval::point(right)).hi();
    sum.is_finite().then_some(sum)
}

fn update_max(slot: &mut Option<f64>, candidate: f64) {
    *slot = Some(slot.map_or(candidate, |current| current.max(candidate)));
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_usize(bytes: &mut Vec<u8>, value: usize) {
    push_u64(bytes, u64::try_from(value).unwrap_or(u64::MAX));
}

fn replay_golden(
    identity: &LoadPathInputIdentity,
    analysis: &LoadPathAnalysis,
    path_weight: Interval,
    active_threshold: Interval,
    active_members: &[usize],
    member_weights: &[Interval],
) -> u64 {
    let mut canonical = Vec::new();
    canonical.extend_from_slice(b"fs-truss-e2e-load-path-replay-golden-v1");
    canonical.extend_from_slice(&identity.version.to_le_bytes());
    canonical.extend_from_slice(identity.problem.as_bytes());
    canonical.extend_from_slice(identity.solver_input.as_bytes());
    canonical.extend_from_slice(identity.solver_certificate.as_bytes());
    push_usize(&mut canonical, identity.nodes.len());
    for node in &identity.nodes {
        push_u64(&mut canonical, node[0]);
        push_u64(&mut canonical, node[1]);
    }
    push_usize(&mut canonical, identity.members.len());
    for &(a, b) in &identity.members {
        push_usize(&mut canonical, a);
        push_usize(&mut canonical, b);
    }
    push_usize(&mut canonical, identity.load_node);
    push_usize(&mut canonical, identity.support_nodes.len());
    for &support in &identity.support_nodes {
        push_usize(&mut canonical, support);
    }
    push_u64(&mut canonical, identity.relative_threshold);
    push_u64(&mut canonical, identity.force_floor);
    push_u64(&mut canonical, active_threshold.lo().to_bits());
    push_u64(&mut canonical, active_threshold.hi().to_bits());
    push_usize(&mut canonical, active_members.len());
    for &member in active_members {
        push_usize(&mut canonical, member);
    }
    push_usize(&mut canonical, member_weights.len());
    for &weight in member_weights {
        push_u64(&mut canonical, weight.lo().to_bits());
        push_u64(&mut canonical, weight.hi().to_bits());
    }
    push_usize(&mut canonical, analysis.members.len());
    for &member in &analysis.members {
        push_usize(&mut canonical, member);
    }
    push_u64(&mut canonical, analysis.weight.to_bits());
    push_usize(
        &mut canonical,
        analysis.bottleneck_member.unwrap_or(usize::MAX),
    );
    canonical.push(u8::from(analysis.path_is_unique));
    push_u64(&mut canonical, path_weight.lo().to_bits());
    push_u64(&mut canonical, path_weight.hi().to_bits());
    ProvenanceHash::of_bytes(&canonical).0
}

fn validate_certified_path_inputs(
    problem: &LayoutCertificateProblem<'_>,
    nodes: &[[f64; 2]],
    members: &[(usize, usize)],
    load_node: usize,
    support_nodes: &[usize],
) -> Result<BTreeSet<usize>, TrussError> {
    if nodes.is_empty() || nodes.len() > MAX_TRUSS_CAMPAIGN_NODES {
        return Err(TrussError::InvalidLoadPath {
            reason: "node count must be within the campaign bound",
        });
    }
    if members.is_empty()
        || members.len() > MAX_TASK_DAG_NODES
        || problem.c().len() != 2 * members.len()
    {
        return Err(TrussError::InvalidLoadPath {
            reason: "member count must match the split-variable certificate problem",
        });
    }
    if load_node >= nodes.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "load node is out of range",
        });
    }
    if support_nodes.is_empty() || support_nodes.len() > nodes.len() {
        return Err(TrussError::InvalidLoadPath {
            reason: "support count must be within 1..=node count",
        });
    }
    if nodes
        .iter()
        .flatten()
        .any(|coordinate| !coordinate.is_finite())
    {
        return Err(TrussError::InvalidLoadPath {
            reason: "node coordinates must be finite",
        });
    }
    for &(a, b) in members {
        if a >= nodes.len() || b >= nodes.len() || a == b {
            return Err(TrussError::InvalidLoadPath {
                reason: "member endpoints must be distinct and in range",
            });
        }
    }
    let supports: BTreeSet<usize> = support_nodes.iter().copied().collect();
    if supports.len() != support_nodes.len()
        || supports.iter().any(|&node| node >= nodes.len())
        || supports.contains(&load_node)
    {
        return Err(TrussError::InvalidLoadPath {
            reason: "supports must be unique, in range, and exclude the load node",
        });
    }
    Ok(supports)
}

/// Attempt the positive interval certificate for a load-to-support path.
///
/// The lower-layer optimum receipt supplies outward signed and split member-
/// force boxes. This gate intervalizes the relative active threshold and every
/// cost-times-split-force weight, requires every included and excluded member
/// to separate, orients every active member only across disjoint support-
/// distance boxes, removes disconnected/interior tasks, and proves one path
/// and one bottleneck by strict lower-versus-rival-upper comparisons.
///
/// # Errors
/// Returns a structured malformed-input, allocation, cancellation, or tropical
/// error. Sound inability to separate is [`LoadPathCertificateStatus::Unavailable`].
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // Full exact receipt and proof story.
pub fn certify_load_path(
    problem: &LayoutCertificateProblem<'_>,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    optimum_status: &LayoutCertificateStatus,
    nodes: &[[f64; 2]],
    members: &[(usize, usize)],
    load_node: usize,
    support_nodes: &[usize],
    cx: &Cx<'_>,
) -> Result<LoadPathCertificateStatus, TrussError> {
    let supports =
        validate_certified_path_inputs(problem, nodes, members, load_node, support_nodes)?;
    path_poll(cx, "load-path certificate admission")?;
    let LayoutCertificateStatus::Certified(optimum) = optimum_status else {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::OptimumCertificateUnavailable,
        ));
    };
    if !optimum.verifies_for_problem(problem, x, y, settings, cx)? {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::SolverIdentityMismatch,
        ));
    }

    let member_count = members.len();
    if optimum.repaired_member_forces().len() != member_count
        || optimum.repaired_split_forces().len() != 2 * member_count
    {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::SolverIdentityMismatch,
        ));
    }

    let mut absolute_forces = Vec::new();
    let mut member_weights = Vec::new();
    absolute_forces
        .try_reserve_exact(member_count)
        .map_err(|_| {
            TrussError::Certificate(LayoutCertificateError::AllocationFailed {
                resource: "load-path absolute-force intervals",
                requested: member_count,
            })
        })?;
    member_weights
        .try_reserve_exact(member_count)
        .map_err(|_| {
            TrussError::Certificate(LayoutCertificateError::AllocationFailed {
                resource: "load-path member-volume intervals",
                requested: member_count,
            })
        })?;
    let mut max_force_lo = 0.0f64;
    let mut max_force_hi = 0.0f64;
    for member in 0..member_count {
        if member.is_multiple_of(64) {
            path_poll(cx, "load-path member interval construction")?;
        }
        let cost = problem.c()[member];
        if !cost.is_finite()
            || cost <= 0.0
            || cost.to_bits() != problem.c()[member_count + member].to_bits()
        {
            return Err(TrussError::InvalidLoadPath {
                reason: "member costs must be finite, positive, and symmetric",
            });
        }
        let force = optimum.repaired_member_forces()[member].abs();
        if !finite_interval(force) || force.lo() < 0.0 {
            return Ok(LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::NonFiniteArithmetic {
                    stage: "absolute member force",
                },
            ));
        }
        max_force_lo = max_force_lo.max(force.lo());
        max_force_hi = max_force_hi.max(force.hi());
        absolute_forces.push(force);

        let split_sum = optimum.repaired_split_forces()[member]
            + optimum.repaired_split_forces()[member_count + member];
        let raw_weight = Interval::point(cost) * split_sum;
        if !finite_interval(raw_weight) || raw_weight.hi() < 0.0 {
            return Ok(LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::NonFiniteArithmetic {
                    stage: "member-volume product",
                },
            ));
        }
        member_weights.push(Interval::new(raw_weight.lo().max(0.0), raw_weight.hi()));
    }

    let force_scale = Interval::new(
        max_force_lo.max(LOAD_PATH_ACTIVE_FORCE_FLOOR),
        max_force_hi.max(LOAD_PATH_ACTIVE_FORCE_FLOOR),
    );
    let raw_threshold = Interval::point(LOAD_PATH_ACTIVE_RELATIVE_THRESHOLD) * force_scale;
    if !finite_interval(raw_threshold) || raw_threshold.hi() <= 0.0 {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NonFiniteArithmetic {
                stage: "active-set threshold",
            },
        ));
    }
    let active_threshold = Interval::new(raw_threshold.lo().max(0.0), raw_threshold.hi());
    let mut active = Vec::new();
    active.try_reserve_exact(member_count).map_err(|_| {
        TrussError::Certificate(LayoutCertificateError::AllocationFailed {
            resource: "load-path active set",
            requested: member_count,
        })
    })?;
    for (member, force) in absolute_forces.iter().copied().enumerate() {
        if force.lo() > active_threshold.hi() {
            active.push(member);
        } else if force.hi() > active_threshold.lo() {
            return Ok(LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::ActiveSetUnseparated { member },
            ));
        }
    }
    if active.len() < 2 {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NoCompleteLoadPath,
        ));
    }

    let distances = nearest_support_distances(nodes, &supports, cx)?;
    let mut oriented = Vec::new();
    oriented.try_reserve_exact(active.len()).map_err(|_| {
        TrussError::Certificate(LayoutCertificateError::AllocationFailed {
            resource: "load-path oriented active members",
            requested: active.len(),
        })
    })?;
    for (active_position, &member) in active.iter().enumerate() {
        if active_position.is_multiple_of(64) {
            path_poll(cx, "load-path active-member orientation")?;
        }
        let (a, b) = members[member];
        let (from, to) = if distances[a].lo() > distances[b].hi() {
            (a, b)
        } else if distances[b].lo() > distances[a].hi() {
            (b, a)
        } else {
            return Ok(LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::OrientationUnseparated { member },
            ));
        };
        let weight = member_weights[member];
        if !finite_interval(weight) || weight.lo() <= 0.0 {
            return Ok(LoadPathCertificateStatus::Unavailable(
                LoadPathCertificateRefusal::NonFiniteArithmetic {
                    stage: "active member-volume interval",
                },
            ));
        }
        oriented.push(CertifiedOrientedMember {
            original: member,
            from,
            to,
            weight,
        });
    }

    let mut reachable = vec![false; nodes.len()];
    reachable[load_node] = true;
    loop {
        path_poll(cx, "load-path forward reachability")?;
        let mut changed = false;
        for member in &oriented {
            if reachable[member.from] && !reachable[member.to] {
                reachable[member.to] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let mut reaches_support = vec![false; nodes.len()];
    for &support in &supports {
        reaches_support[support] = true;
    }
    loop {
        path_poll(cx, "load-path reverse reachability")?;
        let mut changed = false;
        for member in oriented.iter().rev() {
            if reaches_support[member.to] && !reaches_support[member.from] {
                reaches_support[member.from] = true;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    oriented.retain(|member| reachable[member.from] && reaches_support[member.to]);
    oriented.sort_by(|a, b| {
        distances[b.from]
            .midpoint()
            .total_cmp(&distances[a.from].midpoint())
            .then(a.original.cmp(&b.original))
    });
    if oriented.len() < 2 {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NoCompleteLoadPath,
        ));
    }

    let mut starts_at = vec![Vec::new(); nodes.len()];
    for (index, member) in oriented.iter().enumerate() {
        starts_at[member.from].push(index);
    }
    let mut predecessors = vec![Vec::new(); oriented.len()];
    let mut dag = TaskDag::new(
        oriented
            .iter()
            .map(|member| member.weight.midpoint())
            .collect(),
    );
    let mut edge_count = 0usize;
    for (index, member) in oriented.iter().enumerate() {
        if index.is_multiple_of(64) {
            path_poll(cx, "load-path task-edge construction")?;
        }
        for &successor in &starts_at[member.to] {
            if successor <= index {
                return Ok(LoadPathCertificateStatus::Unavailable(
                    LoadPathCertificateRefusal::OrientationUnseparated {
                        member: member.original,
                    },
                ));
            }
            edge_count = edge_count.checked_add(1).ok_or(TrussError::WorkBudget {
                resource: "load-path edge count",
                limit: MAX_TASK_DAG_EDGES,
                observed: usize::MAX,
            })?;
            if edge_count > MAX_TASK_DAG_EDGES {
                return Err(TrussError::WorkBudget {
                    resource: "load-path edge count",
                    limit: MAX_TASK_DAG_EDGES,
                    observed: edge_count,
                });
            }
            predecessors[successor].push(index);
            dag = dag.with_edge(index, successor);
        }
    }

    let critical = dag.critical_path()?;
    if critical.path.len() < 2 {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NoCompleteLoadPath,
        ));
    }
    let path: Vec<CertifiedOrientedMember> =
        critical.path.iter().map(|&index| oriented[index]).collect();
    if path.first().is_none_or(|member| member.from != load_node)
        || path
            .last()
            .is_none_or(|member| !supports.contains(&member.to))
        || path.windows(2).any(|pair| pair[0].to != pair[1].from)
    {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NoCompleteLoadPath,
        ));
    }

    let mut path_weight = Interval::point(0.0);
    for (position, member) in path.iter().enumerate() {
        if position.is_multiple_of(64) {
            path_poll(cx, "load-path selected-weight sum")?;
        }
        path_weight = path_weight + member.weight;
    }
    if !finite_interval(path_weight) || path_weight.lo() <= 0.0 {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NonFiniteArithmetic {
                stage: "critical-path weight sum",
            },
        ));
    }

    let mut witness_position = vec![None; oriented.len()];
    for (position, &member) in critical.path.iter().enumerate() {
        witness_position[member] = Some(position);
    }
    let mut exact_prefix_upper = vec![None; oriented.len()];
    let mut alternative_upper = vec![None; oriented.len()];
    for index in 0..oriented.len() {
        if index.is_multiple_of(64) {
            path_poll(cx, "load-path interval rival comparison")?;
        }
        let member = oriented[index];
        if member.from == load_node {
            if witness_position[index] == Some(0) {
                exact_prefix_upper[index] = Some(member.weight.hi());
            } else {
                alternative_upper[index] = Some(member.weight.hi());
            }
        }
        for &predecessor in &predecessors[index] {
            if let Some(prefix) = alternative_upper[predecessor] {
                let Some(candidate) = checked_upper_sum(prefix, member.weight.hi()) else {
                    return Ok(LoadPathCertificateStatus::Unavailable(
                        LoadPathCertificateRefusal::NonFiniteArithmetic {
                            stage: "rival path upper sum",
                        },
                    ));
                };
                update_max(&mut alternative_upper[index], candidate);
            }
            if let (Some(prefix), Some(previous_position)) = (
                exact_prefix_upper[predecessor],
                witness_position[predecessor],
            ) {
                let Some(candidate) = checked_upper_sum(prefix, member.weight.hi()) else {
                    return Ok(LoadPathCertificateStatus::Unavailable(
                        LoadPathCertificateRefusal::NonFiniteArithmetic {
                            stage: "witness-prefix upper sum",
                        },
                    ));
                };
                if witness_position[index] == Some(previous_position + 1) {
                    exact_prefix_upper[index] = Some(candidate);
                } else {
                    update_max(&mut alternative_upper[index], candidate);
                }
            }
        }
    }
    let witness_last = critical.path.last().copied();
    let mut rival_upper = None;
    for (index, member) in oriented.iter().enumerate() {
        if supports.contains(&member.to) {
            if let Some(candidate) = alternative_upper[index] {
                update_max(&mut rival_upper, candidate);
            }
            if Some(index) != witness_last
                && let Some(candidate) = exact_prefix_upper[index]
            {
                update_max(&mut rival_upper, candidate);
            }
        }
    }
    if rival_upper.is_some_and(|rival| path_weight.lo() <= rival) {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::CriticalPathUnseparated,
        ));
    }

    let Some(bottleneck_task) = critical.path.iter().copied().max_by(|&left, &right| {
        oriented[left]
            .weight
            .midpoint()
            .total_cmp(&oriented[right].weight.midpoint())
            .then(oriented[right].original.cmp(&oriented[left].original))
    }) else {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NoCompleteLoadPath,
        ));
    };
    let bottleneck = oriented[bottleneck_task];
    if critical
        .path
        .iter()
        .copied()
        .any(|task| task != bottleneck_task && oriented[task].weight.hi() >= bottleneck.weight.lo())
    {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::BottleneckUnseparated,
        ));
    }

    let analysis = LoadPathAnalysis {
        members: path.iter().map(|member| member.original).collect(),
        weight: critical.makespan,
        bottleneck_member: Some(bottleneck.original),
        path_is_unique: true,
    };
    if !path_weight.contains(analysis.weight) {
        return Ok(LoadPathCertificateStatus::Unavailable(
            LoadPathCertificateRefusal::NonFiniteArithmetic {
                stage: "representative path weight containment",
            },
        ));
    }
    path_poll(cx, "load-path receipt identity")?;
    let identity = LoadPathInputIdentity {
        version: LOAD_PATH_CERTIFICATE_VERSION,
        problem: optimum.problem_identity(),
        solver_input: optimum.input_identity(),
        solver_certificate: optimum.certificate_identity(),
        nodes: nodes
            .iter()
            .map(|node| [node[0].to_bits(), node[1].to_bits()])
            .collect(),
        members: members.to_vec(),
        load_node,
        support_nodes: support_nodes.to_vec(),
        relative_threshold: LOAD_PATH_ACTIVE_RELATIVE_THRESHOLD.to_bits(),
        force_floor: LOAD_PATH_ACTIVE_FORCE_FLOOR.to_bits(),
    };
    let replay_golden = replay_golden(
        &identity,
        &analysis,
        path_weight,
        active_threshold,
        &active,
        &member_weights,
    );
    path_poll(cx, "load-path receipt publication")?;
    Ok(LoadPathCertificateStatus::Certified(LoadPathCertificate {
        identity,
        analysis,
        path_weight,
        active_threshold,
        active_members: active,
        member_weights,
        replay_golden,
    }))
}

/// Run the TrussPath campaign: a cantilever on an `nx×ny` grid over `[0,w]×[0,h]`,
/// left edge supported, a unit downward load at the free bottom corner.
///
/// # Errors
/// Returns a structured refusal for invalid/unbounded grid parameters, an empty
/// candidate set, an excessive construction/solver budget, or invalid
/// solver-derived path data. Ground and LP construction also return structured
/// allocation or cancellation refusals through the supplied context. The
/// outward certificate stage can additionally return [`TrussError::Certificate`]
/// for malformed retained state, allocation failure, or cancellation.
#[allow(clippy::too_many_lines)] // One bounded campaign, diagnostics, and evidence report pipeline.
pub fn run_campaign(
    nx: usize,
    ny: usize,
    w: f64,
    h: f64,
    gap_tol: f64,
    cx: &Cx<'_>,
) -> Result<TrussReport, TrussError> {
    if nx < 2 || ny < 2 {
        return Err(TrussError::InvalidInput {
            field: "grid dimensions",
            requirement: "must each be at least two",
        });
    }
    let node_count = nx.checked_mul(ny).ok_or(TrussError::InvalidInput {
        field: "grid node count",
        requirement: "must fit usize and the deterministic node budget",
    })?;
    if node_count > MAX_TRUSS_CAMPAIGN_NODES {
        return Err(TrussError::InvalidInput {
            field: "grid node count",
            requirement: "exceeds the 256-node deterministic work budget",
        });
    }
    let ground_checks = node_count
        .checked_mul(node_count)
        .and_then(|square| square.checked_mul(node_count))
        .ok_or(TrussError::WorkBudget {
            resource: "ground-structure triplet checks",
            limit: MAX_TRUSS_GROUND_CHECKS,
            observed: usize::MAX,
        })?;
    if ground_checks > MAX_TRUSS_GROUND_CHECKS {
        return Err(TrussError::WorkBudget {
            resource: "ground-structure triplet checks",
            limit: MAX_TRUSS_GROUND_CHECKS,
            observed: ground_checks,
        });
    }
    let max_extent = f64::MAX.sqrt() * 0.5;
    if !w.is_finite() || w <= 0.0 || w > max_extent {
        return Err(TrussError::InvalidInput {
            field: "width",
            requirement: "must be finite, positive, and safe for squared geometry",
        });
    }
    if !h.is_finite() || h <= 0.0 || h > max_extent {
        return Err(TrussError::InvalidInput {
            field: "height",
            requirement: "must be finite, positive, and safe for squared geometry",
        });
    }
    if !gap_tol.is_finite() || gap_tol <= 0.0 || gap_tol > 1.0 {
        return Err(TrussError::InvalidInput {
            field: "gap tolerance",
            requirement: "must be finite and in 0 < gap_tol <= 1",
        });
    }
    let max_member_length = w.hypot(h) / 1.5;
    if max_member_length < 0.1 {
        return Err(TrussError::NoCandidateMembers);
    }
    let rules = GroundRules::try_new(0.1, max_member_length, Vec::new(), 1e-6)?;
    let ground_limits = GroundLimits::try_new(
        MAX_TRUSS_CAMPAIGN_NODES,
        MAX_TRUSS_CAMPAIGN_NODES * (MAX_TRUSS_CAMPAIGN_NODES - 1) / 2,
        MAX_TRUSS_GROUND_CHECKS,
        MAX_TRUSS_CANDIDATE_MEMBERS,
        1 << 20,
    )?;
    let gs = GroundStructure::try_grid(nx, ny, w, h, &rules, ground_limits, cx)?;
    let m = gs.members().len();
    if m == 0 {
        return Err(TrussError::NoCandidateMembers);
    }
    if m > MAX_TRUSS_CANDIDATE_MEMBERS {
        return Err(TrussError::WorkBudget {
            resource: "candidate member count",
            limit: MAX_TRUSS_CANDIDATE_MEMBERS,
            observed: m,
        });
    }

    // Left edge supported; unit downward load at the free bottom-right node.
    let support_nodes: Vec<usize> = (0..ny).map(|row| row * nx).collect();
    let support_set: BTreeSet<usize> = support_nodes.iter().copied().collect();
    let load_node = nx - 1;
    let supported: Vec<[bool; 2]> = (0..node_count)
        .map(|node| [support_set.contains(&node); 2])
        .collect();
    let loads: Vec<[f64; 2]> = (0..node_count)
        .map(|node| {
            if node == load_node {
                [0.0, -1.0]
            } else {
                [0.0, 0.0]
            }
        })
        .collect();
    let case = LayoutCase::try_new(supported, loads, node_count)?;
    let lp = LayoutLp::try_assemble(&gs, &case, 1.0, LayoutLimits::default(), cx)?;
    // Two sparse multiply-add passes (4*nnz scalar arithmetic), the projected
    // primal update plus extrapolation (6*nvar), and the dual update (3*nrow).
    // Diagnostic checkpoints add two more SpMVs and bounded reductions.
    let per_iteration = lp
        .a()
        .nnz()
        .checked_mul(4)
        .and_then(|steps| {
            lp.c()
                .len()
                .checked_mul(10)
                .and_then(|vector_steps| steps.checked_add(vector_steps))
        })
        .and_then(|steps| {
            lp.b()
                .len()
                .checked_mul(3)
                .and_then(|row_steps| steps.checked_add(row_steps))
        })
        .ok_or(TrussError::WorkBudget {
            resource: "PDHG scalar steps",
            limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
            observed: usize::MAX,
        })?;
    let per_diagnostic = lp
        .a()
        .nnz()
        .checked_mul(4)
        .and_then(|steps| {
            lp.c()
                .len()
                .checked_mul(6)
                .and_then(|vector_steps| steps.checked_add(vector_steps))
        })
        .and_then(|steps| {
            lp.b()
                .len()
                .checked_mul(7)
                .and_then(|row_steps| steps.checked_add(row_steps))
        })
        .and_then(|steps| steps.checked_add(16))
        .ok_or(TrussError::WorkBudget {
            resource: "PDHG scalar steps",
            limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
            observed: usize::MAX,
        })?;
    let iteration_steps =
        per_iteration
            .checked_mul(TRUSS_PDHG_MAX_ITERS)
            .ok_or(TrussError::WorkBudget {
                resource: "PDHG scalar steps",
                limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
                observed: usize::MAX,
            })?;
    let diagnostic_steps = per_diagnostic
        .checked_mul(TRUSS_PDHG_MAX_ITERS.div_ceil(TRUSS_PDHG_CHECK_EVERY))
        .ok_or(TrussError::WorkBudget {
            resource: "PDHG scalar steps",
            limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
            observed: usize::MAX,
        })?;
    let solver_steps =
        iteration_steps
            .checked_add(diagnostic_steps)
            .ok_or(TrussError::WorkBudget {
                resource: "PDHG scalar steps",
                limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
                observed: usize::MAX,
            })?;
    if solver_steps > MAX_TRUSS_PDHG_SCALAR_STEPS {
        return Err(TrussError::WorkBudget {
            resource: "PDHG scalar steps",
            limit: MAX_TRUSS_PDHG_SCALAR_STEPS,
            observed: solver_steps,
        });
    }
    let settings = PdhgSettings {
        max_iters: TRUSS_PDHG_MAX_ITERS,
        gap_tol,
        check_every: TRUSS_PDHG_CHECK_EVERY,
    };
    let (x, y, mut report) = lp.solve(None, None, settings)?;
    let certificate_status = lp.certify_optimum_for_report(
        &x,
        &y,
        settings,
        &mut report,
        LayoutCertificateLimits::default(),
        cx,
    )?;
    let certificate_problem = LayoutCertificateProblem::try_new(lp.a(), lp.c(), lp.b())?;

    // Member force (q⁺ − q⁻) and material volume (both split costs).
    let force = |k: usize| x[k] - x[m + k];
    let volume = |k: usize| lp.c()[k] * x[k] + lp.c()[m + k] * x[m + k];
    let max_force = (0..m).map(|k| force(k).abs()).fold(0.0, f64::max);
    let active_tol =
        LOAD_PATH_ACTIVE_RELATIVE_THRESHOLD * max_force.max(LOAD_PATH_ACTIVE_FORCE_FLOOR);

    let active: Vec<usize> = (0..m).filter(|&k| force(k).abs() > active_tol).collect();
    let volumes: Vec<f64> = (0..m).map(volume).collect();
    let advisory_load_path = analyze_load_path(
        gs.nodes(),
        gs.members(),
        &active,
        &volumes,
        load_node,
        &support_nodes,
    )?;
    let load_path_status = certify_load_path(
        &certificate_problem,
        &x,
        &y,
        settings,
        &certificate_status,
        gs.nodes(),
        gs.members(),
        load_node,
        &support_nodes,
        cx,
    )?;
    let load_path = match &load_path_status {
        LoadPathCertificateStatus::Certified(certificate) => certificate.analysis().clone(),
        LoadPathCertificateStatus::Unavailable(_) => advisory_load_path,
    };
    let num_active = match &load_path_status {
        LoadPathCertificateStatus::Certified(certificate) => certificate.active_members().len(),
        LoadPathCertificateStatus::Unavailable(_) => active.len(),
    };
    let load_path_color = load_path_color_from_certificate(&load_path_status);

    let solver_converged = report.gap.is_finite()
        && report.eq_residual.is_finite()
        && report.gap >= 0.0
        && report.eq_residual >= 0.0
        && report.gap < gap_tol
        && report.eq_residual < gap_tol;
    let optimality_color = optimality_color_from_certificate(
        &certificate_problem,
        &x,
        &y,
        settings,
        &certificate_status,
        report.gap,
        report.eq_residual,
        cx,
    )?;

    Ok(TrussReport {
        num_members: m,
        num_active,
        total_volume: report.volume,
        gap: report.gap,
        eq_residual: report.eq_residual,
        iters: report.iters,
        solver_converged,
        critical_path: load_path.members,
        critical_path_volume: load_path.weight,
        bottleneck_member: load_path.bottleneck_member,
        optimality_color,
        load_path_color,
        load_path_status,
    })
}
