//! Static structural execution planner: block DAG, parallel waves, class &
//! time-scale grouping, tearing inspection, cost estimates, content-addressed
//! plan identity, human explanation, canonical JSON, and causalization audit
//! export (bead `frankensim-rc-1109-root-3gc12.5.7`).
//!
//! All execution workspaces are allocated once at [`crate::SystemBuilder::seal`];
//! [`System::execution_plan`] and [`System::causalization_audit`] are pure
//! read-only inspections that never mutate accepted or candidate state.

use core::fmt::{self, Write as _};

use fs_ports::component::Class;

use super::System;
use super::line::{Line, LineTear};
use crate::diagnostics::{BlockSolveMode, BlockSolverPolicy, FixedPointMethod};
use crate::wiring::Driver;

/// Domain tag for content-addressed [`ExecutionPlan`] hashes.
const EXECUTION_PLAN_DOMAIN: &str = "org.frankensim.fs-system.execution-plan.v1";

/// Schema version tag emitted in canonical JSON exports of [`ExecutionPlan`].
pub const EXECUTION_PLAN_SCHEMA: &str = "frankensim.execution_plan.v1";

/// Schema version tag emitted in canonical JSON exports of [`CausalizationAudit`].
pub const CAUSALIZATION_AUDIT_SCHEMA: &str = "frankensim.causalization_audit.v1";

/// Aggregate mathematical class of a strongly connected block, derived from
/// its member components' [`Class`] declarations and load-line presence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockClass {
    /// Every component in the block declares [`Class::ExactLinear`]. When all
    /// components also offer [`fs_ports::component::PortTangent`] and all
    /// drivers are linear, the block's condensed Jacobian factors depend on
    /// `dt_s` alone and are cached across steps.
    ExactLinear,
    /// Components are [`Class::ExactLinear`] or [`Class::Smooth`] (with at
    /// least one [`Class::Smooth`]), solved by Newton / Krylov / fixed-point.
    Smooth,
    /// Components are [`Class::ExactLinear`] or [`Class::Monotone1d`] (with at
    /// least one [`Class::Monotone1d`]), eligible for safeguarded monotone
    /// scalar solving.
    Monotone1d,
    /// Block contains a [`Class::NonSmooth`] component or a load-line
    /// complementarity/contact component.
    NonSmooth,
    /// Block mixes [`Class::Smooth`] and [`Class::Monotone1d`] components
    /// without any [`Class::NonSmooth`] or load-line component.
    Mixed,
}

impl BlockClass {
    /// Canonical lowercase identifier for reports, JSON, and plan hashes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ExactLinear => "exact_linear",
            Self::Smooth => "smooth",
            Self::Monotone1d => "monotone_1d",
            Self::NonSmooth => "non_smooth",
            Self::Mixed => "mixed",
        }
    }
}

impl fmt::Display for BlockClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// How the structural planner tears and solves a block's interface equations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TearStrategy {
    /// Acyclic block with zero algebraic loop unknowns; evaluated once in
    /// topological order with no linear solve.
    Acyclic,
    /// Copy/follower unknowns are eliminated in topological substitution order
    /// via [`crate::solve::Condensed`]; only the loop-closing tear unknowns are
    /// factored densely with partial pivoting.
    EliminationCondensed,
    /// Woodbury tear on a load-line component's active effort outputs
    /// (`J = M - B D_L P_L`), combining condensed neighbour elimination with a
    /// tiny dense solve of width `active_outputs`.
    LoadLineWoodbury,
    /// Monolithic Newton–Krylov solve over internal [`fs_ports::component::StageResidual`]
    /// unknowns and interface inputs via `fs-solver`.
    MonolithicStageResidual,
    /// Derivative-free fixed-point iteration with dynamic Aitken or vector
    /// IQN-ILS acceleration over the block's algebraic unknowns.
    FixedPointAccelerated,
}

impl TearStrategy {
    /// Canonical lowercase identifier for reports, JSON, and plan hashes.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Acyclic => "acyclic",
            Self::EliminationCondensed => "elimination_condensed",
            Self::LoadLineWoodbury => "load_line_woodbury",
            Self::MonolithicStageResidual => "monolithic_stage_residual",
            Self::FixedPointAccelerated => "fixed_point_accelerated",
        }
    }
}

impl fmt::Display for TearStrategy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Time-scale tier assigned to a block from its components' `limit_dt_s`
/// bounds relative to the slowest block at the probe horizon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum TimeScaleTier {
    /// Step-size limit at least 100x faster than the slowest block.
    Fast,
    /// Step-size limit between 10x and 100x faster than the slowest block.
    Medium,
    /// Step-size limit within 10x of the slowest block, yet bounded below the
    /// probe horizon.
    Slow,
    /// No component in the block restricted the probe horizon.
    Unconstrained,
}

impl TimeScaleTier {
    /// Canonical lowercase identifier for reports and JSON.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Fast => "fast",
            Self::Medium => "medium",
            Self::Slow => "slow",
            Self::Unconstrained => "unconstrained",
        }
    }
}

impl fmt::Display for TimeScaleTier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Group of blocks sharing one [`TimeScaleTier`] at the probe horizon.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeScaleGroup {
    /// Assigned time-scale tier.
    pub tier: TimeScaleTier,
    /// Block indices belonging to this tier, in execution order.
    pub blocks: Vec<usize>,
    /// Smallest finite `limit_dt_s` among blocks in this tier [s].
    pub min_dt_s: f64,
    /// Largest finite `limit_dt_s` among blocks in this tier [s].
    pub max_dt_s: f64,
    /// Ratio of the system's slowest block step limit to `min_dt_s` (`>= 1.0`).
    pub stiffness_ratio: f64,
}

/// Deterministic analytical cost and static-memory estimate for one block,
/// formatted for `frankensim explain` and `fs-plan` cost attribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockCostEstimate {
    /// Estimated component evaluations (`propose` + tangent calls) per warm step.
    pub estimated_evals_per_step: u64,
    /// Estimated dense LU factorization FLOPs per factor refresh (`0` when cached
    /// or matrix-free).
    pub estimated_dense_flops_per_factor: u64,
    /// Estimated substitution + triangular solve FLOPs per linear solve.
    pub estimated_substitution_flops_per_solve: u64,
    /// Pre-allocated static workspace bytes held for this block (zero per-step
    /// heap allocation after warm-up).
    pub preallocated_workspace_bytes: usize,
}

/// Static structural plan for one strongly connected block.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockPlan {
    /// Zero-based block index in deterministic Tarjan topological order.
    pub index: usize,
    /// Zero-based parallel execution wave (`0` has no upstream block
    /// dependencies; blocks in the same wave share no dependency edge).
    pub wave: usize,
    /// Upstream block indices whose outputs feed inputs of this block.
    pub dependencies: Vec<usize>,
    /// Downstream block indices whose inputs read outputs of this block.
    pub dependents: Vec<usize>,
    /// Member component names, in index order.
    pub nodes: Vec<String>,
    /// Declared [`Class`] of each member component, parallel to `nodes`.
    pub node_classes: Vec<Class>,
    /// Aggregate classification of the block.
    pub class: BlockClass,
    /// Active coupled-step solve mode selected under the current solver policy.
    pub solve_mode: BlockSolveMode,
    /// Tearing and linear-algebra strategy for the block.
    pub tear_strategy: TearStrategy,
    /// Whether the block contains an algebraic coupling loop.
    pub looped: bool,
    /// Whether the block's full condensed Jacobian factors are cached across
    /// identical step lengths (`Class::ExactLinear` + declared tangents + linear
    /// drivers).
    pub affine_cached: bool,
    /// Whether the block's load-line neighbour factors, compliance matrix, and
    /// Woodbury basis `Z = M^-1 B` are cached across identical step lengths.
    pub affine_line_cached: bool,
    /// Name of the load-line component in this block, if any.
    pub load_line_node: Option<String>,
    /// Total scalar interface unknowns inside the loop (`0` when acyclic).
    pub total_unknowns: usize,
    /// Number of torn unknowns solved simultaneously (for interface tearing
    /// strategies `<= total_unknowns`; for [`TearStrategy::MonolithicStageResidual`]
    /// equals `stage_unknowns + total_unknowns`).
    pub kept_unknowns: usize,
    /// Number of copy/follower unknowns eliminated by forward substitution.
    pub substituted_count: usize,
    /// Total internal [`fs_ports::component::StageResidual`] unknowns across
    /// components in this block.
    pub stage_unknowns: usize,
    /// Explicit labels (`node.port` or `node.port[k]`) of the torn unknowns
    /// solved simultaneously (for [`TearStrategy::LoadLineWoodbury`], neighbour
    /// kept unknowns plus the load-line component's active tear efforts).
    pub tear_variables: Vec<String>,
    /// Explicit labels of follower unknowns in exact substitution order.
    pub substituted_variables: Vec<String>,
    /// Time-scale tier at the probe horizon.
    pub time_scale_tier: TimeScaleTier,
    /// Tightest `limit_dt_s` returned by any component in this block at the
    /// probe horizon [s], or `None` if a component refused negotiation.
    pub limit_dt_s: Option<f64>,
    /// Analytical per-step cost and static workspace size.
    pub cost: BlockCostEstimate,
}

/// Content-addressed, explainable static execution plan for a sealed [`System`].
///
/// [`Self::plan_hash`] is a domain-separated BLAKE3 digest over the static
/// structural plan (problem identity, solver policy, block ordering, parallel
/// waves, classifications, tear strategies, tear/follower sets, caching flags,
/// and analytical costs). Transient runtime `limit_dt_s` probes are reported in
/// [`Self::time_scale_groups`] and [`BlockPlan::limit_dt_s`] without perturbing
/// [`Self::plan_hash`], so stepping a system preserves its structural plan hash.
#[derive(Debug, Clone, PartialEq)]
pub struct ExecutionPlan {
    /// Domain-separated BLAKE3 digest of the static structural execution plan.
    pub plan_hash: [u8; 32],
    /// Problem identity hash ([`System::problem_id`]).
    pub problem_id: [u8; 32],
    /// Number of components in the sealed system.
    pub node_count: usize,
    /// Number of connections in the sealed system.
    pub connection_count: usize,
    /// Number of prescribed external boundaries.
    pub boundary_count: usize,
    /// Number of strongly connected blocks.
    pub block_count: usize,
    /// Number of blocks containing an algebraic loop.
    pub loop_count: usize,
    /// Whether all hot-path step workspaces are pre-allocated statically at seal.
    pub static_zero_alloc_workspaces: bool,
    /// Probe horizon used when querying `limit_dt_s` for time-scale grouping [s].
    pub probe_dt_s: f64,
    /// Static block plans in deterministic execution order.
    pub blocks: Vec<BlockPlan>,
    /// Parallel execution waves: each wave lists independent block indices that
    /// have zero inter-block dependency edges between them.
    pub waves: Vec<Vec<usize>>,
    /// Blocks grouped by relative time scale at `probe_dt_s`.
    pub time_scale_groups: Vec<TimeScaleGroup>,
    /// Sum of `estimated_evals_per_step` across all blocks.
    pub total_estimated_evals_per_step: u64,
    /// Sum of factorization + substitution FLOPs across all blocks.
    pub total_estimated_flops_per_step: u64,
    /// Total static workspace bytes pre-allocated across all blocks.
    pub total_workspace_bytes: usize,
}

impl ExecutionPlan {
    /// Lowercase 64-character hexadecimal encoding of [`Self::plan_hash`].
    #[must_use]
    pub fn plan_hash_hex(&self) -> String {
        hex_encode(&self.plan_hash)
    }

    /// Lowercase 64-character hexadecimal encoding of [`Self::problem_id`].
    #[must_use]
    pub fn problem_id_hex(&self) -> String {
        hex_encode(&self.problem_id)
    }

    /// Total kept (torn) unknowns summed across all looped blocks.
    #[must_use]
    pub fn total_kept_unknowns(&self) -> usize {
        self.blocks.iter().map(|b| b.kept_unknowns).sum()
    }

    /// Human-readable engineering explanation of the plan (suitable for
    /// `frankensim explain` and diagnostic reports).
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn explain(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "ExecutionPlan ({EXECUTION_PLAN_SCHEMA})");
        let _ = writeln!(out, "  plan_hash:  {}", self.plan_hash_hex());
        let _ = writeln!(out, "  problem_id: {}", self.problem_id_hex());
        let _ = writeln!(
            out,
            "  graph:      {} nodes, {} connections, {} boundaries -> {} blocks ({} looped) in {} waves",
            self.node_count,
            self.connection_count,
            self.boundary_count,
            self.block_count,
            self.loop_count,
            self.waves.len()
        );
        let _ = writeln!(
            out,
            "  runtime:    static_zero_alloc_workspaces={}, est_evals/step={}, est_flops/step={}, workspace={} B",
            self.static_zero_alloc_workspaces,
            self.total_estimated_evals_per_step,
            self.total_estimated_flops_per_step,
            self.total_workspace_bytes
        );

        let _ = writeln!(out, "Waves (topological Tarjan DAG):");
        for (w, members) in self.waves.iter().enumerate() {
            let note = if members.len() > 1 {
                " [parallel-independent]"
            } else {
                ""
            };
            let _ = writeln!(out, "  wave {w}: blocks {members:?}{note}");
        }

        let _ = writeln!(out, "Blocks:");
        for block in &self.blocks {
            let _ = writeln!(
                out,
                "  block {} [wave {}, deps {:?} -> {:?}]: nodes={:?} class={} mode={} tear={}",
                block.index,
                block.wave,
                block.dependencies,
                block.dependents,
                block.nodes,
                block.class,
                block.solve_mode,
                block.tear_strategy,
            );
            let dt_str = match block.limit_dt_s {
                Some(dt) => format!("{dt:.6e} s"),
                None => "refused".to_owned(),
            };
            let _ = writeln!(
                out,
                "    unknowns={} (kept={}, substituted={}, stage_z={}) caches(affine={}, affine_line={}) timescale={} (limit_dt={})",
                block.total_unknowns,
                block.kept_unknowns,
                block.substituted_count,
                block.stage_unknowns,
                block.affine_cached,
                block.affine_line_cached,
                block.time_scale_tier,
                dt_str,
            );
            if let Some(ln) = &block.load_line_node {
                let _ = writeln!(out, "    load_line_node={ln}");
            }
            if block.looped {
                let _ = writeln!(out, "    tear_variables={:?}", block.tear_variables);
            }
            let _ = writeln!(
                out,
                "    cost: evals/step={}, factor_flops={}, solve_flops={}, workspace={} B",
                block.cost.estimated_evals_per_step,
                block.cost.estimated_dense_flops_per_factor,
                block.cost.estimated_substitution_flops_per_solve,
                block.cost.preallocated_workspace_bytes,
            );
        }

        let _ = writeln!(out, "Algebraic Loops:");
        if self.loop_count == 0 {
            let _ = writeln!(out, "  none (acyclic graph)");
        } else {
            for block in self.blocks.iter().filter(|b| b.looped) {
                let _ = writeln!(
                    out,
                    "  block {}: {} of {} unknowns kept via {}; tear={:?}",
                    block.index,
                    block.kept_unknowns,
                    block.total_unknowns,
                    block.tear_strategy,
                    block.tear_variables,
                );
            }
        }

        let _ = writeln!(
            out,
            "Time-Scale Groups (probe_dt={:.6e} s):",
            self.probe_dt_s
        );
        for group in &self.time_scale_groups {
            let _ = writeln!(
                out,
                "  {}: blocks {:?}, dt in [{:.6e}, {:.6e}] s, stiffness_ratio={:.3}x",
                group.tier, group.blocks, group.min_dt_s, group.max_dt_s, group.stiffness_ratio,
            );
        }
        out
    }

    /// Deterministic canonical JSON representation suitable for storing in
    /// `fs-ledger` (`EXECUTION_PLAN_ARTIFACT_KIND`) or emitting from CLI tools.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::new();
        out.push('{');
        push_json_kv_str(&mut out, "schema", EXECUTION_PLAN_SCHEMA, true);
        push_json_kv_str(&mut out, "plan_hash", &self.plan_hash_hex(), false);
        push_json_kv_str(&mut out, "problem_id", &self.problem_id_hex(), false);
        let _ = write!(
            out,
            ",\"node_count\":{},\"connection_count\":{},\"boundary_count\":{},\"block_count\":{},\"loop_count\":{},\"static_zero_alloc_workspaces\":{},\"probe_dt_s\":",
            self.node_count,
            self.connection_count,
            self.boundary_count,
            self.block_count,
            self.loop_count,
            self.static_zero_alloc_workspaces,
        );
        push_json_f64(&mut out, self.probe_dt_s);
        let _ = write!(
            out,
            ",\"total_estimated_evals_per_step\":{},\"total_estimated_flops_per_step\":{},\"total_workspace_bytes\":{},\"waves\":[",
            self.total_estimated_evals_per_step,
            self.total_estimated_flops_per_step,
            self.total_workspace_bytes,
        );
        for (wi, wave) in self.waves.iter().enumerate() {
            if wi > 0 {
                out.push(',');
            }
            push_json_usize_array(&mut out, wave);
        }
        out.push_str("],\"blocks\":[");
        for (bi, b) in self.blocks.iter().enumerate() {
            if bi > 0 {
                out.push(',');
            }
            out.push('{');
            let _ = write!(
                out,
                "\"index\":{},\"wave\":{},\"dependencies\":",
                b.index, b.wave
            );
            push_json_usize_array(&mut out, &b.dependencies);
            out.push_str(",\"dependents\":");
            push_json_usize_array(&mut out, &b.dependents);
            out.push_str(",\"nodes\":");
            push_json_str_array(&mut out, &b.nodes);
            out.push_str(",\"node_classes\":[");
            for (ci, class) in b.node_classes.iter().enumerate() {
                if ci > 0 {
                    out.push(',');
                }
                push_json_str(&mut out, class_str(*class));
            }
            out.push(']');
            push_json_kv_str(&mut out, "class", b.class.as_str(), false);
            push_json_kv_str(&mut out, "solve_mode", b.solve_mode.as_str(), false);
            push_json_kv_str(&mut out, "tear_strategy", b.tear_strategy.as_str(), false);
            let _ = write!(
                out,
                ",\"looped\":{},\"affine_cached\":{},\"affine_line_cached\":{},\"load_line_node\":",
                b.looped, b.affine_cached, b.affine_line_cached,
            );
            match &b.load_line_node {
                Some(name) => push_json_str(&mut out, name),
                None => out.push_str("null"),
            }
            let _ = write!(
                out,
                ",\"total_unknowns\":{},\"kept_unknowns\":{},\"substituted_count\":{},\"stage_unknowns\":{},\"tear_variables\":",
                b.total_unknowns, b.kept_unknowns, b.substituted_count, b.stage_unknowns,
            );
            push_json_str_array(&mut out, &b.tear_variables);
            out.push_str(",\"substituted_variables\":");
            push_json_str_array(&mut out, &b.substituted_variables);
            push_json_kv_str(
                &mut out,
                "time_scale_tier",
                b.time_scale_tier.as_str(),
                false,
            );
            out.push_str(",\"limit_dt_s\":");
            match b.limit_dt_s {
                Some(dt) => push_json_f64(&mut out, dt),
                None => out.push_str("null"),
            }
            let _ = write!(
                out,
                ",\"cost\":{{\"estimated_evals_per_step\":{},\"estimated_dense_flops_per_factor\":{},\"estimated_substitution_flops_per_solve\":{},\"preallocated_workspace_bytes\":{}}}}}",
                b.cost.estimated_evals_per_step,
                b.cost.estimated_dense_flops_per_factor,
                b.cost.estimated_substitution_flops_per_solve,
                b.cost.preallocated_workspace_bytes,
            );
        }
        out.push_str("],\"time_scale_groups\":[");
        for (gi, g) in self.time_scale_groups.iter().enumerate() {
            if gi > 0 {
                out.push(',');
            }
            out.push('{');
            push_json_kv_str(&mut out, "tier", g.tier.as_str(), true);
            out.push_str(",\"blocks\":");
            push_json_usize_array(&mut out, &g.blocks);
            out.push_str(",\"min_dt_s\":");
            push_json_f64(&mut out, g.min_dt_s);
            out.push_str(",\"max_dt_s\":");
            push_json_f64(&mut out, g.max_dt_s);
            out.push_str(",\"stiffness_ratio\":");
            push_json_f64(&mut out, g.stiffness_ratio);
            out.push('}');
        }
        out.push_str("]}");
        out
    }
}

/// Determination status for a causalization audit export (aligned with
/// Machine-IR `DeterminationClass` without making `fs-ir` a runtime dependency).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CausalDetermination {
    /// Sealed graph has zero scalar input equations.
    Empty,
    /// Every scalar input has exactly one admitted driving equation and zero
    /// unmatched equations or variables.
    WellDetermined,
}

impl CausalDetermination {
    /// Canonical lowercase identifier.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Empty => "empty",
            Self::WellDetermined => "well_determined",
        }
    }
}

/// One equation-to-variable structural matching entry in [`CausalizationAudit`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalMatchingEntry {
    /// Global scalar input equation index.
    pub equation_index: usize,
    /// Matched scalar variable (`node.port` or `node.port[k]`).
    pub variable: String,
    /// Block index owning this equation.
    pub block_index: usize,
    /// Driver equation family (`"linear"`, `"on_line"`,
    /// `"equivalent_temperature"`, `"boundary"`, `"fixed"`).
    pub driver_kind: &'static str,
    /// Whether this variable participates in an algebraic loop.
    pub is_algebraic_unknown: bool,
    /// Whether this variable is kept in the tear set (versus substituted).
    pub is_tear_variable: bool,
    /// Source output variables read by this equation.
    pub sources: Vec<String>,
}

/// Summary of one algebraic coupling loop in [`CausalizationAudit`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlgebraicLoopAudit {
    /// Block index in execution order.
    pub block_index: usize,
    /// Member component names.
    pub nodes: Vec<String>,
    /// Block mathematical classification.
    pub class: BlockClass,
    /// Selected tearing strategy.
    pub tear_strategy: TearStrategy,
    /// Total scalar unknowns in the loop.
    pub unknown_count: usize,
    /// Kept tear unknowns solved simultaneously.
    pub kept_count: usize,
    /// Torn variable labels.
    pub tear_variables: Vec<String>,
    /// Substituted follower labels in evaluation order.
    pub substituted_variables: Vec<String>,
}

/// Authority-free structural causalization audit export aligned with Machine-IR
/// `CausalizationReceiptDraft` (`fs-ir/src/machine/causalization.rs`).
///
/// Machine-IR is deliberately kept out of the `fs-system` runtime dependency
/// graph (layer L3 must not depend on L6 `fs-ir`); higher layers can ingest
/// this deterministic audit struct or its canonical JSON directly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalizationAudit {
    /// Audit schema version (`1`).
    pub schema_version: u32,
    /// Problem identity hex digest.
    pub problem_id_hex: String,
    /// Execution plan hex digest.
    pub plan_hash_hex: String,
    /// Determination axis.
    pub determination: CausalDetermination,
    /// Total admitted scalar input equations.
    pub equation_count: usize,
    /// Total admitted scalar input variables.
    pub variable_count: usize,
    /// Unmatched equation count (`0` for an admitted sealed [`System`]).
    pub unmatched_equations: usize,
    /// Unmatched variable count (`0` for an admitted sealed [`System`]).
    pub unmatched_variables: usize,
    /// Complete equation-to-variable matching in deterministic index order.
    pub matching: Vec<CausalMatchingEntry>,
    /// Reported algebraic loops and their tear sets.
    pub algebraic_loops: Vec<AlgebraicLoopAudit>,
}

impl CausalizationAudit {
    /// Deterministic canonical JSON representation suitable for ledger storage
    /// or external Machine-IR receipt bridging.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::new();
        out.push('{');
        push_json_kv_str(&mut out, "schema", CAUSALIZATION_AUDIT_SCHEMA, true);
        let _ = write!(out, ",\"schema_version\":{}", self.schema_version);
        push_json_kv_str(&mut out, "problem_id", &self.problem_id_hex, false);
        push_json_kv_str(&mut out, "plan_hash", &self.plan_hash_hex, false);
        push_json_kv_str(
            &mut out,
            "determination",
            self.determination.as_str(),
            false,
        );
        let _ = write!(
            out,
            ",\"equation_count\":{},\"variable_count\":{},\"unmatched_equations\":{},\"unmatched_variables\":{},\"matching\":[",
            self.equation_count,
            self.variable_count,
            self.unmatched_equations,
            self.unmatched_variables,
        );
        for (i, m) in self.matching.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push('{');
            let _ = write!(out, "\"equation_index\":{}", m.equation_index);
            push_json_kv_str(&mut out, "variable", &m.variable, false);
            let _ = write!(out, ",\"block_index\":{}", m.block_index);
            push_json_kv_str(&mut out, "driver_kind", m.driver_kind, false);
            let _ = write!(
                out,
                ",\"is_algebraic_unknown\":{},\"is_tear_variable\":{},\"sources\":",
                m.is_algebraic_unknown, m.is_tear_variable,
            );
            push_json_str_array(&mut out, &m.sources);
            out.push('}');
        }
        out.push_str("],\"algebraic_loops\":[");
        for (i, lp) in self.algebraic_loops.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push('{');
            let _ = write!(out, "\"block_index\":{},\"nodes\":", lp.block_index);
            push_json_str_array(&mut out, &lp.nodes);
            push_json_kv_str(&mut out, "class", lp.class.as_str(), false);
            push_json_kv_str(&mut out, "tear_strategy", lp.tear_strategy.as_str(), false);
            let _ = write!(
                out,
                ",\"unknown_count\":{},\"kept_count\":{},\"tear_variables\":",
                lp.unknown_count, lp.kept_count,
            );
            push_json_str_array(&mut out, &lp.tear_variables);
            out.push_str(",\"substituted_variables\":");
            push_json_str_array(&mut out, &lp.substituted_variables);
            out.push('}');
        }
        out.push_str("]}");
        out
    }
}

impl System {
    /// Emit the content-addressed, explainable static [`ExecutionPlan`] using a
    /// default 1.0 s horizon for time-scale grouping.
    #[must_use]
    pub fn execution_plan(&self) -> ExecutionPlan {
        self.execution_plan_for_dt(1.0)
    }

    /// Emit the content-addressed, explainable static [`ExecutionPlan`],
    /// querying each component's `limit_dt_s` against `probe_dt_s` [s] to
    /// classify blocks into time-scale groups.
    ///
    /// Non-positive or non-finite `probe_dt_s` falls back to `1.0` s.
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn execution_plan_for_dt(&self, probe_dt_s: f64) -> ExecutionPlan {
        let probe_dt = if probe_dt_s.is_finite() && probe_dt_s > 0.0 {
            probe_dt_s
        } else {
            1.0
        };
        let block_count = self.blocks.len();
        let mut block_of = vec![0_usize; self.nodes.len()];
        for (b_idx, block) in self.blocks.iter().enumerate() {
            for &node in &block.nodes {
                block_of[node] = b_idx;
            }
        }
        let output_owner = self.build_output_owners();

        // Build inter-block DAG edges and parallel execution waves.
        let mut dependencies: Vec<Vec<usize>> = vec![Vec::new(); block_count];
        let mut dependents: Vec<Vec<usize>> = vec![Vec::new(); block_count];
        let mut sources = Vec::new();
        for (input_idx, driver) in self.drivers.iter().enumerate() {
            let dst_block = block_of[self.input_owner[input_idx]];
            sources.clear();
            driver.sources(&mut sources);
            for &src in &sources {
                let src_block = block_of[output_owner[src]];
                if src_block != dst_block {
                    dependencies[dst_block].push(src_block);
                    dependents[src_block].push(dst_block);
                }
            }
        }
        for list in dependencies.iter_mut().chain(dependents.iter_mut()) {
            list.sort_unstable();
            list.dedup();
        }
        let mut block_wave = vec![0_usize; block_count];
        let mut max_wave = 0_usize;
        for b in 0..block_count {
            let w = dependencies[b]
                .iter()
                .map(|&dep| block_wave[dep] + 1)
                .max()
                .unwrap_or(0);
            block_wave[b] = w;
            max_wave = max_wave.max(w);
        }
        let mut waves = if block_count == 0 {
            Vec::new()
        } else {
            vec![Vec::new(); max_wave + 1]
        };
        for (b, &w) in block_wave.iter().enumerate() {
            waves[w].push(b);
        }

        // Probe time-scale bounds per block at the accepted state.
        let mut block_dt_limits: Vec<Option<f64>> = Vec::with_capacity(block_count);
        for block in &self.blocks {
            let mut min_dt = probe_dt;
            let mut valid = true;
            for &node in &block.nodes {
                match self.nodes[node].limit_dt_s(self.accepted[node].as_ref(), probe_dt) {
                    Ok(dt) if dt.is_finite() && dt > 0.0 => {
                        min_dt = min_dt.min(dt);
                    }
                    _ => {
                        valid = false;
                        break;
                    }
                }
            }
            block_dt_limits.push(valid.then_some(min_dt));
        }
        let slowest_dt = block_dt_limits
            .iter()
            .flatten()
            .copied()
            .reduce(f64::max)
            .unwrap_or(probe_dt);

        let mut block_plans = Vec::with_capacity(block_count);
        for (b_idx, block) in self.blocks.iter().enumerate() {
            let work = &self.block_work[b_idx];
            let line = (self.lines[b_idx].len() == 1).then(|| &self.lines[b_idx][0]);
            let node_names: Vec<String> =
                block.nodes.iter().map(|&n| self.names[n].clone()).collect();
            let node_classes: Vec<Class> =
                block.nodes.iter().map(|&n| self.nodes[n].class()).collect();
            let class = classify_block(&node_classes, !self.lines[b_idx].is_empty());
            let solve_mode = self.select_block_mode(b_idx);
            let looped = block.looped();
            let total_unknowns = block.unknowns.len();
            let stage_unknowns: usize = block
                .nodes
                .iter()
                .map(|&n| self.nodes[n].stage_unknowns())
                .sum();
            let load_line_node = self.lines[b_idx]
                .first()
                .map(|l| self.names[l.node].clone());

            let tear_strategy = if !looped {
                TearStrategy::Acyclic
            } else if matches!(solve_mode, BlockSolveMode::MonolithicNewtonKrylov) {
                TearStrategy::MonolithicStageResidual
            } else if matches!(
                solve_mode,
                BlockSolveMode::FixedPointAitken | BlockSolveMode::FixedPointIqnIls
            ) {
                TearStrategy::FixedPointAccelerated
            } else if work.line_tear.is_some() {
                TearStrategy::LoadLineWoodbury
            } else {
                TearStrategy::EliminationCondensed
            };

            let (kept_unknowns, substituted_count, tear_variables, substituted_variables) =
                match tear_strategy {
                    TearStrategy::Acyclic => (0, 0, Vec::new(), Vec::new()),
                    TearStrategy::LoadLineWoodbury => {
                        let tear = work
                            .line_tear
                            .as_ref()
                            .expect("LoadLineWoodbury requires line_tear");
                        let line_ref = line.expect("LoadLineWoodbury requires load line");
                        let mut tear_vars: Vec<String> = work
                            .line_condensed
                            .kept_unknowns()
                            .into_iter()
                            .map(|u| self.scalar_input_label(block.unknowns[u]))
                            .collect();
                        let out_base = self.output_start[line_ref.node];
                        for &local_out in tear.active_outputs() {
                            tear_vars.push(self.scalar_output_label(out_base + local_out));
                        }
                        let sub_vars: Vec<String> = work
                            .line_condensed
                            .follower_order()
                            .into_iter()
                            .map(|u| self.scalar_input_label(block.unknowns[u]))
                            .collect();
                        let kept = work.line_condensed.kept() + tear.width();
                        let sub = work.line_condensed.follow_count();
                        (kept, sub, tear_vars, sub_vars)
                    }
                    TearStrategy::EliminationCondensed => {
                        let tear_vars: Vec<String> = work
                            .condensed
                            .kept_unknowns()
                            .into_iter()
                            .map(|u| self.scalar_input_label(block.unknowns[u]))
                            .collect();
                        let sub_vars: Vec<String> = work
                            .condensed
                            .follower_order()
                            .into_iter()
                            .map(|u| self.scalar_input_label(block.unknowns[u]))
                            .collect();
                        (
                            work.condensed.kept(),
                            work.condensed.follow_count(),
                            tear_vars,
                            sub_vars,
                        )
                    }
                    TearStrategy::MonolithicStageResidual => {
                        let mut tear_vars = Vec::with_capacity(stage_unknowns + total_unknowns);
                        for &n in &block.nodes {
                            let nz = self.nodes[n].stage_unknowns();
                            for k in 0..nz {
                                tear_vars.push(format!("{}.z[{k}]", self.names[n]));
                            }
                        }
                        tear_vars.extend(
                            block
                                .unknowns
                                .iter()
                                .map(|&idx| self.scalar_input_label(idx)),
                        );
                        (stage_unknowns + total_unknowns, 0, tear_vars, Vec::new())
                    }
                    TearStrategy::FixedPointAccelerated => {
                        let tear_vars: Vec<String> = block
                            .unknowns
                            .iter()
                            .map(|&idx| self.scalar_input_label(idx))
                            .collect();
                        (total_unknowns, 0, tear_vars, Vec::new())
                    }
                };

            let ws_bytes = work.workspace_bytes(self.group_count)
                + self.lines[b_idx]
                    .iter()
                    .map(Line::workspace_bytes)
                    .sum::<usize>();
            let cost = estimate_block_cost(
                block.nodes.len(),
                total_unknowns,
                kept_unknowns,
                stage_unknowns,
                tear_strategy,
                work.affine,
                work.affine_line,
                work.condensed.nnz(),
                work.line_condensed.kept(),
                work.line_condensed.nnz(),
                work.line_tear.as_ref().map_or(0, LineTear::width),
                work.line_tear
                    .as_ref()
                    .map_or(0, LineTear::line_column_count),
                ws_bytes,
            );

            let limit_dt = block_dt_limits[b_idx];
            let time_scale_tier = classify_time_scale(limit_dt, probe_dt, slowest_dt);

            block_plans.push(BlockPlan {
                index: b_idx,
                wave: block_wave[b_idx],
                dependencies: dependencies[b_idx].clone(),
                dependents: dependents[b_idx].clone(),
                nodes: node_names,
                node_classes,
                class,
                solve_mode,
                tear_strategy,
                looped,
                affine_cached: looped && work.affine,
                affine_line_cached: looped && work.affine_line,
                load_line_node,
                total_unknowns,
                kept_unknowns,
                substituted_count,
                stage_unknowns,
                tear_variables,
                substituted_variables,
                time_scale_tier,
                limit_dt_s: limit_dt,
                cost,
            });
        }

        let time_scale_groups = build_time_scale_groups(&block_plans, probe_dt, slowest_dt);
        let loop_count = block_plans.iter().filter(|b| b.looped).count();
        let total_estimated_evals_per_step = block_plans
            .iter()
            .map(|b| b.cost.estimated_evals_per_step)
            .sum();
        let total_estimated_flops_per_step = block_plans
            .iter()
            .map(|b| {
                b.cost
                    .estimated_dense_flops_per_factor
                    .saturating_add(b.cost.estimated_substitution_flops_per_solve)
            })
            .sum();
        let total_workspace_bytes = block_plans
            .iter()
            .map(|b| b.cost.preallocated_workspace_bytes)
            .sum();

        let problem_id = *self.problem_id().as_bytes();
        let plan_hash = compute_plan_hash(&problem_id, self.solver_policy, &waves, &block_plans);

        ExecutionPlan {
            plan_hash,
            problem_id,
            node_count: self.nodes.len(),
            connection_count: self.connections.len(),
            boundary_count: self.boundaries.len(),
            block_count,
            loop_count,
            static_zero_alloc_workspaces: true,
            probe_dt_s: probe_dt,
            blocks: block_plans,
            waves,
            time_scale_groups,
            total_estimated_evals_per_step,
            total_estimated_flops_per_step,
            total_workspace_bytes,
        }
    }

    /// Emit a deterministic, authority-free [`CausalizationAudit`] describing
    /// every admitted scalar equation, its causal matching, and every algebraic
    /// loop with its tear set.
    #[must_use]
    pub fn causalization_audit(&self) -> CausalizationAudit {
        let plan = self.execution_plan();
        let mut block_of = vec![0_usize; self.nodes.len()];
        for (b_idx, block) in self.blocks.iter().enumerate() {
            for &node in &block.nodes {
                block_of[node] = b_idx;
            }
        }
        // Track which global scalar inputs are algebraic unknowns and tear vars.
        let total_inputs = self.drivers.len();
        let mut is_unknown = vec![false; total_inputs];
        let mut is_tear = vec![false; total_inputs];
        for (b_idx, block) in self.blocks.iter().enumerate() {
            let work = &self.block_work[b_idx];
            let bp = &plan.blocks[b_idx];
            for &idx in &block.unknowns {
                is_unknown[idx] = true;
            }
            match bp.tear_strategy {
                TearStrategy::EliminationCondensed => {
                    for u in work.condensed.kept_unknowns() {
                        is_tear[block.unknowns[u]] = true;
                    }
                }
                TearStrategy::LoadLineWoodbury => {
                    for u in work.line_condensed.kept_unknowns() {
                        is_tear[block.unknowns[u]] = true;
                    }
                }
                TearStrategy::FixedPointAccelerated | TearStrategy::MonolithicStageResidual => {
                    for &idx in &block.unknowns {
                        is_tear[idx] = true;
                    }
                }
                TearStrategy::Acyclic => {}
            }
        }

        let mut sources_buf = Vec::new();
        let mut matching = Vec::with_capacity(total_inputs);
        for (eq_idx, driver) in self.drivers.iter().enumerate() {
            sources_buf.clear();
            driver.sources(&mut sources_buf);
            let sources = sources_buf
                .iter()
                .map(|&src| self.scalar_output_label(src))
                .collect();
            matching.push(CausalMatchingEntry {
                equation_index: eq_idx,
                variable: self.scalar_input_label(eq_idx),
                block_index: block_of[self.input_owner[eq_idx]],
                driver_kind: driver_kind_str(driver),
                is_algebraic_unknown: is_unknown[eq_idx],
                is_tear_variable: is_tear[eq_idx],
                sources,
            });
        }

        let algebraic_loops = plan
            .blocks
            .iter()
            .filter(|b| b.looped)
            .map(|b| AlgebraicLoopAudit {
                block_index: b.index,
                nodes: b.nodes.clone(),
                class: b.class,
                tear_strategy: b.tear_strategy,
                unknown_count: b.total_unknowns,
                kept_count: b.kept_unknowns,
                tear_variables: b.tear_variables.clone(),
                substituted_variables: b.substituted_variables.clone(),
            })
            .collect();

        CausalizationAudit {
            schema_version: 1,
            problem_id_hex: plan.problem_id_hex(),
            plan_hash_hex: plan.plan_hash_hex(),
            determination: if total_inputs == 0 {
                CausalDetermination::Empty
            } else {
                CausalDetermination::WellDetermined
            },
            equation_count: total_inputs,
            variable_count: total_inputs,
            unmatched_equations: 0,
            unmatched_variables: 0,
            matching,
            algebraic_loops,
        }
    }

    fn build_output_owners(&self) -> Vec<usize> {
        let total_outputs = self.outputs.len();
        let mut owners = Vec::with_capacity(total_outputs);
        for node in 0..self.nodes.len() {
            let len = self.output_start[node + 1] - self.output_start[node];
            owners.extend(core::iter::repeat_n(node, len));
        }
        owners
    }

    fn scalar_input_label(&self, global_input: usize) -> String {
        let node = self.input_owner[global_input];
        let local = global_input - self.input_start[node];
        let ports = self.nodes[node].ports();
        for (port_idx, spec) in ports.iter().enumerate() {
            let range = self.layouts[node].input_range(port_idx);
            if range.contains(&local) {
                let offset = local - range.start;
                if range.len() == 1 {
                    return format!("{}.{}", self.names[node], spec.id().as_str());
                }
                return format!("{}.{}[{offset}]", self.names[node], spec.id().as_str());
            }
        }
        self.input_label[global_input].clone()
    }

    fn scalar_output_label(&self, global_output: usize) -> String {
        let node = self
            .output_start
            .partition_point(|&start| start <= global_output)
            .saturating_sub(1);
        let local = global_output - self.output_start[node];
        let ports = self.nodes[node].ports();
        for (port_idx, spec) in ports.iter().enumerate() {
            let range = self.layouts[node].output_range(port_idx);
            if range.contains(&local) {
                let offset = local - range.start;
                if range.len() == 1 {
                    return format!("{}.{}", self.names[node], spec.id().as_str());
                }
                return format!("{}.{}[{offset}]", self.names[node], spec.id().as_str());
            }
        }
        format!("{}.out[{local}]", self.names[node])
    }
}

fn classify_block(classes: &[Class], has_load_line: bool) -> BlockClass {
    if has_load_line || classes.contains(&Class::NonSmooth) {
        return BlockClass::NonSmooth;
    }
    if classes.iter().all(|&c| c == Class::ExactLinear) {
        return BlockClass::ExactLinear;
    }
    let has_smooth = classes.contains(&Class::Smooth);
    let has_monotone = classes.contains(&Class::Monotone1d);
    match (has_smooth, has_monotone) {
        (true, false) => BlockClass::Smooth,
        (false, true) => BlockClass::Monotone1d,
        _ => BlockClass::Mixed,
    }
}

fn classify_time_scale(limit_dt: Option<f64>, probe_dt: f64, slowest_dt: f64) -> TimeScaleTier {
    let Some(dt) = limit_dt else {
        return TimeScaleTier::Unconstrained;
    };
    if dt >= probe_dt * (1.0 - 1.0e-12) {
        return TimeScaleTier::Unconstrained;
    }
    let ratio = slowest_dt / dt;
    if ratio >= 100.0 {
        TimeScaleTier::Fast
    } else if ratio >= 10.0 {
        TimeScaleTier::Medium
    } else {
        TimeScaleTier::Slow
    }
}

fn build_time_scale_groups(
    blocks: &[BlockPlan],
    probe_dt: f64,
    slowest_dt: f64,
) -> Vec<TimeScaleGroup> {
    let tiers = [
        TimeScaleTier::Fast,
        TimeScaleTier::Medium,
        TimeScaleTier::Slow,
        TimeScaleTier::Unconstrained,
    ];
    let mut groups = Vec::new();
    for tier in tiers {
        let mut indices = Vec::new();
        let mut min_dt = f64::INFINITY;
        let mut max_dt = 0.0_f64;
        for b in blocks {
            if b.time_scale_tier == tier {
                indices.push(b.index);
                let dt = b.limit_dt_s.unwrap_or(probe_dt);
                min_dt = min_dt.min(dt);
                max_dt = max_dt.max(dt);
            }
        }
        if !indices.is_empty() {
            let ratio = if min_dt > 0.0 && min_dt.is_finite() {
                (slowest_dt / min_dt).max(1.0)
            } else {
                1.0
            };
            groups.push(TimeScaleGroup {
                tier,
                blocks: indices,
                min_dt_s: min_dt,
                max_dt_s: max_dt,
                stiffness_ratio: ratio,
            });
        }
    }
    groups
}

#[allow(clippy::too_many_arguments)]
fn estimate_block_cost(
    node_count: usize,
    total_unknowns: usize,
    kept_unknowns: usize,
    stage_unknowns: usize,
    tear_strategy: TearStrategy,
    affine_cached: bool,
    affine_line_cached: bool,
    condensed_nnz: usize,
    line_kept: usize,
    line_nnz: usize,
    tear_width: usize,
    line_columns: usize,
    preallocated_workspace_bytes: usize,
) -> BlockCostEstimate {
    let nodes = node_count as u64;
    match tear_strategy {
        TearStrategy::Acyclic => BlockCostEstimate {
            estimated_evals_per_step: nodes,
            estimated_dense_flops_per_factor: 0,
            estimated_substitution_flops_per_solve: 0,
            preallocated_workspace_bytes,
        },
        TearStrategy::LoadLineWoodbury => {
            let k_n = line_kept as u64;
            let w = tear_width as u64;
            let nnz_n = line_nnz as u64;
            let lc = line_columns as u64;
            let w2 = w.saturating_mul(w);
            let w_lu = (2_u64.saturating_mul(w2).saturating_mul(w) / 3).max(w2);
            let k_n2 = k_n.saturating_mul(k_n);
            let n_lu = 2_u64.saturating_mul(k_n2).saturating_mul(k_n) / 3;
            let factor_flops = if affine_line_cached {
                w_lu
            } else {
                n_lu.saturating_add(w_lu)
            };
            let solve_flops = 2_u64
                .saturating_mul(nnz_n)
                .saturating_add(2_u64.saturating_mul(k_n2))
                .saturating_add(2_u64.saturating_mul(w).saturating_mul(lc))
                .saturating_add(2_u64.saturating_mul(w2));
            let evals = if affine_line_cached {
                nodes.saturating_mul(2)
            } else {
                nodes
                    .saturating_mul(3)
                    .saturating_add(total_unknowns as u64)
            };
            BlockCostEstimate {
                estimated_evals_per_step: evals,
                estimated_dense_flops_per_factor: factor_flops,
                estimated_substitution_flops_per_solve: solve_flops,
                preallocated_workspace_bytes,
            }
        }
        TearStrategy::EliminationCondensed => {
            let k = kept_unknowns as u64;
            let nnz = condensed_nnz as u64;
            let k2 = k.saturating_mul(k);
            let factor_flops = if affine_cached {
                0
            } else {
                (2_u64.saturating_mul(k2).saturating_mul(k) / 3).max(k)
            };
            let solve_flops = 2_u64
                .saturating_mul(nnz)
                .saturating_add(2_u64.saturating_mul(k2));
            let evals = if affine_cached {
                nodes.saturating_mul(2)
            } else {
                nodes
                    .saturating_mul(3)
                    .saturating_add(total_unknowns as u64)
            };
            BlockCostEstimate {
                estimated_evals_per_step: evals,
                estimated_dense_flops_per_factor: factor_flops,
                estimated_substitution_flops_per_solve: solve_flops,
                preallocated_workspace_bytes,
            }
        }
        TearStrategy::MonolithicStageResidual => {
            let n = (stage_unknowns as u64).saturating_add(total_unknowns as u64);
            BlockCostEstimate {
                estimated_evals_per_step: nodes.saturating_mul(6),
                estimated_dense_flops_per_factor: 0,
                estimated_substitution_flops_per_solve: 4_u64.saturating_mul(n),
                preallocated_workspace_bytes,
            }
        }
        TearStrategy::FixedPointAccelerated => {
            let m = total_unknowns as u64;
            BlockCostEstimate {
                estimated_evals_per_step: nodes.saturating_mul(6),
                estimated_dense_flops_per_factor: 0,
                estimated_substitution_flops_per_solve: 4_u64.saturating_mul(m),
                preallocated_workspace_bytes,
            }
        }
    }
}

fn push_field(preimage: &mut Vec<u8>, bytes: &[u8]) {
    preimage.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    preimage.extend_from_slice(bytes);
}

fn push_u64(preimage: &mut Vec<u8>, value: u64) {
    push_field(preimage, &value.to_le_bytes());
}

fn compute_plan_hash(
    problem_id: &[u8; 32],
    policy: BlockSolverPolicy,
    waves: &[Vec<usize>],
    blocks: &[BlockPlan],
) -> [u8; 32] {
    let mut preimage = Vec::new();
    push_field(&mut preimage, problem_id);
    encode_policy(&mut preimage, policy);
    push_u64(&mut preimage, waves.len() as u64);
    for wave in waves {
        push_u64(&mut preimage, wave.len() as u64);
        for &b in wave {
            push_u64(&mut preimage, b as u64);
        }
    }
    push_u64(&mut preimage, blocks.len() as u64);
    for b in blocks {
        push_u64(&mut preimage, b.index as u64);
        push_u64(&mut preimage, b.wave as u64);
        push_u64(&mut preimage, b.dependencies.len() as u64);
        for &d in &b.dependencies {
            push_u64(&mut preimage, d as u64);
        }
        push_u64(&mut preimage, b.nodes.len() as u64);
        for name in &b.nodes {
            push_field(&mut preimage, name.as_bytes());
        }
        for &c in &b.node_classes {
            push_field(&mut preimage, class_str(c).as_bytes());
        }
        push_field(&mut preimage, b.class.as_str().as_bytes());
        push_field(&mut preimage, b.solve_mode.as_str().as_bytes());
        push_field(&mut preimage, b.tear_strategy.as_str().as_bytes());
        push_field(&mut preimage, &[u8::from(b.looped)]);
        push_field(&mut preimage, &[u8::from(b.affine_cached)]);
        push_field(&mut preimage, &[u8::from(b.affine_line_cached)]);
        match &b.load_line_node {
            Some(name) => push_field(&mut preimage, name.as_bytes()),
            None => push_field(&mut preimage, &[]),
        }
        push_u64(&mut preimage, b.total_unknowns as u64);
        push_u64(&mut preimage, b.kept_unknowns as u64);
        push_u64(&mut preimage, b.substituted_count as u64);
        push_u64(&mut preimage, b.stage_unknowns as u64);
        push_u64(&mut preimage, b.tear_variables.len() as u64);
        for v in &b.tear_variables {
            push_field(&mut preimage, v.as_bytes());
        }
        push_u64(&mut preimage, b.substituted_variables.len() as u64);
        for v in &b.substituted_variables {
            push_field(&mut preimage, v.as_bytes());
        }
        push_u64(&mut preimage, b.cost.estimated_evals_per_step);
        push_u64(&mut preimage, b.cost.estimated_dense_flops_per_factor);
        push_u64(&mut preimage, b.cost.estimated_substitution_flops_per_solve);
        push_u64(&mut preimage, b.cost.preallocated_workspace_bytes as u64);
    }
    *fs_blake3::hash_domain(EXECUTION_PLAN_DOMAIN, &preimage).as_bytes()
}

fn encode_policy(preimage: &mut Vec<u8>, policy: BlockSolverPolicy) {
    match policy {
        BlockSolverPolicy::Auto => push_field(preimage, b"auto"),
        BlockSolverPolicy::InterfaceNewton => push_field(preimage, b"interface_newton"),
        BlockSolverPolicy::FixedPointAitken(cfg) => {
            push_field(preimage, b"fixed_point_aitken");
            push_u64(preimage, cfg.initial_relaxation.to_bits());
            push_u64(preimage, cfg.min_relaxation.to_bits());
            push_u64(preimage, cfg.max_relaxation.to_bits());
        }
        BlockSolverPolicy::FixedPointIqnIls(cfg) => {
            push_field(preimage, b"fixed_point_iqn_ils");
            push_u64(preimage, cfg.max_history as u64);
            push_u64(preimage, cfg.relative_rank_tolerance.to_bits());
            push_field(preimage, &[u8::from(cfg.reuse_history_across_steps)]);
            push_u64(preimage, cfg.fallback.initial_relaxation.to_bits());
            push_u64(preimage, cfg.fallback.min_relaxation.to_bits());
            push_u64(preimage, cfg.fallback.max_relaxation.to_bits());
        }
        BlockSolverPolicy::MonolithicNewtonKrylov(cfg) => {
            push_field(preimage, b"monolithic_newton_krylov");
            push_u64(preimage, cfg.absolute_tolerance.to_bits());
            push_u64(preimage, cfg.relative_tolerance.to_bits());
            push_u64(preimage, cfg.max_iterations as u64);
            push_u64(preimage, cfg.linear_restart as u64);
            push_u64(preimage, cfg.max_linear_cycles as u64);
        }
        BlockSolverPolicy::ByTier {
            fixed_point,
            monolithic,
        } => {
            push_field(preimage, b"by_tier");
            match fixed_point {
                FixedPointMethod::Aitken(cfg) => {
                    push_field(preimage, b"aitken");
                    push_u64(preimage, cfg.initial_relaxation.to_bits());
                    push_u64(preimage, cfg.min_relaxation.to_bits());
                    push_u64(preimage, cfg.max_relaxation.to_bits());
                }
                FixedPointMethod::IqnIls(cfg) => {
                    push_field(preimage, b"iqn_ils");
                    push_u64(preimage, cfg.max_history as u64);
                    push_u64(preimage, cfg.relative_rank_tolerance.to_bits());
                    push_field(preimage, &[u8::from(cfg.reuse_history_across_steps)]);
                    push_u64(preimage, cfg.fallback.initial_relaxation.to_bits());
                    push_u64(preimage, cfg.fallback.min_relaxation.to_bits());
                    push_u64(preimage, cfg.fallback.max_relaxation.to_bits());
                }
            }
            push_u64(preimage, monolithic.absolute_tolerance.to_bits());
            push_u64(preimage, monolithic.relative_tolerance.to_bits());
            push_u64(preimage, monolithic.max_iterations as u64);
            push_u64(preimage, monolithic.linear_restart as u64);
            push_u64(preimage, monolithic.max_linear_cycles as u64);
        }
    }
}

fn class_str(class: Class) -> &'static str {
    match class {
        Class::ExactLinear => "exact_linear",
        Class::Smooth => "smooth",
        Class::Monotone1d => "monotone_1d",
        Class::NonSmooth => "non_smooth",
    }
}

fn driver_kind_str(driver: &Driver) -> &'static str {
    match driver {
        Driver::Unset => "unset",
        Driver::Linear(_) => "linear",
        Driver::OnLine { .. } => "on_line",
        Driver::EquivalentTemperature { .. } => "equivalent_temperature",
        Driver::Boundary { .. } => "boundary",
        Driver::Fixed(_) => "fixed",
    }
}

fn hex_encode(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        let _ = write!(s, "{b:02x}");
    }
    s
}

fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn push_json_kv_str(out: &mut String, key: &str, val: &str, first: bool) {
    if !first {
        out.push(',');
    }
    push_json_str(out, key);
    out.push(':');
    push_json_str(out, val);
}

fn push_json_f64(out: &mut String, value: f64) {
    if value.is_finite() {
        let _ = write!(out, "{value:.17e}");
    } else {
        out.push_str("null");
    }
}

fn push_json_usize_array(out: &mut String, values: &[usize]) {
    out.push('[');
    for (i, &v) in values.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{v}");
    }
    out.push(']');
}

fn push_json_str_array(out: &mut String, values: &[String]) {
    out.push('[');
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_json_str(out, v);
    }
    out.push(']');
}
