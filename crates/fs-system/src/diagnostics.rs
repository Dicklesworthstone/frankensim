//! Solve-mode selection, explicit step budgets, typed partial outcomes, and
//! deterministic per-step JSONL diagnostics for the coupled-step driver.

#![deny(unsafe_code)]

use core::fmt;

use fs_ports::component::Receipt;

use crate::accelerator::{AitkenConfig, IqnIlsConfig};

/// Solve mode chosen for one block of the coupled evaluation plan.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockSolveMode {
    /// Acyclic block evaluated once in topological order.
    OrderedSweep,
    /// Fixed-point iteration with bounded Aitken dynamic under-relaxation on
    /// torn port variables (mandatory `Component` tier only).
    FixedPointAitken,
    /// Fixed-point iteration with bounded vector IQN-ILS (Anderson acceleration
    /// with twice-orthogonalized MGS rank filtering) on torn port variables.
    FixedPointIqnIls,
    /// Condensed interface Newton using `PortTangent` or resolved finite
    /// differences.
    InterfaceNewton,
    /// Load-line Newton / Woodbury tear for a block containing a `LoadLine`
    /// component.
    LoadLineNewton,
    /// Monolithic Newton--Krylov through `fs-solver`'s fallible adapter when
    /// every component in the looped block offers `StageResidual`.
    MonolithicNewtonKrylov,
}

impl BlockSolveMode {
    /// Canonical stable mode label for diagnostics and JSONL output.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::OrderedSweep => "ordered_sweep",
            Self::FixedPointAitken => "fixed_point_aitken",
            Self::FixedPointIqnIls => "fixed_point_iqn_ils",
            Self::InterfaceNewton => "interface_newton",
            Self::LoadLineNewton => "load_line_newton",
            Self::MonolithicNewtonKrylov => "monolithic_newton_krylov",
        }
    }
}

impl fmt::Display for BlockSolveMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Configuration for monolithic Newton--Krylov solves of `StageResidual` blocks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MonolithicConfig {
    /// Absolute Euclidean residual tolerance for the monolithic vector `[r_z, r_u]`.
    pub absolute_tolerance: f64,
    /// Relative residual tolerance against initial residual norm (`0.0` disables).
    pub relative_tolerance: f64,
    /// Maximum outer Newton iterations per block trial.
    pub max_iterations: usize,
    /// Inner FGMRES restart length.
    pub linear_restart: usize,
    /// Maximum FGMRES cycles per outer Newton iteration.
    pub max_linear_cycles: usize,
}

impl Default for MonolithicConfig {
    fn default() -> Self {
        Self {
            absolute_tolerance: 1.0e-14,
            relative_tolerance: 0.0,
            max_iterations: 25,
            linear_restart: 30,
            max_linear_cycles: 4,
        }
    }
}

/// Fixed-point acceleration method selected by [`BlockSolverPolicy::ByTier`]
/// for mandatory-tier-only looped blocks without `PortTangent` or `StageResidual`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FixedPointMethod {
    /// Bounded Aitken dynamic under-relaxation.
    Aitken(AitkenConfig),
    /// Bounded vector IQN-ILS with rank-filtered MGS history.
    IqnIls(IqnIlsConfig),
}

/// Policy controlling how looped blocks select their solve mode.
///
/// [`Self::Auto`] selects the highest applicable tier per block while keeping
/// full compatibility with existing `InterfaceNewton` and `LoadLineNewton`
/// graphs. [`Self::ByTier`] additionally routes mandatory-tier-only looped
/// blocks (no `PortTangent`, no `StageResidual`, no `LoadLine`) through
/// derivative-free fixed-point acceleration (`Aitken` or `IQN-ILS`). Explicit
/// variants force a specific mode for benchmarking or cross-mode convergence
/// testing.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BlockSolverPolicy {
    /// Automatic tier-based selection preserving v0 default behavior:
    /// - acyclic block -> `OrderedSweep`
    /// - block with `LoadLine` -> `LoadLineNewton`
    /// - looped block where every node offers `StageResidual` -> `MonolithicNewtonKrylov`
    /// - looped block otherwise -> `InterfaceNewton`
    #[default]
    Auto,
    /// Full tier-based dispatch:
    /// - acyclic block -> `OrderedSweep`
    /// - block with `LoadLine` -> `LoadLineNewton`
    /// - looped block where every node offers `StageResidual` -> `MonolithicNewtonKrylov`
    /// - looped block where at least one node offers `PortTangent` -> `InterfaceNewton`
    /// - looped block with mandatory `Component` tier only -> fixed point (`fixed_point`)
    ByTier {
        /// Accelerator used for mandatory-tier-only looped blocks.
        fixed_point: FixedPointMethod,
        /// Monolithic Newton--Krylov settings for `StageResidual` blocks.
        monolithic: MonolithicConfig,
    },
    /// Force bounded Aitken fixed-point iteration on looped blocks without a load line.
    FixedPointAitken(AitkenConfig),
    /// Force bounded vector IQN-ILS fixed-point iteration on looped blocks without a load line.
    FixedPointIqnIls(IqnIlsConfig),
    /// Force condensed interface Newton on looped blocks without a load line.
    InterfaceNewton,
    /// Force monolithic Newton--Krylov via `StageResidual` on looped blocks.
    MonolithicNewtonKrylov(MonolithicConfig),
}

/// Explicit resource budgets for one coupled step.
///
/// Every field defaults to `None` (unlimited) so warm unconstrained steps incur
/// zero wall-clock overhead. When a budget is reached, the step refuses
/// atomically with [`crate::SystemError::BudgetExhausted`] carrying a
/// [`PartialStepOutcome`], leaving accepted physics, ledger, and controller
/// state bitwise untouched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct StepBudgets {
    /// Maximum component evaluations (`propose`, `jvp`, `stage_residual`, `jvp_z`).
    pub max_evaluations: Option<u64>,
    /// Maximum weighted work units (evaluations weighted by component/block size).
    pub max_work_units: Option<u64>,
    /// Maximum wall-clock duration for the step, in nanoseconds.
    pub max_wall_ns: Option<u64>,
}

/// Which resource budget was exhausted during a coupled step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BudgetKind {
    /// Component evaluation count limit reached.
    Evaluations,
    /// Weighted work-unit limit reached.
    Work,
    /// Wall-clock time limit reached.
    WallTime,
}

impl BudgetKind {
    /// Canonical stable label for JSONL diagnostics.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Evaluations => "evaluations",
            Self::Work => "work",
            Self::WallTime => "wall_time",
        }
    }
}

impl fmt::Display for BudgetKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Typed partial outcome reported when a coupled step exhausts an explicit
/// budget before convergence or publication.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartialStepOutcome {
    /// Which budget limit triggered the stop.
    pub budget: BudgetKind,
    /// Total component evaluations consumed during the attempted step.
    pub spent_evaluations: u64,
    /// Total weighted work units consumed during the attempted step.
    pub spent_work_units: u64,
    /// Elapsed wall-clock time in nanoseconds (zero when wall budget is disabled).
    pub elapsed_wall_ns: u64,
    /// Index of the block executing when the budget ran out.
    pub block_index: usize,
    /// Solve mode active on that block.
    pub mode: BlockSolveMode,
    /// Solver iterations completed on that block before halting.
    pub iterations: u32,
    /// Infinity norm of the last observed interface residual (`NaN` if none yet).
    pub worst_residual: f64,
}

/// Per-block diagnostic entry for one coupled step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlockStepDiagnostic {
    /// Block index in topological execution order.
    pub block: usize,
    /// Solve mode used for this block.
    pub mode: BlockSolveMode,
    /// Solver iterations performed on this block.
    pub iterations: u32,
    /// Component evaluations performed on this block.
    pub evaluations: u64,
    /// Final Euclidean residual norm of the block's interface equations.
    pub residual_l2: f64,
    /// Final infinity residual norm of the block's interface equations.
    pub residual_linf: f64,
}

/// Per-component receipt summary in declared component order.
#[derive(Debug, Clone, PartialEq)]
pub struct ReceiptDiagnostic {
    /// Declared component name.
    pub node: String,
    /// Stored-energy change reported by the component [J].
    pub stored_change_j: f64,
    /// Physical dissipation reported by the component [J].
    pub dissipated_j: f64,
    /// Numerical loss reported by the component [J].
    pub numerical_loss_j: f64,
    /// Entropy generated reported by the component [J/K].
    pub entropy_generated_j_per_k: f64,
}

impl ReceiptDiagnostic {
    /// Build from a node name and its step [`Receipt`].
    #[must_use]
    pub fn from_receipt(node: String, receipt: Receipt) -> Self {
        Self {
            node,
            stored_change_j: receipt.stored_change_j,
            dissipated_j: receipt.dissipated_j,
            numerical_loss_j: receipt.numerical_loss_j,
            entropy_generated_j_per_k: receipt.entropy_generated_j_per_k,
        }
    }
}

/// Final verdict of a diagnosed coupled step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepVerdict {
    /// All blocks converged, audits and admissions passed, and states published.
    Accepted,
    /// A step budget ran out before publication.
    BudgetExhausted {
        /// Which budget ran out.
        budget: BudgetKind,
    },
    /// Cooperative `Cx` cancellation interrupted the step.
    Cancelled,
    /// Another typed driver, audit, or component refusal occurred.
    Refused {
        /// Human-readable refusal description.
        reason: String,
    },
}

/// Complete per-step diagnostic record with deterministic JSONL serialization.
#[derive(Debug, Clone, PartialEq)]
pub struct StepDiagnostic {
    /// Accepted step number after publication (or last accepted step on refusal).
    pub step: u64,
    /// Accepted physical time [s].
    pub time_s: f64,
    /// Attempted or accepted step duration `dt` [s].
    pub dt_s: f64,
    /// Summary solve mode for the step (single block's mode, or `"mixed"`).
    pub mode: &'static str,
    /// Total nonlinear/fixed-point iterations summed across blocks.
    pub iterations: u32,
    /// Number of component-requested shortenings before this trial.
    pub shortenings: u32,
    /// Total component evaluations consumed.
    pub evaluations: u64,
    /// Total weighted work units consumed.
    pub work_units: u64,
    /// Aggregate Euclidean residual norm across blocks.
    pub residual_l2: f64,
    /// Aggregate infinity residual norm across blocks.
    pub residual_linf: f64,
    /// Per-block solve telemetry in execution order.
    pub blocks: Vec<BlockStepDiagnostic>,
    /// Per-component receipts in declared order.
    pub receipts: Vec<ReceiptDiagnostic>,
    /// Step verdict.
    pub verdict: StepVerdict,
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
                use core::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn push_json_f64(out: &mut String, v: f64) {
    use core::fmt::Write as _;
    if v.is_finite() {
        let normalized = if v == 0.0 { 0.0 } else { v };
        let _ = write!(out, "{normalized:.17e}");
    } else {
        out.push_str("null");
    }
}

impl StepDiagnostic {
    /// Render this diagnostic record as one deterministic single-line JSONL object.
    ///
    /// Includes `dt`, `mode`, `iterations`, `residual_norms` (`l2`, `linf`),
    /// `blocks`, `receipts` (in declared component order), and `verdict`.
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::with_capacity(256);
        out.push_str("{\"step\":");
        let _ = write!(out, "{}", self.step);
        out.push_str(",\"time_s\":");
        push_json_f64(&mut out, self.time_s);
        out.push_str(",\"dt\":");
        push_json_f64(&mut out, self.dt_s);
        out.push_str(",\"mode\":");
        push_json_str(&mut out, self.mode);
        out.push_str(",\"iterations\":");
        let _ = write!(out, "{}", self.iterations);
        out.push_str(",\"shortenings\":");
        let _ = write!(out, "{}", self.shortenings);
        out.push_str(",\"evaluations\":");
        let _ = write!(out, "{}", self.evaluations);
        out.push_str(",\"work_units\":");
        let _ = write!(out, "{}", self.work_units);
        out.push_str(",\"residual_norms\":{\"l2\":");
        push_json_f64(&mut out, self.residual_l2);
        out.push_str(",\"linf\":");
        push_json_f64(&mut out, self.residual_linf);
        out.push_str("},\"blocks\":[");
        for (i, b) in self.blocks.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str("{\"block\":");
            let _ = write!(out, "{}", b.block);
            out.push_str(",\"mode\":");
            push_json_str(&mut out, b.mode.as_str());
            out.push_str(",\"iterations\":");
            let _ = write!(out, "{}", b.iterations);
            out.push_str(",\"evaluations\":");
            let _ = write!(out, "{}", b.evaluations);
            out.push_str(",\"residual_l2\":");
            push_json_f64(&mut out, b.residual_l2);
            out.push_str(",\"residual_linf\":");
            push_json_f64(&mut out, b.residual_linf);
            out.push('}');
        }
        out.push_str("],\"receipts\":[");
        for (i, r) in self.receipts.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str("{\"node\":");
            push_json_str(&mut out, &r.node);
            out.push_str(",\"stored_change_j\":");
            push_json_f64(&mut out, r.stored_change_j);
            out.push_str(",\"dissipated_j\":");
            push_json_f64(&mut out, r.dissipated_j);
            out.push_str(",\"numerical_loss_j\":");
            push_json_f64(&mut out, r.numerical_loss_j);
            out.push_str(",\"entropy_generated_j_per_k\":");
            push_json_f64(&mut out, r.entropy_generated_j_per_k);
            out.push('}');
        }
        out.push_str("],\"verdict\":");
        match &self.verdict {
            StepVerdict::Accepted => out.push_str("\"accepted\""),
            StepVerdict::Cancelled => out.push_str("\"cancelled\""),
            StepVerdict::BudgetExhausted { budget } => {
                out.push_str("{\"budget_exhausted\":");
                push_json_str(&mut out, budget.as_str());
                out.push('}');
            }
            StepVerdict::Refused { reason } => {
                out.push_str("{\"refused\":");
                push_json_str(&mut out, reason);
                out.push('}');
            }
        }
        out.push('}');
        out
    }
}
