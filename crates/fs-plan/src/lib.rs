//! fs-plan: per-operator error + cost models and the Error/Time Ledgers
//! (plan §11.4, Bet 12; Decalogue P4). "How accurate is this number and
//! where did the error come from" — and "where did the seconds go" —
//! become QUERIES over attribution trees, not research projects.
//!
//! v1 scope (models; the ALLOCATOR that optimizes over them is gp3.9):
//! - [`cost::CostModel`] — quantile cost predictions (P10/P50/P90) from
//!   observed history via log-log power-law fits with empirical residual
//!   bands; online-updated; refuses instead of guessing on thin data.
//! - [`ledgers::ErrorLedger`] / [`ledgers::TimeLedger`] — attribution
//!   trees with per-source subtotals, dominant-source queries, completeness
//!   lint (no silent error mass), and `explain()` JSON payloads.
//! - [`oracle::PlanCostOracle`] — the Rep Router's `CostOracle` backed by
//!   these models; [`oracle::cost_model_from_tune`] rebuilds models
//!   deterministically from fs-ledger `tune` snapshots.
//!
//! Layer: L6 (HELM). Runtime deps: `std`, fs-geom, fs-ledger; feature-gated
//! VoI adds asupersync cancellation and fs-eproc audit authority.

pub mod alloc;
pub mod cost;
pub mod ledgers;
#[cfg(feature = "moonshot-planner")]
pub mod moonshot;
pub mod oracle;
#[cfg(feature = "voi-queries")]
pub mod voi;

pub use alloc::{
    AllocProblem, AllocationError, Allocator, BudgetInfeasible, Knob, KnobSetting,
    MAX_ALLOCATION_KNOBS, MAX_EXECUTION_TRACKS, MAX_ORACLE_COMBINATIONS, MAX_SETTINGS_PER_KNOB,
    MAX_TOTAL_SETTINGS, Plan, PlanInputError, allocate, oracle_min_error,
};
pub use cost::{
    CostModel, CostObservation, CostPrediction, CostRefusal, MAX_COST_EVALUATION_OBSERVATIONS,
    MAX_COST_OBSERVATIONS, MIN_OBS,
};
pub use ledgers::{
    Contribution, ErrorLedger, ErrorSource, LedgerDefect, Rigor, TimeLedger, TimeLedgerDefect,
    TimeStage,
};
pub use oracle::{
    MAX_PLAN_ORACLE_EDGE_BYTES, MAX_PLAN_ORACLE_EDGES, MAX_PLAN_ORACLE_ERROR_OBSERVATIONS,
    MAX_ROOFLINE_RECEIPT_BYTES, PlanCostOracle, PlanOracleError, ROOFLINE_MACHINE_KEY_BYTES,
    ROOFLINE_RECEIPT_VERSION, ROOFLINE_ROW_SCHEMA, ROOFLINE_TUNE_SHAPE_PREFIX, TuneModelError,
    cost_model_from_tune,
};
#[cfg(feature = "voi-queries")]
pub use voi::{
    MAX_VOI_AUDIT_RECORDS, MAX_VOI_EVALUATIONS, MAX_VOI_GRID, MAX_VOI_NAME_BYTES, MAX_VOI_NODES,
    MAX_VOI_PROBES, MAX_VOI_SCHEDULED_CONTEXTS, MAX_VOI_WORK_UNITS, VOI_AUDIT_ALPHA, VoiError,
};

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
