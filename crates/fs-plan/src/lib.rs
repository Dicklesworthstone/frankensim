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
//! Layer: L6 (HELM). Runtime deps: `std`, fs-geom, fs-ledger.

pub mod alloc;
pub mod cost;
pub mod ledgers;
#[cfg(feature = "moonshot-planner")]
pub mod moonshot;
pub mod oracle;
#[cfg(feature = "voi-queries")]
pub mod voi;

pub use alloc::{
    AllocProblem, Allocator, BudgetInfeasible, Knob, KnobSetting, Plan, allocate, oracle_min_error,
};
pub use cost::{CostModel, CostObservation, CostPrediction, CostRefusal, MIN_OBS};
pub use ledgers::{
    Contribution, ErrorLedger, ErrorSource, LedgerDefect, Rigor, TimeLedger, TimeStage,
};
pub use oracle::{PlanCostOracle, cost_model_from_tune, json_f64_field};

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
