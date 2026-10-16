//! Source-owned Runner V2 foundational work-package declarations.
//!
//! The public `run_*` wrappers are intentionally absent at Stage A. Their
//! implementation, binding, invocation, and attempt evidence belong to the
//! final integration owner.

pub mod base_values;

pub use base_values::{
    RunnerV2BaseValuesStageADeclarationV1, RunnerV2StageADeclarationRootV1, declare_24_1_1_1_1_v1,
};

#[allow(
    unused_imports,
    reason = "the final integration work package will invoke this crate-private evaluator"
)]
pub(crate) use base_values::evaluate_24_1_1_1_1_cell_v1;
