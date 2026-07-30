//! Runner V2 local work-package declarations and rootless evaluator handoffs.
//!
//! Stage-A declarations in this module are source-authoritative, canonical
//! contract data. Evaluator handoffs are deliberately different: they are
//! bounded Rust values with no canonical root, attempt identity, AC57
//! disposition, retained artifact, receipt, telemetry, or authority.

pub mod handoff;
pub mod work_packages;

pub use handoff::{
    RunnerV2LocalWorkPackageHandoffV1, RunnerV2RawCellObservationV1, RunnerV2RawDiagnosticV1,
    RunnerV2RawOutcomeKindV1, RunnerV2RawReasonV1, RunnerV2SafeNumericObservationV1,
    RunnerV2SafeNumericUnitV1, RunnerV2SafeNumericValueV1,
};
