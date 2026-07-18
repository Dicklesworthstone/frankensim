//! fs-opt — the optimization problem IR (plan §9.1). Layer: L4.
//!
//! Optimization problems ARE DATA: typed objective/constraint graphs
//! over MANIFOLD-valued variables — storable, hashable, replayable,
//! and constructible INCREMENTALLY with validation at every step (the
//! agent-ergonomics property):
//!
//! - shapes and fs-qty DIMENSIONS check at build time (adding meters
//!   to square meters is refused with the nodes named);
//! - differentiability CLASS propagates bottom-up, so "this objective
//!   is non-smooth through that min()" is known before any optimizer
//!   runs, and [`Problem::route`] refuses the wrong family with the
//!   poisoning node named;
//! - PDE constraints (`physics(u, θ) = 0`) and stochastic operators
//!   (expectation / CVaR / quantile over UQ configs) are FIRST-CLASS
//!   nodes carrying adjoint/config metadata — the IR represents them,
//!   FLUX/UQ execute them;
//! - manifolds carry retraction metadata ([`Manifold::retract`]) that
//!   a toy Riemannian descent consumes — orientations optimize as
//!   orientations;
//! - canonical serialization round-trips bitwise (floats travel as bit
//!   patterns); admitted semantic identity and exact wire identity are
//!   domain-separated BLAKE3 values, while the legacy FNV-64 body hash
//!   remains a quarantined correlation/corruption tripwire only;
//!   parsing REBUILDS through the validating builder, so tampered
//!   files cannot smuggle ill-typed graphs.
//! - quantified reach-avoid and viability games retain exact quantifier
//!   order, information timing, causal strategy class, model/set context,
//!   horizon, stopping, proof polarity, composition, and no-claim scope in a
//!   separate versioned admission identity ([game]).

mod admission;
mod eval;
pub mod game;
mod guard;
mod ir;
mod serial;

pub use admission::{
    ADMISSION_SCHEMA_VERSION, AdmissionCaps, AdmissionReport, AdmissionViolation, ProblemAdmission,
};
pub use eval::{
    BindingFrame, DescentOptions, DescentReport, DescentStop, RetractionCurve, Value, descend_fn,
    descend_ir, eval, eval_keyed,
};
pub use guard::{
    DeltaPerturbationStep, Endpoint, EscalationKind, EscalationStep, GoodhartGuard, GuardFinding,
    GuardReport, GuardStatus, StepOutcome, StepReport, converged_and_guard_cleared,
};
pub use ir::{
    BilevelRef, Class, Constraint, ConstraintKind, EvalBudget, EvalLimit, Expr,
    MANIFOLD_LAYOUT_SCHEMA_VERSION, Manifold, ManifoldLayout, ManifoldLayoutError, NodeId,
    Objective, ObjectiveEvalSite, OptError, OptimizerFamily, ParamDim, PointDim, ProbeDirection,
    Problem, ProblemBuilder, ProblemTag, Sense, Shape, TangentDim, VarId, Variable, children,
};
pub use serial::{
    ContentHash, DimensionCrosswalkReceipt, FiveToSixRule, LegacyProblemHash, ParsedProblem,
    ProblemSemanticId, WireContentId, WireVersion, canonical_v2_migration_target, parse,
    parse_with_version, problem_hash, serialize, serialize_with_id,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
