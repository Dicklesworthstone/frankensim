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
//!   patterns) and the FNV-64 body hash is the study identity;
//!   parsing REBUILDS through the validating builder, so tampered
//!   files cannot smuggle ill-typed graphs.

mod eval;
mod ir;
mod serial;

pub use eval::{DescentOptions, DescentReport, Value, descend_fn, descend_ir, eval};
pub use ir::{
    Class, Constraint, ConstraintKind, EvalBudget, Expr, Manifold, NodeId, Objective,
    OptError, OptimizerFamily, Problem, ProblemBuilder, ProblemTag, Sense, Shape, VarId,
    Variable,
};
pub use serial::{parse, problem_hash, serialize};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
