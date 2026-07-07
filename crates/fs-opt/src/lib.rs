//! fs-opt — the optimization problem IR (plan §9.1): problems ARE DATA.
//! Typed objective/constraint expression graphs over manifold-valued
//! variables, built incrementally with validation at every step,
//! serializable (deterministic s-expr), hashable, and exactly
//! differentiable (reverse-mode on the DAG).
//!
//! Layer: L4 ASCENT. The graph substrate only — constraint KIND semantics
//! (kinds, repair) belong to fs-constraint; optimizer drivers to later
//! fs-opt beads.
//!
//! Agent ergonomics (the load-bearing property): every added node is
//! type-checked (dimensions via fs-qty, smoothness class propagation), so
//! "this objective is non-smooth through that min()" is KNOWN at build
//! time and routed to the right optimizer family with a diagnostic that
//! names the offending node.

pub mod graph;
pub mod manifold;
pub mod riemann;
pub mod sexpr;

pub use graph::{DiffClass, NodeId, OptimizerFamily, Problem, ValidationError, VarId};
pub use manifold::Manifold;
pub use riemann::descend;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
