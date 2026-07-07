//! PARKED (agent coordination, bead 7tv.1): this battery targeted the
//! numerics-spine agent's IR draft, which lost the crate-structure race
//! to the fuller in-flight implementation (ir/eval/serial modules). The
//! draft modules remain in src/ (graph.rs, manifold.rs, riemann.rs,
//! sexpr.rs — unreferenced by lib.rs, so not compiled) for HARVEST:
//! - manifold.rs: working Sphere/So3 tangent-projection + normalize
//!   retraction (unit-invariant tested to 1e-12 in the draft battery),
//! - riemann.rs: the toy Riemannian descent driver,
//! - sexpr.rs: a RE-VALIDATING parser (rebuilds through typed
//!   constructors so tampered files with dimension violations are
//!   rejected at parse time — a property worth keeping),
//! - graph.rs: hash-consed builder (CSE by construction) + exact
//!   reverse-mode gradient on the DAG with deterministic accumulation.
//! See the bead comment trail for the full coordination note. Delete
//! these parked files only with explicit user permission (RULE 1).
