//! fs-solver — the solver stack (plan §8.9). Layer: L3 FLUX.
//!
//! Matrix-free Krylov methods (CG, MINRES, GMRES(m)) bound by the four
//! workspace contract obligations from day one:
//!
//! - RESUMABLE: solver state is plain data — `clone()` is a
//!   checkpoint, and split runs are bitwise-equal to straight runs
//!   (the fs-time `AdaptiveState` pattern, tested the same way).
//! - CANCELLABLE: iteration granularity — every state is complete
//!   between iterations, so drivers interrupt by simply not calling
//!   `step` again (fs-exec Cx wiring is driver scope).
//! - DETERMINISTIC: all inner products go through the fixed-shape
//!   chunked reduction (fs-tilelang's combiner — shape depends on
//!   length only, never on threads or tiers).
//! - ADJOINT-EQUIPPED: the operator trait carries `apply_transpose`,
//!   and transposed solves run through the same machinery (tested to
//!   converge comparably to primal — the IFT contract's enabler).
//!
//! Error transparency: every solve returns a residual HISTORY and a
//! structured stall diagnosis instead of a timeout mystery.

pub mod krylov;
pub mod op;
pub mod pmg;

pub use krylov::{CgState, GmresState, MinresState, SolveReport, StallDiagnosis};
pub use op::{CsrOp, LinearOp};
pub use pmg::{MaskedTensorOp, PMultigrid};

/// Deterministic inner product: elementwise products folded through
/// the fixed-shape chunked combiner (shape = f(length) only).
#[must_use]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "dot length mismatch");
    let prods: Vec<f64> = a.iter().zip(b).map(|(x, y)| x * y).collect();
    fs_tilelang::deterministic_sum(&prods)
}

/// Deterministic 2-norm.
#[must_use]
pub fn norm2(a: &[f64]) -> f64 {
    fs_math::det::sqrt(dot(a, a))
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
