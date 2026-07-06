//! fs-ivl — certified arithmetic (plan §6.4): outward-rounded intervals and
//! affine forms. This crate is WHAT THE WORD "CERTIFIED" MEANS everywhere
//! else in the system: every operation's postcondition is ENCLOSURE (the
//! result contains the true value set), tested as the G0 containment law.
//!
//! - [`Interval`]: directed rounding via fs-math's `next_up`/`next_down`
//!   nudging — no global rounding-mode state anywhere (grep-lintable:
//!   this workspace never touches the FPU control word). Elementary
//!   functions inherit fs-math's DECLARED ULP budgets.
//! - [`Affine`]/[`AffineCtx`]: noise-symbol forms that kill the dependency
//!   problem on correlated expressions (x − x, deep F-rep DAGs).
//! - The high-precision oracle rungs live in `fs_math::{eft, dd}` (L0;
//!   single implementation shared with fs-la's iterative refinement —
//!   recorded relocation, beads 6ys.8/6ys.12). Quad-double and Taylor
//!   models are recorded follow-up scope.
//!
//! Determinism: everything here is straight-line IEEE arithmetic on
//! fs-math strict functions — cross-ISA bit-deterministic BY CONSTRUCTION
//! (golden-hashed in tests/conformance.rs, verified on both reference
//! ISAs).

pub mod affine;
pub mod interval;

pub use affine::{Affine, AffineCtx};
pub use interval::Interval;

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
