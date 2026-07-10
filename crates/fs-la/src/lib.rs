//! fs-la — Dense linear algebra: GEMM, batched small dense, factorizations, eigensolvers.
//!
//! Layer: L1. See CONTRACT.md for invariants, error model, determinism
//! class, cancellation behavior, and no-claim boundaries. This crate is part
//! of the FrankenSim workspace; the layer dependency direction is enforced by
//! `cargo run -p xtask -- check-all`.

pub mod batched;
pub mod eigen;
pub mod eigen_complex;
pub mod factor;
pub mod gemm;
pub mod mixed;
pub mod rand_nla;

pub use gemm::{gemm_f32, gemm_f64, gemm_f64_parallel, gemm_mixed};

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
