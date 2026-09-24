//! fs-verify — the certified-speculation VERIFIER (bead lmp4.1, the
//! addendum's SINGLE research bet). Layer: L3.
//!
//! The verifier is what lets untrusted proposers be maximally
//! aggressive: the accept test is CERTIFIED, so correctness never
//! depends on the learned component. For the elliptic class,
//! equilibrated-flux a-posteriori estimators (Prager–Synge) give
//! GUARANTEED, CONSTANT-FREE upper bounds on the energy-norm error —
//! the hard analysis is purchased from the literature, and the
//! remaining risk (floating-point rounding) is retired by
//! outward-rounded interval evaluation over mathematically exact
//! quadrature. An accepted candidate carries a VERIFIED color; a
//! rejected or unbounded evaluation carries NOTHING (fail closed).
//!
//! The original 1D polynomial-manufactured class remains available. `tet`
//! adds domain-conditional 3D RT0 energy majorants for linear scalar diffusion
//! with mixed Dirichlet, Neumann and Robin data. It does not promote a
//! geometry, material or product capability's evidence authority. The nonlinear
//! class gets the honest fallback: candidates are WARM STARTS with
//! measured iteration savings, never certificates.

/// Actual steady FEM primal/dual solves with a scoped continuum mean bound.
#[cfg(feature = "thermal-conduction")]
pub mod conduction;
#[cfg(feature = "certified-speculation")]
pub mod economics;
#[cfg(feature = "certified-speculation")]
pub mod estimator;
#[cfg(feature = "certified-speculation")]
pub mod fem1d;
#[cfg(feature = "certified-speculation")]
pub mod interval;
/// Equilibrated tetrahedral flux reconstruction and outward energy majorants.
#[cfg(feature = "certified-speculation")]
pub mod tet;
#[cfg(feature = "certified-speculation")]
pub mod zoo;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
