//! fs-topols — level-set topology optimization (plan §9.5 [S/F], bead
//! 7tv.12): shape-gradient velocity advection with TOPOLOGICAL
//! DERIVATIVES for hole nucleation — genuine topology changes with
//! mathematical justification rather than heuristic hole-punching.
//!
//! Layer: L4 (ASCENT). The level set IS the geometry: [`GridSdf`]
//! implements fs-cutfem's `CutSdf` with an EXACT bilinear enclosure,
//! so fs-solid's CutFEM elasticity evaluates directly on the evolving
//! field — no mesh anywhere in the loop (the marquee coupling).
//!
//! - [`gridsdf`]: the discrete level set (nodal φ, bilinear, certified
//!   per-cell corner enclosures).
//! - [`weno`]: WENO5 + TVD-RK3 narrow-band advection — linear
//!   (order-battery) and Godunov normal-flow (optimizer) Hamiltonians.
//! - [`fim`]: fast-iterative-method redistancing with frozen
//!   interface reconstruction and drift AUDITS (Hausdorff of the zero
//!   set, |∇φ|−1 statistics) — the redistancing-frequency policy's
//!   inputs.
//! - [`veloext`]: interface-normal velocity extension (ascending-|φ|
//!   upwind sweeps) and H¹ smoothing through fs-adjoint's Sobolev
//!   Riesz step.
//! - [`topder`]: the topological derivative of compliance for hole
//!   insertion, with NUMERICALLY GATED constants (nucleation events
//!   predict the compliance change of an actually-punched hole).
//! - [`optimize`]: the compliance descent loop — CutFEM solve on the
//!   SDF, energy-density shape velocity with augmented-Lagrangian
//!   volume control, extend → advect → redistance → audit → (maybe)
//!   nucleate, everything ledgered.

pub mod fim;
pub mod gridsdf;
pub mod optimize;
pub mod topder;
pub mod veloext;
pub mod weno;

pub use fim::{RedistanceAudit, hausdorff, redistance, zero_crossings};
pub use gridsdf::GridSdf;
pub use optimize::{OptimizeReport, OptimizeSettings, optimize_compliance};
pub use topder::{NucleationEvent, nucleate, topological_derivative};
pub use veloext::extend_velocity;
pub use weno::{Velocity, advect, build_band};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
