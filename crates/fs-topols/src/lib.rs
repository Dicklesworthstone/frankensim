//! fs-topols — level-set topology optimization (plan §9.5 [S/F], bead
//! 7tv.12): shape-gradient velocity advection with TOPOLOGICAL
//! DERIVATIVES for hole nucleation — genuine topology changes with
//! mathematical justification rather than heuristic hole-punching.
//!
//! Layer: L4 (ASCENT). The level set IS the geometry: [`GridSdf`]
//! implements fs-cutfem's `CutSdf` with an EXACT bilinear enclosure,
//! so fs-cutfem's canonical elasticity operator evaluates directly on the evolving
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
//! - [`optimize`]: the single-load compliance descent loop — CutFEM solve on the
//!   SDF, energy-density shape velocity with augmented-Lagrangian
//!   volume control, extend → advect → redistance → audit → (maybe)
//!   nucleate, everything ledgered.
//! - [`checkpoint`]: exact durable single-load continuation from retained level
//!   set, global iteration ordinal, and augmented-Lagrange multiplier without
//!   replaying earlier geometry updates.
//! - [`robust_descent`]: evaluated-state simultaneous independent-load compliance descent;
//!   equilibrium is solved per scenario before shape fields are aggregated, and
//!   complete evolved candidates are re-solved before trajectory publication.
//! - [`evaluated`]: transactional publication boundary that independently
//!   re-solves the exact returned geometry so final compliance, area and
//!   snapshot are bound to one actually evaluated design.
//! - [`guarded`]: bounded whole-trajectory candidate acceptance using only
//!   independently evaluated final objectives and material areas.
//! - [`robust`]: independent multi-load final evaluation and transactional
//!   robust selection across bounded candidate trajectories.
//! - [`stress`]: independent sampled plane-strain von Mises evaluation plus
//!   transactional volume/stress-limited candidate publication.
//! - [`robust_stress`]: worst-scenario sampled-stress admission across
//!   simultaneous independent load cases with robust final replay.
//! - [`volume`]: bounded hard-area projection with prescribed fixed nodes and
//!   transactional cancellation; uses the optimizer's numerical area functional.
//! - [`projected`]: feasible-baseline hard-area descent, independently re-solved
//!   accepted geometry and exact accepted-state continuation.
//! - [`projected_stress`]: sampled-stress admission at every same-area update,
//!   bounded candidate retries and exact continuation of a feasible design.

pub mod checkpoint;
pub mod design_regions;
pub mod evaluated;
pub mod fim;
pub mod gridsdf;
pub mod guarded;
pub mod optimize;
pub mod projected;
pub mod projected_stress;
pub mod refinement;
pub mod robust;
#[path = "robust_descent_v2.rs"]
pub mod robust_descent;
pub mod robust_stress;
pub mod stress;
pub mod topder;
pub mod veloext;
pub mod volume;
pub mod weno;

pub use checkpoint::OptimizeCheckpoint;
pub use evaluated::{
    EvaluatedFinalState, EvaluatedOptimizeReport, evaluate_compliance_design,
    optimize_compliance_evaluated,
};
pub use fim::{RedistanceAudit, hausdorff, redistance, zero_crossings};
pub use gridsdf::GridSdf;
pub use guarded::{
    GuardedCandidate, GuardedOptimizeReport, GuardedSettings, GuardedStop,
    optimize_compliance_guarded,
};
pub use optimize::{Cantilever, OptimizeReport, OptimizeSettings, optimize_compliance};
pub use projected_stress::{ProjectedStressCheck, ProjectedStressOptimizer, ProjectedStressUpdate};
pub use robust::{
    RobustAggregate, RobustCandidate, RobustEvaluation, RobustLoadCase, RobustOptimizeReport,
    RobustStop, evaluate_robust_design, optimize_compliance_robust_guarded,
};
pub use robust_descent::{
    RobustDescentReport, optimize_compliance_multi_load, optimize_compliance_multi_load_guarded,
};
pub use robust_stress::{
    RobustSampledStressEvaluation, RobustStressCandidate, RobustStressReport, RobustStressStop,
    evaluate_robust_sampled_stress, optimize_compliance_multi_load_stress_guarded,
};
pub use stress::{
    SampledStressEvaluation, SampledStressLimit, StressGuardedCandidate, StressGuardedReport,
    StressGuardedStop, evaluate_sampled_stress, optimize_compliance_stress_guarded,
};
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