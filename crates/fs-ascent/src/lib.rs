//! fs-ascent — the gradient-based optimizer stack (plan §9.2). Layer:
//! L4 ASCENT.
//!
//! The pipeline made explicit: raw adjoint gradient (fs-adjoint) →
//! Sobolev/Riesz smoothing → optional Riemannian projection
//! (fs-opt's manifold metadata) → optimizer step. Optimizers here are
//! ENGINES; problem structure (typed objective/constraint graphs,
//! manifold metadata) lives in fs-opt — this crate consumes it.
//!
//! House obligations, inherited from the FLUX stack and tested the
//! same way: RESUMABLE states (clone = checkpoint, split runs bitwise
//! equal), deterministic trajectories from seeds (G5), budget-aware
//! stopping through the condition algebra, and a certificate attached
//! to every returned optimum (KKT residuals for constrained solves;
//! gradient norms + stall diagnoses for unconstrained) so converged
//! and stalled are DISTINGUISHABLE outcomes.

pub mod auglag;
pub mod lbfgs;
pub mod riemann;
pub mod stop;
pub mod trust;
pub mod wolfe;

pub use auglag::{AugLagReport, KktResidual, augmented_lagrangian};
pub use lbfgs::{LbfgsReport, LbfgsState};
pub use riemann::{RiemannianLbfgs, RiemannianReport, retract, tangent_project};
pub use stop::{StopReason, StopRule};
pub use trust::{TrustRegionReport, trust_region_newton};
pub use wolfe::{WolfeOutcome, strong_wolfe};

/// The objective callback shape every engine consumes:
/// x ↦ (f(x), ∇f(x)).
pub type FnGrad<'a> = &'a mut dyn FnMut(&[f64]) -> (f64, Vec<f64>);

/// Hessian-vector product callback: (x, v) ↦ H(x)·v.
pub type FnHv<'a> = &'a mut dyn FnMut(&[f64], &[f64]) -> Vec<f64>;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
