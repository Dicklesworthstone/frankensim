//! fs-time — structure-preserving time integration (plan §8.5). Layer:
//! L3 FLUX.
//!
//! Integrators that preserve what the physics preserves: symplectic
//! (Störmer–Verlet, with its discrete-Lagrangian equivalence documented
//! and tested), Lie-group updates over canonical `fs-ga` `So3`/`Se3` states,
//! generalized-α with CONTROLLABLE dissipation, operator-
//! backed first/second-order implicit lanes, IMEX and exponential integrators
//! for stiffness, and embedded-pair adaptivity with a PI controller.
//!
//! The two universal obligations (P7 + §8.7): RESUMABLE state machines
//! (checkpoint = clone; split runs bitwise-equal to straight runs) and
//! DISCRETE ADJOINTS of the stepper (Verlet's ships here, checkpointed
//! through fs-ad's revolve; it is the template for the rest).
//! Hybrid-time/event-order/Zeno object semantics are isolated in [`hybrid`];
//! their validated receipts prove schema consistency, never a Zeno theorem.

pub mod adaptive;
pub mod galpha;
pub mod hybrid;
pub mod lie;
pub mod se3;
#[cfg(feature = "time-slabs")]
pub mod slabs;
pub mod stiff;
pub mod symplectic;

pub use adaptive::{AdaptiveState, PiController, rk45_adaptive};
pub use fs_ga::{Se3, So3, Twist, Vec3};
pub use galpha::{
    FirstOrderGeneralizedAlpha, FirstOrderProblem, FirstOrderState, GeneralizedAlpha,
    ImplicitSolveConfig, ImplicitStepTelemetry, LinearFirstOrderSystem, LinearSecondOrderSystem,
    OperatorFirstOrderGeneralizedAlpha, OperatorGeneralizedAlpha, SecondOrderProblem,
    SecondOrderState, TimeSolveError, first_order_galpha_step, galpha_step,
};
pub use lie::{LieStepError, rigid_body_step, so3_body_exp_step, so3_space_exp_step};
pub use se3::{
    BalanceReceipt, DepSolveParams, DepStepReceipt, RattleProjection, RenormPolicy, RenormReceipt,
    Se3ClaimClass, Se3Error, Se3FixtureDeclaration, Unconstrained, claim_for, dep_free_step,
    dep_momentum_adjoint, run_dep_free, se3_exp_step, se3_exp_step_renorm, se3_rigid_body_step,
    se3_space_exp_step,
};
pub use stiff::{
    ExpEuler, IdentityPreconditioner, Imex2, ImexSolveConfig, ImexSolveError, ImexStage, ImexState,
    ImexStepTelemetry, OperatorImex2, imex2_step,
};
pub use symplectic::{verlet_adjoint, verlet_step};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
