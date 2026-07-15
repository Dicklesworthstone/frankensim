//! fs-time — structure-preserving time integration (plan §8.5). Layer:
//! L3 FLUX.
//!
//! Integrators that preserve what the physics preserves: symplectic
//! (Störmer–Verlet, with its discrete-Lagrangian equivalence documented
//! and tested), Lie-group SO(3) quaternion and SE(3) motor updates via the
//! exponential map, generalized-α with CONTROLLABLE dissipation, operator-
//! backed first/second-order implicit lanes, IMEX and exponential integrators
//! for stiffness, and embedded-pair adaptivity with a PI controller.
//!
//! The two universal obligations (P7 + §8.7): RESUMABLE state machines
//! (checkpoint = clone; split runs bitwise-equal to straight runs) and
//! DISCRETE ADJOINTS of the stepper (Verlet's ships here, checkpointed
//! through fs-ad's revolve; it is the template for the rest).

pub mod adaptive;
pub mod galpha;
pub mod lie;
pub mod se3;
#[cfg(feature = "time-slabs")]
pub mod slabs;
pub mod stiff;
pub mod symplectic;

pub use adaptive::{AdaptiveState, PiController, rk45_adaptive};
pub use galpha::{
    FirstOrderGeneralizedAlpha, FirstOrderProblem, FirstOrderState, GeneralizedAlpha,
    ImplicitSolveConfig, ImplicitStepTelemetry, LinearFirstOrderSystem, LinearSecondOrderSystem,
    OperatorFirstOrderGeneralizedAlpha, OperatorGeneralizedAlpha, SecondOrderProblem,
    SecondOrderState, TimeSolveError, first_order_galpha_step, galpha_step,
};
pub use lie::{quat_exp, quat_exp_step, quat_mul, quat_rotate, rigid_body_step};
pub use se3::{
    BalanceReceipt, DepSolveParams, DepStepReceipt, RattleProjection, RenormPolicy, RenormReceipt,
    Se3ClaimClass, Se3Error, Se3FixtureDeclaration, Twist, Unconstrained, canonicalize_motor,
    claim_for, dep_free_step, dep_momentum_adjoint, run_dep_free, se3_exp_step,
    se3_exp_step_renorm, se3_rigid_body_step,
};
pub use stiff::{
    ExpEuler, IdentityPreconditioner, Imex2, ImexSolveConfig, ImexSolveError, ImexStage, ImexState,
    ImexStepTelemetry, OperatorImex2, imex2_step,
};
pub use symplectic::{verlet_adjoint, verlet_step};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
