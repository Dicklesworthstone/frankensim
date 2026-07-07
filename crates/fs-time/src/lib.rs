//! fs-time — structure-preserving time integration (plan §8.5). Layer:
//! L3 FLUX.
//!
//! Integrators that preserve what the physics preserves: symplectic
//! (Störmer–Verlet, with its discrete-Lagrangian equivalence documented
//! and tested), Lie-group SE(3)/SO(3) via exponential-map updates (no
//! renormalization hacks), generalized-α with CONTROLLABLE dissipation,
//! IMEX and exponential integrators for stiffness, and embedded-pair
//! adaptivity with a PI controller.
//!
//! The two universal obligations (P7 + §8.7): RESUMABLE state machines
//! (checkpoint = clone; split runs bitwise-equal to straight runs) and
//! DISCRETE ADJOINTS of the stepper (Verlet's ships here, checkpointed
//! through fs-ad's revolve; it is the template for the rest).

pub mod adaptive;
pub mod galpha;
pub mod lie;
pub mod stiff;
pub mod symplectic;

pub use adaptive::{AdaptiveState, PiController, rk45_adaptive};
pub use galpha::{GeneralizedAlpha, galpha_step};
pub use lie::{quat_exp_step, rigid_body_step};
pub use stiff::{ExpEuler, imex2_step};
pub use symplectic::{verlet_adjoint, verlet_step};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
