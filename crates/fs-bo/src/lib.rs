//! fs-bo — Bayesian optimization (plan §9.4 [F]). Layer: L4 ASCENT.
//!
//! In-house Gaussian processes (Matérn ½/3⁄2/5⁄2 with ARD, exact
//! Cholesky inference, QMC-multistart marginal-likelihood fitting),
//! deterministic acquisitions (closed-form EI; q-EI through the
//! Cholesky reparameterization over a FIXED scrambled-Sobol normal
//! bank), and BO loops whose every run replays bitwise. Inner
//! optimizers come from the landed stack: fs-ascent L-BFGS for
//! hyperparameters, fs-dfo CMA-ES for acquisition surfaces.
//!
//! TuRBO trust-region BO, multi-fidelity cost-aware acquisition, and
//! inducing-point sparse GPs are included; tape-differentiated
//! acquisition gradients remain a recorded follow-up lane.

pub mod acq;
#[cfg(feature = "tape-acq")]
pub mod acq_grad;
pub mod bo;
/// Joint latent feasibility and reference-capped noisy design improvement.
pub mod constrained;
pub mod gp;
/// Bounded kernel learning that preserves declared per-observation noise.
pub mod hyper;
/// Noisy Bayesian optimization with explicitly budgeted kernel refitting.
pub mod learning;
pub mod mf;
/// Joint-posterior noisy expected improvement for stochastic objectives.
pub mod noisy;
pub mod sparse;
pub mod turbo;

pub use acq::{expected_improvement, normal_bank, phi_cdf, phi_inv, q_expected_improvement};
pub use bo::{BoConfig, BoReport, minimize};
pub use gp::{Gp, Kernel, Matern, fit_hyperparams};
pub use mf::{MfConfig, MfGp, MfKernel, MfReport, fit_mf, mf_minimize};
pub use noisy::{joint_normal_bank, q_noisy_expected_improvement};
pub use sparse::{SparseGp, farthest_point_inducing};
pub use turbo::{TurboConfig, TurboReport, turbo_minimize};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
