//! fs-adjoint — GRADIENT TRUTH (plan §8.7, Bet 4's FLUX realization).
//! Layer: L3 FLUX.
//!
//! Discrete adjoints EVERYWHERE, never differentiated through Krylov
//! iterations: the implicit-function theorem turns "differentiate the
//! solver" into ONE transposed solve sharing the primal's
//! preconditioner infrastructure (fs-solver's `apply_transpose` is in
//! the operator trait for exactly this). Time-dependent adjoints run
//! reverse sweeps under fs-ad's binomial revolve checkpointing. Shape
//! derivatives come in both honest forms: Hadamard boundary integrals
//! on FEEC traces (the mathematically clean form) and density/SIMP
//! volumetric chain rules (where boundary forms are awkward). Sobolev
//! (H¹) smoothing is the Riesz-representation step that turns raw
//! mesh-noise-amplifying gradients into smooth descent directions —
//! "the single most important practical trick in shape optimization".
//!
//! The verification doctrine: every gradient this crate produces is
//! checked by [`verify::verify_gradient`] against central finite
//! differences (and, where forward duals reach, against fs-ad's
//! `Dual64`) — the gate ci-gauntlet wires so that a solver without a
//! passing gradient check cannot merge.

pub mod hadamard;
pub mod ift;
pub mod sobolev;
pub mod timedep;
pub mod verify;

pub use hadamard::{compliance_shape_gradient, volume_shape_gradient};
pub use ift::{AdjointReport, DensityPoisson, ift_gradient_matfree};
pub use sobolev::sobolev_smooth;
pub use timedep::{HeatAdjoint, heat_initial_gradient};
pub use verify::{GradientVerdict, verify_gradient};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
