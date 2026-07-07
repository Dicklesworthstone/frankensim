//! fs-dfo — derivative-free optimization engines (plan §9.3): CMA-ES in
//! its information-geometric form (natural-gradient flow on the Gaussian
//! family — the framing that DICTATES the step-size/covariance couplings
//! and buys invariance properties by construction), BIPOP restarts, and
//! a Nelder–Mead polish baseline.
//!
//! Layer: L4 ASCENT. Engines are IR-agnostic (closure objectives) by
//! design — routing from the fs-opt problem IR is a small wiring bead
//! once that crate stabilizes (deliberate collision avoidance; see bead
//! 7tv.4's comment trail).
//!
//! DETERMINISM: all sampling flows from keyed Philox streams; ranking
//! uses `total_cmp` with lowest-index tie-breaks — the whole evolution
//! is a pure function of the seed (bitwise rerun-tested, cross-ISA
//! golden-hashed).

pub mod cma;
pub mod neldermead;

pub use cma::{BipopReport, CmaParams, CmaReport, bipop_cmaes, cmaes};
pub use neldermead::nelder_mead;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
