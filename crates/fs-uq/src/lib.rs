//! fs-uq — uncertainty quantification (plan §8.8). Layer: L4 ASCENT
//! (propagation WRAPS solvers the way ASCENT optimizers do, and its
//! risk outputs feed ASCENT's robust formulations; sitting at L4 also
//! reuses fs-bo's deterministic Φ⁻¹ for QMC normal germs instead of
//! duplicating the polynomial).
//!
//! Slice 1: Karhunen–Loève random fields with CAPTURED-VARIANCE
//! evidence, polynomial chaos by regression with known-answer gates,
//! the QMC workhorse with its MC advantage MEASURED, and multilevel
//! Monte Carlo with the telescoping sum AUDITED. Seismic machinery
//! (Kanai–Tajimi, CQC, IDA fragility) and e-process anytime-valid
//! stopping are the bead's split lanes.

pub mod kl;
pub mod mlmc;
pub mod pce;

pub use kl::{CovarianceKind, KlExpansion};
pub use mlmc::{MlmcReport, mlmc_estimate};
pub use pce::{PceModel, fit_pce};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
