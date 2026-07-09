//! fs-vessel — Flagship 3 (plan §15.3, bead mye.4): the LAMINAR-POUR
//! VESSEL, smoke tier. Layer: L6 (HELM).
//!
//! Five stages over battle-tested crates: PARAMETERIZE (Chebyshev
//! vessel-of-revolution profile + scalar lip channel), the STABILITY
//! OBJECTIVE (quasi-steady thin-film Reynolds numbers along the pour
//! path feeding fs-cheb Orr–Sommerfeld modal growth — min-max over
//! modes and stations, the differentiable laminarity proxy of
//! Appendix C), VALIDATION (fs-lbm free-surface pours over the lip
//! under a rotating-gravity tilt schedule, strict mass ledger,
//! Carreau viscosity band, Plateau–Rayleigh fragment scoring, and the
//! CONTACT-LINE BRACKET — the genuinely open problem handled by
//! reporting a sensitivity band, never pretended certainty),
//! ROBUSTIFICATION (CVaR over the fluid band, e-raced candidates
//! through fs-race), and the DELIVERABLE (the pour rendered by
//! fs-render's Woodcock tracker bound zero-copy to the simulation's
//! own mass buffer — the marketing shot and the physics are the same
//! bytes).
//!
//! Recorded successors (named, not pretended): level-set lip topology
//! change (fs-topols wiring), cumulant collision, recorded-fluid
//! calibration, the full Appendix C study-program runner.

pub mod pour;
pub mod robust;
pub mod stability;

pub use pour::{PourOutcome, PourRig, render_pour};
pub use robust::{RobustReport, robustify};
pub use stability::{VesselProfile, growth_objective};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
