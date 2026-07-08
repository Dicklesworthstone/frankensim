//! fs-frame — Flagship 2 (plan §15.2, bead mye.3): the SEISMIC-MINIMAL
//! building frame, smoke tier. Layer: L6 HELM.
//!
//! Five stages, each riding a crate that already carries its own
//! battery: LAYOUT (fs-truss ground structure + PDHG LP with the
//! duality gap as certificate), SIZING (Euler floors, catalog
//! up-snapping, code rows), TIME HISTORY (fiber-hinge story model —
//! Mander concrete + Menegotto–Pinto steel through fs-solid sections —
//! under fs-scenario Kanai–Tajimi ensembles, Newmark average
//! acceleration), FRAGILITY (anytime-valid e-stopped exceedance
//! probability via fs-eproc confidence sequences + an fs-uq MLMC
//! report), and CVaR MASS MINIMIZATION (Rockafellar–Uryasev over the
//! section scale, then catalog snap with independent re-check).
//!
//! SMOKE TIER, honestly: one story, two fiber-hinge columns, synthetic
//! motions only. The full-resolution lanes (distributed-plasticity
//! frames, recorded-motion suites, million-member ground structures,
//! variational integrators) are recorded successors in the CONTRACT —
//! named, not pretended.

pub mod cvar;
pub mod fragility;
pub mod history;
pub mod layout;

pub use cvar::{CvarDesign, cvar_mass_min, ensemble_cvar};
pub use fragility::{FragilityReport, e_stopped_fragility};
pub use history::{StoryFrame, StoryParams, peak_drift};
pub use layout::{LayoutReport, layout_and_size};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
