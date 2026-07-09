//! fs-ornith — Flagship 1 (plan §15.1, bead mye.2): the ORNITHOID
//! MULTI-INLET AIRCRAFT, smoke tier. Layer: L6 (HELM).
//!
//! Five stages over battle-tested crates: PARAMETERIZE (sectional
//! candidate with wing thickness, trim angle, inlet position and a
//! flapping gait — every lever exposing a Jacobian action, the BEM
//! adjoint where it exists), SCREEN WIDE (fs-bem panel L/D with a
//! documented drag proxy + fs-vpm flapping wake metric, generations
//! E-RACED through fs-race so dominated candidates die early with the
//! payoff MEASURED), REFINE (fs-lbm channel flow around the rasterized
//! section; forces by CONTROL-VOLUME momentum balance over the public
//! moments — panel-vs-LBM agreement gated within model-form evidence),
//! TRIM & CERTIFY (2-state pitch model per candidate; fs-sos Lyapunov
//! certificate → a CERTIFIED region-of-attraction proxy volume;
//! conformal e-band around the surrogate screen), and the PARETO
//! ATLAS (fs-dfo NSGA-II over L/D × ROA × maneuver × inlet mass-flow,
//! gradient polish via the BEM adjoint; every atlas row carries its
//! certificates and lineage).
//!
//! Recorded successors (named, not pretended): F-rep fuselage blends +
//! IGA-shell surfaces + manifold-harmonics global form, cumulant
//! LBM/LES on sparse VDB lattices with DWR adaptivity, Koopman/DMD
//! trim models with conformal e-bands per trim state, NSGA-III
//! reference-point selection, LBM adjoints [M].

pub mod atlas;
pub mod certify;
pub mod param;
pub mod refine;
pub mod screen;

pub use atlas::{Atlas, AtlasRow, build_atlas};
pub use certify::{CertifyReport, LdSurrogate, certify};
pub use param::{GENE_DIM, OrnithCandidate};
pub use refine::{RefineReport, refine};
pub use screen::{ScreenReport, screen_generation};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
