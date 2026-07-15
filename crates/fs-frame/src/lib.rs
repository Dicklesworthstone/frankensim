//! fs-frame — Flagship 2 (plan §15.2, bead mye.3): the SEISMIC-MINIMAL
//! building frame, smoke tier. Layer: L6 HELM.
//!
//! Five stages, each riding a crate that already carries its own
//! battery: LAYOUT (fs-truss ground structure + PDHG diagnostics + outward
//! optimum certificate),
//! SIZING (Euler floors, catalog
//! up-snapping, code rows), TIME HISTORY (fiber-hinge story model —
//! Mander concrete + Menegotto–Pinto steel through fs-solid sections —
//! under fs-scenario Kanai–Tajimi ensembles, Newmark average
//! acceleration), FRAGILITY (anytime-valid e-stopped exceedance
//! probability via fs-eproc confidence sequences + an fs-uq MLMC
//! report), and CVaR MASS MINIMIZATION (Rockafellar–Uryasev over the
//! section scale, then catalog snap with independent re-check).
//!
//! SMOKE TIER, honestly: one story, two fiber-hinge columns, synthetic
//! motions for fragility/CVaR plus one provenance-pinned El Centro diagnostic
//! fixture. The full-resolution lanes (distributed-plasticity
//! frames, recorded-motion suites, million-member ground structures,
//! variational integrators) are recorded successors in the CONTRACT —
//! named, not pretended.

pub mod cvar;
pub mod fragility;
pub mod history;
pub mod layout;

pub use cvar::{CvarDesign, cvar_mass_min, ensemble_cvar};
pub use fragility::{FragilityReport, e_stopped_fragility};
pub use history::{
    GroundMotion, HistoryError, HistoryLimits, HistoryResponse, StoryFrame, StoryParams, peak_drift,
};
pub use layout::{LayoutError, LayoutReport, layout_and_size};

use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};

/// Refuse ensemble payloads whose scalar vector cannot honestly be interpreted
/// as a non-empty Kanai-Tajimi ground-acceleration history.
///
/// The smoke-tier public APIs still use programmer-contract panics for malformed
/// study specifications. Keeping this check here gives every fragility/CVaR
/// entry point the same fail-closed physics boundary until fs-scenario exposes
/// the typed realization artifact tracked by `frankensim-sj31i.39`.
pub(crate) fn assert_ground_motion_ensemble(ensemble: &StochasticEnsemble) {
    assert!(
        ensemble.members > 0,
        "frame studies require at least one ground-motion ensemble member"
    );
    assert!(
        matches!(ensemble.model, SpectrumModel::KanaiTajimi { .. }),
        "frame studies require a Kanai-Tajimi ground-acceleration ensemble; \
         wind spectra and material-parameter bands are not structural motions"
    );
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
