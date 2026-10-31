//! Stage 5: CVaR-constrained MASS MINIMIZATION in the
//! Rockafellar–Uryasev form: CVaR_β(L) = min_t t + E[(L−t)₊]/(1−β),
//! evaluated over the motion ensemble with the section scale s as the
//! design variable. The smoke-tier bisection assumes a monotone feasible
//! interval in section scale; resonance can violate that assumption.
//! Every nonlinear history must converge, and the final catalog-snapped
//! candidate must independently satisfy the risk limit or be refused.
//! Use [building] for material-bearing finite-catalog design with actual
//! physical mass and no size-response monotonicity assumption.

use crate::history::{GroundMotion, HistoryError, HistoryLimits, StoryFrame, StoryParams};
pub use fs_robust::{EmpiricalCvarReport, RobustError, cvar, empirical_cvar};
use fs_scenario::ensemble::{SpectrumModel, StochasticEnsemble};

/// Material-bearing distributed-plasticity building catalog design.
pub mod building;

/// The CVaR design record.
#[derive(Debug, Clone, PartialEq)]
pub struct CvarDesign {
    /// Feasible bisection endpoint under the smoke-tier monotonicity assumption.
    pub scale_star: f64,
    /// Catalog-snapped scale (≥ scale_star).
    pub scale_snapped: f64,
    /// CVaR at the snapped design (re-checked, must pass).
    pub cvar_snapped: f64,
    /// CVaR at the continuous optimum.
    pub cvar_star: f64,
    /// Mass proxy at the snapped design (scale × member count — the
    /// smoke-tier stand-in for Σ ρAL).
    pub mass: f64,
    /// Bisection iterations.
    pub iters: u32,
}

/// Errors that can arise during fallible frame CVaR evaluation and mass minimization.
#[derive(Debug, Clone, PartialEq)]
pub enum FrameCvarError {
    /// The ensemble was empty or did not declare a Kanai-Tajimi spectrum model.
    EnsembleMalformed(&'static str),
    /// Realization of an ensemble member failed.
    RealizationFailed {
        /// Member index that failed realization.
        member: u32,
    },
    /// The checked nonlinear history refused its inputs or failed to converge.
    History {
        /// Ensemble member whose history was refused.
        member: u32,
        /// Original history-owner diagnostic.
        source: HistoryError,
    },
    /// Non-finite drift loss encountered.
    NonFiniteLoss {
        /// Member index where non-finite loss was observed.
        member: u32,
        /// Non-finite loss value observed.
        value: f64,
    },
    /// Canonical CVaR calculation refused.
    Robust(RobustError),
    /// The declared peak-drift CVaR limit was non-finite or negative.
    InvalidLimit {
        /// Offending limit.
        limit: f64,
    },
    /// Even the maximum scale in the bisection range violates the CVaR limit.
    InfeasibleLimit {
        /// Maximum scale in the bisection range.
        hi_scale: f64,
        /// Observed CVaR at the maximum scale.
        cvar_observed: f64,
        /// Declared CVaR limit.
        limit: f64,
    },
    /// The independently evaluated snapped candidate violates the risk target.
    CatalogRiskExceeded {
        /// Catalog scale evaluated.
        scale: f64,
        /// Actual empirical CVaR at that scale.
        cvar_observed: f64,
        /// Declared CVaR limit.
        limit: f64,
    },
    /// The catalog has no section scale greater than or equal to the continuous optimum.
    NoFeasibleCatalogScale {
        /// Continuous optimal scale.
        scale_star: f64,
    },
}

impl core::fmt::Display for FrameCvarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EnsembleMalformed(msg) => write!(f, "malformed ensemble: {msg}"),
            Self::RealizationFailed { member } => write!(f, "failed to realize member {member}"),
            Self::History { member, source } => write!(f, "history member {member}: {source}"),
            Self::NonFiniteLoss { member, value } => {
                write!(f, "non-finite drift loss {value} at member {member}")
            }
            Self::Robust(err) => write!(f, "canonical CVaR error: {err:?}"),
            Self::InvalidLimit { limit } => {
                write!(
                    f,
                    "CVaR drift limit must be finite and nonnegative, got {limit}"
                )
            }
            Self::InfeasibleLimit {
                hi_scale,
                cvar_observed,
                limit,
            } => write!(
                f,
                "even scale {hi_scale} violates CVaR limit {limit} (observed {cvar_observed})"
            ),
            Self::CatalogRiskExceeded { scale, cvar_observed, limit } => write!(
                f,
                "catalog scale {scale} violates CVaR limit {limit} (observed {cvar_observed})"
            ),
            Self::NoFeasibleCatalogScale { scale_star } => {
                write!(f, "catalog has no section above scale {scale_star}")
            }
        }
    }
}

impl std::error::Error for FrameCvarError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::History { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl From<RobustError> for FrameCvarError {
    fn from(err: RobustError) -> Self {
        Self::Robust(err)
    }
}

/// Validate ground motion ensemble fallibly.
pub(crate) fn try_assert_ground_motion_ensemble(
    ensemble: &StochasticEnsemble,
) -> Result<(), FrameCvarError> {
    if ensemble.members == 0 {
        return Err(FrameCvarError::EnsembleMalformed(
            "frame studies require at least one ground-motion ensemble member",
        ));
    }
    if !matches!(ensemble.model, SpectrumModel::KanaiTajimi { .. }) {
        return Err(FrameCvarError::EnsembleMalformed(
            "frame studies require a Kanai-Tajimi ground-acceleration ensemble; \
             wind spectra and material-parameter bands are not structural motions",
        ));
    }
    Ok(())
}

/// Peak-drift losses over the whole ensemble at section scale `s` (fallible).
///
/// Each smoke-tier history uses the existing checked owner: at most 30 Newton
/// corrections per sample, displacement tolerance 1e-12 m and equilibrium
/// tolerance 1 N. The native realization budget already bounds sample count.
/// Input, nonfinite-state and convergence failures retain their member index.
pub fn try_losses(
    ensemble: &StochasticEnsemble,
    base: StoryParams,
    s: f64,
) -> Result<Vec<f64>, FrameCvarError> {
    try_assert_ground_motion_ensemble(ensemble)?;
    let dt = ensemble.dt.value;
    let mut out = Vec::with_capacity(ensemble.members as usize);
    for member in 0..ensemble.members {
        let real = ensemble
            .realize(member)
            .map_err(|_| FrameCvarError::RealizationFailed { member })?;
        let params = StoryParams { scale: s, ..base };
        let mut frame = StoryFrame::new(params);
        let history = frame
            .run_checked(
                GroundMotion::new(&real.values, dt),
                HistoryLimits::new(real.values.len(), 30, 1e-12, 1.0),
            )
            .map_err(|source| FrameCvarError::History { member, source })?;
        let drift = history.peak_abs_displacement_m / base.h;
        if !drift.is_finite() {
            return Err(FrameCvarError::NonFiniteLoss {
                member,
                value: drift,
            });
        }
        out.push(drift);
    }
    Ok(out)
}

/// Fallible calculation of CVaR of the peak-drift loss over the ensemble at section scale `s`.
pub fn try_ensemble_cvar(
    ensemble: &StochasticEnsemble,
    base: StoryParams,
    s: f64,
    beta: f64,
) -> Result<f64, FrameCvarError> {
    let loss_vec = try_losses(ensemble, base, s)?;
    let rep = empirical_cvar(&loss_vec, beta)?;
    Ok(rep.cvar())
}

/// Minimize mass (∝ scale) subject to CVaR_β(peak drift) ≤ `limit` by
/// bisection on the scale, then snap UP to `catalog` and re-check (fallible).
pub fn try_cvar_mass_min(
    ensemble: &StochasticEnsemble,
    base: StoryParams,
    beta: f64,
    limit: f64,
    catalog: &[f64],
) -> Result<CvarDesign, FrameCvarError> {
    if !limit.is_finite() || limit < 0.0 {
        return Err(FrameCvarError::InvalidLimit { limit });
    }
    let (mut lo, mut hi) = (0.25f64, 4.0f64);
    let cvar_hi = try_ensemble_cvar(ensemble, base, hi, beta)?;
    if cvar_hi > limit {
        return Err(FrameCvarError::InfeasibleLimit {
            hi_scale: hi,
            cvar_observed: cvar_hi,
            limit,
        });
    }
    let mut iters = 0u32;
    let cvar_lo = try_ensemble_cvar(ensemble, base, lo, beta)?;
    if cvar_lo <= limit {
        hi = lo;
    }
    while hi - lo > 0.02 {
        let mid = f64::midpoint(lo, hi);
        let cvar_mid = try_ensemble_cvar(ensemble, base, mid, beta)?;
        if cvar_mid <= limit {
            hi = mid;
        } else {
            lo = mid;
        }
        iters += 1;
    }
    let scale_star = hi;
    let cvar_star = try_ensemble_cvar(ensemble, base, scale_star, beta)?;
    let scale_snapped = catalog
        .iter()
        .copied()
        .filter(|&c| c >= scale_star)
        .fold(f64::INFINITY, f64::min);
    if !scale_snapped.is_finite() {
        return Err(FrameCvarError::NoFeasibleCatalogScale { scale_star });
    }
    let cvar_snapped = try_ensemble_cvar(ensemble, base, scale_snapped, beta)?;
    if cvar_snapped > limit {
        return Err(FrameCvarError::CatalogRiskExceeded {
            scale: scale_snapped,
            cvar_observed: cvar_snapped,
            limit,
        });
    }
    Ok(CvarDesign {
        scale_star,
        scale_snapped,
        cvar_snapped,
        cvar_star,
        mass: scale_snapped * 2.0,
        iters,
    })
}

/// CVaR of the peak-drift loss over the ensemble at section scale
/// `s` — the battery's monotonicity probe and limit-bracketing tool.
#[must_use]
pub fn ensemble_cvar(ensemble: &StochasticEnsemble, base: StoryParams, s: f64, beta: f64) -> f64 {
    try_ensemble_cvar(ensemble, base, s, beta)
        .expect("frame-generated losses and beta must satisfy the canonical CVaR contract")
}

/// Minimize mass (∝ scale) subject to CVaR_β(peak drift) ≤ `limit` by
/// bisection on the scale, then snap UP to `catalog` and re-check.
///
/// # Panics
/// On invalid inputs, failed nonlinear histories, or a refused final catalog
/// candidate. The fallible entry point retains the corresponding diagnostic.
#[must_use]
pub fn cvar_mass_min(
    ensemble: &StochasticEnsemble,
    base: StoryParams,
    beta: f64,
    limit: f64,
    catalog: &[f64],
) -> CvarDesign {
    try_cvar_mass_min(ensemble, base, beta, limit, catalog)
        .expect("cvar_mass_min must find a feasible design under the declared limit and catalog")
}
