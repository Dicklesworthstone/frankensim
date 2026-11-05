//! Anytime-valid fragility under a physical ground-acceleration ensemble.
//!
//! The confidence sequence is the estimator of record. A bounded diagnostic
//! summarizes the same consumed peak drifts. An explicit spectrum additionally
//! enables a paired dt/2dt comparison on one admitted physical band and horizon;
//! the legacy record-periodic spectrum is never silently low-pass filtered.

use crate::history::{GroundMotion, HistoryError, HistoryLimits, StoryFrame, StoryParams};
use fs_eproc::GaussianMixtureCs;
use fs_exec::Cx;
use fs_scenario::ScenarioError;
use fs_scenario::ensemble::{
    RealizationBudget, SpectralDiscretization, SpectrumModel, StochasticEnsemble,
};
use fs_uq::mlmc::{MlmcReport, coupled_mlmc_summary};

/// Physical, statistical and per-realization work controls.
#[derive(Debug, Clone, Copy)]
pub struct FragilityControls {
    /// Nonnegative peak drift-ratio threshold.
    pub drift_limit: f64,
    /// Confidence-sequence error probability, strictly between zero and one.
    pub alpha: f64,
    /// Nonnegative stopping radius; zero requests the complete member budget.
    pub margin: f64,
    /// Native scenario sample and spectral-synthesis work limits per member.
    pub realization_budget: RealizationBudget,
    /// Native nonlinear-history limits per model and member.
    pub history_limits: HistoryLimits,
}

/// The fragility record.
#[derive(Debug, Clone)]
pub struct FragilityReport {
    /// Members actually consumed (the e-stop decides).
    pub members_used: u32,
    /// Exceedance estimate (CS center).
    pub p_hat: f64,
    /// CS radius at stop.
    pub radius: f64,
    /// Did the study stop before exhausting the budget?
    pub stopped_early: bool,
    /// Exceedances observed.
    pub exceedances: u32,
    /// The exact explicit source grid, or the original record-periodic source.
    pub spectrum: Option<SpectralDiscretization>,
    /// Last physical observation time, shared by every diagnostic level [s].
    pub end_time_s: f64,
    /// Observed-member peak-drift diagnostics, without adaptive resampling.
    ///
    /// With an explicit spectrum, level 0 is the 2dt history and level 1 is
    /// fine minus coarse on the same member. A shorter final coarse interval
    /// retains the original endpoint. Otherwise the report has only the dt
    /// level: the legacy spectrum reaches its Nyquist limit and cannot be
    /// decimated honestly. Costs count actual advancing Newmark steps.
    ///
    /// The estimate and variance use the consumed fine values, accounting for
    /// covariance between paired levels. They are descriptive sample moments,
    /// not an anytime confidence interval or a variance-targeted MLMC result.
    pub mlmc: MlmcReport,
}

/// A source, nonlinear solve, resource, or cancellation refusal.
#[derive(Debug)]
pub enum FragilityError {
    /// Malformed ensemble or probability controls.
    InvalidInput(&'static str),
    /// The declared physical band does not fit the coarse observation clock.
    CoarseBandAliased {
        /// Requested physical cutoff [rad/s].
        omega_max_rad_s: f64,
        /// Nominal coarse sample interval [s].
        coarse_dt_s: f64,
    },
    /// Original scenario-owner refusal.
    Realization {
        /// Original ensemble member.
        member: u32,
        /// Original diagnostic.
        source: ScenarioError,
    },
    /// Original nonlinear-history refusal.
    History {
        /// Original ensemble member.
        member: u32,
        /// True for the coarse history; false for the fine history.
        coarse: bool,
        /// Original diagnostic.
        source: HistoryError,
    },
    /// Peak-drift diagnostic storage could not be reserved.
    AllocationFailed,
    /// Cancellation was observed at a realization/history boundary.
    Cancelled,
}

impl core::fmt::Display for FragilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(what) => write!(f, "invalid fragility input: {what}"),
            Self::CoarseBandAliased {
                omega_max_rad_s,
                coarse_dt_s,
            } => write!(
                f,
                "physical cutoff {omega_max_rad_s} rad/s exceeds the Nyquist limit at {coarse_dt_s} s"
            ),
            Self::Realization { member, source } => {
                write!(f, "fragility member {member}: {source}")
            }
            Self::History {
                member,
                coarse,
                source,
            } => {
                write!(f, "fragility member {member}, coarse={coarse}: {source}")
            }
            Self::AllocationFailed => write!(f, "fragility diagnostic allocation refused"),
            Self::Cancelled => write!(f, "fragility study cancelled"),
        }
    }
}

impl std::error::Error for FragilityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Realization { source, .. } => Some(source),
            Self::History { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Run the e-stopped fragility study with the original ensemble spectrum.
///
/// Source row zero establishes initial force equilibrium; N rows advance
/// N-1 steps. The diagnostic has one level because the default source reaches
/// the fine clock's Nyquist frequency. Use
/// [`e_stopped_fragility_with_spectrum`] for a physical paired dt/2dt study.
///
/// # Panics
/// On invalid study controls, a refused realization, or a failed checked
/// nonlinear history. Native realization limits and 30 Newton corrections
/// per step bound each member; the member count bounds the complete study.
#[must_use]
pub fn e_stopped_fragility(
    ensemble: &StochasticEnsemble,
    params: StoryParams,
    drift_limit: f64,
    alpha: f64,
    margin: f64,
) -> FragilityReport {
    let realization_budget = RealizationBudget::default();
    study(
        ensemble,
        params,
        None,
        FragilityControls {
            drift_limit,
            alpha,
            margin,
            realization_budget,
            history_limits: HistoryLimits::new(realization_budget.max_samples, 30, 1e-12, 1.0),
        },
        || Ok(()),
    )
    .expect("fragility source, controls, and nonlinear histories must be admissible")
}

/// Run the original-spectrum study with explicit limits and cancellation.
///
/// This has the same initial-time and single-level diagnostic semantics as
/// [`e_stopped_fragility`], returning the original source/history refusals.
///
/// # Errors
/// Invalid controls, original source/history refusal, allocation failure or
/// observed cancellation. No partial report is returned.
pub fn try_e_stopped_fragility(
    ensemble: &StochasticEnsemble,
    params: StoryParams,
    controls: FragilityControls,
    cx: &Cx<'_>,
) -> Result<FragilityReport, FragilityError> {
    study(ensemble, params, None, controls, || {
        cx.checkpoint().map_err(|_| FragilityError::Cancelled)
    })
}

/// Run a checked fragility study on an explicit physical frequency grid.
///
/// The original scenario owner realizes each consumed member exactly once.
/// The same values drive dt and 2dt histories, with a final dt interval on the
/// coarse level when needed. The cutoff must fit the coarse Nyquist limit:
/// `omega_max_rad_s * 2dt <= pi`. No high-frequency tail is folded or
/// redistributed. Both histories start from the same initial force balance.
///
/// The CS is valid at its stopping time for the admitted finite-band nonlinear
/// model. Bandwidth, frequency quadrature and structural time discretization
/// remain separate approximation choices; the diagnostic is not a convergence
/// certificate or population qualification. Cancellation checkpoints surround
/// each bounded realization/history and publication.
///
/// # Errors
/// Invalid controls/band, original source or history refusal, allocation
/// failure, or observed cancellation. No partial report is returned.
pub fn e_stopped_fragility_with_spectrum(
    ensemble: &StochasticEnsemble,
    params: StoryParams,
    spectrum: SpectralDiscretization,
    controls: FragilityControls,
    cx: &Cx<'_>,
) -> Result<FragilityReport, FragilityError> {
    study(ensemble, params, Some(spectrum), controls, || {
        cx.checkpoint().map_err(|_| FragilityError::Cancelled)
    })
}

#[allow(clippy::too_many_lines)] // One member transaction, then bounded paired diagnostics.
fn study(
    ensemble: &StochasticEnsemble,
    params: StoryParams,
    spectrum: Option<SpectralDiscretization>,
    controls: FragilityControls,
    mut checkpoint: impl FnMut() -> Result<(), FragilityError>,
) -> Result<FragilityReport, FragilityError> {
    checkpoint()?;
    if ensemble.members == 0 || !matches!(ensemble.model, SpectrumModel::KanaiTajimi { .. }) {
        return Err(FragilityError::InvalidInput(
            "a nonempty Kanai-Tajimi ground-acceleration ensemble is required",
        ));
    }
    if !controls.drift_limit.is_finite()
        || controls.drift_limit < 0.0
        || !controls.alpha.is_finite()
        || !(0.0..1.0).contains(&controls.alpha)
        || controls.alpha == 0.0
        || !controls.margin.is_finite()
        || controls.margin < 0.0
    {
        return Err(FragilityError::InvalidInput(
            "finite nonnegative drift limit/margin and alpha strictly in (0,1) are required",
        ));
    }
    let dt = ensemble.dt.value;
    if let Some(grid) = spectrum
        && !(grid.omega_max_rad_s <= core::f64::consts::PI / (2.0 * dt))
    {
        return Err(FragilityError::CoarseBandAliased {
            omega_max_rad_s: grid.omega_max_rad_s,
            coarse_dt_s: 2.0 * dt,
        });
    }
    let mut cs = GaussianMixtureCs::new(0.5, 8.0, controls.alpha);
    let (mut exceedances, mut used, mut stopped_early) = (0u32, 0u32, false);
    let mut peaks = Vec::new();
    let mut end_time_s = 0.0;
    let mut advancing_steps = 0usize;
    for member in 0..ensemble.members {
        checkpoint()?;
        let real = match spectrum {
            Some(grid) => ensemble.realize_with_spectrum(member, grid, controls.realization_budget),
            None => ensemble.realize_with_budget(member, controls.realization_budget),
        }
        .map_err(|source| FragilityError::Realization { member, source })?;
        checkpoint()?;
        let motion = GroundMotion::new(&real.values, dt);
        let mut frame = StoryFrame::new(params);
        let fine = frame
            .run_from_time_zero(motion, controls.history_limits)
            .map_err(|source| FragilityError::History {
                member,
                coarse: false,
                source,
            })?;
        let peak = fine.peak_abs_displacement_m / params.h;
        checkpoint()?;
        let coarse = if spectrum.is_some() {
            StoryFrame::new(params)
                .run_sampled_checked(motion, controls.history_limits, 2)
                .map_err(|source| FragilityError::History {
                    member,
                    coarse: true,
                    source,
                })?
                .peak_abs_displacement_m
                / params.h
        } else {
            peak
        };
        checkpoint()?;
        if !peak.is_finite() || !coarse.is_finite() {
            return Err(FragilityError::InvalidInput(
                "peak drift ratio is not finite",
            ));
        }
        peaks
            .try_reserve(1)
            .map_err(|_| FragilityError::AllocationFailed)?;
        peaks.push([coarse, peak]);
        advancing_steps = real.values.len() - 1;
        end_time_s = *real
            .times
            .last()
            .expect("admitted spectral history has times");
        let x = if peak > controls.drift_limit {
            1.0
        } else {
            0.0
        };
        if x > 0.5 {
            exceedances += 1;
        }
        cs.observe(x);
        used = member + 1;
        if let Some((_, radius)) = cs.interval()
            && used >= 8
            && radius <= controls.margin
        {
            stopped_early = used < ensemble.members;
            break;
        }
    }
    let (p_hat, radius) = cs.interval().expect("observed at least one member");
    checkpoint()?;
    let mlmc = diagnostic(&peaks, advancing_steps, spectrum.is_some())?;
    checkpoint()?;
    Ok(FragilityReport {
        members_used: used,
        p_hat,
        radius,
        stopped_early,
        exceedances,
        spectrum,
        end_time_s,
        mlmc,
    })
}

fn diagnostic(
    peaks: &[[f64; 2]],
    advancing_steps: usize,
    paired: bool,
) -> Result<MlmcReport, FragilityError> {
    // Actual advancing-step work, with coarse histories reused by the paired
    // correction. The native owner includes cross-level covariance through the
    // raw fine moments and never requests another physical realization.
    #[allow(clippy::cast_precision_loss)]
    let fine_cost = advancing_steps as f64;
    let report = if paired {
        #[allow(clippy::cast_precision_loss)]
        let coarse_cost = advancing_steps.div_ceil(2) as f64;
        coupled_mlmc_summary(peaks.iter().copied(), [coarse_cost, fine_cost])
    } else {
        coupled_mlmc_summary(peaks.iter().map(|peak| [peak[1]]), [fine_cost])
    };
    report.map_err(|_| FragilityError::InvalidInput("diagnostic moments are not finite"))
}
