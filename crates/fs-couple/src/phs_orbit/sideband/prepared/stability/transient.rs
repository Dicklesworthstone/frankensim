//! Complete stroboscopic growth in the original disclosed physical state chart.

use super::{C64, FloquetClassification, PhsBaseStabilityReport};
use fs_math::dd::Dd;
use fs_spectral::floquet::transient::{
    StroboscopicGrowthControls, StroboscopicGrowthError, StroboscopicGrowthReport,
    stroboscopic_growth,
};

/// Independent gain convergence and TOTAL three-map resource admission.
#[derive(Clone, Debug)]
pub struct PhsBaseTransientControls {
    /// Shared per-map service limits, including every requested horizon and
    /// its complete optimal initial/final scaled perturbation directions.
    pub growth: StroboscopicGrowthControls,
    /// Absolute natural-log amplitude-gain allowance between adjacent grids.
    /// Both adjacent comparisons must pass at EVERY requested horizon.
    pub log_gain_change_tolerance: f64,
    /// TOTAL new numerical payload for all three retained growth reports,
    /// shared workspaces, one complex map conversion, scales and comparisons.
    /// The borrowed original stability report and caller horizons are excluded.
    /// Their previously admitted storage adds to this operation's live bound.
    /// Allocator bookkeeping is excluded. The conservative bound reserves
    /// THREE shared payload allowances even though their workspaces run serially.
    pub maximum_workspace_bytes: usize,
    /// TOTAL arithmetic admission: three shared work allowances plus all new
    /// conversions/comparisons. This is a conservative cap, not measured FLOPs;
    /// no allowance resets while processing the three temporal resolutions.
    pub maximum_work: u64,
}

/// Every requested horizon on one retained, complete original monodromy.
#[derive(Clone, Debug)]
pub struct PhsBaseTransientLevel {
    pub steps_per_period: usize,
    /// Shared full-dimension samples and their original power/SVD residuals.
    /// Directions are in the disclosed scaled state chart, including memory
    /// and autonomous phase. No direction or neutral coordinate is discarded.
    pub growth: StroboscopicGrowthReport,
    /// Absolute adjacent log-gain changes divided by the independent allowance,
    /// in the same order as growth.samples. None for the first temporal grid.
    /// Two exactly zero gains compare as zero; zero versus nonzero is infinity.
    pub log_gain_change_ratios: Option<Vec<f64>>,
}

/// Finite-horizon amplification separate from asymptotic Floquet stability.
#[derive(Clone, Debug)]
pub struct PhsBaseTransientReport {
    pub omega_rad_s: f64,
    pub period_s: f64,
    pub autonomous: bool,
    /// The squared norm is EXACTLY sum_i |delta x_i / state_scales_i|^2 at the original
    /// phase-zero section. Amplitude gain is the square root of its ratio.
    /// This authored chart is not automatically physical energy, and changing
    /// relative coordinate scales changes the question being answered.
    pub state_scales: Vec<f64>,
    /// Previously measured parent classification, independent of finite gain.
    /// Stable maps may strongly amplify perturbations before asymptotic decay.
    pub asymptotic_classification: FloquetClassification,
    /// Both previously measured complete-map temporal comparisons passed.
    pub full_maps_resolved: bool,
    /// Requires full_maps_resolved AND both adjacent gain comparisons at every
    /// sampled horizon. Singular directions are not compared: their optimum
    /// may be nonunique even when the amplitude gain is accurately resolved.
    pub temporally_resolved: bool,
    /// All three complete shared reports; sample.periods counts full ORIGINAL
    /// physical periods from the same section, including zero when requested.
    pub levels: Vec<PhsBaseTransientLevel>,
    pub admitted_workspace_bytes: usize,
    pub admitted_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PhsBaseTransientError {
    Invalid(&'static str),
    Budget(&'static str),
    Growth(StroboscopicGrowthError),
    Cancelled,
}
impl core::fmt::Display for PhsBaseTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical stroboscopic growth refused: {self:?}")
    }
}
impl std::error::Error for PhsBaseTransientError {}
impl From<StroboscopicGrowthError> for PhsBaseTransientError {
    fn from(error: StroboscopicGrowthError) -> Self {
        match error {
            StroboscopicGrowthError::Cancelled => Self::Cancelled,
            other => Self::Growth(other),
        }
    }
}
type Result<T> = core::result::Result<T, PhsBaseTransientError>;

impl PhsBaseStabilityReport {
    /// Measure optimal complete-state growth using the THREE maps already here.
    ///
    /// The shared service computes ||(S^-1*M*S)^k||_2 and complete optimal
    /// initial/final unit directions for each explicitly requested period count.
    /// Repeating this query calls no physical source, integrator, QR or
    /// eigensolver. Original maps, all original modes, phase evidence and the
    /// parent stability report remain unchanged and available for reuse.
    ///
    /// The measured squared norm is sum_i |delta x_i/state_scales_i|^2 on the ORIGINAL
    /// phase-zero section. Physical perturbations are obtained by multiplying
    /// each direction component by its disclosed state scale. An ordinary gain
    /// may be unrepresentable while its natural log remains available. A zero
    /// power has zero gain, no finite log and no final unit direction; comparing
    /// zero to any nonzero power never passes temporal gain convergence.
    ///
    /// Both adjacent original complete-map checks and all adjacent requested
    /// log-gain checks must pass for temporal resolution. This criterion is
    /// independent of asymptotic stability. Phase and memory remain in the
    /// optimization even for an autonomous cycle. There is no transverse
    /// quotient, within-period maximum, unrequested-horizon maximum, general
    /// physical-energy interpretation or global nonlinear basin conclusion.
    /// These finite-map residuals and comparisons are numerical evidence, not
    /// forward-error or continuous-evolution enclosures.
    ///
    /// # Errors
    /// Invalid complete retained maps/scales, unsorted or repeated horizons,
    /// exhausted aggregate/shared budgets, shared numerical refusal, or
    /// cancellation. No partially evaluated report is returned.
    pub fn analyze_stroboscopic_growth(
        &self,
        periods: &[u64],
        controls: &PhsBaseTransientControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsBaseTransientReport> {
        poll(&mut cancelled)?;
        let (dimension, bytes, work) = admit(self, periods, controls)?;
        let full_maps_resolved = self.levels.iter().skip(1).all(|level| {
            level
                .monodromy_change_ratio
                .is_some_and(|value| value.is_finite() && (0.0..=1.0).contains(&value))
        });
        let mut levels: Vec<PhsBaseTransientLevel> = Vec::with_capacity(3);
        let mut resolved = full_maps_resolved;
        for original in &self.levels {
            poll(&mut cancelled)?;
            let matrix: Vec<_> = original
                .monodromy
                .iter()
                .copied()
                .map(C64::from_re)
                .collect();
            let growth =
                stroboscopic_growth(&matrix, dimension, periods, controls.growth, &mut cancelled)?;
            let changes: Option<Vec<_>> = levels.last().map(|previous| {
                previous
                    .growth
                    .samples
                    .iter()
                    .zip(&growth.samples)
                    .map(|(old, new)| {
                        log_change(
                            old.log_gain,
                            new.log_gain,
                            controls.log_gain_change_tolerance,
                        )
                    })
                    .collect()
            });
            if let Some(changes) = &changes {
                resolved &= changes.iter().all(|&value| value <= 1.0);
            }
            levels.push(PhsBaseTransientLevel {
                steps_per_period: original.steps_per_period,
                growth,
                log_gain_change_ratios: changes,
            });
        }
        poll(&mut cancelled)?;
        Ok(PhsBaseTransientReport {
            omega_rad_s: self.omega_rad_s,
            period_s: self.period_s,
            autonomous: self.autonomous,
            state_scales: self.state_scales.clone(),
            asymptotic_classification: self.classification,
            full_maps_resolved,
            temporally_resolved: resolved,
            levels,
            admitted_workspace_bytes: bytes,
            admitted_work: work,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PhsBaseTransientError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit(
    report: &PhsBaseStabilityReport,
    periods: &[u64],
    controls: &PhsBaseTransientControls,
) -> Result<(usize, usize, u64)> {
    let n = report.state_scales.len();
    let g = controls.growth;
    if n == 0
        || n > 1024
        || n > g.maximum_dimension
        || report.levels.len() != 3
        || !report.period_s.is_finite()
        || report.period_s <= 0.0
        || !report.omega_rad_s.is_finite()
        || report.omega_rad_s <= 0.0
        || report
            .state_scales
            .iter()
            .any(|&value| !value.is_finite() || value <= 0.0)
        || !controls.log_gain_change_tolerance.is_finite()
        || controls.log_gain_change_tolerance <= 0.0
        || periods.is_empty()
        || periods.len() > 4096
        || periods.len() > g.maximum_horizons
        || periods.windows(2).any(|pair| pair[0] >= pair[1])
        || periods.iter().any(|&count| count > g.maximum_periods)
        || g.svd_max_sweeps == 0
        || !g.svd_tolerance.is_finite()
        || g.svd_tolerance <= f64::EPSILON
        || g.svd_tolerance >= 1.0
        || !g.maximum_relative_residual.is_finite()
        || g.maximum_relative_residual <= 0.0
        || g.maximum_relative_residual >= 1.0
    {
        return Err(PhsBaseTransientError::Invalid(
            "complete finite original maps/scales, sorted horizons and positive growth controls required",
        ));
    }
    let first = report.levels[0].steps_per_period;
    if first == 0
        || first.checked_mul(2) != Some(report.levels[1].steps_per_period)
        || first.checked_mul(4) != Some(report.levels[2].steps_per_period)
        || report.levels.iter().any(|level| {
            level.monodromy.len() != n * n || level.monodromy.iter().any(|value| !value.is_finite())
        })
    {
        return Err(PhsBaseTransientError::Invalid(
            "three complete finite original maps on N,2N,4N grids required",
        ));
    }
    let (d, h) = (n as u128, periods.len() as u128);
    // THREE shared full-call allowances conservatively retain prior reports
    // while the next call owns its workspace. Borrowed parent maps are excluded.
    let bytes = 3 * g.maximum_payload_bytes as u128 + 32 * d * d + 256 * d + 128 * h + 8192;
    let work = 3 * g.maximum_work as u128 + 256 * d * d + 128 * d + 128 * h;
    if bytes > controls.maximum_workspace_bytes as u128 || work > controls.maximum_work as u128 {
        return Err(PhsBaseTransientError::Budget(
            "three growth reports and comparisons exceed the TOTAL operation allowance",
        ));
    }
    Ok((n, bytes as usize, work as u64))
}

fn log_change(old: Option<f64>, new: Option<f64>, tolerance: f64) -> f64 {
    match (old, new) {
        (None, None) => 0.0,
        (Some(a), Some(b)) => {
            let change = ((Dd::from_f64(a) - Dd::from_f64(b)).to_f64() / tolerance).abs();
            if change.is_finite() {
                change
            } else {
                f64::INFINITY
            }
        }
        _ => f64::INFINITY,
    }
}
