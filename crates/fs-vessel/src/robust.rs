//! Stage 4: ROBUSTIFICATION. The stability objective is cheap (an
//! Orr–Sommerfeld solve per station), so the fluid band — Carreau
//! parameters spanning the target liquids × pour rates — is swept
//! directly: the NOMINAL design minimizes growth at the band center,
//! the ROBUST design minimizes the CVaR of growth over the band
//! (Rockafellar–Uryasev empirical form from fs-robust). The
//! battery gates the flagship claim: the robust lip beats the nominal
//! lip on off-nominal fluids. Candidate screening under the (noisy)
//! validator runs through fs-race — dominated lips die early with
//! anytime validity, eliminations ledgered.

use crate::stability::{VesselProfile, growth_objective};
pub use fs_robust::{EmpiricalCvarReport, RobustError, cvar, empirical_cvar};

/// The robustification record.
#[derive(Debug, Clone)]
pub struct RobustReport {
    /// Lip width chosen by nominal-only optimization.
    pub nominal_lip: f64,
    /// Lip width chosen by CVaR optimization over the band.
    pub robust_lip: f64,
    /// Worst off-nominal growth of the nominal design.
    pub nominal_offband_growth: f64,
    /// Worst off-nominal growth of the robust design.
    pub robust_offband_growth: f64,
    /// CVaR level used.
    pub beta: f64,
}

/// The (rate, viscosity) band: nominal center plus off-nominal
/// corners (the family of liquids the spout must serve).
#[must_use]
pub fn fluid_band() -> Vec<(f64, f64)> {
    let rates = [0.8f64, 1.0, 1.3];
    let viscs = [0.6f64, 1.0, 1.8];
    let mut band = Vec::new();
    for &r in &rates {
        for &v in &viscs {
            band.push((r, v));
        }
    }
    band
}

/// Optimize the lip width two ways (nominal vs CVaR over the band) by
/// a deterministic golden-section-free grid refinement (the objective
/// is cheap and smooth in the scalar knob), then evaluate both on the
/// OFF-NOMINAL corners.
#[must_use]
pub fn robustify(beta: f64) -> RobustReport {
    let stations = 4;
    let modes = 4;
    let nominal_obj = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        growth_objective(&p, 1.0, 1.0, stations, modes)
    };
    let cvar_obj = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        let losses: Vec<f64> = fluid_band()
            .iter()
            .map(|&(r, v)| growth_objective(&p, r, v, stations, modes))
            .collect();
        cvar(&losses, beta)
            .expect("fluid-band losses and beta must satisfy the canonical CVaR contract")
    };
    let minimize = |f: &dyn Fn(f64) -> f64| -> f64 {
        // Two-stage deterministic grid refinement over the lip range.
        let (mut lo, mut hi) = (0.5f64, 3.0f64);
        let mut best = lo;
        for _ in 0..2 {
            let mut best_v = f64::INFINITY;
            for k in 0..=12 {
                let lip = lo + (hi - lo) * f64::from(k) / 12.0;
                let v = f(lip);
                if v < best_v {
                    best_v = v;
                    best = lip;
                }
            }
            let third = (hi - lo) / 3.0;
            lo = (best - third).max(0.5);
            hi = (best + third).min(3.0);
        }
        best
    };
    let nominal_lip = minimize(&nominal_obj);
    let robust_lip = minimize(&cvar_obj);
    // Off-nominal evaluation: worst growth over the band corners
    // EXCLUDING the nominal center.
    let offband = |lip: f64| -> f64 {
        let p = VesselProfile::carafe(lip);
        fluid_band()
            .iter()
            .filter(|&&(r, v)| (r - 1.0).abs() > 1e-12 || (v - 1.0).abs() > 1e-12)
            .map(|&(r, v)| growth_objective(&p, r, v, stations, modes))
            .fold(f64::NEG_INFINITY, f64::max)
    };
    RobustReport {
        nominal_lip,
        robust_lip,
        nominal_offband_growth: offband(nominal_lip),
        robust_offband_growth: offband(robust_lip),
        beta,
    }
}

/// Physical film-width selection over a finite authored candidate set.
///
/// Widths are dimensionless multipliers of the study's SI section widths.
/// Both nominal and robust choices may coincide or lie on the candidate
/// boundary. No interior optimum or strictly different robust choice is forced.
#[derive(Debug, Clone, PartialEq)]
pub struct LaminarRobustReport {
    /// Every admitted candidate in authored order; consumers can plot the
    /// actual evaluated objectives without solving the entire band again.
    pub candidates: Vec<LaminarCandidateReport>,
    /// Candidate minimizing nominal growth.
    pub nominal_candidate: usize,
    /// Candidate minimizing the fluid band's CVaR.
    pub robust_candidate: usize,
    /// Nominal candidate's width multiplier.
    pub nominal_width_scale: f64,
    /// Robust candidate's width multiplier.
    pub robust_width_scale: f64,
    /// Fluid-band CVaR at the nominal choice, in 1/s.
    pub nominal_band_cvar_per_s: f64,
    /// Fluid-band CVaR at the robust choice, in 1/s.
    pub robust_band_cvar_per_s: f64,
    /// Largest sampled band growth at the nominal choice, in 1/s.
    pub nominal_worst_growth_per_s: f64,
    /// Largest sampled band growth at the robust choice, in 1/s.
    pub robust_worst_growth_per_s: f64,
    /// CVaR level.
    pub beta: f64,
}

/// Physical objective and risk values for one admitted width candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarCandidateReport {
    /// Dimensionless multiplier of every declared section width.
    pub width_scale: f64,
    /// Nominal worst growth, in 1/s.
    pub nominal_growth_per_s: f64,
    /// Canonical empirical band CVaR, in 1/s.
    pub band_cvar_per_s: f64,
    /// Worst growth among the supplied band operating points, in 1/s.
    pub band_worst_growth_per_s: f64,
}

/// A refused physical candidate or risk calculation; no optimum is published.
#[derive(Debug, Clone, PartialEq)]
pub enum LaminarRobustError {
    /// At least one authored width candidate is required.
    EmptyCandidates,
    /// At least one physical band operating point is required.
    EmptyBand,
    /// Native film admission failed.
    Film {
        /// Candidate index.
        candidate: usize,
        /// None for the nominal point; otherwise the band index.
        fluid: Option<usize>,
        /// Original physical refusal.
        source: crate::stability::laminar::LaminarFilmError,
    },
    /// Canonical empirical CVaR refused.
    Risk(RobustError),
}

impl core::fmt::Display for LaminarRobustError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyCandidates => {
                write!(f, "physical film optimization needs a width candidate")
            }
            Self::EmptyBand => write!(f, "physical film optimization needs a fluid-band point"),
            Self::Film {
                candidate,
                fluid,
                source,
            } => {
                write!(f, "film candidate {candidate}, fluid {fluid:?}: {source}")
            }
            Self::Risk(source) => write!(f, "film-band CVaR: {source:?}"),
        }
    }
}

impl std::error::Error for LaminarRobustError {}

/// Optimize the actual free-surface long-wave screen, with explicit SI inputs.
///
/// Every candidate must be admitted at the nominal point AND every band point.
/// A rejected model is returned with its candidate/fluid indices; it is never
/// dropped from the band or assigned an artificial penalty. Candidate order
/// deterministically breaks exact ties. CVaR minimization guarantees no worse
/// band CVaR than the nominal choice, but need not minimize the band's maximum.
pub fn robustify_laminar<S: crate::stability::laminar::LaminarGrowthObjective + ?Sized>(
    study: &S,
    width_scales: &[f64],
    nominal: crate::stability::laminar::LaminarPourConditions,
    band: &[crate::stability::laminar::LaminarPourConditions],
    beta: f64,
) -> Result<LaminarRobustReport, LaminarRobustError> {
    if width_scales.is_empty() {
        return Err(LaminarRobustError::EmptyCandidates);
    }
    if band.is_empty() {
        return Err(LaminarRobustError::EmptyBand);
    }
    // Validate the requested risk level before evaluating any candidate.
    let _ = cvar(&[0.0], beta).map_err(LaminarRobustError::Risk)?;
    let mut nominal_candidate = 0;
    let mut robust_candidate = 0;
    let mut best_nominal = f64::INFINITY;
    let mut best_cvar = f64::INFINITY;
    let mut candidate_risks = Vec::with_capacity(width_scales.len());
    let mut losses = Vec::with_capacity(band.len());
    for (candidate, &width) in width_scales.iter().enumerate() {
        let nominal_growth =
            study
                .growth_per_s(width, nominal)
                .map_err(|source| LaminarRobustError::Film {
                    candidate,
                    fluid: None,
                    source,
                })?;
        if !nominal_growth.is_finite() {
            return Err(LaminarRobustError::Risk(RobustError::BadSample {
                value: nominal_growth,
            }));
        }
        if nominal_growth < best_nominal {
            best_nominal = nominal_growth;
            nominal_candidate = candidate;
        }
        losses.clear();
        let mut worst = f64::NEG_INFINITY;
        for (fluid, &conditions) in band.iter().enumerate() {
            let growth = study.growth_per_s(width, conditions).map_err(|source| {
                LaminarRobustError::Film {
                    candidate,
                    fluid: Some(fluid),
                    source,
                }
            })?;
            losses.push(growth);
            worst = worst.max(growth);
        }
        let risk = cvar(&losses, beta).map_err(LaminarRobustError::Risk)?;
        candidate_risks.push(LaminarCandidateReport {
            width_scale: width,
            nominal_growth_per_s: nominal_growth,
            band_cvar_per_s: risk,
            band_worst_growth_per_s: worst,
        });
        if risk < best_cvar {
            best_cvar = risk;
            robust_candidate = candidate;
        }
    }
    Ok(LaminarRobustReport {
        nominal_candidate,
        robust_candidate,
        nominal_width_scale: width_scales[nominal_candidate],
        robust_width_scale: width_scales[robust_candidate],
        nominal_band_cvar_per_s: candidate_risks[nominal_candidate].band_cvar_per_s,
        robust_band_cvar_per_s: candidate_risks[robust_candidate].band_cvar_per_s,
        nominal_worst_growth_per_s: candidate_risks[nominal_candidate].band_worst_growth_per_s,
        robust_worst_growth_per_s: candidate_risks[robust_candidate].band_worst_growth_per_s,
        beta,
        candidates: candidate_risks,
    })
}
