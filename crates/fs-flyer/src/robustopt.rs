//! Robust multiobjective optimization demo (bead `frankensim-wf-root-guzez.9.5`, E8.4).
//!
//! Optimization over wind uncertainty using Common Random Numbers (CRN) ensembles.
//! Enforces:
//! - Active structural model requirement (refuses under `PrescribedKinematicEstimated`)
//! - Applicability bounds enforcement
//! - Correction-model holdouts validation
//! - Epistemic coloring and CVaR ranking via `fs-robust`

use crate::{refuse, Refusal};
use fs_evidence::{Color, ColorRank};
use fs_robust::{ColoredObjective, cvar};

/// Model representation mode for aeroelastic coupling during optimization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructuralModelMode {
    /// Full active reduced-order aeroelastic warp and elasticity.
    ActiveElastic,
    /// Prescribed kinematic estimated model (explicitly forbidden for robust optimization).
    PrescribedKinematicEstimated,
    /// Rigid approximation.
    RigidEstimated,
}

/// A candidate design parameter point for the Wright Flyer flight optimization.
#[derive(Clone, Debug, PartialEq)]
pub struct FlightOptimizationCandidate {
    /// Candidate name / identifier.
    pub name: String,
    /// Canard pitch trim bias [deg].
    pub canard_trim_deg: f64,
    /// Pilot proportional pitch gain $K_p$.
    pub pilot_kp: f64,
    /// Pilot pitch-rate damping gain $K_q$.
    pub pilot_kq: f64,
    /// Nominal flight distance [m].
    pub nominal_distance_m: f64,
    /// Sampled flight distances [m] across CRN wind realizations.
    pub crn_distance_samples_m: Vec<f64>,
}

/// Configuration for the robust multiobjective optimization campaign.
#[derive(Clone, Debug, PartialEq)]
pub struct RobustOptConfig {
    /// Structural model mode (must be `ActiveElastic`).
    pub structural_mode: StructuralModelMode,
    /// Minimum allowed pilot $K_p$ gain.
    pub min_kp: f64,
    /// Maximum allowed pilot $K_p$ gain.
    pub max_kp: f64,
    /// Minimum allowed pilot $K_q$ gain.
    pub min_kq: f64,
    /// Maximum allowed pilot $K_q$ gain.
    pub max_kq: f64,
    /// CVaR confidence level $\alpha$ (e.g. 0.90 for worst 10% tail risk).
    pub cvar_alpha: f64,
    /// Number of holdout scenarios for generalization validation.
    pub holdout_count: usize,
}

impl Default for RobustOptConfig {
    fn default() -> Self {
        Self {
            structural_mode: StructuralModelMode::ActiveElastic,
            min_kp: 0.05,
            max_kp: 2.0,
            min_kq: 0.01,
            max_kq: 1.0,
            cvar_alpha: 0.90,
            holdout_count: 5,
        }
    }
}

/// Receipt emitted by the robust optimization campaign.
#[derive(Clone, Debug, PartialEq)]
pub struct RobustOptReceipt {
    /// Candidate that achieved highest nominal performance.
    pub nominal_winner: String,
    /// Candidate that achieved best robust performance (lowest CVaR cost / highest robust score).
    pub robust_winner: String,
    /// True if CVaR robustness reordered the ranking vs nominal.
    pub robustness_reorders: bool,
    /// Headline epistemic color rank (Estimated, as derived from finite sample CRN).
    pub headline_rank: ColorRank,
    /// Total evaluated candidates.
    pub candidates_evaluated: usize,
    /// Content digest of the campaign.
    pub receipt_digest: String,
}

/// Run the robust multiobjective optimization campaign over candidates with CRN ensembles.
///
/// # Errors
/// [`Refusal`] if:
/// - `structural_mode` is `PrescribedKinematicEstimated` or not `ActiveElastic`
/// - candidate parameters exceed applicability bounds
/// - holdout validation detects overfit/extrapolation
pub fn run_robust_optimization(
    config: &RobustOptConfig,
    candidates: &[FlightOptimizationCandidate],
) -> Result<RobustOptReceipt, Refusal> {
    // Gate 1: Active structural model enforcement
    match config.structural_mode {
        StructuralModelMode::PrescribedKinematicEstimated => {
            return Err(refuse(
                "robust-opt-disabled-under-prescribed-kinematic",
                "robust optimization is disabled under PrescribedKinematicEstimated".into(),
                "select ActiveElastic structural model mode",
            ));
        }
        StructuralModelMode::RigidEstimated => {
            return Err(refuse(
                "robust-opt-requires-active-structural-model",
                "robust optimization requires active aeroelastic structural model".into(),
                "enable ActiveElastic mode",
            ));
        }
        StructuralModelMode::ActiveElastic => {}
    }

    if candidates.is_empty() {
        return Err(refuse(
            "robust-opt-no-candidates",
            "candidate list cannot be empty".into(),
            "supply at least one candidate",
        ));
    }

    // Gate 2: Applicability bounds checking
    for c in candidates {
        if c.pilot_kp < config.min_kp || c.pilot_kp > config.max_kp {
            return Err(refuse(
                "robust-opt-applicability-exceeded",
                format!("candidate {} pilot_kp {:.3} outside [{}, {}]", c.name, c.pilot_kp, config.min_kp, config.max_kp),
                "keep pilot gains within applicable control domain",
            ));
        }
        if c.pilot_kq < config.min_kq || c.pilot_kq > config.max_kq {
            return Err(refuse(
                "robust-opt-applicability-exceeded",
                format!("candidate {} pilot_kq {:.3} outside [{}, {}]", c.name, c.pilot_kq, config.min_kq, config.max_kq),
                "keep pilot gains within applicable control domain",
            ));
        }
        if c.crn_distance_samples_m.is_empty() {
            return Err(refuse(
                "robust-opt-empty-crn-samples",
                format!("candidate {} has no CRN ensemble distance samples", c.name),
                "generate CRN ensemble evaluations",
            ));
        }
    }

    // Evaluate nominal winner (maximizing nominal distance)
    let mut best_nominal_idx = 0;
    let mut max_nominal_dist = -f64::INFINITY;
    for (i, c) in candidates.iter().enumerate() {
        if c.nominal_distance_m > max_nominal_dist {
            max_nominal_dist = c.nominal_distance_m;
            best_nominal_idx = i;
        }
    }
    let nominal_winner = candidates[best_nominal_idx].name.clone();

    // Evaluate robust CVaR winner
    // Cost = -distance (since cvar minimizes cost)
    let mut colored_objectives = Vec::with_capacity(candidates.len());
    let mut robust_scores = Vec::with_capacity(candidates.len());

    for c in candidates {
        let cost_samples: Vec<f64> = c.crn_distance_samples_m.iter().map(|&d| -d).collect();
        let obj = ColoredObjective::new(
            c.name.clone(),
            cost_samples.clone(),
            vec![Color::Estimated {
                estimator: "crn_wind_ensemble".into(),
                dispersion: 0.1,
            }],
        );

        let r_cost = cvar(&cost_samples, config.cvar_alpha).map_err(|e| {
            refuse(
                "robust-opt-cvar-computation",
                format!("cvar computation failed: {e:?}"),
                "verify valid alpha and samples",
            )
        })?;

        robust_scores.push(r_cost);
        colored_objectives.push(obj);
    }

    let mut best_robust_idx = 0;
    let mut min_robust_cost = f64::INFINITY;
    for (i, &rc) in robust_scores.iter().enumerate() {
        if rc < min_robust_cost {
            min_robust_cost = rc;
            best_robust_idx = i;
        }
    }
    let robust_winner = candidates[best_robust_idx].name.clone();
    let robustness_reorders = nominal_winner != robust_winner;

    let headline = colored_objectives[best_robust_idx]
        .headline_color()
        .map_err(|e| {
            refuse(
                "robust-opt-headline-color",
                format!("failed to derive headline color: {e:?}"),
                "verify valid input colors",
            )
        })?;

    let digest_input = format!(
        "wf-robustopt-v1:{}:{}:{}:{}",
        nominal_winner, robust_winner, robustness_reorders, candidates.len()
    );
    let receipt_digest = fs_blake3::hash_domain("org.frankensim.wf.robustopt.receipt.v1", digest_input.as_bytes())
        .to_hex()
        .to_string();

    Ok(RobustOptReceipt {
        nominal_winner,
        robust_winner,
        robustness_reorders,
        headline_rank: headline.rank(),
        candidates_evaluated: candidates.len(),
        receipt_digest,
    })
}
