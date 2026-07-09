//! fs-oed-e2e — SensorForge: optimal experimental design that knows when to
//! stop. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! You must pick the best of several designs, but their performances are only
//! estimated; you can spend sensors to sharpen them. Which do you measure, and
//! when have you measured enough? This answers both with certificates, composing
//! crates never designed to meet:
//!
//! - **Kalman fusion** ([`fs_assimilate`]): each candidate is a Gaussian belief;
//!   a sensor reading is fused with the exact scalar Kalman update, shrinking that
//!   candidate's posterior variance.
//! - **Value of information** ([`fs_voi`]): at each step the Expected Value of
//!   Perfect Information scores the decision's ambiguity; `recommend` places the
//!   next sensor on the candidate whose measurement most sharpens the DECISION
//!   (not the most-uncertain candidate), and says STOP the instant EVPI falls
//!   below threshold — the design choice is already robust.
//! - **Budget allocation** ([`fs_toleralloc`]): the measurement-precision budget
//!   is then distributed cost-optimally across candidates by sensitivity.
//! - **Honest colors** ([`fs_evidence`]): the posterior variance is `Verified`
//!   (an exact Kalman computation); the EVPI-driven stop is `Estimated`.
//!
//! Deterministic (sensor readings hit each candidate's true value; the Kalman
//! variance update is observation-independent). No dependencies beyond the
//! composed crates.

use fs_assimilate::{Belief, assimilate, point_sensor};
use fs_evidence::{Color, ColorRank};
use fs_toleralloc::{Feature, allocate};
use fs_voi::{Action, ActionKind, DesignEstimate, Recommendation, Uncertainty, evpi, recommend};

/// A candidate design under measurement.
#[derive(Debug, Clone)]
pub struct Candidate {
    /// Name.
    pub name: String,
    /// The (unknown to the planner) true performance a sensor would read.
    pub truth: f64,
    /// Prior mean estimate.
    pub prior_mean: f64,
    /// Prior variance.
    pub prior_var: f64,
    /// This candidate's sensor noise variance.
    pub sensor_noise: f64,
    /// Cost of one measurement.
    pub sensor_cost: f64,
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct OedReport {
    /// Candidate names in the order sensors were placed.
    pub placements: Vec<String>,
    /// Number of sensors placed.
    pub sensors_placed: usize,
    /// Total prior variance across candidates.
    pub prior_total_variance: f64,
    /// Total posterior variance across candidates.
    pub posterior_total_variance: f64,
    /// Fractional variance reduction.
    pub variance_reduction: f64,
    /// EVPI before any sensor.
    pub initial_evpi: f64,
    /// EVPI after the campaign stopped.
    pub final_evpi: f64,
    /// Did the decision become robust (planner chose to STOP)?
    pub decision_robust: bool,
    /// The finally-chosen (lowest-cost posterior) design.
    pub chosen_design: String,
    /// The cost-optimal tolerance allocation `(name, tolerance)`.
    pub allocation: Vec<(String, f64)>,
    /// The posterior-variance color (`Verified` — exact Kalman).
    pub variance_color: Color,
    /// The EVPI color (`Estimated` — decision-theoretic).
    pub evpi_color: Color,
}

fn to_estimates(candidates: &[Candidate], beliefs: &[Belief]) -> Vec<DesignEstimate> {
    candidates
        .iter()
        .zip(beliefs)
        .map(|(c, b)| {
            DesignEstimate::new(
                c.name.clone(),
                b.mean[0],
                Uncertainty {
                    numerical: b.variance(0).sqrt(),
                    statistical: 0.0,
                    model: 0.0,
                },
            )
        })
        .collect()
}

fn total_variance(beliefs: &[Belief]) -> f64 {
    beliefs.iter().map(|b| b.variance(0)).sum()
}

/// Run the SensorForge campaign; stop when EVPI ≤ `threshold` or after
/// `max_sensors` placements.
///
/// # Panics
/// If `candidates` is empty.
#[must_use]
pub fn run_campaign(candidates: &[Candidate], threshold: f64, max_sensors: usize) -> OedReport {
    assert!(!candidates.is_empty(), "need at least one candidate");
    let mut beliefs: Vec<Belief> = candidates
        .iter()
        .map(|c| Belief::scalar(c.prior_mean, c.prior_var))
        .collect();
    let prior_total_variance = total_variance(&beliefs);
    let initial_evpi = evpi(&to_estimates(candidates, &beliefs));

    let mut placements = Vec::new();
    let mut decision_robust = false;
    for _ in 0..max_sensors {
        let estimates = to_estimates(candidates, &beliefs);
        // One candidate measurement per design; measuring reduces its numerical
        // (discretization/estimation) uncertainty.
        let actions: Vec<Action> = candidates
            .iter()
            .map(|c| Action {
                name: format!("measure-{}", c.name),
                kind: ActionKind::Simulate,
                target_design: c.name.clone(),
                reduction: 0.85,
                cost: c.sensor_cost,
            })
            .collect();
        match recommend(&estimates, &actions, threshold) {
            Recommendation::Act { action, .. } => {
                let idx = actions
                    .iter()
                    .position(|a| a.name == action)
                    .expect("recommended action is in the list");
                let obs = point_sensor(
                    0,
                    1,
                    candidates[idx].truth,
                    candidates[idx].sensor_noise,
                    format!("sensor-{}", candidates[idx].name),
                );
                beliefs[idx] = assimilate(&beliefs[idx], &obs).expect("scalar assimilation");
                placements.push(candidates[idx].name.clone());
            }
            Recommendation::Stop { .. } => {
                decision_robust = true;
                break;
            }
        }
    }

    let estimates = to_estimates(candidates, &beliefs);
    let final_evpi = evpi(&estimates);
    let posterior_total_variance = total_variance(&beliefs);
    // fs-voi minimizes: the chosen design is the lowest-cost (best) posterior.
    let chosen = estimates
        .iter()
        .min_by(|a, b| a.mean.total_cmp(&b.mean))
        .map(|d| d.name.clone())
        .unwrap_or_default();

    // Cost-optimal precision allocation: sensitivity = prior std (the leverage
    // each sensor has), colored Verified (a modeled, exact quantity).
    let features: Vec<Feature> = candidates
        .iter()
        .map(|c| Feature {
            name: c.name.clone(),
            sensitivity: c.prior_var.sqrt(),
            sensitivity_color: ColorRank::Verified,
            cost_coeff: c.sensor_cost,
            baseline_tolerance: 0.1,
        })
        .collect();
    let allocation = allocate(&features, 0.02, 3.0)
        .map(|a| {
            a.items
                .into_iter()
                .map(|it| (it.name, it.tolerance))
                .collect()
        })
        .unwrap_or_default();

    OedReport {
        sensors_placed: placements.len(),
        placements,
        prior_total_variance,
        posterior_total_variance,
        variance_reduction: 1.0 - posterior_total_variance / prior_total_variance,
        initial_evpi,
        final_evpi,
        decision_robust,
        chosen_design: chosen,
        allocation,
        variance_color: Color::Verified {
            lo: posterior_total_variance,
            hi: posterior_total_variance,
        },
        evpi_color: Color::Estimated {
            estimator: "evpi-gaussian".to_string(),
            dispersion: final_evpi,
        },
    }
}

/// The worked scenario: four designs with uncertain COST (lower is better). The
/// two cheapest (A, B) are close and uncertain — the decision hinges on
/// measuring THEM, not the clearly-costlier C or D.
#[must_use]
pub fn demo_candidates() -> Vec<Candidate> {
    vec![
        Candidate {
            name: "A".into(),
            truth: 0.60,
            prior_mean: 0.60,
            prior_var: 0.10,
            sensor_noise: 0.01,
            sensor_cost: 1.0,
        },
        Candidate {
            name: "B".into(),
            truth: 0.65,
            prior_mean: 0.65,
            prior_var: 0.12,
            sensor_noise: 0.01,
            sensor_cost: 1.0,
        },
        Candidate {
            name: "C".into(),
            truth: 0.85,
            prior_mean: 0.85,
            prior_var: 0.06,
            sensor_noise: 0.01,
            sensor_cost: 1.0,
        },
        Candidate {
            name: "D".into(),
            truth: 1.10,
            prior_mean: 1.10,
            prior_var: 0.04,
            sensor_noise: 0.01,
            sensor_cost: 1.0,
        },
    ]
}
