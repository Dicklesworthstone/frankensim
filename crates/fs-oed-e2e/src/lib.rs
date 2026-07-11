//! fs-oed-e2e — SensorForge: optimal experimental design that knows when to
//! stop. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! You must pick the best of several designs, but their performances are only
//! estimated; you can spend sensors to sharpen them. Which do you measure, and
//! when have you measured enough? This answers both with evidence, composing
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
//! - **Honest colors** ([`fs_evidence`]): posterior variance and EVPI remain
//!   `Estimated`; their bounded identities commit to every campaign input and
//!   every instrument-bound assimilation candidate.
//!
//! Deterministic (sensor readings hit each candidate's true value; the Kalman
//! variance update is observation-independent). No dependencies beyond the
//! composed crates.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_assimilate::{AssimError, Belief, assimilate_colored, point_sensor};
use fs_evidence::{Color, ColorRank, color_leaf_identity_reason};
use fs_toleralloc::{Feature, allocate};
use fs_voi::{Action, ActionKind, DesignEstimate, Recommendation, Uncertainty, evpi, recommend};

/// Maximum accepted candidate-name length.
pub const MAX_CANDIDATE_NAME_BYTES: usize = 128;
/// Maximum number of candidates in one synchronous campaign.
pub const MAX_CAMPAIGN_CANDIDATES: usize = 256;
/// Maximum number of sensor placements in one synchronous campaign.
pub const MAX_CAMPAIGN_SENSORS: usize = 4_096;
/// Maximum `candidate_count^2 * max_sensors` action-design evaluations.
pub const MAX_CAMPAIGN_EVALUATIONS: usize = 1_000_000;

const REPORT_ID_DOMAIN: &str = "org.frankensim.fs-oed-e2e.report.v2";

fn canonicalize_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

/// A rejected candidate declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CandidateError {
    /// The candidate name cannot serve as a bounded provenance identity.
    InvalidName {
        /// Structural rejection reason.
        reason: &'static str,
    },
    /// A numeric field violates its declared domain.
    InvalidNumber {
        /// Offending field.
        field: &'static str,
        /// Required domain.
        requirement: &'static str,
    },
}

impl fmt::Display for CandidateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName { reason } => {
                write!(f, "candidate name is not an admissible identity: {reason}")
            }
            Self::InvalidNumber { field, requirement } => {
                write!(f, "candidate `{field}` must be {requirement}")
            }
        }
    }
}

impl std::error::Error for CandidateError {}

/// A candidate design under measurement.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    name: String,
    truth: f64,
    prior_mean: f64,
    prior_var: f64,
    sensor_noise: f64,
    sensor_cost: f64,
}

impl Candidate {
    /// Construct a checked candidate.
    ///
    /// # Errors
    /// Returns [`CandidateError`] for an unusable name, a non-finite numeric
    /// field, negative prior variance, or non-positive sensor noise/cost.
    pub fn new(
        name: impl Into<String>,
        truth: f64,
        prior_mean: f64,
        prior_var: f64,
        sensor_noise: f64,
        sensor_cost: f64,
    ) -> Result<Self, CandidateError> {
        let name = name.into();
        let name_reason = if name.len() > MAX_CANDIDATE_NAME_BYTES {
            Some("too-long")
        } else {
            color_leaf_identity_reason(&name)
        };
        if let Some(reason) = name_reason {
            return Err(CandidateError::InvalidName { reason });
        }
        for (field, value) in [("truth", truth), ("prior_mean", prior_mean)] {
            if !value.is_finite() {
                return Err(CandidateError::InvalidNumber {
                    field,
                    requirement: "finite",
                });
            }
        }
        if !prior_var.is_finite() || prior_var < 0.0 {
            return Err(CandidateError::InvalidNumber {
                field: "prior_var",
                requirement: "finite and non-negative",
            });
        }
        for (field, value) in [("sensor_noise", sensor_noise), ("sensor_cost", sensor_cost)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(CandidateError::InvalidNumber {
                    field,
                    requirement: "finite and positive",
                });
            }
        }
        Ok(Self {
            name,
            truth: canonicalize_zero(truth),
            prior_mean: canonicalize_zero(prior_mean),
            prior_var: canonicalize_zero(prior_var),
            sensor_noise: canonicalize_zero(sensor_noise),
            sensor_cost: canonicalize_zero(sensor_cost),
        })
    }

    /// Candidate identity.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Sensor reading used by this deterministic worked campaign.
    #[must_use]
    pub fn truth(&self) -> f64 {
        self.truth
    }

    /// Prior objective mean.
    #[must_use]
    pub fn prior_mean(&self) -> f64 {
        self.prior_mean
    }

    /// Prior objective variance.
    #[must_use]
    pub fn prior_variance(&self) -> f64 {
        self.prior_var
    }

    /// Sensor noise variance.
    #[must_use]
    pub fn sensor_noise(&self) -> f64 {
        self.sensor_noise
    }

    /// Cost of one measurement.
    #[must_use]
    pub fn sensor_cost(&self) -> f64 {
        self.sensor_cost
    }
}

/// A rejected campaign or failed campaign computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OedError {
    /// At least one candidate is required.
    NoCandidates,
    /// The synchronous candidate cap was exceeded.
    TooManyCandidates {
        /// Requested count.
        count: usize,
        /// Accepted maximum.
        max: usize,
    },
    /// The synchronous placement cap was exceeded.
    TooManySensors {
        /// Requested count.
        count: usize,
        /// Accepted maximum.
        max: usize,
    },
    /// The requested planning work exceeds the synchronous campaign budget.
    WorkBudgetExceeded {
        /// Candidate count.
        candidates: usize,
        /// Requested placement cap.
        max_sensors: usize,
        /// Requested action-design evaluations.
        evaluations: usize,
        /// Accepted maximum product.
        max_evaluations: usize,
    },
    /// The EVPI stop threshold must be finite and non-negative.
    InvalidThreshold,
    /// Candidate identities must be unique because actions address them by name.
    DuplicateCandidate {
        /// Repeated identity.
        name: String,
    },
    /// A checked scalar belief unexpectedly rejected an internal access.
    BeliefInvariant(AssimError),
    /// An observation or posterior update failed.
    Assimilation {
        /// Candidate being measured.
        candidate: String,
        /// Structured lower-layer failure.
        source: AssimError,
    },
    /// `fs-voi` returned an action outside the menu it was given.
    UnknownRecommendation {
        /// Returned action identity.
        action: String,
    },
    /// A deterministic derived quantity overflowed or became NaN.
    NonFiniteComputation {
        /// Quantity whose contract failed.
        quantity: &'static str,
    },
    /// The tolerance allocator rejected checked positive-sensitivity inputs.
    AllocationFailed,
    /// The allocator omitted a checked positive-sensitivity candidate.
    MissingAllocation {
        /// Missing candidate identity.
        candidate: String,
    },
}

impl fmt::Display for OedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoCandidates => write!(f, "SensorForge needs at least one candidate"),
            Self::TooManyCandidates { count, max } => {
                write!(f, "candidate count {count} exceeds synchronous cap {max}")
            }
            Self::TooManySensors { count, max } => {
                write!(f, "sensor cap {count} exceeds synchronous cap {max}")
            }
            Self::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations,
                max_evaluations,
            } => write!(
                f,
                "campaign work {candidates}^2 x {max_sensors} = {evaluations} exceeds \
                 {max_evaluations} action-design evaluations"
            ),
            Self::InvalidThreshold => {
                write!(f, "EVPI threshold must be finite and non-negative")
            }
            Self::DuplicateCandidate { name } => {
                write!(f, "candidate identity `{name}` is duplicated")
            }
            Self::BeliefInvariant(source) => write!(f, "scalar belief invariant failed: {source}"),
            Self::Assimilation { candidate, source } => {
                write!(
                    f,
                    "assimilation failed for candidate `{candidate}`: {source}"
                )
            }
            Self::UnknownRecommendation { action } => {
                write!(f, "fs-voi returned unknown action `{action}`")
            }
            Self::NonFiniteComputation { quantity } => {
                write!(f, "campaign produced non-finite `{quantity}`")
            }
            Self::AllocationFailed => write!(f, "precision allocation failed"),
            Self::MissingAllocation { candidate } => {
                write!(f, "precision allocation omitted candidate `{candidate}`")
            }
        }
    }
}

impl std::error::Error for OedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::BeliefInvariant(source) | Self::Assimilation { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Final scalar posterior for one candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct PosteriorSummary {
    /// Candidate identity.
    pub name: String,
    /// Posterior mean.
    pub mean: f64,
    /// Posterior variance.
    pub variance: f64,
}

/// The campaign report.
#[derive(Debug, Clone, PartialEq)]
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
    /// A zero-sensitivity candidate receives `+infinity`, the exact unconstrained
    /// optimum under the first-order allocation model.
    pub allocation: Vec<(String, f64)>,
    /// EVPI before sensing and after every completed placement.
    pub evpi_trace: Vec<f64>,
    /// Final scalar posterior in candidate order.
    pub posteriors: Vec<PosteriorSummary>,
    /// Instrument-bound estimated candidate emitted by each assimilation.
    pub assimilation_colors: Vec<Color>,
    /// The posterior-variance color (`Estimated` until independently certified).
    pub variance_color: Color,
    /// The EVPI color (`Estimated` — decision-theoretic).
    pub evpi_color: Color,
}

fn to_estimates(
    candidates: &[Candidate],
    beliefs: &[Belief],
) -> Result<Vec<DesignEstimate>, OedError> {
    if candidates.len() != beliefs.len() {
        return Err(OedError::NonFiniteComputation {
            quantity: "candidate/belief cardinality",
        });
    }
    candidates
        .iter()
        .zip(beliefs)
        .map(|(c, b)| {
            let mean = b.component_mean(0).map_err(OedError::BeliefInvariant)?;
            let variance = b.variance(0).map_err(OedError::BeliefInvariant)?;
            Ok(DesignEstimate::new(
                c.name.clone(),
                mean,
                Uncertainty {
                    numerical: variance.sqrt(),
                    statistical: 0.0,
                    model: 0.0,
                },
            ))
        })
        .collect()
}

fn total_variance(beliefs: &[Belief]) -> Result<f64, OedError> {
    beliefs.iter().try_fold(0.0, |total, belief| {
        let variance = belief.variance(0).map_err(OedError::BeliefInvariant)?;
        let next = total + variance;
        if next.is_finite() {
            Ok(next)
        } else {
            Err(OedError::NonFiniteComputation {
                quantity: "total variance",
            })
        }
    })
}

fn checked_evpi(estimates: &[DesignEstimate]) -> Result<f64, OedError> {
    let value = evpi(estimates);
    if value.is_finite() && value >= 0.0 {
        Ok(canonicalize_zero(value))
    } else {
        Err(OedError::NonFiniteComputation { quantity: "EVPI" })
    }
}

fn precision_allocation(candidates: &[Candidate]) -> Result<Vec<(String, f64)>, OedError> {
    let features: Vec<Feature> = candidates
        .iter()
        .filter(|candidate| candidate.prior_var > 0.0)
        .map(|candidate| Feature {
            name: candidate.name.clone(),
            sensitivity: candidate.prior_var.sqrt(),
            sensitivity_color: ColorRank::Estimated,
            cost_coeff: candidate.sensor_cost,
            baseline_tolerance: 0.1,
        })
        .collect();
    let allocated: BTreeMap<String, f64> = if features.is_empty() {
        BTreeMap::new()
    } else {
        allocate(&features, 0.02, 3.0)
            .map_err(|_| OedError::AllocationFailed)?
            .items
            .into_iter()
            .map(|item| (item.name, item.tolerance))
            .collect()
    };

    candidates
        .iter()
        .map(|candidate| {
            if candidate.prior_var == 0.0 {
                Ok((candidate.name.clone(), f64::INFINITY))
            } else {
                let tolerance = allocated.get(&candidate.name).copied().ok_or_else(|| {
                    OedError::MissingAllocation {
                        candidate: candidate.name.clone(),
                    }
                })?;
                if !tolerance.is_finite() || tolerance <= 0.0 {
                    return Err(OedError::NonFiniteComputation {
                        quantity: "allocated tolerance",
                    });
                }
                Ok((candidate.name.clone(), tolerance))
            }
        })
        .collect()
}

fn push_bytes(output: &mut Vec<u8>, value: &[u8]) {
    output.extend_from_slice(&(value.len() as u64).to_le_bytes());
    output.extend_from_slice(value);
}

fn push_str(output: &mut Vec<u8>, value: &str) {
    push_bytes(output, value.as_bytes());
}

fn report_identity(
    quantity: &str,
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    placements: &[String],
    posteriors: &[PosteriorSummary],
    assimilation_colors: &[Color],
) -> String {
    let mut canonical = Vec::new();
    push_str(&mut canonical, quantity);
    canonical.extend_from_slice(&(candidates.len() as u64).to_le_bytes());
    for candidate in candidates {
        push_str(&mut canonical, &candidate.name);
        for value in [
            candidate.truth,
            candidate.prior_mean,
            candidate.prior_var,
            candidate.sensor_noise,
            candidate.sensor_cost,
        ] {
            canonical.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    canonical.extend_from_slice(&threshold.to_bits().to_le_bytes());
    canonical.extend_from_slice(&(max_sensors as u64).to_le_bytes());
    canonical.extend_from_slice(&(placements.len() as u64).to_le_bytes());
    for placement in placements {
        push_str(&mut canonical, placement);
    }
    canonical.extend_from_slice(&(posteriors.len() as u64).to_le_bytes());
    for posterior in posteriors {
        push_str(&mut canonical, &posterior.name);
        canonical.extend_from_slice(&posterior.mean.to_bits().to_le_bytes());
        canonical.extend_from_slice(&posterior.variance.to_bits().to_le_bytes());
    }
    canonical.extend_from_slice(&(assimilation_colors.len() as u64).to_le_bytes());
    for color in assimilation_colors {
        push_bytes(&mut canonical, &color.canonical_bytes());
    }
    let identity = format!(
        "sensorforge-{quantity}:v2:{}",
        fs_blake3::hash_domain(REPORT_ID_DOMAIN, &canonical)
    );
    debug_assert!(color_leaf_identity_reason(&identity).is_none());
    identity
}

fn validate_campaign(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
) -> Result<f64, OedError> {
    if candidates.is_empty() {
        return Err(OedError::NoCandidates);
    }
    if candidates.len() > MAX_CAMPAIGN_CANDIDATES {
        return Err(OedError::TooManyCandidates {
            count: candidates.len(),
            max: MAX_CAMPAIGN_CANDIDATES,
        });
    }
    if max_sensors > MAX_CAMPAIGN_SENSORS {
        return Err(OedError::TooManySensors {
            count: max_sensors,
            max: MAX_CAMPAIGN_SENSORS,
        });
    }
    let action_design_pairs =
        candidates
            .len()
            .checked_mul(candidates.len())
            .ok_or(OedError::WorkBudgetExceeded {
                candidates: candidates.len(),
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
    let work =
        action_design_pairs
            .checked_mul(max_sensors)
            .ok_or(OedError::WorkBudgetExceeded {
                candidates: candidates.len(),
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
    if work > MAX_CAMPAIGN_EVALUATIONS {
        return Err(OedError::WorkBudgetExceeded {
            candidates: candidates.len(),
            max_sensors,
            evaluations: work,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
        });
    }
    if !threshold.is_finite() || threshold < 0.0 {
        return Err(OedError::InvalidThreshold);
    }
    let mut names = BTreeSet::new();
    for candidate in candidates {
        if !names.insert(candidate.name.as_str()) {
            return Err(OedError::DuplicateCandidate {
                name: candidate.name.clone(),
            });
        }
    }
    Ok(canonicalize_zero(threshold))
}

struct CampaignState {
    beliefs: Vec<Belief>,
    placements: Vec<String>,
    assimilation_colors: Vec<Color>,
    evpi_trace: Vec<f64>,
    decision_robust: bool,
}

fn execute_placements(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    mut beliefs: Vec<Belief>,
    initial_evpi: f64,
) -> Result<CampaignState, OedError> {
    let actions: Vec<Action> = candidates
        .iter()
        .map(|candidate| Action {
            name: format!("measure-{}", candidate.name),
            kind: ActionKind::Simulate,
            target_design: candidate.name.clone(),
            reduction: 0.85,
            cost: candidate.sensor_cost,
        })
        .collect();
    let mut placements = Vec::new();
    let mut assimilation_colors = Vec::new();
    let mut evpi_trace = vec![initial_evpi];
    let mut decision_robust = false;

    loop {
        let estimates = to_estimates(candidates, &beliefs)?;
        if checked_evpi(&estimates)? <= threshold {
            decision_robust = true;
            break;
        }
        if placements.len() >= max_sensors {
            break;
        }
        let Recommendation::Act { action, .. } = recommend(&estimates, &actions, threshold) else {
            break;
        };
        let idx = actions
            .iter()
            .position(|candidate| candidate.name == action)
            .ok_or(OedError::UnknownRecommendation { action })?;
        let observation = point_sensor(
            0,
            1,
            candidates[idx].truth,
            candidates[idx].sensor_noise,
            format!("sensor-{}", candidates[idx].name),
        )
        .map_err(|source| OedError::Assimilation {
            candidate: candidates[idx].name.clone(),
            source,
        })?;
        let next_count = placements.len() + 1;
        let posterior = assimilate_colored(
            &beliefs[idx],
            std::slice::from_ref(&observation),
            "sensor_count",
            0.0,
            next_count as f64,
        )
        .map_err(|source| OedError::Assimilation {
            candidate: candidates[idx].name.clone(),
            source,
        })?;
        beliefs[idx] = posterior.belief().clone();
        assimilation_colors.push(posterior.color().clone());
        placements.push(candidates[idx].name.clone());
        evpi_trace.push(checked_evpi(&to_estimates(candidates, &beliefs)?)?);
    }

    Ok(CampaignState {
        beliefs,
        placements,
        assimilation_colors,
        evpi_trace,
        decision_robust,
    })
}

fn finish_report(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    prior_total_variance: f64,
    initial_evpi: f64,
    state: CampaignState,
) -> Result<OedReport, OedError> {
    let estimates = to_estimates(candidates, &state.beliefs)?;
    let final_evpi = checked_evpi(&estimates)?;
    let posterior_total_variance = total_variance(&state.beliefs)?;
    let chosen_design = estimates
        .iter()
        .min_by(|a, b| a.mean.total_cmp(&b.mean))
        .map(|design| design.name.clone())
        .ok_or(OedError::NonFiniteComputation {
            quantity: "chosen design",
        })?;
    let allocation = precision_allocation(candidates)?;
    let posteriors = candidates
        .iter()
        .zip(&state.beliefs)
        .map(|(candidate, belief)| {
            Ok(PosteriorSummary {
                name: candidate.name.clone(),
                mean: belief
                    .component_mean(0)
                    .map_err(OedError::BeliefInvariant)?,
                variance: belief.variance(0).map_err(OedError::BeliefInvariant)?,
            })
        })
        .collect::<Result<Vec<_>, OedError>>()?;
    let variance_reduction = if prior_total_variance == 0.0 {
        0.0
    } else {
        let reduction = (prior_total_variance - posterior_total_variance) / prior_total_variance;
        if !reduction.is_finite() {
            return Err(OedError::NonFiniteComputation {
                quantity: "variance reduction",
            });
        }
        canonicalize_zero(reduction)
    };
    let variance_identity = report_identity(
        "posterior-variance",
        candidates,
        threshold,
        max_sensors,
        &state.placements,
        &posteriors,
        &state.assimilation_colors,
    );
    let evpi_identity = report_identity(
        "evpi",
        candidates,
        threshold,
        max_sensors,
        &state.placements,
        &posteriors,
        &state.assimilation_colors,
    );

    Ok(OedReport {
        sensors_placed: state.placements.len(),
        placements: state.placements,
        prior_total_variance,
        posterior_total_variance,
        variance_reduction,
        initial_evpi,
        final_evpi,
        decision_robust: state.decision_robust,
        chosen_design,
        allocation,
        evpi_trace: state.evpi_trace,
        posteriors,
        assimilation_colors: state.assimilation_colors,
        variance_color: Color::Estimated {
            estimator: variance_identity,
            dispersion: f64::INFINITY,
        },
        evpi_color: Color::Estimated {
            estimator: evpi_identity,
            dispersion: final_evpi,
        },
    })
}

/// Run the SensorForge campaign; stop when EVPI <= `threshold` or after
/// `max_sensors` placements.
///
/// The initial STOP condition is evaluated even when `max_sensors == 0`.
///
/// # Errors
/// Returns [`OedError`] for invalid campaign bounds, duplicate candidate names,
/// a lower-layer assimilation/allocation failure, or a non-finite derived value.
pub fn run_campaign(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
) -> Result<OedReport, OedError> {
    let threshold = validate_campaign(candidates, threshold, max_sensors)?;
    let beliefs: Vec<Belief> = candidates
        .iter()
        .map(|c| Belief::scalar(c.prior_mean, c.prior_var))
        .collect::<Result<Vec<_>, _>>()
        .map_err(OedError::BeliefInvariant)?;
    let prior_total_variance = total_variance(&beliefs)?;
    let initial_evpi = checked_evpi(&to_estimates(candidates, &beliefs)?)?;
    let state = execute_placements(candidates, threshold, max_sensors, beliefs, initial_evpi)?;
    finish_report(
        candidates,
        threshold,
        max_sensors,
        prior_total_variance,
        initial_evpi,
        state,
    )
}

/// The worked scenario: four designs with uncertain COST (lower is better). The
/// two cheapest (A, B) are close and uncertain — the decision hinges on
/// measuring THEM, not the clearly-costlier C or D.
pub fn demo_candidates() -> Result<Vec<Candidate>, CandidateError> {
    [
        ("A", 0.60, 0.60, 0.10, 0.01, 1.0),
        ("B", 0.65, 0.65, 0.12, 0.01, 1.0),
        ("C", 0.85, 0.85, 0.06, 0.01, 1.0),
        ("D", 1.10, 1.10, 0.04, 0.01, 1.0),
    ]
    .into_iter()
    .map(
        |(name, truth, prior_mean, prior_var, sensor_noise, sensor_cost)| {
            Candidate::new(
                name,
                truth,
                prior_mean,
                prior_var,
                sensor_noise,
                sensor_cost,
            )
        },
    )
    .collect()
}
