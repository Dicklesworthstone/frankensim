//! fs-schedule-e2e — CampaignSchedule: bounded makespan analysis plus advisory
//! value-of-information scheduling. Layer: L6 (HELM).
//!
//! # The campaign
//!
//! A design campaign asks two orthogonal questions, and this answers both while
//! preserving their different evidence strengths:
//!
//! - **WHEN does it finish?** ([`fs_tropical`]) — the refinement studies form a
//!   precedence DAG with per-study latencies. The completion time is the
//!   longest weighted path — a MAX-PLUS (tropical) critical-path computation
//!   with directed-rounding bounds. It names a bottleneck only when those
//!   bounds prove the critical path is unique.
//! - **WHETHER to keep going?** ([`fs_voi`]) — the candidate designs have
//!   uncertain COST (lower is better) from numerical / statistical / model
//!   sources. The
//!   Expected Value of Perfect Information (EVPI) measures how much the current
//!   ranking ambiguity is worth resolving; `recommend` picks the highest
//!   value-per-cost study, OR says STOP when the decision is already robust
//!   (EVPI below threshold) — the anytime "don't gather data you don't need".
//! - **Honest colors** ([`fs_evidence`]) — the makespan is `Verified` by a
//!   rigorous enclosure; the EVPI-driven recommendation is `Estimated` (a
//!   decision-theoretic model).
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::{Color, NumericalCertificate, validate_color_payload, verified_from};
use fs_tropical::{MAX_TASK_DAG_EDGES, MAX_TASK_DAG_NODES, TaskDag, TropicalError};
use fs_voi::{Action, DesignEstimate, Recommendation, evpi, ranking_flip_probability, recommend};
use std::collections::BTreeSet;

/// Maximum designs or actions admitted to one campaign decision.
pub const MAX_CAMPAIGN_DECISION_ITEMS: usize = 4_096;
/// Maximum action-by-design evaluations admitted to one recommendation.
pub const MAX_CAMPAIGN_DECISION_WORK: usize = 65_536;
/// Maximum ASCII bytes in one campaign identity.
pub const MAX_CAMPAIGN_NAME_BYTES: usize = 128;

/// One refinement study: a node in the schedule DAG.
#[derive(Debug, Clone)]
pub struct Study {
    /// Study name.
    pub name: String,
    /// Latency (compute cost / wall-clock units).
    pub latency: f64,
    /// Prerequisite study indices (must finish first).
    pub deps: Vec<usize>,
}

impl Study {
    /// A study.
    #[must_use]
    pub fn new(name: impl Into<String>, latency: f64, deps: Vec<usize>) -> Study {
        Study {
            name: name.into(),
            latency,
            deps,
        }
    }
}

/// The campaign report.
#[derive(Debug, Clone, PartialEq)]
pub struct ScheduleReport {
    /// The campaign makespan (tropical critical-path length).
    pub makespan: f64,
    /// The makespan's color (`Verified` — outward-rounded tropical enclosure).
    pub makespan_color: Color,
    /// The critical path (study indices, in order).
    pub critical_path: Vec<usize>,
    /// Whether directed bounds prove the returned critical path is unique.
    pub critical_path_is_unique: bool,
    /// Nominal per-study slack, in input order.
    pub slack: Vec<f64>,
    /// Bottleneck index when the critical path is certified unique.
    pub bottleneck_index: Option<usize>,
    /// Bottleneck study name when the critical path is certified unique.
    pub bottleneck: Option<String>,
    /// Studies with positive nominal slack. This is scheduling guidance, not a
    /// certified deferability claim.
    pub slack_studies: Vec<String>,
    /// The current Expected Value of Perfect Information (decision ambiguity).
    pub evpi: f64,
    /// The leading (lowest-cost) design — `fs-voi` minimizes.
    pub leading_design: String,
    /// Probability the ranking flips between the top two designs.
    pub flip_risk: f64,
    /// The VoI recommendation: `Act: <study>` or `Stop: <reason>`.
    pub recommendation: String,
    /// Typed disposition; a deficient action menu is not a robust stop.
    pub disposition: ScheduleDisposition,
    /// Decision-model evidence for the recommendation.
    pub recommendation_color: Color,
    /// Value per cost for an [`ScheduleDisposition::Act`] recommendation.
    pub recommendation_value_per_cost: Option<f64>,
    /// Whether the campaign should STOP (decision already robust).
    pub should_stop: bool,
}

/// Typed campaign-decision disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScheduleDisposition {
    /// Execute the selected information-acquisition action.
    Act,
    /// Stop because EVPI is at or below the declared threshold.
    RobustStop,
    /// EVPI remains material, but the supplied menu has no effective action.
    NoEffectiveAction,
}

/// Structured campaign-scheduling refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScheduleError {
    /// At least one scheduled study is required.
    NoStudies,
    /// At least one candidate design is required.
    NoDesigns,
    /// A deterministic collection limit was exceeded.
    ResourceLimit {
        /// Bounded resource.
        resource: &'static str,
        /// Maximum admitted count.
        limit: usize,
        /// Observed count or conservative lower bound.
        observed: usize,
    },
    /// One indexed field is outside its finite/canonical domain.
    InvalidField {
        /// Record family.
        record: &'static str,
        /// Record index.
        index: usize,
        /// Field name.
        field: &'static str,
        /// Stable requirement.
        requirement: &'static str,
    },
    /// A same-family name is duplicated and therefore ambiguous.
    DuplicateName {
        /// Record family.
        record: &'static str,
        /// Duplicate bounded name.
        name: String,
    },
    /// An action names no admitted design.
    UnknownActionTarget {
        /// Action index.
        action: usize,
    },
    /// Lower max-plus analysis refused the graph.
    Tropical(TropicalError),
    /// Decision algebra produced a non-finite or out-of-range result.
    InvalidDecisionResult {
        /// Result field.
        field: &'static str,
    },
    /// Makespan or recommendation evidence could not be admitted.
    EvidenceRefused,
}

impl core::fmt::Display for ScheduleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoStudies => f.write_str("campaign schedule requires at least one study"),
            Self::NoDesigns => f.write_str("campaign schedule requires at least one design"),
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "campaign {resource} limit {limit} exceeded (observed {observed})"
            ),
            Self::InvalidField {
                record,
                index,
                field,
                requirement,
            } => write!(f, "campaign {record} {index} field {field} {requirement}"),
            Self::DuplicateName { record, name } => {
                write!(f, "campaign {record} name {name:?} is duplicated")
            }
            Self::UnknownActionTarget { action } => {
                write!(f, "campaign action {action} targets no admitted design")
            }
            Self::Tropical(error) => write!(f, "campaign schedule graph refused: {error}"),
            Self::InvalidDecisionResult { field } => {
                write!(f, "campaign decision produced invalid {field}")
            }
            Self::EvidenceRefused => f.write_str("campaign evidence payload was refused"),
        }
    }
}

impl std::error::Error for ScheduleError {}

impl From<TropicalError> for ScheduleError {
    fn from(value: TropicalError) -> Self {
        Self::Tropical(value)
    }
}

fn validate_name(record: &'static str, index: usize, name: &str) -> Result<(), ScheduleError> {
    if name.is_empty()
        || name.len() > MAX_CAMPAIGN_NAME_BYTES
        || !name.bytes().all(|byte| byte.is_ascii_graphic())
    {
        return Err(ScheduleError::InvalidField {
            record,
            index,
            field: "name",
            requirement: "must contain 1..=128 ASCII graphic bytes",
        });
    }
    Ok(())
}

/// Run the CampaignSchedule campaign.
///
/// # Errors
/// Refuses empty, oversized, non-canonical, non-finite, ambiguous, cyclic, or
/// numerically overflowing campaign inputs before minting positive evidence.
#[allow(clippy::too_many_lines)] // One ordered cross-model admission boundary.
pub fn run_campaign(
    studies: &[Study],
    designs: &[DesignEstimate],
    actions: &[Action],
    stop_threshold: f64,
) -> Result<ScheduleReport, ScheduleError> {
    if studies.is_empty() {
        return Err(ScheduleError::NoStudies);
    }
    if designs.is_empty() {
        return Err(ScheduleError::NoDesigns);
    }
    for (resource, observed, limit) in [
        ("studies", studies.len(), MAX_TASK_DAG_NODES),
        ("designs", designs.len(), MAX_CAMPAIGN_DECISION_ITEMS),
        ("actions", actions.len(), MAX_CAMPAIGN_DECISION_ITEMS),
    ] {
        if observed > limit {
            return Err(ScheduleError::ResourceLimit {
                resource,
                limit,
                observed,
            });
        }
    }
    let decision_work =
        designs
            .len()
            .checked_mul(actions.len())
            .ok_or(ScheduleError::ResourceLimit {
                resource: "action-design evaluations",
                limit: MAX_CAMPAIGN_DECISION_WORK,
                observed: usize::MAX,
            })?;
    if decision_work > MAX_CAMPAIGN_DECISION_WORK {
        return Err(ScheduleError::ResourceLimit {
            resource: "action-design evaluations",
            limit: MAX_CAMPAIGN_DECISION_WORK,
            observed: decision_work,
        });
    }
    if !stop_threshold.is_finite() || stop_threshold < 0.0 {
        return Err(ScheduleError::InvalidField {
            record: "campaign",
            index: 0,
            field: "stop_threshold",
            requirement: "must be finite and non-negative",
        });
    }
    let mut study_names = BTreeSet::new();
    let mut dependency_count = 0usize;
    for (index, study) in studies.iter().enumerate() {
        validate_name("study", index, &study.name)?;
        if !study_names.insert(study.name.as_str()) {
            return Err(ScheduleError::DuplicateName {
                record: "study",
                name: study.name.clone(),
            });
        }
        if !study.latency.is_finite() || study.latency < 0.0 {
            return Err(ScheduleError::InvalidField {
                record: "study",
                index,
                field: "latency",
                requirement: "must be finite and non-negative",
            });
        }
        dependency_count =
            dependency_count
                .checked_add(study.deps.len())
                .ok_or(ScheduleError::ResourceLimit {
                    resource: "dependencies",
                    limit: MAX_TASK_DAG_EDGES,
                    observed: usize::MAX,
                })?;
    }
    if dependency_count > MAX_TASK_DAG_EDGES {
        return Err(ScheduleError::ResourceLimit {
            resource: "dependencies",
            limit: MAX_TASK_DAG_EDGES,
            observed: dependency_count,
        });
    }
    let mut design_names = BTreeSet::new();
    for (index, design) in designs.iter().enumerate() {
        validate_name("design", index, &design.name)?;
        if !design_names.insert(design.name.as_str()) {
            return Err(ScheduleError::DuplicateName {
                record: "design",
                name: design.name.clone(),
            });
        }
        if !design.mean.is_finite() {
            return Err(ScheduleError::InvalidField {
                record: "design",
                index,
                field: "mean",
                requirement: "must be finite",
            });
        }
        for (field, value) in [
            ("uncertainty.numerical", design.uncertainty.numerical),
            ("uncertainty.statistical", design.uncertainty.statistical),
            ("uncertainty.model", design.uncertainty.model),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(ScheduleError::InvalidField {
                    record: "design",
                    index,
                    field,
                    requirement: "must be finite and non-negative",
                });
            }
        }
        if !design.uncertainty.total_std().is_finite() {
            return Err(ScheduleError::InvalidField {
                record: "design",
                index,
                field: "uncertainty",
                requirement: "must have a finite quadrature total",
            });
        }
    }
    let mut action_names = BTreeSet::new();
    for (index, action) in actions.iter().enumerate() {
        validate_name("action", index, &action.name)?;
        validate_name("action target", index, &action.target_design)?;
        if !action_names.insert(action.name.as_str()) {
            return Err(ScheduleError::DuplicateName {
                record: "action",
                name: action.name.clone(),
            });
        }
        if !design_names.contains(action.target_design.as_str()) {
            return Err(ScheduleError::UnknownActionTarget { action: index });
        }
        if !action.reduction.is_finite() || !(0.0..=1.0).contains(&action.reduction) {
            return Err(ScheduleError::InvalidField {
                record: "action",
                index,
                field: "reduction",
                requirement: "must be finite and in 0..=1",
            });
        }
        if !action.cost.is_finite() || action.cost <= 0.0 {
            return Err(ScheduleError::InvalidField {
                record: "action",
                index,
                field: "cost",
                requirement: "must be finite and positive",
            });
        }
    }

    // --- WHEN: the tropical critical path over the study precedence DAG. ---
    let latencies: Vec<f64> = studies.iter().map(|s| s.latency).collect();
    let mut dag = TaskDag::new(latencies);
    for (v, s) in studies.iter().enumerate() {
        for &u in &s.deps {
            dag = dag.with_edge(u, v);
        }
    }
    let cp = dag.critical_path()?;
    let bottleneck_index = dag.bottleneck()?;
    let bottleneck = bottleneck_index.map(|i| studies[i].name.clone());
    let slack_studies: Vec<String> = cp
        .slack
        .iter()
        .enumerate()
        .filter(|(_, slack)| **slack > 0.0)
        .map(|(i, _)| studies[i].name.clone())
        .collect();
    let makespan_color = verified_from(&NumericalCertificate::enclosure(
        cp.makespan_lo,
        cp.makespan_hi,
    ))
    .map_err(|_| ScheduleError::EvidenceRefused)?;

    // --- WHETHER: value of information over the candidate designs. ---
    let current_evpi = evpi(designs);
    if !current_evpi.is_finite() || current_evpi < 0.0 {
        return Err(ScheduleError::InvalidDecisionResult { field: "EVPI" });
    }
    // Leading design + the top-two flip risk (the decision's fragility).
    // fs-voi MINIMIZES, so the leader is the lowest-cost design.
    let mut ranked: Vec<&DesignEstimate> = designs.iter().collect();
    ranked.sort_by(|a, b| a.mean.total_cmp(&b.mean));
    let leading_design = ranked.first().map(|d| d.name.clone()).unwrap_or_default();
    let flip_risk = if ranked.len() >= 2 {
        ranking_flip_probability(ranked[0], ranked[1])
    } else {
        0.0
    };
    if !flip_risk.is_finite() || !(0.0..=1.0).contains(&flip_risk) {
        return Err(ScheduleError::InvalidDecisionResult {
            field: "ranking flip probability",
        });
    }
    let (recommendation, disposition, recommendation_value_per_cost) =
        match recommend(designs, actions, stop_threshold) {
            Recommendation::Act {
                action,
                value_per_cost,
            } if value_per_cost.is_finite() && value_per_cost > 0.0 => (
                format!("Act: {action} (value/cost {value_per_cost:.3})"),
                ScheduleDisposition::Act,
                Some(value_per_cost),
            ),
            Recommendation::Act { .. } => {
                return Err(ScheduleError::InvalidDecisionResult {
                    field: "action value per cost",
                });
            }
            Recommendation::Stop { reason } if current_evpi <= stop_threshold => (
                format!("Stop: {reason}"),
                ScheduleDisposition::RobustStop,
                None,
            ),
            Recommendation::Stop { reason } => (
                format!("Expand menu: {reason}"),
                ScheduleDisposition::NoEffectiveAction,
                None,
            ),
        };
    let recommendation_color = Color::Estimated {
        estimator: "fs-voi-gaussian-recommendation-v1".to_string(),
        dispersion: current_evpi,
    };
    validate_color_payload(&recommendation_color).map_err(|_| ScheduleError::EvidenceRefused)?;
    let should_stop = disposition == ScheduleDisposition::RobustStop;

    Ok(ScheduleReport {
        makespan: cp.makespan,
        makespan_color,
        critical_path: cp.path,
        critical_path_is_unique: cp.path_is_unique,
        slack: cp.slack,
        bottleneck_index,
        bottleneck,
        slack_studies,
        evpi: current_evpi,
        leading_design,
        flip_risk,
        recommendation,
        disposition,
        recommendation_color,
        recommendation_value_per_cost,
        should_stop,
    })
}
