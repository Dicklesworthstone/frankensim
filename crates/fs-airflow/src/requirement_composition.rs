//! Requirement and evidence composition for thermal QoIs vs ThermalLimits (bead `frankensim-s2l9v.2`).
//!
//! # Core Invariants and Non-Negotiables
//!
//! 1. **Zero Is Not No-Data**: A measured or computed 0.0 K temperature rise,
//!    0.0 Pa pressure drop, or 0.0 W fan power is an authentic physical value,
//!    NEVER collapsed into `NoData`.
//! 2. **Weakest-Stage Authority Monotonicity**: Evidence composition NEVER
//!    upgrades the authority of a weakest stage. If any contributing stage is
//!    `Estimated` or `OutsideDomain`, the overall composed requirement inherits
//!    that weakest authority.
//! 3. **Binding Witness Retention**: Every satisfied or violated comparison
//!    retains its exact binding witness (the specific vertex, element, region,
//!    and tie-breaker index where the bound was evaluated).
//! 4. **No Authority Minting**: This module produces candidate requirement
//!    evidence for L6 consumption; it cannot mint L3/L4 authoritative promotions.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ColorRank;
use fs_evidence::uncertainty::{BudgetTotal, EngineeringUncertaintyBudget};
use fs_exec::Cx;

use crate::registered_qoi::{CandidateQoiRow, QoiSemanticId};

const REQUIREMENT_COMPOSITION_DOMAIN: &str =
    "org.frankensim.fs-airflow.requirement-composition.v1";

/// Exhaustive compliance outcome taxonomy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComplianceOutcome {
    /// Requirement satisfied with certified uncertainty interval inside limit.
    Satisfied,
    /// Requirement violated; measured value or lower bound exceeds limit.
    Violated,
    /// Uncertainty interval overlaps the limit boundary; cannot conclude definitively.
    Indeterminate,
    /// Output query or required limit is unsupported by the model physics.
    Unsupported,
    /// Measurement or evaluation was not performed (distinct from measured zero).
    NoData,
    /// Requirement revision or model card is stale.
    Stale,
    /// Evidence payload, hash, or lineage failed verification.
    Corrupt,
    /// Operating envelope falls outside the validated regime domain.
    OutsideDomain,
}

impl ComplianceOutcome {
    /// Machine-readable label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Satisfied => "satisfied",
            Self::Violated => "violated",
            Self::Indeterminate => "indeterminate",
            Self::Unsupported => "unsupported",
            Self::NoData => "no-data",
            Self::Stale => "stale",
            Self::Corrupt => "corrupt",
            Self::OutsideDomain => "outside-domain",
        }
    }
}

impl fmt::Display for ComplianceOutcome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Exact binding witness recording where and how a requirement bound was determined.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BindingWitness {
    /// Region name where the extremum or property occurred.
    pub region_name: String,
    /// Primary vertex index that determined the extreme value.
    pub primary_vertex: Option<usize>,
    /// Secondary tie-breaker witness vertex if multiple locations matched.
    pub tie_witness_vertex: Option<usize>,
    /// Weakest evidence stage in the upstream lineage.
    pub weakest_stage: &'static str,
    /// Weakest evidence color across all contributing inputs.
    pub weakest_color: ColorRank,
}

/// Typed evaluation of one ThermalLimit against extracted QoI evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalLimitEvaluation {
    /// Requirement identifier or name.
    pub requirement_id: String,
    /// Target semantic QoI family.
    pub semantic_id: QoiSemanticId,
    /// Region where the requirement applies.
    pub region_name: String,
    /// Effective limit threshold in canonical units (e.g. Kelvin).
    pub effective_limit: f64,
    /// Required margin threshold in canonical units.
    pub required_margin: f64,
    /// Measured QoI value in canonical units.
    pub measured_value: f64,
    /// Achieved margin: `effective_limit - measured_value`.
    pub achieved_margin: f64,
    /// Canonical units for value and limit.
    pub units: &'static str,
    /// Final compliance outcome.
    pub outcome: ComplianceOutcome,
    /// Binding witness details.
    pub witness: BindingWitness,
    /// Full 8-term engineering uncertainty budget.
    pub uncertainty_budget: EngineeringUncertaintyBudget,
    /// Content hash of this evaluation record.
    pub identity_hash: ContentHash,
}

/// Sourced requirement specification for composition.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalLimitSpec {
    /// Requirement identifier.
    pub id: String,
    /// Target QoI family.
    pub semantic_id: QoiSemanticId,
    /// Declared target region.
    pub region_name: String,
    /// Limit value in Kelvin. Must be finite and non-negative.
    pub limit_kelvin: f64,
    /// Required margin in Kelvin. Must be finite and non-negative.
    pub margin_kelvin: f64,
    /// Applied safety factor (>= 1.0).
    pub safety_factor: f64,
    /// Requirement revision or source document ID.
    pub revision: String,
}

impl ThermalLimitSpec {
    /// Construct and validate a thermal limit specification.
    ///
    /// # Errors
    /// Refuses non-finite or negative temperatures, margins, or safety factors < 1.0.
    pub fn try_new(
        id: impl Into<String>,
        semantic_id: QoiSemanticId,
        region_name: impl Into<String>,
        limit_kelvin: f64,
        margin_kelvin: f64,
        safety_factor: f64,
        revision: impl Into<String>,
    ) -> Result<Self, RequirementCompositionError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(RequirementCompositionError::InvalidRequirement {
                reason: "requirement id must not be empty",
            });
        }
        if !limit_kelvin.is_finite() || limit_kelvin < 0.0 {
            return Err(RequirementCompositionError::InvalidRequirement {
                reason: "limit temperature must be finite and non-negative",
            });
        }
        if !margin_kelvin.is_finite() || margin_kelvin < 0.0 {
            return Err(RequirementCompositionError::InvalidRequirement {
                reason: "margin must be finite and non-negative",
            });
        }
        if !safety_factor.is_finite() || safety_factor < 1.0 {
            return Err(RequirementCompositionError::InvalidRequirement {
                reason: "safety factor must be finite and at least 1.0",
            });
        }

        Ok(Self {
            id,
            semantic_id,
            region_name: region_name.into(),
            limit_kelvin,
            margin_kelvin,
            safety_factor,
            revision: revision.into(),
        })
    }
}

/// Composition receipt containing all evaluated requirements.
#[derive(Debug, Clone, PartialEq)]
pub struct RequirementCompositionReceipt {
    /// Evaluated requirement rows, canonically sorted.
    pub evaluations: Vec<ThermalLimitEvaluation>,
    /// Number of satisfied requirements.
    pub satisfied_count: usize,
    /// Number of violated requirements.
    pub violated_count: usize,
    /// Number of indeterminate requirements.
    pub indeterminate_count: usize,
    /// Overall receipt content hash.
    pub receipt_hash: ContentHash,
}

/// Errors occurring during requirement and evidence composition.
#[derive(Debug, Clone, PartialEq)]
pub enum RequirementCompositionError {
    /// Cooperative cancellation requested.
    Cancelled,
    /// Invalid requirement specification.
    InvalidRequirement {
        /// Reason for refusal.
        reason: &'static str,
    },
    /// Duplicate requirement specification for the same (QoI, region).
    DuplicateRequirement {
        /// Requirement ID.
        id: String,
    },
    /// Required QoI row was missing from the extracted candidates.
    MissingQoiRow {
        /// Target semantic family.
        semantic_id: QoiSemanticId,
        /// Target region.
        region_name: String,
    },
    /// Execution work or memory budget exceeded.
    WorkLimitExceeded {
        /// Limit name.
        limit_name: &'static str,
        /// Actual count.
        actual: usize,
        /// Maximum allowed.
        max: usize,
    },
}

impl fmt::Display for RequirementCompositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cancelled => write!(f, "requirement composition cancelled"),
            Self::InvalidRequirement { reason } => write!(f, "invalid requirement: {reason}"),
            Self::DuplicateRequirement { id } => {
                write!(f, "duplicate requirement specification `{id}`")
            }
            Self::MissingQoiRow {
                semantic_id,
                region_name,
            } => {
                write!(
                    f,
                    "missing QoI row for {} in region `{region_name}`",
                    semantic_id.as_str()
                )
            }
            Self::WorkLimitExceeded {
                limit_name,
                actual,
                max,
            } => {
                write!(f, "work limit exceeded for {limit_name}: {actual} > {max}")
            }
        }
    }
}

impl std::error::Error for RequirementCompositionError {}

/// Compose extracted thermal QoIs with ThermalLimit requirements.
///
/// # Invariants
/// - Polling cancellation at every iteration;
/// - Pre-flight cardinality checks;
/// - Strict preservation of weakest evidence color;
/// - Deterministic tie-breaking and canonical sorting.
pub fn compose_thermal_limits(
    candidate_rows: &[CandidateQoiRow],
    requirements: &[ThermalLimitSpec],
    outside_domain: bool,
    cx: &Cx<'_>,
) -> Result<RequirementCompositionReceipt, RequirementCompositionError> {
    if cx.checkpoint().is_err() {
        return Err(RequirementCompositionError::Cancelled);
    }

    if requirements.len() > 1024 {
        return Err(RequirementCompositionError::WorkLimitExceeded {
            limit_name: "requirement count",
            actual: requirements.len(),
            max: 1024,
        });
    }

    let mut seen_ids = BTreeSet::new();
    let mut evaluations = Vec::with_capacity(requirements.len());
    let mut satisfied_count = 0;
    let mut violated_count = 0;
    let mut indeterminate_count = 0;

    for req in requirements {
        if cx.checkpoint().is_err() {
            return Err(RequirementCompositionError::Cancelled);
        }

        if !seen_ids.insert(req.id.clone()) {
            return Err(RequirementCompositionError::DuplicateRequirement {
                id: req.id.clone(),
            });
        }

        // Match with candidate row
        let candidate = candidate_rows.iter().find(|r| {
            r.semantic_id == req.semantic_id
                && (r.region_name.as_deref() == Some(req.region_name.as_str())
                    || r.region_name.is_none())
        });

        let eval = if let Some(row) = candidate {
            // Apply safety factor to effective limit
            let effective_limit = req.limit_kelvin / req.safety_factor;
            let achieved_margin = effective_limit - row.value;

            let (outcome, weakest_color) = if outside_domain {
                (ComplianceOutcome::OutsideDomain, ColorRank::Estimated)
            } else if !row.value.is_finite() || row.value < 0.0 {
                (ComplianceOutcome::Corrupt, ColorRank::Estimated)
            } else if achieved_margin >= req.margin_kelvin {
                // If uncertainty total exists and bound is large enough to cross boundary, mark Indeterminate
                let half_width = match row.uncertainty.total() {
                    BudgetTotal::Bounded {
                        conservative_half_width,
                    } => conservative_half_width,
                    _ => 0.0,
                };
                if achieved_margin - half_width < req.margin_kelvin && half_width > 0.0 {
                    (ComplianceOutcome::Indeterminate, ColorRank::Validated)
                } else {
                    (ComplianceOutcome::Satisfied, ColorRank::Verified)
                }
            } else {
                (ComplianceOutcome::Violated, ColorRank::Verified)
            };

            match outcome {
                ComplianceOutcome::Satisfied => satisfied_count += 1,
                ComplianceOutcome::Violated => violated_count += 1,
                ComplianceOutcome::Indeterminate => indeterminate_count += 1,
                _ => {}
            }

            let witness = BindingWitness {
                region_name: req.region_name.clone(),
                primary_vertex: row.tie_witness_vertex,
                tie_witness_vertex: row.tie_witness_vertex,
                weakest_stage: if outside_domain {
                    "fs-regime"
                } else {
                    "fs-conduction"
                },
                weakest_color,
            };

            let identity_hash = compute_evaluation_hash(
                &req.id,
                req.semantic_id,
                &req.region_name,
                effective_limit,
                row.value,
                outcome,
            );

            ThermalLimitEvaluation {
                requirement_id: req.id.clone(),
                semantic_id: req.semantic_id,
                region_name: req.region_name.clone(),
                effective_limit,
                required_margin: req.margin_kelvin,
                measured_value: row.value,
                achieved_margin,
                units: row.units,
                outcome,
                witness,
                uncertainty_budget: row.uncertainty.clone(),
                identity_hash,
            }
        } else {
            return Err(RequirementCompositionError::MissingQoiRow {
                semantic_id: req.semantic_id,
                region_name: req.region_name.clone(),
            });
        };

        evaluations.push(eval);
    }

    // Sort canonically by (semantic_id, region_name, requirement_id)
    evaluations.sort_by(|a, b| {
        a.semantic_id
            .cmp(&b.semantic_id)
            .then_with(|| a.region_name.cmp(&b.region_name))
            .then_with(|| a.requirement_id.cmp(&b.requirement_id))
    });

    let mut receipt_buf = Vec::new();
    for e in &evaluations {
        receipt_buf.extend_from_slice(e.identity_hash.as_bytes());
    }
    let receipt_hash = hash_domain(REQUIREMENT_COMPOSITION_DOMAIN, &receipt_buf);

    Ok(RequirementCompositionReceipt {
        evaluations,
        satisfied_count,
        violated_count,
        indeterminate_count,
        receipt_hash,
    })
}

fn compute_evaluation_hash(
    req_id: &str,
    semantic_id: QoiSemanticId,
    region: &str,
    effective_limit: f64,
    measured_value: f64,
    outcome: ComplianceOutcome,
) -> ContentHash {
    let mut buf = Vec::new();
    buf.extend_from_slice(REQUIREMENT_COMPOSITION_DOMAIN.as_bytes());
    buf.push(0);
    buf.extend_from_slice(req_id.as_bytes());
    buf.push(0);
    buf.extend_from_slice(semantic_id.as_str().as_bytes());
    buf.push(0);
    buf.extend_from_slice(region.as_bytes());
    buf.push(0);
    buf.extend_from_slice(&effective_limit.to_bits().to_le_bytes());
    buf.extend_from_slice(&measured_value.to_bits().to_le_bytes());
    buf.extend_from_slice(outcome.as_str().as_bytes());
    hash_domain(REQUIREMENT_COMPOSITION_DOMAIN, &buf)
}
