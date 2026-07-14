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
//!   Perfect Information scores the decision's ambiguity; the campaign's
//!   cancellation-aware action-value reduction places the next sensor on the
//!   candidate whose measurement most sharpens the DECISION (not the
//!   most-uncertain candidate), and says STOP the instant EVPI falls below
//!   threshold — the design choice is already robust.
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

use fs_assimilate::{AssimError, Belief, assimilate_colored_with_shared_poll_quota, point_sensor};
use fs_evidence::{Color, ColorRank, color_leaf_identity_reason};
use fs_exec::Cx;
use fs_toleralloc::{Feature, allocate};
use fs_voi::{Action, ActionKind, ActionValue, DesignEstimate, Recommendation, Uncertainty, evpi};

/// Maximum accepted candidate-name length.
pub const MAX_CANDIDATE_NAME_BYTES: usize = 128;
/// Maximum number of candidates in one synchronous campaign.
pub const MAX_CAMPAIGN_CANDIDATES: usize = 256;
/// Maximum number of sensor placements in one synchronous campaign.
pub const MAX_CAMPAIGN_SENSORS: usize = 4_096;
/// Maximum admitted action-design work units. One sensor-action score evaluates
/// the decision model at every retained normal-quadrature point.
pub const MAX_CAMPAIGN_EVALUATIONS: usize = 10_500_000;

/// Semantic version of the sealed SensorForge report estimator identities.
pub const OED_REPORT_IDENTITY_VERSION: u64 = 5;

const REPORT_ID_DOMAIN: &str = "org.frankensim.fs-oed-e2e.report.v5";
const CAMPAIGN_PLANNING_POLICY_VERSION: u64 = 3;
const CAMPAIGN_POLL_POLICY_VERSION: u64 = 2;
const CAMPAIGN_RECORD_POLL_STRIDE: usize = 256;
const CAMPAIGN_ACTION_POLL_STRIDE: usize = 1;

// Nine-point Gauss-Hermite rule, transformed and normalized for expectations
// under N(0, 1). The policy is deterministic and substantially more faithful
// than evaluating only at the unchanged posterior mean. It remains an
// Estimated decision model: no quadrature-remainder certificate is claimed.
const NORMAL_EXPECTATION_RULE: [(f64, f64); 9] = [
    (-4.512_745_863_399_783_5, 2.234_584_400_774_658_3e-5),
    (-3.205_429_002_856_470_3, 0.002_789_141_321_231_769),
    (-2.076_847_978_677_83, 0.049_916_406_765_217_88),
    (-1.023_255_663_789_132_6, 0.244_097_502_894_939_45),
    (0.0, 0.406_349_206_349_206_35),
    (1.023_255_663_789_132_6, 0.244_097_502_894_939_45),
    (2.076_847_978_677_83, 0.049_916_406_765_217_88),
    (3.205_429_002_856_470_3, 0.002_789_141_321_231_769),
    (4.512_745_863_399_783_5, 2.234_584_400_774_658_3e-5),
];
const ACTION_EVALUATION_FACTOR: usize = NORMAL_EXPECTATION_RULE.len() + 2;

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
    /// Cancellation or poll-quota exhaustion was observed at a deterministic
    /// campaign boundary.
    Cancelled {
        /// Phase whose boundary observed the request.
        phase: &'static str,
        /// Placements committed before the request was observed.
        completed_placements: usize,
        /// Logical work units completed before the request was observed.
        completed_work_units: u128,
        /// Admitted worst-case work bound for the requested campaign cap.
        admitted_work_units: u128,
    },
    /// A lower-layer assimilation observed cancellation after this campaign
    /// had completed the selected sensor observation but before committing a
    /// posterior or placement.
    AssimilationCancelled {
        /// Candidate whose posterior update was cancelled.
        candidate: String,
        /// Placements committed before the lower-layer request was observed.
        completed_placements: usize,
        /// Campaign work completed before entering the cancelled update.
        completed_work_units: u128,
        /// Admitted worst-case campaign work bound.
        admitted_work_units: u128,
        /// Structured lower-layer phase and progress evidence.
        source: AssimError,
    },
    /// Executed logical work did not match the exact realized shape or exceeded
    /// the admitted worst-case bound.
    WorkPlanMismatch {
        /// Work credited by the execution ledger.
        completed_work_units: u128,
        /// Exact work implied by the realized early-stop path.
        realized_work_units: u128,
        /// Worst-case work admitted before scientific execution.
        admitted_work_units: u128,
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
    /// The bounded VoI reduction returned an action outside its own menu.
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
                "campaign work {candidates}^2 x {max_sensors} x \
                 {ACTION_EVALUATION_FACTOR} = {evaluations} exceeds \
                {max_evaluations} action-design evaluations"
            ),
            Self::Cancelled {
                phase,
                completed_placements,
                completed_work_units,
                admitted_work_units,
            } => write!(
                f,
                "campaign cancelled or poll budget exhausted during {phase} after \
                 {completed_placements} placements and \
                {completed_work_units}/{admitted_work_units} admitted logical work units"
            ),
            Self::AssimilationCancelled {
                candidate,
                completed_placements,
                completed_work_units,
                admitted_work_units,
                source,
            } => write!(
                f,
                "assimilation for candidate `{candidate}` cancelled after \
                 {completed_placements} committed placements and \
                 {completed_work_units}/{admitted_work_units} admitted campaign work units: {source}"
            ),
            Self::WorkPlanMismatch {
                completed_work_units,
                realized_work_units,
                admitted_work_units,
            } => write!(
                f,
                "campaign work ledger mismatch: completed {completed_work_units}, \
                 realized {realized_work_units}, admitted {admitted_work_units}"
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
                write!(f, "VoI reduction returned unknown action `{action}`")
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
            Self::BeliefInvariant(source)
            | Self::Assimilation { source, .. }
            | Self::AssimilationCancelled { source, .. } => Some(source),
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
///
/// Reports are sealed outputs of [`run_campaign`]. Read-only accessors expose
/// the complete result without permitting callers to replace science fields or
/// evidence identities independently.
#[derive(Debug, Clone, PartialEq)]
pub struct OedReport {
    /// Candidate names in the order sensors were placed.
    placements: Vec<String>,
    /// Number of sensors placed.
    sensors_placed: usize,
    /// Total prior variance across candidates.
    prior_total_variance: f64,
    /// Total posterior variance across candidates.
    posterior_total_variance: f64,
    /// Fractional variance reduction.
    variance_reduction: f64,
    /// EVPI before any sensor.
    initial_evpi: f64,
    /// EVPI after the campaign stopped.
    final_evpi: f64,
    /// Did the decision become robust (planner chose to STOP)?
    decision_robust: bool,
    /// The finally-chosen (lowest-cost posterior) design.
    chosen_design: String,
    /// The cost-optimal tolerance allocation `(name, tolerance)`.
    /// A zero-sensitivity candidate receives `+infinity`, the exact unconstrained
    /// optimum under the first-order allocation model.
    allocation: Vec<(String, f64)>,
    /// EVPI before sensing and after every completed placement.
    evpi_trace: Vec<f64>,
    /// Final scalar posterior in candidate order.
    posteriors: Vec<PosteriorSummary>,
    /// Instrument-bound estimated candidate emitted by each assimilation.
    assimilation_colors: Vec<Color>,
    /// The posterior-variance color (`Estimated` until independently certified).
    variance_color: Color,
    /// The EVPI color (`Estimated` — decision-theoretic).
    evpi_color: Color,
}

impl OedReport {
    /// Candidate names in placement order.
    #[must_use]
    pub fn placements(&self) -> &[String] {
        &self.placements
    }

    /// Number of completed sensor placements.
    #[must_use]
    pub const fn sensors_placed(&self) -> usize {
        self.sensors_placed
    }

    /// Total variance before sensing.
    #[must_use]
    pub const fn prior_total_variance(&self) -> f64 {
        self.prior_total_variance
    }

    /// Total variance after sensing.
    #[must_use]
    pub const fn posterior_total_variance(&self) -> f64 {
        self.posterior_total_variance
    }

    /// Fractional reduction in total variance.
    #[must_use]
    pub const fn variance_reduction(&self) -> f64 {
        self.variance_reduction
    }

    /// EVPI before the first placement.
    #[must_use]
    pub const fn initial_evpi(&self) -> f64 {
        self.initial_evpi
    }

    /// EVPI when the campaign stopped.
    #[must_use]
    pub const fn final_evpi(&self) -> f64 {
        self.final_evpi
    }

    /// Whether the modeled EVPI met the requested stop threshold.
    #[must_use]
    pub const fn decision_robust(&self) -> bool {
        self.decision_robust
    }

    /// Finally chosen design identity.
    #[must_use]
    pub fn chosen_design(&self) -> &str {
        &self.chosen_design
    }

    /// Cost-optimal tolerance allocation in candidate order.
    #[must_use]
    pub fn allocation(&self) -> &[(String, f64)] {
        &self.allocation
    }

    /// EVPI before sensing and after each completed placement.
    #[must_use]
    pub fn evpi_trace(&self) -> &[f64] {
        &self.evpi_trace
    }

    /// Final scalar posterior summaries in candidate order.
    #[must_use]
    pub fn posteriors(&self) -> &[PosteriorSummary] {
        &self.posteriors
    }

    /// Instrument-bound colors emitted by completed assimilations.
    #[must_use]
    pub fn assimilation_colors(&self) -> &[Color] {
        &self.assimilation_colors
    }

    /// Sealed posterior-variance evidence color.
    #[must_use]
    pub const fn variance_color(&self) -> &Color {
        &self.variance_color
    }

    /// Sealed EVPI evidence color.
    #[must_use]
    pub const fn evpi_color(&self) -> &Color {
        &self.evpi_color
    }
}

/// The preflighted worst-case logical work bound. A unit is one bounded
/// candidate/color record visit, one scalar assimilation transaction, or one
/// retained hash; it is a deterministic scheduling/accounting unit, not an
/// instruction count. Early STOP paths realize fewer units and are checked
/// separately rather than padded with phantom work.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CampaignWorkPlan {
    candidates: usize,
    max_sensors: usize,
    action_design_evaluations: usize,
    setup_work_units: u128,
    per_placement_work_units: u128,
    maximum_finalization_work_units: u128,
    admitted_work_units: u128,
}

/// Exact logical work implied by one realized early-stop path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CampaignRealizedWorkPlan {
    completed_placements: u128,
    action_rounds: u128,
    positive_prior_candidates: u128,
    setup_work_units: u128,
    placement_work_units: u128,
    incomplete_action_work_units: u128,
    finalization_work_units: u128,
    realized_work_units: u128,
}

impl CampaignRealizedWorkPlan {
    const fn identity_fields(self) -> [u128; 8] {
        [
            self.completed_placements,
            self.action_rounds,
            self.positive_prior_candidates,
            self.setup_work_units,
            self.placement_work_units,
            self.incomplete_action_work_units,
            self.finalization_work_units,
            self.realized_work_units,
        ]
    }
}

impl CampaignWorkPlan {
    fn checked(candidates: usize, max_sensors: usize) -> Result<Self, OedError> {
        let action_design_pairs =
            candidates
                .checked_mul(candidates)
                .ok_or(OedError::WorkBudgetExceeded {
                    candidates,
                    max_sensors,
                    evaluations: usize::MAX,
                    max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
                })?;
        let action_design_evaluations = action_design_pairs
            .checked_mul(max_sensors)
            .and_then(|work| work.checked_mul(ACTION_EVALUATION_FACTOR))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        if action_design_evaluations > MAX_CAMPAIGN_EVALUATIONS {
            return Err(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: action_design_evaluations,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            });
        }

        let n = candidates as u128;
        let placements = max_sensors as u128;
        // Setup: validation, belief construction, prior variance, and initial
        // estimates/EVPI (five candidate scans). Sensor actions are rebuilt
        // after every posterior update because their effect depends on P and R.
        let setup_work_units = n.checked_mul(5).ok_or(OedError::WorkBudgetExceeded {
            candidates,
            max_sensors,
            evaluations: usize::MAX,
            max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
        })?;
        // Each placement: one action construction scan; for every action, a
        // target lookup, posterior-template construction, and one EVPI scan per
        // normal-expectation node; then the action-menu record, chosen-action
        // lookup, sensor+assimilation transaction, and refreshed estimates/EVPI.
        let per_placement_work_units = n
            .checked_mul(n)
            .and_then(|work| work.checked_mul(ACTION_EVALUATION_FACTOR as u128))
            .and_then(|work| work.checked_add(n.checked_mul(5)?))
            .and_then(|work| work.checked_add(2))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        // Worst-case final science consumes twelve candidate scans. Each full
        // report identity reserves three candidate-sized and three
        // max-placement-sized sequences, the trace's initial value, and one
        // bounded hash. The realized plan later substitutes the actual positive
        // priors and placements. The last unit is publication.
        let maximum_finalization_work_units = n
            .checked_mul(18)
            .and_then(|work| work.checked_add(placements.checked_mul(6)?))
            .and_then(|work| work.checked_add(5))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;
        let admitted_work_units = per_placement_work_units
            .checked_mul(placements)
            .and_then(|work| work.checked_add(setup_work_units))
            .and_then(|work| work.checked_add(maximum_finalization_work_units))
            .ok_or(OedError::WorkBudgetExceeded {
                candidates,
                max_sensors,
                evaluations: usize::MAX,
                max_evaluations: MAX_CAMPAIGN_EVALUATIONS,
            })?;

        Ok(Self {
            candidates,
            max_sensors,
            action_design_evaluations,
            setup_work_units,
            per_placement_work_units,
            maximum_finalization_work_units,
            admitted_work_units,
        })
    }

    fn realized(
        self,
        completed_placements: usize,
        action_rounds: usize,
        positive_prior_candidates: usize,
    ) -> Option<CampaignRealizedWorkPlan> {
        if completed_placements > self.max_sensors
            || action_rounds < completed_placements
            || action_rounds.checked_sub(completed_placements)? > 1
            || positive_prior_candidates > self.candidates
        {
            return None;
        }
        let n = self.candidates as u128;
        let completed_placements = completed_placements as u128;
        let action_rounds = action_rounds as u128;
        let positive_prior_candidates = positive_prior_candidates as u128;
        let placement_work_units = self
            .per_placement_work_units
            .checked_mul(completed_placements)?;
        let incomplete_rounds = action_rounds.checked_sub(completed_placements)?;
        let incomplete_action_work_units = n
            .checked_mul(n)?
            .checked_mul(ACTION_EVALUATION_FACTOR as u128)?
            .checked_add(n.checked_mul(2)?)?
            .checked_mul(incomplete_rounds)?;
        // Final science is 7n + 5m, where m is the number of positive-prior
        // candidates accepted by fs-toleralloc. Two report identities add
        // 6n + 6s + 4, and publication adds one.
        let finalization_work_units = n
            .checked_mul(13)?
            .checked_add(positive_prior_candidates.checked_mul(5)?)?
            .checked_add(completed_placements.checked_mul(6)?)?
            .checked_add(5)?;
        let realized_work_units = self
            .setup_work_units
            .checked_add(placement_work_units)?
            .checked_add(incomplete_action_work_units)?
            .checked_add(finalization_work_units)?;
        if realized_work_units > self.admitted_work_units {
            return None;
        }
        Some(CampaignRealizedWorkPlan {
            completed_placements,
            action_rounds,
            positive_prior_candidates,
            setup_work_units: self.setup_work_units,
            placement_work_units,
            incomplete_action_work_units,
            finalization_work_units,
            realized_work_units,
        })
    }
}

#[derive(Debug)]
struct CampaignProgress {
    completed_placements: usize,
    completed_work_units: u128,
    // Private invocation-global ledger: only `checkpoint` and the nested
    // assimilation transaction can decrease this value, and no caller can
    // replace it between campaign phases.
    polls_remaining: u32,
}

impl CampaignProgress {
    fn new(cx: &Cx<'_>, completed_work_units: u128) -> Self {
        Self {
            completed_placements: 0,
            completed_work_units,
            polls_remaining: cx.budget().poll_quota,
        }
    }

    fn advance(&mut self, units: u128) {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(units)
            .expect("admitted campaign progress cannot exceed u128");
    }

    fn checkpoint(
        &mut self,
        cx: &Cx<'_>,
        plan: CampaignWorkPlan,
        phase: &'static str,
    ) -> Result<(), OedError> {
        if self.polls_remaining == 0 {
            return Err(self.cancelled(plan, phase));
        }
        if self.polls_remaining != u32::MAX {
            self.polls_remaining -= 1;
        }
        cx.checkpoint().map_err(|_| self.cancelled(plan, phase))
    }

    fn finish(
        &self,
        plan: CampaignWorkPlan,
        realized: CampaignRealizedWorkPlan,
    ) -> Result<(), OedError> {
        if self.completed_work_units == realized.realized_work_units
            && realized.realized_work_units <= plan.admitted_work_units
        {
            Ok(())
        } else {
            Err(OedError::WorkPlanMismatch {
                completed_work_units: self.completed_work_units,
                realized_work_units: realized.realized_work_units,
                admitted_work_units: plan.admitted_work_units,
            })
        }
    }

    fn cancelled(&self, plan: CampaignWorkPlan, phase: &'static str) -> OedError {
        OedError::Cancelled {
            phase,
            completed_placements: self.completed_placements,
            completed_work_units: self.completed_work_units,
            admitted_work_units: plan.admitted_work_units,
        }
    }
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
                    numerical: 0.0,
                    statistical: variance.sqrt(),
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
    // `fs_voi` preserves input order as its final equal-mean tie-break.
    // Candidate declarations are an unordered menu here, so impose the
    // campaign's stable identity order before calling the current top-two
    // approximation. The full multi-alternative EVPI replacement remains a
    // separate scientific upgrade, but equal means must not make this policy
    // depend on caller vector order in the meantime.
    let mut canonical = estimates.to_vec();
    canonical.sort_by(|left, right| left.name.cmp(&right.name));
    let value = evpi(&canonical);
    if value.is_finite() && value >= 0.0 {
        Ok(canonicalize_zero(value))
    } else {
        Err(OedError::NonFiniteComputation { quantity: "EVPI" })
    }
}

/// Stable scalar Kalman variance update `P' = P R / (P + R)` without the
/// overflowing intermediate `P * R`. Candidate construction and Belief enforce
/// the input domains; this check protects the independently callable planner
/// path from a derived floating-point failure.
fn predicted_posterior_variance(prior: f64, noise: f64) -> Result<f64, OedError> {
    if prior == 0.0 {
        return Ok(0.0);
    }
    if !prior.is_finite() || prior < 0.0 || !noise.is_finite() || noise <= 0.0 {
        return Err(OedError::NonFiniteComputation {
            quantity: "sensor posterior variance inputs",
        });
    }
    // Divide the smaller operand only by a value at least as large. This is
    // algebraically identical to `P R / (P + R)`, avoids the overflowing
    // product, and—unlike scaling both inputs by `max(P, R)`—does not erase a
    // representable subnormal posterior when P and R span the full exponent
    // range.
    let posterior = if prior <= noise {
        prior / (1.0 + prior / noise)
    } else {
        noise / (1.0 + noise / prior)
    };
    if !posterior.is_finite() || posterior < 0.0 || posterior > prior {
        return Err(OedError::NonFiniteComputation {
            quantity: "predicted sensor posterior variance",
        });
    }
    Ok(canonicalize_zero(posterior))
}

fn sensor_actions(candidates: &[Candidate], beliefs: &[Belief]) -> Result<Vec<Action>, OedError> {
    if candidates.len() != beliefs.len() {
        return Err(OedError::NonFiniteComputation {
            quantity: "candidate/belief cardinality",
        });
    }
    candidates
        .iter()
        .zip(beliefs)
        .map(|(candidate, belief)| {
            let prior = belief.variance(0).map_err(OedError::BeliefInvariant)?;
            let posterior = predicted_posterior_variance(prior, candidate.sensor_noise)?;
            let reduction = if prior == 0.0 {
                0.0
            } else {
                let std_ratio = (posterior / prior).clamp(0.0, 1.0).sqrt();
                canonicalize_zero((1.0 - std_ratio).clamp(0.0, 1.0))
            };
            Ok(Action {
                name: format!("measure-{}", candidate.name),
                kind: ActionKind::Sample,
                target_design: candidate.name.clone(),
                reduction,
                cost: candidate.sensor_cost,
            })
        })
        .collect()
}

/// Outcome-integrated value of a scalar Gaussian sensor action. The posterior
/// variance is the exact declared Kalman-model update. Posterior-mean movement
/// is integrated under its pre-posterior Gaussian distribution with the fixed
/// rule above, so a noisy sensor cannot inherit a fictitious universal effect.
fn expected_sensor_action_value(
    estimates: &[DesignEstimate],
    action: &Action,
    before: f64,
) -> Result<ActionValue, OedError> {
    let target = estimates
        .iter()
        .position(|estimate| estimate.name == action.target_design)
        .ok_or_else(|| OedError::UnknownRecommendation {
            action: action.name.clone(),
        })?;
    let prior_mean = estimates[target].mean;
    if action.kind != ActionKind::Sample {
        return Err(OedError::NonFiniteComputation {
            quantity: "non-sensor action in sensor planner",
        });
    }
    let prior_std = estimates[target].uncertainty.total_std();
    let mut posterior_estimates = estimates.to_vec();
    posterior_estimates[target].uncertainty.statistical *= (1.0 - action.reduction).clamp(0.0, 1.0);
    let posterior_std = posterior_estimates[target].uncertainty.total_std();
    let mean_shift_variance = (prior_std * prior_std - posterior_std * posterior_std).max(0.0);
    let mean_shift_std = mean_shift_variance.sqrt();

    let mut expected_remaining_evpi = 0.0;
    for (normal_node, probability_weight) in NORMAL_EXPECTATION_RULE {
        let posterior_mean = prior_mean + normal_node * mean_shift_std;
        if !posterior_mean.is_finite() {
            return Err(OedError::NonFiniteComputation {
                quantity: "predictive posterior mean",
            });
        }
        posterior_estimates[target].mean = canonicalize_zero(posterior_mean);
        expected_remaining_evpi += probability_weight * checked_evpi(&posterior_estimates)?;
    }
    // Preserve the exact identity map after executing the declared fixed-shape
    // quadrature work. Summing nine identical weighted EVPI values can land a
    // single ulp below `before`; that rounding artifact is not sensor value.
    if action.reduction == 0.0 {
        expected_remaining_evpi = before;
    }
    if !expected_remaining_evpi.is_finite() || expected_remaining_evpi < 0.0 {
        return Err(OedError::NonFiniteComputation {
            quantity: "expected posterior EVPI",
        });
    }
    let value = canonicalize_zero((before - expected_remaining_evpi).max(0.0));
    let value_per_cost = if action.cost.is_finite() && action.cost > 0.0 {
        value / action.cost
    } else {
        0.0
    };
    if !value.is_finite() || !value_per_cost.is_finite() {
        return Err(OedError::NonFiniteComputation {
            quantity: "sensor action value",
        });
    }
    Ok(ActionValue {
        action: action.name.clone(),
        value,
        cost: action.cost,
        value_per_cost,
    })
}

fn precision_allocation(candidates: &[Candidate]) -> Result<(Vec<(String, f64)>, usize), OedError> {
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
    let positive_prior_candidates = features.len();
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

    let allocation = candidates
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
        .collect::<Result<Vec<_>, OedError>>()?;
    Ok((allocation, positive_prior_candidates))
}

fn push_bytes(output: &mut Vec<u8>, value: &[u8]) {
    output.extend_from_slice(&(value.len() as u64).to_le_bytes());
    output.extend_from_slice(value);
}

fn push_str(output: &mut Vec<u8>, value: &str) {
    push_bytes(output, value.as_bytes());
}

#[derive(Debug, Clone, Copy)]
struct ReportIdentityOutputs<'a> {
    placements: &'a [String],
    sensors_placed: usize,
    prior_total_variance: f64,
    posterior_total_variance: f64,
    variance_reduction: f64,
    initial_evpi: f64,
    final_evpi: f64,
    decision_robust: bool,
    chosen_design: &'a str,
    allocation: &'a [(String, f64)],
    evpi_trace: &'a [f64],
    posteriors: &'a [PosteriorSummary],
    assimilation_colors: &'a [Color],
    variance_color_dispersion: f64,
    evpi_color_dispersion: f64,
}

#[derive(Debug, Clone, Copy)]
struct ReportIdentitySource<'a> {
    candidates: &'a [Candidate],
    threshold: f64,
    max_sensors: usize,
    outputs: ReportIdentityOutputs<'a>,
    plan: CampaignWorkPlan,
    realized: CampaignRealizedWorkPlan,
}

fn push_estimated_color_descriptor(output: &mut Vec<u8>, quantity: &str, dispersion: f64) {
    push_str(output, quantity);
    output.extend_from_slice(&fs_evidence::COLOR_ALGEBRA_VERSION.to_le_bytes());
    push_str(output, "Estimated");
    output.extend_from_slice(&dispersion.to_bits().to_le_bytes());
}

fn report_identity(
    quantity: &str,
    source: &ReportIdentitySource<'_>,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<String, OedError> {
    report_identity_with_rule(quantity, source, &NORMAL_EXPECTATION_RULE, progress, cx)
}

#[allow(clippy::too_many_lines)] // One canonical manifest keeps field order auditable.
fn report_identity_with_rule(
    quantity: &str,
    source: &ReportIdentitySource<'_>,
    expectation_rule: &[(f64, f64)],
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<String, OedError> {
    let candidates = source.candidates;
    let threshold = source.threshold;
    let max_sensors = source.max_sensors;
    let outputs = source.outputs;
    let plan = source.plan;
    let realized = source.realized;
    let mut canonical = Vec::new();
    canonical.extend_from_slice(&OED_REPORT_IDENTITY_VERSION.to_le_bytes());
    push_str(&mut canonical, cx.mode().name());
    let stream = cx.stream_key();
    for value in [stream.seed, stream.kernel_id, stream.tile, stream.iteration] {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    let budget = cx.budget();
    match budget.deadline {
        Some(deadline) => {
            canonical.push(1);
            canonical.extend_from_slice(&deadline.as_nanos().to_le_bytes());
        }
        None => canonical.push(0),
    }
    canonical.extend_from_slice(&budget.poll_quota.to_le_bytes());
    match budget.cost_quota {
        Some(cost_quota) => {
            canonical.push(1);
            canonical.extend_from_slice(&cost_quota.to_le_bytes());
        }
        None => canonical.push(0),
    }
    canonical.push(budget.priority);
    for value in [
        plan.candidates as u128,
        plan.max_sensors as u128,
        plan.action_design_evaluations as u128,
        plan.setup_work_units,
        plan.per_placement_work_units,
        plan.maximum_finalization_work_units,
        plan.admitted_work_units,
    ] {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    for value in realized.identity_fields() {
        canonical.extend_from_slice(&value.to_le_bytes());
    }
    canonical.extend_from_slice(&CAMPAIGN_PLANNING_POLICY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&CAMPAIGN_POLL_POLICY_VERSION.to_le_bytes());
    canonical.extend_from_slice(&fs_voi::EVPI_SEMANTICS_VERSION.to_le_bytes());
    canonical.extend_from_slice(&(expectation_rule.len() as u64).to_le_bytes());
    for &(node, weight) in expectation_rule {
        canonical.extend_from_slice(&node.to_bits().to_le_bytes());
        canonical.extend_from_slice(&weight.to_bits().to_le_bytes());
    }
    canonical.extend_from_slice(&(CAMPAIGN_RECORD_POLL_STRIDE as u64).to_le_bytes());
    canonical.extend_from_slice(&(CAMPAIGN_ACTION_POLL_STRIDE as u64).to_le_bytes());
    push_str(&mut canonical, quantity);
    canonical.extend_from_slice(&(candidates.len() as u64).to_le_bytes());
    for (index, candidate) in candidates.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity candidates")?;
        }
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
        progress.advance(1);
    }
    canonical.extend_from_slice(&threshold.to_bits().to_le_bytes());
    canonical.extend_from_slice(&(max_sensors as u64).to_le_bytes());
    canonical.extend_from_slice(&(outputs.placements.len() as u64).to_le_bytes());
    for (index, placement) in outputs.placements.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity placements")?;
        }
        push_str(&mut canonical, placement);
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.sensors_placed as u64).to_le_bytes());
    for value in [
        outputs.prior_total_variance,
        outputs.posterior_total_variance,
        outputs.variance_reduction,
        outputs.initial_evpi,
        outputs.final_evpi,
    ] {
        canonical.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    canonical.push(u8::from(outputs.decision_robust));
    push_str(&mut canonical, outputs.chosen_design);
    canonical.extend_from_slice(&(outputs.allocation.len() as u64).to_le_bytes());
    for (index, (name, tolerance)) in outputs.allocation.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity allocation")?;
        }
        push_str(&mut canonical, name);
        canonical.extend_from_slice(&tolerance.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.evpi_trace.len() as u64).to_le_bytes());
    for (index, value) in outputs.evpi_trace.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity EVPI trace")?;
        }
        canonical.extend_from_slice(&value.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.posteriors.len() as u64).to_le_bytes());
    for (index, posterior) in outputs.posteriors.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity posteriors")?;
        }
        push_str(&mut canonical, &posterior.name);
        canonical.extend_from_slice(&posterior.mean.to_bits().to_le_bytes());
        canonical.extend_from_slice(&posterior.variance.to_bits().to_le_bytes());
        progress.advance(1);
    }
    canonical.extend_from_slice(&(outputs.assimilation_colors.len() as u64).to_le_bytes());
    for (index, color) in outputs.assimilation_colors.iter().enumerate() {
        if index.is_multiple_of(CAMPAIGN_RECORD_POLL_STRIDE) {
            progress.checkpoint(cx, plan, "report identity assimilation colors")?;
        }
        push_bytes(&mut canonical, &color.canonical_bytes());
        progress.advance(1);
    }
    // The estimator strings are the hashes being derived and therefore cannot
    // appear in their own preimage. Bind the stable color algebra, variant, and
    // both dispersions; the sealed report is constructed only after both
    // estimator strings have been derived from this complete source.
    push_estimated_color_descriptor(
        &mut canonical,
        "posterior-variance",
        outputs.variance_color_dispersion,
    );
    push_estimated_color_descriptor(&mut canonical, "evpi", outputs.evpi_color_dispersion);
    progress.checkpoint(cx, plan, "report identity hash")?;
    let identity = format!(
        "sensorforge-{quantity}:v{OED_REPORT_IDENTITY_VERSION}:{}",
        fs_blake3::hash_domain(REPORT_ID_DOMAIN, &canonical)
    );
    progress.advance(1);
    debug_assert!(color_leaf_identity_reason(&identity).is_none());
    Ok(identity)
}

fn validate_campaign(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
) -> Result<(f64, CampaignWorkPlan), OedError> {
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
    let plan = CampaignWorkPlan::checked(candidates.len(), max_sensors)?;
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
    Ok((canonicalize_zero(threshold), plan))
}

struct CampaignState {
    beliefs: Vec<Belief>,
    placements: Vec<String>,
    assimilation_colors: Vec<Color>,
    evpi_trace: Vec<f64>,
    decision_robust: bool,
    action_rounds: usize,
}

fn recommend_with_cancellation(
    estimates: &[DesignEstimate],
    actions: &[Action],
    current_evpi: f64,
    threshold: f64,
    plan: CampaignWorkPlan,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<Recommendation, OedError> {
    if current_evpi <= threshold {
        return Ok(Recommendation::Stop {
            reason: format!("decision robust: EVPI {current_evpi:.3e} <= {threshold:.3e}"),
        });
    }

    let mut best = None;
    for action in actions {
        progress.checkpoint(cx, plan, "action-value tile")?;
        let value = expected_sensor_action_value(estimates, action, current_evpi)?;
        progress.advance((estimates.len() as u128) * (ACTION_EVALUATION_FACTOR as u128) + 1);
        if value.value <= 0.0 || value.value_per_cost <= 0.0 {
            continue;
        }
        let replace = best.as_ref().is_none_or(|current: &ActionValue| {
            match value.value_per_cost.total_cmp(&current.value_per_cost) {
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Equal => value.action < current.action,
                std::cmp::Ordering::Less => false,
            }
        });
        if replace {
            best = Some(value);
        }
    }
    progress.checkpoint(cx, plan, "action-value drain")?;

    Ok(match best {
        Some(value) => Recommendation::Act {
            action: value.action,
            value_per_cost: value.value_per_cost,
        },
        None => Recommendation::Stop {
            reason: "no action changes the decision".to_string(),
        },
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn execute_placements(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    mut beliefs: Vec<Belief>,
    mut estimates: Vec<DesignEstimate>,
    initial_evpi: f64,
    plan: CampaignWorkPlan,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<CampaignState, OedError> {
    let mut placements = Vec::new();
    let mut assimilation_colors = Vec::new();
    let mut evpi_trace = vec![initial_evpi];
    let mut decision_robust = false;
    let mut action_rounds = 0;
    let mut current_evpi = initial_evpi;

    loop {
        if current_evpi <= threshold {
            decision_robust = true;
            break;
        }
        if placements.len() >= max_sensors {
            break;
        }
        progress.checkpoint(cx, plan, "action construction")?;
        let actions = sensor_actions(candidates, &beliefs)?;
        progress.advance(candidates.len() as u128);
        progress.checkpoint(cx, plan, "action construction drain")?;
        let recommendation = recommend_with_cancellation(
            &estimates,
            &actions,
            current_evpi,
            threshold,
            plan,
            progress,
            cx,
        )?;
        action_rounds += 1;
        let Recommendation::Act { action, .. } = recommendation else {
            break;
        };
        progress.checkpoint(cx, plan, "chosen-action lookup")?;
        let idx = actions
            .iter()
            .position(|candidate| candidate.name == action)
            .ok_or(OedError::UnknownRecommendation { action })?;
        progress.advance(candidates.len() as u128);
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
        progress.advance(1);
        let next_count = placements.len() + 1;
        let posterior = assimilate_colored_with_shared_poll_quota(
            &beliefs[idx],
            std::slice::from_ref(&observation),
            "sensor_count",
            0.0,
            next_count as f64,
            cx,
            &mut progress.polls_remaining,
        )
        .map_err(|source| {
            if matches!(source, AssimError::Cancelled { .. }) {
                OedError::AssimilationCancelled {
                    candidate: candidates[idx].name.clone(),
                    completed_placements: progress.completed_placements,
                    completed_work_units: progress.completed_work_units,
                    admitted_work_units: plan.admitted_work_units,
                    source,
                }
            } else {
                OedError::Assimilation {
                    candidate: candidates[idx].name.clone(),
                    source,
                }
            }
        })?;
        progress.advance(1);
        // Request -> drain -> finalize: do not publish the scratch posterior
        // into campaign state until the lower-layer transaction has drained
        // and this deterministic commit boundary is still live.
        progress.checkpoint(cx, plan, "placement commit")?;
        beliefs[idx] = posterior.belief().clone();
        assimilation_colors.push(posterior.color().clone());
        placements.push(candidates[idx].name.clone());
        progress.completed_placements = placements.len();

        progress.checkpoint(cx, plan, "posterior estimate refresh")?;
        estimates = to_estimates(candidates, &beliefs)?;
        progress.advance(candidates.len() as u128);
        progress.checkpoint(cx, plan, "posterior EVPI refresh")?;
        current_evpi = checked_evpi(&estimates)?;
        progress.advance(candidates.len() as u128);
        evpi_trace.push(current_evpi);
    }

    Ok(CampaignState {
        beliefs,
        placements,
        assimilation_colors,
        evpi_trace,
        decision_robust,
        action_rounds,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn finish_report(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    prior_total_variance: f64,
    initial_evpi: f64,
    state: CampaignState,
    plan: CampaignWorkPlan,
    progress: &mut CampaignProgress,
    cx: &Cx<'_>,
) -> Result<OedReport, OedError> {
    progress.checkpoint(cx, plan, "final estimate summary")?;
    let estimates = to_estimates(candidates, &state.beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "final EVPI")?;
    let final_evpi = checked_evpi(&estimates)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "final variance")?;
    let posterior_total_variance = total_variance(&state.beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "chosen-design reduction")?;
    let chosen_design = estimates
        .iter()
        .min_by(|a, b| a.mean.total_cmp(&b.mean).then_with(|| a.name.cmp(&b.name)))
        .map(|design| design.name.clone())
        .ok_or(OedError::NonFiniteComputation {
            quantity: "chosen design",
        })?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "precision allocation")?;
    let (allocation, positive_prior_candidates) = precision_allocation(candidates)?;
    progress.advance((candidates.len() as u128) * 2 + (positive_prior_candidates as u128) * 5);
    progress.checkpoint(cx, plan, "posterior summaries")?;
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
    progress.advance(candidates.len() as u128);
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
    let sensors_placed = state.placements.len();
    let realized = plan
        .realized(
            sensors_placed,
            state.action_rounds,
            positive_prior_candidates,
        )
        .ok_or(OedError::WorkPlanMismatch {
            completed_work_units: progress.completed_work_units,
            realized_work_units: u128::MAX,
            admitted_work_units: plan.admitted_work_units,
        })?;
    let variance_color_dispersion = f64::INFINITY;
    let evpi_color_dispersion = final_evpi;
    let (variance_identity, evpi_identity) = {
        let source = ReportIdentitySource {
            candidates,
            threshold,
            max_sensors,
            outputs: ReportIdentityOutputs {
                placements: &state.placements,
                sensors_placed,
                prior_total_variance,
                posterior_total_variance,
                variance_reduction,
                initial_evpi,
                final_evpi,
                decision_robust: state.decision_robust,
                chosen_design: &chosen_design,
                allocation: &allocation,
                evpi_trace: &state.evpi_trace,
                posteriors: &posteriors,
                assimilation_colors: &state.assimilation_colors,
                variance_color_dispersion,
                evpi_color_dispersion,
            },
            plan,
            realized,
        };
        let variance_identity = report_identity("posterior-variance", &source, progress, cx)?;
        let evpi_identity = report_identity("evpi", &source, progress, cx)?;
        (variance_identity, evpi_identity)
    };

    progress.checkpoint(cx, plan, "report publication")?;
    progress.advance(1);
    progress.finish(plan, realized)?;

    Ok(OedReport {
        sensors_placed,
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
            dispersion: variance_color_dispersion,
        },
        evpi_color: Color::Estimated {
            estimator: evpi_identity,
            dispersion: evpi_color_dispersion,
        },
    })
}

/// Run the SensorForge campaign under an explicit execution context; stop when
/// EVPI <= `threshold` or after `max_sensors` placements.
///
/// The complete worst-case work bound is checked before scientific work starts,
/// and the exact realized early-stop shape is checked before publication. The
/// initial STOP condition is evaluated even when `max_sensors == 0`, and
/// cancellation is polled at deterministic action/record boundaries.
///
/// # Errors
/// Returns [`OedError`] for invalid campaign bounds, duplicate candidate names,
/// observed cancellation, a lower-layer assimilation/allocation failure, or a
/// non-finite derived value. A cancellation never returns a partial report.
pub fn run_campaign(
    candidates: &[Candidate],
    threshold: f64,
    max_sensors: usize,
    cx: &Cx<'_>,
) -> Result<OedReport, OedError> {
    let (threshold, plan) = validate_campaign(candidates, threshold, max_sensors)?;
    let mut progress = CampaignProgress::new(cx, candidates.len() as u128);
    progress.checkpoint(cx, plan, "campaign admission")?;
    let beliefs: Vec<Belief> = candidates
        .iter()
        .map(|c| Belief::scalar(c.prior_mean, c.prior_var))
        .collect::<Result<Vec<_>, _>>()
        .map_err(OedError::BeliefInvariant)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "prior variance")?;
    let prior_total_variance = total_variance(&beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "initial estimates")?;
    let estimates = to_estimates(candidates, &beliefs)?;
    progress.advance(candidates.len() as u128);
    progress.checkpoint(cx, plan, "initial EVPI")?;
    let initial_evpi = checked_evpi(&estimates)?;
    progress.advance(candidates.len() as u128);
    let state = execute_placements(
        candidates,
        threshold,
        max_sensors,
        beliefs,
        estimates,
        initial_evpi,
        plan,
        &mut progress,
        cx,
    )?;
    finish_report(
        candidates,
        threshold,
        max_sensors,
        prior_total_variance,
        initial_evpi,
        state,
        plan,
        &mut progress,
        cx,
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

#[cfg(test)]
mod tests {
    use super::{
        CAMPAIGN_PLANNING_POLICY_VERSION, CAMPAIGN_POLL_POLICY_VERSION, CampaignProgress,
        CampaignRealizedWorkPlan, CampaignWorkPlan, Candidate, NORMAL_EXPECTATION_RULE,
        OED_REPORT_IDENTITY_VERSION, PosteriorSummary, REPORT_ID_DOMAIN, ReportIdentityOutputs,
        ReportIdentitySource, canonicalize_zero, demo_candidates, predicted_posterior_variance,
        report_identity_with_rule,
    };
    use fs_evidence::Color;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

    #[derive(Clone)]
    struct IdentityFixture {
        candidates: Vec<Candidate>,
        threshold: f64,
        max_sensors: usize,
        placements: Vec<String>,
        sensors_placed: usize,
        prior_total_variance: f64,
        posterior_total_variance: f64,
        variance_reduction: f64,
        initial_evpi: f64,
        final_evpi: f64,
        decision_robust: bool,
        chosen_design: String,
        allocation: Vec<(String, f64)>,
        evpi_trace: Vec<f64>,
        posteriors: Vec<PosteriorSummary>,
        assimilation_colors: Vec<Color>,
        variance_color_dispersion: f64,
        evpi_color_dispersion: f64,
        realized: CampaignRealizedWorkPlan,
    }

    impl IdentityFixture {
        fn new() -> Self {
            let candidates = demo_candidates().expect("demo candidates");
            let allocation = candidates
                .iter()
                .enumerate()
                .map(|(index, candidate)| (candidate.name().to_string(), 0.1 + index as f64 * 0.01))
                .collect();
            let posteriors = candidates
                .iter()
                .map(|candidate| PosteriorSummary {
                    name: candidate.name().to_string(),
                    mean: candidate.prior_mean(),
                    variance: candidate.prior_variance(),
                })
                .collect();
            let plan = CampaignWorkPlan::checked(candidates.len(), 1).expect("fixture work plan");
            let realized = plan
                .realized(1, 1, candidates.len())
                .expect("fixture realized work plan");
            Self {
                candidates,
                threshold: 0.1,
                max_sensors: 1,
                placements: vec!["A".to_string()],
                sensors_placed: 1,
                prior_total_variance: 0.32,
                posterior_total_variance: 0.20,
                variance_reduction: 0.375,
                initial_evpi: 0.4,
                final_evpi: 0.2,
                decision_robust: false,
                chosen_design: "A".to_string(),
                allocation,
                evpi_trace: vec![0.4, 0.2],
                posteriors,
                assimilation_colors: vec![Color::Estimated {
                    estimator: "sensor-A-v1".to_string(),
                    dispersion: 0.01,
                }],
                variance_color_dispersion: f64::INFINITY,
                evpi_color_dispersion: 0.2,
                realized,
            }
        }

        fn source(&self) -> ReportIdentitySource<'_> {
            ReportIdentitySource {
                candidates: &self.candidates,
                threshold: self.threshold,
                max_sensors: self.max_sensors,
                outputs: ReportIdentityOutputs {
                    placements: &self.placements,
                    sensors_placed: self.sensors_placed,
                    prior_total_variance: self.prior_total_variance,
                    posterior_total_variance: self.posterior_total_variance,
                    variance_reduction: self.variance_reduction,
                    initial_evpi: self.initial_evpi,
                    final_evpi: self.final_evpi,
                    decision_robust: self.decision_robust,
                    chosen_design: &self.chosen_design,
                    allocation: &self.allocation,
                    evpi_trace: &self.evpi_trace,
                    posteriors: &self.posteriors,
                    assimilation_colors: &self.assimilation_colors,
                    variance_color_dispersion: self.variance_color_dispersion,
                    evpi_color_dispersion: self.evpi_color_dispersion,
                },
                plan: CampaignWorkPlan::checked(self.candidates.len(), self.max_sensors)
                    .expect("fixture work plan"),
                realized: self.realized,
            }
        }
    }

    fn with_test_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 3,
                    iteration: 4,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn report_identities(fixture: &IdentityFixture, rule: &[(f64, f64)]) -> (String, String) {
        with_test_cx(|cx| {
            let source = fixture.source();
            let mut progress = CampaignProgress::new(cx, 0);
            let variance =
                report_identity_with_rule("posterior-variance", &source, rule, &mut progress, cx)
                    .expect("variance identity");
            let evpi = report_identity_with_rule("evpi", &source, rule, &mut progress, cx)
                .expect("EVPI identity");
            (variance, evpi)
        })
    }

    #[test]
    fn scalar_variance_update_preserves_representable_extreme_posteriors() {
        for (prior, noise) in [
            (1.0e300, 1.0e-300),
            (1.0e-300, 1.0e300),
            (f64::MAX, f64::from_bits(1)),
        ] {
            let posterior = predicted_posterior_variance(prior, noise)
                .expect("positive finite scalar variances have a posterior");
            assert!(posterior > 0.0, "a representable posterior was erased");
            assert!(posterior <= prior.min(noise));
        }
    }

    #[test]
    fn normal_expectation_rule_is_positive_normalized_and_symmetric() {
        let weight_sum: f64 = NORMAL_EXPECTATION_RULE
            .iter()
            .map(|(_, weight)| weight)
            .sum();
        assert!((weight_sum - 1.0).abs() <= 8.0 * f64::EPSILON);
        for (left, right) in NORMAL_EXPECTATION_RULE
            .iter()
            .zip(NORMAL_EXPECTATION_RULE.iter().rev())
        {
            assert!(left.1 > 0.0);
            assert_eq!(left.0.to_bits(), canonicalize_zero(-right.0).to_bits());
            assert_eq!(left.1.to_bits(), right.1.to_bits());
        }
    }

    #[test]
    fn report_identity_versions_and_final_work_shape_are_locked() {
        assert_eq!(OED_REPORT_IDENTITY_VERSION, 5);
        assert_eq!(REPORT_ID_DOMAIN, "org.frankensim.fs-oed-e2e.report.v5");
        assert_eq!(CAMPAIGN_PLANNING_POLICY_VERSION, 3);
        assert_eq!(CAMPAIGN_POLL_POLICY_VERSION, 2);
        assert_eq!(fs_voi::EVPI_SEMANTICS_VERSION, 1);

        let plan = CampaignWorkPlan::checked(4, 12).expect("admitted work plan");
        assert_eq!(plan.maximum_finalization_work_units, 18 * 4 + 6 * 12 + 5);
        assert_eq!(plan.admitted_work_units, 2_545);
        assert_eq!(
            plan.realized(0, 0, 4)
                .expect("immediate STOP shape")
                .realized_work_units,
            97
        );
        assert_eq!(
            plan.realized(0, 1, 4)
                .expect("one completed zero-value action round")
                .realized_work_units,
            281
        );
        assert_eq!(
            plan.realized(12, 12, 4)
                .expect("full placement shape")
                .realized_work_units,
            plan.admitted_work_units
        );

        let identities = report_identities(&IdentityFixture::new(), &NORMAL_EXPECTATION_RULE);
        assert!(
            identities
                .0
                .starts_with("sensorforge-posterior-variance:v5:")
        );
        assert!(identities.1.starts_with("sensorforge-evpi:v5:"));
    }

    #[test]
    fn exact_normal_expectation_rule_bits_are_identity_semantic() {
        let fixture = IdentityFixture::new();
        let baseline = report_identities(&fixture, &NORMAL_EXPECTATION_RULE);

        let mut changed_node = NORMAL_EXPECTATION_RULE;
        changed_node[0].0 = f64::from_bits(changed_node[0].0.to_bits() ^ 1);
        let node_identity = report_identities(&fixture, &changed_node);
        assert_ne!(baseline.0, node_identity.0);
        assert_ne!(baseline.1, node_identity.1);

        let mut changed_weight = NORMAL_EXPECTATION_RULE;
        changed_weight[0].1 = f64::from_bits(changed_weight[0].1.to_bits() ^ 1);
        let weight_identity = report_identities(&fixture, &changed_weight);
        assert_ne!(baseline.0, weight_identity.0);
        assert_ne!(baseline.1, weight_identity.1);
    }

    #[test]
    fn every_sealed_report_output_moves_both_identities() {
        let baseline_fixture = IdentityFixture::new();
        let baseline = report_identities(&baseline_fixture, &NORMAL_EXPECTATION_RULE);
        let mut mutations = Vec::new();

        let mut changed = baseline_fixture.clone();
        changed.placements[0] = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.sensors_placed = 2;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.prior_total_variance = 0.33;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posterior_total_variance = 0.21;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.variance_reduction = 0.376;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.initial_evpi = 0.41;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.final_evpi = 0.21;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.decision_robust = true;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.chosen_design = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.allocation[0].0 = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.allocation[0].1 = 0.11;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.evpi_trace[0] = 0.41;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].name = "B".to_string();
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].mean = 0.61;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.posteriors[0].variance = 0.11;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.assimilation_colors[0] = Color::Estimated {
            estimator: "sensor-B-v1".to_string(),
            dispersion: 0.01,
        };
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.assimilation_colors[0] = Color::Estimated {
            estimator: "sensor-A-v1".to_string(),
            dispersion: 0.02,
        };
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.variance_color_dispersion = 1.0e300;
        mutations.push(changed);
        let mut changed = baseline_fixture.clone();
        changed.evpi_color_dispersion = 0.21;
        mutations.push(changed);

        for (index, mutation) in mutations.iter().enumerate() {
            let identity = report_identities(mutation, &NORMAL_EXPECTATION_RULE);
            assert_ne!(
                baseline.0, identity.0,
                "variance identity ignored output mutation {index}"
            );
            assert_ne!(
                baseline.1, identity.1,
                "EVPI identity ignored output mutation {index}"
            );
        }

        for index in 0..8 {
            let mut changed = baseline_fixture.clone();
            match index {
                0 => changed.realized.completed_placements += 1,
                1 => changed.realized.action_rounds += 1,
                2 => changed.realized.positive_prior_candidates += 1,
                3 => changed.realized.setup_work_units += 1,
                4 => changed.realized.placement_work_units += 1,
                5 => changed.realized.incomplete_action_work_units += 1,
                6 => changed.realized.finalization_work_units += 1,
                7 => changed.realized.realized_work_units += 1,
                _ => unreachable!("retained realized-work field count"),
            }
            let identity = report_identities(&changed, &NORMAL_EXPECTATION_RULE);
            assert_ne!(
                baseline.0, identity.0,
                "variance identity ignored realized-work field {index}"
            );
            assert_ne!(
                baseline.1, identity.1,
                "EVPI identity ignored realized-work field {index}"
            );
        }
    }
}
