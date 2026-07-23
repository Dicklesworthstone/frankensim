//! Graph-aware certify-or-escalate routing.
//!
//! This module translates a model-form escalation request into an exact
//! [`fs_ladder::FidelityGraph`] route. It does not own graph adequacy: the
//! graph's `cheapest_adequate` query remains the authority. The returned plan,
//! gap, and calibration records are content-addressed so HELM can persist
//! them, but this L4 crate does not depend upward on a ledger implementation.

use crate::{ConformalBand, Decision, certify_or_escalate};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::EscalationAdvice;
use fs_evidence::uncertainty::{
    BudgetTotal, DominantEngineeringTerm, EngineeringUncertaintyBudget, EngineeringUncertaintyKind,
    TermValue,
};
use fs_ladder::{
    Adequacy, EdgeEvidenceResolver, EdgeId, FidelityGraph, FidelityGraphId, ModelId,
    ModelRecommendation, QoiId, QueryContext, QueryEvidenceRef, QueryRefusal, SelectionBasis,
};
use std::collections::BTreeSet;
use std::fmt;

const PLAN_IDENTITY_DOMAIN: &str = "frankensim.fs-surrogate.escalation-plan.v1";
const GAP_IDENTITY_DOMAIN: &str = "frankensim.fs-surrogate.escalation-gap.v1";
const CALIBRATION_IDENTITY_DOMAIN: &str = "frankensim.fs-surrogate.escalation-calibration.v1";

/// Exact evidence artifact produced by the escalated run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EscalationRunEvidenceRef(ContentHash);

impl EscalationRunEvidenceRef {
    /// Wrap a content-addressed run receipt.
    #[must_use]
    pub const fn new(value: ContentHash) -> Self {
        Self(value)
    }

    /// Underlying content identity.
    #[must_use]
    pub const fn hash(self) -> ContentHash {
        self.0
    }
}

/// Why the graph planner was invoked.
#[derive(Debug, Clone, PartialEq)]
pub enum EscalationCause {
    /// A finite conformal band exceeded the decision tolerance.
    BandTooWide {
        /// Observed band half-width.
        half_width: f64,
        /// Decision-relevant half-width.
        tolerance: f64,
    },
    /// The query was outside the surrogate's admitted validity domain.
    OutsideValidityDomain,
    /// The conformal policy state was malformed or unbounded.
    InvalidConformalPolicy,
    /// A caller invoked the planner directly after model-form dominance.
    ModelFormDominant,
}

/// One ordered edge in an escalation walk.
#[derive(Debug, Clone, PartialEq)]
pub struct EscalationStep {
    /// Zero-based execution order.
    pub rank: u32,
    /// Exact graph edge.
    pub edge: EdgeId,
    /// Model being left.
    pub source: ModelId,
    /// More-informative model being entered.
    pub target: ModelId,
    /// Predicted total target-model cost in seconds.
    pub predicted_target_cost_s: f64,
    /// Pairwise discrepancy assessment retained by the graph query.
    pub assessed_relative_discrepancy: Option<f64>,
    /// Exact resolver receipt used by the graph query.
    pub evidence: QueryEvidenceRef,
}

/// Estimated post-route uncertainty.
///
/// This is planning data, not a numerical enclosure. The estimate replaces
/// the independent model-form term with the selected model's maximum adequate
/// pairwise discrepancy. Unknown totals and correlated model-form terms remain
/// explicitly indeterminate.
#[derive(Debug, Clone, PartialEq)]
pub enum PredictedPostUncertainty {
    /// A finite eight-term estimate can be formed.
    Bounded {
        /// Pre-escalation conservative half-width in the budget unit.
        prior_half_width: f64,
        /// Pre-escalation relative half-width.
        prior_relative: f64,
        /// Predicted model-form half-width after escalation.
        predicted_model_form_half_width: f64,
        /// Predicted total half-width after escalation.
        predicted_total_half_width: f64,
        /// Predicted relative total half-width after escalation.
        predicted_total_relative: f64,
    },
    /// No finite combined prediction is justified.
    Indeterminate {
        /// Stable reason suitable for an agent diagnostic.
        reason: &'static str,
    },
}

impl PredictedPostUncertainty {
    /// Predicted total relative uncertainty when finite.
    #[must_use]
    pub const fn total_relative(&self) -> Option<f64> {
        match self {
            Self::Bounded {
                predicted_total_relative,
                ..
            } => Some(*predicted_total_relative),
            Self::Indeterminate { .. } => None,
        }
    }

    /// Pre-escalation relative uncertainty when finite.
    #[must_use]
    pub const fn prior_relative(&self) -> Option<f64> {
        match self {
            Self::Bounded { prior_relative, .. } => Some(*prior_relative),
            Self::Indeterminate { .. } => None,
        }
    }
}

/// A deterministic, content-addressed graph escalation plan.
#[derive(Debug, Clone, PartialEq)]
pub struct EscalationPlan {
    cause: EscalationCause,
    advice: EscalationAdvice,
    budget: ContentHash,
    unit: String,
    reference_magnitude: f64,
    recommendation: ModelRecommendation,
    steps: Vec<EscalationStep>,
    predicted: PredictedPostUncertainty,
}

impl EscalationPlan {
    /// Trigger that led to this route.
    #[must_use]
    pub const fn cause(&self) -> &EscalationCause {
        &self.cause
    }

    /// Directional evidence advice that authorized graph routing.
    #[must_use]
    pub const fn advice(&self) -> EscalationAdvice {
        self.advice
    }

    /// Exact eight-term budget identity.
    #[must_use]
    pub const fn budget(&self) -> ContentHash {
        self.budget
    }

    /// Unit shared by the budget and prediction.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// Absolute QoI magnitude used to convert absolute and relative widths.
    #[must_use]
    pub const fn reference_magnitude(&self) -> f64 {
        self.reference_magnitude
    }

    /// Graph-owned recommendation and complete replay explanation.
    #[must_use]
    pub const fn recommendation(&self) -> &ModelRecommendation {
        &self.recommendation
    }

    /// Ordered multi-edge walk.
    #[must_use]
    pub fn steps(&self) -> &[EscalationStep] {
        &self.steps
    }

    /// First edge to execute.
    #[must_use]
    pub fn next_edge(&self) -> EdgeId {
        self.steps[0].edge
    }

    /// Estimated post-route uncertainty.
    #[must_use]
    pub const fn predicted_uncertainty(&self) -> &PredictedPostUncertainty {
        &self.predicted
    }

    /// Canonical semantic bytes for ledger persistence.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        put_u16(&mut out, 1);
        encode_cause(&mut out, &self.cause);
        put_u8(&mut out, advice_tag(self.advice));
        put_hash(&mut out, self.budget);
        put_string(&mut out, &self.unit);
        put_f64(&mut out, self.reference_magnitude);
        encode_recommendation(&mut out, &self.recommendation);
        put_u32(&mut out, self.steps.len() as u32);
        for step in &self.steps {
            put_u32(&mut out, step.rank);
            put_edge(&mut out, step.edge);
            put_model(&mut out, step.source);
            put_model(&mut out, step.target);
            put_f64(&mut out, step.predicted_target_cost_s);
            put_optional_f64(&mut out, step.assessed_relative_discrepancy);
            put_hash(&mut out, step.evidence.hash());
        }
        encode_prediction(&mut out, &self.predicted);
        out
    }

    /// Domain-separated identity over the complete plan and explanation.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(PLAN_IDENTITY_DOMAIN, &self.canonical_bytes())
    }
}

/// Evidence work suggested when no route can be authorized.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProbeSuggestion {
    /// Resolve or refresh the cost/discrepancy evidence on an existing edge.
    ResolveEdge {
        /// Candidate edge.
        edge: EdgeId,
        /// Less-informative endpoint.
        source: ModelId,
        /// More-informative endpoint.
        target: ModelId,
    },
    /// Add and probe an outgoing edge for a model that currently has none.
    AddOutgoingEdge {
        /// Model with no graph exit.
        source: ModelId,
        /// Failing QoI.
        qoi: QoiId,
    },
    /// Register a current model that is absent from this graph.
    AddModelNode {
        /// Missing model identity.
        model: ModelId,
    },
}

/// Stable reason why escalation was not authorized.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EscalationGapReason {
    /// The directional advice did not request model fidelity.
    AdviceDoesNotPermitModelEscalation,
    /// The eight-term budget does not identify model form as dominant.
    BudgetDoesNotSupportModelEscalation,
    /// Budget and query name different QoIs.
    QoiMismatch,
    /// The graph tolerance differs from the certify-or-escalate tolerance.
    ToleranceMismatch,
    /// The reference magnitude is non-positive or non-finite.
    InvalidReferenceMagnitude,
    /// The graph judges the current model already adequate.
    CurrentModelAlreadyAdequate,
    /// The current model is not registered in the graph.
    UnknownCurrentModel,
    /// No evidenced route reaches an adequate model.
    NoAdequateRoute,
    /// A route exists only outside the session cost budget.
    SessionBudgetBlocked,
    /// The failure signal has not crossed or rearmed the hysteresis band.
    HysteresisHold,
    /// The route would exceed the run's escalation-edge cap.
    EscalationLimitReached,
    /// The route would reuse an edge already issued in this run.
    RepeatedEdgeBlocked,
}

/// Content-addressed refusal with concrete evidence-acquisition suggestions.
#[derive(Debug, Clone, PartialEq)]
pub struct EscalationGap {
    reason: EscalationGapReason,
    advice: EscalationAdvice,
    graph: FidelityGraphId,
    budget: ContentHash,
    current: ModelId,
    qoi: QoiId,
    budget_s: f64,
    decision_tolerance: f64,
    reference_magnitude: f64,
    candidate_plan: Option<ContentHash>,
    candidate_edges: Vec<EdgeId>,
    suggested_probes: Vec<ProbeSuggestion>,
    required_cost_s: Option<f64>,
}

impl EscalationGap {
    /// Why no route was authorized.
    #[must_use]
    pub const fn reason(&self) -> EscalationGapReason {
        self.reason
    }

    /// Exact directional advice presented to the planner.
    #[must_use]
    pub const fn advice(&self) -> EscalationAdvice {
        self.advice
    }

    /// Exact graph that refused the route.
    #[must_use]
    pub const fn graph(&self) -> FidelityGraphId {
        self.graph
    }

    /// Exact eight-term budget that triggered the refusal.
    #[must_use]
    pub const fn budget(&self) -> ContentHash {
        self.budget
    }

    /// Current model supplied to the planner.
    #[must_use]
    pub const fn current_model(&self) -> ModelId {
        self.current
    }

    /// Candidate plan computed before a policy-only refusal, if any.
    #[must_use]
    pub const fn candidate_plan(&self) -> Option<ContentHash> {
        self.candidate_plan
    }

    /// Existing graph edges relevant for evidence work.
    #[must_use]
    pub fn candidate_edges(&self) -> &[EdgeId] {
        &self.candidate_edges
    }

    /// Canonically ordered evidence-work suggestions.
    #[must_use]
    pub fn suggested_probes(&self) -> &[ProbeSuggestion] {
        &self.suggested_probes
    }

    /// Predicted cost of the otherwise adequate route when budget blocked it.
    #[must_use]
    pub const fn required_cost_s(&self) -> Option<f64> {
        self.required_cost_s
    }

    /// Canonical semantic bytes for ledger persistence.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        put_u16(&mut out, 1);
        put_u8(&mut out, gap_reason_tag(self.reason));
        put_u8(&mut out, advice_tag(self.advice));
        put_hash(&mut out, self.graph.hash());
        put_hash(&mut out, self.budget);
        put_model(&mut out, self.current);
        put_string(&mut out, self.qoi.as_str());
        put_f64(&mut out, self.budget_s);
        put_f64(&mut out, self.decision_tolerance);
        put_f64(&mut out, self.reference_magnitude);
        put_optional_hash(&mut out, self.candidate_plan);
        put_u32(&mut out, self.candidate_edges.len() as u32);
        for edge in &self.candidate_edges {
            put_edge(&mut out, *edge);
        }
        put_u32(&mut out, self.suggested_probes.len() as u32);
        for suggestion in &self.suggested_probes {
            encode_probe(&mut out, suggestion);
        }
        put_optional_f64(&mut out, self.required_cost_s);
        out
    }

    /// Domain-separated identity over the complete gap report.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(GAP_IDENTITY_DOMAIN, &self.canonical_bytes())
    }
}

/// Strict hysteresis and per-run escalation limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EscalationPolicy {
    trigger_factor: f64,
    rearm_factor: f64,
    max_edges_per_run: u32,
}

impl EscalationPolicy {
    /// Admit a strict hysteresis band around the decision tolerance.
    ///
    /// `trigger_factor` must be at least one, `rearm_factor` must be in
    /// `[0, 1)`, and the trigger must strictly exceed the rearm boundary.
    pub fn try_new(
        trigger_factor: f64,
        rearm_factor: f64,
        max_edges_per_run: u32,
    ) -> Result<Self, EscalationPolicyError> {
        if !trigger_factor.is_finite()
            || !rearm_factor.is_finite()
            || trigger_factor < 1.0
            || !(0.0..1.0).contains(&rearm_factor)
            || trigger_factor <= rearm_factor
        {
            return Err(EscalationPolicyError::InvalidHysteresis);
        }
        if max_edges_per_run == 0 {
            return Err(EscalationPolicyError::ZeroEscalationLimit);
        }
        Ok(Self {
            trigger_factor,
            rearm_factor,
            max_edges_per_run,
        })
    }

    /// Width multiplier that triggers escalation.
    #[must_use]
    pub const fn trigger_factor(self) -> f64 {
        self.trigger_factor
    }

    /// Width multiplier below which the latch rearms.
    #[must_use]
    pub const fn rearm_factor(self) -> f64 {
        self.rearm_factor
    }

    /// Maximum graph edges that one run may issue.
    #[must_use]
    pub const fn max_edges_per_run(self) -> u32 {
        self.max_edges_per_run
    }
}

impl Default for EscalationPolicy {
    fn default() -> Self {
        Self {
            trigger_factor: 1.05,
            rearm_factor: 0.90,
            max_edges_per_run: 3,
        }
    }
}

/// Policy admission failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EscalationPolicyError {
    /// The hysteresis band is non-finite, inverted, or degenerate.
    InvalidHysteresis,
    /// A zero edge cap would make every policy unusable.
    ZeroEscalationLimit,
}

impl fmt::Display for EscalationPolicyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidHysteresis => f.write_str(
                "escalation hysteresis needs finite trigger >= 1, rearm in [0,1), and trigger > rearm",
            ),
            Self::ZeroEscalationLimit => {
                f.write_str("escalation policy needs at least one edge per run")
            }
        }
    }
}

impl std::error::Error for EscalationPolicyError {}

/// Mutable, run-local anti-thrash state.
#[derive(Debug, Clone)]
pub struct EscalationSession {
    policy: EscalationPolicy,
    latched: bool,
    consumed_edges: u32,
    issued_edges: BTreeSet<EdgeId>,
}

impl EscalationSession {
    /// Start an unlatched run.
    #[must_use]
    pub fn new(policy: EscalationPolicy) -> Self {
        Self {
            policy,
            latched: false,
            consumed_edges: 0,
            issued_edges: BTreeSet::new(),
        }
    }

    /// Number of graph edges already authorized.
    #[must_use]
    pub const fn consumed_edges(&self) -> u32 {
        self.consumed_edges
    }

    /// Whether a prior escalation is waiting for recovery below the rearm
    /// boundary.
    #[must_use]
    pub const fn is_latched(&self) -> bool {
        self.latched
    }

    fn observe_accepted_surrogate(&mut self, width: f64, tolerance: f64) {
        if relative_to_tolerance(width, tolerance) <= self.policy.rearm_factor {
            self.latched = false;
        }
    }

    fn authorize(
        &mut self,
        plan: &EscalationPlan,
        failure_width: f64,
        tolerance: f64,
    ) -> Result<(), EscalationGapReason> {
        if relative_to_tolerance(failure_width, tolerance) < self.policy.trigger_factor
            || self.latched
        {
            return Err(EscalationGapReason::HysteresisHold);
        }
        let edge_count = u32::try_from(plan.steps.len()).unwrap_or(u32::MAX);
        if self
            .consumed_edges
            .checked_add(edge_count)
            .is_none_or(|total| total > self.policy.max_edges_per_run)
        {
            return Err(EscalationGapReason::EscalationLimitReached);
        }
        if plan
            .steps
            .iter()
            .any(|step| self.issued_edges.contains(&step.edge))
        {
            return Err(EscalationGapReason::RepeatedEdgeBlocked);
        }
        for step in &plan.steps {
            self.issued_edges.insert(step.edge);
        }
        self.consumed_edges += edge_count;
        self.latched = true;
        Ok(())
    }
}

/// Graph-aware certify-or-escalate result.
#[derive(Debug, Clone, PartialEq)]
pub enum GraphDecision {
    /// Existing compatibility policy admitted the surrogate.
    UseSurrogate {
        /// Band half-width backing the decision.
        band_half_width: f64,
    },
    /// Execute the first edge of this retained route.
    Escalate {
        /// Exact multi-edge plan.
        plan: Box<EscalationPlan>,
    },
    /// Refuse blind escalation and acquire the named evidence instead.
    Gap {
        /// Exact refusal.
        gap: Box<EscalationGap>,
    },
}

/// Plan the cheapest graph route whose selected model has adequate pairwise
/// discrepancy evidence for the exact query.
pub fn plan_graph_escalation(
    graph: &FidelityGraph,
    resolver: &impl EdgeEvidenceResolver,
    current: ModelId,
    context: &QueryContext,
    budget: &EngineeringUncertaintyBudget,
    advice: EscalationAdvice,
    reference_magnitude: f64,
) -> Result<EscalationPlan, Box<EscalationGap>> {
    plan_graph_escalation_with_cause(
        graph,
        resolver,
        current,
        context,
        budget,
        advice,
        reference_magnitude,
        EscalationCause::ModelFormDominant,
    )
}

/// Compatibility-preserving graph-aware policy.
///
/// The existing [`crate::certify_or_escalate`] verdict is evaluated first.
/// Only its `Escalate` branch can invoke graph routing. The legacy band and
/// `decision_tolerance` are in the budget's QoI unit; this adapter divides
/// both by `reference_magnitude` before comparing them to the graph's explicit
/// relative-discrepancy tolerance.
#[allow(clippy::too_many_arguments)]
pub fn certify_or_escalate_with_graph(
    band: &ConformalBand,
    in_validity_domain: bool,
    decision_tolerance: f64,
    graph: &FidelityGraph,
    resolver: &impl EdgeEvidenceResolver,
    current: ModelId,
    context: &QueryContext,
    budget: &EngineeringUncertaintyBudget,
    advice: EscalationAdvice,
    reference_magnitude: f64,
    session: &mut EscalationSession,
) -> GraphDecision {
    match certify_or_escalate(band, in_validity_domain, decision_tolerance) {
        Decision::UseSurrogate { band_half_width } => {
            session.observe_accepted_surrogate(band_half_width, decision_tolerance);
            GraphDecision::UseSurrogate { band_half_width }
        }
        Decision::Escalate { .. } => {
            let cause = escalation_cause(band, in_validity_domain, decision_tolerance);
            let plan = match plan_graph_escalation_with_cause(
                graph,
                resolver,
                current,
                context,
                budget,
                advice,
                reference_magnitude,
                cause,
            ) {
                Ok(plan) => plan,
                Err(gap) => return GraphDecision::Gap { gap },
            };
            let relative_tolerance = decision_tolerance / reference_magnitude;
            if !relative_tolerance.is_finite()
                || relative_tolerance < 0.0
                || context.max_relative_discrepancy().to_bits()
                    != canonical_zero(relative_tolerance).to_bits()
            {
                return GraphDecision::Gap {
                    gap: Box::new(gap_from_plan(&plan, EscalationGapReason::ToleranceMismatch)),
                };
            }
            let failure_width =
                failure_width(band, budget, reference_magnitude).unwrap_or(f64::INFINITY);
            match session.authorize(&plan, failure_width, relative_tolerance) {
                Ok(()) => GraphDecision::Escalate {
                    plan: Box::new(plan),
                },
                Err(reason) => GraphDecision::Gap {
                    gap: Box::new(gap_from_plan(&plan, reason)),
                },
            }
        }
    }
}

/// Actual-run observation used to calibrate the planner itself.
#[derive(Debug, Clone, PartialEq)]
pub struct EscalationCalibrationRecord {
    plan: ContentHash,
    actual_budget: ContentHash,
    run_evidence: EscalationRunEvidenceRef,
    predicted_cost_s: f64,
    actual_cost_s: f64,
    prior_relative: Option<f64>,
    predicted_post_relative: Option<f64>,
    actual_post_relative: Option<f64>,
}

impl EscalationCalibrationRecord {
    /// Bind an actual escalated run to its prediction.
    pub fn try_new(
        plan: &EscalationPlan,
        actual_budget: &EngineeringUncertaintyBudget,
        actual_cost_s: f64,
        run_evidence: EscalationRunEvidenceRef,
    ) -> Result<Self, CalibrationError> {
        if actual_budget.qoi() != plan.recommendation.explanation.qoi.as_str() {
            return Err(CalibrationError::QoiMismatch);
        }
        if actual_budget.unit() != plan.unit {
            return Err(CalibrationError::UnitMismatch);
        }
        if !actual_cost_s.is_finite() || actual_cost_s < 0.0 {
            return Err(CalibrationError::InvalidActualCost);
        }
        let scale = plan.reference_magnitude;
        let actual_post_relative = bounded_total(actual_budget).and_then(|value| {
            let relative = value / scale;
            relative.is_finite().then_some(relative)
        });
        Ok(Self {
            plan: plan.content_id(),
            actual_budget: actual_budget.content_id(),
            run_evidence,
            predicted_cost_s: plan.recommendation.predicted_cost_s,
            actual_cost_s: canonical_zero(actual_cost_s),
            prior_relative: plan.predicted.prior_relative(),
            predicted_post_relative: plan.predicted.total_relative(),
            actual_post_relative,
        })
    }

    /// Exact plan whose prediction is being checked.
    #[must_use]
    pub const fn plan(&self) -> ContentHash {
        self.plan
    }

    /// Exact post-run eight-term budget.
    #[must_use]
    pub const fn actual_budget(&self) -> ContentHash {
        self.actual_budget
    }

    /// Evidence receipt for the escalated run.
    #[must_use]
    pub const fn run_evidence(&self) -> EscalationRunEvidenceRef {
        self.run_evidence
    }

    /// Predicted total execution cost in seconds.
    #[must_use]
    pub const fn predicted_cost_s(&self) -> f64 {
        self.predicted_cost_s
    }

    /// Observed total execution cost in seconds.
    #[must_use]
    pub const fn actual_cost_s(&self) -> f64 {
        self.actual_cost_s
    }

    /// Predicted post-escalation relative uncertainty, when finite.
    #[must_use]
    pub const fn predicted_post_relative(&self) -> Option<f64> {
        self.predicted_post_relative
    }

    /// Observed post-escalation relative uncertainty, when finite.
    #[must_use]
    pub const fn actual_post_relative(&self) -> Option<f64> {
        self.actual_post_relative
    }

    /// Predicted minus actual execution cost.
    #[must_use]
    pub fn cost_error_s(&self) -> f64 {
        self.actual_cost_s - self.predicted_cost_s
    }

    /// Predicted reduction in relative total uncertainty.
    #[must_use]
    pub fn predicted_improvement(&self) -> Option<f64> {
        Some(self.prior_relative? - self.predicted_post_relative?)
    }

    /// Observed reduction in relative total uncertainty.
    #[must_use]
    pub fn actual_improvement(&self) -> Option<f64> {
        Some(self.prior_relative? - self.actual_post_relative?)
    }

    /// Prediction error for post-escalation relative uncertainty.
    #[must_use]
    pub fn uncertainty_error(&self) -> Option<f64> {
        Some(self.actual_post_relative? - self.predicted_post_relative?)
    }

    /// Canonical semantic bytes for ledger persistence.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        put_u16(&mut out, 1);
        put_hash(&mut out, self.plan);
        put_hash(&mut out, self.actual_budget);
        put_hash(&mut out, self.run_evidence.hash());
        put_f64(&mut out, self.predicted_cost_s);
        put_f64(&mut out, self.actual_cost_s);
        put_optional_f64(&mut out, self.prior_relative);
        put_optional_f64(&mut out, self.predicted_post_relative);
        put_optional_f64(&mut out, self.actual_post_relative);
        out
    }

    /// Domain-separated identity over prediction, outcome, and exact evidence.
    #[must_use]
    pub fn content_id(&self) -> ContentHash {
        hash_domain(CALIBRATION_IDENTITY_DOMAIN, &self.canonical_bytes())
    }
}

/// Calibration record admission failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationError {
    /// Actual budget names a different QoI.
    QoiMismatch,
    /// Actual budget uses a different unit.
    UnitMismatch,
    /// Actual cost is negative or non-finite.
    InvalidActualCost,
}

impl fmt::Display for CalibrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::QoiMismatch => f.write_str("actual escalation budget names a different QoI"),
            Self::UnitMismatch => f.write_str("actual escalation budget uses a different unit"),
            Self::InvalidActualCost => {
                f.write_str("actual escalation cost must be finite and non-negative")
            }
        }
    }
}

impl std::error::Error for CalibrationError {}

#[allow(clippy::too_many_arguments)]
fn plan_graph_escalation_with_cause(
    graph: &FidelityGraph,
    resolver: &impl EdgeEvidenceResolver,
    current: ModelId,
    context: &QueryContext,
    budget: &EngineeringUncertaintyBudget,
    advice: EscalationAdvice,
    reference_magnitude: f64,
    cause: EscalationCause,
) -> Result<EscalationPlan, Box<EscalationGap>> {
    let gap_context = GapContext {
        graph,
        current,
        query: context,
        budget: budget.content_id(),
        advice,
        reference_magnitude,
    };
    let base_gap = |reason| Box::new(make_gap(gap_context, reason, None));
    if advice != EscalationAdvice::EscalateModelFidelity {
        return Err(base_gap(
            EscalationGapReason::AdviceDoesNotPermitModelEscalation,
        ));
    }
    if !matches!(
        budget.dominant(),
        DominantEngineeringTerm::Known {
            kind: EngineeringUncertaintyKind::ModelForm,
            ..
        }
    ) {
        return Err(base_gap(
            EscalationGapReason::BudgetDoesNotSupportModelEscalation,
        ));
    }
    if budget.qoi() != context.qoi().as_str() {
        return Err(base_gap(EscalationGapReason::QoiMismatch));
    }
    if !reference_magnitude.is_finite() || reference_magnitude <= 0.0 {
        return Err(base_gap(EscalationGapReason::InvalidReferenceMagnitude));
    }

    let recommendation = match graph.cheapest_adequate(current, context, resolver) {
        Ok(recommendation) => recommendation,
        Err(refusal) => {
            return Err(Box::new(query_gap(resolver, gap_context, &refusal)));
        }
    };
    if recommendation.path.is_empty() {
        return Err(base_gap(EscalationGapReason::CurrentModelAlreadyAdequate));
    }

    let mut steps = Vec::with_capacity(recommendation.path.len());
    for (rank, edge) in recommendation.path.iter().enumerate() {
        let Some(row) = recommendation
            .explanation
            .considered
            .iter()
            .find(|row| row.edge == *edge)
        else {
            return Err(base_gap(EscalationGapReason::NoAdequateRoute));
        };
        steps.push(EscalationStep {
            rank: u32::try_from(rank).unwrap_or(u32::MAX),
            edge: row.edge,
            source: row.source,
            target: row.target,
            predicted_target_cost_s: row.target_cost_s,
            assessed_relative_discrepancy: row.assessed_relative_discrepancy,
            evidence: row.evidence,
        });
    }
    let selected_discrepancy = recommendation
        .explanation
        .considered
        .iter()
        .filter(|row| {
            row.source == recommendation.model && row.source_adequacy == Adequacy::Adequate
        })
        .filter_map(|row| row.assessed_relative_discrepancy)
        .max_by(f64::total_cmp);
    let predicted = predict_post_uncertainty(budget, reference_magnitude, selected_discrepancy);
    Ok(EscalationPlan {
        cause,
        advice,
        budget: budget.content_id(),
        unit: budget.unit().to_string(),
        reference_magnitude,
        recommendation,
        steps,
        predicted,
    })
}

#[derive(Clone, Copy)]
struct GapContext<'a> {
    graph: &'a FidelityGraph,
    current: ModelId,
    query: &'a QueryContext,
    budget: ContentHash,
    advice: EscalationAdvice,
    reference_magnitude: f64,
}

fn query_gap(
    resolver: &impl EdgeEvidenceResolver,
    request: GapContext<'_>,
    refusal: &QueryRefusal,
) -> EscalationGap {
    if matches!(refusal, QueryRefusal::UnknownStart(_)) {
        return make_gap(request, EscalationGapReason::UnknownCurrentModel, None);
    }
    let unlimited = QueryContext::new(
        request.query.qoi().clone(),
        request
            .query
            .regime()
            .iter()
            .map(|(axis, value)| (axis.clone(), *value)),
        request.query.problem_size(),
        f64::MAX,
        request.query.max_relative_discrepancy(),
    )
    .expect("an admitted context remains admitted with a finite maximum budget");
    if let Ok(candidate) = request
        .graph
        .cheapest_adequate(request.current, &unlimited, resolver)
        && !candidate.path.is_empty()
        && candidate.predicted_cost_s > request.query.budget_s()
    {
        return make_gap(
            request,
            EscalationGapReason::SessionBudgetBlocked,
            Some(candidate.predicted_cost_s),
        );
    }
    make_gap(request, EscalationGapReason::NoAdequateRoute, None)
}

fn make_gap(
    request: GapContext<'_>,
    reason: EscalationGapReason,
    required_cost_s: Option<f64>,
) -> EscalationGap {
    let outgoing = request
        .graph
        .edges()
        .values()
        .filter(|edge| edge.source() == request.current)
        .collect::<Vec<_>>();
    let candidate_edges = outgoing.iter().map(|edge| edge.id()).collect();
    let suggested_probes = if reason == EscalationGapReason::UnknownCurrentModel {
        vec![ProbeSuggestion::AddModelNode {
            model: request.current,
        }]
    } else if outgoing.is_empty() {
        vec![ProbeSuggestion::AddOutgoingEdge {
            source: request.current,
            qoi: request.query.qoi().clone(),
        }]
    } else {
        outgoing
            .into_iter()
            .map(|edge| ProbeSuggestion::ResolveEdge {
                edge: edge.id(),
                source: edge.source(),
                target: edge.target(),
            })
            .collect()
    };
    EscalationGap {
        reason,
        advice: request.advice,
        graph: request.graph.identity(),
        budget: request.budget,
        current: request.current,
        qoi: request.query.qoi().clone(),
        budget_s: request.query.budget_s(),
        decision_tolerance: request.query.max_relative_discrepancy(),
        reference_magnitude: request.reference_magnitude,
        candidate_plan: None,
        candidate_edges,
        suggested_probes,
        required_cost_s,
    }
}

fn gap_from_plan(plan: &EscalationPlan, reason: EscalationGapReason) -> EscalationGap {
    EscalationGap {
        reason,
        advice: plan.advice,
        graph: plan.recommendation.explanation.graph,
        budget: plan.budget,
        current: plan.recommendation.explanation.start,
        qoi: plan.recommendation.explanation.qoi.clone(),
        budget_s: plan.recommendation.explanation.budget_s,
        decision_tolerance: plan.recommendation.explanation.max_relative_discrepancy,
        reference_magnitude: plan.reference_magnitude,
        candidate_plan: Some(plan.content_id()),
        candidate_edges: plan.steps.iter().map(|step| step.edge).collect(),
        suggested_probes: plan
            .steps
            .iter()
            .map(|step| ProbeSuggestion::ResolveEdge {
                edge: step.edge,
                source: step.source,
                target: step.target,
            })
            .collect(),
        required_cost_s: Some(plan.recommendation.predicted_cost_s),
    }
}

fn escalation_cause(
    band: &ConformalBand,
    in_validity_domain: bool,
    tolerance: f64,
) -> EscalationCause {
    if !in_validity_domain {
        EscalationCause::OutsideValidityDomain
    } else if band.half_width.is_finite()
        && band.half_width >= 0.0
        && band.alpha.is_finite()
        && band.alpha > 0.0
        && band.alpha < 1.0
        && tolerance.is_finite()
        && tolerance >= 0.0
    {
        EscalationCause::BandTooWide {
            half_width: band.half_width,
            tolerance,
        }
    } else {
        EscalationCause::InvalidConformalPolicy
    }
}

fn failure_width(
    band: &ConformalBand,
    budget: &EngineeringUncertaintyBudget,
    reference_magnitude: f64,
) -> Option<f64> {
    let model_relative = term_half_width(
        budget.term(EngineeringUncertaintyKind::ModelForm).value(),
        EngineeringUncertaintyKind::ModelForm,
    )
    .map(|value| value / reference_magnitude);
    if band.half_width.is_finite() && band.half_width >= 0.0 {
        let band_relative = band.half_width / reference_magnitude;
        Some(model_relative.map_or(band_relative, |value| value.max(band_relative)))
    } else {
        model_relative
    }
}

fn predict_post_uncertainty(
    budget: &EngineeringUncertaintyBudget,
    scale: f64,
    selected_discrepancy: Option<f64>,
) -> PredictedPostUncertainty {
    let Some(selected_discrepancy) = selected_discrepancy else {
        return PredictedPostUncertainty::Indeterminate {
            reason: "selected model has no finite discrepancy prediction",
        };
    };
    let Some(prior_half_width) = bounded_total(budget) else {
        return PredictedPostUncertainty::Indeterminate {
            reason: "eight-term budget total is unknown or unbounded",
        };
    };
    let model_term = budget.term(EngineeringUncertaintyKind::ModelForm).value();
    if matches!(model_term, TermValue::CorrelatedBlock(_)) {
        return PredictedPostUncertainty::Indeterminate {
            reason: "model-form term belongs to a covariance block",
        };
    }
    let Some(prior_model_form) = term_half_width(model_term, EngineeringUncertaintyKind::ModelForm)
    else {
        return PredictedPostUncertainty::Indeterminate {
            reason: "model-form half-width is unknown",
        };
    };
    let predicted_model_form_half_width = selected_discrepancy * scale;
    let predicted_total_half_width =
        (prior_half_width - prior_model_form).max(0.0) + predicted_model_form_half_width;
    let prior_relative = prior_half_width / scale;
    let predicted_total_relative = predicted_total_half_width / scale;
    if !predicted_model_form_half_width.is_finite()
        || !predicted_total_half_width.is_finite()
        || !prior_relative.is_finite()
        || !predicted_total_relative.is_finite()
        || (selected_discrepancy > 0.0 && predicted_model_form_half_width == 0.0)
    {
        return PredictedPostUncertainty::Indeterminate {
            reason: "finite prediction arithmetic overflowed or underflowed",
        };
    }
    PredictedPostUncertainty::Bounded {
        prior_half_width,
        prior_relative,
        predicted_model_form_half_width,
        predicted_total_half_width,
        predicted_total_relative,
    }
}

fn bounded_total(budget: &EngineeringUncertaintyBudget) -> Option<f64> {
    match budget.total() {
        BudgetTotal::Bounded {
            conservative_half_width,
        } => Some(conservative_half_width),
        BudgetTotal::Unknown { .. } | BudgetTotal::Unbounded { .. } => None,
    }
}

fn term_half_width(value: &TermValue, kind: EngineeringUncertaintyKind) -> Option<f64> {
    match value {
        TermValue::IntervalBound { upper, .. } => Some(*upper),
        TermValue::Distribution(summary) => Some(summary.conservative_half_width),
        TermValue::Ensemble(summary) => Some(summary.conservative_half_width),
        TermValue::CorrelatedBlock(block) => {
            let index = block.members().iter().position(|member| *member == kind)?;
            Some(block.covariance()[index * block.members().len() + index].sqrt())
        }
        TermValue::Negligible { .. } => Some(0.0),
        _ => None,
    }
}

fn encode_recommendation(out: &mut Vec<u8>, value: &ModelRecommendation) {
    put_model(out, value.model);
    put_f64(out, value.predicted_cost_s);
    put_u32(out, value.path.len() as u32);
    for edge in &value.path {
        put_edge(out, *edge);
    }
    let explanation = &value.explanation;
    put_hash(out, explanation.graph.hash());
    put_string(out, explanation.qoi.as_str());
    put_u64(out, explanation.problem_size);
    put_f64(out, explanation.budget_s);
    put_f64(out, explanation.max_relative_discrepancy);
    put_model(out, explanation.start);
    put_u32(out, explanation.regime.len() as u32);
    for (axis, coordinate) in &explanation.regime {
        put_string(out, axis.as_str());
        put_f64(out, *coordinate);
    }
    put_u8(
        out,
        match explanation.basis {
            SelectionBasis::UniqueGraphMaximum => 1,
            SelectionBasis::IncomparableMaximaCostThenIdentity => 2,
            SelectionBasis::CheapestAdequate => 3,
        },
    );
    put_u32(out, explanation.considered.len() as u32);
    for row in &explanation.considered {
        put_edge(out, row.edge);
        put_model(out, row.source);
        put_model(out, row.target);
        put_f64(out, row.source_cost_s);
        put_f64(out, row.target_cost_s);
        put_u8(
            out,
            match row.source_adequacy {
                Adequacy::Adequate => 1,
                Adequacy::Inadequate => 2,
                Adequacy::Unknown => 3,
            },
        );
        put_optional_f64(out, row.assessed_relative_discrepancy);
        put_hash(out, row.evidence.hash());
        put_u32(out, row.validity_specificity);
        put_u32(out, row.informativeness_specificity);
    }
    put_u32(out, explanation.unresolved.len() as u32);
    for edge in &explanation.unresolved {
        put_edge(out, *edge);
    }
    put_u32(out, explanation.incomparable_maxima.len() as u32);
    for model in &explanation.incomparable_maxima {
        put_model(out, *model);
    }
    put_u32(out, explanation.matched_specificity);
}

fn encode_cause(out: &mut Vec<u8>, cause: &EscalationCause) {
    match cause {
        EscalationCause::BandTooWide {
            half_width,
            tolerance,
        } => {
            put_u8(out, 1);
            put_f64(out, *half_width);
            put_f64(out, *tolerance);
        }
        EscalationCause::OutsideValidityDomain => put_u8(out, 2),
        EscalationCause::InvalidConformalPolicy => put_u8(out, 3),
        EscalationCause::ModelFormDominant => put_u8(out, 4),
    }
}

fn encode_prediction(out: &mut Vec<u8>, prediction: &PredictedPostUncertainty) {
    match prediction {
        PredictedPostUncertainty::Bounded {
            prior_half_width,
            prior_relative,
            predicted_model_form_half_width,
            predicted_total_half_width,
            predicted_total_relative,
        } => {
            put_u8(out, 1);
            for value in [
                *prior_half_width,
                *prior_relative,
                *predicted_model_form_half_width,
                *predicted_total_half_width,
                *predicted_total_relative,
            ] {
                put_f64(out, value);
            }
        }
        PredictedPostUncertainty::Indeterminate { reason } => {
            put_u8(out, 2);
            put_string(out, reason);
        }
    }
}

fn encode_probe(out: &mut Vec<u8>, suggestion: &ProbeSuggestion) {
    match suggestion {
        ProbeSuggestion::ResolveEdge {
            edge,
            source,
            target,
        } => {
            put_u8(out, 1);
            put_edge(out, *edge);
            put_model(out, *source);
            put_model(out, *target);
        }
        ProbeSuggestion::AddOutgoingEdge { source, qoi } => {
            put_u8(out, 2);
            put_model(out, *source);
            put_string(out, qoi.as_str());
        }
        ProbeSuggestion::AddModelNode { model } => {
            put_u8(out, 3);
            put_model(out, *model);
        }
    }
}

const fn advice_tag(advice: EscalationAdvice) -> u8 {
    match advice {
        EscalationAdvice::NoneNeeded => 1,
        EscalationAdvice::RefineNumerics => 2,
        EscalationAdvice::GatherMoreSamples => 3,
        EscalationAdvice::EscalateModelFidelity => 4,
    }
}

const fn gap_reason_tag(reason: EscalationGapReason) -> u8 {
    match reason {
        EscalationGapReason::AdviceDoesNotPermitModelEscalation => 1,
        EscalationGapReason::BudgetDoesNotSupportModelEscalation => 2,
        EscalationGapReason::QoiMismatch => 3,
        EscalationGapReason::ToleranceMismatch => 4,
        EscalationGapReason::InvalidReferenceMagnitude => 5,
        EscalationGapReason::CurrentModelAlreadyAdequate => 6,
        EscalationGapReason::UnknownCurrentModel => 12,
        EscalationGapReason::NoAdequateRoute => 7,
        EscalationGapReason::SessionBudgetBlocked => 8,
        EscalationGapReason::HysteresisHold => 9,
        EscalationGapReason::EscalationLimitReached => 10,
        EscalationGapReason::RepeatedEdgeBlocked => 11,
    }
}

fn put_u8(out: &mut Vec<u8>, value: u8) {
    out.push(value);
}

fn put_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&canonical_zero(value).to_bits().to_le_bytes());
}

fn put_optional_f64(out: &mut Vec<u8>, value: Option<f64>) {
    match value {
        Some(value) => {
            put_u8(out, 1);
            put_f64(out, value);
        }
        None => put_u8(out, 0),
    }
}

fn put_string(out: &mut Vec<u8>, value: &str) {
    put_u32(out, value.len() as u32);
    out.extend_from_slice(value.as_bytes());
}

fn put_hash(out: &mut Vec<u8>, value: ContentHash) {
    out.extend_from_slice(value.as_bytes());
}

fn put_optional_hash(out: &mut Vec<u8>, value: Option<ContentHash>) {
    match value {
        Some(value) => {
            put_u8(out, 1);
            put_hash(out, value);
        }
        None => put_u8(out, 0),
    }
}

fn put_edge(out: &mut Vec<u8>, value: EdgeId) {
    put_hash(out, value.hash());
}

fn put_model(out: &mut Vec<u8>, value: ModelId) {
    put_hash(out, value.hash());
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn relative_to_tolerance(width: f64, tolerance: f64) -> f64 {
    if tolerance == 0.0 {
        if width == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        width / tolerance
    }
}
