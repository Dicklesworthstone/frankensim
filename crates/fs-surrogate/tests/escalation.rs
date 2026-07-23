//! G0/G3 graph-aware certify-or-escalate battery.

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::EscalationAdvice;
use fs_evidence::uncertainty::{
    EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
    TermValue, UncertaintyArtifactRef,
};
use fs_ladder::{
    ContextPredicateSet, CostModelRef, CostRelationRef, DiscrepancyModelRef, DiscrepancyReference,
    EdgeEvidenceResolver, EdgeId, FidelityEdge, FidelityGraph, FidelityNode, Informativeness,
    ModelCardRef, ModelId, QoiId, QueryContext, QueryEvidenceRef, RegimeAxis, ResolvedEdgeEvidence,
    TransferRef, ValidityDomain,
};
use fs_surrogate::ConformalBand;
use fs_surrogate::escalation::{
    EscalationCalibrationRecord, EscalationGapReason, EscalationPolicy, EscalationRunEvidenceRef,
    EscalationSession, GraphDecision, PredictedPostUncertainty, ProbeSuggestion,
    certify_or_escalate_with_graph, plan_graph_escalation,
};
use std::collections::BTreeMap;

fn hash(domain: &str, label: &str) -> ContentHash {
    hash_domain(domain, label.as_bytes())
}

fn assert_near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-12,
        "expected {expected:.17e}, got {actual:.17e}"
    );
}

fn model(label: &str) -> ModelId {
    ModelId::new(hash("test.fs-surrogate.model.v1", label))
}

fn card(label: &str) -> ModelCardRef {
    ModelCardRef::new(hash("test.fs-surrogate.card.v1", label))
}

fn cost(label: &str) -> CostModelRef {
    CostModelRef::new(hash("test.fs-surrogate.cost.v1", label))
}

fn discrepancy(label: &str) -> DiscrepancyModelRef {
    DiscrepancyModelRef::new(hash("test.fs-surrogate.discrepancy.v1", label))
}

fn transfer(label: &str) -> TransferRef {
    TransferRef::new(hash("test.fs-surrogate.transfer.v1", label))
}

fn receipt(label: &str) -> QueryEvidenceRef {
    QueryEvidenceRef::new(hash("test.fs-surrogate.query.v1", label))
}

fn add_node(graph: &mut FidelityGraph, label: &str) -> ModelId {
    let id = model(label);
    graph
        .add_node(FidelityNode::new(id, card(label), label).unwrap())
        .unwrap();
    id
}

fn edge(source: ModelId, target: ModelId, label: &str) -> FidelityEdge {
    FidelityEdge::new(
        source,
        target,
        cost(label),
        discrepancy(label),
        transfer(label),
        ValidityDomain::universal(),
        Informativeness::new(ContextPredicateSet::universal()),
    )
    .unwrap()
}

#[derive(Default)]
struct Resolver {
    rows: BTreeMap<EdgeId, ResolvedEdgeEvidence>,
}

impl Resolver {
    fn insert(
        &mut self,
        edge: &FidelityEdge,
        source_cost_s: f64,
        target_cost_s: f64,
        assessed_relative_discrepancy: Option<f64>,
        label: &str,
    ) {
        let cost_model = match edge.cost() {
            CostRelationRef::Model(reference) => Some(reference),
            CostRelationRef::LegacyRelativeCost(_) => None,
        };
        let discrepancy_model = match edge.discrepancy() {
            DiscrepancyReference::Model(reference) => Some(reference),
            DiscrepancyReference::UnknownLegacy => None,
        };
        self.rows.insert(
            edge.id(),
            ResolvedEdgeEvidence::new(
                cost_model,
                discrepancy_model,
                source_cost_s,
                target_cost_s,
                assessed_relative_discrepancy,
                receipt(label),
            )
            .unwrap(),
        );
    }
}

impl EdgeEvidenceResolver for Resolver {
    fn resolve(
        &self,
        edge: &FidelityEdge,
        _context: &QueryContext,
    ) -> Option<ResolvedEdgeEvidence> {
        self.rows.get(&edge.id()).copied()
    }
}

struct Fixture {
    graph: FidelityGraph,
    resolver: Resolver,
    correlation: ModelId,
    rans: ModelId,
}

fn fixture() -> Fixture {
    let mut graph = FidelityGraph::new("recirculation").unwrap();
    let correlation = add_node(&mut graph, "correlation");
    let rans = add_node(&mut graph, "rans");
    let les = add_node(&mut graph, "les");
    let correlation_rans = edge(correlation, rans, "correlation-rans");
    let rans_les = edge(rans, les, "rans-les");
    graph.add_edge(correlation_rans.clone()).unwrap();
    graph.add_edge(rans_les.clone()).unwrap();
    let mut resolver = Resolver::default();
    resolver.insert(&correlation_rans, 1.0, 40.0, Some(0.20), "corr-rans");
    resolver.insert(&rans_les, 40.0, 400.0, Some(0.02), "rans-les");
    Fixture {
        graph,
        resolver,
        correlation,
        rans,
    }
}

fn context(budget_s: f64, tolerance: f64) -> QueryContext {
    QueryContext::new(
        QoiId::new("recirculation.drag").unwrap(),
        [(RegimeAxis::new("reynolds").unwrap(), 50_000.0)],
        1_000_000,
        budget_s,
        tolerance,
    )
    .unwrap()
}

fn budget(model_form_half_width: f64) -> EngineeringUncertaintyBudget {
    let terms = EngineeringUncertaintyKind::ALL
        .into_iter()
        .map(|kind| {
            let half_width = if kind == EngineeringUncertaintyKind::ModelForm {
                model_form_half_width
            } else {
                1.0
            };
            EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::interval(half_width, half_width).unwrap(),
                UncertaintyArtifactRef::new(
                    &format!("test:{}", kind.name()),
                    hash("test.fs-surrogate.uncertainty.v1", kind.name()),
                )
                .unwrap(),
            )
            .unwrap()
        })
        .collect();
    EngineeringUncertaintyBudget::try_new("recirculation.drag", "newton", terms).unwrap()
}

fn wide_band() -> ConformalBand {
    ConformalBand {
        half_width: 20.0,
        alpha: 0.1,
    }
}

#[test]
fn helpful_edge_routes_to_the_cheapest_adequate_model() {
    let fixture = fixture();
    let budget = budget(20.0);
    let plan = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget,
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap();
    assert_eq!(plan.recommendation().model, fixture.rans);
    assert_eq!(plan.steps().len(), 1);
    assert_eq!(plan.next_edge(), plan.steps()[0].edge);
    assert_eq!(plan.steps()[0].source, fixture.correlation);
    assert_eq!(plan.steps()[0].target, fixture.rans);
    assert_eq!(
        plan.recommendation().predicted_cost_s.to_bits(),
        40.0f64.to_bits()
    );
    match plan.predicted_uncertainty() {
        PredictedPostUncertainty::Bounded {
            prior_half_width,
            predicted_model_form_half_width,
            predicted_total_half_width,
            predicted_total_relative,
            ..
        } => {
            // The eight-term budget deliberately aggregates upward, so the
            // total can sit a few ulps above the mathematical sum.
            assert_near(*prior_half_width, 27.0);
            assert_near(*predicted_model_form_half_width, 2.0);
            assert_near(*predicted_total_half_width, 9.0);
            assert_near(*predicted_total_relative, 0.09);
        }
        PredictedPostUncertainty::Indeterminate { reason } => {
            panic!("unexpected indeterminate prediction: {reason}")
        }
    }
    assert_eq!(plan.canonical_bytes(), plan.canonical_bytes());
    assert_eq!(plan.content_id(), plan.content_id());
}

#[test]
fn graph_aware_policy_preserves_use_surrogate_and_plans_escalation() {
    let fixture = fixture();
    let budget = budget(20.0);
    let query = context(100.0, 0.05);
    let mut session = EscalationSession::new(EscalationPolicy::default());
    let narrow = ConformalBand {
        half_width: 1.0,
        alpha: 0.1,
    };
    assert!(matches!(
        certify_or_escalate_with_graph(
            &narrow,
            true,
            5.0,
            &fixture.graph,
            &fixture.resolver,
            fixture.correlation,
            &query,
            &budget,
            EscalationAdvice::EscalateModelFidelity,
            100.0,
            &mut session,
        ),
        GraphDecision::UseSurrogate {
            band_half_width: 1.0
        }
    ));
    let decision = certify_or_escalate_with_graph(
        &wide_band(),
        true,
        5.0,
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &query,
        &budget,
        EscalationAdvice::EscalateModelFidelity,
        100.0,
        &mut session,
    );
    assert!(matches!(decision, GraphDecision::Escalate { .. }));
    assert_eq!(session.consumed_edges(), 1);
    assert!(session.is_latched());
}

#[test]
fn hysteresis_holds_repeated_noise_until_a_real_recovery() {
    let fixture = fixture();
    let budget = budget(20.0);
    let query = context(100.0, 0.05);
    let mut session = EscalationSession::new(EscalationPolicy::default());
    assert!(matches!(
        certify_or_escalate_with_graph(
            &wide_band(),
            true,
            5.0,
            &fixture.graph,
            &fixture.resolver,
            fixture.correlation,
            &query,
            &budget,
            EscalationAdvice::EscalateModelFidelity,
            100.0,
            &mut session,
        ),
        GraphDecision::Escalate { .. }
    ));
    let held = certify_or_escalate_with_graph(
        &wide_band(),
        true,
        5.0,
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &query,
        &budget,
        EscalationAdvice::EscalateModelFidelity,
        100.0,
        &mut session,
    );
    assert!(matches!(
        held,
        GraphDecision::Gap { ref gap }
            if gap.reason() == EscalationGapReason::HysteresisHold
    ));

    let recovered = ConformalBand {
        half_width: 1.0,
        alpha: 0.1,
    };
    assert!(matches!(
        certify_or_escalate_with_graph(
            &recovered,
            true,
            5.0,
            &fixture.graph,
            &fixture.resolver,
            fixture.correlation,
            &query,
            &budget,
            EscalationAdvice::EscalateModelFidelity,
            100.0,
            &mut session,
        ),
        GraphDecision::UseSurrogate { .. }
    ));
    assert!(!session.is_latched());
    let repeated = certify_or_escalate_with_graph(
        &wide_band(),
        true,
        5.0,
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &query,
        &budget,
        EscalationAdvice::EscalateModelFidelity,
        100.0,
        &mut session,
    );
    assert!(matches!(
        repeated,
        GraphDecision::Gap { ref gap }
            if gap.reason() == EscalationGapReason::RepeatedEdgeBlocked
    ));
}

#[test]
fn no_edge_returns_an_honest_probe_gap() {
    let mut graph = FidelityGraph::new("isolated").unwrap();
    let current = add_node(&mut graph, "correlation");
    let gap = plan_graph_escalation(
        &graph,
        &Resolver::default(),
        current,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap_err();
    assert_eq!(gap.reason(), EscalationGapReason::NoAdequateRoute);
    assert!(gap.candidate_edges().is_empty());
    assert!(matches!(
        gap.suggested_probes(),
        [ProbeSuggestion::AddOutgoingEdge { source, .. }] if *source == current
    ));
    assert_eq!(gap.content_id(), gap.content_id());
}

#[test]
fn unknown_current_model_requests_graph_registration() {
    let fixture = fixture();
    let missing = model("unregistered");
    let gap = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        missing,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap_err();
    assert_eq!(gap.reason(), EscalationGapReason::UnknownCurrentModel);
    assert!(matches!(
        gap.suggested_probes(),
        [ProbeSuggestion::AddModelNode { model }] if *model == missing
    ));
}

#[test]
fn session_budget_block_names_the_required_cost() {
    let fixture = fixture();
    let gap = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(10.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap_err();
    assert_eq!(gap.reason(), EscalationGapReason::SessionBudgetBlocked);
    assert_eq!(gap.required_cost_s().unwrap().to_bits(), 40.0f64.to_bits());
}

#[test]
fn non_model_form_advice_and_budget_fail_closed() {
    let fixture = fixture();
    let wrong_advice = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::GatherMoreSamples,
        100.0,
    )
    .unwrap_err();
    assert_eq!(
        wrong_advice.reason(),
        EscalationGapReason::AdviceDoesNotPermitModelEscalation
    );
    let different_wrong_advice = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::RefineNumerics,
        100.0,
    )
    .unwrap_err();
    assert_ne!(
        wrong_advice.content_id(),
        different_wrong_advice.content_id()
    );

    let numerics_dominate = budget(0.5);
    let wrong_budget = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &numerics_dominate,
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap_err();
    assert_eq!(
        wrong_budget.reason(),
        EscalationGapReason::BudgetDoesNotSupportModelEscalation
    );

    let invalid_reference = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        -100.0,
    )
    .unwrap_err();
    assert_eq!(
        invalid_reference.reason(),
        EscalationGapReason::InvalidReferenceMagnitude
    );
    let zero_reference = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        0.0,
    )
    .unwrap_err();
    assert_ne!(invalid_reference.content_id(), zero_reference.content_id());
}

#[test]
fn per_run_edge_cap_refuses_a_multi_step_route() {
    let mut graph = FidelityGraph::new("deep").unwrap();
    let correlation = add_node(&mut graph, "correlation");
    let rans = add_node(&mut graph, "rans");
    let les = add_node(&mut graph, "les");
    let dns = add_node(&mut graph, "dns");
    let first = edge(correlation, rans, "correlation-rans");
    let second = edge(rans, les, "rans-les");
    let third = edge(les, dns, "les-dns");
    for value in [&first, &second, &third] {
        graph.add_edge(value.clone()).unwrap();
    }
    let mut resolver = Resolver::default();
    resolver.insert(&first, 1.0, 10.0, Some(0.20), "first");
    resolver.insert(&second, 10.0, 40.0, Some(0.10), "second");
    resolver.insert(&third, 40.0, 400.0, Some(0.02), "third");
    let policy = EscalationPolicy::try_new(1.05, 0.90, 1).unwrap();
    let mut session = EscalationSession::new(policy);
    let decision = certify_or_escalate_with_graph(
        &wide_band(),
        true,
        5.0,
        &graph,
        &resolver,
        correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
        &mut session,
    );
    assert!(matches!(
        decision,
        GraphDecision::Gap { ref gap }
            if gap.reason() == EscalationGapReason::EscalationLimitReached
    ));
    assert_eq!(session.consumed_edges(), 0);
}

#[test]
fn calibration_record_binds_prediction_actuals_and_run_evidence() {
    let fixture = fixture();
    let plan = plan_graph_escalation(
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
    )
    .unwrap();
    let actual = budget(1.0);
    let record = EscalationCalibrationRecord::try_new(
        &plan,
        &actual,
        45.0,
        EscalationRunEvidenceRef::new(hash("test.fs-surrogate.run-evidence.v1", "rans-run")),
    )
    .unwrap();
    assert_eq!(record.cost_error_s().to_bits(), 5.0f64.to_bits());
    assert_near(record.predicted_improvement().unwrap(), 0.18);
    assert_near(record.actual_improvement().unwrap(), 0.19);
    assert_near(record.uncertainty_error().unwrap(), -0.01);
    assert_eq!(record.content_id(), record.content_id());
    let changed_run = EscalationCalibrationRecord::try_new(
        &plan,
        &actual,
        45.0,
        EscalationRunEvidenceRef::new(hash("test.fs-surrogate.run-evidence.v1", "different-run")),
    )
    .unwrap();
    assert_ne!(record.content_id(), changed_run.content_id());
}

#[test]
fn policy_rejects_degenerate_hysteresis_and_zero_caps() {
    assert!(EscalationPolicy::try_new(0.9, 0.8, 1).is_err());
    assert!(EscalationPolicy::try_new(1.0, 1.0, 1).is_err());
    assert!(EscalationPolicy::try_new(1.1, 0.9, 0).is_err());
}

#[test]
fn context_tolerance_mismatch_cannot_route() {
    let fixture = fixture();
    let mut session = EscalationSession::new(EscalationPolicy::default());
    let decision = certify_or_escalate_with_graph(
        &wide_band(),
        true,
        4.0,
        &fixture.graph,
        &fixture.resolver,
        fixture.correlation,
        &context(100.0, 0.05),
        &budget(20.0),
        EscalationAdvice::EscalateModelFidelity,
        100.0,
        &mut session,
    );
    assert!(matches!(
        decision,
        GraphDecision::Gap { ref gap }
            if gap.reason() == EscalationGapReason::ToleranceMismatch
                && gap.candidate_plan().is_some()
    ));
}
