//! The cooling fidelity-graph instance battery (bead f85xj.10.4).
//!
//! Assembles the production instance from real carded model surfaces and
//! proves the three flagship demonstrations plus the instance-validation
//! battery: cards resolve, edges stay fresh, the legacy `cht()` embedding
//! is intact, and the structural digest pins the graph shape under the
//! golden-bump discipline.

use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ladder::{
    EdgeEvidenceResolver, FidelityEdge, LadderRegistry, QoiId, QueryContext, QueryEvidenceRef,
    QueryRefusal, RegimeAxis, ResolvedEdgeEvidence, SelectionBasis,
};
use fs_ledger::Ledger;
use fs_plan::{
    COOLING_CARD_DOMAIN, COOLING_INSTANCE_NAME, CoolingFidelityInstance,
    RADIATION_REGIME_SPLIT_TOLERANCE, assemble_cooling_instance, card_ref_for,
    cooling_instance::labels, record_fidelity_campaign,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C001_0104_0000,
                kernel_id: 104,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn instance() -> CoolingFidelityInstance {
    with_cx(|cx| assemble_cooling_instance(cx).expect("the production instance assembles"))
}

/// Resolver backed by the instance's own fitted evidence: it echoes the
/// exact cost/discrepancy references stored on each edge together with the
/// retained numeric summary, which is precisely the adapter contract a
/// product consumer implements.
struct InstanceResolver<'a> {
    instance: &'a CoolingFidelityInstance,
}

impl EdgeEvidenceResolver for InstanceResolver<'_> {
    fn resolve(
        &self,
        edge: &FidelityEdge,
        _context: &QueryContext,
    ) -> Option<ResolvedEdgeEvidence> {
        let summary = self
            .instance
            .summaries
            .iter()
            .find(|summary| summary.edge == edge.id())?;
        ResolvedEdgeEvidence::new(
            Some(summary.cost_ref),
            Some(summary.discrepancy_ref),
            summary.mean_source_cost_s,
            summary.mean_target_cost_s,
            Some(summary.fitted_max_rel),
            QueryEvidenceRef::new(hash_domain(
                "test.fs-plan.cooling-instance.query-receipt.v1",
                summary.edge.to_string().as_bytes(),
            )),
        )
        .ok()
    }
}

fn model(instance: &CoolingFidelityInstance, label: &str) -> fs_ladder::ModelId {
    instance
        .cards
        .iter()
        .find(|card| card.label == label)
        .expect("every named node is carded")
        .model
}

fn radiation_query(surface_temperature_k: f64) -> QueryContext {
    QueryContext::new(
        QoiId::new("outward-radiative-flux").expect("qoi"),
        [
            (
                RegimeAxis::new("surface_temperature_k").expect("axis"),
                surface_temperature_k,
            ),
            (
                RegimeAxis::new("ambient_temperature_k").expect("axis"),
                300.0,
            ),
            (RegimeAxis::new("emissivity").expect("axis"), 0.8),
        ],
        1,
        1.0,
        RADIATION_REGIME_SPLIT_TOLERANCE,
    )
    .expect("query context")
}

#[test]
fn instance_assembles_with_carded_nodes_and_evidenced_edges() {
    let instance = instance();

    assert!(
        instance.cards.len() >= 5,
        "DONE-WHEN floor: at least 5 carded nodes, got {}",
        instance.cards.len()
    );
    assert_eq!(instance.campaign.graph.nodes().len(), instance.cards.len());
    assert!(
        instance.campaign.edges.len() >= 3,
        "DONE-WHEN floor: at least 3 evidenced edges, got {}",
        instance.campaign.edges.len()
    );
    assert_eq!(instance.campaign.gaps.len(), 3);
    assert_eq!(instance.summaries.len(), instance.campaign.edges.len());

    // Cards resolve: every node's stored reference is recomputable from the
    // retained card's canonical bytes under the published domain.
    for card in &instance.cards {
        assert_eq!(
            card.card_ref,
            card_ref_for(&card.card),
            "card binding for {} must be recomputable",
            card.label
        );
        let node = instance
            .campaign
            .graph
            .node(card.model)
            .expect("every retained card names a graph node");
        assert_eq!(node.card(), card.card_ref);
        // The domain separation is part of the promise.
        assert_eq!(
            card.card_ref.as_bytes(),
            hash_domain(
                COOLING_CARD_DOMAIN,
                card.card.to_ledger_row_json().as_bytes()
            )
            .as_bytes()
        );
    }

    // Both radiation edges carry an independent closed-form reference, so
    // informativeness is evidenced there and honestly unknown elsewhere.
    let radiation_target = model(&instance, labels::RADIATION_GRAY_DIFFUSE);
    for fitted in &instance.campaign.edges {
        if fitted.edge.target() == radiation_target {
            assert!(fitted.informativeness_supported);
        } else {
            assert!(!fitted.informativeness_supported);
            assert!(fitted.edge.informativeness().predicates().is_unknown());
        }
    }
}

#[test]
fn structural_digest_is_cost_independent_and_pinned() {
    let first = instance();
    let second = instance();
    // Wall-clock probe costs differ between assemblies; the structural
    // digest must not see them, while the full graph identity legitimately
    // does (it binds the fitted cost/discrepancy artifacts).
    assert_eq!(first.structural_digest, second.structural_digest);

    // Golden-bump discipline: this constant pins the instance SHAPE (node
    // set and bindings, edge contexts, gap list). Bumping it requires a
    // semantic reason in the commit that changes it, per the repository's
    // golden protocol.
    let digest = first.structural_digest.to_string();
    assert_eq!(
        digest, GOLDEN_STRUCTURAL_DIGEST,
        "cooling-instance shape changed; if intentional, bump the golden with cause"
    );
}

/// See `structural_digest_is_cost_independent_and_pinned` before touching.
/// Frozen 2026-07-27 from the first assembled instance (7 nodes, 4 fitted
/// edge contexts, 3 gaps); any change to that shape must re-derive this in
/// the same commit with the semantic cause stated.
const GOLDEN_STRUCTURAL_DIGEST: &str =
    "d83b94768948dbf0f124f650aa7b6301da02fb9c97cb2123f38924af6b6b6713";

#[test]
fn demo_cost_is_not_authority_in_the_validated_regime() {
    let instance = instance();
    let resolver = InstanceResolver {
        instance: &instance,
    };
    let start = model(&instance, labels::RADIATION_LINEARIZED);
    let context = radiation_query(324.0);

    // The graph maximum is the gray-diffuse rung...
    let best = instance
        .campaign
        .graph
        .best_model_for(start, &context, &resolver)
        .expect("evidence exists near the anchor");
    assert_eq!(best.model, model(&instance, labels::RADIATION_GRAY_DIFFUSE));

    // ...but inside the near-anchor validated regime the calibrated cheap
    // rung is ADEQUATE at the demo tolerance, and it is selected precisely
    // because its fitted discrepancy band is tight where it was validated:
    // cost never created that authority, evidence did.
    let cheapest = instance
        .campaign
        .graph
        .cheapest_adequate(start, &context, &resolver)
        .expect("the calibrated cheap rung is adequate near the anchor");
    assert_eq!(cheapest.model, start);
    assert_eq!(cheapest.explanation.basis, SelectionBasis::CheapestAdequate);
    let near_summary = instance
        .summaries
        .iter()
        .find(|summary| {
            summary.source == start && summary.fitted_max_rel < RADIATION_REGIME_SPLIT_TOLERANCE
        })
        .expect("the near-anchor edge is fitted");
    assert!(
        near_summary.mean_source_cost_s <= near_summary.mean_target_cost_s * 10.0,
        "sanity: the demo does not depend on any particular cost ratio"
    );
    assert!(
        !cheapest.explanation.considered.is_empty(),
        "the recommendation must carry its reasoning"
    );
}

#[test]
fn demo_same_qoi_routes_to_different_models_by_regime() {
    let instance = instance();
    let resolver = InstanceResolver {
        instance: &instance,
    };
    let start = model(&instance, labels::RADIATION_LINEARIZED);
    let fine = model(&instance, labels::RADIATION_GRAY_DIFFUSE);

    // Near the linearization anchor the cheap rung's validated band is
    // inside tolerance: the query keeps the cheap model.
    let near = instance
        .campaign
        .graph
        .cheapest_adequate(start, &radiation_query(324.0), &resolver)
        .expect("near-anchor adequacy");
    assert_eq!(near.model, start);

    // Far below the anchor the SAME QoI under the SAME tolerance is out of
    // the cheap rung's evidenced adequacy (fitted band ≈ 15–21%), so the
    // stand-in query REFUSES — cheapest_adequate never lets an inadequate
    // source impersonate the maximum — and the certify-or-escalate policy
    // routes the regime to the gray-diffuse maximum instead.
    let far_context = radiation_query(306.0);
    let refusal = instance
        .campaign
        .graph
        .cheapest_adequate(start, &far_context, &resolver)
        .expect_err("the cheap rung's far-cold band exceeds the tolerance");
    assert!(matches!(refusal, QueryRefusal::NoAdequateModel { .. }));
    let far = instance
        .campaign
        .graph
        .best_model_for(start, &far_context, &resolver)
        .expect("the fine rung remains the evidenced maximum far from the anchor");
    assert_eq!(far.model, fine);
    assert_ne!(
        near.model, far.model,
        "same QoI, different regime, different model"
    );
}

#[test]
fn demo_honest_gap_is_a_refusal_plus_acquisition_demand() {
    let instance = instance();
    let resolver = InstanceResolver {
        instance: &instance,
    };
    // Natural convection: the node exists and is carded, but no edge
    // touches it, so a query from it must refuse rather than improvise.
    let start = model(&instance, labels::NATURAL_CHURCHILL_CHU);
    let context = QueryContext::new(
        QoiId::new("nusselt-number").expect("qoi"),
        [
            (RegimeAxis::new("Ra").expect("axis"), 1.0e7),
            (RegimeAxis::new("Pr").expect("axis"), 0.7),
        ],
        1,
        1.0,
        0.2,
    )
    .expect("query context");
    let refusal = instance
        .campaign
        .graph
        .best_model_for(start, &context, &resolver)
        .expect_err("no evidenced comparison touches the natural-convection node");
    assert_eq!(
        refusal,
        QueryRefusal::NoApplicableEvidence {
            start,
            qoi: QoiId::new("nusselt-number").expect("qoi"),
        }
    );

    // The refusal is paired with an explicit, actionable acquisition
    // demand retained in the campaign itself.
    let gap = instance
        .campaign
        .gaps
        .iter()
        .find(|gap| gap.source == start)
        .expect("the natural-convection gap is recorded");
    assert!(
        gap.reason.contains("no shared-validity geometry"),
        "the gap names what evidence must be acquired: {}",
        gap.reason
    );
}

#[test]
fn freshness_and_ledger_retention_hold() {
    let instance = instance();

    // Fresh against its own authority; stale the moment the corpus moves.
    let authority = instance.campaign.authority.clone();
    assert!(!instance.campaign.assess_freshness(&authority).is_stale());
    let mut moved = authority;
    moved.corpus_version += 1;
    assert!(instance.campaign.assess_freshness(&moved).is_stale());

    // Atomic ledger retention of graph, campaign, and per-edge artifacts.
    let ledger = Ledger::open(":memory:").expect("ledger");
    let receipt = record_fidelity_campaign(&ledger, &instance.campaign, 1_000, 2_000)
        .expect("atomic retention");
    assert_eq!(receipt.edge_artifacts.len(), instance.campaign.edges.len());
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn legacy_cht_embedding_remains_intact() {
    // The instance-validation battery keeps the legacy ladder path alive:
    // the Proposal-7 CHT declaration still registers, resolves, and
    // migrates into the graph model while this instance exists beside it.
    let registry = LadderRegistry::cht();
    let ladder = registry.ladder("cht").expect("cht ladder resolves");
    assert_eq!(ladder.top().index, 2);
    let embedded = LadderRegistry::cht_graph().expect("cht graph migration");
    assert!(!embedded.graph().nodes().is_empty());

    let instance = instance();
    assert_eq!(
        instance.campaign.graph.name(),
        COOLING_INSTANCE_NAME,
        "the production instance is its own graph, not a cht mutation"
    );
}

#[test]
fn probe_metadata_stays_joinable() {
    let instance = instance();
    // Every summary joins a fitted edge and carries the exact references
    // stored on that edge; a resolver echoing them is accepted, which is
    // what the query demos rely on.
    let mut seen = BTreeMap::new();
    for summary in &instance.summaries {
        let fitted = instance
            .campaign
            .edges
            .iter()
            .find(|fitted| fitted.edge.id() == summary.edge)
            .expect("summary joins a fitted edge");
        assert_eq!(fitted.edge.source(), summary.source);
        assert_eq!(fitted.edge.target(), summary.target);
        assert!(summary.fitted_max_rel.is_finite());
        assert!(summary.mean_source_cost_s > 0.0);
        assert!(summary.mean_target_cost_s > 0.0);
        seen.insert(summary.edge, ());
    }
    assert_eq!(seen.len(), instance.summaries.len(), "edge ids are unique");
}
