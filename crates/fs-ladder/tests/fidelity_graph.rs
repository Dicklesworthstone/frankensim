//! G0/G3 battery for the context- and QoI-specific fidelity graph.

use fs_blake3::{ContentHash, hash_domain};
use fs_ladder::{
    Adequacy, ClosedInterval, ContextClause, ContextPredicateSet, CostModelRef,
    EdgeEvidenceResolver, EdgeId, FIDELITY_GRAPH_SCHEMA_VERSION, FidelityEdge, FidelityGraph,
    FidelityNode, GraphError, Informativeness, Ladder, LadderRegistry, ModelCardRef, ModelId,
    QoiId, QoiSelector, QueryContext, QueryEvidenceRef, QueryRefusal, Refine1d, RegimeAxis,
    ResolvedEdgeEvidence, SelectionBasis, TransferRef, ValidityDomain,
};
use std::collections::BTreeMap;

fn hash(domain: &str, label: &str) -> ContentHash {
    hash_domain(domain, label.as_bytes())
}

fn model(label: &str) -> ModelId {
    ModelId::new(hash("test.fs-ladder.model.v1", label))
}

fn card(label: &str) -> ModelCardRef {
    ModelCardRef::new(hash("test.fs-ladder.card.v1", label))
}

fn cost(label: &str) -> CostModelRef {
    CostModelRef::new(hash("test.fs-ladder.cost.v1", label))
}

fn discrepancy(label: &str) -> fs_ladder::DiscrepancyModelRef {
    fs_ladder::DiscrepancyModelRef::new(hash("test.fs-ladder.discrepancy.v1", label))
}

fn transfer(label: &str) -> TransferRef {
    TransferRef::new(hash("test.fs-ladder.transfer.v1", label))
}

fn query_receipt(label: &str) -> QueryEvidenceRef {
    QueryEvidenceRef::new(hash("test.fs-ladder.query-evidence.v1", label))
}

fn exact_context(qoi: &str, axis: &str, lower: f64, upper: f64) -> ContextPredicateSet {
    ContextPredicateSet::new([ContextClause::new(
        QoiSelector::Exact(QoiId::new(qoi).unwrap()),
        [(
            RegimeAxis::new(axis).unwrap(),
            ClosedInterval::new(lower, upper).unwrap(),
        )],
    )
    .unwrap()])
    .unwrap()
}

fn native_edge(
    source: ModelId,
    target: ModelId,
    label: &str,
    predicates: ContextPredicateSet,
) -> FidelityEdge {
    FidelityEdge::new(
        source,
        target,
        cost(label),
        discrepancy(label),
        transfer(label),
        ValidityDomain::new(predicates.clone()),
        Informativeness::new(predicates),
    )
    .unwrap()
}

fn add_node(graph: &mut FidelityGraph, label: &str) -> ModelId {
    let id = model(label);
    graph
        .add_node(FidelityNode::new(id, card(label), label).unwrap())
        .unwrap();
    id
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
        adequacy: Adequacy,
        receipt: &str,
    ) {
        let cost_model = match edge.cost() {
            fs_ladder::CostRelationRef::Model(reference) => Some(reference),
            fs_ladder::CostRelationRef::LegacyRelativeCost(_) => None,
        };
        let discrepancy_model = match edge.discrepancy() {
            fs_ladder::DiscrepancyReference::Model(reference) => Some(reference),
            fs_ladder::DiscrepancyReference::UnknownLegacy => None,
        };
        self.rows.insert(
            edge.id(),
            ResolvedEdgeEvidence::new(
                cost_model,
                discrepancy_model,
                source_cost_s,
                target_cost_s,
                match adequacy {
                    Adequacy::Adequate => Some(0.01),
                    Adequacy::Inadequate => Some(0.10),
                    Adequacy::Unknown => None,
                },
                query_receipt(receipt),
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

struct Synthetic {
    graph: FidelityGraph,
    resolver: Resolver,
    correlation: ModelId,
    rans: ModelId,
    les: ModelId,
    calibrated: ModelId,
}

fn synthetic_graph() -> Synthetic {
    let mut graph = FidelityGraph::new("cooling.synthetic").unwrap();
    let correlation = add_node(&mut graph, "correlation");
    let rans = add_node(&mut graph, "rans");
    let les = add_node(&mut graph, "les");
    let calibrated = add_node(&mut graph, "calibrated");

    let mean = exact_context("temperature.mean", "reynolds", 1_000.0, 100_000.0);
    let maximum = exact_context("temperature.max", "reynolds", 1_000.0, 100_000.0);
    let corr_to_rans_mean = native_edge(correlation, rans, "corr-rans-mean", mean.clone());
    let rans_to_calibrated = native_edge(rans, calibrated, "rans-calibrated-mean", mean);
    let corr_to_rans_max = native_edge(correlation, rans, "corr-rans-max", maximum.clone());
    let rans_to_les = native_edge(rans, les, "rans-les-max", maximum);

    for edge in [
        corr_to_rans_mean.clone(),
        rans_to_calibrated.clone(),
        corr_to_rans_max.clone(),
        rans_to_les.clone(),
    ] {
        graph.add_edge(edge).unwrap();
    }

    let mut resolver = Resolver::default();
    resolver.insert(
        &corr_to_rans_mean,
        1.0,
        40.0,
        Adequacy::Adequate,
        "corr-rans-mean",
    );
    // Cost and epistemic direction deliberately disagree: the calibrated
    // correlation is more informative here while much cheaper than RANS.
    resolver.insert(
        &rans_to_calibrated,
        40.0,
        2.0,
        Adequacy::Inadequate,
        "rans-calibrated-mean",
    );
    resolver.insert(
        &corr_to_rans_max,
        1.0,
        40.0,
        Adequacy::Inadequate,
        "corr-rans-max",
    );
    resolver.insert(
        &rans_to_les,
        40.0,
        400.0,
        Adequacy::Adequate,
        "rans-les-max",
    );
    Synthetic {
        graph,
        resolver,
        correlation,
        rans,
        les,
        calibrated,
    }
}

fn context(qoi: &str, reynolds: f64, budget_s: f64) -> QueryContext {
    QueryContext::new(
        QoiId::new(qoi).unwrap(),
        [(RegimeAxis::new("reynolds").unwrap(), reynolds)],
        1_000_000,
        budget_s,
        0.05,
    )
    .unwrap()
}

#[test]
fn g0_graph_refuses_self_loops_missing_cards_and_duplicates() {
    let a = model("a");
    let predicates = ContextPredicateSet::universal();
    assert!(matches!(
        FidelityEdge::new(
            a,
            a,
            cost("self"),
            discrepancy("self"),
            transfer("self"),
            ValidityDomain::new(predicates.clone()),
            Informativeness::new(predicates.clone()),
        ),
        Err(GraphError::SelfLoop(found)) if found == a
    ));

    let mut graph = FidelityGraph::new("validation").unwrap();
    graph
        .add_node(FidelityNode::new(a, card("a"), "a").unwrap())
        .unwrap();
    assert!(matches!(
        graph.add_node(FidelityNode::new(a, card("a2"), "a2").unwrap()),
        Err(GraphError::DuplicateNode(found)) if found == a
    ));
    let b = model("b");
    let edge = native_edge(a, b, "missing-b", predicates);
    assert!(matches!(
        graph.add_edge(edge),
        Err(GraphError::MissingNode(found)) if found == b
    ));
}

#[test]
fn g3_clause_and_graph_identity_ignore_declaration_order() {
    let qoi = QoiId::new("temperature.mean").unwrap();
    let reynolds = RegimeAxis::new("reynolds").unwrap();
    let mach = RegimeAxis::new("mach").unwrap();
    let clause_a = ContextClause::new(
        QoiSelector::Exact(qoi.clone()),
        [
            (reynolds.clone(), ClosedInterval::new(1.0, 2.0).unwrap()),
            (mach.clone(), ClosedInterval::new(0.0, 0.3).unwrap()),
        ],
    )
    .unwrap();
    let clause_b = ContextClause::new(
        QoiSelector::Exact(qoi),
        [
            (mach, ClosedInterval::new(0.0, 0.3).unwrap()),
            (reynolds, ClosedInterval::new(1.0, 2.0).unwrap()),
        ],
    )
    .unwrap();
    let broad = ContextClause::universal();
    let predicates_a = ContextPredicateSet::new([clause_a, broad.clone(), broad.clone()]).unwrap();
    let predicates_b = ContextPredicateSet::new([broad, clause_b]).unwrap();
    assert_eq!(predicates_a, predicates_b);

    let a = model("a");
    let b = model("b");
    let node_a = FidelityNode::new(a, card("a"), "a").unwrap();
    let node_b = FidelityNode::new(b, card("b"), "b").unwrap();
    let edge_a = native_edge(a, b, "ab", predicates_a);
    let edge_b = native_edge(a, b, "ab", predicates_b);
    assert_eq!(edge_a.id(), edge_b.id());

    let mut first = FidelityGraph::new("order").unwrap();
    first.add_node(node_a.clone()).unwrap();
    first.add_node(node_b.clone()).unwrap();
    first.add_edge(edge_a).unwrap();
    let mut second = FidelityGraph::new("order").unwrap();
    second.add_node(node_b).unwrap();
    second.add_node(node_a).unwrap();
    second.add_edge(edge_b).unwrap();
    assert_eq!(first.canonical_bytes(), second.canonical_bytes());
    assert_eq!(first.identity(), second.identity());
}

#[test]
fn g3_signed_zero_is_canonical_in_domains_and_queries() {
    assert_eq!(
        ClosedInterval::new(-0.0, 0.0).unwrap(),
        ClosedInterval::new(0.0, 0.0).unwrap()
    );
    let negative = QueryContext::new(
        QoiId::new("temperature.mean").unwrap(),
        [(RegimeAxis::new("offset").unwrap(), -0.0)],
        1,
        -0.0,
        -0.0,
    )
    .unwrap();
    let positive = QueryContext::new(
        QoiId::new("temperature.mean").unwrap(),
        [(RegimeAxis::new("offset").unwrap(), 0.0)],
        1,
        0.0,
        0.0,
    )
    .unwrap();
    assert_eq!(negative, positive);
}

#[test]
fn semantic_mutations_move_edge_and_graph_identity() {
    let a = model("a");
    let b = model("b");
    let base = exact_context("temperature.mean", "reynolds", 1.0, 2.0);
    let qoi_mutant = exact_context("temperature.max", "reynolds", 1.0, 2.0);
    let domain_mutant = exact_context("temperature.mean", "reynolds", 1.0, 3.0);
    let edge = native_edge(a, b, "base", base);
    assert_ne!(edge.id(), native_edge(a, b, "base", qoi_mutant).id());
    assert_ne!(edge.id(), native_edge(a, b, "base", domain_mutant).id());
    assert_ne!(
        edge.id(),
        FidelityEdge::new(
            a,
            b,
            cost("cost-mutant"),
            discrepancy("base"),
            transfer("base"),
            ValidityDomain::universal(),
            Informativeness::legacy_total_order(),
        )
        .unwrap()
        .id()
    );
    assert_ne!(
        edge.id(),
        FidelityEdge::new(
            a,
            b,
            cost("base"),
            discrepancy("discrepancy-mutant"),
            transfer("base"),
            ValidityDomain::universal(),
            Informativeness::legacy_total_order(),
        )
        .unwrap()
        .id()
    );
}

#[test]
fn canonical_transport_round_trips_and_refuses_mutants() {
    let synthetic = synthetic_graph();
    let bytes = synthetic.graph.canonical_bytes();
    let decoded = FidelityGraph::from_canonical_bytes(&bytes).unwrap();
    assert_eq!(decoded, synthetic.graph);
    assert_eq!(decoded.identity(), synthetic.graph.identity());

    let mut version = bytes.clone();
    version[0] = 99;
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&version),
        Err(GraphError::UnsupportedSchema(99))
    ));

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&trailing),
        Err(GraphError::Decode("trailing canonical bytes"))
    ));

    // The first edge identity begins after the graph header and four
    // fixed-size node rows. Find it by decoding a one-edge graph where the
    // final 32-byte root can be mutated without relying on this graph's
    // variable-length offsets.
    let mut tiny = FidelityGraph::new("tiny").unwrap();
    let a = add_node(&mut tiny, "a");
    let b = add_node(&mut tiny, "b");
    tiny.add_edge(native_edge(a, b, "ab", ContextPredicateSet::universal()))
        .unwrap();
    let tiny_bytes = tiny.canonical_bytes();
    let edge_id = tiny.edges().values().next().unwrap().id();
    let offset = tiny_bytes
        .windows(32)
        .position(|window| window == edge_id.as_bytes())
        .unwrap();
    let mut identity_mutant = tiny_bytes;
    identity_mutant[offset] ^= 1;
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&identity_mutant),
        Err(GraphError::IdentityMismatch {
            what: "edge identity"
        })
    ));
}

#[test]
fn qoi_specific_queries_choose_different_models_and_explain_cost_not_authority() {
    let synthetic = synthetic_graph();
    let mean = synthetic
        .graph
        .best_model_for(
            synthetic.correlation,
            &context("temperature.mean", 10_000.0, 50.0),
            &synthetic.resolver,
        )
        .unwrap();
    assert_eq!(mean.model, synthetic.calibrated);
    assert_eq!(mean.predicted_cost_s.to_bits(), 2.0f64.to_bits());
    assert_eq!(mean.path.len(), 2);
    assert_eq!(mean.explanation.basis, SelectionBasis::UniqueGraphMaximum);
    assert!(mean.explanation.unresolved.is_empty());
    assert_eq!(mean.explanation.start, synthetic.correlation);
    assert_eq!(
        mean.explanation.max_relative_discrepancy.to_bits(),
        0.05f64.to_bits()
    );
    assert_eq!(
        mean.explanation.regime[&RegimeAxis::new("reynolds").unwrap()].to_bits(),
        10_000.0f64.to_bits()
    );

    let maximum = synthetic
        .graph
        .best_model_for(
            synthetic.correlation,
            &context("temperature.max", 10_000.0, 500.0),
            &synthetic.resolver,
        )
        .unwrap();
    assert_eq!(maximum.model, synthetic.les);
    assert_eq!(maximum.predicted_cost_s.to_bits(), 400.0f64.to_bits());
    assert_eq!(maximum.path.len(), 2);
}

#[test]
fn cheapest_adequate_uses_pairwise_discrepancy_and_not_rung_index() {
    let synthetic = synthetic_graph();
    let mean = synthetic
        .graph
        .cheapest_adequate(
            synthetic.correlation,
            &context("temperature.mean", 10_000.0, 50.0),
            &synthetic.resolver,
        )
        .unwrap();
    assert_eq!(mean.model, synthetic.correlation);
    assert_eq!(mean.predicted_cost_s.to_bits(), 1.0f64.to_bits());
    assert_eq!(mean.explanation.basis, SelectionBasis::CheapestAdequate);

    let maximum = synthetic
        .graph
        .cheapest_adequate(
            synthetic.correlation,
            &context("temperature.max", 10_000.0, 500.0),
            &synthetic.resolver,
        )
        .unwrap();
    assert_eq!(maximum.model, synthetic.rans);
    assert_eq!(maximum.predicted_cost_s.to_bits(), 40.0f64.to_bits());
}

#[test]
fn contradictory_or_unresolved_pairwise_evidence_blocks_adequacy() {
    let mut graph = FidelityGraph::new("adequacy-conjunction").unwrap();
    let a = add_node(&mut graph, "a");
    let b = add_node(&mut graph, "b");
    let c = add_node(&mut graph, "c");
    let ab = native_edge(a, b, "ab", ContextPredicateSet::universal());
    let ac = native_edge(a, c, "ac", ContextPredicateSet::universal());
    graph.add_edge(ab.clone()).unwrap();
    graph.add_edge(ac.clone()).unwrap();
    let mut contradictory = Resolver::default();
    contradictory.insert(&ab, 1.0, 2.0, Adequacy::Adequate, "ab");
    contradictory.insert(&ac, 1.0, 3.0, Adequacy::Inadequate, "ac");
    assert!(matches!(
        graph.cheapest_adequate(a, &context("anything", 1.0, 10.0), &contradictory),
        Err(QueryRefusal::NoAdequateModel { .. })
    ));

    let mut unresolved = Resolver::default();
    unresolved.insert(&ab, 1.0, 2.0, Adequacy::Adequate, "ab-only");
    assert!(matches!(
        graph.cheapest_adequate(a, &context("anything", 1.0, 10.0), &unresolved),
        Err(QueryRefusal::NoAdequateModel { .. })
    ));
}

#[test]
fn empty_informativeness_and_out_of_domain_are_honest_unknowns() {
    let mut graph = FidelityGraph::new("unknown").unwrap();
    let a = add_node(&mut graph, "a");
    let b = add_node(&mut graph, "b");
    let edge = FidelityEdge::new(
        a,
        b,
        cost("ab"),
        discrepancy("ab"),
        transfer("ab"),
        ValidityDomain::universal(),
        Informativeness::unknown(),
    )
    .unwrap();
    let mut resolver = Resolver::default();
    resolver.insert(&edge, 1.0, 2.0, Adequacy::Adequate, "ab");
    graph.add_edge(edge).unwrap();
    assert!(matches!(
        graph.best_model_for(a, &context("temperature.mean", 10_000.0, 5.0), &resolver),
        Err(QueryRefusal::NoApplicableEvidence { .. })
    ));

    let synthetic = synthetic_graph();
    assert!(matches!(
        synthetic.graph.best_model_for(
            synthetic.correlation,
            &context("temperature.mean", 1_000_000.0, 500.0),
            &synthetic.resolver
        ),
        Err(QueryRefusal::NoApplicableEvidence { .. })
    ));
}

#[test]
fn mismatched_resolution_refs_fail_closed_and_are_named() {
    let mut graph = FidelityGraph::new("mismatch").unwrap();
    let a = add_node(&mut graph, "a");
    let b = add_node(&mut graph, "b");
    let edge = native_edge(a, b, "ab", ContextPredicateSet::universal());
    graph.add_edge(edge.clone()).unwrap();
    let mut resolver = Resolver::default();
    resolver.rows.insert(
        edge.id(),
        ResolvedEdgeEvidence::new(
            Some(cost("wrong-cost")),
            Some(discrepancy("ab")),
            1.0,
            2.0,
            Some(0.01),
            query_receipt("mismatch"),
        )
        .unwrap(),
    );
    let refusal = graph
        .best_model_for(a, &context("anything", 1.0, 10.0), &resolver)
        .unwrap_err();
    assert!(matches!(refusal, QueryRefusal::NoApplicableEvidence { .. }));
}

#[test]
fn an_unaffordable_source_cannot_teleport_to_a_cheap_target() {
    let mut graph = FidelityGraph::new("budget-path").unwrap();
    let a = add_node(&mut graph, "a");
    let b = add_node(&mut graph, "b");
    let edge = native_edge(a, b, "ab", ContextPredicateSet::universal());
    graph.add_edge(edge.clone()).unwrap();
    let mut resolver = Resolver::default();
    resolver.insert(&edge, 100.0, 1.0, Adequacy::Unknown, "ab");
    assert!(matches!(
        graph.best_model_for(a, &context("anything", 1.0, 10.0), &resolver),
        Err(QueryRefusal::NoApplicableEvidence { .. })
    ));
}

#[test]
fn incomparable_maxima_use_replay_tie_break_without_epistemic_claim() {
    let mut graph = FidelityGraph::new("fork").unwrap();
    let root = add_node(&mut graph, "root");
    let left = add_node(&mut graph, "left");
    let right = add_node(&mut graph, "right");
    let predicates = ContextPredicateSet::universal();
    let root_left = native_edge(root, left, "root-left", predicates.clone());
    let root_right = native_edge(root, right, "root-right", predicates);
    graph.add_edge(root_right.clone()).unwrap();
    graph.add_edge(root_left.clone()).unwrap();
    let mut resolver = Resolver::default();
    resolver.insert(&root_left, 1.0, 10.0, Adequacy::Unknown, "left");
    resolver.insert(&root_right, 1.0, 10.0, Adequacy::Unknown, "right");
    let result = graph
        .best_model_for(root, &context("anything", 1.0, 20.0), &resolver)
        .unwrap();
    assert_eq!(
        result.explanation.basis,
        SelectionBasis::IncomparableMaximaCostThenIdentity
    );
    assert_eq!(result.model, left.min(right));
    assert_eq!(result.explanation.incomparable_maxima.len(), 1);
}

#[test]
fn adding_specific_evidence_never_makes_explanation_less_specific() {
    fn one_edge(predicates: ContextPredicateSet) -> (FidelityGraph, Resolver, ModelId) {
        let mut graph = FidelityGraph::new("specificity").unwrap();
        let a = add_node(&mut graph, "a");
        let b = add_node(&mut graph, "b");
        let edge = native_edge(a, b, "ab", predicates);
        let mut resolver = Resolver::default();
        resolver.insert(&edge, 1.0, 2.0, Adequacy::Unknown, "ab");
        graph.add_edge(edge).unwrap();
        (graph, resolver, a)
    }
    let (broad_graph, broad_resolver, broad_start) = one_edge(ContextPredicateSet::universal());
    let (specific_graph, specific_resolver, specific_start) =
        one_edge(exact_context("temperature.mean", "reynolds", 1.0, 2.0));
    let query = context("temperature.mean", 1.5, 10.0);
    let broad = broad_graph
        .best_model_for(broad_start, &query, &broad_resolver)
        .unwrap();
    let specific = specific_graph
        .best_model_for(specific_start, &query, &specific_resolver)
        .unwrap();
    assert_eq!(broad.model, specific.model);
    assert!(specific.explanation.matched_specificity >= broad.explanation.matched_specificity);
}

#[test]
fn removing_an_edge_cannot_increase_evidenced_path_depth() {
    let synthetic = synthetic_graph();
    let full = synthetic
        .graph
        .best_model_for(
            synthetic.correlation,
            &context("temperature.max", 10_000.0, 500.0),
            &synthetic.resolver,
        )
        .unwrap();

    let mut reduced = FidelityGraph::new("reduced").unwrap();
    let correlation = add_node(&mut reduced, "correlation");
    let rans = add_node(&mut reduced, "rans");
    let edge = native_edge(
        correlation,
        rans,
        "corr-rans-max",
        exact_context("temperature.max", "reynolds", 1_000.0, 100_000.0),
    );
    reduced.add_edge(edge.clone()).unwrap();
    let mut resolver = Resolver::default();
    resolver.insert(&edge, 1.0, 40.0, Adequacy::Inadequate, "reduced");
    let after_removal = reduced
        .best_model_for(
            correlation,
            &context("temperature.max", 10_000.0, 500.0),
            &resolver,
        )
        .unwrap();
    assert!(after_removal.path.len() <= full.path.len());
}

#[test]
fn legacy_ladder_embedding_is_lossless_and_moves_transfers_to_edges() {
    let ladder = Ladder::new("demo", "coarse", 1.0, "coarse note").then(
        Box::new(Refine1d),
        "fine",
        5.0,
        "fine note",
    );
    let descriptor = ladder.descriptor();
    let embedded = ladder.into_fidelity_graph().unwrap();
    assert_eq!(
        embedded.graph().embedded_ladder_descriptor().unwrap(),
        descriptor
    );
    assert_eq!(embedded.graph().nodes().len(), 2);
    assert_eq!(embedded.graph().edges().len(), 1);
    let edge = embedded.graph().edges().values().next().unwrap().id();
    let coarse = [0.0, 2.0, 4.0];
    let fine = embedded
        .transfers()
        .prolongate(embedded.graph(), edge, &coarse)
        .unwrap();
    let back = embedded
        .transfers()
        .restrict(embedded.graph(), edge, &fine)
        .unwrap();
    assert_eq!(
        back.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
        coarse
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    let bytes = embedded.graph().canonical_bytes();
    let decoded = FidelityGraph::from_canonical_bytes(&bytes).unwrap();
    assert_eq!(decoded.embedded_ladder_descriptor().unwrap(), descriptor);

    let legacy_kernel_offset = 2 + 4 + "demo".len() + 1 + 4;
    let mut mismatched_kernel = bytes.clone();
    mismatched_kernel[legacy_kernel_offset] = b'x';
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&mismatched_kernel),
        Err(GraphError::IdentityMismatch {
            what: "legacy kernel and graph name"
        })
    ));

    let note_offset = bytes
        .windows(b"coarse note".len())
        .position(|window| window == b"coarse note")
        .unwrap();
    let mut control_note = bytes.clone();
    control_note[note_offset] = b'\n';
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&control_note),
        Err(GraphError::InvalidName {
            field: "legacy rung note",
            ..
        })
    ));

    let mut empty_legacy = Vec::new();
    empty_legacy.extend_from_slice(&FIDELITY_GRAPH_SCHEMA_VERSION.to_le_bytes());
    empty_legacy.extend_from_slice(&(5_u32).to_le_bytes());
    empty_legacy.extend_from_slice(b"empty");
    empty_legacy.push(1);
    empty_legacy.extend_from_slice(&(5_u32).to_le_bytes());
    empty_legacy.extend_from_slice(b"empty");
    empty_legacy.extend_from_slice(&0_u32.to_le_bytes());
    empty_legacy.extend_from_slice(&0_u32.to_le_bytes());
    assert!(matches!(
        FidelityGraph::from_canonical_bytes(&empty_legacy),
        Err(GraphError::IdentityMismatch {
            what: "legacy ladder must contain a rung"
        })
    ));

    assert!(matches!(
        decoded
            .clone()
            .add_node(FidelityNode::new(model("extra"), card("extra"), "extra").unwrap()),
        Err(GraphError::LegacyEmbeddingImmutable)
    ));
}

#[test]
fn cht_instance_migrates_without_changing_the_legacy_registry() {
    let legacy = LadderRegistry::cht();
    let graph = LadderRegistry::cht_graph().unwrap();
    let descriptor = graph.graph().embedded_ladder_descriptor().unwrap();
    assert_eq!(descriptor.kernel(), legacy.ladder("cht").unwrap().kernel());
    assert_eq!(descriptor.rungs().len(), 3);
    assert_eq!(descriptor.rungs()[0].name, "correlation-Nu");
    assert_eq!(descriptor.rungs()[1].name, "RANS");
    assert_eq!(descriptor.rungs()[2].name, "LES");
}
