//! G0/G2/G3 conformance for the Phase 0B-C content-addressed evidence graph and
//! anytime consequence-times-doubt allocator.
//!
//! These tests intentionally exercise only descriptive graph snapshots and
//! allocation receipts. They do not treat a graph edge, score, or allocation
//! decision as authenticated evidence or statement truth.

use fs_blake3::ContentHash;
use fs_govern::{
    LaneCharter,
    evidence_contract::{
        AUTHORITY_ALGEBRA_VERSION, AssumptionSet, AttackEdge, AuthorityBudget, AuthorityState,
        CapabilityBinding, ClaimInstance, ClaimLaneBinding, ClaimStatement,
        CounterexampleCandidate, DomainVariable, EvidenceKind, EvidenceRef, ExactInstanceAdmission,
        FiveExplicits, InferenceRule, InvalidationState, KernelState, NoClaimBoundary,
        NonvacuityState, QuantifiedDomain, Quantifier, QuantifierBlock, ReproductionState,
        SatisfiabilityState, ScaleState, SupportEdge, TruthState, UnitSystem, VersionBinding,
    },
    evidence_graph::{
        ALLOCATION_CANDIDATE_IDENTITY_DOMAIN, ALLOCATION_DECISION_IDENTITY_DOMAIN,
        ALLOCATION_POLICY_IDENTITY_DOMAIN, ANYTIME_ACCOUNTING_IDENTITY_DOMAIN, AllocationCandidate,
        AllocationFloors, AllocationPolicy, AnytimeAccountingCandidate, DoubtProfile,
        EVIDENCE_GRAPH_AUTHORITY_ALGEBRA_TAG, EVIDENCE_GRAPH_VERSION, FIXED_RATE_SCALE, FixedRate,
        GRAPH_NODE_IDENTITY_DOMAIN, GRAPH_SNAPSHOT_IDENTITY_DOMAIN, GraphError, GraphNode,
        GraphNodeKind, GraphSnapshot, MAX_ALLOCATION_CANDIDATES, MAX_ALLOCATION_SELECTIONS,
        MAX_GRAPH_NODES, MAX_GRAPH_TEXT_BYTES, PlanningPass, WorkKind, plan_allocations,
        plan_allocations_with_cancel,
    },
};
use std::collections::BTreeSet;

/// Fixture-input namespace only; production graph identities use the public
/// v2 domains asserted independently below.
const TEST_EVIDENCE_GRAPH_FIXTURE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.test-evidence-graph.v2";

fn hash(label: &str) -> ContentHash {
    fs_blake3::hash_domain(
        TEST_EVIDENCE_GRAPH_FIXTURE_IDENTITY_DOMAIN,
        label.as_bytes(),
    )
}

fn push_graph_identity_field(out: &mut Vec<u8>, tag: u8, bytes: &[u8]) {
    out.push(tag);
    out.extend_from_slice(
        &u64::try_from(bytes.len())
            .expect("bounded graph identity fixture fits u64")
            .to_le_bytes(),
    );
    out.extend_from_slice(bytes);
}

fn lane(statement: &str, independence_class: &str) -> LaneCharter {
    LaneCharter::new(
        statement,
        "bounded design envelope",
        &["deterministic seed", "registered falsifier"],
        "candidate decision authority",
        "boring reference mechanism",
        "blind holdout and adversarial counterexample family",
        independence_class,
    )
    .expect("valid proof lane")
}

fn claim(statement_text: &str, independence_class: &str) -> ClaimInstance {
    let statement = ClaimStatement::new(&[statement_text]).expect("claim statement");
    let domain = QuantifiedDomain::new(
        vec![
            QuantifierBlock::commutative(
                Quantifier::ForAll,
                vec![
                    DomainVariable::new("design", "registered envelope").expect("domain variable"),
                ],
            )
            .expect("quantifier block"),
        ],
        &["finite candidate set"],
    )
    .expect("quantified domain");
    let assumptions = AssumptionSet::new(&["deterministic seed"]).expect("assumptions");
    let charter = lane(statement_text, independence_class);
    let statement_key = statement.identity().to_string();
    let binding = ClaimLaneBinding::new(
        &statement,
        &domain,
        &assumptions,
        &charter,
        hash(&format!("{statement_key}/lane-binding-artifact")),
        hash(&format!("{statement_key}/lane-binding-reviewer")),
    )
    .expect("claim/lane binding");
    let explicits = FiveExplicits::new(
        UnitSystem::dimensionless(),
        7,
        AuthorityBudget {
            work_units: 100,
            memory_bytes: 4096,
            wall_time_millis: 1_000,
            reviewer_slots: 1,
        },
        vec![VersionBinding::new("graph-fixture", "2").expect("version")],
        vec![CapabilityBinding::new("descriptive-planning", 1).expect("capability")],
    )
    .expect("Five Explicits");
    ClaimInstance::new(
        statement,
        domain,
        assumptions,
        binding,
        explicits,
        NoClaimBoundary::new(&[
            "not authenticated authority",
            "not calibrated probability",
            "not durable budget reservation",
        ])
        .expect("no-claim boundary"),
    )
    .expect("claim instance")
}

fn unknown_state(claim: &ClaimInstance) -> AuthorityState {
    AuthorityState::unknown(claim.clone()).expect("unknown authority state")
}

fn evidence(claim: &ClaimInstance, kind: EvidenceKind, label: &str) -> EvidenceRef {
    EvidenceRef::new(
        kind,
        claim.identity(),
        hash(&format!("{label}/artifact")),
        hash(&format!("{label}/checker")),
        AUTHORITY_ALGEBRA_VERSION,
    )
    .expect("evidence reference")
}

fn support(source: &AuthorityState, target: &ClaimInstance, label: &str) -> SupportEdge {
    let rule = InferenceRule::new(
        &format!("{label} support rule"),
        1,
        hash(&format!("{label}/rule-definition")),
    )
    .expect("inference rule");
    SupportEdge::new(source, target, &rule, support_evidence(target, label)).expect("support edge")
}

fn support_evidence(target: &ClaimInstance, label: &str) -> EvidenceRef {
    evidence(target, EvidenceKind::Support, label)
}

fn support_evidence_node(target: &ClaimInstance, label: &str) -> GraphNode {
    GraphNode::evidence(support_evidence(target, label))
}

fn attack(target: &ClaimInstance, label: &str) -> (CounterexampleCandidate, AttackEdge) {
    let candidate = CounterexampleCandidate::new(
        target,
        evidence(target, EvidenceKind::Counterexample, label),
    )
    .expect("counterexample candidate");
    let edge = AttackEdge::new(
        &candidate,
        target,
        evidence(target, EvidenceKind::Attack, &format!("{label}/attack")),
    )
    .expect("attack edge");
    (candidate, edge)
}

fn attack_evidence_nodes(target: &ClaimInstance, label: &str) -> [GraphNode; 2] {
    [
        GraphNode::evidence(evidence(target, EvidenceKind::Counterexample, label)),
        GraphNode::evidence(evidence(
            target,
            EvidenceKind::Attack,
            &format!("{label}/attack"),
        )),
    ]
}

fn doubt(
    calibrated_uncertainty: u32,
    attack_coverage: u32,
    independent_support: u32,
    assumption_resolution: u32,
) -> DoubtProfile {
    DoubtProfile::new(
        FixedRate::new(calibrated_uncertainty).expect("calibrated uncertainty"),
        FixedRate::new(attack_coverage).expect("attack coverage"),
        FixedRate::new(independent_support).expect("independent support"),
        FixedRate::new(assumption_resolution).expect("assumption resolution"),
    )
}

fn anytime(label: &str) -> AnytimeAccountingCandidate {
    AnytimeAccountingCandidate::new(
        "descriptive wsr e-process candidate",
        1,
        17,
        hash(&format!("{label}/anytime-state")),
        hash(&format!("{label}/anytime-evidence")),
    )
    .expect("anytime-accounting candidate")
}

#[derive(Clone, Copy)]
struct CandidateSpec<'a> {
    statement: &'a str,
    independence_class: &'a str,
    label: &'a str,
    work_artifact: &'a str,
    kind: WorkKind,
    cost: u64,
    utility_weight: u32,
    doubt: DoubtProfile,
    prior: &'a str,
    correlation_class: &'a str,
    include_anytime: bool,
}

impl<'a> CandidateSpec<'a> {
    fn new(
        statement: &'a str,
        independence_class: &'a str,
        label: &'a str,
        kind: WorkKind,
    ) -> Self {
        Self {
            statement,
            independence_class,
            label,
            work_artifact: label,
            kind,
            cost: 1,
            utility_weight: 1,
            doubt: doubt(
                500_000,
                FIXED_RATE_SCALE,
                FIXED_RATE_SCALE,
                FIXED_RATE_SCALE,
            ),
            prior: label,
            correlation_class: label,
            include_anytime: true,
        }
    }

    fn cost(mut self, cost: u64) -> Self {
        self.cost = cost;
        self
    }

    fn work(mut self, work_artifact: &'a str) -> Self {
        self.work_artifact = work_artifact;
        self
    }

    fn utility(mut self, utility_weight: u32) -> Self {
        self.utility_weight = utility_weight;
        self
    }

    fn doubt(mut self, doubt: DoubtProfile) -> Self {
        self.doubt = doubt;
        self
    }

    fn prior(mut self, prior: &'a str) -> Self {
        self.prior = prior;
        self
    }

    fn correlation(mut self, correlation_class: &'a str) -> Self {
        self.correlation_class = correlation_class;
        self
    }

    fn without_anytime(mut self) -> Self {
        self.include_anytime = false;
        self
    }
}

fn allocation_candidate(
    graph: &GraphSnapshot,
    claim: &ClaimInstance,
    spec: CandidateSpec<'_>,
) -> AllocationCandidate {
    AllocationCandidate::new(
        graph,
        claim,
        &lane(spec.statement, spec.independence_class),
        spec.kind,
        hash(&format!("{}/work", spec.work_artifact)),
        spec.cost,
        spec.utility_weight,
        spec.doubt,
        hash(&format!("{}/prior", spec.prior)),
        hash(&format!("{}/correlation", spec.correlation_class)),
        spec.include_anytime.then(|| anytime(spec.label)),
    )
    .expect("allocation candidate")
}

fn allocation_graph(specs: &[(&ClaimInstance, &str, u64)]) -> GraphSnapshot {
    let mut nodes = Vec::with_capacity(specs.len() * 2);
    for (claim, consumer_label, consequence) in specs {
        nodes.push(GraphNode::claim(claim));
        nodes.push(
            GraphNode::consumer(
                claim,
                hash(&format!("{consumer_label}/consumer")),
                *consequence,
            )
            .expect("allocation consumer"),
        );
    }
    GraphSnapshot::new(nodes, vec![], vec![]).expect("allocation graph")
}

#[allow(clippy::too_many_arguments)] // The policy contract deliberately exposes each guard.
fn allocation_policy(
    total_budget: u64,
    no_action_reserve: u64,
    max_selections: u32,
    min_correlation_classes: u32,
    require_anytime: bool,
    floors: AllocationFloors,
    label: &str,
) -> AllocationPolicy {
    AllocationPolicy::new(
        total_budget,
        no_action_reserve,
        max_selections,
        min_correlation_classes,
        require_anytime,
        floors,
        hash(&format!("{label}/utility-model")),
        hash(&format!("{label}/sensitivity")),
    )
    .expect("allocation policy")
}

fn permutations<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
    fn visit<T: Clone>(prefix: &mut Vec<T>, rest: &mut Vec<T>, out: &mut Vec<Vec<T>>) {
        if rest.is_empty() {
            out.push(prefix.clone());
            return;
        }
        for index in 0..rest.len() {
            let item = rest.remove(index);
            prefix.push(item.clone());
            visit(prefix, rest, out);
            prefix.pop();
            rest.insert(index, item);
        }
    }

    let mut out = Vec::new();
    visit(&mut Vec::new(), &mut items.to_vec(), &mut out);
    out
}

fn representative_orders<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
    let mut orders = vec![items.to_vec()];
    let mut reversed = items.to_vec();
    reversed.reverse();
    orders.push(reversed);
    for offset in 1..items.len() {
        let mut rotated = items.to_vec();
        rotated.rotate_left(offset);
        orders.push(rotated);
    }
    orders
}

#[test]
#[allow(clippy::too_many_lines)] // Independent v2 node and snapshot preimages are intentionally explicit.
fn g0_graph_v2_domains_and_node_snapshot_preimages_are_schema_locked() {
    assert_eq!(EVIDENCE_GRAPH_VERSION, 2);
    assert_eq!(EVIDENCE_GRAPH_AUTHORITY_ALGEBRA_TAG, 250);
    assert_eq!(
        GRAPH_NODE_IDENTITY_DOMAIN,
        "frankensim.fs-govern.evidence-graph-node.v2"
    );
    assert_eq!(
        GRAPH_SNAPSHOT_IDENTITY_DOMAIN,
        "frankensim.fs-govern.evidence-graph-snapshot.v2"
    );
    assert_eq!(
        ANYTIME_ACCOUNTING_IDENTITY_DOMAIN,
        "frankensim.fs-govern.anytime-accounting-candidate.v2"
    );
    assert_eq!(
        ALLOCATION_CANDIDATE_IDENTITY_DOMAIN,
        "frankensim.fs-govern.evidence-allocation-candidate.v2"
    );
    assert_eq!(
        ALLOCATION_POLICY_IDENTITY_DOMAIN,
        "frankensim.fs-govern.evidence-allocation-policy.v2"
    );
    assert_eq!(
        ALLOCATION_DECISION_IDENTITY_DOMAIN,
        "frankensim.fs-govern.evidence-allocation-decision.v2"
    );
    assert_eq!(
        TEST_EVIDENCE_GRAPH_FIXTURE_IDENTITY_DOMAIN,
        "frankensim.fs-govern.test-evidence-graph.v2"
    );

    let claim = claim("graph v2 preimage lock", "graph-v2-preimage/class");
    let state = unknown_state(&claim);
    let counterexample_evidence = evidence(
        &claim,
        EvidenceKind::Counterexample,
        "graph-v2-preimage/counterexample",
    );
    let counterexample = CounterexampleCandidate::new(&claim, counterexample_evidence)
        .expect("counterexample candidate");
    let checker = hash("graph-v2-preimage/checker");
    let falsifier = hash("graph-v2-preimage/falsifier");
    let consumer = hash("graph-v2-preimage/consumer");
    let nodes = [
        GraphNode::claim(&claim),
        GraphNode::authority_with_consequence(&state, 7).expect("nonzero authority consequence"),
        GraphNode::assumptions(&claim),
        GraphNode::evidence(counterexample_evidence),
        GraphNode::checker(&claim, checker).expect("checker node"),
        GraphNode::falsifier(&claim, falsifier).expect("falsifier node"),
        GraphNode::consumer(&claim, consumer, 11).expect("consumer node"),
        GraphNode::counterexample(&counterexample),
    ];

    // Claim: tag 1, exact claim and lane payload.
    let mut preimage = Vec::new();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[1]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 3, claim.proof_lane().as_hash().as_bytes());
    assert_eq!(
        nodes[0].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // AuthorityState: tag 2, state and consequence payload (not a duplicate Claim fixture).
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[2]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 3, state.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 4, &7_u64.to_le_bytes());
    assert_eq!(
        nodes[1].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // AssumptionSet: tag 3 and exact assumption-set root.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[3]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(
        &mut preimage,
        3,
        claim.assumptions().identity().as_hash().as_bytes(),
    );
    assert_eq!(
        nodes[2].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // Evidence: tag 4 and every retained EvidenceRef field.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[4]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(
        &mut preimage,
        3,
        counterexample_evidence.identity().as_hash().as_bytes(),
    );
    push_graph_identity_field(&mut preimage, 4, b"counterexample");
    push_graph_identity_field(
        &mut preimage,
        5,
        counterexample_evidence.artifact().as_bytes(),
    );
    push_graph_identity_field(
        &mut preimage,
        6,
        counterexample_evidence.checker().as_bytes(),
    );
    push_graph_identity_field(&mut preimage, 7, &2_u32.to_le_bytes());
    assert_eq!(
        nodes[3].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // Checker: tag 5.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[5]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 3, checker.as_bytes());
    assert_eq!(
        nodes[4].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // Falsifier: tag 6.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[6]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 3, falsifier.as_bytes());
    assert_eq!(
        nodes[5].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // Consumer: tag 7, semantic consumer identity and consequence.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[7]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(&mut preimage, 3, consumer.as_bytes());
    push_graph_identity_field(&mut preimage, 4, &11_u64.to_le_bytes());
    assert_eq!(
        nodes[6].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    // Counterexample: tag 8, candidate and exact counterexample-evidence roots.
    preimage.clear();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 1, &[8]);
    push_graph_identity_field(&mut preimage, 2, claim.identity().as_hash().as_bytes());
    push_graph_identity_field(
        &mut preimage,
        3,
        counterexample.identity().as_hash().as_bytes(),
    );
    push_graph_identity_field(
        &mut preimage,
        4,
        counterexample_evidence.identity().as_hash().as_bytes(),
    );
    assert_eq!(
        nodes[7].identity().as_hash(),
        &fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-node.v2", &preimage)
    );

    let snapshot = GraphSnapshot::new(nodes.into_iter().rev().collect(), vec![], vec![])
        .expect("all-node-kinds schema-lock snapshot");
    let mut canonical_node_identities = nodes.map(|node| node.identity());
    canonical_node_identities.sort_unstable();
    assert_eq!(
        snapshot
            .nodes()
            .iter()
            .map(GraphNode::identity)
            .collect::<Vec<_>>(),
        canonical_node_identities.to_vec()
    );
    let mut snapshot_preimage = Vec::new();
    push_graph_identity_field(&mut snapshot_preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut snapshot_preimage, 250, &2_u32.to_le_bytes());
    for node_identity in canonical_node_identities {
        push_graph_identity_field(
            &mut snapshot_preimage,
            1,
            node_identity.as_hash().as_bytes(),
        );
    }
    assert_eq!(
        snapshot.identity(),
        fs_blake3::hash_domain(
            "frankensim.fs-govern.evidence-graph-snapshot.v2",
            &snapshot_preimage,
        )
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Edge-bearing snapshot preimage pins three sorted collection tags.
fn g0_edge_bearing_snapshot_literal_preimage_pins_node_support_and_attack_order() {
    let upstream_a = claim("snapshot upstream A", "snapshot-edge/a");
    let upstream_b = claim("snapshot upstream B", "snapshot-edge/b");
    let downstream = claim("snapshot downstream", "snapshot-edge/downstream");
    let state_a = unknown_state(&upstream_a);
    let state_b = unknown_state(&upstream_b);
    let support_a = support(&state_a, &downstream, "snapshot-edge/support-a");
    let support_b = support(&state_b, &downstream, "snapshot-edge/support-b");
    let (counterexample_a, attack_a) = attack(&upstream_a, "snapshot-edge/attack-a");
    let (counterexample_b, attack_b) = attack(&upstream_b, "snapshot-edge/attack-b");

    let nodes = vec![
        GraphNode::claim(&upstream_a),
        GraphNode::authority(&state_a),
        GraphNode::claim(&upstream_b),
        GraphNode::authority(&state_b),
        GraphNode::claim(&downstream),
        GraphNode::counterexample(&counterexample_a),
        GraphNode::counterexample(&counterexample_b),
        support_evidence_node(&downstream, "snapshot-edge/support-a"),
        support_evidence_node(&downstream, "snapshot-edge/support-b"),
        attack_evidence_nodes(&upstream_a, "snapshot-edge/attack-a")[0],
        attack_evidence_nodes(&upstream_a, "snapshot-edge/attack-a")[1],
        attack_evidence_nodes(&upstream_b, "snapshot-edge/attack-b")[0],
        attack_evidence_nodes(&upstream_b, "snapshot-edge/attack-b")[1],
    ];
    let snapshot = GraphSnapshot::new(
        nodes.iter().copied().rev().collect(),
        vec![support_b, support_a],
        vec![attack_b, attack_a],
    )
    .expect("edge-bearing schema-lock snapshot");

    let mut expected_nodes = nodes.iter().map(GraphNode::identity).collect::<Vec<_>>();
    expected_nodes.sort_unstable();
    let mut expected_support = [support_a.identity(), support_b.identity()];
    expected_support.sort_unstable();
    let mut expected_attacks = [attack_a.identity(), attack_b.identity()];
    expected_attacks.sort_unstable();
    assert_eq!(
        snapshot
            .nodes()
            .iter()
            .map(GraphNode::identity)
            .collect::<Vec<_>>(),
        expected_nodes
    );
    assert_eq!(
        snapshot
            .support_edges()
            .iter()
            .map(SupportEdge::identity)
            .collect::<Vec<_>>(),
        expected_support.to_vec()
    );
    assert_eq!(
        snapshot
            .attack_edges()
            .iter()
            .map(AttackEdge::identity)
            .collect::<Vec<_>>(),
        expected_attacks.to_vec()
    );

    let mut preimage = Vec::new();
    push_graph_identity_field(&mut preimage, 0, &2_u32.to_le_bytes());
    push_graph_identity_field(&mut preimage, 250, &2_u32.to_le_bytes());
    for node in &expected_nodes {
        push_graph_identity_field(&mut preimage, 1, node.as_hash().as_bytes());
    }
    for edge in expected_support {
        push_graph_identity_field(&mut preimage, 2, edge.as_hash().as_bytes());
    }
    for edge in expected_attacks {
        push_graph_identity_field(&mut preimage, 3, edge.as_hash().as_bytes());
    }
    assert_eq!(
        snapshot.identity(),
        fs_blake3::hash_domain("frankensim.fs-govern.evidence-graph-snapshot.v2", &preimage,)
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Node, support-edge, attack-edge, and payload mutations move one root.
fn g0_graph_identity_is_permutation_invariant_and_semantic_mutations_move_the_root() {
    let upstream = claim("upstream residual is bounded", "identity/upstream");
    let downstream = claim("downstream objective is certified", "identity/downstream");
    let upstream_state = unknown_state(&upstream);
    let support = support(&upstream_state, &downstream, "identity");
    let alternate_support = support(&upstream_state, &downstream, "identity/alternate");
    let (counterexample, attack) = attack(&upstream, "identity");
    let (alternate_counterexample, alternate_attack) = attack(&upstream, "identity/alternate");
    let nodes = vec![
        GraphNode::claim(&upstream),
        GraphNode::authority(&upstream_state),
        GraphNode::claim(&downstream),
        GraphNode::consumer(&downstream, hash("identity/consumer"), 17).expect("consumer"),
        GraphNode::counterexample(&counterexample),
        GraphNode::counterexample(&alternate_counterexample),
        support_evidence_node(&downstream, "identity"),
        support_evidence_node(&downstream, "identity/alternate"),
        attack_evidence_nodes(&upstream, "identity")[0],
        attack_evidence_nodes(&upstream, "identity")[1],
        attack_evidence_nodes(&upstream, "identity/alternate")[0],
        attack_evidence_nodes(&upstream, "identity/alternate")[1],
    ];
    let expected = GraphSnapshot::new(
        nodes.clone(),
        vec![support, alternate_support],
        vec![attack, alternate_attack],
    )
    .expect("reference snapshot");

    for permutation in representative_orders(&nodes) {
        let actual = GraphSnapshot::new(
            permutation,
            vec![alternate_support, support],
            vec![alternate_attack, attack],
        )
        .expect("permuted snapshot");
        assert_eq!(actual.identity(), expected.identity());
        assert_eq!(actual.nodes(), expected.nodes());
    }
    let reversed_edges = GraphSnapshot::new(
        nodes.clone(),
        vec![alternate_support, support],
        vec![alternate_attack, attack],
    )
    .expect("reversed edge sets");
    assert_eq!(reversed_edges.identity(), expected.identity());

    let changed_weight = GraphSnapshot::new(
        vec![
            GraphNode::claim(&upstream),
            GraphNode::authority(&upstream_state),
            GraphNode::claim(&downstream),
            GraphNode::consumer(&downstream, hash("identity/consumer"), 18)
                .expect("changed consumer"),
            GraphNode::counterexample(&counterexample),
            GraphNode::counterexample(&alternate_counterexample),
            support_evidence_node(&downstream, "identity"),
            support_evidence_node(&downstream, "identity/alternate"),
            attack_evidence_nodes(&upstream, "identity")[0],
            attack_evidence_nodes(&upstream, "identity")[1],
            attack_evidence_nodes(&upstream, "identity/alternate")[0],
            attack_evidence_nodes(&upstream, "identity/alternate")[1],
        ],
        vec![support, alternate_support],
        vec![attack, alternate_attack],
    )
    .expect("weight-mutated snapshot");
    assert_ne!(changed_weight.identity(), expected.identity());

    let changed_consumer = GraphSnapshot::new(
        vec![
            GraphNode::claim(&upstream),
            GraphNode::authority(&upstream_state),
            GraphNode::claim(&downstream),
            GraphNode::consumer(&downstream, hash("identity/other-consumer"), 17)
                .expect("changed consumer identity"),
            GraphNode::counterexample(&counterexample),
            GraphNode::counterexample(&alternate_counterexample),
            support_evidence_node(&downstream, "identity"),
            support_evidence_node(&downstream, "identity/alternate"),
            attack_evidence_nodes(&upstream, "identity")[0],
            attack_evidence_nodes(&upstream, "identity")[1],
            attack_evidence_nodes(&upstream, "identity/alternate")[0],
            attack_evidence_nodes(&upstream, "identity/alternate")[1],
        ],
        vec![support, alternate_support],
        vec![attack, alternate_attack],
    )
    .expect("consumer-mutated snapshot");
    assert_ne!(changed_consumer.identity(), expected.identity());

    let removed_support =
        GraphSnapshot::new(nodes.clone(), vec![support], vec![attack, alternate_attack])
            .expect("support-edge removal");
    assert_ne!(removed_support.identity(), expected.identity());
    let removed_attack = GraphSnapshot::new(nodes, vec![support, alternate_support], vec![attack])
        .expect("attack-edge removal");
    assert_ne!(removed_attack.identity(), expected.identity());
}

#[test]
fn g3_cosmetic_claim_relabeling_preserves_graph_identity() {
    let canonical = claim("cosmetic graph relabel claim", "relabel/class");
    let relabeled = claim("  cosmetic   graph relabel   claim  ", "relabel/class");
    assert_eq!(canonical.identity(), relabeled.identity());
    assert_eq!(GraphNode::claim(&canonical), GraphNode::claim(&relabeled));

    let first = GraphSnapshot::new(
        vec![
            GraphNode::claim(&canonical),
            GraphNode::consumer(&canonical, hash("relabel/consumer"), 7)
                .expect("canonical consumer"),
        ],
        vec![],
        vec![],
    )
    .expect("canonical graph");
    let second = GraphSnapshot::new(
        vec![
            GraphNode::claim(&relabeled),
            GraphNode::consumer(&relabeled, hash("relabel/consumer"), 7)
                .expect("relabeled consumer"),
        ],
        vec![],
        vec![],
    )
    .expect("cosmetically relabeled graph");
    assert_eq!(first, second);
    assert_eq!(first.identity(), second.identity());
}

#[test]
fn g3_each_graph_node_kind_and_authority_weight_moves_snapshot_identity() {
    let claim = claim("graph node identity claim", "node-identity/class");
    let state = unknown_state(&claim);
    let base = GraphNode::claim(&claim);
    let evidence_ref = evidence(&claim, EvidenceKind::KernelProof, "node-identity/evidence");
    let evidence_node = GraphNode::evidence(evidence_ref);
    assert_eq!(
        evidence_node.kind(),
        GraphNodeKind::Evidence {
            evidence: evidence_ref.identity(),
            kind: EvidenceKind::KernelProof,
            claim: claim.identity(),
            artifact: evidence_ref.artifact(),
            checker: evidence_ref.checker(),
            schema_version: AUTHORITY_ALGEBRA_VERSION,
        }
    );
    let node_variants = vec![
        vec![base],
        vec![base, GraphNode::assumptions(&claim)],
        vec![base, evidence_node],
        vec![
            base,
            GraphNode::checker(&claim, hash("node-identity/checker")).expect("checker node"),
        ],
        vec![
            base,
            GraphNode::falsifier(&claim, hash("node-identity/falsifier")).expect("falsifier node"),
        ],
        vec![base, GraphNode::authority(&state)],
        vec![
            base,
            GraphNode::authority_with_consequence(&state, 2).expect("weighted authority node"),
        ],
        vec![
            base,
            GraphNode::consumer(&claim, hash("node-identity/consumer"), 1).expect("consumer node"),
        ],
    ];
    let roots = node_variants
        .into_iter()
        .map(|nodes| {
            GraphSnapshot::new(nodes, vec![], vec![])
                .expect("node-kind snapshot")
                .identity()
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(roots.len(), 8);
}

#[test]
fn g0_duplicate_semantics_and_unknown_claim_endpoints_are_refused() {
    let first = claim("first duplicate guard claim", "duplicates/first");
    let second = claim("second duplicate guard claim", "duplicates/second");
    let first_node = GraphNode::claim(&first);
    assert_eq!(
        GraphSnapshot::new(vec![first_node, first_node], vec![], vec![]),
        Err(GraphError::DuplicateNode {
            node: first_node.identity(),
        })
    );

    let shared_consumer = hash("duplicates/shared-consumer");
    let duplicate_semantic = GraphSnapshot::new(
        vec![
            GraphNode::claim(&first),
            GraphNode::claim(&second),
            GraphNode::consumer(&first, shared_consumer, 3).expect("first consumer"),
            GraphNode::consumer(&second, shared_consumer, 5).expect("second consumer"),
        ],
        vec![],
        vec![],
    );
    assert_eq!(
        duplicate_semantic,
        Err(GraphError::DuplicateSemanticIdentity {
            kind: "consumer",
            identity: shared_consumer,
        })
    );

    let orphan_consumer =
        GraphNode::consumer(&first, hash("duplicates/orphan"), 1).expect("orphan consumer node");
    assert_eq!(
        GraphSnapshot::new(
            vec![GraphNode::claim(&second), orphan_consumer],
            vec![],
            vec![]
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "claim",
            identity: *first.identity().as_hash(),
        })
    );

    assert!(matches!(
        GraphNode::consumer(&first, hash("duplicates/zero-consequence"), 0),
        Err(GraphError::InvalidValue {
            what: "consumer consequence must be nonzero"
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // One refusal table pins every support-DAG guard.
fn g0_support_refuses_duplicate_unknown_self_and_circular_edges() {
    let first = claim("support source one", "support/first");
    let second = claim("support source two", "support/second");
    let third = claim("support source three", "support/third");
    let first_state = unknown_state(&first);
    let second_state = unknown_state(&second);
    let third_state = unknown_state(&third);
    let first_to_second = support(&first_state, &second, "support/first-second");

    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&first),
                GraphNode::authority(&first_state),
                GraphNode::claim(&second),
                support_evidence_node(&second, "support/first-second"),
            ],
            vec![first_to_second, first_to_second],
            vec![],
        ),
        Err(GraphError::DuplicateEdge {
            kind: "support",
            identity: *first_to_second.identity().as_hash(),
        })
    );

    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&first),
                GraphNode::claim(&second),
                support_evidence_node(&second, "support/first-second"),
            ],
            vec![first_to_second],
            vec![],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "support source authority",
            identity: *first_state.identity().as_hash(),
        })
    );

    assert_eq!(
        GraphSnapshot::new(
            vec![GraphNode::claim(&first), GraphNode::authority(&first_state)],
            vec![first_to_second],
            vec![],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "support target claim",
            identity: *second.identity().as_hash(),
        })
    );

    let self_edge = support(&first_state, &first, "support/self");
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&first),
                GraphNode::authority(&first_state),
                support_evidence_node(&first, "support/self"),
            ],
            vec![self_edge],
            vec![],
        ),
        Err(GraphError::SelfSupport {
            claim: first.identity(),
        })
    );

    let second_to_third = support(&second_state, &third, "support/second-third");
    let third_to_first = support(&third_state, &first, "support/third-first");
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&first),
                GraphNode::authority(&first_state),
                GraphNode::claim(&second),
                GraphNode::authority(&second_state),
                GraphNode::claim(&third),
                GraphNode::authority(&third_state),
                support_evidence_node(&second, "support/first-second"),
                support_evidence_node(&third, "support/second-third"),
                support_evidence_node(&first, "support/third-first"),
            ],
            vec![first_to_second, second_to_third, third_to_first],
            vec![],
        ),
        Err(GraphError::SupportCycle)
    );
}

#[test]
fn g0_support_and_attack_edges_require_exact_evidence_nodes() {
    let source = claim("edge evidence source", "edge-evidence/source");
    let target = claim("edge evidence target", "edge-evidence/target");
    let source_state = unknown_state(&source);
    let support = support(&source_state, &target, "edge-evidence/support");
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&source),
                GraphNode::authority(&source_state),
                GraphNode::claim(&target),
            ],
            vec![support],
            vec![],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "support evidence",
            identity: *support.evidence().as_hash(),
        })
    );

    let (counterexample, attack) = attack(&target, "edge-evidence/attack");
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&target),
                GraphNode::counterexample(&counterexample),
            ],
            vec![],
            vec![attack],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "counterexample evidence",
            identity: *counterexample.evidence().as_hash(),
        })
    );
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&target),
                GraphNode::counterexample(&counterexample),
                attack_evidence_nodes(&target, "edge-evidence/attack")[0],
            ],
            vec![],
            vec![attack],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "attack evidence",
            identity: *attack.evidence().as_hash(),
        })
    );

    let complete = GraphSnapshot::new(
        vec![
            GraphNode::claim(&source),
            GraphNode::authority(&source_state),
            GraphNode::claim(&target),
            support_evidence_node(&target, "edge-evidence/support"),
            GraphNode::counterexample(&counterexample),
            attack_evidence_nodes(&target, "edge-evidence/attack")[0],
            attack_evidence_nodes(&target, "edge-evidence/attack")[1],
        ],
        vec![support],
        vec![attack],
    )
    .expect("both exact edge evidence nodes are present");
    assert_eq!(complete.support_edges(), &[support]);
    assert_eq!(complete.attack_edges(), &[attack]);
}

#[test]
fn g0_support_reachability_is_monotone_and_shared_consumers_are_counted_once() {
    let root = claim("root approximation is adequate", "reach/root");
    let left = claim("left derived bound holds", "reach/left");
    let right = claim("right derived bound holds", "reach/right");
    let leaf = claim("leaf decision is safe", "reach/leaf");
    let root_state = unknown_state(&root);
    let left_state = unknown_state(&left);
    let right_state = unknown_state(&right);
    let leaf_state = unknown_state(&leaf);

    let nodes = vec![
        GraphNode::claim(&root),
        GraphNode::authority(&root_state),
        GraphNode::consumer(&root, hash("reach/root-consumer"), 2).expect("root consumer"),
        GraphNode::claim(&left),
        GraphNode::authority(&left_state),
        GraphNode::consumer(&left, hash("reach/left-consumer"), 3).expect("left consumer"),
        GraphNode::claim(&right),
        GraphNode::authority(&right_state),
        GraphNode::consumer(&right, hash("reach/right-consumer"), 5).expect("right consumer"),
        GraphNode::claim(&leaf),
        GraphNode::authority_with_consequence(&leaf_state, 13).expect("leaf authority consequence"),
        GraphNode::consumer(&leaf, hash("reach/shared-leaf-consumer"), 11).expect("leaf consumer"),
        support_evidence_node(&left, "reach/root-left"),
        support_evidence_node(&right, "reach/root-right"),
        support_evidence_node(&leaf, "reach/left-leaf"),
        support_evidence_node(&leaf, "reach/right-leaf"),
    ];
    let decomposed = GraphSnapshot::new(
        nodes.clone(),
        vec![
            support(&root_state, &left, "reach/root-left"),
            support(&root_state, &right, "reach/root-right"),
            support(&left_state, &leaf, "reach/left-leaf"),
            support(&right_state, &leaf, "reach/right-leaf"),
        ],
        vec![],
    )
    .expect("diamond support graph");
    assert_eq!(decomposed.consequence(root.identity()), Ok(37));
    assert_eq!(decomposed.consequence(left.identity()), Ok(28));
    assert_eq!(decomposed.consequence(right.identity()), Ok(30));
    assert_eq!(decomposed.consequence(leaf.identity()), Ok(24));

    let without_right_branch = GraphSnapshot::new(
        nodes,
        vec![
            support(&root_state, &left, "reach/root-left"),
            support(&left_state, &leaf, "reach/left-leaf"),
        ],
        vec![],
    )
    .expect("strict subgraph");
    assert_eq!(without_right_branch.consequence(root.identity()), Ok(31));
    assert!(
        decomposed
            .consequence(root.identity())
            .expect("consequence")
            >= without_right_branch
                .consequence(root.identity())
                .expect("subgraph consequence")
    );

    let absent = claim("absent consequence query", "reach/absent");
    assert_eq!(
        decomposed.consequence(absent.identity()),
        Err(GraphError::UnknownEndpoint {
            kind: "consequence claim",
            identity: *absent.identity().as_hash(),
        })
    );
}

#[test]
fn g0_reachable_authority_consequence_is_explicit_deduplicated_and_nonzero() {
    let root = claim("authority consequence root", "authority-consequence/root");
    let downstream = claim(
        "authority consequence downstream",
        "authority-consequence/downstream",
    );
    let root_state = unknown_state(&root);
    let downstream_state = unknown_state(&downstream);
    let downstream_refuted = AuthorityState::new(
        downstream.clone(),
        TruthState::Refuted,
        SatisfiabilityState::Unknown,
        NonvacuityState::Unknown,
        ExactInstanceAdmission::NotEvaluated,
        KernelState::NotChecked,
        ScaleState::NotQualified,
        ReproductionState::NotAttempted,
        InvalidationState::Clear,
    )
    .expect("distinct downstream authority state");
    let edge = support(&root_state, &downstream, "authority-consequence/support");
    let snapshot = GraphSnapshot::new(
        vec![
            GraphNode::claim(&root),
            GraphNode::authority(&root_state),
            GraphNode::claim(&downstream),
            GraphNode::authority_with_consequence(&downstream_state, 50)
                .expect("weighted downstream authority"),
            GraphNode::authority_with_consequence(&downstream_refuted, 10)
                .expect("lower historical downstream authority"),
            support_evidence_node(&downstream, "authority-consequence/support"),
        ],
        vec![edge],
        vec![],
    )
    .expect("authority consequence graph");
    assert_eq!(snapshot.consequence(root.identity()), Ok(51));
    assert_eq!(snapshot.consequence(downstream.identity()), Ok(50));
    assert!(matches!(
        GraphNode::authority_with_consequence(&downstream_state, 0),
        Err(GraphError::InvalidValue {
            what: "authority-state consequence must be nonzero"
        })
    ));
}

#[test]
fn g0_attacks_are_exact_counted_and_do_not_create_support_cycles() {
    let first = claim("first mutually constrained claim", "attacks/first");
    let second = claim("second mutually constrained claim", "attacks/second");
    let (first_counterexample, first_attack) = attack(&first, "attacks/first");
    let (second_counterexample, second_attack) = attack(&second, "attacks/second");
    let nodes = vec![
        GraphNode::claim(&first),
        GraphNode::counterexample(&first_counterexample),
        GraphNode::claim(&second),
        GraphNode::counterexample(&second_counterexample),
        attack_evidence_nodes(&first, "attacks/first")[0],
        attack_evidence_nodes(&first, "attacks/first")[1],
        attack_evidence_nodes(&second, "attacks/second")[0],
        attack_evidence_nodes(&second, "attacks/second")[1],
    ];
    let snapshot = GraphSnapshot::new(nodes.clone(), vec![], vec![second_attack, first_attack])
        .expect("mutual threat descriptions are not circular support");
    assert_eq!(snapshot.attack_count(first.identity()), 1);
    assert_eq!(snapshot.attack_count(second.identity()), 1);
    assert_eq!(
        snapshot.attack_count(claim("absent attack", "attacks/absent").identity()),
        0
    );

    assert_eq!(
        GraphSnapshot::new(nodes.clone(), vec![], vec![first_attack, first_attack],),
        Err(GraphError::DuplicateEdge {
            kind: "attack",
            identity: *first_attack.identity().as_hash(),
        })
    );
    assert_eq!(
        GraphSnapshot::new(
            vec![
                GraphNode::claim(&first),
                GraphNode::claim(&second),
                attack_evidence_nodes(&first, "attacks/first")[1],
            ],
            vec![],
            vec![first_attack],
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "attack counterexample",
            identity: *first_counterexample.identity().as_hash(),
        })
    );
}

#[test]
fn g0_checked_consequence_arithmetic_refuses_overflow() {
    let claim = claim("overflow guard claim", "overflow");
    let snapshot = GraphSnapshot::new(
        vec![
            GraphNode::claim(&claim),
            GraphNode::consumer(&claim, hash("overflow/consumer-max"), u64::MAX)
                .expect("maximum consumer"),
            GraphNode::consumer(&claim, hash("overflow/consumer-one"), 1).expect("unit consumer"),
        ],
        vec![],
        vec![],
    )
    .expect("snapshot construction does not pre-sum consequences");
    assert_eq!(
        snapshot.consequence(claim.identity()),
        Err(GraphError::ArithmeticOverflow {
            what: "reachable consequence",
        })
    );
}

#[test]
fn g0_doubt_components_are_named_bounded_and_conservatively_combined() {
    let profile = DoubtProfile::new(
        FixedRate::new(200_000).expect("uncertainty"),
        FixedRate::new(750_000).expect("attack coverage"),
        FixedRate::new(500_000).expect("independent support"),
        FixedRate::new(800_000).expect("assumption resolution"),
    );
    assert_eq!(profile.combined().parts_per_million(), 760_000);
    assert_eq!(
        profile.calibrated_uncertainty().parts_per_million(),
        200_000
    );
    assert_eq!(profile.attack_coverage().parts_per_million(), 750_000);
    assert_eq!(profile.independent_support().parts_per_million(), 500_000);
    assert_eq!(profile.assumption_resolution().parts_per_million(), 800_000);

    let fully_resolved = DoubtProfile::new(
        FixedRate::ZERO,
        FixedRate::ONE,
        FixedRate::ONE,
        FixedRate::ONE,
    );
    assert_eq!(fully_resolved.combined(), FixedRate::ZERO);

    let uncovered_attacks = DoubtProfile::new(
        FixedRate::ZERO,
        FixedRate::ZERO,
        FixedRate::ONE,
        FixedRate::ONE,
    );
    assert_eq!(uncovered_attacks.combined(), FixedRate::ONE);
    assert!(uncovered_attacks.combined() > profile.combined());
    assert_eq!(
        FixedRate::new(FIXED_RATE_SCALE + 1),
        Err(GraphError::InvalidRate {
            observed: FIXED_RATE_SCALE + 1,
        })
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One table owns tie order and both budget caps.
fn g0_allocator_conserves_budget_and_breaks_ties_by_canonical_identity() {
    let statements = [
        "tie candidate alpha",
        "tie candidate beta",
        "tie candidate gamma",
    ];
    let classes = ["tie/class-alpha", "tie/class-beta", "tie/class-gamma"];
    let claims = statements
        .iter()
        .zip(classes)
        .map(|(statement, class)| claim(statement, class))
        .collect::<Vec<_>>();
    let graph = allocation_graph(&[
        (&claims[0], "tie/alpha", 100),
        (&claims[1], "tie/beta", 100),
        (&claims[2], "tie/gamma", 100),
    ]);
    let candidates = claims
        .iter()
        .zip(statements)
        .zip(classes)
        .enumerate()
        .map(|(index, ((claim, statement), class))| {
            allocation_candidate(
                &graph,
                claim,
                CandidateSpec::new(
                    statement,
                    class,
                    ["tie/alpha", "tie/beta", "tie/gamma"][index],
                    WorkKind::IndependentCheck,
                )
                .cost(3)
                .correlation(["tie/corr-alpha", "tie/corr-beta", "tie/corr-gamma"][index]),
            )
        })
        .collect::<Vec<_>>();
    let policy = allocation_policy(10, 1, 2, 2, true, AllocationFloors::zero(), "tie");
    let expected = plan_allocations(&graph, &policy, &candidates).expect("tie decision");
    assert_eq!(expected.used_budget(), 6);
    assert_eq!(expected.unallocated_budget(), 4);
    assert_eq!(
        expected.used_budget() + expected.unallocated_budget(),
        policy.total_budget()
    );
    assert!(expected.unallocated_budget() >= policy.no_action_reserve());
    assert_eq!(expected.selected().len(), 2);
    assert!(!expected.no_action_selected());
    assert_eq!(expected.graph(), graph.identity());
    assert_eq!(expected.policy(), policy.identity());
    assert_eq!(expected.utility_model(), policy.utility_model());
    assert_eq!(
        expected.sensitivity_artifact(),
        policy.sensitivity_artifact()
    );
    assert_eq!(expected.no_action_reserve(), 1);
    assert_eq!(expected.unused_allocatable_budget(), 3);
    assert_eq!(
        expected.no_action_reserve() + expected.unused_allocatable_budget(),
        expected.unallocated_budget()
    );
    assert_eq!(
        expected
            .selected()
            .iter()
            .map(|reservation| reservation.cost())
            .sum::<u64>(),
        expected.used_budget()
    );

    let ranked_ids = expected
        .ranked()
        .iter()
        .map(|row| row.candidate())
        .collect::<Vec<_>>();
    let mut sorted_ids = ranked_ids.clone();
    sorted_ids.sort_unstable();
    assert_eq!(ranked_ids, sorted_ids);
    for row in expected.ranked() {
        let candidate = candidates
            .iter()
            .find(|candidate| candidate.identity() == row.candidate())
            .expect("ranked candidate retained in input set");
        assert_eq!(row.claim(), candidate.claim());
        assert_eq!(row.work_artifact(), candidate.work_artifact());
        assert_eq!(row.prior(), candidate.prior());
        assert_eq!(row.doubt_profile(), candidate.doubt());
        assert_eq!(row.anytime_state(), candidate.anytime());
        assert_eq!(
            row.anytime_accounting(),
            candidate
                .anytime()
                .map(AnytimeAccountingCandidate::identity)
        );
    }
    for reservation in expected.selected() {
        let candidate = candidates
            .iter()
            .find(|candidate| candidate.identity() == reservation.candidate())
            .expect("selected candidate retained in input set");
        assert_eq!(
            reservation.correlation_class(),
            candidate.correlation_class()
        );
        assert_eq!(reservation.prior(), candidate.prior());
        assert_eq!(reservation.doubt_profile(), candidate.doubt());
        assert_eq!(reservation.anytime_state(), candidate.anytime());
        assert_eq!(
            reservation.anytime_accounting(),
            candidate
                .anytime()
                .map(AnytimeAccountingCandidate::identity)
        );
    }
    assert!(
        expected
            .ranked()
            .windows(2)
            .all(|pair| pair[0].score().weighted_value() == pair[1].score().weighted_value())
    );
    assert!(expected.ranked()[0].selected());
    assert!(expected.ranked()[1].selected());
    assert!(!expected.ranked()[2].selected());

    for permutation in permutations(&candidates) {
        let replay = plan_allocations(&graph, &policy, &permutation).expect("permuted replay");
        assert_eq!(replay, expected);
        assert_eq!(replay.identity(), expected.identity());
    }

    let budget_bound = allocation_policy(
        8,
        2,
        3,
        0,
        true,
        AllocationFloors::zero(),
        "tie/budget-bound",
    );
    let bounded = plan_allocations(&graph, &budget_bound, &candidates).expect("budget-bound plan");
    assert_eq!(bounded.selected().len(), 2);
    assert_eq!(bounded.used_budget(), 6);
    assert_eq!(bounded.unallocated_budget(), 2);
}

#[test]
#[allow(clippy::too_many_lines)] // Each orthogonal anti-gaming cap needs an isolated witness.
fn g0_lane_independence_correlation_scalar_budget_and_selection_caps_each_bind() {
    let first_statement = "cap guard first claim";
    let second_statement = "cap guard second claim";
    let first = claim(first_statement, "caps/shared-independent-class");
    let second = claim(second_statement, "caps/shared-independent-class");
    let graph = allocation_graph(&[(&first, "caps/first", 100), (&second, "caps/second", 90)]);

    let same_lane = [
        allocation_candidate(
            &graph,
            &first,
            CandidateSpec::new(
                first_statement,
                "caps/shared-independent-class",
                "caps/same-lane-a",
                WorkKind::Falsification,
            )
            .correlation("caps/same-lane-corr-a"),
        ),
        allocation_candidate(
            &graph,
            &first,
            CandidateSpec::new(
                first_statement,
                "caps/shared-independent-class",
                "caps/same-lane-b",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/same-lane-corr-b"),
        ),
    ];
    let open_policy = allocation_policy(10, 0, 10, 0, true, AllocationFloors::zero(), "caps/open");
    assert_eq!(
        plan_allocations(&graph, &open_policy, &same_lane)
            .expect("same-lane plan")
            .selected()
            .len(),
        1
    );

    let split_lane_same_class = [
        allocation_candidate(
            &graph,
            &first,
            CandidateSpec::new(
                first_statement,
                "caps/shared-independent-class",
                "caps/split-a",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/split-corr-a"),
        ),
        allocation_candidate(
            &graph,
            &second,
            CandidateSpec::new(
                second_statement,
                "caps/shared-independent-class",
                "caps/split-b",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/split-corr-b"),
        ),
    ];
    assert_ne!(
        split_lane_same_class[0].lane(),
        split_lane_same_class[1].lane()
    );
    assert_eq!(
        split_lane_same_class[0].independence_class(),
        split_lane_same_class[1].independence_class()
    );
    assert_eq!(
        plan_allocations(&graph, &open_policy, &split_lane_same_class)
            .expect("split-lane same-class plan")
            .selected()
            .len(),
        1
    );

    let second_distinct_statement = "cap guard distinct-class claim";
    let second_distinct = claim(second_distinct_statement, "caps/distinct-independent-class");
    let correlation_graph = allocation_graph(&[
        (&first, "caps/correlation-first", 100),
        (&second_distinct, "caps/correlation-second", 90),
    ]);
    let same_correlation = [
        allocation_candidate(
            &correlation_graph,
            &first,
            CandidateSpec::new(
                first_statement,
                "caps/shared-independent-class",
                "caps/correlation-a",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/shared-correlation"),
        ),
        allocation_candidate(
            &correlation_graph,
            &second_distinct,
            CandidateSpec::new(
                second_distinct_statement,
                "caps/distinct-independent-class",
                "caps/correlation-b",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/shared-correlation"),
        ),
    ];
    assert_eq!(
        plan_allocations(&correlation_graph, &open_policy, &same_correlation)
            .expect("correlated plan")
            .selected()
            .len(),
        1
    );

    let one_selection =
        allocation_policy(10, 0, 1, 0, true, AllocationFloors::zero(), "caps/global");
    let distinct = [
        same_correlation[0].clone(),
        allocation_candidate(
            &correlation_graph,
            &second_distinct,
            CandidateSpec::new(
                second_distinct_statement,
                "caps/distinct-independent-class",
                "caps/global-b",
                WorkKind::IndependentCheck,
            )
            .correlation("caps/global-corr-b"),
        ),
    ];
    assert_eq!(
        plan_allocations(&correlation_graph, &one_selection, &distinct)
            .expect("global selection cap")
            .selected()
            .len(),
        1
    );
}

#[test]
#[allow(clippy::too_many_lines)] // All four floors and no-action are one coupled policy law.
fn g0_floors_empty_selection_and_no_action_reserve_are_explicit_and_budget_bounded() {
    let inputs = [
        (
            "floor falsification claim",
            "floors/falsification-class",
            "floors/falsification",
            WorkKind::Falsification,
        ),
        (
            "floor independent-check claim",
            "floors/independent-class",
            "floors/independent",
            WorkKind::IndependentCheck,
        ),
        (
            "floor holdout claim",
            "floors/holdout-class",
            "floors/holdout",
            WorkKind::Holdout,
        ),
        (
            "floor exploration claim",
            "floors/exploration-class",
            "floors/exploration",
            WorkKind::Exploration,
        ),
    ];
    let claims = inputs
        .iter()
        .map(|(statement, class, _, _)| claim(statement, class))
        .collect::<Vec<_>>();
    let graph = allocation_graph(&[
        (&claims[0], "floors/falsification", 1),
        (&claims[1], "floors/independent", 1),
        (&claims[2], "floors/holdout", 1),
        (&claims[3], "floors/exploration", 1),
    ]);
    let resolved_doubt = doubt(0, FIXED_RATE_SCALE, FIXED_RATE_SCALE, FIXED_RATE_SCALE);
    let candidates = inputs
        .iter()
        .enumerate()
        .map(|(index, (statement, class, label, kind))| {
            allocation_candidate(
                &graph,
                &claims[index],
                CandidateSpec::new(statement, class, label, *kind)
                    .cost(2)
                    .doubt(resolved_doubt)
                    .correlation(label),
            )
        })
        .collect::<Vec<_>>();
    let floors = AllocationFloors {
        falsification: 2,
        independent_check: 2,
        holdout: 2,
        exploration: 2,
    };
    let policy = allocation_policy(10, 2, 4, 4, true, floors, "floors");
    let decision = plan_allocations(&graph, &policy, &candidates).expect("floor plan");
    assert_eq!(decision.used_budget(), 8);
    assert_eq!(decision.unallocated_budget(), 2);
    assert_eq!(decision.selected().len(), 4);
    assert!(!decision.no_action_selected());
    assert!(
        decision
            .selected()
            .iter()
            .all(|reservation| reservation.score().weighted_value() == 0)
    );
    let selected_kinds = decision
        .selected()
        .iter()
        .map(|reservation| reservation.kind())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        selected_kinds,
        BTreeSet::from([
            WorkKind::Falsification,
            WorkKind::IndependentCheck,
            WorkKind::Holdout,
            WorkKind::Exploration,
        ])
    );

    let no_floor = allocation_policy(
        10,
        2,
        4,
        0,
        true,
        AllocationFloors::zero(),
        "floors/no-action",
    );
    let no_action = plan_allocations(&graph, &no_floor, &candidates).expect("no-action plan");
    assert!(no_action.no_action_selected());
    assert!(no_action.selected().is_empty());
    assert_eq!(no_action.used_budget(), 0);
    assert_eq!(no_action.unallocated_budget(), 10);

    assert_eq!(
        plan_allocations(&graph, &policy, &candidates[..3]),
        Err(GraphError::FloorUnsatisfied {
            kind: WorkKind::Exploration,
            required: 2,
            selected: 0,
        })
    );

    let positive_candidates = candidates
        .iter()
        .enumerate()
        .map(|(index, _)| {
            let (statement, class, label, kind) = inputs[index];
            allocation_candidate(
                &graph,
                &claims[index],
                CandidateSpec::new(statement, class, label, kind)
                    .cost(2)
                    .correlation(label),
            )
        })
        .collect::<Vec<_>>();
    let reserve_bound = allocation_policy(
        10,
        4,
        4,
        0,
        true,
        AllocationFloors::zero(),
        "floors/reserve-bound",
    );
    let reserved = plan_allocations(&graph, &reserve_bound, &positive_candidates)
        .expect("no-action reserve bound");
    assert_eq!(reserved.selected().len(), 3);
    assert_eq!(reserved.used_budget(), 6);
    assert_eq!(reserved.unallocated_budget(), 4);
}

#[test]
fn g0_snapshot_duplicates_and_anytime_accounting_are_refused_before_selection() {
    let statement = "input validation candidate claim";
    let class = "input-validation/class";
    let claim = claim(statement, class);
    let first_graph = allocation_graph(&[(&claim, "input-validation/first", 10)]);
    let second_graph = allocation_graph(&[(&claim, "input-validation/second", 11)]);
    let complete = allocation_candidate(
        &first_graph,
        &claim,
        CandidateSpec::new(
            statement,
            class,
            "input-validation/complete",
            WorkKind::IndependentCheck,
        ),
    );
    let missing_anytime = allocation_candidate(
        &first_graph,
        &claim,
        CandidateSpec::new(
            statement,
            class,
            "input-validation/missing-anytime",
            WorkKind::IndependentCheck,
        )
        .without_anytime(),
    );
    let required = allocation_policy(
        10,
        0,
        2,
        0,
        true,
        AllocationFloors::zero(),
        "input-validation/required",
    );

    assert_eq!(
        plan_allocations(&second_graph, &required, &[complete.clone()]),
        Err(GraphError::SnapshotMismatch)
    );
    assert_eq!(
        plan_allocations(
            &first_graph,
            &required,
            &[complete.clone(), complete.clone()]
        ),
        Err(GraphError::DuplicateCandidate {
            identity: complete.identity(),
        })
    );
    assert_eq!(
        plan_allocations(&first_graph, &required, &[missing_anytime.clone()]),
        Err(GraphError::AnytimeAccountingRequired {
            candidate: missing_anytime.identity(),
        })
    );

    let optional = allocation_policy(
        10,
        0,
        2,
        0,
        false,
        AllocationFloors::zero(),
        "input-validation/optional",
    );
    let accepted = plan_allocations(&first_graph, &optional, &[missing_anytime])
        .expect("policy explicitly permits absent descriptive anytime accounting");
    assert_eq!(accepted.selected().len(), 1);

    assert_eq!(
        AnytimeAccountingCandidate::new(
            "descriptive wsr e-process candidate",
            1,
            0,
            hash("input-validation/zero-observation-state"),
            hash("input-validation/zero-observation-evidence"),
        ),
        Err(GraphError::InvalidValue {
            what: "anytime-accounting observation count must be nonzero",
        })
    );
}

#[test]
fn g0_allocation_candidates_require_registered_claim_exact_lane_cost_and_utility() {
    let statement = "allocation constructor claim";
    let class = "allocation-constructor/class";
    let registered = claim(statement, class);
    let foreign_statement = "allocation constructor foreign claim";
    let foreign_class = "allocation-constructor/foreign-class";
    let foreign = claim(foreign_statement, foreign_class);
    let graph = allocation_graph(&[(&registered, "allocation-constructor", 10)]);
    let profile = doubt(
        500_000,
        FIXED_RATE_SCALE,
        FIXED_RATE_SCALE,
        FIXED_RATE_SCALE,
    );
    let accounting = Some(anytime("allocation-constructor"));

    assert_eq!(
        AllocationCandidate::new(
            &graph,
            &foreign,
            &lane(foreign_statement, foreign_class),
            WorkKind::IndependentCheck,
            hash("allocation-constructor/foreign-work"),
            1,
            1,
            profile,
            hash("allocation-constructor/foreign-prior"),
            hash("allocation-constructor/foreign-correlation"),
            accounting.clone(),
        ),
        Err(GraphError::UnknownEndpoint {
            kind: "allocation claim",
            identity: *foreign.identity().as_hash(),
        })
    );
    assert_eq!(
        AllocationCandidate::new(
            &graph,
            &registered,
            &lane(foreign_statement, foreign_class),
            WorkKind::IndependentCheck,
            hash("allocation-constructor/wrong-lane-work"),
            1,
            1,
            profile,
            hash("allocation-constructor/wrong-lane-prior"),
            hash("allocation-constructor/wrong-lane-correlation"),
            accounting.clone(),
        ),
        Err(GraphError::ClaimLaneMismatch)
    );
    assert!(matches!(
        AllocationCandidate::new(
            &graph,
            &registered,
            &lane(statement, class),
            WorkKind::IndependentCheck,
            hash("allocation-constructor/zero-cost-work"),
            0,
            1,
            profile,
            hash("allocation-constructor/zero-cost-prior"),
            hash("allocation-constructor/zero-cost-correlation"),
            accounting.clone(),
        ),
        Err(GraphError::InvalidValue {
            what: "allocation candidate cost must be nonzero"
        })
    ));
    assert!(matches!(
        AllocationCandidate::new(
            &graph,
            &registered,
            &lane(statement, class),
            WorkKind::IndependentCheck,
            hash("allocation-constructor/zero-utility-work"),
            1,
            0,
            profile,
            hash("allocation-constructor/zero-utility-prior"),
            hash("allocation-constructor/zero-utility-correlation"),
            accounting,
        ),
        Err(GraphError::InvalidValue {
            what: "allocation utility weight must be nonzero"
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // Exact-limit and plus-one mutations share bounded fixtures.
fn g0_public_limits_and_checked_policy_arithmetic_refuse_before_partial_planning() {
    let statement = "bounded limit claim";
    let class = "limits/class";
    let claim = claim(statement, class);
    let node = GraphNode::claim(&claim);
    assert_eq!(
        GraphSnapshot::new(vec![node; MAX_GRAPH_NODES + 1], vec![], vec![]),
        Err(GraphError::TooLarge {
            what: "nodes",
            observed: MAX_GRAPH_NODES + 1,
            cap: MAX_GRAPH_NODES,
        })
    );

    let oversized_method = "x".repeat(MAX_GRAPH_TEXT_BYTES + 1);
    assert_eq!(
        AnytimeAccountingCandidate::new(
            &oversized_method,
            1,
            1,
            hash("limits/oversized-state"),
            hash("limits/oversized-evidence"),
        ),
        Err(GraphError::TooLarge {
            what: "anytime-accounting method",
            observed: MAX_GRAPH_TEXT_BYTES + 1,
            cap: MAX_GRAPH_TEXT_BYTES,
        })
    );

    assert!(matches!(
        AllocationPolicy::new(
            0,
            0,
            1,
            0,
            false,
            AllocationFloors::zero(),
            hash("limits/zero-budget-model"),
            hash("limits/zero-budget-sensitivity"),
        ),
        Err(GraphError::InvalidValue {
            what: "allocation budget must be nonzero"
        })
    ));
    assert!(matches!(
        AllocationPolicy::new(
            2,
            0,
            1,
            2,
            false,
            AllocationFloors::zero(),
            hash("limits/diversity-model"),
            hash("limits/diversity-sensitivity"),
        ),
        Err(GraphError::InvalidValue {
            what: "correlation diversity floor exceeds selection cap"
        })
    ));
    assert!(matches!(
        AllocationPolicy::new(
            2,
            1,
            2,
            0,
            false,
            AllocationFloors {
                falsification: 2,
                independent_check: 0,
                holdout: 0,
                exploration: 0,
            },
            hash("limits/floor-budget-model"),
            hash("limits/floor-budget-sensitivity"),
        ),
        Err(GraphError::InvalidValue {
            what: "allocation floors exceed budget after no-action reserve"
        })
    ));
    assert!(matches!(
        AllocationPolicy::new(
            1,
            2,
            1,
            0,
            false,
            AllocationFloors::zero(),
            hash("limits/reserve-model"),
            hash("limits/reserve-sensitivity"),
        ),
        Err(GraphError::InvalidValue {
            what: "no-action reserve exceeds total budget"
        })
    ));
    let maximum_selection_cap =
        u32::try_from(MAX_ALLOCATION_SELECTIONS).expect("selection cap fits u32");
    assert!(
        AllocationPolicy::new(
            1,
            0,
            maximum_selection_cap,
            0,
            false,
            AllocationFloors::zero(),
            hash("limits/max-selection-model"),
            hash("limits/max-selection-sensitivity"),
        )
        .is_ok()
    );
    assert!(matches!(
        AllocationPolicy::new(
            1,
            0,
            maximum_selection_cap + 1,
            0,
            false,
            AllocationFloors::zero(),
            hash("limits/over-selection-model"),
            hash("limits/over-selection-sensitivity"),
        ),
        Err(GraphError::InvalidValue {
            what: "allocation selection cap is outside the supported range"
        })
    ));
    assert_eq!(
        AllocationPolicy::new(
            u64::MAX,
            0,
            2,
            0,
            false,
            AllocationFloors {
                falsification: u64::MAX,
                independent_check: 1,
                holdout: 0,
                exploration: 0,
            },
            hash("limits/floor-overflow-model"),
            hash("limits/floor-overflow-sensitivity"),
        ),
        Err(GraphError::ArithmeticOverflow {
            what: "floor budget",
        })
    );

    let graph = allocation_graph(&[(&claim, "limits", 1)]);
    let candidate = allocation_candidate(
        &graph,
        &claim,
        CandidateSpec::new(
            statement,
            class,
            "limits/candidate",
            WorkKind::IndependentCheck,
        ),
    );
    let policy = allocation_policy(
        1,
        0,
        1,
        0,
        true,
        AllocationFloors::zero(),
        "limits/candidate-count",
    );
    assert_eq!(
        plan_allocations(
            &graph,
            &policy,
            &vec![candidate; MAX_ALLOCATION_CANDIDATES + 1],
        ),
        Err(GraphError::TooLarge {
            what: "allocation candidates",
            observed: MAX_ALLOCATION_CANDIDATES + 1,
            cap: MAX_ALLOCATION_CANDIDATES,
        })
    );
}

#[test]
fn g0_zero_floors_do_not_sum_irrelevant_huge_candidate_costs() {
    let first_statement = "huge no-action candidate one";
    let second_statement = "huge no-action candidate two";
    let first = claim(first_statement, "huge-cost/first-class");
    let second = claim(second_statement, "huge-cost/second-class");
    let graph = allocation_graph(&[
        (&first, "huge-cost/first", 1),
        (&second, "huge-cost/second", 1),
    ]);
    let resolved = doubt(0, FIXED_RATE_SCALE, FIXED_RATE_SCALE, FIXED_RATE_SCALE);
    let candidates = [
        allocation_candidate(
            &graph,
            &first,
            CandidateSpec::new(
                first_statement,
                "huge-cost/first-class",
                "huge-cost/first",
                WorkKind::IndependentCheck,
            )
            .cost(u64::MAX)
            .doubt(resolved),
        ),
        allocation_candidate(
            &graph,
            &second,
            CandidateSpec::new(
                second_statement,
                "huge-cost/second-class",
                "huge-cost/second",
                WorkKind::IndependentCheck,
            )
            .cost(u64::MAX)
            .doubt(resolved),
        ),
    ];
    let policy = allocation_policy(1, 0, 2, 0, true, AllocationFloors::zero(), "huge-cost");
    let decision = plan_allocations(&graph, &policy, &candidates)
        .expect("zero floors do not inspect irrelevant aggregate cost");
    assert!(decision.no_action_selected());
    assert_eq!(decision.used_budget(), 0);
    assert_eq!(decision.unallocated_budget(), 1);
}

#[test]
fn g0_cancellation_at_every_pass_returns_no_partial_decision_and_replay_is_stable() {
    let statement = "cancellable falsification claim";
    let class = "cancellation/class";
    let claim = claim(statement, class);
    let graph = allocation_graph(&[(&claim, "cancellation", 100)]);
    let candidate = allocation_candidate(
        &graph,
        &claim,
        CandidateSpec::new(
            statement,
            class,
            "cancellation/candidate",
            WorkKind::Falsification,
        )
        .cost(2),
    );
    let candidates = [candidate];
    let policy = allocation_policy(
        4,
        1,
        1,
        1,
        true,
        AllocationFloors {
            falsification: 2,
            independent_check: 0,
            holdout: 0,
            exploration: 0,
        },
        "cancellation",
    );
    let baseline = plan_allocations(&graph, &policy, &candidates).expect("baseline decision");

    for target in [
        PlanningPass::InputValidation,
        PlanningPass::Feasibility,
        PlanningPass::Floor(WorkKind::Falsification),
        PlanningPass::Ranking,
        PlanningPass::Finalization,
    ] {
        let mut observed_target = false;
        let cancelled = plan_allocations_with_cancel(&graph, &policy, &candidates, |pass| {
            if pass == target {
                observed_target = true;
                true
            } else {
                false
            }
        });
        assert!(observed_target, "planner never polled {target:?}");
        assert_eq!(cancelled, Err(GraphError::Cancelled { pass: target }));
        assert_eq!(
            plan_allocations(&graph, &policy, &candidates),
            Ok(baseline.clone()),
            "a cancelled pass must not mutate a later replay"
        );
    }

    let mut finalization_polls = 0u8;
    let after_identity_assembly =
        plan_allocations_with_cancel(&graph, &policy, &candidates, |pass| {
            if pass == PlanningPass::Finalization {
                finalization_polls += 1;
                finalization_polls == 2
            } else {
                false
            }
        });
    assert_eq!(finalization_polls, 2);
    assert_eq!(
        after_identity_assembly,
        Err(GraphError::Cancelled {
            pass: PlanningPass::Finalization,
        })
    );
    assert_eq!(
        plan_allocations(&graph, &policy, &candidates),
        Ok(baseline),
        "cancellation after receipt assembly must still publish no partial decision"
    );
}

#[test]
fn g2_high_consequence_weak_claim_wins_without_starving_blind_falsification() {
    let high_statement = "high consequence weak structural claim";
    let blind_statement = "blind falsification sentinel";
    let high = claim(high_statement, "priority/high-class");
    let blind = claim(blind_statement, "priority/blind-class");
    let graph = allocation_graph(&[
        (&high, "priority/high", 1_000),
        (&blind, "priority/blind", 1),
    ]);
    let high_candidate = allocation_candidate(
        &graph,
        &high,
        CandidateSpec::new(
            high_statement,
            "priority/high-class",
            "priority/high-candidate",
            WorkKind::IndependentCheck,
        )
        .cost(4)
        .utility(2)
        .doubt(doubt(
            900_000,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
        ))
        .correlation("priority/high-correlation"),
    );
    let blind_candidate = allocation_candidate(
        &graph,
        &blind,
        CandidateSpec::new(
            blind_statement,
            "priority/blind-class",
            "priority/blind-candidate",
            WorkKind::Falsification,
        )
        .cost(2)
        .doubt(doubt(
            100_000,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
        ))
        .correlation("priority/blind-correlation"),
    );
    let policy = allocation_policy(
        6,
        0,
        2,
        2,
        true,
        AllocationFloors {
            falsification: 2,
            independent_check: 0,
            holdout: 0,
            exploration: 0,
        },
        "priority",
    );
    let decision = plan_allocations(
        &graph,
        &policy,
        &[blind_candidate.clone(), high_candidate.clone()],
    )
    .expect("priority plan");
    assert_eq!(decision.ranked()[0].candidate(), high_candidate.identity());
    assert!(
        decision.ranked()[0].score().weighted_value()
            > decision.ranked()[1].score().weighted_value()
    );
    let selected = decision
        .selected()
        .iter()
        .map(|reservation| reservation.candidate())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        selected,
        BTreeSet::from([high_candidate.identity(), blind_candidate.identity()])
    );
    assert_eq!(decision.used_budget(), 6);
    assert_eq!(decision.unallocated_budget(), 0);
}

#[test]
#[allow(clippy::too_many_lines)] // Duplicate, diversity, and de-correlation mutations are coupled.
fn g3_correlated_duplicates_cannot_multiply_budget_or_fake_diversity() {
    let first_statement = "correlated duplicate first claim";
    let second_statement = "correlated duplicate second claim";
    let sentinel_statement = "independent sentinel claim";
    let first = claim(first_statement, "correlation/first-class");
    let second = claim(second_statement, "correlation/second-class");
    let sentinel = claim(sentinel_statement, "correlation/sentinel-class");
    let graph = allocation_graph(&[
        (&first, "correlation/first", 100),
        (&second, "correlation/second", 100),
        (&sentinel, "correlation/sentinel", 10),
    ]);
    let first_duplicate = allocation_candidate(
        &graph,
        &first,
        CandidateSpec::new(
            first_statement,
            "correlation/first-class",
            "correlation/duplicate-a",
            WorkKind::IndependentCheck,
        )
        .cost(2)
        .correlation("correlation/shared-evidence"),
    );
    let second_duplicate = allocation_candidate(
        &graph,
        &second,
        CandidateSpec::new(
            second_statement,
            "correlation/second-class",
            "correlation/duplicate-b",
            WorkKind::IndependentCheck,
        )
        .cost(2)
        .correlation("correlation/shared-evidence"),
    );
    let sentinel_candidate = allocation_candidate(
        &graph,
        &sentinel,
        CandidateSpec::new(
            sentinel_statement,
            "correlation/sentinel-class",
            "correlation/sentinel-candidate",
            WorkKind::Exploration,
        )
        .cost(2)
        .correlation("correlation/independent-sentinel"),
    );
    let candidates = [
        first_duplicate.clone(),
        second_duplicate.clone(),
        sentinel_candidate.clone(),
    ];
    let policy = allocation_policy(6, 0, 3, 2, true, AllocationFloors::zero(), "correlation");
    let decision = plan_allocations(&graph, &policy, &candidates).expect("correlation plan");
    let selected = decision
        .selected()
        .iter()
        .map(|reservation| reservation.candidate())
        .collect::<BTreeSet<_>>();
    assert_eq!(selected.len(), 2);
    assert!(selected.contains(&sentinel_candidate.identity()));
    assert_ne!(
        selected.contains(&first_duplicate.identity()),
        selected.contains(&second_duplicate.identity()),
        "exactly one representative of a correlation class may consume budget"
    );
    assert_eq!(decision.used_budget(), 4);
    assert_eq!(decision.unallocated_budget(), 2);

    let impossible_diversity = allocation_policy(
        6,
        0,
        3,
        3,
        true,
        AllocationFloors::zero(),
        "correlation/impossible-diversity",
    );
    assert_eq!(
        plan_allocations(&graph, &impossible_diversity, &candidates),
        Err(GraphError::DiversityFloorUnsatisfied {
            required: 3,
            selected: 2,
        })
    );

    let independent_second = allocation_candidate(
        &graph,
        &second,
        CandidateSpec::new(
            second_statement,
            "correlation/second-class",
            "correlation/independent-b",
            WorkKind::IndependentCheck,
        )
        .cost(2)
        .correlation("correlation/truly-independent-b"),
    );
    let independent_candidates = [first_duplicate, independent_second, sentinel_candidate];
    let independent = plan_allocations(&graph, &impossible_diversity, &independent_candidates)
        .expect("three genuinely distinct correlation classes");
    assert_eq!(independent.selected().len(), 3);
    assert_eq!(independent.used_budget(), 6);
}

#[test]
fn g3_support_decomposition_preserves_reachable_consequence_and_score() {
    let root_statement = "decomposition root claim";
    let middle_statement = "decomposition intermediate lemma";
    let leaf_statement = "decomposition consumer claim";
    let root = claim(root_statement, "decomposition/root-class");
    let middle = claim(middle_statement, "decomposition/middle-class");
    let leaf = claim(leaf_statement, "decomposition/leaf-class");
    let root_state = unknown_state(&root);
    let middle_state = unknown_state(&middle);
    let nodes = vec![
        GraphNode::claim(&root),
        GraphNode::authority(&root_state),
        GraphNode::claim(&middle),
        GraphNode::authority(&middle_state),
        GraphNode::claim(&leaf),
        GraphNode::consumer(&leaf, hash("decomposition/leaf-consumer"), 100)
            .expect("leaf consumer"),
        support_evidence_node(&leaf, "decomposition/direct"),
        support_evidence_node(&middle, "decomposition/root-middle"),
        support_evidence_node(&leaf, "decomposition/middle-leaf"),
    ];
    let direct = GraphSnapshot::new(
        nodes.clone(),
        vec![
            support(&root_state, &middle, "decomposition/root-middle"),
            support(&root_state, &leaf, "decomposition/direct"),
        ],
        vec![],
    )
    .expect("direct support graph");
    let decomposed = GraphSnapshot::new(
        nodes,
        vec![
            support(&root_state, &middle, "decomposition/root-middle"),
            support(&middle_state, &leaf, "decomposition/middle-leaf"),
        ],
        vec![],
    )
    .expect("decomposed support graph");
    assert_ne!(direct.identity(), decomposed.identity());
    assert_eq!(direct.consequence(root.identity()), Ok(102));
    assert_eq!(decomposed.consequence(root.identity()), Ok(102));

    let spec = CandidateSpec::new(
        root_statement,
        "decomposition/root-class",
        "decomposition/candidate",
        WorkKind::IndependentCheck,
    );
    let direct_candidate = allocation_candidate(&direct, &root, spec);
    let decomposed_candidate = allocation_candidate(&decomposed, &root, spec);
    assert_ne!(direct_candidate.identity(), decomposed_candidate.identity());
    let policy = allocation_policy(1, 0, 1, 1, true, AllocationFloors::zero(), "decomposition");
    let direct_decision =
        plan_allocations(&direct, &policy, &[direct_candidate]).expect("direct decision");
    let decomposed_decision = plan_allocations(&decomposed, &policy, &[decomposed_candidate])
        .expect("decomposed decision");
    assert_eq!(
        direct_decision.ranked()[0].score().weighted_value(),
        decomposed_decision.ranked()[0].score().weighted_value()
    );
    assert_ne!(direct_decision.identity(), decomposed_decision.identity());
}

#[test]
#[allow(clippy::too_many_lines)] // The diamond, duplicate attacks, and clean control are one law.
fn g3_attack_coverage_and_transitive_attack_count_change_descriptive_score() {
    let root_statement = "attacked upstream claim";
    let left_statement = "left attack propagation lemma";
    let right_statement = "right attack propagation lemma";
    let child_statement = "downstream allocation claim";
    let root = claim(root_statement, "attack-score/root-class");
    let left = claim(left_statement, "attack-score/left-class");
    let right = claim(right_statement, "attack-score/right-class");
    let child = claim(child_statement, "attack-score/child-class");
    let root_state = unknown_state(&root);
    let left_state = unknown_state(&left);
    let right_state = unknown_state(&right);
    let support = [
        support(&root_state, &left, "attack-score/root-left"),
        support(&root_state, &right, "attack-score/root-right"),
        support(&left_state, &child, "attack-score/left-child"),
        support(&right_state, &child, "attack-score/right-child"),
    ];
    let (counterexample_a, attack_a) = attack(&root, "attack-score/a");
    let (counterexample_b, attack_b) = attack(&root, "attack-score/b");
    let attacked_graph = GraphSnapshot::new(
        vec![
            GraphNode::claim(&root),
            GraphNode::authority(&root_state),
            GraphNode::claim(&left),
            GraphNode::authority(&left_state),
            GraphNode::claim(&right),
            GraphNode::authority(&right_state),
            GraphNode::claim(&child),
            GraphNode::consumer(&child, hash("attack-score/consumer"), 100).expect("consumer"),
            support_evidence_node(&left, "attack-score/root-left"),
            support_evidence_node(&right, "attack-score/root-right"),
            support_evidence_node(&child, "attack-score/left-child"),
            support_evidence_node(&child, "attack-score/right-child"),
            GraphNode::counterexample(&counterexample_a),
            GraphNode::counterexample(&counterexample_b),
            attack_evidence_nodes(&root, "attack-score/a")[0],
            attack_evidence_nodes(&root, "attack-score/a")[1],
            attack_evidence_nodes(&root, "attack-score/b")[0],
            attack_evidence_nodes(&root, "attack-score/b")[1],
        ],
        support.to_vec(),
        vec![attack_a, attack_b],
    )
    .expect("attacked graph");
    let clean_graph = GraphSnapshot::new(
        vec![
            GraphNode::claim(&root),
            GraphNode::authority(&root_state),
            GraphNode::claim(&left),
            GraphNode::authority(&left_state),
            GraphNode::claim(&right),
            GraphNode::authority(&right_state),
            GraphNode::claim(&child),
            GraphNode::consumer(&child, hash("attack-score/consumer"), 100).expect("consumer"),
            support_evidence_node(&left, "attack-score/root-left"),
            support_evidence_node(&right, "attack-score/root-right"),
            support_evidence_node(&child, "attack-score/left-child"),
            support_evidence_node(&child, "attack-score/right-child"),
        ],
        support.to_vec(),
        vec![],
    )
    .expect("clean comparison graph");
    let covered_profile = doubt(
        200_000,
        FIXED_RATE_SCALE,
        FIXED_RATE_SCALE,
        FIXED_RATE_SCALE,
    );
    let attacked = allocation_candidate(
        &attacked_graph,
        &child,
        CandidateSpec::new(
            child_statement,
            "attack-score/child-class",
            "attack-score/attacked",
            WorkKind::Falsification,
        )
        .doubt(covered_profile),
    );
    let clean = allocation_candidate(
        &clean_graph,
        &child,
        CandidateSpec::new(
            child_statement,
            "attack-score/child-class",
            "attack-score/clean",
            WorkKind::Falsification,
        )
        .doubt(covered_profile),
    );
    let policy = allocation_policy(1, 0, 1, 1, true, AllocationFloors::zero(), "attack-score");
    let attacked_decision = plan_allocations(&attacked_graph, &policy, &[attacked])
        .expect("attacked downstream decision");
    let clean_decision =
        plan_allocations(&clean_graph, &policy, &[clean]).expect("clean downstream decision");
    let attacked_score = attacked_decision.ranked()[0].score();
    let clean_score = clean_decision.ranked()[0].score();
    assert_eq!(attacked_graph.attack_count(root.identity()), 2);
    assert_eq!(attacked_graph.attack_count(left.identity()), 2);
    assert_eq!(attacked_graph.attack_count(right.identity()), 2);
    assert_eq!(attacked_graph.attack_count(child.identity()), 2);
    assert_eq!(attacked_graph.direct_attack_count(root.identity()), 2);
    assert_eq!(attacked_graph.direct_attack_count(child.identity()), 0);
    assert_eq!(clean_graph.attack_count(child.identity()), 0);
    assert_eq!(attacked_score.attack_count(), 2);
    assert_eq!(clean_score.attack_count(), 0);
    assert_eq!(attacked_score.consequence(), 100);
    assert_eq!(attacked_score.doubt(), FixedRate::ONE);
    assert_eq!(attacked_score.weighted_value(), 100);
    assert_eq!(clean_score.doubt().parts_per_million(), 200_000);
    assert_eq!(clean_score.weighted_value(), 20);
    assert!(attacked_score.weighted_value() > clean_score.weighted_value());
}

#[test]
fn g3_utility_perturbation_changes_rank_while_prior_perturbation_is_retained() {
    let target_statement = "utility perturbation target";
    let reference_statement = "utility perturbation reference";
    let target = claim(target_statement, "perturbation/target-class");
    let reference = claim(reference_statement, "perturbation/reference-class");
    let graph = allocation_graph(&[
        (&target, "perturbation/target", 100),
        (&reference, "perturbation/reference", 100),
    ]);
    let target_base = allocation_candidate(
        &graph,
        &target,
        CandidateSpec::new(
            target_statement,
            "perturbation/target-class",
            "perturbation/target",
            WorkKind::IndependentCheck,
        )
        .utility(1)
        .correlation("perturbation/target-correlation"),
    );
    let reference_candidate = allocation_candidate(
        &graph,
        &reference,
        CandidateSpec::new(
            reference_statement,
            "perturbation/reference-class",
            "perturbation/reference",
            WorkKind::IndependentCheck,
        )
        .utility(2)
        .correlation("perturbation/reference-correlation"),
    );
    let policy = allocation_policy(1, 0, 1, 1, true, AllocationFloors::zero(), "perturbation");
    let baseline = plan_allocations(
        &graph,
        &policy,
        &[target_base.clone(), reference_candidate.clone()],
    )
    .expect("baseline ranking");
    assert_eq!(
        baseline.ranked()[0].candidate(),
        reference_candidate.identity()
    );

    let target_boosted = allocation_candidate(
        &graph,
        &target,
        CandidateSpec::new(
            target_statement,
            "perturbation/target-class",
            "perturbation/target",
            WorkKind::IndependentCheck,
        )
        .utility(3)
        .correlation("perturbation/target-correlation"),
    );
    let boosted = plan_allocations(
        &graph,
        &policy,
        &[reference_candidate, target_boosted.clone()],
    )
    .expect("boosted ranking");
    assert_eq!(boosted.ranked()[0].candidate(), target_boosted.identity());
    assert!(
        boosted.ranked()[0].score().weighted_value()
            > baseline.ranked()[0].score().weighted_value()
    );

    let prior_mutated = allocation_candidate(
        &graph,
        &target,
        CandidateSpec::new(
            target_statement,
            "perturbation/target-class",
            "perturbation/target",
            WorkKind::IndependentCheck,
        )
        .utility(1)
        .prior("perturbation/alternative-prior")
        .correlation("perturbation/target-correlation"),
    );
    assert_ne!(prior_mutated.identity(), target_base.identity());
    let base_prior_decision =
        plan_allocations(&graph, &policy, &[target_base]).expect("base-prior decision");
    let changed_prior_decision =
        plan_allocations(&graph, &policy, &[prior_mutated]).expect("changed-prior decision");
    assert_eq!(
        base_prior_decision.ranked()[0].score(),
        changed_prior_decision.ranked()[0].score(),
        "the opaque prior artifact is retained but is not invented into a numeric prior"
    );
    assert_ne!(
        base_prior_decision.identity(),
        changed_prior_decision.identity(),
        "a different retained prior must move the decision receipt identity"
    );
}

#[test]
fn g3_candidate_identity_is_sensitive_to_every_planning_input() {
    let statement = "candidate identity mutation claim";
    let class = "candidate-identity/class";
    let claim = claim(statement, class);
    let graph = allocation_graph(&[(&claim, "candidate-identity", 100)]);
    let base = CandidateSpec::new(
        statement,
        class,
        "candidate-identity/base",
        WorkKind::IndependentCheck,
    )
    .work("candidate-identity/work")
    .prior("candidate-identity/prior")
    .correlation("candidate-identity/correlation");
    let variants = [
        base,
        base.cost(2),
        base.utility(2),
        base.doubt(doubt(
            600_000,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
            FIXED_RATE_SCALE,
        )),
        base.doubt(doubt(500_000, 900_000, FIXED_RATE_SCALE, FIXED_RATE_SCALE)),
        base.doubt(doubt(500_000, FIXED_RATE_SCALE, 900_000, FIXED_RATE_SCALE)),
        base.doubt(doubt(500_000, FIXED_RATE_SCALE, FIXED_RATE_SCALE, 900_000)),
        base.prior("candidate-identity/other-prior"),
        base.correlation("candidate-identity/other-correlation"),
        base.work("candidate-identity/other-work"),
        base.without_anytime(),
        CandidateSpec {
            kind: WorkKind::Falsification,
            ..base
        },
        CandidateSpec {
            label: "candidate-identity/other-anytime",
            ..base
        },
    ];
    let identities = variants
        .into_iter()
        .map(|spec| allocation_candidate(&graph, &claim, spec).identity())
        .collect::<BTreeSet<_>>();
    assert_eq!(identities.len(), variants.len());

    let changed_graph = allocation_graph(&[(&claim, "candidate-identity", 101)]);
    let original = allocation_candidate(&graph, &claim, base);
    let rebound = allocation_candidate(&changed_graph, &claim, base);
    assert_ne!(original.graph(), rebound.graph());
    assert_ne!(original.identity(), rebound.identity());
}

#[test]
#[allow(clippy::too_many_lines)] // Each canonical field has one change-only mutation.
fn g3_anytime_accounting_and_policy_identities_bind_every_declared_field() {
    let make_accounting = |method, schema, observations, state, evidence| {
        AnytimeAccountingCandidate::new(method, schema, observations, hash(state), hash(evidence))
            .expect("anytime identity fixture")
    };
    let accounting = [
        make_accounting("wsr e-process", 1, 7, "anytime/state", "anytime/evidence"),
        make_accounting(
            "safe wsr e-process",
            1,
            7,
            "anytime/state",
            "anytime/evidence",
        ),
        make_accounting("wsr e-process", 2, 7, "anytime/state", "anytime/evidence"),
        make_accounting("wsr e-process", 1, 8, "anytime/state", "anytime/evidence"),
        make_accounting(
            "wsr e-process",
            1,
            7,
            "anytime/other-state",
            "anytime/evidence",
        ),
        make_accounting(
            "wsr e-process",
            1,
            7,
            "anytime/state",
            "anytime/other-evidence",
        ),
    ];
    assert_eq!(
        accounting
            .iter()
            .map(AnytimeAccountingCandidate::identity)
            .collect::<BTreeSet<_>>()
            .len(),
        accounting.len()
    );
    assert_eq!(
        accounting[0].identity(),
        make_accounting(
            "  wsr   e-process ",
            1,
            7,
            "anytime/state",
            "anytime/evidence",
        )
        .identity()
    );

    let base_floors = AllocationFloors {
        falsification: 1,
        independent_check: 0,
        holdout: 0,
        exploration: 0,
    };
    let make_policy = |total,
                       reserve,
                       max_selections,
                       min_classes,
                       require_anytime,
                       floors,
                       utility,
                       sensitivity| {
        AllocationPolicy::new(
            total,
            reserve,
            max_selections,
            min_classes,
            require_anytime,
            floors,
            hash(utility),
            hash(sensitivity),
        )
        .expect("policy identity fixture")
    };
    let policies = [
        make_policy(
            10,
            1,
            4,
            1,
            true,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            11,
            1,
            4,
            1,
            true,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            2,
            4,
            1,
            true,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            5,
            1,
            true,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            4,
            2,
            true,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            4,
            1,
            false,
            base_floors,
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            4,
            1,
            true,
            AllocationFloors {
                falsification: 0,
                independent_check: 1,
                holdout: 0,
                exploration: 0,
            },
            "policy/utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            4,
            1,
            true,
            base_floors,
            "policy/other-utility",
            "policy/sensitivity",
        ),
        make_policy(
            10,
            1,
            4,
            1,
            true,
            base_floors,
            "policy/utility",
            "policy/other-sensitivity",
        ),
    ];
    assert_eq!(
        policies
            .iter()
            .map(AllocationPolicy::identity)
            .collect::<BTreeSet<_>>()
            .len(),
        policies.len()
    );
}

#[test]
fn g3_crossing_lane_and_correlation_constraints_choose_the_feasible_portfolio() {
    let first_statement = "crossing constraint first claim";
    let second_statement = "crossing constraint second claim";
    let first = claim(first_statement, "crossing/first-class");
    let second = claim(second_statement, "crossing/second-class");
    let graph = allocation_graph(&[
        (&first, "crossing/first", 100),
        (&second, "crossing/second", 100),
    ]);
    let tempting = allocation_candidate(
        &graph,
        &first,
        CandidateSpec::new(
            first_statement,
            "crossing/first-class",
            "crossing/tempting",
            WorkKind::Falsification,
        )
        .utility(3)
        .correlation("crossing/correlation-x"),
    );
    let lane_alternative = allocation_candidate(
        &graph,
        &first,
        CandidateSpec::new(
            first_statement,
            "crossing/first-class",
            "crossing/lane-alternative",
            WorkKind::Falsification,
        )
        .utility(2)
        .correlation("crossing/correlation-y"),
    );
    let correlation_alternative = allocation_candidate(
        &graph,
        &second,
        CandidateSpec::new(
            second_statement,
            "crossing/second-class",
            "crossing/correlation-alternative",
            WorkKind::Holdout,
        )
        .utility(1)
        .correlation("crossing/correlation-x"),
    );
    let policy = allocation_policy(
        2,
        0,
        2,
        2,
        true,
        AllocationFloors {
            falsification: 1,
            independent_check: 0,
            holdout: 1,
            exploration: 0,
        },
        "crossing",
    );
    let candidates = vec![
        tempting.clone(),
        correlation_alternative.clone(),
        lane_alternative.clone(),
    ];
    let decision =
        plan_allocations(&graph, &policy, &candidates).expect("a feasible floor portfolio exists");
    let selected = decision
        .selected()
        .iter()
        .map(|reservation| reservation.candidate())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        selected,
        BTreeSet::from([
            lane_alternative.identity(),
            correlation_alternative.identity(),
        ])
    );
    assert!(!selected.contains(&tempting.identity()));
    assert_eq!(decision.used_budget(), 2);
    for permutation in permutations(&candidates) {
        assert_eq!(
            plan_allocations(&graph, &policy, &permutation),
            Ok(decision.clone())
        );
    }
}
