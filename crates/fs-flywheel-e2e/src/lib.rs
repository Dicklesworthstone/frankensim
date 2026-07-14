//! FLYWHEEL CLOSES (bead lmp4.18; [F] — behind the `flywheel-e2e`
//! feature until its Gauntlet tier is green): the whole-loop harness
//! testing the addendum's CENTRAL CLAIM — that speculation (9),
//! incremental recompute (2), the sheaf-adjudicated merge (10), and
//! tombstones (E) COMPOUND, not merely work in isolation.
//!
//! The workload is the benchmark corpus's design-iteration model on the
//! CHT wedge: two concurrent agents iterate a design whose per-edit DAG
//! sizes and certifiable skip sets come from the recorded edit traces
//! (`fs-benchmark`), with a seeded fraction of candidate designs being
//! π-equivalent re-visits of tombstoned dead ends. COSTS ARE MODELED
//! UNITS from the corpus (real API calls, modeled physics): the loop
//! MECHANICS are measured; wall-clock physics lands with the wedge's
//! real solvers.
//!
//! The measurement (review round 3): isolated speedups per proposal and
//! the composed loop, over N seeded replays, asserting composed >
//! max(isolated) by a stated margin with across-replay variance
//! reported — plus laundering-across-the-loop (an estimated speculation
//! result is never upgraded or re-rooted anywhere downstream), whole-loop
//! determinism (G5: identical input-bound report and evidence commitments),
//! and a cancellation storm (G4: a mid-loop cancel leaves consistent event
//! and evidence prefixes with no residue).
#![cfg(feature = "flywheel-e2e")]

use std::collections::BTreeMap;

use fs_evidence::{Color, IntervalOp, NumericalCertificate, NumericalKind, verified_from};
use fs_geom::sheaf_merge::{BranchState, MergeOutcome, three_way_merge};
use fs_geom::sheaf_repair::SheafSkeleton;
use fs_ledger::tombstone::{Descriptor, ExplorationVerdict, TombstoneIndex};
use fs_ledger::{
    ColorGraph, ColorNode, ContentHash, PolicyDecision, SourceOrigin, SourceOriginRequest,
    SourceOriginVerifier, WaiverGrant, hash_bytes,
};
use fs_qty::{Dims, QtyAny};
use fs_recompute::{NodeRecord, ParamValue, SkipDecision, Store};
use fs_spececo::{Decision, ProposerTelemetry, SolveRecord, decide};

/// Which proposals are switched ON for a run (a feature-toggle matrix
/// by design: the harness measures every on/off combination).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct LoopConfig {
    /// Proposal 9: certified speculation (warm starts).
    pub speculation: bool,
    /// Proposal 2: incremental recompute (skips).
    pub recompute: bool,
    /// Proposal 10: sheaf-adjudicated merge (parallel credit).
    pub merge: bool,
    /// Proposal E: tombstone gate (dead candidates blocked).
    pub tombstones: bool,
    /// Cancel after this many stage transitions (G4 storm), if any.
    pub cancel_after_stages: Option<usize>,
}

impl LoopConfig {
    /// Everything off: the baseline.
    #[must_use]
    pub fn baseline() -> LoopConfig {
        LoopConfig {
            speculation: false,
            recompute: false,
            merge: false,
            tombstones: false,
            cancel_after_stages: None,
        }
    }

    /// Everything on: the closed loop.
    #[must_use]
    pub fn composed() -> LoopConfig {
        LoopConfig {
            speculation: true,
            recompute: true,
            merge: true,
            tombstones: true,
            cancel_after_stages: None,
        }
    }
}

/// One run's report — the whole flywheel's telemetry in one trace.
#[derive(Debug)]
pub struct LoopReport {
    /// Exact proposal/cancellation configuration used for this replay.
    pub config: LoopConfig,
    /// Number of design iterations requested before any cancellation.
    pub requested_iterations: usize,
    /// Logical workload seed (independent of worker identity).
    pub seed: u64,
    /// Total modeled cost (wall-analog units).
    pub total_cost: f64,
    /// Iterations completed, including terminal block/dead outcomes
    /// (`== requested_iterations` unless cancelled).
    pub iterations: usize,
    /// True when a G4 cancel fired mid-loop.
    pub cancelled: bool,
    /// Structured stage events (the trace; hashable for G5).
    pub events: Vec<String>,
    /// Speculation accept rate (0 when off).
    pub accept_rate: f64,
    /// Ops skipped by recompute.
    pub skips: usize,
    /// Merges resolved / conflicted.
    pub merges: (usize, usize),
    /// Candidates blocked by the tombstone gate.
    pub tombstone_blocks: usize,
    /// Append-only evidence lineage retained by this run.
    pub color_graph: ColorGraph,
    /// Current headline node in [`Self::color_graph`].
    pub headline_node: u64,
}

impl LoopReport {
    /// The current scientifically admissible headline color, resolved from
    /// its retained lineage. A waived headline is intentionally unavailable
    /// through this accessor.
    #[must_use]
    pub fn headline(&self) -> Option<&Color> {
        self.color_graph
            .node(self.headline_node)
            .and_then(ColorNode::scientific_color)
    }

    /// The G5 trace hash over every report field and every retained evidence
    /// field. The encoding is versioned, domain-separated, and length-prefixed.
    #[must_use]
    pub fn trace_hash(&self) -> String {
        let mut bytes = vec![2_u8];
        push_field(&mut bytes, b"frankensim/fs-flywheel-e2e/loop-report/v2");
        bytes.push(u8::from(self.config.speculation));
        bytes.push(u8::from(self.config.recompute));
        bytes.push(u8::from(self.config.merge));
        bytes.push(u8::from(self.config.tombstones));
        match self.config.cancel_after_stages {
            Some(stage) => {
                bytes.push(1);
                push_usize(&mut bytes, stage);
            }
            None => bytes.push(0),
        }
        push_usize(&mut bytes, self.requested_iterations);
        bytes.extend_from_slice(&self.seed.to_le_bytes());
        bytes.extend_from_slice(&self.total_cost.to_bits().to_le_bytes());
        push_usize(&mut bytes, self.iterations);
        bytes.push(u8::from(self.cancelled));
        push_len(&mut bytes, self.events.len());
        for event in &self.events {
            push_field(&mut bytes, event.as_bytes());
        }
        bytes.extend_from_slice(&self.accept_rate.to_bits().to_le_bytes());
        push_usize(&mut bytes, self.skips);
        push_usize(&mut bytes, self.merges.0);
        push_usize(&mut bytes, self.merges.1);
        push_usize(&mut bytes, self.tombstone_blocks);
        bytes.extend_from_slice(&self.headline_node.to_le_bytes());

        push_len(&mut bytes, self.color_graph.nodes().len());
        for node in self.color_graph.nodes() {
            push_color_node(&mut bytes, node);
        }
        push_len(&mut bytes, self.color_graph.rows().len());
        for row in self.color_graph.rows() {
            push_field(&mut bytes, row.as_bytes());
        }
        hash_bytes(&bytes).to_hex()
    }
}

impl PartialEq for LoopReport {
    fn eq(&self, other: &Self) -> bool {
        self.config == other.config
            && self.requested_iterations == other.requested_iterations
            && self.seed == other.seed
            && self.total_cost.to_bits() == other.total_cost.to_bits()
            && self.iterations == other.iterations
            && self.cancelled == other.cancelled
            && self.events == other.events
            && self.accept_rate.to_bits() == other.accept_rate.to_bits()
            && self.skips == other.skips
            && self.merges == other.merges
            && self.tombstone_blocks == other.tombstone_blocks
            && self.headline_node == other.headline_node
            && color_graphs_equal(&self.color_graph, &other.color_graph)
    }
}

const MODELED_SOURCE_NODE: &str = "modeled-baseline-qoi";
const MODELED_SOURCE_PRODUCER: &str = "fs-flywheel-e2e/modeled-baseline-v1";

#[derive(Debug, Clone, Copy)]
struct ModeledSourceVerifier;

impl SourceOriginVerifier for ModeledSourceVerifier {
    fn verify(&self, request: &SourceOriginRequest<'_>) -> PolicyDecision {
        let certificate = modeled_source_certificate();
        let color = verified_from(&certificate)
            .expect("the fixed modeled baseline certificate is a valid enclosure");
        let origin = modeled_source_origin(certificate);
        let expected = SourceOriginRequest::new(MODELED_SOURCE_NODE, &color, &origin);
        let policy = hash_bytes(b"frankensim/fs-flywheel-e2e/modeled-source-policy/v1");
        if request.canonical_bytes() == expected.canonical_bytes() {
            PolicyDecision::accept(policy)
        } else {
            PolicyDecision::reject(policy)
        }
    }
}

fn modeled_source_certificate() -> NumericalCertificate {
    NumericalCertificate::enclosure(0.0, 1e-9)
}

fn modeled_source_origin(certificate: NumericalCertificate) -> SourceOrigin {
    let mut artifact = Vec::new();
    push_field(
        &mut artifact,
        b"frankensim/fs-flywheel-e2e/modeled-source-certificate/v1",
    );
    push_numerical_certificate(&mut artifact, certificate);
    SourceOrigin::Certificate {
        producer: MODELED_SOURCE_PRODUCER.to_string(),
        certificate_hash: hash_bytes(&artifact),
        certificate,
    }
}

fn push_len(out: &mut Vec<u8>, len: usize) {
    let len = u64::try_from(len).expect("a Rust allocation length fits in u64");
    out.extend_from_slice(&len.to_le_bytes());
}

fn push_usize(out: &mut Vec<u8>, value: usize) {
    let value = u64::try_from(value).expect("a Rust usize fits in u64");
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_field(out: &mut Vec<u8>, field: &[u8]) {
    push_len(out, field.len());
    out.extend_from_slice(field);
}

fn push_optional_hash(out: &mut Vec<u8>, hash: Option<ContentHash>) {
    match hash {
        Some(hash) => {
            out.push(1);
            out.extend_from_slice(hash.as_bytes());
        }
        None => out.push(0),
    }
}

fn push_interval_op(out: &mut Vec<u8>, operation: Option<IntervalOp>) {
    out.push(match operation {
        None => 0,
        Some(IntervalOp::Add) => 1,
        Some(IntervalOp::Mul) => 2,
        Some(IntervalOp::Hull) => 3,
    });
}

fn push_numerical_certificate(out: &mut Vec<u8>, certificate: NumericalCertificate) {
    out.push(match certificate.kind {
        NumericalKind::Exact => 0,
        NumericalKind::Enclosure => 1,
        NumericalKind::Estimate => 2,
        NumericalKind::NoClaim => 3,
    });
    out.extend_from_slice(&certificate.lo.to_bits().to_le_bytes());
    out.extend_from_slice(&certificate.hi.to_bits().to_le_bytes());
}

fn push_source_origin(out: &mut Vec<u8>, origin: Option<&SourceOrigin>) {
    match origin {
        None => out.push(0),
        Some(SourceOrigin::Certificate {
            producer,
            certificate_hash,
            certificate,
        }) => {
            out.push(1);
            push_field(out, producer.as_bytes());
            out.extend_from_slice(certificate_hash.as_bytes());
            push_numerical_certificate(out, *certificate);
        }
        Some(SourceOrigin::Anchoring {
            dataset_id,
            content_hash,
            regime,
        }) => {
            out.push(2);
            push_field(out, dataset_id.as_bytes());
            out.extend_from_slice(content_hash.as_bytes());
            push_len(out, regime.bounds().len());
            for (axis, (lo, hi)) in regime.bounds() {
                push_field(out, axis.as_bytes());
                out.extend_from_slice(&lo.to_bits().to_le_bytes());
                out.extend_from_slice(&hi.to_bits().to_le_bytes());
            }
        }
    }
}

fn push_waiver(out: &mut Vec<u8>, waiver: Option<&fs_ledger::Waiver>) {
    match waiver {
        None => out.push(0),
        Some(waiver) => {
            out.push(1);
            push_field(out, waiver.id.as_bytes());
            push_field(out, waiver.signer.as_bytes());
            push_field(out, waiver.reason.as_bytes());
        }
    }
}

fn push_grant(out: &mut Vec<u8>, grant: &WaiverGrant) {
    push_waiver(out, Some(&grant.annotation));
    push_field(out, grant.key_id.as_bytes());
    push_field(out, grant.scope.as_bytes());
    push_field(out, grant.node_name.as_bytes());
    push_field(out, &grant.claimed_color);
    push_len(out, grant.parent_hashes.len());
    for hash in &grant.parent_hashes {
        out.extend_from_slice(hash.as_bytes());
    }
    out.extend_from_slice(&grant.expires_day.to_le_bytes());
    push_field(out, &grant.signature);
}

fn push_optional_grant(out: &mut Vec<u8>, grant: Option<&WaiverGrant>) {
    match grant {
        None => out.push(0),
        Some(grant) => {
            out.push(1);
            push_grant(out, grant);
        }
    }
}

fn push_color_node(out: &mut Vec<u8>, node: &ColorNode) {
    out.extend_from_slice(&node.id().to_le_bytes());
    push_field(out, node.name().as_bytes());
    push_field(out, &node.declared_color_unverified().canonical_bytes());
    push_len(out, node.parents().len());
    for parent in node.parents() {
        out.extend_from_slice(&parent.to_le_bytes());
    }
    push_interval_op(out, node.operation());
    push_len(out, node.demotions().len());
    for demotion in node.demotions() {
        push_usize(out, demotion.parent_index());
        out.extend_from_slice(&demotion.parent_id().to_le_bytes());
        push_field(out, demotion.reason().dataset.as_bytes());
        push_field(out, demotion.reason().axis.as_bytes());
        out.extend_from_slice(&demotion.reason().value.to_bits().to_le_bytes());
    }
    push_source_origin(out, node.origin());
    push_optional_hash(out, node.origin_policy_fingerprint());
    push_waiver(out, node.waiver());
    push_optional_grant(out, node.grant());
    push_optional_hash(out, node.waiver_policy_fingerprint());
    match node.waiver_admission_day() {
        Some(day) => {
            out.push(1);
            out.extend_from_slice(&day.to_le_bytes());
        }
        None => out.push(0),
    }
    push_len(out, node.waiver_dependencies().len());
    for dependency in node.waiver_dependencies() {
        out.extend_from_slice(&dependency.authorizing_node().to_le_bytes());
        push_interval_op(out, dependency.operation());
        push_grant(out, dependency.grant());
        out.extend_from_slice(dependency.policy_fingerprint().as_bytes());
        out.extend_from_slice(&dependency.admission_day().to_le_bytes());
    }
    out.extend_from_slice(node.hash().as_bytes());
}

fn origins_equal(a: Option<&SourceOrigin>, b: Option<&SourceOrigin>) -> bool {
    let mut a_bytes = Vec::new();
    let mut b_bytes = Vec::new();
    push_source_origin(&mut a_bytes, a);
    push_source_origin(&mut b_bytes, b);
    a_bytes == b_bytes
}

fn demotions_equal(a: &ColorNode, b: &ColorNode) -> bool {
    a.demotions().len() == b.demotions().len()
        && a.demotions().iter().zip(b.demotions()).all(|(a, b)| {
            a.parent_index() == b.parent_index()
                && a.parent_id() == b.parent_id()
                && a.reason().dataset == b.reason().dataset
                && a.reason().axis == b.reason().axis
                && a.reason().value.to_bits() == b.reason().value.to_bits()
        })
}

fn color_nodes_equal(a: &ColorNode, b: &ColorNode) -> bool {
    a.id() == b.id()
        && a.name() == b.name()
        && a.declared_color_unverified().canonical_bytes()
            == b.declared_color_unverified().canonical_bytes()
        && a.parents() == b.parents()
        && a.operation() == b.operation()
        && demotions_equal(a, b)
        && origins_equal(a.origin(), b.origin())
        && a.origin_policy_fingerprint() == b.origin_policy_fingerprint()
        && a.waiver() == b.waiver()
        && a.grant() == b.grant()
        && a.waiver_policy_fingerprint() == b.waiver_policy_fingerprint()
        && a.waiver_admission_day() == b.waiver_admission_day()
        && a.waiver_dependencies() == b.waiver_dependencies()
        && a.hash() == b.hash()
}

fn color_graphs_equal(a: &ColorGraph, b: &ColorGraph) -> bool {
    a.rows() == b.rows()
        && a.nodes().len() == b.nodes().len()
        && a.nodes()
            .iter()
            .zip(b.nodes())
            .all(|(a, b)| color_nodes_equal(a, b))
}

#[derive(Debug, Clone, Copy)]
enum DrawDomain {
    CandidateVelocity,
    SpeculationDecision,
    MergeTaint,
    MergeXGauge,
    MergeYGauge,
}

impl DrawDomain {
    const fn tag(self) -> u64 {
        match self {
            Self::CandidateVelocity => 0x8d58_6c94_1ec8_4a37,
            Self::SpeculationDecision => 0x4f7a_0d35_25bb_f3e1,
            Self::MergeTaint => 0xd316_2a6e_92f1_75cb,
            Self::MergeXGauge => 0x63c9_b804_4a17_e25d,
            Self::MergeYGauge => 0xa279_5f1b_dc34_0986,
        }
    }
}

fn mix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn logical_draw(
    seed: u64,
    domain: DrawDomain,
    iteration: usize,
    agent: usize,
    operation: usize,
) -> f64 {
    let iteration = u64::try_from(iteration).expect("a Rust iteration index fits in u64");
    let agent = u64::try_from(agent).expect("a Rust agent index fits in u64");
    let operation = u64::try_from(operation).expect("a Rust operation index fits in u64");
    let mut value = mix64(seed ^ domain.tag());
    value = mix64(value ^ iteration);
    value = mix64(value ^ agent.wrapping_mul(0xd6e8_feb8_6659_fd93));
    value = mix64(value ^ operation.wrapping_mul(0xa5a3_56e4_7f9c_2d17));
    ((value >> 11) as f64) / (1u64 << 53) as f64
}

fn completed_branch_cost(config: &LoopConfig, branch_costs: &[f64; 2], completed: usize) -> f64 {
    let completed = completed.min(branch_costs.len());
    if config.merge {
        branch_costs[..completed]
            .iter()
            .copied()
            .fold(0.0, f64::max)
    } else {
        branch_costs[..completed].iter().sum()
    }
}

fn wedge_descriptor(name: &str, velocity: f64, scale: f64) -> Descriptor {
    let mut params = BTreeMap::new();
    params.insert(
        "velocity".to_string(),
        QtyAny::new(velocity, Dims([1, 0, -1, 0, 0, 0])),
    );
    params.insert(
        "length".to_string(),
        QtyAny::new(scale, Dims([1, 0, 0, 0, 0, 0])),
    );
    params.insert(
        "viscosity".to_string(),
        QtyAny::new(1.8e-5, Dims([2, 0, -1, 0, 0, 0])),
    );
    Descriptor {
        name: name.to_string(),
        params,
    }
}

/// A 3-patch merge skeleton (the wedge's chart layout stand-in).
fn merge_skeleton() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// Run the loop: `iterations` design steps by two concurrent agents on
/// the corpus's first edit trace, with the configured proposals live.
/// Deterministic in `seed` (G5).
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn run_loop(config: &LoopConfig, iterations: usize, seed: u64) -> LoopReport {
    let trace = fs_benchmark::edit_traces()[0];
    let ops_per_iter = trace.total_ops;
    let skippable = trace.correct_skips;
    let mut events = Vec::new();
    let mut total_cost = 0.0f64;
    let mut stages = 0usize;
    let mut cancelled = false;

    let mut tombstones = TombstoneIndex::new();
    // Pre-seed the graveyard: three dead designs at known velocities.
    for v in [10.0, 20.0, 40.0] {
        tombstones.record_falsification_kill(
            wedge_descriptor("cht-wedge bracket", v, 0.1),
            "{\"kind\":\"tombstone\"}",
            vec!["estimated".to_string()],
            50.0,
            "2026-07-08",
            "agent:corpus",
        );
    }
    let mut store = Store::new();
    let mut telemetry = ProposerTelemetry::new();
    let mut accepts = 0usize;
    let mut proposals = 0usize;
    let mut skips = 0usize;
    let mut merges_ok = 0usize;
    let mut merges_conflict = 0usize;
    let mut blocks = 0usize;
    // The modeled baseline begins at an authenticated, retained certificate
    // leaf. Every accepted estimate is appended to this SAME graph and the
    // headline advances through a derived node; no stage may re-root it.
    let baseline_certificate = modeled_source_certificate();
    let baseline_color = verified_from(&baseline_certificate)
        .expect("the fixed modeled baseline certificate is a valid enclosure");
    let mut color_graph = ColorGraph::new();
    let mut headline_node = color_graph
        .source_with_origin(
            MODELED_SOURCE_NODE,
            &baseline_color,
            modeled_source_origin(baseline_certificate),
            &ModeledSourceVerifier,
        )
        .expect("the modeled-source capability admits its exact fixed request");
    let skeleton = merge_skeleton();
    let mut done = 0usize;

    'outer: for iter in 0..iterations {
        // ---- Stage E: the tombstone gate over this iteration's candidate.
        stages += 1;
        if let Some(limit) = config.cancel_after_stages
            && stages > limit
        {
            cancelled = true;
            break 'outer;
        }
        // Every third candidate is a REVISIT of a dead design (same
        // descriptor family, same pi-neighborhood); fresh explorations
        // carry genuinely distinct descriptors and physics (a fin array
        // at a different scale — decades away in pi-space).
        let revisit = iter % 3 == 2;
        let velocity_draw = logical_draw(seed, DrawDomain::CandidateVelocity, iter, 0, 0);
        let velocity = if revisit {
            20.0 + velocity_draw * 0.4
        } else {
            100.0 + 50.0 * velocity_draw
        };
        let candidate = if revisit {
            wedge_descriptor("cht-wedge bracket", velocity, 0.1)
        } else {
            wedge_descriptor(
                &format!("cht-wedge fin-array rev{}", iter % 5),
                velocity,
                0.5,
            )
        };
        if config.tombstones {
            if let ExplorationVerdict::Blocked { .. } = tombstones.pre_exploration_check(&candidate)
            {
                blocks += 1;
                events.push(format!("iter={iter} stage=tombstone verdict=blocked"));
                done += 1;
                continue; // the whole candidate's cost is saved
            }
            events.push(format!(
                "iter={iter} stage=tombstone verdict=clear v={velocity:.3}"
            ));
        } else if revisit {
            // Without the gate the dead candidate is fully re-solved by
            // BOTH agents (and then re-discovered dead).
            total_cost += 2.0 * ops_per_iter as f64;
            events.push(format!(
                "iter={iter} stage=dead-resolve cost={ops_per_iter}"
            ));
            done += 1;
            continue;
        }

        // ---- Stages 2+9: each agent solves its branch's DAG.
        let mut branch_costs = [0.0f64; 2];
        let mut completed_branches = 0usize;
        for agent in 0..branch_costs.len() {
            stages += 1;
            if let Some(limit) = config.cancel_after_stages
                && stages > limit
            {
                total_cost += completed_branch_cost(config, &branch_costs, completed_branches);
                cancelled = true;
                break 'outer;
            }
            for op in 0..ops_per_iter {
                let record = NodeRecord {
                    op_id: format!("wedge-op-{op}"),
                    input_hashes: Vec::new(),
                    params: vec![(
                        "iter-group".to_string(),
                        // Skippable ops share params across iterations
                        // (the unchanged part of the design); the rest
                        // change every iteration.
                        ParamValue::f(if op < skippable {
                            0.0
                        } else {
                            #[allow(clippy::cast_precision_loss)]
                            {
                                (iter * 2 + agent) as f64
                            }
                        }),
                    )],
                    code_version_hash: hash_bytes(b"wedge-v1"),
                    rng_seed: 7,
                    achieved_error: 1e-8,
                    required_tolerance: 1e-6,
                };
                if config.recompute
                    && matches!(store.can_skip(&record, 1e-6), SkipDecision::Hit { .. })
                {
                    skips += 1;
                    continue; // certified skip: no cost
                }
                // Speculation: a proposer offers a warm start.
                let mut op_cost = 1.0;
                if config.speculation {
                    proposals += 1;
                    let bound =
                        if logical_draw(seed, DrawDomain::SpeculationDecision, iter, agent, op)
                            < 0.7
                        {
                            5e-7
                        } else {
                            1e-3
                        };
                    let decision = decide(bound, 1e-6);
                    let accepted = decision == Decision::AcceptOutright;
                    if accepted {
                        accepts += 1;
                        op_cost = 0.15; // verification-only cost
                        // Each accepted step retains both the proposer estimate
                        // and its modeled holdout-probe estimate. Their distinct
                        // dataset identities force a derived evidence identity
                        // before the result joins the headline, so even the first
                        // accepted step cannot later masquerade as a source leaf.
                        let proposer =
                            format!("wedge-proposer-v1/agent-{agent}/dataset-agent-{agent}");
                        let proposal_node = color_graph
                            .source(
                                &format!("proposal/iter-{iter}/agent-{agent}/op-{op}"),
                                Color::Estimated {
                                    estimator: proposer,
                                    dispersion: 0.05,
                                },
                            )
                            .expect("fixed proposer evidence is a valid estimate leaf");
                        let holdout = format!(
                            "wedge-validation-probe-v1/agent-{agent}/dataset-holdout-{}",
                            1 - agent
                        );
                        let holdout_node = color_graph
                            .source(
                                &format!("holdout/iter-{iter}/agent-{agent}/op-{op}"),
                                Color::Estimated {
                                    estimator: holdout,
                                    dispersion: 0.05,
                                },
                            )
                            .expect("fixed holdout evidence is a valid estimate leaf");
                        let accepted_node = color_graph
                            .derive(
                                &format!("accepted/iter-{iter}/agent-{agent}/op-{op}"),
                                &[proposal_node, holdout_node],
                                IntervalOp::Hull,
                                None,
                                &BTreeMap::new(),
                                None,
                            )
                            .expect("heterogeneous accepted evidence derives conservatively");
                        headline_node = color_graph
                            .derive(
                                &format!("headline/iter-{iter}/agent-{agent}/op-{op}"),
                                &[headline_node, accepted_node],
                                IntervalOp::Hull,
                                None,
                                &BTreeMap::new(),
                                None,
                            )
                            .expect("the write gate derives the conservative headline");
                    }
                    telemetry.record(&SolveRecord::new(
                        "wedge-proposer-v1",
                        "re-1e5",
                        accepted,
                        bound,
                        if accepted { 40 } else { -2 },
                    ));
                }
                branch_costs[agent] += op_cost;
                let _ = store.put(record, b"artifact");
            }
            events.push(format!(
                "iter={iter} stage=solve agent={agent} cost={:.2}",
                branch_costs[agent]
            ));
            completed_branches += 1;
        }

        // ---- Stage 10: merge the two branches.
        stages += 1;
        if let Some(limit) = config.cancel_after_stages
            && stages > limit
        {
            total_cost += completed_branch_cost(config, &branch_costs, completed_branches);
            cancelled = true;
            break 'outer;
        }
        if config.merge {
            let base = vec![0.0; 3];
            // Gauge-style concurrent edits (occasionally a circulation
            // taint that genuinely conflicts).
            let taint = logical_draw(seed, DrawDomain::MergeTaint, iter, 0, 0) < 0.1;
            let x = BranchState {
                provenance: format!("agent-x@{iter}"),
                mismatch: if taint {
                    skeleton.d1t(&[0.05])
                } else {
                    skeleton.d0(&[
                        0.0,
                        0.01 * logical_draw(seed, DrawDomain::MergeXGauge, iter, 0, 0),
                        0.0,
                    ])
                },
                assignments: BTreeMap::new(),
            };
            let y = BranchState {
                provenance: format!("agent-y@{iter}"),
                mismatch: skeleton.d0(&[
                    0.0,
                    0.0,
                    0.01 * logical_draw(seed, DrawDomain::MergeYGauge, iter, 0, 0),
                ]),
                assignments: BTreeMap::new(),
            };
            match three_way_merge(&skeleton, &base, &x, &y, None, 1e-6, 1e-6) {
                MergeOutcome::Resolved { .. } | MergeOutcome::Trivial { .. } => {
                    merges_ok += 1;
                    // Parallel credit: wall time is the max branch.
                    total_cost += branch_costs[0].max(branch_costs[1]);
                    events.push(format!("iter={iter} stage=merge verdict=resolved"));
                }
                _ => {
                    merges_conflict += 1;
                    // Conflict: serialize + redo the cheaper branch.
                    total_cost +=
                        branch_costs[0] + branch_costs[1] + branch_costs[0].min(branch_costs[1]);
                    events.push(format!("iter={iter} stage=merge verdict=conflict"));
                }
            }
        } else {
            // No merge machinery: agents serialize.
            total_cost += branch_costs[0] + branch_costs[1];
            events.push(format!("iter={iter} stage=serialize"));
        }
        done += 1;
    }

    #[allow(clippy::cast_precision_loss)]
    let accept_rate = if proposals == 0 {
        0.0
    } else {
        accepts as f64 / proposals as f64
    };
    LoopReport {
        config: *config,
        requested_iterations: iterations,
        seed,
        total_cost,
        iterations: done,
        cancelled,
        events,
        accept_rate,
        skips,
        merges: (merges_ok, merges_conflict),
        tombstone_blocks: blocks,
        color_graph,
        headline_node,
    }
}

/// Isolated + composed speedups over one seed (baseline_cost / cost).
/// A zero-work comparison is defined as the neutral ratio `1.0` rather than
/// the indeterminate IEEE-754 expression `0.0 / 0.0`.
#[must_use]
pub fn speedups(iterations: usize, seed: u64) -> (BTreeMap<&'static str, f64>, f64) {
    let base = run_loop(&LoopConfig::baseline(), iterations, seed).total_cost;
    let ratio = |cost: f64| {
        if base == 0.0 && cost == 0.0 {
            1.0
        } else {
            base / cost
        }
    };
    let one = |f: fn(&mut LoopConfig)| {
        let mut c = LoopConfig::baseline();
        f(&mut c);
        ratio(run_loop(&c, iterations, seed).total_cost)
    };
    let mut isolated = BTreeMap::new();
    isolated.insert("speculation", one(|c| c.speculation = true));
    isolated.insert("recompute", one(|c| c.recompute = true));
    isolated.insert("merge", one(|c| c.merge = true));
    isolated.insert("tombstones", one(|c| c.tombstones = true));
    let composed = ratio(run_loop(&LoopConfig::composed(), iterations, seed).total_cost);
    (isolated, composed)
}
