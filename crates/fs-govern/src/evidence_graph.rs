//! Phase 0B-C descriptive support/threat graph and allocation planning.
//!
//! This module consumes the exact candidate identities defined by
//! [`crate::evidence_contract`] and the validated lane charters defined by
//! [`crate::lanes`]. It does not authenticate evidence, execute an e-process,
//! reserve resources in [`crate::lanes::PortfolioLedger`], persist a graph, or
//! mint runtime authority. Every public result here is immutable planning data
//! that a later authenticated checker/ledger adapter must verify again.

#![allow(missing_docs)]

use crate::{
    evidence_contract::{
        AUTHORITY_ALGEBRA_VERSION, AssumptionSetId, AttackEdge, AuthorityState, AuthorityStateId,
        ClaimInstance, ClaimInstanceId, CounterexampleCandidate, CounterexampleId, EvidenceId,
        EvidenceKind, EvidenceRef, SupportEdge,
    },
    lanes::{LaneCharter, ProofLaneId},
};
use fs_blake3::ContentHash;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

pub const EVIDENCE_GRAPH_VERSION: u32 = 2;
/// Canonical field tag binding graph nodes and snapshots to the authority
/// algebra they interpret. This tag is part of the v2 identity schema.
pub const EVIDENCE_GRAPH_AUTHORITY_ALGEBRA_TAG: u8 = 250;
pub const FIXED_RATE_SCALE: u32 = 1_000_000;
pub const MAX_GRAPH_NODES: usize = 1_024;
pub const MAX_GRAPH_EDGES: usize = 4_096;
pub const MAX_ALLOCATION_CANDIDATES: usize = 1_024;
pub const MAX_ALLOCATION_SELECTIONS: usize = 256;
pub const MAX_FEASIBILITY_SEARCH_STATES: usize = 100_000;
pub const MAX_GRAPH_TEXT_BYTES: usize = 4_096;

pub const GRAPH_NODE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.evidence-graph-node.v2";
pub const GRAPH_SNAPSHOT_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.evidence-graph-snapshot.v2";
pub const ANYTIME_ACCOUNTING_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.anytime-accounting-candidate.v2";
pub const ALLOCATION_CANDIDATE_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.evidence-allocation-candidate.v2";
pub const ALLOCATION_POLICY_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.evidence-allocation-policy.v2";
pub const ALLOCATION_DECISION_IDENTITY_DOMAIN: &str =
    "frankensim.fs-govern.evidence-allocation-decision.v2";

fn push_field(out: &mut Vec<u8>, tag: u8, bytes: &[u8]) {
    out.push(tag);
    out.extend_from_slice(
        &u64::try_from(bytes.len())
            .expect("bounded evidence-graph field length fits u64")
            .to_le_bytes(),
    );
    out.extend_from_slice(bytes);
}

fn push_hash(out: &mut Vec<u8>, tag: u8, hash: &ContentHash) {
    push_field(out, tag, hash.as_bytes());
}

fn require_hash(what: &'static str, hash: ContentHash) -> Result<ContentHash, GraphError> {
    if hash.as_bytes().iter().all(|byte| *byte == 0) {
        Err(GraphError::MissingIdentity { what })
    } else {
        Ok(hash)
    }
}

fn canonical_text(what: &'static str, raw: &str) -> Result<String, GraphError> {
    if raw.len() > MAX_GRAPH_TEXT_BYTES {
        return Err(GraphError::TooLarge {
            what,
            observed: raw.len(),
            cap: MAX_GRAPH_TEXT_BYTES,
        });
    }
    if raw.split_whitespace().next().is_none() {
        return Err(GraphError::EmptyField { what });
    }
    Ok(raw.split_whitespace().collect::<Vec<_>>().join(" "))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphError {
    EmptyField {
        what: &'static str,
    },
    TooLarge {
        what: &'static str,
        observed: usize,
        cap: usize,
    },
    MissingIdentity {
        what: &'static str,
    },
    InvalidRate {
        observed: u32,
    },
    InvalidValue {
        what: &'static str,
    },
    DuplicateNode {
        node: GraphNodeId,
    },
    DuplicateSemanticIdentity {
        kind: &'static str,
        identity: ContentHash,
    },
    DuplicateEdge {
        kind: &'static str,
        identity: ContentHash,
    },
    UnknownEndpoint {
        kind: &'static str,
        identity: ContentHash,
    },
    EndpointMismatch {
        what: &'static str,
    },
    SelfSupport {
        claim: ClaimInstanceId,
    },
    SupportCycle,
    ArithmeticOverflow {
        what: &'static str,
    },
    SnapshotMismatch,
    ClaimLaneMismatch,
    DuplicateCandidate {
        identity: AllocationCandidateId,
    },
    AnytimeAccountingRequired {
        candidate: AllocationCandidateId,
    },
    FloorUnsatisfied {
        kind: WorkKind,
        required: u64,
        selected: u64,
    },
    DiversityFloorUnsatisfied {
        required: u32,
        selected: u32,
    },
    FeasibilitySearchLimit {
        explored: usize,
        cap: usize,
    },
    Cancelled {
        pass: PlanningPass,
    },
}

impl core::fmt::Display for GraphError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyField { what } => write!(f, "evidence-graph field `{what}` is empty"),
            Self::TooLarge {
                what,
                observed,
                cap,
            } => write!(f, "evidence-graph {what} has size {observed}, cap {cap}"),
            Self::MissingIdentity { what } => write!(f, "missing identity for {what}"),
            Self::InvalidRate { observed } => {
                write!(f, "fixed rate {observed} exceeds scale {FIXED_RATE_SCALE}")
            }
            Self::InvalidValue { what } => write!(f, "invalid evidence-graph value: {what}"),
            Self::DuplicateNode { node } => write!(f, "duplicate graph node {node}"),
            Self::DuplicateSemanticIdentity { kind, identity } => {
                write!(f, "duplicate {kind} semantic identity {identity}")
            }
            Self::DuplicateEdge { kind, identity } => {
                write!(f, "duplicate {kind} edge {identity}")
            }
            Self::UnknownEndpoint { kind, identity } => {
                write!(f, "unknown {kind} endpoint {identity}")
            }
            Self::EndpointMismatch { what } => write!(f, "graph endpoint mismatch: {what}"),
            Self::SelfSupport { claim } => write!(f, "claim {claim} cannot support itself"),
            Self::SupportCycle => f.write_str("support graph contains a circular support path"),
            Self::ArithmeticOverflow { what } => {
                write!(f, "checked evidence-graph arithmetic overflow in {what}")
            }
            Self::SnapshotMismatch => {
                f.write_str("allocation candidate binds another graph snapshot")
            }
            Self::ClaimLaneMismatch => {
                f.write_str("allocation charter does not mint the claim's exact proof lane")
            }
            Self::DuplicateCandidate { identity } => {
                write!(f, "duplicate allocation candidate {identity}")
            }
            Self::AnytimeAccountingRequired { candidate } => write!(
                f,
                "allocation candidate {candidate} lacks required anytime accounting"
            ),
            Self::FloorUnsatisfied {
                kind,
                required,
                selected,
            } => write!(
                f,
                "{} floor {required} cannot be met; selected {selected}",
                kind.code()
            ),
            Self::DiversityFloorUnsatisfied { required, selected } => write!(
                f,
                "correlation-diversity floor {required} cannot be met; selected {selected}"
            ),
            Self::FeasibilitySearchLimit { explored, cap } => write!(
                f,
                "allocation feasibility search explored {explored} states, cap {cap}"
            ),
            Self::Cancelled { pass } => {
                write!(f, "evidence planning cancelled during {}", pass.code())
            }
        }
    }
}

impl std::error::Error for GraphError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GraphNodeId(ContentHash);

impl GraphNodeId {
    #[must_use]
    pub fn as_hash(&self) -> &ContentHash {
        &self.0
    }
}

impl core::fmt::Display for GraphNodeId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphNodeKind {
    Claim {
        claim: ClaimInstanceId,
        lane: ProofLaneId,
    },
    AuthorityState {
        state: AuthorityStateId,
        claim: ClaimInstanceId,
        consequence: u64,
    },
    AssumptionSet {
        assumptions: AssumptionSetId,
        claim: ClaimInstanceId,
    },
    Evidence {
        evidence: EvidenceId,
        kind: EvidenceKind,
        claim: ClaimInstanceId,
        artifact: ContentHash,
        checker: ContentHash,
        schema_version: u32,
    },
    Checker {
        checker: ContentHash,
        claim: ClaimInstanceId,
    },
    Falsifier {
        falsifier: ContentHash,
        claim: ClaimInstanceId,
    },
    Consumer {
        consumer: ContentHash,
        claim: ClaimInstanceId,
        consequence: u64,
    },
    Counterexample {
        candidate: CounterexampleId,
        claim: ClaimInstanceId,
        evidence: EvidenceId,
    },
}

impl GraphNodeKind {
    fn tag(self) -> u8 {
        match self {
            Self::Claim { .. } => 1,
            Self::AuthorityState { .. } => 2,
            Self::AssumptionSet { .. } => 3,
            Self::Evidence { .. } => 4,
            Self::Checker { .. } => 5,
            Self::Falsifier { .. } => 6,
            Self::Consumer { .. } => 7,
            Self::Counterexample { .. } => 8,
        }
    }

    #[must_use]
    pub fn claim(self) -> ClaimInstanceId {
        match self {
            Self::Claim { claim, .. }
            | Self::AuthorityState { claim, .. }
            | Self::AssumptionSet { claim, .. }
            | Self::Evidence { claim, .. }
            | Self::Checker { claim, .. }
            | Self::Falsifier { claim, .. }
            | Self::Consumer { claim, .. }
            | Self::Counterexample { claim, .. } => claim,
        }
    }

    fn encode(self, out: &mut Vec<u8>) {
        push_field(out, 1, &[self.tag()]);
        push_hash(out, 2, self.claim().as_hash());
        match self {
            Self::Claim { lane, .. } => push_hash(out, 3, lane.as_hash()),
            Self::AuthorityState {
                state, consequence, ..
            } => {
                push_hash(out, 3, state.as_hash());
                push_field(out, 4, &consequence.to_le_bytes());
            }
            Self::AssumptionSet { assumptions, .. } => {
                push_hash(out, 3, assumptions.as_hash());
            }
            Self::Evidence {
                evidence,
                kind,
                artifact,
                checker,
                schema_version,
                ..
            } => {
                push_hash(out, 3, evidence.as_hash());
                push_field(out, 4, kind.code().as_bytes());
                push_hash(out, 5, &artifact);
                push_hash(out, 6, &checker);
                push_field(out, 7, &schema_version.to_le_bytes());
            }
            Self::Checker { checker, .. } => push_hash(out, 3, &checker),
            Self::Falsifier { falsifier, .. } => push_hash(out, 3, &falsifier),
            Self::Consumer {
                consumer,
                consequence,
                ..
            } => {
                push_hash(out, 3, &consumer);
                push_field(out, 4, &consequence.to_le_bytes());
            }
            Self::Counterexample {
                candidate,
                evidence,
                ..
            } => {
                push_hash(out, 3, candidate.as_hash());
                push_hash(out, 4, evidence.as_hash());
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GraphNode {
    kind: GraphNodeKind,
    identity: GraphNodeId,
}

impl GraphNode {
    fn from_kind(kind: GraphNodeKind) -> Self {
        let mut canonical = Vec::new();
        push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
        push_field(
            &mut canonical,
            EVIDENCE_GRAPH_AUTHORITY_ALGEBRA_TAG,
            &AUTHORITY_ALGEBRA_VERSION.to_le_bytes(),
        );
        kind.encode(&mut canonical);
        Self {
            kind,
            identity: GraphNodeId(fs_blake3::hash_domain(
                GRAPH_NODE_IDENTITY_DOMAIN,
                &canonical,
            )),
        }
    }

    #[must_use]
    pub fn claim(claim: &ClaimInstance) -> Self {
        Self::from_kind(GraphNodeKind::Claim {
            claim: claim.identity(),
            lane: claim.proof_lane(),
        })
    }

    #[must_use]
    pub fn authority(state: &AuthorityState) -> Self {
        Self::authority_with_consequence(state, 1)
            .expect("the default authority consequence is nonzero")
    }

    /// Add an authority-state node with an explicit downstream consequence.
    /// The weight is descriptive planning input, not authenticated authority.
    pub fn authority_with_consequence(
        state: &AuthorityState,
        consequence: u64,
    ) -> Result<Self, GraphError> {
        if consequence == 0 {
            return Err(GraphError::InvalidValue {
                what: "authority-state consequence must be nonzero",
            });
        }
        Ok(Self::from_kind(GraphNodeKind::AuthorityState {
            state: state.identity(),
            claim: state.claim().identity(),
            consequence,
        }))
    }

    #[must_use]
    pub fn assumptions(claim: &ClaimInstance) -> Self {
        Self::from_kind(GraphNodeKind::AssumptionSet {
            assumptions: claim.assumptions().identity(),
            claim: claim.identity(),
        })
    }

    #[must_use]
    pub fn evidence(evidence: EvidenceRef) -> Self {
        Self::from_kind(GraphNodeKind::Evidence {
            evidence: evidence.identity(),
            kind: evidence.kind(),
            claim: evidence.claim(),
            artifact: evidence.artifact(),
            checker: evidence.checker(),
            schema_version: evidence.schema_version(),
        })
    }

    pub fn checker(claim: &ClaimInstance, checker: ContentHash) -> Result<Self, GraphError> {
        Ok(Self::from_kind(GraphNodeKind::Checker {
            checker: require_hash("checker", checker)?,
            claim: claim.identity(),
        }))
    }

    pub fn falsifier(claim: &ClaimInstance, falsifier: ContentHash) -> Result<Self, GraphError> {
        Ok(Self::from_kind(GraphNodeKind::Falsifier {
            falsifier: require_hash("falsifier", falsifier)?,
            claim: claim.identity(),
        }))
    }

    pub fn consumer(
        claim: &ClaimInstance,
        consumer: ContentHash,
        consequence: u64,
    ) -> Result<Self, GraphError> {
        if consequence == 0 {
            return Err(GraphError::InvalidValue {
                what: "consumer consequence must be nonzero",
            });
        }
        Ok(Self::from_kind(GraphNodeKind::Consumer {
            consumer: require_hash("consumer", consumer)?,
            claim: claim.identity(),
            consequence,
        }))
    }

    #[must_use]
    pub fn counterexample(candidate: &CounterexampleCandidate) -> Self {
        Self::from_kind(GraphNodeKind::Counterexample {
            candidate: candidate.identity(),
            claim: candidate.target(),
            evidence: candidate.evidence(),
        })
    }

    #[must_use]
    pub fn kind(&self) -> GraphNodeKind {
        self.kind
    }

    #[must_use]
    pub fn identity(&self) -> GraphNodeId {
        self.identity
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GraphSnapshot {
    nodes: Vec<GraphNode>,
    support: Vec<SupportEdge>,
    attacks: Vec<AttackEdge>,
    claims: BTreeMap<ClaimInstanceId, ProofLaneId>,
    support_adjacency: BTreeMap<ClaimInstanceId, BTreeSet<ClaimInstanceId>>,
    authority_consequences: BTreeMap<ClaimInstanceId, u64>,
    consumers: BTreeMap<ClaimInstanceId, Vec<(GraphNodeId, u64)>>,
    direct_attack_counts: BTreeMap<ClaimInstanceId, u64>,
    attack_counts: BTreeMap<ClaimInstanceId, u64>,
    identity: ContentHash,
}

impl GraphSnapshot {
    #[allow(clippy::too_many_lines)]
    pub fn new(
        mut nodes: Vec<GraphNode>,
        mut support: Vec<SupportEdge>,
        mut attacks: Vec<AttackEdge>,
    ) -> Result<Self, GraphError> {
        if nodes.is_empty() {
            return Err(GraphError::InvalidValue {
                what: "graph snapshot must contain at least one node",
            });
        }
        if nodes.len() > MAX_GRAPH_NODES {
            return Err(GraphError::TooLarge {
                what: "nodes",
                observed: nodes.len(),
                cap: MAX_GRAPH_NODES,
            });
        }
        let edge_count = support
            .len()
            .checked_add(attacks.len())
            .ok_or(GraphError::ArithmeticOverflow { what: "edge count" })?;
        if edge_count > MAX_GRAPH_EDGES {
            return Err(GraphError::TooLarge {
                what: "edges",
                observed: edge_count,
                cap: MAX_GRAPH_EDGES,
            });
        }

        nodes.sort_by_key(GraphNode::identity);
        for pair in nodes.windows(2) {
            if pair[0].identity == pair[1].identity {
                return Err(GraphError::DuplicateNode {
                    node: pair[0].identity,
                });
            }
        }

        let mut claims = BTreeMap::new();
        let mut authorities = BTreeMap::new();
        let mut authority_consequences = BTreeMap::<ClaimInstanceId, u64>::new();
        let mut evidence_nodes = BTreeMap::new();
        let mut counterexamples = BTreeMap::new();
        let mut consumers = BTreeMap::<ClaimInstanceId, Vec<(GraphNodeId, u64)>>::new();
        let mut consumer_ids = BTreeMap::<ContentHash, ClaimInstanceId>::new();
        for node in &nodes {
            match node.kind {
                GraphNodeKind::Claim { claim, lane } => {
                    if claims.insert(claim, lane).is_some() {
                        return Err(GraphError::DuplicateSemanticIdentity {
                            kind: "claim",
                            identity: *claim.as_hash(),
                        });
                    }
                }
                GraphNodeKind::AuthorityState {
                    state,
                    claim,
                    consequence,
                } => {
                    if authorities.insert(state, claim).is_some() {
                        return Err(GraphError::DuplicateSemanticIdentity {
                            kind: "authority state",
                            identity: *state.as_hash(),
                        });
                    }
                    authority_consequences
                        .entry(claim)
                        .and_modify(|current| *current = (*current).max(consequence))
                        .or_insert(consequence);
                }
                GraphNodeKind::Evidence {
                    evidence, claim, ..
                } => {
                    if evidence_nodes.insert(evidence, claim).is_some() {
                        return Err(GraphError::DuplicateSemanticIdentity {
                            kind: "evidence",
                            identity: *evidence.as_hash(),
                        });
                    }
                }
                GraphNodeKind::Counterexample {
                    candidate,
                    claim,
                    evidence,
                } => {
                    if counterexamples
                        .insert(candidate, (claim, evidence))
                        .is_some()
                    {
                        return Err(GraphError::DuplicateSemanticIdentity {
                            kind: "counterexample",
                            identity: *candidate.as_hash(),
                        });
                    }
                }
                GraphNodeKind::Consumer {
                    consumer,
                    claim,
                    consequence,
                } => {
                    if consumer_ids.insert(consumer, claim).is_some() {
                        return Err(GraphError::DuplicateSemanticIdentity {
                            kind: "consumer",
                            identity: consumer,
                        });
                    }
                    consumers
                        .entry(claim)
                        .or_default()
                        .push((node.identity, consequence));
                }
                GraphNodeKind::AssumptionSet { .. }
                | GraphNodeKind::Checker { .. }
                | GraphNodeKind::Falsifier { .. } => {}
            }
        }

        for node in &nodes {
            if !matches!(node.kind, GraphNodeKind::Claim { .. })
                && !claims.contains_key(&node.kind.claim())
            {
                return Err(GraphError::UnknownEndpoint {
                    kind: "claim",
                    identity: *node.kind.claim().as_hash(),
                });
            }
        }
        for (claim, evidence) in counterexamples.values() {
            match evidence_nodes.get(evidence).copied() {
                None => {
                    return Err(GraphError::UnknownEndpoint {
                        kind: "counterexample evidence",
                        identity: *evidence.as_hash(),
                    });
                }
                Some(evidence_claim) if evidence_claim != *claim => {
                    return Err(GraphError::EndpointMismatch {
                        what: "counterexample evidence is not bound to the target claim",
                    });
                }
                Some(_) => {}
            }
        }

        support.sort_by_key(SupportEdge::identity);
        for pair in support.windows(2) {
            if pair[0].identity() == pair[1].identity() {
                return Err(GraphError::DuplicateEdge {
                    kind: "support",
                    identity: *pair[0].identity().as_hash(),
                });
            }
        }
        let mut support_adjacency = claims
            .keys()
            .copied()
            .map(|claim| (claim, BTreeSet::new()))
            .collect::<BTreeMap<_, _>>();
        for edge in &support {
            let Some(source_claim) = authorities.get(&edge.source()).copied() else {
                return Err(GraphError::UnknownEndpoint {
                    kind: "support source authority",
                    identity: *edge.source().as_hash(),
                });
            };
            if !claims.contains_key(&edge.target()) {
                return Err(GraphError::UnknownEndpoint {
                    kind: "support target claim",
                    identity: *edge.target().as_hash(),
                });
            }
            if claims.get(&edge.target()).copied() != Some(edge.proof_lane()) {
                return Err(GraphError::EndpointMismatch {
                    what: "support edge proof lane does not match target claim",
                });
            }
            match evidence_nodes.get(&edge.evidence()).copied() {
                None => {
                    return Err(GraphError::UnknownEndpoint {
                        kind: "support evidence",
                        identity: *edge.evidence().as_hash(),
                    });
                }
                Some(evidence_claim) if evidence_claim != edge.target() => {
                    return Err(GraphError::EndpointMismatch {
                        what: "support evidence is not bound to the target claim",
                    });
                }
                Some(_) => {}
            }
            if source_claim == edge.target() {
                return Err(GraphError::SelfSupport {
                    claim: source_claim,
                });
            }
            support_adjacency
                .entry(source_claim)
                .or_default()
                .insert(edge.target());
        }
        validate_support_dag(&support_adjacency)?;

        attacks.sort_by_key(AttackEdge::identity);
        for pair in attacks.windows(2) {
            if pair[0].identity() == pair[1].identity() {
                return Err(GraphError::DuplicateEdge {
                    kind: "attack",
                    identity: *pair[0].identity().as_hash(),
                });
            }
        }
        let mut direct_attack_counts = BTreeMap::new();
        let mut attack_counts = BTreeMap::new();
        for edge in &attacks {
            let Some((candidate_target, _)) = counterexamples.get(&edge.candidate()).copied()
            else {
                return Err(GraphError::UnknownEndpoint {
                    kind: "attack counterexample",
                    identity: *edge.candidate().as_hash(),
                });
            };
            if candidate_target != edge.target() {
                return Err(GraphError::EndpointMismatch {
                    what: "attack candidate target differs from attack edge target",
                });
            }
            if claims.get(&edge.target()).copied() != Some(edge.proof_lane()) {
                return Err(GraphError::EndpointMismatch {
                    what: "attack edge proof lane does not match target claim",
                });
            }
            match evidence_nodes.get(&edge.evidence()).copied() {
                None => {
                    return Err(GraphError::UnknownEndpoint {
                        kind: "attack evidence",
                        identity: *edge.evidence().as_hash(),
                    });
                }
                Some(evidence_claim) if evidence_claim != edge.target() => {
                    return Err(GraphError::EndpointMismatch {
                        what: "attack evidence is not bound to the target claim",
                    });
                }
                Some(_) => {}
            }
            let direct = direct_attack_counts.entry(edge.target()).or_insert(0u64);
            *direct = direct
                .checked_add(1)
                .ok_or(GraphError::ArithmeticOverflow {
                    what: "direct attack count",
                })?;
            for affected in reachable_claims(&support_adjacency, edge.target()) {
                let propagated = attack_counts.entry(affected).or_insert(0u64);
                *propagated = propagated
                    .checked_add(1)
                    .ok_or(GraphError::ArithmeticOverflow {
                        what: "propagated attack count",
                    })?;
            }
        }

        for values in consumers.values_mut() {
            values.sort_by_key(|(identity, _)| *identity);
        }
        let mut canonical = Vec::new();
        push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
        push_field(
            &mut canonical,
            EVIDENCE_GRAPH_AUTHORITY_ALGEBRA_TAG,
            &AUTHORITY_ALGEBRA_VERSION.to_le_bytes(),
        );
        for node in &nodes {
            push_hash(&mut canonical, 1, node.identity.as_hash());
        }
        for edge in &support {
            push_hash(&mut canonical, 2, edge.identity().as_hash());
        }
        for edge in &attacks {
            push_hash(&mut canonical, 3, edge.identity().as_hash());
        }
        let identity = fs_blake3::hash_domain(GRAPH_SNAPSHOT_IDENTITY_DOMAIN, &canonical);
        Ok(Self {
            nodes,
            support,
            attacks,
            claims,
            support_adjacency,
            authority_consequences,
            consumers,
            direct_attack_counts,
            attack_counts,
            identity,
        })
    }

    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub fn nodes(&self) -> &[GraphNode] {
        &self.nodes
    }

    #[must_use]
    pub fn support_edges(&self) -> &[SupportEdge] {
        &self.support
    }

    #[must_use]
    pub fn attack_edges(&self) -> &[AttackEdge] {
        &self.attacks
    }

    #[must_use]
    pub fn contains_claim(&self, claim: ClaimInstanceId) -> bool {
        self.claims.contains_key(&claim)
    }

    #[must_use]
    pub fn attack_count(&self, claim: ClaimInstanceId) -> u64 {
        self.attack_counts.get(&claim).copied().unwrap_or(0)
    }

    /// Count attacks whose target is exactly `claim`, before downstream
    /// propagation through support edges.
    #[must_use]
    pub fn direct_attack_count(&self, claim: ClaimInstanceId) -> u64 {
        self.direct_attack_counts.get(&claim).copied().unwrap_or(0)
    }

    pub fn consequence(&self, claim: ClaimInstanceId) -> Result<u64, GraphError> {
        if !self.claims.contains_key(&claim) {
            return Err(GraphError::UnknownEndpoint {
                kind: "consequence claim",
                identity: *claim.as_hash(),
            });
        }
        let reachable = reachable_claims(&self.support_adjacency, claim);
        let mut seen_consumers = BTreeSet::new();
        let mut consequence = 0u64;
        for reachable_claim in reachable {
            if let Some(weight) = self.authority_consequences.get(&reachable_claim) {
                consequence =
                    consequence
                        .checked_add(*weight)
                        .ok_or(GraphError::ArithmeticOverflow {
                            what: "reachable consequence",
                        })?;
            }
            for (consumer, weight) in self
                .consumers
                .get(&reachable_claim)
                .map_or(&[][..], Vec::as_slice)
            {
                if seen_consumers.insert(*consumer) {
                    consequence =
                        consequence
                            .checked_add(*weight)
                            .ok_or(GraphError::ArithmeticOverflow {
                                what: "reachable consequence",
                            })?;
                }
            }
        }
        Ok(consequence)
    }
}

fn reachable_claims(
    adjacency: &BTreeMap<ClaimInstanceId, BTreeSet<ClaimInstanceId>>,
    root: ClaimInstanceId,
) -> BTreeSet<ClaimInstanceId> {
    let mut reachable = BTreeSet::from([root]);
    let mut queue = VecDeque::from([root]);
    while let Some(current) = queue.pop_front() {
        if let Some(targets) = adjacency.get(&current) {
            for target in targets {
                if reachable.insert(*target) {
                    queue.push_back(*target);
                }
            }
        }
    }
    reachable
}

fn validate_support_dag(
    adjacency: &BTreeMap<ClaimInstanceId, BTreeSet<ClaimInstanceId>>,
) -> Result<(), GraphError> {
    let mut indegree = adjacency
        .keys()
        .copied()
        .map(|claim| (claim, 0usize))
        .collect::<BTreeMap<_, _>>();
    for targets in adjacency.values() {
        for target in targets {
            let degree = indegree
                .get_mut(target)
                .ok_or(GraphError::UnknownEndpoint {
                    kind: "support DAG claim",
                    identity: *target.as_hash(),
                })?;
            *degree = degree
                .checked_add(1)
                .ok_or(GraphError::ArithmeticOverflow {
                    what: "support indegree",
                })?;
        }
    }
    let mut ready = indegree
        .iter()
        .filter_map(|(claim, degree)| (*degree == 0).then_some(*claim))
        .collect::<BTreeSet<_>>();
    let mut visited = 0usize;
    while let Some(claim) = ready.pop_first() {
        visited = visited
            .checked_add(1)
            .ok_or(GraphError::ArithmeticOverflow {
                what: "support traversal",
            })?;
        if let Some(targets) = adjacency.get(&claim) {
            for target in targets {
                let degree = indegree
                    .get_mut(target)
                    .ok_or(GraphError::UnknownEndpoint {
                        kind: "support DAG claim",
                        identity: *target.as_hash(),
                    })?;
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(*target);
                }
            }
        }
    }
    if visited == adjacency.len() {
        Ok(())
    } else {
        Err(GraphError::SupportCycle)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FixedRate(u32);

impl FixedRate {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(FIXED_RATE_SCALE);

    pub const fn new(parts_per_million: u32) -> Result<Self, GraphError> {
        if parts_per_million <= FIXED_RATE_SCALE {
            Ok(Self(parts_per_million))
        } else {
            Err(GraphError::InvalidRate {
                observed: parts_per_million,
            })
        }
    }

    #[must_use]
    pub const fn parts_per_million(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn complement(self) -> Self {
        Self(FIXED_RATE_SCALE - self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DoubtProfile {
    calibrated_uncertainty: FixedRate,
    attack_coverage: FixedRate,
    independent_support: FixedRate,
    assumption_resolution: FixedRate,
}

impl DoubtProfile {
    #[must_use]
    pub const fn new(
        calibrated_uncertainty: FixedRate,
        attack_coverage: FixedRate,
        independent_support: FixedRate,
        assumption_resolution: FixedRate,
    ) -> Self {
        Self {
            calibrated_uncertainty,
            attack_coverage,
            independent_support,
            assumption_resolution,
        }
    }

    #[must_use]
    pub const fn calibrated_uncertainty(self) -> FixedRate {
        self.calibrated_uncertainty
    }

    #[must_use]
    pub const fn attack_coverage(self) -> FixedRate {
        self.attack_coverage
    }

    #[must_use]
    pub const fn independent_support(self) -> FixedRate {
        self.independent_support
    }

    #[must_use]
    pub const fn assumption_resolution(self) -> FixedRate {
        self.assumption_resolution
    }

    /// Conservative union of the four named doubt sources. Products of the
    /// remaining-confidence terms round down, so the returned doubt never
    /// rounds inward.
    #[must_use]
    pub fn combined(self) -> FixedRate {
        let risks = [
            self.calibrated_uncertainty,
            self.attack_coverage.complement(),
            self.independent_support.complement(),
            self.assumption_resolution.complement(),
        ];
        let mut remaining = u64::from(FIXED_RATE_SCALE);
        for risk in risks {
            remaining =
                remaining * u64::from(FIXED_RATE_SCALE - risk.0) / u64::from(FIXED_RATE_SCALE);
        }
        let remaining = u32::try_from(remaining).expect("fixed-rate product remains within scale");
        FixedRate(FIXED_RATE_SCALE - remaining)
    }

    /// Combine caller-supplied doubt with graph-visible, unadjudicated attack
    /// candidates. The graph has no authenticated adjudication state, so any
    /// reachable attack is treated as maximum doubt instead of being silently
    /// discounted by a favorable caller-supplied coverage value.
    #[must_use]
    pub fn combined_with_unadjudicated_attacks(self, attack_count: u64) -> FixedRate {
        if attack_count == 0 {
            self.combined()
        } else {
            FixedRate::ONE
        }
    }

    fn encode(self, out: &mut Vec<u8>, tag: u8) {
        let mut value = Vec::new();
        for rate in [
            self.calibrated_uncertainty,
            self.attack_coverage,
            self.independent_support,
            self.assumption_resolution,
        ] {
            value.extend_from_slice(&rate.0.to_le_bytes());
        }
        push_field(out, tag, &value);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnytimeAccountingCandidate {
    method: String,
    schema_version: u32,
    observations: u64,
    state_root: ContentHash,
    evidence_artifact: ContentHash,
    identity: ContentHash,
}

impl AnytimeAccountingCandidate {
    pub fn new(
        method: &str,
        schema_version: u32,
        observations: u64,
        state_root: ContentHash,
        evidence_artifact: ContentHash,
    ) -> Result<Self, GraphError> {
        if schema_version == 0 {
            return Err(GraphError::InvalidValue {
                what: "anytime-accounting schema version must be nonzero",
            });
        }
        if observations == 0 {
            return Err(GraphError::InvalidValue {
                what: "anytime-accounting observation count must be nonzero",
            });
        }
        let method = canonical_text("anytime-accounting method", method)?;
        let state_root = require_hash("anytime-accounting state", state_root)?;
        let evidence_artifact = require_hash("anytime-accounting evidence", evidence_artifact)?;
        let mut canonical = Vec::new();
        push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
        push_field(&mut canonical, 1, method.as_bytes());
        push_field(&mut canonical, 2, &schema_version.to_le_bytes());
        push_field(&mut canonical, 3, &observations.to_le_bytes());
        push_hash(&mut canonical, 4, &state_root);
        push_hash(&mut canonical, 5, &evidence_artifact);
        let identity = fs_blake3::hash_domain(ANYTIME_ACCOUNTING_IDENTITY_DOMAIN, &canonical);
        Ok(Self {
            method,
            schema_version,
            observations,
            state_root,
            evidence_artifact,
            identity,
        })
    }

    #[must_use]
    pub fn method(&self) -> &str {
        &self.method
    }

    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub fn observations(&self) -> u64 {
        self.observations
    }

    #[must_use]
    pub fn state_root(&self) -> ContentHash {
        self.state_root
    }

    #[must_use]
    pub fn evidence_artifact(&self) -> ContentHash {
        self.evidence_artifact
    }

    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum WorkKind {
    Falsification,
    IndependentCheck,
    Holdout,
    Exploration,
    Reproduction,
}

impl WorkKind {
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Falsification => "falsification",
            Self::IndependentCheck => "independent-check",
            Self::Holdout => "holdout",
            Self::Exploration => "exploration",
            Self::Reproduction => "reproduction",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Falsification => 1,
            Self::IndependentCheck => 2,
            Self::Holdout => 3,
            Self::Exploration => 4,
            Self::Reproduction => 5,
        }
    }

    const fn index(self) -> usize {
        match self {
            Self::Falsification => 0,
            Self::IndependentCheck => 1,
            Self::Holdout => 2,
            Self::Exploration => 3,
            Self::Reproduction => 4,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AllocationCandidateId(ContentHash);

impl AllocationCandidateId {
    #[must_use]
    pub fn as_hash(&self) -> &ContentHash {
        &self.0
    }
}

impl core::fmt::Display for AllocationCandidateId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllocationCandidate {
    graph: ContentHash,
    claim: ClaimInstanceId,
    lane: ProofLaneId,
    independence_class: ContentHash,
    kind: WorkKind,
    work_artifact: ContentHash,
    cost: u64,
    utility_weight: u32,
    doubt: DoubtProfile,
    prior: ContentHash,
    correlation_class: ContentHash,
    anytime: Option<AnytimeAccountingCandidate>,
    identity: AllocationCandidateId,
}

impl AllocationCandidate {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        graph: &GraphSnapshot,
        claim: &ClaimInstance,
        charter: &LaneCharter,
        kind: WorkKind,
        work_artifact: ContentHash,
        cost: u64,
        utility_weight: u32,
        doubt: DoubtProfile,
        prior: ContentHash,
        correlation_class: ContentHash,
        anytime: Option<AnytimeAccountingCandidate>,
    ) -> Result<Self, GraphError> {
        if !graph.contains_claim(claim.identity()) {
            return Err(GraphError::UnknownEndpoint {
                kind: "allocation claim",
                identity: *claim.identity().as_hash(),
            });
        }
        let lane = charter.lane_id();
        if lane != claim.proof_lane() {
            return Err(GraphError::ClaimLaneMismatch);
        }
        if cost == 0 {
            return Err(GraphError::InvalidValue {
                what: "allocation candidate cost must be nonzero",
            });
        }
        if utility_weight == 0 {
            return Err(GraphError::InvalidValue {
                what: "allocation utility weight must be nonzero",
            });
        }
        let work_artifact = require_hash("allocation work artifact", work_artifact)?;
        let prior = require_hash("allocation prior", prior)?;
        let correlation_class = require_hash("allocation correlation class", correlation_class)?;
        let independence_class = charter.independence_class_id();
        let mut canonical = Vec::new();
        push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
        push_hash(&mut canonical, 1, &graph.identity);
        push_hash(&mut canonical, 2, claim.identity().as_hash());
        push_hash(&mut canonical, 3, lane.as_hash());
        push_hash(&mut canonical, 4, &independence_class);
        push_field(&mut canonical, 5, &[kind.tag()]);
        push_hash(&mut canonical, 6, &work_artifact);
        push_field(&mut canonical, 7, &cost.to_le_bytes());
        push_field(&mut canonical, 8, &utility_weight.to_le_bytes());
        doubt.encode(&mut canonical, 9);
        push_hash(&mut canonical, 10, &prior);
        push_hash(&mut canonical, 11, &correlation_class);
        match &anytime {
            None => push_field(&mut canonical, 12, &[0]),
            Some(accounting) => {
                push_field(&mut canonical, 12, &[1]);
                push_hash(&mut canonical, 13, &accounting.identity);
            }
        }
        let identity = AllocationCandidateId(fs_blake3::hash_domain(
            ALLOCATION_CANDIDATE_IDENTITY_DOMAIN,
            &canonical,
        ));
        Ok(Self {
            graph: graph.identity,
            claim: claim.identity(),
            lane,
            independence_class,
            kind,
            work_artifact,
            cost,
            utility_weight,
            doubt,
            prior,
            correlation_class,
            anytime,
            identity,
        })
    }

    #[must_use]
    pub fn graph(&self) -> ContentHash {
        self.graph
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
    }

    #[must_use]
    pub fn independence_class(&self) -> ContentHash {
        self.independence_class
    }

    #[must_use]
    pub fn kind(&self) -> WorkKind {
        self.kind
    }

    #[must_use]
    pub fn work_artifact(&self) -> ContentHash {
        self.work_artifact
    }

    #[must_use]
    pub fn cost(&self) -> u64 {
        self.cost
    }

    #[must_use]
    pub fn utility_weight(&self) -> u32 {
        self.utility_weight
    }

    #[must_use]
    pub fn doubt(&self) -> DoubtProfile {
        self.doubt
    }

    #[must_use]
    pub fn prior(&self) -> ContentHash {
        self.prior
    }

    #[must_use]
    pub fn correlation_class(&self) -> ContentHash {
        self.correlation_class
    }

    #[must_use]
    pub fn anytime(&self) -> Option<&AnytimeAccountingCandidate> {
        self.anytime.as_ref()
    }

    #[must_use]
    pub fn identity(&self) -> AllocationCandidateId {
        self.identity
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AllocationFloors {
    pub falsification: u64,
    pub independent_check: u64,
    pub holdout: u64,
    pub exploration: u64,
}

impl AllocationFloors {
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            falsification: 0,
            independent_check: 0,
            holdout: 0,
            exploration: 0,
        }
    }

    fn get(self, kind: WorkKind) -> u64 {
        match kind {
            WorkKind::Falsification => self.falsification,
            WorkKind::IndependentCheck => self.independent_check,
            WorkKind::Holdout => self.holdout,
            WorkKind::Exploration => self.exploration,
            WorkKind::Reproduction => 0,
        }
    }

    fn checked_sum(self) -> Option<u64> {
        self.falsification
            .checked_add(self.independent_check)?
            .checked_add(self.holdout)?
            .checked_add(self.exploration)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllocationPolicy {
    total_budget: u64,
    no_action_reserve: u64,
    max_selections: u32,
    min_correlation_classes: u32,
    require_anytime_accounting: bool,
    floors: AllocationFloors,
    utility_model: ContentHash,
    sensitivity_artifact: ContentHash,
    identity: ContentHash,
}

impl AllocationPolicy {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        total_budget: u64,
        no_action_reserve: u64,
        max_selections: u32,
        min_correlation_classes: u32,
        require_anytime_accounting: bool,
        floors: AllocationFloors,
        utility_model: ContentHash,
        sensitivity_artifact: ContentHash,
    ) -> Result<Self, GraphError> {
        if total_budget == 0 {
            return Err(GraphError::InvalidValue {
                what: "allocation budget must be nonzero",
            });
        }
        if no_action_reserve > total_budget {
            return Err(GraphError::InvalidValue {
                what: "no-action reserve exceeds total budget",
            });
        }
        if max_selections == 0
            || usize::try_from(max_selections).unwrap_or(usize::MAX) > MAX_ALLOCATION_SELECTIONS
        {
            return Err(GraphError::InvalidValue {
                what: "allocation selection cap is outside the supported range",
            });
        }
        if min_correlation_classes > max_selections {
            return Err(GraphError::InvalidValue {
                what: "correlation diversity floor exceeds selection cap",
            });
        }
        let allocatable = total_budget - no_action_reserve;
        let floor_sum = floors.checked_sum().ok_or(GraphError::ArithmeticOverflow {
            what: "floor budget",
        })?;
        if floor_sum > allocatable {
            return Err(GraphError::InvalidValue {
                what: "allocation floors exceed budget after no-action reserve",
            });
        }
        let utility_model = require_hash("allocation utility model", utility_model)?;
        let sensitivity_artifact =
            require_hash("allocation sensitivity artifact", sensitivity_artifact)?;
        let mut canonical = Vec::new();
        push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
        push_field(&mut canonical, 1, &total_budget.to_le_bytes());
        push_field(&mut canonical, 2, &no_action_reserve.to_le_bytes());
        push_field(&mut canonical, 3, &max_selections.to_le_bytes());
        push_field(&mut canonical, 4, &min_correlation_classes.to_le_bytes());
        push_field(&mut canonical, 5, &[u8::from(require_anytime_accounting)]);
        for (tag, value) in [
            (6, floors.falsification),
            (7, floors.independent_check),
            (8, floors.holdout),
            (9, floors.exploration),
        ] {
            push_field(&mut canonical, tag, &value.to_le_bytes());
        }
        push_hash(&mut canonical, 10, &utility_model);
        push_hash(&mut canonical, 11, &sensitivity_artifact);
        let identity = fs_blake3::hash_domain(ALLOCATION_POLICY_IDENTITY_DOMAIN, &canonical);
        Ok(Self {
            total_budget,
            no_action_reserve,
            max_selections,
            min_correlation_classes,
            require_anytime_accounting,
            floors,
            utility_model,
            sensitivity_artifact,
            identity,
        })
    }

    #[must_use]
    pub fn total_budget(&self) -> u64 {
        self.total_budget
    }

    #[must_use]
    pub fn no_action_reserve(&self) -> u64 {
        self.no_action_reserve
    }

    #[must_use]
    pub fn max_selections(&self) -> u32 {
        self.max_selections
    }

    #[must_use]
    pub fn min_correlation_classes(&self) -> u32 {
        self.min_correlation_classes
    }

    #[must_use]
    pub fn require_anytime_accounting(&self) -> bool {
        self.require_anytime_accounting
    }

    #[must_use]
    pub fn floors(&self) -> AllocationFloors {
        self.floors
    }

    #[must_use]
    pub fn utility_model(&self) -> ContentHash {
        self.utility_model
    }

    #[must_use]
    pub fn sensitivity_artifact(&self) -> ContentHash {
        self.sensitivity_artifact
    }

    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }

    fn allocatable(&self) -> u64 {
        self.total_budget - self.no_action_reserve
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AllocationScore {
    consequence: u64,
    doubt: FixedRate,
    attack_count: u64,
    utility_weight: u32,
    weighted_value: u128,
}

impl AllocationScore {
    fn for_candidate(
        graph: &GraphSnapshot,
        candidate: &AllocationCandidate,
    ) -> Result<Self, GraphError> {
        let consequence = graph.consequence(candidate.claim)?;
        let attack_count = graph.attack_count(candidate.claim);
        let doubt = candidate
            .doubt
            .combined_with_unadjudicated_attacks(attack_count);
        let weighted_value = u128::from(consequence)
            .checked_mul(u128::from(doubt.0))
            .and_then(|value| value.checked_mul(u128::from(candidate.utility_weight)))
            .ok_or(GraphError::ArithmeticOverflow {
                what: "consequence-times-doubt score",
            })?
            / u128::from(FIXED_RATE_SCALE);
        Ok(Self {
            consequence,
            doubt,
            attack_count,
            utility_weight: candidate.utility_weight,
            weighted_value,
        })
    }

    #[must_use]
    pub fn consequence(self) -> u64 {
        self.consequence
    }

    #[must_use]
    pub fn doubt(self) -> FixedRate {
        self.doubt
    }

    #[must_use]
    pub fn attack_count(self) -> u64 {
        self.attack_count
    }

    #[must_use]
    pub fn utility_weight(self) -> u32 {
        self.utility_weight
    }

    #[must_use]
    pub fn weighted_value(self) -> u128 {
        self.weighted_value
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanningPass {
    InputValidation,
    Feasibility,
    Floor(WorkKind),
    Ranking,
    Finalization,
}

impl PlanningPass {
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::InputValidation => "input-validation",
            Self::Feasibility => "feasibility",
            Self::Floor(WorkKind::Falsification) => "falsification-floor",
            Self::Floor(WorkKind::IndependentCheck) => "independent-check-floor",
            Self::Floor(WorkKind::Holdout) => "holdout-floor",
            Self::Floor(WorkKind::Exploration) => "exploration-floor",
            Self::Floor(WorkKind::Reproduction) => "reproduction-floor",
            Self::Ranking => "ranking",
            Self::Finalization => "finalization",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RankedAllocationCandidate {
    candidate: AllocationCandidateId,
    claim: ClaimInstanceId,
    kind: WorkKind,
    lane: ProofLaneId,
    independence_class: ContentHash,
    correlation_class: ContentHash,
    work_artifact: ContentHash,
    cost: u64,
    doubt_profile: DoubtProfile,
    prior: ContentHash,
    anytime_state: Option<AnytimeAccountingCandidate>,
    score: AllocationScore,
    selected: bool,
}

impl RankedAllocationCandidate {
    #[must_use]
    pub fn candidate(&self) -> AllocationCandidateId {
        self.candidate
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn kind(&self) -> WorkKind {
        self.kind
    }

    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
    }

    #[must_use]
    pub fn independence_class(&self) -> ContentHash {
        self.independence_class
    }

    #[must_use]
    pub fn correlation_class(&self) -> ContentHash {
        self.correlation_class
    }

    #[must_use]
    pub fn work_artifact(&self) -> ContentHash {
        self.work_artifact
    }

    #[must_use]
    pub fn cost(&self) -> u64 {
        self.cost
    }

    #[must_use]
    pub fn doubt_profile(&self) -> DoubtProfile {
        self.doubt_profile
    }

    #[must_use]
    pub fn prior(&self) -> ContentHash {
        self.prior
    }

    #[must_use]
    pub fn anytime_accounting(&self) -> Option<ContentHash> {
        self.anytime_state
            .as_ref()
            .map(AnytimeAccountingCandidate::identity)
    }

    #[must_use]
    pub fn anytime_state(&self) -> Option<&AnytimeAccountingCandidate> {
        self.anytime_state.as_ref()
    }

    #[must_use]
    pub fn score(&self) -> AllocationScore {
        self.score
    }

    #[must_use]
    pub fn selected(&self) -> bool {
        self.selected
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllocationReservationCandidate {
    candidate: AllocationCandidateId,
    claim: ClaimInstanceId,
    lane: ProofLaneId,
    independence_class: ContentHash,
    correlation_class: ContentHash,
    kind: WorkKind,
    work_artifact: ContentHash,
    cost: u64,
    doubt_profile: DoubtProfile,
    prior: ContentHash,
    anytime_state: Option<AnytimeAccountingCandidate>,
    score: AllocationScore,
}

impl AllocationReservationCandidate {
    #[must_use]
    pub fn candidate(&self) -> AllocationCandidateId {
        self.candidate
    }

    #[must_use]
    pub fn claim(&self) -> ClaimInstanceId {
        self.claim
    }

    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
    }

    #[must_use]
    pub fn independence_class(&self) -> ContentHash {
        self.independence_class
    }

    #[must_use]
    pub fn correlation_class(&self) -> ContentHash {
        self.correlation_class
    }

    #[must_use]
    pub fn kind(&self) -> WorkKind {
        self.kind
    }

    #[must_use]
    pub fn work_artifact(&self) -> ContentHash {
        self.work_artifact
    }

    #[must_use]
    pub fn cost(&self) -> u64 {
        self.cost
    }

    #[must_use]
    pub fn doubt_profile(&self) -> DoubtProfile {
        self.doubt_profile
    }

    #[must_use]
    pub fn prior(&self) -> ContentHash {
        self.prior
    }

    #[must_use]
    pub fn anytime_accounting(&self) -> Option<ContentHash> {
        self.anytime_state
            .as_ref()
            .map(AnytimeAccountingCandidate::identity)
    }

    #[must_use]
    pub fn anytime_state(&self) -> Option<&AnytimeAccountingCandidate> {
        self.anytime_state.as_ref()
    }

    #[must_use]
    pub fn score(&self) -> AllocationScore {
        self.score
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AllocationDecisionCandidate {
    graph: ContentHash,
    policy: ContentHash,
    utility_model: ContentHash,
    sensitivity_artifact: ContentHash,
    ranked: Vec<RankedAllocationCandidate>,
    selected: Vec<AllocationReservationCandidate>,
    used_budget: u64,
    no_action_reserve: u64,
    unused_allocatable_budget: u64,
    unallocated_budget: u64,
    no_action_selected: bool,
    identity: ContentHash,
}

impl AllocationDecisionCandidate {
    #[must_use]
    pub fn graph(&self) -> ContentHash {
        self.graph
    }

    #[must_use]
    pub fn policy(&self) -> ContentHash {
        self.policy
    }

    #[must_use]
    pub fn utility_model(&self) -> ContentHash {
        self.utility_model
    }

    #[must_use]
    pub fn sensitivity_artifact(&self) -> ContentHash {
        self.sensitivity_artifact
    }

    #[must_use]
    pub fn ranked(&self) -> &[RankedAllocationCandidate] {
        &self.ranked
    }

    #[must_use]
    pub fn selected(&self) -> &[AllocationReservationCandidate] {
        &self.selected
    }

    #[must_use]
    pub fn used_budget(&self) -> u64 {
        self.used_budget
    }

    #[must_use]
    pub fn no_action_reserve(&self) -> u64 {
        self.no_action_reserve
    }

    #[must_use]
    pub fn unused_allocatable_budget(&self) -> u64 {
        self.unused_allocatable_budget
    }

    #[must_use]
    pub fn unallocated_budget(&self) -> u64 {
        self.unallocated_budget
    }

    #[must_use]
    pub fn no_action_selected(&self) -> bool {
        self.no_action_selected
    }

    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }
}

#[derive(Clone, Default)]
struct SelectionState {
    selected: BTreeSet<usize>,
    lanes: BTreeSet<ProofLaneId>,
    independence_classes: BTreeSet<ContentHash>,
    correlation_classes: BTreeSet<ContentHash>,
    kind_spend: [u64; 5],
    used: u64,
}

impl SelectionState {
    fn can_select(&self, candidate: &AllocationCandidate, policy: &AllocationPolicy) -> bool {
        if self.selected.len()
            >= usize::try_from(policy.max_selections).unwrap_or(MAX_ALLOCATION_SELECTIONS)
        {
            return false;
        }
        if self.lanes.contains(&candidate.lane)
            || self
                .independence_classes
                .contains(&candidate.independence_class)
            || self
                .correlation_classes
                .contains(&candidate.correlation_class)
        {
            return false;
        }
        self.used
            .checked_add(candidate.cost)
            .is_some_and(|needed| needed <= policy.allocatable())
    }

    fn select(&mut self, index: usize, candidate: &AllocationCandidate) -> Result<(), GraphError> {
        self.used =
            self.used
                .checked_add(candidate.cost)
                .ok_or(GraphError::ArithmeticOverflow {
                    what: "selected allocation budget",
                })?;
        let spend = &mut self.kind_spend[candidate.kind.index()];
        *spend = spend
            .checked_add(candidate.cost)
            .ok_or(GraphError::ArithmeticOverflow {
                what: "allocation floor spend",
            })?;
        self.selected.insert(index);
        self.lanes.insert(candidate.lane);
        self.independence_classes
            .insert(candidate.independence_class);
        self.correlation_classes.insert(candidate.correlation_class);
        Ok(())
    }

    fn selected_for_kind(&self, kind: WorkKind) -> u64 {
        self.kind_spend[kind.index()]
    }

    fn meets_requirements(&self, policy: &AllocationPolicy) -> bool {
        [
            WorkKind::Falsification,
            WorkKind::IndependentCheck,
            WorkKind::Holdout,
            WorkKind::Exploration,
        ]
        .into_iter()
        .all(|kind| self.selected_for_kind(kind) >= policy.floors.get(kind))
            && u32::try_from(self.correlation_classes.len()).unwrap_or(u32::MAX)
                >= policy.min_correlation_classes
    }

    fn progress_key(&self, policy: &AllocationPolicy) -> (u64, u32, usize) {
        let floor_progress = [
            WorkKind::Falsification,
            WorkKind::IndependentCheck,
            WorkKind::Holdout,
            WorkKind::Exploration,
        ]
        .into_iter()
        .map(|kind| self.selected_for_kind(kind).min(policy.floors.get(kind)))
        .fold(0u64, u64::saturating_add);
        let diversity = u32::try_from(self.correlation_classes.len())
            .unwrap_or(u32::MAX)
            .min(policy.min_correlation_classes);
        (floor_progress, diversity, self.selected.len())
    }
}

fn precheck_raw_feasibility(
    policy: &AllocationPolicy,
    candidates: &[AllocationCandidate],
) -> Result<(), GraphError> {
    for kind in [
        WorkKind::Falsification,
        WorkKind::IndependentCheck,
        WorkKind::Holdout,
        WorkKind::Exploration,
    ] {
        let required = policy.floors.get(kind);
        if required == 0 {
            continue;
        }
        let available = candidates
            .iter()
            .filter(|candidate| candidate.kind == kind)
            .fold(0u64, |sum, candidate| {
                sum.saturating_add(candidate.cost).min(required)
            });
        if available < required {
            return Err(GraphError::FloorUnsatisfied {
                kind,
                required,
                selected: available,
            });
        }
    }
    let available_classes = candidates
        .iter()
        .map(|candidate| candidate.correlation_class)
        .collect::<BTreeSet<_>>()
        .len();
    let available_classes = u32::try_from(available_classes).unwrap_or(u32::MAX);
    if available_classes < policy.min_correlation_classes {
        return Err(GraphError::DiversityFloorUnsatisfied {
            required: policy.min_correlation_classes,
            selected: available_classes,
        });
    }
    Ok(())
}

fn select_feasible_state<F>(
    policy: &AllocationPolicy,
    candidates: &[AllocationCandidate],
    scored: &[(usize, AllocationScore)],
    cancelled: &mut F,
) -> Result<SelectionState, GraphError>
where
    F: FnMut(PlanningPass) -> bool,
{
    precheck_raw_feasibility(policy, candidates)?;
    let mut stack = vec![(0usize, SelectionState::default())];
    let mut best = SelectionState::default();
    let mut explored = 0usize;
    while let Some((position, state)) = stack.pop() {
        if explored == MAX_FEASIBILITY_SEARCH_STATES {
            return Err(GraphError::FeasibilitySearchLimit {
                explored,
                cap: MAX_FEASIBILITY_SEARCH_STATES,
            });
        }
        explored += 1;
        if cancelled(PlanningPass::Feasibility) {
            return Err(GraphError::Cancelled {
                pass: PlanningPass::Feasibility,
            });
        }
        if state.progress_key(policy) > best.progress_key(policy) {
            best.clone_from(&state);
        }
        if state.meets_requirements(policy) {
            return Ok(state);
        }
        let Some((index, _)) = scored.get(position) else {
            continue;
        };
        let candidate = &candidates[*index];
        let included = if state.can_select(candidate, policy) {
            let mut included = state.clone();
            included.select(*index, candidate)?;
            Some(included)
        } else {
            None
        };
        stack.push((position + 1, state));
        if let Some(included) = included {
            // LIFO makes the higher-ranked include branch deterministic first.
            stack.push((position + 1, included));
        }
    }

    for kind in [
        WorkKind::Falsification,
        WorkKind::IndependentCheck,
        WorkKind::Holdout,
        WorkKind::Exploration,
    ] {
        let required = policy.floors.get(kind);
        let selected = best.selected_for_kind(kind);
        if selected < required {
            return Err(GraphError::FloorUnsatisfied {
                kind,
                required,
                selected,
            });
        }
    }
    Err(GraphError::DiversityFloorUnsatisfied {
        required: policy.min_correlation_classes,
        selected: u32::try_from(best.correlation_classes.len()).unwrap_or(u32::MAX),
    })
}

pub fn plan_allocations(
    graph: &GraphSnapshot,
    policy: &AllocationPolicy,
    candidates: &[AllocationCandidate],
) -> Result<AllocationDecisionCandidate, GraphError> {
    plan_allocations_with_cancel(graph, policy, candidates, |_| false)
}

#[allow(clippy::too_many_lines)]
pub fn plan_allocations_with_cancel<F>(
    graph: &GraphSnapshot,
    policy: &AllocationPolicy,
    candidates: &[AllocationCandidate],
    mut cancelled: F,
) -> Result<AllocationDecisionCandidate, GraphError>
where
    F: FnMut(PlanningPass) -> bool,
{
    if cancelled(PlanningPass::InputValidation) {
        return Err(GraphError::Cancelled {
            pass: PlanningPass::InputValidation,
        });
    }
    if candidates.len() > MAX_ALLOCATION_CANDIDATES {
        return Err(GraphError::TooLarge {
            what: "allocation candidates",
            observed: candidates.len(),
            cap: MAX_ALLOCATION_CANDIDATES,
        });
    }
    let mut seen = BTreeSet::new();
    let mut scored = Vec::with_capacity(candidates.len());
    for (index, candidate) in candidates.iter().enumerate() {
        if candidate.graph != graph.identity {
            return Err(GraphError::SnapshotMismatch);
        }
        if !seen.insert(candidate.identity) {
            return Err(GraphError::DuplicateCandidate {
                identity: candidate.identity,
            });
        }
        if policy.require_anytime_accounting && candidate.anytime.is_none() {
            return Err(GraphError::AnytimeAccountingRequired {
                candidate: candidate.identity,
            });
        }
        scored.push((index, AllocationScore::for_candidate(graph, candidate)?));
        if cancelled(PlanningPass::InputValidation) {
            return Err(GraphError::Cancelled {
                pass: PlanningPass::InputValidation,
            });
        }
    }
    scored.sort_by(|(left_index, left_score), (right_index, right_score)| {
        right_score
            .weighted_value
            .cmp(&left_score.weighted_value)
            .then_with(|| {
                candidates[*left_index]
                    .identity
                    .cmp(&candidates[*right_index].identity)
            })
    });

    for kind in [
        WorkKind::Falsification,
        WorkKind::IndependentCheck,
        WorkKind::Holdout,
        WorkKind::Exploration,
    ] {
        if policy.floors.get(kind) > 0 && cancelled(PlanningPass::Floor(kind)) {
            return Err(GraphError::Cancelled {
                pass: PlanningPass::Floor(kind),
            });
        }
    }
    let mut state = select_feasible_state(policy, candidates, &scored, &mut cancelled)?;

    for (index, score) in &scored {
        if score.weighted_value == 0 {
            continue;
        }
        let candidate = &candidates[*index];
        if state.can_select(candidate, policy) {
            state.select(*index, candidate)?;
        }
        if cancelled(PlanningPass::Ranking) {
            return Err(GraphError::Cancelled {
                pass: PlanningPass::Ranking,
            });
        }
    }

    if cancelled(PlanningPass::Finalization) {
        return Err(GraphError::Cancelled {
            pass: PlanningPass::Finalization,
        });
    }

    let mut ranked = Vec::with_capacity(scored.len());
    let mut selected = Vec::with_capacity(state.selected.len());
    for (index, score) in scored {
        let candidate = &candidates[index];
        let is_selected = state.selected.contains(&index);
        let anytime_state = candidate.anytime.clone();
        ranked.push(RankedAllocationCandidate {
            candidate: candidate.identity,
            claim: candidate.claim,
            kind: candidate.kind,
            lane: candidate.lane,
            independence_class: candidate.independence_class,
            correlation_class: candidate.correlation_class,
            work_artifact: candidate.work_artifact,
            cost: candidate.cost,
            doubt_profile: candidate.doubt,
            prior: candidate.prior,
            anytime_state: anytime_state.clone(),
            score,
            selected: is_selected,
        });
        if is_selected {
            selected.push(AllocationReservationCandidate {
                candidate: candidate.identity,
                claim: candidate.claim,
                lane: candidate.lane,
                independence_class: candidate.independence_class,
                correlation_class: candidate.correlation_class,
                kind: candidate.kind,
                work_artifact: candidate.work_artifact,
                cost: candidate.cost,
                doubt_profile: candidate.doubt,
                prior: candidate.prior,
                anytime_state,
                score,
            });
        }
    }
    let unused_allocatable_budget =
        policy
            .allocatable()
            .checked_sub(state.used)
            .ok_or(GraphError::ArithmeticOverflow {
                what: "unused allocatable budget",
            })?;
    let unallocated_budget = policy
        .no_action_reserve
        .checked_add(unused_allocatable_budget)
        .ok_or(GraphError::ArithmeticOverflow {
            what: "unallocated budget",
        })?;
    let no_action_selected = selected.is_empty();
    let mut canonical = Vec::new();
    push_field(&mut canonical, 0, &EVIDENCE_GRAPH_VERSION.to_le_bytes());
    push_hash(&mut canonical, 1, &graph.identity);
    push_hash(&mut canonical, 2, &policy.identity);
    for row in &ranked {
        push_hash(&mut canonical, 3, row.candidate.as_hash());
        push_field(&mut canonical, 4, &row.score.weighted_value.to_le_bytes());
        push_field(&mut canonical, 5, &[u8::from(row.selected)]);
        push_hash(&mut canonical, 6, &row.prior);
        match row.anytime_accounting() {
            None => push_field(&mut canonical, 7, &[0]),
            Some(accounting) => {
                push_field(&mut canonical, 7, &[1]);
                push_hash(&mut canonical, 8, &accounting);
            }
        }
    }
    push_field(&mut canonical, 9, &state.used.to_le_bytes());
    push_field(&mut canonical, 10, &policy.no_action_reserve.to_le_bytes());
    push_field(&mut canonical, 11, &unused_allocatable_budget.to_le_bytes());
    push_field(&mut canonical, 12, &unallocated_budget.to_le_bytes());
    push_field(&mut canonical, 13, &[u8::from(no_action_selected)]);
    push_hash(&mut canonical, 14, &policy.utility_model);
    push_hash(&mut canonical, 15, &policy.sensitivity_artifact);
    let identity = fs_blake3::hash_domain(ALLOCATION_DECISION_IDENTITY_DOMAIN, &canonical);
    if cancelled(PlanningPass::Finalization) {
        return Err(GraphError::Cancelled {
            pass: PlanningPass::Finalization,
        });
    }
    Ok(AllocationDecisionCandidate {
        graph: graph.identity,
        policy: policy.identity,
        utility_model: policy.utility_model,
        sensitivity_artifact: policy.sensitivity_artifact,
        ranked,
        selected,
        used_budget: state.used,
        no_action_reserve: policy.no_action_reserve,
        unused_allocatable_budget,
        unallocated_budget,
        no_action_selected,
        identity,
    })
}
