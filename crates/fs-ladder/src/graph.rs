//! Context- and QoI-specific fidelity graphs.
//!
//! A directed edge is an evidenced statement that its target is more
//! informative than its source in a declared context. Cost is deliberately
//! orthogonal to that direction: an informative target may be cheaper, and a
//! costly target acquires no authority merely from expense.

use crate::{Ladder, Rung, Transfer};
use fs_blake3::{ContentHash, hash_domain};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::error::Error;
use std::fmt;

/// Canonical transport and identity schema version.
pub const FIDELITY_GRAPH_SCHEMA_VERSION: u16 = 1;
/// Maximum model nodes admitted to one graph.
pub const MAX_GRAPH_NODES: usize = 4_096;
/// Maximum evidence edges admitted to one graph.
pub const MAX_GRAPH_EDGES: usize = 32_768;
/// Maximum disjunctive clauses admitted to one context predicate set.
pub const MAX_CONTEXT_CLAUSES: usize = 1_024;
/// Maximum regime axes admitted to one predicate or query point.
pub const MAX_REGIME_AXES: usize = 256;
/// Maximum bytes admitted for a semantic name.
pub const MAX_NAME_BYTES: usize = 128;
/// Maximum UTF-8 bytes retained for one legacy rung note.
pub const MAX_NOTE_BYTES: usize = 4_096;
/// Maximum canonical bytes accepted by the decoder.
pub const MAX_CANONICAL_BYTES: usize = 16 * 1024 * 1024;

const GRAPH_IDENTITY_DOMAIN: &str = "frankensim.fs-ladder.fidelity-graph.v1";
const EDGE_IDENTITY_DOMAIN: &str = "frankensim.fs-ladder.fidelity-edge.v1";
const LEGACY_MODEL_DOMAIN: &str = "frankensim.fs-ladder.legacy-model.v1";
const LEGACY_CARD_DOMAIN: &str = "frankensim.fs-ladder.legacy-model-card-ref.v1";
const LEGACY_TRANSFER_DOMAIN: &str = "frankensim.fs-ladder.legacy-transfer-ref.v1";

macro_rules! typed_hash {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ContentHash);

        impl $name {
            /// Construct from the exact external or domain-separated root.
            #[must_use]
            pub const fn new(hash: ContentHash) -> Self {
                Self(hash)
            }

            /// The underlying 32-byte root.
            #[must_use]
            pub const fn hash(self) -> ContentHash {
                self.0
            }

            /// Raw identity bytes.
            #[must_use]
            pub fn as_bytes(&self) -> &[u8; 32] {
                self.0.as_bytes()
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}

typed_hash!(
    /// Exact semantic identity of a model implementation/configuration.
    ModelId
);
typed_hash!(
    /// Exact content identity of the model card governing a node.
    ModelCardRef
);
typed_hash!(
    /// Exact identity of an `fs-plan` cost-model artifact.
    CostModelRef
);
typed_hash!(
    /// Exact identity of an `fs-evidence` discrepancy-model artifact.
    DiscrepancyModelRef
);
typed_hash!(
    /// Exact identity of a transfer implementation/configuration.
    TransferRef
);
typed_hash!(
    /// Exact identity of one graph edge.
    EdgeId
);
typed_hash!(
    /// Exact identity of a complete canonical fidelity graph.
    FidelityGraphId
);
typed_hash!(
    /// Exact receipt for a query-time cost/discrepancy evaluation.
    QueryEvidenceRef
);

macro_rules! semantic_name {
    ($(#[$meta:meta])* $name:ident, $field:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(String);

        impl $name {
            /// Parse one bounded semantic name.
            ///
            /// Names use visible ASCII alphanumerics plus `.`, `_`, `-`, `:`,
            /// and `/`; whitespace and confusable Unicode are refused.
            pub fn new(value: impl Into<String>) -> Result<Self, GraphError> {
                let value = value.into();
                validate_name($field, &value)?;
                Ok(Self(value))
            }

            /// Borrow the canonical name.
            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }
    };
}

semantic_name!(
    /// Typed quantity-of-interest identity.
    QoiId,
    "qoi"
);
semantic_name!(
    /// Typed regime-axis identity.
    RegimeAxis,
    "regime axis"
);

fn validate_name(field: &'static str, value: &str) -> Result<(), GraphError> {
    if value.is_empty()
        || value.len() > MAX_NAME_BYTES
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"._-:/".contains(&byte))
    {
        return Err(GraphError::InvalidName {
            field,
            value: value.to_string(),
        });
    }
    Ok(())
}

fn validate_note(value: &str) -> Result<(), GraphError> {
    if value.len() > MAX_NOTE_BYTES || value.chars().any(char::is_control) {
        return Err(GraphError::InvalidName {
            field: "legacy rung note",
            value: value.to_string(),
        });
    }
    Ok(())
}

fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

/// Closed finite interval on one regime axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedInterval {
    lower: f64,
    upper: f64,
}

impl ClosedInterval {
    /// Construct `[lower, upper]`.
    pub fn new(lower: f64, upper: f64) -> Result<Self, GraphError> {
        if !lower.is_finite() || !upper.is_finite() || lower > upper {
            return Err(GraphError::InvalidInterval { lower, upper });
        }
        Ok(Self {
            lower: canonical_zero(lower),
            upper: canonical_zero(upper),
        })
    }

    /// Inclusive lower endpoint.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Inclusive upper endpoint.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    fn contains(self, value: f64) -> bool {
        value >= self.lower && value <= self.upper
    }
}

/// QoI selector used by a context clause.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum QoiSelector {
    /// Match every QoI.
    Any,
    /// Match exactly one QoI.
    Exact(QoiId),
}

/// One conjunctive context clause.
///
/// Clauses are ORed by [`ContextPredicateSet`]. Inside a clause, the QoI
/// selector and every regime interval must match.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextClause {
    qoi: QoiSelector,
    axes: BTreeMap<RegimeAxis, ClosedInterval>,
}

impl ContextClause {
    /// Construct a clause from a QoI selector and bounded regime intervals.
    pub fn new(
        qoi: QoiSelector,
        axes: impl IntoIterator<Item = (RegimeAxis, ClosedInterval)>,
    ) -> Result<Self, GraphError> {
        let mut canonical = BTreeMap::new();
        for (axis, interval) in axes {
            if canonical.insert(axis.clone(), interval).is_some() {
                return Err(GraphError::DuplicateRegimeAxis(axis));
            }
            if canonical.len() > MAX_REGIME_AXES {
                return Err(GraphError::LimitExceeded {
                    what: "context regime axes",
                    limit: MAX_REGIME_AXES,
                });
            }
        }
        Ok(Self {
            qoi,
            axes: canonical,
        })
    }

    /// A universal clause.
    #[must_use]
    pub fn universal() -> Self {
        Self {
            qoi: QoiSelector::Any,
            axes: BTreeMap::new(),
        }
    }

    /// QoI selector.
    #[must_use]
    pub const fn qoi(&self) -> &QoiSelector {
        &self.qoi
    }

    /// Canonically ordered regime constraints.
    #[must_use]
    pub const fn axes(&self) -> &BTreeMap<RegimeAxis, ClosedInterval> {
        &self.axes
    }

    fn specificity(&self) -> u32 {
        u32::from(matches!(self.qoi, QoiSelector::Exact(_))) + self.axes.len() as u32
    }

    fn matches(&self, context: &QueryContext) -> bool {
        let qoi_matches = match &self.qoi {
            QoiSelector::Any => true,
            QoiSelector::Exact(qoi) => qoi == context.qoi(),
        };
        qoi_matches
            && self.axes.iter().all(|(axis, interval)| {
                context
                    .regime()
                    .get(axis)
                    .is_some_and(|value| interval.contains(*value))
            })
    }
}

/// Canonical OR-set of context clauses.
///
/// Empty means "unknown/no evidenced context", not universal.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextPredicateSet {
    clauses: Vec<ContextClause>,
}

impl ContextPredicateSet {
    /// Canonicalize an OR-set. Input order and exact duplicates are
    /// nonsemantic.
    pub fn new(clauses: impl IntoIterator<Item = ContextClause>) -> Result<Self, GraphError> {
        let mut keyed = Vec::new();
        for clause in clauses {
            if keyed.len() >= MAX_CONTEXT_CLAUSES {
                return Err(GraphError::LimitExceeded {
                    what: "context clauses",
                    limit: MAX_CONTEXT_CLAUSES,
                });
            }
            keyed.push((clause_bytes(&clause), clause));
        }
        keyed.sort_by(|left, right| left.0.cmp(&right.0));
        keyed.dedup_by(|left, right| left.0 == right.0);
        Ok(Self {
            clauses: keyed.into_iter().map(|(_, clause)| clause).collect(),
        })
    }

    /// Explicit unknown predicate set.
    #[must_use]
    pub const fn unknown() -> Self {
        Self {
            clauses: Vec::new(),
        }
    }

    /// A predicate set matching every finite query context.
    #[must_use]
    pub fn universal() -> Self {
        Self {
            clauses: vec![ContextClause::universal()],
        }
    }

    /// Canonical clauses.
    #[must_use]
    pub fn clauses(&self) -> &[ContextClause] {
        &self.clauses
    }

    /// Whether no evidenced context exists.
    #[must_use]
    pub fn is_unknown(&self) -> bool {
        self.clauses.is_empty()
    }

    fn matched_specificity(&self, context: &QueryContext) -> Option<u32> {
        self.clauses
            .iter()
            .filter(|clause| clause.matches(context))
            .map(ContextClause::specificity)
            .max()
    }
}

/// Contexts in which an edge comparison has evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidityDomain(ContextPredicateSet);

impl ValidityDomain {
    /// Wrap a canonical context set.
    #[must_use]
    pub const fn new(predicates: ContextPredicateSet) -> Self {
        Self(predicates)
    }

    /// No validated context.
    #[must_use]
    pub const fn unknown() -> Self {
        Self(ContextPredicateSet::unknown())
    }

    /// Valid for every represented context.
    #[must_use]
    pub fn universal() -> Self {
        Self(ContextPredicateSet::universal())
    }

    /// Borrow the predicates.
    #[must_use]
    pub const fn predicates(&self) -> &ContextPredicateSet {
        &self.0
    }
}

/// Contexts in which the edge target is evidenced as more informative.
#[derive(Debug, Clone, PartialEq)]
pub struct Informativeness(ContextPredicateSet);

impl Informativeness {
    /// Wrap a canonical context set.
    #[must_use]
    pub const fn new(predicates: ContextPredicateSet) -> Self {
        Self(predicates)
    }

    /// Explicitly unknown informativeness.
    #[must_use]
    pub const fn unknown() -> Self {
        Self(ContextPredicateSet::unknown())
    }

    /// Total-order predicate used only by lossless legacy-ladder embedding.
    #[must_use]
    pub fn legacy_total_order() -> Self {
        Self(ContextPredicateSet::universal())
    }

    /// Borrow the predicates.
    #[must_use]
    pub const fn predicates(&self) -> &ContextPredicateSet {
        &self.0
    }
}

/// Exact query context supplied to cost/discrepancy resolvers.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryContext {
    qoi: QoiId,
    regime: BTreeMap<RegimeAxis, f64>,
    problem_size: u64,
    budget_s: f64,
    max_relative_discrepancy: f64,
}

impl QueryContext {
    /// Construct one context. Axis order is nonsemantic; duplicates refuse.
    pub fn new(
        qoi: QoiId,
        regime: impl IntoIterator<Item = (RegimeAxis, f64)>,
        problem_size: u64,
        budget_s: f64,
        max_relative_discrepancy: f64,
    ) -> Result<Self, GraphError> {
        if !budget_s.is_finite() || budget_s < 0.0 {
            return Err(GraphError::InvalidNumber {
                field: "query budget seconds",
                value: budget_s,
            });
        }
        if !max_relative_discrepancy.is_finite() || max_relative_discrepancy < 0.0 {
            return Err(GraphError::InvalidNumber {
                field: "maximum relative discrepancy",
                value: max_relative_discrepancy,
            });
        }
        let mut canonical = BTreeMap::new();
        for (axis, value) in regime {
            if !value.is_finite() {
                return Err(GraphError::InvalidNumber {
                    field: "query regime coordinate",
                    value,
                });
            }
            if canonical
                .insert(axis.clone(), canonical_zero(value))
                .is_some()
            {
                return Err(GraphError::DuplicateRegimeAxis(axis));
            }
            if canonical.len() > MAX_REGIME_AXES {
                return Err(GraphError::LimitExceeded {
                    what: "query regime axes",
                    limit: MAX_REGIME_AXES,
                });
            }
        }
        Ok(Self {
            qoi,
            regime: canonical,
            problem_size,
            budget_s: canonical_zero(budget_s),
            max_relative_discrepancy: canonical_zero(max_relative_discrepancy),
        })
    }

    /// Requested QoI.
    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    /// Exact regime coordinates.
    #[must_use]
    pub const fn regime(&self) -> &BTreeMap<RegimeAxis, f64> {
        &self.regime
    }

    /// Declared problem size interpreted by the referenced cost model.
    #[must_use]
    pub const fn problem_size(&self) -> u64 {
        self.problem_size
    }

    /// Maximum admitted predicted model cost in seconds.
    #[must_use]
    pub const fn budget_s(&self) -> f64 {
        self.budget_s
    }

    /// Largest admitted relative discrepancy for `cheapest_adequate`.
    #[must_use]
    pub const fn max_relative_discrepancy(&self) -> f64 {
        self.max_relative_discrepancy
    }
}

/// Legacy rung metadata retained only for a lossless v1-ladder embedding.
#[derive(Debug, Clone, PartialEq, Eq)]
struct LegacyRung {
    index: u32,
    relative_cost_bits: u64,
    note: String,
}

/// One model node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FidelityNode {
    id: ModelId,
    card: ModelCardRef,
    label: String,
    legacy: Option<LegacyRung>,
}

impl FidelityNode {
    /// Construct a native graph node.
    pub fn new(
        id: ModelId,
        card: ModelCardRef,
        label: impl Into<String>,
    ) -> Result<Self, GraphError> {
        let label = label.into();
        validate_name("model label", &label)?;
        Ok(Self {
            id,
            card,
            label,
            legacy: None,
        })
    }

    /// Model identity.
    #[must_use]
    pub const fn id(&self) -> ModelId {
        self.id
    }

    /// Governing model-card identity.
    #[must_use]
    pub const fn card(&self) -> ModelCardRef {
        self.card
    }

    /// Human/agent-readable bounded label.
    #[must_use]
    pub fn label(&self) -> &str {
        &self.label
    }
}

/// Cost relation carried by an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CostRelationRef {
    /// Native reference to a queryable `fs-plan` cost-model artifact.
    Model(CostModelRef),
    /// Exact legacy scalar retained during v1 ladder embedding.
    ///
    /// This remains advisory and never masquerades as a fitted cost model.
    LegacyRelativeCost(u64),
}

impl CostRelationRef {
    /// Construct an exact legacy scalar reference from a positive finite value.
    fn legacy(value: f64) -> Result<Self, GraphError> {
        if !value.is_finite() || value <= 0.0 {
            return Err(GraphError::InvalidNumber {
                field: "legacy relative cost",
                value,
            });
        }
        Ok(Self::LegacyRelativeCost(value.to_bits()))
    }

    /// Recover a retained legacy scalar when present.
    #[must_use]
    pub fn legacy_relative_cost(self) -> Option<f64> {
        match self {
            Self::Model(_) => None,
            Self::LegacyRelativeCost(bits) => Some(f64::from_bits(bits)),
        }
    }
}

/// Pairwise discrepancy evidence carried by an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscrepancyReference {
    /// Native reference to a retained discrepancy-model artifact.
    Model(DiscrepancyModelRef),
    /// Legacy ladder declared an order but carried no discrepancy evidence.
    UnknownLegacy,
}

/// Directed contextual comparison between two models.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityEdge {
    id: EdgeId,
    source: ModelId,
    target: ModelId,
    cost: CostRelationRef,
    discrepancy: DiscrepancyReference,
    transfer: TransferRef,
    validity: ValidityDomain,
    informativeness: Informativeness,
}

impl FidelityEdge {
    /// Construct a native evidence edge.
    pub fn new(
        source: ModelId,
        target: ModelId,
        cost: CostModelRef,
        discrepancy: DiscrepancyModelRef,
        transfer: TransferRef,
        validity: ValidityDomain,
        informativeness: Informativeness,
    ) -> Result<Self, GraphError> {
        Self::build(
            source,
            target,
            CostRelationRef::Model(cost),
            DiscrepancyReference::Model(discrepancy),
            transfer,
            validity,
            informativeness,
        )
    }

    fn build(
        source: ModelId,
        target: ModelId,
        cost: CostRelationRef,
        discrepancy: DiscrepancyReference,
        transfer: TransferRef,
        validity: ValidityDomain,
        informativeness: Informativeness,
    ) -> Result<Self, GraphError> {
        if source == target {
            return Err(GraphError::SelfLoop(source));
        }
        let mut edge = Self {
            id: EdgeId::new(ContentHash([0; 32])),
            source,
            target,
            cost,
            discrepancy,
            transfer,
            validity,
            informativeness,
        };
        edge.id = EdgeId::new(hash_domain(EDGE_IDENTITY_DOMAIN, &edge_body_bytes(&edge)));
        Ok(edge)
    }

    /// Edge identity.
    #[must_use]
    pub const fn id(&self) -> EdgeId {
        self.id
    }

    /// Less-informative endpoint under matching conditions.
    #[must_use]
    pub const fn source(&self) -> ModelId {
        self.source
    }

    /// More-informative endpoint under matching conditions.
    #[must_use]
    pub const fn target(&self) -> ModelId {
        self.target
    }

    /// Cost relation reference.
    #[must_use]
    pub const fn cost(&self) -> CostRelationRef {
        self.cost
    }

    /// Pairwise discrepancy reference.
    #[must_use]
    pub const fn discrepancy(&self) -> DiscrepancyReference {
        self.discrepancy
    }

    /// Transfer implementation/configuration reference.
    #[must_use]
    pub const fn transfer(&self) -> TransferRef {
        self.transfer
    }

    /// Context in which the comparison has evidence.
    #[must_use]
    pub const fn validity(&self) -> &ValidityDomain {
        &self.validity
    }

    /// Context in which the target is more informative.
    #[must_use]
    pub const fn informativeness(&self) -> &Informativeness {
        &self.informativeness
    }
}

/// Versioned fidelity graph.
#[derive(Debug, Clone, PartialEq)]
pub struct FidelityGraph {
    name: String,
    legacy_kernel: Option<String>,
    nodes: BTreeMap<ModelId, FidelityNode>,
    edges: BTreeMap<EdgeId, FidelityEdge>,
}

impl FidelityGraph {
    /// Construct an empty native graph.
    pub fn new(name: impl Into<String>) -> Result<Self, GraphError> {
        let name = name.into();
        validate_name("graph name", &name)?;
        Ok(Self {
            name,
            legacy_kernel: None,
            nodes: BTreeMap::new(),
            edges: BTreeMap::new(),
        })
    }

    /// Bounded semantic graph name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonically ordered nodes.
    #[must_use]
    pub const fn nodes(&self) -> &BTreeMap<ModelId, FidelityNode> {
        &self.nodes
    }

    /// Canonically ordered edges.
    #[must_use]
    pub const fn edges(&self) -> &BTreeMap<EdgeId, FidelityEdge> {
        &self.edges
    }

    /// Add one unique model node.
    pub fn add_node(&mut self, node: FidelityNode) -> Result<(), GraphError> {
        if self.legacy_kernel.is_some() {
            return Err(GraphError::LegacyEmbeddingImmutable);
        }
        if self.nodes.len() >= MAX_GRAPH_NODES {
            return Err(GraphError::LimitExceeded {
                what: "graph nodes",
                limit: MAX_GRAPH_NODES,
            });
        }
        if self.nodes.contains_key(&node.id) {
            return Err(GraphError::DuplicateNode(node.id));
        }
        let requested = self
            .canonical_size()
            .checked_add(encoded_node_len(&node))
            .ok_or(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            })?;
        if requested > MAX_CANONICAL_BYTES {
            return Err(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            });
        }
        self.nodes.insert(node.id, node);
        Ok(())
    }

    /// Add one edge after resolving both model-card nodes.
    pub fn add_edge(&mut self, edge: FidelityEdge) -> Result<(), GraphError> {
        if self.legacy_kernel.is_some() {
            return Err(GraphError::LegacyEmbeddingImmutable);
        }
        if self.edges.len() >= MAX_GRAPH_EDGES {
            return Err(GraphError::LimitExceeded {
                what: "graph edges",
                limit: MAX_GRAPH_EDGES,
            });
        }
        for model in [edge.source, edge.target] {
            if !self.nodes.contains_key(&model) {
                return Err(GraphError::MissingNode(model));
            }
        }
        if self.edges.contains_key(&edge.id) {
            return Err(GraphError::DuplicateEdge(edge.id));
        }
        let requested = self
            .canonical_size()
            .checked_add(encoded_edge_len(&edge))
            .ok_or(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            })?;
        if requested > MAX_CANONICAL_BYTES {
            return Err(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            });
        }
        self.edges.insert(edge.id, edge);
        Ok(())
    }

    /// Lookup a node.
    #[must_use]
    pub fn node(&self, id: ModelId) -> Option<&FidelityNode> {
        self.nodes.get(&id)
    }

    /// Lookup an edge.
    #[must_use]
    pub fn edge(&self, id: EdgeId) -> Option<&FidelityEdge> {
        self.edges.get(&id)
    }

    /// Deterministic canonical transport bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.canonical_size());
        put_u16(&mut out, FIDELITY_GRAPH_SCHEMA_VERSION);
        put_string(&mut out, &self.name);
        match &self.legacy_kernel {
            None => put_u8(&mut out, 0),
            Some(kernel) => {
                put_u8(&mut out, 1);
                put_string(&mut out, kernel);
            }
        }
        put_u32(&mut out, self.nodes.len() as u32);
        for node in self.nodes.values() {
            out.extend_from_slice(node.id.as_bytes());
            out.extend_from_slice(node.card.as_bytes());
            put_string(&mut out, &node.label);
            match &node.legacy {
                None => put_u8(&mut out, 0),
                Some(legacy) => {
                    put_u8(&mut out, 1);
                    put_u32(&mut out, legacy.index);
                    put_u64(&mut out, legacy.relative_cost_bits);
                    put_string(&mut out, &legacy.note);
                }
            }
        }
        put_u32(&mut out, self.edges.len() as u32);
        for edge in self.edges.values() {
            out.extend_from_slice(edge.id.as_bytes());
            let body = edge_body_bytes(edge);
            put_bytes(&mut out, &body);
        }
        out
    }

    fn canonical_size(&self) -> usize {
        let legacy_kernel = self
            .legacy_kernel
            .as_ref()
            .map_or(0, |kernel| 4 + kernel.len());
        2 + 4
            + self.name.len()
            + 1
            + legacy_kernel
            + 4
            + self.nodes.values().map(encoded_node_len).sum::<usize>()
            + 4
            + self.edges.values().map(encoded_edge_len).sum::<usize>()
    }

    /// Decode bounded canonical transport and reject alternate encodings.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, GraphError> {
        if bytes.len() > MAX_CANONICAL_BYTES {
            return Err(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            });
        }
        let mut cursor = Cursor::new(bytes);
        let version = cursor.u16()?;
        if version != FIDELITY_GRAPH_SCHEMA_VERSION {
            return Err(GraphError::UnsupportedSchema(version));
        }
        let name = cursor.string("graph name")?;
        let mut graph = Self::new(name)?;
        let legacy_kernel = match cursor.u8()? {
            0 => None,
            1 => {
                let kernel = cursor.string("legacy kernel")?;
                validate_name("legacy kernel", &kernel)?;
                Some(kernel)
            }
            _ => return Err(GraphError::Decode("invalid legacy-kernel tag")),
        };
        let node_count = cursor.count("graph nodes", MAX_GRAPH_NODES)?;
        for _ in 0..node_count {
            let id = ModelId::new(ContentHash(cursor.fixed_32()?));
            let card = ModelCardRef::new(ContentHash(cursor.fixed_32()?));
            let label = cursor.string("model label")?;
            let legacy = match cursor.u8()? {
                0 => None,
                1 => {
                    let index = cursor.u32()?;
                    let relative_cost_bits = {
                        let bits = cursor.u64()?;
                        CostRelationRef::legacy(f64::from_bits(bits))?;
                        bits
                    };
                    let note = cursor.string_with_limit("legacy rung note", MAX_NOTE_BYTES)?;
                    validate_note(&note)?;
                    Some(LegacyRung {
                        index,
                        relative_cost_bits,
                        note,
                    })
                }
                _ => return Err(GraphError::Decode("invalid legacy-rung tag")),
            };
            let mut node = FidelityNode::new(id, card, label)?;
            node.legacy = legacy;
            graph.add_node(node)?;
        }
        let edge_count = cursor.count("graph edges", MAX_GRAPH_EDGES)?;
        for _ in 0..edge_count {
            let encoded_id = EdgeId::new(ContentHash(cursor.fixed_32()?));
            let body = cursor.bytes("edge body", MAX_CANONICAL_BYTES)?;
            let edge = decode_edge(body)?;
            if edge.id != encoded_id {
                return Err(GraphError::IdentityMismatch {
                    what: "edge identity",
                });
            }
            graph.add_edge(edge)?;
        }
        if !cursor.is_finished() {
            return Err(GraphError::Decode("trailing canonical bytes"));
        }
        graph.legacy_kernel = legacy_kernel;
        if graph.legacy_kernel.is_some() {
            graph.validate_legacy_embedding()?;
        }
        if graph.canonical_bytes() != bytes {
            return Err(GraphError::NonCanonical);
        }
        Ok(graph)
    }

    fn validate_legacy_embedding(&self) -> Result<(), GraphError> {
        let descriptor = self
            .embedded_ladder_descriptor()
            .ok_or(GraphError::IdentityMismatch {
                what: "legacy ladder descriptor",
            })?;
        if descriptor.kernel != self.name {
            return Err(GraphError::IdentityMismatch {
                what: "legacy kernel and graph name",
            });
        }
        if descriptor.rungs.is_empty() {
            return Err(GraphError::IdentityMismatch {
                what: "legacy ladder must contain a rung",
            });
        }
        if self.edges.len() != descriptor.rungs.len().saturating_sub(1) {
            return Err(GraphError::IdentityMismatch {
                what: "legacy ladder path cardinality",
            });
        }
        let mut models = Vec::with_capacity(descriptor.rungs.len());
        for rung in &descriptor.rungs {
            let preimage = legacy_rung_bytes(&descriptor.kernel, rung);
            let model = ModelId::new(hash_domain(LEGACY_MODEL_DOMAIN, &preimage));
            let node = self.nodes.get(&model).ok_or(GraphError::IdentityMismatch {
                what: "legacy model identity",
            })?;
            let expected_card = ModelCardRef::new(hash_domain(LEGACY_CARD_DOMAIN, &preimage));
            if node.card != expected_card {
                return Err(GraphError::IdentityMismatch {
                    what: "legacy model-card identity",
                });
            }
            models.push(model);
        }
        for index in 0..models.len().saturating_sub(1) {
            let mut transfer_preimage = Vec::new();
            put_string(&mut transfer_preimage, &descriptor.kernel);
            put_u32(&mut transfer_preimage, index as u32);
            let expected_transfer =
                TransferRef::new(hash_domain(LEGACY_TRANSFER_DOMAIN, &transfer_preimage));
            let edge = self
                .edges
                .values()
                .find(|edge| edge.source == models[index] && edge.target == models[index + 1])
                .ok_or(GraphError::IdentityMismatch {
                    what: "legacy ladder path edge",
                })?;
            if edge.cost != CostRelationRef::legacy(descriptor.rungs[index + 1].relative_cost)?
                || edge.discrepancy != DiscrepancyReference::UnknownLegacy
                || edge.transfer != expected_transfer
                || edge.validity != ValidityDomain::universal()
                || edge.informativeness != Informativeness::legacy_total_order()
            {
                return Err(GraphError::IdentityMismatch {
                    what: "legacy ladder edge semantics",
                });
            }
        }
        Ok(())
    }

    /// Domain-separated identity over exact canonical bytes.
    #[must_use]
    pub fn identity(&self) -> FidelityGraphId {
        FidelityGraphId::new(hash_domain(GRAPH_IDENTITY_DOMAIN, &self.canonical_bytes()))
    }

    /// Select a maximal evidenced model reachable inside the cost budget.
    ///
    /// Incomparable maxima are not silently called scientifically equivalent.
    /// The result names them and uses cost then model identity only as an
    /// operational replay tie-break.
    pub fn best_model_for(
        &self,
        start: ModelId,
        context: &QueryContext,
        resolver: &impl EdgeEvidenceResolver,
    ) -> Result<ModelRecommendation, QueryRefusal> {
        let search = self.search(start, context, resolver)?;
        let mut non_maximal = BTreeSet::new();
        for edge in &search.considered {
            if search.traversed.contains(&edge.edge) {
                non_maximal.insert(edge.source);
            }
        }
        let mut maxima: Vec<_> = search
            .reachable
            .iter()
            .filter(|(model, candidate)| {
                !non_maximal.contains(model) && candidate.cost_s.is_finite()
            })
            .map(|(model, candidate)| (*model, candidate.clone()))
            .collect();
        if maxima.is_empty() {
            return Err(QueryRefusal::NoApplicableEvidence {
                start,
                qoi: context.qoi.clone(),
            });
        }
        maxima.sort_by(candidate_order);
        let (selected, candidate) = maxima[0].clone();
        let incomparable_maxima = maxima.iter().skip(1).map(|(model, _)| *model).collect();
        let basis = if maxima.len() == 1 {
            SelectionBasis::UniqueGraphMaximum
        } else {
            SelectionBasis::IncomparableMaximaCostThenIdentity
        };
        Ok(ModelRecommendation {
            model: selected,
            predicted_cost_s: candidate.cost_s,
            path: candidate.path,
            explanation: QueryExplanation {
                graph: self.identity(),
                qoi: context.qoi.clone(),
                problem_size: context.problem_size,
                budget_s: context.budget_s,
                max_relative_discrepancy: context.max_relative_discrepancy,
                start,
                regime: context.regime.clone(),
                basis,
                considered: search.considered,
                unresolved: search.unresolved,
                incomparable_maxima,
                matched_specificity: candidate.specificity,
            },
        })
    }

    /// Choose the cheapest reachable source model carrying an adequate
    /// pairwise discrepancy assessment for this exact context.
    pub fn cheapest_adequate(
        &self,
        start: ModelId,
        context: &QueryContext,
        resolver: &impl EdgeEvidenceResolver,
    ) -> Result<ModelRecommendation, QueryRefusal> {
        let search = self.search(start, context, resolver)?;
        let unresolved_sources: BTreeSet<_> = search
            .unresolved
            .iter()
            .filter_map(|edge| self.edges.get(edge).map(|edge| edge.source))
            .collect();
        let mut by_source: BTreeMap<ModelId, Vec<&ConsideredEdge>> = BTreeMap::new();
        for evidence in &search.considered {
            if search.reachable.contains_key(&evidence.source) {
                by_source.entry(evidence.source).or_default().push(evidence);
            }
        }
        let mut adequate: BTreeMap<ModelId, Candidate> = BTreeMap::new();
        for (source, evidence_rows) in by_source {
            if unresolved_sources.contains(&source)
                || evidence_rows
                    .iter()
                    .any(|evidence| evidence.source_adequacy != Adequacy::Adequate)
            {
                continue;
            }
            let source_cost_s = evidence_rows
                .iter()
                .map(|evidence| evidence.source_cost_s)
                .max_by(f64::total_cmp)
                .unwrap_or(f64::INFINITY);
            if source_cost_s > context.budget_s {
                continue;
            }
            let reachable = &search.reachable[&source];
            let candidate = Candidate {
                cost_s: source_cost_s,
                path: reachable.path.clone(),
                specificity: reachable
                    .specificity
                    .saturating_add(evidence_rows.iter().fold(0, |total, evidence| {
                        total
                            .saturating_add(evidence.validity_specificity)
                            .saturating_add(evidence.informativeness_specificity)
                    })),
            };
            match adequate.get(&source) {
                Some(current)
                    if candidate_order(
                        &(source, current.clone()),
                        &(source, candidate.clone()),
                    ) != Ordering::Greater => {}
                _ => {
                    adequate.insert(source, candidate);
                }
            }
        }
        let mut candidates: Vec<_> = adequate.into_iter().collect();
        if candidates.is_empty() {
            return Err(QueryRefusal::NoAdequateModel {
                start,
                qoi: context.qoi.clone(),
            });
        }
        candidates.sort_by(candidate_order);
        let (selected, candidate) = candidates[0].clone();
        Ok(ModelRecommendation {
            model: selected,
            predicted_cost_s: candidate.cost_s,
            path: candidate.path,
            explanation: QueryExplanation {
                graph: self.identity(),
                qoi: context.qoi.clone(),
                problem_size: context.problem_size,
                budget_s: context.budget_s,
                max_relative_discrepancy: context.max_relative_discrepancy,
                start,
                regime: context.regime.clone(),
                basis: SelectionBasis::CheapestAdequate,
                considered: search.considered,
                unresolved: search.unresolved,
                incomparable_maxima: Vec::new(),
                matched_specificity: candidate.specificity,
            },
        })
    }

    fn search(
        &self,
        start: ModelId,
        context: &QueryContext,
        resolver: &impl EdgeEvidenceResolver,
    ) -> Result<SearchResult, QueryRefusal> {
        if !self.nodes.contains_key(&start) {
            return Err(QueryRefusal::UnknownStart(start));
        }
        let mut considered = Vec::new();
        let mut unresolved = Vec::new();
        for edge in self.edges.values() {
            let Some(validity_specificity) = edge.validity.0.matched_specificity(context) else {
                continue;
            };
            let Some(informativeness_specificity) =
                edge.informativeness.0.matched_specificity(context)
            else {
                continue;
            };
            let Some(resolved) = resolver.resolve(edge, context) else {
                unresolved.push(edge.id);
                continue;
            };
            if !resolved.matches(edge) {
                unresolved.push(edge.id);
                continue;
            }
            considered.push(ConsideredEdge {
                edge: edge.id,
                source: edge.source,
                target: edge.target,
                source_cost_s: resolved.source_cost_s,
                target_cost_s: resolved.target_cost_s,
                source_adequacy: resolved.adequacy(context),
                assessed_relative_discrepancy: resolved.assessed_relative_discrepancy,
                evidence: resolved.evidence,
                validity_specificity,
                informativeness_specificity,
            });
        }
        considered.sort_by_key(|edge| edge.edge);
        unresolved.sort_unstable();
        unresolved.dedup();

        let mut outgoing: BTreeMap<ModelId, Vec<&ConsideredEdge>> = BTreeMap::new();
        for edge in &considered {
            outgoing.entry(edge.source).or_default().push(edge);
        }
        let mut reachable = BTreeMap::new();
        reachable.insert(
            start,
            Candidate {
                cost_s: f64::INFINITY,
                path: Vec::new(),
                specificity: 0,
            },
        );
        let mut traversed = BTreeSet::new();
        let mut queue = VecDeque::from([start]);
        while let Some(source) = queue.pop_front() {
            let source_path = reachable
                .get(&source)
                .map(|candidate| candidate.path.clone())
                .unwrap_or_default();
            let Some(edges) = outgoing.get(&source) else {
                continue;
            };
            for evidence in edges {
                if evidence.source_cost_s > context.budget_s {
                    continue;
                }
                let current = reachable.get(&source).cloned();
                let candidate = Candidate {
                    cost_s: evidence.source_cost_s,
                    path: source_path.clone(),
                    specificity: current.map_or(0, |value| value.specificity),
                };
                update_candidate(&mut reachable, source, candidate);
                if evidence.target_cost_s > context.budget_s {
                    continue;
                }
                traversed.insert(evidence.edge);
                let mut path = source_path.clone();
                path.push(evidence.edge);
                let candidate = Candidate {
                    cost_s: evidence.target_cost_s,
                    path,
                    specificity: reachable
                        .get(&source)
                        .map_or(0, |value| value.specificity)
                        .saturating_add(evidence.validity_specificity)
                        .saturating_add(evidence.informativeness_specificity),
                };
                if update_candidate(&mut reachable, evidence.target, candidate) {
                    queue.push_back(evidence.target);
                }
            }
        }
        Ok(SearchResult {
            reachable,
            considered,
            unresolved,
            traversed,
        })
    }

    /// Recover the exact declaration of an embedded legacy ladder.
    #[must_use]
    pub fn embedded_ladder_descriptor(&self) -> Option<LadderDescriptor> {
        let kernel = self.legacy_kernel.clone()?;
        let mut rungs: Vec<_> = self
            .nodes
            .values()
            .map(|node| {
                let legacy = node.legacy.as_ref()?;
                Some(Rung {
                    index: legacy.index,
                    name: node.label.clone(),
                    relative_cost: f64::from_bits(legacy.relative_cost_bits),
                    note: legacy.note.clone(),
                })
            })
            .collect::<Option<Vec<_>>>()?;
        rungs.sort_by_key(|rung| rung.index);
        if rungs
            .iter()
            .enumerate()
            .any(|(index, rung)| rung.index as usize != index)
        {
            return None;
        }
        Some(LadderDescriptor { kernel, rungs })
    }
}

fn update_candidate(
    reachable: &mut BTreeMap<ModelId, Candidate>,
    model: ModelId,
    candidate: Candidate,
) -> bool {
    let replace = reachable.get(&model).is_none_or(|current| {
        candidate
            .cost_s
            .total_cmp(&current.cost_s)
            .then(candidate.path.cmp(&current.path))
            .then(current.specificity.cmp(&candidate.specificity))
            == Ordering::Less
    });
    if replace {
        reachable.insert(model, candidate);
    }
    replace
}

fn candidate_order(left: &(ModelId, Candidate), right: &(ModelId, Candidate)) -> Ordering {
    left.1
        .cost_s
        .total_cmp(&right.1.cost_s)
        .then(left.0.cmp(&right.0))
        .then(left.1.path.cmp(&right.1.path))
}

#[derive(Debug, Clone)]
struct Candidate {
    cost_s: f64,
    path: Vec<EdgeId>,
    specificity: u32,
}

struct SearchResult {
    reachable: BTreeMap<ModelId, Candidate>,
    considered: Vec<ConsideredEdge>,
    unresolved: Vec<EdgeId>,
    traversed: BTreeSet<EdgeId>,
}

/// Context-specific adequacy of an edge's source relative to its target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Adequacy {
    /// The exact discrepancy assessment meets the caller's tolerance.
    Adequate,
    /// The exact discrepancy assessment exceeds the caller's tolerance.
    Inadequate,
    /// The resolver has no conclusion.
    Unknown,
}

/// Query-time evaluation of the exact cost/discrepancy artifacts on one edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResolvedEdgeEvidence {
    cost_model: Option<CostModelRef>,
    discrepancy_model: Option<DiscrepancyModelRef>,
    source_cost_s: f64,
    target_cost_s: f64,
    assessed_relative_discrepancy: Option<f64>,
    evidence: QueryEvidenceRef,
}

impl ResolvedEdgeEvidence {
    /// Construct one finite, nonnegative resolution.
    pub fn new(
        cost_model: Option<CostModelRef>,
        discrepancy_model: Option<DiscrepancyModelRef>,
        source_cost_s: f64,
        target_cost_s: f64,
        assessed_relative_discrepancy: Option<f64>,
        evidence: QueryEvidenceRef,
    ) -> Result<Self, GraphError> {
        for (field, value) in [
            ("resolved source cost seconds", source_cost_s),
            ("resolved target cost seconds", target_cost_s),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(GraphError::InvalidNumber { field, value });
            }
        }
        if let Some(value) = assessed_relative_discrepancy
            && (!value.is_finite() || value < 0.0)
        {
            return Err(GraphError::InvalidNumber {
                field: "assessed relative discrepancy",
                value,
            });
        }
        Ok(Self {
            cost_model,
            discrepancy_model,
            source_cost_s: canonical_zero(source_cost_s),
            target_cost_s: canonical_zero(target_cost_s),
            assessed_relative_discrepancy: assessed_relative_discrepancy.map(canonical_zero),
            evidence,
        })
    }

    fn adequacy(self, context: &QueryContext) -> Adequacy {
        match self.assessed_relative_discrepancy {
            Some(value) if value <= context.max_relative_discrepancy => Adequacy::Adequate,
            Some(_) => Adequacy::Inadequate,
            None => Adequacy::Unknown,
        }
    }

    fn matches(self, edge: &FidelityEdge) -> bool {
        let cost_matches = match edge.cost {
            CostRelationRef::Model(expected) => self.cost_model == Some(expected),
            CostRelationRef::LegacyRelativeCost(_) => self.cost_model.is_none(),
        };
        let discrepancy_matches = match edge.discrepancy {
            DiscrepancyReference::Model(expected) => self.discrepancy_model == Some(expected),
            DiscrepancyReference::UnknownLegacy => self.discrepancy_model.is_none(),
        };
        cost_matches && discrepancy_matches
    }
}

/// Adapter boundary implemented by `fs-plan`/`fs-evidence` consumers.
///
/// `None` is an honest unresolved edge. Returned artifact references must
/// exactly match the edge or the graph treats the resolution as unresolved.
pub trait EdgeEvidenceResolver {
    /// Resolve cost and discrepancy evidence for one exact edge/context.
    fn resolve(&self, edge: &FidelityEdge, context: &QueryContext) -> Option<ResolvedEdgeEvidence>;
}

/// One exact edge included in a query explanation.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsideredEdge {
    /// Edge identity.
    pub edge: EdgeId,
    /// Source model.
    pub source: ModelId,
    /// Target model.
    pub target: ModelId,
    /// Predicted source cost in seconds.
    pub source_cost_s: f64,
    /// Predicted target cost in seconds.
    pub target_cost_s: f64,
    /// Pairwise source adequacy.
    pub source_adequacy: Adequacy,
    /// Resolver-supplied relative discrepancy used to derive adequacy.
    pub assessed_relative_discrepancy: Option<f64>,
    /// Exact query-evaluation receipt.
    pub evidence: QueryEvidenceRef,
    /// Specificity of the matched validity clause.
    pub validity_specificity: u32,
    /// Specificity of the matched informativeness clause.
    pub informativeness_specificity: u32,
}

/// Why a deterministic recommendation was selected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionBasis {
    /// The selected model is the only reachable graph maximum in budget.
    UniqueGraphMaximum,
    /// Maxima were incomparable; cost then model identity was used only for
    /// deterministic operation, not epistemic ranking.
    IncomparableMaximaCostThenIdentity,
    /// Cheapest reachable model with exact adequate discrepancy evidence.
    CheapestAdequate,
}

/// Recorded reasoning behind one recommendation.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryExplanation {
    /// Exact graph identity.
    pub graph: FidelityGraphId,
    /// Requested QoI.
    pub qoi: QoiId,
    /// Problem size supplied to cost resolvers.
    pub problem_size: u64,
    /// Cost budget in seconds.
    pub budget_s: f64,
    /// Maximum admitted relative discrepancy.
    pub max_relative_discrepancy: f64,
    /// Requested starting model.
    pub start: ModelId,
    /// Exact finite regime coordinates.
    pub regime: BTreeMap<RegimeAxis, f64>,
    /// Selection rule.
    pub basis: SelectionBasis,
    /// Context-matching edges with exact resolution receipts.
    pub considered: Vec<ConsideredEdge>,
    /// Context-matching edges without an exact usable resolution.
    pub unresolved: Vec<EdgeId>,
    /// Other incomparable graph maxima, in deterministic order.
    pub incomparable_maxima: Vec<ModelId>,
    /// Sum of matched clause specificity along the selected path.
    pub matched_specificity: u32,
}

/// Query result plus its replay explanation.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelRecommendation {
    /// Selected model.
    pub model: ModelId,
    /// Resolved predicted total cost in seconds.
    pub predicted_cost_s: f64,
    /// Exact selected path from the requested start model.
    pub path: Vec<EdgeId>,
    /// Complete bounded reasoning record.
    pub explanation: QueryExplanation,
}

/// Structured query refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryRefusal {
    /// The requested start model is absent.
    UnknownStart(ModelId),
    /// No context-matching, resolved maximum exists.
    NoApplicableEvidence {
        /// Requested start.
        start: ModelId,
        /// Requested QoI.
        qoi: QoiId,
    },
    /// No reachable source has an adequate pairwise discrepancy assessment.
    NoAdequateModel {
        /// Requested start.
        start: ModelId,
        /// Requested QoI.
        qoi: QoiId,
    },
}

impl fmt::Display for QueryRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownStart(model) => {
                write!(f, "fidelity query start model {model} is not in the graph")
            }
            Self::NoApplicableEvidence { start, qoi } => write!(
                f,
                "no resolved, context-matching fidelity evidence is reachable from {start} for QoI {qoi}"
            ),
            Self::NoAdequateModel { start, qoi } => write!(
                f,
                "no reachable model from {start} has adequate pairwise discrepancy evidence for QoI {qoi}"
            ),
        }
    }
}

impl Error for QueryRefusal {}

/// Runtime transfer operators keyed by exact graph-edge identity.
#[derive(Default)]
pub struct GraphTransfers {
    transfers: BTreeMap<EdgeId, Box<dyn Transfer>>,
}

impl fmt::Debug for GraphTransfers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GraphTransfers")
            .field("edge_count", &self.transfers.len())
            .finish()
    }
}

impl GraphTransfers {
    /// Construct an empty transfer registry.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            transfers: BTreeMap::new(),
        }
    }

    /// Register an exact edge transfer. Duplicate registrations refuse.
    pub fn register(
        &mut self,
        edge: EdgeId,
        transfer: Box<dyn Transfer>,
    ) -> Result<(), GraphError> {
        if self.transfers.contains_key(&edge) {
            return Err(GraphError::DuplicateTransfer(edge));
        }
        self.transfers.insert(edge, transfer);
        Ok(())
    }

    /// Apply source-to-target transfer on an existing edge.
    pub fn prolongate(
        &self,
        graph: &FidelityGraph,
        edge: EdgeId,
        source: &[f64],
    ) -> Result<Vec<f64>, GraphError> {
        if graph.edge(edge).is_none() {
            return Err(GraphError::MissingEdge(edge));
        }
        let transfer = self
            .transfers
            .get(&edge)
            .ok_or(GraphError::MissingTransfer(edge))?;
        Ok(transfer.prolongate(source))
    }

    /// Apply target-to-source transfer on an existing edge.
    pub fn restrict(
        &self,
        graph: &FidelityGraph,
        edge: EdgeId,
        target: &[f64],
    ) -> Result<Vec<f64>, GraphError> {
        if graph.edge(edge).is_none() {
            return Err(GraphError::MissingEdge(edge));
        }
        let transfer = self
            .transfers
            .get(&edge)
            .ok_or(GraphError::MissingTransfer(edge))?;
        Ok(transfer.restrict(target))
    }
}

/// Transfer-owning lossless embedding of one legacy ladder.
#[derive(Debug)]
pub struct EmbeddedLadderGraph {
    graph: FidelityGraph,
    transfers: GraphTransfers,
}

impl EmbeddedLadderGraph {
    /// Embedded graph.
    #[must_use]
    pub const fn graph(&self) -> &FidelityGraph {
        &self.graph
    }

    /// Edge-keyed transfer registry.
    #[must_use]
    pub const fn transfers(&self) -> &GraphTransfers {
        &self.transfers
    }

    /// Split graph declarations from runtime operators.
    #[must_use]
    pub fn into_parts(self) -> (FidelityGraph, GraphTransfers) {
        (self.graph, self.transfers)
    }
}

/// Exact serializable declaration of a legacy ladder.
#[derive(Debug, Clone, PartialEq)]
pub struct LadderDescriptor {
    kernel: String,
    rungs: Vec<Rung>,
}

impl LadderDescriptor {
    /// Kernel identity.
    #[must_use]
    pub fn kernel(&self) -> &str {
        &self.kernel
    }

    /// Exact ordered rung declarations.
    #[must_use]
    pub fn rungs(&self) -> &[Rung] {
        &self.rungs
    }
}

impl Ladder {
    /// Snapshot the exact serializable legacy declaration.
    #[must_use]
    pub fn descriptor(&self) -> LadderDescriptor {
        LadderDescriptor {
            kernel: self.kernel.clone(),
            rungs: self.rungs.clone(),
        }
    }

    /// Consume a legacy ladder into a path graph while preserving every rung
    /// field and moving every transfer to its exact edge.
    pub fn into_fidelity_graph(self) -> Result<EmbeddedLadderGraph, GraphError> {
        let Ladder {
            kernel,
            rungs,
            transfers,
        } = self;
        let mut graph = FidelityGraph::new(kernel.clone())?;
        let mut model_ids = Vec::with_capacity(rungs.len());
        for rung in &rungs {
            CostRelationRef::legacy(rung.relative_cost)?;
            validate_note(&rung.note)?;
            let preimage = legacy_rung_bytes(&kernel, rung);
            let id = ModelId::new(hash_domain(LEGACY_MODEL_DOMAIN, &preimage));
            let card = ModelCardRef::new(hash_domain(LEGACY_CARD_DOMAIN, &preimage));
            let mut node = FidelityNode::new(id, card, rung.name.clone())?;
            node.legacy = Some(LegacyRung {
                index: rung.index,
                relative_cost_bits: rung.relative_cost.to_bits(),
                note: rung.note.clone(),
            });
            graph.add_node(node)?;
            model_ids.push(id);
        }
        let mut graph_transfers = GraphTransfers::new();
        for (index, transfer) in transfers.into_iter().enumerate() {
            let source = model_ids[index];
            let target = model_ids[index + 1];
            let mut transfer_preimage = Vec::new();
            put_string(&mut transfer_preimage, &kernel);
            put_u32(&mut transfer_preimage, index as u32);
            let transfer_ref =
                TransferRef::new(hash_domain(LEGACY_TRANSFER_DOMAIN, &transfer_preimage));
            let edge = FidelityEdge::build(
                source,
                target,
                CostRelationRef::legacy(rungs[index + 1].relative_cost)?,
                DiscrepancyReference::UnknownLegacy,
                transfer_ref,
                ValidityDomain::universal(),
                Informativeness::legacy_total_order(),
            )?;
            let edge_id = edge.id;
            graph.add_edge(edge)?;
            graph_transfers.register(edge_id, transfer)?;
        }
        graph.legacy_kernel = Some(kernel);
        if graph.canonical_size() > MAX_CANONICAL_BYTES {
            return Err(GraphError::LimitExceeded {
                what: "canonical graph bytes",
                limit: MAX_CANONICAL_BYTES,
            });
        }
        graph.validate_legacy_embedding()?;
        Ok(EmbeddedLadderGraph {
            graph,
            transfers: graph_transfers,
        })
    }
}

/// Graph construction/transport refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum GraphError {
    /// A semantic name is empty, oversized, or outside the ASCII grammar.
    InvalidName {
        /// Field kind.
        field: &'static str,
        /// Refused value.
        value: String,
    },
    /// A finite closed interval was not supplied.
    InvalidInterval {
        /// Lower endpoint.
        lower: f64,
        /// Upper endpoint.
        upper: f64,
    },
    /// A finite number satisfying the field's sign constraint was not supplied.
    InvalidNumber {
        /// Field kind.
        field: &'static str,
        /// Refused value.
        value: f64,
    },
    /// One regime axis appeared twice.
    DuplicateRegimeAxis(RegimeAxis),
    /// A bounded collection exceeded its limit.
    LimitExceeded {
        /// Collection kind.
        what: &'static str,
        /// Maximum admitted count/bytes.
        limit: usize,
    },
    /// A model identity was inserted twice.
    DuplicateNode(ModelId),
    /// An edge identity was inserted twice.
    DuplicateEdge(EdgeId),
    /// An edge referenced an absent model.
    MissingNode(ModelId),
    /// An edge targeted its own source.
    SelfLoop(ModelId),
    /// A runtime transfer was registered twice.
    DuplicateTransfer(EdgeId),
    /// A lossless legacy embedding cannot be mutated in place.
    LegacyEmbeddingImmutable,
    /// A requested graph edge is absent.
    MissingEdge(EdgeId),
    /// A requested edge has no runtime transfer.
    MissingTransfer(EdgeId),
    /// Canonical transport used an unsupported version.
    UnsupportedSchema(u16),
    /// Canonical bytes were malformed or truncated.
    Decode(&'static str),
    /// Parsed semantics re-encoded differently.
    NonCanonical,
    /// Encoded identity disagreed with the semantic body.
    IdentityMismatch {
        /// Identity kind.
        what: &'static str,
    },
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName { field, value } => write!(
                f,
                "invalid {field} `{value}`; use 1..={MAX_NAME_BYTES} visible ASCII identity bytes"
            ),
            Self::InvalidInterval { lower, upper } => write!(
                f,
                "invalid closed regime interval [{lower}, {upper}]; endpoints must be finite and ordered"
            ),
            Self::InvalidNumber { field, value } => {
                write!(f, "invalid {field} value {value}")
            }
            Self::DuplicateRegimeAxis(axis) => {
                write!(f, "regime axis `{axis}` was declared twice")
            }
            Self::LimitExceeded { what, limit } => {
                write!(f, "{what} exceed bounded limit {limit}")
            }
            Self::DuplicateNode(model) => write!(f, "model node {model} is duplicated"),
            Self::DuplicateEdge(edge) => write!(f, "fidelity edge {edge} is duplicated"),
            Self::MissingNode(model) => write!(f, "fidelity edge references missing model {model}"),
            Self::SelfLoop(model) => write!(f, "fidelity self-loop on model {model} is refused"),
            Self::DuplicateTransfer(edge) => {
                write!(f, "runtime transfer for edge {edge} is duplicated")
            }
            Self::LegacyEmbeddingImmutable => write!(
                f,
                "lossless legacy ladder embeddings are immutable; construct a native graph to add nodes or edges"
            ),
            Self::MissingEdge(edge) => write!(f, "fidelity edge {edge} is absent"),
            Self::MissingTransfer(edge) => {
                write!(f, "runtime transfer for fidelity edge {edge} is absent")
            }
            Self::UnsupportedSchema(version) => {
                write!(f, "unsupported fidelity-graph schema version {version}")
            }
            Self::Decode(reason) => write!(f, "fidelity-graph decode refused: {reason}"),
            Self::NonCanonical => write!(f, "fidelity-graph transport is not canonical"),
            Self::IdentityMismatch { what } => {
                write!(f, "fidelity-graph {what} does not match its semantic body")
            }
        }
    }
}

impl Error for GraphError {}

fn legacy_rung_bytes(kernel: &str, rung: &Rung) -> Vec<u8> {
    let mut out = Vec::new();
    put_string(&mut out, kernel);
    put_u32(&mut out, rung.index);
    put_string(&mut out, &rung.name);
    put_u64(&mut out, rung.relative_cost.to_bits());
    put_string(&mut out, &rung.note);
    out
}

fn encoded_node_len(node: &FidelityNode) -> usize {
    32 + 32
        + 4
        + node.label.len()
        + 1
        + node
            .legacy
            .as_ref()
            .map_or(0, |legacy| 4 + 8 + 4 + legacy.note.len())
}

fn encoded_edge_len(edge: &FidelityEdge) -> usize {
    32 + 4 + edge_body_bytes(edge).len()
}

fn clause_bytes(clause: &ContextClause) -> Vec<u8> {
    let mut out = Vec::new();
    encode_clause(&mut out, clause);
    out
}

fn encode_clause(out: &mut Vec<u8>, clause: &ContextClause) {
    match &clause.qoi {
        QoiSelector::Any => put_u8(out, 0),
        QoiSelector::Exact(qoi) => {
            put_u8(out, 1);
            put_string(out, qoi.as_str());
        }
    }
    put_u32(out, clause.axes.len() as u32);
    for (axis, interval) in &clause.axes {
        put_string(out, axis.as_str());
        put_u64(out, interval.lower.to_bits());
        put_u64(out, interval.upper.to_bits());
    }
}

fn encode_predicates(out: &mut Vec<u8>, predicates: &ContextPredicateSet) {
    put_u32(out, predicates.clauses.len() as u32);
    for clause in &predicates.clauses {
        encode_clause(out, clause);
    }
}

fn edge_body_bytes(edge: &FidelityEdge) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(edge.source.as_bytes());
    out.extend_from_slice(edge.target.as_bytes());
    match edge.cost {
        CostRelationRef::Model(reference) => {
            put_u8(&mut out, 0);
            out.extend_from_slice(reference.as_bytes());
        }
        CostRelationRef::LegacyRelativeCost(bits) => {
            put_u8(&mut out, 1);
            put_u64(&mut out, bits);
        }
    }
    match edge.discrepancy {
        DiscrepancyReference::Model(reference) => {
            put_u8(&mut out, 0);
            out.extend_from_slice(reference.as_bytes());
        }
        DiscrepancyReference::UnknownLegacy => put_u8(&mut out, 1),
    }
    out.extend_from_slice(edge.transfer.as_bytes());
    encode_predicates(&mut out, &edge.validity.0);
    encode_predicates(&mut out, &edge.informativeness.0);
    out
}

fn decode_edge(bytes: &[u8]) -> Result<FidelityEdge, GraphError> {
    let mut cursor = Cursor::new(bytes);
    let source = ModelId::new(ContentHash(cursor.fixed_32()?));
    let target = ModelId::new(ContentHash(cursor.fixed_32()?));
    let cost = match cursor.u8()? {
        0 => CostRelationRef::Model(CostModelRef::new(ContentHash(cursor.fixed_32()?))),
        1 => {
            let bits = cursor.u64()?;
            CostRelationRef::legacy(f64::from_bits(bits))?
        }
        _ => return Err(GraphError::Decode("invalid cost-relation tag")),
    };
    let discrepancy = match cursor.u8()? {
        0 => DiscrepancyReference::Model(DiscrepancyModelRef::new(ContentHash(cursor.fixed_32()?))),
        1 => DiscrepancyReference::UnknownLegacy,
        _ => return Err(GraphError::Decode("invalid discrepancy-reference tag")),
    };
    let transfer = TransferRef::new(ContentHash(cursor.fixed_32()?));
    let validity = ValidityDomain::new(decode_predicates(&mut cursor)?);
    let informativeness = Informativeness::new(decode_predicates(&mut cursor)?);
    if !cursor.is_finished() {
        return Err(GraphError::Decode("trailing edge-body bytes"));
    }
    FidelityEdge::build(
        source,
        target,
        cost,
        discrepancy,
        transfer,
        validity,
        informativeness,
    )
}

fn decode_predicates(cursor: &mut Cursor<'_>) -> Result<ContextPredicateSet, GraphError> {
    let count = cursor.count("context clauses", MAX_CONTEXT_CLAUSES)?;
    let mut clauses = Vec::with_capacity(count);
    for _ in 0..count {
        let qoi = match cursor.u8()? {
            0 => QoiSelector::Any,
            1 => QoiSelector::Exact(QoiId::new(cursor.string("qoi")?)?),
            _ => return Err(GraphError::Decode("invalid QoI selector tag")),
        };
        let axis_count = cursor.count("context regime axes", MAX_REGIME_AXES)?;
        let mut axes = Vec::with_capacity(axis_count);
        for _ in 0..axis_count {
            let axis = RegimeAxis::new(cursor.string("regime axis")?)?;
            let lower = f64::from_bits(cursor.u64()?);
            let upper = f64::from_bits(cursor.u64()?);
            axes.push((axis, ClosedInterval::new(lower, upper)?));
        }
        clauses.push(ContextClause::new(qoi, axes)?);
    }
    ContextPredicateSet::new(clauses)
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

fn put_string(out: &mut Vec<u8>, value: &str) {
    put_bytes(out, value.as_bytes());
}

fn put_bytes(out: &mut Vec<u8>, value: &[u8]) {
    put_u32(out, value.len() as u32);
    out.extend_from_slice(value);
}

struct Cursor<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Cursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], GraphError> {
        let end = self
            .position
            .checked_add(count)
            .ok_or(GraphError::Decode("cursor overflow"))?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(GraphError::Decode("truncated canonical bytes"))?;
        self.position = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, GraphError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, GraphError> {
        self.take(2)?
            .try_into()
            .map(u16::from_le_bytes)
            .map_err(|_| GraphError::Decode("invalid u16"))
    }

    fn u32(&mut self) -> Result<u32, GraphError> {
        self.take(4)?
            .try_into()
            .map(u32::from_le_bytes)
            .map_err(|_| GraphError::Decode("invalid u32"))
    }

    fn u64(&mut self) -> Result<u64, GraphError> {
        self.take(8)?
            .try_into()
            .map(u64::from_le_bytes)
            .map_err(|_| GraphError::Decode("invalid u64"))
    }

    fn fixed_32(&mut self) -> Result<[u8; 32], GraphError> {
        self.take(32)?
            .try_into()
            .map_err(|_| GraphError::Decode("invalid 32-byte identity"))
    }

    fn count(&mut self, what: &'static str, maximum: usize) -> Result<usize, GraphError> {
        let count = self.u32()? as usize;
        if count > maximum {
            return Err(GraphError::LimitExceeded {
                what,
                limit: maximum,
            });
        }
        Ok(count)
    }

    fn bytes(&mut self, what: &'static str, maximum: usize) -> Result<&'a [u8], GraphError> {
        let count = self.count(what, maximum)?;
        self.take(count)
    }

    fn string(&mut self, what: &'static str) -> Result<String, GraphError> {
        self.string_with_limit(what, MAX_NAME_BYTES)
    }

    fn string_with_limit(
        &mut self,
        what: &'static str,
        maximum: usize,
    ) -> Result<String, GraphError> {
        let bytes = self.bytes(what, maximum)?;
        let value = std::str::from_utf8(bytes)
            .map_err(|_| GraphError::Decode("semantic string is not UTF-8"))?;
        Ok(value.to_string())
    }

    fn is_finished(&self) -> bool {
        self.position == self.bytes.len()
    }
}
