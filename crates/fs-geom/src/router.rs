//! The Rep Router (plan §7.3, Bet 1): conversions form a directed
//! multigraph — nodes are chart KINDS, edges are converters annotated with
//! (cost model, error model, certificate availability) — and a request
//! ("reach kind K with composed error ≤ ε under cost budget B") is solved
//! as a PARETO shortest-path problem over that graph.
//!
//! - **Planner**: label-correcting multi-objective search over
//!   (cost, composed absolute error, uncertified-edge count), exact on
//!   nonnegative weights (non-simple paths are dominated by their simple
//!   reductions; the planner enumerates simple paths only). Winner rule,
//!   in order: all-certified preferred, then minimal cost, then minimal
//!   error, then lexicographic path — deterministic and explainable
//!   (AGENTS.md invariant; P2).
//! - **Cost oracle**: cost models are machine-specific ON PURPOSE
//!   (§14.1 — kernel-level winners flip between the reference machines).
//!   The router reads measured costs through [`CostOracle`] — an
//!   abstraction because fs-geom (L2) must not depend on fs-ledger (L6);
//!   HELM wires the ledger `tune` table behind it, and
//!   [`MemoryCostOracle`] serves in-process learning and tests.
//! - **Execution**: [`Router::execute`] runs a plan's chain through
//!   [`EdgeRunner`]s, COMPOSES the per-edge Evidence receipts
//!   (fs-evidence `Op::Add` over achieved-error enclosures — total error
//!   of a chain is bounded by the sum), and records actual cost/error
//!   back into the oracle so later plans improve.
//! - **Refusals teach** (P10): no admissible path → a structured refusal
//!   naming the BINDING constraint and the cheapest relaxations.
//!
//! Error composition semantics (declared per edge, composed per path, all
//! against the request's reference `scale`): `Exact` contributes nothing;
//! `AdditiveAbs(a)` adds `a`; `MultiplicativeRel(r)` amplifies incoming
//! error by `(1+r)` and adds `r·scale` — conservative first-order model,
//! CONTRACT.md documents the boundary.

use core::fmt;
use std::collections::BTreeMap;

use fs_evidence::{Certified, Evidence, NumericalCertificate, Op, ProvenanceHash};
use fs_exec::Cx;

/// Per-edge error model with its path-composition rule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ErrorModel {
    /// No error introduced (e.g. exact-predicate mesh→SDF sign paths).
    Exact,
    /// Adds an absolute (Hausdorff-class) error.
    AdditiveAbs(f64),
    /// Relative error against the request scale; amplifies incoming error.
    MultiplicativeRel(f64),
}

impl ErrorModel {
    /// Composed absolute error after crossing this edge.
    fn compose(self, incoming_abs: f64, scale: f64) -> f64 {
        match self {
            ErrorModel::Exact => incoming_abs,
            ErrorModel::AdditiveAbs(a) => incoming_abs + a,
            ErrorModel::MultiplicativeRel(r) => incoming_abs * (1.0 + r) + r * scale,
        }
    }
}

/// A registered converter edge.
#[derive(Debug, Clone, PartialEq)]
pub struct ConverterSpec {
    /// Source chart kind.
    pub from: String,
    /// Destination chart kind.
    pub to: String,
    /// Unique edge name (ledger/tune key; e.g. `"frep->sdf/sampled-v1"`).
    pub name: String,
    /// A-priori cost estimate in seconds (used until measurements exist).
    pub base_cost_s: f64,
    /// Declared error model.
    pub error: ErrorModel,
    /// Whether the edge's error claim is certificate-backed. Certified
    /// edges are PREFERRED by the winner rule; uncertified (estimated)
    /// edges may have their additive error magnitude replaced by learned
    /// measurements — certificates are never "learned" away.
    pub certified: bool,
}

/// Measured-cost source (HELM wires the ledger `tune` table behind this;
/// L2 stays ledger-free).
pub trait CostOracle {
    /// Mean measured wall cost for an edge, if any measurement exists.
    fn measured_cost_s(&self, edge: &str) -> Option<f64>;
    /// Mean measured absolute error for an edge, if any.
    fn measured_error_abs(&self, edge: &str) -> Option<f64>;
    /// Atomically record one validated execution batch.
    ///
    /// # Errors
    /// Invalid, overflowing, or capacity-exceeding evidence is refused before
    /// it can influence later routes. An implementation must apply all rows or
    /// none: a chain cannot leave half of a learning update behind.
    fn record_batch(
        &mut self,
        observations: &[ValidatedEdgeObservation],
    ) -> Result<(), CostOracleError>;
}

/// Strength of the structurally validated edge evidence that authorized one
/// learned observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeEvidenceClass {
    /// Rigorous `Certified<f64>` evidence.
    Certified,
    /// Bounded non-rigorous estimate constructed through the sealed edge API.
    Estimated,
}

/// One execution observation admitted by [`Router::execute`]. Fields are
/// private so raw scalars cannot become future routing authority.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedEdgeObservation {
    edge: String,
    cost_s: f64,
    conservative_error_abs: f64,
    evidence_class: EdgeEvidenceClass,
    provenance: ProvenanceHash,
}

impl ValidatedEdgeObservation {
    /// Exact converter identity.
    #[must_use]
    pub fn edge(&self) -> &str {
        &self.edge
    }

    /// Measured wall cost in seconds.
    #[must_use]
    pub fn cost_s(&self) -> f64 {
        self.cost_s
    }

    /// Conservative learned absolute error (the receipt's upper bound, not
    /// its point QoI).
    #[must_use]
    pub fn conservative_error_abs(&self) -> f64 {
        self.conservative_error_abs
    }

    /// Evidence strength that admitted the observation.
    #[must_use]
    pub fn evidence_class(&self) -> EdgeEvidenceClass {
        self.evidence_class
    }

    /// Producing operation identity retained with the observation.
    #[must_use]
    pub fn provenance(&self) -> ProvenanceHash {
        self.provenance
    }
}

/// Why an oracle refused an executed edge's measurement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CostOracleError {
    /// The edge identity is empty, malformed, or unknown to the oracle.
    InvalidEdge {
        /// Human-readable, agent-actionable diagnosis.
        problem: String,
    },
    /// A measured scalar is outside its declared finite domain.
    InvalidMeasurement {
        /// Rejected field (`cost_s` or `error_abs`).
        field: &'static str,
        /// Human-readable, agent-actionable diagnosis.
        problem: String,
    },
    /// The oracle's bounded evidence budget is exhausted.
    CapacityExceeded {
        /// Bounded collection that is full.
        resource: &'static str,
        /// Maximum admitted entries.
        limit: usize,
    },
    /// A backend/model refused the otherwise well-shaped record.
    Backend {
        /// Human-readable, agent-actionable diagnosis.
        problem: String,
    },
}

impl fmt::Display for CostOracleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidEdge { problem } => write!(f, "invalid oracle edge: {problem}"),
            Self::InvalidMeasurement { field, problem } => {
                write!(f, "invalid oracle measurement {field}: {problem}")
            }
            Self::CapacityExceeded { resource, limit } => {
                write!(f, "oracle {resource} capacity {limit} is exhausted")
            }
            Self::Backend { problem } => write!(f, "oracle backend refused record: {problem}"),
        }
    }
}

impl core::error::Error for CostOracleError {}

/// Maximum distinct edges retained by the in-memory test/learning oracle.
pub const MAX_MEMORY_ORACLE_EDGES: usize = 4_096;
/// Maximum visible-ASCII bytes in a converter/oracle identity.
pub const MAX_ROUTER_ID_BYTES: usize = 256;

fn validate_identity(field: &str, value: &str) -> Result<(), String> {
    if value.is_empty() {
        return Err(format!("{field} identity is empty"));
    }
    if value.len() > MAX_ROUTER_ID_BYTES {
        return Err(format!(
            "{field} identity is {} bytes; maximum is {MAX_ROUTER_ID_BYTES}",
            value.len()
        ));
    }
    if !value.bytes().all(|byte| (b'!'..=b'~').contains(&byte)) {
        return Err(format!("{field} identity must contain visible ASCII only"));
    }
    Ok(())
}

/// In-memory running-mean oracle (tests, in-process learning).
#[derive(Debug, Clone, Default)]
pub struct MemoryCostOracle {
    rows: BTreeMap<String, (f64, f64, u32)>,
}

impl MemoryCostOracle {
    /// An empty oracle.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn record_one(
        &mut self,
        observation: &ValidatedEdgeObservation,
    ) -> Result<(), CostOracleError> {
        let edge = observation.edge();
        let cost_s = observation.cost_s();
        let error_abs = observation.conservative_error_abs();
        validate_identity("edge", edge)
            .map_err(|problem| CostOracleError::InvalidEdge { problem })?;
        if !cost_s.is_finite() || cost_s <= 0.0 {
            return Err(CostOracleError::InvalidMeasurement {
                field: "cost_s",
                problem: "must be positive and finite".to_string(),
            });
        }
        if !error_abs.is_finite() || error_abs < 0.0 {
            return Err(CostOracleError::InvalidMeasurement {
                field: "error_abs",
                problem: "must be nonnegative and finite".to_string(),
            });
        }
        let prior = self.rows.get(edge).copied();
        if prior.is_none() && self.rows.len() >= MAX_MEMORY_ORACLE_EDGES {
            return Err(CostOracleError::CapacityExceeded {
                resource: "edges",
                limit: MAX_MEMORY_ORACLE_EDGES,
            });
        }
        let (cost_sum, error_sum, count) = prior.unwrap_or((0.0, 0.0, 0));
        let next_cost = cost_sum + cost_s;
        let next_error = error_sum + error_abs;
        let next_count = count
            .checked_add(1)
            .ok_or(CostOracleError::CapacityExceeded {
                resource: "observations_per_edge",
                limit: u32::MAX as usize,
            })?;
        if !next_cost.is_finite() {
            return Err(CostOracleError::Backend {
                problem: "cost accumulator overflowed".to_string(),
            });
        }
        if !next_error.is_finite() {
            return Err(CostOracleError::Backend {
                problem: "error accumulator overflowed".to_string(),
            });
        }
        self.rows
            .insert(edge.to_string(), (next_cost, next_error, next_count));
        Ok(())
    }
}

impl CostOracle for MemoryCostOracle {
    fn measured_cost_s(&self, edge: &str) -> Option<f64> {
        self.rows.get(edge).map(|(c, _, n)| c / f64::from(*n))
    }

    fn measured_error_abs(&self, edge: &str) -> Option<f64> {
        self.rows.get(edge).map(|(_, e, n)| e / f64::from(*n))
    }

    fn record_batch(
        &mut self,
        observations: &[ValidatedEdgeObservation],
    ) -> Result<(), CostOracleError> {
        let mut candidate = self.clone();
        for observation in observations {
            candidate.record_one(observation)?;
        }
        *self = candidate;
        Ok(())
    }
}

/// A routing request: reach `to` from `from` with composed absolute error
/// ≤ `max_abs_error` and predicted cost ≤ `max_cost_s`; `scale` is the
/// reference magnitude that grounds relative error models.
#[derive(Debug, Clone, PartialEq)]
pub struct RouteRequest {
    /// Source chart kind.
    pub from: String,
    /// Destination chart kind.
    pub to: String,
    /// Reference magnitude for relative→absolute error grounding.
    pub scale: f64,
    /// Composed absolute error budget.
    pub max_abs_error: f64,
    /// Predicted cost budget, seconds.
    pub max_cost_s: f64,
}

/// A winning conversion chain.
#[derive(Debug, Clone, PartialEq)]
pub struct RoutePlan {
    /// Edge names in execution order.
    edges: Vec<String>,
    /// Exact registered specs that authorized the plan.
    planned_specs: Vec<ConverterSpec>,
    /// Exact request/budgets under which the plan was admitted.
    request: RouteRequest,
    /// Predicted total cost, seconds.
    predicted_cost_s: f64,
    /// Composed absolute error bound.
    composed_abs_error: f64,
    /// True when every edge is certificate-backed.
    all_certified: bool,
}

impl RoutePlan {
    /// Planned converter identities in execution order.
    #[must_use]
    pub fn edges(&self) -> &[String] {
        &self.edges
    }

    /// Predicted total wall cost in seconds.
    #[must_use]
    pub fn predicted_cost_s(&self) -> f64 {
        self.predicted_cost_s
    }

    /// Planned conservative absolute-error bound.
    #[must_use]
    pub fn composed_abs_error(&self) -> f64 {
        self.composed_abs_error
    }

    /// Whether every planned converter requires certified evidence.
    #[must_use]
    pub fn all_certified(&self) -> bool {
        self.all_certified
    }

    /// Exact request that admitted this plan.
    #[must_use]
    pub fn request(&self) -> &RouteRequest {
        &self.request
    }
}

/// One Pareto-optimal candidate the planner considered (explainability).
#[derive(Debug, Clone, PartialEq)]
pub struct RouteCandidate {
    /// Edge names in order.
    pub edges: Vec<String>,
    /// Predicted cost.
    pub cost_s: f64,
    /// Composed error.
    pub abs_error: f64,
    /// Number of uncertified edges (0 = fully certified).
    pub uncertified_edges: u32,
    /// Whether it satisfied both budgets.
    pub admissible: bool,
}

/// The full explanation of a routing decision (agent-facing; the hook
/// fs-ir admission uses for chart-compatibility checking).
#[derive(Debug, Clone)]
pub struct RouteExplanation {
    /// Every Pareto-optimal candidate at the target, deterministic order.
    pub candidates: Vec<RouteCandidate>,
    /// Index into `candidates` of the winner, if any was admissible.
    pub winner: Option<usize>,
    /// Why the winner won (or why nothing did).
    pub reason: String,
}

/// Which budget made a request infeasible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Binding {
    /// No path connects the kinds at all.
    NoPath,
    /// The error budget binds (some path fits the cost budget).
    Error,
    /// The cost budget binds (some path fits the error budget).
    Cost,
    /// Every path violates both budgets.
    Both,
}

/// A refusal that teaches (P10): the binding constraint and the cheapest
/// relaxations that would admit a path.
#[derive(Debug, Clone, PartialEq)]
pub struct RouteRefusal {
    /// The binding constraint.
    pub binding: Binding,
    /// Smallest composed error any path achieves (ignoring cost).
    pub best_abs_error: Option<f64>,
    /// Smallest predicted cost any path achieves (ignoring error).
    pub best_cost_s: Option<f64>,
    /// Ranked textual fixes.
    pub fixes: Vec<String>,
}

impl fmt::Display for RouteRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RouteInfeasible: binding={:?}; best achievable error={:?}, cost={:?}s; fixes: {}",
            self.binding,
            self.best_abs_error,
            self.best_cost_s,
            self.fixes.join("; ")
        )
    }
}

impl core::error::Error for RouteRefusal {}

/// Planning failure: malformed authority is distinct from ordinary
/// infeasibility.
#[derive(Debug, Clone, PartialEq)]
pub enum RoutePlanError {
    /// A request field is malformed or outside its finite domain.
    InvalidRequest {
        /// Rejected field.
        field: &'static str,
        /// Actionable diagnosis.
        problem: String,
    },
    /// An oracle returned a malformed measurement.
    InvalidOracle {
        /// Converter whose learned value was rejected.
        edge: String,
        /// Rejected oracle field.
        field: &'static str,
        /// Actionable diagnosis.
        problem: String,
    },
    /// Finite inputs overflowed while composing a path.
    Arithmetic {
        /// Converter at which composition failed.
        edge: String,
        /// Actionable diagnosis.
        problem: String,
    },
    /// No route satisfies the valid request.
    Infeasible(RouteRefusal),
}

impl fmt::Display for RoutePlanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRequest { field, problem } => {
                write!(f, "invalid route request {field}: {problem}")
            }
            Self::InvalidOracle {
                edge,
                field,
                problem,
            } => write!(f, "invalid oracle {field} for edge {edge:?}: {problem}"),
            Self::Arithmetic { edge, problem } => {
                write!(f, "route arithmetic failed at edge {edge:?}: {problem}")
            }
            Self::Infeasible(refusal) => refusal.fmt(f),
        }
    }
}

impl core::error::Error for RoutePlanError {}

impl From<RouteRefusal> for RoutePlanError {
    fn from(refusal: RouteRefusal) -> Self {
        Self::Infeasible(refusal)
    }
}

/// Errors from registry misuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RouterError {
    /// Duplicate edge name.
    DuplicateEdge(String),
    /// Nonsensical spec (empty names, negative cost/error).
    InvalidSpec(String),
}

impl fmt::Display for RouterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RouterError::DuplicateEdge(name) => {
                write!(f, "converter edge {name:?} is already registered")
            }
            RouterError::InvalidSpec(why) => write!(f, "invalid converter spec: {why}"),
        }
    }
}

/// One Pareto label during search.
#[derive(Debug, Clone)]
struct Label {
    cost: f64,
    err: f64,
    uncertified: u32,
    path: Vec<usize>,
    nodes: Vec<String>,
}

struct ExecutedRoute {
    receipt: Evidence<f64>,
    measured_cost_s: f64,
    observations: Vec<ValidatedEdgeObservation>,
}

impl Label {
    /// Pareto dominance (≤ on all axes, < on at least one).
    fn dominates(&self, other: &Label) -> bool {
        let le = self.cost <= other.cost
            && self.err <= other.err
            && self.uncertified <= other.uncertified;
        let lt =
            self.cost < other.cost || self.err < other.err || self.uncertified < other.uncertified;
        le && lt
    }
}

/// The conversion-graph registry and planner.
#[derive(Debug, Default)]
pub struct Router {
    edges: Vec<ConverterSpec>,
}

impl Router {
    /// An empty router.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a converter edge.
    ///
    /// # Errors
    /// [`RouterError`] on duplicate names or nonsensical specs.
    pub fn register(&mut self, spec: ConverterSpec) -> Result<(), RouterError> {
        for (field, value) in [
            ("name", spec.name.as_str()),
            ("from", spec.from.as_str()),
            ("to", spec.to.as_str()),
        ] {
            validate_identity(field, value).map_err(RouterError::InvalidSpec)?;
        }
        if !spec.base_cost_s.is_finite() || spec.base_cost_s < 0.0 {
            return Err(RouterError::InvalidSpec(format!(
                "base_cost_s must be finite and >= 0, got {}",
                spec.base_cost_s
            )));
        }
        if let ErrorModel::AdditiveAbs(a) | ErrorModel::MultiplicativeRel(a) = spec.error
            && (!a.is_finite() || a < 0.0)
        {
            return Err(RouterError::InvalidSpec(format!(
                "error magnitude must be finite and >= 0, got {a}"
            )));
        }
        if self.edges.iter().any(|e| e.name == spec.name) {
            return Err(RouterError::DuplicateEdge(spec.name));
        }
        self.edges.push(spec);
        // Deterministic expansion order regardless of registration order.
        self.edges.sort_by(|a, b| a.name.cmp(&b.name));
        Ok(())
    }

    /// The registered edges (deterministic order).
    #[must_use]
    pub fn edges(&self) -> &[ConverterSpec] {
        &self.edges
    }

    /// Effective cost of an edge: measured mean when the oracle has one,
    /// else the a-priori base cost.
    fn edge_cost(spec: &ConverterSpec, oracle: &dyn CostOracle) -> Result<f64, RoutePlanError> {
        let cost = oracle
            .measured_cost_s(&spec.name)
            .unwrap_or(spec.base_cost_s);
        if !cost.is_finite() || cost < 0.0 {
            return Err(RoutePlanError::InvalidOracle {
                edge: spec.name.clone(),
                field: "measured_cost_s",
                problem: "must be finite and nonnegative".to_string(),
            });
        }
        Ok(cost)
    }

    /// Effective error model: learned measurements may replace the
    /// declared magnitude ONLY on uncertified additive edges (estimates
    /// learn; certificates are never learned away).
    fn edge_error(
        spec: &ConverterSpec,
        oracle: &dyn CostOracle,
    ) -> Result<ErrorModel, RoutePlanError> {
        if spec.certified {
            return Ok(spec.error);
        }
        match (spec.error, oracle.measured_error_abs(&spec.name)) {
            (ErrorModel::AdditiveAbs(_), Some(measured)) => {
                if !measured.is_finite() || measured < 0.0 {
                    return Err(RoutePlanError::InvalidOracle {
                        edge: spec.name.clone(),
                        field: "measured_error_abs",
                        problem: "must be finite and nonnegative".to_string(),
                    });
                }
                Ok(ErrorModel::AdditiveAbs(measured))
            }
            (declared, _) => Ok(declared),
        }
    }

    fn validate_request(req: &RouteRequest) -> Result<(), RoutePlanError> {
        for (field, value) in [("from", req.from.as_str()), ("to", req.to.as_str())] {
            validate_identity(field, value)
                .map_err(|problem| RoutePlanError::InvalidRequest { field, problem })?;
        }
        for (field, value) in [
            ("scale", req.scale),
            ("max_abs_error", req.max_abs_error),
            ("max_cost_s", req.max_cost_s),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(RoutePlanError::InvalidRequest {
                    field,
                    problem: format!("must be finite and nonnegative, got {value}"),
                });
            }
        }
        Ok(())
    }

    /// All Pareto-optimal simple paths from `req.from` to `req.to`
    /// (label-correcting; deterministic).
    fn pareto_front(
        &self,
        req: &RouteRequest,
        oracle: &dyn CostOracle,
        include: &dyn Fn(&ConverterSpec) -> bool,
    ) -> Result<Vec<Label>, RoutePlanError> {
        Self::validate_request(req)?;
        let start = Label {
            cost: 0.0,
            err: 0.0,
            uncertified: 0,
            path: Vec::new(),
            nodes: vec![req.from.clone()],
        };
        let mut fronts: BTreeMap<String, Vec<Label>> = BTreeMap::new();
        fronts.insert(req.from.clone(), vec![start.clone()]);
        let mut queue = vec![start];
        while let Some(label) = queue.pop() {
            let at = label
                .nodes
                .last()
                .expect("labels always have a node")
                .clone();
            for (idx, spec) in self.edges.iter().enumerate() {
                if !include(spec) || spec.from != at || label.nodes.iter().any(|n| n == &spec.to) {
                    continue; // not outgoing, or would revisit (simple paths)
                }
                let edge_cost = Self::edge_cost(spec, oracle)?;
                let cost = label.cost + edge_cost;
                if !cost.is_finite() {
                    return Err(RoutePlanError::Arithmetic {
                        edge: spec.name.clone(),
                        problem: "path cost overflowed the finite domain".to_string(),
                    });
                }
                let err = Self::edge_error(spec, oracle)?.compose(label.err, req.scale);
                if !err.is_finite() || err < 0.0 {
                    return Err(RoutePlanError::Arithmetic {
                        edge: spec.name.clone(),
                        problem: "path error composition left the finite nonnegative domain"
                            .to_string(),
                    });
                }
                let uncertified = label
                    .uncertified
                    .checked_add(u32::from(!spec.certified))
                    .ok_or_else(|| RoutePlanError::Arithmetic {
                        edge: spec.name.clone(),
                        problem: "uncertified-edge counter overflowed".to_string(),
                    })?;
                let next = Label {
                    cost,
                    err,
                    uncertified,
                    path: {
                        let mut p = label.path.clone();
                        p.push(idx);
                        p
                    },
                    nodes: {
                        let mut n = label.nodes.clone();
                        n.push(spec.to.clone());
                        n
                    },
                };
                let front = fronts.entry(spec.to.clone()).or_default();
                if front.iter().any(|l| l.dominates(&next)) {
                    continue;
                }
                front.retain(|l| !next.dominates(l));
                front.push(next.clone());
                queue.push(next);
            }
        }
        let mut result = fronts.remove(&req.to).unwrap_or_default();
        // Deterministic reporting order: error, then cost, then path names.
        result.sort_by(|a, b| {
            a.err
                .total_cmp(&b.err)
                .then(a.cost.total_cmp(&b.cost))
                .then_with(|| self.path_names(&a.path).cmp(&self.path_names(&b.path)))
        });
        Ok(result)
    }

    fn path_names(&self, path: &[usize]) -> Vec<String> {
        path.iter().map(|&i| self.edges[i].name.clone()).collect()
    }

    /// Solve a request: the cheapest admissible chain, certified edges
    /// preferred (winner rule in module docs).
    ///
    /// # Errors
    /// [`RouteRefusal`] naming the binding constraint and relaxations.
    pub fn plan(
        &self,
        req: &RouteRequest,
        oracle: &dyn CostOracle,
    ) -> Result<RoutePlan, RoutePlanError> {
        self.plan_with_edge_filter(req, oracle, |_| true)
    }

    /// Solve while excluding edges that do not satisfy `include`. The filter
    /// is deterministic policy input; it does not rebuild or clone the graph.
    ///
    /// # Errors
    /// Any malformed request/oracle/arithmetic failure or an infeasible
    /// filtered graph.
    pub fn plan_with_edge_filter(
        &self,
        req: &RouteRequest,
        oracle: &dyn CostOracle,
        include: impl Fn(&ConverterSpec) -> bool,
    ) -> Result<RoutePlan, RoutePlanError> {
        let front = self.pareto_front(req, oracle, &include)?;
        let admissible: Vec<&Label> = front
            .iter()
            .filter(|l| l.err <= req.max_abs_error && l.cost <= req.max_cost_s)
            .collect();
        if let Some(winner) = self.pick_winner(&admissible) {
            return Ok(RoutePlan {
                edges: self.path_names(&winner.path),
                planned_specs: winner
                    .path
                    .iter()
                    .map(|&index| self.edges[index].clone())
                    .collect(),
                request: req.clone(),
                predicted_cost_s: winner.cost,
                composed_abs_error: winner.err,
                all_certified: winner.uncertified == 0,
            });
        }
        Err(Self::refusal(req, &front).into())
    }

    fn pick_winner<'a>(&self, admissible: &[&'a Label]) -> Option<&'a Label> {
        admissible
            .iter()
            .min_by(|a, b| {
                a.uncertified
                    .cmp(&b.uncertified)
                    .then(a.cost.total_cmp(&b.cost))
                    .then(a.err.total_cmp(&b.err))
                    .then_with(|| self.path_names(&a.path).cmp(&self.path_names(&b.path)))
            })
            .copied()
    }

    fn refusal(req: &RouteRequest, front: &[Label]) -> RouteRefusal {
        if front.is_empty() {
            return RouteRefusal {
                binding: Binding::NoPath,
                best_abs_error: None,
                best_cost_s: None,
                fixes: vec![format!(
                    "no converter chain connects {:?} to {:?}; register a converter or \
                     request a reachable chart kind",
                    req.from, req.to
                )],
            };
        }
        let best_err = front.iter().map(|l| l.err).fold(f64::INFINITY, f64::min);
        let best_cost = front.iter().map(|l| l.cost).fold(f64::INFINITY, f64::min);
        let error_feasible = best_err <= req.max_abs_error;
        let cost_feasible = best_cost <= req.max_cost_s;
        let binding = match (error_feasible, cost_feasible) {
            (false, true) => Binding::Error,
            (true, false) => Binding::Cost,
            _ => Binding::Both,
        };
        let mut fixes = Vec::new();
        if !error_feasible {
            fixes.push(format!(
                "relax max_abs_error from {} to {best_err} (the best any chain achieves)",
                req.max_abs_error
            ));
        }
        if !cost_feasible {
            fixes.push(format!(
                "relax max_cost_s from {} to {best_cost} (the cheapest chain)",
                req.max_cost_s
            ));
        }
        if binding == Binding::Both {
            fixes.push(
                "no single relaxation suffices: the low-error and low-cost chains differ; \
                 consider registering a certified direct converter"
                    .to_string(),
            );
        }
        RouteRefusal {
            binding,
            best_abs_error: Some(best_err),
            best_cost_s: Some(best_cost),
            fixes,
        }
    }

    /// Full decision explanation: every Pareto candidate, the winner, and
    /// why (agent-facing; deterministic).
    pub fn explain(
        &self,
        req: &RouteRequest,
        oracle: &dyn CostOracle,
    ) -> Result<RouteExplanation, RoutePlanError> {
        let front = self.pareto_front(req, oracle, &|_| true)?;
        let candidates: Vec<RouteCandidate> = front
            .iter()
            .map(|l| RouteCandidate {
                edges: self.path_names(&l.path),
                cost_s: l.cost,
                abs_error: l.err,
                uncertified_edges: l.uncertified,
                admissible: l.err <= req.max_abs_error && l.cost <= req.max_cost_s,
            })
            .collect();
        let admissible: Vec<&Label> = front
            .iter()
            .filter(|l| l.err <= req.max_abs_error && l.cost <= req.max_cost_s)
            .collect();
        let winner_label = self.pick_winner(&admissible);
        let winner = winner_label.and_then(|w| front.iter().position(|l| l.path == w.path));
        let reason = match winner_label {
            Some(w) if w.uncertified == 0 => format!(
                "fully certified chain at minimal cost {}s within error budget",
                w.cost
            ),
            Some(w) => format!(
                "no fully certified admissible chain; cheapest admissible has {} uncertified \
                 edge(s) at cost {}s",
                w.uncertified, w.cost
            ),
            None => Self::refusal(req, &front).to_string(),
        };
        Ok(RouteExplanation {
            candidates,
            winner,
            reason,
        })
    }

    /// Execute a plan through per-edge runners: composes the edges'
    /// Evidence receipts (error enclosures add along a chain) and records
    /// actual cost/error into the oracle so later plans improve.
    ///
    /// # Errors
    /// The failing edge's name and its runner error.
    pub fn execute(
        &self,
        plan: &RoutePlan,
        runners: &BTreeMap<String, Box<dyn EdgeRunner>>,
        oracle: &mut dyn CostOracle,
        cx: &Cx<'_>,
    ) -> Result<ChainOutcome, ExecuteError> {
        self.validate_plan_binding(plan)?;
        let executed = Self::run_plan(plan, runners, cx)?;
        Self::record_and_validate(plan, executed, oracle)
    }

    fn validate_plan_binding(&self, plan: &RoutePlan) -> Result<(), ExecuteError> {
        if plan.edges.len() != plan.planned_specs.len() {
            return Err(ExecuteError {
                edge: "<plan>".to_string(),
                kind: ExecuteErrorKind::InvalidPlan,
                detail: "sealed plan edge/spec cardinality mismatch".to_string(),
            });
        }
        for (edge, planned_spec) in plan.edges.iter().zip(&plan.planned_specs) {
            if self.edges.iter().find(|spec| spec.name == *edge) != Some(planned_spec) {
                return Err(ExecuteError {
                    edge: edge.clone(),
                    kind: ExecuteErrorKind::InvalidPlan,
                    detail: "plan was minted for a different or changed converter registry"
                        .to_string(),
                });
            }
        }
        Ok(())
    }

    fn run_plan(
        plan: &RoutePlan,
        runners: &BTreeMap<String, Box<dyn EdgeRunner>>,
        cx: &Cx<'_>,
    ) -> Result<ExecutedRoute, ExecuteError> {
        let mut composed: Option<Evidence<f64>> = None;
        let mut total_cost = 0.0;
        let mut observations = Vec::with_capacity(plan.edges.len());
        for (edge, planned_spec) in plan.edges.iter().zip(&plan.planned_specs) {
            let (receipt, observation) = Self::run_edge(edge, planned_spec, runners, cx)?;
            total_cost += observation.cost_s;
            if !total_cost.is_finite() {
                return Err(ExecuteError {
                    edge: edge.clone(),
                    kind: ExecuteErrorKind::Arithmetic,
                    detail: "measured chain cost overflowed the finite domain".to_string(),
                });
            }
            composed = Some(Self::compose_actual(composed, receipt, edge)?);
            observations.push(observation);
        }
        let receipt = composed
            .unwrap_or_else(|| Evidence::exact(0.0, ProvenanceHash::of_bytes(b"identity-route")));
        Ok(ExecutedRoute {
            receipt,
            measured_cost_s: total_cost,
            observations,
        })
    }

    fn run_edge(
        edge: &str,
        planned_spec: &ConverterSpec,
        runners: &BTreeMap<String, Box<dyn EdgeRunner>>,
        cx: &Cx<'_>,
    ) -> Result<(Evidence<f64>, ValidatedEdgeObservation), ExecuteError> {
        let runner = runners.get(edge).ok_or_else(|| ExecuteError {
            edge: edge.to_string(),
            kind: ExecuteErrorKind::MissingRunner,
            detail: "no runner registered for this edge".to_string(),
        })?;
        let outcome = runner.run(cx).map_err(|detail| ExecuteError {
            edge: edge.to_string(),
            kind: ExecuteErrorKind::Runner,
            detail,
        })?;
        if planned_spec.certified && outcome.evidence_class() != EdgeEvidenceClass::Certified {
            return Err(ExecuteError {
                edge: edge.to_string(),
                kind: ExecuteErrorKind::Evidence,
                detail: "certified converter returned only estimated evidence".to_string(),
            });
        }
        let measured_cost_s = outcome.measured_cost_s();
        let evidence_class = outcome.evidence_class();
        let receipt = outcome.into_receipt();
        let observation = ValidatedEdgeObservation {
            edge: edge.to_string(),
            cost_s: measured_cost_s,
            conservative_error_abs: receipt.numerical.hi,
            evidence_class,
            provenance: receipt.provenance,
        };
        Ok((receipt, observation))
    }

    fn compose_actual(
        composed: Option<Evidence<f64>>,
        receipt: Evidence<f64>,
        edge: &str,
    ) -> Result<Evidence<f64>, ExecuteError> {
        let Some(acc) = composed else {
            return Ok(receipt);
        };
        let sum = acc.qoi + receipt.qoi;
        if !sum.is_finite() {
            return Err(ExecuteError {
                edge: edge.to_string(),
                kind: ExecuteErrorKind::Arithmetic,
                detail: "measured error QoI overflowed the finite domain".to_string(),
            });
        }
        let combined = Evidence::combine(Op::Add, &acc, &receipt, sum);
        if !combined.numerical.lo.is_finite()
            || !combined.numerical.hi.is_finite()
            || combined.numerical.lo > combined.numerical.hi
        {
            return Err(ExecuteError {
                edge: edge.to_string(),
                kind: ExecuteErrorKind::Arithmetic,
                detail: "measured error enclosure overflowed the finite domain".to_string(),
            });
        }
        Ok(combined)
    }

    fn record_and_validate(
        plan: &RoutePlan,
        executed: ExecutedRoute,
        oracle: &mut dyn CostOracle,
    ) -> Result<ChainOutcome, ExecuteError> {
        oracle
            .record_batch(&executed.observations)
            .map_err(|error| ExecuteError {
                edge: executed
                    .observations
                    .last()
                    .map_or_else(|| "<identity-route>".to_string(), |row| row.edge.clone()),
                kind: ExecuteErrorKind::OracleRecord,
                detail: error.to_string(),
            })?;
        if executed.receipt.numerical.hi > plan.composed_abs_error {
            return Err(ExecuteError {
                edge: plan
                    .edges
                    .last()
                    .cloned()
                    .unwrap_or_else(|| "<identity-route>".to_string()),
                kind: ExecuteErrorKind::PlanContradiction,
                detail: format!(
                    "actual conservative error {} exceeds planned bound {}",
                    executed.receipt.numerical.hi, plan.composed_abs_error
                ),
            });
        }
        if executed.measured_cost_s > plan.request.max_cost_s {
            return Err(ExecuteError {
                edge: plan
                    .edges
                    .last()
                    .cloned()
                    .unwrap_or_else(|| "<identity-route>".to_string()),
                kind: ExecuteErrorKind::PlanContradiction,
                detail: format!(
                    "actual measured cost {} exceeds request budget {}",
                    executed.measured_cost_s, plan.request.max_cost_s
                ),
            });
        }
        Ok(ChainOutcome {
            receipt: executed.receipt,
            measured_cost_s: executed.measured_cost_s,
        })
    }
}

/// Executes one converter edge; the receipt's QoI is the edge's ACHIEVED
/// absolute error with its enclosure (rigorous for certified edges).
pub trait EdgeRunner {
    /// Run the conversion under `cx`.
    ///
    /// # Errors
    /// A human/agent-readable failure description.
    fn run(&self, cx: &Cx<'_>) -> Result<EdgeOutcome, String>;
}

#[derive(Debug, Clone)]
enum EdgeReceipt {
    Certified(Certified<f64>),
    Estimated(Evidence<f64>),
}

/// One executed edge's actuals.
#[derive(Debug, Clone)]
pub struct EdgeOutcome {
    /// Achieved-error evidence (QoI = absolute error).
    receipt: EdgeReceipt,
    /// Measured wall cost, seconds.
    measured_cost_s: f64,
}

impl EdgeOutcome {
    /// Construct an outcome carrying rigorous certified error evidence.
    ///
    /// # Errors
    /// Refuses non-positive/non-finite cost or evidence outside the
    /// nonnegative absolute-error domain.
    pub fn certified(
        receipt: Certified<f64>,
        measured_cost_s: f64,
    ) -> Result<Self, EdgeOutcomeError> {
        Self::validate_cost(measured_cost_s)?;
        if !receipt.qoi.is_finite() || receipt.qoi < 0.0 || receipt.numerical.lo < 0.0 {
            return Err(EdgeOutcomeError::InvalidEvidence(
                "certified absolute error and its enclosure must be finite and nonnegative"
                    .to_string(),
            ));
        }
        Ok(Self {
            receipt: EdgeReceipt::Certified(receipt),
            measured_cost_s,
        })
    }

    /// Construct a bounded non-rigorous error estimate. This is the only
    /// estimated edge-evidence constructor; arbitrary public `Evidence`
    /// literals cannot enter routing authority.
    ///
    /// # Errors
    /// Refuses malformed bounds, a point outside its band, or invalid cost.
    pub fn estimated(
        error_abs: f64,
        lo: f64,
        hi: f64,
        provenance: ProvenanceHash,
        measured_cost_s: f64,
    ) -> Result<Self, EdgeOutcomeError> {
        Self::validate_cost(measured_cost_s)?;
        if !(error_abs.is_finite()
            && lo.is_finite()
            && hi.is_finite()
            && lo >= 0.0
            && lo <= error_abs
            && error_abs <= hi)
        {
            return Err(EdgeOutcomeError::InvalidEvidence(
                "estimated absolute error requires finite 0 <= lo <= value <= hi".to_string(),
            ));
        }
        let mut receipt = Evidence::exact(error_abs, provenance);
        receipt.numerical = NumericalCertificate::estimate(lo, hi);
        Ok(Self {
            receipt: EdgeReceipt::Estimated(receipt),
            measured_cost_s,
        })
    }

    fn validate_cost(measured_cost_s: f64) -> Result<(), EdgeOutcomeError> {
        if !measured_cost_s.is_finite() || measured_cost_s <= 0.0 {
            return Err(EdgeOutcomeError::InvalidCost);
        }
        Ok(())
    }

    fn evidence_class(&self) -> EdgeEvidenceClass {
        match self.receipt {
            EdgeReceipt::Certified(_) => EdgeEvidenceClass::Certified,
            EdgeReceipt::Estimated(_) => EdgeEvidenceClass::Estimated,
        }
    }

    fn measured_cost_s(&self) -> f64 {
        self.measured_cost_s
    }

    fn into_receipt(self) -> Evidence<f64> {
        match self.receipt {
            EdgeReceipt::Certified(receipt) => receipt.into_evidence(),
            EdgeReceipt::Estimated(receipt) => receipt,
        }
    }
}

/// Why a converter outcome could not enter the sealed execution path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EdgeOutcomeError {
    /// Wall cost is not positive and finite.
    InvalidCost,
    /// Error evidence is malformed for an absolute-error receipt.
    InvalidEvidence(String),
}

impl fmt::Display for EdgeOutcomeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCost => write!(f, "edge wall cost must be positive and finite"),
            Self::InvalidEvidence(problem) => write!(f, "invalid edge evidence: {problem}"),
        }
    }
}

impl core::error::Error for EdgeOutcomeError {}

/// A fully executed chain: the composed receipt and measured cost.
#[derive(Debug, Clone)]
pub struct ChainOutcome {
    /// Composed achieved-error evidence for the whole chain.
    pub receipt: Evidence<f64>,
    /// Total measured cost, seconds.
    pub measured_cost_s: f64,
}

/// An execution failure, attributed to its edge.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExecuteError {
    /// The failing edge name.
    pub edge: String,
    /// Which execution boundary refused.
    pub kind: ExecuteErrorKind,
    /// What went wrong.
    pub detail: String,
}

/// Structured execution failure class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecuteErrorKind {
    /// The opaque plan does not match this router's registered converters.
    InvalidPlan,
    /// No runner exists for the planned edge.
    MissingRunner,
    /// The edge runner itself failed.
    Runner,
    /// The runner's evidence class contradicts the converter contract.
    Evidence,
    /// Finite per-edge values overflowed during chain composition.
    Arithmetic,
    /// The edge ran, but its evidence was invalid and the oracle refused it.
    OracleRecord,
    /// Valid actual evidence contradicts the plan's admitted bound/budget.
    PlanContradiction,
}

impl fmt::Display for ExecuteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "route execution failed at edge {:?} ({:?}): {}",
            self.edge, self.kind, self.detail
        )
    }
}

impl core::error::Error for ExecuteError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::SphereChart;
    use crate::{Convert, ErrBudget, SampledSdf};
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                asupersync::types::Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn edge(
        from: &str,
        to: &str,
        name: &str,
        cost: f64,
        error: ErrorModel,
        certified: bool,
    ) -> ConverterSpec {
        ConverterSpec {
            from: from.to_string(),
            to: to.to_string(),
            name: name.to_string(),
            base_cost_s: cost,
            error,
            certified,
        }
    }

    /// A graph with parallel edges, a cycle, and mixed certificates.
    fn test_router() -> Router {
        let mut r = Router::new();
        for spec in [
            edge(
                "frep",
                "sdf",
                "frep->sdf/coarse",
                1.0,
                ErrorModel::AdditiveAbs(0.02),
                true,
            ),
            edge(
                "frep",
                "sdf",
                "frep->sdf/fine",
                4.0,
                ErrorModel::AdditiveAbs(0.005),
                true,
            ),
            edge(
                "sdf",
                "mesh",
                "sdf->mesh/dc",
                2.0,
                ErrorModel::AdditiveAbs(0.01),
                true,
            ),
            edge(
                "mesh",
                "sdf",
                "mesh->sdf/wind",
                1.5,
                ErrorModel::Exact,
                true,
            ),
            edge(
                "frep",
                "mesh",
                "frep->mesh/direct",
                2.5,
                ErrorModel::AdditiveAbs(0.05),
                false,
            ),
            edge(
                "mesh",
                "spline",
                "mesh->spline/fit",
                3.0,
                ErrorModel::MultiplicativeRel(0.01),
                false,
            ),
            edge(
                "sdf",
                "spline",
                "sdf->spline/fit",
                8.0,
                ErrorModel::AdditiveAbs(0.008),
                true,
            ),
        ] {
            r.register(spec).unwrap();
        }
        r
    }

    fn req(to: &str, max_err: f64, max_cost: f64) -> RouteRequest {
        RouteRequest {
            from: "frep".to_string(),
            to: to.to_string(),
            scale: 1.0,
            max_abs_error: max_err,
            max_cost_s: max_cost,
        }
    }

    /// Brute-force oracle: every simple path by DFS, Pareto-filtered.
    fn oracle_front(r: &Router, request: &RouteRequest) -> Vec<(f64, f64, u32)> {
        let oracle = MemoryCostOracle::new();
        let mut paths: Vec<(f64, f64, u32)> = Vec::new();
        let mut stack = vec![(
            request.from.clone(),
            0.0f64,
            0.0f64,
            0u32,
            vec![request.from.clone()],
        )];
        while let Some((at, cost, err, unc, nodes)) = stack.pop() {
            if at == request.to {
                paths.push((cost, err, unc));
                continue;
            }
            for spec in r.edges() {
                if spec.from == at && !nodes.contains(&spec.to) {
                    let mut n2 = nodes.clone();
                    n2.push(spec.to.clone());
                    stack.push((
                        spec.to.clone(),
                        cost + Router::edge_cost(spec, &oracle).expect("fixture oracle is valid"),
                        Router::edge_error(spec, &oracle)
                            .expect("fixture oracle is valid")
                            .compose(err, request.scale),
                        unc + u32::from(!spec.certified),
                        n2,
                    ));
                }
            }
        }
        // Pareto filter.
        let mut front: Vec<(f64, f64, u32)> = Vec::new();
        for p in &paths {
            if !paths.iter().any(|q| {
                (q.0 <= p.0 && q.1 <= p.1 && q.2 <= p.2) && (q.0 < p.0 || q.1 < p.1 || q.2 < p.2)
            }) {
                front.push(*p);
            }
        }
        front.sort_by(|a, b| {
            a.1.total_cmp(&b.1)
                .then(a.0.total_cmp(&b.0))
                .then(a.2.cmp(&b.2))
        });
        front.dedup();
        front
    }

    #[test]
    fn planner_front_matches_bruteforce_oracle() {
        let r = test_router();
        let oracle = MemoryCostOracle::new();
        for target in ["sdf", "mesh", "spline"] {
            let request = req(target, f64::MAX, f64::MAX);
            let mut got: Vec<(f64, f64, u32)> = r
                .explain(&request, &oracle)
                .expect("valid request and oracle")
                .candidates
                .iter()
                .map(|c| (c.cost_s, c.abs_error, c.uncertified_edges))
                .collect();
            got.sort_by(|a, b| {
                a.1.total_cmp(&b.1)
                    .then(a.0.total_cmp(&b.0))
                    .then(a.2.cmp(&b.2))
            });
            got.dedup();
            let want = oracle_front(&r, &request);
            assert_eq!(got, want, "Pareto front mismatch for target {target}");
            assert!(!want.is_empty());
        }
    }

    #[test]
    fn winner_prefers_certified_and_is_deterministic() {
        let r = test_router();
        let oracle = MemoryCostOracle::new();
        // Both the uncertified direct mesh edge (2.5s, 0.05) and the
        // certified 2-hop (3.0s, 0.03) are admissible: certified must win
        // despite costing more.
        let plan = r.plan(&req("mesh", 0.06, 10.0), &oracle).unwrap();
        assert!(
            plan.all_certified,
            "certified chain must be preferred: {plan:?}"
        );
        assert_eq!(plan.edges, vec!["frep->sdf/coarse", "sdf->mesh/dc"]);
        // Determinism across rebuilds with shuffled registration order.
        let mut r2 = Router::new();
        let mut specs: Vec<ConverterSpec> = r.edges().to_vec();
        specs.reverse();
        for s in specs {
            r2.register(s).unwrap();
        }
        for _ in 0..25 {
            assert_eq!(r2.plan(&req("mesh", 0.06, 10.0), &oracle).unwrap(), plan);
        }
    }

    #[test]
    fn learned_costs_change_the_plan() {
        let r = test_router();
        let mut oracle = MemoryCostOracle::new();
        // A-priori: coarse (1.0s) + dc (2.0s) = 3.0s certified chain wins
        // over fine (4.0s) + dc.
        let before = r.plan(&req("mesh", 0.06, 100.0), &oracle).unwrap();
        assert_eq!(before.edges[0], "frep->sdf/coarse");
        // Measurements reveal the coarse edge is actually slow (10s).
        let observation = ValidatedEdgeObservation {
            edge: "frep->sdf/coarse".to_string(),
            cost_s: 10.0,
            conservative_error_abs: 0.02,
            evidence_class: EdgeEvidenceClass::Certified,
            provenance: ProvenanceHash::of_bytes(b"learned-cost-fixture"),
        };
        oracle
            .record_batch(&[observation.clone(), observation])
            .unwrap();
        let after = r.plan(&req("mesh", 0.06, 100.0), &oracle).unwrap();
        assert_eq!(
            after.edges[0], "frep->sdf/fine",
            "measured costs must reroute the plan: {after:?}"
        );
    }

    #[test]
    fn refusals_name_the_binding_constraint_and_relaxations() {
        let r = test_router();
        let oracle = MemoryCostOracle::new();
        // Error binds: cost budget generous, error impossible.
        let RoutePlanError::Infeasible(e) =
            r.plan(&req("spline", 1e-9, 1000.0), &oracle).unwrap_err()
        else {
            panic!("valid request must fail as infeasible");
        };
        assert_eq!(e.binding, Binding::Error);
        assert!(e.fixes[0].contains("relax max_abs_error"), "{e}");
        // Cost binds.
        let RoutePlanError::Infeasible(e) = r.plan(&req("spline", 1.0, 0.5), &oracle).unwrap_err()
        else {
            panic!("valid request must fail as infeasible");
        };
        assert_eq!(e.binding, Binding::Cost);
        // No path at all.
        let RoutePlanError::Infeasible(e) = r
            .plan(
                &RouteRequest {
                    from: "spline".to_string(),
                    to: "frep".to_string(),
                    scale: 1.0,
                    max_abs_error: 1.0,
                    max_cost_s: 1.0,
                },
                &oracle,
            )
            .unwrap_err()
        else {
            panic!("valid request must fail as infeasible");
        };
        assert_eq!(e.binding, Binding::NoPath);
    }

    /// Synthetic runner with a known achieved error and cost.
    struct FixedRunner {
        err: f64,
        cost: f64,
    }

    impl EdgeRunner for FixedRunner {
        fn run(&self, _cx: &Cx<'_>) -> Result<EdgeOutcome, String> {
            let receipt = Evidence::enclosed(
                self.err,
                0.0,
                self.err,
                ProvenanceHash::of_bytes(b"fixed-runner"),
            )
            .certified()
            .map_err(|error| error.to_string())?;
            EdgeOutcome::certified(receipt, self.cost).map_err(|error| error.to_string())
        }
    }

    struct EstimatedRunner {
        value: f64,
        hi: f64,
        cost: f64,
    }

    impl EdgeRunner for EstimatedRunner {
        fn run(&self, _cx: &Cx<'_>) -> Result<EdgeOutcome, String> {
            EdgeOutcome::estimated(
                self.value,
                0.0,
                self.hi,
                ProvenanceHash::of_bytes(b"estimated-runner"),
                self.cost,
            )
            .map_err(|error| error.to_string())
        }
    }

    #[derive(Default)]
    struct MalformedOracle {
        cost: Option<f64>,
        error: Option<f64>,
    }

    impl CostOracle for MalformedOracle {
        fn measured_cost_s(&self, _edge: &str) -> Option<f64> {
            self.cost
        }

        fn measured_error_abs(&self, _edge: &str) -> Option<f64> {
            self.error
        }

        fn record_batch(
            &mut self,
            _observations: &[ValidatedEdgeObservation],
        ) -> Result<(), CostOracleError> {
            Ok(())
        }
    }

    #[test]
    fn planning_refuses_nonfinite_specs_requests_oracle_and_composition() {
        let mut router = Router::new();
        assert!(matches!(
            router.register(edge(
                "a",
                "b",
                "a->b/inf",
                f64::INFINITY,
                ErrorModel::Exact,
                true,
            )),
            Err(RouterError::InvalidSpec(_))
        ));
        router
            .register(edge(
                "a",
                "b",
                "a->b/max",
                f64::MAX,
                ErrorModel::Exact,
                true,
            ))
            .unwrap();
        router
            .register(edge(
                "b",
                "c",
                "b->c/max",
                f64::MAX,
                ErrorModel::Exact,
                true,
            ))
            .unwrap();
        let request = RouteRequest {
            from: "a".to_string(),
            to: "c".to_string(),
            scale: 1.0,
            max_abs_error: 1.0,
            max_cost_s: f64::MAX,
        };
        assert!(matches!(
            router.plan(&request, &MemoryCostOracle::new()),
            Err(RoutePlanError::Arithmetic { .. })
        ));
        let mut invalid_request = request.clone();
        invalid_request.max_abs_error = f64::INFINITY;
        assert!(matches!(
            router.plan(&invalid_request, &MemoryCostOracle::new()),
            Err(RoutePlanError::InvalidRequest {
                field: "max_abs_error",
                ..
            })
        ));
        let direct = Router {
            edges: vec![edge(
                "a",
                "c",
                "a->c/direct",
                1.0,
                ErrorModel::AdditiveAbs(0.1),
                false,
            )],
        };
        assert!(matches!(
            direct.plan(
                &request,
                &MalformedOracle {
                    cost: Some(f64::NAN),
                    error: None,
                }
            ),
            Err(RoutePlanError::InvalidOracle {
                field: "measured_cost_s",
                ..
            })
        ));
        assert!(matches!(
            direct.plan(
                &request,
                &MalformedOracle {
                    cost: None,
                    error: Some(-1.0),
                }
            ),
            Err(RoutePlanError::InvalidOracle {
                field: "measured_error_abs",
                ..
            })
        ));
    }

    #[test]
    fn sealed_outcomes_refuse_malformed_and_certified_to_estimated_laundering() {
        assert!(matches!(
            EdgeOutcome::estimated(0.1, 0.2, 0.3, ProvenanceHash::of_bytes(b"bad-band"), 1.0,),
            Err(EdgeOutcomeError::InvalidEvidence(_))
        ));
        with_cx(|cx| {
            let mut router = Router::new();
            router
                .register(edge(
                    "a",
                    "b",
                    "a->b/certified",
                    1.0,
                    ErrorModel::AdditiveAbs(0.1),
                    true,
                ))
                .unwrap();
            let request = RouteRequest {
                from: "a".to_string(),
                to: "b".to_string(),
                scale: 1.0,
                max_abs_error: 0.1,
                max_cost_s: 2.0,
            };
            let mut oracle = MemoryCostOracle::new();
            let plan = router.plan(&request, &oracle).unwrap();
            let runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::from([(
                "a->b/certified".to_string(),
                Box::new(EstimatedRunner {
                    value: 0.05,
                    hi: 0.1,
                    cost: 1.0,
                }) as Box<dyn EdgeRunner>,
            )]);
            let error = router
                .execute(&plan, &runners, &mut oracle, cx)
                .unwrap_err();
            assert_eq!(error.kind, ExecuteErrorKind::Evidence);
            assert_eq!(oracle.measured_error_abs("a->b/certified"), None);
        });
    }

    #[test]
    fn actual_contradictions_refuse_but_valid_observations_still_teach() {
        with_cx(|cx| {
            let mut router = Router::new();
            router
                .register(edge(
                    "a",
                    "b",
                    "a->b/estimated",
                    0.5,
                    ErrorModel::AdditiveAbs(0.01),
                    false,
                ))
                .unwrap();
            let request = RouteRequest {
                from: "a".to_string(),
                to: "b".to_string(),
                scale: 1.0,
                max_abs_error: 0.1,
                max_cost_s: 1.0,
            };
            let mut oracle = MemoryCostOracle::new();
            let plan = router.plan(&request, &oracle).unwrap();
            let runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::from([(
                "a->b/estimated".to_string(),
                Box::new(EstimatedRunner {
                    value: 0.005,
                    hi: 0.02,
                    cost: 2.0,
                }) as Box<dyn EdgeRunner>,
            )]);
            let error = router
                .execute(&plan, &runners, &mut oracle, cx)
                .unwrap_err();
            assert_eq!(error.kind, ExecuteErrorKind::PlanContradiction);
            assert_eq!(
                oracle.measured_error_abs("a->b/estimated"),
                Some(0.02),
                "learning must retain the conservative upper bound, not the 0.005 QoI"
            );
            assert_eq!(oracle.measured_cost_s("a->b/estimated"), Some(2.0));
        });
    }

    #[test]
    fn oracle_batches_are_atomic_and_plans_are_registry_bound() {
        let mut oracle = MemoryCostOracle::new();
        let valid = ValidatedEdgeObservation {
            edge: "a->b/valid".to_string(),
            cost_s: 1.0,
            conservative_error_abs: 0.1,
            evidence_class: EdgeEvidenceClass::Estimated,
            provenance: ProvenanceHash::of_bytes(b"valid"),
        };
        let invalid = ValidatedEdgeObservation {
            edge: "b->c/invalid".to_string(),
            cost_s: f64::NAN,
            conservative_error_abs: 0.1,
            evidence_class: EdgeEvidenceClass::Estimated,
            provenance: ProvenanceHash::of_bytes(b"invalid"),
        };
        assert!(oracle.record_batch(&[valid, invalid]).is_err());
        assert_eq!(oracle.measured_cost_s("a->b/valid"), None);

        with_cx(|cx| {
            let mut first = Router::new();
            first
                .register(edge(
                    "a",
                    "b",
                    "a->b/shared",
                    1.0,
                    ErrorModel::AdditiveAbs(0.1),
                    true,
                ))
                .unwrap();
            let request = RouteRequest {
                from: "a".to_string(),
                to: "b".to_string(),
                scale: 1.0,
                max_abs_error: 1.0,
                max_cost_s: 2.0,
            };
            let mut oracle = MemoryCostOracle::new();
            let plan = first.plan(&request, &oracle).unwrap();
            let mut changed = Router::new();
            changed
                .register(edge(
                    "a",
                    "b",
                    "a->b/shared",
                    1.0,
                    ErrorModel::AdditiveAbs(0.2),
                    true,
                ))
                .unwrap();
            let runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::from([(
                "a->b/shared".to_string(),
                Box::new(FixedRunner {
                    err: 0.05,
                    cost: 1.0,
                }) as Box<dyn EdgeRunner>,
            )]);
            let error = changed
                .execute(&plan, &runners, &mut oracle, cx)
                .unwrap_err();
            assert_eq!(error.kind, ExecuteErrorKind::InvalidPlan);
        });
    }

    #[test]
    fn g3_execution_composes_receipts_and_feeds_the_oracle() {
        with_cx(|cx| {
            let r = test_router();
            let mut oracle = MemoryCostOracle::new();
            let plan = r.plan(&req("mesh", 0.06, 10.0), &oracle).unwrap();
            let mut runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::new();
            runners.insert(
                "frep->sdf/coarse".to_string(),
                Box::new(FixedRunner {
                    err: 0.015,
                    cost: 0.9,
                }),
            );
            runners.insert(
                "sdf->mesh/dc".to_string(),
                Box::new(FixedRunner {
                    err: 0.007,
                    cost: 1.8,
                }),
            );
            let out = r.execute(&plan, &runners, &mut oracle, cx).unwrap();
            // Composed receipt: errors add; the chain's measured error must
            // sit inside the plan's composed certificate (G3).
            assert!((out.receipt.qoi - 0.022).abs() < 1e-15);
            assert!(
                out.receipt.qoi <= plan.composed_abs_error,
                "measured composed error {} exceeded the composed certificate {}",
                out.receipt.qoi,
                plan.composed_abs_error
            );
            assert!(out.receipt.numerical.hi >= out.receipt.qoi);
            // Actuals recorded: the oracle learned both edges.
            assert!(oracle.measured_cost_s("frep->sdf/coarse").is_some());
            assert!(oracle.measured_cost_s("sdf->mesh/dc").is_some());
            // Missing runner is a structured, edge-attributed failure.
            let empty: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::new();
            let err = r.execute(&plan, &empty, &mut oracle, cx).unwrap_err();
            assert_eq!(err.edge, "frep->sdf/coarse");
        });
    }

    #[test]
    fn execution_refuses_invalid_outcome_without_laundering_actuals() {
        with_cx(|cx| {
            let r = test_router();
            let oracle = MemoryCostOracle::new();
            let plan = r.plan(&req("sdf", 0.03, 10.0), &oracle).unwrap();
            let edge = plan.edges[0].clone();
            let mut runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::new();
            runners.insert(
                edge.clone(),
                Box::new(FixedRunner {
                    err: 0.01,
                    cost: f64::NAN,
                }),
            );
            let mut oracle = MemoryCostOracle::new();
            let error = r
                .execute(&plan, &runners, &mut oracle, cx)
                .expect_err("nonfinite actual must fail at the sealed outcome boundary");
            assert_eq!(error.edge, edge);
            assert_eq!(error.kind, ExecuteErrorKind::Runner);
            assert_eq!(oracle.measured_cost_s(&error.edge), None);
            assert_eq!(oracle.measured_error_abs(&error.edge), None);
        });
    }

    #[test]
    fn real_convert_edge_runs_under_the_router() {
        with_cx(|cx| {
            struct SphereToSdf;
            impl EdgeRunner for SphereToSdf {
                fn run(&self, cx: &Cx<'_>) -> Result<EdgeOutcome, String> {
                    let sphere = SphereChart {
                        center: crate::Point3::new(0.0, 0.0, 0.0),
                        radius: 1.0,
                    };
                    let start = std::time::Instant::now();
                    let converted: fs_evidence::Certified<SampledSdf> =
                        Convert::<SampledSdf>::convert(
                            &sphere,
                            ErrBudget { abs_sd_error: 0.05 },
                            cx,
                        )
                        .map_err(|d| d.to_string())?;
                    let receipt = Evidence {
                        value: converted.qoi,
                        qoi: converted.qoi,
                        numerical: converted.numerical,
                        statistical: converted.statistical,
                        model: converted.model.clone(),
                        sensitivity: converted.sensitivity.clone(),
                        provenance: converted.provenance,
                        adjoint_ref: converted.adjoint_ref,
                    }
                    .certified()
                    .map_err(|error| error.to_string())?;
                    EdgeOutcome::certified(receipt, start.elapsed().as_secs_f64())
                        .map_err(|error| error.to_string())
                }
            }
            let mut r = Router::new();
            r.register(edge(
                "frep",
                "sdf",
                "frep->sdf/sampled-real",
                1.0,
                ErrorModel::AdditiveAbs(0.05),
                true,
            ))
            .unwrap();
            let mut oracle = MemoryCostOracle::new();
            let plan = r
                .plan(
                    &RouteRequest {
                        from: "frep".to_string(),
                        to: "sdf".to_string(),
                        scale: 1.0,
                        max_abs_error: 0.05,
                        max_cost_s: 60.0,
                    },
                    &oracle,
                )
                .unwrap();
            let mut runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::new();
            runners.insert("frep->sdf/sampled-real".to_string(), Box::new(SphereToSdf));
            let out = r.execute(&plan, &runners, &mut oracle, cx).unwrap();
            assert!(
                out.receipt.qoi <= plan.composed_abs_error,
                "real conversion's achieved error {} exceeded the plan bound {}",
                out.receipt.qoi,
                plan.composed_abs_error
            );
            assert!(oracle.measured_cost_s("frep->sdf/sampled-real").is_some());
        });
    }
}
