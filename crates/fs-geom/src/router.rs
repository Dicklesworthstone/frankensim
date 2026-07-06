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

use fs_evidence::{Evidence, Op, ProvenanceHash};
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
#[derive(Debug, Clone)]
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
    /// Record one executed edge's actuals.
    fn record(&mut self, edge: &str, cost_s: f64, error_abs: f64);
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
}

impl CostOracle for MemoryCostOracle {
    fn measured_cost_s(&self, edge: &str) -> Option<f64> {
        self.rows.get(edge).map(|(c, _, n)| c / f64::from(*n))
    }

    fn measured_error_abs(&self, edge: &str) -> Option<f64> {
        self.rows.get(edge).map(|(_, e, n)| e / f64::from(*n))
    }

    fn record(&mut self, edge: &str, cost_s: f64, error_abs: f64) {
        let row = self.rows.entry(edge.to_string()).or_insert((0.0, 0.0, 0));
        row.0 += cost_s;
        row.1 += error_abs;
        row.2 += 1;
    }
}

/// A routing request: reach `to` from `from` with composed absolute error
/// ≤ `max_abs_error` and predicted cost ≤ `max_cost_s`; `scale` is the
/// reference magnitude that grounds relative error models.
#[derive(Debug, Clone)]
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
    pub edges: Vec<String>,
    /// Predicted total cost, seconds.
    pub predicted_cost_s: f64,
    /// Composed absolute error bound.
    pub composed_abs_error: f64,
    /// True when every edge is certificate-backed.
    pub all_certified: bool,
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
        if spec.name.is_empty() || spec.from.is_empty() || spec.to.is_empty() {
            return Err(RouterError::InvalidSpec("empty name/from/to".to_string()));
        }
        if spec.base_cost_s.is_nan() || spec.base_cost_s < 0.0 {
            return Err(RouterError::InvalidSpec(format!(
                "base_cost_s must be ≥ 0, got {}",
                spec.base_cost_s
            )));
        }
        if let ErrorModel::AdditiveAbs(a) | ErrorModel::MultiplicativeRel(a) = spec.error
            && (a.is_nan() || a < 0.0)
        {
            return Err(RouterError::InvalidSpec(format!(
                "error magnitude must be ≥ 0, got {a}"
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
    fn edge_cost(spec: &ConverterSpec, oracle: &dyn CostOracle) -> f64 {
        oracle
            .measured_cost_s(&spec.name)
            .unwrap_or(spec.base_cost_s)
    }

    /// Effective error model: learned measurements may replace the
    /// declared magnitude ONLY on uncertified additive edges (estimates
    /// learn; certificates are never learned away).
    fn edge_error(spec: &ConverterSpec, oracle: &dyn CostOracle) -> ErrorModel {
        if spec.certified {
            return spec.error;
        }
        match (spec.error, oracle.measured_error_abs(&spec.name)) {
            (ErrorModel::AdditiveAbs(_), Some(measured)) => ErrorModel::AdditiveAbs(measured),
            (declared, _) => declared,
        }
    }

    /// All Pareto-optimal simple paths from `req.from` to `req.to`
    /// (label-correcting; deterministic).
    fn pareto_front(&self, req: &RouteRequest, oracle: &dyn CostOracle) -> Vec<Label> {
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
                if spec.from != at || label.nodes.iter().any(|n| n == &spec.to) {
                    continue; // not outgoing, or would revisit (simple paths)
                }
                let next = Label {
                    cost: label.cost + Self::edge_cost(spec, oracle),
                    err: Self::edge_error(spec, oracle).compose(label.err, req.scale),
                    uncertified: label.uncertified + u32::from(!spec.certified),
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
        result
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
    ) -> Result<RoutePlan, RouteRefusal> {
        let front = self.pareto_front(req, oracle);
        let admissible: Vec<&Label> = front
            .iter()
            .filter(|l| l.err <= req.max_abs_error && l.cost <= req.max_cost_s)
            .collect();
        if let Some(winner) = self.pick_winner(&admissible) {
            return Ok(RoutePlan {
                edges: self.path_names(&winner.path),
                predicted_cost_s: winner.cost,
                composed_abs_error: winner.err,
                all_certified: winner.uncertified == 0,
            });
        }
        Err(Self::refusal(req, &front))
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
    #[must_use]
    pub fn explain(&self, req: &RouteRequest, oracle: &dyn CostOracle) -> RouteExplanation {
        let front = self.pareto_front(req, oracle);
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
        RouteExplanation {
            candidates,
            winner,
            reason,
        }
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
        let mut composed: Option<Evidence<f64>> = None;
        let mut total_cost = 0.0;
        for edge in &plan.edges {
            let runner = runners.get(edge).ok_or_else(|| ExecuteError {
                edge: edge.clone(),
                detail: "no runner registered for this edge".to_string(),
            })?;
            let outcome = runner.run(cx).map_err(|detail| ExecuteError {
                edge: edge.clone(),
                detail,
            })?;
            oracle.record(edge, outcome.measured_cost_s, outcome.receipt.qoi);
            total_cost += outcome.measured_cost_s;
            composed = Some(match composed {
                None => outcome.receipt,
                Some(acc) => {
                    let sum = acc.qoi + outcome.receipt.qoi;
                    Evidence::combine(Op::Add, &acc, &outcome.receipt, sum)
                }
            });
        }
        let receipt = composed
            .unwrap_or_else(|| Evidence::exact(0.0, ProvenanceHash::of_bytes(b"identity-route")));
        Ok(ChainOutcome {
            receipt,
            measured_cost_s: total_cost,
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

/// One executed edge's actuals.
#[derive(Debug, Clone)]
pub struct EdgeOutcome {
    /// Achieved-error evidence (QoI = absolute error).
    pub receipt: Evidence<f64>,
    /// Measured wall cost, seconds.
    pub measured_cost_s: f64,
}

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
    /// What went wrong.
    pub detail: String,
}

impl fmt::Display for ExecuteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "route execution failed at edge {:?}: {}",
            self.edge, self.detail
        )
    }
}

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
                        cost + Router::edge_cost(spec, &oracle),
                        Router::edge_error(spec, &oracle).compose(err, request.scale),
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
            let request = req(target, f64::INFINITY, f64::INFINITY);
            let mut got: Vec<(f64, f64, u32)> = r
                .explain(&request, &oracle)
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
        oracle.record("frep->sdf/coarse", 10.0, 0.02);
        oracle.record("frep->sdf/coarse", 10.0, 0.02);
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
        let e = r.plan(&req("spline", 1e-9, 1000.0), &oracle).unwrap_err();
        assert_eq!(e.binding, Binding::Error);
        assert!(e.fixes[0].contains("relax max_abs_error"), "{e}");
        // Cost binds.
        let e = r.plan(&req("spline", 1.0, 0.5), &oracle).unwrap_err();
        assert_eq!(e.binding, Binding::Cost);
        // No path at all.
        let e = r
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
            .unwrap_err();
        assert_eq!(e.binding, Binding::NoPath);
    }

    /// Synthetic runner with a known achieved error and cost.
    struct FixedRunner {
        err: f64,
        cost: f64,
    }

    impl EdgeRunner for FixedRunner {
        fn run(&self, _cx: &Cx<'_>) -> Result<EdgeOutcome, String> {
            Ok(EdgeOutcome {
                receipt: Evidence::enclosed(
                    self.err,
                    0.0,
                    self.err,
                    ProvenanceHash::of_bytes(b"fixed-runner"),
                ),
                measured_cost_s: self.cost,
            })
        }
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
                    Ok(EdgeOutcome {
                        receipt: Evidence {
                            value: converted.qoi,
                            qoi: converted.qoi,
                            numerical: converted.numerical,
                            statistical: converted.statistical,
                            model: converted.model.clone(),
                            sensitivity: converted.sensitivity.clone(),
                            provenance: converted.provenance,
                            adjoint_ref: converted.adjoint_ref,
                        },
                        measured_cost_s: start.elapsed().as_secs_f64(),
                    })
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
