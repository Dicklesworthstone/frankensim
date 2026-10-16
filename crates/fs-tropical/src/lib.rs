//! fs-tropical — tropical (max-plus) critical-path analytics. Layer: L1.
//!
//! Task-DAG timing lives in the MAX-PLUS SEMIRING (`⊕ = max`, `⊗ = +`): the
//! makespan of an acyclic pipeline is the LONGEST PATH, and the steady-state
//! throughput of a pipelined (cyclic) DAG is the TROPICAL EIGENVALUE — the
//! maximum cycle mean, whose supporting cycle is the CRITICAL CIRCUIT.
//!
//! Scheduling theory as an ONLINE INSTRUMENT, not a whiteboard aphorism:
//! - [`TaskDag::critical_path`] gives the makespan, the critical path, and the
//!   per-task SLACK (relaxing a task off the critical path buys NOTHING);
//! - [`TaskDag::tune_next`] ranks a certified unique critical path so the
//!   autotuner never claims one tied path can move the makespan, and
//!   [`TaskDag::bottleneck`] names the top task when one exists;
//! - [`max_cycle_mean`] (Karp) computes the tropical eigenvalue and
//!   [`critical_circuit`] recovers its cycle, cross-checked against a
//!   brute-force enumeration;
//! - [`RecommendationAudit`] tracks recommendation → realized outcome so the
//!   instrument only earns default-on once its advice is validated ([M]).
//!
//! Deterministic; no dependencies.

/// The max-plus additive identity (`−∞`): "no path / no edge".
pub const NEG_INF: f64 = f64::NEG_INFINITY;

/// Maximum tasks admitted to one critical-path analysis.
pub const MAX_TASK_DAG_NODES: usize = 4_096;
/// Maximum precedence edges admitted to one critical-path analysis.
pub const MAX_TASK_DAG_EDGES: usize = 65_536;
/// Maximum dense matrix dimension admitted to Karp's cubic algorithm.
pub const MAX_TROPICAL_MATRIX_NODES: usize = 256;
/// Maximum dimension admitted to factorial simple-cycle enumeration.
pub const MAX_BRUTE_FORCE_NODES: usize = 10;
/// Maximum realized-outcome rows retained by one recommendation audit.
pub const MAX_RECOMMENDATION_AUDIT_ENTRIES: usize = 65_536;

/// `a ⊕ b = max(a, b)`.
#[must_use]
pub fn oplus(a: f64, b: f64) -> f64 {
    a.max(b)
}

/// `a ⊗ b = a + b` (absorbing at `−∞`).
#[must_use]
pub fn otimes(a: f64, b: f64) -> f64 {
    if a == NEG_INF || b == NEG_INF {
        NEG_INF
    } else {
        a + b
    }
}

/// A structured failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TropicalError {
    /// The DAG has a cycle (no topological order / no finite makespan).
    Cyclic,
    /// An edge references a task index out of range.
    BadEdge {
        /// The offending endpoint.
        node: usize,
    },
    /// A task latency is NaN, infinite, or negative.
    InvalidLatency {
        /// Offending task.
        node: usize,
        /// Stable domain diagnosis.
        reason: &'static str,
    },
    /// One recommendation-audit field is outside its admitted domain.
    InvalidAuditField {
        /// Stable field name.
        field: &'static str,
        /// Stable domain requirement.
        requirement: &'static str,
    },
    /// A deterministic work or storage cap was exceeded.
    ResourceLimit {
        /// Bounded resource.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Observed request.
        observed: usize,
    },
    /// A dense weight matrix is not square.
    NonSquareMatrix {
        /// Offending row.
        row: usize,
        /// Required row length.
        expected: usize,
        /// Actual row length.
        actual: usize,
    },
    /// A matrix weight is neither finite nor the `NEG_INF` no-edge sentinel.
    InvalidWeight {
        /// Matrix row.
        row: usize,
        /// Matrix column.
        col: usize,
    },
    /// Finite inputs overflowed a max-plus path calculation.
    NumericalOverflow {
        /// Stable operation name.
        operation: &'static str,
    },
}

impl core::fmt::Display for TropicalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cyclic => f.write_str("task graph contains a cycle"),
            Self::BadEdge { node } => write!(f, "task edge references missing node {node}"),
            Self::InvalidLatency { node, reason } => {
                write!(f, "task {node} latency {reason}")
            }
            Self::InvalidAuditField { field, requirement } => {
                write!(f, "recommendation audit {field} {requirement}")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "tropical {resource} limit {limit} exceeded by request {observed}"
            ),
            Self::NonSquareMatrix {
                row,
                expected,
                actual,
            } => write!(
                f,
                "tropical matrix row {row} has length {actual}; expected {expected}"
            ),
            Self::InvalidWeight { row, col } => write!(
                f,
                "tropical matrix weight at ({row},{col}) must be finite or NEG_INF"
            ),
            Self::NumericalOverflow { operation } => {
                write!(f, "tropical numerical overflow during {operation}")
            }
        }
    }
}

impl std::error::Error for TropicalError {}

/// A task DAG: per-task latencies + precedence edges `(from → to)`.
#[derive(Debug, Clone, PartialEq)]
pub struct TaskDag {
    latency: Vec<f64>,
    edges: Vec<(usize, usize)>,
}

struct DagTopology {
    adj: Vec<Vec<usize>>,
    preds: Vec<Vec<usize>>,
    order: Vec<usize>,
}

struct EarliestFinishes {
    values: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    best_pred: Vec<Option<usize>>,
}

/// A critical-path analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct CriticalPath {
    /// The total makespan (longest path length).
    pub makespan: f64,
    /// Directed-rounding lower bound for the real-valued max-plus makespan.
    pub makespan_lo: f64,
    /// Directed-rounding upper bound for the real-valued max-plus makespan.
    pub makespan_hi: f64,
    /// Whether directed bounds prove this is the unique critical path.
    pub path_is_unique: bool,
    /// The critical path (task indices, source → sink).
    pub path: Vec<usize>,
    /// Per-task slack (`0` ⇔ on the critical path).
    pub slack: Vec<f64>,
}

impl TaskDag {
    /// A DAG with the given task latencies and no edges.
    #[must_use]
    pub fn new(latency: Vec<f64>) -> TaskDag {
        TaskDag {
            latency,
            edges: Vec::new(),
        }
    }

    /// Add a precedence edge `from → to` (builder).
    #[must_use]
    pub fn with_edge(mut self, from: usize, to: usize) -> TaskDag {
        self.edges.push((from, to));
        self
    }

    /// Number of tasks.
    #[must_use]
    pub fn len(&self) -> usize {
        self.latency.len()
    }

    /// Is the DAG empty?
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.latency.is_empty()
    }

    fn validated_topology(&self) -> Result<DagTopology, TropicalError> {
        let n = self.latency.len();
        if n > MAX_TASK_DAG_NODES {
            return Err(TropicalError::ResourceLimit {
                resource: "task nodes",
                limit: MAX_TASK_DAG_NODES,
                observed: n,
            });
        }
        if self.edges.len() > MAX_TASK_DAG_EDGES {
            return Err(TropicalError::ResourceLimit {
                resource: "task edges",
                limit: MAX_TASK_DAG_EDGES,
                observed: self.edges.len(),
            });
        }
        for (node, latency) in self.latency.iter().copied().enumerate() {
            if !latency.is_finite() {
                return Err(TropicalError::InvalidLatency {
                    node,
                    reason: "must be finite",
                });
            }
            if latency < 0.0 {
                return Err(TropicalError::InvalidLatency {
                    node,
                    reason: "must be non-negative",
                });
            }
        }
        for &(from, to) in &self.edges {
            if from >= n {
                return Err(TropicalError::BadEdge { node: from });
            }
            if to >= n {
                return Err(TropicalError::BadEdge { node: to });
            }
        }

        let mut adj = vec![Vec::new(); n];
        let mut preds = vec![Vec::new(); n];
        let mut indegree = vec![0usize; n];
        for &(from, to) in &self.edges {
            adj[from].push(to);
            preds[to].push(from);
            indegree[to] += 1;
        }
        let mut queue: Vec<usize> = (0..n).filter(|&node| indegree[node] == 0).collect();
        let mut order = Vec::with_capacity(n);
        while let Some(node) = queue.pop() {
            order.push(node);
            for &successor in &adj[node] {
                indegree[successor] -= 1;
                if indegree[successor] == 0 {
                    queue.push(successor);
                }
            }
        }
        if order.len() != n {
            return Err(TropicalError::Cyclic);
        }
        Ok(DagTopology { adj, preds, order })
    }

    fn earliest_finishes(
        &self,
        preds: &[Vec<usize>],
        order: &[usize],
    ) -> Result<EarliestFinishes, TropicalError> {
        let n = self.latency.len();
        let mut values = vec![0.0_f64; n];
        let mut lower = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut best_pred = vec![None; n];
        for &node in order {
            // Avoid a 0.0 bias: a zero-finish predecessor still belongs in the
            // source-to-sink witness even though it does not change makespan.
            let mut best_finish = f64::NEG_INFINITY;
            for &predecessor in &preds[node] {
                if values[predecessor] > best_finish {
                    best_finish = values[predecessor];
                    best_pred[node] = Some(predecessor);
                }
            }
            values[node] = best_finish.max(0.0) + self.latency[node];
            if !values[node].is_finite() {
                return Err(TropicalError::NumericalOverflow {
                    operation: "critical-path accumulation",
                });
            }
            let predecessor_lo = preds[node]
                .iter()
                .map(|&predecessor| lower[predecessor])
                .fold(0.0_f64, f64::max);
            let predecessor_hi = preds[node]
                .iter()
                .map(|&predecessor| upper[predecessor])
                .fold(0.0_f64, f64::max);
            lower[node] = (predecessor_lo + self.latency[node]).next_down();
            upper[node] = (predecessor_hi + self.latency[node]).next_up();
        }
        Ok(EarliestFinishes {
            values,
            lower,
            upper,
            best_pred,
        })
    }

    /// Compute the critical path (CPM longest-path in the max-plus semiring).
    ///
    /// # Errors
    /// Refuses invalid latencies, bounded-resource excess, out-of-range edges,
    /// cycles, and finite path-sum overflow.
    pub fn critical_path(&self) -> Result<CriticalPath, TropicalError> {
        let n = self.latency.len();
        let DagTopology { adj, preds, order } = self.validated_topology()?;
        if n == 0 {
            return Ok(CriticalPath {
                makespan: 0.0,
                makespan_lo: 0.0,
                makespan_hi: 0.0,
                path_is_unique: false,
                path: Vec::new(),
                slack: Vec::new(),
            });
        }
        let EarliestFinishes {
            values: ef,
            lower: ef_lo,
            upper: ef_hi,
            best_pred,
        } = self.earliest_finishes(&preds, &order)?;
        // A complete precedence path ends at a terminal task. With
        // non-negative latencies, the maximum terminal finish equals the
        // maximum finish over all tasks, while selecting among terminals keeps
        // a zero-latency tail in the returned source-to-sink witness.
        let terminals: Vec<usize> = (0..n).filter(|&node| adj[node].is_empty()).collect();
        let makespan = terminals
            .iter()
            .map(|&node| ef[node])
            .fold(0.0_f64, f64::max);
        let makespan_lo = terminals
            .iter()
            .map(|&node| ef_lo[node])
            .fold(0.0_f64, f64::max);
        let makespan_hi = terminals
            .iter()
            .map(|&node| ef_hi[node])
            .fold(0.0_f64, f64::max);
        let sink = terminals
            .iter()
            .copied()
            .max_by(|&a, &b| ef[a].total_cmp(&ef[b]))
            .ok_or(TropicalError::Cyclic)?;
        // backtrack the critical path.
        let mut path = vec![sink];
        let mut cur = sink;
        while let Some(p) = best_pred[cur] {
            path.push(p);
            cur = p;
        }
        path.reverse();
        let mut path_is_unique = terminals
            .iter()
            .copied()
            .filter(|&node| node != sink)
            .all(|node| ef_lo[sink] > ef_hi[node]);
        for pair in path.windows(2) {
            let chosen = pair[0];
            let child = pair[1];
            path_is_unique &= preds[child]
                .iter()
                .filter(|&&candidate| candidate != chosen)
                .all(|&candidate| ef_lo[chosen] > ef_hi[candidate]);
        }
        // latest finish (reverse pass) → slack.
        let mut lf = vec![makespan; n];
        for &u in order.iter().rev() {
            let mut latest = makespan;
            let mut has_succ = false;
            for &v in &adj[u] {
                has_succ = true;
                latest = latest.min(lf[v] - self.latency[v]);
            }
            if has_succ {
                lf[u] = latest;
            }
        }
        let slack: Vec<f64> = (0..n).map(|i| lf[i] - ef[i]).collect();
        Ok(CriticalPath {
            makespan,
            makespan_lo,
            makespan_hi,
            path_is_unique,
            path,
            slack,
        })
    }

    /// Positive-latency tasks on the certified unique critical path, ranked by
    /// latency. When directed bounds cannot prove a unique path, returns an
    /// empty advisory list rather than claiming one task can reduce the
    /// makespan.
    pub fn tune_next(&self) -> Result<Vec<usize>, TropicalError> {
        let cp = self.critical_path()?;
        if !cp.path_is_unique {
            return Ok(Vec::new());
        }
        let mut critical: Vec<usize> = cp
            .path
            .into_iter()
            .filter(|&task| self.latency[task] > 0.0)
            .collect();
        critical.sort_by(|&a, &b| self.latency[b].total_cmp(&self.latency[a]));
        Ok(critical)
    }

    /// The single strictly highest positive-latency bottleneck on a certified
    /// unique critical path. Tied top tasks, all-zero paths, and tied or
    /// rounding-ambiguous critical paths return `None` rather than naming a
    /// non-unique or ineffective target.
    pub fn bottleneck(&self) -> Result<Option<usize>, TropicalError> {
        let ranked = self.tune_next()?;
        let Some(&first) = ranked.first() else {
            return Ok(None);
        };
        if ranked
            .get(1)
            .is_some_and(|&second| self.latency[first].total_cmp(&self.latency[second]).is_eq())
        {
            return Ok(None);
        }
        Ok(Some(first))
    }
}

/// The tropical eigenvalue (maximum cycle mean) of a max-plus weight matrix
/// (`matrix[u][v]` = edge weight `u → v`, [`NEG_INF`] if absent). `Ok(None)`
/// if there is no cycle. Karp's algorithm.
// A max-plus matrix kernel: `d[k][v]` is inherently 2D-indexed by node column,
// so index loops are the correct, readable form here.
#[allow(clippy::needless_range_loop)]
pub fn max_cycle_mean(matrix: &[Vec<f64>]) -> Result<Option<f64>, TropicalError> {
    let n = matrix.len();
    validate_weight_matrix(matrix, MAX_TROPICAL_MATRIX_NODES, "matrix nodes")?;
    if n == 0 {
        return Ok(None);
    }
    // d[k][v] = longest path of EXACTLY k edges ending at v (from any start).
    let mut d = vec![vec![NEG_INF; n]; n + 1];
    for v in 0..n {
        d[0][v] = 0.0;
    }
    for k in 1..=n {
        for v in 0..n {
            let mut best = NEG_INF;
            for (u, row) in matrix.iter().enumerate() {
                if d[k - 1][u] != NEG_INF && row[v] != NEG_INF {
                    let candidate = d[k - 1][u] + row[v];
                    if !candidate.is_finite() {
                        return Err(TropicalError::NumericalOverflow {
                            operation: "Karp path accumulation",
                        });
                    }
                    best = best.max(candidate);
                }
            }
            d[k][v] = best;
        }
    }
    // λ = max_v min_{0<=k<n} (d[n][v] − d[k][v]) / (n − k).
    let mut lambda = NEG_INF;
    for v in 0..n {
        if d[n][v] == NEG_INF {
            continue;
        }
        let mut worst = f64::INFINITY;
        for k in 0..n {
            if d[k][v] == NEG_INF {
                continue;
            }
            let mean = (d[n][v] - d[k][v]) / (n - k) as f64;
            if !mean.is_finite() {
                return Err(TropicalError::NumericalOverflow {
                    operation: "Karp cycle mean",
                });
            }
            worst = worst.min(mean);
        }
        if worst.is_finite() {
            lambda = lambda.max(worst);
        }
    }
    if lambda == NEG_INF {
        Ok(None)
    } else {
        Ok(Some(lambda))
    }
}

fn validate_weight_matrix(
    matrix: &[Vec<f64>],
    limit: usize,
    resource: &'static str,
) -> Result<(), TropicalError> {
    let n = matrix.len();
    if n > limit {
        return Err(TropicalError::ResourceLimit {
            resource,
            limit,
            observed: n,
        });
    }
    for (row_index, row) in matrix.iter().enumerate() {
        if row.len() != n {
            return Err(TropicalError::NonSquareMatrix {
                row: row_index,
                expected: n,
                actual: row.len(),
            });
        }
        for (col, weight) in row.iter().copied().enumerate() {
            if weight != NEG_INF && !weight.is_finite() {
                return Err(TropicalError::InvalidWeight {
                    row: row_index,
                    col,
                });
            }
        }
    }
    Ok(())
}

/// The critical circuit (the cycle achieving the maximum cycle mean) by
/// brute-force simple-cycle enumeration — for the small timing matrices this
/// instrument runs on, and the ground truth [`max_cycle_mean`] is checked
/// against.
pub fn critical_circuit(matrix: &[Vec<f64>]) -> Result<Option<Vec<usize>>, TropicalError> {
    Ok(brute_force_max_cycle(matrix)?.map(|(cycle, _)| cycle))
}

/// Brute-force maximum-mean simple cycle: `(cycle nodes, mean)`.
pub fn brute_force_max_cycle(
    matrix: &[Vec<f64>],
) -> Result<Option<(Vec<usize>, f64)>, TropicalError> {
    let n = matrix.len();
    validate_weight_matrix(matrix, MAX_BRUTE_FORCE_NODES, "brute-force cycle nodes")?;
    let mut best: Option<(Vec<usize>, f64)> = None;
    // enumerate simple cycles whose smallest node is the start (dedup rotations).
    for start in 0..n {
        let mut stack = vec![start];
        let mut on_stack = vec![false; n];
        on_stack[start] = true;
        dfs_cycles(
            matrix,
            start,
            start,
            0.0,
            &mut stack,
            &mut on_stack,
            &mut best,
        )?;
        on_stack[start] = false;
    }
    Ok(best)
}

fn dfs_cycles(
    matrix: &[Vec<f64>],
    start: usize,
    u: usize,
    weight: f64,
    stack: &mut Vec<usize>,
    on_stack: &mut [bool],
    best: &mut Option<(Vec<usize>, f64)>,
) -> Result<(), TropicalError> {
    let n = matrix.len();
    for v in 0..n {
        let w = matrix[u][v];
        if w == NEG_INF {
            continue;
        }
        if v == start {
            // closed a cycle back to the start.
            let len = stack.len();
            let total = weight + w;
            if !total.is_finite() {
                return Err(TropicalError::NumericalOverflow {
                    operation: "brute-force cycle accumulation",
                });
            }
            let mean = total / len as f64;
            let is_better = best
                .as_ref()
                .is_none_or(|(_, old_mean)| mean.total_cmp(old_mean).is_gt());
            if is_better {
                *best = Some((stack.clone(), mean));
            }
        } else if v > start && !on_stack[v] {
            on_stack[v] = true;
            stack.push(v);
            let next_weight = weight + w;
            if !next_weight.is_finite() {
                return Err(TropicalError::NumericalOverflow {
                    operation: "brute-force path accumulation",
                });
            }
            dfs_cycles(matrix, start, v, next_weight, stack, on_stack, best)?;
            stack.pop();
            on_stack[v] = false;
        }
    }
    Ok(())
}

/// One recommendation-outcome record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AuditEntry {
    /// The task the instrument recommended tuning.
    pub recommended: usize,
    /// The makespan improvement it predicted.
    pub predicted_gain: f64,
    /// The makespan improvement actually realized.
    pub realized_gain: f64,
}

/// The promotion self-audit ([M]): tracks whether recommendations actually
/// moved the makespan, gating promotion from advisory-only.
#[derive(Debug, Clone, Default)]
pub struct RecommendationAudit {
    entries: Vec<AuditEntry>,
}

impl RecommendationAudit {
    /// An empty audit.
    #[must_use]
    pub fn new() -> RecommendationAudit {
        RecommendationAudit {
            entries: Vec::new(),
        }
    }

    /// Record a recommendation and its realized outcome.
    pub fn record(&mut self, entry: AuditEntry) -> Result<(), TropicalError> {
        if self.entries.len() >= MAX_RECOMMENDATION_AUDIT_ENTRIES {
            return Err(TropicalError::ResourceLimit {
                resource: "recommendation audit entries",
                limit: MAX_RECOMMENDATION_AUDIT_ENTRIES,
                observed: self.entries.len().saturating_add(1),
            });
        }
        if entry.recommended >= MAX_TASK_DAG_NODES {
            return Err(TropicalError::ResourceLimit {
                resource: "recommended task index",
                limit: MAX_TASK_DAG_NODES,
                observed: entry.recommended.saturating_add(1),
            });
        }
        if !entry.predicted_gain.is_finite() || entry.predicted_gain < 0.0 {
            return Err(TropicalError::InvalidAuditField {
                field: "predicted_gain",
                requirement: "must be finite and non-negative",
            });
        }
        if !entry.realized_gain.is_finite() {
            return Err(TropicalError::InvalidAuditField {
                field: "realized_gain",
                requirement: "must be finite",
            });
        }
        self.entries.push(entry);
        Ok(())
    }

    /// The hit rate: the fraction of recommendations that realized a positive
    /// makespan improvement.
    #[must_use]
    pub fn hit_rate(&self) -> f64 {
        if self.entries.is_empty() {
            return 0.0;
        }
        let hits = self
            .entries
            .iter()
            .filter(|e| e.realized_gain > 0.0)
            .count();
        hits as f64 / self.entries.len() as f64
    }

    /// Promote the instrument from advisory only when it has enough samples and
    /// a high enough hit rate.
    pub fn promoted(&self, min_samples: usize, min_hit_rate: f64) -> Result<bool, TropicalError> {
        if min_samples == 0 {
            return Err(TropicalError::InvalidAuditField {
                field: "min_samples",
                requirement: "must be at least one",
            });
        }
        if !min_hit_rate.is_finite() || !(0.0..=1.0).contains(&min_hit_rate) {
            return Err(TropicalError::InvalidAuditField {
                field: "min_hit_rate",
                requirement: "must be finite and in 0..=1",
            });
        }
        Ok(self.entries.len() >= min_samples && self.hit_rate() >= min_hit_rate)
    }
}
