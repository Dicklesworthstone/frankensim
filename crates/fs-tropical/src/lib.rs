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
//! - [`TaskDag::tune_next`] ranks the critical tasks so the autotuner knows
//!   which kernel to tune next, and [`TaskDag::bottleneck`] names the top one;
//! - [`max_cycle_mean`] (Karp) computes the tropical eigenvalue and
//!   [`critical_circuit`] recovers its cycle, cross-checked against a
//!   brute-force enumeration;
//! - [`RecommendationAudit`] tracks recommendation → realized outcome so the
//!   instrument only earns default-on once its advice is validated ([M]).
//!
//! Deterministic; no dependencies.

/// The max-plus additive identity (`−∞`): "no path / no edge".
pub const NEG_INF: f64 = f64::NEG_INFINITY;

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
}

/// A task DAG: per-task latencies + precedence edges `(from → to)`.
#[derive(Debug, Clone, PartialEq)]
pub struct TaskDag {
    latency: Vec<f64>,
    edges: Vec<(usize, usize)>,
}

/// A critical-path analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct CriticalPath {
    /// The total makespan (longest path length).
    pub makespan: f64,
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

    /// Compute the critical path (CPM longest-path in the max-plus semiring).
    ///
    /// # Errors
    /// [`TropicalError::Cyclic`] on a cycle; [`TropicalError::BadEdge`] on an
    /// out-of-range endpoint.
    pub fn critical_path(&self) -> Result<CriticalPath, TropicalError> {
        let n = self.latency.len();
        if n == 0 {
            return Ok(CriticalPath {
                makespan: 0.0,
                path: Vec::new(),
                slack: Vec::new(),
            });
        }
        let mut adj = vec![Vec::new(); n];
        let mut preds = vec![Vec::new(); n];
        let mut indeg = vec![0usize; n];
        for &(u, v) in &self.edges {
            if u >= n {
                return Err(TropicalError::BadEdge { node: u });
            }
            if v >= n {
                return Err(TropicalError::BadEdge { node: v });
            }
            adj[u].push(v);
            preds[v].push(u);
            indeg[v] += 1;
        }
        // Kahn topological order.
        let mut queue: Vec<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
        let mut order = Vec::with_capacity(n);
        let mut indeg2 = indeg;
        while let Some(u) = queue.pop() {
            order.push(u);
            for &v in &adj[u] {
                indeg2[v] -= 1;
                if indeg2[v] == 0 {
                    queue.push(v);
                }
            }
        }
        if order.len() != n {
            return Err(TropicalError::Cyclic);
        }
        // earliest finish + the max predecessor (for backtracking).
        let mut ef = vec![0.0_f64; n];
        let mut best_pred = vec![None; n];
        for &u in &order {
            // track the max-ef predecessor without the 0.0 bias, so a
            // zero-earliest-finish predecessor is still recorded on the path.
            let mut best_ef = f64::NEG_INFINITY;
            for &p in &preds[u] {
                if ef[p] > best_ef {
                    best_ef = ef[p];
                    best_pred[u] = Some(p);
                }
            }
            ef[u] = best_ef.max(0.0) + self.latency[u];
        }
        let makespan = ef.iter().copied().fold(0.0_f64, f64::max);
        let sink = ef
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map_or(0, |(i, _)| i);
        // backtrack the critical path.
        let mut path = vec![sink];
        let mut cur = sink;
        while let Some(p) = best_pred[cur] {
            path.push(p);
            cur = p;
        }
        path.reverse();
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
            path,
            slack,
        })
    }

    /// The critical tasks ranked by latency (the tune-next list): tuning a task
    /// OFF the critical path (positive slack) cannot move the makespan, so only
    /// zero-slack tasks appear, highest-latency first.
    #[must_use]
    pub fn tune_next(&self, cp: &CriticalPath) -> Vec<usize> {
        let mut critical: Vec<usize> = (0..self.latency.len())
            .filter(|&i| cp.slack[i].abs() < 1e-9)
            .collect();
        critical.sort_by(|&a, &b| {
            self.latency[b]
                .partial_cmp(&self.latency[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        critical
    }

    /// The single bottleneck: the highest-latency task on the critical path.
    #[must_use]
    pub fn bottleneck(&self, cp: &CriticalPath) -> Option<usize> {
        self.tune_next(cp).first().copied()
    }
}

/// The tropical eigenvalue (maximum cycle mean) of a max-plus weight matrix
/// (`matrix[u][v]` = edge weight `u → v`, [`NEG_INF`] if absent). `None` if
/// there is no cycle. Karp's algorithm.
// A max-plus matrix kernel: `d[k][v]` is inherently 2D-indexed by node column,
// so index loops are the correct, readable form here.
#[allow(clippy::needless_range_loop)]
#[must_use]
pub fn max_cycle_mean(matrix: &[Vec<f64>]) -> Option<f64> {
    let n = matrix.len();
    if n == 0 {
        return None;
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
                    best = best.max(d[k - 1][u] + row[v]);
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
            worst = worst.min(mean);
        }
        if worst.is_finite() {
            lambda = lambda.max(worst);
        }
    }
    if lambda == NEG_INF {
        None
    } else {
        Some(lambda)
    }
}

/// The critical circuit (the cycle achieving the maximum cycle mean) by
/// brute-force simple-cycle enumeration — for the small timing matrices this
/// instrument runs on, and the ground truth [`max_cycle_mean`] is checked
/// against.
#[must_use]
pub fn critical_circuit(matrix: &[Vec<f64>]) -> Option<Vec<usize>> {
    brute_force_max_cycle(matrix).map(|(cycle, _)| cycle)
}

/// Brute-force maximum-mean simple cycle: `(cycle nodes, mean)`.
#[must_use]
pub fn brute_force_max_cycle(matrix: &[Vec<f64>]) -> Option<(Vec<usize>, f64)> {
    let n = matrix.len();
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
        );
        on_stack[start] = false;
    }
    best
}

fn dfs_cycles(
    matrix: &[Vec<f64>],
    start: usize,
    u: usize,
    weight: f64,
    stack: &mut Vec<usize>,
    on_stack: &mut [bool],
    best: &mut Option<(Vec<usize>, f64)>,
) {
    let n = matrix.len();
    for v in 0..n {
        let w = matrix[u][v];
        if w == NEG_INF {
            continue;
        }
        if v == start {
            // closed a cycle back to the start.
            let len = stack.len();
            let mean = (weight + w) / len as f64;
            let is_better = best.as_ref().is_none_or(|(_, m)| mean > *m + 1e-12);
            if is_better {
                *best = Some((stack.clone(), mean));
            }
        } else if v > start && !on_stack[v] {
            on_stack[v] = true;
            stack.push(v);
            dfs_cycles(matrix, start, v, weight + w, stack, on_stack, best);
            stack.pop();
            on_stack[v] = false;
        }
    }
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
    pub fn record(&mut self, entry: AuditEntry) {
        self.entries.push(entry);
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
    #[must_use]
    pub fn promoted(&self, min_samples: usize, min_hit_rate: f64) -> bool {
        self.entries.len() >= min_samples && self.hit_rate() >= min_hit_rate
    }
}
