# CONTRACT: fs-tropical

Tropical (max-plus) critical-path analytics: task-DAG timing in the max-plus
semiring as an online scheduling instrument.

## Purpose and layer

Layer L1 (fundamental graph/algebra). No dependencies — pure Rust.

## Public types and semantics

- `oplus(a,b) = max`, `otimes(a,b) = +` (absorbing at `NEG_INF`) — the max-plus
  semiring.
- `TaskDag::new(latencies).with_edge(u,v)` — a precedence DAG.
  `critical_path() -> Result<CriticalPath, TropicalError>` gives the makespan
  (CPM longest path), an outward-rounded makespan enclosure, a deterministic
  witness path, certified-uniqueness flag, and nominal per-task slack; finite
  non-negative latencies, task/edge caps, path sums, and endpoints are checked.
- `tune_next() -> Result<_>` — positive-latency tasks on the certified unique
  critical path, ranked by latency; an ambiguous/tied path returns no
  single-target advice. `bottleneck()` returns a task only when the largest
  positive latency is itself unique.
- `max_cycle_mean(&[Vec<f64>]) -> Result<Option<f64>, TropicalError>` — the
  tropical eigenvalue (Karp), the pipelined steady-state throughput; `None` if
  acyclic. Dense dimensions are capped and matrices must be square with finite
  weights or exact `NEG_INF` no-edge sentinels.
- `critical_circuit` / `brute_force_max_cycle` — the maximum-mean cycle (nodes,
  mean), the ground truth `max_cycle_mean` is checked against.
- `RecommendationAudit` — bounded, validated `record`, `hit_rate`, and fallible
  `promoted(min_samples, min_hit_rate)`: the [M] promotion gate on outcomes.

## Invariants

- `critical_path` returns the longest source-to-sink path length as the
  nominal makespan and a directed-rounding enclosure containing the real-valued
  max-plus result. A zero-latency terminal remains in the witness; an empty DAG
  has the zero identity makespan but no path and therefore no uniqueness claim.
- `max_cycle_mean` (Karp) equals the brute-force maximum simple-cycle mean.
- `bottleneck` is the strictly highest positive-latency task on a certified
  unique path. A tied critical set, tied maximum task latency, or all-zero path
  yields no single-task claim because the named target would not be unique or
  actionable.
- The audit only reports `promoted` with at least one required sample, enough
  retained observations, and a finite threshold in `0..=1`.

## Error model

Structured `TropicalError` for invalid domains/audit fields, bounded-work
refusals, cycles, bad edges, and finite-arithmetic overflow; no panics.

## Determinism class

Fully deterministic: every quantity is a pure function of the DAG / matrix.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/tropical.rs`: the max-plus semiring; directed makespan containment;
source-to-terminal witnesses including zero-latency tails; unique-path
positive-work tune-next and path/task-tie refusal; cyclic rejection; Karp
equality with exact-order brute force; bounded validated promotion audit; and
malformed, ragged, non-finite, overflowing, and oversized inputs fail closed.

## No-claim boundaries

- `critical_circuit` / `brute_force_max_cycle` enumerate SIMPLE cycles — for the
  small timing matrices this instrument runs on and therefore refuse more than
  `MAX_BRUTE_FORCE_NODES`; `max_cycle_mean` (Karp) is the
  efficient `O(n·E)` path for large graphs.
- The LIVE pipeline (windowed re-computation over the ledger's event stream with
  drift detection, and contributing the tropical kernels upstream to
  FrankenNetworkx) consumes this core over `fs-ledger` events — this crate
  provides the batch analytics + the audit.
- The [M] discipline is advisory-until-validated; `RecommendationAudit` is the
  bounded self-audit machinery, while policy identity/authentication and
  durable cross-process promotion remain the planner/ledger's responsibility.
