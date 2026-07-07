# CONTRACT: fs-spececo

Certified-speculation accept/reject economics (plan addendum, Proposal 9): the
decision-and-learning layer over the verifier and the proposer zoo.

## Purpose and layer

Layer L6 (orchestration/decision). No numerical dependencies — pure decision +
telemetry + drift logic; the verifier (fs-verify) and proposers produce the
bounds and candidates this reasons about.

## Public types and semantics

- `decide(certified_bound, query_tolerance) -> Decision` — `AcceptOutright`
  iff `certified_bound` is finite, non-negative, and `<= query_tolerance`;
  otherwise `WarmStart`. FAIL-SAFE: a non-finite or negative bound (or a
  non-finite tolerance) never accepts, so a bad proposer can only trigger a
  warm-start, never a false accept.
- `SolveRecord { proposer_id, regime, accepted, bound, iterations_saved }` —
  the Error-Ledger solve-node telemetry fields (`iterations_saved` may be
  negative: a warm-start worse than a cold start).
- `ProposerTelemetry` — folds records per `(proposer, regime)` into `Stats
  { attempts, accepts, positive_saves, net_iterations_saved }`. `accept_rate`
  and `mean_iterations_saved` return 0 (conservative) when there is no
  telemetry — never a divide-by-zero. A negative warm-start is never counted a
  "positive save" but does lower the net total.
- `DriftDetector::new(min_samples, demote_below, restore_above)` +
  `update(&telem, proposer, regime) -> bool` / `is_demoted(...)` — demotes a
  proposer in a regime when its accept-rate collapses to/under `demote_below`,
  restores it only when the rate clears `restore_above`, with a minimum sample
  count gating both. `new` panics if `restore_above < demote_below`.

## Invariants

- `decide` never accepts on a non-finite/negative bound or non-finite
  tolerance.
- Drift state does not change below `min_samples` (noise-gated), and cannot
  flap: `demote_below < restore_above` (a strict hysteresis band) means a rate
  inside the band holds the previous state.
- Telemetry accessors are total (0 on missing data).

## Error model

Total functions; the only failure is a programmer error (`DriftDetector::new`
with an inverted hysteresis band), which panics.

## Determinism class

Fully deterministic: telemetry folding, decisions, and drift transitions are
pure functions of the recorded stream (no RNG, no I/O); replay reproduces every
decision and demotion.

## Cancellation behavior

None (synchronous pure functions; the true solve a warm-start triggers runs
under the caller's `Cx`).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/economics.rs` (Proposal 9, 8 cases): decision accepts only a
finite/non-negative bound within tolerance (equal accepted); fail-safe on
NaN/inf/negative bounds and non-finite tolerance; telemetry accumulates
accepts + savings with a negative warm-start counted as not-a-win but net-
lowering; unknown-pair rates are 0 (no panic); drift does not demote below the
minimum sample count; drift demotes on collapse and HOLDS through a partial
recovery in the hysteresis band, then restores above the upper threshold;
`new` rejects an inverted band; determinism.

## No-claim boundaries

- This crate owns the DECISION + telemetry accumulation + drift logic. The
  telemetry SCHEMA persistence (writing `SolveRecord` fields onto ledger solve
  nodes) is fs-ledger's; a consumer wires the two.
- `iterations_saved` is supplied by the caller (measured against a cold solve);
  this crate accumulates but does not itself run solvers.
- The speculative RACE (running proposer + target concurrently under one `Cx`,
  cancelling the loser) is the executor's concern; this crate decides accept
  vs warm-start given a bound, it does not schedule the race.
- Drift thresholds are policy inputs; learning them from fleet data is a later
  refinement.
