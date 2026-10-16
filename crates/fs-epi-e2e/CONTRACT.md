# CONTRACT: fs-epi-e2e

The epistemic type-system end-to-end suite (plan addendum, Proposal 3's Layer-2
conformance harness): the runnable battery that exercises the whole type system
and is the artifact of record that it FAILS SAFE.

## Purpose and layer

Layer L6. An integration crate: depends on `fs-evidence` (colors / laundering /
falsifier), `fs-opt` (the Goodhart guard), `fs-robust` (objective epistemics),
`fs-package` + `fs-checker` (evidence round-trip). It composes them; it owns no
new primitive.

## Public types and semantics

- `run_battery() -> EpiE2eReport` — runs all five stages.
- `EpiE2eReport` has sealed construction and private storage. `stages()` exposes
  an ordered read-only view; `complete()` requires the exact five canonical
  stage identities in order plus well-formed evidence; `passed()` is true only
  when that complete schema passes; `stage(name)` provides named lookup.
- `StageLog` has private fields and read-only `stage()`, `passed()`, and
  `events()` accessors. Its pass accessor also requires at least one nonblank
  evidence event. Events are returned as DATA, never printed, and callers
  cannot forge a passing log by struct construction.
- Per-stage entry points (`stage_laundering`, `stage_falsifier`,
  `stage_goodhart_guard`, `stage_objective_epistemics`,
  `stage_evidence_roundtrip`) for granular runs.

## The five stages (each a fail-closed assertion)

1. **Laundering** — `compose(verified, estimated)` yields estimated (min rank,
   no upgrade); a validated claim OUT of its regime auto-demotes to estimated,
   one IN its regime is preserved.
2. **Falsifier** — `catalog_gate` names a class with no falsifier declaration;
   the bounded diagnostic consequence×doubt allocator spends monotonically and
   zero claims → zero spend. This catalog lint does not admit a release claim.
3. **Goodhart guard** — a discretization-exploit endpoint is REFUSED (`Failed`)
   even when the other escalation steps pass; a genuine smooth optimum with the
   full escalation set is honored (`Cleared`); a guard missing steps stays
   `Provisional` (never false-cleared).
4. **Objective epistemics** — `robust_optimum` refuses an un-colored objective;
   the weakest input colors the headline; a colored, monotone fragility curve.
5. **Evidence round-trip** — a package re-verifies through the solver-free
   checker, renders its budget pie, and a tampered package fails with a
   localized `content-address-mismatch` finding.

## Invariants

- The full battery passes and is DETERMINISTIC (`run_battery() ==
  run_battery()`).
- Every stage's load-bearing fail-closed behavior holds; the guard clears only
  when a step ran for every escalation kind (so the suite registers a full set
  to demonstrate an honored optimum).

## Error model

No errors/panics; a stage records `passed = false` with its events on any
failure.

## Determinism class

Fully deterministic: all subsystems it drives are deterministic; no RNG, no I/O.

## Cancellation behavior

None (synchronous).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/e2e.rs` (Layer-2 conformance, 7 cases) plus internal schema-adversary
tests: the full battery passes only with the exact ordered five-stage schema and
nonblank evidence; empty, missing, duplicated, reordered, unexpected, or
blank-event reports fail closed; laundering fails closed; the falsifier-catalog lint names
an unpaired class; the guard refuses exploits but honors genuine optima and stays
provisional when it cannot check; objective epistemics holds the contract; the
evidence package round-trips and tamper is caught; the battery is deterministic.

## No-claim boundaries

- The falsifier stage checks declaration-catalog completeness and allocator
  arithmetic only. It does not prove a checker implementation is independent,
  bind a retained run to the exact claim instance, or authorize release.
- The suite emits its log events as returned DATA; wiring them onto the base
  plan's structured tracing / ledger event sinks is the harness integration.
- The guard's non-δ escalation steps (rung k+1, cross-representation,
  estimator-independence) are represented by trivially-passing stand-ins to
  demonstrate the CLEARED path; the real capability steps live in their
  subsystems.
- This is the Layer-2 (epistemic type system) e2e; the HELM/FLUX e2e suites are
  separate base-plan artifacts.
