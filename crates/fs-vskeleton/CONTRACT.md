# CONTRACT: fs-vskeleton

## Purpose and layer
The PV vertical skeleton (patch Rev R): the tiny end-to-end slice proving the
typed-value semantics — SDF → PDE → objective → adjoint check → optimize →
ledger replay → deterministic rerun → report. Layer: L6 (orchestrates; may
depend on everything). This crate is a PROVING ARTIFACT: real crates (fs-ir,
fs-ledger, fs-exec, fs-geom, fs-opt) supersede its minis; its e2e suite
remains as the continuum's smallest regression test.

## Public types and semantics
`run_study(study_text, db_path) -> StudyOutcome` (objective/radius/grad-check
traces, budget spent, report, artifact hashes); `replay(db_path)` (integrity
scan + re-execute + hash compare); `sexpr` (minimal total s-expr reader);
`model` (StudySpec parse w/ mandatory seed+budget, EdgeLaw one-source-of-truth
stencil, CG w/ cancellation polls, adjoint + central-difference gradients);
`ledger::MiniLedger` (fsqlite ops/artifacts/edges, domain-separated BLAKE3 content addressing, format-versioned: pre-v2 FNV ledgers are version-refused).

## Invariants
- Bitwise deterministic: same study → identical artifact hashes across runs
  (fixed-chunk parallel maps, fixed-order reductions).
- Gradient truth: every optimizer step gates on adjoint-vs-central-difference
  rel err < 1e-4 or the study aborts (plan §8.7 in miniature).
- Budgets are enforced (BudgetExhausted), never advisory (P4).
- Cancellation is request → drain → finalize; ledger never holds torn state.
- Replay refuses tampered ledgers (byte-corruption detection).

## Error model
All errors are teaching strings naming the fix (BudgetExhausted,
GradientCheckFailed, LedgerCorruption, SolverStalled, parse errors with
positions). No panics on any study input (parser garbage-battery-tested).

## Determinism class
Deterministic (single ISA): bit-stable across runs and thread schedules by
construction. Cross-ISA claims deferred to fs-math/G5.

## Cancellation behavior
Cooperative AtomicBool polls at row/iteration granularity; asupersync-scope
integration is fs-exec's bead (Budget vocabulary already smoke-tested there).

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
tests/e2e.rs pv-001..pv-005 (determinism, replay, corruption, optimization +
gradient gates, budget teaching errors) + 7 model/parser unit tests including
the Poisson series-reference check (peak u ≈ 0.0736713 for -Δu=1).

## No-claim boundaries
Performance (unoptimized by design); no RNG
consumption (seed recorded for provenance only); 2D scalar physics only.
