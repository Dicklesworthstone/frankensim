# fs-flywheel-e2e — CONTRACT

## 1. Purpose
FLYWHEEL CLOSES (bead lmp4.18): the whole-loop e2e harness testing the
addendum's central claim — that speculation (9), incremental recompute
(2), the sheaf-adjudicated merge (10), and tombstones (E) COMPOUND.
[F]-gated behind the `flywheel-e2e` feature until its Gauntlet tier +
kill metric are green.

## 2. Public surface
`LoopConfig` (per-proposal toggles + the G4 cancel point),
`run_loop(config, iterations, seed) -> LoopReport` (total modeled cost,
per-stage event trace, accept rate, skips, merge verdicts, tombstone
blocks, the end-to-end headline color, `trace_hash()` for G5),
`speedups(iterations, seed)` (isolated per-proposal + composed).

## 3. Guarantees
- Deterministic in `seed`: bit-equal costs and trace hashes on replay
  (G5, tested).
- A mid-loop cancel yields a CLEAN PREFIX trace: no torn stages, no
  cross-run residue (G4, five cancel points tested).
- Colors: accepted speculation results are Estimated and the headline
  composes by the weakest-input rule; downstream upgrades are refused
  by the ColorGraph write gate (laundering-across-the-loop, tested).
- The compounding measurement follows the review-round-3 protocol:
  isolated speedups AND composed over 5 seeded replays, composed >
  max(isolated) by the stated 1.15x margin, coefficient of variation
  reported and bounded.

## 4. Errors and refusals
The harness is total over its config space; member-crate refusals
(merge conflicts, tombstone blocks, skip misses) are OUTCOMES the
report records, not errors.

## 5. Determinism
One LCG stream, fixed stage order, fixed-order accumulation; the
member crates (store hashing, merge Gauss–Seidel, tombstone pi-groups)
are deterministic per their own contracts.

## 6. Cancellation
`cancel_after_stages` models the G4 storm at every stage boundary; the
loop unwinds before any state mutation for the interrupted stage.

## 7. Performance envelope
Modeled cost units; the battery runs in ~40 s (12-iteration loops, 36
full runs). Not a wall-clock benchmark.

## 8. Dependencies
fs-benchmark (the corpus), fs-spececo (9), fs-recompute (2), fs-geom
`sheaf-merge` (10), fs-ledger tombstones (E), fs-evidence (colors),
fs-qty (dimensions).

## Feature flags
`flywheel-e2e` ([F], default OFF): gates the entire harness per the
Ambition-Tag rule until its Gauntlet tier + kill metric are green.

## Conformance tests
tests/e2e.rs — fw-001 compounding (margin + CV over 5 replays), fw-002
laundering-across-the-loop, fw-003 G5 whole-loop determinism, fw-004 G4
cancellation storm, fw-005 telemetry completeness; tests/dbg.rs —
config-sweep smoke.

## Unsafe boundary
No `unsafe` anywhere in this crate.

## No-claim boundaries
- COSTS ARE MODELED UNITS from the corpus's op counts: the loop
  MECHANICS are measured; wall-clock physics compounding lands when
  the wedge's real solvers (CHT vertical) replace the cost model.
- The Proposal-8 query stage is a soft edge (Phase 2 per the polish
  note); the headline-color tail models its admission check.
- The two-agent concurrency model is synchronous round-based; a live
  multi-process swarm trial is the xpck.3 milestone's territory.
- The G4 storm is a deterministic stage-boundary cancel model, not the
  base plan's thread-storm harness (fs-exec owns that).
