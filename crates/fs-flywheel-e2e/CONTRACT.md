# fs-flywheel-e2e — CONTRACT

FLYWHEEL CLOSES (bead lmp4.18): the whole-loop e2e harness testing the
addendum's central claim — that speculation (9), incremental recompute
(2), the sheaf-adjudicated merge (10), and tombstones (E) COMPOUND.

## Purpose and layer
Layer L6 (integration harness over the flywheel crates). [F]-gated per
the Ambition-Tag rule.

## Public types and semantics
- `LoopConfig`: per-proposal toggles + the G4 cancel point; `baseline()`
  and `composed()` are the measurement anchors.
- `run_loop(config, iterations, seed) -> LoopReport`: the modeled
  two-agent design-iteration workload on the corpus's CHT-wedge edit
  trace — total modeled cost, per-stage event trace, accept rate,
  skips, merge verdicts, tombstone blocks, the end-to-end headline
  color, and `trace_hash()` for G5.
- `speedups(iterations, seed)`: isolated per-proposal speedups + the
  composed loop (baseline cost ratios).

## Invariants
- Colors: accepted speculation is Estimated; the headline composes by
  the weakest-input rule; downstream upgrades are refused by the
  ColorGraph write gate (laundering-across-the-loop).
- The compounding measurement follows the review-round-3 protocol:
  isolated AND composed over 5 seeded replays, composed >
  max(isolated) by the stated 1.15x margin, coefficient of variation
  reported and bounded.
- Member-crate refusals (merge conflicts, tombstone blocks, skip
  misses) are OUTCOMES the report records, not errors.

## Error model
The harness is total over its config space; member-crate errors would
be bugs (panics), not runtime conditions.

## Determinism class
Deterministic in `seed`: one LCG stream, fixed stage order,
fixed-order accumulation — bit-equal costs and trace hashes on replay
(G5, tested). Member crates are deterministic per their contracts.

## Cancellation behavior
`cancel_after_stages` models the G4 storm at every stage boundary; the
loop unwinds BEFORE any state mutation for the interrupted stage, so a
partial trace is a clean prefix (tested at five cancel points).

## Unsafe boundary
No `unsafe` anywhere in this crate.

## Feature flags
`flywheel-e2e` ([F], default OFF): gates the entire harness until its
Gauntlet tier + kill metric are green.

## Conformance tests
tests/e2e.rs — fw-001 compounding (margin + CV over 5 replays), fw-002
laundering-across-the-loop, fw-003 G5 whole-loop determinism, fw-004
G4 cancellation storm, fw-005 telemetry completeness; tests/dbg.rs —
config-sweep smoke. Battery ~40 s.

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
