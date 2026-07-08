# CONTRACT: fs-report

Automatic lab notebooks + semantic design diffs: reproducibility as a side
effect.

## Purpose and layer

Layer L6 (HELM). No dependencies — pure Rust (in-house FNV-1a for content
addressing).

## Public types and semantics

- `Quantity { value, unit }` — a dimensioned value (units on every value).
- `ReproStep { op, args }` — one replayable operation of the reproducibility IR.
- `Block` — `Prose` / `Metric { name, quantity }` / `Step(ReproStep)`.
- `LabNotebook { title, seed, version, blocks }` — builder (`prose`, `metric`,
  `step`); `metrics`, `repro_ir` (the exact reproducing IR), `render_markdown`
  (deterministic), `content_hash` (FNV-1a of the render — content-addressed).
- `FeatureDelta { name, before, after, abs_change, rel_change, unit }` +
  `describe()`; `semantic_diff(before, after)` — a per-feature geometric
  attribution ranked by significance.

## Invariants

- The render is DETERMINISTIC, so `content_hash` is stable across runs and
  changes whenever any content changes (no silent drift) — replaying the same
  study reproduces the same artifact hash (the reproducibility loop closes).
- Every metric renders with its UNIT (P10 extends to reports).
- `repro_ir` returns the study's steps in order (the exact reproducing recipe).
- `semantic_diff` recovers per-feature absolute + relative edits and ranks them
  by `|relative change|` (largest first), with the feature name as tiebreak.

## Error model

Total functions; no panics.

## Determinism class

Fully deterministic: rendering, hashing, and diffing are pure functions of the
notebook / designs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/report.rs` (6 cases): the notebook renders all sections with units;
metrics carry their units; the notebook carries the exact reproducing IR; the
reproducibility loop closes by content hash (stable + change-sensitive);
semantic diff recovers known edits ranked by significance; determinism.

## No-claim boundaries

- v0 is the notebook DATA MODEL + deterministic Markdown render + content hash +
  the reproducing IR, and a scalar-feature semantic diff. The fuller deliverable
  — FrankenPandas frames over the Design Ledger, HTML with embedded data tables
  and LUMEN renders, convergence tables, Error/Time-Ledger attributions, and
  report generation being itself a ledgered op — is staged.
- `semantic_diff` compares scalar feature maps; the GEOMETRIC diff proper
  (varifold / optimal-transport distance with transport-plan visualization and
  per-region attribution across chart types) is the fuller deliverable.
- The ≤100 ms latency-lane serving guarantee is an fs-exec two-lane integration,
  measured there — out of scope for this pure crate.
