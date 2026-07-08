# CONTRACT: fs-benchmark

The wedge-vertical benchmark & trace corpus (plan addendum, Proposal 7): the
single shared, versioned, deterministic artifact many kill criteria measure
against.

## Purpose and layer

Layer UTIL (versioned data + measurement helpers). Depends only on `fs-evidence`
(the `ColorRank` for reference-answer colors). Governance Rule 2 says an
un-instrumented kill measurement counts as killed; this corpus instruments
six+ proposals so they are not killed-by-default.

## Public types and semantics

- Datasets (all `const`): `query_set()` (conjugate-heat-transfer `QueryCase`s
  with `reference_answer`, `reference_cost`, `reference_color`),
  `design_tasks()` (`DesignTask`s with known optima), `edit_traces()`
  (`EditTrace`s with known-correct skip sets), `mms_battery()` (`MmsCase`
  elliptic references), `merge_trials()` (`MergeTrial`s with harmonic-conflict
  counts).
- Measurement helpers: `speedup(baseline, candidate)` (Prop 8 / 2, `>= 2×`),
  `win_rate(&[bool])` (Prop 1, `>= 0.70`), `rate(count, total)` /
  `conflict_rate(trial)` (Prop 10, `< 0.25`), `accept_rate(accepts, attempts)`
  (Prop 9). All guard divide-by-zero.
- `instrumented_proposals()` — the `InstrumentedProposal`s (proposal, dataset,
  kill metric) this corpus discharges.
- `corpus_digest() -> u64` — an FNV-1a content digest over the whole corpus
  (bit-stable → replayable).
- `audit() -> CorpusAudit` — every dataset non-empty and every instrumented
  proposal references a real dataset.

## Invariants

- DETERMINISM: `corpus_digest` is bit-stable across runs (const data) — the
  replayability the acceptance criteria demand.
- Every reference answer carries a color class; all three classes appear.
- GOVERNANCE RULE 2: every instrumented proposal references a real, non-empty
  dataset and declares a kill metric — no proposal is killed-by-default for
  lack of an instrument.
- Measurement helpers never divide by zero (return `0.0`).

## Error model

Total functions; no panics. Completeness gaps surface as `CorpusAudit::gaps`.

## Determinism class

Fully deterministic: all data is `const`; helpers are pure functions.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/benchmark.rs` (Proposal 7, 8 cases): every dataset populated; reference
answers carry their color (all three classes); the measurement helpers compute
the kill numbers (speedup / win-rate / accept-rate + divide-by-zero guards);
the harmonic-conflict rate from merge trials; the skippable fraction from edit
traces; Governance-Rule-2 discharge for all eight instrumented proposals; the
deterministic digest; the complete audit.

## No-claim boundaries

- The datasets are SMALL, representative fixtures encoding the corpus SHAPE and
  the measurement API; populating them with the full high-precision reference
  solves (and the real recorded traces) is the vertical-kernel work that
  consumes this contract.
- `reference_color` records the color CLASS (`ColorRank`); the full `Color`
  certificate payload is materialized by the reference solver.
- The corpus provides the DATA + the measurement helpers; each proposal's bead
  computes its own kill number by feeding its results through them.
- Coupling to the base-plan Gauntlet G1/G2 registries + fs-roofline for the
  cost model is a downstream integration.
