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
  elliptic references), `merge_trials()` (synthetic `MergeTrial` fixtures with
  candidate-remainder counts for exercising the corpus shape and rate API).
- Measurement helpers: `speedup`, `win_rate`, `rate`, `conflict_rate`, and
  `accept_rate` return typed failures for absent or invalid denominators instead
  of manufacturing a zero-valued measurement.
- `resolve_query_reference` applies the safe deny-all admission policy;
  `resolve_query_reference_with_verifier` can attach positive authority only
  after an injected verifier authenticates the exact retained query context.
- `instrumented_proposals()` declares evaluator schemas. A proposal counts as
  instrumented only when `evaluate_proposal` reconstructs every typed role from
  independently referenced retained evidence.
- `corpus_digest() -> ContentHash` is a schema-versioned, length-framed BLAKE3
  identity over the complete corpus semantics.
- `audit_corpus` applies the safe deny-all admission policy.
  `audit_corpus_with_verifier` allows a composition root to audit valid positive
  query evidence with its real admission capability. `audit()` remains the
  deny-all audit of the built-in corpus.

## Invariants

- DETERMINISM: `corpus_digest` is bit-stable across runs (const data) — the
  replayability the acceptance criteria demand.
- Positive query authority is impossible without exact retained-evidence
  resolution, a context-correct receipt, and an accepting injected verifier.
- GOVERNANCE RULE 2: declaring a non-empty dataset and kill metric is
  insufficient; all evaluator roles must resolve before instrumentation is
  available.
- Measurement helpers never divide by zero; they return a typed refusal.

## Error model

Evidence and metric failures are typed. Completeness and admission gaps surface
as `CorpusAudit::gaps`; the default audit deliberately refuses positive rank.

## Determinism class

Fully deterministic: all data is `const`; helpers are pure functions.

## Cancellation behavior

None.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/benchmark.rs`: retained-evidence resolution, positive-query admission,
default deny-all versus injected-authority audit behavior, typed denominator
receipts, proposal evaluator reconstruction/refusal, semantic-field mutation
sensitivity, deterministic identity, and fail-closed completeness auditing.

## No-claim boundaries

- The datasets are SMALL, representative fixtures encoding the corpus SHAPE and
  the measurement API; populating them with the full high-precision reference
  solves (and the real recorded traces) is the vertical-kernel work that
  consumes this contract.
- Built-in query references are Estimated declarations. Verified or Validated
  references require an admission verifier supplied by the composition root.
- The corpus provides the DATA + the measurement helpers; each proposal's bead
  computes its own kill number by feeding its results through them.
- Merge-trial counts are synthetic fixtures for the guarded candidate-remainder
  path; they are neither retained realistic trace evidence nor certified H¹ or
  topology counts. The full Proposal 10 gate must additionally count
  escalations, refusals, and type conflicts on retained trials.
- Coupling to the base-plan Gauntlet G1/G2 registries + fs-roofline for the
  cost model is a downstream integration.
