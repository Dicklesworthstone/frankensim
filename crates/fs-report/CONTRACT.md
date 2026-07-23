# CONTRACT: fs-report

Automatic lab notebooks + semantic design diffs: reproducibility as a side
effect.

## Purpose and layer

Layer L6 (HELM). Pure deterministic projection over L6 session decisions and
lower-layer evidence; notebook identities use `fs-obs::IdentityBuilder`.

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
- `decision_headline_markdown(assessment)` — compact tri-state decision
  headline with the effective requirement, declared safety factor and policy,
  exact evidence/context/replay identities, flip conditions, paired
  attribution headlines, and the complete indented audit projection.
- `regime_no_claims_markdown(audit)` — deterministic final operating-envelope
  no-claim section. It omits fully in-domain QoIs and renders every demoted
  receipt with the human diagnosis, strong receipt identity, and exact canonical
  JSON. An override remains visible but cannot restore color.
- `project_regime_audit_outputs(package, audit)` — atomic product projection
  returning the reviewer-facing no-claim section and the evidence package with
  every demotion retained from that same immutable audit. It refuses as a unit;
  callers cannot receive a partial projection.

## Invariants

- The render is DETERMINISTIC, so `content_hash` is stable across runs and
  changes whenever any content changes (no silent drift) — replaying the same
  study reproduces the same artifact hash (the reproducibility loop closes).
- Every metric renders with its UNIT (P10 extends to reports).
- `repro_ir` returns the study's steps in order (the exact reproducing recipe).
- `semantic_diff` recovers per-feature absolute + relative edits and ranks them
  by `|relative change|` (largest first), with the feature name as tiebreak.
- A decision headline never maps `indeterminate` to either binary verdict and
  always retains the exact `DecisionAssessment` and replay-package identities.
- Safety factors are reported as already applied to the effective sourced
  limit; `fs-report` does not invent a multiplication or division convention.
- A regime no-claim section is absent only when every final receipt is fully
  in-domain. Rendering is sorted by QoI and receipt identity, includes exact
  split-sweep partitions, and never averages coverage or drops a demotion.
- The coupled regime projection derives both artifacts from one audit. A
  demoted receipt therefore cannot appear in only the human or only the machine
  projection returned by that API; an in-domain audit yields no section and an
  unchanged package.

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

`tests/report.rs`: the notebook renders all sections with units;
metrics carry their units; the notebook carries the exact reproducing IR; the
reproducibility loop closes by content hash (stable + change-sensitive);
semantic diff recovers known edits ranked by significance; determinism; final
regime demotions render with exact receipts while fully in-domain audits emit no
no-claim section; the coupled projection proves every demoted receipt appears
in both artifacts and that in-domain projection is an exact no-op.

`tests/decision.rs` covers deterministic indeterminate and binary decision
headlines, units and authority identities, explicit flip actions, and the
projection-only no-claim boundary.

## No-claim boundaries

- v0 is the notebook DATA MODEL + deterministic Markdown render + content hash +
  the reproducing IR, and a scalar-feature semantic diff. The fuller deliverable
  — FrankenPandas frames over the Design Ledger, HTML with embedded data tables
  and LUMEN renders, convergence tables, Error/Time-Ledger attributions, and
  report generation being itself a ledgered op — is staged.
- The decision headline is a presentation of an already-admitted
  `fs-session::DecisionAssessment`; it does not resolve content hashes,
  authenticate requirement or policy sources, recompute compliance or
  attribution, price evidence actions, or certify scientific evidence.
- The regime no-claim section is a presentation of an `fs-regime` audit. Its
  receipt identities bind exact bytes but do not authenticate card ownership,
  calibration sources, or the completeness of the orchestrator-supplied card
  set and operating envelope.
- `retain_regime_demotions_in_package` adds one deterministic declaration per
  demoted receipt and is an exact-retry-idempotent no-op for fully in-domain
  receipts. Each declaration is permanently `Estimated` with infinite
  dispersion and package-root-binds the collection provenance, receipt
  identity, and exact canonical JSON. It deliberately does not use
  `SemanticWitness`, which is reserved for `Verified` source certificates; the
  projection authenticates no card authority and cannot restore color.
- Reference-run assembly for every requirement-bearing QoI, deterministic HTML,
  report routing through the CLI, and `frankensim explain` remain staged behind
  their prerequisite work.
- `semantic_diff` compares scalar feature maps; the GEOMETRIC diff proper
  (varifold / optimal-transport distance with transport-plan visualization and
  per-region attribution across chart types) is the fuller deliverable.
- The ≤100 ms latency-lane serving guarantee is an fs-exec two-lane integration,
  measured there — out of scope for this pure crate.
