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
  headline with the effective requirement, complete requirement and
  safety-factor document/version/locator lineage,
  exact evidence/context/replay identities, flip conditions, paired
  attribution headlines, and the complete indented audit projection.
- `project_decision_gate_markdown(authority, compliance)` — deterministic
  project-intent walkthrough retaining the context, intended decision, gate,
  consequence, both full source lineages, content identities, lower-layer
  tri-state verdict, and explicit `admitted` versus
  `refused: this context requires a determinate assessment` outcome.
- `no_useful_bound_markdown(refusal)` — a distinct `NoUsefulBound` visual
  class showing achieved enclosure/width, failed caller threshold, unit,
  decision context, closed cause, exact E09 suggested reformulation, and an
  explicit no-certificate/no-compliance boundary. Decision headlines invoke
  it automatically for typed useful-bound indeterminacy.
- `regime_no_claims_markdown(audit)` — deterministic final operating-envelope
  no-claim section. It omits fully in-domain QoIs and renders every demoted
  receipt with the human diagnosis, strong receipt identity, and exact canonical
  JSON. An override remains visible but cannot restore color.
- `project_regime_audit_outputs(package, audit)` — atomic product projection
  returning the reviewer-facing no-claim section and the evidence package with
  every demotion retained from that same immutable audit. It refuses as a unit;
  callers cannot receive a partial projection.
- `AsBuiltQoiDelta::try_new(qoi, unit, nominal, as_built, geometry)` — one
  nominal-versus-as-built comparison row. It borrows the geometry term rather
  than copying it, and refuses at construction when the cited term is not the
  budget's geometry slot (`AsBuiltDeltaError::NotGeometryTerm`) or a value is
  not finite (`NonFiniteValue`). `shift()` is `as_built - nominal`;
  `geometry_half_width()` reads the conservative half-width, taking the UPPER
  endpoint of an `IntervalBound` because that variant brackets the half-width
  itself; `shift_resolved()` is a TRI-STATE — `Some(true)` when the shift
  strictly exceeds that half-width, `Some(false)` when it does not, and `None`
  when the term publishes no half-width at all, so an unanswerable question is
  never rendered as a negative answer.
- `nominal_vs_as_built_markdown(rows)` — deterministic design-build-measure-
  update table rendering each QoI's shift, the attributed geometry term with its
  `role@digest`, the resolution tri-state, and a correlation line derived from
  the rows themselves: one shared propagation record across several QoIs renders
  as a single correlated block, distinct records render as an explicit refusal
  of that claim, and one row claims no cross-QoI correlation at all. Row order
  is caller order, not a sort.

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
- Project gate rendering never changes the supplied verdict. Only
  non-safety-critical scoping displays indeterminate as admitted;
  design/sign-off and safety-critical contexts display it as refused.
- A `NoUsefulBound` headline never uses compliant, non-compliant, verified, or
  certificate status. Its cause and reformulation remain visible.
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
headlines, units and full authority lineages, explicit flip actions, the
projection-only no-claim boundary, and the same indeterminate physics rendered
as admitted for advisory scoping versus refused for compliance sign-off.
`tests/useful_bound.rs` covers the distinct refusal visual class, cause,
suggested E09 reformulation, and explicit no-certificate boundary.

`tests/as_built.rs` covers the nominal-versus-as-built projection: fail-closed
row construction (non-geometry term, non-finite value), the resolution
tri-state including the strict boundary case and the `Unknown` term a rejected
upstream linearization produces, the `IntervalBound` upper-endpoint reading,
the shared-record correlated-block sentence, the explicit refusal when terms
cite distinct records, single-row and empty tables making no correlation claim,
render determinism, and caller-order row rendering.

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
- The project-gate block presents declared project intent beside a supplied
  lower-layer verdict. It does not itself enforce the gate (the
  `fs-project::ProjectDecisionAuthority` adapter does), and its scoping
  `admitted` label does not authenticate requirement sources or authorize
  design release/compliance sign-off.
- The `NoUsefulBound` block projects a lower-layer refusal. It does not prove
  that the achieved enclosure is rigorous, select the threshold, or establish
  that the suggested reformulation will yield useful evidence.
- The nominal-versus-as-built table is a PRESENTATION of two already-computed
  solves and one already-computed geometry term. It runs no solve, performs no
  covariance propagation, and does not verify that the two values came from the
  same study, the same mesh, the same physics, or the placements they claim.
  The rendered shift is a difference between two computed numbers, never a
  measurement of the physical part. The geometry term bounds only the
  pose-propagated contribution and covers no discretization, parameter,
  boundary-condition, or model-form error, so `shift_resolved() == Some(true)`
  means the shift exceeds the POSE band alone and is not a claim that the shift
  is physically real. `Some(false)` is likewise not evidence that nothing
  changed. The correlation line reports only whether the supplied terms cite one
  propagation record; agreement on an identity is not proof that the upstream
  block is correctly calibrated, and this crate never resolves that identity.
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
