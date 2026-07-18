# CONTRACT: fs-grammar-e2e

GrammarForge — certified-fabricable geometric program discovery. Layer L4
(ASCENT).

## Purpose and layer

Composes `fs-shapeprog` (CSG programs + certified rewrites), `fs-archive`
(MAP-Elites), `fs-fab` (manufacturability), `fs-evidence` (Verified). Deps
point downward.

## Public types and semantics

- `target() -> Geom` — the peanut target (two unit spheres at `x=±0.8`).
- `build_program(r1, r2, d, o) -> Geom` — a candidate CSG program.
- `assess_simplification(program, radius_threshold, samples) ->
  SimplificationAssessment` — binds the exact local threshold, before/after
  sizes, optional finite certificate, exact bits of any rejected invalid
  certificate, optional finite outward sampled check, typed
  `SimplificationCheckStatus`, and exact underlying `SimplifyRefusal`.
- `SimplificationCheckStatus` — stable typed outcomes and wire codes:
  `Certified=0`, `StructuralEmptyAgreement=1`, `SimplifierRefused=2`,
  `NonFiniteCertificate=3`, `NegativeCertificate=4`,
  `DiscrepancyEvidenceRefused=5`, `CertificateCheckExceeded=6`, and
  `ThresholdMismatch=7`.
- `SimplificationSummary` — shared native/WASM aggregate accounting. It records
  the local radius threshold separately from the maximum compositional
  certificate and maximum admitted conservative outward finite-sample check,
  plus typed counts. Both the individual assessment and aggregate summary are
  sealed construction-only evidence types with read-only accessors; callers
  cannot manufacture or mutate a `Certified` state by setting public fields.
  Equality follows exact IEEE-754 provenance bits, so NaN-payload refusal replay
  and signed-zero threshold identity agree with aggregation semantics.
- `run_campaign(match_tol, simplify_radius_threshold) -> GrammarReport` —
  sweeps a program grid, illuminates (size × fab-margin), simplifies each elite,
  and re-verifies every simplification through the shared summary. Historical
  scalar report fields (`simplified_count`, sizes, `max_certified_error`, and
  `simplification_sound`) remain snapshot mirrors for compatibility; the typed
  `simplification` summary is authoritative.

## Invariants

- ILLUMINATION: a MAP-Elites archive of the best-matching program per
  (size × fab-margin) niche.
- CERTIFICATE-PRESERVING SIMPLIFICATION: some elites shrink; for EVERY admitted
  elite the independently re-measured outward SDF discrepancy between the
  original and simplified program is no larger than its finite compositional
  certificate. No comparison epsilon can turn a conservative check exceedance
  into an admitted result.
- THRESHOLD IS NOT ERROR: `simplify_radius_threshold` is the strict local
  `|offset radius|` admission threshold. It is not a global error budget. An
  admitted `0.02` offset at threshold `0.03` has the exact context-free
  certificate `0.04`; sequential effects and retained rounded parents can make
  the aggregate envelope larger still.
- FAIL-CLOSED STATUS: simplifier refusal, non-finite or finite-negative
  certificate, discrepancy-evidence refusal, conservative sample-check
  exceedance, or threshold mismatch makes `SimplificationSummary::is_sound()`
  false and prevents a `Verified` headline. Structural-empty `+∞/+∞` agreement
  is typed and counted separately rather than confused with invalid non-finite
  evidence.
- COMPLETENESS: campaign soundness additionally requires exactly one assessment
  per elite through `is_complete_and_sound(expected)`. A certified strict subset
  cannot authorize an “all elites checked” headline.
- The headline is `Verified` iff the best program matches within `match_tol`, is
  fab-satisfied, and the full simplification summary is sound.
- Deterministic (fixed program grid + sample grid; no RNG).

## Error model

For representable in-memory programs whose recursive ShapeProg operations
complete, `assess_simplification` propagates structured core refusal and maps
invalid certificate/evidence states to typed status rather than treating
transactional rollback (`program == original`, `max_error == 0`) as proof.
`run_campaign` expects a non-empty archive (the fixed grid guarantees it).

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

The crate-local unit suite injects otherwise unreachable defensive certificate
states to verify sealed aggregation, counter accounting, severity, and
fail-closed completeness. `tests/grammar.rs` provides G0/G3 coverage for a
fabricable illuminated family;
the exact `0.02/0.03 -> 0.04` threshold/envelope distinction; zero/drop/keep and
sequential cases; translation, sign, and constructor-nesting covariance; the
mandatory `0.006/0.01` sequential-drop witness; typed refusal, missing evidence,
and structural-empty agreement; sampled-witness containment; target identity;
and deterministic summary replay. Assertions and JSON diagnostics report both
the local threshold and returned certificate.

## No-claim boundaries

- The grammar is a fixed two-primitive CSG family; independent discrepancy is
  sampled on a finite grid and rounded outward. It is a conservative executable
  admission check, not a continuum proof or a mathematically downward-rounded
  lower bound. Global authority comes only from the compositional ShapeProg
  certificate under its admitted finite-evaluation hypotheses.
- `radius_threshold` controls local rewrite admission only. No consumer may
  infer `max_certified_error <= radius_threshold` or silently replace the sound
  `2*|radius|` envelope with the historical real-arithmetic `|radius|` claim.
- `StructuralEmptyAgreement` is authorized only when the fail-closed core
  discrepancy checker admits the pair and every supplied direct value is the
  structural `+∞` sentinel; arbitrary non-finite evidence remains refused.
- Fabricability uses a single minimum-feature proxy. A full generative grammar
  plus adjoint fitness is the fuller `fs-shapeprog`/`fs-xform` deliverable.
