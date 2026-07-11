# CONTRACT: fs-ir

> Status: ACTIVE (FrankenScript core, IR language v2). Owns the typed AST,
> both concrete syntaxes, study recognition, and verb lowering. Admission
> (dimensional/chart/budget/capability checks) is the gp3.5 bead's;
> the operator catalog is gp3.6's.

## Purpose and layer

FrankenScript — the system's one true interface (plan §11.1, Decalogue
P10): a typed, versioned IR with two isomorphic concrete syntaxes
(canonical s-expressions; lossless JSON mapping), both parsing to the same
typed AST. Layer: L6 (HELM). Runtime deps: `std` + fs-qty.

## Public types and semantics

- `Node`/`NodeKind`/`Span` — every node carries a byte span. Atoms are the
  real nouns: `Int`, `Float` (finite only), `Qty` (fs-qty SI value + dims +
  the ORIGINAL literal text — fs-qty normalizes 65deg → rad, so the
  literal is preserved verbatim for lossless printing; equality is
  value+dims), `Count` (information/core grants: B/KiB/MiB/GiB/cores —
  deliberately outside fs-qty's SI domain), `Seed` (0x… u64), `Str`,
  `Symbol`, `Keyword`, `List`.
- `CountValue` preserves bare integers as exact `u128` and decimal/exponent
  spellings as a bounded exact decimal (`u128` significand + base-10
  exponent). Whole-byte/core enforcement uses checked integer arithmetic;
  binary floating point is reporting-only. Mixed syntax classes remain
  distinct identities (`2B` differs from `2.0B`), while each class has one
  deterministic canonical spelling.
- `Node::same_shape` — semantic equality ignoring spans and Qty
  presentation; the isomorphism property is stated in terms of it.
- `sexpr::parse/print` — total reader with spans, comments (`;`),
  string escapes, deterministic atom classification (numeric-leading
  tokens MUST fully parse — a number with a garbage suffix is a structured
  error, never silently a symbol), depth cap (adversarial nesting refuses
  structurally). Printer output reparses to the same shape.
- `json::parse/print` — the lossless mapping (single-key tagged objects:
  i/f/q/c/seed/s/sym/kw; arrays = lists). Qty/Count/Seed reuse the s-expr
  literal grammar inside strings so ONE classifier owns numeric semantics
  for both syntaxes. Unknown tags and tag/literal mismatches refuse with
  spans.
- `Study::from_node` — recognizes Appendix C study forms: name, seed,
  versions/budget/capability clauses, `(let …)` bindings, body;
  `constellation_lock()` extracts the versions pin. Duplicate Five-Explicit
  pillars and duplicate let names refuse as ambiguous instead of replacing an
  earlier declaration. Extraction only — validity POLICY lives in `admission`
  (below).
- `lower::lower` — high-level verbs (`optimize-shape`, `simulate-pour`)
  expand to explicit IR with an inspectable trace naming every injected
  default (progressive disclosure with nothing hidden); idempotent;
  malformed verb usage refuses with the verb's span.
- `IrError` — span + stable kind code + detail + fix hint (refusals
  teach). `IR_VERSION` — the language version this build reads/writes.
- `query` (addendum Proposal 8 — declarative query language v0): a query is
  `(QoI, Target, budget_usd, deadline_s)` where `Qoi` is a fixed MENU —
  `MaxOverRegion`, `Integral` (linear), `Exceedance` (probabilistic, needs a
  named environment) — each advertising `QoiMeta { linear, adjoint_available,
  ladder_applicable }` for the planner and a `value_dims(field_dims)`
  (max→field dims; integral→field·m³; exceedance→dimensionless). `Target` is
  `Tolerance{value,dims}`, `Confidence(f64)`, or `ToleranceAndConfidence`.
  `Query::admit(&FieldRegistry) -> QueryAdmission` type-checks a query in
  constant time (no solves) over six fixed-order checks — `query.field`
  (the QoI's field must exist), `query.budget` (finite positive $),
  `query.deadline` (finite positive s), `query.confidence` (strictly in
  `(0,1)` — 100% is uncertifiable), `query.target` (finite positive
  tolerance), `query.dimensions` (tolerance dims == QoI value dims; exceedance
  threshold dims == field dims) — emitting the admission bead's `Finding`s
  with ranked teaching `RankedFix`es. `Query::from_node`/`to_node` give the
  `(query …)` IR surface, round-tripping under `same_shape`. This is the
  addendum's declarative surface; imperative solver access is the
  internal/expert path.

- `admission` (the gp3.5 bead): `admit(node, &AdmissionContext) ->
  AdmissionReport` runs six timed dimensions — Five Explicits structure,
  dimensional analysis (fs-qty dims inferred bottom-up through `+ - * /
  min max` and comparisons; unknown verbs never false-reject), budget
  feasibility (fs-plan cost models over exact numeric `:dof`/`:size`/`:modes`
  features, with malformed or duplicate explicit features refused rather than
  priced as unit size; p90
  totals vs the `(budget (wall …))` bound, with RANKED cost-model-derived
  fixes: coarsen / surrogate-screen / relax), capability sufficiency
  (finite non-negative session grants, session-token and self-contained
  explicit globs vs namespaced verbs, and finite declared asks). Capability
  fields are exact keyword/value pairs; operator grants are exact names or
  namespace wildcards of the form `foo.*`. Wall/memory budget clauses have
  exact arity; structured operator-specific budget clauses remain extensible
  until the catalog lands. Chart
  routability (fs-geom Router as an admission predicate with the
  RouteRefusal's own fixes attached), and regime gating (explicit
  `(assert (regime.allows …))` plus `flux.*` verbs checked against an
  fs-regime report; policy-graded Reject/Warn). Findings carry spans,
  diagnoses, and `RankedFix { action, predicted_wall_s, qoi_impact }`.

- `planner` module (addendum Proposal 8, bead lmp4.16; [F], behind
  `ladder-planner` → optional fs-verify dep): the GREEDY LADDER WALK —
  not a general planner (Governance Rule 1). The operator menu
  `{cache, speculate, solve-rung, refine-to-target, climb}` runs
  greedily with costs LEARNED from telemetry (`CostTable`; cold
  entries fall back to the conservative default). Speculation
  verifies a prolongated coarse answer WITHOUT solving; refinement is
  the textbook equidistribution criterion (split every element above
  the per-element target `tol²/n` with per-element depth from its own
  gap). A discharged answer's bound is a REAL equilibrated enclosure
  (VERIFIED color) and never violates the query's certified-accuracy
  contract; the cannot-discharge boundary refuses with the best
  achieved certified interval (never a false in-budget answer);
  operator choice tie-breaks deterministically (G5 replay).
  `baseline_uniform` is the fixed control the kill criterion measures
  against.

- `anytime` module (addendum Proposal 8, bead lmp4.17; ships behind
  `ladder-planner` but its CONTRACT survives even a frozen planner —
  the product win): `run_anytime` drives the planner over an
  increasing budget ladder — the first rung is the IMMEDIATE wide
  certified interval, tightening is MONOTONE, every step carries its
  Proposal-3 color and a PRICED "what would tighten this" hint
  (`tighten_hint`: gap extrapolation naming the next menu move and the
  hot region where refinement concentrated; cold telemetry degrades to
  the generic priced form). REFUSAL semantics: an undischargeable
  query returns the achieved certified interval, the price of the gap,
  and the explicit no-point-estimate clause — never a silent
  best-effort number. Replays reproduce trajectories bit-for-bit (G5).

## Invariants

1. Isomorphism: `parse(print(x))` has the same shape as `x`, per syntax
   and across syntaxes (property-tested on generated programs and the
   Appendix C fixtures).
2. Both parsers are total: any input yields a value or a structured error
   with an in-bounds span; recursion is depth-capped (no stack overflow).
3. No silent reinterpretation: numeric-leading tokens either fully parse
   as int/float/quantity/count or refuse; non-finite literals refuse.
4. Count authority is exact end to end: decimal text cannot round into a
   different byte/core claim, checked unit scaling precedes admission, and
   `SessionCapability` carries integer memory/core grants without an `f64`
   projection.
5. Lowering is explicit, inspectable, and idempotent; the trace names
   every injected default.

- Admission determinism: same study + context → byte-identical
  `diagnosis()`; findings sorted (check, span).
- Admission latency is milliseconds-class on Appendix C studies (six
  checks timed individually; conformance logs and bounds the total).
- Zero false admits on the violation zoo; missing verifiers (no Router,
  no RegimeReport) degrade to WARN verification-gap findings, never to
  silent admits of violations they could not check.

## Error model

All fallible APIs return `IrError` (span, stable `IrErrorKind::code()`,
detail, hint). Never panics across the crate boundary (fuzz-tested).

## Determinism class

Parsing, printing, and lowering are pure functions of their input text —
bit-deterministic across runs, thread counts, and ISAs.

## Cancellation behavior

No compute loops beyond input length; parsing is bounded by source size
and the depth cap. No `Cx` integration needed at this layer.

## Unsafe boundary

None. Safe Rust only.

## Feature flags

- `ladder-planner` [F] (default OFF) — the greedy ladder-walk planner
  (`dep:fs-verify`); disabled until its Gauntlet tier + kill metric are
  green. Gates the `planner`, `plancal`, and `anytime` targets. All
  other v1 behavior is `[S]` default-path.

## Conformance tests

`tests/conformance.rs`: Appendix C spout + frame studies as verbatim
fixtures (names, seeds, locks, lets, and typed-noun counts asserted);
isomorphism property over 200 generated programs plus the fixtures
(s-expr, JSON, and cross-syntax cycles); 8000-parse garbage battery with
in-bounds spans and non-empty hints plus 100k-deep nesting rejections;
span-accuracy cases (bad seed, bad quantity); verb lowering explicitness,
trace content, idempotence, and structured refusal; version-pin
round-trip through both syntaxes.

`tests/query.rs` (suite `fs-ir/query`, addendum Proposal 8): the wedge QoI
menu is expressible with correct metadata; `value_dims` follows the
functional; well-posed queries admit; the FIVE ill-posed classes each reject
on a distinct check with a teaching fix (zero budget, past deadline, 100%
confidence, field-absent-from-design, self-contradictory dimensions), plus
off-dimension exceedance thresholds and integral-tolerance-needs-volume-dims;
multiple faults are reported together; admission is deterministic (identical
verdict on replay); and the `(query …)` IR form round-trips (tolerance,
exceedance+confidence) with a teaching error on a non-query form.

`tests/admission.rs` (suite `fs-ir/admission`): ad-001 Appendix C admits
cleanly + ms-class latency + determinism; ad-002 five-study violation zoo
(all rejected on the right dimension, fixes attached); ad-002b malformed,
negative, non-finite, empty, and duplicate resource grants/pillars fail closed,
and self-contained operator grants constrain the study; ad-003 dimensional
spans pinpoint the offending operand, products stay legal; ad-004
BudgetInfeasible with ranked cost-derived fixes + fix-quality harness
(applying fixes admits); ad-005 Router-backed feasibility; ad-006 regime
gating with alternatives + policy grading; ad-007 2000 mutants + all
truncation prefixes never panic (a fuzz-found scanner panic became a
structured refusal).

## No-claim boundaries

- No operator catalog or per-operator semantic versions — gp3.6; the
  `IR_VERSION` constant covers the language only.
- JSON `\uXXXX` escapes cover Unicode scalar values only (surrogate
  pairs are rejected with a structured error).
- The verb table is v1-small (optimize-shape, simulate-pour); verbs are
  data to extend, not a framework.
- Qty literals must be written in units fs-qty accepts; information
  units are Counts, not quantities, by design.
- Admission's dimensional pass covers arithmetic/comparison heads;
  verb-signature dimension contracts (per-operator expected dims) land
  with the operator registry.
- Chart requirements are supplied by lowering/callers; admission does not
  yet derive them from raw study text.
- `SessionCapability` is admission's view of a token; issuance,
  revocation, and idempotency keys are fs-session's bead (gp3.7).
  A self-contained `(capability ...)` clause supports static planning and
  source-level admission only; it does not mint runtime authority. Plan §11.3's
  session token remains mandatory before execution.
- IR v2 changes Count identity from binary-float-backed count atoms to exact
  integer/decimal atoms. V1 canonical artifacts must be reparsed and re-emitted
  under v2 before their new identity is recorded; no silent v1 hash migration
  is claimed.
- The query language is v0: a FIXED QoI menu (max/integral/exceedance), not
  a general program surface. `Query::admit` type-checks well-posedness and
  dimensions ONLY — it does NOT plan, cost, or execute a query (the greedy
  fidelity-ladder planner and the anytime/refusal result semantics are
  separate addendum beads). Field dimensions come from a caller-supplied
  `FieldRegistry` (the design's typed fields, Proposal 13); this module does
  not itself derive fields from geometry. `budget_usd` is a priced dollar
  budget distinct from the wall/memory/core grants of the `(budget …)` study
  clause. The returned answer's COLOR (verified/validated/estimated) is
  attached by the query result, not here.

## No-claim boundaries (planner)

- v0 discharges the verifier's 1-D elliptic kernel class; the 2-D
  cutfem DWR (fs-dwr) and real physics kernels plug into the same walk
  as the ladder registry grows rungs.
- Cost units are solved cells (the flywheel's telemetry currency);
  wall-clock costs arrive with the perf-CI lane.
- Confidence targets (`Target::Confidence`) are the e-process beads'
  contract; v0 discharges tolerance targets.
- The kill measurement (>=2x vs mid-rung+uniform; measured 4.31x on the
  steep-feature fixture) is per-fixture evidence, not a universal claim
  — the wedge query set re-measures it as kernels land.

## No-claim boundaries (anytime)

- The hint's price is an O(h) extrapolation from the achieved bound —
  an estimate for teaching, not a certified cost bound; Proposal C's
  full value-of-information ranking replaces it when C lands (the soft
  dependency the bead names).
- Tile-boundary interruptibility rides the planner's operator
  granularity (each budget rung is a clean cut); sub-operator
  cancellation lands with the fs-exec tile integration.
