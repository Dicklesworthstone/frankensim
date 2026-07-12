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
- `VersionedProgram` — the persisted/replayed artifact boundary shared by both
  syntaxes. It canonically wraps a program as
  `(frankensim-ir :version 2 :program <node>)`, binds the language version into
  serialized identity, and refuses older/newer semantics unless a caller first
  performs an explicit audited migration. Bare parsers remain syntax-only.
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
  `RouteRefusal`'s own fixes attached; malformed/spec-mismatched oracle
  authority and bounded-search exhaustion reject distinctly; a route containing
  any converter not declared certificate-backed is estimated and cannot
  authorize admission), and regime gating (explicit
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
  contract. `ProblemFamily`, `CachedAnswer`, and `CostTable` have checked
  constructors and private state; `plan -> Result<PlanOutcome, PlanError>`
  validates finite theta/tolerance/budget, a non-empty strictly increasing
  non-zero rung ladder, family coefficients/boundaries, meshes, candidates,
  telemetry, verifier enclosures/refusals, and every arithmetic result. Bounded
  FEM refusals retain their structured `Fem1dError` source in `PlanError`
  instead of being flattened into a plausible result. Cache-declared
  bounds are never trusted: hits are independently re-verified against a
  lower-layer canonical MMS class identity before discharge. The family owns
  an admitted immutable `fs-verify::MmsClass`; theta scaling constructs a new
  admitted class, and cache keys bind that class's versioned canonical bytes
  instead of independently serializing planner coefficients. Cache records
  normalize signed-zero mesh and nodal values before retention. Ladder transfers
  use deterministic P1 interpolation over the actual coarse coordinates,
  including non-dyadic ladders and adaptive-to-uniform moves. Actual solve/speculation
  costs are admitted before execution, so spend never exceeds budget. Learned
  costs rank the zero-cost refine/climb transition only; they never veto exact
  affordable work. Transition telemetry is pending until its first downstream
  verification/solve executes, records that actual compute cost before the
  resulting certificate checkpoint, and is dropped when admission aborts before
  downstream work. A rejected climb speculation completes that transition;
  any subsequent fine-rung solve is separately charged as `SolveRung`. The
  family degree is capped at five (six coefficients), matching the exactness
  envelope of the verifier's five-point squared-residual quadrature; its work
  preflight counts all five quadrature evaluations per coefficient/cell. The
  homogeneous trace, six-coefficient cap, signed/trailing-zero
  normalization, derived forcing, and stable identity all come from the single
  lower-layer fs-verify admission type. Parameter scaling is refused if
  independently rounded products no longer sum exactly to zero at `x=1`. The
  cannot-discharge boundary refuses with the best achieved certified interval;
  if no solve is affordable it returns `RefusedWithoutAnswer` with no interval
  or color rather than fabricating evidence. Operator choice tie-breaks
  deterministically (G5 replay).
  Every finite audit-log bound retains a private-constructor
  `VerifierCertificate`: the guarded `fs-evidence` color, stable verifier-family
  identity, and reconstructed-flux hash travel together instead of reminting a
  bare `Color::Verified` from a discarded scalar. The cell cap is derived as
  one less than the lower-layer `fs-verify` node cap. Resource admission caps cells,
  family coefficients, fidelity-rung count, and their coefficient-by-cell-by-
  quadrature-point work product. Uniform/adaptive mesh, indicator,
  prolongation, trajectory, and
  family-scaling allocations use fallible reservation; exact downstream solve
  cost and the combined compute envelope are admitted before mesh allocation.
  `baseline_uniform` is the fixed control the kill criterion measures
  against and is fallible under the same numerical/family boundary.

- `anytime` module (addendum Proposal 8, bead lmp4.17; ships behind
  `ladder-planner` but its CONTRACT survives even a frozen planner —
  the product win): `run_anytime_observed` drives one cumulative planner
  execution through operational certification/pre-work checkpoints. It emits
  each affordable budget-rung certificate before work, allocation, cache
  insertion, or telemetry for a later rung. An observer can return `Stop`; the
  returned report is then the exact deterministic prefix and no later side
  effect occurs. `run_anytime` is the collector wrapper whose observer always
  continues. Total work never exceeds the final budget, no certificate appears
  before its cumulative cost is affordable, and retaining the best checked
  certificate makes tightening MONOTONE. The first affordable rung is the
  IMMEDIATE wide certified interval, and every step carries its guarded
  Proposal-3 color plus verifier family/flux identity and a PRICED "what would
  tighten this" hint
  (`tighten_hint`: gap extrapolation naming the next menu move and the
  hot region where refinement concentrated; cold telemetry degrades to
  the generic priced form). REFUSAL semantics: an undischargeable
  query returns the achieved certified interval, the price of the gap,
  and the explicit no-point-estimate clause — never a silent
  best-effort number. `run_anytime -> Result<AnytimeReport, PlanError>` rejects
  empty, non-finite, non-positive, or non-increasing budget ladders before work.
  A valid rung too small for the initial solve contributes no trajectory point;
  if it is the final rung the report explicitly says that no certified interval
  or color exists. `tighten_hint` is likewise fallible and cannot emit NaN/∞
  gap prices. Hints use the cost table at emission time and cannot consume
  telemetry from future work. Budget-ladder length is resource-capped and its
  trajectory allocation is fallible. Replays reproduce trajectories,
  certificate identities, and observer-selected prefixes bit-for-bit (G5).

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

Syntax/study/lowering APIs return `IrError` (span, stable
`IrErrorKind::code()`, detail, hint). Feature-gated planner/anytime APIs return
`PlanError`, and valid but under-budget queries return structured
`PlanOutcome` refusals. Neither boundary panics on malformed caller data.

## Determinism class

Parsing, printing, and lowering are pure functions of their input text.
Planner replay is deterministic for the same family, query, ladders, cache
contents, and learned cost table: fixed operator ordering, exact cache-key
framing, coordinate-ordered prolongation, and deterministic tie-breaking.

## Cancellation behavior

Parsing is bounded by source size and the depth cap. The feature-gated planner
does perform numerical work. Its operational anytime API is synchronously
stoppable between certified budget rungs: `Stop` prevents later planner work,
allocation, cache insertion, and telemetry. Sub-operator cancellation still
lands with fs-exec integration; no claim is made that a running solve or
verification can be interrupted inside its admitted coefficient-by-cell work
envelope.

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

`tests/planner.rs` + planner unit tests (`ladder-planner`, G0/G3/G5): existing
accuracy, kill-ratio, cache, refusal, calibration, canonical family/cache
identity equivalence and semantic-mutation separation, and replay checks plus
empty/zero/non-monotone rungs; non-finite theta/tolerance/budget; malformed
family/mesh/candidates; poisoned cost samples; unaffordable initial solves;
independent replay of a falsely certified cache answer; non-dyadic
prolongation; adaptive-to-uniform coordinate interpolation; pessimistic learned
costs that cannot veto affordable exact work; bounded family/rung/cell and
combined coefficient-by-cell resource drivers; pre-allocation budget refusal;
verifier authority retained on every finite audit bound; and aborted
transitions that do not enter observed telemetry.

`tests/anytime.rs` (`ladder-planner`, G0/G5): monotone verified trajectories,
priced refusal/hints, cache termination, empty/zero/non-monotone budget/rung
ladders, malformed scalar/hint inputs, and explicit no-interval/no-color output
when the final budget cannot fund one solve. A counting-cache regression proves
that an entire budget ladder executes the planner once. Operational observer
regressions prove callback order, actual rung/spend receipts, contemporaneous
hints without future telemetry, verifier-family/flux identity retention, and
that `Stop` prevents later work telemetry and cache insertion while retaining
telemetry for a completed speculative transition.

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
- Router certification is currently a validated declaration on
  `ConverterSpec`, not an authenticated checker/ledger receipt. Admission
  refuses explicitly estimated routes, but full opaque admitted-converter
  authority remains part of the scientific-evidence migration; callers must not
  interpret the declaration Boolean as independent proof.
- `SessionCapability` is admission's view of a token; issuance,
  revocation, and idempotency keys are fs-session's bead (gp3.7).
  A self-contained `(capability ...)` clause supports static planning and
  source-level admission only; it does not mint runtime authority. Plan §11.3's
  session token remains mandatory before execution.
- IR v2 changes Count identity from binary-float-backed count atoms to exact
  integer/decimal atoms. V1 canonical artifacts must be reparsed and re-emitted
  under v2 before their new identity is recorded; no silent v1 hash migration
  is claimed.
- Bare `sexpr::parse`/`json::parse` intentionally do not infer an artifact
  version. Persisted or replayed programs must use `VersionedProgram`; callers
  that ledger a bare AST have no version-binding claim.
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
- The v0 synchronous numerical envelope additionally refuses when polynomial
  coefficients times mesh cells exceeds `MAX_POLYNOMIAL_CELL_WORK`; this is a
  deterministic resource guard, not a wall-time certificate.
- Cache storage/transport authentication remains the content-addressed store's
  responsibility. This planner treats cache data as untrusted and re-verifies
  its numerical claim, but does not authenticate who wrote the entry.
- The v0 family boundary checks finite polynomial structure and homogeneous
  endpoints; it does not prove that arbitrary caller-supplied polynomial
  semantics represent the intended physical model.
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
- Operational interruption is rung-granular. A callback can stop before the
  next operator and receive a clean deterministic prefix, but sub-operator
  cancellation lands with the fs-exec tile integration.
