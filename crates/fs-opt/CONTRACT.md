# fs-opt — CONTRACT

## Purpose and layer

L4 (ASCENT). The optimization problem IR (plan §9.1): optimization
problems ARE DATA — typed objective/constraint graphs over
manifold-valued variables, storable, hashable, replayable, and
constructible INCREMENTALLY with validation at every step (the
agent-ergonomics property). The IR REPRESENTS physics and stochastic
structure; FLUX/UQ execute it.

## Public types and semantics

- `ProblemBuilder` → `Problem`: hash-consed expression arena (repeated
  subexpressions return the SAME `NodeId` — CSE by construction).
  `Problem` is SEALED: fields are crate-private; the public surface is
  read-only accessors (`vars()`, `exprs()`, `objectives()`,
  `constraints()`, `tags()`, `budget()`) plus CHECKED id-indexed
  accessors (`expr`, `variable`, `shape`, `node_dims`, `class`,
  `node_depth`, `reachable`) plus the sealed graph's
  `total_admission_work`; unknown ids refuse instead of panicking. Every
  builder constructor validates through the SAME versioned leaf rules
  the admission validator uses (`derive_expr` + leaf policies):
  shapes (`Scalar`/`Vector(n)`), fs-qty DIMENSIONS (add/compare need
  equal dims; mul/div/dot/norm_sq combine exponents with CHECKED
  refusal — no silent saturation; powi scales with checked refusal;
  sqrt halves even exponents; transcendentals demand dimensionless),
  node/variable existence, parameter ranges, FINITE constants (exact
  bit pattern retained in the refusal), finite NONNEGATIVE objective
  weights (`-0.0` refused; `Sense` carries direction), manifold
  policy (`Rn` dim ≥ 1, `Sphere` ambient ≥ 2, `Stiefel` 1 ≤ p ≤ n,
  CHECKED point/tangent formulas), validated tags (nonzero capped
  fidelity levels, finite open-interval chance probabilities, typed
  bilevel references), scalar-only objective/constraint roots, and
  versioned per-item/aggregate caps (`AdmissionCaps`,
  `ProblemBuilder::with_caps`), including graph depth, conservative
  retained/canonical bytes, exact input-wire bytes, deterministic
  admission work (one item visit plus expression child edges), and
  target-`usize` packed point storage. External strings are checked
  against per-field and aggregate byte/work/depth limits while borrowed;
  PDE/UQ ownership and insertion of the fixed-size, collision-checked
  fingerprint happen only after every cap passes. Cap+1 rejection leaves
  vector capacities, owned-string capacity, intern tables, ids,
  byte/work/storage totals, budget, and ordering unchanged.
- `Problem::admit` / `admit_with_caps` → `ProblemAdmission`: the
  single versioned re-validation chokepoint (schema
  `ADMISSION_SCHEMA_VERSION`). Re-derives every node's
  shape/dims/class from the expression list (cache agreement),
  proves reference validity and acyclicity from arena ordering,
  re-checks every leaf policy and cap. Cheap aggregate count/alignment
  failures return their complete deterministic preflight section before
  proportional graph work. Within that count envelope, validation work
  and retained bytes accumulate only until their first cap crossing, and
  depth is re-derived with a max+1 early exit before shape/class table
  allocation. The remaining admitted-size sections then scan
  deterministically — or mint
  the `ProblemSemanticId` and lists quarantined legacy identities on
  success. Builder output always admits (same rules, pinned by test).
- Identity is DOMAIN-SEPARATED with no implicit conversion:
  `ProblemSemanticId` (BLAKE3 over the domain-tagged canonical v3
  body; minted by admission; publicly constructible from a full-width
  hash only as a bilevel REFERENCE), `WireContentId` (BLAKE3 over the
  domain-tagged exact artifact bytes; minted ONLY by
  `serialize_with_id`/strict parsing — programmatic construction
  never manufactures one), and `LegacyProblemHash` (the quarantined
  FNV-1a 64; correlation and corruption tripwire only, NO authority).
- Node kinds: arithmetic (`add/sub/mul/div/neg/powi/sqrt/exp/ln/
  tanh`), vector reductions (`dot/norm_sq/component`), kinks
  (`min/max` are scalar-only, `abs` is scalar; all are C0),
  `pde_residual` (FLUX study reference with
  ADJOINT AVAILABILITY metadata), `expectation`/`cvar`/`quantile`
  (UQ config references; CVaR/quantile are C0).
- `Class` propagation: bottom-up minimum of children and each node's
  own contribution — "this objective is non-smooth through that
  min()" is knowable at BUILD time. `Problem::route(family)` refuses
  L-BFGS/Newton on C0 graphs NAMING the poisoning node, refuses
  gradient families on adjoint-less PDE nodes NAMING the study, and
  admits subgradient/gradient-free families. `class_trace()` names
  every node's class.
- `Manifold` (`Rn`, `Sphere`, `So3` as unit quaternions, `Stiefel`)
  with `point_dim`/`tangent_dim`/`param_dim` and `retract` (Rn
  translation, Sphere normalize, SO(3) quaternion exponential,
  Stiefel Gram-Schmidt/QR) — the metadata the gradient stack
  consumes. The public retraction boundary validates the manifold and
  requires exact point and parameter storage before any indexing or
  zip-based arithmetic; malformed slices refuse with a typed input and
  expected/actual lengths. Point, step, candidate, and output components
  must remain finite. Sphere/SO(3) bases must have unit squared norm and
  Stiefel bases must have an identity Gram matrix within `1e-10`; candidate
  squared norms at or below `1e-24` refuse as singular instead of being
  clamped and normalized into fabricated points. `descend_fn`/`descend_ir`
  are the TOY consumers proving iterates stay ON their manifolds. Both
  descent entry points are LEAF-GATED before the first objective
  evaluation or any descent arithmetic (bead j3vb5 / review High #6):
  the manifold validates through the admission policy, the start
  point must match `point_dim` with FINITE components (offending bit
  pattern retained in the refusal), and `DescentOptions` must carry
  finite positive `fd_h` and `lr` (a NaN/zero step would divide
  through the FD quotient; a negative rate would silently ascend).
  The internal objective seam is fallible: non-finite closure results
  and typed IR evaluation failures propagate without publishing a
  report or converting the refusal into a panic. Each finite-difference
  step reserves its full two-probes-per-parameter gradient plus terminal
  valuation before the first probe, so budget exhaustion cannot leave a
  partially spent gradient. Cancellation is polled before/after f0,
  probe evaluations, retraction boundaries, and final publication.
- Structure: multi-objective (weights), constraint KINDS (`EqZero`,
  `LeZero` — semantics/repair are fs-constraint's), `ProblemTag`
  (multi-fidelity, chance-constrained, bilevel via typed
  `BilevelRef`: `Semantic(ProblemSemanticId)` full-width or
  `LegacyFnv(LegacyProblemHash)` QUARANTINED — never interchangeable,
  never widened), `EvalBudget` (P4, enforced by consumers).
- Serialization: `serialize` writes the canonical six-base `fsopt v3`
  line-based text form and `parse` round-trips it BITWISE (floats travel as
  bit patterns); `serialize_with_id` additionally mints the artifact's
  `WireContentId`. In v3, `tag bilevel <64-hex>` carries a semantic id and
  `tag bilevel_legacy <16-hex>` keeps a quarantined legacy identity
  EXPLICIT. v2 input (`tag bilevel <16-hex>`) remains readable with the
  identity quarantined in the type; the v2 → v3 step is a pure
  identity-typing re-encoding with no semantic receipt.
  `parse_with_version` also accepts strict explicit `fsopt v1`
  bytes emitted by either known historical v1 token writer, maps the absent
  amount exponent to `mol = 0`, and returns an immutable
  `DimensionCrosswalkReceipt`. Headerless input is refused because no
  historical writer emitted it and it has no authoritative schema identity.
  The receipt binds BLAKE3 hashes of the complete old artifact and the
  complete canonical V2 target artifact (where the five-to-six dimension
  semantics land; re-derivable via the fallible
  `canonical_v2_migration_target`) under
  the sole `AppendMoleZero` rule; `parse` refuses v1 because it
  cannot return that mandatory evidence. `ParsedProblem` keeps its fields
  private and exposes read-only provenance accessors (`source_version`,
  typed `source_hash`, `wire_content_id`, `migration`) plus `into_parts`,
  so an inconsistent provenance tuple cannot be constructed through
  the public API. Every admitted artifact has exactly one terminal `hash`
  directive; `problem_hash` (in-house FNV-1a 64 over the canonical body,
  returned as the quarantined `LegacyProblemHash` type) is correlation
  only — semantic identity is admission's `ProblemSemanticId`. Legacy FNV
  integrity is verified over the historical v1 body
  without normalization. The reader rebuilds an exact canonical v1
  artifact using both known historical v1 token encodings and requires complete
  byte equality with one of them before it may issue the sole-rule receipt;
  v2/v3 receive the same complete-byte comparison against their canonical
  writers. CRLF,
  blank/missing-final-newline forms, noncanonical token spellings, malformed
  escapes, wrong IDs, extra fields, and duplicate/missing budget directives
  therefore refuse even when their recomputed FNV is internally consistent.
  Semantic bilevel references have no v2 spelling, so that migration
  writer refuses with `WireIncompatible` rather than emitting a
  self-invalid downconversion. Parsing preflights exact artifact bytes
  and the minimum of structural/work directive envelopes before allocating
  its line table. Token decoding is two-pass: malformed percent escapes
  and any decoded cap+1 token refuse before its output buffer exists;
  decoded UTF-8 validity is checked when that bounded buffer becomes a
  `String`. The parser then REBUILDS through the validating builder and
  verifies the integrity hash —
  tampered, oversized, or ill-typed files refuse as `Parse` with line
  numbers.
- `eval`: memoized evaluation of algebraic subgraphs; the sealed root
  depth and aggregate-work receipts are checked against the default
  admission schedule before memo allocation. Supplied runtime bindings
  are checked for exact manifold point length and finite components first;
  a refusal identifies the exact variable, component, and IEEE-754 bits.
  Every computed scalar/vector node result is then checked for finiteness
  with node/component attribution before it enters the memo or becomes a
  public result. The walk itself is EXPLICIT-STACK (reachability worklist
  + bottom-up arena-order sweep;
  bead frankensim-xf8v7) so no admitted graph — at the depth cap or
  otherwise — can overflow the call stack. `ir::children` is public so
  downstream evaluators can drive the same iterative traversal.
  PDE/stochastic nodes refuse with `Unevaluable` NAMING their executor.
- `GoodhartGuard` (addendum Proposal D): treats an optimizer `Endpoint`
  (`design`, `objective`, `label`; `from_descent` bridges `DescentReport`)
  as an adversarial example. A FIXED four-step escalation ladder
  (`EscalationKind::ORDER` = rung-k+1, cross-representation, δ-perturbation,
  estimator-independence) runs pluggable `EscalationStep`s; each yields a
  `StepOutcome` (`Passed` / `Vetoed{reason}` / `NotPerformed{reason}`).
  Aggregation: any veto → `GuardStatus::Failed` (+ a `GuardFinding` the
  caller files as a tombstone/bug report); else any unregistered step →
  `Provisional`; else `Cleared`. `is_honored()` is true ONLY on `Cleared`
  — the endpoint certificate stays provisional on any skipped check (never
  a false clear). `converged_and_guard_cleared(converged, &report)` is the
  amended contract ("converged AND guard-cleared"). One concrete step ships:
  `DeltaPerturbationStep` re-evaluates a supplied objective at deterministic
  `±δ` coordinate probes and vetoes a found-better point (not a true optimum)
  or a sharp crack (optimum not in a smooth basin), failing closed on
  non-finite values.
- The game module (RE.Q1) is the versioned object language for quantified
  reach-avoid and viability games. Initial, target, unsafe, control,
  disturbance, parameter, hybrid-mode, and hybrid-event domains have
  non-interchangeable identity types and retain exact model-version,
  state-space, frame, unit, and time-unit context. The ordered quantifier
  prefix prints each exists/forall clause and full typed domain identity.
  Information patterns state hidden, initial, current, history, delayed, or
  forbidden-future observations; a dependency-free open-loop strategy may
  carry the empty pattern. Strategy descriptions distinguish open-loop,
  state-feedback, nonanticipative, hybrid-mode, set-valued, and unresolved
  policies and bind their exact temporal dependencies. Finite/infinite
  horizons, stopping rules, deterministic, differential-game, finite-hybrid,
  and admitted-DAE model classes, proof polarity, composition, and analysis
  budgets are identity-bearing semantics. Admission canonicalizes unordered
  grants, dependencies, strategies, and parallel/product components while
  retaining quantifier order and sequential order plus multiplicity. Its claim
  availability is only later-checker eligibility: unsupported infinite
  horizons, unresolved Zeno behavior, regularized hybrid models, unresolved or
  set-valued strategies, and Unknown polarity remain explicitly Unknown.

## Invariants

1. Seeded ill-typed constructions refuse with teaching text naming
   ops/nodes, and a 600-op fuzz storm matches an independent validity
   model exactly (opt-001).
2. build→serialize→parse through canonical `fsopt v3` yields an IDENTICAL
   problem; hashes are stable across identical builds, differ across edits,
   and guard integrity. Exact explicit five-dimension v1 inputs from both known
   historical writers remain readable through the receipt-bearing API with
   `mol = 0`; their
   historical FNV is unchanged, their complete old/new artifacts are bound by
   BLAKE3, dimension arities are strict, and missing, duplicate, nonterminal,
   or malformed hash directives fail closed. Recomputed-FNV adversarial tests
   also lock rejection of extra fields/operands, wrong variable IDs, non-Boolean
   flags, malformed or invalid-UTF-8 escapes, CRLF/blank/missing-final-newline
   forms, missing/duplicate budgets, and noncanonical v2 spellings (opt-002 plus
   focused serializer unit tests).
3. Hash-consing gives CSE identity; substitution commutes with
   evaluation BITWISE; `neg∘neg` and `min(x,x)` are bitwise identities
   (opt-003).
4. Class propagation + routing: kinks poison smooth families with the
   node named; adjoint-less PDE nodes refuse gradient families with
   the study named; the class trace covers every node (opt-004).
5. The toy Riemannian descent consumes manifold metadata: Sphere
   reaches the analytic minimizer staying unit, SO(3) aligns with a
   unit quaternion throughout, Stiefel stays orthonormal to 1e-10 and
   finds the top invariant subspace (opt-005).
6. P4/P7: the attached budget stops descent with a RECEIPT (not an
   error), never exceeds its cap, reuses the already-counted initial
   value when no step lands, and reserves the complete FD gradient plus
   terminal evaluation before spending any probe. Raw manifold descriptors
   and `x0` lengths refuse before the objective closure is called;
   cancellation returns the teaching error; PDE/stochastic nodes name
   their executor when asked
   to evaluate (opt-005/006).
7. G3 unit rescaling: the live `descend_fn` step is equivariant when a
   one-dimensional quadratic's start, target, and finite-difference step are
   coherently rescaled by a nonidentity power of two. The final coordinate
   scales by `s`, both objective receipts scale by `s²`, and evaluation,
   step, and budget receipts remain exact (`tests/metamorphic.rs`).

8. RE.Q1 admission never rewrites quantifier order. In particular, a universal
   disturbance followed by an existential open-loop control refuses as a
   trajectory-clairvoyant lowering, while an explicitly nonanticipative policy
   may use only observations and delays granted by its information pattern.
   A positive-lag grant never authorizes time-zero `InitialOnly` or `Current`
   access; a zero-lag grant authorizes both, but not retained history. Future
   access always refuses. Finite viability uses fixed-horizon stopping so its
   throughout-horizon obligation cannot be truncated, and custom/hybrid stops
   without an encoded terminal outcome refuse. A nonanticipative policy's
   finite memory bound must fit the retained-strategy-node budget. V1 fixes the
   player roles to existential control and universal disturbance; opposite
   polarities refuse instead of being sealed with contradictory strategy
   ownership. Hidden, delayed, mode, event, model-version, frame, unit, and
   typed-domain mismatches fail closed before publication. Equivalent ordering
   of semantic sets replays to one receipt, but swapped quantifier order or
   strategy classes produce distinct identities.

## Error model

`OptError` teaching errors throughout: unknown ids, shape/dimension
mismatches and dimension overflow (with exponent vectors shown; both
scaling `DimOverflow` and combining `DimSumOverflow`), non-dimensionless
transcendentals, odd-sqrt dims, bad parameters/indices, non-scalar
roots, `ManifoldInvalid` (violated policy named), `NonFinite` (payload
name + exact bit pattern), `BindingNonFinite`/`EvalNonFinite` (runtime
location + exact bits), `RetractionLen`/`RetractionNonFinite`/
`RetractionDomain` (input, manifold rule, location, and measurement),
`CapExceeded` (cap name + count + limit), `BindingCount`/`BindingLen`
(declared vs supplied),
`WireIncompatible` (historical version + unrepresentable typed construct),
`NonsmoothForFamily` (node + kind + class),
`NoAdjoint` (node + study), `Unevaluable` (node + executor), `Parse`
(line + what), `Cancelled`, `BudgetExhausted` (spent count receipt).
Whole-problem re-validation refuses with a deterministically ordered
`AdmissionReport`: the cheap count/alignment preflight gathers all of
its findings before early refusal; work, retained bytes, and graph depth
then fail at the first aggregate crossing, and admitted-size section
scans are index ordered.

Game admission adds deterministic typed refusals for hard pre-work caps,
decoded schema drift, set/model context, quantifier domains, anticipative or
unavailable observations, strategy causality, horizon/stopping conflicts,
composition, cancellation, and canonical identity.

## Determinism class

Fully deterministic: `BTreeMap` interning, index-ordered ids, bitwise
float serialization, in-house FNV hashing, domain-separated BLAKE3
identity minting, no time or randomness.
Identical build sequences give identical problems, semantic ids,
hashes, and bytes; identical rejections give identical reports
(opt-002/003 and adm-004/005/006 are the trip-wires).

RE.Q1 game identities are bit-deterministic: semantic sets are canonically
ordered, quantifier and sequential-composition order remain significant,
floating signed zero is normalized, and model, units, information, strategy,
polarity, horizon, stopping, composition, and budget are domain separated.

## Cancellation behavior

`descend_fn`/`descend_ir` poll `cx.checkpoint()` before and after f0,
each finite-difference probe and retraction boundary, the landed step,
terminal evaluation, and report publication. Cancellation returns
`OptError::Cancelled` without publishing a partial report. Budget
exhaustion is a RECEIPT (`budget_stopped` in the report), not an error;
the iterate remains valid, no partial gradient is spent, and `evals`
never exceeds the positive cap (P4).

Game admission polls before proportional scans, during bounded information,
strategy-dependency, and composition traversal, and at identity publication.
Cancellation returns a typed Cancelled issue and no receipt.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`parked-ir-battery` — compiles the PARKED numerics-spine draft
battery (`tests/ir_battery.rs`, see its header): it targeted a
parallel draft of this crate that lost the crate-structure race; the
draft modules (`graph.rs`, `manifold.rs`, `riemann.rs`, `sexpr.rs`,
`expr.rs`) remain in `src/` UNREFERENCED for harvest (notably: exact
reverse-mode gradients and the s-expression re-validating parser).
Off by default; nothing else is gated.

## Conformance tests

`tests/conformance.rs`, cases opt-001..opt-006 — JSON-line verdicts,
seeded LCG randomness, fs-obs Custom event carrying the fixture
problem hash and routing refusal. Any reimplementation must pass the
suite unchanged.

The RE.Q1 game battery covers canonical schema replay and semantic-set
permutation, quantifier-order and open-loop/feedback identity separation,
clairvoyance refusal, empty open-loop information, hidden and delayed
disturbance information, strategy-memory budgets, finite/infinite horizons,
objective/stopping coherence, hybrid mode information, unresolved/excluded
Zeno scope, unit/model/domain mismatch, future access, DAE/stopping failures,
parallel-component uniqueness, sequential multiplicity, schema/cap mutation,
canonical diagnostics, and cancellation.

`tests/admission.rs` (beads sj31i.48 / xf8v7) — G0/G4/G5
leaf-policy tables (manifold boundaries incl. checked `Stiefel`
overflow, non-finite payloads with bit retention, weight policy incl.
`-0.0`, tag domains, checked dimension combining), builder-rollback
identity, builder/admission agreement, mutation-sensitive semantic
identity, deterministic bounded admission reports under explicit caps,
domain-separated identity round trips, v3 bilevel + legacy
quarantine, checked id accessors, binding validation, saturated component
diagnostics, scalar-only Min/Max, fail-closed graph-depth / aggregate
retained-byte builder rollback, exact aggregate-work receipts, measured
cap+1 storage atomicity, and the max-depth builder/parser/admission/eval
boundary. Unit fixtures additionally pin string-bearing and bitwise-float
fingerprint identity plus exact fallback inside a simulated hash-collision
bucket.

`tests/guard.rs` (Proposal D, 15 cases): no-steps→provisional-not-honored;
all-pass→cleared→honored; a veto→failed with a finding; an unregistered
step keeps the endpoint provisional (never cleared on a skipped check);
fixed step order; the amended contract needs BOTH converged and cleared;
determinism; first-registered-step-of-a-kind wins; `from_descent` bridge;
and δ-perturbation passes a smooth optimum, vetoes found-better and
sharp-crack exploits, fails closed on non-finite, and treats an empty
design as vacuously robust — plus the realistic v0 state (δ-only →
provisional).

`tests/metamorphic.rs` (bead frankensim-2uce): 256 shrinkable G3 cases (seed
`0x2ACE_0402`) apply the shared `fs-propcheck` unit-rescaling declaration
`quadratic-descent-power-of-two-units` to the live `descend_fn` quadratic step
at `2e-12` absolute-relative numeric tolerance. Every generated transform is a
nonidentity power of two, and every generated start is separated from its
target. The comparator checks coordinate and objective equivariance plus exact
discrete receipts. Existing opt-005/006 fixed pins remain unchanged.

## No-claim boundaries

- The G3 descent adopter covers coherent power-of-two rescaling of a bounded,
  one-dimensional quadratic under the toy fixed-step finite-difference
  consumer. It does not claim arbitrary-unit conditioning, general optimizer
  convergence, manifold-coordinate invariance, or dimensional correctness for
  caller-defined objectives.
- `max_total_work` is a deterministic structural admission envelope
  (retained items plus expression edges), not a wall-clock or cycle-count
  performance model. Per-field byte caps separately bound string hashing
  and decoding work.
- The cap+1 allocation fixture snapshots all builder vector capacities,
  intern entries, owned-string capacity, and accounting totals. It does
  not claim that Rust's process-global allocator performs no temporary
  bookkeeping inside error formatting.
- Gradients here are FD-through-retraction toys; exact adjoints and
  reverse-mode graph gradients are the gradient-stack bead (the
  parked draft's `graph.rs` already prototypes reverse-mode — harvest
  it there).
- PDE and stochastic nodes are REPRESENTED and validated, not
  executed; FLUX studies and UQ runners bind to them in their beads.
- Constraint semantics (kinds, repair, feasibility restoration) are
  fs-constraint's; this crate carries kind + name only.
- RE.Q1 admits game meaning but computes no winning, reachable, or viable set.
  Eligible means only that a later theorem module may attempt the requested
  polarity. Infinite-horizon results and hybrid games without retained Zeno
  exclusion remain Unknown. A regularization witness identifies a declared
  transformation but RE.Q1 does not resolve or prove distinct source/result
  lineage. Inner, outer, exact, and Unknown sets are never silently
  interchanged.
- Component, witness, strategy, interface, observation, stopping-rule, and
  cell-decomposition identities are unresolved replay references at this
  layer. Admission preserves their exact bytes and local context but does not
  prove that referenced artifacts exist or are scientifically authoritative;
  later theorem checkers must resolve and validate them before making claims.
- FrankenScript `ascent.optimize` lowering binds to this IR when the
  HELM surface lands.
- Bilevel tags reference inner problems by TYPED identity; admission
  does not verify that a referenced inner problem exists, was
  admitted, or matches its id — inner-problem storage/resolution is a
  later bead. Legacy FNV references are quarantined provenance, never
  upgraded, and confer no authority.
- `ProblemSemanticId` normalization is exactly the canonical v3 body:
  hash-consing dedupes structurally identical subexpressions, but no
  deeper equivalence (variable renaming, algebraic identities,
  objective reordering) is folded — two problems can mean the same
  mathematics and carry different semantic ids.
- Direct `Problem` construction/mutation is prevented by sealed
  (crate-private) fields — the compiler enforces the seal. A
  trybuild-style compile-fail harness pinning that property as a test
  artifact is tracked follow-up work, as is a G4 adversarial
  storm/cancellation lane for admission itself.
- `ProblemAdmission` proves the checks RAN; it is not an authenticity
  or provenance anchor (no signature, no ledger binding here — HELM
  owns that).
- `Stiefel` descent uses ambient FD directions (overcomplete but
  convergent with the QR retraction); proper tangent bases join with
  the gradient stack.
- The Goodhart guard is the POLICY ENGINE only. Three of its four steps
  (rung-k+1, cross-representation, estimator-independence) need machinery
  that does not exist yet (the fidelity-ladder registry, a live Rep Router
  re-solve, ≥2 estimator families) and are `NotPerformed` until callers
  inject them — so a v0 endpoint clears to `Provisional`, never `Cleared`,
  by design. `GuardFinding`s are PRODUCED here (L4); writing them to the
  ledger as tombstones/bug reports is HELM's job (no upward dependency).
  The endpoint-vs-random catch-rate kill measurement (G4/statistical) is a
  Gauntlet harness bead, not this crate.
