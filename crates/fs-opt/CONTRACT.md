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
  Every constructor validates: shapes (`Scalar`/`Vector(n)`), fs-qty
  DIMENSIONS (add/compare need equal dims; mul/dot add exponents; div
  subtracts; powi scales with checked refusal rather than clamping; sqrt halves even exponents; transcendentals
  demand dimensionless), node/variable existence, parameter ranges,
  and scalar-only objective/constraint roots.
- Node kinds: arithmetic (`add/sub/mul/div/neg/powi/sqrt/exp/ln/
  tanh`), vector reductions (`dot/norm_sq/component`), kinks
  (`min/max/abs` — C0), `pde_residual` (FLUX study reference with
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
  consumes; `descend_fn`/`descend_ir` are the TOY consumers proving
  iterates stay ON their manifolds.
- Structure: multi-objective (weights), constraint KINDS (`EqZero`,
  `LeZero` — semantics/repair are fs-constraint's), `ProblemTag`
  (multi-fidelity, chance-constrained, bilevel-by-hash), `EvalBudget`
  (P4, enforced by consumers).
- Serialization: `serialize` writes the canonical six-base `fsopt v2`
  line-based text form and `parse` round-trips it BITWISE (floats travel as
  bit patterns). `parse_with_version` also accepts strict explicit `fsopt v1`
  bytes emitted by either known historical v1 token writer, maps the absent
  amount exponent to `mol = 0`, and returns an immutable
  `DimensionCrosswalkReceipt`. Headerless input is refused because no
  historical writer emitted it and it has no authoritative schema identity.
  The receipt
  binds BLAKE3 hashes of the complete old artifact and complete canonical v2
  artifact under the sole `AppendMoleZero` rule; `parse` refuses v1 because it
  cannot return that mandatory evidence. `ParsedProblem` keeps its fields
  private and exposes read-only provenance accessors plus `into_parts`, so an
  inconsistent source-version/hash/receipt tuple cannot be constructed through
  the public API. Every admitted artifact has exactly one terminal `hash`
  directive;
  `problem_hash` (in-house FNV-1a 64 over the canonical body) is the
  study identity. Legacy FNV identity is verified over the historical v1 body
  without normalization into v2. The reader rebuilds an exact canonical v1
  artifact using both known historical v1 token encodings and requires complete
  byte equality with one of them before it may issue the sole-rule receipt;
  v2 receives the same complete-byte comparison against `serialize`. CRLF,
  blank/missing-final-newline forms, noncanonical token spellings, malformed
  escapes, wrong IDs, extra fields, and duplicate/missing budget directives
  therefore refuse even when their recomputed FNV is internally consistent.
  Parsing REBUILDS through the validating builder and verifies the integrity
  hash — tampered or ill-typed files refuse with line numbers.
- `eval`: memoized evaluation of algebraic subgraphs; PDE/stochastic
  nodes refuse with `Unevaluable` NAMING their executor.
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

## Invariants

1. Seeded ill-typed constructions refuse with teaching text naming
   ops/nodes, and a 600-op fuzz storm matches an independent validity
   model exactly (opt-001).
2. build→serialize→parse through canonical `fsopt v2` yields an IDENTICAL
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
   error); cancellation returns the teaching error; PDE/stochastic
   nodes name their executor when asked to evaluate (opt-006).

## Error model

`OptError` teaching errors throughout: unknown ids, shape/dimension
mismatches and dimension overflow (with exponent vectors shown), non-dimensionless
transcendentals, odd-sqrt dims, bad parameters/indices, non-scalar
roots, `NonsmoothForFamily` (node + kind + class),
`NoAdjoint` (node + study), `Unevaluable` (node + executor), `Parse`
(line + what), `Cancelled`, `BudgetExhausted` (spent count receipt).

## Determinism class

Fully deterministic: `BTreeMap` interning, index-ordered ids, bitwise
float serialization, in-house FNV hashing, no time or randomness.
Identical build sequences give identical problems, hashes, and bytes
(opt-002/003 are the trip-wires).

## Cancellation behavior

`descend_fn`/`descend_ir` poll `cx.checkpoint()` every step and
return `OptError::Cancelled` between steps. Budget exhaustion is a
RECEIPT (`budget_stopped` in the report), not an error — the iterate
remains valid (P4).

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

`tests/guard.rs` (Proposal D, 15 cases): no-steps→provisional-not-honored;
all-pass→cleared→honored; a veto→failed with a finding; an unregistered
step keeps the endpoint provisional (never cleared on a skipped check);
fixed step order; the amended contract needs BOTH converged and cleared;
determinism; first-registered-step-of-a-kind wins; `from_descent` bridge;
and δ-perturbation passes a smooth optimum, vetoes found-better and
sharp-crack exploits, fails closed on non-finite, and treats an empty
design as vacuously robust — plus the realistic v0 state (δ-only →
provisional).

## No-claim boundaries

- Gradients here are FD-through-retraction toys; exact adjoints and
  reverse-mode graph gradients are the gradient-stack bead (the
  parked draft's `graph.rs` already prototypes reverse-mode — harvest
  it there).
- PDE and stochastic nodes are REPRESENTED and validated, not
  executed; FLUX studies and UQ runners bind to them in their beads.
- Constraint semantics (kinds, repair, feasibility restoration) are
  fs-constraint's; this crate carries kind + name only.
- FrankenScript `ascent.optimize` lowering binds to this IR when the
  HELM surface lands.
- Bilevel tags reference inner problems by hash; inner-problem
  storage/resolution is a later bead.
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
