# CONTRACT: fs-diffreal-e2e

The differentiation & reality end-to-end suite (plan addendum, Proposal 11 /
Layer-3 conformance): a runnable battery over selected differentiation,
synthetic as-built, and tolerance fixtures.

## Purpose and layer

Layer L6. This integration crate depends on `fs-ad`, `fs-adjoint` (with
`ledger-transpose`), `fs-asbuilt`, `fs-assimilate`, `fs-blake3`, `fs-evidence`,
`fs-exec`, and `fs-toleralloc`. It composes those crates and owns the local
scalar fixture VJPs, an operator-name-bounded registry facade, the
sensitivity-seal schema, typed stage events, and the battery's fail-closed
report policy. Registry cardinality and declaration diagnostics are not
bounded. The crate does not own a general differentiation or numerical
primitive.

## Public types and semantics

- `run_battery(&Cx) -> Result<DiffRealReport, DiffRealError>` runs all four
  stages. Cancellation or ambient-budget refusal returns a typed error and
  publishes no partial report.
- `DiffRealReport { stages }` — `complete()`, `all_required_passed()`,
  `promotion_ready()`, fail-closed compatibility alias `passed()`, and
  `stage(name)`.
- `StageLog { stage, requirement, status, evidence_identity, events }` is a
  diagnostic record. Its `StageEvent` values are typed deterministic DATA;
  floating-point payloads are stored as exact bits and are never printed here.
- `StageRequirement::{Required, Optional}` — makes the decision policy
  explicit. The four current battery stages are all required.
- `StageStatus::{Passed, Failed(StageReason), Gated(StageReason),
  Refused(StageReason)}` — distinguishes an evaluated failure from work that
  was not validly evaluated.
- `StageReason { code, detail }` — a stable programmatic code plus deterministic
  human-readable detail. `Display` for requirements, statuses, reasons, and
  logs is deterministic and preserves this distinction.
- `DifferentiationRegistry` bounds operator names used through its facade over
  the shared `VjpRegistry`; `production_vjp_registry()` constructs the fixed
  fixture.
- `differentiate_path(ops, registry, x, cx) -> Result<PathDerivative,
  DifferentiationError>` records the admitted scalar path on the shared
  `fs-adjoint` tape and executes its registered reverse sweep. The first
  missing VJP blocks the gradient; no silent zero is substituted.
- `verify_sensitivity(ops, registry, x, cx) -> Result<SealedSensitivity,
  DifferentiationError>` is restricted to the canonical fixture path and mints
  opaque sensitivity evidence only after reverse, dual, and two-step FD
  agreement. The receipt supports integrity recomputation and checked scalar
  input-unit rescaling.
- `stage_differentiation_with_registry` is the injected-VJP falsifier seam.
- `stage_tolerance_allocation_with_samples` is the injected-sample
  linearization falsifier seam; it is not a probability API.
- All four normal stage entry points take `&Cx` and return
  `Result<StageLog, DiffRealError>`.
- `DifferentiationError` and `DiffRealError` preserve typed structural,
  representability, cancellation, and budget refusals.

## The four stages (each a fail-closed assertion)

1. **Differentiation** — the stage records the local
   affine→square→identity scalar fixture on `fs_adjoint::transpose::Tape`,
   executes the registered VJP transpose, and seals the result only when it
   agrees with the independent `fs-ad` dual result and two-step
   `fd_falsifier`. A perturbed VJP produces `Failed`. Missing VJP coverage is
   checked before numerical evaluation, so it cannot be hidden by NaN or
   infinity. Non-finite inputs, forward values, gradients, and rescaled values
   are typed refusals.
2. **As-built loop** — a synthetic three-point rigid transform registers
   (residual carried forward), the as-built delta is an Estimated candidate
   carrying a calibration-label candidate, a seeded defect is LOCALIZED
   (argmax deviation), and point-sensor assimilation reduces model-data
   misfit. No scan ingestion, custody, measurement covariance, metrology
   validation, or calibration authority is inferred from the label.
3. **Tolerance allocation** — both feature sensitivities come from
   `verify_sensitivity`; `ColorRank::Verified` is assigned only after receipt
   agreement and integrity recomputation. Allocation direction and GD&T
   carriage are checked separately. `robustness_check` evaluates only the
   supplied scalar samples against the first-order bound, and its typed event
   always records `probability_claimed=false`.
4. **Gated spacetime** — `fs-time` temporal-complex support and its owning bead
   are shipped, but the coupled spacetime fixture is not integrated and
   activated in this battery. This required stage is honestly `Gated`, not
   silently passed.

## Status, completeness, and promotion policy

| Required-stage status | Assertion evaluated? | Can report be complete? | Can promote? |
| --- | --- | --- | --- |
| `Passed` | yes | yes | yes, if every required stage passed |
| `Failed(reason)` | yes | yes | no |
| `Gated(reason)` | no | no | no |
| `Refused(reason)` | no | no | no |

Completeness is schema-aware, not a vacuous `all()` over whatever records a
caller supplied. The four fixed required stage names must each appear exactly
once in fixed relative order, be marked `Required`, carry their exact versioned
fixture/evidence identity, contain non-empty diagnostic events, and have an
evaluated status. Missing, reordered, or duplicate stages, unexpected required
stages, blank diagnostics, identity drift, gates, and refusals all fail closed.

`all_required_passed()` separately requires the same valid schema and a
`Passed` status for every required stage. `promotion_ready()` is the conjunction
of completeness and all-required-passed. `passed()` is retained as a
compatibility alias for `promotion_ready()`; it does not restore the old
boolean semantics.

Additional stages must declare `Optional`. A well-formed optional gated,
refused, or failed diagnostic does not block the decision over the fixed
required set. An optional record cannot replace, rename, or downgrade one of
the four required stages.

## Invariants

- The current full battery is DETERMINISTIC for equal `Cx` provenance and
  inputs, but is intentionally **not complete or promotion-ready** while the
  required spacetime integration stage is gated.
- Differentiation paths contain 1–16 operators; each name is at most 64 bytes.
  The first missing VJP wins in forward order and is checked before non-finite
  input rejection.
- Published primal, gradient, oracle, and rescaled values are finite.
- `SealedSensitivity` binds the policy version, length-framed operator names,
  input, primal, production gradient, dual gradient, both FD values, and the FD
  tolerance with a domain-separated content identity.
- Scientific derivative disagreement produces `StageStatus::Failed`; it can
  never mint a sealed sensitivity or a Verified tolerance input.
- The as-built defect is localized to the seeded index. The tolerance stage
  derives numeric sensitivity and `Verified` rank from integrity-checked
  receipts, tightens high sensitivity, loosens low sensitivity, and requires at
  least one loosened GD&T row carrying those derived fields. GD&T rows do not
  retain the receipt or its identity.
- Crate-produced tolerance-stage sampled-linearization events always carry
  `probability_claimed=false`; caller-authored diagnostics carry no authority.
- Differentiation fixture identity v2 and tolerance fixture identity v3 are not
  interchangeable with their earlier schemas.
- No required `Gated` or `Refused` stage can make `complete()`,
  `all_required_passed()`, `promotion_ready()`, or `passed()` return true.
- A required `Failed` stage is an evaluated result, so it may be complete, but
  it can never be all-required-passed or promotion-ready.

## Error model

Scientific assertion failures are `StageStatus::Failed` with structured reason
codes and retained typed events. Deliberate unavailability is `Gated`; an
inability to evaluate admissibly is `Refused`. Fixed-fixture allocation or
sample-check refusals become a refused stage instead of panicking. Production,
dual, or FD disagreement becomes a failed stage. Cancellation, insufficient
ambient stage quota, and lower-layer as-built/assimilation errors remain typed
`DiffRealError` values and suppress the partial battery report.

## Determinism class

The fixed crate-authored battery and `production_vjp_registry` are fully
deterministic for equal inputs and `Cx` provenance: no RNG and no I/O. Stage
order, exact-bit numeric event fields, status/reason codes, versioned fixture
identities, sealed content identities, and `Display` output are stable.
Caller-injected `Vjp` implementations and caller-authored events are outside
this determinism claim.

## Cancellation behavior

Every stage accepts `&Cx`, checks cancellation before fixed-work admission, and
has a nonzero `cost_quota` threshold. Differentiation additionally checks at
bounded forward-operator boundaries, around transpose/oracle work, and before
publication. The as-built stage passes `Cx` through its bounded lower-layer
operations. Tolerance checks between sealed sensitivity, allocation, reporting,
and sampled-linearization phases. Spacetime checks before recording its gate.

Cancellation or insufficient ambient quota suppresses the partial battery
report. These fixed threshold checks do not consume quota and are not a
monotonic, non-reissuable, whole-battery budget ledger.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

This crate exposes no feature flags. It enables `fs-adjoint/ledger-transpose`
on its dependency to exercise the shared tape/VJP path.

## Conformance tests

`tests/e2e.rs` covers the shared tape/VJP result, a perturbed-VJP kill test,
dual-number gradient agreement, the coarse/fine FD study, sealed-receipt replay/integrity,
valid/invalid unit rescaling, NaN, both infinities, finite-input intermediate
overflow, a non-finite VJP result, missing-VJP precedence, typed as-built events, allocation direction,
adverse supplied samples, explicit probability no-claim, the spacetime gate,
one exact sampled-linearization event display golden, status/log displays,
deterministic replay, pre- and injected mid-stage cancellation, and zero-cost
admission for all four stages and the full battery.

The private-field unit test in `src/lib.rs` mutates one bound sensitivity field
and proves integrity recomputation rejects it. Report-policy unit tests cover
all-passed, failed, gated, refused, optional, missing, duplicated,
misidentified, reordered, malformed, and unlogged schema cases.

## No-claim boundaries

- "Production" here means the shared `fs-adjoint` `Tape`/`VjpRegistry`
  implementation. The names `sdf`, `spline`, and `solve` denote crate-local
  affine, square, and identity scalar fixtures. This battery does not execute
  production geometry, meshing, spline, or solver kernels and does not certify
  an SDF→mesh/spline→solve derivative.
- The independent dual and two-step finite-difference computations falsify the
  local reverse sweep only at selected scalar inputs. They are not external
  validation, a global derivative bound, a convergence-order certificate, or
  proof over an operating domain. `verify_sensitivity` is fixture-specific;
  its oracles encode `(2x + 1)²` and reject every other path.
- `SealedSensitivity::identity` is an unkeyed content-integrity binding, not a
  signature, authenticated provenance, ledger receipt, authorization, or
  independent scientific certificate. Battery-local `ColorRank::Verified`
  means only that this fixed scalar fixture passed the declared reverse/dual/FD
  policy.
- `gradient_in_input_units` applies a caller-supplied positive scalar through
  the chain rule. It does not validate dimensions or unit provenance.
- The sampled-linearization check compares only caller-supplied scalar QoI
  values with a first-order bound. It neither proves those values came from
  actual tolerance corners nor enumerates corners, characterizes dependence or
  tails, or establishes `P(performance ∈ spec)`. `variance_budget` expresses a
  probability target only under its normal/first-order model assumptions.
- GD&T rows carry sensitivity and color, but not the full metrology domain,
  residual uncertainty, custody, or manufacturing-process evidence.
- Fixed cost thresholds are not consumed, monotonic, or non-reissuable budget
  accounting. Public raw `Vjp` implementations remain subject to the shared
  registry's arity contract; this crate does not turn an arbitrary panicking
  implementation into a typed refusal.
- Stage 4 is GATED because this battery has no integrated, activated coupled
  spacetime fixture. This is not a claim that `fs-time` or its temporal-complex
  bead is unbuilt.
- `StageLog::evidence_identity` is a versioned fixture/schema label, not the
  `SealedSensitivity` hash or external proof. `StageEvent` and `StageLog` are
  returned data, not persisted ledger events, and caller-constructible
  diagnostics carry no promotion authority.
