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
report policy. It also owns canonical fixed-fixture transcripts, versioned
stage/result receipts, the as-built stage's fixed invocation plan, the ordered
report root, and the DiffReal-specific external-authority seam. The invocation
plan composes the typed affine accounting implemented by `fs-exec`; this crate
does not mint a second spend authority. Registry cardinality and declaration
diagnostics are not bounded. The crate does not own a general differentiation
or numerical primitive.

`fs-report` and `fs-scenario` are DEV-dependencies only, added for the as-built
product-chain battery (bead `frankensim-extreal-program-f85xj.12.2`). The
library composes neither, and neither appears in this crate's runtime graph.

## Public types and semantics

- `run_battery(&Cx) -> Result<DiffRealReport, DiffRealError>` runs all four
  stages. Cancellation or ambient-budget refusal returns a typed error and
  publishes no partial report.
- `run_battery_with_clock(&Cx, &dyn TimeSource)` runs the same sealed path with
  an injected monotonic clock. It exists for deterministic absolute-deadline
  testing; the clock is not a caller-controlled scientific input.
- `DiffRealReport` exposes read-only diagnostics, stage receipts, the ordered
  report root, exact execution identity, `complete()`,
  `all_required_passed()`, `structurally_ready()`, integrity/replay checks, and
  `authenticate(...)`. It exposes no raw promotion predicate.
- `StageReceipt` is opaque crate-minted DATA. Its root binds the schema version,
  exact fixed-fixture input root, complete private result root, canonical typed
  log fields, full `Cx` stream/mode/budget identity, and policy versions. The
  as-built receipt additionally retains and binds the complete, integrity-
  checked `fs_exec::InvocationReceipt`; every other fixed stage must carry no
  invocation receipt under the current schema.
- `PromotionReceiptVerifier` authenticates one immutable, domain-separated
  ordered report root. `NoPromotionReceiptVerifier` is the default deny-all
  authority. Attestation, request, and atomic decision must carry the exact
  same policy fingerprint.
- `AuthenticatedDiffRealReport` is minted only after version, integrity,
  expected-`Cx`, semantic, attestation-policy, verifier-policy, and authority
  checks pass. It is the only type exposing `promotion_ready()`.
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
  linearization falsifier seam; it is not a probability API. Neither injected
  seam can mint a stage receipt or enter a sealed report.
- All four normal stage entry points take `&Cx` and return
  `Result<StageLog, DiffRealError>`. `stage_as_built_loop_with_clock` is the
  deterministic-clock form of the same fixed as-built path; neither diagnostic
  entry point can mint a `StageReceipt` outside `run_battery[_with_clock]`.
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
2. **As-built loop** — one preflight admits the complete synthetic transaction
   before scientific work starts. One non-cloneable root transfers affine
   child leases to setup, registration, comparison, prior construction,
   colored assimilation, and publication. Work, polls, cost, evaluations,
   concurrent memory, and retained output remain distinct typed dimensions;
   immutable accuracy and capability identities travel in the root envelope.
   Preflight derives assimilation work from the actual execution mode. Its
   292-poll/64-KiB/16-KiB child grant remains admission headroom; the accepted
   fixed receipt must carry the lower planner's exact memory/output usage and
   the fixture's exact 46-poll mixed-stride spend.
   The three-point rigid transform registers (residual carried forward), the
   as-built delta is an Estimated candidate carrying a calibration-label
   candidate, a seeded defect is LOCALIZED using the retained deterministic
   argmax, and point-sensor assimilation reduces model-data misfit. The stage
   uses the `misfit_before`/`misfit_after` values retained by that single
   assimilation result; it does not issue two redundant standalone misfit
   evaluations. No scan ingestion, custody, measurement covariance, metrology
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
`Passed` status for every required stage. `structurally_ready()` is the
unauthenticated conjunction of completeness and all-required-passed; its name
is an explicit no-authority boundary. A report becomes promotion-eligible only
after `authenticate` validates every receipt and semantic transcript against an
exact replay `Cx`, then one injected authority authenticates the ordered report
root through a domain-separated subject that also binds the purpose, full
execution identity, both receipt versions, and exact local policy fingerprint.
Only the returned opaque wrapper has `promotion_ready()`.

Stage roots are domain-separated, canonical, length-framed encodings with
explicit numeric tags for every enum variant and exact IEEE-754 bits for every
floating value. They do not use `Display`, `Debug`, JSON, sorting, or an
unobserved ambient wall clock as identity. A deadline-bearing as-built receipt
binds the injected monotonic clock observation used to enforce that absolute
deadline; a deadline-free invocation performs and binds no clock observation.
The ordered report root binds stage name/root pairs without reordering;
omission, duplication, or reordering changes or invalidates it. Unknown
stage/report/invocation receipt versions fail before authority verification.

Additional stages must declare `Optional`. A well-formed optional gated,
refused, or failed diagnostic does not block the decision over the fixed
required set. An optional record cannot replace, rename, or downgrade one of
the four required stages.

## Invariants

- The current full battery is DETERMINISTIC for equal `Cx` provenance, inputs,
  and injected logical-clock observations. Deadline-free runs do not observe a
  clock. It may be externally authenticated as an exact run, but is
  intentionally **not complete or promotion-ready** while the required
  spacetime integration stage is gated.
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
- The as-built transaction has exactly one root admission and one ordered child
  topology. Capacity cannot be cloned or reissued between nested scientific
  stages. Admission checks the complete fixed plan before work; one-below poll
  or cost envelopes refuse before a stage log exists. Unused capacity returns
  exactly once, while consumed work/cost/evaluations and retained output do not.
- The as-built invocation receipt is acceptable only with `Completed`
  disposition, no latched failure, valid typed conservation/topology/memory
  semantics, and a root that recomputes. Every leaf's direct work, polls, cost,
  evaluations, memory request/release/peak, and retained output must equal the
  fixed implementation plan; the transaction/root aggregates must equal those
  leaves plus setup. A cap-only receipt with omitted polls, reservations, or
  publication is rejected even if its root is recomputed. All admitted
  temporary memory is released before sealing. The completed receipt is part of
  the as-built `StageReceipt` root and therefore of the ordered report root.
- The fixed as-built plan binds its accuracy obligation, capability scope,
  absolute deadline, and full `Cx` identity. Accuracy or capability obligations
  cannot be weakened by a nested stage. Work, poll, cost, evaluation, memory,
  and output quantities are never converted into one another.
- Crate-produced tolerance-stage sampled-linearization events always carry
  `probability_claimed=false`; caller-authored diagnostics carry no authority.
- Differentiation fixture identity v2, as-built fixture identity v2, and
  tolerance fixture identity v3 are not interchangeable with earlier schemas.
- Every fixed stage receipt binds the logical stream (`seed`, `kernel_id`,
  `tile`, `iteration`), execution mode, deadline presence/value, poll quota,
  cost-quota presence/value, priority, exact fixed inputs, complete private
  results, diagnostic semantics, and all named policy versions. Any independent
  mutation fails integrity or moves the root.
- One report-level authority decision authenticates the ordered stage-root
  sequence atomically. A stage hash alone is never described as authenticated.
- No required `Gated` or `Refused` stage can make `complete()`,
  `all_required_passed()`, `structurally_ready()`, or authenticated
  `promotion_ready()` return true.
- A required `Failed` stage is an evaluated result, so it may be complete, but
  it can never be all-required-passed or promotion-ready.

## Error model

Scientific assertion failures are `StageStatus::Failed` with structured reason
codes and retained typed events. Deliberate unavailability is `Gated`; an
inability to evaluate admissibly is `Refused`. Fixed-fixture allocation or
sample-check refusals become a refused stage instead of panicking. Production,
dual, or FD disagreement becomes a failed stage. As-built preflight and runtime
admission failures remain typed `DiffRealError::Invocation(InvocationError)`.
After root admission, every non-completed terminal path returns
`InvocationDidNotComplete` with the complete immutable failure receipt, its
disposition and redundant root, plus the original typed cause when one exists.
Cancellation, expired deadline, insufficient typed capacity, and lower-layer
as-built/assimilation errors suppress the partial stage log and battery report.

## Determinism class

The fixed crate-authored battery and `production_vjp_registry` are fully
deterministic for equal inputs, `Cx` provenance, and injected logical-clock
observations: no RNG and no I/O. Deadline-free execution does not consult the
clock, so its invocation receipt replays bit-for-bit without wall-time input.
Stage order, child issue order and IDs, typed resource consumption, exact-bit
numeric event fields, status/reason codes, versioned fixture identities,
complete result roots, stage roots, ordered report roots, sealed content
identities, and `Display` output are stable. Caller-injected `Vjp`
implementations and caller-authored events are outside this determinism claim.

## Cancellation behavior

Every stage accepts `&Cx`, checks cancellation before fixed-work admission, and
has a nonzero `cost_quota` threshold. Differentiation additionally checks at
bounded forward-operator boundaries, around transpose/oracle work, and before
publication. The as-built stage admits one invocation ledger before work, then
polls the injected absolute deadline before spending each typed poll and before
observing cancellation. Its lower-layer registration, comparison, belief, and
assimilation entry points spend the same affine child authority. A deadline hit
requests cancellation, and every failing path drains child leases and memory
reservations before terminal handling. Tolerance checks between sealed
sensitivity, allocation, reporting, and sampled-linearization phases. Spacetime
checks before recording its gate.

Cancellation, expired deadline, insufficient ambient quota, or any affine
accounting refusal suppresses the partial as-built `StageLog` and battery
report. A post-admission failure returns its drained receipt as error evidence;
only a completed, integrity-valid invocation receipt may be attached to the
as-built stage receipt.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

This crate exposes no feature flags. It enables `fs-adjoint/ledger-transpose`
on its dependency to exercise the shared tape/VJP path.

## Conformance tests

`tests/as_built_chain.rs` covers the as-built product chain end to end —
register (`fs-asbuilt`) → bind (`fs-scenario`) → resolve → propagate + budget
(`fs-asbuilt`/`fs-evidence`) → render (`fs-report`) — asserting that the
scenario cites exactly the registration it was bound to, that a foreign
citation does not resolve to it, that the two gap QoIs sharing the pose
translation are genuinely correlated, that every geometry term cites the one
shared propagation record, that the render reaches the same shared-record
conclusion and names that record, that an unchanged QoI still renders, and that
per-stage forensic logging retains every stage plus the correlation structure.
A second case pins determinism: two independent runs produce the same
registration identity and the same propagation record, and the binding
transports that identity unchanged.

`tests/e2e.rs` covers the shared tape/VJP result, a perturbed-VJP kill test,
dual-number gradient agreement, the coarse/fine FD study, sealed-receipt
replay/integrity, valid/invalid unit rescaling, NaN, both infinities,
finite-input intermediate overflow, a non-finite VJP result, missing-VJP
precedence, typed as-built events, allocation direction, adverse supplied
samples, explicit probability no-claim, the spacetime gate, one exact
sampled-linearization event display golden, status/log displays, deterministic
replay, pre- and injected mid-stage cancellation, and zero-cost admission for
all four stages and the full battery. The as-built invocation battery adds G0
completed topology/conservation/memory checks, exact-envelope success, G3
one-below poll and cost refusals, exact 46/69 assimilation/transaction poll
spend, mode-aware Fast/Deterministic work plans, and receipt-root binding; G4
cancellation and absolute-deadline suppression; and G5 bit-for-bit replay under
an injected clock. Its exact phase topology also locks out redundant standalone
misfit child work. Other G3 receipt batteries cover independent input/result/event/root
mutations, omission, duplication, reordering, wrong signatures,
unknown/revoked keys, policy disagreement, and unknown versions. Other G5
batteries cover independent mismatch of every `Cx` stream, mode, deadline,
poll, cost, and priority field.

The private-field unit test in `src/lib.rs` mutates one bound sensitivity field
and proves integrity recomputation rejects it. Report-policy unit tests cover
all-passed, failed, gated, refused, optional, missing, duplicated,
misidentified, reordered, malformed, and unlogged schema cases. Private result
transcripts bind values omitted from diagnostic events: registration transform
and all deviations, estimator/color/regime, posterior mean/covariance/color and
misfits, full allocation/GD&T rows and totals, and the robustness verdict.

## No-claim boundaries

- `tests/as_built_chain.rs` proves that ONE identity survives the five-stage
  as-built product chain: the registration's `model_identity()` is what the
  scenario cites, and the propagation record's identity is what every geometry
  term and the rendered correlation line agree on. Its "solve" is a LINEAR QoI
  evaluator, not a physics solve. The battery therefore establishes that the
  seams compose and the identities flow; it does not establish that any thermal
  or structural result is correct, does not exercise a real solver consuming
  `PlacementBasis::AsBuilt`, and does not validate the registration against any
  measured artifact. The fiducial scan is synthetic with fixed literal
  perturbations.
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
- The as-built invocation receipt proves local affine accounting integrity, not
  optimal resource use, scheduler fairness, physical energy cost, or scientific
  validity. Its conservative poll, memory, and output capacities are fixed
  policy envelopes; unused returned capacity is not evidence that the policy is
  globally tight. Public raw `Vjp` implementations remain subject to the shared
  registry's arity contract; this crate does not turn an arbitrary panicking
  implementation into a typed refusal.
- Stage 4 is GATED because this battery has no integrated, activated coupled
  spacetime fixture. This is not a claim that `fs-time` or its temporal-complex
  bead is unbuilt.
- `StageLog::evidence_identity` is a versioned fixture/schema label, not the
  `SealedSensitivity` hash or external proof. `StageEvent` and `StageLog` are
  returned data, not persisted ledger events, and caller-constructible
  diagnostics carry no promotion authority.
- `StageReceipt` and `DiffRealReport::receipt_root` are unkeyed
  content-integrity bindings, not signatures or authority. Authentication
  begins only when an injected verifier accepts a detached attestation over the
  exact ordered root and policy fingerprint.
- Successful authentication proves only that the injected authority approved
  that exact root under its declared policy. It does not independently prove
  scientific correctness, trusted-hardware execution, external validation,
  calibration custody, ledger persistence, release admission, or permission to
  strengthen any evidence color.
