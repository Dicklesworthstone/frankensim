# CONTRACT: fs-diffreal-e2e

The differentiation & reality end-to-end suite (plan addendum, Proposal 11 /
Layer-3 conformance): a runnable battery exercising adjoints + reality-as-a-chart.

## Purpose and layer

Layer L6. An integration crate: depends on `fs-evidence`, `fs-asbuilt`,
`fs-assimilate`, `fs-exec`, and `fs-toleralloc`. It composes them; it owns no
new numerical primitive. It does own the battery's fail-closed report policy.

## Public types and semantics

- `run_battery(&Cx) -> Result<DiffRealReport, DiffRealError>` — runs all four
  stages; lower-layer cancellation or scientific refusal returns a typed error
  and publishes no partial report.
- `DiffRealReport { stages }` — `complete()`, `all_required_passed()`,
  `promotion_ready()`, fail-closed compatibility alias `passed()`, and
  `stage(name)`.
- `StageLog { stage, requirement, status, evidence_identity, events }` — the
  structured per-stage log. Events are DATA, never printed.
- `StageRequirement::{Required, Optional}` — makes the decision policy
  explicit. The four current battery stages are all required.
- `StageStatus::{Passed, Failed(StageReason), Gated(StageReason),
  Refused(StageReason)}` — distinguishes an evaluated failure from work that
  was not validly evaluated.
- `StageReason { code, detail }` — a stable programmatic code plus deterministic
  human-readable detail. `Display` for requirements, statuses, reasons, and
  logs is deterministic and preserves this distinction.
- Per-stage entry points (`stage_differentiation`, `stage_as_built_loop`,
  `stage_tolerance_allocation`, `stage_spacetime_gated`).
- `differentiate_path(ops, has_vjp, x)` — a gradient over an op pipeline that
  returns `Err` naming the first op with a missing VJP (blocks; never a silent
  zero).

## The four stages (each a fail-closed assertion)

1. **Differentiation** — the adjoint (reverse-mode) gradient agrees with finite
   differences within tolerance; a full-VJP-coverage path differentiates; a
   forced-remesh path with a missing VJP is BLOCKED (structured error).
2. **As-built loop** — a scanned fixture registers (residual carried forward),
   the as-built delta is an Estimated candidate carrying calibration provenance,
   a seeded defect is LOCALIZED (argmax deviation), and registration-free
   point-sensor assimilation reduces the model-data misfit. No calibration
   authority is inferred from a caller-supplied string.
3. **Tolerance allocation** — the high-sensitivity feature is tightened, the low
   one loosened, every loosened tolerance is justified by a certified
   sensitivity, and the band-extremes check confirms `P(in-spec)`.
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
- The differentiation, as-built, and tolerance stages currently return
  `Passed`; the spacetime stage returns `Gated` with a stable reason code and
  versioned evidence identity.
- A missing VJP blocks the gradient; the as-built defect is localized to the
  seeded index; the tolerance allocation tightens-high / loosens-low.
- No required `Gated` or `Refused` stage can make `complete()`,
  `all_required_passed()`, `promotion_ready()`, or `passed()` return true.
- A required `Failed` stage is an evaluated result, so it may be complete, but
  it can never be all-required-passed or promotion-ready.

## Error model

Scientific assertion failures are `StageStatus::Failed` with structured reason
codes and retained events. Deliberate unavailability is `Gated`; an inability
to evaluate admissibly is `Refused`. The fixed tolerance fixture converts an
unexpected budget/allocation refusal into a `Refused` stage instead of
panicking. Lower-layer as-built/assimilation errors, including cancellation,
remain typed `DiffRealError` values and suppress the partial battery report.

## Determinism class

Fully deterministic for equal inputs and `Cx` provenance: every subsystem it
drives is deterministic; no RNG and no I/O. Stage order, status codes, reason
codes/details, versioned identities, events, and `Display` output are stable.

## Cancellation behavior

The as-built/assimilation stage accepts an explicit `Cx` and polls through its
lower-layer operations. A cancellation becomes a typed `DiffRealError`; no
partial `DiffRealReport` is returned. The analytic differentiation and fixed
tolerance stages are bounded synchronous fixture work and currently do not poll
independently.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/e2e.rs` (Layer-3 conformance): the current required gate makes the real
battery incomplete and non-promotable; differentiation agrees with FD and
blocks a missing VJP; the as-built loop localizes a defect and reduces misfit;
tolerance allocation tightens-high / loosens-low and confirms the sampled
robustness fixture; the spacetime stage is honestly gated; synthetic fixed-
schema reports exercise all-passed, failed, gated, and refused truth-table
rows; optional-stage policy is explicit; missing/duplicate/misidentified or
reordered/unlogged stages fail closed; status/log display is deterministic and
distinguishes failure from unavailability; replay is deterministic; and
cancellation propagates without a partial report.

## No-claim boundaries

- Stage 1 uses a SELF-CONTAINED analytic adjoint + finite-difference check and a
  VJP-coverage gate to demonstrate adjoint-vs-FD agreement and missing-VJP
  blocking; the production seam-crossing gradient (SDF→mesh→solve) runs on
  fs-adjoint's certified adjoints.
- Stage 4 is GATED because this battery has no integrated, activated coupled
  spacetime fixture. This is not a claim that `fs-time` or its temporal-complex
  bead is unbuilt.
- `evidence_identity` is a versioned binding to the stage fixture/schema. It is
  not a content hash, certificate, independent verification receipt, or claim
  that the stage's scientific result is externally validated.
- The suite emits log events as returned DATA; wiring them to structured
  tracing / ledger sinks is the harness integration. Event payloads remain
  human-readable strings; typed event schemas are outside this bead.
