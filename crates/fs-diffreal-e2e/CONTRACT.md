# CONTRACT: fs-diffreal-e2e

The differentiation & reality end-to-end suite (plan addendum, Proposal 11 /
Layer-3 conformance): a runnable battery exercising adjoints + reality-as-a-chart.

## Purpose and layer

Layer L6. An integration crate: depends on `fs-evidence`, `fs-asbuilt`,
`fs-assimilate`, `fs-toleralloc`. It composes them; it owns no new primitive.

## Public types and semantics

- `run_battery() -> DiffRealReport` — runs all four stages.
- `DiffRealReport { stages }` — `passed()`, `stage(name)`.
- `StageLog { stage, passed, events }` — the structured per-stage log (events
  are DATA, never printed).
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
   the as-built δ is validated-color + anchored to a calibration cert, a seeded
   defect is LOCALIZED (argmax deviation), and registration-free point-sensor
   assimilation reduces the model-data misfit.
3. **Tolerance allocation** — the high-sensitivity feature is tightened, the low
   one loosened, every loosened tolerance is justified by a certified
   sensitivity, and the band-extremes check confirms `P(in-spec)`.
4. **Gated spacetime** — the temporal-complex stage is honestly reported GATED
   (its bead is unbuilt), not silently passed.

## Invariants

- The full battery passes and is DETERMINISTIC (`run_battery() ==
  run_battery()`).
- A missing VJP blocks the gradient; the as-built defect is localized to the
  seeded index; the tolerance allocation tightens-high / loosens-low.

## Error model

No panics; a stage records `passed = false` with its events on any failure.

## Determinism class

Fully deterministic: every subsystem it drives is deterministic; no RNG, no I/O.

## Cancellation behavior

None (synchronous).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/e2e.rs` (Layer-3 conformance, 6 cases): the full battery passes with all
four stages logged; differentiation agrees with FD + blocks a missing VJP; the
as-built loop localizes a defect + reduces misfit; tolerance allocation
tightens-high / loosens-low + confirms robustness; the spacetime stage is
honestly gated; determinism.

## No-claim boundaries

- Stage 1 uses a SELF-CONTAINED analytic adjoint + finite-difference check and a
  VJP-coverage gate to demonstrate adjoint-vs-FD agreement and missing-VJP
  blocking; the production seam-crossing gradient (SDF→mesh→solve) runs on
  fs-adjoint's certified adjoints.
- Stage 4 (the spacetime complex) is GATED — its bead is unbuilt; the stage is
  reported as skipped, not asserted.
- The suite emits log events as returned DATA; wiring them to structured
  tracing / ledger sinks is the harness integration.
