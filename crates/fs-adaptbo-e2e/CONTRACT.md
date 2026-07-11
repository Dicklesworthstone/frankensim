# CONTRACT: fs-adaptbo-e2e

AnytimeBO — Bayesian optimization with anytime-valid evidence against a
declared stall-rate null. Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-bo` (GP + Expected Improvement), `fs-eproc` (betting e-process +
Gaussian-mixture confidence sequence), `fs-evidence` (Verified/Estimated). Deps
point downward.

## Public types and semantics

- `objective(x) -> f64` — a tilted double-well polynomial on `[0, 4]`.
- `run_campaign(max_iters, delta, alpha) -> AdaptBoReport` — runs GP/EI BO,
  stops when a betting e-process on the per-step stall indicator rejects at
  `alpha`, and reports an anytime-valid confidence sequence on the best-value
  trace.

## Invariants

- CONVERGENCE: the loop finds the global well (`x ≈ 3`, value `≈ −0.45`), beating
  the shallower well.
- ANYTIME-VALID STALL DECISION: the e-process stops the search when its log
  e-value crosses `ln(1/α)` (Ville's threshold), so peeking after every
  iteration keeps type-I error `≤ α` under the declared conditional stall-rate
  null. The report carries that e-value as a statistical certificate candidate.
- The observed incumbent and GP surrogate remain `Estimated`; neither the
  e-process decision nor the running trace diagnostic is a certified enclosure
  of the global optimum.
- A shrinking mixture-boundary trace diagnostic is always reported, with no
  fixed-mean coverage claim for the adaptive incumbent sequence.
- With an impossibly small `delta`, no step stalls and the search runs to the cap.
- Deterministic (fixed grid + polynomial objective; no RNG, no libm).

## Error model

Total on the default grid; a singular GP fit (`try_fit → None`) ends the loop.
Malformed statistical controls and iteration counts above 64 panic before
work begins.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/adaptbo.rs`: BO converges and stops with an anytime-valid e-value
candidate (log-e past the Ville threshold); the incumbent stays Estimated; a
tiny delta never declares a stall; determinism; malformed/unbounded controls
refuse before work.

## No-claim boundaries

1-D grid acquisition (not a continuous inner optimizer); the stall indicator is a
binary improvement test; the confidence sequence uses a fixed sub-Gaussian σ.
The trace CS is not a global-optimum confidence interval, and rejecting the
stall null is not a convergence proof. Multi-fidelity / TuRBO acquisition are
`fs-bo`'s fuller deliverables.
