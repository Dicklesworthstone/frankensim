# CONTRACT: fs-adaptbo-e2e

AnytimeBO — Bayesian optimization that provably knows when to stop. Layer L4
(ASCENT).

## Purpose and layer

Composes `fs-bo` (GP + Expected Improvement), `fs-eproc` (betting e-process +
Gaussian-mixture confidence sequence), `fs-evidence` (Verified/Estimated). Deps
point downward.

## Public types and semantics

- `objective(x) -> f64` — a tilted double-well polynomial on `[0, 4]`.
- `run_campaign(max_iters, delta, alpha) -> AdaptBoReport` — runs GP/EI BO,
  stops when a betting e-process on the per-step stall indicator rejects at
  `alpha`, and reports the anytime-valid confidence sequence on the optimum.

## Invariants

- CONVERGENCE: the loop finds the global well (`x ≈ 3`, value `≈ −0.45`), beating
  the shallower well.
- ANYTIME-VALID STOP: the e-process stops the search when its log e-value crosses
  `ln(1/α)` (Ville's threshold), so peeking after every iteration keeps the
  false-stop rate `≤ α`. `stopped_early` ⇒ `Verified` stop color; the GP
  surrogate is `Estimated`.
- A shrinking anytime-valid interval on the optimum is always reported.
- With an impossibly small `delta`, no step stalls and the search runs to the cap.
- Deterministic (fixed grid + polynomial objective; no RNG, no libm).

## Error model

Total on the default grid; a singular GP fit (`try_fit → None`) ends the loop.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/adaptbo.rs` (3): BO converges and stops with an anytime-valid certificate
(log-e past the Ville threshold); a tiny delta never declares a stall;
determinism.

## No-claim boundaries

1-D grid acquisition (not a continuous inner optimizer); the stall indicator is a
binary improvement test; the confidence sequence uses a fixed sub-Gaussian σ.
Multi-fidelity / TuRBO acquisition are `fs-bo`'s fuller deliverables.
