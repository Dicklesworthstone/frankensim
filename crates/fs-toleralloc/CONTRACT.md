# CONTRACT: fs-toleralloc

Adjoint-driven tolerance allocation (plan addendum, Proposal 11's commercial
kicker): spend tight tolerances only where sensitivity is large; loosen the
rest with a certified justification.

## Purpose and layer

Layer L4 (optimization). Depends only on `fs-evidence` (`ColorRank` for the
sensitivity's color). Pure, deterministic.

## Public types and semantics

- `Feature { name, sensitivity, sensitivity_color, cost_coeff, baseline_tolerance }`.
- `allocate(&[Feature], variance_budget, k) -> Result<Allocation, ToleranceError>`
  — cost-optimal tolerances `tᵢ ∝ (cᵢ / sᵢ²)^{1/3}`, normalized so the QoI
  variance `Σ sᵢ²(tᵢ/k)²` exactly meets the budget. Each `TolItem` records the
  tolerance, its certified sensitivity + color, and an `Action` (Tighten /
  Loosen / Unchanged vs baseline).
- `robustness_check(&Allocation, extreme_qois, nominal_qoi, k, margin) ->
  RobustnessVerdict` — compares the first-order `linearized_std` against the QoI
  at sampled tolerance-band extremes; `confirmed` iff the extremes stay within
  `k · linearized_std · (1 + margin)`.
- `gdt_report(&Allocation) -> Vec<Suggestion>` — every entry (and every loosened
  tolerance) carries the certified sensitivity + color that justifies it.
- `variance_budget(spec_margin, target) -> Result<f64, ToleranceError>` — the
  budget for `P(|QoI − nom| ≤ spec_margin) ≥ target`, via the inverse normal.
- `ToleranceError` — `NoFeatures` / `NonPositive` / `BadBudget`.

## Invariants

- The allocation TIGHTENS high-sensitivity features and LOOSENS low-sensitivity
  ones, and meets the variance budget exactly (`achieved_variance == budget`).
- `robustness_check` flags where the first-order linearization is exceeded at
  the band extremes (it does not silently trust the linearization).
- Every GD&T suggestion carries a certified sensitivity (with its color) — no
  unjustified tolerance change.

## Error model

Structured `ToleranceError`; no panics.

## Determinism class

Fully deterministic: the allocation, robustness check, and budget are pure
functions of the inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/toleralloc.rs` (Proposal 11, 6 cases): tolerance is spent where
sensitivity is large (tighten high / loosen low, budget met); bad-input
rejection; the band-extremes robustness check (confirm + flag); the GD&T report
attaches a certified sensitivity to every loosened tolerance; the variance
budget follows the in-spec probability; determinism.

## No-claim boundaries

- Sensitivities are SUPPLIED (from Proposal 1 adjoint `∂QoI/∂geometry` fields);
  this crate consumes them and their color, it does not compute them.
- First-order variance propagation assumes independent features and a
  linearization; `robustness_check` is the guard for where that fails, and
  full correlated / higher-order propagation is a refinement.
- The cost model `cᵢ / tᵢ` is a convex placeholder; a real manufacturing cost
  curve is a drop-in.
- Emitting the report into a GD&T/CAD annotation format is a downstream
  integration.
