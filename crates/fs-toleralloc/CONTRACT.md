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
  Loosen / Unchanged vs baseline). Normalization is evaluated in log space so
  finite positive inputs do not overflow merely because a sensitivity is
  squared; a mathematically unrepresentable public result is refused.
- `robustness_check(&Allocation, extreme_qois, nominal_qoi, k, margin) ->
  Result<RobustnessVerdict, ToleranceError>` — compares the first-order
  `linearized_std` against the QoI at sampled tolerance-band extremes;
  `confirmed` iff the extremes stay within `k · linearized_std · (1 + margin)`.
  An empty extreme set has no evidentiary meaning and is a structured refusal,
  never a vacuous confirmation.
- `gdt_report(&Allocation) -> Result<Vec<Suggestion>, ToleranceError>` — every
  entry (and every loosened tolerance) carries the certified sensitivity + color
  that justifies it. Forged/deserialized items with unsafe fields or ambiguous
  names are refused before report publication.
- `variance_budget(spec_margin, target) -> Result<f64, ToleranceError>` — the
  budget for `P(|QoI − nom| ≤ spec_margin) ≥ target`, via the inverse normal.
  The quantile is evaluated from the central probability or upper-tail mass
  directly, so representable targets adjacent to zero and one do not first
  round to the singular CDF endpoints.
- `ToleranceError` identifies the exact invalid feature field, public argument,
  sampled extreme, canonical-name collision, or derived quantity. Numeric
  reasons are stable `ScalarIssue` values rather than formatted floating-point
  text.

## Invariants

- The allocation TIGHTENS high-sensitivity features and LOOSENS low-sensitivity
  ones, and meets the variance budget exactly (`achieved_variance == budget`).
- Every admitted scalar is finite and in its declared domain. Every published
  tolerance, cost, variance, standard deviation, deviation, and bound is finite;
  positive quantities remain strictly positive.
- Feature names are non-empty, have no surrounding whitespace or control
  characters, and are unique under locale-independent Unicode lowercase
  comparison. Output order is input order, which is also the stable tie-break.
- `robustness_check` flags where the first-order linearization is exceeded at
  the band extremes. It refuses empty, non-finite, negative-domain, or
  unrepresentable evidence rather than silently trusting the linearization.
- Every GD&T suggestion carries a certified sensitivity (with its color) — no
  unjustified tolerance change.

## Error model

Structured `ToleranceError`; no panics. NaN never reaches `f64::max` or a
comparison: all scalar inputs are admitted before arithmetic, each derived
quantity is checked before publication, and sampled maxima use an explicit
ordered comparison over finite values.

## Determinism class

Fully deterministic: the allocation, robustness check, and budget are pure
functions of the inputs. Accumulation and output use input order; canonical-name
collision reporting always identifies the first and colliding input positions.
Bitwise reproducibility holds CROSS-ISA: every transcendental routes through
`fs_math::det` (bead frankensim-lyms; platform libm is not correctly rounded
and differs across ISAs), and the crate is registered in the `check-libm`
doctrine lint. `sqrt` stays primitive (IEEE-754 correct rounding).

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/toleralloc.rs` covers tolerance direction and budget adherence;
field-specific zero/negative/NaN/infinity rejection; empty, unstable, duplicate,
and case-colliding names; finite boundary behavior and derived overflow refusal;
empty/poisoned/unrepresentable robustness evidence; GD&T sensitivity carriage;
probability-to-variance conversion; G3 common-sensitivity rescaling; and G5
repeatability plus input-order tie-breaking.

The stricter robustness/admission policy is evidence-semantic. The consuming
`fs-diffreal-e2e` tolerance fixture binds it as
`fs-diffreal-e2e/tolerance-allocation-fixture/v3`; v1/v2 evidence must not be
silently reinterpreted under the sealed-sensitivity, typed-event, sampled-only
policy.

## No-claim boundaries

- Sensitivities are SUPPLIED (from Proposal 1 adjoint `∂QoI/∂geometry` fields);
  this crate consumes them and their color, it does not compute them.
- First-order variance propagation assumes independent features and a
  linearization; `robustness_check` is the guard for where that fails, and
  full correlated / higher-order propagation is a refinement.
- The cost model `cᵢ / tᵢ` is a convex placeholder; a real manufacturing cost
  curve is a drop-in.
- Canonical ambiguity detection uses deterministic Unicode lowercase comparison,
  not full Unicode normalization or locale-sensitive case folding. Callers that
  need a narrower naming grammar must enforce it before allocation.
- Emitting the report into a GD&T/CAD annotation format is a downstream
  integration.
