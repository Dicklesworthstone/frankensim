# CONTRACT: fs-toleralloc

Adjoint-driven tolerance allocation (plan addendum, Proposal 11's commercial
kicker): spend tight tolerances only where sensitivity is large; loosen the
rest with a certified justification.

## Purpose and layer

Layer L4 (optimization). Depends only on `fs-evidence` (`ColorRank` for the
sensitivity's color) and `fs-math` (deterministic scalar kernels). Pure,
deterministic.

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
- `AdmittedCorrelationModel::try_new(namespace, version, digest, dimension,
  lower_factor)` — admits at most 128 ordered axes and a row-major
  lower-triangular factor `L` with canonical `+0.0` above the diagonal,
  nonnegative diagonal, finite coefficients, and binary64-computed row norms
  near one. Therefore the implied `C = L Lᵀ` is positive semidefinite by
  construction; its computed diagonal is near one, but admission is not an
  exact-real diagonal enclosure. The namespace, nonzero version, nonzero
  caller-supplied digest, exact factor, and maximum measured row-norm defect are
  retained; they are provenance, not population validation.
- `propagate_correlated_stack(&AdmittedCorrelationModel, &[CorrelatedStackTerm])
  -> Result<CorrelatedStackReceipt, CorrelatedStackError>` — binds the model's
  positional factor order to bounded unique term names, signed sensitivities
  and their colors, and strictly positive standard deviations supplied by the
  caller. It evaluates
  `aᵀ C a = ||Lᵀa||²` with scaled/compensated arithmetic, where
  `aᵢ = sensitivityᵢ σᵢ`, and retains independent and correlated standard
  deviations/variances plus their signed variance delta. Correlation may
  increase or decrease variance. A zero correlated standard deviation/variance
  is published only when every supplied sensitivity is exactly zero; a
  numerical zero reached from nonzero inputs is refused because this binary64
  lane cannot certify exact cancellation. The signed delta is the subtraction
  of the two published binary64 variances, so its zero is diagnostic rather
  than a certificate of no exact-real correlation effect.
- `ToleranceError` identifies the exact invalid feature field, public argument,
  sampled extreme, canonical-name collision, or derived quantity. Numeric
  reasons are stable `ScalarIssue` values rather than formatted floating-point
  text.
- `CorrelationAdmissionError` and `CorrelatedStackError` identify malformed
  model identity/factor/axis data and non-representable first-order results.

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
- A correlated factor is bounded and lower triangular, and its
  binary64-computed row norms are near one; it is never silently normalized.
  PSD follows from the factor construction. The exact admitted factor and
  caller-supplied external identity survive in every receipt.
- Correlated propagation uses signed sensitivities, preserves caller axis
  order, and reports the counterfactual independent result. It never assumes
  that correlation inflates variance.

## Error model

Structured `ToleranceError`; no panics. NaN never reaches `f64::max` or a
comparison: all scalar inputs are admitted before arithmetic, each derived
quantity is checked before publication, and sampled maxima use an explicit
ordered comparison over finite values.

Correlation-model and stack refusals are separately typed. Dimension, factor
length, term count, and each term-name byte length are checked before retained
work; a receipt retains at most `128 * 256` original term-name bytes, and every
retained lowercase comparison key has the same per-name cap. Oversized
namespace errors retain only a bounded UTF-8 prefix. Normalized projection sums
are compensated in positional order; nonzero normalized terms, products, and
variances that underflow, ambiguous projection zeros, and finite inputs whose
result overflows are refused instead of becoming clean-looking zero or
infinity. Negative zero is refused where an exact-zero sensitivity or factor
coefficient would otherwise acquire two semantic encodings.

## Determinism class

Fully deterministic: the allocation, robustness check, and budget are pure
functions of the inputs. Accumulation and output use input order; canonical-name
collision reporting always identifies the first and colliding input positions.
Correlation admission scans row-major order, and correlated propagation uses
positional factor order for both compensated projections and scale-safe norms.
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
repeatability plus input-order tie-breaking. The additive correlated-stack
battery covers a manufactured finite population where `ρ = 0.8` changes the
variance from the false independent value `2.0` to `3.6`; positive/negative
correlation with signed sensitivities; exhaustive population enumeration and a
fixed-seed 200,000-sample Monte Carlo cross-check; exact identity and caller
positional-term retention; factor, dimension, namespace, digest, measured row
normalization, bounded-name, overflow, and underflow refusals; fail-closed
singular-factor cancellation; all-zero sensitivity binding; and deterministic
replay.

The stricter robustness/admission policy is evidence-semantic. The consuming
`fs-diffreal-e2e` tolerance fixture binds it as
`fs-diffreal-e2e/tolerance-allocation-fixture/v3`; v1/v2 evidence must not be
silently reinterpreted under the sealed-sensitivity, typed-event, sampled-only
policy.

## Perfect-dependence grouped allocation v1

`allocate_grouped(&GroupedDependenceModel, variance_budget, k)` is the bounded
dependency-aware complement to the legacy independent allocator. The admitted
dependence class is deliberately fixed: members of one nonempty group share a
perfectly positively correlated standardized shock, while different groups
are mutually uncorrelated at second order. Sensitivities are strictly positive
magnitudes, so coherent member loadings add and never silently cancel.

For feature sensitivity `s_i`, cost coefficient `c_i`, tolerance `t_i`, and
group `g`, this lane minimizes `sum(c_i / t_i)` subject to

`sum_g (sum_(i in g) s_i t_i / k)^2 <= variance_budget`.

Define `J_g = sum_(i in g) sqrt(c_i s_i)`,
`D = sum_g J_g^(4/3)`, and `alpha = sqrt(variance_budget / D)`. The unique
strictly convex optimum is

`t_i = k alpha sqrt(c_i / s_i) J_g^(-1/3)`.

The implementation evaluates this formula with deterministic max-shifted
log-sum-exp, stable ordinal order, compensated public sums, and one common
binary64 scale correction that preserves every KKT ratio while tightening the
published budget residual. It refuses a vanished positive log-sum-exp term,
nonrepresentable tolerance/loading/cost/variance, empty or oversized group and
feature tables (128 each), empty groups, unstable or colliding names, invalid
group references, and nonpositive/non-finite inputs. No iterative optimizer or
callback is used.

The privately constructed receipt retains the exact caller model and identity,
budget, `k`, feature order and membership, tolerances, baseline actions, costs,
coherent loadings, per-group log shapes/standard deviations/variances, and the
counterfactual independent variances. Each published group standard deviation
is the compensated sum of that receipt's published member loadings, so group
and budget diagnostics audit through one numeric path. The receipt publishes
coherent-minus-independent variance deltas, budget and closed-form cost
residuals, plus the largest log-domain KKT stationarity mismatch. Singleton
groups reduce to the legacy independent allocation algebra.

`tests/dependency_allocation.rs` supplies G0/G3/G5 evidence. Its three-feature
manufactured fixture has exact optimal tolerances `1/2`, group variances `1`
and `1`, total variance `2`, and total cost `8`. Re-evaluating the legacy
independent allocation under the declared coherent groups exceeds the true
budget by a large margin. An unequal coherent pair independently checks the
within-group sensitivity/cost ratio and receipt closure. The battery also
covers common sensitivity and cost rescaling, singleton reduction, exact
retained replay, group/name/reference admission, positive-domain refusals, and
lost log-domain contribution refusal.

Group membership and perfect `+1` dependence are supplied assertions, not
inferred or calibrated facts. This lane makes no negative/partial correlation,
cross-group covariance, signed cancellation, nonlinear response,
distributional, tail/reliability, confidence, causality, or manufacturing-fit
claim. It does not prove the declared block correlation describes the physical
population. Sensitivity and tolerance scalars carry no typed units; compatible
units remain a caller responsibility. The semantic digest is retained but not
computed from, verified against, or authenticated for the group/feature bytes.
Binary64 residuals are diagnostics, not exact-real certificates.

## Structured finite-population propagation v1

`propagate_structured_population(&StructuredPopulationModel)` is the bounded
nonlinear, hierarchical, and mode-switching complement to the analytic
correlated-stack lane. It atomically admits and evaluates one externally
identified finite population with these semantics:

- The model is a parent-before-child tree of at most 8,192 globally keyed
  nodes, 4,096 weighted terminal leaves, and depth 16. Ordinal zero is the one
  root; every other node names an earlier branch parent; leaves cannot own
  children; empty branches are refused. Leaf `NonZeroU64` weights are exact
  finite-population multiplicities. Their checked sum may not exceed `2^53`,
  so every weight is exact when converted to binary64 for moments.
- Each leaf supplies one finite raw clearance and selects one of at most 64
  stable response laws. A law clamps to finite ordered bounds and contains at
  most 64 quadratic pieces over strictly increasing knots. Every interior knot
  explicitly names its lower or upper piece as the equality owner. This makes
  branch selection total and disjoint without callbacks or predicate-order
  ambiguity.
- The selected piece is evaluated in fixed Horner order as
  `(a * x + b) * x + c`, with two explicit multiplies and two explicit adds.
  Nonzero-input products that round to zero and any non-finite intermediate are
  typed refusals; the deadband result is canonical `+0.0`.
- Every receipt retains the caller-supplied identity alongside the complete
  hierarchy, leaf multiplicities, raw and clamped clearances, clamp
  dispositions, response laws, selected mode ordinals, and evaluated outputs.
  Per-node receipts publish weighted mean, immediate-child within variance,
  immediate-child between variance, total population variance, standard
  deviation, and the binary64 decomposition residual. Per-mode receipts publish
  the same observed population partition; unobserved declared modes remain
  explicit with zero weight and absent moments.
- Moment accumulation is deterministic and population-weighted (denominator
  `W`, never `W - 1`). Outputs are normalized by one global maximum magnitude;
  child and mode contributions are accumulated in stable ordinal order with
  compensated sums, then rescaled with checked arithmetic. The root also
  retains an independently accumulated two-pass mean/variance and residuals
  against the hierarchical and mode decompositions.

`tests/structured_propagation.rs` supplies G0/G3/G5 evidence. Its exhaustive
64-occurrence process -> lot -> part gear-backlash population has the exact
oracle `sum(w*y)=14`, `sum(w*y^2)=76`, mean `7/32`, and variance `1167/1024`.
It exercises inclusive deadband boundaries, saturated negative/positive drive
modes, quadratic response, process/lot and mode laws of total variance,
integer-weight replication, affine-response metamorphism, bitwise replay,
topology/key/law admission, and overflow/underflow refusal. The fixture also
shows why an exact zero nominal deadband derivative cannot claim zero variance
for a finite population that switches contact modes.

This lane is descriptive evidence for the supplied finite population. Integer
weights are multiplicities, not calibrated probabilities. The hierarchy does
not imply random-effects independence or causality. Supplied clearances and
quadratic laws are not certified distributions, geometry, contact mechanics,
fits, monotonicity, convexity, continuity, or physical validity. Domain
clamping is a response rule, not proof of interference freedom. There is no
hysteresis, time-dependent switching, arbitrary nonlinear callback,
confidence interval, tail/reliability claim, dependency-aware allocation, or
automatic inference of a correlation model. Correlation must already be
embodied by the finite population or propagated explicitly upstream; it is
never silently inferred or double-counted here. Binary64 residuals are
diagnostics, not exact-real equality certificates. The semantic digest is
caller-supplied and retained but is not computed from, verified against, or
authenticated for the hierarchy/law bytes. Raw clearances and response scalars
carry no typed units; dimensional compatibility and coefficient units remain a
caller admission responsibility.

## No-claim boundaries

- Sensitivities are SUPPLIED (from Proposal 1 adjoint `∂QoI/∂geometry` fields);
  this crate consumes them and their color, it does not compute them.
- `allocate` remains the legacy independent-feature optimizer and old receipts
  are never silently reinterpreted. `allocate_grouped` solves only its admitted
  block-diagonal perfect-`+1` dependence class; neither allocator solves general
  covariance-constrained or nonlinear dependency-aware allocation. The
  additive correlated lane evaluates a supplied fixed first-order stack rather
  than allocating it.
- Factor admission proves PSD structure from `C = L Lᵀ`. Its row-norm test is a
  binary64-computed near-unit check, not an exact-real unit-diagonal proof. It
  does not prove that the caller-supplied digest authenticates the factor bytes,
  that an equivalent matrix has a unique factor (especially when singular), or
  that the model represents a named manufacturing population/process.
- The seed model carries a positional dimension, not semantic axis identifiers.
  A receipt retains the caller's positional term association, but does not
  prove that those positions match the external model's declared axes/order.
- Correlated propagation proves only binary64 first-order moments for supplied
  sensitivities and standard deviations. It makes no Gaussian, quantile,
  reliability, nonlinear, hierarchical, mode-switching, tail, calibration, or
  causal-process claim. `robustness_check` remains the sampled guard for the
  legacy `Allocation`; a correlated nonlinear guard is future work.
- The signed variance delta compares the two published binary64 variances. A
  zero delta neither proves independent axes nor excludes an exact-real effect
  hidden by rounding.
- Correlation coefficients are dimensionless, but this seed API does not yet
  carry typed QoI/axis units. A caller must bind compatible units before
  constructing terms; the receipt alone does not prove dimensional closure.
- The receipt is gear-consumable but no gear flagship currently consumes it.
  Machine-IR lowering, datum/GD&T/surface-texture/fits schemas, assembly/process
  lineage, and Monte Carlo or experimental population validation remain
  downstream work.
- The cost model `cᵢ / tᵢ` is a convex placeholder; a real manufacturing cost
  curve is a drop-in.
- Canonical ambiguity detection uses deterministic Unicode lowercase comparison,
  not full Unicode normalization or locale-sensitive case folding. Callers that
  need a narrower naming grammar must enforce it before allocation.
- Emitting the report into a GD&T/CAD annotation format is a downstream
  integration.
