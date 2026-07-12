# CONTRACT: fs-robust

Objective epistemics (plan addendum, Proposal F): the three colors applied to
the GOAL itself.

## Purpose and layer

Layer L4 (optimization). Depends only on `fs-evidence` (UTIL, the `Color`
lattice). Pure risk + coloring algebra; the environment sample paths (common
random numbers) are supplied by the caller (fs-scenario).

## Public types and semantics

- `cvar(samples, alpha)` — finite-sample Conditional Value at Risk: the mean
  of exactly the worst `(1 - alpha)` empirical mass, with fractional weight on
  the boundary order statistic when the tail mass is non-integral.
- `weakest_color(&[Color]) -> Option<Color>` — the lowest-rank color.
- `ColoredObjective { design, cost_samples, input_colors }` — `robust_value`
  (CVaR), `nominal_value` (mean), `headline_color` (weakest input color; errors
  if un-colored).
- `robust_optimum(&[ColoredObjective], alpha) -> RobustReport` — the design
  minimizing CVaR; ENFORCES the amended optimization contract (every candidate
  must be colored) and returns the weakest-input headline color.
- `dominated_by_nominal(robust_cost, nominal_plus_safety_cost) -> Result<bool,
  RobustError>` — the Proposal-F kill-criterion test; non-finite costs are
  refused rather than silently suppressing domination.
- `fragility_curve(capacity_samples, intensities, color) -> ColoredFragility` —
  `P(failure)` = fraction of capacities below each intensity, with a color band;
  output points are canonically sorted by finite intensity.
- `RobustError` — `EmptySamples` / `BadAlpha` / `BadSample` /
  `UncoloredObjective` / `NoCandidates`.

## Invariants

- CONTRACT: `robust_optimum` and `headline_color` REFUSE an un-colored
  objective (no optimizing a fiction with certified precision).
- WEAKEST-INPUT RULE: a headline's color is the minimum-rank input color — a
  verified solve under an estimated hazard is an estimated answer.
- CVaR of the tail is `>=` the mean; `alpha ∈ (0, 1)` enforced; non-finite
  samples are refused instead of sorted into risk or headline values. CVaR and
  nominal means use bounded convex combinations, so finite constant samples
  near `f64::MAX` do not overflow through an intermediate sum. Both statistics
  canonically total-order samples, making them permutation-invariant.
- A fragility curve is sorted and nondecreasing in intensity (`P(failure)` =
  CDF of demand exceeding capacity), independent of caller input order.

## Error model

Structured `RobustError` values; no panics.

## Determinism class

Fully deterministic: CVaR, the optimum, and fragility curves are pure functions
of the supplied samples (sampling / common random numbers are the caller's).

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/robust.rs` (Proposal F, 11 cases): CVaR weights integral and non-integral
worst tails + rejects bad inputs; the weakest-input color rule; robust vs
nominal optima diverge; the un-colored-objective contract (+ no-candidates);
extreme finite means remain finite; the kill-criterion dominance test;
non-finite dominance costs refuse; unsorted fragility inputs produce monotone
canonical curves; determinism.

## No-claim boundaries

- The default robust measure is CVaR over the supplied samples; general
  distributionally-robust optimization over an explicit ambiguity set, and
  adjoint hazard-sensitivity (Proposal 1) via sample-average approximation with
  common random numbers, are the fuller Proposal-F deliverable — this crate
  provides the coloring + CVaR + contract core.
- Enforcing "no optimization against an un-colored objective" at the actual
  `fs-opt` boundary is that crate's integration; here it is enforced in
  `robust_optimum`.
- Environment/hazard/cost artifacts as colored, versioned, provenance-carrying
  distribution objects are fs-scenario's upgrade; this crate consumes their
  samples + colors.
