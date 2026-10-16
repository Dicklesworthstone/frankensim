# CONTRACT: fs-robustopt-e2e

ProofRobust — pick a design that is certified globally optimal (SOS proof) AND
robust to worst-case perturbation (CVaR). Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-sos` (SOS global-optimality certificates), `fs-robust` (CVaR +
colored objectives), `fs-evidence` (Verified/Estimated). Deps point downward.

## Public types and semantics

- `Family { name, a, b, c }` — nominal cost `p(x)=ax²+bx+c` (`a>0`); `x_star`,
  `perturbed_costs(sigma, n)`.
- `FamilyVerdict { name, nominal_cost, x_star, nominal_color, robust_cost }`.
- `run_campaign(&[Family], alpha, sigma, n) -> RobustOptReport` — proves each
  global optimum by composing a radius-scoped SOS lower bound (on a radius
  containing the vertex), convex-quadratic monotonicity outside that radius,
  and an interval-evaluated attained upper bound; ranks by CVaR robustness.
- `demo_families()` — three families with `x*=2` but different curvature.

## Invariants

- Every finite convex family's nominal optimum is SOS-`Verified` by a true
  enclosure containing `c − b²/4a`; a non-convex, non-finite, or numerically
  unrepresentable family is `Estimated` (no global certificate).
- ROBUSTNESS REORDERS: the lowest-nominal family ("champion") is not the CVaR
  winner ("flat") — a flatter family wins under worst-case perturbation.
- NO LAUNDERING: the robust (CVaR) headline is `Estimated` (a sample statistic).
- Deterministic (a fixed tolerance grid; no RNG).

## Error model

`run_campaign` panics only on an empty family list.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/robustopt.rs` (3): global optima SOS-proven + robustness reorders; a
downward family is not certified; determinism.

## No-claim boundaries

Nominal cost is a scalar convex quadratic (univariate SOS is complete); the
multivariate moment/Positivstellensatz SDP and a full NSGA Pareto front are the
fuller `fs-sos`/`fs-dfo` deliverables. Perturbation is a 1-D tolerance grid.
