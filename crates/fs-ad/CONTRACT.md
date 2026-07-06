# CONTRACT: fs-ad

## Purpose and layer
Forward-mode automatic differentiation (plan §6.6 regime 1): generic,
nestable dual numbers and the `Real` generic-scalar contract that lets
kernels run unchanged on values or derivatives. Layer: L1. (Reverse mode /
IFT adjoint infrastructure = the fs-ad-adjoint-infra bead.)

## Public types and semantics
- `ift::{ift_gradient, IftReport}` — implicit-function-theorem adjoints:
  dJ/dp at a solution of F(u,p)=0 via one adjoint solve
  ((∂F/∂u)ᵀλ = ∂J/∂u through fs-la LU `solve_transpose`); Jacobians built
  densely column-by-column with single-lane duals (deterministic seeding
  order). `IftReport` carries the PRIMAL residual (the gradient formula is
  exact only at F = 0 — callers get the honesty number) and the adjoint
  residual. Singular ∂F/∂u surfaces as `FactorError` (the IFT hypothesis
  failed), never a wrong gradient.
- `revolve::{checkpointed_adjoint, full_adjoint, min_budget,
  RevolveStats}` — binary-treeverse checkpointed reverse sweeps: peak
  snapshots ≤ ⌈log₂L⌉+1 (asserted via instrumentation), forward
  re-evaluations ≤ L·⌈log₂L⌉ (asserted). HEADLINE INVARIANT: the
  checkpointed adjoint is BITWISE equal to the full-storage adjoint
  (deterministic recomputation reproduces identical states) — tested.
  Insufficient budget is a structured panic, not a silent overrun.
- `gradcheck::{gradcheck, GradCheckReport}` — the CI gradient-gate
  primitive: dual gradient vs central FD with scale-aware relative error;
  JSON-line Display. Catches the derivative-killing bug class (tested on
  a value()/from_f64 round-trip specimen: O(1) error detected).

- `Real` — the scalar contract (zero/one/from_f64/value, arithmetic ops,
  mul_add, recip, sqrt, abs, exp, ln, sin, cos, tanh, powi). `f64`'s impl
  routes elementary functions through fs-math STRICT det — genericity
  preserves cross-ISA determinism.
- `Dual<T: Real, const N: usize> { re, eps: [T; N] }` — implements `Real`,
  so NESTED duals give higher-order derivatives from one implementation.
  `Dual64<N>` alias.
- Helpers: `gradient` (N-lane seeding), `jvp` (directional), and
  `second_directional` (nested duals → exact vᵀHv).

## Invariants
- PRIMAL FIDELITY: evaluating through Dual is bit-identical to the scalar
  path (same strict functions, same order, FUSED mul_add primal — tested
  bitwise on 2000 random composite evaluations). A gradient check can never
  be confounded by primal drift.
- Packed lanes ≡ single lanes bitwise (Dual<4> vs 4×Dual<1>, tested).
- Comparison convention: PartialEq/PartialOrd compare the primal ONLY
  (branching-on-values; kinks give per-branch one-sided derivatives —
  documented forward-AD semantics).
- Conventions at non-smooth points: abs'(0) = 0 (subgradient choice);
  sqrt'(0) = +inf (honestly unbounded, never clamped).

## Error model
Total functions; derivative singularities produce inf/NaN honestly.

## Determinism class
Deterministic CROSS-ISA (inherits fs-math strict + pure IEEE arithmetic).

## Cancellation behavior
Straight-line arithmetic; no poll points.

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
Gradient-vs-central-FD on 500 random points of a 3-deep composite (rel
< 5e-9); primal bitwise fidelity battery; analytic first+second derivatives
(sin x²); JVP ≡ grad·v; GENERIC NEWTON differentiated through convergence
(d c^(1/3)/dc to 1e-10); kink/singularity conventions; lane-packing
equivalence.

## No-claim boundaries
- Qty-typed duals: requires fs-qty generalization to Qty<S: Real> (recorded
  follow-up; until then dimension discipline at kernel boundaries).
- Sparsity-aware Jacobian seeding (graph coloring) — consumer-driven.
- Explicit SIMD for eps arrays (autovectorized today; measured lanes when a
  consumer profiles it).
- powf/general pow (needs fs-math extensions).
- Reverse mode & tape bridge (fs-ad-adjoint-infra bead).
