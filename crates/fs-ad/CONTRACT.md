# CONTRACT: fs-ad

## Purpose and layer
Forward-mode automatic differentiation (plan §6.6 regime 1): generic,
nestable dual numbers and the `Real` generic-scalar contract that lets
kernels run unchanged on values or derivatives; adjoint infrastructure
(IFT, revolve checkpointing incl. binomial + spill, matrix-free tangent
route) and the OPT-IN FrankenTorch reverse-mode tape bridge (feature
`torch-bridge` — the workspace's first ft-* constellation dependency).
Layer: L1.

## Public types and semantics
- `ift::{ift_gradient_matrix_free, MatrixFreeIftReport}` (bead o3ui) —
  the TANGENT route for large N: (∂F/∂u) is only APPLIED (one
  directional-dual pass per application; nothing N×N formed) and the
  LINEAR SOLVER IS CALLER-SUPPLIED (`solve(apply, b)`) — fs-solver's
  Krylov stack plugs in at L3, fs-ad stays solver-agnostic. One solve
  per parameter (right shape for few parameters); solve quality is
  MEASURED with a fresh operator application (`tangent_residual`),
  never trusted from the solver.
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
- `revolve::{checkpointed_adjoint_binomial, beta, binomial_reevals}`
  (bead o3ui) — the Griewank–Walther binomially-OPTIMAL schedule for
  FIXED budgets down to s = 1: measured forward re-evaluations EQUAL
  the theorem minimum r·L − β(s+1, r−1) (gated as equality, budgets
  1..8), bitwise-equal to full storage, and never worse than treeverse
  at equal RAM (gated: 255 vs 356 at L = 100, RAM 8). Budget counts
  PARKED checkpoints (β semantics); RAM peak is budget + 1 (the live
  state), reported truthfully in `RevolveStats`. The optimal split is
  found by scanning the O(r) β-breakpoints of the piecewise-linear DP
  cost — a closed-form "l̂ = β(s, r−1)" split was MEASURED 10 re-evals
  above optimal at s = 2 (kept as the documented rejection).
- `revolve::{checkpointed_adjoint_spilling, SnapshotStore, SpillStats}`
  (bead o3ui) — the ledger-spill escape valve: snapshots beyond the RAM
  budget go to a caller-provided store instead of refusing. CONTRACT on
  the store: byte-exact round-trip (that is what preserves bitwise
  equality — gated at RAM budget 2 on a 100-step chain: 92 spills, 163
  restores, store drained, still bitwise). Keys are written once, read
  many, evicted when dead. The fsqlite-backed adapter belongs to
  fs-ledger (L6); this trait is the L1 seam. The NO-SPILL default
  (`checkpointed_adjoint`) keeps the structured panic.
- `gradcheck::{gradcheck, GradCheckReport}` — the CI gradient-gate
  primitive: dual gradient vs central FD with scale-aware relative error;
  JSON-line Display. Catches the derivative-killing bug class (tested on
  a value()/from_f64 round-trip specimen: O(1) error detected).

- `Real` — the scalar contract (zero/one/from_f64/value, arithmetic ops,
  mul_add, recip, sqrt, abs, exp, ln, sin, cos, tanh, asin, acos, atan,
  atan2, powi). `f64`'s impl routes elementary functions through fs-math
  STRICT det — genericity preserves cross-ISA determinism.
- `bridge::{TapeReal, reverse_gradient, taped_vjp}` (feature
  `torch-bridge`, bead o3ui) — REVERSE MODE via the FrankenTorch scalar
  tape: `TapeReal` is a Copy handle onto a thread-local ft-autograd
  Tape implementing the FULL `Real` surface (Strict execution mode), so
  kernels generic over `Real` get O(cost(f)) full gradients unchanged.
  `taped_vjp(f, x, bar)` = Jᵀ·bar in one backward pass, shaped exactly
  as revolve's `reverse(i, state, bar)` — checkpointing composes with
  taped segments (gated). One tape scope per thread at a time (nested
  scopes panic loudly); use outside a scope panics loudly.
- `Dual<T: Real, const N: usize> { re, eps: [T; N] }` — implements `Real`,
  so NESTED duals give higher-order derivatives from one implementation.
  `Dual64<N>` alias. It also implements fs-la's `GemmScalar` recursively when
  `T` does, which admits packed and nested dual values to the public
  `gemm_scalar_checked` reference kernel without making fs-la depend on fs-ad.
  Structural exact-zero checks include the primal and every derivative lane;
  a zero primal with a live sensitivity is never short-circuited.
- Helpers: `gradient` (N-lane seeding), `jvp` (directional), and
  `second_directional` (nested duals → exact vᵀHv).

## Invariants
- PRIMAL FIDELITY: evaluating through Dual is bit-identical to the scalar
  path (same strict functions, same order, FUSED mul_add primal — tested
  bitwise on 2000 random composite evaluations). A gradient check can never
  be confounded by primal drift.
- Packed lanes ≡ single lanes bitwise (Dual<4> vs 4×Dual<1>, tested).
- In fs-la's scalar-generic GEMM, packed `Dual64<2>` lanes are bitwise equal to
  two `Dual64<1>` executions and the packed primal is bitwise equal to the
  optimized f64 result on the fixed dyadic bridge fixture. Nested duals execute
  through the same recursive `GemmScalar` implementation.
- Comparison convention: PartialEq/PartialOrd compare the primal ONLY
  (branching-on-values; kinks give per-branch one-sided derivatives —
  documented forward-AD semantics).
- Conventions at non-smooth points: abs'(0) = 0 (subgradient choice);
  sqrt'(0) = +inf (honestly unbounded, never clamped).
- Endpoint and full-exponent conventions: asin'(1) = +inf; acos'(-1) =
  -inf; and `powi(i32::MIN)` at x = 1 returns primal 1 with derivative
  -2147483648. These singularities and the full i32 exponent domain are
  represented honestly, never clamped or wrapped.

## Error model
Dual arithmetic is total; derivative singularities produce inf/NaN honestly.
The fs-la bridge retains fs-la's typed `GemmShapeError` and transactional
output guarantee for invalid extents or lengths.

## Determinism class
Deterministic CROSS-ISA (inherits fs-math strict + pure IEEE arithmetic).
EXCEPTION (declared): the `torch-bridge` feature computes through
FrankenTorch's Strict mode — deterministic per ft's own contract but NOT
fs-math det (and `mul_add` is composed, unfused, on the tape). Bridge
results are cross-checked against forward duals to tolerances (1e-10),
never bitwise; kernels needing the bit-contract stay on f64/Dual.

## Cancellation behavior
Straight-line arithmetic; no poll points.

## Unsafe boundary
None.

## Feature flags
`torch-bridge` (default OFF): enables `bridge` and the ft-autograd /
ft-core constellation path dependencies — the workspace's first ft-*
wiring. Default builds stay forward-dual only with zero new deps.

## Conformance tests
The shared `fs-casebook` runner records three bounded G0 cases: exact value,
gradient, JVP, and second-directional known answers for a two-variable
polynomial; selected documented abs/sqrt/inverse-trig/powi boundary conventions;
and treeverse bit equality with full storage, resource bounds, and insufficient-
budget refusal. Canonical frames bind operation identities, formulas, inputs,
directions, expected outputs, schedules, budgets, and refusal policy. Disclosed
seed `0xF5AD_0001` corrupts one exact gradient bit and must produce one
structured red record plus `assert_green` refusal. This portable tranche adds
no torch-bridge, spill-store, binomial, IFT, or gradcheck claim, no performance
claim, and no fresh dual-ISA execution proof; it remains central-package-proof
pending.

`tests/la_dual_bridge_casebook.rs` is the literal fs-ad/fs-la cross-crate G0
entrypoint. It drives `Dual64<2>` and nested Dual values through
`gemm_scalar_checked`, pins exact primal/derivative bits against single-lane and
optimized-f64 references, and exercises full-scalar α/β zero semantics plus
typed preflight refusal with unchanged C. Seed `0xF5A1_0000` flips one expected
derivative bit and must yield one replay-identical red record rejected by
`assert_green`. This proves the reference integration seam, not optimized
dual-SIMD performance or fresh cross-ISA execution.

Gradient-vs-central-FD on 500 random points of a 3-deep composite (rel
< 5e-9); primal bitwise fidelity battery; analytic first+second derivatives
(sin x²); JVP ≡ grad·v; GENERIC NEWTON differentiated through convergence
(d c^(1/3)/dc to 1e-10); kink/singularity conventions; lane-packing
equivalence. INVERSE TRIG (bead t88x): asin/acos/atan/atan2 on `Real`
(f64 → det::*, Dual chain rules incl. binary atan2 partials
(x·dy − y·dx)/(x²+y²)); gradcheck lanes — inverse gauntlet vs central
FD (500 pts, rel < 5e-9), primal BITWISE vs scalar, analytic first +
second derivatives through nested duals, honest endpoints (asin′(1) =
+∞, acos′(−1) = −∞ since acos decreases; never clamped).

## No-claim boundaries
- Qty-typed duals: requires fs-qty generalization to Qty<S: Real> (recorded
  follow-up; until then dimension discipline at kernel boundaries).
- Sparsity-aware Jacobian seeding (graph coloring) — consumer-driven.
- Explicit SIMD for eps arrays (autovectorized today; measured lanes when a
  consumer profiles it).
- fs-la's scalar-generic GEMM is the correctness/reference seam; packed
  microkernels, parallel scheduling, cancellation polling, and performance
  evidence for Dual values remain unclaimed.
- powf/general pow (needs fs-math extensions).
- Matrix-free ADJOINT solves (one solve for MANY parameters): needs
  Jᵀ·v over VECTOR residuals, i.e. reverse mode through F itself —
  the scalar tape bridge covers scalar objectives; the tensor-tape
  vector bridge is the recorded follow-up.
- The fsqlite-backed `SnapshotStore` adapter: fs-ledger (L6) territory;
  the trait seam + byte-exact contract are ready (bead comment posted
  to the fs-ledger owner).
