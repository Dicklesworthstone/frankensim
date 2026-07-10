# fs-ascent CONTRACT

## Purpose and layer

Layer: **L4 ASCENT** (deps: fs-adjoint/fs-solver L3, fs-opt L4,
fs-la L1, fs-math L0). The gradient-based optimizer stack (plan
§9.2): the pipeline raw adjoint gradient → Sobolev smoothing →
optional Riemannian projection → optimizer, with ENGINES here and
problem STRUCTURE (typed graphs, manifold metadata) in fs-opt. Every
returned optimum carries a certificate (gradient norm, KKT residuals)
so converged and stalled are distinguishable outcomes.

## Public types and semantics

- `FnGrad`/`FnHv` — the callback shapes every engine consumes.
- `LbfgsState` — resumable L-BFGS (two-loop recursion, m ≈ 17-class
  memory, strong-Wolfe line search, curvature-pair skip guard,
  steepest-descent fallback when stale memory yields a non-descent
  direction after resume on nonconvex terrain). Checkpoint = clone;
  split runs bitwise-equal to straight runs.
- `wolfe::strong_wolfe` — bracket + bisection zoom, deterministic
  control flow, evaluation counts returned.
- `trust::trust_region_newton` — Steihaug-CG on the quadratic model
  with NEGATIVE-CURVATURE boundary steps (counted in the report),
  classical radius laws (¼/¾ thresholds); `hv_fd_of_gradients` is the
  interim Hessian-vector product with its O(√ε) accuracy in the name
  (second-order adjoints are recorded follow-up).
- `auglag::augmented_lagrangian` — PHR augmented Lagrangian
  (equalities + inequalities) over L-BFGS inner solves, classical
  penalty schedule, multipliers returned, and a `KktResidual`
  certificate (stationarity, primal feasibility, inequality dual
  feasibility, complementarity) on every outcome.
- `riemann::{tangent_project, retract, RiemannianLbfgs}` — the
  manifold OPERATIONS for fs-opt's metadata (Rn/Sphere/So3;
  Stiefel is metadata-only and panics loudly), plus Riemannian L-BFGS
  with projection-based vector transport and Armijo curve search;
  reports carry the worst manifold violation along the path.
- `pareto::{weighted_sum_sweep, epsilon_constraint_sweep}` —
  gradient-based Pareto TRACING: warm-started L-BFGS continuation
  along a weight schedule (exact on convex fronts) and warm-started
  augmented-Lagrangian ε-constraint continuation (covers CONCAVE
  fronts, where weighted sums provably collapse to extremes —
  exhibited in the battery); every ε-constraint point carries its
  KKT certificate.
- `stop::{StopRule, StopReason}` — the stopping-condition algebra:
  grad-norm / objective / budget / stall leaves under Any/All
  combinators, with REASON attribution in every report.
- `interior::interior_point` (bead ijil) — the LOG-BARRIER
  interior-point option: barrier subproblems minimized by the
  resumable L-BFGS, equalities on the AL term, softmax phase-1 for
  infeasible starts, barrier multiplier estimates ν = μ/(−cᵢ) feeding
  the SAME KktResidual certificate. Fixture parity with AL gated
  (same optimum, same active multiplier on the landed KKT fixture;
  circle fixture to the analytic ν = 0.5).
- `sqp::sqp` (bead ijil) — active-set SQP for tightly-constrained
  small-dimension POLISH: damped-BFGS Lagrangian Hessian, QP
  subproblems via the dense fs-la KKT factorization, ℓ1-merit
  backtracking; parity with AL to 5e-8 on x AND multipliers; measured
  warm-start polish envelope 10 iterations to 5e-11 KKT
  (identity-seeded BFGS needs curvature pairs — the gate is the
  measurement, not a wish). Large-scale SQP (sparse KKT, TR
  globalization) is recorded follow-up.
- `runner::{Study, Packing, StudyReport}` (bead ijil) — the
  Problem-IR study runner: ALL variables pack across the manifold
  product (per-block riemann ops keep the sphere spherical to 1e-12
  along the whole path), central-difference tangent gradients through
  `fs_opt::eval` (the live IR is evaluation-only — documented),
  `Problem.budget.max_evals` threads into the stop algebra with
  CORRECT attribution (`StopReason::Budget`, not a composite shrug),
  constraints route to AL/IP/SQP through packed FD adapters, and
  studies are RESUMABLE (clone = checkpoint; split runs bitwise equal
  INCLUDING eval accounting — the cached current-objective is what
  makes segment boundaries invisible).

## Invariants

- Deterministic trajectories from identical inputs (G5-tested);
  resumable states bitwise across split points.
- L-BFGS curvature pairs are admitted only with sᵀy above the
  roundoff floor (SPD memory preserved on nonconvex problems).
- Riemannian iterates remain ON the manifold to roundoff (the
  retraction is the constraint — no renormalization hacks; violation
  tracked and reported).
- KKT certificates are computed from the PROBLEM's callbacks, not
  the inner solver's internal state.
- KKT certification is fail-closed: decision vectors, objective and
  constraint values, gradients, multiplier vectors, and Jacobian-
  transpose actions must have exact dimensions and finite entries.
  Inequality dual feasibility is an explicit residual, so a negative
  multiplier cannot cancel stationarity into a false certificate.

## Error model

Structured panics on dimension mismatches, non-finite callback data,
non-positive or non-finite tolerances, non-descent directions handed
to the line search, and metadata-only manifolds used in descent
(modeling errors surfaced loudly). Pareto sweeps reject non-finite
decisions/objective values/gradients, invalid weights or epsilons,
non-positive tolerances, and gradient dimension mismatches at the API
boundary. Optimizer non-convergence with a well-formed problem is a
REPORTED outcome (reason + certificate), never a panic.

## Determinism class

Bit-deterministic across runs; golden FNV-64 over L-BFGS, trust-
region, and Riemannian trajectories: `0xb28d_3cf4_99e8_9071`,
recorded on Apple M4 Pro, verified on Threadripper (x86_64).

## Cancellation behavior

Iteration-granular (the fs-solver pattern): states are complete
between `run` calls; budgets flow through the stop algebra. Cx
wiring is driver scope.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

`tests/ascent_battery.rs` (8 cases): L-BFGS on Rosenbrock(10) inside
an evaluation envelope (< 600 evals) certifying by gradient norm;
bitwise resume at 3 cut points + G5 repeat; THE FLAGSHIP full
pipeline — PDE-constrained density misfit (fs-adjoint DensityPoisson)
with the IFT adjoint gradient verify_gradient-GATED before the
optimizer sees it, driven to ≤ 1e−12 misfit by L-BFGS; trust-region
Newton–Krylov solving Rosenbrock with negative-curvature steps
counted and gradient certification; augmented Lagrangian on an
equality+active-inequality fixture recovering the analytic optimum
with all four KKT residuals ≤ 1e−6 and a positive active multiplier;
Riemannian L-BFGS minimizing a Rayleigh quotient on S¹¹ to the
Jacobi-verified smallest eigenvalue with manifold violation ≤ 1e−14
along the whole path; stop-rule attribution (Budget and Stall
distinguished); TR-Newton exact-Hv vs FD-Hv head-to-head (3 iters /
43 solves vs 5 iters / 1625 evaluations); cross-ISA golden hash.
`tests/pareto_battery.rs` (4 cases): the convex quadratic pair's
closed-form front matched to 2.2e−16 (machine precision — the
weighted-sum optimum is analytic); Fonseca–Fleming CONCAVE front —
weighted sums collapse to exactly 2 clusters while ε-constraint
covers f₁ ∈ [0.15, 0.85] on the true Pareto set (x₁ = x₂ to 1e−4)
with all KKT residuals ≤ 1e−7; bitwise replay; Pareto golden
`0x301b_04df_db91_3965`; fail-fast guards for invalid weights,
epsilons, tolerances, decision vectors, objective values, and gradient
dimensions.

## No-claim boundaries

- Second-order adjoints (adjoint-of-adjoint Hv) are follow-up scope;
  `hv_fd_of_gradients` is the documented O(√ε) interim, and exact
  duals cover small parameter counts.
- No Adam-family stochastic optimizer yet (surrogate-training-
  adjacent; lands with its consumer). Interior-point + SQP: RESOLVED
  (bead ijil; AL remains the constrained default). FrankenScipy
  cross-validation: RESOLVED (bead ijil) with the API-check outcome on
  record — fsci-opt 0.1.0's `minimize`/`slsqp` accept NO general
  constraint callbacks (bounds/penalty only), so the oracle pairing is
  unconstrained parity vs `minimize(Bfgs/LBfgsB)` (agree within 1e-4
  on fsci's own Rosenbrock from a shared global-basin start) and
  constrained parity vs `differential_evolution_constrained` (seeded;
  AL/IP/SQP all within 1e-6 of the DE oracle). MEASURED FINDING kept:
  Rosenbrock n ≥ 4 is BIMODAL — from the classic start our L-BFGS and
  fsci's BFGS landed in DIFFERENT basins (both genuinely stationary;
  basin choice is not a parity criterion), and from another start the
  roles flipped.
- Riemannian line search is Armijo (strong Wolfe on manifolds needs
  transported-derivative bookkeeping — follow-up); Stiefel and
  fixed-volume level sets are metadata-only until their consumer
  beads supply retractions.
- Multi-variable manifold products + the Problem-IR driver: RESOLVED
  (bead ijil, `runner` module). Remaining runner scope: reverse-mode
  IR gradients (the live IR is evaluation-only), L-BFGS/TR engines
  behind the runner (projected gradient is the v1 driver), and
  ledgered study artifacts (fs-ledger wiring) — recorded follow-ups.
