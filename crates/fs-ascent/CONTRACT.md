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
  certificate (stationarity, feasibility, complementarity) on every
  outcome.
- `riemann::{tangent_project, retract, RiemannianLbfgs}` — the
  manifold OPERATIONS for fs-opt's metadata (Rn/Sphere/So3;
  Stiefel is metadata-only and panics loudly), plus Riemannian L-BFGS
  with projection-based vector transport and Armijo curve search;
  reports carry the worst manifold violation along the path.
- `stop::{StopRule, StopReason}` — the stopping-condition algebra:
  grad-norm / objective / budget / stall leaves under Any/All
  combinators, with REASON attribution in every report.

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

## Error model

Structured panics on dimension mismatches, non-descent directions
handed to the line search, and metadata-only manifolds used in
descent (modeling errors surfaced loudly). Optimizer non-convergence
is a REPORTED outcome (reason + certificate), never a panic.

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
with all three KKT residuals ≤ 1e−6 and a positive active multiplier;
Riemannian L-BFGS minimizing a Rayleigh quotient on S¹¹ to the
Jacobi-verified smallest eigenvalue with manifold violation ≤ 1e−14
along the whole path; stop-rule attribution (Budget and Stall
distinguished); cross-ISA golden hash.

## No-claim boundaries

- Second-order adjoints (adjoint-of-adjoint Hv) are follow-up scope;
  `hv_fd_of_gradients` is the documented O(√ε) interim, and exact
  duals cover small parameter counts.
- No Adam-family stochastic optimizer yet (surrogate-training-
  adjacent; lands with its consumer). No interior-point or SQP
  (augmented Lagrangian is the constrained default; IP/SQP join with
  a consumer needing tightly-constrained polish). No FrankenScipy
  cross-validation battery yet (the §12 oracle contract lane).
- Riemannian line search is Armijo (strong Wolfe on manifolds needs
  transported-derivative bookkeeping — follow-up); Stiefel and
  fixed-volume level sets are metadata-only until their consumer
  beads supply retractions.
- Multi-variable manifold products and the fs-opt Problem-IR driver
  (optimize a `fs_opt::Problem` end-to-end) land with the study
  runner; engines here consume callbacks.
