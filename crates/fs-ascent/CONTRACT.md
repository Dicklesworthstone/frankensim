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
- `pareto::nondominated_front` — exact two-objective minimization filter for
  tracing and benchmark candidate streams. Dominance requires no worse in
  either objective and strictly better in at least one; bit-identical objective
  rows retain the first occurrence. The output preserves input order, so
  disconnected front segments keep the schedule's stable presentation while
  dominated bridge rows disappear.
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
- `runner::{Study, Packing, StudyReport, StudyError, StudyForkReceipt}`
  (beads ijil and 7tv.21.15) — the
  Problem-IR study runner: ALL variables pack across the manifold
  product (per-block riemann ops keep the sphere spherical to 1e-12
  along the whole path), central-difference tangent gradients through
  `fs_opt::eval` (the live IR is evaluation-only — documented),
  `Problem.budget.limit` threads into the stop algebra as explicit
  `Unlimited` or `Limited(NonZeroU64)`, with
  CORRECT attribution (`StopReason::Budget`, not a composite shrug),
  constraints route to AL/IP/SQP through packed FD adapters, and
  studies are RESUMABLE (clone = checkpoint; split runs bitwise equal
  INCLUDING eval accounting — the cached current-objective is what
  makes segment boundaries invisible). Every checkpoint is bound to the
  full-width semantic identity of that cache's admitted `Problem`;
  `try_run` refuses a different identity before mutation. Objective steering
  uses `fork_for`: an exact variable schema is mandatory, the parent is
  untouched, and the child keeps the branch point/driver configuration while
  clearing objective history, step/evaluation counters, and problem-dependent
  caches. An unchanged semantic problem is refused so a fork cannot reset the
  same problem's budget accounting. The returned receipt retains both problem
  IDs, parent accounting, branch-point float bits, and driver bits.

## Invariants

- Deterministic trajectories from identical inputs (G5-tested);
  resumable states bitwise across split points.
- Resumption cannot reinterpret cached state under a different semantic
  problem. Reweighting is an explicit world fork with immutable parent state
  and fresh branch-local accounting.
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
boundary. Nondominance filtering also refuses malformed public candidates,
including non-finite objective, gradient-residual, or KKT evidence. Optimizer
non-convergence with a well-formed problem is a
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

`tests/wfg_moo_battery.rs` covers WFG-scale front geometry without claiming
the full deceptive or multimodal transformation stack. Concave, convex, and
WFG1-class mixed fixtures retain their measured sweep behavior as linted,
wire-validated `fs-obs` `ConformanceCase` rows. The WFG2-class disconnected
fixture samples the deterministic five-basin shape on a fixed 1001-point
parameter grid, runs the public stable nondominance filter, and requires six
separated record-low segments with both objective seams retained and every
dominated bridge removed. Repeat and reversed-input runs must reproduce the
same complete front identity after presentation order is restored. An exact
duplicate of the upper seam catches a test-local comparator mutant that treats
equality as dominance and incorrectly lets the two copies erase each other.
The structured receipt binds crate versions, grid, objective scales, WFG2
basin count, every retained parameter/objective bit, segment/front counts, and
the caught mutation with fixed input seed zero.

`tests/bbob_budget_ledger.rs` and `tests/gradient_budget_ledger.rs` emit the
four-metric observation envelope for every optimizer fixture: observed nfev/ERT,
pinned nfev ceiling, observed success rate, and minimum success rate, all as
wire-validated `fs-obs::BenchmarkResult` events with machine id zero. Their
gate metadata comes only from `tests/support/budget_trend.rs`, not duplicated
tuple literals. `tests/budget_trend_manifest.rs` (bead 7tv.21.12) audits the
closed fourteen-row `(suite, kernel)` inventory, refuses missing/duplicate or
metadata-drifted rows, and deterministically renders schema
`frankensim-ascent-budget-trend-v1` for central trend ingestion.

`tests/cutest_scale_battery.rs` (bead 7tv.21.13) extends the small-dimensional
CUTEst-class oracle tranche with extended Rosenbrock, extended Powell singular,
and variably-dimensioned families at 16 and 64 variables. Before an optimizer
result is admitted, each analytic gradient must match FrankenScipy's independent
finite-difference gate and the same gate must reject a one-coordinate mutant.
L-BFGS must then reach the known zero objective within a pinned evaluation
ceiling and reproduce the complete evaluated-point trace, public state, history,
accounting, and report bit for bit. Every case is a wire-validated `fs-obs`
`ConformanceCase` row with fixed-input seed zero and a versioned identity over
the fixture, engine configuration, initial point, and every evaluated point.

`tests/runner_battery.rs` (4 cases): Problem-IR product-manifold packing,
problem-owned budget attribution, clone-checkpoint replay, and packed
constraint routing. All aggregate outcomes are linted and wire-validated
`fs-obs` `ConformanceCase` events with fixed-input seed zero. The G5 replay
case admits the problem and binds its full-width semantic ID, study
parameters, initial point, exact final point, complete public objective
history, step/evaluation accounting, and terminal report into canonical
identities. An object-shaped companion receipt retains the exact initial,
final, and objective-history float bits plus the reference, independent
repeat, and three cut-bound checkpoint-resume identities. Each per-cut row
retains the checkpoint state, stop reason, step/history counts, and study/report
evaluation accounting. The case passes only if all three cuts are genuine
iteration-cap splits and every canonical final state is byte-identical to the
uninterrupted reference.

`tests/runner_fork_battery.rs` (bead 7tv.21.15) covers gradient-study P9
world-forking. A two-step equal-weight trunk holds a live cached objective;
direct continuation under a different semantic problem must return a typed,
mutation-free mismatch, while an unchanged-problem fork must refuse a hidden
accounting reset. Opposite 0.95/0.05 objective reweightings then fork from the
exact same point, start with empty branch-local history/accounting, converge
toward opposite analytic optima, and reproduce both fork receipts and complete
public terminal states bit for bit on independent repeats. A renamed variable
schema is refused before a child exists. The wire-validated `fs-obs` receipt
binds parent/child problem identities, branch-point/driver bits, parent
accounting, terminal states, and fixed input seed zero; the aggregate
`ConformanceCase` is emitted only after all direct assertions pass.

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
- The runner G5 receipts cover one fixed, sequential projected-gradient
  Problem-IR fixture with in-process clone checkpoints and one finite
  in-process two-objective reweight-fork fixture. They do not claim persisted
  checkpoint serialization, generalized lineage replay, other optimizer-family
  studies, `Cx` or worker-count replay, cancellation-storm recovery, cross-ISA
  equality, fs-ledger artifact replay, or performance. The fork receipt is
  ledger-ready data, not evidence that fs-ledger persistence is wired. Those
  claims require their own retained study/host evidence. Expectations that fail
  before an aggregate event remain ordinary Rust test diagnostics.
- Constraint Jacobian-transpose callbacks remain a mathematical trust
  boundary: fs-ascent checks exact dimensions, finite values, and the
  mandatory `J^T 0 = 0` linearity identity, but independent derivative
  verification is required before stronger correctness claims.
- The budget-trend manifest is a machine-independent regression-gate
  declaration over objective-evaluation counts. It is not a wall-clock or
  throughput claim, does not establish cross-ISA optimizer equivalence, and
  does not itself persist or compare historical runs. Central CI/ledger tooling
  owns retention, history selection, alert policy, and machine evidence.
- The CUTEst-scale battery embeds three smooth analytic problem families; it
  does not parse CUTEst SIF files, link a CUTEst runtime, cover the full CUTEst
  taxonomy, or establish wall-clock or cross-ISA performance. Its 16/64-variable
  rows prove finite gradient-gate, evaluation-budget, and same-process replay
  behavior for the named starts only. Wider dimensions, constrained/nonsmooth
  families, retained two-ISA execution, and historical trend persistence remain
  separate evidence obligations. The versioned `fs-obs` trace identity is a
  deterministic legacy-FNV drift fingerprint, not a cryptographic authenticity
  or scientific-authority anchor.
- The WFG2 disconnected battery certifies exact nondominance semantics and the
  sampled five-basin front geometry only. It does not implement the complete
  WFG2 reduction/transform pipeline, prove coverage between grid samples,
  establish optimizer convergence across dominated gaps, cover all nine WFG
  families, or produce cross-ISA/performance evidence. Its retained identity
  is replay evidence, not cryptographic scientific authority.
