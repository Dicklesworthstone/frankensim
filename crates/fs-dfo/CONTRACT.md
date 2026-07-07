# CONTRACT: fs-dfo

> Status: PARTIAL — CMA-ES (IGO form), BIPOP restarts, and Nelder–Mead
> are in force; sep/low-rank CMA, NES parameterizations, DE, DIRECT,
> TR-DFO, and fs-exec population waves are recorded follow-up scope.

## Purpose and layer
Derivative-free optimization engines (plan §9.3). Layer: **L4 ASCENT**.
Deps: fs-rand (keyed sampling), fs-la (eigendecomposition), fs-math.
Engines are IR-agnostic (closure objectives) BY DESIGN — routing from
the fs-opt problem IR is a wiring bead once that crate stabilizes
(deliberate collision avoidance, bead 7tv.4 trail).

## Public types and semantics
- `cmaes(f, x0, CmaParams, seed) -> CmaReport` — full-covariance CMA-ES
  with the standard Hansen couplings (log-rank weights, rank-µ + rank-1
  covariance updates, cumulative step-size adaptation): the natural-
  gradient/IGO parameterization. Eigendecomposition via the landed
  cyclic Jacobi with symmetrization + eigenvalue flooring (SPD
  maintenance). Stagnation stops: TolX (σ·√λmax < 1e-12·σ₀) OR TolFun
  (no f_best improvement > 1e-12 relative for 120 generations) — the
  TolFun rule is what frees budget for restart ladders (measured).
- `bipop_cmaes(...) -> BipopReport` — large/small alternation with
  population doubling; per-run budgets scale with λ (≈250 generations);
  the ladder cap counts LARGE runs only; restarts launch from
  deterministic Philox-perturbed starts. The report carries the
  schedule (evidence).
- `nelder_mead(...)` — deterministic simplex polish (no randomness).

## Invariants
1. DETERMINISM: the full evolution is a pure function of the seed
   (keyed Philox sampling, `total_cmp` ranking, lowest-index
   tie-breaks) — bitwise rerun-tested and cross-ISA golden-hashed.
2. IGO INVARIANCE: strictly monotone transforms of the objective give
   BITWISE-identical trajectories (tested with exp and cube transforms;
   precondition documented: monotonicity at the resolution of sampled
   values — x³ underflow and exp saturation break injectivity near
   machine-precision convergence, measured during bring-up).
3. Translation equivariance (behavioral, tested).
4. Covariance stays SPD (symmetrize + floor at each refresh).
5. BIPOP schedule = doublings of the base population (shape-tested).

## Error model
Structured panics for modeling errors (empty dimension). Convergence
failure is DATA: `converged: false` with best-found + diagnostics.

## Determinism class
Bit-deterministic per seed, cross-ISA (golden hash
`0x5441_10a6_afb1_70a1`, bumped once when the TolFun stagnation
criterion was added — semantic justification recorded; verified
identical on both reference ISAs).

## Cancellation behavior
Single-threaded v1; population waves as fs-exec sibling scopes with
cancellation draining are the recorded follow-up (G4 scope there).

## Unsafe boundary
None.

## Feature flags
None.

## Conformance tests
tests/dfo_battery.rs (benchmarks incl. condition-1e6 ellipsoid,
determinism, invariance, BIPOP schedule, NM polish, golden hash);
tests/probe_tmp.rs (success-rate + stagnation-stop regression; filename
is bring-up history).

## No-claim boundaries
- No published-ERT-table parity claims yet (in-repo BBOB-class fixtures
  only; the external COCO battery is follow-up).
- Sep-CMA/low-rank (dim > ~200), NES, DE, DIRECT, TR-DFO: not built.
- No constraint handling (fs-constraint owns kinds; integration later).
- No parallel evaluation waves yet (fs-exec bead).
