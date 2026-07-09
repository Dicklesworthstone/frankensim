# CONTRACT: fs-dfo

> Status: PARTIAL — CMA-ES (IGO form), BIPOP restarts, Nelder–Mead,
> balanced entropic OT, MOO helpers, and discrete-support
> Wasserstein-DRO inner sup are in force; sep/low-rank CMA, NES
> parameterizations, DE, DIRECT, TR-DFO, and fs-exec population waves
> are recorded follow-up scope.

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
- `wasserstein_worst_case(losses, costs, n, rho) -> DroReport` —
  exact dual evaluation for a discrete-support Wasserstein-DRO inner
  supremum, with deterministic dual minimization and a recovered
  support distribution. Kink cases split mass fractionally across
  active supports instead of pretending a single argmax distribution
  realizes the dual value.

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
is bring-up history); tests/ot_battery.rs (4 cases): Sinkhorn
marginal feasibility ≤ 1e−8 and cost symmetry ≤ 1e−8; the ε-ladder
approaching the 1D monotone-coupling CLOSED FORM monotonically
(0.11 → 0.005 at ε = 0.05 → 0.002); translation covariance
(W₂² of equal translates = t² within 0.2%); OT golden
`0x58eb_8443_224c_a689`; tests/moo_battery.rs (7 cases): hypervolume vs
hand-computed 2D/3D values including dominated/out-of-reference
degenerate cases; non-dominated-sort front assignment exact;
NSGA-II on ZDT1/ZDT2 at standard budgets (pop 80 × 200 generations —
short runs measurably leave the f2-minimal arm unexplored, documented
in the test) with mean front gap ≤ 0.05 (measured ≤ 0.0008), full f1
spread, hypervolume beating scrambled-Sobol random at MATCHED
evaluations (0.87 vs 0.19 / 0.54 vs 0.00), and bitwise replay; knee
detection hitting a synthetic elbow exactly; CVaR Rockafellar–Uryasev
on 2·10⁵ Gaussian samples vs the closed form μ + σφ(z_β)/(1−β) within
0.02 (and the RU minimizer matching the VaR); MOO golden hash
`0xaf70_6167_593f_51cc`; MC hypervolume vs exact at m = 2/3 within
0.01 absolute and m = 6 CLOSED FORMS beyond exact reach (single point
∏(ref − p) and two-point inclusion-exclusion, within 8% at 1.2·10⁵
Sobol samples, deterministic per seed with the hit count returned as
the honesty knob); bounded archive — dominated inserts are no-ops,
dominating inserts evict, over-capacity eviction removes the TRUE
least contributor (verified against brute-force keep-(k−1) subset
enumeration); tests/nsga3_battery.rs (4 cases): Das–Dennis counts
(C(p+m−1, m−1) exact: 91 at (3,12), 70 at (5,4)) and on-simplex
membership; DTLZ2(m=3) at standard budgets — worst | ‖f‖−1 | =
0.0238 against the unit-sphere-octant front with 98%
reference-direction coverage; the MANY-OBJECTIVE claim at m = 5:
NSGA-III beats NSGA-II 6.51 vs 4.39 on MC-estimated hypervolume at
matched budget (the classic motivation, measured by composition with
mc_hypervolume) plus bitwise replay; NSGA-III golden
`0xd912_6c49_f1b1_6897`.
tests/dro_battery.rs (4 cases): the one-sample kink LP recovers the
fractional q = [0.5, 0.5] and worst-case value 0.5; tiny-scale kinks
use scale-relative recovery rather than an absolute lambda cutoff;
large radius saturates at max loss; public guards reject
empty/non-finite losses, zero sample counts, invalid radii/costs, and
rows without a zero-cost stay-put support. tests/dro_oracle_battery.rs
(4 cases): closed-form endpoints and monotonicity in rho, exact
tiny-LP strong-duality oracle, robust-decision shift demo, and frozen
DRO golden hash `0xd21c_d092_b4a5_ba98`.

## No-claim boundaries
- No published-ERT-table parity claims yet (in-repo BBOB-class fixtures
  only; the external COCO battery is follow-up).
- NSGA-III normalization uses the ideal point with FIRST-front
  per-objective maxima as the nadir estimate; the full ASF
  extreme-point construction is the recorded refinement. MOEA/D is
  the bead's remaining decomposition-based sibling.
- `mc_hypervolume` is the m > 4 path; its accuracy knob is the
  sample count (standard error √(p(1−p)/n)) — no silent precision
  claim. `HvArchive` eviction uses EXACT contributions (m ≤ 4);
  MC-contribution eviction joins with its many-objective consumer.
- Module `ot`: BALANCED entropic OT only (equal masses asserted);
  unbalanced/partial transport and Sinkhorn divergences (debiasing)
  are follow-up lanes.
- Module `dro`: discrete candidate support only; continuous ambiguity
  sets, adaptive support generation, and coupled decision-dependent
  losses are follow-up scope.
- MOO slice-1 scope (module `moo`): NSGA-II, exact hypervolume m ≤ 4,
  knee, sample-CVaR. NSGA-III reference directions / MOEA/D
  (many-objective), MC hypervolume beyond m = 4,
  hypervolume-contribution archiving, gradient-based Pareto tracing
  (fs-ascent continuation), ledger world-forking steering, and chance
  constraints are the bead's recorded split lanes.
- Sep-CMA/low-rank (dim > ~200), NES, DE, DIRECT, TR-DFO: not built.
- No constraint handling (fs-constraint owns kinds; integration later).
- No parallel evaluation waves yet (fs-exec bead).
