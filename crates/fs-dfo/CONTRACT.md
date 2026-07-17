# CONTRACT: fs-dfo

> Status: PARTIAL — CMA-ES (IGO form), BIPOP restarts, Nelder–Mead,
> balanced entropic OT, NSGA-II/III, MOEA/D, hypervolume helpers and
> archives, canonical sample-CVaR delegated to `fs-robust`, and
> discrete-support Wasserstein-DRO inner sup are in force; sep/low-rank CMA, NES
> parameterizations, DE, DIRECT, TR-DFO, and fs-exec population waves
> are recorded follow-up scope.

## Purpose and layer
Derivative-free optimization engines (plan §9.3). Layer: **L4 ASCENT**.
Deps: fs-rand (keyed sampling), fs-la (eigendecomposition), fs-math, fs-robust
(canonical empirical risk algebra).
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
- `empirical_cvar(samples, alpha) -> Result<EmpiricalCvarReport, RobustError>`
  — a direct re-export of the canonical `fs-robust` finite-sample risk algebra.
  The report carries overflow-safe CVaR, the deterministic lower empirical
  VaR/Rockafellar–Uryasev minimizer, boundary rank, and fractional boundary
  weight; DFO does not maintain a second order-statistic implementation.

- `steer` module (bead qlvf, lane a): WORLD-FORKING steering (P9).
  `SteeredStudy` = the deterministic (population, stream-index,
  weights) triple + a base seed + the steering LINEAGE; `fork` never
  mutates the parent and records a ledger-ready `SteerEvent`;
  `advance` is a pure function of (seed, stream-index, weights) so
  every branch replays bitwise from its lineage (`fingerprint` is the
  witness). Chance constraints live in `fs_uq::chance` — fs-uq sits
  above fs-dfo through fs-bo, so the integration points that way (the
  dependency cycle the first draft hit is the layer diagram talking).

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
failure is DATA: `converged: false` with best-found + diagnostics. Empirical
CVaR instead returns the canonical structured `RobustError` and does not panic.

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
`0x58eb_8443_224c_a689`; tests/moo_battery.rs: hypervolume vs
hand-computed 2D/3D values including dominated/out-of-reference
degenerate cases; non-dominated-sort front assignment exact;
NSGA-II on ZDT1/ZDT2 at standard budgets (pop 80 × 200 generations —
short runs measurably leave the f2-minimal arm unexplored, documented
in the test) with mean front gap ≤ 0.05 (measured ≤ 0.0008), full f1
spread, hypervolume beating scrambled-Sobol random at MATCHED
evaluations (0.87 vs 0.19 / 0.54 vs 0.00), and bitwise replay; helper
edges for crowding; knee detection hitting a synthetic elbow exactly;
canonical empirical CVaR on 2·10⁵ Gaussian samples vs the closed form
μ + σφ(z_β)/(1−β) within 0.02 (and the lower RU minimizer matching the
VaR); fractional boundaries, tied integer boundaries, mixed-sign extremes,
constant `f64::MAX`, permutations, and direct `fs-robust` parity; MOO golden
hash `0x606f_35d4_bfb8_822a` remains frozen until the canonical accumulation is
run and any bit change is deliberately justified; MC hypervolume vs exact at
m = 2/3 within
0.01 absolute and m = 6 CLOSED FORMS beyond exact reach (single point
∏(ref − p) and two-point inclusion-exclusion, within 8% at 1.2·10⁵
Sobol samples, deterministic per seed with the hit count returned as
the honesty knob); public MC-hypervolume edges match exact
hypervolume's malformed-point ignore policy, empty dimension returns
zero, and zero samples are rejected; bounded archive — dominated
inserts are no-ops, dominating inserts evict, over-capacity eviction
removes the TRUE least contributor (verified against brute-force
keep-(k−1) subset enumeration); tests/nsga3_battery.rs (11 cases):
Das–Dennis counts
(C(p+m−1, m−1) exact: 91 at (3,12), 70 at (5,4)) and on-simplex
membership; public guards reject zero-division directions, invalid
NSGA-III reference directions, empty MOEA/D weights, and zero
neighborhoods at the API boundary; DTLZ2(m=3) at standard budgets —
worst | ‖f‖−1 | = 0.0238 against the unit-sphere-octant front with
98% reference-direction coverage; the MANY-OBJECTIVE claim at m = 5:
NSGA-III beats NSGA-II 6.51 vs 4.39 on MC-estimated hypervolume at
matched budget (the classic motivation, measured by composition with
mc_hypervolume) plus bitwise replay; MOEA/D converges on ZDT1 and is
competitive with NSGA-III on DTLZ2(m=3) at matched budget (2.7457 vs
2.7775 hypervolume in the current deterministic battery) plus bitwise
replay; NSGA-III golden
`0xd912_6c49_f1b1_6897`.
tests/dro_battery.rs (4 cases): the one-sample kink LP recovers the
fractional q = [0.5, 0.5] and worst-case value 0.5; tiny-scale kinks
use scale-relative recovery rather than an absolute lambda cutoff;
large radius saturates at max loss; public guards reject
empty/non-finite losses, zero sample counts, invalid radii/costs, and
rows without a zero-cost stay-put support. tests/dro_oracle_battery.rs
(4 cases): closed-form endpoints and monotonicity in rho, exact
tiny-LP strong-duality oracle, robust-decision shift demo, and frozen
DRO golden hash `0xd21c_d092_b4a5_ba98`. tests/steer.rs has two logged
conformance cases plus one empty-population guard regression: a fork leaves
its parent untouched while both branches replay bitwise, opposite weights
steer sibling forks toward opposite Pareto extremes, and zero population is
rejected at construction. Aggregate results use canonical fs-obs
`ConformanceCase` events with their exact input seeds; fork and steering
measurements use separate validated fs-obs `Custom` events that also retain
those seeds. Assertions precede each passing aggregate verdict, so earlier
failures remain ordinary Rust test diagnostics.

## No-claim boundaries
- No published-ERT-table parity claims yet (in-repo BBOB-class fixtures
  only; the external COCO battery is follow-up).
- NSGA-III normalization uses the ideal point with FIRST-front
  per-objective maxima as the nadir estimate; the full ASF
  extreme-point construction is the recorded refinement.
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
- MOO scope (module `moo`): NSGA-II/III, MOEA/D, exact hypervolume
  m ≤ 4, MC hypervolume beyond m = 4, bounded hypervolume archive,
  and knee detection. Canonical sample-CVaR is re-exported from `fs-robust` at
  the crate root. Gradient-based Pareto tracing (fs-ascent
  continuation), ledger world-forking steering, and chance constraints
  are still split lanes.
- Sep-CMA/low-rank (dim > ~200), NES, DE, DIRECT, TR-DFO: not built.
- No constraint handling (fs-constraint owns kinds; integration later).
- No parallel evaluation waves yet (fs-exec bead).
