# fs-uq CONTRACT

## Purpose and layer

Layer: **L4 ASCENT** (deps: fs-bo/fs-robust L4, fs-la/fs-rand L1, fs-math L0).
Uncertainty quantification (plan §8.8): random-field inputs,
spectral and sampling propagation, and multilevel Monte Carlo.
Layered at L4 deliberately: propagation WRAPS solvers the way
ASCENT optimizers do, its risk outputs feed ASCENT's robust
formulations, and the placement reuses fs-bo's deterministic Φ⁻¹
for QMC Gaussian germs instead of duplicating the polynomial.

## Public types and semantics

- `KlExpansion` — Karhunen–Loève by dense Jacobi eigendecomposition
  of the covariance at a point set (exponential and
  squared-exponential families), truncated at a captured-variance
  target with the fraction ACTUALLY captured reported (the evidence
  the bead requires); `realize` (germ → field), `qmc_germs`
  (scrambled Sobol through Φ⁻¹ on the leading ≤ 10 variance-dominant
  modes — the embedded Joe–Kuo table's cap, larger tables being
  fs-rand's recorded follow-up — with Philox normal tails),
  `covariance_reconstruction_error` (the truncation-quality audit).
- `pce::{fit_pce, PceModel, hermite_orthonormal}` — polynomial chaos
  by regression: orthonormalized probabilists' Hermite basis,
  graded-lex total-degree truncation, least squares via Cholesky
  normal equations (ridge 1e−12) with an oversampling assertion
  (n ≥ 2·basis); mean/variance drop out of coefficients.
- `mlmc::mlmc_estimate` — multilevel Monte Carlo with pilot variance
  estimation and the Giles allocation N_ℓ ∝ √(V_ℓ/C_ℓ); the sampler
  contract makes COUPLING explicit (level ℓ returns the correction
  Y_ℓ driven by one germ — that coupling is what makes variance
  decay); per-level statistics ledgered in the report.

- Slice 2 (bead o5kc): `seismic` — `KanaiTajimi` PSD + spectral-
  representation synthesis with counter-based phases (bit-replayable
  per seed), `sdof_peak` (Newmark average-acceleration response
  ordinates), `cqc` (Der Kiureghian coefficients; collapses to `srss`
  for separated modes, exceeds it for close modes),
  `bilinear_peak_ductility` + `ida_fragility` (the nonlinear IDA
  workhorse; monotone exceedance curves). `anytime` —
  `estimate_probability_anytime`: a sub-Gaussian (σ = 1/2) mixture
  confidence sequence stops the moment its RADIUS meets the decision's
  half-width — valid AT the stopping time by construction (fs-eproc's
  `interval()` returns (center, radius); mis-reading it as (lo, hi)
  produced 58% miss rates in development — the semantics are
  load-bearing). `cvar` and `empirical_cvar` are direct re-exports of
  `fs-robust`, the canonical ASCENT risk algebra; the report includes CVaR,
  lower empirical VaR/Rockafellar–Uryasev minimizer, boundary rank, and
  fractional boundary weight. `adaptive` — `adaptive_mlmc`: levels added while
  the extrapolated bias `|mean_L|/(2^α − 1)` exceeds tol/2, with weak/
  strong rates FITTED from level statistics.

- `chance` module (bead qlvf, lane b): `chance_constrained_min` —
  `P(g(x,ξ) ≤ 0) ≥ 1−α` enforced through probability estimates that are
  ANYTIME-STOPPED (each query ends the
  moment the CS is decision-grade). The complete bounded query family is
  predeclared and receives a Bonferroni allocation of the familywise budget;
  adaptive selection may return only a candidate whose simultaneous CS lower
  bound is feasible. If none is admitted, the API returns
  `NoFeasibleChanceCandidate` rather than a penalized but uncertified design.

## Invariants

- The KL truncation reports its captured-variance fraction; the
  reconstruction error audit bounds what was dropped.
- PCE bases are orthonormal (mean = c₀, variance = Σc²_{≥1} are
  exact identities, not approximations).
- The MLMC estimate is exactly the sum of level means (telescoping
  bookkeeping is an identity, tested to 1e−14).
- Everything is deterministic per seed.

## Error model

Structured panics on dimension mismatches, under-sampled PCE
regressions, and non-SPD normal equations. Statistical quality is
REPORTED (captured variance, reconstruction error, per-level
variances), never assumed. The canonical CVaR re-exports instead return
structured `RobustError` values for empty samples, invalid tail levels, and
non-finite losses; they do not panic.

## Determinism class

Bit-deterministic per seed; golden FNV-64 over KL spectra/fields and
PCE coefficients: `0x0ed2_4974_dc37_bbc6`, recorded on Apple M4 Pro,
verified on Threadripper (x86_64). CROSS-ISA by construction: every
transcendental routes through `fs_math::det` (bead frankensim-lyms;
the golden's KL/PCE flow was already libm-free, and the adaptive-MLMC
rate fit and seismic paths now route det:: too, so no path depends on
platform libm). Registered in the `check-libm` doctrine lint. log2 is
computed as `det::ln(x) / LN_2` (det has no log2); `sqrt`/`powi(±3)`
stay primitive (IEEE correct rounding / pinned lowering).

## Cancellation behavior

Sampling loops are bounded (pilot + allocation); MLMC state is the
per-level running sums — iteration-granular by construction. Cx
wiring is driver scope.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

`tests/uq_battery.rs` (5 cases): KL on a 49-point grid — captured
variance ≥ 0.99 with hard truncation (< half the points for the
smooth family), covariance reconstruction ≤ 5% Frobenius, and
QMC-sampled field variance at a probe matching the retained-diagonal
target within 5%; PCE known-answer — Y = exp(a·ξ) coefficients match
the closed form e^{a²/2}aᵏ/√k! (six coefficients), mean/variance to
1%/5%, surrogate max error < 0.02 on fresh points; QMC-vs-MC
MEASURED — scrambled Sobol beats plain MC by > 3× RMSE at n = 2048,
d = 5 over 20 replicates (the workhorse claim, ledgered); MLMC on a
synthetic ladder with classic decay rates — telescoping identity to
1e−14, estimate within 5% of the closed-form target, allocation
favoring coarse levels > 10×, and cost win > 3× vs single-level MC
at MATCHED estimator variance; cross-ISA golden hash.

`tests/slice2.rs`: UQ-side CVaR re-exports match the canonical `fs-robust`
report for fractional and tied integer boundaries, mixed-sign extremes,
constant `f64::MAX`, and input permutations; invalid risk inputs return
structured errors rather than panicking.

## No-claim boundaries

- Slice-1 scope. The bead's split lanes: seismic machinery
  (Kanai–Tajimi spectra via fs-scenario, response-spectrum CQC fast
  path, IDA fragility curves — needs fs-solid-advanced), e-process
  anytime-valid stopping (fs-eproc integration — every estimate
  under confidence sequences), CVaR/quantile risk functionals
  feeding ASCENT robust formulations, and adaptive MLMC (level
  addition on bias estimates; the ladder here is caller-fixed).
- KL is dense-eigen (fine to ~10³ points); Nyström/randomized
  scaling rides fs-la-randomized when a consumer needs it.
- PCE is regression-only (no sparse grids/projection quadrature);
  sparse-grid construction joins with its consumer.
- Sobol' sensitivity indices are derivable from the PCE coefficients
  but not yet exposed as API.

## No-claim boundaries (slice 2)

- The IDA workhorse is a bilinear SDOF: the fiber-frame nonlinear
  time-history path arrives with fs-solid-advanced (tfz.14, in
  flight) — the fragility MACHINERY is what this slice ships.
- Ground motions are Kanai–Tajimi synthetic; recorded-suite ingestion
  rides fs-scenario.
- The CS applies to BOUNDED [0,1] outcomes via Hoeffding
  sub-Gaussianity; release and debug builds both refuse non-finite or
  out-of-range samples before updating the process. The target half-width must
  be finite and non-negative. Unbounded losses need a declared sigma.
- `cvar` is the standard finite-sample empirical CVaR (Acerbi–Tasche): the mean
  over the worst `n·(1−β)` fraction of losses, with a FRACTIONAL weight on the
  boundary order statistic when `n·β` is not an integer. A plain top-`⌈n(1−β)⌉`
  mean under-reports the shortfall (anti-conservative — the wrong direction for a
  risk measure); the fractional weight makes the estimator honest. The
  canonical `fs-robust` implementation uses bounded convex combinations so
  finite constant and mixed-sign extreme samples do not overflow intermediate
  sums. Empty samples, non-finite samples, or a tail level outside `(0,1)`
  return `RobustError`; returning zero risk for an empty loss set would be a
  false certificate.
- Adaptive-MLMC rate fits assume geometric level decay; oscillatory
  correction means defeat the extrapolation (documented, as in the
  classical literature).
