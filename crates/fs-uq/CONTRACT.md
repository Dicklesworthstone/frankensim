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
- `product_pce::{PceConfig, PceExecution, PceReport, PceSummary}` — fixed
  Hermite/shifted-Legendre product regression from physical-model callbacks,
  with disjoint held-out physical validation. `ProductPceModel` retains each
  independent germ law alongside its coefficients, sharing Hermite and ANOVA
  algebra with the unchanged standalone `PceModel`. Execution retains
  product admission and bounded, scaled `fs-la` QR least squares. Numerical
  rank/conditioning are screened before the first physical solve; no ridge
  changes the declared least-squares problem.
- `mlmc::mlmc_estimate` — multilevel Monte Carlo with pilot variance
  estimation and the Giles allocation N_ℓ ∝ √(V_ℓ/C_ℓ); the sampler
  contract makes COUPLING explicit (level ℓ returns the correction
  Y_ℓ driven by one germ — that coupling is what makes variance
  decay); per-level statistics ledgered in the report.
- `product_lhs::LhsExecution` executes independent replicated Latin hypercubes.
  Its `freeze_linear_control_variate` and `assess_linear_control_variate`
  methods share the existing physical marginal means and centered-control
  arithmetic. The complete LHS plan and layout bind the frozen coefficients.
  Assessment replays only sampler coordinates, uses complete hypercubes as the
  independent statistical units, and leaves raw observations, compliance and
  checkpoint bytes unchanged. Partial replicates are retained work, excluded
  from both means. Errors are descriptive between-replicate estimates, not
  confidence intervals, optional-stopping guarantees or model-error bounds;
  variance reduction is measured rather than guaranteed. Assessment polls
  cancellation during permutation replay, centering and reductions.

- Slice 2 (bead o5kc): `seismic` — `KanaiTajimi` PSD + spectral-
  representation synthesis with counter-based phases (bit-replayable
  per seed), `sdof_peak` (Newmark average-acceleration response
  ordinates), `cqc` (Der Kiureghian coefficients; collapses to `srss`
  for separated modes, exceeds it for close modes),
  `bilinear_peak_ductility` + `ida_fragility` (the nonlinear IDA
  workhorse). Both oscillator paths use implicit Newmark average acceleration
  from initial rest; the first acceleration sample initializes equilibrium,
  so N samples span N−1 intervals and empty/single-sample records return zero.
  The bilinear path solves its monotone piecewise-linear equilibrium in at
  most two scalar solves per interval, without allocating or iterating Newton;
  plastic history changes only after the accepted branch is solved. It uses
  a parallel elastic and elastic-perfectly-plastic spring with post-yield
  stiffness ratio in `[0,1]`. CQC preserves signed modal contributions and
  exact unit diagonal correlation; distinct undamped frequencies are
  uncorrelated, while identical undamped frequencies combine coherently.
  Pair symmetry halves coefficient work and ordered frequency ratios avoid
  overflow across wide bands. `anytime` —
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
The scalar seismic response/CQC fronts also assert their finite physical
parameters, array dimensions and representable arithmetic; NaNs must not
silently disappear into a seemingly valid peak response.

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

`seismic` unit tests: implicit bilinear dynamics match the linear limit for
undamped, critical and overdamped oscillators; a short-period oscillator stays
within its analytic constant-load amplitude bound at `omega*dt = 5`; a yielded
step solves its accepted restoring force; cyclic plastic peaks converge to an
independent 600,000-step explicit reference. CQC checks zero damping, signed
degeneracy, unequal damping, frequency-unit changes and extreme peak scales.

### Ground-motion bandwidth and sample clock

`GroundMotionSampling` fixes the physical cutoff in rad/s, midpoint frequency
count, sample interval and record length. `synthesize_with_sampling` realizes
the one-sided angular-frequency density: each line carries variance
`S(omega_k)*delta_omega`, so its random-phase cosine amplitude is
`sqrt(2*S(omega_k)*delta_omega)`. A per-Hz density is `2*pi*S(2*pi*f)`.
Every line is strictly below the sampling Nyquist frequency. An explicit band
that exceeds Nyquist is refused; high-frequency power is neither folded into
lower structural modes nor renormalized into the retained band.

The legacy `synthesize` convenience call uses the visible `default_sampling`
choice `omega_max=min(8*wg,pi/dt)`. This is a finite-band approximation, not
the full untruncated Kanai-Tajimi variance. The former dimensional `wg.max(1)`
floor is absent. Keeping the explicit bandwidth, frequency count and seed fixed
preserves common samples under clock refinement and preserves prefixes when
extending a record. Frequency resolution remains a quadrature choice: narrow
soil resonances require sufficient lines. The model requires positive finite
ground frequency/damping and finite nonnegative intensity.

`ida_fragility_with_sampling` routes this exact record clock and bandwidth into
the existing bilinear time-history consumer, preparing spectral amplitudes once
for the suite and reusing records across intensity levels. The original IDA
entry point remains the 0.01-s, 1500-sample, 96-line convenience fixture.
Regression tests independently reconstruct Fourier-line powers, finite-period
variance/covariance and rad/s-to-Hz normalization; check fixed-band clock and
record refinement, physical unit/intensity scaling, and the actual IDA route.

## No-claim boundaries

- Newmark's linear stability does not remove temporal accuracy requirements:
  callers must resolve ground-motion bandwidth, structural periods and yield
  reversals. `ida_fragility` retains its explicit fixture clock of 0.01 s;
  `ida_fragility_with_sampling` accepts a caller-refined clock and fixed physical
  bandwidth. Observed monotonicity is not a theorem for every nonlinear suite.
  The zero-damping CQC convention is the coherent/uncorrelated modal limit,
  not a finite stationary variance for undamped white-noise excitation.

- `product_plan::UqPropagator` currently executes only fixed-count Philox Monte
  Carlo with 1..=256 parameters and 2..=1,000,000 samples. Independent Gaussian
  and uniform marginals require explicit independence when more than one input
  has nonzero uncertainty; the default dependence is Unknown. Deterministic
  inputs do not create a dependence requirement. `JointGaussian` explicitly
  declares the joint law; Gaussian marginals plus a generic `Correlated` matrix
  do not suffice and refuse. The joint law is sampled through a numerical PSD
  Cholesky factor, including consistent singular
  pivots (roundoff tolerance `64*EPSILON*dimension`). Finite symmetry, unit
  diagonal and correlation range are checked before sampling. This numerical
  factorization is not a certified PSD proof or authority for a general copula.
  Parameters follow declaration order; each logical sample uses Philox kernel
  `0x0517`, tile equal to the sample index. This replaces the earlier single
  sequential stream, so old product-plan sample sequences are not preserved.
  Result hashes use the `org.frankensim.uq.result.v2` domain and include every
  result field without the former six-decimal truncation of observed ranges.
- The single-level product-plan driver refuses QMC/PCE/MLMC/epistemic bounding,
  confidence-band sampling, unstated inputs and non-Gaussian correlated models
  before model calls.
  A confidence half-width is not a probability measure. Units are required
  labels; dimensional compatibility and physical applicability remain caller
  obligations. Non-finite model values stop with an explicit refusal and the
  actual evaluation count. Empirical ranges, percentiles, compliance frequency
  and mean standard error are always Estimated; none is a guaranteed interval,
  anytime confidence sequence, PDE error bound or model validation. There is
  no cancellation/resume or Cooling project driver in this API yet.
  `tests/product_plan.rs` checks analytic correlated-Gaussian variance,
  singular dependence, pre-callback refusals, finite statistics, and authority.
- `product_mlmc::MlmcExecution` is the method-specific MLMC producer for
  1..=64 caller-declared resolution levels with at least two corrections per
  level. It reuses the existing Gaussian/uniform probability-law admission,
  numerical joint-Gaussian factor and physical Philox sampler. The immutable
  plan budget equals `N0 + 2*sum(Nl, l>0)` physical model calls and remains
  capped at one million. Level-major corrections have disjoint sample ordinals;
  every fine/coarse pair receives the identical physical parameter vector.
  A completed fine half-pair survives chunk boundaries, cancellation and clone
  resumption without another fine solve. An interrupted unfinished call retries
  its same level and parameters; its internal cost belongs to the solver.
  Errors and nonfinite physical values/corrections are terminal, with no sample
  skipping or replacement. Completed corrections and call counts remain visible.
  The global estimate appears only after the entire fixed allocation and is
  `sum(mean(Yl))`, with descriptive sampling variance `sum(s_l^2/Nl)`.
  Scaled moment arithmetic is shared with the independent-replicate QMC owner.
  A singleton partial level has undefined dispersion; unrepresentable variance
  remains absent even when its standard error is representable. Optional
  compliance applies the indicator to each physical level before differencing.
  Its finite-sample telescoping probability estimate can leave `[0,1]` and is
  not clipped. Correction samples are not a physical-QoI distribution and
  produce no pooled quantiles or CVaR. All reports remain Estimated. There is
  no automatic level/allocation adaptation, bias/continuum/model-error bound,
  anytime stopping rule, hidden-randomness coupling or durable checkpoint format.
  `tests/product_mlmc.rs` covers a backward-Euler cooling ladder with an exact
  finest-level expectation, independent reconstruction of correction statistics,
  one-level MC equivalence, correlated inputs, paired-call continuation,
  compliance, admission failures and terminal numerical/model refusals.
- `product_pce::PceExecution` admits independent finite Gaussian/uniform
  mixtures or an admitted joint Gaussian law, preserving unknown dependence
  and non-Gaussian correlated-law refusals. Standard-normal and U(0,1) germs
  use the existing Philox seed/kernel/ordinal sampler in declaration order.
  The callback receives the declared physical transformation, also exposed by
  `physical_parameters`; singular covariance and constant marginals require no
  inverse. Uniform values use bounded convex endpoint interpolation, not a
  potentially overflowing interval width. Germs outside [0,1] refuse.
  The envelope is 1..=16 declared parameters, degree 0..=8, <=128 total-degree
  terms, <=65,536 physical calls, and <=16,777,216 training row/term/term
  products. At least twice the basis count is reserved for training; at least
  two disjoint later ordinals remain for held-out validation. Numerical pivot
  checks and the existing `fs-la` R-factor 1-norm condition estimate screen the
  scaled design before callbacks. This floating-point screen is not a certified
  condition-number bound. Uniform inputs do not silently become Gaussian.
  Every paid finite observation survives chunks, cancellation and clone resume.
  Interrupted unfinished callbacks retry identical inputs without consuming an
  ordinal. Model/sampling/fit failures are terminal; no data is replaced.
  Only the full allocation exposes the fitted `ProductPceModel` and summary.
  Its mean, variance and gradients belong to the product surrogate; coefficient variance
  is response variance, not an estimator standard error or confidence interval.
  Unrepresentable variance is absent. Training and held-out RMS residuals and
  held-out mean error are descriptive checks against physics, not error bounds.
  Optional compliance uses held-out physical indicators only, with an empirical
  frequency and descriptive standard error; it does not evaluate the surrogate
  at the threshold or reuse the fitting data. Germ Sobol indices cannot be
  relabeled physical-coordinate indices under correlated inputs. Evidence stays
  Estimated with infinite dispersion because no total error bound is supplied.
  Cancellation polls bracket physical calls and final fitting; bounded design
  construction, QR and diagnostics remain synchronous. No durable checkpoint,
  degree adaptation, automatic method selection or anytime stopping is claimed.
  `tests/product_pce.rs` covers analytical thermal-product and cooling-ODE
  moments, usable surrogate derivatives, joint/singular Gaussian MC sample
  identity, held-out separation, zero uncertainty, chunk/interruption/fit-boundary
  continuation, design admission and terminal physical/numerical refusals.
  `tests/product_pce_model.rs` adds shifted-Legendre orthogonality, endpoint
  derivatives, mixed-law coefficient sensitivity and Gaussian arithmetic
  identity. `tests/product_pce_mixed.rs` adds uniform/mixed physical moments,
  ordinary-MC draw identity, held-out separation and paid-work continuation.
  `MIXED_PCE.md` documents use and the return-type/accessor migration. These
  new tests require native execution; their presence is not a runtime pass.
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
- `PceModel::sobol_indices` exposes algebraic independent-germ sensitivity
  indices of the supplied surrogate, with no physical-model accuracy guarantee.

## No-claim boundaries (slice 2)

- Direct-model `SobolExecution` evaluates fixed Philox A/B/hybrid designs under
  explicitly independent Gaussian/uniform physical inputs. Main and total
  effects are Jansen variance-ratio estimates and include no confidence interval
  or optional-stopping guarantee. Only a complete predeclared design exposes
  estimates; finite-sample indices are not clipped, and zero sampled variance
  remains undefined rather than implying zero sensitivity.
- `product_shapley::ShapleyExecution` attributes physical QoI variance to
  1..=8 nonconstant Gaussian inputs under explicit Independent or JointGaussian
  dependence. The original Monte Carlo plan budget must equal
  `N*2*(2^d-1)`, with N>=2 and at most one million physical calls. Every proper
  coalition S receives a permuted copy of the original law, admitted under the
  existing numerical PSD and singular-pivot policy before physical work. The
  conditioned coordinates come first in that factor; sharing their normal
  prefix between two otherwise independent draws freezes the physical X_S
  exactly. Callback inputs return to original parameter order and physical
  units. No covariance inversion, physical-coordinate swap, jitter or hidden
  independence assumption substitutes for conditional sampling. Constant
  physical parameters remain in the caller's fixed model so conditioning on
  them cannot reveal an unobserved latent variable.
  Half the mean squared conditional-pair difference estimates
  `c(S)=E[Var(Y|X_S)]`; the empty subset estimates output variance and the full
  subset has zero conditional variance under the fixed deterministic-model
  contract. Each physical effect sums `c(S)-c(S union {i})` with weights
  `1/(d*choose(d-1,|S|))`, normalized by c(empty). These Shapley effects share
  both dependence and interaction contributions; they are not causal effects
  or independent-input Sobol indices. Conditional variance fractions and
  finite-sample effects remain unclipped even when sample noise violates
  monotonicity or [0,1]. The effect sum is an algebraic identity up to rounding,
  not evidence of estimator accuracy. Zero/unrepresentable sampled variance
  leaves effects unavailable after completion, never a false zero sensitivity.
  Shared Sobol compensated reductions and scaled/centered normalization retain
  small changes on large offsets and handle opposite finite extremes.
  Paid calls, including one half of a pair, survive chunks, cancellation and
  clone resumption. Unfinished callbacks retry identical inputs; model/sampling
  failures are terminal and never skip a draw. Only the full allocation exposes
  effects. Evidence stays Estimated with infinite dispersion; no confidence
  interval, optional-stopping rule, durable checkpoint or adaptive subset search
  is claimed. `tests/product_shapley.rs` checks analytical correlated and
  interaction effects, physical conditioning and base-MC identity, independently
  reconstructed variance costs, permutation/affine behavior, paid-call replay,
  admission, undefined variance, unclipped effects and terminal failures.
- `product_importance::ImportanceExecution` executes physical rare-event
  failure probabilities under explicit Independent/JointGaussian Gaussian laws,
  including singular covariance. A fixed, prechosen independent-germ shift mu
  defines q(z)=epsilon*phi(z)+(1-epsilon)*phi(z-mu), with 0<epsilon<=1 and finite
  reciprocal. The separate Philox component stream preserves original normal
  addressing. Zero shift or epsilon=1 retains ordinary-MC input bits and unit
  weights. Physical units, admitted factor and the fixed QoI ceiling are retained.
  Failure means QoI>ceiling. Log-sum-exp yields target/proposal weights; invalid
  configurations and unrepresentable draw weights refuse before model calls.
  Ordinary weighted means are neither self-normalized nor clipped to [0,1].
  Scaled reductions expose descriptive standard error, weighted-failure ESS and
  mean likelihood weight; these diagnostics are not confidence guarantees.
  `assess_failure_probability` sends epsilon*weight*failure in [0,1] to the
  existing sigma=1/2 Gaussian-mixture confidence sequence and rescales by
  epsilon; proposals equal to the target retain the ordinary unit-indicator CS.
  This range-based anytime bound can be much broader than empirical
  rare-event error. Clipping probability bounds never creates convergence.
  Fix proposal/model/threshold/alpha independently of these outcomes; population
  transfer, finite-grid/rounding and solver/model error are not covered. No
  outward certificate, adaptive proposal fitting or durable checkpoint is minted.
  Budgeted chunks, cancellation and cloning retain paid observations; unfinished
  calls retry identical inputs/weights. Numerical/model failures are terminal
  and cannot expose a selected-prefix estimate or confidence. Evidence remains
  Estimated. `tests/product_importance.rs` covers analytic five-sigma events,
  physical correlations/singularity, independently reconstructed likelihood
  moments and CS, MC identity, paid-call resumption and terminal refusals.
- `SobolExecution::checkpoint` / `restore` retain exact ordered observations,
  including partial rows, in bounded `FSSOB001` bytes. Identity binds every
  `UqPlan` field, original evaluation cap, Philox addressing, A/B/hybrid pairing
  and caller-supplied evaluator/model identity. Restore regenerates only sampler
  inputs, never re-evaluates completed physical calls, and reconstructs all
  statistics from retained observation bits. Errors/nonfinite observations are
  terminal and cannot yield resumable prefixes. Completed checkpoints remain
  terminal even when sensitivity normalization is undefined. Checksums detect
  corruption but do not authenticate or validate the producing model; callers
  must bind all fixed inputs and evaluator semantics and trust the source.

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
- `docs/CERTIFICATE_REGIMES.md` routes long-horizon means, broadband
  observables, and duty-cycle reliability to statistical evidence with the
  exact sampling/dependence/stopping and model-form context. Current UQ
  machinery is thin or staged for those complete routes; neither one simulated
  trajectory nor an empirical estimator alone validates the model, regime,
  population transfer, or rare-event tail.
