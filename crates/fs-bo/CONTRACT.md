# fs-bo CONTRACT

## Purpose and layer

Layer: **L4 ASCENT** (deps: fs-ascent/fs-dfo L4, fs-la/fs-rand L1,
fs-math L0). Bayesian optimization (plan §9.4 [F]): in-house
Gaussian processes with exact Cholesky inference, deterministic
acquisitions, and BO loops that beat baselines with ledgered
evidence and replay bitwise. Inner optimizers come from the landed
stack: fs-ascent L-BFGS for hyperparameters, fs-dfo CMA-ES for
acquisition surfaces.

## Public types and semantics

- `Kernel`/`Matern` — Matérn ½, 3⁄2, 5⁄2 with ARD lengthscales and
  signal variance; the r → ∞ limit is guarded (degenerate
  lengthscales during hyperparameter search hit inf·0 = NaN in the
  polynomial×exp forms — the guard returns the true limit 0; the
  first draft crashed as a fake NotSpd mid-L-BFGS).
- `Gp` — exact zero-mean GP: `fit`/`try_fit` (Cholesky of
  K + σ_n²I; try_fit makes NotSpd a REJECTABLE hyperparameter
  candidate, not a crash — near-duplicate exploitation points with
  tiny noise are a normal BO condition), `predict` (mean/variance
  via triangular solves), `predict_joint` (posterior mean + jittered
  covariance Cholesky — the q-EI reparameterization substrate),
  `lml` (log marginal likelihood).
- `fit_hyperparams` — maximize LML over log-parameters by fs-ascent
  L-BFGS with FD gradients (≤ D+2 params — cost trivial at fixture
  scale; analytic LML gradients recorded follow-up), QMC-multistarted
  from scrambled-Sobol points (the plan's named detail); clamped
  exponentiation with a 1e−8 noise floor.
- `acq` — closed-form EI; Φ via Abramowitz–Stegun 7.1.26 (~1.5e−7
  absolute) and Φ⁻¹ via Acklam (~1e−9 relative), both DETERMINISTIC
  polynomials with documented accuracy — no platform libm;
  `normal_bank` (scrambled Sobol through Φ⁻¹ — fixed common random
  numbers); `q_expected_improvement` through the Cholesky
  reparameterization f = μ + L·z over the fixed bank.
- `turbo::turbo_minimize` — TuRBO-class trust-region BO (the honest
  answer to BO's dimensionality ceiling): local GP inside an adaptive
  hyperrectangle with ARD-weighted sides, Thompson sampling over
  Sobol candidates through the joint-posterior Cholesky (fixed Philox
  noise — common random numbers, bitwise-replayable), success/failure
  counters (double/halve), restarts on collapse keeping the global
  best; local data capped at `max_local` nearest points and
  hyperparameters refit every `refit_every` iterations (unbounded
  local sets measurably stalled the 30-d battery).
- `mf` — two-fidelity joint GP via the ICM kernel
  K((x,m),(x',m')) = B[m][m']·k_x(x,x'), B = LLᵀ Cholesky-
  parameterized (PSD by construction; the between-fidelity
  correlation is LEARNED and reported); `mf_minimize` = EI on the
  HIGH-fidelity posterior with the MFEI-class fidelity rule
  (evaluate cheap when corr²·cost-ratio > 1, corr from the joint
  posterior at [(x,lo),(x,hi)]); cost-indexed traces ledgered.
- `sparse::{SparseGp, farthest_point_inducing}` — inducing-point
  DTC/SoR GPs with the TITSIAS ELBO as the honesty instrument (the
  trace slack (1/2σ²)·tr(K_XX − Q_XX) reports what the approximation
  discards); inversion-lemma identities through fs-la Cholesky,
  O(n·m²); deterministic farthest-point selection.
- `acq_grad` (feature `tape-acq`): EXACT q-EI gradients by one
  reverse pass through the taped chain (kernels → posterior →
  Cholesky → reparameterized hinge) via fs-ad's FrankenTorch scalar
  bridge — no matrix-level Cholesky backward needed; `qei_ascent` =
  probe-then-polish L-BFGS argmax. DETERMINISM CLASS inherited from
  the bridge: tolerance-verified vs FD, never bitwise, EXCLUDED from
  cross-ISA goldens (the production f64 q-EI stays on det kernels).
- `Gp::try_fit_diag` — HETEROSCEDASTIC fits: per-point noise
  variances K + diag(σᵢ²); the Cholesky/LML/predict paths are
  unchanged once the matrix is built (predictions are latent-f).
- `bo::minimize` — Sobol initialization, per-iteration y
  STANDARDIZATION (EI is affine-invariant when applied consistently;
  without it the signal box cannot span arbitrary objective scales —
  measured: raw-y BO LOST to random on Branin), hyperparameter refit,
  acquisition argmax by CMA-ES restarts, greedy q-EI batch growth
  under common random numbers. Deterministic per seed.

## Invariants

- Kernel matrices are PSD on distinct points (G0-gated across all
  three families); posterior consistency: fitting the union of
  datasets in any order yields the same posterior.
- Acquisition surfaces are deterministic functions (fixed z-banks);
  whole BO runs replay bitwise.
- q-EI is monotone under batch growth (fixed bank) and matches
  closed-form EI at q = 1 to MC tolerance.

## Error model

Structured panics on dimension mismatches, invalid fidelity indices,
invalid TuRBO configuration knobs, and direct `fit` of non-SPD
systems; the hyperparameter search path uses `try_fit` (rejection,
not crash). `phi_inv` asserts p ∈ (0,1).

## Determinism class

Bit-deterministic per seed; golden FNV-64 over GP posteriors,
acquisitions, and a short BO run: `0x5db8_f433_fd71_f738` (bumped
once from `0x4f5a_0601_3cd1_6f46` when predict_joint's flat jitter
became adaptive escalation — required for TuRBO's collapsed trust
regions; semantic justification recorded in the battery); TuRBO
golden `0xe671_9eef_01a1_b960`. Recorded on Apple M4 Pro, verified
on Threadripper (x86_64).

## Cancellation behavior

Iteration-granular: the BO loop is resumable between batches (all
state is the (X, y) history); inner optimizers are the landed
resumable/deterministic engines. Cx wiring is driver scope.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

- `tape-acq` (default OFF) — acquisition gradients through the
  FrankenTorch tape bridge (`dep:fs-ad` + `fs-ad/torch-bridge`); gates
  the `acq_grad_battery` and `probe_grad` integration targets.

## Conformance tests

`tests/bo_battery.rs` (8 cases): kernel PSD across the family × 40
random points; posterior union-order consistency ≤ 1e−8; known-answer
posteriors (noiseless interpolation mean == y and var → 0 at data;
one-point closed form to 1e−12); kernel/GP dimension mismatches fail
fast; Φ table values + Φ⁻¹ round-trip ≤ 1e−6; EI ≥ 0, EI ≈ 0 at
data, q-EI batch dominance + q-EI(1) vs closed-form EI within 5% MC
tolerance; EI-BO beats scrambled-Sobol random search on Branin at
matched budget over a fixed seed set (median 0.80 vs 1.44, optimum
0.398 — the ledgered comparison) with whole-run bitwise replay;
cross-ISA golden hash.
`tests/turbo_battery.rs` (4 cases): trust-region mechanics G0 (a
smooth bowl tracked to 1e−6; a needle objective collapses the TR
into ≥ 1 restart); invalid new TuRBO knobs fail fast; Ackley-30 at
matched 300-eval budget — TuRBO (11.09 median) beats QMC-random
(12.22) with CMA-ES REPORTED alongside (10.22; both are legitimate
high-d optimizers — the ledger records the numbers rather than
cherry-picking the comparison) plus bitwise replay; TuRBO golden.
`tests/mf_battery.rs` (3 cases + dims guard): correlation recovery
(learned ρ ≈ 1.0 on the linear-bias Branin pair) with high-fidelity
variance reduced at 10/10 held-out probes when cheap data is added;
allocation + cost-to-target at 10:1 costs — low-fidelity evaluations
dominate (78% share) and MF-BO reaches 0.3995 median vs
single-fidelity EI-BO's 1.7179 at MATCHED total cost (optimum
0.3979 — the documented win) with bitwise replay; MF golden
`0x6411_f077_1d5e_9f88`.
`tests/sparse_battery.rs` (5 cases): EXACTNESS RECOVERY at Z = X —
predictions match the exact GP to 2.8e−10 (mean) / 1.5e−9 (variance)
and the ELBO is tight to the exact LML at relative 1e−6 (the two
sides come from DIFFERENT factorization paths, so absolute 1e−6 is
not an honest expectation — measured and documented); the ELBO
LOWER-BOUNDS the exact LML at every m (instance-checked, monotone
under nested farthest-point selection); the accuracy ladder —
RMSE-vs-exact 0.114 → 0.019 → 0.003 at m = 10/30/80 on n = 300;
duplicate-row farthest-point ties never reselect an already chosen
row; sparse golden `0x0138_e24a_db84_4bec`.
`tests/hetero_battery.rs` (4 cases): the COINCIDENT-POINT closed
form (2×2 posterior algebra by hand; the precision-weighted mean
sits at the low-noise observation); declared-noisy clusters do not
drag predictions (hetero err 0.011 vs homoscedastic 1.14 on the same
data — measured); E-RACING candidate elimination via fs-eproc
confidence sequences — the optimum found with 49% fewer samples than
uniform max replication AND every stopped CS covering its true mean
at the stopping time (the Bet-5 validity claim, instance-checked);
two fixture bugs documented in-test (half-width stopping is
mean-independent for fixed-σ CSs, making cross-candidate comparisons
vacuous; Hoeffding σ = ½ was 3× conservative vs the true noise —
clamping is a contraction so the actual σ is valid); hetero golden
`0xe9b3_f6b5_69ee_258b`.
`tests/acq_grad_battery.rs` + `tests/probe_grad.rs` (feature
`tape-acq`, aarch64 evidence; cross-ISA not claimed per the bridge's
determinism class): taped q-EI primal parity with the production
path at 1.3e−7 relative; gradient vs FD worst 2.0e−8 WITH ALIVENESS
GUARDS (the first gate passed vacuously on a flat region while the
Σ-diagonal's √0 backward poisoned gradients with NaN — Matérn-5⁄2 is
C² so the r = 0 kernel value is a CONSTANT with exactly zero
gradient, and the probe regression pins an active surface forever);
probe-then-polish ascent improves its seed 0.314 → 0.374 in 47
reverse passes (FD-equivalent 423) with CMA-ES REPORTED (0.488 — a
Monte-Carlo needle spike at near-duplicate candidate blocks; gating
local polish on a needle hunt would be dishonest in both
directions).
`tests/mf_battery.rs`: multi-fidelity correlation recovery and
variance reduction, dimension/fidelity mismatch fail-fast behavior,
cost-aware allocation and replay, and MF golden.
`tests/bo_study_replay.rs` (3 Casebook cases plus a seeded red
self-test): runs the production sequential-EI loop on the short Branin
golden configuration (`n_init = 6`, two one-point iterations), records
all eight objective callbacks and every public `BoReport` bit, and
independently reconstructs the callback/report correspondence,
objective values, and three-point best-so-far trace. The complete
canonical output frame replays exactly under the same configuration,
seed, crate-version set, RNG stream-semantics version, and Casebook
record version. A disclosed `StreamKey` selects one initial reported
coordinate and one low mantissa bit, changes only that result bit, and
reseals the payload: payload validation still succeeds, while the
reference-output identity gate returns a typed mismatch. The red
Casebook record is byte-stable across regeneration and its
`assert_green` refusal is caught and retained before the final green
three-case report.

## No-claim boundaries

- All bead lanes landed (TuRBO, two-fidelity ICM, DTC/Titsias
  sparse GPs, heteroscedastic fits, e-racing stopping, tape
  acquisition gradients). The taped path is feature-gated and
  fixture-scale (dense training view); production-scale tape argmax
  joins with a consumer needing it. Heteroscedastic noise is
  CALLER-SUPPLIED (from fs-uq estimates); joint noise-model learning
  is out of scope. Sparse inducing LOCATIONS are fixed
  (farthest-point); variational location optimization joins with a
  consumer at genuine 10⁴-point scale. The MF module is
  TWO-fidelity; deeper ladders and per-evidence-model-ledger
  discrepancy models join with their consumers.
- TuRBO's joint fallback drops cross-correlations only where the
  posterior is numerically void (degenerate candidate sets); the
  30-d battery budget is sized for the debug profile — the claim is
  comparative at matched budget, not budget-specific.
- Acquisition gradients are derivative-free (CMA-ES); the
  FrankenTorch reparameterized-gradient tape through q-EI is the
  named follow-up (the fixed-bank surfaces are already
  differentiable-by-construction when the tape lands).
- q-NEI (noisy EI) is not implemented — EI with the standardized
  noise floor covers the deterministic-objective regime shipped
  here.
- No hyperparameter marginalization (point estimates by LML); no
  input warping.
- The short-study replay receipt covers one finite Branin
  configuration and same-process same-ISA execution. It makes no
  optimizer-quality, all-objective, all-configuration, all-seed,
  cross-ISA, cancellation/`Cx`, persistence, authenticated-ledger, or
  performance claim. Its FNV identity and local reference gate are
  evidence-fixture plumbing, not a cryptographic artifact seal or a
  production admission service.
