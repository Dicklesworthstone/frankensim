# Optimizing an objective with declared observation noise

`fs_bo::noisy::minimize_noisy` runs sequential or batched Bayesian optimization
without treating the smallest noisy observation as an exact incumbent. Each
callback returns a `NoisyObservation { value, noise_variance }`. The variance
belongs to that observation: when `value` is an average of independent draws,
use the variance of the average rather than the raw sample variance.

```bash
cargo run -p fs-bo --example noisy_design
cargo test -p fs-bo --lib noisy
cargo test -p fs-bo
```

The example supplies a synthetic one-dimensional objective with explicitly
seeded, heteroscedastic uniform observation noise. It requests six initial
points and six two-point batches, hence 18 objective evaluations on successful
completion. It is a usage example, not a measured optimization benchmark or a
physics-validation claim.

## Model and acquisition

The configuration declares a fixed Matérn kernel, prior mean, search interval,
seed, batch size, Monte Carlo sample count, and acquisition-search work limits.
All input coordinates remain in their original units; the kernel's signal
variance and the observation variances are in squared objective units. The
fit subtracts the declared prior mean and uses `Gp::try_fit_diag` without
rescaling, learning, replacing, or flooring the supplied noise variances.
Zero variance is permitted. A singular training covariance refuses rather
than silently inventing measurement noise.

For each joint posterior draw, q-NEI computes:

```text
max(0, minimum latent value at evaluated inputs
       - minimum latent value at proposed inputs)
```

The baseline and candidates are sampled together, so their cross-covariance is
retained. Equal coordinates share the same latent draw, including signed zero;
a candidate already in the baseline therefore has exactly zero improvement.
This objective does not estimate the information value of replicate sampling.
The implementation inherits the existing joint posterior's adaptive numerical
jitter and documented degenerate-covariance fallback; it introduces no stronger
claim about that factorization.

Batch selection is sequential-greedy: already chosen points stay fixed while
the next point maximizes the joint acquisition. Each batch is completely
selected before invoking its objective callbacks. All searches for that batch
use row-wise prefixes of one common normal bank. The leading columns are the
existing scrambled-Sobol normals; widths beyond the embedded Sobol table use
an independently keyed Philox Monte Carlo tail, not fabricated Sobol dimensions.

## Results and limits

`NoisyBoReport` retains every callback input and the original observation and
variance. Its `incumbent_trace` recommends the lowest posterior-mean evaluated
design after initialization and after each batch, with latent posterior
variance. It is deliberately not a monotone raw-minimum trace: later data can
revise an earlier recommendation. Posterior means and variances are model
estimates, not bounds or certificates.

The fixed-kernel driver is synchronous and dense, with the full observation
history as its baseline. Work limits cover initialization, batch count,
sample-bank size, and CMA-ES restart/evaluation counts; no within-factorization
cancellation, statistical stopping, or sparse-baseline approximation is added.
Optional kernel learning is described below. Replay requires the same
configuration and callback outcomes, not just the optimizer seed. No cross-ISA
or optimizer quality claim is inferred from the deterministic BO goldens.

Fourteen focused unit tests are included across the acquisition and fixed-kernel
driver: correlated-Gaussian and noiseless-limit comparisons; a noisy-outlier
case; near-coincident covariance; duplicate-coordinate identity; normal-bank
replay; input refusals; exact callback/history accounting; per-observation
variance retention; posterior recommendations; and sequential/batched execution
paths. These test sources were added without a local Rust toolchain: compilation
and Rust test execution remain unverified. Independent Python Gaussian
calculations check selected equations, not the compiled Rust implementation.

## Learning the kernel without replacing the observation noise

`fs_bo::learning::minimize_noisy_with_learning` is an opt-in consumer of the
same acquisition/history engine. The existing fixed-kernel function and its
configuration remain unchanged. Enable the example's learning lane with:

```bash
cargo run -p fs-bo --example noisy_design -- --learn-kernel
cargo test -p fs-bo --lib hyper
cargo test -p fs-bo --lib learning
```

Supply a `NoisyLearningConfig` containing a positive `refit_every` and a
`HeteroFitConfig`. Declare one lengthscale interval per input dimension and a
signal-VARIANCE interval in squared objective units. Endpoints are positive
physical values, not logarithms; equal endpoints freeze a parameter. The
initial kernel must lie inside these intervals. Declare the number of starts,
maximum local iterations per start, total likelihood evaluations per refit,
projected-gradient tolerance, and seed explicitly. The first start is always
the previous selected kernel; subsequent starts use independently keyed Philox
samples in the log box. Neither noise nor the prior mean is learned.

Learning occurs after initialization and each `refit_every` completed batches,
including the final batch when it lands on that cadence. Other stages still
condition the latest kernel on ALL observations. A refit reuses its winning GP
instead of spending an extra factorization outside the declared allowance.
The example learns after 0, 2, 4 and 6 completed batches, with at most 400
likelihood probes total, and has three ordinary posterior-only conditioning
fits. It still makes exactly 18 objective calls on success.

The reusable `hyper::fit_heteroscedastic` entry point accepts centered values
and the original per-observation variances. It computes analytic Matérn
half/three-halves/five-halves likelihood gradients using the existing training
Cholesky. Bounded projected Armijo steps and seeded restarts retain the best
finite evaluated likelihood. Rejected likelihood probes count against the
same allowance. Exhausted budgets return the retained model and its projected
log-gradient residual, not a convergence claim. A singular initial model
refuses without inventing noise; exact duplicated noiseless constraints are
rejected before a rounded Cholesky pivot can admit them. Noisy replicates are
not deduplicated in the likelihood.

`NoisyLearnedReport` preserves the original objective report plus each refit's
observation count, selected kernel, warm-start and selected likelihoods,
resolved seed, spent likelihood probes, and projected residual. Likelihoods
are compared only within the SAME data prefix. The total likelihood-probe
count is separate from ordinary posterior-only fits and acquisition work.

Nine learner tests cover finite-difference gradients, a hand-derived signal
optimum, heteroscedastic posterior parity, frozen bounds, budget limits,
best-restart retention, invalid/singular input, cancellation and replay. Six
consumer tests cover fixed-kernel trajectory parity at a one-probe allowance,
cadence and warm starts, complete-data conditioning, replay, admission before
callbacks, and initialization-only learning. All 15 new Rust tests remain
unexecuted in the authoring environment: DSR, RCH and Cargo are absent.
Independent Python checks matched the three kernel-gradient fixtures to central
differences within 4.9e-11 and recovered signal variance 1.9999999534 for the
one-point analytic optimum 2. These are equation checks, not Rust test receipts.

The learner's controlled API accepts a continuation hook before fits, between
inverse-column solves and between search steps; a dense Cholesky or one column
solve is not preemptible. The BO driver remains synchronous. Point-estimate
kernel fitting does not add hyperparameter marginalization, replicate
value-of-information, statistical certification, sparse history or an
all-objective performance guarantee.
