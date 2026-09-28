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

The driver is synchronous and dense, with the full observation history as its
baseline. Work limits cover initialization, batch count, sample-bank size, and
CMA-ES restart/evaluation counts; no within-factorization cancellation,
statistical stopping, sparse-baseline approximation, or automatic noisy-model
hyperparameter training is added. Replay requires the same configuration and
callback outcomes, not just the optimizer seed. No cross-ISA or optimizer
quality claim is inferred from the existing deterministic BO goldens.

Fourteen focused unit tests are included across the acquisition and driver:
correlated-Gaussian and noiseless-limit comparisons; a noisy-outlier case;
near-coincident covariance; duplicate-coordinate identity; normal-bank replay;
input refusals; exact callback/history accounting; per-observation variance
retention; posterior recommendations; and sequential/batched execution paths.
These test sources were added without a local Rust toolchain: compilation and
Rust test execution remain unverified. Independent Python Gaussian calculations
check selected equations, not the compiled Rust implementation.
