# Mean estimation with one nominal coupled adjoint

```sh
cargo run --release -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-adjoint-mean.json \
  --mean-control adjoint
```

This opt-in path first executes the real steady cooling model at the declared
input means with its coupled adjoint enabled. It fixes the resulting coefficients
before drawing any Monte Carlo observations, then runs every ordinary cooling
sample through the existing same-binary child evaluator.

For actual observations `Y_i` and original sampled inputs `X_i`, the additional
mean estimate is the average of `Y_i - g^T (X_i - E[X])`. The declared Gaussian
means or uniform midpoints determine `E[X]`; the base request and sample means
cannot substitute for those expectations. Joint Gaussian inputs use their
original joint sampler and marginal expectations. No coefficient is fitted to
these observations, and no control is selected afterwards for a favorable result.
A poor but fixed linearization can increase variance; that increase is reported.

Supported uncertain inputs are external inlet temperature, fan-speed ratio,
declared scalar surface heat-transfer coefficients and named contact resistance.
Fan coefficients use the total fan/flow/convection derivative. Log-coordinate
surface, speed and resistance derivatives are converted to the linear parameter
coordinates actually sampled. Missing derivatives refuse; they are never zero.
Other inputs continue to use ordinary UQ without this option.

All original result fields stay unadjusted, including `mean_k`, `std_dev_k`,
quantiles, bounds, `sampling_standard_error_k` and empirical compliance. The new
`mean_control_variate` object reports the separate adjusted mean and standard
error, actual adjusted/raw variance ratio, parameter means, frozen gradient,
nominal objective and adjoint residual. The example requests 32 actual samples
and one additional nominal forward/adjoint evaluation: 33 model calls. They share
one wall-time allowance; the nominal solve is not a random observation. The total
model-call cap is 10000, including the nominal call.

This initial command adapter requires one complete fixed-count Monte Carlo run.
It refuses sequential decision, candidate-selection, QMC, sensitivity and
checkpoint/resume options, as well as transient, radiation, recirculation and
mesh-study bases, before any nominal solve. No partial adjusted distribution is
published. The reusable `fs_uq` assessment itself supports raw-checkpoint recovery
with the original frozen control and cancellation without rerunning physics;
a durable control binding for this command remains separate work.

These are descriptive fixed-sample mean estimates for the declared numerical
model, not confidence sequences or physical safety bounds. Adjusted samples are
not temperatures to use for compliance, quantiles, maxima or CVaR. The control
has no effect on the original Bernoulli compliance event. A nonsmooth maximum's
selected local adjoint can be a poor predictor; no differentiability or guaranteed
variance-reduction claim is added.

This implements the reusable mean-control part of `frankensim-rc-root-q61wp.72`
and its first real cooling consumer. It does **not** close that bead: distribution
sections, native `.fsim` solve/report/ledger integration and its durable nominal
adjoint binding are still required. The canonical product remains `.fsim`; this
command is the existing experimental cooling-network surface, not a new product.

Focused native checks:

```sh
cargo test --release -p fs-uq --lib product_control_variate::
cargo test --release -p fs-cli --bin frankensim uq_command::execute::mean_control::
cargo test --release -p fs-cli --bin frankensim mean_control_requires_explicit
cargo test --release -p fs-cli --test cooling_uq_mean_control
```
