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

The command requires a fixed-count Monte Carlo **plan**, but it can now execute
that plan in durable chunks. It still refuses sequential decision,
candidate-selection, QMC and sensitivity combinations, and transient, radiation,
recirculation and mesh-study bases before any nominal solve. Partial prefixes
are progress records, never completed adjusted distributions.

## Durable interruption and recovery

The existing `--checkpoint`, `--resume` and `--max-new-samples` options retain the
actual nominal objective, adjoint residual, frozen coefficients and ordered raw
observations in **one** versioned file. The nominal is published before the first
random sample. In particular, `--max-new-samples 0` pays for and saves the nominal
adjoint without starting any Monte Carlo solve. On recovery, the original adjoint
is read, not recomputed or refitted, and only the missing sample ordinals run.

```sh
# These names must be fresh: an existing checkpoint is never overwritten.
mkdir -p target/uq-adjoint
cargo run --release -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-adjoint-mean.json \
  --mean-control adjoint \
  --checkpoint target/uq-adjoint/first.bin --max-new-samples 8
# Exit 6 means the valid prefix is retained, not that a completed estimate exists.

cargo run --release -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-adjoint-mean.json \
  --mean-control adjoint \
  --resume target/uq-adjoint/first.bin --checkpoint target/uq-adjoint/complete.bin
```

Use the same executable, fixed base bytes, parameter bindings, seed, full sample
budget and deterministic runtime profile. A different per-invocation
`wall_seconds` or chunk length is allowed; it does not change the sampled law or
reset the lifetime sample budget. Within each invocation the nominal (when new),
samples and final mean assessment share one evaluation-time allowance. Complete
checkpoints export the same result without nominal or sample physics.

Ordinary raw-only checkpoints cannot acquire a control after sampling, and a
controlled checkpoint cannot silently become an ordinary run. Changes to nominal
diagnostics, coefficients, model identity or plan are refused before reserving a
new output. Actual sample failures invalidate the current output rather than
filtering inconvenient observations. Storage errors stop before another sample.
The existing atomic write/sync path publishes the entire state together; there
is no separately updated coefficient sidecar.

Progress records include whether the control was restored and completed model
counts both over the retained lifetime and for the current invocation. Complete
results retain lifetime counts: 32 samples plus one nominal remain 33 completed
model calls regardless of chunking. Interrupted attempts, or work lost in a
crash before a checkpoint is durable, are not counted and may need repeating.
This is not an exactly-once or cross-invocation wall-budget guarantee. If a final
assessment exhausts its allowance, a complete saved checkpoint can be resumed
with a new allowance without repeating physics.

The checksum detects corruption and binds all these values; it is not an
authentication mechanism, proof of physical validity or proof of statistical
independence. Resume only checkpoints from a trusted source.

These are descriptive fixed-sample mean estimates for the declared numerical
model, not confidence sequences or physical safety bounds. Adjusted samples are
not temperatures to use for compliance, quantiles, maxima or CVaR. The control
has no effect on the original Bernoulli compliance event. A nonsmooth maximum's
selected local adjoint can be a poor predictor; no differentiability or guaranteed
variance-reduction claim is added.

This implements the reusable mean-control part of `frankensim-rc-root-q61wp.72`
and its first real cooling consumer. It does **not** close that bead: distribution
sections and native `.fsim` solve/report/ledger integration are still required.
The durable nominal binding described here belongs to the experimental command. The canonical product remains `.fsim`; this
command is the existing experimental cooling-network surface, not a new product.

Focused native checks:

```sh
cargo test --release -p fs-uq --lib product_control_variate::
cargo test --release -p fs-cli --bin frankensim uq_command::execute::mean_control::
cargo test --release -p fs-cli --bin frankensim mean_control_requires_explicit
cargo test --release -p fs-cli --test cooling_uq_mean_control
```
