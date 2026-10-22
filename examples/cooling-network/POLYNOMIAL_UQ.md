# Polynomial uncertainty through the cooling product

`cooling-network-uq BASE UQ --pce DEGREE,TRAINING` runs the existing physical
cooling command at every training and held-out point. The UQ document's
`samples` is the TOTAL count, including held-out observations. Degree is 0..8;
training must contain at least twice the total-degree basis size and leave at
least two held-out calls. Existing library size/conditioning limits and the
command's 10,000-call ceiling apply before the first physical evaluation.

For the four-parameter `uq-fan-hotspot.json` example, first set `samples` to 32
in a separate request file, then run:

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  my-uq-32.json --pce 1,24
```

This fixed design has five polynomial terms, 24 training solves and eight
held-out solves. It is an illustrative allocation, not a claim that degree one
resolves the particular model. An eight-sample one-parameter study can use
`--pce 1,6`; the existing four-parameter eight-sample request cannot.

Gaussian scatter uses the existing Hermite basis; independent uniform tolerances
use shifted Legendre polynomials. Explicitly joint Gaussian laws retain their
admitted correlation map, including singular laws. Unknown multivariate
dependence and unsupported probability models are refused, never replaced.
The parser, material/geometry inputs, nonlinear solve, boundary lowering and
same-binary child watchdog are the original cooling product's owners.

Steady objectives and explicit `qoi.kind=transient-sampled-peak` use the same
path. A transient observation is the initial/accepted-endpoint maximum over a
COMPLETE declared trajectory, including fixed repeated cycles, not a final
value or an independent timestep. Nonlinear or switching peaks need not be
well represented by a low-degree polynomial; held-out residuals expose some
fit error but do not bound unsampled errors or continuous-time peaks.

The `frankensim.cooling-network-uq.pce.v1` report separates:

* `surrogate_mean_k`, `surrogate_variance_k2`, `surrogate_std_dev_k`: coefficient
  moments of the fitted polynomial, NOT an unbiased physical-model estimate.
* Training/validation counts, R-factor condition estimate, `training_rmse_k`,
  `validation_rmse_k` and `validation_mean_error_k`: numerical/physical checks.
* `held_out_probability_of_compliance` and its descriptive standard error:
  actual held-out cooling responses below the optional temperature limit.
* `germ_basis`, `indices`, `coefficients_k` and original parameter declarations:
  enough to reconstruct the bounded `fs_uq::ProductPceModel` in germ coordinates.

No surrogate threshold indicators or training observations enter compliance.
The evidence remains Estimated; no confidence sequence, physical qualification,
optional stopping or automatic degree selection is claimed. In particular,
coefficient response variance is not sampling variance of an estimated mean.

The same soft `wall_seconds` allowance covers polynomial setup, every cooling
child, fitting and result preparation. Children are killed/reaped at their
watchdog deadline; bounded synchronous QR/reductions are not preemptible. A
model refusal or expired allowance emits a diagnostic and NO partial fit.
This command mode has no durable codec: `--checkpoint`, `--resume` and explicit
chunking refuse before input reads or output reservation. It also refuses QMC,
Sobol sensitivity, adjoint mean control, sequential compliance and candidate
selection rather than ignoring incompatible options. Existing MC/QMC/recovery
behavior without `--pce` is unchanged.

Focused native checks:

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::polynomial
cargo test -p fs-cli --test cooling_network_pce
```

The integration tests invoke the actual binary, independently reconstruct a
fit using real cooling responses, retain complete transient observations, and
exercise incompatible flags, timeout and physical refusal. They must execute
natively before any Rust test-pass or cooling-accuracy claim is made.
