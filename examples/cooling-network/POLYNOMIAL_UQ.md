# Polynomial uncertainty through the cooling product

`cooling-network-uq BASE UQ --pce DEGREE,TRAINING` runs the existing physical
cooling command at every training and held-out point. The UQ document's
`samples` is the TOTAL count, including held-out observations. Degree is 0..8;
training must contain at least twice the total-degree basis size and leave at
least two held-out calls. Existing library size/conditioning limits and the
command's 10,000-call ceiling apply before the first physical evaluation.

The supplied two-parameter example mixes uniform inlet temperature with Gaussian
component-power scatter, keeping the original fan, geometry and material model:

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-polynomial-hotspot.json --pce 1,20
```

This fixed design has three polynomial terms, 20 training solves and twelve
held-out solves. It is an illustrative allocation, not qualification of a cooling
device. Gaussian power is not silently clipped: an out-of-domain draw refuses.
For four parameters, degree one needs five terms and at least ten training
solves. The existing four-parameter eight-sample request therefore cannot run
a degree-one fit; an eight-sample ONE-parameter study can use `--pce 1,6`.

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

## Independently corrected physical means

Use **`--pce-mean DEGREE,TRAINING,FIT`** instead of `--pce` to estimate the
original numerical-model expectation using a frozen nonlinear mean control.
`FIT` includes its training prefix AND held-out validation. The remaining
`samples - FIT` calls are independent production; at least two must remain.
The original 10,000-call and wall allowances cover BOTH phases, not just
production. The fixed control weight is explicitly one; it is never fitted to
production outcomes or selected from validation results.

Using the same supplied 32-sample request:

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-polynomial-hotspot.json --pce-mean 1,20,24
```

This pays for 20 training, four held-out and eight production cooling evaluations.
The library `PceMeanExecution` owns both phases, freezes the original fit, and
uses the distinct production key `seed XOR 0x5043454d45414e31`. No second
regression, sampler, physical model or control-reduction implementation is
introduced by the command. Avoid reusing either key in separate studies.

The reported physical-mean estimator averages
`Y - (P(germ) - E[P])` on production responses alone. With a fixed fit and
independent draws from the ideal declared law, the correction has zero
expectation. An inaccurate surrogate can INCREASE variance rather than replace
the physical mean by a biased surrogate mean. Finite sampling grids, arithmetic,
model validity and solver error remain outside that ideal-law statement.

`frankensim.cooling-network-uq.pce-mean.v1` reports `mean_k`, `raw_mean_k`, both
sampling standard errors and `observed_variance_ratio`. It reports fitting and
production counts and both seeds explicitly. Its nested `fit` retains the
polynomial and held-out diagnostics, separately from production. The
`raw_production_*` variance, ranges, percentiles and compliance frequency use
actual UNADJUSTED production cooling outputs; adjusted samples are not a
physical temperature distribution. No probability is computed by thresholding
the polynomial or its adjusted values.

Steady and complete-transient observations share the same adapter and watchdog.
Timeout at either phase, fitting, assessment or publication returns BUDGET with
no partial controlled estimate. Physical failures return a terminal diagnostic.
No fallback silently replaces a failed control by the raw mean. Recovery/chunk,
other sampling modes and decision flags remain incompatible, and `--pce` and
`--pce-mean` are mutually exclusive. The library can resume in memory, but this
command currently has no durable combined checkpoint. Completion is neither an
accuracy guarantee nor a demonstrated variance or runtime improvement.

The same focused test file also independently reconstructs the combined study
using actual cooling solves, checks total work, preserves raw production
compliance, exercises insufficient-production/timeout admission and verifies
that a degree-zero control leaves complete-transient MC statistics unchanged.
