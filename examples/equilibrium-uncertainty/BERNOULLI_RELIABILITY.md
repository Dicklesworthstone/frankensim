# High-reliability decisions with Bernoulli mixture bounds

`equilibrium_uq` now offers an explicit binary-outcome confidence sequence:

```text
--require-probability 0.99 --confidence-alpha 0.05 \
--confidence-method bernoulli-mixture
```

The method is available for both a selected displacement event and
`--all-constraints`. In joint mode, a success still means **every authored
requirement passed on the same physical parameter draw**. It does not multiply
marginal probabilities. The inference choice changes neither the sampler nor
physical equality bands, force/penetration/energy limits, or solver admission.
Failed draws terminate the run; they are not removed or replaced.

## Run the existing contact model

This bounded force range lies strictly within all five example requirements,
although each draw still runs the actual two-case equilibrium model:

```sh
cargo run -p fs-cli --bin equilibrium_uq -- \
  examples/equilibrium-uncertainty/joint-reliability.model \
  examples/equilibrium-uncertainty/joint-reliability.fit \
  --method mc --samples 2048 --seed 73 --all-constraints --independent \
  --uniform-x force-N -0.1 0.1 \
  --equality-tolerance settled-band 0.00390625 \
  --require-probability 0.99 --confidence-alpha 0.05 \
  --confidence-method bernoulli-mixture
```

Use `--case load-a --target 0 --limit-m 0.02` instead of the all-constraints
and equality-tolerance options to assess just the first displacement. This
retains the original displacement output schema and its explicit count of
unassessed requirements.

Omitting `--confidence-method` retains the existing Gaussian-mixture method;
`--confidence-method gaussian-mixture` selects it explicitly with the same
output. Either named method requires the paired probability/alpha policy.
RQMC refuses these confidence options before reading model files because its
within-net points are not independent Bernoulli observations. Its existing
between-scramble errors are unchanged.

## Meaning of the new bounds

`fs-eproc::bernoulli::BernoulliMixtureCs` mixes the *ordered binary sequence
likelihood* against a fixed Beta(1/2,1/2) tuning distribution. With s successes
and f failures it uses the likelihood-ratio process

```text
E(p) = B(s+1/2, f+1/2) / (B(1/2, 1/2) p^s (1-p)^f).
```

Under a fixed conditional Bernoulli success probability p, this is a test
martingale (a supermartingale at the null endpoints). Inverting the threshold
`E(p) < 1/alpha` gives a time-uniform confidence sequence via Ville's inequality.
This is a likelihood-mixture confidence set, **not a Bayesian posterior
credible interval**; the mixing distribution does not assert physical prior
knowledge. The construction follows the conjugate-mixture method in Howard,
Ramdas, McAuliffe and Sekhon, *Time-uniform, nonparametric, nonasymptotic
confidence sequences* (2021), section 3.2, arXiv:1810.08240.

The numerical implementation uses Boolean observations, compensated log
predictive updates, and bounded bisection retaining outer numerical brackets.
It has an explicit one-million-observation envelope, matching this command's
sample cap. These floating-point computations are **not an outward-rounded
interval certificate**. Zero failures give nonzero uncertainty; small budgets
can and do return `inconclusive`.

The UQ adapter `UqExecution::assess_bernoulli_compliance` reconstructs the same
process from every retained threshold indicator. It changes no work counters,
observations, report fields, or checkpoint bytes. Its interval is asymmetric;
`converged` requires both distances from the empirical mean to be within the
caller's target half-width. The command instead compares lower/upper bounds
to the fixed probability requirement at the original completed prefixes
2, 4, 8, ... and the final sample cap.

## Efficiency and limits

An independent no-failure numerical reference at alpha=0.05 resolves a 99%
requirement after 1,024 samples under this method, versus 65,536 under the
existing generic Gaussian bound, using the same doubling schedule. This
illustrates rare-event efficiency, **not a universal speedup or a measured
native mechanical run**. Other data can favor the original method.

Choose the inference method, event, alpha and model **before inspecting the
samples**. Computing both bounds and selecting the tighter after seeing data,
changing equality bands, or searching across models/seeds requires separate
multiplicity control. Joint confidence concerns one complete event; per-row
statistics remain descriptive. No physical-model, discretization, finite-grid,
or material-calibration guarantee follows from these sampling bounds.

Output identifies `beta-half-bernoulli-mixture-confidence-sequence` and its
sampling-only scope. A budget stop stays inconclusive unless the interval
actually resolves the requirement. Empirical mean/dispersion errors at a
data-dependent stopping time remain descriptive, not new confidence bounds.

Focused native checks:

```sh
rch exec -- cargo test -p fs-eproc --lib bernoulli::
rch exec -- cargo test -p fs-eproc --test bernoulli_mixture
rch exec -- cargo test -p fs-uq --test bernoulli_compliance
rch exec -- cargo test -p fs-cli --bin equilibrium_uq --test equilibrium_bernoulli
```

The tests are authored but not executed in the implementation environment,
which has no rch, cargo or rustc. Independent SciPy beta-function/root checks
and finite-horizon crossing calculations do not execute Rust or fs-math.
