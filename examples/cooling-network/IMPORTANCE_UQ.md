# Rare-event probabilities through actual cooling solves

`cooling-network-uq BASE UQ --importance EPSILON,SHIFT...` uses the existing
`fs_uq::ImportanceExecution` to evaluate the strict event
`QoI > temperature_limit_k`. The limit must be present in the UQ document;
equality is compliant. Every observation calls the original cooling product,
with its original geometry, materials, target lowering and child watchdog.

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-importance-hotspot.json --importance 0.2,3
```

The example declares Gaussian inlet-temperature scatter with mean 300 K and
standard deviation 0.25 K. It spends 64 cooling solves. The illustrative shift
of three is in independent standard-normal coordinates, NOT three kelvin.
Choose the proposal, threshold, sample count and model before seeing these
production responses. This example does not establish that the selected limit
is rare or that the proposal is effective for this particular cooling device.

Provide exactly one finite shift per declared parameter in declaration order.
The defensive target mass EPSILON must be in (0,1] with finite reciprocal.
The proposal is `epsilon*phi(z) + (1-epsilon)*phi(z-shift)` in independent germs;
the existing admitted Gaussian factor maps these to physical parameters. Explicit
independence or a joint-Gaussian law is required, including admitted singular
laws. Uniform marginals, unknown dependence and invalid shifts refuse before
physical evaluation. Input-domain failures are not clipped, replaced or redrawn.
Zero shift or epsilon one retains the ordinary target Monte Carlo sampling path.

## What the result means

`frankensim.cooling-network-uq.importance.v1` reports the original parameter
law, sampling key, actual proposal, paid calls and the strict failure threshold.
`estimated_failure_probability` is the mean of the target/proposal-weighted
failure indicators. It is NOT the raw frequency of failures under the proposal,
not normalized by the observed total weight, and not clipped into [0,1].
`failure_standard_error` and `failure_effective_samples` are descriptive
sampling diagnostics, not confidence bounds or valid stopping rules.
`mean_likelihood_weight` diagnoses importance weighting; its population value
under the ideal sampling model is one. `proposal_failure_count` is explicitly
a diagnostic of the proposal, not a target-law probability.

No physical temperature mean, quantiles or distribution are inferred from the
biased proposal responses. No observed failures or a zero empirical standard
error does not establish zero risk. A single shift can miss disconnected modes
or increase variance. Finite-generator grids, numerical likelihood evaluation,
solver error and physical-model validity are not bounded by this report.

Steady objectives and explicit `qoi.kind=transient-sampled-peak` share the same
adapter. One transient observation is the initial/accepted-endpoint maximum of
one COMPLETE declared trajectory, including fixed repeated cycles. It is not a
single time step or a continuous-time peak certificate.

The existing 10,000-call ceiling applies. One soft wall allowance includes
setup, all physical solves and report construction. Bounded synchronous work
is not preemptible, but a deadline reached before publication produces BUDGET,
with no partial probability. Physical refusal remains a terminal failure even
when a deadline also expires. There is no fallback to ordinary Monte Carlo.

This command mode has no durable checkpoint. Recovery, explicit chunks, QMC,
PCE, Sobol/Shapley sensitivity, mean controls, sequential compliance and design
selection refuse before reading input files or reserving outputs. Existing
modes without the new option retain their behavior.

Focused native checks:

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::importance
cargo test -p fs-cli --test cooling_importance
```

These tests exercise actual cooling calls, independent reconstruction through
the library, target-MC identity, complete transient observations, admission,
deadlines and terminal physics failures. They require native execution before
claiming a Rust test pass or a physical-accuracy result.
