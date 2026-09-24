# Which physical uncertainties drive cooling temperature?

Run variance-based sensitivity against the actual cooling solver:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-sensitivity-hotspot.json \
  --sensitivity sobol
```

The example declares four independent uncertain inputs: air density, inlet
temperature, fan speed, and chip power. It requests 128 pick-freeze rows,
768 actual cooling evaluations in total. The 300-second wall allowance is a
user-adjustable cap, not a promise that this design completes on every machine.
The input intervals are illustrative assumptions, not measured tolerances.

The `frankensim.cooling-network-uq.sensitivity.v1` result names every physical
input and reports `first_order` and `total_order`. First-order effects estimate
the part of temperature variance associated with that input alone. Total-order
effects include interactions involving that input; a substantial gap can point
to combinations worth investigating. They are estimates under the supplied
probability model, not local derivatives or causal attributions. No surrogate
fitting or alternative thermal model is inserted into the physical evaluation
path.

## Design and cost

For d uncertain inputs and N rows, set `samples = N * (d + 2)`, with N at least
2 and the total at most 10,000. Each row evaluates two independent base vectors
A and B, followed by d hybrids: A with coordinate i replaced by B_i. The
existing Philox product sampler addresses the bases at ordinals 2*r and 2*r+1.
There is no silent rounding, dropped tail, budget increase or QMC fallback.

Only explicit `correlation.kind=independent` Gaussian and uniform marginals
are admitted, each with nonzero uncertainty. Keep deterministic parameters in
the base request instead. Even a joint-Gaussian identity matrix is refused by
this narrow adapter: declare independence explicitly. Ordinary coordinate
swaps would not preserve a dependent physical-input law. The existing MC/QMC
modes continue to support their own broader joint-law admission.

The estimators are the Jansen pick-freeze ratios:

```
V    = unbiased sample variance of the 2*N base outputs (A and B only)
S_i  = 1 - mean((Y_B - Y_ABi)^2) / (2*V)
ST_i =     mean((Y_A - Y_ABi)^2) / (2*V)
```

Finite-sample estimates may be negative or above one and are not clipped.
Total effects need not sum to one because interactions overlap. Tiny designs
are suitable for checking wiring, not for stable rankings. `base_output_std_dev_k`
is a descriptive base-sample dispersion, not an error bar on the indices.
Sampling error is not quantified here. Zero sampled variance or an
unrepresentable normalization produces an explicit refusal, never all-zero
indices implying that nothing matters.

## Physical observables and interruptions

All existing uncertainty-target validation and physical domain restrictions
remain active. A rejected material state, fan input, contact solve or child
failure aborts rather than filtering or replacing the observation. Nonlinear
responses and interactions are evaluated by the existing producer.

Steady requests analyze the declared temperature objective. Transient requests
must declare `qoi.kind=transient-sampled-peak`; each base and hybrid call then
executes one complete trajectory, including the configured fixed repeat count.
The response is its initial/accepted-endpoint peak, not the final temperature
or an individual timestep. Use the existing transient targets, including
initial temperature and interval loads, as described in `TRANSIENT_UQ.md`.
A `temperature_limit_k` does not turn sensitivity into a compliance-indicator
analysis and does not produce an approval or probability decision.

A wall interruption returns `BUDGET` with no partial sensitivity distribution.
This CLI mode currently refuses checkpoint/resume/chunk flags, QMC flags,
sequential compliance and candidate-selection flags. `SobolExecution` supports
sample-boundary cancellation, partial-row retention and in-memory clone/resume,
but has no durable checkpoint transport. Its estimates are unavailable until
the complete predeclared design succeeds.

No confidence interval, adaptive stopping guarantee, continuous-time peak
bound, mesh convergence, physical-model validation or safety claim is made.
The existing PCE `sobol_indices()` remains a separate algebraic decomposition
of a supplied independent-germ surrogate; it is not used here.

## Focused native checks

```bash
cargo test -p fs-uq --test product_sensitivity
cargo test -p fs-cli --bin frankensim uq_command::execute::sensitivity
cargo test -p fs-cli --test cooling_sensitivity
```

The tests include analytical additive/interaction models, marginal sampling,
exact input pairing, cancellation/retry, numerical range, and real steady and
transient cooling commands. Their presence is not a claim that native tests
have run on this checkout.

Reference: M. J. W. Jansen, *Analysis of variance designs for model output*,
Computer Physics Communications 117 (1999), 35-43;
doi:10.1016/S0010-4655(98)00154-4.
