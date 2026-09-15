# Empirical UQ over actual coupled Cooling solves

`cooling-network-uq` runs fixed-count Monte Carlo by executing the existing
steady `cooling-network` product once for every sampled parameter vector. It does
not substitute `fs-uq` unit-test callbacks for the cooling model.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json
```

The first argument is an ordinary **steady** cooling-network request. It must
have `objective.gradient=false` and may not contain transient or design controls.
The second document uses schema `frankensim.cooling-network-uq.v1` and declares a
fixed sample count, deterministic seed, whole-run wall budget, optional
temperature limit, explicit dependence, and uncertain parameters.

Supported parameter targets are:

* `air-density`
* `air-specific-heat`
* `inlet-temperature` with its external-inlet index
* `fan-speed-ratio`
* `surface-htc` for a surface that declares scalar `htc_w_m2_k`
* `component-power` for a named `solid.component_power` row

A sampled scalar `h` cannot override a correlation-derived convection surface.
Changing component power also recomputes the declared component total before the
sample is admitted, so the existing power-map audit remains active. Fan speed,
inlet temperature and all other sampled values pass through the ordinary cooling
parser and producer gates on every sample.

Distributions are `gaussian` (`mean`, `std_dev`) or `uniform` (`lo`, `hi`). A
uniform support must remain inside the target's basic physical domain. Gaussian
sampling is not silently truncated: if a sampled value leaves the physical or
model domain, that model evaluation refuses the **entire** UQ execution rather
than dropping the sample and biasing the distribution.

Dependence is explicit:

```json
{"kind":"independent"}
{"kind":"unknown"}
{"kind":"joint-gaussian","matrix":[[1,0.4],[0.4,1]]}
```

`unknown` with multiple random dimensions refuses before a cooling sample is
run. `joint-gaussian` uses the existing `fs-uq` numerical PSD admission and
requires Gaussian marginals. Correlations alone are not interpreted as a copula.

The result reports empirical mean, standard deviation, p05/p50/p95, observed
min/max, descriptive standard error, and optional empirical probability
`temperature <= temperature_limit_k`. All are labeled **Estimated empirical
Monte Carlo**. No confidence sequence, optional-stopping guarantee, model-form
uncertainty bound, mesh-convergence certificate, experimental validation, or
native `.fsim`/ledger package is inferred.

Each sample is handed to a child invocation of the same `frankensim
--json cooling-network /dev/stdin` binary. This deliberately spends process
startup time to ensure the parser and physics product boundary are identical to a
normal cooling run. The present handoff therefore targets the repository's Unix
platforms (Linux/macOS). The parent drains bounded child output concurrently and
kills the active child if the declared whole-UQ wall budget expires. No partial
distribution is published after budget exhaustion or a failed model sample.

For deterministic replay, the Philox sample ordinal is keyed to the UQ seed and
fixed plan. Running the same base model and UQ document produces the same sample
vectors; deterministic cooling mode then provides the same empirical reduction
on the same admitted software/hardware profile.
