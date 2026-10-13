# One cooling design across declared operating cases

`cooling-network resistance-design` accepts an optional `operating_envelope`
inside the existing design policy. Fan speed and selected resistances are ONE
shared design across all cases. The nominal model is always case zero; up to
seven additional named cases may change source magnitude and existing external
supply temperatures. Geometry, material laws, fluid properties, fan curve and
hydraulic topology remain those of the input model.

```json
"operating_envelope": {
  "smoothing_k": 0.02,
  "power_relative_tolerance": 1e-12,
  "cases": [
    {"name": "hot-load", "source_multiplier": 1.3,
     "inlet_temperatures": [{"node": 0, "temperature_k": 301}]},
    {"name": "light-load", "source_multiplier": 0.8}
  ]
}
```

Case names are unique; `nominal` is reserved. Omitted source multipliers mean
1, and omitted inlet entries retain the original model's temperature. Source
multipliers must be finite and nonnegative. They scale uniform/nodal generation
or every component's watts without changing P1 footprints. Component sources
are reprojected with the existing audited power-map producer using the explicit
`power_relative_tolerance` in [0,1). Temperature overrides are positive absolute
kelvin and may name only existing supplies, once each. Every case is constructed
from the ORIGINAL model, never by accumulating earlier case/trial changes.

For case constraint temperatures T_j and their own temperature limits L_j,
the optimizer uses an envelope of excesses m_j = T_j - L_j:

    E = max(m_j) + tau * ln(sum_j exp((m_j - max(m_j))/tau))

with `tau = smoothing_k`. This is an unnormalized upper envelope, NOT an average
or a probability. In exact arithmetic, max(m_j) <= E <= max(m_j)+tau*ln(n). The inequality is E <= 0.
The implementation also checks each case's actual selected mean/maximum against
its own limit; neither another passing case nor a softmax weight near zero can mask
a failure. Existing explicit hotspot smoothing applies inside EACH peak case;
the separate case smoothing never changes the selected vertex set. Both biases
are disclosed and are formula properties, not floating-point/FEM certificates.

Each full optimizer sample evaluates every case and contracts the existing
complete coupled hydraulic/thermal gradients with the envelope's weights.
Design-change and optional hydraulic air-power costs are charged ONCE: sources
and supply temperatures do not change the declared frozen-property hydraulics.
No finite differences drive the search and no extra solid solves create seeds.

`max_evaluations` counts PHYSICAL CASE evaluations, not case sweeps. It must
reserve two complete sweeps (initialization and final replay). A limit not
divisible by the case count leaves the remainder unused; no partial sweep is
accepted. The optimizer's `evaluations` counts case sweeps, explicitly labelled
in `optimizer_evaluation_unit`. All cases share the original wall-time budget;
individual physics budgets remain per case. Final replay rechecks every case's
value and full sample gradient, not just their possibly cancelling combination.

The ordinary top-level cooling result describes the NOMINAL case at the selected
design. `resistance_design.operating_envelope.cases` contains a complete cooling
result and inlet declarations for every replayed case. The design's
`selected_temperature_k` is the largest actual case temperature, and
`selected_constraint_temperature_k` is the nominal limit plus E. With unequal
limits this is a reference-equivalent value, NOT an actual case temperature.
`operating_envelope.constraint_excess_k` reports E directly. Without an operating envelope the
single-case behavior, cost, and budget units remain unchanged.

This admits a finite declared operating set, not all possible loads, an
uncertainty distribution, a safety guarantee, or a continuum maximum. A failed
case/derivative/cancelled sweep publishes no selected design. Exhausting a budget
without a feasible sample is not a proof that the declared design is infeasible.

`operating_envelope.scaled_constraint_gradient` is the envelope gradient in
optimizer coordinate order (selected log-resistances, then optional log-speed).
The ordinary top-level thermal/hydraulic gradients still describe NOMINAL,
not the aggregated constraint.

Run the included policy with the existing whole-solid hotspot model:

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-hotspot-model.json \
  examples/cooling-network/operating-envelope-design.json
```

## Different component duty patterns

On a component-power model a case can instead declare absolute watt overrides:

```json
{"name": "processor-busy",
 "component_powers_w": {"chip": 1.2, "memory": 0.2}}
```

Unlisted components retain their ORIGINAL power. Names must already exist and
watts must be finite/nonnegative, including zero. Existing vertex footprints
are retained; the existing power-map producer recomputes total and delivered
power. Overrides are declarations with unstated uncertainty, not estimates
inheriting the original component's band. `source_multiplier` and
`component_powers_w` cannot coexist in one case. Uniform/nodal-source models
reject component overrides rather than inventing footprints. This allows the
same design to address competing workloads at different locations, not merely
scaled copies of one spatial heating pattern.

## Different component requirements

Each additional case may provide its own `objective` selector and positive
`temperature_limit_k`. Omitting either inherits the original model's selector
or the design's nominal limit. The selectors are `mean_wall_region`,
`max_wall_region`, `max_solid_temperature`, and `max_vertices`; exactly one is
required when overriding. The nested selector does not accept a `gradient`
switch: complete gradients remain mandatory for the shared optimization.

```json
{"name": "chip-duty", "source_multiplier": 1.2,
 "objective": {"max_vertices": [4, 7]},
 "temperature_limit_k": 302.1, "maximum_smoothing_k": 0.03}
```

A peak case requires its own positive `maximum_smoothing_k` or inherits the
nominal design's explicit peak smoothing. A mean case cannot supply a peak
smoothing value and never becomes a maximum merely because nominal is a peak.
The nominal requirement is still always included. To enforce multiple
requirements under the same operating load, declare multiple named cases with
the same source/inlet settings and different selectors or limits. All count
against the same eight-case cap and physical evaluation budget.

Every case's own selector seeds the existing total thermal/hydraulic adjoint.
The envelope combines excesses over the individual limits, not absolute
maximum temperatures: a cooler component with a tighter allowed temperature
can be the controlling requirement. The raw and smoothed selected temperatures
must both pass their OWN case limit before a candidate can be retained.

Every replayed case row includes `temperature_limit_k` and the full result for
that case's objective. `selected_temperature_k` is the largest raw selected
mean/maximum across cases; it is not the maximum normalized violation when
requirements differ. The nominal top-level result and gradients are restored
after both successful and failed sweeps. The per-case objective definitions,
limits and gradients remain fixed through the entire optimization.

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-model.json \
  examples/cooling-network/component-requirements-design.json
```
