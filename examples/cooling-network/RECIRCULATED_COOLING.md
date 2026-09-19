# Cooling with imposed return air

The `cooling-network` request can bind an adiabatic return from a solved
exhaust boundary to one or more fresh-supply boundaries. The actual intake
mixtures are solved along with the solid/air feedback, not substituted once
using a previous exhaust temperature.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/recirculated-slab.json
```

Add an optional root object:

```json
"recirculation": {
  "model": "prescribed-adiabatic-return",
  "source": "Caller declaration identifying the imposed fractions",
  "temperature_tolerance_k": 1e-9,
  "links": [
    {"supply_node": 0, "return_node": 3, "fraction": 0.6},
    {"supply_node": 1, "return_node": 3, "fraction": 0.2}
  ]
}
```

A fraction is relative to the RECEIVING supply's heat-capacity flow, not to
the donating exhaust. Each supply must retain positive fresh makeup; an
exhaust cannot be drawn down beyond its actual available capacity. Duplicate
pairs, nonfinite quantities, unknown fields, invalid endpoints and fully
closed supply loops refuse. Zero fractions retain the original once-through
numerics, and the output explicitly identifies the zero-return case.

Boundary `temperature_k` values denote **fresh makeup**. The implicit return
solve determines mixed intake temperatures. Declared fractions remain fixed
for temperature, coefficient, contact and fan-speed derivatives. Fresh-inlet
sensitivities include the return feedback; they are not derivatives with
respect to a held-fixed mixed intake. The fan's existing common-flow scaling
identity continues to apply when return fractions remain fixed.

The model is bound through the shared transport constructor, so changes of
fan speed, convection, mesh, duty-cycle endpoint or thermal design do not drop
the returns. Admission is repeated at the actual candidate flow. Existing
steady radiation, matching contact, material and transient restrictions still
apply; these features do not acquire broader physical validity from mixing.

A requested result adds a `recirculation` object with the exact configured
links, actual mixed/fresh temperatures, returned and fresh capacity rates,
returned stream temperatures, the independent fresh/undrawn-exhaust heat
gain, its watt imbalance, and its mixing residual. Existing internal network
heat diagnostics remain separate. A positive requested return without its
numerical report refuses publication rather than silently reverting to
once-through behavior.

**Scope:** this is an imposed, instantaneous, adiabatic return with an
external pressure reset. It does not solve return-duct pressure losses, fan
heat, humidity, transit delay, or fluid storage, and it is not a validated
hardware prediction or an uncertainty certificate. Native `.fsim`/ledger
workflows are unchanged.

## Workload transients

Accepted transient endpoints report evolving mixed supplies and outer
fresh/undrawn-exhaust heat in each history row. The complete cycle reports
`fresh_exhaust_energy_gain_j` and `fresh_exhaust_energy_residual_j`, alongside
its original internal-network energy diagnostics. The outer window balances
solid storage against source input, exhaust and any modeled radiation. A
repeated run's `transient` object still describes only the final cycle in
local time; the new fresh/exhaust integral is not a cumulative all-cycle field.

## Sensitivity to the return fraction

With `objective.gradient=true`, `recirculation_sensitivity.links` reports
`dobjective_dfraction_k` for each configured link, including zero-fraction
links. It is a derivative per unit fraction: multiply by 0.01 for the local
linearized effect of one percentage point. Both signs are possible: a return
can heat a colder supply or cool a hotter one.

These derivatives reuse the total fresh-temperature adjoint. Writing the
mixed intake as `x_i = (1 - sum_j r_ij) F_i + sum_j r_ij T_j`, the corresponding
supply multiplier is the fresh-temperature gradient divided by the fresh
fraction. Multiplying that multiplier by `T_j - F_i` gives the return-fraction
gradient, retaining the same implicit solid/air and modeled radiation
feedback. No perturbed physics solve is used to form this derivative.

All other model inputs and hydraulic flows stay fixed. Perturbations must
respect available exhaust and positive makeup; a zero-fraction derivative is
a right-hand model derivative, not a promise that increasing that link is
admissible. Peak objectives retain the existing active-vertex limitation.
There is no finite-change, monotonicity or physical-validation claim.

## Find a fraction meeting a temperature target

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/design-recirculation.json
```

Inside `recirculation`, add:

```json
"design": {
  "supply_node": 0,
  "return_node": 3,
  "min_fraction": 0,
  "max_fraction": 0.9,
  "temperature_limit_k": 318,
  "fraction_tolerance": 1e-5,
  "temperature_tolerance_k": 1e-5,
  "max_evaluations": 80
}
```

The selected pair must already exist in `links`. Only its fraction changes;
all other returns, materials, heating, fan settings and boundary declarations
are retained. This is a steady scalar search, including optional modeled
radiation, not a nested transient, mesh, effective-h or fan-speed optimization.

The search evaluates both bounds using the complete coupled model. One
passing and one failing endpoint form a local threshold bracket regardless
of which end passes. Available derivatives suggest safeguarded Newton trials;
otherwise bisection is used. Every proposed fraction undergoes full admission
and a real numerical solve. Success returns an actually evaluated passing
field satisfying both bracket-width and temperature tolerances. If both
bounds pass, it returns the declared maximum with `both-bounds-feasible`;
if both fail, it refuses without asserting that every interior point fails.
Cancellation, invalid physics or exhausted evaluation budgets publish no
partial design.

The output adds `recirculation_design` with the selected fraction, failed
endpoint, evaluation history, work counts and a `resolved_request`. Saving
that object as JSON and running `cooling-network` on it independently replays
the selected fixed-fraction model without repeating the search. This is not
a maximum-return, unique-root, global-optimum or hardware-compliance
certificate.
