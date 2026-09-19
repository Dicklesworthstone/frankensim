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
