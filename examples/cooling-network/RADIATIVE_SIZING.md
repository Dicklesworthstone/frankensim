# Effective-convection sizing with radiative heat transfer

The `cooling-network` effective-h target search accepts a `radiation` policy.
It evaluates the complete solid/air/radiation problem at every candidate, with
fixed geometry and hydraulics. Both supported steady radiation producers work:
area-mean gray exchange with isothermal surroundings, and a closed gray-diffuse
enclosure with caller-supplied, admitted view factors.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/size-radiative-slab.json
```

The example reuses the two-face mixed-air slab geometry, with a 323 K limit on
the first face. The search changes only the last face's effective convective
coefficient, between 10 and 1000 W/(m² K). The two declared surroundings are
300 K and 310 K, with emissivities 0.85 and 0.6. These are synthetic model
inputs, not measured or calibrated hardware.

## Search and result contract

Keep `objective.gradient=true`. The existing `design` fields and both stopping
tolerances are unchanged. The chosen surface must use a declared effective h:
a flow-derived convection law cannot be overridden by this search. A steady
h search cannot be combined with another steady design search or a transient
schedule. Other request admission rules remain in force.

Each proposal runs the full selected radiative producer, including its nonlinear
heat gates, air feedback, material/contact handling and total steady adjoint.
The search does not freeze radiative coefficients from a previous candidate,
add radiation to the air heat, or fall back to a nonradiating solve on failure.
A derivative only proposes an interior trial; its freshly evaluated objective
decides whether the trial passes.

The final `solid_temperatures_k`, `walls`, gradients and `radiation` report all
belong to the selected candidate. `design.total_solid_solves` sums the actual
producer-reported work across candidates, including radiative inner and
reconstruction solves; it is not just the count of outer air iterations.
The radiation report retains its own selected-candidate work and heat accounting.

A successful result is an evaluated local target bracket or a feasible declared
minimum. It is not proof of global monotonicity, global optimality, interior
infeasibility, physical validation or hardware compliance. Exhausted search or
radiative work remains a budget refusal with no partial result published.

## Independent numerical references

Eliminating the air references gives two nonlinear slab-face equations. Their
independent roots for the same 323 K first-face target are:

| Declared radiative model | Target effective h, W/(m² K) |
| --- | ---: |
| No radiation | 194.5014417544 |
| Example's isothermal surroundings | 54.9941404052 |
| Synthetic equal-area two-patch enclosure, unit mutual view factors | 190.7116021084 |

These are independently calculated mathematical reference values, not recorded
executions of the Rust binary. In the enclosure reference, the supplied view
factors are an idealized test input, not a geometric visibility claim about the
solid mesh. The tests also check the independently differentiated total thermal
response rather than treating the search's own adjoint as its oracle.

```bash
cargo test -p fs-cli --bin frankensim network_command::design::tests
cargo test -p fs-cli --test cooling_radiative_design
cargo test -p fs-cli --test cooling_radiation
```

The new end-to-end tests cover both radiative models, the selected nodal field
and total derivative, exact selected-candidate replay, separate heat ownership,
true work counts, a feasible minimum, missing brackets, exhausted budgets,
admission refusals, and deterministic patch-order replay. They require a Rust
workspace environment; they were not executed in the patch-authoring environment.
