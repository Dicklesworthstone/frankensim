# Joint fan-speed and resistance design

The existing `cooling-network resistance-design` command also admits an optional
`fan_speed_control` in its design policy. It optimizes that speed simultaneously
with the selected branch resistances, not through a nested fan-sizing search.

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-hotspot-model.json \
  examples/cooling-network/joint-fan-resistance-design.json
```

The policy addition is:

```json
"fan_speed_control": {
  "minimum_speed_ratio": 0.6,
  "maximum_speed_ratio": 1.8,
  "weight": 0.2
}
```

These bounds are absolute fan-law speed ratios, not multipliers of the initial
speed. They must contain the model's original speed, have nonzero width, and be
admitted by the retained fan curve. Only the original single fan bank varies;
its count, arrangement and curve samples stay fixed. Existing request limits
and derivative refusals still apply at every candidate.

For speed `s`, the final optimizer coordinate is `y = ln(s / s_original)`.
It follows all selected resistance coordinates in declaration order. By default the cost
adds `0.5 * weight * y^2` to the existing weighted resistance-change penalty.
This is a declared design-change preference, not electrical energy or money.
Every trial re-solves the hydraulic and thermal systems at the SAME new speed
and resistances; the total thermal adjoint supplies all constraint derivatives.
No trial compounds the preceding trial's speed or resistance multipliers.

An empty `controls` array is permitted only with `fan_speed_control`, allowing
a speed-only optimization through the same engine. Omitting `fan_speed_control`
keeps the old fixed-speed resistance problem. Up to 64 resistance controls plus
one speed coordinate are supported, with one shared optimizer evaluation budget
and one reserved full selected-candidate replay.

The ordinary top-level `fan` report belongs to the replayed selected speed.
`resistance_design.fan_speed_control` additionally records the original and
selected speeds, bounds, weight, and zero-based optimizer coordinate index.
It is null when speed is not optimized. Terminal KKT diagnostics still belong
only to `optimizer_terminal.point_log_multipliers`, not a different incumbent.

Mean and explicit smooth-upper-envelope peak requirements retain the semantics
in `RESISTANCE_DESIGN.md`. A successful result is a replayed nominal design,
not a global optimum, certified continuum temperature limit, electrical-power model,
manufacturing specification, or validated hardware recommendation.

## Add an operating air-power preference

Declare `"air_power_cost": { "weight_per_w": 10 }` to add
`weight_per_w * Q * (p_inlet - p_outlet)` to the change penalty. The weight must
be finite and positive; omission preserves the old objective and does not run
a power adjoint. It requires `hydraulics.fan`, but can be used with either a
fixed fan or `fan_speed_control`. All selected resistance and speed weights
remain positive: this is a weighted tradeoff, not an unregularized pure-power
minimizer. The thermal requirement and its gradient are unchanged.

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-hotspot-model.json \
  examples/cooling-network/air-power-design.json
```

The power term uses the same accepted pressure/flow state as the thermal trial.
Its adjoint includes both product-rule terms and movement of the operating
point; it does not infer power changes from speed alone when resistances vary.
This opt-in term adds one hydraulic linearization and power adjoint per complete
model evaluation, not an extra solid solve. They use the declared graph and
linear budgets. As before, the overall wall-time budget covers the entire run.

`resistance_design.air_power_cost` reports the selected air power, weighting and
cost contribution; `electrical_power_w` stays null. The top-level `fan` report
and all selected controls belong to the independently replayed candidate. A
smaller weighted cost need not imply smaller power in every problem, and the
constraint need not be active. Neither endpoint KKT diagnostics nor an observed
power reduction establishes a global optimum or hardware efficiency claim.
