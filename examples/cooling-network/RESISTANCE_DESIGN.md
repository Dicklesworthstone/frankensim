# Bounded resistance design against a mean-temperature limit

Build with `cooling-design` and run the two-file command from the repository root:

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-model.json \
  examples/cooling-network/resistance-design.json
```

The model is an ordinary cooling-network v1 request. The separate design file
selects branches, resistance multiplier bounds, a temperature requirement and
optimizer work limits. Neither input file is changed. A build without the feature
returns the explicit `cooling-design-feature` diagnostic for this command.

## What is optimized

For selected branch i, let `x_i = ln(R_i / R_i_original)`. The solver minimizes
`0.5 * sum(weight_i * x_i^2)` subject to the multiplier bounds and
`(mean_wall_temperature - temperature_limit_k) / temperature_scale_k <= 0`.
The cost is a caller-weighted penalty for changing the original design, not a
manufacturing, fan-power or money model. Bounds must be positive, finite, have
nonzero width and contain multiplier 1; the starting point is the original model.
At most 64 unique branches are admitted. Unselected resistances remain fixed.

Every sample recomputes the fan/prescribed-pressure equilibrium, junction mixing,
exchanger capacity, any admitted flow-derived convection, and the shared solid.
The existing total hydraulic/thermal adjoint supplies the constraint gradient;
no finite differences or frozen-wall/flow approximations drive the search. The
optimizer is fs-ascent's existing box-constrained projected augmented Lagrangian,
not a new solver. It cannot silently substitute a zero for an unavailable slope.

The model must select `mean_wall_region` and set `objective.gradient` to `true`.
Maximum-temperature objectives are refused because this engine requires a smooth
inequality. The limit concerns the selected area mean, NOT every solid vertex.
Prescribed-pressure and single fan-bank drives are supported. Existing thermal
contacts and fixed heterogeneous materials keep their original model meaning.
Transient, radiation, recirculation, mesh-convergence and nested sizing requests
are refused rather than solving a simpler model under the same request.

## Budgets and results

`max_evaluations` includes initialization, rejected line-search trials, and one
reserved final physical replay; it must be at least two. `max_iterations`,
`max_multiplier_updates` and `max_backtracks` bound the other optimizer work.
`kkt_tolerance` acts on the dimensionless, temperature-scaled optimization
problem. The fixed numerical policy is recorded in the result's `policy` field.
The model's `wall_seconds` covers the entire search and replay. Its graph,
coupling, linear and derivative budgets still apply to each physical evaluation.

A successful response retains the ordinary complete cooling-network result and
adds `resistance_design`, including selected SI resistances, multipliers, bounds,
cost and work counts. `total_solid_solves` counts primal coupled-solid solves;
it is not a count of adjoint solves or Krylov iterations. The selected model is
independently recomputed and must reproduce its feasible value and gradient
before any result is published.

`locally-converged` means the selected point is also the optimizer's converged
terminal point. `feasible-incumbent` means a fully evaluated passing design was
found, but convergence is not claimed for it. `optimizer_terminal` reports the
terminal coordinates and first-order diagnostics separately: those diagnostics
must not be attributed to a different selected incumbent. A spent budget without
an evaluated passing sample returns `cooling-network-design-no-feasible-sample`,
not a proof that the design space is infeasible. Cancellation or a producer
failure publishes no partial design.

These are nominal fixed-geometry designs, not guaranteed realizable duct shapes,
global optima, uncertainty bounds, continuum-error certificates or validated
hardware requirements. A mean-temperature design does not certify a peak limit.
