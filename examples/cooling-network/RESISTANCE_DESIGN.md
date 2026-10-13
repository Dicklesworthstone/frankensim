# Bounded resistance design against a temperature limit

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
`(constraint_temperature - temperature_limit_k) / temperature_scale_k <= 0`.
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

The model must set `objective.gradient` to `true`. With `mean_wall_region`, the
constraint temperature is that surface's area mean; omit `maximum_smoothing_k`.
That requirement does NOT constrain every solid vertex. A maximum selector is
also supported, under the explicit upper-envelope policy below.
Prescribed-pressure and single fan-bank drives are supported. Existing thermal
contacts and fixed heterogeneous materials keep their original model meaning.
Transient, radiation, recirculation, mesh-convergence and nested sizing requests
are refused rather than solving a simpler model under the same request.

## Hotspot constraints

A complete whole-solid example is included:

```sh
cargo run --release -p fs-cli --features cooling-design -- --json \
  cooling-network resistance-design \
  examples/cooling-network/resistance-hotspot-model.json \
  examples/cooling-network/resistance-hotspot-design.json
```

Select `max_solid_temperature`, `max_wall_region`, or `max_vertices` in the model.
The design file must additionally declare a finite positive `maximum_smoothing_k`
(tau). The constraint is the UNNORMALIZED log-sum-exp over the selected nodes:

`S(T) = m + tau * ln(sum_i exp((T_i - m)/tau))`, where `m = max_i T_i`.

In exact arithmetic, `m <= S <= m + tau*ln(n)`. Thus the smooth constraint is
conservative relative to the selected discrete maximum. No division by n or
subtraction of `tau*ln(n)` is applied. Softmax nodal weights seed the existing
coupled adjoint; tied or changing hottest vertices do not require choosing one
arbitrary active derivative. No extra primal solid solve is needed for the seed.

Smaller tau reduces conservatism but makes the optimization problem sharper.
For the example's 12 vertices and tau=0.04 K, the worst-case mathematical bias is
about 0.0994 K; the actual bias depends on the temperature field. This is a
formula property, not an interval bound on floating-point evaluation or FEM error.
A very large tau can make the declared smooth requirement unattainable even
when the raw peak would pass; that is not proof the original peak problem is
infeasible. The fixed tau never changes silently during a search.

Both S and the actual selected maximum must meet the limit when retaining and
independently replaying a feasible candidate. `objective.kind` is then
`smooth_upper_envelope`: its `value_k` and reported derivatives concern S.
`maximum_selection`, `actual_maximum_k`, `maximum_smoothing_k`, and
`analytic_bias_bound_k` separately disclose the original selector and raw peak.
The active vertex, position and tie count still describe that raw maximum.
The complete temperature field remains in `solid_temperatures_k`.
Ordinary single-model maximum evaluations retain their original unsmoothed
value and active-vertex derivative; only this explicit design policy enables S.

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
cost and work counts. `selected_temperature_k` is the actual mean or selected
peak; `selected_constraint_temperature_k` is the mean or smooth upper envelope.
`initial_temperature_k` reports the initial actual mean or peak.
`constraint_form` distinguishes `area-mean` from `smooth-maximum-upper-envelope`.
`total_solid_solves` counts primal coupled-solid solves; it is not a count of
adjoint solves or Krylov iterations. The selected model is independently
recomputed and must reproduce its feasible value and gradient before publication.

`locally-converged` means the selected point is also the optimizer's converged
terminal point. `feasible-incumbent` means a fully evaluated passing design was
found, but convergence is not claimed for it. `optimizer_terminal` reports the
terminal coordinates and first-order diagnostics separately: those diagnostics
must not be attributed to a different selected incumbent. For a peak request,
they concern the fixed smooth upper-envelope constraint, not the nonsmooth raw
maximum. A spent budget without an evaluated passing sample returns
`cooling-network-design-no-feasible-sample`, not proof of design-space
infeasibility. Cancellation or a producer failure publishes no partial design.

These are nominal fixed-geometry designs, not guaranteed realizable duct shapes,
global optima, uncertainty bounds, continuum-error certificates or validated
hardware requirements. A selected-wall or selected-node limit does not constrain
unselected vertices, and no discrete peak check certifies the continuum maximum.
