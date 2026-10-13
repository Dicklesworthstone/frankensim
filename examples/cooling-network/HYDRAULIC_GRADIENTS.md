# Temperature gradients for hydraulic controls

Set `objective.gradient` to `true` in a steady once-through cooling-network
request, then run:

```sh
frankensim --json cooling-network request.json
```

The existing mean-wall and selected-active-vertex temperature objectives now
report an additional `hydraulic_sensitivity` object. The original conditional
thermal and contact derivatives retain their meanings.

When `status` is `available`, the result contains:

- `branch_order` and `dobjective_dlog_resistance_k`: each branch's derivative
  with respect to a relative resistance change, in input branch order. Divide
  by the nominal resistance for a derivative with respect to resistance itself.
- `boundary_nodes` and `dobjective_dboundary_pressure_k_pa`: prescribed-pressure
  derivatives in the original boundary-input order, not sorted node order.
- For fan-driven requests, `dobjective_dlog_fan_speed_ratio_k` and
  `dobjective_dfan_curve_pressure_k_pa`: total derivatives including movement
  of the operating point. Curve entries use original sample order, with fixed
  flow coordinates. Prescribed-pressure derivatives are null in this mode.
- `dobjective_dflow_k_s_m3`: intermediate signed-flow weights, including
  shared-solid feedback and flow-derived convection. These are chain-rule
  weights, not a promise that individual branch flows can vary independently.

The calculation differentiates the accepted equilibrium rather than solver
iteration history. It includes changing exchanger capacities, junction mixing,
external supply/exhaust exchange, and smooth-card Reynolds response. The
hydraulic state must exactly match the state used for the thermal derivative.
The hydraulic residual report is a floating-point linear-system check, not a
bound on temperature-gradient error. Hydraulic linear work uses
`budgets.linear_iterations`; the interface uses `budgets.derivative_iterations`.

Stagnant/near-stagnant flows, zero reservoir exchange, curve knots/domain
boundaries and convection cards without smooth derivatives report this extra
control family as unavailable. Existing conditional gradients remain available.
Recirculation, radiation, transient and nested fan-speed-design paths retain
their existing specialized gradients; this addition does not freeze their
feedback to manufacture an individual-resistance gradient.

All gradients are local nominal estimates at fixed geometry, material laws,
loads and fluid properties. Maximum-temperature objectives retain the reported
active-vertex/tie semantics. This is not a finite design-change guarantee,
hardware validation, uncertainty enclosure or continuum-error certificate.
