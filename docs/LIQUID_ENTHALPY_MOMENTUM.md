# Live liquid material feedback on the enthalpy timeline

The source-bound `LiquidFlowEnthalpyInput` now has a `feedback` entry for a
temperature-dependent liquid viscosity. It uses the same pinned
`LiquidTemperatureLaw`, sourced heat-capacity chart, pressure/composition
context and actual borrowed caloric law as the existing frozen-flow binding.

```rust,ignore
// binding: LiquidFlowEnthalpyInput, already naming the original liquid/caloric
// sources, actual voxel-to-thermal cells and incoming reservoir enthalpies.
let observer = binding.feedback(
    &accepted_thermal, &voxel_domain, &simple_config,
    existing_thermal_observer, &cancel_gate, max_momentum_solves,
)?;
let (next_thermal, next_observer, energy) = accepted_thermal.advance_observed(
    &observer, target_time_s, thermal_config, checkpoint,
)?;
```

Every actual thermal substep follows one physical sequence: obtain T from that
cell's current h and original caloric law; resolve the pinned liquid material;
assemble its local molecular viscosity in the [existing full-stress momentum
operator](MATERIAL_FLOW_FEEDBACK.md); pass the actual conservative face fluxes
to the original `Advected` owner. That owner supplies the same-step mass and h
fluxes, boundary-energy account and combined conductive/advective CFL. The
wrapped source/contact/history observer remains on those same substeps. This
adds no alternate heat integrator or fabricated Cp-based temperature update.

`applied()` retains the last solved flow and its exact starting epoch and source
samples. `endpoint()` retains material evaluations at the accepted thermal
endpoint. `audit()` sums the original advective accounting over the advance.
Boundary enthalpy is already included in the thermal audit and must not be
added a second time. `inner()` preserves the original physical history.

This source-bound entry admits fully liquid, validity-wide constant density.
Viscosity can vary locally and between substeps, and source-domain failures
refuse instead of switching claims or extrapolating. The existing Cp/rho claim
pins, complete non-temperature context and actual caloric-object identity are
still checked. Both the initial and proposed endpoint must be admitted. Any
material, momentum, budget, history or cancellation failure discards the whole
tentative thermal/flow/history candidate.

Momentum is quasisteady: its equilibration time must be small relative to the
physical thermal evolution time. The returned flow belongs to the explicitly
reported starting temperature, not an implicit endpoint solve. Time and space
refinement remain necessary. Density-change continuity, phase mobility, moving
surfaces, inertia history, pressure-work heat and viscous dissipation are not
supplied. Thermal conduction remains with the original thermal/contact owner;
local k is retained in the same material samples for that owner's use.

`tests/liquid_flow_enthalpy/feedback.rs` exercises data-only viscosity curves,
native Darcy flow and donor-cell enthalpy, fixed-viscosity parity, pinned source
support, borrowed law/context consistency and atomic refusal. The independent
Darcy-series energy recurrence is a reference check, not native Rust test
execution or experimental material qualification.
