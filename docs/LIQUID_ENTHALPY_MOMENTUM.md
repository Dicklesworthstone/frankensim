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
supplied. Unassigned thermal/contact links retain their original owner.
Existing full fluid-face links can opt into the same-sample conduction below.

`tests/liquid_flow_enthalpy/feedback.rs` exercises data-only viscosity curves,
native Darcy flow and donor-cell enthalpy, fixed-viscosity parity, pinned source
support, borrowed law/context consistency and atomic refusal. The independent
Darcy-series energy recurrence is a reference check, not native Rust test
execution or experimental material qualification.

## Molecular conduction from the same material epoch

Call `observer.with_fluid_conduction(&accepted_thermal, &fluid_face_link_ids)?`
to assign existing full fluid-face thermal links. Each actual substep then
uses the **same resolved transport sample** for molecular viscosity and thermal
conductivity. Changing conductivity data changes the next temperature field,
which changes local viscosity and the next momentum solve.

For neighboring cells with normal widths `dx_a, dx_b`, the applied link law is

```
G = A / (dx_a / (2 k_a) + dx_b / (2 k_b))
```

The shared `fs_thermal::enthalpy::series_conductance_w_k` also supplies the
existing material circuit network. The native graded voxel geometry supplies
area and both half-cell distances; no material query, enthalpy law or interface
resistance is added. The original thermal link owner applies opposite heat
rates and uses these same conductances in its combined conductive/advective
timestep bound. Authored links and both covered areas remain unchanged.

This opt-in assignment requires an existing link across adjacent mapped fluid
voxels, their native total surface areas, and full face coverage on **both**
sides of the authoritative contact patch. Partial/solid interfaces, changed
endpoints, duplicate or parallel links on an owned face, and a wrapped observer
supplying that same link's conductance refuse. Unassigned contacts retain their
original owner. Source time and every applied cell enthalpy must still match
the actual thermal state before the stored transport bundle can supply k.
The assignment is bounded by the existing connection budget.

The added native tests exercise a conductivity-only card swap affecting heat
spread, outlet energy, subsequent viscosity-driven flow and the actual CFL
substeps; a graded-cell series-resistance/mass balance; and atomic refusal for
partial coverage, parallel/competing ownership and stale source epochs. These
remain numerical coupling fixtures, not experimentally qualified material laws.
