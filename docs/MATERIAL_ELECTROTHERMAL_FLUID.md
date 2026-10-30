# Conducting fluids in material-driven electrothermal scenes

`frankensim electrothermal scene.json` admits a fluid material alongside
solid materials without changing the electrical or thermal solver. Add
`--transient time.json` to advance supported scenes in time. The existing
`frankensim.electrothermal.v1` scene remains all-solid when `fluid` is omitted.

Add this optional object to the ordinary scene:

```json
"fluid": {
  "material_index": 0,
  "reference_temperature_k": 300,
  "velocity_m_s": [0.1, 0, 0],
  "flow_faces": ["open", "open", "wall", "wall", "wall", "wall"]
}
```

Pair nonzero boundary flow with open thermal faces in the same x-, x+, y-,
y+, z-, z+ order. For the positive x velocity above, use:

```json
"thermal_faces": [
  {"inflow_temperature_k": 300},
  {"backflow_temperature_k": 300},
  "adiabatic", "adiabatic", "adiabatic", "adiabatic"
]
```

The inlet supplies the entering fluid temperature to the original flux
assembly. The outlet advects the adjacent fluid cell's temperature with zero
diffusive outflow; if flow reverses, its positive `backflow_temperature_k`
supplies entering enthalpy. Native open-face rules treat solid portions of
these boundaries as adiabatic. A temperature, heat-flux, convective or adiabatic
thermal face remains closed to mass flow; prescribing velocity through it
refuses. Both new temperature variants require positive kelvin and apply to
steady and transient scenes.

The index refers to the existing `materials` array, not another registry. Every
cell assigned that material becomes a fluid cell. Other assignments remain
solid. The material's existing `conducting` boolean independently chooses bulk
DC conduction or explicit insulation. Missing resistivity for a conducting
material refuses; missing source data never becomes insulation.

The same canonical FSMCDPK card supplies the reference density [kg/m3], specific
heat capacity [J/(kg K)], dynamic viscosity [Pa s] and thermal conductivity
[W/(m K)]. Canonical property names are `density`, `specific_heat_capacity`,
`dynamic_viscosity` and `thermal_conductivity`. The conductivity requirement
retains the source's exact existing dimensional or static ThermalConductivity
schema. The other three use their dimensional scalar requirements. All four
resolve atomically at `reference_temperature_k` and the material definition's
other state coordinates, with the common claim-selection and validity checks.
The shared `FluidTransportState` and `MaterialChtFluid` owners produce the
numerical transport values. No air or representative-liquid preset is used.

## Physics that actually runs

`ElectricSetup::solve_with_cell_resistivity` admits conducting fluid and solid
cells into the **same** current-balance, electrode, sparse-solve and conservative
Joule-deposition implementation. `solve_cell_material_electrothermal` adds that
heat to authored sources and rechecks the electrical material law at the final
thermal temperature. The CLI's existing nonlinear heat iteration also resolves
conductivity at local temperature. It refuses unsupported states rather than
clamping or extending a material table.

The existing energy solver advects heat with the prescribed face volume flux
and the reference `density * specific_heat_capacity`. Changing either source
property therefore changes advective heat transport without changing code or
the prescribed velocity. Electrical and thermal conductivity swaps propagate
through current, Joule heating and the temperature field. The result contains a
`fluid_reference` object with its source-state identity, temperature, density,
heat capacity, kinematic viscosity, velocity and measured volume divergence;
local thermal/electrical state identities remain in `source_states`.

## Explicit model limits

This is **prescribed flow**, not a pressure or momentum solve. The three velocity
components are sampled uniformly on eligible faces. `flow_faces` are in
x-, x+, y-, y+, z-, z+ order. `open` permits the prescribed flux; `wall` blocks
normal flux and is not evidence of a solved no-slip velocity profile. Interior
solid faces block flow. A field whose local volume divergence exceeds floating-
point roundoff refuses, so incompatible obstacles or closed ends cannot silently
create/destroy transported fluid. Use a compatible straight passage, or zero
velocity with six walls for a quiescent conductor.

There is one fluid material because the current heat owner has one reference
fluid heat-capacity flux. In a steady run, density, heat capacity and viscosity
remain frozen at the explicitly requested reference state. That approximation
does not establish the validity of constant properties across arbitrary
temperature changes. Transient admission checks the full density source interval
and uses the original heat-capacity chart for both storage and transport.
Viscosity does not alter a prescribed velocity. Local conductivity and electrical
resistivity do update. There is no variable-density enthalpy advection, melting,
phase-front movement, ion/species transport, electrochemistry, motional EMF,
Lorentz forcing, AC response or experimental qualification.

An electrical electrode box is an ideal equipotential, **not an implicit
mechanical wall**. Assign a solid material to electrode cells when the geometry
should block fluid flow. Solid-solid contact resistance cannot be assigned to
the fluid material.

## Transient fluid and solid storage

Use the ordinary time request described in `MATERIAL_ELECTROTHERMAL_TRANSIENT.md`:

```sh
frankensim --json electrothermal scene.json --transient time.json
```

The fluid's density must be constant throughout the declared temperature interval
and match its original reference claim. Specific heat may vary with temperature.
The existing source-line query checks continuous support and every crossed
breakpoint, including supported sampled surfaces at fixed pressure and
composition. Missing support, source switching, nonpositive Cp, changing density
or inlet/backflow temperatures outside the interval refuse.

The same original bounded enthalpy chart owns both cell storage
`m * [h(T_new) - h(T_old)]` and conservative face transport
`rho * Q * h(T_upstream)`. Each nonlinear pass supplies the existing native
finite-volume operator with per-cell caloric tangents. Both sides of an
interior face use the SAME upstream value and slope, including the enthalpy
offset. They therefore exchange equal and opposite energy even when local Cp
differs. Incoming and backflow enthalpy are evaluated directly on that same
chart at the boundary temperature; outgoing enthalpy follows the fluid cell.

Fluid transients use **first-order upwind enthalpy transport**. The steady
command retains its existing reference-Cp Patankar scheme. A transient result
identifies its scheme as `"fluid_advection":"source-enthalpy-first-order-upwind"`
(`null` for all-solid scenes). The native field API refuses a power-law scheme
for nonlinear enthalpy because its reference-Cp blend is a different numerical
model. No source data or physics implementation changes are needed when an
admitted material's Cp curve changes.

The original backward-Euler operator handles conduction, storage, prescribed
advection, local Joule heat, contacts and exterior radiation together. The outer
material iteration rechecks true chart storage and sums absolute face-local
advective linearization errors, so compensating local errors cannot hide in a
global balance. Each accepted time record independently uses true boundary
enthalpy, rather than the frozen tangent's flux, along with native conduction
and physical heat exchange. Source-chart interpolation error remains separate
from numerical energy closure and material uncertainty.

Cell masses remain fixed while equal incoming and outgoing mass transports
enthalpy through each control volume. A common enthalpy datum cancels because
flow conserves volume locally and density is constant. Storage diagonals are
not external heat. Changing source Cp changes both heat storage and advected
enthalpy; changing a neighboring solid's Cp changes fluid temperature through
shared conduction. Local electrical resistivity and thermal conductivity are
re-evaluated at the accepted temperature, and the existing voltage/current and
authored-power schedules remain active.

The result retains `fluid_reference` as reference-state evidence alongside
source/chart identities and per-step energy accounts. Failed nonlinear,
material-domain, conservation or cancellation checks publish no partial
trajectory and leave accepted input state intact.

This remains a constant-density incompressible thermal approximation on fixed
geometry. The fluid role assumes one stable phase; it does not discover phase
boundaries. Reference viscosity remains frozen and flow remains prescribed.
Pressure/momentum evolution, density-driven volume change and phase transition
require their corresponding owners. Spatial and timestep refinement are needed
to establish upwind/backward-Euler accuracy.

## Focused checks

The new Rust tests cover source-only changes in current and temperature,
source heat-capacity changes in advective cooling with unchanged prescribed
flow, dimensional/semantic conductivity preservation, invalid/reference-domain
and nonconservative-flow refusals, and pack/scene-to-command execution with
replay. Three transient regressions additionally compare the original coupled
solver with an independent three-control-volume balance through drive switch-off,
exercise fluid/solid capacity coupling, and reject interior density variation,
unsupported inlet states, late storage escape and cancellation. Additional native
field regressions check donor-face conservation, reversed flow and enthalpy-datum
invariance. A source-varying-Cp trajectory is compared with independent analytic
three-cell balances through drive switch-off, including exact chart admission
and unchanged-input retry. The file-driven
regression reaches both steady and transient paths. These are authored synthetic
fixtures, not measured liquid-metal or electrolyte calibration; the new Rust
regressions have not been executed in this implementation environment.
