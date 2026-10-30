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
temperature changes. Transient admission additionally checks the full density
and heat-capacity source interval as described below.
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

The fluid's density and specific heat must be constant throughout the declared
temperature interval and exactly match the original reference claims. The
existing material source-line query checks continuous support and **every**
crossed breakpoint, including supported sampled surfaces at the fixed pressure
and composition. Equal endpoint values cannot hide an interior change. Missing
support, claim switching, varying fluid density or Cp, and inlet/backflow
temperatures outside the interval refuse before time evolution.

This admission makes the native face energy flux `rho * Cp * Q * T` consistent
with fluid cell storage `rho * V * Cp * (T_new - T_old)`. The original enthalpy
chart owns storage, using the scene's explicit single-fluid-phase assumption;
solid cells retain their existing potentially variable Cp charts. The same
native backward-Euler energy solver handles conduction, conservative prescribed
advection, fluid and solid storage, local Joule heating, thermal contacts and
admitted exterior radiation. No second advection or storage solver is added.

Cell masses remain fixed while equal incoming and outgoing fluid mass carries
enthalpy across the control volumes. The arbitrary enthalpy reference cancels
because flow is volume-conservative and density is constant. Physical outflow
in each time record includes the native boundary advection and conduction,
without counting the storage diagonal as external heat. The accepted endpoint
rechecks electrical resistivity and thermal conductivity at local temperature.

Changing only source Cp changes both transient heat storage and the transported
heat-capacity flux. Changing a neighboring solid's Cp changes the fluid's
temperature through their shared conduction solve. Voltage/current and authored
power schedules use the existing time request. The result retains
`fluid_reference` alongside final local source identities and per-step energy
accounts. A failed step publishes no partial trajectory and leaves inputs intact.

This is a constant-density incompressible thermal approximation on fixed
geometry. The fluid role assumes one stable phase; it does not discover phase
boundaries. Prescribed flow, frozen reference viscosity, absence of momentum or
pressure evolution, and the other limits above still apply. Varying fluid Cp
requires a conservative enthalpy-advection owner and is deliberately refused.
Backward Euler requires timestep refinement to establish temporal accuracy.

## Focused checks

The new Rust tests cover source-only changes in current and temperature,
source heat-capacity changes in advective cooling with unchanged prescribed
flow, dimensional/semantic conductivity preservation, invalid/reference-domain
and nonconservative-flow refusals, and pack/scene-to-command execution with
replay. Three transient regressions additionally compare the original coupled
solver with an independent three-control-volume balance through drive switch-off,
exercise fluid/solid capacity coupling, and reject interior rho/Cp variation,
unsupported inlet states, late storage escape and cancellation. The file-driven
regression reaches both steady and transient paths. These are authored synthetic
fixtures, not measured liquid-metal or electrolyte calibration; the new Rust
regressions have not been executed in this implementation environment.
