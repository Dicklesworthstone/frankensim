# Conducting fluids in material-driven electrothermal scenes

`frankensim electrothermal scene.json` can now admit a fluid material alongside
solid materials without changing the electrical or thermal solver. The existing
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
fluid heat-capacity flux. Density, heat capacity and viscosity remain frozen at
the explicitly requested reference state. This approximation does not establish
the validity of constant properties across arbitrary temperature changes.
Viscosity does not alter a prescribed velocity. Local conductivity and electrical
resistivity do update. There is no variable-density enthalpy advection, melting,
phase-front movement, ion/species transport, electrochemistry, motional EMF,
Lorentz forcing, AC response or experimental qualification.

An electrical electrode box is an ideal equipotential, **not an implicit
mechanical wall**. Assign a solid material to electrode cells when the geometry
should block fluid flow. Solid-solid contact resistance cannot be assigned to
the fluid material. `--transient` with `fluid` currently refuses: the existing
solid-only storage path must not be mislabeled as a liquid advection model.

## Focused checks

The new Rust tests cover source-only changes in current and temperature,
source heat-capacity changes in advective cooling with unchanged prescribed
flow, dimensional/semantic conductivity preservation, invalid/reference-domain
and nonconservative-flow refusals, and pack/scene-to-command execution with
replay and rejection of unsupported transient storage. They use authored
synthetic data, not measured liquid-metal or electrolyte calibration.
