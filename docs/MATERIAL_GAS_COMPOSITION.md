# Conservative runtime gas composition replacement

`ClosedGasParcel::replace_model_conserving_mass_energy` is the uniform-state
composition boundary for material-driven gas simulations.

A caller supplies a complete `ThermochemicalGasModel`: species thermodynamics,
mole composition, and source-backed transport cards. The operation holds parcel
mass, volume, and total internal energy fixed, solves the new model for
temperature/pressure, then resolves viscosity, thermal conductivity, sound speed,
impedance, heat capacities, and the rest of the common gas state from that same
accepted material state.

This means a simulation can replace ordinary air data with a helium-rich,
CO2-rich, or expanded species model as a **data/model change**, without changing
piston, duct, radiation, thermal, or acoustic equations. Consumers that already
read `ResolvedThermochemicalGasState` automatically see the new properties.

The operation is atomic: failed source support, an unbracketed conserved-energy
state, numerical budget exhaustion, or cancellation leaves the old model and
state intact. Extra query axes such as specimen/source conditions are retained.

It is intentionally **not** a mass-transfer, diffusion, valve, combustion, or
reaction model. When composition changes because matter crossed a boundary, the
open-system owner must account for species mass, incoming/outgoing enthalpy,
kinetic energy, and any chemical/source work. This API provides the conservative
material-state replacement needed by such an owner; it does not invent those
fluxes.

The accepted balance is

```text
m_new = m_old
V_new = V_old
m * u_new(T_new, composition_new) = m * u_old
```

and the new EOS supplies pressure at the conserved density. The energy residual
is bounded by the caller's existing `GasEnergySolve` tolerance.
