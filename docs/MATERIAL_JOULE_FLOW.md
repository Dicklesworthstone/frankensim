# Material-bound electrical heating and transient flow

`MaterialSlab::step_electrothermal` connects the existing voxel current solver
to the existing source-coupled thermal LBM in both dimensions. Its counterpart,
`step_electrothermal_midpoint`, also closes bulk-resistivity feedback within a
step. No additional electrical solver, thermal solver or material database is
introduced.

The input is `material::electrothermal::ElectricHeatInput`: existing electrodes,
optional existing `ElectricalContact` data, and explicit cell/temperature/energy
limits. A callback receives `(cell_index, temperature_kelvin, retained_transport)`
and returns `Some(resistivity_ohm_m)` or explicit insulation (`None`). A missing
material property must be an error, not insulation. The callback can capture
material-card assignments, pressure, composition and source receipts. It must
represent the same constitutive law on repeated queries. Numerical coefficients
and the returned balances are not evidence of experimental calibration.

## Geometry and energy

The electrical grid is derived from the slab's retained spatial units. In 3D
both owners use the same x-fastest cubic cells. In 2D the electrical grid is one
extruded cell thick, with caller-specified depth; electrode currents and energies
scale with that depth under a fixed voltage while temperature evolution does not.

Thermal wall cells are electrical insulators, and electrodes may not overlap them.
Electrodes constrain electrical potential/current in fluid cells; they do not
create solid obstacles or model the actual mechanics/chemistry of an electrode.
Electrical exterior faces are insulated. This is intentionally independent of
the thermal slab's periodic lateral fluid boundaries. A periodic array of imposed
electrical constraints with insulating electrical seams is not a periodic current
solve. Contacts must join electrically conducting interior cells; their bulk and
surface heating retain the current solver's original partition and are deposited
exactly once. There is no hidden heat source for an electrode/wall reservoir.

Every accepted step reports electrode energy, spatial Joule powers, source energy,
stored sensible-heat change and independently measured wall exchange. Electrical
work/deposition, source transfer and thermal balance are checked separately before
publishing. The combined energy tolerance is the sum of the thermal tolerance and
twice the coupling tolerance. This accounts for electrical-to-sensible-heat energy,
not kinetic energy, viscous dissipation, magnetic energy or pressure-volume work.

## Temperature feedback

The explicit method resolves resistance at the accepted initial temperature.
The midpoint method solves

    T_next = thermal_step(T_old, Joule((T_old + T_next)/2))

with a bounded, optionally under-relaxed fixed-point iteration. Every trial starts
from the same accepted populations, so nonlinear iterations do not advance time
repeatedly. Acceptance checks the largest cellwise midpoint-temperature mismatch,
then re-solves the current at the actual old/new midpoint and checks the L1 sum of
cell source discrepancies times dt. Opposing local errors cannot cancel. The
returned electrical field is the one that actually supplied heat; its evaluation
Kelvin field and both consistency residuals accompany it.

Controls are `material::electrothermal::midpoint::MidpointOptions`. Invalid inputs,
unsupported material states, solver failures, thermal limits, failed balances,
nonconvergence and cancellation leave the slab intact. Captured external callback
receipts remain caller-owned, including trial records. At peak, the original,
trial and thermal-owner tentative grids coexist. Cancellation is polled around
whole thermal steps and within the existing electrical solve, not inside every
thermal-cell update. This path makes no performance or unconditional stability
claim. Midpoint source evaluation alone is not a second-order claim for the split
flow/heat integrator; nonconvergence is not proof of physical runaway.

Density, heat capacity and thermal transport remain the slab's retained snapshot.
Permitted existing transport rebinds may be applied between steps. Fixed contact
values and drives apply throughout each midpoint solve. No phase change,
electrochemistry, ionic transport, magnetic force, skin effect or motional EMF is
implied.

## Data-only example

The two files differ only in heat capacity. Their parameters are explicitly
synthetic. The same current heats for the first half of the trajectory and is
then switched off without resetting the thermal or mechanical state. Temperature
changes also change resistance and hence subsequent electrical power.

```sh
cargo run -p fs-lbm --example material_joule_flow -- \
  crates/fs-lbm/examples/material_joule_flow/low_capacity.txt 100
cargo run -p fs-lbm --example material_joule_flow -- \
  crates/fs-lbm/examples/material_joule_flow/high_capacity.txt 100
cargo test -p fs-lbm --test material_electrothermal_flow
```

The example emits CSV temperature/current/power records and a cumulative sensible
energy balance. Rust tests and example compilation were not executed in the
implementation environment because Cargo/rustc were unavailable; independent
circuit/thermal equations do not replace native execution or measured validation.
