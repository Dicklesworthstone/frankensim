# Source-backed solid and coolant storage

`recirculating_cooling --solid SOLID.txt` advances the actual spatial solid
enthalpy field together with the finite coolant reservoir. The original command
without this option retains its quasi-steady solid calculation and CSV columns.
No second conduction, caloric, pump, or coupling solver is introduced.

```sh
cargo run -p fs-airflow --example recirculating_cooling -- \
  --solid crates/fs-airflow/examples/recirculating_cooling/solid_low_capacity.txt \
  crates/fs-airflow/examples/pressure_cooling/constant.txt \
  crates/fs-airflow/examples/pressure_cooling/pump.txt \
  0.1 1 4 20 0 299 302 20 80
```

The two existing file arguments and nine numeric arguments keep their original
meanings: tank mass [kg], step [s], count, body heating [W], reservoir heating
[W], supported inlet-temperature bounds [K], and network-pressure bounds [Pa].
The initial solid and reservoir are both at 300 K. All nested source, pump,
flow, and numerical brackets still have to support every trial.

The optional solid table uses the SAME bounded source/license/pressure header
reader as the coolant table, but its rows are `T_K rho_kg_m3 Cp_J_kg_K k_W_m_K`.
There is no solid viscosity column or fabricated liquid claim. Three ordinary
material-card properties supply density, heat capacity and conductivity. The
source pressure is a material-query condition, not the pump pressure rise.
Two to 256 positive, temperature-ordered rows are supported. The source-backed
conductivity table retains the source knots and refuses extrapolation. The
original material caloric owner creates the solid's enthalpy chart with a
0.01 J/kg interpolation allowance, separate from source-data uncertainty.

Replace only `solid_low_capacity.txt` with `solid_high_capacity.txt` to double
the authored solid heat capacity without changing fluid, pump, geometry or
solver code. Both files have the same density and conductivity. All supplied
material and equipment data are synthetic examples, not measurements.

## What the time step actually solves

The library's `StoredCoolingState` binds the accepted nodal solid enthalpy and
reservoir to one clock and the original `EnthalpyBackwardEuler` reference-mass
owner. Every nested trial calls that existing FEM solver from the SAME accepted
old solid field over the SAME backward-Euler duration as the reservoir. It does
not repeatedly advance the solid during nonlinear iterations. The exact solved
Robin fluxes enter the existing liquid circuit; the body power is not injected
directly into the reservoir as well. The source state at 300 K supplies invariant
solid reference density. Evolving chart density never creates or destroys mass.

All original local nodal, source-domain, solid-energy, regional heat, pump,
flow and reservoir-energy gates remain active. Two additional checks bound
aggregate interface mismatch and combined solid-plus-reservoir stored energy
minus external heat. Both histories are published together after a final
cancellation check. On error neither accepted history changes. Original typed
boundary/conduction errors survive the coupling callback boundary.

CSV appends `solid_stored_change_j`, `solid_external_heat_j`,
`interface_balance_j`, and `combined_balance_j` only when `--solid` is supplied.
The original `stored_change_j` remains reservoir storage; it is NOT the total.
`solid_heat_w` is actual heat crossing the coolant boundary and need not equal
body generation while the solid stores or releases energy. All reported values
come from the returned solvers, not analytical answer substitution.

The example retains the existing one-million branch-evaluation and 8192 solid
callback caps per physical step. Each spatial response uses up to 32 Newton
attempts, restart 24 and 16 linear cycles, with a 1e-8 J solid energy allowance.
Reservoir closure stays at 1e-4 J; the additional interface/combined gates are
5e-4 J each. None of these permits material extrapolation. The original 0.01 W
and 0.01*dt J hydraulic-heating omission limits are unchanged.

## Physical and verification boundary

The solid is stationary on fixed geometry with invariant reference mass. The
example uses a single solid-phase caloric chart and isotropic k(T); the core
adapter also preserves original tensor/element conductivity assignments and
matching thermal interfaces. It does not add material motion, expansion,
phase-dependent mechanics, fluid residence time, boiler physics or hydraulic
heat deposition. The reservoir is fixed-pressure/free-volume, not a sealed
rigid tank. Closing the thermal accounts is not a full mechanical-energy proof.

```sh
cargo test -p fs-airflow --test stored_cooling
cargo test -p fs-airflow --example recirculating_cooling
```

Five library regressions and the example's material/argument regressions include
actual FEM calls. They were authored but not run in the authoring environment:
Cargo/rustc are unavailable. Independent NumPy equations check pulse storage,
continued tank warming after power removal, solid-data substitution, combined
energy and temporal refinement. They do not execute the Rust/IQN/FEM code.
