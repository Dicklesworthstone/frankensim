# Parallel material-driven cooling

The existing `pressure_cooling` executable now accepts `--parallel`. The same
material table, physical solid, heat source and FEM boundary callback are used
in either topology. Omitting the flag preserves the original series path.

```sh
cargo run -p fs-airflow --example pressure_cooling -- --parallel \
  crates/fs-airflow/examples/pressure_cooling/constant.txt 50 20 0.0002 0.008
cargo run -p fs-airflow --example pressure_cooling -- --parallel \
  crates/fs-airflow/examples/pressure_cooling/temperature_dependent.txt 50 20 0.0002 0.008
```

Both 0.003 m2 wall regions now receive separate streams from the common 300 K
inlet and share the supplied 50 Pa pressure difference. The two positive
mass-flow bounds are **per branch**, not a bound on their combined flow. Each
branch retains the original authored 5e9 /m3 linear hydraulic geometry factor.
The solid mesh, conductivity and total 20 W heat source are unchanged. Negative
heat input still extracts heat.

Every wall iteration recomputes both branch operating points from their current
bulk material states. One common solid solve receives the actual convection
boundary rows. The returned outlet is the adiabatic enthalpy mixture of both
streams on their shared source chart; temperatures are not simply averaged.
Changing viscosity, density, heat capacity or conductivity data therefore
reaches the flow split, heat transfer and solid field through the library APIs.
The supplied tables are synthetic examples, not measured coolants.

Output adds `topology`, `branch_mass_flows_kg_s`, `branch_outlets_k` and
`branch_pressure_drops_pa`. Branch arrays follow the `first`, `last` order.
`mass_flow_kg_s` is their total and `outlet_k` is their mixed outlet. For series
mode, those arrays contain one complete two-region stream. The legacy
`pressure_drop_pa` field reports the accepted first complete inlet-to-outlet
path; **parallel pressure drops are never summed**. Both actual path drops
remain visible in the branch array.

For the constant table, the reduced model predicts 0.005 kg/s in each branch at
50 Pa: 0.01 kg/s total. Absorbing 20 W then raises the mixed outlet to 301 K.
These are analytic reference checks, not claimed output from native execution.
The parallel solve uses at most 64 common-solid callbacks and a single cap of
32768 full branch thermal evaluations across all wall and pressure trials.
Pressure roots are resolved more tightly than the outer thermal gates. Local
source/correlation admission and per-region power checks remain mandatory.

The parallel model assumes one stable, low-compressibility liquid source,
positive branch flows and ideal common manifolds. It does not implement pump
curves, arbitrary graph topology, recirculation, pressure-dependent material
state, boiling, deformation or kinetic-pressure recovery. Hydraulic dissipation
is reported but not deposited into heat; the **whole network** must remain
below the same explicit 0.01 W omission allowance. This is not a complete
mechanical/thermal energy certificate or experimental validation.

```sh
cargo test -p fs-airflow --test liquid_parallel
cargo test -p fs-airflow --example pressure_cooling
```

The integration tests include actual tetrahedral conduction, flow redistribution,
source-backed enthalpy mixing and conservation gates. Native compilation and
these tests were not run in the authoring environment: Cargo/rustc are absent.
