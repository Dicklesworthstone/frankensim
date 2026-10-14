# Recirculating material-driven coolant

`recirculating_cooling` closes the previous pump/coolant/solid calculation around
finite coolant storage. Returned fluid is no longer discarded into an unlimited
cold supply. The same source-backed enthalpy chart owns both transported heat
and the reservoir state, and each inlet-temperature trial reruns the existing
pump, branch, convection and common-solid calculations.

From the workspace root, with the pinned dependencies materialized:

```sh
cargo run -p fs-airflow --example recirculating_cooling -- \
  crates/fs-airflow/examples/pressure_cooling/temperature_dependent.txt \
  crates/fs-airflow/examples/pressure_cooling/pump.txt \
  0.1 5 4 20 -5 290 315 20 80
```

The arguments after the two files are reservoir mass [kg], time step [s], step
count, solid heat input [W], separate reservoir heat input [W], the supported
supply-temperature bracket [K], and the network pressure bracket [Pa]. Negative
heat inputs remove energy. There are at most 128 steps. This example starts at
300 K and uses two parallel coolant paths; each has the existing explicit
0.0002--0.012 kg/s mass-flow bracket. It does not widen any bracket on failure.

Only replace the coolant file with `constant.txt` to compare material data.
Both files and the supplied pump curve are synthetic examples, not measured
materials or validated equipment. Both source readers are imported from the
existing `pressure_cooling` example; no new file format, material database,
head law, transport equation or solid solver is introduced. The original
once-through command is unchanged.

## Conservation and retained state

For constant reservoir mass M at fixed source pressure, backward Euler solves

    M * (h_new - h_old) = dt * sum(mdot * (h_return - h_new)) + external_heat_J

Every trial uses h_new to set every supply inlet and reevaluates all material
properties and the pump operating point. The right-hand side uses NET circuit
heat uptake. Counting raw returning mdot*h would spuriously create energy just
by circulating fluid. No fixed Cp or temperature averaging replaces the source
chart. The nonlinear solve requires an absolute joule residual; a small change
in temperature is not acceptance. Source validity remains independent of the
chosen numerical tolerances.

`LiquidReservoir` retains its original source chart, mass, enthalpy state and
accepted time. `advance_reservoir_with_pump` uses prescribed wall temperatures;
`advance_reservoir_with_pump_coupled_iqn` uses the existing common-solid callback.
Both are inherent methods of `ParallelLiquidNetwork` and use one time-step
implementation. A one-branch network represents a serial circuit.

Only a fully accepted step changes reservoir state/time. The returned value
retains the actual pump and network result; for common-solid coupling it also
retains the original solid field/report and the exact applied and recomputed
thermal boundary states. Trial callbacks must not modify accepted solid history
or perform irreversible external writes. A transient solid callback must itself
use the same backward-Euler interval; no solid history scheme is inferred.

The executable uses actual tetrahedral FEM with the same 0.06 x 0.10 x 0.01 m
block, declared 5 W/(m K) solid conductivity and two 0.003 m2 cooled patches as
`pressure_cooling`. The other faces are adiabatic. Only the reservoir is dynamic
in this example: the solid is solved quasi-steadily at each accepted endpoint.
CSV output contains reservoir/return temperatures, total mass flow, both solid
wall temperatures, pressure, stored/returned/external joules, the balance defect,
omitted hydraulic joules, independent solid heat and actual solid-call counts.
It contains computed results, never substituted reference answers.

## Budgets and physical boundary

Each physical step permits at most 64 inlet trials, one million full branch
thermal evaluations and 8192 actual common-solid solves. These are SHARED across
all inlet and pump-pressure trials, not fresh allowances per candidate. Branch
and pump controls remain active. The guarded scalar secant includes periodic
midpoint progress without claiming a unique or stable nonlinear solution.

The stable liquid is perfectly mixed in a constant-mass reservoir at the
original thermodynamic material pressure. Volume is free to change; this is
not a sealed rigid tank. Pipes and manifolds have no residence delay or storage.
There is no evaporation, mixing of unrelated liquids, phase evolution, pressure
EOS update, pump-speed dynamics, cavitation or automatic pump-heat deposition.

Hydraulic dissipation is still reported, NOT deposited into thermal energy.
The example enforces the original aggregate 0.01 W omission limit and an
additional per-step 0.01*dt J limit. These bound all nested trial calculations;
neither the joule residual nor pressure closure certifies full mechanical plus
thermal energy. A real cooler may be represented through explicit cold-wall
exchanges or a separately supplied, justified reservoir heat input.

## Verification status

```sh
cargo test -p fs-airflow --test liquid_recirculation
cargo test -p fs-airflow --example recirculating_cooling
```

Seven new library tests and two new example tests were authored, not executed:
the authoring environment has no Cargo/rustc. Native compilation, Rust tests and
FEM example execution remain unverified. Existing pump-reader tests are reused.
Independent Python checks exercised 60 nonlinear pump/return/reservoir cases,
three constant-property balances and first-order time refinement. They do not
execute the repository's Rust, IQN implementation or FEM solver.

With the constant table, a 0.1 kg reservoir, 20 W solid input and -5 W reservoir
input imply 0.075 K/s: after 5, 10, 15 and 20 s the tank should be 300.375,
300.750, 301.125 and 301.500 K. This follows from energy balance and is a reference
check, not claimed program output. The variable-Cp table changes the enthalpy-to-
temperature relation as well as the flow and solid field.

The well-mixed-storage distinction between enthalpy flow, external heat and
mechanical work follows the Modelica Fluid lumped-volume convention. The fixed-
pressure enthalpy form here explicitly allows changing volume; it does not copy
that library's fixed-volume internal-energy model:
https://build.openmodelica.org/Documentation/Modelica.Fluid.Interfaces.PartialLumpedVolume.html
