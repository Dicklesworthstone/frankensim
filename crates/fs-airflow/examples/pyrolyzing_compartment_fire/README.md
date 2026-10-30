# Heated solid feeding a reacting, ventilated compartment

This example couples the existing spatial solid-decomposition model to the
existing reacting gas and pressure-driven ventilation models. Heating determines
the fuel release; neither a fuel-feed history nor combustion power is supplied.
It exercises oxygen starvation, unburned fuel accumulation, natural ventilation
and pressure-dependent mechanical exhaust with combined conservation accounts.
The default solid now shrinks as material volatilizes, with implicit conduction
and reaction followed by conservative remapping of retained enthalpy.

The material is an **explicitly manufactured methane-releasing solid**. It is
not PMMA, wood, or a calibrated experimental fire. The source-backed PMMA/TGA
examples remain separate: their effective outgoing enthalpy cannot simply be
substituted for the original gas's physical formation-enthalpy reference.

## Run

```sh
cargo run --release -p fs-airflow --features compartment-fire \
  --example pyrolyzing_compartment_fire -- all 60 0.02 24 50000 receding
```

Arguments are case (`all`, `closed`, `restricted`, `open`, `fan-exhaust`), end
time in seconds, maximum timestep in seconds, solid cell count, external
irradiance in W/m², and geometry (`receding` or `fixed`). The values above are
the defaults. `fixed` selects the earlier fixed-geometry explicit solid solver
for comparison. CSV goes to stdout;
source declarations and final balance summaries go to stderr. The face
temperature in each row is the lagged face temperature used over the accepted
interval, while the cell and gas temperatures are the accepted endpoint values.

The example bounds work and admits at most 96 cells, 120 s, 100 kW/m² external
irradiance, and 200,000 accepted steps. The minimum network step is 1e-5 s.
Those are numerical input ceilings, not a guarantee that every such input stays
within the model's physical domain.

## Existing owners and energy reference

The solid uses the original `fs-thermal::enthalpy::pyrolysis::PyrolysisSimulation`
on the original deposited-cell thermal graph, wrapped by the existing
`pyrolysis::slab::RecedingSlab` in the default mode. Its ordinary links conduct
heat through the slab; its Arrhenius reaction changes retained mass and exports
integrated volatile mass and enthalpy. The coupled implicit solver reuses that
same local reaction/energy inversion. After each accepted solid interval, actual
retained mass and the supplied constant density determine the slab thickness.
Limited enthalpy reconstruction and overlap integrals rezone the retained mass
onto equal-mass cells, and the original caloric inverse recovers temperatures.
Geometry and series conductances are rebuilt with the back face fixed. This
transfers extensive enthalpy, rather than averaging temperatures or inventing
residual mass. The existing `VolatileMixture` uses the original
NASA gas calorics to supply that outgoing enthalpy and map release onto the
canonical gas species axis.

The manufactured solid's reference is explicit:

\[
h_s(300\,\mathrm{K}) = h_{\mathrm{CH_4}}(300\,\mathrm{K})
                        - 10^6\,\mathrm{J/kg}.
\]

Its constant solid heat capacity and the original temperature-dependent gas
heat capacity then determine the decomposition enthalpy at other temperatures.
No additional empirical heat of pyrolysis is subtracted from that same account.
The original NASA9 CH4, CO2, H2O, N2 and O2 data, continuity convention and source
revision are reused unchanged from
[`ventilated_compartment_fire`](../ventilated_compartment_fire/README.md).

For each attempted network interval:

1. Solve a bounded face energy balance using the accepted gas and front-cell
   temperatures and current slab thickness. The front half-cell conductance
   balances convection, absorbed external irradiation, and radiation to a 300 K
   blackbody environment.
2. Advance an unpublished solid candidate with that boundary power. Transfer
   every cell's actual integrated volatile release to the gas using the original
   formation-energy convention. Apply the opposite gas convection power once.
3. Transport gas through the actual pressure-driven openings. Inside every
   pressure trial, apply the original inventory-limited reaction
   `CH4 + 2 O2 -> CO2 + 2 H2O` to the transported inventory and recompute pressure.
4. Install the accepted gas and matching solid candidate together. Rejected
   intervals recompute release at the shorter attempted duration from the same
   accepted solid; they never accumulate discarded sources.

Chemical conversion preserves physical chemical-plus-thermal energy through
the original species references. The reported `reference300K_hrr_w` is a
300 K reaction-enthalpy diagnostic, not an additional heat source. The finite-rate
`ZoneKinetics` operator is available separately; this example deliberately
states its instantaneous complete-oxidation closure instead of inventing a
calibrated methane rate law.

## Authored configuration

| Input | Value |
|---|---:|
| Rigid room floor area × height | 0.1 m² × 0.2 m = 20 litres |
| Initial gas/ambient state | 300 K, 101325 Pa, dry 0.79/0.21 N2/O2 mole mixture |
| Solid area × thickness | 0.01 m² × 0.006 m |
| Solid initial density / mass | 850 kg/m³ / 0.051 kg |
| Solid heat capacity | 1500 J/(kg K) |
| Solid conductivity | 0.2 W/(m K) |
| Solid temperature domain | 250–1000 K |
| Arrhenius prefactor / activation energy | 1e8 s⁻¹ / 100000 J/mol |
| Char mass yield | 0 |
| Face convection coefficient / gray emissivity | 10 W/(m² K) / 0.9 |
| External irradiation / radiative environment | 50000 W/m² / 300 K |
| Door height / discharge coefficient | 0.2 m / 0.7 |
| Closed / restricted / open door width | 0 / 0.001 / 0.02 m |

The gas occupies one active layer, the well-mixed limit of the original two-zone
model. The slab is an external fuel source thermally coupled to a fixed-volume
gas vessel. Its lost volume is not added to the vessel volume; an immersed
moving gas/solid interface would need a separate volume and pressure-work model.
Other room walls are adiabatic. The small sealed room deliberately becomes
highly pressurized; this is a rigid control-volume example, not a structural
model of a room or vessel.

The `fan-exhaust` case adds an explicitly authored fan and complete duct/outlet
loss to the restricted door. It uses the original static-to-total pressure
conversion, actual donor-density scaling and pressure-dependent operating-point
solve. Shaft heat is deposited downstream into ambient, so it is reported
separately from the finite solid/room energy balance. The supplied curve is a
fixture, not a measured product rating. Out-of-curve operation refuses; the
solver does not extrapolate or clamp it.

Its source curve runs from 100 Pa static rise at zero flow to zero static rise
at 0.01 m³/s, with one fan at unit speed. The reference density is the initial
room density, per-fan outlet area is 0.002 m², and the complete quadratic loss
resistance is 1e8 Pa s²/m⁶. Total efficiency is 0.5 and rejected shaft-loss
fraction is zero. The admitted donor/reference density ratio is 0.05–2, and
the low-pressure domain limits the pressure fraction to 0.02. These are explicit
authored fixture choices; no uncertain fan data are inferred.

## Completed 60-second receding-slab runs

The [retained actual outputs](../../../../data/reference/compartment_fire/receding_pyrolyzing_compartment_example.json)
contain all four default cases, a fan-case timestep refinement and a separate
fan-case mesh refinement. At maximum timestep 0.02 s and 24 solid cells:

| Quantity | Closed | Restricted door | Open door | Fan exhaust |
|---|---:|---:|---:|---:|
| Released fuel (g) | 12.260347 | 9.822454 | 11.424330 | 10.149217 |
| Burned fuel (g) | 1.368521 | 0.378477 | 2.300467 | 0.645526 |
| Fuel remaining in gas (g) | 10.891826 | 5.355061 | 0.922996 | 3.615231 |
| Final gas temperature (K) | 1452.101 | 654.250 | 1449.497 | 826.878 |
| Final pressure (Pa) | 951799.424 | 101327.134 | 101324.411 | 101324.008 |
| Final slab thickness (mm) | 4.557606 | 4.844417 | 4.655961 | 4.805974 |
| Surface recession (mm) | 1.442394 | 1.155583 | 1.344039 | 1.194026 |
| First oxygen-limited accepted time (s) | 17.92 | 12.30 | 12.29 | 12.30 |

All four cases pass the executable mass, canonical-species, elemental,
physical-energy, clock and pressure checks. The largest combined mass residual
was 1.74e-16 kg, elemental residual 5.13e-14 mol, combined energy residual
9.62e-7 J, and pressure residual 1.00e-5 Pa. The open case required eight
whole-step halvings; the other three required none. Every accepted interval
remapped the solid, with no geometry retries in these runs. Accumulated absolute
remap energy residuals were below 6.31e-9 J. Those residuals measure conservation
closure, not physical accuracy or the error introduced by redistribution.

The exhaust fan ended at 9.67963e-4 m³/s and 0.0413375 W shaft power. It
exhausted 0.0250541 kg of gas and delivered 3.69815 J of shaft energy downstream
into ambient. Open cases export unburned fuel, so burned plus remaining fuel
does not by itself equal released fuel.

The two selected refinements keep all other inputs unchanged:

```sh
cargo run --release -p fs-airflow --features compartment-fire \
  --example pyrolyzing_compartment_fire -- fan-exhaust 60 0.01 24 50000 receding
cargo run --release -p fs-airflow --features compartment-fire \
  --example pyrolyzing_compartment_fire -- fan-exhaust 60 0.02 48 50000 receding
```

| Change from the 24-cell, 0.02 s fan case | 0.01 s timestep | 48 solid cells |
|---|---:|---:|
| Released fuel | −0.002116% | +0.602647% |
| Burned fuel | +0.125202% | +0.064552% |
| Final gas temperature | +0.033386 K | −0.561744 K |
| Surface recession | −0.0000253 mm | +0.007196 mm |

Both runs pass the same checks. Their largest combined energy residuals were
1.53e-6 J and 9.21e-7 J, respectively. These are two sensitivity checks, not a
convergence-order or mesh-independence claim. The remap's independent quadratic
enthalpy refinement and datum-invariance checks live with the thermal owner.

CSV retains the original fields and appends `solid_thickness_m` and
`surface_recession_m`. In receding mode, `minimum_retained_cell_fraction` compares
current equal-mass spatial slots to their initial masses; slot identifiers no
longer label individual material parcels. Actual remap and geometry-work counts
appear in the final stderr summary.

## Retained 30-second fixed-geometry comparison

The [retained actual outputs](../../../../data/reference/compartment_fire/pyrolyzing_compartment_example.json)
contain the earlier four fixed-geometry cases and the fan case with a halved
timestep. Their reproduction command now explicitly selects the old mode:

```sh
cargo run --release -p fs-airflow --features compartment-fire \
  --example pyrolyzing_compartment_fire -- all 30 0.02 24 50000 fixed
```

At maximum timestep 0.02 s and 24 solid cells:

| Quantity | Closed | Restricted door | Open door | Fan exhaust |
|---|---:|---:|---:|---:|
| Released fuel (g) | 4.388407 | 3.638615 | 3.972363 | 3.702770 |
| Burned fuel (g) | 1.368521 | 0.377853 | 1.104748 | 0.492334 |
| Fuel remaining in gas (g) | 3.019886 | 2.740425 | 0.681380 | 2.104442 |
| Final gas temperature (K) | 2093.939 | 938.699 | 1574.367 | 1063.739 |
| Final pressure (Pa) | 945352.438 | 101325.557 | 101324.393 | 101323.803 |
| First oxygen-limited accepted time (s) | 17.74 | 12.24 | 12.225 | 12.22 |
| Minimum final cell mass / initial cell mass | 0.001461 | 0.055015 | 0.017762 | 0.046164 |

All cases pass the executable mass, species, elemental, physical-energy, clock
and pressure checks. The largest combined mass residual was 9.37e-17 kg, largest
elemental residual 1.65e-14 mol, largest combined energy residual 2.90e-6 J, and
largest pressure-equation residual 1.00e-5 Pa. The open-door case required 11
whole-step halvings; the other three completed without retries. Burned fuel plus
remaining fuel need not equal released fuel in open cases because the actual
vents also export unburned gas.

The exhaust fan's final volume flow was 9.67144e-4 m³/s and shaft power
0.0359141 W. Over the run it exhausted 0.0173896 kg of gas and delivered 2.51170 J
of shaft energy downstream into ambient. Those powers are not added again to
the room balance.

Halving only the fan case's maximum timestep to 0.01 s changed released fuel
by about 0.0173%, burned fuel by 0.172%, and final gas temperature by 0.1942 K.
All checks still passed; the maximum combined energy residual was 1.50e-5 J.
This is a single temporal refinement on the same mesh, not a convergence-order
or mesh-independence claim. A `closed 33 0.02 24 50000 fixed` scope check deliberately
exits with the named retained-mass/recession refusal; it is not a completed
33-second prediction.

## Numerical and physical scope

Every accepted step checks combined solid/gas mass, canonical gas species,
C/H/N/O elements, physical energy, common clock and pressure residual. The closed
case also checks the independent oxygen-capacity limit from the initial O2
inventory. Tiny balance residuals alone do not establish time or mesh accuracy.

The default receding model admits a homogeneous, noncharring slab with explicit
constant density and conductivity. Geometry retries bound the fractional mass
loss per solid interval. It preserves actual volatile mass and enthalpy while
remapping retained material; it does not model char porosity, swelling, stress,
thermal expansion, solid oxidation or complete burnout. The fixed comparison
keeps its original 1e-4 initial-cell-mass refusal. Neither mode inserts residual
mass, a heat-capacity floor or release clipping, and neither promises that every
admitted duration will stay within the supplied caloric/model domain.

Implicit solid stepping and geometry remapping remain first order in time;
the boundary and gas-source coupling is also first order. Symmetric composition
would not by itself establish a second-order method or an entropy theorem.
The instantaneous,
well-mixed oxidation law does not predict ignition, extinction, flame structure,
CO or soot yields, or high-temperature dissociation. The external optical
boundary does not model smoke/flame radiative feedback. The separate
`radiating_compartment` example exercises participating radiation coupled to
finite walls. Sprays and real-gas phases remain outside this example.
