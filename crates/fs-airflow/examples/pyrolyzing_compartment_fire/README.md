# Heated solid feeding a reacting, ventilated compartment

This example couples the existing spatial solid-decomposition model to the
existing reacting gas and pressure-driven ventilation models. Heating determines
the fuel release; neither a fuel-feed history nor combustion power is supplied.
It exercises oxygen starvation, unburned fuel accumulation, natural ventilation
and pressure-dependent mechanical exhaust with combined conservation accounts.

The material is an **explicitly manufactured methane-releasing solid**. It is
not PMMA, wood, or a calibrated experimental fire. The source-backed PMMA/TGA
examples remain separate: their effective outgoing enthalpy cannot simply be
substituted for the original gas's physical formation-enthalpy reference.

## Run

```sh
cargo run --release -p fs-airflow --features compartment-fire \
  --example pyrolyzing_compartment_fire -- all 30 0.02 24 50000
```

Arguments are case (`all`, `closed`, `restricted`, `open`, `fan-exhaust`), end
time in seconds, maximum timestep in seconds, solid cell count, and external
irradiance in W/m². The values above are the defaults. CSV goes to stdout;
source declarations and final balance summaries go to stderr. The face
temperature in each row is the lagged face temperature used over the accepted
interval, while the cell and gas temperatures are the accepted endpoint values.

The example bounds work and admits at most 96 cells, 120 s, 100 kW/m² external
irradiance, and 200,000 accepted steps. The minimum network step is 1e-5 s.
Those are numerical input ceilings, not a guarantee that every such input stays
within the model's physical domain.

## Existing owners and energy reference

The solid uses the original `fs-thermal::enthalpy::pyrolysis::PyrolysisSimulation` on the
original deposited-cell thermal graph. Its ordinary links conduct heat through
the slab; its Arrhenius reaction changes retained mass and exports integrated
volatile mass and enthalpy. The existing `VolatileMixture` uses the original
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
   temperatures. The front half-cell conductance balances convection, absorbed
   external irradiation, and radiation to a 300 K blackbody environment.
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
model. Other room walls are adiabatic. The small sealed room deliberately becomes
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

## Completed 30-second runs

The [retained actual outputs](../../../../data/reference/compartment_fire/pyrolyzing_compartment_example.json)
contain all four default cases and the fan case with a halved timestep. At
maximum timestep 0.02 s and 24 solid cells:

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
or mesh-independence claim. A `closed 33 0.02 24 50000` scope check deliberately
exits with the named retained-mass/recession refusal; it is not a completed
33-second prediction.

## Numerical and physical scope

Every accepted step checks combined solid/gas mass, canonical gas species,
C/H/N/O elements, physical energy, common clock and pressure residual. The closed
case also checks the independent oxygen-capacity limit from the initial O2
inventory. Tiny balance residuals alone do not establish time or mesh accuracy.

Solid geometry and conductance stay fixed as material decomposes. The example
refuses when any cell retains less than 1e-4 of its initial mass. A receding
surface/remeshing model is required to continue through depletion; no residual
mass, heat-capacity floor, or release clipping is introduced. This scope limit
can be reached earlier on a finer mesh, so the maximum permitted input duration
must not be read as a guaranteed run duration.

The boundary and source composition is first order. The instantaneous,
well-mixed oxidation law does not predict ignition, extinction, flame structure,
CO or soot yields, or high-temperature dissociation. The external optical
boundary does not model smoke/flame radiative feedback. No pores, oxidation of
the solid, char layer, shrinkage, sprays, or real-gas phases are inferred.
