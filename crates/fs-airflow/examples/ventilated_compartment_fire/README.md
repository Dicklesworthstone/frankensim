# Ventilated reacting compartment

This example connects real species inflow, pressure-driven doorway transport,
and oxygen-limited chemical conversion using the existing `fs-gas` and
`fs-airflow` owners. Its purpose is to exercise conserved reacting ventilation,
including starvation and unburned fuel accumulation. It is a G0/G1 numerical
demonstration, not a calibrated or experimentally validated compartment-fire
prediction.

Run from the repository root:

```sh
cargo run -p fs-airflow --features compartment-fire --example ventilated_compartment_fire -- 60 0.02
```

The optional arguments are end time in seconds (default 60), maximum timestep
in seconds (default 0.02), and source coupling (`coupled`, the default, or
`split`). The latter retains the earlier post-pressure-solve source composition
for a direct numerical comparison. End time is bounded by 120 s and the maximum
timestep must lie in 0.001–0.25 s. Gas states outside the declared 200–3500 K
caloric domain refuse; they are never clipped.

## What is coupled

All three cases have a rigid 1 m³ room initially containing dry air at 300 K and
101325 Pa, with an explicitly authored 0.79/0.21 N2/O2 mole mixture. A pure methane
feed enters at 0.002 kg/s and 300 K. All gas occupies one active layer, the
well-mixed limit of the existing two-zone model. This isolates reacting
ventilation from an additional unvalidated interlayer mixing law.

The closed case has exactly zero opening width. Restricted ventilation uses a
0.02 m wide doorway, and the adequately ventilated case uses a 0.2 m wide
doorway. Both are 1 m high and use the existing piecewise hydrostatic
bidirectional flow and pressure solve, with nominal discharge coefficient 0.7.
The restricted case has no imposed combustion-power or air-supply cap: its
actual transported oxygen inventory limits burning through the same reaction.

The default pressure-coupled timestep is:

1. For a trial opening pressure, transport the original donor species and
   partial enthalpy, including methane feed and doorway inflow/outflow.
2. Apply `InventoryLimitedReaction::apply_split` to the resulting actual
   inventory with `CH4 + 2 O2 -> CO2 + 2 H2O`.
3. Recompute the opening-pressure residual from the reacted gas. Every Newton,
   finite-difference and line-search evaluation repeats these operations from
   the accepted inventories; sources never accumulate across rejected trials.
4. Publish the converged network and its matching source result together.

`ZoneNetwork::step_with_pressure_sources` reuses the original dense Newton/LU
solver, source admission and balance/publication path. Its preparation is
immutable during the pressure solve; accepted diagnostics are returned from the
actual converged source candidate. The optional `split` method closes pressure
before applying chemistry and separately reports the resulting pressure jump.

The reaction requests conversion of the methane actually remaining in the
room; the original reaction adapter limits progress by the available reactants.
Oxygen that has already left through a vent cannot also burn fuel. The source
step preserves the physical clock, so transport and chemistry do not count the
same interval twice.

The methane feed brings its original **enthalpy**, including its formation
reference. Chemical conversion changes the sensible-energy gauge consistently
with the changed species inventory and preserves physical chemical-plus-thermal
energy. The reported released enthalpy is evaluated at an explicit 300 K
reference temperature and is a diagnostic. It is never added as another heat
source.

The reaction adapter retains a mechanism object to reuse exact elemental/charge
and molecular-mass admission. Its required unit Arrhenius-law carrier is explicitly
unused; this example does not evaluate or claim a methane kinetic rate.

## Original thermochemical data

Five original NASA9 species are retained: CH4, CO2, H2O, N2, and O2. The original
molar masses and the 200–1000 K and 1000–6000 K coefficients come from
[NASA CEA thermo.inp](https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp),
commit `2f79a647e85737742fb08a127da00d1ea1d069da`, blob
`5ddeb5b213701d0e8b79dbda7b48ede14e31404a`.

The standard pressure is **1 bar**, as specified in
[NASA TP-2002-211556](https://ntrs.nasa.gov/citations/20020085330). The source
convention is distinct from the room's initial pressure.

The existing explicit `with_caloric_continuity` owner anchors the low-temperature
region at 300 K. It retains every source Cp coefficient and all low-temperature
enthalpy/entropy values. Only bounded constant h/s offsets are applied to the
high-temperature region. The admitted budgets are 0.01 J/mol and
0.00001 J/(mol K), and the JSON output prints the actual offsets for each species.
These budgets bound the authored seam adjustment, not experimental uncertainty.

## Numerical accounting

The output retains canonical species order and reports:

- Final species inventories, injected and reacted fuel, oxygen-limit onset,
  temperature, pressure, and the diagnostic heat-release rate.
- Actual net external species transfer, separated from internal chemical
  conversion, and species-balance residuals.
- C/H/N/O elemental balances including inlet and outlet transfer.
- Physical energy balance and the chemical source's energy neutrality.
- Maximum mass/species residuals, mass-fraction sum error, pressure-equation
  residual, and the separate pressure change caused by the chemical split.
- Accepted step and retry counts, accepted pressure-iteration totals, actual
  source-evaluation counts, and the terminal accepted duration. The source count
  includes discarded trials; the iteration total counts accepted solves.

The closed case also checks a separate stoichiometric oracle: reacted methane
must equal the smaller of total supplied methane and
`initial_O2_mass * M_CH4 / (2 * M_O2)`. For the default initial state this oxygen
capacity is 0.068426044707 kg CH4; at 0.002 kg/s it is reached after
34.213022354 s. At 60 s, 0.051573955293 kg CH4 must remain unburned in the closed
room. This follows from the initial inventory, independently of the caloric
solver or the timestep.

## Physical and numerical limits

The walls are adiabatic and rigid. The sealed case deliberately becomes hot and
pressurized; it is a control-volume limit, not a building structural prediction.
The single-step complete-oxidation closure assumes instantaneous mixing and
burning whenever both reactants exist. It does not describe ignition,
extinction, finite-rate chemistry, flame structure, high-temperature
dissociation, CO/soot formation, pyrolysis, radiation, or condensation of product
water. Temperature-dependent ideal-gas heat capacities alone do not supply
those missing mechanisms.

The transport/chemistry composition remains first order. Coupling the source
into pressure closure removes the subsequent chemical pressure jump, while
donor temperature, composition and interface geometry in the transport flux
remain frozen over the interval. Small algebraic conservation residuals do not
establish temporal accuracy or agreement with a real fire experiment.


## Executed 60-second comparison

The retained [actual numerical outputs](../../../../data/reference/compartment_fire/ventilation_limited_methane_example.json)
contain all three pressure-coupled cases at maximum timesteps 0.02 and 0.01 s,
plus explicitly labeled two-case split baselines. Every coupled case completed
without timestep retries and passed the executable species, element,
physical-energy, mass-fraction, and closed-room stoichiometric checks.

Results at a maximum timestep of 0.02 s:

| Quantity | Closed | Restricted ventilation | Adequate ventilation |
|---|---:|---:|---:|
| Opening width (m) | 0 | 0.02 | 0.2 |
| Final temperature (K) | 2493.122425 | 1660.925889 | 1160.910333 |
| Final pressure (Pa) | 997107.918804 | 101321.905955 | 101321.702066 |
| Reacted methane (kg) | 0.068426044707 | 0.042634673926 | 0.120000000000 |
| Remaining methane (kg) | 0.051573955293 | 0.026355515959 | 0 |
| Methane carried out (kg) | 0 | 0.051009810116 | 0 |
| Remaining oxygen (kg) | 0 | 0 | 0.044724870270 |
| First oxygen-limited step end (s) | 34.22 | 8.96 | None |
| Final diagnostic release (W) | 0 | 23429.525433 | 100051.883345 |
| Physical energy imbalance (J) | 4.15e-7 | 7.28e-10 | -1.06e-8 |

The restricted case demonstrates continued burning sustained by incoming
oxygen after the resident oxygen is exhausted. Its final release is about
23.43 kW against the approximately 100.05 kW complete-oxidation potential of the
fuel feed. The 0.12 kg supplied methane is accounted for as reacted fuel, fuel
remaining in the room, and unburned fuel leaving through the opening. It burns
less cumulative fuel than the sealed case because ventilation also carries
oxygen out during the early transient. The sealed case eventually stops
reacting, while the adequately ventilated case burns all supplied methane.
These outcomes come from the actual inventories and opening flow, with no
additional combustion-power limit.

Halving the maximum timestep changes the restricted case's final temperature
by 0.040590 K, pressure by -0.0000083 Pa, and release by 0.06743 W. Its reacted
fuel becomes 0.042651636247 kg and oxygen limitation is first reported at
8.97 s; onset is resolved only to the accepted step. The adequately ventilated
case's final temperature changes by less than 3e-10 K and pressure by less than
7e-9 Pa. These are late-time results for this example, not a global convergence
order claim. The independently computed closed-volume NASA energy balance
gives 2493.122425259 K, consistent with both executed closed temperatures.

Across all six coupled runs, physical-energy imbalance is below 9.24e-7 J,
final species imbalance below 3.09e-14 kg, and elemental imbalance below
4.31e-12 mol. The full post-source pressure residual is below 1e-5 Pa, and
there is no additional chemical pressure jump after acceptance. Numerical
precision here does not measure physical model accuracy.

### Pressure coupling compared with the retained split method

The adequately ventilated case exposes the pressure error of applying a hot
chemical source after closing the vent pressure. Both methods use the same
reaction, feed, caloric data and opening law:

| Method and maximum timestep | Final pressure (Pa) | Maximum post-solve pressure change (Pa) | Accepted pressure iterations | Source evaluations | Timestep halvings |
|---|---:|---:|---:|---:|---:|
| Split, 0.02 s | 101939.152732 | 797.187903 | 24417 | 3009 | 18 |
| Coupled, 0.02 s | 101321.702066 | 0 | 5804 | 11412 | 0 |
| Split, 0.01 s | 101475.989971 | 398.657180 | 43032 | 6007 | 12 |
| Coupled, 0.01 s | 101321.702066 | 0 | 11025 | 21073 | 0 |

The split residual measures transport-only pressure closure; a small value
does not include the subsequent chemistry. The coupled method evaluates more
chemical candidates, but closes the pressure equation on the reacted state
and uses fewer accepted pressure iterations in this example. These are actual
algorithmic counts, not a controlled runtime benchmark. The retained output
also records terminal accepted durations, since instantaneous split pressure
depends on that duration.

Verification used a native Cargo build of the example and all 66 `fs-airflow`
library tests, including eight focused pressure-source tests, against the
original-source gas, airflow, thermochemistry, quantity, material and evidence
owners. The selected-owner native manifest excludes unrelated workspace
dependencies; this is not a full-workspace Cargo or DSR pass.
