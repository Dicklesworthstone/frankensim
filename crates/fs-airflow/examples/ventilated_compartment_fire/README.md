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

The optional arguments are end time in seconds (default 60) and maximum timestep
in seconds (default 0.02). End time is bounded by 120 s and the maximum timestep
must lie in 0.001–0.25 s. Gas states outside the declared 200–3500 K caloric domain
refuse; they are never clipped.

## What is coupled

Both cases have a rigid 1 m³ room initially containing dry air at 300 K and
101325 Pa, with an explicitly authored 0.79/0.21 N2/O2 mole mixture. A pure methane
feed enters at 0.002 kg/s and 300 K. All gas occupies one active layer, the
well-mixed limit of the existing two-zone model. This isolates reacting
ventilation from an additional unvalidated interlayer mixing law.

The closed case has exactly zero opening width. The ventilated case has a
0.2 m wide, 1 m high wall opening using the existing piecewise hydrostatic
bidirectional flow and pressure solve, with nominal discharge coefficient 0.7.

Each accepted timestep is:

1. Transport the original donor species and partial enthalpy, including methane
   feed and any incoming or outgoing doorway gas.
2. Apply `InventoryLimitedReaction::apply_split` to the resulting actual
   inventory with `CH4 + 2 O2 -> CO2 + 2 H2O`.
3. Publish the accepted network and source result together.

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
- Accepted step and retry counts.

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

The transport/chemistry composition is first order. Its separate post-transport
chemical pressure change is reported and should decrease with timestep
refinement. Small algebraic conservation residuals do not establish temporal
accuracy or agreement with a real fire experiment.


## Executed 60-second comparison

The retained [actual numerical outputs](../../../../data/reference/compartment_fire/ventilation_limited_methane_example.json)
contain both cases at maximum timesteps 0.02 and 0.01 s. Both completed and passed
the executable species, element, physical-energy, mass-fraction, and closed-room
stoichiometric checks.

| Quantity | Closed, 0.02 s | Closed, 0.01 s | Ventilated, 0.02 s | Ventilated, 0.01 s |
|---|---:|---:|---:|---:|
| Final temperature (K) | 2493.122425259 | 2493.122425258 | 1160.910333075 | 1161.280178600 |
| Reacted methane (kg) | 0.068426044707 | 0.068426044707 | 0.120000000000 | 0.120000000000 |
| Remaining methane (kg) | 0.051573955293 | 0.051573955293 | 0 | 0 |
| Remaining oxygen (kg) | 0 | 0 | 0.044997421957 | 0.044752258791 |
| Final diagnostic release (W) | 0 | 0 | 100051.883345 | 100051.883345 |
| Physical energy imbalance (J) | 4.15e-7 | -9.24e-7 | -1.04e-8 | -1.97e-8 |

The maximum elemental inventory defect across these runs is 4.39e-12 mol;
the largest final species defect is below 3.0e-14 kg. The independently computed
closed-volume NASA energy balance gives 2493.122425259 K, consistent with both
executed temperatures. Numerical precision here does not measure physical
model accuracy.

The ventilated final temperature changes by 0.36985 K on timestep halving,
and its remaining oxygen changes by 0.00024516 kg. Maximum chemical split
pressure changes fall from approximately 797.2 Pa to 398.7 Pa. The pressure
solver itself has residual below 1e-5 Pa, but that is the **transport** closure;
it does not remove the subsequent chemical pressure increment. Final reported
pressures also depend on the accepted terminal step. This example does not
establish converged instantaneous compartment pressure or a second-order global
time integrator.

Verification used a narrow native build of this example against the existing
original-source gas, airflow, thermochemistry, quantity, material, and evidence
libraries because the environment lacked full-workspace build capacity. It is
not a full-workspace Cargo or DSR result.
