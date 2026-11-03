# Speciated gold leaching and finite carbon contact

`fs_material::cyanide_leach` connects the existing NBS acid/base data to a
finite, isothermal gold-leach batch. It solves HCN/CN, water, and any other
declared acid families at every kinetic stage. The evolving dicyanoaurate
inventory participates in electroneutrality, so free cyanide and pH change
with dissolution. Carbonate buffering can use the matched profiles in
[MINING_CARBONATE_BUFFERING.md](MINING_CARBONATE_BUFFERING.md).

## Chemistry and rate are separate inputs

The existing exact atom/charge owner supplies the dissolved-oxygen Elsner
equation:

```text
4 Au(s) + 8 CN[-](aq) + O2(aq) + 2 H2O(l)
    -> 4 Au(CN)2[-](aq) + 4 OH[-](aq)
```

If `x` is moles of gold dissolved, the original equation extent is `x/4`.
The closed batch derives stocks from their initial values:

| Inventory | Current amount |
| --- | --- |
| Accessible solid gold | `A0 - x` |
| Uncomplexed HCN + CN | `T0 - 2*x` |
| Dissolved oxygen | `O0 - x/4` |
| Dicyanoaurate | `G0 + x` |

At each stage, the original `IdealAcidBaseSystem` determines how much of the
remaining uncomplexed pool is free CN. Counting only the initial free CN as
available reagent would incorrectly discard its HCN reservoir. The cyanide
balance includes two bound CN groups per dicyanoaurate; it does not describe
copper complexes, WAD-CN, oxidation losses or other cyanicide minerals.

A distinct source card declares the empirical law
`r_Au = k * n_accessible_Au * c_free_CN * c_dissolved_O2`, with
`k` in m6/(mol2 s). The actual pH, free-CN concentration, dissolved-O2
concentration, gold inventory, elapsed time, volume and temperature must be
inside that card's applicability box. Gold availability is an empirical
inventory term, not a solid activity or an inferred reactive surface.

The NBS equilibrium data do **not** supply `k` or establish these reaction
orders. Primary ore studies report their own fitted orders and source
conditions: [de Andrade Lima and Hodouin (2005)](https://doi.org/10.1016/j.hydromet.2005.06.001)
and [Ling et al. (1996)](https://doi.org/10.1016/0008-4433(96)00003-1).
Those studies motivate separate reagent inputs; their fitted laws are not
silently substituted for this explicitly declared first-order law.

Activities use the caller's common `amount_at_unit_activity`, whereas rate
concentrations use independently supplied liquid volume. For the NBS molal
standard state, 10 kg solvent gives `n_ref = 10 mol`; it does not follow from
having 0.01 m3 of solution. Water retains fixed unit activity. Stoichiometric
water consumption and OH production are distinct from equilibrium OH and
do not close a finite H/O solvent balance or a bath heat calculation.

## Oxygen contact while gold continues to dissolve

`cyanide_leach::aeration::AeratedGoldCyanidationBatch` adds a maintained,
prescribed oxygen partial pressure to the same chemical batch. It reuses the
existing `IdealGasLiquidPartition` and the original NBS gas/aqueous oxygen
profiles. For `O2(g) -> O2(aq)`, the ideal activity laws give

```text
n_DO_sat = n_ref * K * p_O2 / p_ref
K = exp(-Delta_r_G_standard / (R*T))
```

The gas species' partial pressure is supplied explicitly. The standard
reference pressure stays with the original thermochemical source; it is not
a total-air pressure, a humidity correction or a measured pressure envelope.
The solvent-based amount at unit activity remains distinct from liquid volume.
These Henry-law conventions follow the
[IUPAC recommendations](https://doi.org/10.1515/pac-2020-0302).

At 298.15 K, the matched NBS profiles give a standard dissolution Gibbs energy
of 16,400 J/mol and `K = 0.0013391965031285026`. For `n_ref = 10 mol` and
`p_O2 = 21,000 Pa`, the target is approximately 0.002812312657 mol dissolved
oxygen. With a separately declared 0.01 m3 liquid volume, this is
0.281231266 mol/m3. Source compilation, original phase ownership, equilibrium
activities and pressure/basis scaling are exercised by
`xtask/tests/mining_gas_solubility_sources.rs`.

A **separate** source card supplies `kLa` in inverse seconds. At fixed volume,
the explicitly adopted well-mixed transfer closure is

```text
J = kLa * (n_DO_sat - n_DO)
dx/dt = r_Au
d(n_DO)/dt = J - r_Au/4
```

The transfer structure follows
[US EPA, Aeration Equipment Evaluation, section 4, equation 16](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100GGR5.TXT).
That clean-water report does not calibrate a slurry or ore-specific aerator.
No transfer coefficient is inferred from Gibbs energy, source names or
stoichiometric reaction coefficients.

Dissolved oxygen is an independently evolved amount. Starting at zero DO
therefore permits subsequent leaching after oxygen enters. The chemical
ceiling remains `min(initial_accessible_Au, initial_uncomplexed_CN / 2)`;
the initial oxygen inventory no longer imposes a ceiling when transfer is
active. Free CN and pH still come from simultaneous speciation at every
kinetic stage, with evolving dicyanoaurate in the charge balance.

Gross oxygen uptake and gross stripping have independent cumulative and
interval accounts. The original Elsner equation supplies consumption, so
the independently checked oxygen balance is

```text
current_DO + gold_dissolved/4 + cumulative_O2_out
    = initial_DO + cumulative_O2_in
```

Pressure is fixed during one advance and may change on the next advance.
At `kLa = 0`, exchange stops. At zero oxygen partial pressure and positive
`kLa`, the boundary strips oxygen; zero pressure is not an exchange-off
switch. The prescribed reservoir has no finite headspace or gas-feed cap.
For a finite gas inventory and metered oxygen additions, the separate
`cyanide_leach::gas_transfer` owner retains those explicit stocks.

The actual pH, free CN, DO, time, volume and temperature remain subject to
the original chemical rate card. The transfer card separately admits actual
partial pressure, DO, time, volume and temperature, with optional saturation
and activity-basis bounds. Caller-supplied extras cannot overwrite owned
conditions. The external gas species cannot also be a liquid acid-family
member or spectator.

Every advance is atomic and bounded by work/attempt limits and cancellation.
Full/two-half RK4 trials independently control raw local differences in gold,
DO, gross uptake and gross stripping. These estimates are not certified
global error bounds. Positive inventories are not clipped, and the separately
reported oxygen balance allowance does not modify any stock.

## Run the connected calculation

The `gold_leach_to_carbon` example loads eleven actual NBS source profiles,
retains their receipt-backed equilibrium constants and gold molar mass, and
reports the direct standard Gibbs energy of the gold reaction. It compares
lower/higher alkalinity, limited dissolved oxygen and carbonate buffering.
Each batch advances through three intervals while retaining its accepted
stocks and clock.

```bash
gold_process_packs=$(mktemp -d)
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry gold-solid-nbs1989-thermochemistry dicyanoaurate-ion-aqueous-nbs1989 oxygen-aqueous-nbs1989-thermochemistry carbonic-component-aqueous-nbs1989 bicarbonate-ion-aqueous-nbs1989 carbonate-ion-aqueous-nbs1989; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$gold_process_packs/$pack.fsmatpk"
done
cargo run -p fs-material --example gold_leach_to_carbon -- \
  --materials "$gold_process_packs"
cargo test -p fs-material cyanide_leach
```

The oxygen-limited case starts with 0.001 mol accessible Au but only
0.00005 mol dissolved O2. Its gold-dissolution ceiling is therefore 0.0002 mol,
or 20% of the accessible gold, irrespective of how much cyanide is present.
No atmospheric oxygen replenishment is assumed.

The resulting dissolved analytical Au then enters the existing finite-time
Langmuir/LDF carbon owner. Activated and natural carbon compete for that
finite amount. The output separates their gold uptake and checks the
combined unleached-solid, liquid and adsorbed gold balance. All example
inventories, rates, capacities and applicability boxes are clearly authored
illustrations; no ore or plant calibration is claimed.

This carbon handoff is deliberately on the existing **analytical Au** basis.
The carbon stage treats pH as an externally held condition and does not
return a fully speciated cyanide liquor: no ligand/counterion adsorption law
has been supplied. Simultaneous CIL chemistry, cyanide/alkalinity dosing, through-flow, heat,
ore-specific liberation and measured kinetic calibration remain separate
capabilities.

Numerical advances have bounded work and cancellation checks and publish
only after the full requested interval succeeds. Local and accumulated error
estimates are numerical diagnostics, not certified global bounds or source
uncertainty. The exact coupled fixture in the unit tests holds pH at 7 while
gold, free CN, HCN and oxygen decline; it has the independent solution
`n_Au(t) = 1/sqrt(1 + k*t/2)` under the fixture's declared unit scales.

## Run the aerated calculation

The optional mode preserves the original eleven-pack invocation. Compile the
additional existing oxygen gas profile into the same directory, then run:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/oxygen-gas-nbs1989-thermochemistry/manifest.tsv \
  --out "$gold_process_packs/oxygen-gas-nbs1989-thermochemistry.fsmatpk"
cargo run -p fs-material --example gold_leach_to_carbon -- \
  --materials "$gold_process_packs" --aerated
```

This mode runs three explicitly authored cases at 21,000 Pa oxygen partial
pressure and 298.15 K:

| Case | Initial dissolved O2 | Authored kLa |
| --- | ---: | ---: |
| Sealed control | 0.00005 mol | 0 /s |
| Weak oxygen contact | 0 mol | 0.0001 /s |
| Stronger oxygen contact | 0 mol | 0.01 /s |

Each starts with 0.001 mol accessible gold and 0.01 mol uncomplexed cyanide
in the independently specified 0.01 m3 liquid and 10 kg solvent basis.
The chemical source box explicitly covers DO up to 0.4 mol/m3 in this mode.
The cases differ in their declared initial oxygen and transfer inputs; they
are not a fitted comparison of commercial aerators.

The example prints the source-derived saturation, actual pH/free CN/DO,
oxygen entering and leaving, reaction consumption, conservation residuals,
numerical estimates and bounded-work counts. Its three accepted intervals
end at 10, 50 and 200 seconds. The resulting dissolved analytical Au then
enters the same finite carbon-contact calculation, with separate activated
and natural carbon uptake and a complete analytical-gold balance.

The aeration tests use independent transfer-only exponentials, the closed
leach's analytic solution, initially oxygen-free Taylor coefficients, pressure
changes and uptake/stripping reversal. They also exercise amount/volume/basis
scaling, actual-source refusals and rollback after cancellation or exhausted
budgets. The example's rate and transfer cards remain synthetic; use measured
source cards with applicable ore, slurry and equipment conditions before
treating the results as a process prediction.
