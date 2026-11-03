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
The batch derives stocks from their initial values:

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
has been supplied. Simultaneous CIL chemistry, oxygen mass transfer, reagent
dosing, through-flow, heat, ore-specific liberation and measured kinetic
calibration remain separate capabilities.

Numerical advances have bounded work and cancellation checks and publish
only after the full requested interval succeeds. Local and accumulated error
estimates are numerical diagnostics, not certified global bounds or source
uncertainty. The exact coupled fixture in the unit tests holds pH at 7 while
gold, free CN, HCN and oxygen decline; it has the independent solution
`n_Au(t) = 1/sqrt(1 + k*t/2)` under the fixture's declared unit scales.
