# Finite hydrated lime and shared cyanide pH

The NBS calcium profiles support a finite crystalline Ca(OH)2 inventory in the
same reference convention as the existing HCN/CN− and water sources.
[The `lime_buffered_cyanide` example](../../../crates/fs-thermochem/examples/lime_buffered_cyanide.rs)
solves dissolved calcium, residual solid, pH and free cyanide together. Small
lime inventories can dissolve completely; an inventory exceeding the selected
solid's equilibrium capacity retains solid. Source receipts and the calcium,
cyanide, charge and saturation residuals remain visible.

## Primary values and source identity

The values are literal rows from [NIST's corrected NBS tables](https://doi.org/10.18434/M32124),
[official CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv),
Table 94:Ca, *Calcium (Prepared 1968)*. The
[original NBS publication](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf), sections
5.5–5.8, gives the standard-state conventions; the
[digitization paper](https://doi.org/10.6028/jres.125.007) explains the corrected
dataset. All claims apply at **298.15 K and 100000 Pa**.

| Profile | Literal state | M, kg/mol | ΔfH°, J/mol | ΔfG°, J/mol | S°, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `calcium-oxide-solid-nbs1989-thermochemistry` | CaO, `cr` | 0.0560794 | −635090 | −604030 | 39.75 |
| `calcium-hydroxide-solid-nbs1989-thermochemistry` | Ca(OH)2, `cr` | 0.0740948 | −986090 | −898490 | 83.39 |
| `calcium-ion-aqueous-nbs1989` | Ca2+, `ao` | 0.0400800 | −542830 | −553580 | −53.1 |

Each profile contains these four claims with `Unstated` uncertainty. Calcium's
negative entropy is **conventional ionic entropy**, on the H+(aq)-zero
reference. The solid rows use standard molar entropy. The crystalline Ca(OH)2
row has no named polymorph suffix, so the pack does not invent a
`source-phase-portlandite` qualifier.

All three require typed static temperature and pressure,
`source-elemental-reference-nbs1989=1`, `source-standard-state=1`,
`engineering-reference-model=1` and `source-specimen-qualified=0`.
The solids additionally require `source-phase-solid=1` and
`source-pure-species=1`. Calcium requires `source-phase-aqueous=1`,
`source-aqueous-unit-activity-infinite-dilution=1`,
`source-individual-ion=1` and `source-ionic-reference-hydrogen-zero=1`.

The existing JANAF/NIST lime and USGS calcium/portlandite profiles retain their
own evaluations. Their values cannot be relabeled as NBS participant sources.
For example, NBS CaO entropy is 39.75 J/(mol K), and NBS Ca2+ conventional
entropy is −53.1 J/(mol K); both differ from those separate evaluations.
The literal molar masses also retain the NBS table's atomic-weight basis.

## Conserved slaking and acid/base endpoints

The existing `LimeHydration` and `HydratedLimeNeutralization` reactions can
now use this complete NBS chain. `CalciumHydroxidePrecipitation` supplies
the reverse neutralization in the four-species hydrolysis form required by
`HydrolysisPrecipitation`.

| Reaction | Exact equation | ΔH°, J/mol extent | Direct ΔG°, J/mol extent | ΔS°, J/(mol K) |
| --- | --- | ---: | ---: | ---: |
| Lime hydration | CaO(s) + H2O(l) → Ca(OH)2(s) | −65170 | −57331 | −26.27 |
| Hydrated-lime neutralization | Ca(OH)2(s) + 2 H+ → Ca2+ + 2 H2O(l) | −128400 | −129348 | 3.33 |
| Calcium-hydroxide precipitation | Ca2+ + 2 H2O(l) → Ca(OH)2(s) + 2 H+ | 128400 | 129348 | −3.33 |

The water and H+ partners are the existing NBS profiles. Direct ΔG° uses each
source's independent formation-Gibbs column. The caloric route uses ΔH°−TΔS°
and remains separately reported; the tabulated routes need not agree exactly.
Slaking heat is reference bookkeeping for a declared reacted amount. It does
not calculate a slurry temperature or dissolution rate.

The finite-phase solver uses pure-water and solid unit activities and ideal
aqueous activities `a_i=n_i/n_ref`. It re-solves water ionization and HCN/CN−
distribution at the same pH as calcium/solid partition. For the dissolution
direction Ca(OH)2(s) → Ca2+ + 2 OH−, the literal NBS formation energies give
ΔG° = **30422 J/mol**, hence `Ksp=exp(-30422/(R*298.15))`. This is a derived
reference equilibrium constant, not a fitted solubility measurement.

A present solid requires `a_Ca*a_OH^2/Ksp=1` within the declared numerical
tolerance. A completely dissolved finite inventory must be undersaturated or
at saturation. The model reports the natural log of that ratio. An empty
calcium component has no defined saturation ratio.

## Run the actual-source example

Compile the nine example profiles and the additional CaO profile used by the
slaking CLI:

```bash
mkdir -p target/lime-packs
for pack in \
  hydrogen-cyanide-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  cyanide-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry \
  calcium-ion-aqueous-nbs1989 \
  calcium-hydroxide-solid-nbs1989-thermochemistry \
  sodium-ion-aqueous-nbs1989 \
  sodium-cyanide-cubic-solid-nbs1989-thermochemistry \
  calcium-oxide-solid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "target/lime-packs/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example lime_buffered_cyanide -- \
  --materials target/lime-packs
```

The default authored scenarios use 1 kg solvent and 0.01 mol already-dissolved
NaCN, with calcium-hydroxide component totals of 0, 0.001 and 0.1 mol.
The supplied lime is a finite solid-plus-dissolved component; complete
dissolution is an output of the declared model rather than an input assumption.
The already-dissolved NaCN supplies both analytical cyanide and its sodium
counterion. Actual-source molar masses convert the explicitly declared pure
reagent masses.

A single case can instead be requested by mass:

```bash
cargo run -p fs-thermochem --example lime_buffered_cyanide -- \
  --materials target/lime-packs \
  --lime-mass-kg 0.00740948 \
  --cyanide-mass-kg 0.000490077 \
  --solvent-mass-kg 1
```

The solvent mass fixes the molal activity basis: 1 kg gives `n_ref=1 mol`.
It is not a physical solution volume or an assumption that one mole of solvent
is present.

For source-bound reaction amounts and both Gibbs routes:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction lime-hydration-nbs1989 --equilibrium --gibbs \
  target/lime-packs/calcium-oxide-solid-nbs1989-thermochemistry.fsmatpk \
  target/lime-packs/water-liquid-nbs1989-thermochemistry.fsmatpk \
  target/lime-packs/calcium-hydroxide-solid-nbs1989-thermochemistry.fsmatpk \
  0.0560794
```

The final argument is kilograms of the first species **actually reacted**;
this invocation declares one mole of CaO reacted. The additional CLI modes are:

| Mode | Pack/species order |
| --- | --- |
| `hydrated-lime-neutralization-nbs1989` | Ca(OH)2(s), H+(aq), Ca2+(aq), H2O(l), all NBS |
| `calcium-hydroxide-precipitation` | Ca2+(aq), H2O(l), Ca(OH)2(s), H+(aq), all NBS |

Their final argument is respectively reacted pure Ca(OH)2 mass or reacted
calcium-ion mass, not a bulk-solution mass. `mineral_processing --help`
prints the exact pack names. These endpoint diagnostics and the simultaneous
finite-phase example serve different purposes.

## Scope

This is a fixed-temperature, fixed-solvent, ideal-activity calculation with
one explicitly selected crystalline solid. Calcium ion pairs, carbonate and
sulfate precipitation, competing solids and nonideal activity coefficients are
not included. The example does not model metal-cyanide competition or a gold
leach trajectory. Commercial purity, wetting, slaking/dissolution kinetics,
mixing heat and heat transfer require separate evidence. The finite equilibrium
does not supply a plant dose, operating pH recommendation or calibrated
recovery law.
