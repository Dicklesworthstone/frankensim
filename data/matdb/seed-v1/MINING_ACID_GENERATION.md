# Aqueous pyrite oxidation and acid-generation references

Compatible profiles from [USGS Bulletin 2131](https://doi.org/10.3133/b2131)
support pyrite oxidation, oxidation of dissolved ferrous iron, and attack on
pyrite by ferric ions. Complete formation-enthalpy and formation-Gibbs source
sets allow reagent accounting, standard reaction heat, and thermodynamic
driving-force calculations at a declared activity state. They use the normal
`matdb-pack` compiler and typed `NormalizedPack` queries. The selected rows were
checked in the original report, including its
[government-document mirror](https://digital.library.unt.edu/ark:/67531/metadc1114954/).

## Compatible source values

All quantities below apply at **298.15 K and 100000 Pa** with
`source-elemental-reference-usgs-b2131=1`. H is standard molar formation
enthalpy; G is standard molar formation Gibbs energy; S is standard molar
entropy for pure neutral species and conventional ionic entropy for individual
ions. The Gibbs property is `standard-molar-formation-gibbs-energy`, with
dimensions of energy per amount.

| Profile | Source page | Molar mass, kg/mol | H, J/mol | G, J/mol | S, J/(mol K) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `pyrite-solid-usgs-b2131-thermochemistry` | [12](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/18/) | 0.119979 | -171500 | -160100 | 52.9 |
| `oxygen-gas-usgs-b2131` | [8](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/14/) | 0.031999 | 0 | 0 | 205.15 |
| `water-liquid-usgs-b2131` | [16](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/22/) | 0.018015 | -285800 | -237100 | 70.0 |
| `iron-ii-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.055847 | -91100 | -90000 | -107.1 |
| `iron-iii-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.055847 | -49900 | -16700 | -280.0 |
| `sulfate-ion-aqueous-usgs-b2131` | [18](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/24/) | 0.096064 | -909300 | -744000 | 18.5 |
| `hydrogen-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.001008 | 0 | 0 | 0 |

These profiles retain the USGS evaluation and rounded formula
weights. Existing JANAF pyrite, oxygen and water profiles keep their separate
values and reference identity. In particular, USGS pyrite is 119.979 g/mol and
-171.5 kJ/mol, while the existing JANAF profile is 119.967 g/mol and
-171.544 kJ/mol. Neither is silently relabeled to combine reference sets.
The selected pyrite row uses the standard elemental sulfur reference;
the report's alternative S2-gas-reference tables are excluded.

Every profile requires typed temperature and pressure, `source-standard-state=1`,
`engineering-reference-model=1`, and `source-specimen-qualified=0`. Neutral
species require their gas/liquid/solid phase and `source-pure-species=1`;
pyrite additionally requires `source-phase-pyrite=1`.

The ions require `source-phase-aqueous=1`, `source-individual-ion=1`,
`source-aqueous-unit-activity-infinite-dilution=1` and
`source-ionic-reference-hydrogen-zero=1`. They have no pure-species assertion.
USGS's aqueous `m=1` denotes a hypothetical ideal state at unit activity and
infinite dilution, not a measured finite-concentration 1 mol/kg solution.
The entropy property is `conventional-standard-molar-entropy`, preserving the
negative ferrous- and ferric-ion values; it cannot be queried as absolute ionic
entropy. Ferric ions are individual Fe3+(aq), distinct from total dissolved iron,
ferric-sulfate complexes and iron-hydroxide solids.

The original rows print these uncertainty magnitudes:

| Source | H magnitude, J/mol | G magnitude, J/mol | S magnitude, J/(mol K) |
| --- | ---: | ---: | ---: |
| Pyrite | 1700 | 1700 | 0.1 |
| Oxygen | Not printed for the reference zero | Not printed for the reference zero | 0.02 |
| Liquid water | 100 | 100 | 0.1 |
| Ferrous ion | 3000 | 2000 | 2.0 |
| Ferric ion | 5000 | 2000 | 13.0 |
| Sulfate ion | 400 | 400 | 0.4 |

The H and S magnitudes are exposed as dimensioned metadata for pyrite, water,
ferrous ion, ferric ion and sulfate; oxygen has an S metadata claim. Gibbs
magnitudes are preserved in their source observations. The ferric profile
retains its separately tabulated formation-Gibbs claim alongside the H and S
metadata. The source's calorimetric convention is twice the standard deviation
of the mean.
It supplies no universal coverage probability or cross-species covariance.
All claims retain `Unstated` uncertainty semantics; the metadata are not
independent Gaussian standard deviations or certified reaction-error bounds.

## Ferrous endpoint, reagent inventory and heat

`MineralReaction::PyriteOxidationFerrous` declares:

```text
2 FeS2(pyrite) + 7 O2(g) + 2 H2O(l)
    -> 2 Fe2+(aq) + 4 SO4(2-)(aq) + 4 H+(aq)
```

The existing exact atom/charge checker admits that integer equation. Its
six compatible formation sources give **-2904800 J per mole of the written
reaction**, equivalent to -1452400 J per mole of pyrite consumed. Oxygen and
hydrogen ions remain required sources even though their formation enthalpies
are zero by convention.

At **119.979 kg of pure pyrite actually reacted**, the source mass represents
1000 mol pyrite and 500 mol equation extent. The balance consumes 3500 mol O2
and 1000 mol water, produces 1000 mol Fe2+, 2000 mol sulfate and 2000 mol H+,
and gives **-1452.4 MJ** standard enthalpy change. The independently rounded
source masses retain a +0.0005 kg products-minus-reactants residual. This does
not alter the exact atom/charge balance or justify adjusting source values.

For a finite inventory of 10 mol pyrite, 14 mol oxygen and 20 mol water,
oxygen limits the possible extent to 2 mol of the written equation. Advancing
that declared extent consumes 4 mol pyrite and generates 8 mol H+; any prior
proton inventory remains accounted for separately. This is a material-availability
ceiling. It does not predict the extent reached in a particular residence time.

## Ferrous oxidation and ferric leaching

`MineralReaction::FerrousOxidationFerric` declares:

```text
4 Fe2+(aq) + O2(g) + 4 H+(aq) -> 4 Fe3+(aq) + 2 H2O(l)
```

`MineralReaction::FerricPyriteOxidation` declares:

```text
FeS2(pyrite) + 14 Fe3+(aq) + 8 H2O(l)
    -> 15 Fe2+(aq) + 2 SO4(2-)(aq) + 16 H+(aq)
```

These stoichiometries follow equations 3–4 of
[USGS Open-File Report 99-348](https://pubs.usgs.gov/of/1999/0348/report.pdf),
with the ferrous equation multiplied by four to retain integer coefficients.
They also appear as equations [2] and [4] on printed page 2-3 of the
EPA-hosted [Handbook for Constructed Wetlands Receiving Acid Mine Drainage](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=94001BR7.TXT)
(December 1991 draft). The
implementation explicitly retains O2(g), liquid water and dissolved individual
ions as the source phases. The ferric-leach equation consumes ferric inventory;
it has no direct oxygen term. Ferrous oxidation can replenish that ferric
inventory only while its own oxygen and proton reactants remain available.

The compatible source sums give:

| Integer reaction | Standard enthalpy, J/mol extent | Standard Gibbs energy, J/mol extent | ln K° at 298.15 K |
| --- | ---: | ---: | ---: |
| Pyrite oxidation to dissolved ferrous iron | -2904800 | -2361600 | 952.6587 |
| Ferrous oxidation | -406800 | -181000 | 73.0146 |
| Ferric leaching of pyrite | -28600 | -547300 | 220.7783 |

Two ferric-leach extents plus seven ferrous-oxidation extents recover one
integer pyrite-ferrous equation after cancelling the intermediate iron ions.
Both the H and G sums above reproduce that net balance. An executable sequence
still needs each stage's reactants to be present when the stage advances;
cancellation on paper supplies no initial ferric inventory or reaction rate.

The finite-inventory regression starts with 28 mol Fe3+, attacks two mol
pyrite, then oxidizes 28 mol of the resulting Fe2+ in a separate explicit
stage. It restores the original ferric stock and retains the net ferrous,
sulfate, proton, water and oxygen changes. Asking the same initial inventory
to fund both gross reaction demands simultaneously refuses; products from a
future stage cannot supply a present reagent deficit. Extents remain explicit
inputs, with no inferred microbial rate, dissolved-oxygen transfer or pH.

## Standard Gibbs energy and supplied activities

`ConservedReaction::standard_gibbs_change` consumes a complete set of
`FormationGibbsSource` values. Each selected claim must carry its compatible
phase, elemental reference and, for ions, aqueous activity and ionic-reference
conventions. Oxygen and H+ terms in an equation still require their source
claims even though their formation Gibbs energies are zero by convention.
The result retains the original material answers and receipts.

Using the IUPAC definitions of
[standard equilibrium constant](https://goldbook.iupac.org/terms/view/S05915)
and [activity](https://goldbook.iupac.org/terms/view/A00115), the consumer evaluates
the following at the common source-admitted temperature, with positive
stoichiometric coefficients for products and negative ones for reactants:

$$
\Delta_r G^\circ = \sum_i \nu_i\Delta_f G_i^\circ,
\qquad \ln K^\circ = -\frac{\Delta_r G^\circ}{RT}.
$$

`StandardReactionGibbs::at_activities` then evaluates:

$$
\ln Q = \sum_i \nu_i\ln a_i,
\qquad \Delta_r G = RT(\ln Q-\ln K^\circ).
$$

Every active species needs one finite, strictly positive dimensionless activity,
including water and pyrite. Pure condensed phases may be assigned activity one
when that assumption applies; the caller declares it explicitly. Gas activity
is a fugacity ratio relative to the retained standard state. Identifying it
with a pressure ratio requires an additional ideal-gas assumption. Ion
activities must use the source's aqueous convention; raw molar concentrations
are not accepted as an activity model.

Negative `ΔrG` identifies a forward thermodynamic driving force at the supplied
state; positive values identify a reverse driving force. It is the local
derivative `dG/dξ`, so multiplying it by a finite extent does not integrate a
process whose activities change. The API does not solve for equilibrium
composition, pH, activity coefficients or kinetics. Missing, duplicate,
unknown-species, zero, negative and non-finite activities refuse. These point
claims admit 298.15 K and 100000 Pa only; Cp data elsewhere do not extend their
Gibbs domain.

For pyrite oxidation, direct `exp(ln K°)` overflows a finite `f64`; the result
keeps the finite `ln K°` and reports `LogOnlyOverflow`. The log-domain activity
calculation remains usable. At unit activities `ln Q=0` and `ΔrG=ΔrG°`.

The printed entropy rows also give reaction entropy -1822.05 J/(mol K).
At 298.15 K, reaction H-T*S is -2361555.7925 J/mol, while the independently
tabulated formation-Gibbs sum is -2361600 J/mol. Their 44.2075 J/mol rounding
difference is retained. Subtracting temperature times an individual species'
entropy from its formation enthalpy is not a formation-Gibbs calculation.

## Run the source-backed balance

```bash
acid_pack_dir="$(mktemp -d)"
for pack in pyrite-solid-usgs-b2131-thermochemistry oxygen-gas-usgs-b2131 water-liquid-usgs-b2131 iron-ii-ion-aqueous-usgs-b2131 iron-iii-ion-aqueous-usgs-b2131 sulfate-ion-aqueous-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$acid_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction pyrite-oxidation-ferrous --gibbs \
  --activity 'FeS2(s)=1' --activity 'O2(g)=0.21' --activity 'H2O(l)=1' \
  --activity 'Fe[2+](aq)=0.001' --activity 'SO4[2-](aq)=0.001' \
  --activity 'H+(aq)=0.0001' \
  "$acid_pack_dir/pyrite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$acid_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/sulfate-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 119.979
cargo test -p xtask --test mining_acid_generation --test mining_iron_redox --test mining_aqueous_redox_sources
```

All options precede the pack list; species identifiers are quoted literally.
The supplied numbers are illustrative dimensionless activities, with no
concentration-to-activity or oxygen-pressure conversion inferred. This state
gives `ln Q` about -67.36336 and `ΔrG` about -2528590.87 J/mol extent. The final
119.979 is kilograms of pure pyrite actually reacted for the separate
stoichiometric and enthalpy calculation. Omit the `--activity` options to
request only the standard Gibbs result with `--gibbs`.

The focused tests compile actual packs, verify literal values and receipts,
check phase/reference admission, exercise finite reagent inventories and
reject a caller's attempt to relabel JANAF pyrite as USGS. The iron-redox
regressions verify the ferric source, exact cycle cancellation, source-backed
Hess heat and finite staged inventories. The Gibbs regressions verify direct
formation-Gibbs sums and explicit activity calculations.

The same CLI supports each redox stage, including optional standard-equilibrium
output from the original enthalpy and entropy claims:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --reaction ferrous-oxidation-ferric \
  "$acid_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/water-liquid-usgs-b2131.fsmatpk" 55.847
cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --reaction ferric-pyrite-oxidation \
  "$acid_pack_dir/pyrite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$acid_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/sulfate-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 119.979
```

The first basis is 55.847 kg chemical-equivalent mass of reacted dissolved
Fe2+, not solution mass or an isolated pure-ion phase. It represents 1000 mol
Fe2+ and 250 mol written extent, with -101.7 MJ standard enthalpy change. The
second consumes 1000 mol pyrite at 1000 mol written extent, with -28.6 MJ heat;
ferric regeneration is a separate explicit stage. Neither CLI invocation
implicitly supplies a finite inventory or predicts conversion.

Add `--gibbs` before the pack list to report the separately tabulated
formation-Gibbs route alongside `--equilibrium`. The CLI labels both routes
and retains their source-rounding differences. Supplying `--activity` requires
the direct formation-Gibbs claims even when `--equilibrium` is also requested.

The selected reactions retain the declared dissolved iron endpoint. Hydroxide
or sulfate precipitation, dissolved-oxygen transfer, a concentration-to-activity
model, microbial kinetics, process pH and heat of mixing require additional
models. Standard heat does not predict a pile or tank temperature, and the
ionic acid balance does not specify a commercial neutralizing-reagent dose.
The compatible calcite/bicarbonate neutralization source set is documented in
[Reactive mineral thermochemistry](MINERAL_THERMOCHEMISTRY.md#aqueous-references-for-calcite-neutralization).
