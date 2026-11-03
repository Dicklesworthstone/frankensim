# Aqueous pyrite oxidation and acid-generation references

These six profiles use [USGS Bulletin 2131](https://doi.org/10.3133/b2131)
and complete the ferrous-pyrite endpoint and ferrous/ferric redox cycle with the
existing hydrogen-ion profile. They use the normal `matdb-pack` compiler and typed
`NormalizedPack` queries. The selected rows were checked in the original report,
including its [government-document mirror](https://digital.library.unt.edu/ark:/67531/metadc1114954/).

The ferric profile also retains its separately tabulated formation-Gibbs claim.
Its formation-enthalpy and entropy uncertainty magnitudes are exposed as dimensional
metadata without replacing that existing Gibbs datum.

## Compatible source values

All quantities below apply at **298.15 K and 100000 Pa** with
`source-elemental-reference-usgs-b2131=1`. H is standard molar formation
enthalpy; S is standard molar entropy for pure neutral species and conventional
ionic entropy for individual ions.

| Profile | Source page | Molar mass, kg/mol | H, J/mol | S, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `pyrite-solid-usgs-b2131-thermochemistry` | [12](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/18/) | 0.119979 | -171500 | 52.9 |
| `oxygen-gas-usgs-b2131` | [8](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/14/) | 0.031999 | 0 | 205.15 |
| `water-liquid-usgs-b2131` | [16](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/22/) | 0.018015 | -285800 | 70.0 |
| `iron-ii-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.055847 | -91100 | -107.1 |
| `iron-iii-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.055847 | -49900 | -280.0 |
| `sulfate-ion-aqueous-usgs-b2131` | [18](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/24/) | 0.096064 | -909300 | 18.5 |
| Existing `hydrogen-ion-aqueous-usgs-b2131` | [7](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/13/) | 0.001008 | 0 | 0 |

The new profiles deliberately retain the USGS evaluation and rounded formula
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
negative ferrous-ion value; it cannot be queried as absolute ionic entropy.

Printed uncertainty magnitudes are retained as dimensioned metadata:

| Source | H magnitude, J/mol | S magnitude, J/(mol K) |
| --- | ---: | ---: |
| Pyrite | 1700 | 0.1 |
| Oxygen | Not printed for the reference zero | 0.02 |
| Liquid water | 100 | 0.1 |
| Ferrous ion | 3000 | 2.0 |
| Ferric ion | 5000 | 13.0 |
| Sulfate ion | 400 | 0.4 |

The source's calorimetric convention is twice the standard deviation of the
mean. It supplies no universal coverage probability or cross-species covariance.
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

The same printed entropy rows give a reaction entropy of -1822.05 J/(mol K).
At 298.15 K, H-T*S yields -2361555.7925 J/mol of reaction. Independently
summing the report's rounded formation-Gibbs entries gives -2361600 J/mol;
their 44.2075 J/mol difference remains visible rather than forcing either
source quantity to agree by changing its published rounding.

## Ferrous oxidation and ferric attack on pyrite

The two explicit endpoints follow equations 3–4 of
[USGS Open-File Report 99-348](https://pubs.usgs.gov/of/1999/0348/report.pdf),
with the ferrous equation multiplied by four to retain integer coefficients:

```text
FerrousOxidationFerric:
4 Fe2+(aq) + O2(g) + 4 H+(aq) -> 4 Fe3+(aq) + 2 H2O(l)

FerricPyriteOxidation:
FeS2(pyrite) + 14 Fe3+(aq) + 8 H2O(l)
    -> 15 Fe2+(aq) + 2 SO4(2-)(aq) + 16 H+(aq)
```

The compatible sources give -406800 J/mol of the written ferrous-oxidation
equation and -28600 J/mol of the ferric-pyrite equation. Seven oxidation
extents plus two ferric-attack extents cancel the ferric inventory and recover
the existing two-pyrite ferrous endpoint exactly, both in atom/charge counts
and in its -2904800 J/mol standard heat.

The finite-inventory regression starts with 28 mol Fe3+, attacks two mol
pyrite, then oxidizes 28 mol of the resulting Fe2+ in a separate explicit
stage. It restores the original ferric stock and retains the net ferrous,
sulfate, proton, water and oxygen changes. Asking the same initial inventory
to fund both gross reaction demands simultaneously refuses; products from a
future stage cannot supply a present reagent deficit. Extents remain explicit
inputs, with no inferred microbial rate, dissolved-oxygen transfer or pH.

## Run the source-backed balance

```bash
acid_pack_dir="$(mktemp -d)"
for pack in pyrite-solid-usgs-b2131-thermochemistry oxygen-gas-usgs-b2131 water-liquid-usgs-b2131 iron-ii-ion-aqueous-usgs-b2131 iron-iii-ion-aqueous-usgs-b2131 sulfate-ion-aqueous-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$acid_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction pyrite-oxidation-ferrous \
  "$acid_pack_dir/pyrite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$acid_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/sulfate-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 119.979
cargo test -p xtask --test mining_acid_generation
cargo test -p xtask --test mining_iron_redox
```

The original three focused tests compile the actual packs, verify literal values and
receipts, check phase/reference admission, exercise the oxygen-limited
inventory and reject a caller's attempt to relabel JANAF pyrite as USGS.
Three additional redox tests verify the ferric source, exact cycle cancellation,
source-backed Hess heat and finite staged inventories.

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

This endpoint retains ferrous iron in solution and uses gaseous oxygen as its
source reference. Hydroxide or sulfate precipitation,
dissolved-oxygen transfer, finite-concentration activities, microbial kinetics,
process pH and heat of mixing need separate models. The standard heat does not
predict a pile or tank temperature, and the ionic acid balance does not specify
a commercial neutralizing-reagent dose.
