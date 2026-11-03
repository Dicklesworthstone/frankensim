# Gold cyanidation from a complete NBS source set

Six profiles supply all participants in the existing
`MineralReaction::GoldCyanidation` Elsner equation:

```text
4 Au(s) + 8 CN[-](aq) + O2(g) + 2 H2O(l)
    -> 4 Au(CN)2[-](aq) + 4 OH[-](aq)
```

The existing exact atom/charge checker, finite-inventory owner, standard
enthalpy/entropy and direct formation-Gibbs consumers use these profiles
without a new reaction or material schema. They add **24 claims**, four per
profile: molar mass, formation enthalpy, formation Gibbs energy and entropy.

## Primary numerical source and conventions

Values come from [NIST's NBS tables with 1989 errata](https://doi.org/10.18434/M32124),
specifically the [official corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
The rows were also checked against the
[original NBS tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf):
printed pages 2-37, 2-38, 2-89 and 2-164. NIST describes the digitization and
errata handling in [Reed (2020)](https://doi.org/10.6028/jres.125.007).

| Profile | Source formula/state | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | S, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `gold-solid-nbs1989-thermochemistry` | Au, cr | 196.9670 | 0 | 0 | 47.40 |
| `cyanide-ion-aqueous-nbs1989` | CN-, ao | 26.0179 | 150.6 | 172.4 | 94.1 |
| `oxygen-gas-nbs1989-thermochemistry` | O2, g | 31.9988 | 0 | 0 | 205.138 |
| `water-liquid-nbs1989-thermochemistry` | H2O, l | 18.0154 | -285.830 | -237.129 | 69.91 |
| `dicyanoaurate-ion-aqueous-nbs1989` | Au(CN)2-, ao | 249.0028 | 242.3 | 285.8 | 172. |
| `hydroxide-ion-aqueous-nbs1989` | OH-, ao | 17.0074 | -229.994 | -157.244 | -10.75 |

These are the source's printed units; the normal compiler converts mass to
kg/mol and energy to J/mol. All claims require **298.15 K and 100000 Pa** and
the distinct `source-elemental-reference-nbs1989=1` selector. The name means
NIST's NBS/errata evaluation, not that all entries were measured or revised in
1989: the underlying gold table was prepared in 1967.

The source's sections 5.5-5.8, printed pages 2-20 to 2-22, define the selected
states. Gold is pure crystalline Au; oxygen is molecular O2 in its pure
ideal-gas reference; water is pure liquid. The aqueous `ao` rows are individual
ions as written, in a hypothetical ideal solution at unit molality with the
activity coefficient referred to infinite dilution. Their H/G/S convention
sets the corresponding H+(ao) quantities to zero. Their entropies therefore use
`conventional-standard-molar-entropy`, including the negative OH- value.
The three neutral species use `standard-molar-entropy`. Charge conservation
allows the common ionic reference to cancel in the reaction sum.

These ionic profiles do not represent total cyanide, total dissolved gold,
HCN, solid NaCN, sodium/potassium dicyanoaurate salts or commercial solutions.
Aqueous queries require the existing individual-ion, unit-activity/infinite-
dilution and H+-zero selectors. Pure phases require their phase and pure-species
selectors. Every claim also retains the engineering-reference/specimen
selectors and `Unstated` uncertainty; no per-row probability model or
cross-species covariance is invented.

## Two retained thermodynamic routes

For one mole of the equation as written, which consumes four moles of gold:

| Quantity | Result |
| --- | ---: |
| Standard reaction enthalpy | -583916 J/mol extent |
| Standard reaction entropy | -642.358 J/(mol extent K) |
| Direct formation-Gibbs sum | -390718 J/mol extent |
| Caloric route, H - T S | -392396.9623 J/mol extent |
| Caloric minus direct route | -1678.9623 J/mol extent |

The source columns remain independent observations. The discrepancy is exposed,
not repaired by changing a number or attributed solely to rounding.
`--gibbs` uses the direct formation-Gibbs sum; `--equilibrium` uses H and S.
Both calculate their own dimensionless standard equilibrium constant through
ln K = -G/(R T). Activity calculations use the direct-G route.

The later [US Bureau of Mines IC 9429](https://stacks.cdc.gov/view/cdc/206360/cdc_206360_DS1.pdf),
Table D-9, printed page 154 (PDF page 172), gives **+297.6 kJ/mol** for
Au(CN)2- formation Gibbs energy. Its column is headed *minus* formation Gibbs
energy and prints -297.6. That later evaluation is distinct from the NBS
+285.8 kJ/mol row. This source set preserves the complete historical NBS
evaluation; it does not combine the later G value with NBS H/S or advertise
a reconciled contemporary gold-complex dataset.

## Run the source-backed consumer

The CLI requires packs in Au, CN-, O2, H2O, Au(CN)2-, OH- order:

```bash
gold_cyanidation_pack_dir=/tmp/mining-gold-cyanidation
mkdir -p "$gold_cyanidation_pack_dir"
gold_cyanidation_packs=(
  gold-solid-nbs1989-thermochemistry
  cyanide-ion-aqueous-nbs1989
  oxygen-gas-nbs1989-thermochemistry
  water-liquid-nbs1989-thermochemistry
  dicyanoaurate-ion-aqueous-nbs1989
  hydroxide-ion-aqueous-nbs1989
)
gold_cyanidation_args=()
for pack in "${gold_cyanidation_packs[@]}"; do
  pack_file="$gold_cyanidation_pack_dir/$pack.fsmatpk"
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" --out "$pack_file"
  gold_cyanidation_args+=("$pack_file")
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction gold-cyanidation --gibbs --equilibrium \
  "${gold_cyanidation_args[@]}" 0.787868
cargo test -p xtask --test mining_gold_cyanidation_sources
```

The final argument is **actually reacted pure-gold mass**, not ore mass.
Here 0.787868 kg is 4 mol Au and one mole of reaction extent. The reference
balance consumes 0.2081432 kg CN-, 0.0319988 kg O2 and 0.0360308 kg water,
forms 0.9960112 kg Au(CN)2- and 0.0680296 kg OH-, and returns -583916 J
standard enthalpy change.

To evaluate a local driving force, add one `--activity 'SPECIES=VALUE'` for
each of the six species, before the pack arguments. These are dimensionless
activities under the retained reference conventions, not molar concentrations.
With all other activities one, lowering CN- activity from one to 0.1 raises
the direct reaction Gibbs energy by 8 R T ln(10); lowering OH- activity to
0.1 lowers it by 4 R T ln(10). The activity result is dG/dextent at that
supplied state, not the total free-energy change over an arbitrary finite batch.

## Finite material availability

The integration tests exercise the real packs and existing inventory owner.
Starting with 4 mol Au, 4 mol CN-, 1 mol O2 and 2 mol water, absent products
represented explicitly by zero, cyanide permits at most 0.5 mol extent.
The endpoint contains 2 mol Au, no CN-, 0.5 mol O2, 1 mol water and 2 mol each
of Au(CN)2- and OH-. The corresponding standard enthalpy change is -291958 J.

With 4 mol Au, 8 mol CN-, 0.25 mol O2 and 2 mol water, oxygen instead limits
extent to 0.25 mol. Only 1 mol gold can dissolve, producing 1 mol each of the
complex and hydroxide. Requests above either material ceiling refuse.
The G0/G1/G3 tests also verify the 24 literal source values, uncertainty and
receipts, both thermodynamic routes, activity response, mass-to-extent scaling,
and refusal of temperature/pressure/reference/ionic-entropy substitutions.

These are standard-reference thermodynamic and material balances. They do not
predict a leach rate, mineral accessibility, HCN/CN- speciation, dissolved oxygen
supply, copper-cyanide competition, carbon loading or plant recovery. Reaction
extent and activities still require the corresponding physical models or
supported inputs.
