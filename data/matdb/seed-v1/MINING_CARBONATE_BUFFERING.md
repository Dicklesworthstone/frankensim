# Matched carbonate buffering of aqueous cyanide

Three NBS profiles complete a two-step carbonate family that the existing
`IdealAcidBaseSystem` can solve together with HCN/CN- and water. This makes
dissolved carbonate alkalinity affect the free-cyanide fraction through a
shared pH and explicit counterion inventories.

The profiles add **12 claims** and preserve the same NBS/1989 evaluation as
the existing proton, hydroxide, water and cyanide packs. The existing USGS
carbonate profiles retain their separate evaluation; they are not relabelled
to satisfy the coupled owner's requirement for identical shared sources.

## Literal primary-source values

| Profile | Source formula/state | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | S, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `carbonic-component-aqueous-nbs1989` | H2CO3, ao, formal component | 62.0254 | -699.65 | -623.08 | 187.4 |
| `bicarbonate-ion-aqueous-nbs1989` | HCO3-, ao | 61.0174 | -691.99 | -586.77 | 91.2 |
| `carbonate-ion-aqueous-nbs1989` | CO3(2-), ao | 60.0094 | -677.14 | -527.81 | -56.9 |

These are literal entries in the
[NIST corrected NBS dataset](https://doi.org/10.18434/M32124), available as an
[official CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
They were also checked visually against the
[original NBS tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf), Table 23:C,
printed page 2-83 for the two ions and 2-84 for the carbonic component.
NIST documents the digitization in [Reed (2020)](https://doi.org/10.6028/jres.125.007).

All claims require **298.15 K, 100000 Pa**, and
`source-elemental-reference-nbs1989=1`. The underlying carbon table was
prepared in 1966; the pack name identifies NIST's evaluation with 1989 errata,
not a new measurement date. Uncertainty remains `Unstated`; no probability
distribution, cross-species covariance or propagated confidence is invented.

### The formal carbonic component is essential

The source explicitly defines its H2CO3 entry as **equivalent to CO2(ao) plus
H2O(l)**. The 62.0254 g/mol formula weight includes that water formula unit.
It does not describe independently resolved hydrated molecular H2CO3.

Every claim in that profile therefore also requires
`source-carbonic-component-co2-plus-water=1`. Its `H2CO3(aq)` identity in
the reaction and family means this formal source component. A generic neutral
aqueous query without this qualifier refuses. No CO2 hydration constant or
fraction of molecular carbonic acid follows from this profile; adding an
independent hydration step to the same component would double-count its water.

The neutral component uses `standard-molar-entropy`. Bicarbonate and carbonate
use `conventional-standard-molar-entropy` on the H+-zero reference, including
the negative carbonate entropy. All aqueous profiles use the hypothetical
unit-molality reference with activity referred to infinite dilution.
Neither molarity nor a measured 1 mol/kg process solution is inferred.

## Reactions and numerical constants

The new conserved variants are:

```text
CarbonicComponentDissociation:
    H2CO3(aq, formal CO2+H2O component) -> H+(aq) + HCO3[-](aq)

BicarbonateDissociation:
    HCO3[-](aq) -> H+(aq) + CO3[2-](aq)
```

| Step | Standard H, J/mol extent | Standard S, J/(mol extent K) | Direct G, J/mol extent | H - T S, J/mol extent | Direct-G pK |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carbonic component to bicarbonate | 7660 | -96.2 | 36310 | 36342.03 | 6.361236782081848 |
| Bicarbonate to carbonate | 14850 | -148.1 | 58960 | 59006.015 | 10.329345102493686 |

The dimensionless constants follow ln K = -G/(RT). Direct-G and caloric
values remain independent source routes; their discrepancies are not removed
by changing a table entry. The coupled acid/base owner uses the direct-G route.

## Compile and run the original consumers

Compile these eight profiles into a new directory:

```bash
carbonate_pack_dir=$(mktemp -d)
for pack in carbonic-component-aqueous-nbs1989 bicarbonate-ion-aqueous-nbs1989 carbonate-ion-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry hydrogen-cyanide-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$carbonate_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction carbonic-component-dissociation --gibbs --equilibrium \
  "$carbonate_pack_dir/carbonic-component-aqueous-nbs1989.fsmatpk" \
  "$carbonate_pack_dir/hydrogen-ion-aqueous-nbs1989.fsmatpk" \
  "$carbonate_pack_dir/bicarbonate-ion-aqueous-nbs1989.fsmatpk" 0.0620254

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction bicarbonate-dissociation --gibbs --equilibrium \
  "$carbonate_pack_dir/bicarbonate-ion-aqueous-nbs1989.fsmatpk" \
  "$carbonate_pack_dir/hydrogen-ion-aqueous-nbs1989.fsmatpk" \
  "$carbonate_pack_dir/carbonate-ion-aqueous-nbs1989.fsmatpk" 0.0610174

cargo test -p xtask --test mining_carbonate_speciation_sources
```

The two mass arguments declare one mole of the first formal component or
bicarbonate actually reacted. The CLI reports standard reaction quantities;
its `--equilibrium` flag does not solve a composition.

For simultaneous composition, build an `AcidFamily` named, for example,
`inorganic-carbon` from the two source-bound `AcidDissociation` steps, in
the order shown above. Join it with the existing `uncomplexed-cyanide` family
and `WaterAutoprotolysis` in `IdealAcidBaseSystem::new`. Supply both analytical
totals and explicit spectator ions to `speciate`.

The carbon total is the sum of formal H2CO3, HCO3- and CO3(2-) amounts,
each containing one carbon atom. The separate HCN+CN- pool excludes metal
complexes and WAD cyanide. The solver preserves both pools while sharing
proton and hydroxide activities. Its fixed solvent basis does not assert a
closed solvent H/O balance.

## Demonstrated buffering response

The actual-source integration test compares two declared ideal solutions on
the same **10 kg solvent basis**, giving n_ref = 10 mol at unit molal activity.
Both start at pH 10 with 0.01 mol of uncomplexed cyanide; each has the explicit
Na+ amount required by its own initial charge balance. One also contains
0.1 mol of inorganic carbon. The same 0.005 mol strong-acid addition is
represented by added Cl- spectator inventory.

| Inorganic carbon | Initial free CN-, mol | Final pH | Final free CN-, mol |
| --- | ---: | ---: | ---: |
| 0 mol | 0.00854073031671 | 9.13057963584 | 0.00441517915449 |
| 0.1 mol | 0.00854073031671 | 9.90602505634 | 0.00824989563011 |

Carbonate reduces the pH change and preserves more free cyanide under these
declared ideal conditions. These are source-backed model results with synthetic
inventories, not measured plant recoveries or a recommended reagent dose.

The three focused integration tests compile the real packs and check all
12 new values, source/component admission, both thermodynamic routes,
two-family speciation against independently evaluated polynomial fractions,
charge and analytical-pool closure, the acid-pulse response, and invariance
when all stocks and solvent basis scale together.

This increment does not resolve true CO2 hydration, gas exchange, calcium
complexes or carbonate precipitation. Lime dissolution, finite oxygen
replenishment, activity corrections and leach kinetics require their own
source-bound owners; an inorganic-carbon inventory does not imply a lime dose.
