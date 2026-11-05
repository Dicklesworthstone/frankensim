# Sulfite reagents and finite SO2/air cyanide-treatment accounts

Four material profiles add **15 literal source claims** from the
[NBS tables with the 1989 errata](https://doi.org/10.18434/M32124).
The original corrected CSV is pinned in each manifest by URL and SHA-256.
The source states are Table 14:S, HSO3- and SO3(2-) individual aqueous ions,
and Table 99:Na, crystalline anhydrous Na2SO3 and Na2S2O5.

| Pack prefix | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | S, J/(mol K) | Cp, J/(mol K) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `hydrogen-sulfite-ion-aqueous` | 81.0702 | -626.22 | -527.73 | 139.7, conventional | unavailable |
| `sulfite-ion-aqueous` | 80.0622 | -635.5 | -486.5 | -29., conventional | unavailable |
| `sodium-metabisulfite-solid` | 190.1046 | -1478.2 | unavailable | unavailable | unavailable |
| `sodium-sulfite-solid` | 126.0418 | -1100.8 | -1012.5 | 145.94 | 120.25 |

The aqueous prefixes end in `-nbs1989`; solid prefixes end in
`-nbs1989-thermochemistry`. Every claim requires 298.15 K, 100000 Pa,
the NBS elemental reference and its declared phase/source selectors.
Individual-ion entropy uses the H+(aq)-zero convention; negative sulfite
entropy is not changed or reinterpreted as an absolute entropy. All statistical
uncertainties remain unstated. The single solid Cp point is not a temperature law.
Metabisulfite's missing Gibbs energy/entropy remain query refusals, including
when its known enthalpy can be used in a reaction-heat calculation.

## Conserved endpoints

`fs_thermochem::stoichiometry::minerals::MineralReaction` provides seven
exact atom/charge-balanced equations. The productive endpoints are the
free-cyanide, free-sulfate forms of the SO2/air balances in
[SGS T3 SGS 018, May 2005, equations 13–14](https://www.sgs.com/pt/-/media/sgscorp/documents/corporate/brochures/sgs-min-wa017-cyanide-destruction-en-11.cdn.pt.pdf).
Choosing free CN- is explicit: WAD-CN is an analytical measurand, not one
species whose atoms, metal inventory and charge can be inferred.

| Variant | Declared equation |
| --- | --- |
| `SodiumMetabisulfiteDissolution` | Na2S2O5(s) + H2O(l) → 2 Na+ + 2 HSO3- |
| `SodiumSulfiteDissolution` | Na2SO3(s) → 2 Na+ + SO3(2-) |
| `HydrogenSulfiteDissociation` | HSO3- → H+ + SO3(2-) |
| `CyanideBisulfiteOxidation` | CN- + HSO3- + O2(aq) → OCN- + SO4(2-) + H+ |
| `CyanideSulfiteOxidation` | CN- + SO3(2-) + O2(aq) → OCN- + SO4(2-) |
| `SulfiteOxygenOxidation` | 2 SO3(2-) + O2(aq) → 2 SO4(2-) |
| `BisulfiteOxygenOxidation` | 2 HSO3- + O2(aq) → 2 SO4(2-) + 2 H+ |

All ionic species in this table are aqueous. O2(aq) is a finite dissolved
stock; an air flow cannot be substituted for it. Dissolving metabisulfite
retains its sodium and consumes solvent water. Productive oxidation and
unproductive reagent oxidation share the same admitted oxygen/reagent stocks.
They have different cyanide conversion and proton accounts.

Matched NBS references give 31350 J/mol for metabisulfite preparation,
-567950 J/mol for productive bisulfite oxidation and -554400 J/mol for the
written two-bisulfite side reaction. These are standard chemical enthalpy
changes, not heat duty at arbitrary process conditions. Separately, the
HSO3-/SO3(2-) family uses the original acid/base solver and independently
tabulated formation Gibbs energies; its dissociation delta G is 41230 J/mol.
That two-state family does not represent low-pH dissolved SO2 speciation.

## Run the actual source-pack consumer

Compile the nine participating packs into an explicit output directory:

```bash
mkdir -p /tmp/mining-sulfite-packs
for pack in \
  sodium-metabisulfite-solid-nbs1989-thermochemistry \
  hydrogen-sulfite-ion-aqueous-nbs1989 \
  sodium-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry \
  cyanide-ion-aqueous-nbs1989 \
  cyanate-ion-aqueous-nbs1989 \
  sulfate-ion-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  oxygen-aqueous-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "/tmp/mining-sulfite-packs/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example cyanide_treatment -- \
  sulfite /tmp/mining-sulfite-packs 10 4 100 7 4 6 1
```

These **authored accounting inputs** are 10 mol free cyanide, 4 mol crystalline
metabisulfite, 100 mol water and 7 mol dissolved oxygen. Declared dissolution
extent is 4 mol; productive and side-reaction extents are 6 and 1 mol of their
written equations. After preparation and oxidation, the result retains
4 mol CN-, 6 mol cyanate, 8 mol sulfate, 8 mol sodium and 8 mol proton
equivalents, with 96 mol water and exhausted dissolved oxygen. Nitrogen is
still present. Combined standard enthalpy change is -3836700 J.

Both stages must succeed before accepted output is printed. A request for
8 mol productive conversion has only 7 mol dissolved oxygen and refuses.
Preparation and oxidation cannot borrow newly formed bisulfite within one
simultaneous stage; the completed preparation stage must supply it.

`cargo test -p xtask --test mining_sulfite_sources` compiles the real manifests,
round-trips source packs and receipts, checks all 15 literals and seven reaction
heats, exercises finite staged inventories and shared oxygen refusal, and
checks one charge-balanced sulfite acid/base solution. The chemistry library's
original reference-reaction test also checks exact conservation and replay.
Neither those tests nor these inputs establish copper-catalyst kinetics,
gas transfer, commercial dose, complete dissolution, metal-complex conversion,
an automatic pH setpoint, WAD/total-CN removal, or final water quality.
