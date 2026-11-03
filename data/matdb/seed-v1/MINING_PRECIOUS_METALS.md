# Silver cyanidation and shared-zinc precious-metal recovery

Four NBS source profiles complete the references for native-silver cyanidation
and productive gold/silver zinc cementation. `MineralReaction` supplies the
conserved endpoints, and `precious_metal_recovery` consumes compiled material
packs through the existing finite-inventory `ReactionNetwork` owner.

## Material references

| Pack | Species | Original NBS location |
| --- | --- | --- |
| `silver-solid-nbs1989-thermochemistry` | Ag(cr) | Table 37:Ag, p. 2-160 |
| `dicyanoargentate-ion-aqueous-nbs1989` | Ag(CN)2[-](aq) | Table 37:Ag, p. 2-163 |
| `zinc-solid-nbs1989-thermochemistry` | Zn(cr) | Table 33:Zn, p. 2-138 |
| `tetracyanozincate-ion-aqueous-nbs1989` | Zn(CN)4[2-](aq) | Table 33:Zn, p. 2-143 |

Each profile supplies molar mass, formation enthalpy, independently tabulated
formation Gibbs energy, and standard entropy at **298.15 K and 100000 Pa**.
The [original NBS tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf) supply the
literal numerical rows. Historical molar masses are retained. Ionic entropy
is conventional, referenced to H+(aq) = 0; uncertainty remains unstated.

These profiles use the existing `source-elemental-reference-nbs1989` NBS
reference-family selector. All four imported properties and the source phases
were independently checked against [NIST's corrected CSV with the 1989 errata](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv)
and match the cited original table rows, including printed numerical precision.
Each pack records both sources. The selector does not imply that the underlying
rows were re-evaluated in 1989. No JANAF, USGS or USBM evaluation is relabelled.

## Conserved equations and finite stocks

```text
4 Ag(s) + 8 CN[-](aq) + O2(aq) + 2 H2O(l)
    -> 4 Ag(CN)2[-](aq) + 4 OH[-](aq)
2 Au(CN)2[-](aq) + Zn(s) -> 2 Au(s) + Zn(CN)4[2-](aq)
2 Ag(CN)2[-](aq) + Zn(s) -> 2 Ag(s) + Zn(CN)4[2-](aq)
```

The library variants are `SilverCyanidation` (gas oxygen),
`SilverCyanidationDissolvedOxygen`, `GoldCyanideCementationWithZinc`, and
`SilverCyanideCementationWithZinc`. The silver cementation balance follows
section 13.3.8 of [Aya's Zgounder technical report](https://www.sec.gov/Archives/edgar/data/1826836/000162828026024934/exhibit9986-technicalrepor.htm).
The [EPA process description](https://semspub.epa.gov/work/HQ/189773.pdf)
identifies zinc precipitation for both metals.

Both cementation reactions consume one shared zinc inventory. Their total
gross demand must fit the initial stock before any material or heat result is
returned. Extents are explicit caller decisions; no recovery priority or
kinetic selectivity is inferred. Leach products fund a subsequent explicit
cementation stage, not another reaction within the same simultaneous stage.
Bound cyanide remains in the zinc complex; it is not returned as free CN-.

## Run with the actual source packs

```bash
precious_pack_dir=/tmp/frankensim-precious-metal-packs
mkdir -p "$precious_pack_dir"
for pack in silver-solid-nbs1989-thermochemistry \
  dicyanoargentate-ion-aqueous-nbs1989 \
  zinc-solid-nbs1989-thermochemistry \
  tetracyanozincate-ion-aqueous-nbs1989 \
  gold-solid-nbs1989-thermochemistry dicyanoaurate-ion-aqueous-nbs1989 \
  cyanide-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 \
  oxygen-aqueous-nbs1989-thermochemistry water-liquid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$precious_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example precious_metal_recovery -- \
  silver-leach "$precious_pack_dir" 4 8 1 2 1
cargo run -p fs-thermochem --example precious_metal_recovery -- \
  cementation "$precious_pack_dir" 10 20 8 3 5
```

All amounts are moles; an extent is moles of the entire integer equation.
`silver-leach` takes initial silver, free cyanide, dissolved oxygen and water,
then the requested extent. This case produces 4 mol silver complex and 4 mol
hydroxide and exhausts the four reactant stocks. Its standard chemical
enthalpy change is -460216 J.

`cementation` takes initial gold complex, silver complex and zinc, followed
by separate gold/silver extents. The shown case consumes 8 mol zinc, produces
6 mol Au and 10 mol Ag, leaves 4 mol gold complex and 10 mol silver complex,
and forms 8 mol zinc complex. Its standard chemical enthalpy change is
-1418400 J. Changing the final extent from 5 to 6 requires 9 mol zinc and
refuses the whole stage without publishing a partial material or heat result.

Output includes initial, consumed, formed and final amounts, final species
masses, source pack identities/hashes, standard enthalpy, and both direct
Gibbs and enthalpy-minus-temperature-times-entropy diagnostics. The original
source disagreement between the two Gibbs routes is retained.

These are productive chemistry endpoints. Zinc side consumption, solution
activities, clarification, deaeration, filtration, smelting, heat of mixing,
and operating-temperature duty require separate inputs/models. Native Ag(s)
does not represent silver sulfides or an electrum activity model. Favorable
standard Gibbs energy does not determine a recovery target or rate.

`xtask/tests/mining_precious_metal_sources.rs` compiles the actual source packs
and checks all 16 new claims, all four reaction endpoints, source-basis mass
closure, shared zinc, explicit leach/cementation stages, receipt replay and
refusal of unsupported phase, reference and temperature queries.
