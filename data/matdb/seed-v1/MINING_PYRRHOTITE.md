# Pyrrhotite composition, oxidant demand and declared reaction endpoints

Three matched NBS material packs now drive four executable pyrrhotite
reactions through `fs-thermochem` and the ordinary `mineral_processing`
example. They distinguish the source's crystalline alpha FeS reference from
its sulfur-rich Fe₇S₈ reference. That distinction changes mineral mass,
dissolved-oxygen demand, dissolved-iron production and direct acid generation.

## Source phases and formula bases

The values come from the [NIST digitization with 1989 errata](https://doi.org/10.18434/M32124)
of Wagman et al.'s *NBS Tables of Chemical Thermodynamic Properties*,
[original tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf), Tables 41:Fe
and 14:S. The exact rows are available in the
[corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
All claims refer to **298.15 K and 100000 Pa**.

| Pack | Exact source phase/formula basis | M, g/mol | ΔfH°, kJ/mol | ΔfG°, kJ/mol | S°, J/(mol K) | Cp°, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `pyrrhotite-fes-alpha-solid-nbs1989` | Fe1.000S, iron-rich pyrrhotite, alpha, crystalline | 87.9110 | -100.0 | -100.4 | 60.29 | 50.54 |
| `pyrrhotite-fe7s8-solid-nbs1989` | Fe₇S₈, sulfur-rich pyrrhotite, crystalline | 647.4410 | -736.4 | -748.5 | 485.8 | 398.57 |
| `sulfur-orthorhombic-solid-nbs1989` | Rhombic crystalline S, per mole of S atoms | 32.0640 | 0 | 0 | 31.80 | 22.64 |

Each row supplies five ordinary claims with original source observations,
typed quantities and unstated statistical uncertainty. Fe₇S₈ quantities are
per formula unit containing **seven Fe and eight S atoms**. Dividing only its
molar mass by eight, or treating its formation energies as values per Fe0.875S,
would corrupt the reaction balance. Sulfur's zero formation energies are
explicit elemental-reference claims; their source is still required.

The pyrrhotite packs require their own phase flags and exact Fe:S ratios
(`source-iron-atoms-per-sulfur-atom=1` or `0.875`). Neither pack describes an
arbitrary Fe1-xS solid solution, mackinawite, pyrite, mineral mixtures or an ore
lot. No continuous composition interpolation is introduced. Point Cp values
do not extend the formation-enthalpy/Gibbs domain or provide thermal feedback.

## Oxygen and ferric pathways

The dissolved-oxygen balance follows equation 3 in
[USGS Open-File Report 02-195, chapter L](https://pubs.usgs.gov/of/2002/of02-195/OF02-195L.pdf).
Its oxygen coefficient is **2 - x/2** for Fe1-xS. The implementation admits
only the two selected source compositions and retains exact integer equations:

```text
FeS(alpha-pyrrhotite) + 2 O2(aq)
    -> Fe2+(aq) + SO4(2-)(aq)

2 Fe7S8(pyrrhotite) + 31 O2(aq) + 2 H2O(l)
    -> 14 Fe2+(aq) + 16 SO4(2-)(aq) + 4 H+(aq)
```

At the same eight-mole sulfur basis, eight moles of FeS consume 16 mol O₂ and
produce eight mol Fe²⁺ with no direct proton term. One mole of Fe₇S₈ consumes
15.5 mol O₂ and produces seven mol Fe²⁺ and two mol H⁺. Subsequent ferrous
oxidation, ferric hydrolysis and neutralization can change net acid demand;
those are separate reactions. The zero direct proton term for FeS does not
establish zero long-term acid generation from an ore containing it.

The ferric-ion pathways follow equation 3 in
[the same report, chapter B](https://pubs.usgs.gov/of/2002/of02-195/OF02-195B.pdf):

```text
FeS(alpha-pyrrhotite) + 2 Fe3+(aq) -> 3 Fe2+(aq) + S(s)
Fe7S8(pyrrhotite) + 14 Fe3+(aq) -> 21 Fe2+(aq) + 8 S(s)
```

Ferric ions are consumed inventory. The additional ferrous iron from the
mineral is accounted for separately from reduced ferric reagent. These
endpoints select orthorhombic elemental sulfur as the standard reference
product; they do not assert that an actual leach forms that polymorph or that
a sulfur film has a particular passivating effect. Sulfate formation, oxygen
uptake and proton production are not silently added to these ferric endpoints.

The four public `MineralReaction` variants are
`IronRichPyrrhotiteDissolvedOxygenOxidation`,
`SulfurRichPyrrhotiteDissolvedOxygenOxidation`,
`FerricIronRichPyrrhotiteLeach` and `FerricSulfurRichPyrrhotiteLeach`.
Each goes through the existing exact elemental and charge checker. A caller
selects an endpoint and a reacted extent; the records supply no rate,
residence-time conversion, reaction priority, gas transfer or recovery model.

## Source-driven reaction heat and Gibbs energy

The new packs combine with the existing **same-evaluation NBS** references
for dissolved O₂, liquid water, H⁺, sulfate and ferrous/ferric ions. They do
not substitute the separately tabulated USGS or JANAF values.

| Example reaction option | ΔrH°, J/mol written equation | ΔrG°, J/mol written equation |
| --- | ---: | ---: |
| `pyrrhotite-fes-oxidation-dissolved-oxygen` | -874970 | -755830 |
| `pyrrhotite-fe7s8-oxidation-dissolved-oxygen` | -13388560 | -11554222 |
| `ferric-pyrrhotite-fes-leach` | -70300 | -126900 |
| `ferric-pyrrhotite-fe7s8-leach` | -455700 | -842600 |

The Fe₇S₈ oxygen equation consumes **two** mineral formula units per equation
extent. Thus 647.441 kg of Fe₇S₈ actually reacted is 1000 mol mineral and
500 mol equation extent. The ferric Fe₇S₈ equation uses one mineral unit per
extent. Standard Gibbs energy is summed from the independent formation-Gibbs
column; it is not calculated as a species' formation enthalpy minus T times
its absolute entropy. Optional supplied activities remain explicit through
the existing `--activity` interface.

The shared-inventory integration test starts with three mol FeS, two mol
Fe₇S₈, 20 mol dissolved O₂ and 100 mol water. Requested extents 2 and 0.5
consume two mol FeS and one mol Fe₇S₈, leaving 0.5 mol O₂ and producing nine
mol Fe²⁺, ten mol sulfate and two mol H⁺. The matched standard enthalpy
increment is -8.44422 MJ. Increasing the FeS extent to three would demand
21.5 mol O₂ and refuses atomically. No atmospheric reservoir is borrowed to
fund a dissolved-oxygen deficit.

## Run the material-driven calculation

Compile the ordinary source packs, then use the existing example:

```bash
pyrrhotite_pack_dir="$(mktemp -d)"
for pack in pyrrhotite-fes-alpha-solid-nbs1989 pyrrhotite-fe7s8-solid-nbs1989 sulfur-orthorhombic-solid-nbs1989 oxygen-aqueous-nbs1989-thermochemistry iron-ii-ion-aqueous-nbs1989 iron-iii-ion-aqueous-nbs1989 sulfate-ion-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$pyrrhotite_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction pyrrhotite-fe7s8-oxidation-dissolved-oxygen --gibbs \
  "$pyrrhotite_pack_dir/pyrrhotite-fe7s8-solid-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/oxygen-aqueous-nbs1989-thermochemistry.fsmatpk" \
  "$pyrrhotite_pack_dir/water-liquid-nbs1989-thermochemistry.fsmatpk" \
  "$pyrrhotite_pack_dir/iron-ii-ion-aqueous-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/sulfate-ion-aqueous-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/hydrogen-ion-aqueous-nbs1989.fsmatpk" 647.441

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction ferric-pyrrhotite-fes-leach --gibbs \
  "$pyrrhotite_pack_dir/pyrrhotite-fes-alpha-solid-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/iron-iii-ion-aqueous-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/iron-ii-ion-aqueous-nbs1989.fsmatpk" \
  "$pyrrhotite_pack_dir/sulfur-orthorhombic-solid-nbs1989.fsmatpk" 87.911
```

`--help` lists all four options and exact pack order. The final argument is
mass of the declared mineral actually reacted, not total ore throughput.
`cargo test -p xtask --test mining_pyrrhotite_sources` compiles the original
manifests, verifies all fifteen literal facts and source receipts, checks
source/phase/composition refusal, and exercises all four reaction sums plus
the finite shared-oxygen stage. The ordinary stoichiometry tests also prove
exact conservation and finite ferric/oxygen reagent limits.
