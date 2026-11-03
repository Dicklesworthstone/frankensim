# Matched copper-cyanide thermochemistry for finite ligand accounting

These six profiles supply a compatible NBS source set for Cu(I), Cu(II),
crystalline CuCN and the di-, tri- and tetracyanocuprate(I) ions. They address
FrankenMining plan §17.3's cyanide-soluble-copper input and §17.7's
feed-dependent reagent consumption: complexed cyanide cannot also be counted
as the free CN activity used by a gold dissolution law.

## Original source values

All new rows are in **Table 36:Cu, prepared 1968**, in Wagman et al.,
*The NBS Tables of Chemical Thermodynamic Properties*, as distributed by
[NIST with the 1989 errata](https://doi.org/10.18434/M32124).
The [corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv)
was checked against its [official SHA256](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv.sha256):
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.
The original publication is [NBS/JPCRD Supplement 2](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf).

| Profile | Species and source state | Molar mass, g/mol | Formation H, kJ/mol | Direct formation G, kJ/mol | S, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `copper-i-ion-aqueous-nbs1989` | Cu+, ao | 63.5400 | 71.67 | 49.98 | 40.6 |
| `copper-ii-ion-aqueous-nbs1989` | Cu2+, ao | 63.5400 | 64.77 | 65.49 | -99.6 |
| `copper-i-cyanide-solid-nbs1989` | CuCN, cr | 89.5579 | 96.2 | 111.3 | 84.5 |
| `dicyanocuprate-i-ion-aqueous-nbs1989` | Cu(CN)2-, ao | 115.5758 | unavailable | 257.8 | unavailable |
| `tricyanocuprate-i-ion-aqueous-nbs1989` | Cu(CN)3(2-), ao | 141.5937 | unavailable | 403.8 | unavailable |
| `tetracyanocuprate-i-ion-aqueous-nbs1989` | Cu(CN)4(3-), ao | 167.6116 | unavailable | 566.6 | unavailable |

The six profiles contain **18 literal claims**. Missing complex enthalpies and
entropies remain unavailable; no zero or synthetic value fills those columns.
The source's old rounded molar masses remain unchanged. All claims require
298.15 K, 100000 Pa and `source-elemental-reference-nbs1989`, with `Unstated`
uncertainty. Aqueous ions retain the unit-molality/infinite-dilution and
H+(aq)-zero conventions; their entropy is conventional. CuCN uses the pure
crystalline solid reference, without an inferred polymorph. The existing
USGS Cu(II) profile is a separate evaluation.

## Conserved formation and ligand demand

With the existing `cyanide-ion-aqueous-nbs1989` value
G(CN-)=172.4 kJ/mol, the source-backed cumulative formation equations are:

| Exactly conserved reaction | Standard reaction G, J/mol | log10 of formation constant |
| --- | ---: | ---: |
| Cu+ + 2 CN- -> Cu(CN)2- | -136980 | 23.997857736424443 |
| Cu+ + 3 CN- -> Cu(CN)3(2-) | -163380 | 28.622937633063407 |
| Cu+ + 4 CN- -> Cu(CN)4(3-) | -172980 | 30.304784868204848 |

Constants are calculated as `ln beta_n=-DeltaG/(R*T)` from the direct Gibbs
column. Each binds one Cu(I), consumes exactly n free cyanide ligands and
produces charge `1-n`. The stepwise di-to-tri and tri-to-tetra additions have
reaction G of -26400 and -9600 J/mol, respectively. The equations are built
through the original exact atom/charge owner; no new mineral-name lookup or
implicit redox conversion is required.

The solid dissolution equation `CuCN(s) -> Cu+ + CN-` has
reaction G=111080 J/mol and log10 K=-19.460374049949095. Adding the first
complex-formation equation gives `CuCN(s) + CN- -> Cu(CN)2-`, with reaction
G=-25900 J/mol. This Gibbs cycle preserves the same original participant
values and provides a thermodynamic diagnostic, not a precipitation or
dissolution-rate prediction.

For an explicitly ideal dilute model, `a_complex = beta_n*a_Cu*a_CN^n`.
The ligand variable here is **free CN-**, not HCN+CN, total cyanide including
complexes, or WAD-CN from an analytical method. A finite solution calculation
must close at least the Cu inventory, all cyanide ligands and charge balance
with water and HCN/CN equilibrium. The aqueous activity basis is independent
of solution volume. These data do not themselves predict that composition.

Cu(II) cannot silently enter the Cu(I)-cyanide family. A reduction reaction
and its redox partners would need their own exact conservation and sources.
The reference data also do not supply copper-mineral accessibility, cyanide
oxidation, thiocyanate/cyanate production, nonideal activities, temperature
dependence, carbon adsorption, kinetics or a calibrated reagent dose.

## Compile and exercise the actual data

```bash
copper_cyanide_pack_dir="$(mktemp -d)"
for pack in copper-i-ion-aqueous-nbs1989 copper-ii-ion-aqueous-nbs1989 copper-i-cyanide-solid-nbs1989 dicyanocuprate-i-ion-aqueous-nbs1989 tricyanocuprate-i-ion-aqueous-nbs1989 tetracyanocuprate-i-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$copper_cyanide_pack_dir/$pack.fsmatpk"
done
cargo test -p xtask --test mining_copper_cyanide_sources
```

The focused integration compiles the real TSVs, checks every literal row and
receipt, verifies formation constants and finite ligand inventories, checks
the stepwise/solid Gibbs cycles and free-ligand activity dependence, and
refuses missing heat data, unbalanced Cu(II) substitution and unsupported
source conditions. Its explicit-extent checks exercise conservation. The
additional coupled tests below solve finite solution composition.

## Coupled dissolved Cu(I), cyanide and pH

`fs_thermochem::stoichiometry::thermodynamics::acid_base::complexation`
provides the simultaneous ideal closure. Bind the three consecutive
`LigandAssociation` steps to the original conserved reactions and Gibbs
sources: Cu+ plus **two** CN- forms the dicyano ion; the next two steps each
add one CN-. Join them with `MetalLigandFamily::new`, then pass that family
and the original HCN/water `IdealAcidBaseSystem` into `IdealLigandSystem::new`.
The CN- source selected by the acid and metal steps must be identical.

Call `speciate` with the total cyanide ligand amount, every dissolved metal
total, any other acid-family totals, explicit spectator counterions, the
common amount at unit activity, numerical controls and cancellation callback.
`LigandSpeciation` retains free CN, HCN+CN, Cu-bound ligands, every copper
member, pH, original sources and independently recomputed balances. It solves
the finite ligand and charge equations together; it does not prescribe a
fixed cyanide consumption per mole of copper.

For an illustrative 0.010 mol cyanide pool on `n_ref=10 mol`, with 0.010 mol
Na+ and one mol chloride counterion per mol dissolved Cu(I), the actual NBS
rows give these **calculated ideal-model** results:

| Dissolved Cu(I), mol | Free CN-, mol | Cu-bound CN, mol | pH |
| ---: | ---: | ---: | ---: |
| 0 | 0.00876922811862 | 0 | 10.08542435335 |
| 0.002 | 0.00334128696158 | 0.00589899308257 | 9.87590079290 |

The inventories are authored comparison inputs, not an ore assay or operating
recommendation. The 0.002 mol copper case contains approximately
0.00013055666116 mol dicyano, 0.00183989359512 mol tricyano and
0.00002954974372 mol tetracyano complexes, plus a positive free-Cu+ tail.
The [actual-source coupled tests](../../../xtask/tests/mining_copper_cyanide_sources/speciation.rs)
compare against an independent high-precision balance reduction, check
ligand scarcity and empty pools, and retain tiny positive species amounts.
Near ligand exhaustion, pH is sensitive to the allowed ligand-balance error;
the reference tests tighten the ligand and logarithmic bracket tolerances.
Numerical residuals and retained brackets should accompany the calculated pH.

For the source's ideal unit-molality convention, `n_ref=10 mol` represents
10 kg solvent times the 1 mol/kg standard molality. It is not a 10 m3 solution
volume. Scaling all inventories and this basis together preserves activities;
changing the basis alone changes pH and the copper coordination distribution.
The tests exercise both cases. No solvent inventory or temperature-dependent
equilibrium, CuCN precipitation, Cu(I)/Cu(II) redox or dissolution of an ore
copper mineral is inferred. All supplied dissolved Cu(I) remains in the
explicit aqueous family until a separate supported process transfers it.
