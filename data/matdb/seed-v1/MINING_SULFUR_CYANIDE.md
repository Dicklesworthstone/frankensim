# Silver-sulfide leaching and sulfur cyanide demand

These three NBS source packs connect silver-sulfide dissolution to the sulfur
species that consume additional cyanide. The runnable finite-inventory example
shares cyanide and dissolved oxygen with the existing gold endpoint, then
accounts for a separately requested thiosulfate reaction. This supports
FrankenMining plan §§17.3 and 17.7's mineralogy-dependent reagent accounting.
Native silver, silver sulfide, free CN-, metal-bound cyanide and thiocyanate
remain distinct inventories throughout.

## Matched original reference data

The literal rows come from [NIST's NBS tables with the 1989 errata](https://doi.org/10.18434/M32124)
and its [corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
The CSV SHA-256 is
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.
The [original publication](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf)
defines the reference conventions. All fourteen numerical claims require
**298.15 K and 100000 Pa** and `source-elemental-reference-nbs1989=1`.

| Source pack | Literal NBS identity | M, g/mol | Formation H, kJ/mol | Direct formation G, kJ/mol | S, J/(mol K) | Cp, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `silver-sulfide-alpha-solid-nbs1989` | Table 37:Ag, Ag2S, alpha orthorhombic, cr | 247.8040 | -32.59 | -40.67 | 144.01 | 76.53 |
| `thiosulfate-ion-aqueous-nbs1989` | Table 14:S, S2O3(2-), ao | 112.1262 | -652.3* | -522.5 | 67. | unavailable |
| `thiocyanate-ion-aqueous-nbs1989` | Table 23:C, SCN-, ao | 58.0819 | 76.44 | 92.71 | 144.3 | -40.2 |

The source phase is deliberately named
`Ag2S(s-nbs-alpha-orthorhombic)` in the reaction. It retains the historical
NBS description and its required phase selector. The pack does not relabel
that source as acanthite, argentite, the separate NBS beta row or a measured
ore specimen. Its 298.15 K formation enthalpy is distinct from the adjacent
0 K value.

The thiosulfate enthalpy's asterisk marks an erratum correction, as explained
in [Reed's digitization methods](https://doi.org/10.6028/jres.125.007).
The numerical value is -652.3 kJ/mol; the marker remains in the source notes.
It is not an uncertainty magnitude. Uncertainty for every imported value
remains `Unstated`; no joint distribution or covariance is invented.

Both aqueous ions use hypothetical unit molality at infinite dilution and the
H+(aq)-zero convention. Their entropies use
`conventional-standard-molar-entropy`. The negative SCN- partial molar Cp uses
`conventional-standard-molar-isobaric-heat-capacity`; it must not be used as a
bulk-solution heat capacity or clamped to zero. The two available Cp values
are reference points, not temperature laws. The absent thiosulfate Cp remains
unavailable.

## Three conserved endpoints

```text
2 Ag2S(s-nbs-alpha-orthorhombic) + 8 CN-(aq) + 2 O2(aq) + H2O(l)
    -> 4 Ag(CN)2-(aq) + S2O3(2-)(aq) + 2 OH-(aq)

S2O3(2-)(aq) + CN-(aq) -> SCN-(aq) + SO3(2-)(aq)

S(s) + CN-(aq) -> SCN-(aq)
```

The first selects the sulfur product observed in
[Luna and Lapidus, *Cyanidation kinetics of silver sulfide* (2000)](https://doi.org/10.1016/S0304-386X(00)00072-4).
Their synthetic-silver-sulfide experiments identified thiosulfate and the
dicyanoargentate/hydroxide endpoint. This supplies the pathway choice; it
does not establish a universal product distribution for an ore.

The sulfur-transfer pathways are discussed in
[Luthy and Bruce, DOE report FE-2496-32](https://www.osti.gov/servlets/purl/6166154)
and their [1979 research paper](https://doi.org/10.1021/es60160a016).
Their work distinguishes the behavior of different reduced-sulfur species.
No empirical rate coefficient from that study is imported here. The
elemental-sulfur equation is normalized per mole S to match the existing
`sulfur-orthorhombic-solid-nbs1989` reference, rather than per mole S8.

| Library variant | Reaction ID | Standard reaction H, J/mol equation | Direct standard reaction G, J/mol equation |
| --- | --- | ---: | ---: |
| `SilverSulfideCyanidationThiosulfate` | `silver-sulfide-cyanidation-thiosulfate` | -861478 | -708519 |
| `ThiosulfateCyanideSulfurTransfer` | `thiosulfate-cyanide-sulfur-transfer` | -57360 | -43690 |
| `ElementalSulfurCyanideConsumption` | `elemental-sulfur-cyanide-consumption` | -74160 | -79690 |

All participants use the same NBS reference family. Existing gold,
dicyanoaurate, dicyanoargentate, free cyanide, hydroxide, dissolved oxygen,
water, sulfite and orthorhombic sulfur packs are reused unchanged. The `S(s)`
and `SO3[2-](aq)` identifiers connect directly to the earlier pyrrhotite and
sulfite endpoints. Exact atom/charge conservation is checked by the original
`ConservedReaction` and `ReactionNetwork` owners.

The independently tabulated Gibbs column is retained alongside the separate
H-minus-T-times-S calculation. Those results differ because the original
reference columns are not adjusted to force agreement. Negative standard
Gibbs energy does not select reaction progress or determine a rate.

## Run shared cyanide and oxygen accounting

```bash
sulfur_pack_dir=$(mktemp -d)
for pack in silver-sulfide-alpha-solid-nbs1989 \
  thiosulfate-ion-aqueous-nbs1989 thiocyanate-ion-aqueous-nbs1989 \
  sulfur-orthorhombic-solid-nbs1989 sulfite-ion-aqueous-nbs1989 \
  gold-solid-nbs1989-thermochemistry dicyanoaurate-ion-aqueous-nbs1989 \
  dicyanoargentate-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 oxygen-aqueous-nbs1989-thermochemistry \
  water-liquid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$sulfur_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example sulfur_cyanide_demand -- \
  "$sulfur_pack_dir" 4 2 1 18 3 10 1 1 1 1
```

The six initial amounts are Au, Ag2S, elemental S, free CN-, dissolved O2 and
water, in moles. The final four arguments are the requested extents of gold
cyanidation, silver-sulfide cyanidation, elemental-sulfur transfer and
thiosulfate transfer. An extent is per complete integer equation, including
the two Ag2S formula units in the silver reaction. All initial products are
explicitly zero in this example.

The first three reactions share one initial inventory. Their combined
withdrawals must fit that inventory. The second stage consumes the first
stage's thiosulfate and remaining CN-. This is the existing conservative
gross-inventory admission policy with an explicitly declared stage ordering;
it is not a physical claim that intermediate products cannot be consumed
during simultaneous chemistry. A kinetic or coupled-equilibrium treatment
would need its own supported model.

For these **constructed extents**, the calculation gives:

| Quantity | Initial | After gold/silver/sulfur stage | After thiosulfate stage |
| --- | ---: | ---: | ---: |
| Free CN-, mol | 18 | 1 | 0 |
| Dissolved O2, mol | 3 | 0 | 0 |
| Au(CN)2-, mol | 0 | 4 | 4 |
| Ag(CN)2-, mol | 0 | 4 | 4 |
| Thiosulfate, mol | 0 | 1 | 0 |
| Thiocyanate, mol | 0 | 1 | 2 |
| Sulfite, mol | 0 | 0 | 1 |
| Water, mol | 10 | 7 | 7 |
| Hydroxide, mol | 0 | 6 | 6 |

The final dissolved metal complexes contain 16 CN units; thiocyanate contains
the other two. Those two units are not available as free cyanide. All four Au
atoms, four Ag atoms and three S atoms remain accounted for in moles.
The source-basis total mass stays **2.0600126 kg**. Standard chemical enthalpy
changes are **-1507854 J** and **-57360 J**, totaling **-1565214 J**.
The example prints each stage's gross consumption/formation, final species
masses, original pack identities/hashes and both Gibbs diagnostics.

Changing the initial cyanide argument from 18 to 17 permits the first stage
but refuses the requested second stage. Changing oxygen from 3 to 2 refuses
the combined first-stage demand even though either metal endpoint is
individually feasible. The example admits both stages before printing an
accepted result, so either failure produces no partial success report.

## Scope and focused verification

`cargo test -p xtask --test mining_sulfur_cyanide_sources` compiles the real
source files with the existing pack compiler. It checks all fourteen facts,
reference and phase refusals, preserved conventional negative Cp, source
receipt replay, three independent mass/heat/Gibbs balances, staged shared
stocks, and independently reconstructed Au/Ag/S/C/N/H/O and charge totals.
The stoichiometry tests also join the new sulfur pathway to the existing
pyrrhotite, gold and sulfite equations to check shared species definitions.

These are reference chemical balances for supplied extents. They do not
provide mineral accessibility, kinetics, a sulfur-product equilibrium,
polysulfide distributions, free-cyanide pH speciation, mass transfer,
carbon adsorption, metal recovery forecasts, heat of mixing, dilution or
operating-temperature duty. Net reactive-species charge is conserved;
spectator counterions and the full electroneutral solution remain outside the
standalone inventory. No generic ore assay or plant calibration is inferred.
