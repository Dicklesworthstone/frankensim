# Zinc ferrite, acid demand and iron precipitation

Iron-bearing zinc calcine can contain zinc ferrite as well as zincite. Those
phases have different leach requirements: dissolving pure ZnFe2O4 releases two
ferric ions per zinc and consumes eight protons before subsequent iron
precipitation. This source set makes that zinc/iron/acid account explicit in
the existing conserved reactions and finite-inventory consumers.

The industrial context is described in sections 2.3.2–2.3.3 of the
[Australian NPI zinc processing manual](https://www.dcceew.gov.au/sites/default/files/documents/fzinc.pdf):
roasting can form zinc ferrite, which requires more reactive leach conditions
than zincite, followed by iron removal from the leach solution. The pure
goethite reference below is a declared phase, not an identification of every
industrial iron precipitate.

## Compatible, phase-specific source data

Both packs transcribe Robie and Hemingway (1995),
[USGS Bulletin 2131](https://doi.org/10.3133/b2131),
[primary report](https://pubs.usgs.gov/bul/2131/report.pdf).
The manifests link the original page images.

| Pack | Printed page | Formula | M [g/mol] | Hf° [kJ/mol] | Gf° [kJ/mol] | S° [J/(mol K)] |
|---|---:|---|---:|---:|---:|---:|
| `franklinite-solid-usgs-b2131-thermochemistry` | 22 | ZnFe2O4 | 241.082 | -1188.1 | -1082.1 | 150.7 |
| `hematite-solid-usgs-b2131-thermochemistry` | 16 | Fe2O3 | 159.692 | -826.2 | -744.4 | 87.4 |

The franklinite row supplies the **pure ZnFe2O4 endmember**. It does not
represent Mn-bearing natural franklinite, substituted spinels or an assayed
calcine mixture. The hematite values retain the USGS evaluation; substituting
the separate JANAF hematite pack would mix source conventions.

Each pack has seven claims: mass, Hf°, Gf°, S° and the three printed uncertainty
magnitudes. The Hf°/Gf°/S° magnitudes are 4.0 kJ/mol, 4.0 kJ/mol and
0.3 J/(mol K) for franklinite, and 1.3 kJ/mol, 1.3 kJ/mol and 0.2 J/(mol K)
for hematite. They remain dimensional metadata with `Unstated` uncertainty;
no coverage probability or covariance is invented.

All claims require **298.15 K, 100000 Pa**, the named pure solid phase and
`source-elemental-reference-usgs-b2131=1`. They are engineering reference
data, not specimen-qualified measurements. Independent source Gf° is kept
separate from the caloric `ΔrH° - T ΔrS°` route; rounded source columns need
not produce identical results.

## Three exact endpoints

The following values use the written integer equations and the same USGS
source family for every participant. Sulfate spectators are omitted.

| CLI reaction | Integer equation | ΔrH° [J/mol equation] | Direct ΔrG° [J/mol equation] |
|---|---|---:|---:|
| `zinc-ferrite-formation` | ZnO(s) + Fe2O3(s) → ZnFe2O4(s) | -11400 | -17300 |
| `zinc-ferrite-acid-leach` | ZnFe2O4(s) + 8 H+(aq) → Zn2+(aq) + 2 Fe3+(aq) + 4 H2O(l) | -208300 | -47000 |
| `zinc-ferrite-goethite-leach` | ZnFe2O4(s) + 2 H+(aq) → Zn2+(aq) + 2 FeO(OH)(s) | -90500 | -48800 |

Formation is the declared oxide-combination balance; it does not predict a
roaster mechanism or the fraction of zinc converted to ferrite. Dissolution
uses the existing USGS conventional aqueous-ion reference states.

The existing goethite endpoint is
`Fe3+(aq) + 2 H2O(l) -> FeO(OH)(s) + 3 H+(aq)`.
Adding two such extents to one dissolution extent gives the combined
equation, with enthalpy `-208300 + 2 × 58900 = -90500 J` and direct Gibbs
energy `-47000 + 2 × (-900) = -48800 J`.

This distinguishes **eight gross protons admitted to dissolution** from
**two net protons consumed after iron precipitation**. A staged
`ReactionNetwork` must have the acid and intermediates available at each
step. Its simultaneous operation cannot spend the six protons produced
later by hydrolysis. The combined equation is an aggregate account; it does
not choose the precipitate, determine pH or prescribe a plant acid dose.

## Run the source-backed cases

From the repository root, compile the two new packs and six existing partners:

```bash
ferrite_pack_dir=$(mktemp -d)
for pack in franklinite-solid-usgs-b2131-thermochemistry hematite-solid-usgs-b2131-thermochemistry zincite-solid-usgs-b2131-thermochemistry hydrogen-ion-aqueous-usgs-b2131 zinc-ion-aqueous-usgs-b2131 iron-iii-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131 goethite-solid-usgs-b2131-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$ferrite_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction zinc-ferrite-formation --equilibrium --gibbs \
  "$ferrite_pack_dir/zincite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$ferrite_pack_dir/hematite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$ferrite_pack_dir/franklinite-solid-usgs-b2131-thermochemistry.fsmatpk" 0.081389

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction zinc-ferrite-acid-leach --equilibrium --gibbs \
  "$ferrite_pack_dir/franklinite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$ferrite_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ferrite_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ferrite_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ferrite_pack_dir/water-liquid-usgs-b2131.fsmatpk" 0.241082

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction zinc-ferrite-goethite-leach --equilibrium --gibbs \
  "$ferrite_pack_dir/franklinite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$ferrite_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ferrite_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ferrite_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" 0.241082
```

The last argument is kilograms of the first species actually consumed,
not ore or slurry mass. These inputs each represent one reaction extent.
The CLI prints source identities and receipts, exact conservation, source-mass
rounding residuals and both thermodynamic routes.

`cargo test -p xtask --test mining_zinc_ferrite_sources` checks all 14 new
literal claims, incompatible phase/state refusals, actual-pack Hess closure,
finite ferrite formation, gross-acid admission and sequential iron removal.

These endpoints supply no high-temperature caloric extrapolation, kinetics,
accessible-zinc fraction, solid-solution activity, phase proportions,
impurity selectivity or process recovery. Those require additional measured
inputs or separately validated models.
