# Zinc extraction and refining reference chemistry

The zinc source set connects sphalerite ore to zincite calcine, dissolved Zn(II)
and solid zinc. It also accounts for zinc dissolution by ferric oxidant in a
mixed sulfide leach and copper removal by zinc-dust cementation. All five
endpoints use the existing exact atom/charge,
finite-inventory, source-enthalpy and source-Gibbs consumers in `fs-thermochem`.

## Source values and validity

The five packs transcribe the indicated rows from Robie and Hemingway (1995),
[USGS Bulletin 2131](https://doi.org/10.3133/b2131),
[primary report](https://pubs.usgs.gov/bul/2131/report.pdf).
Each manifest cites an original page image. Values below retain source precision.

| Pack | Printed page | M [g/mol] | Hf° [kJ/mol] | Gf° [kJ/mol] | S° [J/(mol K)] |
|---|---:|---:|---:|---:|---:|
| `sphalerite-solid-usgs-b2131-thermochemistry` | 13 | 97.456 | -204.1 | -199.6 | 58.7 |
| `zincite-solid-usgs-b2131-thermochemistry` | 20 | 81.389 | -350.5 | -320.4 | 43.2 |
| `zinc-ion-aqueous-usgs-b2131` | 10 | 65.39 | -153.4 | -147.3 | -109.8 |
| `zinc-solid-usgs-b2131` | 10 | 65.39 | 0.0 | 0.0 | 41.63 |
| `sulfur-dioxide-gas-usgs-b2131` | 18 | 64.065 | -296.8 | -300.1 | 248.2 |

Every claim is restricted to **298.15 K and 100000 Pa**, the named phase and
`source-elemental-reference-usgs-b2131=1`. Sphalerite and zincite require their
own phase flags. The Zn(II) ion additionally requires the individual-ion,
infinite-dilution unit-activity and H+(aq)-zero conventions. Its negative
entropy is conventional ionic entropy. Formation Gibbs values come directly
from the independent source column; they are not `Hf - T * species entropy`.

Sphalerite uses the standard-state row on printed page 13, with orthorhombic
elemental sulfur. The different formation values in the high-temperature
table on printed page 170 use an S2(g) reference and are not substituted.
The SO2 pack keeps the roasting calculation on the same USGS convention as
its mineral participants; the separate JANAF evaluation retains its identity.

The 33 claims include dimensional metadata for the printed Hf/Gf/S uncertainty
magnitudes. Solid zinc supplies only its printed entropy magnitude; the two
elemental formation zeros are not assigned invented error bars. All claims
retain `Unstated` uncertainty, with no fabricated probability or covariance.

## Five explicit processing endpoints

The roasting, leaching and electrolysis stages follow section 2.3 of the
[Australian NPI zinc processing manual](https://www.dcceew.gov.au/sites/default/files/documents/fzinc.pdf).
The leach and productive cell equations below omit sulfate spectators. The
ferric endpoint is equation 1 in [Fowler and Crundwell (1999)](https://doi.org/10.1128/AEM.65.12.5285-5292.1999).

| CLI reaction | Integer equation | ΔrH° [J/mol equation] | ΔrG° [J/mol equation] |
|---|---|---:|---:|
| `sphalerite-roasting` | 2 ZnS(s) + 3 O2(g) → 2 ZnO(s) + 2 SO2(g) | -886400 | -841800 |
| `zincite-acid-leach` | ZnO(s) + 2 H+(aq) → Zn2+(aq) + H2O(l) | -88700 | -64000 |
| `zinc-electrowinning` | 2 Zn2+(aq) + 2 H2O(l) → 2 Zn(s) + O2(g) + 4 H+(aq) | 878400 | 768800 |
| `ferric-sphalerite-leach` | ZnS(s) + 2 Fe3+(aq) → Zn2+(aq) + 2 Fe2+(aq) + S(s) | -31700 | -94300 |
| `copper-cementation-with-zinc` | Cu2+(aq) + Zn(s) → Cu(s) + Zn2+(aq) | -218300 | -212400 |

Copper cementation uses the existing USGS copper-metal and Cu(II) packs.
The NPI manual section 2.3.4 describes zinc-dust displacement during solution
purification. The declared endpoint consumes one mole of metallic zinc and
adds one mole of dissolved Zn(II) for each mole of copper removed. Abundant
dissolved zinc cannot replace the metallic reductant. This single-species
balance supplies neither impurity selectivity nor complete electrolyte
purification.

One roasting extent, two leach extents and one electrowinning extent cancel
ZnO, Zn(II), protons and water, leaving `2 ZnS + 2 O2 -> 2 Zn + 2 SO2`.
The standard Hess sums are **-185400 J/mol equation** for enthalpy and
**-201000 J/mol equation** for Gibbs energy. This cancellation is across
sequential stages: a simultaneous `ReactionNetwork` request cannot spend
calcine, dissolved zinc or acid that a later operation has not produced yet.

Zinc electrowinning transfers two electrons per deposited zinc atom and four
per written equation. The source-based signed reversible potential is
**-1.99201262793359 V**, from `-ΔrG°/(4F)` in the written direction. This is
not an operating cell voltage. `ElectrolysisReaction` separately admits passed
charge, caller-supplied current efficiency and finite zinc/water inventories;
unallocated charge has no inferred side-reaction chemistry.

## Run the source-backed cases

From the repository root, compile the five new packs and their eight existing
reference partners into a new output directory:

```bash
zinc_pack_dir=$(mktemp -d)
for pack in sphalerite-solid-usgs-b2131-thermochemistry zincite-solid-usgs-b2131-thermochemistry zinc-ion-aqueous-usgs-b2131 zinc-solid-usgs-b2131 sulfur-dioxide-gas-usgs-b2131 oxygen-gas-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131 iron-iii-ion-aqueous-usgs-b2131 iron-ii-ion-aqueous-usgs-b2131 sulfur-orthorhombic-solid-usgs-b2131 copper-ii-ion-aqueous-usgs-b2131 copper-solid-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$zinc_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction sphalerite-roasting --gibbs \
  "$zinc_pack_dir/sphalerite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$zinc_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zincite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$zinc_pack_dir/sulfur-dioxide-gas-usgs-b2131.fsmatpk" 0.194912

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction zincite-acid-leach --gibbs \
  "$zinc_pack_dir/zincite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$zinc_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/water-liquid-usgs-b2131.fsmatpk" 0.081389

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction zinc-electrowinning --gibbs \
  "$zinc_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zinc-solid-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 0.13078

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction ferric-sphalerite-leach --gibbs \
  "$zinc_pack_dir/sphalerite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$zinc_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/sulfur-orthorhombic-solid-usgs-b2131.fsmatpk" 0.097456

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction copper-cementation-with-zinc --gibbs \
  "$zinc_pack_dir/copper-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zinc-solid-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/copper-solid-usgs-b2131.fsmatpk" \
  "$zinc_pack_dir/zinc-ion-aqueous-usgs-b2131.fsmatpk" 0.063546
```

The final input is kilograms of the first species actually consumed. For the
electrowinning case this is dissolved Zn(II) mass, not electrolyte mass. The CLI
prints each source identity, exact conservation identity, masses, source-mass
rounding residual and standard reaction energies. `--activity 'SPECIES=VALUE'`
accepts a complete explicit set of dimensionless activities for a local Gibbs
and reversible-potential evaluation.

`cargo test -p xtask --test mining_zinc_sources` checks all literal values and
state/phase refusals, the actual-pack Hess cycle, finite sequential inventory,
zinc-limited copper cementation and source-mass Faraday yield with excess-charge
refusal.

These reference endpoints do not supply a roast temperature or sensible heat,
conversion, zinc-ferrite formation, solid-solution properties, selective leach
kinetics, solution purification, passivation or plant recovery. Iron-bearing
calcines and zinc electrolytes need additional phase and impurity models;
ordinary point thermochemistry cannot infer those process outcomes.
