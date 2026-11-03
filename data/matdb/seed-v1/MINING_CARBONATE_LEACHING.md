# Copper-carbonate leaching and acidic calcite demand

The malachite and neutral carbonic profiles complete two existing conserved
reaction endpoints with material-backed mass, enthalpy and Gibbs calculations:

```text
Cu2CO3(OH)2(s) + 4 H+(aq) -> 2 Cu[2+](aq) + H2CO3(aq) + 2 H2O(l)
CaCO3(s, calcite) + 2 H+(aq) -> Ca[2+](aq) + H2CO3(aq)
```

## Source and component basis

Both new profiles use [USGS Bulletin 2131](https://pubs.usgs.gov/bul/2131/report.pdf),
printed pages 25 (malachite) and 15 (un-ionized carbonic component), at
298.15 K and 100000 Pa. Each adds six claims: formula weight, formation
enthalpy, direct formation Gibbs energy, entropy, and separately dimensioned
reported enthalpy/entropy uncertainty magnitudes. All retain `Unstated`
confidence semantics. Existing proton, calcium, copper(II), calcite and water
profiles supply the other participants under the same USGS reference convention.

| New profile | Formula weight, g/mol | Formation enthalpy, kJ/mol | Formation Gibbs energy, kJ/mol | Entropy, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `malachite-solid-usgs-b2131-thermochemistry` | 221.116 | -1054.0 | -890.2 | 166.3 |
| `carbonic-acid-aqueous-usgs-b2131` | 62.025 | -699.7 | -623.2 | 184.7 |

**The carbonic component includes a water formula unit.** The original NBS
compilation cited by USGS defines its H2CO3 entry as equivalent to aqueous CO2
plus liquid water. This is explicit in the official
[NIST corrected NBS dataset](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
The pack preserves the USGS numerical evaluation on this formal component basis.
It does not supply properties or activity for independently resolved hydrated
molecular H2CO3. Its required `source-carbonic-reference-usgs-b2131=1` selector
also prevents accidentally substituting an ion or an unqualified CO2 profile.
Its entropy is a neutral-component entropy, not individual-ion entropy.

The equations therefore retain two explicit product waters for malachite.
Changing the carbon product to CO2(aq) requires three waters and appropriately
converted component data. Adding a separate hydration step to this H2CO3 basis
would double-count that component conversion.

## Compile and run

Use a fresh output directory; the compiler refuses to overwrite existing packs.

```bash
carbonate_pack_dir=$(mktemp -d)
for pack in malachite-solid-usgs-b2131-thermochemistry carbonic-acid-aqueous-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131 copper-ii-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131 calcite-solid-usgs-b2131-thermochemistry calcium-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$carbonate_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --gibbs --equilibrium --reaction malachite-carbonic-acid-leach \
  "$carbonate_pack_dir/malachite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$carbonate_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$carbonate_pack_dir/copper-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$carbonate_pack_dir/carbonic-acid-aqueous-usgs-b2131.fsmatpk" \
  "$carbonate_pack_dir/water-liquid-usgs-b2131.fsmatpk" 221.116

cargo run -p fs-thermochem --example mineral_processing -- \
  --gibbs --equilibrium --reaction calcite-carbonic-acid-neutralization \
  "$carbonate_pack_dir/calcite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$carbonate_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$carbonate_pack_dir/calcium-ion-aqueous-usgs-b2131.fsmatpk" \
  "$carbonate_pack_dir/carbonic-acid-aqueous-usgs-b2131.fsmatpk" 100.087

cargo test -p xtask --test mining_carbonate_leach_sources
```

## Outputs and limits

For 1000 mol reacted malachite, the consumer reports 221.116 kg mineral and
4.032 kg proton equivalents consumed; 127.092 kg dissolved Cu(II), 62.025 kg
formal carbonic component and 36.030 kg explicit water produced. The standard
enthalpy change is **-87.5 MJ** and direct standard Gibbs energy is
**-77.0 kJ/mol extent**. With 10 mol malachite and 6 mol H+, only 1.5 mol extent
is available: 8.5 mol mineral remain and 3 mol dissolved Cu(II) form.

Exact atoms and charge close. The literal printed formula weights produce a
**-0.001 kg** product-minus-reactant residual at 1000 mol extent. The consumer
reports this source rounding; it does not change the malachite formula weight
to force numerical mass equality.

For calcite's carbonic endpoint, the reaction consumes **two proton equivalents**
per formula unit, with **-35.3 kJ/mol** standard reaction enthalpy and
**-48.3 kJ/mol** direct standard Gibbs energy. The existing bicarbonate endpoint
consumes one proton equivalent. Endpoint choice must come from an appropriate
solution model; these data do not automatically select pH or speciation.

The `--equilibrium` caloric route independently combines the retained enthalpy
and entropy values: -76.28956 kJ/mol for malachite and -46.2689385 kJ/mol for
calcite. The direct-Gibbs route uses separate tabulated formation values. The
two results are retained explicitly; no fitted offset or propagated confidence
interval is invented to reconcile the compilation's independent properties.

These are declared reacted amounts and finite inventory ceilings. They do not
predict ore recovery, kinetics, gas evolution, heat of acid dilution, plant
temperature, or commercial sulfuric-acid solution dose. A kinetic model must
also identify this same component basis and supply its own calibrated rates.
