# Copper-oxide acid leaching

Two compatible profiles complete the source set for a declared tenorite
acid-leaching endpoint:

`CuO(tenorite) + 2 H+(aq) -> Cu[2+](aq) + H2O(l)`

`MineralReaction::CupricOxideAcidLeach` builds this equation through the
existing exact atom/charge checker. The new profiles use
[USGS Bulletin 2131](https://pubs.usgs.gov/bul/2131/report.pdf), printed pages
6 (Cu²⁺) and 16 (tenorite), at 298.15 K and 100000 Pa. They combine with the
existing hydrogen-ion and liquid-water USGS profiles. The separate JANAF
CuO profile retains its own evaluation and reference selector.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Direct formation Gibbs energy, J/mol | Entropy, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `copper-ii-oxide-solid-usgs-b2131-thermochemistry` | 0.079545 | −156100 | −128300 | 42.6, absolute |
| `copper-ii-ion-aqueous-usgs-b2131` | 0.063546 | 64900 | 65100 | −98.0, conventional ionic |

The new tenorite pack adds six claims and reuses the existing six-claim Cu²⁺
profile. In addition to the four quantities above,
each retains separately dimensioned reported uncertainty magnitudes for
formation enthalpy and entropy: 2000 J/mol and 0.2 J/(mol K) for tenorite;
1000 J/mol and 4.0 J/(mol K) for Cu²⁺. The reported Gibbs magnitudes remain
in source observations. Every claim retains `Unstated` confidence semantics;
no independent error distributions or propagated intervals are invented.

Tenorite requires `source-phase-tenorite=1`, solid and pure-species selectors.
Cu²⁺ requires the aqueous unit-activity/infinite-dilution standard state,
individual-ion identity and H⁺ zero convention. Its negative entropy remains
`conventional-standard-molar-entropy`. All claims require the USGS elemental
reference, the typed reference temperature/pressure and the existing
engineering-reference/specimen selectors. No finite-concentration activity,
complexation, heat-capacity extrapolation or ore calibration is supplied.

## Run the complete balance

```bash
copper_leach_pack_dir=/tmp/mining-copper-leach
mkdir -p "$copper_leach_pack_dir"
for pack in copper-ii-oxide-solid-usgs-b2131-thermochemistry hydrogen-ion-aqueous-usgs-b2131 copper-ii-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$copper_leach_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --gibbs --equilibrium --reaction cupric-oxide-acid-leach \
  "$copper_leach_pack_dir/copper-ii-oxide-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$copper_leach_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$copper_leach_pack_dir/copper-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$copper_leach_pack_dir/water-liquid-usgs-b2131.fsmatpk" 79.545
cargo test -p xtask --test mining_copper_leach_sources
```

The consumer reports 1000 mol extent: **79.545 kg CuO + 2.016 kg H⁺ →
63.546 kg Cu²⁺ + 18.015 kg water**, with **−64.8 MJ** standard enthalpy change.
The molar heat is −64800 J/mol extent. An inventory of 10 mol CuO and 6 mol H⁺
permits at most 3 mol extent, leaving 7 mol CuO and producing 3 mol each of
Cu²⁺ and water. It gives −194400 J at that ceiling; 4 mol extent refuses.
These proton equivalents are net ionic demand, not commercial sulfuric-acid
solution mass or total acid consumption by an ore blend.

`ConservedReaction::standard_gibbs_change` and the CLI `--gibbs` route consume
the direct formation-Gibbs claims and return **−43700 J/mol extent**.
The existing CLI `--equilibrium`
option instead uses matched enthalpy and entropy: ΔS° = −70.6 J/(mol K),
giving ΔG° = **−43750.61 J/mol**. The 50.61 J/mol difference preserves the
independently rounded source columns; neither route silently adjusts a datum.
The integration tests exercise the direct-Gibbs route and verify every
load-bearing material receipt.

This two-reactant equation can also be passed to
`fs_material::reaction_batch::BilinearBatchReaction` with explicit CuO/H⁺ roles
and a separately supported kinetic card. The reference data themselves do not
predict leaching rate, liberation, equilibrium composition or actual recovery.
