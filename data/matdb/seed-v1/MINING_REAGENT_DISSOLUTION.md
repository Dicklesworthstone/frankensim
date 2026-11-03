# Sodium reagent dissolution and counterion accounting

The NBS sodium set connects declared dissolved masses of NaCN and NaOH to
their cyanide/base equivalents **and their sodium counterions**. It supplies
compatible standard reaction energies for the existing NBS cyanide, hydroxide
and shared-pH consumers. A dissolution extent is an explicit input; these
source rows do not determine how much of a commercial reagent dissolves.

## Source values and phase identity

The primary source is the [official corrected NIST NBS
CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv),
[dataset DOI](https://doi.org/10.18434/M32124), **Table 99:Na, Sodium (Prepared
1980)**. Its header fixes these values at 298.15 K and 0.1 MPa (1 bar).
The [original NBS publication](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf),
sections 5.5–5.8, supplies the reference conventions; the
[digitization paper](https://doi.org/10.6028/jres.125.007) describes the edition
with 1989 errata.

| Pack | Source state | M [g/mol] | Hf° [kJ/mol] | Gf° [kJ/mol] | S° [J/(mol K)] |
|---|---|---:|---:|---:|---:|
| `sodium-ion-aqueous-nbs1989` | Na+, `ao` | 22.9898 | -240.12 | -261.905 | 59.0 |
| `sodium-hydroxide-solid-nbs1989-thermochemistry` | NaOH, `cr` | 39.9972 | -425.609 | -379.494 | 64.455 |
| `sodium-cyanide-cubic-solid-nbs1989-thermochemistry` | NaCN, `cr I, cubic` | 49.0077 | -87.49 | -76.43 | 115.60 |

There are four claims per pack, **12 total**, with `Unstated` uncertainty.
All retain the original source precision, observation receipts, phase,
**298.15 K**, **100000 Pa** and `source-elemental-reference-nbs1989=1`.
The NBS/1989 label names the distributed evaluation, not a measurement date.

Na+ uses conventional individual-ion entropy and the H+(aq)=0 reference.
Its `ao` state is hypothetical unit molality referred to infinite dilution;
it is not a measured 1 mol/kg solution or a molarity convention. The solids
are pure reference constituents, not hydrated or commercial products.

The NaCN pack requires `source-phase-sodium-cyanide-cubic=1` on every claim.
The adjacent orthorhombic `cr II` row has a different enthalpy and lacks G/S;
its enthalpy is not combined with the cubic row. Existing JANAF/NIST NaOH and
NaCN solid packs retain their separate evaluations and caloric uses. They are
not relabeled as NBS to complete this aqueous reaction set.

## Two conserved dissolution endpoints

| CLI reaction | Written equation | ΔrH° [J/mol equation] | Direct ΔrG° [J/mol equation] | ΔrS° [J/(mol K)] |
|---|---|---:|---:|---:|
| `sodium-hydroxide-dissolution` | NaOH(s) → Na+(aq) + OH−(aq) | -44505 | -39655 | -16.205 |
| `sodium-cyanide-dissolution` | NaCN(s, cubic) → Na+(aq) + CN−(aq) | -2030 | -13075 | 37.5 |

Each written extent consumes one mole of the declared solid and produces one
mole of each ion. The existing exact atom/charge checker and finite-inventory
owner admit these extents; a request exceeding the available solid refuses.
The source molar masses also close independently: Na+ plus OH− is 39.9972
g/mol, and Na+ plus CN− is 49.0077 g/mol.

`ΔrH° - T ΔrS°` gives -39673.47925 J/mol for NaOH and -13210.625 J/mol for
NaCN. These caloric results remain distinct from the directly tabulated
formation-G route. The consumer does not adjust source data to make them
coincide. Standard dissolution enthalpy is not a temperature-rise prediction
or a finite-concentration heat-of-mixing model.

## Execute the actual source calculations

From the repository root:

```bash
reagent_pack_dir=$(mktemp -d)
for pack in sodium-ion-aqueous-nbs1989 sodium-hydroxide-solid-nbs1989-thermochemistry sodium-cyanide-cubic-solid-nbs1989-thermochemistry hydroxide-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$reagent_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction sodium-hydroxide-dissolution --equilibrium --gibbs \
  "$reagent_pack_dir/sodium-hydroxide-solid-nbs1989-thermochemistry.fsmatpk" \
  "$reagent_pack_dir/sodium-ion-aqueous-nbs1989.fsmatpk" \
  "$reagent_pack_dir/hydroxide-ion-aqueous-nbs1989.fsmatpk" 0.0399972

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction sodium-cyanide-dissolution --equilibrium --gibbs \
  "$reagent_pack_dir/sodium-cyanide-cubic-solid-nbs1989-thermochemistry.fsmatpk" \
  "$reagent_pack_dir/sodium-ion-aqueous-nbs1989.fsmatpk" \
  "$reagent_pack_dir/cyanide-ion-aqueous-nbs1989.fsmatpk" 0.0490077

cargo test -p xtask --test mining_reagent_dissolution_sources
```

The final CLI argument is kilograms of pure solid **declared actually
consumed**, not slurry mass or an unqualified commercial dose. Each example
represents one mole of reaction extent. Source identity, receipts, exact
conservation, mass accounting and both energy routes remain visible.

## Connect the dissolved amounts to pH

The integration test compiles the three new sources and the five original NBS
CN−, OH−, HCN, H+ and water packs. Declared dissolved masses of 4.90077 g cubic
NaCN and 0.799944 g NaOH give 0.1 mol cyanide/sodium and 0.02 mol
hydroxide/sodium respectively. The coupled acid/base owner then uses the
uncomplexed cyanide total, explicit sodium spectator amount and an explicit
10 kg solvent basis, corresponding to 10 mol at unit molal activity.

Adding the declared caustic increases pH and free cyanide while conserving
the HCN+CN− analytical pool and charge. Hydroxide is re-solved by the shared
water/protonation owner; adding it again as a fixed spectator would count the
base twice. Dissolved cyanide is not automatically all free CN−, and neither
the uncomplexed pool nor this calculation represents total or WAD cyanide.

The tests also verify the 12 literal source claims, phase/state/reference
refusals, separate direct/caloric energies, finite solid availability and
the source-mass-to-counterion calculation. These checks do not establish
commercial purity, dissolution kinetics, concentrated-electrolyte activity,
plant pH calibration, cyanide destruction or an optimal reagent dose.

## Add reagents to an active leach bath

`cyanide_leach::dosing::DissolvedSodiumReagent` is accepted by the closed
`GoldCyanidationBatch`, finite-headspace `gas_transfer::AeratedGoldCyanidation`
and prescribed-pressure `aeration::AeratedGoldCyanidationBatch` owners:

```rust
use fs_material::cyanide_leach::dosing::DissolvedSodiumReagent;
use fs_qty::Amount;

let added = bath.dose_reagent(
    DissolvedSodiumReagent::SodiumCyanide,
    Amount::new(0.001),
    1_000_000,
    || false,
)?;
let totals = bath.reagent_inputs();
assert_eq!(added.inputs, totals);
```

The argument is an amount in moles **already dissolved into the bath**. NaCN
adds the same amount of analytical cyanide and sodium. NaOH adds sodium/base
equivalents. The shared equilibrium owner determines the new pH, HCN, free CN
and OH; optional competing-metal inventories remain fixed while their bound
cyanide repartitions. The initial charge stays immutable. `ReagentInputs`
tracks cumulative NaCN, NaOH and sodium additions separately, and dose and
advance receipts report independent cyanide and sodium balance residuals.

Dosing preserves the accepted clock, gold dissolution, actual oxygen stocks
and all existing gas input/transfer accounts. It can replenish a cyanide-limited
bath without resetting its reaction history. A closed bath still has only its
original dissolved oxygen available. Each complete candidate is admitted by
the original chemical source; finite-headspace dosing also checks its retained
transfer source at the changed chemistry. The prescribed-pressure owner checks
the next call's pressure boundary and transfer source on the next advance.
Invalid or unresolved additions, source refusal, work exhaustion and
cancellation leave the accepted bath and input accounts unchanged.

The actual-source `aerated_gold_leach` example now connects the mass calculation
above to this active-bath API. Add the four following packs to the same
`reagent_pack_dir` created above, then run the measured-Henry route:

```bash
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry oxygen-water-henry-sander2023-rettich2000; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$reagent_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-material --example aerated_gold_leach -- \
  --materials "$reagent_pack_dir" \
  --henry-pack "$reagent_pack_dir/oxygen-water-henry-sander2023-rettich2000.fsmatpk" \
  --dose-reagents
```

At the existing 120-second interruption, the example declares 0.0490077 g
NaCN and 0.0399972 g NaOH actually dissolved, **0.001 mol of each**. It resolves
the compatible NBS reaction and molar-mass receipts, adds both reagents, prints
the new pH/free CN and cumulative reagent accounts, and resumes the original
contact-stop, oxygen-dose and contact-restart scenario. The existing Gibbs
oxygen-partition route also accepts `--dose-reagents`; its oxygen source
compilation is in [MINING_OXYGEN_TRANSFER.md](MINING_OXYGEN_TRANSFER.md).

These masses, source boxes and kinetic parameters are numerical illustrations,
not operating recommendations or measured ore calibration. Neither this API
nor the example predicts solid dissolution completion, stock-solution carrier
volume, commercial purity, finite-concentration activities or heat feedback.
Temperature and volume remain fixed. The pure-water Henry datum remains an
explicit approximation when used for cyanide liquor.
