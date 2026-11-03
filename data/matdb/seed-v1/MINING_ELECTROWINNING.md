# Copper electrowinning: compatible sources and reversible potential

The USGS reference set now closes both the chemical energy and Faraday basis
of the existing oxygen-evolving copper electrowinning endpoint:

```text
2 Cu2+(aq) + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+(aq)
```

`MineralReaction::CopperElectrowinning` conserves atoms and charge exactly.
`ElectrolysisReaction` associates two transferred electrons with each deposited
copper atom, hence **four per written equation**. It separately supports
charge-driven productive inventories with explicit current efficiency.

## Source data

All entries use [USGS Bulletin 2131](https://doi.org/10.3133/b2131),
298.15 K, 100000 Pa and `source-elemental-reference-usgs-b2131=1`.
The two copper rows are on [printed page 6](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/12/).
Water, oxygen and H+ are compatible existing profiles; the elemental-reference
zeros remain required sources.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Formation Gibbs energy, J/mol | Standard entropy, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `copper-ii-ion-aqueous-usgs-b2131` | 0.063546 | 64900 | 65100 | -98.0, conventional ionic |
| `copper-solid-usgs-b2131` | 0.063546 | 0 | 0 | 33.14 |
| `water-liquid-usgs-b2131` | 0.018015 | -285800 | -237100 | 70.0 |
| `oxygen-gas-usgs-b2131` | 0.031999 | 0 | 0 | 205.15 |
| `hydrogen-ion-aqueous-usgs-b2131` | 0.001008 | 0 | 0 | 0, conventional ionic |

Cu2+ is the individual aqueous ion on the source's hypothetical unit-activity,
infinite-dilution and H+-zero reference. It is not total dissolved copper,
copper sulfate solution or a copper complex. Copper solid is the elemental
reference phase. Their phase and reference selectors belong to the source
claims and cannot be supplied by a caller to legitimize an incompatible pack.

The source's printed uncertainty magnitudes remain distinct from probability
models. Cu2+ enthalpy/entropy and solid-copper entropy magnitudes are separately
dimensioned metadata; the 0.1 kJ/mol Cu2+ Gibbs magnitude is retained in its
observation. All claims retain `Unstated` uncertainty semantics. Source
rounding does not establish an error interval for a derived voltage.

## Chemical driving force and voltage

The directly tabulated formation values give, per written equation:

- Standard reaction enthalpy: **+441800 J/mol extent**.
- Standard reaction Gibbs energy: **+344000 J/mol extent**.
- Natural log of the dimensionless equilibrium constant: **-138.76803667515603**.

The [IUPAC standard cell-potential relation](https://goldbook.iupac.org/terms/view/S05913)
is `E° = -Delta_r G°/(n F)`. `ElectrolysisReaction::standard_potential`
evaluates the original reaction's source set and returns the original model,
Gibbs evidence, receipts and signed typed voltage. Its `at_activities` method
retains the same basis and uses `Delta_r G = Delta_r G° + R T ln Q`.

The standard potential for the written electrowinning direction is
**-0.891327190438547 V**. The negative sign identifies the externally driven
direction. Its opposing magnitude is about **0.8913 V** at these standard
states. Setting only proton activity to 0.01, with every other activity
explicitly one, changes the signed reversible potential to
**-0.7730084910689824 V**. These digits make the numerical regression
reproducible; they do not add precision to the underlying source values.

This is the [zero-current reversible potential](https://goldbook.iupac.org/terms/view/E01934).
It is not the operating voltage at finite current. No overpotential, ohmic
loss, concentration-to-activity model or current-efficiency factor is inferred.
The supplied activities define a local state, so its Gibbs derivative is not
multiplied by a finite extent whose activities would change.

## Run the source-backed case

```bash
ew_pack_dir="$(mktemp -d)"
for pack in copper-ii-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131 copper-solid-usgs-b2131 oxygen-gas-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$ew_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction copper-electrowinning --gibbs \
  --activity 'Cu[2+](aq)=1' --activity 'H2O(l)=1' \
  --activity 'Cu(s)=1' --activity 'O2(g)=1' --activity 'H+(aq)=0.01' \
  "$ew_pack_dir/copper-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$ew_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$ew_pack_dir/copper-solid-usgs-b2131.fsmatpk" \
  "$ew_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$ew_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 63.546
cargo test -p xtask --test mining_copper_gibbs_sources
```

The final mass is 63.546 kg of actually reacted Cu2+ chemical equivalent,
not electrolyte mass. All active species need explicit activities, including
the metal and water. The source integration tests also connect a productive
4F charge to one equation extent and 0.127092 kg copper, verify the retained
source receipts and reject missing reference participants or excess charge.
