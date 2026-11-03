# Cassiterite reduction and carbon-gas reference chemistry

Five material packs supply a coherent source family for tin reduction and
the accompanying graphite/CO/CO2 balance. Three `MineralReaction` endpoints
consume those packs through the existing exact conservation, finite inventory,
formation-enthalpy and formation-Gibbs APIs. The original mineral-density
pack remains separate: mineral density alone supplies none of these energies.

## Original reference data

The following literal rows were checked against original page images in
Robie and Hemingway (1995), [USGS Bulletin 2131](https://doi.org/10.3133/b2131),
available in the [USGS report](https://pubs.usgs.gov/bul/2131/report.pdf) and
its [UNT archival copy](https://digital.library.unt.edu/ark:/67531/metadc1114954/).
All values describe **298.15 K and 100000 Pa**.

| Pack | Printed page | M, g/mol | Hf°, kJ/mol | Gf°, kJ/mol | S°, J/(mol K) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `cassiterite-solid-usgs-b2131-thermochemistry` | 19 | 150.709 | -577.6 | -515.8 | 49.0 |
| `tin-white-solid-usgs-b2131` | 9 | 118.710 | 0.0 | 0.0 | 51.18 |
| `graphite-solid-usgs-b2131` | 5 | 12.011 | 0.0 | 0.0 | 5.74 |
| `carbon-monoxide-gas-usgs-b2131` | 15 | 28.010 | -110.5 | -137.1 | 197.3 |
| `carbon-dioxide-gas-usgs-b2131` | 15 | 44.010 | -393.5 | -394.4 | 213.8 |

The 31 claims retain these four properties per species and eleven printed
uncertainty magnitudes as dimensional metadata. All claims keep `Unstated`
uncertainty: source magnitudes do not establish a joint probability model or
certified reaction-energy interval. The gas entropy magnitudes print as 0.0;
that rounding is preserved without asserting exact entropy. Formation zeros
for elemental graphite and white tin follow the source reference convention.

Cassiterite requires `source-phase-cassiterite=1`, tin requires
`source-phase-tin-white-tetragonal=1`, and carbon requires
`source-phase-graphite=1`, in addition to the common source-state selectors.
Liquid tin, grey tin, diamond, coke and coal cannot substitute for these named
reference phases. The packs use `source-elemental-reference-usgs-b2131=1`;
the existing JANAF CO/CO2 data and NASA tin caloric curves retain their separate
evaluations. No reference family is silently relabelled.

The standard-state cassiterite row gives S° = 49.0 with magnitude 0.2. Its
separate fitted high-temperature table on printed page 219 gives 49.02 at
298.15 K. This pack selects the original standard-state row consistently
instead of combining a fitted ordinate with the other table's uncertainty.

## Three explicit equations

| CLI mode | Equation | ΔrH°, J/mol equation | ΔrG°, J/mol equation |
| --- | --- | ---: | ---: |
| `cassiterite-co-reduction` | SnO2(s) + 2 CO(g) → Sn(s) + 2 CO2(g) | 11600 | 1200 |
| `graphite-gasification` | C(s) + CO2(g) → 2 CO(g) | 172500 | 120200 |
| `cassiterite-graphite-reduction` | SnO2(s) + 2 C(s) → Sn(s) + 2 CO(g) | 356600 | 241600 |

The enum variants are `CassiteriteCarbonMonoxideReduction`,
`GraphiteCarbonDioxideGasification` and `CassiteriteGraphiteReduction`.
The [primary graphite-reduction study](https://doi.org/10.1016/0040-6031(88)87012-6)
distinguishes the SnO2/Sn gas equilibrium and Boudouard carbon gasification.
The library equations provide aggregate reference endpoints; they do not
import reaction rates or assume these are elementary mechanisms.

Adding the first equation and twice the second cancels internal CO/CO2 to
recover the third. Their independent source enthalpy and Gibbs sums satisfy
the same Hess identity. The separately evaluated H − TΔS route is retained
and differs slightly from the direct Gibbs tabulation because the original
values are rounded independently. No property is adjusted to erase that gap.

`ReactionNetwork` admits gross consumption against each stage's original
stocks. Newly produced CO can fund a subsequent reduction call; it cannot
fund another requested reaction in the same simultaneous call. An aggregate
carbon endpoint also needs two available moles of graphite per mole of SnO2.
An insufficient inventory refuses the whole stage without clipping recovery.

## Run the source-backed calculation

```bash
tin_pack_dir=$(mktemp -d)
for pack in cassiterite-solid-usgs-b2131-thermochemistry tin-white-solid-usgs-b2131 graphite-solid-usgs-b2131 carbon-monoxide-gas-usgs-b2131 carbon-dioxide-gas-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$tin_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction cassiterite-co-reduction --equilibrium --gibbs \
  "$tin_pack_dir/cassiterite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$tin_pack_dir/carbon-monoxide-gas-usgs-b2131.fsmatpk" \
  "$tin_pack_dir/tin-white-solid-usgs-b2131.fsmatpk" \
  "$tin_pack_dir/carbon-dioxide-gas-usgs-b2131.fsmatpk" 0.150709

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction cassiterite-graphite-reduction --gibbs \
  "$tin_pack_dir/cassiterite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$tin_pack_dir/graphite-solid-usgs-b2131.fsmatpk" \
  "$tin_pack_dir/tin-white-solid-usgs-b2131.fsmatpk" \
  "$tin_pack_dir/carbon-monoxide-gas-usgs-b2131.fsmatpk" 0.150709
```

Each example consumes one mole (0.150709 kg) of pure cassiterite and forms
one mole (0.118710 kg) of white solid tin. The CO route consumes 0.056020 kg
CO and forms 0.088020 kg CO2. The graphite route consumes 0.024022 kg graphite
and forms 0.056020 kg CO. Rounded source molar masses leave respective
formed-minus-consumed residuals of +0.000001 kg and -0.000001 kg per written
equation; exact elemental and charge conservation is unaffected.

The focused regression is `cargo test -p xtask --test mining_tin_sources`.
It compiles the real source packs, verifies all 31 literal claims and canonical
pack identities, checks unsupported states and receipt replay, verifies the
Hess cycle and source-mass residuals, and exercises finite CO/graphite stocks
across simultaneous and sequential reaction stages.

These point references supply neither a furnace operating temperature nor
sensible heat, melting duty, reaction kinetics, gas-transfer rate, equilibrium
conversion, coke quality, SnO volatilization, Fe-Sn hardhead losses or slag
partitioning. Furnace conditions require their own supported phases and
temperature-dependent properties; a 298.15 K reference heat is not a smelting
energy requirement or a recovered-metal forecast.
