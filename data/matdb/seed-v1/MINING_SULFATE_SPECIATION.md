# Sulfate and bisulfate in acidic process solutions

This source closure adds the individual aqueous ions **HSO4−** and **SO4²−**
from the same NBS evaluation as the existing proton, hydroxide and liquid-water
packs. The existing coupled acid/base solver can now calculate their shared
proton activity, free-ion partition and water autoprotolysis from an explicit
analytical total and spectator-ion inventory. This supports the acid/base part
of dilute sulfate leach-solution and neutralization calculations.

## Original source values and conventions

The primary numerical source is the [NIST corrected CSV of the NBS
tables](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv),
[dataset DOI 10.18434/M32124](https://doi.org/10.18434/M32124),
Table 14:S (prepared 1964), the `HSO4<sup>-</sup>` and
`SO4<sup>2-</sup>` rows with state `ao`. The CSV header fixes the state at
298.15 K and 0.1 MPa (1 bar). These rows were read directly from the official
dataset. The original publication is Wagman et al. (1982), *The NBS Tables of
Chemical Thermodynamic Properties*, JPCRD 11 Supplement 2; the
[digitization paper](https://doi.org/10.6028/jres.125.007) describes the
machine-readable edition with 1989 errata.

| Profile | Molar mass, kg/mol | Formation H, J/mol | Direct formation G, J/mol | Conventional ionic S, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `hydrogen-sulfate-ion-aqueous-nbs1989` | 0.0970696 | −887340 | −755910 | 131.8 |
| `sulfate-ion-aqueous-nbs1989` | 0.0960616 | −909270 | −744530 | 20.1 |

Every claim requires `source-elemental-reference-nbs1989`, the aqueous phase,
the hypothetical unit-molality/infinite-dilution standard state, the individual
ion identity and the H+(aq)=0 convention. The pressure is **100000 Pa**.
The individual-ion entropy is conventional, not an absolute species entropy.
No confidence distribution or cross-species covariance is invented;
uncertainty remains `Unstated`.

The pre-existing `sulfate-ion-aqueous-usgs-b2131` remains a separate evaluation.
It is not silently substituted into this NBS source set, even though both
name the sulfate ion and use a 1 bar reference pressure.

## Conserved acid step and compatible thermochemistry

`MineralReaction::HydrogenSulfateDissociation` provides the exact conserved
equation:

```text
HSO4−(aq) -> H+(aq) + SO4²−(aq)
```

One written reaction extent removes one mole of hydrogen sulfate and forms one
mole each of proton and sulfate. The integer equation conserves H, S, O and
charge. It is the hydrogen-sulfate acid step; its existence does not select
the molecular-H2SO4 fraction in a real acid solution.

With the original matched NBS sources:

| Quantity per mole of the written equation | Value |
| --- | ---: |
| Standard reaction H | −21930 J/mol |
| Standard reaction S | −111.7 J/(mol K) |
| Direct standard reaction G | +11380 J/mol |
| Reaction G derived from H − T S | +11373.355 J/mol |
| Ka from the direct G values | 0.0101463597511199 |
| pKa from the direct G values | 1.993689743323917 |

The directly tabulated and caloric routes remain distinct. Their 6.645 J/mol
difference is retained; source data are not adjusted to make the two routes
identical. Both routes retain the original source receipts.

## Coupled ideal solution calculation

The analytical pool is **HSO4− + SO4²−**, excluding metal-sulfate complexes,
gypsum or other solids, sulfide, sulfite and molecular H2SO4. The caller supplies
the total pool amount, all fixed spectator ions, and one amount at unit relative
activity for the fixed solvent basis. For the demonstrated 10 kg solvent basis,
that amount is 10 mol under the unit-molality convention; it is not inferred
from solution volume.

The solver uses `a = n/n_ref`, activity coefficients of one, and unit water
activity. At one shared proton activity, it solves sulfate-family mass action,
water autoprotolysis and the full declared charge balance together. A sulfate
family with no positive spectator inventory acquires its balancing charge from
the free proton pool. A bisulfate salt requires one fixed positive charge per
sulfur; a sulfate salt requires two. Added positive spectator equivalents
represent the separately declared effect of strong-base neutralization; the
solver does not calculate reagent purity, salt dissolution or commercial dose.

The executable prints separate acid and base amounts, pH, free proton and
hydroxide amounts, family-total residual, charge residual and logarithmic
mass-action residual:

```bash
sulfate_pack_dir=/path/to/compiled-materials
mkdir -p "$sulfate_pack_dir"
for profile in hydrogen-sulfate-ion-aqueous-nbs1989 \
  sulfate-ion-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$profile/manifest.tsv" \
    --out "$sulfate_pack_dir/$profile.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_speciation -- \
  --materials "$sulfate_pack_dir" --family sulfate
```

The existing cyanide example remains available with `--family cyanide` or by
omitting the family option. The sulfate mode demonstrates pure water, an acidic
uncomplexed sulfate pool, the one- and two-counterion salt limits, added-base
equivalents and acidification by explicitly declared spectator anions.

## Focused checks and limits

`cargo test -p xtask --test mining_sulfate_speciation` compiles the actual TSVs
and verifies the literal source values, direct/caloric reaction sums, source
receipts and pKa. The acidic-pool and bisulfate-salt cases are compared with
independent quadratic roots with an explicit bound for the water correction.
Other checks retain total sulfur and charge, verify the response to declared
neutralization and acidification, preserve pH under consistent amount/basis
scaling, and reject missing sources, out-of-reference temperature and a
caller attempt to relabel the USGS sulfate evaluation as NBS.

This is a reference-temperature, ideal, uncomplexed acid/base model. It does
not supply a concentrated-electrolyte activity model, molecular sulfuric-acid
dissociation, metal complexation, mineral precipitation, temperature evolution,
heat of dilution, corrosion, reaction kinetics or a measured plant pH. The
fixed-solvent calculation conserves the analytical family and charge; it does
not assert a closed finite-solvent hydrogen/oxygen balance.
