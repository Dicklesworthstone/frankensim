# Mining mineral reference data

Five source profiles add 19 claims to the existing material compiler and catalog.
They extend the existing pyrite, chalcopyrite, chalcocite, cassiterite, dolomite
and magnetite reference profiles with bornite, sphalerite, galena, arsenopyrite
and malachite. Each addition retains the publisher's mineral identity and
measured-density role. It supplies reference data for explicitly selected grains;
it does not characterize a mine's ore.

| Bundle prefix | Formula / phase | Source D(meas.), g/cm³ | Mohs endpoints |
| --- | --- | --- | --- |
| `bornite` | Cu₅FeS₄, orthorhombic | 5.06–5.08 | 3–3.25 |
| `sphalerite` | (Zn,Fe)S, cubic | 3.9–4.1 | 3.5–4 |
| `galena` | PbS, cubic | 7.58 | 2.5–2.75 |
| `arsenopyrite` | FeAsS, monoclinic | 6.07(15) | 5.5–6 |
| `malachite` | Cu₂(CO₃)(OH)₂, monoclinic | 4.05(2) | 3.5–4 |

Append `-handbook-mineralogy` to the prefix for the pack ID. Each manifest cites
the corresponding original *Handbook of Mineralogy* sheet, Physical Properties,
version 1 (2001–2005), hosted by the Mineralogical Society of America. It also
pins the downloaded PDF's SHA-256. The publisher's
[page guide](https://handbookofmineralogy.org/) identifies the density unit as
g/cm³. Source PDFs, illustrations, descriptive text and table layouts are not
redistributed; the bundles contain selected attributed facts and original notes.

For bornite and sphalerite, `density-reported-range-minimum` and
`density-reported-range-maximum` retain the source endpoints using the existing
mineral-range vocabulary. Neither profile supplies a fabricated point `density`.
For the two parenthetical values, `density` carries the printed central value
and `density-source-parenthetical-spread` carries 0.15 or 0.02 g/cm³. No
confidence level, distribution or prediction interval is supplied. All 19 claims
therefore retain `UncertaintyModel::Unstated`. Calculated crystallographic
densities are not silently substituted for the measured values.

Mohs values use distinct dimensionless endpoint properties:
`hardness-mohs-source-minimum` and `hardness-mohs-source-maximum`. Mohs hardness
is an ordinal scratch scale. It supplies no conversion to Vickers hardness,
abrasion index, fracture energy or Bond work index. Mineral mixing and ore
grindability require their own measured inputs and model assumptions.

Every query must explicitly select all six source conditions:

- `source-mineral-<name>=1`;
- `mineral-reference-model=1`;
- `source-test-temperature-known=0`;
- `source-test-pressure-known=0`;
- `bulk-ore-qualified=0`;
- `source-specimen-qualified=0`.

The source rows do not state a test temperature, pressure or common specimen.
Those missing conditions do not become room-temperature measurements or a
temperature response. Mineral formula labels identify the reported phase; they
do not promise end-member purity. Sphalerite explicitly retains iron substitution.

For example, the lower bornite mineral-density endpoint is available through
the ordinary catalog:

```bash
cargo run -p xtask -- matdb-catalog \
  --pack bornite-handbook-mineralogy \
  --property density-reported-range-minimum \
  --at source-mineral-bornite=1 --at mineral-reference-model=1 \
  --at source-test-temperature-known=0 --at source-test-pressure-known=0 \
  --at bulk-ore-qualified=0 --at source-specimen-qualified=0
```

The result is 5060 kg/m³ with the source receipt. It is not a bulk density for a
copper-bearing rock or a recovery prediction. The focused integration test
`cargo test -p xtask --test mining_mineral_coverage` compiles all five manifests,
round-trips their canonical packs, queries every source value, and checks the
missing-context, incompatible-context, wrong-unit and unsupported-property cases.
