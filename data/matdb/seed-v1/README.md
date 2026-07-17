# FrankenSim material seed data v1

This directory retains the human-reviewable source inputs compiled by
`xtask matdb-pack`. Generated binary packs are deliberately not committed:
the compiler must reproduce them from the pinned manifest and source record.

## Current species-association tranches

| Directory | Species | NASA molecular weight | Program role |
| --- | --- | ---: | --- |
| `methane/` | `CH4` | `16.04246 g/mol` | fuel constituent |
| `nitrogen/` | `N2` | `28.01340 g/mol` | dry-air/exhaust constituent |
| `oxygen/` | `O2` | `31.99880 g/mol` | dry-air/exhaust constituent |
| `argon/` | `Ar` | `39.94800 g/mol` | dry-air/exhaust constituent |
| `carbon-dioxide/` | `CO2` | `44.00950 g/mol` | exhaust constituent |
| `water-vapor/` | `H2O` | `18.01528 g/mol` | exhaust/humidity constituent |
| `carbon-monoxide/` | `CO` | `28.01010 g/mol` | exhaust constituent |

Each directory seeds one immutable species association with ideal-gas phase
and EOS, reference pressure `100 kPa` (the source report's `1 bar` gas
standard state), and elemental-reference convention
`NASA-TP-2002-211556-reference-elements-298.15K-1bar`.

The primary source is McBride, Zehe, and Gordon, *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*,
NASA/TP-2002-211556 (2002), NTRS document `20020085330`. Appendix B reports
the seven gas molecular weights, while the Standard States section defines the
ideal-gas standard pressure as `1 bar`. The NTRS record marks the report
publicly distributable and as a work of the U.S. Government whose public use
is permitted. These seeds copy only the factual associations above, retain
NASA attribution, and do not copy third-party figures or tables.

As independent spot checks, the NIST Chemistry WebBook SRD 69 pages report
`16.0425` for methane, `28.0134` for nitrogen, `31.9988` for oxygen, `39.948`
for argon, `44.0095` for carbon dioxide, `18.0153` for water, and `28.0101`
for carbon monoxide. Each displayed value agrees with the corresponding NASA
value within one half-unit at NIST's displayed precision. NIST values are
comparison oracles only; they are not pack sources and do not replace the
retained NASA values.

Primary and comparison references:

- <https://ntrs.nasa.gov/citations/20020085330>
- <https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf>
- <https://sti.nasa.gov/disclaimers/>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080>

## Aluminum 6061-T6 cryogenic exact-point tranche

`aluminum-6061-t6-cryogenic/` is the first committed bulk-material tranche.
It identifies the condition as Aluminum 6061-T6, UNS AA96061, and retains
thermal conductivity, specific heat capacity, and Young's modulus at exactly
`77 K` and `293 K`. Each value is an evaluation of the polynomial printed on
NIST's Aluminum 6061-T6 material-property page, using the page's displayed
coefficients and equation. The runtime claims have degenerate temperature
validity intervals (`T_min = T_max`), so they do not imply interpolation or a
continuous curve.

NIST reports curve-fit errors relative to the underlying data of `0.5%` for
thermal conductivity, `5%` for specific heat, and `1%` for Young's modulus.
Those figures do not state a confidence level or degrees of freedom, so the
source records retain them as observation caveats and explicitly encode the
runtime uncertainty as `Unstated`. They are not laundered into statistical
confidence intervals.

As independent G3 comparison evidence, the 1966 NASA thermophysical-property
compilation tabulates 6061-T6 thermal conductivity as `82 W/(m K)` at `75 K`
and `155 W/(m K)` at `300 K`. Those nearby-temperature values differ from the
NIST-derived `77 K` and `293 K` points by less than `3%`; they are comparison
oracles only and do not overwrite the NIST-derived claims. NASA's NTRS record
for its 1969 Aluminum 6061 handbook independently confirms that the handbook
covers cryogenic, ambient, and elevated-temperature properties, but this
tranche does not silently import values from that separate work.

The NIST page is not marked as copyrighted. NIST's Copyrights & Disclaimers
page says unmarked NIST-site information is public information that may be
distributed or copied and requests appropriate credit. The manifest records
that specific redistribution decision and retains NIST attribution; it does
not generalize the decision to separately licensed NIST research products.

Material references:

- <https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061>
- <https://www.nist.gov/copyrights-disclaimers>
- <https://ntrs.nasa.gov/api/citations/19660014513/downloads/19660014513.pdf>
- <https://ntrs.nasa.gov/citations/19690000065>

To compile the sources into canonical runtime packs:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/methane/manifest.tsv \
  --out /path/to/CH4.fsspcpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aluminum-6061-t6-cryogenic/manifest.tsv \
  --out /path/to/aluminum-6061-t6-cryogenic.fsmatpk
```

## No-claim boundary

The gas tranches are species identity/standard-state associations, not complete
material cards. They do not specify an air or exhaust mixture, composition
basis, humidity level, combustion state, or validity domain. In particular,
the H2O association does not select a wet-air composition and the CO
association does not claim incomplete combustion. They supply no heat-capacity
coefficients, equation evaluator, uncertainty model for thermodynamic
properties, reaction mechanism, equilibrium result, or transport data. A
decimal agreement check is not an uncertainty estimate.

The Aluminum 6061-T6 tranche is not a general-purpose design card. It contains
six polynomial-derived scalar points for one named temper, not density,
composition tolerances, anisotropy, batch/process variation, plasticity,
strength, fatigue, fracture, corrosion, joining, or a continuous-temperature
constitutive model. NIST's curve-fit error is retained but not represented as a
confidence interval. The NASA comparison uses nearby temperatures and is only
a coarse independent agreement check. These missing claims and the remaining
named materials and interface systems keep bead
`frankensim-ext-matdb-seed-dataset-1sxe` open.
