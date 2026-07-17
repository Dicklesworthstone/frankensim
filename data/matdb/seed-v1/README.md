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

Each directory seeds one immutable species association with ideal-gas phase
and EOS, reference pressure `100 kPa` (the source report's `1 bar` gas
standard state), and elemental-reference convention
`NASA-TP-2002-211556-reference-elements-298.15K-1bar`.

The primary source is McBride, Zehe, and Gordon, *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*,
NASA/TP-2002-211556 (2002), NTRS document `20020085330`. Appendix B reports
the five gas molecular weights, while the Standard States section defines the
ideal-gas standard pressure as `1 bar`. The NTRS record marks the report
publicly distributable and as a work of the U.S. Government whose public use
is permitted. These seeds copy only the factual associations above, retain
NASA attribution, and do not copy third-party figures or tables.

As independent spot checks, the NIST Chemistry WebBook SRD 69 pages report
`16.0425` for methane, `28.0134` for nitrogen, `31.9988` for oxygen, `39.948`
for argon, and `44.0095` for carbon dioxide. Each displayed value agrees with
the corresponding NASA value within one half-unit at NIST's displayed
precision. NIST values are comparison oracles only; they are not pack sources
and do not replace the retained NASA values.

Primary and comparison references:

- <https://ntrs.nasa.gov/citations/20020085330>
- <https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf>
- <https://sti.nasa.gov/disclaimers/>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440371>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389>

To compile the source into a canonical runtime pack:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/methane/manifest.tsv \
  --out /path/to/CH4.fsspcpk
```

## No-claim boundary

These tranches are species identity/standard-state associations, not complete
material cards. In particular, the constituent associations do not specify an
air or exhaust mixture, composition basis, humidity, combustion state, or
validity domain. They supply no heat-capacity coefficients, equation evaluator,
uncertainty model for thermodynamic properties, reaction mechanism, equilibrium
result, or transport data. A decimal agreement check is not an uncertainty
estimate. Those claims require later, separately sourced seed records and keep
bead `frankensim-ext-matdb-seed-dataset-1sxe` open.
