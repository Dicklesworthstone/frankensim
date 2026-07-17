# FrankenSim material seed data v1

This directory retains the human-reviewable source inputs compiled by
`xtask matdb-pack`. Generated binary packs are deliberately not committed:
the compiler must reproduce them from the pinned manifest and source record.

## Current tranche

`methane/` seeds one immutable `CH4` species association:

- molar mass: `16.04246 g/mol`;
- standard-state phase and EOS: ideal gas;
- reference pressure: `100 kPa` (the source report's `1 bar` gas standard
  state); and
- elemental-reference convention:
  `NASA-TP-2002-211556-reference-elements-298.15K-1bar`.

The primary source is McBride, Zehe, and Gordon, *NASA Glenn Coefficients for
Calculating Thermodynamic Properties of Individual Species*,
NASA/TP-2002-211556 (2002), NTRS document `20020085330`. Appendix B reports
the `CH4` gas molecular weight, while the Standard States section defines the
ideal-gas standard pressure as `1 bar`. The NTRS record marks the report
publicly distributable and as a work of the U.S. Government whose public use
is permitted. This seed copies only the factual association above, retains
NASA attribution, and does not copy third-party figures or tables.

As an independent spot check, the NIST Chemistry WebBook SRD 69 methane page
reports molecular weight `16.0425`. That displayed value agrees with the NASA
value within `0.00005 g/mol`, one half-unit at NIST's displayed precision. The
NIST value is a comparison oracle only; it is not the pack's source and is not
substituted for the retained NASA value.

Primary and comparison references:

- <https://ntrs.nasa.gov/citations/20020085330>
- <https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf>
- <https://sti.nasa.gov/disclaimers/>
- <https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828>

To compile the source into a canonical runtime pack:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/methane/manifest.tsv \
  --out /path/to/CH4.fsspcpk
```

## No-claim boundary

This first tranche is a species identity/standard-state association, not a
complete methane material card. It supplies no heat-capacity coefficients,
equation evaluator, uncertainty model for thermodynamic properties, validity
interval for such properties, reaction mechanism, equilibrium result, or
transport data. The decimal agreement check is not an uncertainty estimate.
Those claims require later, separately sourced seed records and keep bead
`frankensim-ext-matdb-seed-dataset-1sxe` open.
