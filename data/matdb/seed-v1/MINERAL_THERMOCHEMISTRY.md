# Reactive mineral thermochemistry

These three profiles add 11 claims through the existing `matdb-pack` compiler
and catalog. They complement mineral density and hardness references with
phase-specific molar quantities needed for thermal and chemical accounting.

| Profile | Molar mass, kg/mol | Formation enthalpy at 298.15 K, J/mol | Standard entropy, J/(mol K) | Cp domain |
| --- | ---: | ---: | ---: | --- |
| `calcite-solid-usgs-b2131-thermochemistry` | 0.100087 | −1207400 | 91.71 | 298.15–500 K |
| `pyrite-solid-janaf-thermochemistry` | 0.119967 | −171544 | 52.916 | 298.15–500 K |
| `quartz-alpha-janaf-standard-reference` | 0.0600843 | −910857 | 41.463 | No new Cp claim |

Formation enthalpy and entropy are point claims at 298.15 K and 100000 Pa.
Calcite and pyrite retain four literal molar Cp observations, converted to
J/(kg K) with the selected source's formula weight. Their piecewise-linear
interpolation is an authored approximation with unstated uncertainty. Quartz
already has a NIST alpha-quartz caloric profile; this increment supplies only
its missing standard reference quantities instead of duplicating that curve.

## Source identities and conventions

Calcite follows [USGS Bulletin 2131](https://doi.org/10.3133/b2131),
[printed page 268](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/274/).
Its reference convention is `source-elemental-reference-usgs-b2131=1`.
The source specifies one bar and standard elemental states, while its elemental
evaluations combine JANAF, CODATA and separate USGS calcium measurements.
This profile is not relabeled as a JANAF evaluation for convenience.

Pyrite follows [JANAF Fe-028](https://janaf.nist.gov/tables/Fe-028.html),
cross-checked against printed page 1199 of the
[Fourth Edition iron tables](https://janaf.nist.gov/pdf/JANAF-FourthEd-1998-Iron.pdf).
The selected formation enthalpy is −171.544 kJ/mol. The WebBook's alternative
−167.36 kJ/mol is documented in the source notes and is not averaged into this
record. The selected table's formula weight, 119.967 g/mol, is used consistently
for its mass and Cp normalization. Printed uncertainty magnitudes remain in the
source notes without being converted into an invented probability distribution.

Quartz follows [JANAF O-037](https://janaf.nist.gov/tables/O-037.html) for
enthalpy and entropy, with explicitly attributed
[WebBook molar mass](https://webbook.nist.gov/cgi/cbook.cgi?ID=C14808607&Mask=2).
Both JANAF profiles require `source-elemental-reference-janaf=1`.

All profiles require typed temperature and pressure, `source-phase-solid=1`,
`source-pure-species=1`, `source-standard-state=1`,
`engineering-reference-model=1`, and `source-specimen-qualified=0`. The
additional polymorph selectors are `source-phase-calcite=1`,
`source-phase-pyrite=1`, or `source-phase-quartz-alpha=1` respectively.
Each profile's own elemental-reference selector is mandatory.

These records do not supply aqueous-ion properties, dissolution/oxidation
rates, decomposition boundaries, equilibrium conversion or plant recovery.
Reaction heat needs a complete, compatible set of participant formation data;
an isolated mineral enthalpy is not a heat-of-leaching prediction.

## Executable checks

```bash
cargo test -p xtask --test mining_mineral_thermochemistry
```

The three tests compile the actual manifests, verify canonical packs and query
receipts, check molar units and printed source values, and refuse incompatible
phases, conventions, pressures or temperatures. Heating tests integrate the
compiled Cp curves and compare against independently printed enthalpy columns.
For 1.5 kg heated from 300 to 400 K, the piecewise-linear estimates are
135497.12 J for calcite and 81983.17 J for pyrite. Their differences from the
printed heat-content changes are about 0.823% and 0.503%, respectively. These
regression tolerances describe coarse interpolation, not certified uncertainty.
The calcite check correctly interprets its `(H(T)-H(298.15))/T` source column.
