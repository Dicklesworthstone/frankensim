# Reactive mineral thermochemistry

The five solid-phase profiles below add 18 claims through the existing `matdb-pack` compiler
and catalog. They complement mineral density and hardness references with
phase-specific molar quantities needed for thermal and chemical accounting.

| Profile | Molar mass, kg/mol | Formation enthalpy at 298.15 K, J/mol | Standard entropy, J/(mol K) | Cp domain |
| --- | ---: | ---: | ---: | --- |
| `calcite-solid-usgs-b2131-thermochemistry` | 0.100087 | −1207400 | 91.71 | 298.15–500 K |
| `pyrite-solid-janaf-thermochemistry` | 0.119967 | −171544 | 52.916 | 298.15–500 K |
| `quartz-alpha-janaf-standard-reference` | 0.0600843 | −910857 | 41.463 | No new Cp claim |
| `gypsum-solid-usgs-b2131-thermochemistry` | 0.172172 | −2023000 | 193.8 | Not supplied |
| `anhydrite-solid-usgs-b2131-thermochemistry` | 0.136142 | −1434400 | 107.4 | 298.15–1000 K |

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

Gypsum (CaSO4·2H2O) and orthorhombic anhydrite (CaSO4) follow the same USGS
reference convention. Their standard quantities come from
[printed page 28](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/34/).
Anhydrite's [page 285 table](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/291/)
supplies nine literal molar Cp entries: 101.23, 101.57, 116.68, 127.18, 135.03,
141.54, 147.55, 153.59 and 160.04 J/(mol K) at 298.15, 300, 400, 500, 600,
700, 800, 900 and 1000 K. These become mass-specific knots using the same
source formula weight; linear interpolation is an authored approximation.
Gypsum receives no Cp curve because its selected reference table does not
provide one. Wallboard properties are not substituted for pure hydrated mineral
thermochemistry. Printed uncertainty magnitudes remain source notes without
an invented probability distribution or a hard-bound interpretation.

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

All solid-phase profiles above require typed temperature and pressure, `source-phase-solid=1`,
`source-pure-species=1`, `source-standard-state=1`,
`engineering-reference-model=1`, and `source-specimen-qualified=0`. The
additional phase selectors are `source-phase-calcite=1`,
`source-phase-pyrite=1`, `source-phase-quartz-alpha=1`,
`source-phase-gypsum=1`, or `source-phase-anhydrite=1` respectively.
Each profile's own elemental-reference selector is mandatory.

These solid-phase records do not supply aqueous-ion properties, dissolution/oxidation
rates, decomposition boundaries, equilibrium conversion or plant recovery.
Reaction heat needs a complete, compatible set of participant formation data;
an isolated mineral enthalpy is not a heat-of-leaching prediction.

## Executable checks

```bash
cargo test -p xtask --test mining_mineral_thermochemistry
cargo test -p xtask --test mining_sulfate_thermochemistry
```

The six tests compile the actual manifests, verify canonical packs and query
receipts, check molar units and printed source values, and refuse incompatible
phases, conventions, pressures or temperatures. Heating tests integrate the
compiled Cp curves and compare against independently printed enthalpy columns.
For 1.5 kg heated from 300 to 400 K, the piecewise-linear estimates are
135497.12 J for calcite and 81983.17 J for pyrite. Their differences from the
printed heat-content changes are about 0.823% and 0.503%, respectively. These
regression tolerances describe coarse interpolation, not certified uncertainty.
The calcite check correctly interprets its `(H(T)-H(298.15))/T` source column.

For 2 kg of anhydrite heated from 300 to 400 K, the compiled curve gives
160310.56 J, differing by approximately 0.424% from the independently printed
heat-content change. This comparison uses the same divided-by-temperature
column convention. All nine Cp knots and the exact 1000 K endpoint are checked;
1000.01 K refuses. Formation enthalpy and entropy still apply only at 298.15 K.
The new sulfate records alone do not predict scale growth, hydrate stability,
dehydration rate or a complete reaction heat without compatible participants.

## Aqueous references for calcite neutralization

Three additional profiles supply 13 claims from the same USGS Bulletin 2131
reference convention as calcite. The primary rows are on printed pages 7, 6,
and 15 of the [report](https://pubs.usgs.gov/bul/2131/report.pdf), respectively.
All values below apply at 298.15 K and 100000 Pa.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Conventional ionic entropy, J/(mol K) |
| --- | ---: | ---: | ---: |
| `hydrogen-ion-aqueous-usgs-b2131` | 0.001008 | 0 | 0 |
| `calcium-ion-aqueous-usgs-b2131` | 0.040078 | -543000 | -56.2 |
| `bicarbonate-ion-aqueous-usgs-b2131` | 0.061017 | -689900 | 98.4 |

The source defines aqueous `m=1` as the hypothetical ideal state at unit
activity with ions at infinite dilution. It is not a measured finite-
concentration solution. Claims require `source-phase-aqueous=1`,
`source-aqueous-unit-activity-infinite-dilution=1`,
`source-ionic-reference-hydrogen-zero=1`, `source-individual-ion=1`, and
`source-elemental-reference-usgs-b2131=1`, in addition to the typed reference
temperature, pressure, standard-state and engineering-reference selectors.
Individual ions do not carry a `source-pure-species` assertion.

Entropy is named `conventional-standard-molar-entropy`, distinctly from the
solid-phase `standard-molar-entropy` property. The H+ zeros define the ionic
reference; calcium's negative conventional entropy is retained. Calcium and
bicarbonate also expose the source's reported uncertainty magnitudes as
dimensioned metadata: 1000 and 200 J/mol for formation enthalpy, and 1.0 and
0.5 J/(mol K) for conventional entropy. The compilation's calorimetric
convention is twice the standard deviation of the mean; it supplies no
universal confidence probability or covariance for these rows. Main claims
therefore retain `Unstated` uncertainty. No certified reaction interval is
invented from these magnitudes.

Together with the existing calcite pack, these sources close the declared
bicarbonate endpoint:

`CaCO3(calcite) + H+(aq) -> Ca[2+](aq) + HCO3[-](aq)`

The existing exact atom/charge balance and source-backed formation-enthalpy
consumer give **-25500 J/mol of reaction extent**. For 100.087 kg of actually
reacted pure calcite, this is 1000 mol extent and **-25.5 MJ**. It consumes
1.008 kg H+ and forms 40.078 kg calcium ions and 61.017 kg bicarbonate ions.
These are net ionic quantities, not the mass of a commercial acid solution.

Compile each profile using the existing `cargo run -p xtask -- matdb-pack
--manifest .../manifest.tsv --out ...fsmatpk` command, then run:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction calcite-bicarbonate-neutralization \
  /tmp/mining-materials/calcite-solid-usgs-b2131-thermochemistry.fsmatpk \
  /tmp/mining-materials/hydrogen-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/calcium-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/bicarbonate-ion-aqueous-usgs-b2131.fsmatpk 100.087
cargo test -p xtask --test mining_aqueous_sources
```

The integration checks compile all four actual packs, verify source receipts,
preserve the ionic convention, and exercise an acid-limited inventory: 10 mol
calcite plus 6 mol H+ permits at most 6 mol extent, leaving 4 mol calcite and
giving -153000 J at that ceiling. Overdrawn inventory refuses. This is an
availability bound, not a predicted conversion or dissolution rate. The
two-proton dissolved-CO2 endpoint remains a different reaction requiring its
own complete source set; bicarbonate data are not relabeled as CO2(aq).
Finite-concentration activities, heat of mixing, pH selection, gas transfer
and thermal plant behavior require additional models.
