# Reactive mineral thermochemistry

The solid-phase profiles below use the existing `matdb-pack` compiler and
catalog. They complement mineral density and hardness references with
phase-specific molar quantities needed for thermal and chemical accounting.
Calcite also has a directly tabulated standard formation Gibbs energy for
reaction driving-force calculations with the compatible aqueous source set.

| Profile | Molar mass, kg/mol | Formation enthalpy at 298.15 K, J/mol | Formation Gibbs energy at 298.15 K, J/mol | Standard entropy, J/(mol K) | Cp domain |
| --- | ---: | ---: | ---: | ---: | --- |
| `calcite-solid-usgs-b2131-thermochemistry` | 0.100087 | −1207400 | −1128500 | 91.71 | 298.15–500 K |
| `dolomite-solid-usgs-b2131-thermochemistry` | 0.184401 | −2324500 | Not supplied | 155.2 | 298.15–500 K |
| `pyrite-solid-janaf-thermochemistry` | 0.119967 | −171544 | Not supplied | 52.916 | 298.15–500 K |
| `quartz-alpha-janaf-standard-reference` | 0.0600843 | −910857 | Not supplied | 41.463 | No new Cp claim |
| `gypsum-solid-usgs-b2131-thermochemistry` | 0.172172 | −2023000 | Not supplied | 193.8 | Not supplied |
| `anhydrite-solid-usgs-b2131-thermochemistry` | 0.136142 | −1434400 | Not supplied | 107.4 | 298.15–1000 K |

Formation enthalpy, formation Gibbs energy and entropy are point claims at
298.15 K and 100000 Pa. A Cp curve does not extend their reference-state domain.
Calcite and pyrite retain four literal molar Cp observations, converted to
J/(kg K) with the selected source's formula weight. Their piecewise-linear
interpolation is an authored approximation with unstated uncertainty. Quartz
already has a NIST alpha-quartz caloric profile; this increment supplies only
its missing standard reference quantities instead of duplicating that curve.
Dolomite similarly retains four Cp observations and separately dimensioned
metadata for the source's reported enthalpy and entropy uncertainty magnitudes.

## Source identities and conventions

Calcite follows [USGS Bulletin 2131](https://doi.org/10.3133/b2131),
[printed page 268](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/274/).
Its reference convention is `source-elemental-reference-usgs-b2131=1`.
The source specifies one bar and standard elemental states, while its elemental
evaluations combine JANAF, CODATA and separate USGS calcium measurements.
This profile is not relabeled as a JANAF evaluation for convenience.
The selected page's literal formation Gibbs energy is **−1128.5 kJ/mol**;
nearby values from other compilations are not substituted. Its
`standard-molar-formation-gibbs-energy` claim has separate source evidence and
retains the same calcite phase and USGS reference selectors.

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

The separate `pyrite-solid-usgs-b2131-thermochemistry` profile supplies a
compatible USGS pyrite reference for aqueous redox accounting. Together with
O2(g), H2O(l), Fe2+(aq), the new Fe3+(aq) profile, sulfate and H+, it supplies
complete enthalpy and Gibbs sources for pyrite oxidation, ferrous oxidation
and ferric leaching. Their primary values, reference conventions, reaction
equations and executable activity example are in
[Aqueous pyrite oxidation and acid-generation references](MINING_ACID_GENERATION.md).

All solid-phase profiles above require typed temperature and pressure, `source-phase-solid=1`,
`source-pure-species=1`, `source-standard-state=1`,
`engineering-reference-model=1`, and `source-specimen-qualified=0`. The
additional phase selectors are `source-phase-calcite=1`,
`source-phase-dolomite=1`,
`source-phase-pyrite=1`, `source-phase-quartz-alpha=1`,
`source-phase-gypsum=1`, or `source-phase-anhydrite=1` respectively.
Each profile's own elemental-reference selector is mandatory.

The solid-phase records supply no aqueous-ion properties, dissolution/oxidation
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

Three aqueous profiles use the same USGS Bulletin 2131 reference convention
as calcite. Their primary rows are on printed pages 7, 6,
and 15 of the [report](https://pubs.usgs.gov/bul/2131/report.pdf), respectively.
All values below apply at 298.15 K and 100000 Pa.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Formation Gibbs energy, J/mol | Conventional ionic entropy, J/(mol K) |
| --- | ---: | ---: | ---: | ---: |
| `hydrogen-ion-aqueous-usgs-b2131` | 0.001008 | 0 | 0 | 0 |
| `calcium-ion-aqueous-usgs-b2131` | 0.040078 | -543000 | -553600 | -56.2 |
| `bicarbonate-ion-aqueous-usgs-b2131` | 0.061017 | -689900 | -586800 | 98.4 |

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

The independently tabulated Gibbs sources give **−11900 J/mol extent**,
`ln K° = 4.8004059` and dimensionless `K° ≈ 121.55975` at 298.15 K.
`ConservedReaction::standard_gibbs_change` retains each selected claim and
receipt, including the zero-reference H+ claim. With explicitly supplied
dimensionless activities, `StandardReactionGibbs::at_activities` uses

$$
Q = \frac{a_{\mathrm{Ca^{2+}}}a_{\mathrm{HCO_3^-}}}
         {a_{\mathrm{CaCO_3}}a_{\mathrm{H^+}}},
\qquad \Delta_r G = RT\ln(Q/K^\circ).
$$

This follows the IUPAC definitions of
[activity](https://goldbook.iupac.org/terms/view/A00115) and the
[standard equilibrium constant](https://goldbook.iupac.org/terms/view/S05915).
If pure calcite is assigned activity one, that value is still an explicit
input. `Q<K°` gives a forward thermodynamic driving force, and `Q>K°` a reverse
one. The calculation evaluates the declared state; it does not select pH,
compute solution speciation or solve for a finite equilibrium extent. Every
active species needs a finite positive activity in the retained source
convention. The returned `ΔrG` is a local derivative per mole of equation
extent, so a changing-activity process requires more than multiplying it by
reacted mass. The formation-Gibbs column remains separate from formation
enthalpy and conventional ionic entropy; no species-level `Hf−T*S` shortcut
or temperature extrapolation is used.

Compile each profile using the existing `cargo run -p xtask -- matdb-pack
--manifest .../manifest.tsv --out ...fsmatpk` command, then run:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction calcite-bicarbonate-neutralization --gibbs \
  /tmp/mining-materials/calcite-solid-usgs-b2131-thermochemistry.fsmatpk \
  /tmp/mining-materials/hydrogen-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/calcium-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/bicarbonate-ion-aqueous-usgs-b2131.fsmatpk 100.087
cargo test -p xtask --test mining_aqueous_sources
cargo test -p xtask --test mining_aqueous_redox_sources
```

The integration checks compile all four actual packs, verify source receipts,
preserve the ionic convention, and exercise an acid-limited inventory: 10 mol
calcite plus 6 mol H+ permits at most 6 mol extent, leaving 4 mol calcite and
giving -153000 J at that ceiling. Overdrawn inventory refuses. This is an
availability bound, not a predicted conversion or dissolution rate. The
two-proton dissolved-CO2 endpoint remains a different reaction requiring its
own complete source set; bicarbonate data are not relabeled as CO2(aq).
The CLI's `--gibbs` option reports the standard reaction Gibbs energy and
equilibrium constant when every pack provides the required source claim.
Repeated `--activity 'SPECIES=VALUE'` options before the pack list additionally
evaluate the supplied activity state; all participants, including calcite, are
required. A full source-pack example appears in the
[acid-generation reference](MINING_ACID_GENERATION.md#run-the-source-backed-balance).
Concentration-to-activity models, heat of mixing, pH selection, gas transfer
and thermal plant behavior remain separate inputs or models.

## Dolomite gangue neutralization and heating

`dolomite-solid-usgs-b2131-thermochemistry` and
`magnesium-ion-aqueous-usgs-b2131` add 11 claims from the same
[USGS Bulletin 2131](https://pubs.usgs.gov/bul/2131/report.pdf) convention.
The dolomite point quantities come from printed page 25; Mg²⁺ comes from page 8.
Dolomite's formula mass is 184.401 g/mol. The magnesium profile supplies
24.305 g/mol, formation enthalpy −467.0 kJ/mol, and conventional ionic entropy
−137.0 J/(mol K), all at 298.15 K and 100000 Pa. It requires the same aqueous
unit-activity, infinite-dilution and H⁺ zero selectors as calcium. Its negative
ionic entropy is not relabeled as an absolute entropy.

The source's reported uncertainty magnitudes are retained as dimensional
metadata: 1.5 kJ/mol and 0.3 J/(mol K) for dolomite, and 0.6 kJ/mol and
4.0 J/(mol K) for Mg²⁺. The claims retain `Unstated` confidence semantics.
No independent Gaussian uncertainties or propagated reaction intervals are
inferred from these magnitudes.

Together with the existing H⁺, Ca²⁺ and bicarbonate profiles, these records
complete the declared gangue acid-consumption endpoint:

`CaMg(CO3)2(dolomite) + 2 H+(aq) -> Ca[2+](aq) + Mg[2+](aq) + 2 HCO3[-](aq)`

The source-backed reaction enthalpy is **−65300 J/mol of reaction extent**.
For 184.401 kg of actually reacted pure dolomite, the consumer reports
**−65.3 MJ**, consumes 2.016 kg H⁺ and forms 40.078 kg Ca²⁺, 24.305 kg Mg²⁺
and 122.034 kg bicarbonate. For an inventory of 10 mol dolomite and 6 mol H⁺,
the existing exact stoichiometry allows at most 3 mol extent, leaves 7 mol
dolomite and gives −195900 J at that ceiling. This acid-limited amount is
an availability bound, not a predicted reaction rate or equilibrium endpoint.

Compile the five profiles with `xtask matdb-pack`, then run:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction dolomite-bicarbonate-neutralization \
  /tmp/mining-materials/dolomite-solid-usgs-b2131-thermochemistry.fsmatpk \
  /tmp/mining-materials/hydrogen-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/calcium-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/magnesium-ion-aqueous-usgs-b2131.fsmatpk \
  /tmp/mining-materials/bicarbonate-ion-aqueous-usgs-b2131.fsmatpk 184.401
cargo test -p xtask --test mining_dolomite_sources
```

Dolomite also supplies Cp at 298.15, 300, 400 and 500 K from printed page 269:
157.51, 158.06, 183.54 and 201.89 J/(mol K), normalized by 0.184401 kg/mol.
The authored linear interpolation gives 185248.45 J for 2 kg heated from
300 to 400 K. The independently printed heat-content column gives 186040.21 J,
a difference of about 0.426%. The 0.5% regression tolerance describes the
coarse interpolation; it is not certified uncertainty. Page 25's point
entropy of 155.2 is retained separately from page 269's more finely printed
155.18 J/(mol K). The Cp domain ends at 500 K; decomposition, dissolution
kinetics and off-reference reaction heat require separate models.
