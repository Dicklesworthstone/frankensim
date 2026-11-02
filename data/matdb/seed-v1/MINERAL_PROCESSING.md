# Mineral processing: public sources and callable physics

These **31 source profiles and 97 property claims** provide reusable mineral
and reagent inputs for FrankenMining's material-corpus and Crucible work.
They use the existing `matdb-pack` compiler, `NormalizedPack` queries and
`matdb-catalog` discovery. No separate material registry or schema is required.

The physical consumers are `fs_material::comminution`, `fs_material::slurry`,
`fs_material::recovery` and `fs_thermochem::stoichiometry`.
These are upstream building blocks: they do
not implement a FrankenMining flowsheet or complete its site-calibration
beads. In particular, FrankenMining's plan §13.10 keeps located ore tests
separate from public non-spatial material references, and §17.3 requires
declared validity for processing correlations.

## Ore and gangue constituents

The six `*-handbook-mineralogy-reference` profiles retain the Mineralogical
Society of America's [Handbook of Mineralogy](https://handbookofmineralogy.org/)
source identities and measured-density entries. Each manifest links the
individual mineral record. Formula metadata identify an ideal mineral class;
they are not a composite assay or a machine-authenticated elemental matrix.

| Profile prefix | Mineral formula | Density data, kg/m³ |
| --- | --- | ---: |
| `pyrite` | FeS2 | 5018 |
| `chalcopyrite` | CuFeS2 | Reported endpoints 4100 and 4300 |
| `chalcocite` | Cu2S | Reported endpoints 5500 and 5800 |
| `cassiterite` | SnO2 | Reported endpoints 6980 and 7010 |
| `dolomite` | CaMg(CO3)2 | 2860 |
| `magnetite` | Fe3O4 | 5175 |
| `bornite` | Cu5FeS4 | Reported endpoints 5060 and 5080 |
| `sphalerite` | (Zn,Fe)S | Reported endpoints 3900 and 4100 |
| `galena` | PbS | 7580 |
| `arsenopyrite` | FeAsS | 6070; source parenthetical magnitude 150 |
| `malachite` | Cu2CO3(OH)2 | 4050; source parenthetical magnitude 20 |

The last five profiles use the `-handbook-mineralogy` suffix. Their
[mineral reference guide](MINING_MINERALS.md) details the additional 19 claims,
including ordinal Mohs endpoints and the source's parenthetical density
magnitudes. Those magnitudes carry no inferred probability or confidence level.
These profiles also require their named mineral and unqualified specimen flags.

Single values use `density`. Ranges use the separate properties
`density-reported-range-minimum` and `density-reported-range-maximum`.
A reported range has neither an invented midpoint nor a probability law.
These records leave measurement temperature, pressure and statistical
uncertainty unstated; they do not support a temperature sweep.

For example, inspect the two chalcocite endpoints:

```bash
cargo run -p xtask -- matdb-catalog \
  --pack chalcocite-handbook-mineralogy-reference \
  --property density-reported-range-minimum \
  --property density-reported-range-maximum \
  --at mineral-reference-model=1 \
  --at source-test-temperature-known=0 \
  --at source-test-pressure-known=0 \
  --at bulk-ore-qualified=0
```

The packs do not substitute for bulk-ore density, porosity, mineral proportions,
liberation, Bond work index, flotation response or leach recovery. Existing
quartz, calcite, hematite, water and oxygen profiles remain independently
identified sources with their own conditions.

## Processing reagents and standard thermochemistry

Every reagent profile contains four claims: `molar-mass`,
`standard-molar-formation-enthalpy`, `standard-molar-entropy` and
`specific-heat-capacity`. Formation enthalpy and entropy apply only at
**298.15 K and 100000 Pa** with the declared JANAF elemental reference.

| Profile | Material/phase | Formation enthalpy, J/mol | Cp support, K |
| --- | --- | ---: | --- |
| `calcium-oxide-solid-nist-mining` | CaO(s), quicklime constituent | −635090 | 298.15–1000 |
| `calcium-hydroxide-solid-nist-mining` | Ca(OH)2(s), hydrated lime constituent | −986090 | 298.15–500 |
| `sodium-hydroxide-solid-nist-mining` | NaOH(s) | −425930 | 298.15–500 |
| `sodium-cyanide-solid-nist-mining` | NaCN(s) | −90710 | 298.15–500 |
| `sodium-carbonate-solid-nist-shomate` | Anhydrous Na2CO3(s), soda-ash constituent | −1130770 | 298.15–700 |
| `sulfuric-acid-liquid-janaf-mining` | Pure H2SO4(l) | −813989 | 298.15–400 |
| `water-liquid-janaf-mining` | Pure H2O(l) | −285830 | 298.15–360 |

The five solid profiles use [NIST Chemistry WebBook](https://webbook.nist.gov/)
data with explicitly sampled Shomate Cp curves. The liquid profiles retain
selected literal [JANAF H-075](https://janaf.nist.gov/tables/H-075.html) and
[H-063](https://janaf.nist.gov/tables/H-063.html) Cp points. Linear interpolation
is an authored reference approximation; source/interpolation uncertainties
remain unstated. Formation enthalpies are not extrapolated with temperature.
The chosen WebBook values are retained at their published rounding, separately
from more precise fit coefficients or alternate evaluations.

The [sodium carbonate source](https://webbook.nist.gov/cgi/cbook.cgi?ID=C497198&Mask=2)
provides the pure anhydrous solid's first Shomate branch, 298–723 K. The pack
retains 298.15–700 K with 42 Cp knots; it excludes the higher-temperature
branch, the melt, hydrates and dissolved carbonate/bicarbonate species.

All queries require `source-standard-state=1`,
`source-elemental-reference-janaf=1`, `source-pure-species=1`,
`engineering-reference-model=1`, `source-specimen-qualified=0`, the correct
`source-phase-solid=1` or `source-phase-liquid=1`, and typed temperature/pressure.
Pure liquid sulfuric acid does not represent an aqueous acid concentration;
solid NaOH/NaCN do not supply aqueous ion formation properties. Commercial
purity, slurry composition, solvation and heat of dilution remain separate.

The [reactive mineral thermochemistry guide](MINERAL_THERMOCHEMISTRY.md)
documents 11 additional claims for calcite, pyrite and alpha-quartz, including
bounded calcite/pyrite Cp curves. It preserves their separate USGS and JANAF
reference conventions and reuses the existing quartz caloric profile.

### Carbonate calcination and iron-oxide oxidation

Six additional profiles supply `molar-mass`,
`standard-molar-formation-enthalpy` and `standard-molar-entropy` at the same
298.15 K / 100000 Pa JANAF reference. Their selected crystalline and gaseous
species close two complete mixed-phase heat balances:

| Profile | Formation enthalpy, J/mol | Entropy, J/(mol K) | Original JANAF table |
| --- | ---: | ---: | --- |
| `magnesite-solid-janaf-mining` | −1111689 | 65.854 | [C-079](https://janaf.nist.gov/tables/C-079.html) |
| `magnesium-oxide-solid-janaf-mining` | −601241 | 26.924 | [Mg-008](https://janaf.nist.gov/tables/Mg-008.html) |
| `carbon-dioxide-gas-janaf-mining` | −393522 | 213.795 | [C-095](https://janaf.nist.gov/tables/C-095.html) |
| `magnetite-solid-janaf-mining` | −1120894 | 145.266 | [Fe-032](https://janaf.nist.gov/tables/Fe-032.html) |
| `hematite-solid-janaf-mining` | −825503 | 87.400 | [Fe-030](https://janaf.nist.gov/tables/Fe-030.html) |
| `oxygen-gas-janaf-mining` | 0 | 205.147 | [O-029](https://janaf.nist.gov/tables/O-029.html) |

The gas profiles require `source-phase-gas=1`; the solids require
`source-phase-solid=1`. Other reference flags are the same as above. Oxygen's
zero formation enthalpy still has an original source claim and receipt.
Existing heat-capacity and gas-EOS models remain separate sources; these six
profiles supply reference-point quantities without inventing temperature curves.

`MgCO3(s) -> MgO(s) + CO2(g)` gives **+116926 J/mol of reaction**. At
84.3139 kg of pure carbonate actually reacted, the source molar masses give
approximately 1000 mol extent, 44.0095 kg CO2 formed and **+116.926 MJ**.

`4 Fe3O4(s) + O2(g) -> 6 Fe2O3(s)` gives **−469442 J/mol of the written
reaction**, or −117360.5 J per mole of magnetite consumed. Fe2O3 denotes the
hematite endpoint. These are standard chemical-energy changes; sensible heat,
phase stability, conversion, kinetics and kiln fuel duty need further inputs.

### Copper oxide reduction and oxidation

Four more profiles use the same typed JANAF standard-state convention. Each
contains molar mass, formation enthalpy and absolute molar entropy; all are
reference points at **298.15 K and 100000 Pa**.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Entropy, J/(mol K) | Source |
| --- | ---: | ---: | ---: | --- |
| `copper-ii-oxide-solid-janaf-mining` | 0.079545 | −156063 | 42.594 | [Cu-015](https://janaf.nist.gov/tables/Cu-015.html) |
| `copper-i-oxide-solid-janaf-mining` | 0.143091 | −170707 | 92.360 | [Cu-019](https://janaf.nist.gov/tables/Cu-019.html) |
| `copper-solid-janaf-mining` | 0.063546 | 0 | 33.164 | [Cu-002](https://janaf.nist.gov/tables/Cu-002.html) |
| `carbon-monoxide-gas-janaf-mining` | 0.0280101 | −110527 | 197.653 | [C-093](https://janaf.nist.gov/tables/C-093.html) |

Molar masses retain the NIST WebBook display precision cited in each manifest.
Metallic copper has zero standard formation enthalpy and nonzero absolute
entropy; both remain source-backed claims. Together with the existing O2/CO2
profiles, these close three explicit reactions:

| Written reaction | Standard enthalpy, J/mol of equation extent |
| --- | ---: |
| `2 Cu2O(s) + O2(g) -> 4 CuO(s)` | −282838 |
| `CuO(s) + CO(g) -> Cu(s) + CO2(g)` | −126932 |
| `Cu2O(s) + CO(g) -> 2 Cu(s) + CO2(g)` | −112288 |

For 79.545 kg of pure CuO actually reacted, the reduction consumes 28.0101 kg
CO and forms 63.546 kg Cu plus 44.0095 kg CO2, with approximately −126.932 MJ
standard enthalpy change. The independently rounded source molar masses give
a **+0.0004 kg products-minus-reactants residual** at that extent. The example
prints this discrepancy; it never changes an attributed source to force a
zero. Exact atom/charge conservation applies to the integer equation.
This supplies reductant and metal bookkeeping for an explicit conversion,
not a smelter operating temperature, reaction rate, conversion forecast or
furnace fuel duty. Hot, molten or aqueous states refuse these source points.

```bash
copper_pack_dir="$(mktemp -d)"
for pack in copper-ii-oxide-solid-janaf-mining copper-i-oxide-solid-janaf-mining copper-solid-janaf-mining carbon-monoxide-gas-janaf-mining carbon-dioxide-gas-janaf-mining oxygen-gas-janaf-mining; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$copper_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction cupric-oxide-reduction \
  "$copper_pack_dir/copper-ii-oxide-solid-janaf-mining.fsmatpk" \
  "$copper_pack_dir/carbon-monoxide-gas-janaf-mining.fsmatpk" \
  "$copper_pack_dir/copper-solid-janaf-mining.fsmatpk" \
  "$copper_pack_dir/carbon-dioxide-gas-janaf-mining.fsmatpk" \
  79.545
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction cuprous-oxide-reduction \
  "$copper_pack_dir/copper-i-oxide-solid-janaf-mining.fsmatpk" \
  "$copper_pack_dir/carbon-monoxide-gas-janaf-mining.fsmatpk" \
  "$copper_pack_dir/copper-solid-janaf-mining.fsmatpk" \
  "$copper_pack_dir/carbon-dioxide-gas-janaf-mining.fsmatpk" \
  143.091
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction cuprous-oxide-oxidation \
  "$copper_pack_dir/copper-i-oxide-solid-janaf-mining.fsmatpk" \
  "$copper_pack_dir/oxygen-gas-janaf-mining.fsmatpk" \
  "$copper_pack_dir/copper-ii-oxide-solid-janaf-mining.fsmatpk" \
  286.182
```

## Run a source-backed reaction balance

Compile the three native source packs and run the checked example:

```bash
mining_pack_dir="$(mktemp -d)"
for pack in calcium-oxide-solid-nist-mining water-liquid-janaf-mining calcium-hydroxide-solid-nist-mining; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$mining_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  "$mining_pack_dir/calcium-oxide-solid-nist-mining.fsmatpk" \
  "$mining_pack_dir/water-liquid-janaf-mining.fsmatpk" \
  "$mining_pack_dir/calcium-hydroxide-solid-nist-mining.fsmatpk" \
  56.077
```

The final argument is **kilograms of pure CaO actually reacted**, not commercial
feed or a predicted conversion. The source molar mass yields approximately
1000 mol extent and 1000 mol stoichiometric water consumption. The selected
formation values give

`CaO(s) + H2O(l) -> Ca(OH)2(s)`

`Delta H = -986090 - (-635090) - (-285830) = -65170 J/mol`.

At that extent the standard enthalpy change is approximately **−65.17 MJ**.
The example prints source pack identities/hashes and retains query receipts.
It does not compute plant temperature, heat transfer, dissolution, conversion
rate or an operating dose. Original uncertainty is retained without producing
a fabricated error band. Independently rounded source molar masses need not
be bit-exactly additive; exact elemental conservation is a separate proof.

The same executable also supports the two new mineral reactions. Compile their
sources, then select the desired reaction explicitly:

```bash
for pack in magnesite-solid-janaf-mining magnesium-oxide-solid-janaf-mining carbon-dioxide-gas-janaf-mining magnetite-solid-janaf-mining oxygen-gas-janaf-mining hematite-solid-janaf-mining; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$mining_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction magnesite-calcination \
  "$mining_pack_dir/magnesite-solid-janaf-mining.fsmatpk" \
  "$mining_pack_dir/magnesium-oxide-solid-janaf-mining.fsmatpk" \
  "$mining_pack_dir/carbon-dioxide-gas-janaf-mining.fsmatpk" \
  84.3139
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction magnetite-oxidation \
  "$mining_pack_dir/magnetite-solid-janaf-mining.fsmatpk" \
  "$mining_pack_dir/oxygen-gas-janaf-mining.fsmatpk" \
  "$mining_pack_dir/hematite-solid-janaf-mining.fsmatpk" \
  926.132
```

The oxidation invocation consumes 4000 mol magnetite at 1000 mol extent,
requiring 1000 mol O2 and forming 6000 mol hematite. The executable verifies
the expected pack identities, prints every consumed/formed amount and mass,
and reports the standard enthalpy change with its sign.

### Other callable reaction foundations

`ConservedReaction` accepts existing `ElementalMatrix`, `StoichiometricMatrix`
and `ChargeVector` artifacts. Admission checks exact element and charge
balance. `amounts` scales the equation by a supplied extent;
`extent_from_reacted_mass` makes the mass/mole basis explicit;
`limiting_extent` bounds possible extent from all reactant inventories, with a
conservative rounding check against reagent overdraw. `advance_inventory`
applies an explicitly supplied extent to every active species, including
initially absent products entered as zero. It returns canonical inventories
without mutating the inputs and refuses extents exceeding available reagents.

The same example binary can run inventory accounting without property packs:

```bash
cargo run -p fs-thermochem --example mineral_processing -- --malachite-inventory 10 12 100 3
```

The inputs are 10 mol malachite, 12 mol protons, 100 mol initial water and a
caller-declared extent of 3 mol. The final inventory contains 7 mol malachite,
0 mol protons, 6 mol dissolved Cu(II), 3 mol dissolved CO2 and 109 mol water.
An extent above 3 mol refuses. This is a declared material-flow calculation;
the availability bound does not predict reaction progress or a commercial dose.

`MineralReaction::build` provides twelve explicit endpoint choices:

| Reference | What the equation accounts for |
| --- | --- |
| `LimeHydration` | One mole water per mole CaO, forming solid Ca(OH)2 |
| `CalciteAcidNeutralization` | Two mol H+ per mol calcite, dissolved CO2 endpoint |
| `DolomiteAcidNeutralization` | Four mol H+ per mol dolomite, dissolved CO2 endpoint |
| `CalciteBicarbonateNeutralization` | One mol H+ per mol calcite, retained bicarbonate endpoint |
| `DolomiteBicarbonateNeutralization` | Two mol H+ per mol dolomite, retained bicarbonate endpoint |
| `HydratedLimeNeutralization` | Two mol H+ consumed and two mol water formed per mol Ca(OH)2 |
| `MalachiteAcidLeach` | Four mol H+ per mol malachite, forming two mol dissolved Cu(II) |
| `GoldCyanidation` | Two mol cyanide and 0.25 mol oxygen per mol gold, net ionic Elsner balance |
| `PyriteOxidationFerrous` | 3.5 mol O2 consumed and 2 mol H+ generated per mol pyrite |
| `PyriteOxidationFerricHydroxide` | 3.75 mol O2 consumed and 4 mol H+ generated per mol pyrite |
| `MagnesiteCalcination` | Forms one mol MgO and one mol gaseous CO2 per mol magnesite |
| `MagnetiteOxidationToHematite` | Consumes 0.25 mol O2 and forms 1.5 mol hematite per mol magnetite |

Pyrite equations use integer multiples of two/four FeS2; an extent is per
equation as written. These endpoints follow the primary references cited in
`stoichiometry/minerals.rs` and do not choose pH/Eh, speciation or precipitates.
Lime hydration, magnesite
calcination and magnetite oxidation have complete formation-enthalpy inputs
in this seed slice. Missing aqueous-ion data for other reaction heats remain
missing.

`faraday_yield(charge, molar_mass, electrons_per_formula, efficiency)` supplies
`n = Q * efficiency / (z * F)` and a typed mass/amount record. With explicit
Cu(II) reduction, 1000 Ah, 0.063546 kg/mol and efficiency 1, the independent
reference is 1.185494183238305 kg copper. Efficiency is a caller input;
electrode kinetics, current distribution and electrical energy are not inferred.

## Ball-mill energy and throughput

`fs_material::comminution` implements the uncorrected
[GMG 2021 Bond reference](https://gmggroup.org/wp-content/uploads/2024/07/GUIDELINE__Determining-the-Bond-Efficiency-of-Industrial-Grinding-Circuits_2021.pdf).
`BondWorkIndex` has separate metric-tonne and US-short-ton constructors.
`BondGrindingCard` requires an explicit specimen, test/preparation method,
closing screen, provenance and F80/P80 validity domain. Sizes passed to the
runtime are SI metres; energy is J/kg of dry solids. No mineral-name default
work index is supplied. The uncorrected scope requires P80 at least 70 µm.

The familiar micrometre form is
`W = 10 Wi (1/sqrt(P80) - 1/sqrt(F80))`.
For the **authored test scenario** Wi=15 kWh/dry metric tonne,
F80=2500 µm and P80=100 µm, energy is 12 kWh/t (43200 J/kg).
At 12 MW net pinion power and utilization 0.9, the corresponding arithmetic
throughput is 250 kg/s, or 900 dry t/h. A real site needs its own measured work
index and applicable circuit corrections; this example is not a plant forecast.

### Grinding a declared ore blend

`AdditiveGrindingBlend::new` accepts `GrindingFeed` records containing each
feed's original `BondGrindingCard`, dry mass and F80. It explicitly selects
an additive-specific-energy approximation at a common P80:
`E_blend = sum(m_i * E_i) / sum(m_i)`. Each source's feed and product bounds
remain enforced. The resulting energy is used to calculate throughput;
averaging work indices, feed sizes or individual throughputs is not the rule.

An authored two-feed example at P80=100 µm makes the difference concrete:

| Feed | Dry mass share | Wi, kWh/t | F80, µm | Energy, kWh/t |
| --- | ---: | ---: | ---: | ---: |
| A | 25% | 10 | 2500 | 8 |
| B | 75% | 20 | 400 | 10 |

The blend requires 9.5 kWh/t, or 34200 J/kg. At 1 MW net pinion power and
0.9 utilization, this is approximately 26.316 kg/s, or 94.737 dry t/h.
The approximation preserves each measured or authored test context but does
not prove that real ores grind independently when mixed. A site must establish
the additive regime or supply a different supported blend model.

### Solving for grind size from energy or throughput

The original `BondGrindingCard` also exposes
`product_size_for_specific_energy(feed, energy)` and
`product_size_for_throughput(feed, pinion_power, utilization, dry_mass_flow)`.
They invert the same relation with the original test context and validity
window. For Wi=15 kWh/t and F80=2500 µm, either 43200 J/kg or 12 MW at 0.9
utilization and 250 kg/s returns P80=100 µm. Higher flow at fixed power implies
coarser grinding within this model. A target outside the card's supported
size window, or positive flow during zero utilization, refuses.

`ball_mill_critical_speed` returns the ideal single-ball centrifuging frequency
from inside-liner diameter, ball diameter and explicit gravity:
`n = sqrt(2*g/(D-d))/(2*pi)` revolutions/second.
It is not an operating-speed recommendation, power-draw model or DEM solver.

## From dry mill throughput to slurry and process water

`fs_material::slurry` combines dry constituent mass flows with one liquid
carrier. `DensitySource::resolve` selects a real pack claim at its original
typed conditions and retains the material evidence and verified usage receipt.
The explicitly selected additive-volume model computes
`Q = sum(mass_flow_i / density_i)` and `pulp_density = total_mass_flow / Q`.
It exposes solids mass and volume fractions separately. Unknown mineral test
temperatures stay unknown; this engineering reference mixture does not become
a measured bulk ore or a qualified common-temperature material.

`dilute_to_mass_fraction` computes carrier addition to reach an explicit mass
fraction. Requests requiring removal return `WaterRemovalRequired` with the
removal rate. `split` is nonselective and `mix` conserves each named constituent;
conflicting source/query contexts cannot silently merge. Zero flow is valid,
with density and composition absent instead of NaN or invented defaults.

The runnable example uses the Bond scenario above, a declared 50% magnetite /
20% pyrite / 30% dolomite dry blend, an initial 50 kg/s water flow, and a target
of 65% solids by mass. Its outputs are:

| Output | Authored reference scenario |
| --- | ---: |
| Dry solids | 250 kg/s |
| Added water | 84.615384615 kg/s |
| Total water | 134.615384615 kg/s |
| Total slurry volume | 0.195356574114 m³/s |
| Pulp density | 1968.786493931 kg/m³ |
| Solids mass fraction | 0.65 |
| Solids volume fraction | 0.308883870454 |

The water density comes from the existing IAPWS liquid-water pack at 298.15 K
and 100000 Pa. This is an explicit reference calculation, not a site forecast.

```bash
for pack in magnetite-handbook-mineralogy-reference pyrite-handbook-mineralogy-reference dolomite-handbook-mineralogy-reference water-liquid-iapws-sr6-08; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$mining_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-material --example mineral_slurry -- \
  "$mining_pack_dir/magnetite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/pyrite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/dolomite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/water-liquid-iapws-sr6-08.fsmatpk" \
  0.65
```

The optional final argument is a solids mass fraction, not a percentage.
The executable checks pack identities and prints every source hash. It does
not infer packing, settling, entrained gas, viscosity or a separation recovery.

### Yield stress and pipeline pressure requirements

`fs_material::slurry::rheology` adds a source-qualified Herschel–Bulkley law
and its laminar circular-pipe integral. The fit declares yield stress, reference
excess stress at a reference shear rate and a dimensionless flow index. Its
temperature, solids mass fraction and shear-rate limits remain mandatory.
`SlurryConditions::from_stream` carries the existing slurry's mass fraction
into the law while requiring an explicit temperature.

`LaminarCircularPipe::pressure_for_flow` solves for frictional pressure at a
requested signed volume flow. It reports the unsheared core, wall stress/rate,
hydraulic power and actual numerical residual. The source domain and an
explicit caller-declared laminar ceiling bound the solve.

```bash
cargo run -p fs-material --example mineral_pipe -- 0.005 0.6
```

The authored example uses a 0.1 m inner radius, 100 m length, 298.15 K,
25 Pa yield stress, 8 Pa excess stress at 10/s, and flow index 0.6.
At 0.005 m³/s it requires approximately 71680.085 Pa, with a 0.069754 m
unsheared-core radius and 358.400 W hydraulic dissipation. Reversing flow
reverses pressure while retaining positive dissipation. Solids mass fraction
0.9 and volume flow 100 m³/s refuse the declared support. These are arithmetic
scenario results, not measured mineral rheology. Turbulence, settling, entrance
losses, fittings, elevation and pump efficiency require separate closures.

## Particle classification and component recovery

`fs_material::separation` supplies size-resolved `ParticleStream` values with
explicit dry solids, water and solid-phase element inventories. Every class
declares the same element axes, so an unassayed element cannot become zero.
Splits carry each class and its composition into both outlets. Water uses its
own partition, and product grade is calculated from actual product dry mass.
A zero-solids outlet has no defined grade.

`ClassifierCard` uses an original source card with corrected cut size d50,
sharpness m and water-bypass fraction b. The corrected underflow fraction is
`1-exp(-ln(2)*(d/d50)^m)`; actual underflow is `b+(1-b)*corrected`.
It evaluates the stream's actual solids mass fraction and each class diameter
against the source domain. Geometry, liberation or missing test conditions are
not inferred. Independent product/reject evaluation retains small tails.

`fs_material::recovery::distributed` adds single, rectangular, gamma and
explicit discrete rate distributions. Batch and ideal plug flow integrate
the rate distribution; a steady train of N equal mixed tanks instead uses
`(1+k*t/N)^(-N)` for each population's remaining fraction, where t is total
mean residence. Mixed tanks accept single or discrete rates; continuous
distributions refuse that unsupported closure. All rates, ultimate recoveries,
specimen contexts and typed conditions must be supplied explicitly.

```bash
cargo run -p fs-material --example mineral_separation -- 120
```

This authored single-pass example combines Bond grinding with four explicit
valuable/gangue and fine/coarse classes, a classifier and two mixed recovery
tanks. Valuable particles have two rate populations; gangue has its own rate.
No ore name supplies these parameters. With the example's default inputs:

| Stream | Dry solids, kg/s | Cu, kg/s | Water, kg/s | Cu grade, dry mass |
| --- | ---: | ---: | ---: | ---: |
| Mill product | 32.866187724 | 3.286618772 | 32.866187724 | 10.000000% |
| Classifier underflow | 27.388489770 | 2.662161206 | 6.573237545 | 9.720000% |
| Classifier overflow / recovery feed | 5.477697954 | 0.624457567 | 26.292950179 | 11.400000% |
| Concentrate | 1.171076732 | 0.330523973 | 2.629295018 | 28.223938% |
| Tailings | 4.306621222 | 0.293933593 | 23.663655161 | 6.825155% |

Underflow plus concentrate plus tailings close the feed balances: displayed
residuals are zero for dry solids and water and 4.44e-16 kg/s for copper.
The throughput uses Wi=15 kWh/dry metric tonne, F80=2500 µm, P80=200 µm,
1 MW net pinion power and 0.9 utilization. The full PSD is an explicit authored
input consistent with that discrete P80; the Bond equation does not predict it.
Zero residence produces no recovered solids and an undefined concentrate grade.
Residence beyond the declared 3600 s refuses. This example leaves underflow as
an outlet; a closed recycle circuit, liberation/entrainment and reagent response
still need their own models and site evidence.

## Recovery curves and material held in process tanks

`fs_material::recovery::RecoveryKinetics` evaluates two explicitly selected
contact-time models from an original `ConstitutiveModelCard`:

| Law ID | Declared rate population | Recovery at contact time t |
| --- | --- | --- |
| `single-rate-recovery` | All recoverable material has rate k | `Rinf * (1 - exp(-k*t))` |
| `uniform-rate-recovery` | Rates uniformly distributed over 0 to kmax | `Rinf * (1 - (1 - exp(-kmax*t))/(kmax*t))` |

Both have a defined zero-time limit of zero. The card supplies dimensionless
`ultimate-recovery` and an inverse-second `rate-constant` or
`maximum-rate-constant`, with a typed `elapsed-time` validity window and every
applicable grind, chemistry or temperature coordinate. Missing sources,
unsupported model versions, incorrect dimensions and out-of-box conditions
refuse. The full card, original source hashes and provenance stay available.
These formulas do not select a rate distribution or calibrate parameters for
an ore. For the authored Rinf=0.9 and rate endpoint 1/(6 h), the two curves
give respectively 0.568909 and 0.331091 at six hours; the model choice matters.

### Retained inventory and changing feed grades

`recovery::residence::StirredTankCascade` supplies the separate forward
hydraulic calculation for equal, full, perfectly mixed tanks. It advances a
conserved tracer/metal concentration in **kg per m³ of slurry**, retaining
each tank's state between intervals. Ore grade in g/dry tonne must first be
converted using dry-solids mass per slurry volume. The slurry module's actual
volume-flow result can be supplied directly as the train flow.

An interval reports inlet and independently integrated outlet mass, mean
outlet concentration, initial/final inventory and a conservation residual.
Changing flow changes exchanged volume; stopped flow retains inventory and
has no outlet mean. Failed intervals leave the prior state intact. The caller
supplies active mixed volume and actual train flow, explicitly excluding any
bypass/dead space. Mean residence is active volume divided by train flow.

The example reproduces the plan's authored eight-tank, 24-hour mean-residence
scenario: feed grade is normally 3 g/t, rises to 9 g/t for one 12-hour interval,
then returns to 3 g/t. Applying a separately declared 88% extraction at the
transport outlet yields these reporting-window results:

| 12-hour interval | Mean tail grade, g/t | Naive same-interval recovery |
| --- | ---: | ---: |
| Feed pulse | 0.366053 | 95.932746% |
| First interval after pulse | 0.548899 | 81.703370% |
| Second interval after pulse | 0.713954 | 76.201535% |
| Third interval after pulse | 0.504046 | 83.198457% |

The apparent changes arise from material held in the tanks while feed grade
changes. The extraction fraction remains 88%. The example applies that
fraction downstream of transport; it does not couple a kinetic law into the
hydraulic state or infer in-tank recovery. Its independent contact-time curve
table illustrates the two recovery laws separately.

```bash
cargo run -p fs-material --example mineral_residence
```

The output includes source/card hashes, per-interval inventory changes and
independent balance residuals, followed by stopped-flow and half-flow cases.
The implementation admits at most 64 tanks and dimensionless exposure 4096
per call; longer intervals can be split into successive calls. These are
forward model calculations, not deconvolution, measured plant calibration or
a fitted leach/adsorption model.

## Focused verification

Normal workspace commands for these slices:

```bash
cargo test -p xtask --test mining_mineral_sources --test mining_mineral_coverage --test mining_reagent_sources --test mining_reagent_references --test mining_oxide_sources --test mining_copper_sources --test mining_reaction_balances
cargo test -p fs-material comminution
cargo test -p fs-material slurry::tests
cargo test -p fs-material recovery::
cargo test -p xtask --test mining_slurry_sources
cargo test -p fs-thermochem stoichiometry::tests
```

The source/consumer tests exercise the actual compiler, canonical decoding,
typed query receipts, unsupported-condition refusals, known-answer physics,
deterministic ordering and source-reference mutations. They establish these
bounded data/formula capabilities, not site calibration or full-workspace
qualification. The source-backed command-line example is an additional real
consumer of the shipped packs.
