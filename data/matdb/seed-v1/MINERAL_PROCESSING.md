# Mineral processing: public sources and callable physics

Native-silver cyanidation and shared-zinc gold/silver cementation now have
matched material references, finite-inventory reactions and a runnable
source-pack consumer; see [precious-metal recovery](MINING_PRECIOUS_METALS.md).

These source profiles provide reusable mineral
and reagent inputs for FrankenMining's material-corpus and Crucible work.
They use the existing `matdb-pack` compiler, `NormalizedPack` queries and
`matdb-catalog` discovery. No separate material registry or schema is required.

The physical consumers include `fs_material::comminution`, `fs_material::separation`,
`fs_material::slurry`, `fs_material::recovery`, `fs_material::solvent_extraction`,
`fs_material::aqueous_adsorption` and `fs_thermochem::stoichiometry`.
These are upstream building blocks: they do
not implement a FrankenMining flowsheet or complete its site-calibration
beads. In particular, FrankenMining's plan §13.10 keeps located ore tests
separate from public non-spatial material references, and §17.3 requires
declared validity for processing correlations.

## Solution recovery and ferric treatment

The [solution-recovery guide](MINING_SOLUTION_RECOVERY.md) exposes finite
solvent extraction and stripping from original equilibrium isotherms, a runnable
extraction-to-electrowinning chain, and finite-bath aqueous Langmuir adsorption
with explicit carbon capacity and loaded-carbon desorption. The examples label
their authored parameters; measured process cards keep their own source support.

The [iron-precipitation guide](MINING_IRON_PRECIPITATION.md) adds compatible
USGS goethite and portlandite profiles. Their source heat and exact conserved
endpoints connect ferric hydrolysis, proton generation and hydrated-lime
neutralization. The processing CLI runs each endpoint from the actual compiled
packs, including source-matched enthalpy/entropy equilibrium diagnostics.

## Copper-carbonate leaching and acidic gangue demand

The [carbonate-leaching guide](MINING_CARBONATE_LEACHING.md) completes malachite
leaching and calcite's two-proton carbonic endpoint with compatible USGS
material packs, actual mass and acid accounting, reaction heat, and separate
direct/caloric Gibbs calculations. The carbonic component retains the source's
formal aqueous-CO2-plus-water basis and requires its explicit reference selector.

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
| `stannite` | Cu2FeSnS4 | Reported endpoints 4300 and 4500 |
| `cuprite` | Cu2O | 6140 |
| `tenorite` | CuO | 6450, explicitly synthetic source material |
| `kaolinite` | Al2Si2O5(OH)4 | Reported endpoints 2610 and 2680 |

The last nine profiles use the `-handbook-mineralogy` suffix. Their
[mineral reference guide](MINING_MINERALS.md) details the additional 33 claims,
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
documents solid-phase reference data for calcite, dolomite, pyrite, alpha-quartz,
gypsum and anhydrite, including bounded carbonate/pyrite Cp curves and nine literal anhydrite
Cp points through 1000 K. It preserves their separate USGS and JANAF reference
conventions and reuses the existing quartz caloric profile. Gypsum's hydrated
reference is distinct from anhydrite and from commercial wallboard records;
no missing gypsum Cp, aqueous properties or dehydration kinetics are inferred.
Compatible USGS H⁺, Ca²⁺, Mg²⁺ and bicarbonate profiles supply aqueous data
and close the calcite and dolomite bicarbonate-neutralization heat balances.
The executable `mineral_processing` consumer reports both net ionic acid demand
and standard reaction heat from those actual source packs; see the guide for
pack order, the 298.15 K reference state and the separate Cp domains.

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

### Copper-oxide acid leaching

The USGS tenorite profile adds six claims and completes
`CuO(s) + 2 H+(aq) -> Cu[2+](aq) + H2O(l)` with the existing USGS copper-ion,
proton and water sources. `MineralReaction::CupricOxideAcidLeach` supplies the exact
conserved endpoint; `mineral_processing --reaction cupric-oxide-acid-leach`
reports its four-pack balance. The source heat is −64800 J/mol extent and
the direct formation-Gibbs route gives −43700 J/mol. See the
[copper-leaching source guide](MINING_COPPER_LEACHING.md) for source values,
pack commands, inventory limits and the distinct entropy-based Gibbs route.

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

### Sulfide roasting and sulfur-to-acid balances

Two new gas profiles close the sulfur reaction chain with the existing pyrite,
oxygen, hematite, water and sulfuric-acid sources. Each adds `molar-mass`,
`standard-molar-formation-enthalpy` and `standard-molar-entropy`, at
298.15 K / 100000 Pa with the same explicit JANAF elemental reference.

| Profile | Molar mass, kg/mol | Formation enthalpy, J/mol | Entropy, J/(mol K) | JANAF table |
| --- | ---: | ---: | ---: | --- |
| `sulfur-dioxide-gas-janaf-mining` | 0.064064 | −296842 | 248.212 | [O-034](https://janaf.nist.gov/tables/O-034.html) |
| `sulfur-trioxide-gas-janaf-mining` | 0.080063 | −395765 | 256.769 | [O-058](https://janaf.nist.gov/tables/O-058.html) |

These retain gas standard states, not an assertion that a gas is stable under
arbitrary plant conditions. No temperature curve or uncertainty band is invented.
The existing `pyrite-solid-janaf-thermochemistry` pack retains JANAF's own
119.967 g/mol mass and −171544 J/mol formation enthalpy. Its named pyrite phase
condition is required; another FeS2 phase or a newer molar-mass evaluation is
not silently substituted. The example selects that phase explicitly.

| Endpoint | Integer equation | Standard heat, J/mol of written equation |
| --- | --- | ---: |
| `pyrite-roasting` | 4 FeS2 + 11 O2 → 2 Fe2O3 + 8 SO2 | −3339566 |
| `sulfur-dioxide-oxidation` | 2 SO2 + O2 → 2 SO3 | −197846 |
| `sulfur-trioxide-hydration` | SO3(g) + H2O(l) → H2SO4(l) | −132394 |
| `pyrite-to-acid` | 4 FeS2 + 15 O2 + 8 H2O → 2 Fe2O3 + 8 H2SO4(l) | −5190102 |

The last row is independently element/charge balanced and queries all five
sources. Its heat equals one roasting extent plus four SO2-oxidation extents
plus eight SO3-hydration extents. This Hess-cycle check also catches a wrong
reaction-mole basis. All sulfur is assigned to pure liquid acid in that endpoint;
conversion, capture efficiency, catalytic rates, aqueous acid concentration,
heat of dilution and thermal losses require additional inputs.

Compile the five packs and run the net balance for 119.967 kg of pure pyrite
actually reacted:

```bash
sulfur_pack_dir="$(mktemp -d)"
for pack in pyrite-solid-janaf-thermochemistry oxygen-gas-janaf-mining water-liquid-janaf-mining hematite-solid-janaf-mining sulfuric-acid-liquid-janaf-mining; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$sulfur_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction pyrite-to-acid \
  "$sulfur_pack_dir/pyrite-solid-janaf-thermochemistry.fsmatpk" \
  "$sulfur_pack_dir/oxygen-gas-janaf-mining.fsmatpk" \
  "$sulfur_pack_dir/water-liquid-janaf-mining.fsmatpk" \
  "$sulfur_pack_dir/hematite-solid-janaf-mining.fsmatpk" \
  "$sulfur_pack_dir/sulfuric-acid-liquid-janaf-mining.fsmatpk" \
  119.967
```

That source mass represents approximately 1000 mol pyrite and 250 mol of the
written equation. The declared endpoint consumes 3750 mol O2 and 2000 mol water,
forms 500 mol hematite and 2000 mol pure acid, and gives approximately
−1297.5255 MJ standard enthalpy change. The program reports source masses and
their mass-sum residual, retaining differences in source atomic-weight
evaluations and rounding rather than adjusting attributed values. These
quantities are not an operating-dose or furnace temperature prediction.
`--help` lists pack order for each individual stage.

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

`MineralReaction::build` provides explicit endpoint choices, including:

| Reference | What the equation accounts for |
| --- | --- |
| `LimeHydration` | One mole water per mole CaO, forming solid Ca(OH)2 |
| `CalciteAcidNeutralization` | Two mol H+ per mol calcite, dissolved CO2 endpoint |
| `CalciteCarbonicAcidNeutralization` | Two mol H+ per mol calcite, source-defined neutral carbonic component endpoint |
| `DolomiteAcidNeutralization` | Four mol H+ per mol dolomite, dissolved CO2 endpoint |
| `CalciteBicarbonateNeutralization` | One mol H+ per mol calcite, retained bicarbonate endpoint |
| `DolomiteBicarbonateNeutralization` | Two mol H+ per mol dolomite, retained bicarbonate endpoint |
| `HydratedLimeNeutralization` | Two mol H+ consumed and two mol water formed per mol Ca(OH)2 |
| `MalachiteAcidLeach` | Four mol H+ per mol malachite, forming two mol dissolved Cu(II) |
| `MalachiteCarbonicAcidLeach` | Four mol H+ per mol malachite, two mol Cu2+, one formal H2CO3 component and two explicit waters |
| `GoldCyanidation` | Two mol cyanide and 0.25 mol oxygen per mol gold, net ionic Elsner balance |
| `PyriteOxidationFerrous` | 3.5 mol O2 consumed and 2 mol H+ generated per mol pyrite |
| `PyriteOxidationFerricHydroxide` | 3.75 mol O2 consumed and 4 mol H+ generated per mol pyrite |
| `MagnesiteCalcination` | Forms one mol MgO and one mol gaseous CO2 per mol magnesite |
| `MagnetiteOxidationToHematite` | Consumes 0.25 mol O2 and forms 1.5 mol hematite per mol magnetite |
| `CuprousOxideOxidation` | Consumes 0.5 mol O2 and forms two mol CuO per mol Cu2O |
| `CupricOxideCarbonMonoxideReduction` | One mol CO per mol CuO, forming copper and CO2 |
| `CupricOxideAcidLeach` | Two mol H+ per mol tenorite CuO, forming one mol Cu2+ and one mol water |
| `CuprousOxideCarbonMonoxideReduction` | One mol CO per mol Cu2O, forming two mol copper and one mol CO2 |
| `PyriteRoastingToHematite` | Consumes 2.75 mol O2 and forms 0.5 mol hematite and two mol SO2 per mol pyrite |
| `SulfurDioxideOxidation` | Consumes 0.5 mol O2 and forms one mol SO3 per mol SO2 |
| `SulfurTrioxideHydration` | One mol water per mol SO3, forming pure liquid acid |
| `PyriteToHematiteAndSulfuricAcid` | Net declared endpoint: 3.75 mol O2 and two mol water per mol pyrite, forming 0.5 mol hematite and two mol pure acid |
| `CopperElectrowinning` | Productive pathway consumes one mol Cu2+ and one mol water per mol deposited Cu, forming 0.5 mol O2 and two mol H+; two electrons per deposited Cu |

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

## Finite-time reaction with a depleting reagent inventory

`fs_material::reaction_batch::BilinearBatchReaction` connects a conserved
two-reactant equation to an explicit source-bound kinetic law. For the written
equation `a A + b B -> products`, it solves
`dξ/dt = (k/V) (nA0-aξ) (nB0-bξ)` at constant volume and temperature. The
coefficient has units m³/(mol s); both empirical orders are one, independently
of the stoichiometric coefficients. Lower acid inventory reduces the rate as
well as the maximum possible conversion. The runtime returns complete initial
and final species inventories, gross consumed/formed amounts, actual extent,
the availability bound, and the original equation and source card.

The original card must bound typed elapsed time, volume, temperature and both
reactant concentrations. The runtime supplies these coordinates from actual
inputs and admits the complete monotone path through its endpoints. Additional
declared conditions remain mandatory; a caller cannot substitute a different
concentration in metadata. Zero reagent/time/rate preserves the complete
inventory after source admission. Under-resolved progress and loss of a
positive finite-time remaining stock refuse rather than becoming artificial
complete conversion.

The executable example uses the dolomite bicarbonate endpoint and an explicitly
**synthetic** `k=0.005 m³/(mol s)`, 1 m³ reference volume and 298.15 K. At 10 s,
10 mol dolomite plus 20 mol H⁺ gives 5 mol extent and consumes 10 mol H⁺. With
only 2 mol initial H⁺, progress is slower and stays below the 1 mol availability
ceiling. An ore-specific coefficient must come from compatible measurements;
this example is not a published dolomite rate. Mineral inventory/volume is
an empirical availability variable, not a solid activity or a surface law.

The optional material-pack mode computes standard heat at that kinetic extent
from the complete compatible USGS dolomite/H⁺/Ca²⁺/Mg²⁺/bicarbonate source set:

```bash
leach_pack_dir=/tmp/mining-materials
mkdir -p "$leach_pack_dir"
for pack in dolomite-solid-usgs-b2131-thermochemistry hydrogen-ion-aqueous-usgs-b2131 calcium-ion-aqueous-usgs-b2131 magnesium-ion-aqueous-usgs-b2131 bicarbonate-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$leach_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-material --example mineral_leaching -- \
  --thermochemistry "$leach_pack_dir"
cargo test -p fs-material reaction_batch
```

The 5 mol case gives −326500 J from the source-backed −65300 J/mol reaction
enthalpy. This combines a declared kinetic model with real standard-state
thermochemistry; it does not infer calibration, mixing heat, finite-solution
activities, coupled reaction competition or plant temperature.

### Simultaneous mineral competition for a shared reagent

`fs_material::reaction_batch::network::BilinearBatchNetwork` advances several
original bilinear reaction models against one actual shared reagent stock.
This makes gangue consumption reduce the acid available to valuable minerals
throughout the calculation. Models retain their own kinetic cards and validity;
each numerical stage uses the same inventory and absolute elapsed time.

```bash
cargo run -p fs-material --example competing_mineral_leaching
```

The example starts with 10 mol CuO, 10 mol calcite and 12 mol H+ in 1 m3 at
298.15 K, and advances 30 seconds. Both equations consume two protons per
mineral formula unit; calcite uses the explicit neutral carbonic component
basis described in [carbonate leaching](MINING_CARBONATE_LEACHING.md). Rates
are authored examples, not measured ore constants. The returned `NetworkAdvance`
retains each reaction's extent, gross consumption/production and complete
remaining inventory for subsequent source-backed heat or process calculations.

The caller controls maximum time step, work and estimated extent error.
Unsupported source states or exhausted budgets return no partial result.
The bounded model requires one shared consumed reagent, a distinct other
reactant per equation, and no product-to-reactant feedback. It holds volume,
temperature and coefficients fixed; it does not infer pH, speciation, thermal
feedback or a calibrated recovery. Step-doubling estimates remain numerical
estimates, not certified global bounds.

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

### Continuous grinding with residence-dependent particle sizes

`fs_material::comminution::population_balance::ContinuousPopulationBalanceMill`
solves a steady, perfectly mixed mill from original selection frequencies
[1/s] and daughter mass fractions for each size class. `at_residence_time`
uses a supplied residence, while `at_dry_holdup` derives residence from dry
hold-up / actual dry feed. All rates and daughters retain their original
specimen, method, source and typed operating limits.

```bash
cargo run -p fs-material --example continuous_ball_mill
```

The synthetic example holds 6000 kg dry solids in the mill. Changing dry feed
from 100 to 200 kg/s shortens residence from 60 to 30 seconds and changes the
fine-class product fraction from 52.25% to 35.60%. Dry mass, every assayed
element and outlet water remain accounted for. The inferred 3000 kg water
hold-up assumes water shares the same perfect-mixer residence as solids.

The coarse-to-fine population balance uses strictly finer daughter classes and
retains a terminal nonbreaking fine class. Selection rates are explicit
source inputs; Bond work index and per-pass redistribution do not supply them.
This is a local steady first-order grinding model, with no inferred plant
calibration, liberation, transient response, size-dependent discharge or media
trajectory. Its [model-family reference](https://doi.org/10.3390/min12060679)
does not supply the example's authored parameters.

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

### Closing the mill/classifier recycle loop

`comminution::breakage::SizeRedistributionCard` supplies a source-declared
per-pass size-transfer matrix over at most 128 named classes. Columns retain
unit mass within a stated binary64 allowance without normalization. Material
either keeps its class or moves to a strictly smaller diameter; transferred
fragments carry the parent's explicit elemental composition. Water passes
unchanged. The card is a congruent size-redistribution approximation, not an
inferred liberation or first-principles kinetic model.

`comminution::circuit::ClosedGrindingCircuit` combines that card with the
existing classifier. With fresh flow f, mill feed x, mill matrix B and classifier
underflow fractions r, it solves `x = f + diag(r) B x` in decreasing-size order.
No iterative circulating-load guess is required. Water closes independently.
The source cards then evaluate the actual solved mill feed and discharge;
source limits on solids fraction and dry mill throughput cannot be satisfied
with caller-supplied substitute values. Every local and external dry-solid,
element and water account is checked against an explicit numerical tolerance.
The caller separately declares the maximum acceptable circulating load.

```bash
cargo run -p fs-material --example mineral_circuit
cargo run -p fs-material --example mineral_circuit -- 0.75
```

The optional argument is the fraction of coarse material retained in its size
class on each pass, default 0.5. The remaining coarse material becomes fine.
The authored two-size example supplies fresh coarse/fine masses and copper
assays, a 100 µm corrected classifier cut, sharpness 2 and water bypass 0.2.
Its default closed balance is:

| Stream | Dry solids, kg/s | Cu, kg/s | Water, kg/s | Discrete P80, µm |
| --- | ---: | ---: | ---: | ---: |
| Fresh feed | 12.000000 | 2.200000 | 12.000000 | 200 |
| Mill feed | 38.333333 | 7.166667 | 15.000000 | 200 |
| Mill discharge | 38.333333 | 7.166667 | 15.000000 | 200 |
| Recycle underflow | 26.333333 | 4.966667 | 3.000000 | 200 |
| Final overflow product | 12.000000 | 2.200000 | 12.000000 | 100 |

Circulating load is 219.444% of fresh dry solids. The largest checked balance
discrepancy is approximately 3.55e-15 kg/s. Raising coarse retention to 0.75
increases mill feed to 52.826087 kg/s and circulating load to 340.217%, while
the external solids, copper and water balances remain closed.

The separate authored Bond estimate uses circuit F80=200 µm, P80=100 µm and
Wi=15 kWh/dry metric tonne: 15816.234 J/kg of FRESH solids, requiring 189.795 kW
net pinion power at 12 kg/s on-stream fresh throughput. It is not multiplied by mill
circulating tonnage, and no additional efficiency effect is inferred from the
redistribution matrix. At retained fraction 1, discrete product P80 equals F80;
the circuit remains reportable but this strict-size-reduction Bond estimate is
unavailable. Discrete passing sizes come only from the declared bins, without
within-bin interpolation. Site fits must justify the fixed coefficients and
their operating domains; this example supplies no plant calibration.

#### Explicit classifier-underflow purge

The same `comminution::circuit::ClosedGrindingCircuit` also supports
`solve_with_underflow_recycle(fresh, mill_conditions, classifier_conditions,
settings, underflow_recycle)`. The final argument is the existing `Partition`:
its product fraction returns to the mill and its reject fraction leaves as a
nonselective purge. Independent fractions preserve a very small positive
bleed even when the returned fraction rounds to one. `solve` remains the
full-return case.

The solution retains the raw classifier underflow, returned stream, purge,
overflow product and original split. It checks the underflow split, the
fresh-plus-recycle mixer, and the combined external boundary independently
for every dry class, declared element and water account. The existing
`product_minus_fresh_*` fields retain their literal overflow-minus-fresh
meaning; the new `external_minus_fresh_*` fields test overflow plus purge
against fresh feed. No purged material is hidden in an artificial balance
adjustment or counted as saleable fine product.

For recycle `s`, purge `t`, classifier underflow `r` and overflow `q`, each
triangular denominator is `t + s*(q + r*fragmented_fraction)`. Water uses
`t + s*water_overflow`. Consequently a positive purge can release solids or
water that otherwise have no external exit, while zero purge retains the
existing trap refusal. Circulating-load limits apply to the actually returned
dry flow, and the original cards still admit the solved operating state.

The existing example accepts a second optional argument, the underflow purge
fraction, default zero. Its declared authored support includes 12..210 kg/s
dry mill flow and solids mass fractions 0.5..0.95 so the complete purge range
can be demonstrated without implying empirical calibration:

```bash
cargo run -p fs-material --example mineral_circuit -- 0.5 0.1
cargo run -p fs-material --example mineral_circuit -- 0.5 1.0
```

Full purge produces one mill pass followed by classification, with no recycle.
The example prints purge and combined external streams separately; a zero
dry-flow stream has an undefined P80. Its optional Bond estimate now uses the
combined overflow-plus-purge P80 with fresh throughput. Using only the fine
overflow would omit coarse material deliberately discharged from the circuit.
For full purge in the authored case, the combined external P80 equals F80 and
the strict-size-reduction Bond estimate is unavailable, even though the
separately displayed fine overflow is smaller.

#### Coupling available power to the admitted circuit

`comminution::circuit::power::BondCircuitPrototype` joins the existing circuit
and Bond owners. Its constructor independently solves a reference feed with
the original mill card, classifier, conditions, numerical settings and exact
underflow return/purge partition. `at_net_pinion_power(Power)` then calculates
fresh capacity, while `at_fresh_throughput(MassFlowRate)` calculates the required
net pinion power. Both use fresh-feed F80 and combined overflow-plus-purge P80.

Every request scales the fresh size classes, elemental accounts and water
together and re-solves the original circuit at the actual operating state.
The original mill throughput and solids-fraction limits still apply; an
admitted reference does not grant permission to extrapolate the source cards.
The result retains the original prototype and exposes the complete scaled
circuit, specific energy, size basis and net pinion power. A changed discrete
size quantile, absent size reduction or unrepresentable positive flow refuses.

The example accepts net pinion power in kW as a third optional argument:

```bash
cargo run -p fs-material --example mineral_circuit -- 0.5 0 100
cargo run -p fs-material --example mineral_circuit -- 0.5 0.01 100
```

The first command uses the authored two-size case above and prints the scaled
fresh, mill-feed, recycle, overflow and purge accounts at 100 kW. The second
retains an explicit 1% underflow purge. Requests above the mill card's supported
flow window fail instead of returning an extrapolated throughput.

These are steady **on-stream** rates. Reporting-period utilization may be
applied separately to production totals; it must not lower the actual running
flow used to admit the mill source. The work index, transfer matrix and
classifier remain caller-supplied empirical inputs. This coupling does not
infer motor or auxiliary loads, ore liberation, wear or power-dependent
breakage kinetics.

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

### Selective rougher/cleaner separation with retained density evidence

`slurry::separation::SelectiveSeparator` accepts a concentrate partition for
every exact dry constituent and an independent carrier-liquid partition.
Both outlet `SlurryStream` values retain their original material sources,
queries and receipts. The result exposes per-component recovery, dry-component
fractions, separate water recovery and balance residuals. Its fractions are
prescribed inputs with provenance; the API does not manufacture a calibrated
flotation, gravity or cyclone correlation from a mineral name.

The runnable `mineral_partition_circuit` example uses the same four compiled
packs above. It sends rougher concentrate to a cleaner and combines rougher
and cleaner tailings through the existing stream mixer. All partitions are
explicitly authored. Its source-backed volume/density calculation gives:

| Stream | Dry solids, kg/s | Water, kg/s | Pulp density, kg/m³ |
| --- | ---: | ---: | ---: |
| Feed | 200 | 100 | 1979.185247 |
| Rougher concentrate | 105 | 20 | 3047.961295 |
| Final concentrate | 87.1 | 5 | 4202.496057 |
| Combined tailings | 112.9 | 95 | 1603.399888 |

The final concentrate contains 85.5 kg/s magnetite, 1.2 kg/s pyrite and
0.4 kg/s dolomite: **98.163031% magnetite by dry constituent mass** at
**85.5% overall magnetite recovery**. Each solid and water balance closes
in the example. These are constituent fractions, not elemental metal assays.
Mineral-by-size IDs can receive distinct partitions without being collapsed.
Missing/extra partitions refuse; zero-solids grade and zero-feed recovery
remain undefined rather than becoming zero or NaN.

```bash
cargo run -p fs-material --example mineral_partition_circuit -- \
  "$mining_pack_dir/magnetite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/pyrite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/dolomite-handbook-mineralogy-reference.fsmatpk" \
  "$mining_pack_dir/water-liquid-iapws-sr6-08.fsmatpk"
```

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

## Continuous treatment driven by actual slurry flow

`recovery::continuous::ContinuousRecoveryCard` connects actual flow, active
mixed volume and a separately declared stationary local kinetic law. For a
single rate and N equal full ideal tanks, the recovered feed fraction is
`Rinf * (1 - (1 + k*V/(N*Q))^(-N))`. The uniform-rate law analytically
averages survival over `[0,kmax]`, preserving the same particle rate through
all tanks. It never substitutes batch recovery evaluated at the mean residence.

The two law IDs, `cstr-single-rate-recovery` and `cstr-uniform-rate-recovery`,
require their own original source hashes and provenance, typed finite
`mean-residence-time` and `tank-count` domains, and all other conditions.
They reject a finite batch `elapsed-time` domain. A batch fit is therefore not
silently extrapolated over an infinite residence-time distribution. These
cards do not replace the separate transient inventory owner or the explicit
discrete-population card above. Zero flow refuses a stationary treatment result.

The API consumes the slurry model's existing typed result directly:

```rust
let treatment = model.evaluate_train(
    8,
    active_mixed_volume,
    slurry.state().total_volume_flow,
    contained_metal_feed,
    &conditions,
)?;
```

The response retains the complete source card, admitted conditions, actual
volume/flow, derived mean residence, recovered and unrecovered fractions,
feed/recovered/tail contained-metal rates, and their numerical closure.
An assay is an independent input; a mineral's reference density does not
supply a metal grade or predict its recovery.

The existing `mineral_slurry` command now continues its source-backed Bond,
water-addition and volume calculation into an authored eight-tank treatment
example. Active volume is fixed at 17000 m³, pH is 10, and contained-metal
grade is an independently authored 3 g/dry metric tonne. Both laws have
Rinf=0.9 and mean local rate 1/(6 h); the uniform population spans 0..2/(6 h).

| Solids mass fraction | Actual slurry flow, m³/s | Mean residence, h | Single-rate recovery | Uniform-rate recovery |
| --- | ---: | ---: | ---: | ---: |
| 0.65 | 0.195356574114 | 24.172323064 | 86.554859543% | 77.331774645% |
| 0.50 | 0.311082927870 | 15.179946565 | 80.010623676% | 70.330265324% |

This holds dry throughput and the authored contained-metal feed at 250 kg/s
and 2.7 kg/h. Additional water increases actual volume flow and reduces the
residence/recovery predicted by these declared laws. At 10% solids the slurry
calculation remains valid, but treatment reports `OutsideValidity`: the
residence has left the separately declared 12..48 h support. No value is clipped
or extrapolated. Carrier-flow changes from dissolution, adsorption, reagent
limits, transients and site calibration require their own models.

A separate small example accepts actual flow as an optional m³/s argument:

```bash
cargo run -p fs-material --example continuous_recovery
cargo run -p fs-material --example continuous_recovery -- 0.097678285
```

At its authored 24 h nominal residence, eight tanks give 86.488340192% recovery
for the single-rate card and 77.243303571% for the equal-mean uniform card.
Half flow gives 48 h and respectively 89.648437500% and 83.574368019%.
These source-bound arithmetic scenarios are not calibrated plant forecasts.

## Closed rougher-cleaner recycle and purge

`fs_material::slurry::recycle::RougherCleanerCircuit` closes the cleaner-tail
return to the rougher feed. It retains both existing `SelectiveSeparator`
prescriptions and every original constituent density query and receipt. Fresh
feed, mixed rougher feed, both stage results, returned cleaner tail, purge,
final concentrate, and combined external tailings remain inspectable.

For each exact dry-component identity, and separately for the liquid carrier,
let `r` and `c` be rougher and cleaner fractions to concentrate and `s` the
nonselective cleaner-tail return. The steady mixed feed is
`F / ((1-r) + r*(c + (1-c)*(1-s)))`. Summing positive exit paths preserves
small product exits that subtracting a rounded loop gain from one would lose.
The existing separators, splitter and mixer then replay the solved circuit.
Both `fresh = final product + rougher tails + purge` and
`mixed = fresh + recycle` are checked independently; no balance plug is added.

The real-pack example uses the same four density packs as the open circuit:

```bash
cargo run -p fs-material --example mineral_recycle -- \
  /tmp/mining-materials/magnetite-handbook-mineralogy-reference.fsmatpk \
  /tmp/mining-materials/pyrite-handbook-mineralogy-reference.fsmatpk \
  /tmp/mining-materials/dolomite-handbook-mineralogy-reference.fsmatpk \
  /tmp/mining-materials/water-liquid-iapws-sr6-08.fsmatpk 0.8
```

Its authored fresh feed is 100 kg/s magnetite, 20 kg/s pyrite, 80 kg/s
dolomite and 100 kg/s water. Rougher fractions are 0.95/0.30/0.05 for the
minerals and 0.20 for water; cleaner fractions are 0.90/0.20/0.10 and 0.25.
With 80% cleaner-tail return, final product contains respectively
92.5324675325, 1.4851485149 and 0.4149377593 kg/s of dry constituents plus
5.6818181818 kg/s water. Magnetite recovery is 92.5324675325%, compared with
85.5% for the same two stages without return. This is a consequence of the
declared constant partitions, not a calibrated recovery improvement claim.

The solver reports each component's circulating load, overall recovery, and
both residuals. Missing partitions, fully trapped identities, nonfinite
flows and erased positive paths refuse. A fully trapped zero-feed identity
also refuses because its stationary circulating inventory is undetermined.
External closure is scaled by fresh flow, not the inflated internal flow;
near-trapped ill-conditioned circuits can therefore refuse. These are steady
prescribed partitions: changing liberation, capacity, reagent conditions or
residence with recycle requires separately qualified models.

## Source-backed aqueous calcite acid demand and heat

The three aqueous USGS profiles for H+, Ca2+ and HCO3- complete the existing
calcite bicarbonate reaction with a compatible reference set. The executable
`mineral_processing --reaction calcite-bicarbonate-neutralization` accepts
calcite, hydrogen-ion, calcium-ion and bicarbonate-ion packs in that order,
followed by kilograms of actually reacted pure calcite. The standard result
is -25.5 kJ/mol extent at 298.15 K and 1 bar. Exact atoms/charge, source
molar-mass conversions and all four formation-enthalpy receipts remain
inspectable. The associated integration test also closes an explicitly
acid-limited inventory and refuses reagent overdraw.

See [the aqueous reference table and runnable example](MINERAL_THERMOCHEMISTRY.md#aqueous-references-for-calcite-neutralization)
for the 13 new claims, conventional ionic entropy, source uncertainty
magnitudes, and net-ionic versus commercial-reagent distinction. These
reference-state quantities do not select process pH, predict dissolution
kinetics, or substitute for the separate dissolved-CO2 endpoint.

## Aqueous pyrite oxidation and the iron redox cycle

The compatible USGS pyrite, oxygen, liquid-water, ferrous-ion and sulfate
profiles complete `PyriteOxidationFerrous` with the existing hydrogen-ion
source. `mineral_processing --reaction pyrite-oxidation-ferrous` reports the
six-source amount, mass and heat balance. One mole of pyrite at this declared
endpoint consumes 3.5 mol oxygen and one mol water, and produces one mol Fe2+,
two mol sulfate and two mol H+. The standard heat is -1452.4 kJ/mol pyrite.
See [the source values, uncertainty conventions and runnable pack commands](MINING_ACID_GENERATION.md).

The ferric-ion profile completes `FerrousOxidationFerric` and
`FerricPyriteOxidation`, also exposed as `ferrous-oxidation-ferric` and
`ferric-pyrite-oxidation` CLI modes. Their written-equation heats are -406.8
and -28.6 kJ/mol respectively. Seven ferrous-oxidation extents plus two
ferric-attack extents recover the original two-pyrite net equation and heat.
Finite-stage tests restore an explicit ferric inventory and reject attempts
to spend future reaction products before they exist. Actual rates, aqueous
speciation and precipitation still require separately admitted models.

## Source-bound equilibrium constants and chemical driving force

The `stoichiometry::thermodynamics::caloric` source route exposes
`ConservedReaction::standard_equilibrium`, which uses the same explicit
`FormationEnthalpySource` associations as the reaction-heat calculation and
also resolves entropy from each original pack. It returns standard reaction
enthalpy, entropy, Gibbs energy and the natural logarithm of the dimensionless
equilibrium constant at the common source temperature and pressure:

`Delta G° = Delta H° - T Delta S°`, `ln K = -Delta G°/(R T)`.

These relations follow the [IUPAC Gibbs-energy definition](https://goldbook.iupac.org/terms/view/G02629)
and [standard equilibrium constant](https://goldbook.iupac.org/terms/view/S05915).
The original enthalpy and entropy rounding is retained. Independently rounded
formation-Gibbs tables need not yield an identical number. The adjacent
`standard_gibbs_change` method consumes direct formation-Gibbs claims; neither
method substitutes a species formation enthalpy minus its absolute entropy
for a formation-Gibbs datum. The caloric route combines exactly conserved
reaction sums, so elemental and common ionic reference contributions cancel.

Charged aqueous species require `conventional-standard-molar-entropy` and
source-declared hydrogen-zero, individual-ion and unit-activity/infinite-dilution
qualifiers. Negative conventional ionic entropies are retained. The selected
enthalpy and entropy claims must agree on phase, polymorph and solution
convention. Every input receipt and uncertainty remains inspectable; no
propagated confidence interval or temperature extrapolation is added.

`ReactionEquilibrium::at_log_activities` accepts one explicit natural-log
activity per active species and evaluates `Delta G = Delta G° + R T ln Q`.
The coefficients come from the original conserved equation. Unit-activity
phases are supplied explicitly with log activity zero. The API refuses missing,
duplicate, nonfinite or unrepresentable inputs and does not infer activities
from concentrations, solve speciation, choose pH or predict reaction rates.
Logarithmic K remains available when exponentiating it would exceed binary64;
the optional `equilibrium_constant()` conversion refuses that range loss.

Prefix any supported source-pack CLI invocation with `--equilibrium`:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --reaction pyrite-oxidation-ferrous \
  "$acid_pack_dir/pyrite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$acid_pack_dir/oxygen-gas-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/iron-ii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/sulfate-ion-aqueous-usgs-b2131.fsmatpk" \
  "$acid_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 119.979
```

Use the pack-building loop in the acid-generation guide first. The selected
rows give -1822.05 J/(mol K) reaction entropy and -2361555.7925 J/mol standard
Gibbs energy per written two-pyrite equation. These are reference-state
quantities; process-water activities and finite conversion remain separate.

## Finite-inventory electrowinning

`stoichiometry::electrolysis::ElectrolysisReaction` connects the existing
Faraday amount identity to an original conserved reaction and finite reagent
stocks. `advance_charge` accepts integrated charge; `advance_current` accepts
constant current and elapsed time. Both require an explicit current efficiency,
deposited product, molar mass and electron count. They advance every modeled
species atomically and refuse charge-driven reagent overdraw without clipping
the requested efficiency.

`MineralReaction::CopperElectrowinning` supplies the net productive endpoint
`2 Cu2+(aq) + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+(aq)`: two electrons per
deposited copper atom and four per written equation. The result retains the
model, passed/productive/unallocated charge, extent, product amount/mass and
modeled final inventory. Unallocated charge has unspecified side chemistry;
its water, oxygen and proton effects are not inferred.

```bash
cargo run -p fs-thermochem --example electrowinning
cargo run -p fs-thermochem --example electrowinning -- \
  --charge 3600000 0.95 1000 10000
```

Both authored cases deposit approximately 1.126219474 kg copper using an
explicit 0.063546 kg/mol mass basis. The example also reports the remaining
dissolved copper, modeled water/oxygen/proton inventory and unallocated
charge. Cell voltage, electrical energy, transport limitations, efficiency
calibration and side-reaction balances require their own inputs and models.

## Focused verification

Normal workspace commands for these slices:

```bash
cargo test -p xtask --test mining_mineral_sources --test mining_mineral_coverage --test mining_reagent_sources --test mining_reagent_references --test mining_oxide_sources --test mining_copper_sources --test mining_reaction_balances
cargo test -p fs-material comminution
cargo test -p fs-material slurry::
cargo test -p fs-material recovery::
cargo test -p xtask --test mining_slurry_sources
cargo test -p xtask --test mining_aqueous_sources
cargo test -p xtask --test mining_acid_generation --test mining_iron_redox --test mining_equilibrium
cargo test -p fs-thermochem stoichiometry::tests
cargo test -p fs-thermochem stoichiometry::electrolysis
```

The source/consumer tests exercise the actual compiler, canonical decoding,
typed query receipts, unsupported-condition refusals, known-answer physics,
deterministic ordering and source-reference mutations. They establish these
bounded data/formula capabilities, not site calibration or full-workspace
qualification. The source-backed command-line example is an additional real
consumer of the shipped packs.

