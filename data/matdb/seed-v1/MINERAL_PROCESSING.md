# Mineral processing: public sources and callable physics

These **18 source profiles and 56 property claims** provide reusable mineral
and reagent inputs for FrankenMining's material-corpus and Crucible work.
They use the existing `matdb-pack` compiler, `NormalizedPack` queries and
`matdb-catalog` discovery. No separate material registry or schema is required.

The physical consumers are `fs_material::comminution` and
`fs_thermochem::stoichiometry`. These are upstream building blocks: they do
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

`MineralReaction::build` provides ten explicit endpoint choices:

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

Pyrite equations use integer multiples of two/four FeS2; an extent is per
equation as written. These endpoints follow the primary references cited in
`stoichiometry/minerals.rs`
and do not choose pH/Eh, speciation or precipitates. Only lime hydration has
complete formation-enthalpy inputs in this new seed slice. Missing aqueous-ion
data for other reaction heats remain missing.

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

`ball_mill_critical_speed` returns the ideal single-ball centrifuging frequency
from inside-liner diameter, ball diameter and explicit gravity:
`n = sqrt(2*g/(D-d))/(2*pi)` revolutions/second.
It is not an operating-speed recommendation, power-draw model or DEM solver.

## Focused verification

Normal workspace commands for these slices:

```bash
cargo test -p xtask --test mining_mineral_sources --test mining_mineral_coverage --test mining_reagent_sources --test mining_reagent_references --test mining_reaction_balances
cargo test -p fs-material comminution::tests
cargo test -p fs-thermochem stoichiometry::tests
```

The source/consumer tests exercise the actual compiler, canonical decoding,
typed query receipts, unsupported-condition refusals, known-answer physics,
deterministic ordering and source-reference mutations. They establish these
bounded data/formula capabilities, not site calibration or full-workspace
qualification. The source-backed command-line example is an additional real
consumer of the shipped packs.
