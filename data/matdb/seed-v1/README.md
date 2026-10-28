# FrankenSim material seed data v1

This directory retains the human-reviewable source inputs compiled by
`xtask matdb-pack`. Generated binary packs are deliberately not committed:
the compiler must reproduce them from the pinned manifest and source record.

## Aerogel and refractory-fibre insulation

Four packs add Pyrogel XTE revision 4.2 and Superwool Plus blanket at nominal
96, 128 and 160 kg/m³. Each holds density, a source conductivity table, one
mechanical-strength observation and a separate authored conductivity model.

| Pack | Source temperature range | Nominal density | Mechanical source value |
| --- | --- | ---: | --- |
| `aspen-pyrogel-xte-rev42` | Mean 0–600 °C | 200 kg/m³ | At least 3 psi compressive resistance at 10% deformation, with 2 psi preload |
| `morgan-superwool-plus-96-2025` | Mean 200–1000 °C | 96 kg/m³ | Average tensile strength 50 kPa |
| `morgan-superwool-plus-128-2025` | Mean 200–1000 °C | 128 kg/m³ | Average tensile strength 75 kPa |
| `morgan-superwool-plus-160-2025` | Mean 200–1000 °C | 160 kg/m³ | Average tensile strength 90 kPa |

The Pyrogel source is the
[Aspen producer sheet, revision 4.2](https://www.aerogel.com/wp-content/uploads/2021/06/Pyrogel-XTE-Datasheet-English-v4.2.pdf).
The Superwool source is
[Morgan BL.01, 18 April 2025](https://www.morganthermalceramics.com/media/l52g3r3j/superwool-plus-blankets_eng.pdf).
Both source tables describe effective blanket transport at mean test temperature.
ASTM's [C201 method summary](https://store.astm.org/standards/c201) identifies the
large-gradient, steady mean-temperature measurement used by Superwool.

The measured table role is `apparent-thermal-conductivity`, with typed
`mean-test-temperature` and tabulated-only evaluation. It refuses between-point
queries. The separate `thermal-conductivity` curve uses linear interpolation
and the local `temperature` coordinate only when
`engineering-reference-model=1` and
`mean-test-temperature-as-local-model=1` are supplied. This is an explicit
engineering approximation, with no measured-local-law or interpolation-error
claim. Its effective transport already aggregates the source's heat-transfer
mechanisms; no separate pore-radiation contribution is inferred.

All thermal queries retain the exact grade, unknown moisture and gas pressure,
unknown shared-specimen history, direction, test temperature difference and
specimen thickness. Superwool also retains unknown compressive test load.
Pyrogel conductivity instead requires the source's **2 psi compressive stress**,
stored as `test-compressive-stress=13789.514586336722 Pa`; this is distinct
from ambient gas pressure and from the mechanical test's
`test-preload-pressure`.

No heat capacity is supplied, so these packs support bounded steady-transport
calculations and do not complete a transient-heating profile. Nominal density
does not supply a density-temperature law. Service-temperature or melting
ratings are not promoted to conductivity-domain extensions or latent heat.

### Run the insulation data through a nonlinear wall

The existing `material_wall` example accepts arbitrary source packs in a
`nonlinear` mode. It retains every source knot and query receipt, then calls
the existing `LayeredWall::solve_temperature_dependent` owner. The original
three-path CFAST invocation remains available.

For these four insulation packs, choosing the engineering curve is explicit:
both `engineering-reference-model=1` and
`mean-test-temperature-as-local-model=1` must be supplied. The producer's
mean-test-temperature measurements remain separate tabulated-only claims.
The following are authored steady planar calculations using that local-k
approximation; they are not measurements of a 50 mm product assembly.

From the repository root, in Bash:

```bash
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/aspen-pyrogel-xte-rev42/manifest.tsv --out /tmp/pyrogel-xte.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/morgan-superwool-plus-128-2025/manifest.tsv --out /tmp/superwool-128.fsmatpk

wall_axes=(
  --axis engineering-reference-model=1
  --axis mean-test-temperature-as-local-model=1
  --axis source-moisture-condition-known=0
  --axis source-gas-pressure-known=0
  --axis source-specimens-matched=0
  --axis source-thermal-direction-known=0
  --axis source-test-temperature-difference-known=0
  --axis source-test-thickness-known=0
)

cargo run -p fs-conduction --example material_wall -- nonlinear \
  --left-k 873.15 --right-k 473.15 --area-m2 0.01 \
  --layer /tmp/pyrogel-xte.fsmatpk 0.05 "${wall_axes[@]}" \
  --axis source-grade-aspen-pyrogel-xte-rev42=1 \
  --axis test-compressive-stress=13789.514586336722

cargo run -p fs-conduction --example material_wall -- nonlinear \
  --left-k 873.15 --right-k 473.15 --area-m2 0.01 \
  --layer /tmp/superwool-128.fsmatpk 0.05 "${wall_axes[@]}" \
  --axis source-grade-morgan-superwool-plus-128-2025=1 \
  --axis source-compressive-load-known=0
```

The CLI types each source axis from its claim declaration. In particular,
Pyrogel's `test-compressive-stress` is pressure in Pa, preserving the source's
2 psi compressive test load; it is not ambient gas pressure. The Superwool
source leaves compressive test load unknown.

For the declared 873.15 K / 473.15 K face pair, 0.05 m thickness and
0.01 m² area, independent piecewise integration gives:

| Local-k reference | Integral of k from 473.15 to 873.15 K, W/m | Heat flux, W/m² | Heat rate, W | Midplane, K |
| --- | ---: | ---: | ---: | ---: |
| Pyrogel XTE | 20.35 | 407 | 4.07 | 731.2251150283 |
| Superwool Plus 96 | 37 | 740 | 7.40 | 720.0897989875 |
| Superwool Plus 128 | 33 | 660 | 6.60 | 714.7380433164 |
| Superwool Plus 160 | 29 | 580 | 5.80 | 720.0126966597 |

These values follow from `q = integral(k(T), Tc..Th) / L`, using the stated
linear engineering interpolation. They do not replace source uncertainty.
To use another Superwool density, compile its own manifest and supply the
matching grade axis; density classes are not interchangeable.

The model's nonlinear midplane can also be exposed as two identical 25 mm
layers with an explicitly zero contact resistance:

```bash
cargo run -p fs-conduction --example material_wall -- nonlinear \
  --left-k 873.15 --right-k 473.15 --area-m2 0.01 \
  --layer /tmp/pyrogel-xte.fsmatpk 0.025 "${wall_axes[@]}" \
  --axis source-grade-aspen-pyrogel-xte-rev42=1 \
  --axis test-compressive-stress=13789.514586336722 \
  --contact 0 \
  --layer /tmp/pyrogel-xte.fsmatpk 0.025 "${wall_axes[@]}" \
  --axis source-grade-aspen-pyrogel-xte-rev42=1 \
  --axis test-compressive-stress=13789.514586336722
```

The shared layer face is approximately 731.225115 K, with unchanged heat
rate. Optional `--films LEFT RIGHT` and every `--contact R` are explicit
scenario resistances in m² K/W. Omitted films prescribe the material face
temperatures. Each layer owns the `--axis` entries following it; the next
layer requires a contact declaration, including zero for ideal contact.

Output includes pack/card identities, typed source-query receipts, layer
face temperatures, signed left-to-right flux, total heat rate, secant R/U
for that boundary pair, and numerical Fourier/boundary residuals. These
residuals diagnose the solve and are not an error bound for the source data.
Every layer must support the full reservoir interval. Omitted model
selectors, wrong compressive load, undeclared axes and out-of-domain
temperatures are refused before a successful result is printed. This
example accepts distinct positive reservoir temperatures, up to 16 layers,
64 axes per layer, 512 arguments, 128 KiB argument bytes, and 16 MiB per
source pack. No heat capacity or transient behavior is inferred.

## Cured electronics-potting resins: MG Chemicals 2026

Five producer grades add 34 claims for cured thermosets: 834HTC and 832TC thermally conductive epoxies, 832C translucent epoxy, 832FX flexible epoxy, and 8810 rigid polyurethane. Each pack contains cured density, thermal conductivity, specific heat capacity, independently reported thermal diffusivity, tensile strength, and electrical resistivity. Four also contain compressive strength; 832FX does not report it.

| Pack | Claims | Producer sheet |
| --- | ---: | --- |
| `mg-chemicals-834htc-cured-2026` | 7 | [834HTC, 23 January 2026, Ver. 6.4](https://mgchemicals.com/downloads/tds/tds-834htc-2parts.pdf) |
| `mg-chemicals-832tc-cured-2026` | 7 | [832TC, 04 February 2026, Ver. 4.3](https://mgchemicals.com/downloads/tds/tds-832tc-2parts.pdf) |
| `mg-chemicals-832c-cured-2026` | 7 | [832C, 20 February 2026, Ver. 4.2](https://mgchemicals.com/downloads/tds/tds-832c-2parts.pdf) |
| `mg-chemicals-832fx-cured-2026` | 6 | [832FX, 23 January 2026, Ver. 4.3](https://mgchemicals.com/downloads/tds/tds-832fx-2parts.pdf) |
| `mg-chemicals-8810-cured-2026` | 7 | [8810, 23 January 2026, Ver. 5.2](https://mgchemicals.com/downloads/tds/tds-8810-2parts.pdf) |

The three thermal properties have point validity at **298.15 K**. Cured density is explicitly separate from liquid-mixture density and retains an unknown measurement-temperature flag. The sources specify product mix ratios and several possible cure schedules; the actual cure used for the retained property tests is unstated. Grade, cured state, specified mass mix ratio, unknown cure schedule, moisture condition, pressure, and specimen matching remain typed validity conditions. Unknown test temperatures, measurement directions, and loading rates remain explicit where relevant. Every imported claim has unstated uncertainty. Service-temperature ratings and curing temperatures do not extend thermal validity.

The source diffusivity remains an independent observation. Rounded source conductivity, density, and heat capacity do not exactly reproduce it, particularly for 832C and 8810. Consumers must identify a diffusivity computed from `k/(rho*Cp)` as a derived value and preserve its source dependencies; it must not overwrite the reported claim or imply matched specimens. These packs provide measured strength roles without inventing elastic moduli, Poisson ratios, viscoelasticity, curing kinetics, or properties over an operating-temperature window.

## Printed polymers, compliant materials, ordinary steels and ceramics

This expansion adds **23 source packs and 193 property claims**, including
additive improvements to five existing alloy/mineral profiles. New data are
automatically searchable and inspectable through `matdb-catalog`.

| Family | New packs | Added claims | Coverage |
| --- | ---: | ---: | --- |
| Printed plastics | 3 | 32 | Prusament PLA/PETG MK3 profiles; Ultrafuse TPU 95A, with orientation, print settings and conditioning |
| Named rubbers and casting resins | 5 | 24 | Parker E3609 EPDM and N0674 NBR; Ecoflex 00-30, Dragon Skin 20 and VytaFlex 40 |
| Flexible polyurethane foam | 3 | 27 | PORON 4701-40 at nominal 240, 320 and 480 kg/m³ |
| Ordinary steels | 2 | 20 | A36 hot-rolled and 1018 cold-finished supplier references with density, stiffness, heat capacity, conductivity, strength and ductility |
| Stainless steels and elements | 4 | 20 | 430/1.4016 and 316L/1.4404; bismuth and vanadium |
| Ceramics and mineral crystals | 6 | 57 | AD-96/AD-995 alumina, AlN 180, Hexoloy SA SiC, MACOR's 2013 source, sylvite; frequency-qualified permittivity added to existing halite and periclase |
| Existing alloy completion | 0 | 13 | INCONEL 625 dynamic elastic curves; MONEL 400 nominal moduli; TIMETAL 6-4 endpoint-mean expansion |

The printed profiles distinguish orientation and actual test protocol. The TPU
sheet's filament density remains separate from printed bulk density. A rubber
value labeled “100% modulus” is retained as **tensile stress at strain 1**;
it does not supply a small-strain elastic or hyperelastic law. Silicone and
urethane castings retain seven-day cure and mixing conditions. Tear strengths
keep the original Die B/Die C geometry. Compression-set observations keep the
temperature, dwell and imposed deformation distinct from a recoverable strain.

PORON retains its 25% compression point and reported crosshead speed. Typical
values, minimum/maximum limits and nominal density remain separate claims.
No foam receives an invented heat capacity, crushing curve or hysteresis law.

A36 and 1018 are explicitly acknowledged supplier-guide engineering references.
The source supplies no numeric test temperature, matched coupon or statistical
uncertainty. This supports a declared axial or conductive-slab reference, but
does not admit a temperature-span heat-storage law or structural design
allowable. Imperial units are converted explicitly; the unspecified Btu
convention is recorded as an International Table conversion assumption.
Independently reported E, G and Poisson ratio are not forced into consistency.

For example, resolve the A36 reference's available mechanical and thermal
properties together, retaining all unknown-condition acknowledgements:

```sh
cargo run -p xtask -- matdb-catalog --pack steel-a36-onlinemetals-reference \
  --property density --property young_modulus --property poisson_ratio \
  --property specific-heat-capacity --property thermal-conductivity \
  --at engineering-reference-model=1 --at source-specimens-matched=0 \
  --at source-test-temperature-known=0 --at source-test-method-known=0 \
  --at source-hot-rolled=1
```

The ceramic and elemental packs preserve differences between heat-capacity,
conductivity and mechanical test temperatures. Unspecified “room temperature”
does not acquire an invented Kelvin value. Dielectric observations retain their
frequencies, sparse frequency/temperature samples remain discrete where the
source supplies no interpolation law, and mean expansion retains both interval
endpoints. The existing halite and periclase thermal claims remain intact.

INCONEL 625 provides separate annealed and solution-treated dynamic E/G/nu
curves from 294.15 to 1144.15 K, preserving the 0.75-inch rod and unknown
frequency. MONEL tension, compression, torsion and Poisson data remain
independent nominal constants with an explicitly authored calculation reference.
TIMETAL adds three interval means; no instantaneous expansion curve is inferred.

Focused source/query and physical-consumer regressions live in the
`material_catalog_expansion` target's `printed_elastomers`, `ceramics_elements`
and `alloy_completion` modules. They exercise source values, canonical hashes,
receipt replay, rejected conditions, axial force/mass/work, slab resistance and
endpoint expansion using the production material owners.


## Compression observations and engineering ceramics

This tranche adds **12 source packs, 100 property claims and 18 tabulated
datasets containing 63 source points**. Source observations, specification
limits and authored engineering approximations keep separate meanings.

| New source family | Packs | Coverage |
| --- | --- | --- |
| Plastazote LD33/LD45 and Evazote EV30/EV50 | 4 | Density, discrete compression stresses, tension/tear properties and compression-set recovery; LD33/LD45 also carry conductivity at the stated mean temperature |
| PORON 4701-30, 320 kg/m3; BISCO HT-800 | 2 | Named polyurethane and silicone foams; typical values and specification limits remain distinct |
| ATI Zircadyne 702/705 | 2 | Temperature-dependent modulus and tensile/proof strength, with process and direction-average selectors |
| Corning MACOR, Kyocera SN241O/SN240O, Morgan CVD SiC | 4 | Grade-specific thermal/mechanical properties, interval mean expansion and explicit frozen-reference assumptions |

The Zotefoams June 2026 sheets supply compression stresses at 25% and 50%,
plus compression-set observations after 30 minutes and 24 hours of recovery.
The compression-set records retain the preceding compression, hold time and
temperature. These are discrete observations, not a recovery differential
equation or a calibrated hysteresis/relaxation law. Density is explicitly
skin-to-skin. Unstated specimen geometry and loading rate remain unknown.

PORON selects the 320 kg/m3 column and its reported commercial thickness band.
Its 25% compression point preserves the source's 0.51 cm/min crosshead speed.
BISCO's August 2026 sheet leaves numeric compression deflection unstated, so its
pressure-valued `compression-force-deflection` is not relabeled as a 25% stress
point. Its less-than-5% compression-set specification remains an exclusive
limit. No conductivity is inferred for the Evazote grades and no foam receives
an invented heat capacity.

Zircadyne's numbered temperature rows are separate from its unspecified
room-temperature observations; directional averages and cold-worked/annealed
strength conditions do not become a universal isotropic law. The four ceramic
thermal references require explicit `engineering-reference-model=1` and
`source-temperatures-matched=0`. Their one-kelvin frozen intervals are authored
approximations, not measured temperature dependence. Mean expansion retains
both interval endpoints. Morgan's high-resistivity SiC is separate from its ELR
grade, and its greater-than resistivity bound is not an exact resistivity.

### Use compression observations in a bounded specimen calculation

`fs-material::compression_reference` and the
`material_compression_reference` example turn a selected nominal
compression-stress table into force, shortening, segment stiffness and model
loading work for explicit original geometry. Every original source knot and
receipt survives. By default, the source must permit continuous interpolation.
`--interpolate-observations` explicitly authors a piecewise-linear reference
between discrete source points; ordinary material queries still reject
unreported intermediate strains.

A specimen begins at the **first source strain**, with zero subsequent loading
work. No unloaded origin or missing loading segment is invented. Original
density, length and area determine fixed reference mass. All loading must
increase compression inside the original strain span; a known source crosshead
speed must match each ramp. Unloading, holds, recovery, creep, failure, lateral
contraction, dissipated heat and recoverable energy are outside this reference.
The reported tangent is the right segment at an interior knot and the left
segment at the final endpoint.

For a 10 mm long specimen with 100 mm2 original area, the new foam references
give the following independent expectations at 50% compression. Work is only
the integral from the first supported 25% compression:

| Grade | Reference mass (kg) | Force at 50% (N) | Model work from 25% to 50% (J) |
| --- | ---: | ---: | ---: |
| Plastazote LD33 | 0.000033 | 13.4 | 0.025 |
| Plastazote LD45 | 0.000045 | 16.4 | 0.031875 |
| Evazote EV30 | 0.000030 | 10.0 | 0.017375 |
| Evazote EV50 | 0.000050 | 11.6 | 0.02075 |

Compile and run LD33:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/zotefoams-plastazote-ld33-2026/manifest.tsv \
  --out /tmp/ld33.fsmatpk
cargo run -p fs-material --example material_compression_reference -- \
  /tmp/ld33.fsmatpk density compressive-stress compressive-strain 0.01 0.0001 \
  --axis source-grade-plastazote-ld33-2026 1 \
  --axis density-includes-skins 1 \
  --axis source-temperature-known 0 \
  --axis source-moisture-condition-known 0 \
  --axis source-pressure-known 0 \
  --axis source-orientation-known 0 \
  --axis source-test-rate-known 0 \
  --axis source-specimen-geometry-known 0 \
  --axis compressive-strain 0.25 \
  --interpolate-observations \
  --step 0.00375 1 --step 0.005 1
```

The one-second foam ramps are authored durations because the source rate is
unknown. The existing `ensinger-tecapeek-gf30-natural` pack supplies three
compression points at 1%, 2% and 5% and a known test speed. With the same original
geometry, its 5% reference gives 12,000 N and 2.985 J of model loading work
**since 1% compression**, with mass 0.00153 kg. Select
`test-speed=0.00008333333333333333`, `--speed-axis test-speed`, the grade and
unknown-condition selectors, and `compressive-strain=0.01`. The monotonic
steps `--step 0.0002 1.2 --step 0.0005 3.6` reproduce the source speed.

Focused regression targets are `cargo test -p fs-material --test
compression_reference` and the existing `xtask`
`material_catalog_expansion` target's `cushioning`, `thermal_engineering` and
`compression_reference` modules. They cover original source values and
conditions, independent piecewise work, exact source receipts, claim pins,
step subdivision, deterministic replay and refusal before history mutation.


## Magnesium and added elastic data for nickel alloys

This expansion adds five AZ31B claims and completes two existing nickel-alloy
packs with three independently reported moduli. Their existing thermal,
electrical and expansion records remain intact.

| Pack | Addition | Source scope |
| --- | --- | --- |
| `magnesium-az31b-luxfer-2018` | Density, Young modulus, Poisson ratio, heat capacity, conductivity | AZ31B rolled sheet/plate/coil; physical-table temper, temperature, method and direction unstated |
| `inconel-718-aged-special-metals-2007` | Dynamic Young and torsional moduli | SMC-045 Table 4; hot-rolled-flat specimens and their own test/treatment context |
| `nilo-36-special-metals-2004` | Young modulus | SMC-031 Table 16; temperature, method, direction and temper unstated |

AZ31B retains the printed SI modulus of 44 GPa and records the inconsistent
6500-ksi companion without silently substituting it. The inconsistent expansion
entry is omitted. Conditions from the separate strength table are not imported.

INCONEL 718 uses `dynamic-young-modulus` and `dynamic-torsional-modulus`.
They remain dynamic bending/torsion results at the native 70 °F point, with
individual test frequencies unknown. They do not become static axial inputs
or an independent three-constant elastic calibration. The printed derived
Poisson ratio is not added as a new measurement. The existing
`engineering-reference-model=1` and `source-specimens-matched=0` requirements
also apply to the new moduli, preserving the distinction between the elastic
and transport specimens.

NILO 36 gains only its source-reported 140 GPa Young modulus. Its existing
mean-expansion curve still uses endpoint temperature relative to the required
20 °C reference; it is not an instantaneous expansion law. No independent shear
modulus or Poisson ratio is supplied.

The axial reference accepts `young-modulus` or `young_modulus` for **tension**,
so an E-only metal can supply an authored one-dimensional reference without
inventing a shear modulus or Poisson ratio. The original key, source receipt
and unknown protocol/direction remain intact. Compression still requires an
independently named compression modulus; dynamic, flexural, strength and
point-stress substitutions refuse.

For example, use NILO 36 with authored 50 mm / 100 mm² geometry:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/nilo-36-special-metals-2004/manifest.tsv --out /tmp/nilo-36.fsmatpk
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/nilo-36.fsmatpk tension density young-modulus 0.05 0.0001 0.001 \
  --axis engineering-reference-model 1 --axis source-specimens-matched 0 \
  --axis temperature 293.15 --axis source-grade-nilo-36-special-metals-2004 1 \
  --axis source-pressure-known 0 --axis source-elastic-temperature-known 0 \
  --axis source-elastic-method-known 0 --axis source-elastic-direction-known 0 \
  --axis source-elastic-temper-known 0 --step 0.00005 1
```

The analytic reference is 14,000 N, 0.04055 kg and 0.35 J of model loading work
for E = 140 GPa and density = 8110 kg/m³. These follow from the authored linear
reference, not measured coupon results. The density's 293.15 K reference,
duration and strain ceiling are scenario choices; the modulus's test
temperature remains unknown. No tensile-test protocol, matched-specimen
calibration, yield margin, unloading response or three-dimensional isotropy is
established.

## Engineering grades, adhesives and building materials

This expansion adds **26 source packs and 210 property claims**, including
14 tabulated datasets with 78 source knots. All packs are discovered by
`matdb-catalog`; there is no separate registry to update. Each manifest pins its
producer/source edition and each claim retains typed units, source observations,
selection conditions and unstated uncertainty where the source supplies none.

| Family | New materials | Packs | Claims |
| --- | --- | ---: | ---: |
| Elements | Indium wire, palladium, platinum, iridium | 4 | 33 |
| Metal alloys | MONEL 400, NILO 36, TIMETAL 6-4, aged INCONEL 718 | 4 | 16 |
| Optical crystals | Lithium fluoride; CVD FLIR-grade zinc sulfide | 2 | 8 |
| Engineering plastics | TECAPEEK GF30; TECATRON GF40 and SX; TECAPAI CM XP403 | 4 | 30 |
| Structural foams | Divinycell H60, H80, H130 and H200 PVC | 4 | 56 |
| Resins and adhesives | LOCTITE EA 9466; EpoxAcast 655/103; Resbond 940 zircon ceramic adhesive | 3 | 19 |
| Reinforced laminate | Norplex NP509 melamine/woven glass | 1 | 11 |
| Concrete and cement | QUIKRETE 5000; CTS CEMENT ALL | 2 | 16 |
| Wood panels | MEDITE MR 18 mm MDF; WISA UPM007CPR 18 mm, 13-ply birch plywood | 2 | 21 |

### Find an exact grade and inspect its conditions

```sh
cargo run -p xtask -- matdb-catalog --pack nilo-36-special-metals-2004
cargo run -p xtask -- matdb-catalog --pack diab-divinycell-h80-2021
cargo run -p xtask -- matdb-catalog --pack henkel-loctite-ea9466-seven-day-film-2025
cargo run -p xtask -- matdb-catalog --pack norplex-np509-melamine-glass-coupons-2022
cargo run -p xtask -- matdb-catalog --pack wisa-birch-plywood-18mm-upm007-2026
cargo run -p xtask -- matdb-catalog --search ensinger --property tensile-modulus
```

The distinctions in these packs affect which queries are admissible:

- Platinum annealed and work-hardened strength/elongation are separate states.
  Palladium annealed-wire strength requires its own specimen selector.
- NILO 36 expansion is an **interval mean** from a required reference temperature;
  it is not instantaneous alpha. Alloy transport curves have explicit interpolation
  policies and bounded support; MONEL's source-extrapolated conductivity rows
  are excluded.
- The LiF and ZnS thermal references use explicitly authored, narrow frozen
  calculation intervals with `engineering-reference-model=1` and
  `source-temperatures-matched=0`. They do not assert that the source measured
  all properties at one temperature. ZnS scalar elasticity separately requires
  `isotropic-engineering-model=1`.
- Ensinger tensile and compression test speeds are actual crosshead speeds.
  Compression stress samples remain `TabulatedOnly`; an unreported strain
  does not become a measured or interpolated stress.
- Divinycell nominal mechanical values and declared minima use separate keys
  and selectors. Density/mechanics are specified at 296.15 K, whereas conductivity
  is specified at 293.15 K. No foam heat capacity or complete crushing law is added.
- Henkel's 295.15 K cure temperature is distinct from its unspecified test
  temperature. EpoxAcast selects hardener 103 and the stated mass ratio; its
  liquid-mixture gravity does not become cured density. Resbond conductivity
  exposes its assumed Btu convention.
- NP509 test coupon thicknesses vary between bending, tension and compression.
  Its thermal expansion X/Y labels are preserved without inventing a mapping
  onto the mechanical lengthwise/crosswise directions.
- QUIKRETE's fresh-mixture density cannot be selected as hardened density.
  Cement strengths are discrete age observations. CEMENT ALL's first-hour result
  is explicitly measured after final set; later source-labelled ages retain their
  unspecified zero-time definition. Its pinned 2026 PDF method takes precedence
  over the differing live HTML method.
- MEDITE values marked as minima remain lower-limit properties; wet-cycle
  results require a different state. WISA's mean axial and bending moduli remain
  distinct, with explicit face-grain direction, thickness and ply count.
  Characteristic strengths and racking density do not become mean properties.

### Run the new plastic and foam through the axial reference example

These examples resolve material values from the compiled packs. Geometry,
monotonic displacement and the 0.1% reference strain ceiling are authored
scenario inputs. The ceiling is not a measured elasticity or failure limit.

Compile the selected packs:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/ensinger-tecapeek-gf30-natural/manifest.tsv --out /tmp/peek-gf30.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/diab-divinycell-h80-2021/manifest.tsv --out /tmp/divinycell-h80.fsmatpk
```

TECAPEEK GF30 tension, at the producer's tensile-modulus crosshead speed:

```sh
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/peek-gf30.fsmatpk tension density tensile-modulus 0.05 0.0001 0.001 \
  --axis source-grade-tecapeek-gf30-natural 1 --axis source-temperature-known 0 \
  --axis source-moisture-condition-known 0 --axis source-pressure-known 0 \
  --axis source-orientation-known 0 --axis test-speed 1.6666666666666667e-5 \
  --speed-axis test-speed --step 0.00005 3
```

Divinycell H80 nominal compression perpendicular to the sheet:

```sh
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/divinycell-h80.fsmatpk compression density compression-modulus 0.05 0.0001 0.001 \
  --axis source-grade-divinycell-h80-2021 1 --axis temperature 296.15 \
  --axis source-moisture-condition-known 0 --axis source-pressure-known 0 \
  --axis source-test-rate-known 0 --axis normal-to-sheet 1 \
  --axis source-is-nominal-value 1 --step -0.00005 1
```

Independent expectations for the authored linear references:

| Reference | Nominal force (N) | Mass (kg) | Model loading work (J) |
| --- | ---: | ---: | ---: |
| TECAPEEK GF30 tension | 630 | 0.00765 | 0.01575 |
| Divinycell H80 compression | -9 | 0.0004 | 0.000225 |

The foam's one-second ramp duration is an authored choice because its source
rate is unstated. Neither reference supplies unloading, relaxation, recoverable
energy, fracture, crushing, Poisson ratio or full three-dimensional elasticity.
The source unknown-orientation selector on the reinforced PEEK also remains
explicit; it is not a laminate-direction calibration.

### Focused native regressions

The existing `material_catalog_expansion` integration target includes 14 added
tests across the `engineering`, `alloys_minerals`, `polymers_foams` and
`adhesives` child modules. They compile all 26 packs, verify canonical decoding
and evidence receipts, check SI values and actual interpolation, and exercise
missing/wrong grade, age, temperature, rate, direction, thickness, cure and
source-assumption refusals. The iridium case also feeds the original bounded
slab-resistance consumer and refuses an unsupported hotter span.

## Reinforced engineering plastics and industrial PPSU

Four additional Ensinger packs contribute 36 claims for PPS, glass-filled PBT,
carbon-filled PEEK and industrial PPSU. They retain each producer sheet's printed
revision; a newer-looking download filename is not treated as a new revision.

| Pack | Material | Printed revision | Claims |
| --- | --- | --- | ---: |
| `ensinger-tecatron-natural` | PPS | 2023-07-19 AD | 9 |
| `ensinger-tecadur-pbt-gf30-natural` | PBT GF30 | 2023-07-19 AE | 9 |
| `ensinger-tecapeek-cf30-black` | PEEK CF30 | 2023-07-19 AI | 8 |
| `ensinger-tecason-p-black-radel-5500` | PPSU, Radel 5500 stock shapes | 2024-12-09 A0 | 10 |

The first three packs include density, heat capacity, conductivity, independently
named tensile/flexural properties and discrete compression-stress points. PPS
and PBT also have reported compression moduli. The current CF30 sheet supplies
no compression modulus; neither an older revision nor a stress/strain quotient
fills that gap. Reinforcement percentage basis and measurement orientation
remain unspecified. These scalar values do not supply an isotropic or
orthotropic stiffness tensor.

For these three ISO-based sheets, each mechanical claim carries its source
crosshead speed; flexure and compression also retain the printed 10 N setting
as `source-force-setting`. Its role is unspecified, so it is not interpreted as
preload or the calculated specimen force. Compression points are tabulated
only. Thermal properties have no numerical test temperature or temperature law.

The industrial PPSU pack uses the separate ASTM sheet. Its mechanical
properties are point-bound to 73 °F (295.92777777777775 K), with test rate unknown.
The tensile strength is explicitly a yield-strength property. Compression at
1% and 10% consists of two separate scalar point claims: the latter carries
the producer's public-source footnote, as do heat capacity and conductivity.
Intermediate strains and continuous compression-curve queries refuse.
The source's unspecified Btu variant is made explicit through the
`source-btu-conversion-it=1` choice for thermal conversions. Density and
thermal test temperatures remain unknown.

For a PPSU tensile reference using the same authored 50 mm / 100 mm² geometry:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/ensinger-tecason-p-black-radel-5500/manifest.tsv --out /tmp/ppsu-radel-5500.fsmatpk
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/ppsu-radel-5500.fsmatpk tension density tensile-modulus 0.05 0.0001 0.001 \
  --axis source-grade-tecason-p-black-radel-5500 1 \
  --axis source-moisture-condition-known 0 --axis source-pressure-known 0 \
  --axis source-measurement-orientation-known 0 \
  --axis source-density-temperature-known 0 --axis source-test-rate-known 0 \
  --axis temperature 295.92777777777775 --step 0.00005 1
```

The duration and strain ceiling are authored model choices. Each manifest
links its exact producer PDF, and each claim retains the relevant protocol
and secondary-source observations. No medical-grade equivalence or
three-dimensional constitutive law is inferred.

## Additional structural foams: ROHACELL and AIREX

Four producer-attributed packs add 50 claims for low-density PMI and structural
PVC foams. The ROHACELL 31/71 IG-F packs use the
[April 2022 Evonik sheet](https://products.evonik.com/assets/35/34/ROHACELL_IG_F_2022_April_EN_243534.pdf);
the AIREX C70.55/C70.130 packs use the
[August 2022 metric sheet](https://www.3accorematerials.com/uploads/pdf/TDS-AIREX-C70-E-08.2022-ex-Europe.pdf).

| Pack | Density, kg/m³ | Compression modulus, MPa | Tensile modulus, MPa |
| --- | ---: | --- | ---: |
| `evonik-rohacell-31-igf-2022` | 32 | 17 | 36 |
| `evonik-rohacell-71-igf-2022` | 75 | 73 | 92 |
| `airex-c70-55-2022` | 60 | 69 ISO/DIN; 58 ASTM | 45 |
| `airex-c70-130-2022` | 130 | 170 ISO/DIN; 145 ASTM | 115 |

These are typical/average producer values. AIREX minimum mechanical values use
separate `-minimum` property keys; the source's DNV-GL minimum definition is not
a confidence interval. Its compression is normal to the sheet, tension is
in-plane, and the shear plane is unknown. Density and compression retain the
block-center, at-least-40-mm specimen condition. ISO/DIN and ASTM compression
and shear modulus rows remain distinct selectable claims. Conductivity retains
the source's unspecified numerical room temperature and unknown direction;
no heat capacity or Poisson ratio is inferred.

For example, compile C70.55 and use its ISO/DIN compression modulus in an
authored 50-mm-long, 100-mm², 0.1%-strain axial reference:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/airex-c70-55-2022/manifest.tsv --out /tmp/airex-c70-55.fsmatpk
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/airex-c70-55.fsmatpk compression density compression-modulus 0.05 0.0001 0.001 \
  --axis source-grade-airex-c70-55-2022 1 \
  --axis source-temperature-known 0 --axis source-moisture-condition-known 0 \
  --axis source-pressure-known 0 --axis source-test-rate-known 0 \
  --axis block-center-specimen 1 --axis source-specimen-at-least-40mm 1 \
  --axis loading-in-sheet-plane 0 --axis loading-normal-to-sheet 1 \
  --axis compression-modulus-method-astm-c365 0 --step -0.00005 1
```

Selecting `compression-modulus-method-astm-c365=1` chooses the separate ASTM
row. Omitting that selector leaves the modulus query incomplete. The one-second
ramp and strain ceiling are scenario choices because the source supplies no
test rate or validated linear range. This uses the axial modulus reference;
it does not establish foam crushing, recovery, unloading, or three-dimensional
isotropy.

ROHACELL additionally requires the exact grade flag, unknown-temperature,
moisture and pressure flags, `source-test-rate-known=0`,
`source-loading-direction-known=0`, and `full-sheet-at-least-10mm=1` for
density. Its stated density spreads are manufacturing variation. The reported
ISO/DIN methods are retained; ASTM comparisons are not relabeled as performed
tests.

## Run axial reference calculations for plastics and foams

`material_axial_reference` resolves density and a named tensile or compression
modulus from a compiled source pack. It reports nominal axial force, displacement,
mass and model loading work for a declared geometry and monotonic loading path.

Compile PEEK and ROHACELL 51 IG-F:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/ensinger-tecapeek-natural/manifest.tsv --out /tmp/peek-axial.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/evonik-rohacell-51-igf-2022/manifest.tsv --out /tmp/rohacell-axial.fsmatpk
```

The following PEEK reference uses a 50 mm length, 100 mm² original area, 0.1%
strain ceiling, and a 50 µm displacement over three seconds. The ramp matches
the producer's declared tensile-modulus crosshead speed.

```sh
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/peek-axial.fsmatpk tension density tensile-modulus 0.05 0.0001 0.001 \
  --axis source-grade-tecapeek-natural 1 --axis source-temperature-known 0 \
  --axis source-moisture-condition-known 0 --axis source-pressure-known 0 \
  --axis test-speed 1.6666666666666667e-5 --speed-axis test-speed \
  --step 0.00005 3
```

This ROHACELL compression reference uses the same geometry and displacement
magnitude. The source test rate is unstated, so the one-second duration is an
explicit scenario choice.

```sh
cargo run -p fs-material --example material_axial_reference -- \
  /tmp/rohacell-axial.fsmatpk compression density compression-modulus 0.05 0.0001 0.001 \
  --axis source-grade-rohacell-51-igf-2022 1 --axis source-temperature-known 0 \
  --axis source-moisture-condition-known 0 --axis source-test-rate-known 0 \
  --axis full-sheet-at-least-10mm 1 --step -0.00005 1
```

Independent linear-reference expectations are:

| Case | Force (N) | Mass (kg) | Model loading work (J) |
| --- | ---: | ---: | ---: |
| PEEK tension | 420 | 0.00655 | 0.0105 |
| ROHACELL compression | -4.3 | 0.00026 | 0.0001075 |

These are authored zero-origin linear references using the reported moduli.
The strain ceiling is not a measured linearity or failure limit. The calculation
does not infer Poisson ratio, lateral contraction, unloading, relaxation,
recovery or strength. Tensile, flexural and compression moduli remain distinct,
and discrete compression-stress points are not consumed as moduli.

Supply repeated `--step` arguments for a monotonic history. At most 64 source
axes and 64 steps are accepted; a rejected step prevents the example from
emitting a successful path. The module documentation describes the complete
argument format and required source conditions.

## Run spatial transient heating from catalog curves

The `material_transient` example connects source Cp(T), rho(T) and k(T)
to the existing spatial enthalpy solver. It accepts one material or two layers
with a declared finite contact resistance. Every material retains its source
conditions and receipts, and mass uses density frozen at an explicit 300 K
reference state.

Compile the desired source packs:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/stainless-316-nasa-tp216435/manifest.tsv --out /tmp/nasa316.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/aluminum-2024-t3-first-heating-km2024/manifest.tsv --out /tmp/aa2024.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/peek-nasa-thermic-plate-linearized-model/manifest.tsv --out /tmp/nasapeek.fsmatpk
```

Then run steel, the newly added PEEK model, or an aluminum/steel stack:

```sh
cargo run -p fs-conduction --example material_transient -- nasa316 /tmp/nasa316.fsmatpk 0.01 8 120 120
cargo run -p fs-conduction --example material_transient -- nasapeek /tmp/nasapeek.fsmatpk 0.01 8 1200 600
cargo run -p fs-conduction --example material_transient -- aa2024 /tmp/aa2024.fsmatpk 0.01 8 120 120 nasa316 /tmp/nasa316.fsmatpk 0.01 0.0001
```

Arguments are profile, pack, thickness in metres, cells per layer, duration in
seconds and time steps, optionally followed by a second profile, pack,
thickness and area-specific contact resistance in m² K/W. The experiment
declares a 20×20 mm section, 300 K initial temperature, 360/300 K end reservoirs,
1000 W/(m² K) end-film coefficients and insulated sides. Each accepted step
prints spatial temperature profiles, reservoir/contact heat rates and a joule
balance. The final report compares total stored-energy change with accumulated
external heat.

The constructor explicitly selects isotropic scalar transport. The source
must support the whole 300–360 K path: the AA2024 profile retains first-heating
conditions, and PEEK explicitly selects the NASA constant-density model.
Unknown-temperature scalar packs do not acquire a thermal trajectory through
this adapter. Reference geometry and mass remain fixed; this calculation does
not supply thermal expansion, phase change or material aging.

## Search and inspect the material catalog

`matdb-catalog` searches the shipped source bundles using the actual offline
compiler. It supports case-insensitive text search and repeated required-property
filters; all supplied filters must match.

```bash
cargo run -p xtask -- matdb-catalog --search melamine
cargo run -p xtask -- matdb-catalog --search ensinger \
  --property thermal-conductivity --property specific-heat-capacity
cargo run -p xtask -- matdb-catalog --pack g10cr-nist-cryogenic
cargo run -p xtask -- matdb-catalog --pack smartply-strongdeck-osb4-30mm-reference
```

The list is deterministic TSV with source ID, internal pack ID, family, claim
count, canonical property names, manifest path, compilation status and diagnostic.
`--pack` selects the exact **source directory name**: internal pack IDs can be
shared by independent sources, including argon species and thermochemical packs.
Material/interface inspection includes coherent-SI values, curve knots, surface
samples and triangles, typed validity axes, interpolation, uncertainty, source
citations, observation methods and caveats. Species/model inspection is explicitly
metadata-only. No generated binary pack is published by inspection.

A source that the compiler refuses remains visible with its diagnostic in an
unfiltered listing; it cannot satisfy a required-property filter. Exact inspection
of that source fails. An unrelated refused source does not hide usable materials.
Property presence is a discovery aid: consumers still need a compatible state,
the required quantities and their own physical admission checks.

### Evaluate required properties at an explicit state

Repeat `--at AXIS=VALUE` to evaluate every requested property at the same
source-supported point. Coordinates use canonical SI, including kelvin for
absolute temperature; axis quantity types come from the source claims.
For example, inspect density and conductivity from the NASA 316 reference:

```sh
cargo run -p xtask -- matdb-catalog --pack stainless-316-nasa-tp216435 \
  --property density --property thermal-conductivity \
  --at temperature=350 --at source-generic-316-identity=1 \
  --at source-uns-designation-known=0 --at source-heat-treatment-known=0 \
  --at source-product-form-known=0 --at source-pressure-known=0
```

This evaluates 7938.6 kg/m³ and 14.3 W/(m K) using the source's declared linear
interpolation. Replace `--pack` with `--search` to retain only matching source
bundles whose requested properties evaluate successfully. Output includes SI
values, quantity types, selected claim hashes, interpolation decisions and
uncertainty. The results remain `Estimate` evaluations.

Every supplied coordinate must be declared by every selected claim. A density
entry with an unknown test temperature cannot silently satisfy a query at a
numerical temperature simply because the same pack has a temperature-dependent
conductivity curve. Missing conditions, unsupported temperatures or orientations,
discrete-knot gaps, ambiguous observations, conflicting quantity types and
extra axes refuse; the command supplies no implicit condition defaults.
Inspect the pack without `--at` to see its exact required axes and caveats.
Successful samples do not establish a matched specimen or replace a physical
consumer's additional admission checks.

## Everyday materials expansion: boards, polymers, alloys and composites

This increment adds **38 source-specific packs and 195 property claims**.
Together with the existing wood, plywood, metal, insulation and resin sources,
these additions provide substantially broader everyday construction and fabrication
coverage. Counts refer to source conditions and claims, rather than complete
constitutive models or independent experimental measurements.

| Source group | Packs | Claims | Named additions |
|---|---:|---:|---|
| Elements and minerals | 8 | 31 | Chromium, niobium, sapphire, barium fluoride, zinc selenide, magnesium fluoride, crystalline quartz and calcite |
| Alloys | 3 | 12 | Hydro 6063-T6 and 6082-T6 extrusions; Special Metals INCONEL 625 with thermal/electrical temperature curves |
| Fibers and composites | 4 | 20 | NIST G10-CR fiberglass/epoxy; Norplex NP510A FR-4 coupons; Toray T700S fiber and separately identified epoxy-composite coupons |
| Plastics, films, resins, glues and foam | 12 | 71 | Ensinger PA66/POM-H/PEI; PLEXIGLAS GS; Mylar 92A and 1000A; 3M DP420 epoxy; DOWSIL 732 silicone; Basotect G+ melamine foam; SABIC LDPE; Titebond Original; Smooth-Cast 300 polyurethane |
| Boards and cementitious products | 11 | 61 | Finsa MDF/HDF/particleboard; Kronospan melamine-faced P2; Fermacell gypsum-fiber and cement boards; Sika repair mortar; two Bauroc AAC grades; SMARTPLY OSB/4; Metsa Kerto S LVL |

Every new source directory retains its primary citation, attributed numerical
facts, original condition notes, source units or explicit conversion notes and
unstated uncertainty where no uncertainty was supplied. The following distinctions
are enforced by property names and query coordinates:

- Mylar machine/transverse directions and film gauges remain separate.
- Crystal/G10-CR conduction and OSB/LVL stiffness retain their published directions.
- Tensile, flexural and compression moduli retain their measurement meanings.
- Product specifications, tolerances and fifth-percentile values remain explicit;
  density ranges and strength limits do not become invented point measurements.
- Fresh mortar, uncured wood glue, cured silicone and seven-day casting resin
  retain different states. Discrete cure ages do not imply a continuous curing law.
- Bare carbon fiber, FR-4 coupons and normalized epoxy composites retain their
  own constructions, thicknesses, cure conditions and available property sets.
- INCONEL transport retains its anneal condition; its nominal reference density
  does not supply an elevated-temperature density curve.

The focused test targets are:

```bash
cargo test -p xtask --test elements_alloys_composites_expansion
cargo test -p xtask --test building_products_expansion
cargo test -p xtask --test polymers_resins_expansion
cargo test -p xtask --test material_catalog_search
```

The 20 material source/query integration cases and three catalog CLI cases passed
in a narrow native build using the pinned `nightly-2026-07-06` compiler and the
actual repository compiler/runtime modules. They cover compilation, canonical
serialization, hash-verified decoding, SI queries, receipt replay and unsupported
state refusals. Full workspace/DSR validation and solver-level or experimental
qualification were not performed for this increment.

## Continuous NASA PEEK thermal-model curves

`peek-nasa-thermic-plate-linearized-model/` makes the published NASA THERMIC
plate correlations usable across **300–525 K**. It is a separately selected
model pack; the original `peek-nasa-thermic-plate/` exact-point observations
remain unchanged. The NASA plate's grade, supplier and processing history are
unspecified, so it is distinct from the Ensinger PEEK stock-shape pack.

Two curves retain 46 exact-decimal evaluations apiece, spaced every 5 K, of
the conductivity and Cp cubics in NASA/TM-20210014330. Linear interpolation
introduces at most 0.000010446578125 W/(m K) in k and 0.06480046875 J/(kg K)
in Cp relative to those cubics, in real arithmetic. These bounds exclude
floating-point roundoff and any subsequent enthalpy-chart interpolation.
The integrated Cp approximation over the whole interval is
319647.6954140625 J/kg; the polynomial integral is 319641.29033203125 J/kg.
This numerical difference is not experimental uncertainty, which stays unstated.

The density claim retains the paper's **1264 kg/m³ constant-density model
input** over the same interval. It does not assert measured density throughout
that interval. Queries must explicitly supply `source-constant-density-model=1`
and `source-pressure-atmospheric=1`; density also requires
`source-test-temperature-known=0`. Temperature uses the absolute-temperature
quantity kind. Outside-range and missing-model-condition requests refuse.

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/peek-nasa-thermic-plate-linearized-model/manifest.tsv \
  --out /tmp/peek-nasa-thermic-plate-linearized-model.fsmatpk
cargo test -p xtask --test material_catalog_expansion g0_g3_peek_model_curves
frankensim discover examples/material-discovery/peek-nasa-thermal-model.json \
  /tmp/peek-nasa-thermic-plate-linearized-model.fsmatpk
```

The source report and original fit limitations are linked in the manifest and
observations. The two numerical curves are three model claims together with
density, not 92 independent property measurements.

## Engineering polymers with separate loading modes

Four producer-grade stock-shape packs add **32 claims**. Each retains density,
heat capacity, conductivity, tensile modulus and strength, flexural modulus,
compression modulus, and a discrete compression-stress table.

| Directory | Polymer | Density (kg/m³) | Cp (J/kg/K) | k (W/m/K) | Tensile / flexural / compression modulus (MPa) |
| --- | --- | ---: | ---: | ---: | --- |
| `ensinger-tecapeek-natural` | PEEK | 1310 | 1100 | 0.27 | 4200 / 4200 / 3400 |
| `ensinger-tecapei-natural` | PEI | 1280 | 1200 | 0.21 | 3200 / 3300 / 2800 |
| `ensinger-tecason-s-natural` | PSU | 1240 | 1200 | 0.21 | 2700 / 2600 / 2300 |
| `ensinger-tecamid-66-natural` | PA66 | 1150 | 1500 | 0.36 | 3500 / 3100 / 2700 |

These are producer comparison data with unstated statistical uncertainty.
Numerical test temperature and moisture conditioning remain unknown; PA66
retains the producer's immediately-after-machining condition. Each modulus
keeps its test method, specimen and crosshead-speed context. A tensile modulus
is not a flexural or compression modulus, and none supplies a Poisson ratio.
Compression stresses at 1%, 2% and 5% strain are tabulated-only claims: they
do not establish an interpolated stress-strain or unloading law. Temperature
service limits do not establish a thermal-property trajectory.

Compile any pack with the existing compiler:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/ensinger-tecapeek-natural/manifest.tsv \
  --out /tmp/ensinger-tecapeek-natural.fsmatpk
```

## Engineering alloys with independent elastic parameters

Three added copper-alloy packs contain density, Young's modulus, shear modulus,
specific heat and thermal conductivity. Their original source data are sufficient
for the explicit `YoungShear` resolver in `fs-material::state_point`; no
Poisson ratio is fabricated in the database.

| Pack | Primary source | Mechanical use |
|---|---|---|
| `phosphor-bronze-c52100-cda-reference` | [Copper Development Association C52100](https://alloys.copper.org/alloy/C52100) | Nominal copper-tin-phosphorus alloy elasticity from the published E/G pair |
| `beryllium-copper-c17200-ngk-reference` | [NGK physical properties](https://www.ngkmetals.com/technical/properties/physical-properties/) | Supplier-specific copper-beryllium E/G reference |
| `copper-nickel-beryllium-c17510-ngk-reference` | [NGK physical properties](https://www.ngkmetals.com/technical/properties/physical-properties/) | Supplier-specific copper-nickel-beryllium E/G reference |

The canonical source property keys are `density`, `young-modulus`,
`shear-modulus`, `specific-heat-capacity`, and `thermal-conductivity`.
Density and the two thermal properties retain the printed 68 °F / 293.15 K state.
The moduli have no reported test temperature and require
`source-test-temperature-known=0`; they have no invented temperature interval.
All five properties retain unknown temper and pressure through
`source-temper-known=0` and `source-pressure-known=0`.

A combined density/E/G query uses those three explicit unknown-condition flags
and `temperature=293.15`. Using the elastic values in this nominal isotropic
model is a caller assumption; it does not establish co-measured compatibility,
a particular heat treatment, strength allowables or a temperature-independent
constitutive law. Source observations retain original imperial values, exact SI
conversion factors and the explicit Btu_IT convention. Conversion digits do not
increase the source precision.

Compile any of these manifests through the same `matdb-pack` command used below.
All three were admitted by the actual compiler, including its repeat-compilation
and hash-pinned decode checks.

## Run a sourced mechanical coupon

`fs-solid` now includes `material_coupon`: a source-pack-driven static tensile
or compressive coupon with explicit geometry, load, supports and constitutive
assumptions. Material values come from the compiled pack. Isotropic modes admit
E/nu, E/G or K/G; the full FPL reference-ratio mode retains orthotropy and an
explicit L/R/T orientation. No source temperature is invented for wood.

```bash
coupon_dir="$(mktemp -d /tmp/frankensim-coupon.XXXXXX)"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/phosphor-bronze-c52100-cda-reference/manifest.tsv \
  --out "$coupon_dir/bronze.fsmatpk"
cargo run -p fs-solid --example material_coupon -- \
  --template phosphor-bronze-c52100 > "$coupon_dir/bronze.tsv"
cargo run -p fs-solid --example material_coupon -- \
  "$coupon_dir/bronze.fsmatpk" "$coupon_dir/bronze.tsv"
```

Templates also exist for `aluminum-5005` and `loblolly-pine`; compile their
`aluminum-5005-austral-wright-reference` and `pine-loblolly-fpl-gtr282` manifests
respectively. The TSV declares exact property-role mappings and every source
coordinate. Values are coherent SI, so 12% dry-mass moisture is `0.12`.
Axis quantity types are read from the actual selected source declarations.
Conflicting claims, incomplete source coordinates, missing 3-D wood data,
unsupported response axes and excessive recovered strain refuse.

The default coupon is 100 x 20 x 20 mm with four axial cells and 1 MPa traction
(400 N). Left-face normal rollers and minimal transverse pins allow Poisson
contraction. Native runs of the actual source compiler and physical consumer gave:

| Source and orientation | Mass, kg | Mean axial displacement, micrometres |
|---|---:|---:|
| Aluminum 5005 | 0.108 | 1.466275660 |
| C52100 phosphor bronze | 0.352088388 | 0.906485861 |
| Loblolly pine, longitudinal | 0.022848 | 7.390983001 |
| Loblolly pine, radial | 0.022848 | 65.406929210 |

For the radial case, replace the wood template's orientation row with
`orientation<TAB>0<TAB>1<TAB>0<TAB>-1<TAB>0<TAB>0<TAB>0<TAB>0<TAB>1` using actual
tab separators. The radial/longitudinal displacement ratio agrees with the
retained source's `1 / er_over_el` (8.849557522), with unchanged mass and load.
These are nominal fixed-state, small-strain model results; general rotated
anisotropy remains a coarse finite-element calculation, and no strength
allowable, loading-history, fracture or mesh-error claim is added.

Output includes applied force, support reactions, displacement, mass, stored
elastic energy, strain and force residuals, alongside the selected source
coordinates, claims and physical identities. `load-displacement-product-j` is
the final F dot u pairing, twice stored elastic energy in this zero-initial-strain
linear case; it is not work integrated over a loading history.

```bash
cargo test -p fs-material state_point::isotropic_parameters
cargo test -p fs-solid --lib nodal_force_tests
cargo test -p fs-solid --example material_coupon
```

## Run a sourced layered heat calculation

The steady planar wall consumer is in
`fs_conduction::layered`. Each layer retains its own `ConductivityTable`,
selected source receipts, thickness and optional material-card identity.
Film and contact resistances are explicit inputs. The result includes heat
flux, total heat rate, U-value, and the temperature on both faces of every
layer. Density and heat capacity are not required for a steady calculation.

This example uses three existing NIST CFAST nominal model-input records:

```sh
mkdir -p /tmp/frankensim-wall
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/gypsum-board-5-8-cfast/manifest.tsv --out /tmp/frankensim-wall/gypsum.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/glass-fiber-insulation-cfast/manifest.tsv --out /tmp/frankensim-wall/glass-fiber.fsmatpk
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/brick-common-cfast-sp1041/manifest.tsv --out /tmp/frankensim-wall/brick.fsmatpk
cargo run -p fs-conduction --example material_wall -- /tmp/frankensim-wall/gypsum.fsmatpk /tmp/frankensim-wall/glass-fiber.fsmatpk /tmp/frankensim-wall/brick.fsmatpk
```

The scenario declares 16 mm gypsum, 88 mm glass-fiber insulation and 76 mm brick,
293.15/273.15 K reservoirs, 10 m² area, film resistances of 0.13/0.04 m² K/W
and contact resistances of 0.01/0.002 m² K/W. The analytical reference is
R'' = 2.5875555556 m² K/W, U = 0.3864651323 W/(m² K), q'' = 7.7293026451 W/m²
and heat rate = 77.2930264514 W.

These source records identify CFAST inputs rather than measured specimens.
Unknown temperature, moisture, manufacturing and measurement conditions stay
explicitly unknown in the query receipts. The example's temperature span is a
declared nominal-model operating interval; it does not establish experimental
coverage. Film/contact values are scenario assumptions. This planar steady
calculation includes no moisture transport, lateral thermal bridging, thermal
radiation, heat sources or transient storage.

The focused catalog integration regression also compiles these sources, stores
and reopens their immutable cards, discovers the exact conductivity keys, and
passes those selections to the physical wall owner:

```sh
cargo test -p fs-conduction --test layered
cargo test -p xtask --test catalog_layered_wall -- --nocapture
```

## October 2026 material expansion

The October 28 expansion adds **34 source bundles and 272 property claims**:
14 elemental/mineral/glass/alloy references, eight wood species, four engineering
plastics, six insulation/structural-foam products, and two cured epoxy systems.
These are source-condition bundles, not 34 universally complete material laws.
Each manifest identifies its original source and redistribution terms; each
property retains units, observations, uncertainty and required query conditions.

| Family | New source bundles | Useful coverage |
| --- | --- | --- |
| Elements | `diamond-crystran-thermal-reference`, `germanium-crystran-thermal-reference`, `molybdenum-attl-reference`, `tantalum-efinea-reference`, `magnesium-espi-reference` | Density, heat capacity and thermal conductivity; selected electrical resistivity and modulus values |
| Minerals and glass | `fluorite-caf2-crystran-reference`, `halite-nacl-crystran-reference`, `periclase-mgo-crystran-reference`, `fused-silica-crystran-reference` | Complete thermal triplets; amorphous fused silica also supplies Young modulus and Poisson ratio |
| Alloys | `aluminum-5005-austral-wright-reference`, `aluminum-5052-austral-wright-reference`, `bronze-c93200-cda-reference`, `bronze-c95400-cda-reference`, `cupronickel-c70600-cda-reference` | Thermal triplets and stiffness, with shear modulus, resistivity and explicitly derived Poisson ratios where supported |
| Woods | `pine-{loblolly,longleaf,slash,ponderosa}-fpl-gtr282`, `oak-{northern-red,white}-fpl-gtr282`, `cedar-northern-white-fpl-gtr282`, `fir-subalpine-fpl-gtr282` | Density, bending/longitudinal stiffness and strength at 12% moisture; six conifers also retain directional stiffness and Poisson ratios |
| Engineering plastics | `ensinger-tecaform-ah-natural`, `ensinger-tecamid-6-natural`, `ensinger-tecapet-white`, `ensinger-tecaflon-pvdf-natural` | POM-C, PA6, PET and PVDF density/Cp/conductivity, distinct test moduli and discrete compression observations |
| Insulation | `eps-nist-srm1453-2013`, `rockwool-prorox-sl960-en-2020`, `austrotherm-xps-top30-sf-2026`, `austrotherm-xps-top70-sf-2026` | EPS conductivity as a density/mean-temperature surface; stone-wool conductivity over 50–660 °C; product/thickness-specific XPS declared conductivity and mechanical data |
| Structural foam | `evonik-rohacell-51-igf-2022`, `evonik-rohacell-110-igf-2022` | Separate closed-cell PMI grades with density, compression/tensile/shear response |
| Resins | `west-system-105-205-neat-cured`, `west-system-gflex-650-neat-cured` | Separate neat epoxies after 14-day room-temperature cure: tensile/flexural/compressive properties and specific gravity |

### Compile and discover the additions

Run from the repository root; choose a previously absent output filename:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/ensinger-tecaform-ah-natural/manifest.tsv \
  --out /tmp/tecaform-pom-c.fsmatpk

cargo run -p fs-cli -- --json discover \
  examples/material-discovery/tecaform-pom-c-thermal.json \
  /tmp/tecaform-pom-c.fsmatpk
```

Discovery consumes the generated pack paths explicitly. The following requests
select useful property sets and preserve the required source context:

| Discovery request under `examples/material-discovery/` | Applicable new packs |
| --- | --- |
| `engineering-reference-thermal-293k.json` | Five alloys, molybdenum, tantalum, germanium and fluorite at the declared 293.15 K reference |
| `minerals-reference-thermal-273k.json` | Diamond, halite and periclase frozen references over 273–274 K |
| `glass-magnesium-reference-thermal-298k.json` | Fused silica and magnesium frozen references over 298.15–299.15 K |
| `clear-wood-12pct.json` | All eight new wood species; longitudinal mechanics, not an isotropic elastic request |
| `tecaform-pom-c-thermal.json`, `tecamid-pa6-thermal.json`, `tecapet-pet-thermal.json`, `tecaflon-pvdf-thermal.json` | Exact named plastic grades and their stated unknown-condition flags |
| `eps-srm1453-conductivity.json` | SRM1453 lot, 281–313 K mean temperature, 37–46 kg/m3, measured 20 K temperature difference |
| `rockwool-sl960-conductivity.json` | European 2020 SL960 conductivity curve with its actual source range |
| `rohacell-51-mechanical.json` | Full-sheet 51 IG-F density and compression modulus; 110 IG-F remains a distinct grade |
| `west-105-205-mechanical.json` | Neat 105/205 epoxy at 14-day cure; G/flex remains a distinct formulation |
| `austrotherm-xps-top30-conductivity.json` | TOP30 SF at a supported declared board thickness |

### Conditions that remain explicit

The 14 elemental/mineral/glass/alloy packs are **frozen engineering references**.
Their one-kelvin intervals are author-declared calculation ranges, not measured
temperature curves. Queries require `engineering-reference-model=1` and
`source-temperatures-matched=0`; per-property observations retain actual source
temperatures. Derived Poisson ratios are labeled as derivations. Cubic-crystal
thermal scalars do not establish isotropic elasticity.

USDA values describe clear-specimen species averages at 12% moisture, not
graded-lumber allowables. The six conifers retain the input names used by
`fs_material::state_point::resolve_orthotropic_elastic_state_point_from_reference_ratios`.
The two oaks deliberately lack unsupported directional ratios. Density derives
from the table's dry-mass specific gravity and stated moisture convention.

Plastic compression points remain `TabulatedOnly`: for example, the POM data
support the measured 1%, 2% and 5% points, not an invented 3% stress or unloading
law. Test moduli remain distinct from isotropic Young modulus. PA6 retains its
immediately-after-machining condition. Unknown source temperatures stay unknown.

EPS uses the published affine surface, not a density-independent constant.
Its 1.5% uncertainty retains NIST's stated posterior-predictive interpretation.
ROCKWOOL and XPS data do not supply missing heat capacity or porous-acoustic laws.
XPS declared conductivity is not silently relabeled as a national design value;
unsupported board thicknesses and thermally bonded variants remain separate.
ROHACELL manufacturing spreads are not statistical confidence intervals.
Epoxy specific gravity remains dimensionless because no reference-water
temperature is supplied; no bulk density or joint adhesion law is inferred.

Focused regression command:

```sh
cargo test -p xtask --test material_catalog_expansion
```

The regression exercises the real compiler, canonical pack decoding, source
literals and SI conversions, typed queries and receipt replay, interpolation,
and refusals for unsupported material conditions. It verifies the database and
evaluator path; it does not qualify every source as a measured specimen or run
a full thermo-mechanical simulation for every new material.

## Compiler profiles

Bulk scalar/curve sources use `material-tsv-v1`. Ordered interface-system
sources use `interface-tsv-v1`: they retain that same provenance, uncertainty,
validity, frame, and joint-statistics grammar, and additionally require exactly
one record for each ordered surface and one complete system context:

```text
surface_a<TAB>chemistry<TAB>phase<TAB>process<TAB>revision<TAB>texture_frame
surface_b<TAB>chemistry<TAB>phase<TAB>process<TAB>revision<TAB>texture_frame
context<TAB>medium<TAB>third_body-or--<TAB>environment<TAB>named_history
```

The compiler emits `FSINTPK` v1, preserving surface order and named history as
identity-bearing state. This profile does not carry constitutive model cards,
infer symmetry, or turn a fixture claim into authority for either bulk material.
The committed directories below remain bulk/species sources until separately
identified, redistributable interface evidence is added.

## Common material acquisition, September 2026

`shell-heat-transfer-oil-s2-2011` retains four common-temperature rows from
Shell's May 2011 technical data sheet for a named refined mineral heat-transfer
oil: density, heat capacity, conductivity and kinematic viscosity at 0, 40,
100 and 200 °C. The kinematic adapter in `fs-material::liquid` resolves all four
at one point and derives dynamic viscosity using that point's density. This is
not a hydraulic viscosity-grade specification. Pressure, precise additive
composition and measurement uncertainty are unstated; aged-oil use is excluded.
The source's printed Prandtl row disagrees with its other property tables by
more than rounding, and that unresolved discrepancy is retained in the source
observations. Derived Prandtl values are engineering calculations, not another
measurement. Sparse linear interpolation has no accuracy bound. The pack does
not support friction, wear, oxidation or boiling, nor extrapolation beyond its
retained rows. The source document itself is not redistributed; this pack contains
only selected attributed numerical facts and original explanatory notes.
Discovery request: `examples/material-discovery/shell-heat-transfer-oil-s2.json`.

Scalar/curve source packs use manifest v3 with explicit property-name mappings and typed
temperature axes. Names such as `specific-heat-capacity` and
`thermal-conductivity` therefore match current discovery requests. A directory
is a source bundle with its own conditions; it is not automatically a complete
material card or a qualified simulation.

`copper-annealed-iacs-nbs-hb100` and `aluminum-ec-h19-nbs-hb109` supply
temperature-dependent **bulk electrical resistivity** in ohm metres, for
standard annealed 100% IACS copper and EC-H19 aluminum respectively. Copper
retains the handbook's rounded 20 °C reference and published volume-resistivity
slope over a bounded 10–30 °C range. Aluminum retains five literal tabulated
values over 0–30 °C with linear interpolation and the source's approximation
neglecting volume expansion. Neither pack substitutes a wire-resistance
temperature coefficient for a volume-resistivity coefficient, or infers a
modern UNS/alloy specification from a historical conductivity designation.
Both retain unknown pressure and unstated statistical uncertainty.

The `fs-material::conductor` adapter resolves the selected claim and its receipt
to uniform DC resistance using caller-supplied **actual geometry at the queried
state**, `R = rho_e L/A`. It supplies resistance to the existing circuit owner;
it does not infer thermal expansion, contact resistance, stranding, AC skin
effects, ampacity or a self-heating temperature trajectory. Discovery requests
are in `examples/material-discovery/{copper,aluminum}-dc-conductor.json`.

`stainless-316-20c-engineering-reference` supplies a complete six-property
isotropic thermoelastic input set at exactly 293.15 K: density, Young modulus,
Poisson ratio, specific heat, thermal conductivity and instantaneous expansion.
It explicitly joins Metalcor's 316 density at 20 °C and Austral Wright's typical
annealed 316 Poisson ratio at 20 °C with direct NIST 316 equation evaluations
at the same temperature. This is a nominal bulk austenitic-grade engineering
approximation; only the Poisson source specifies annealing. It does not establish
a common specimen, product form, heat treatment, pressure or joint uncertainty.
The manifest distributes attributed numerical facts and derived values; it
does not license or redistribute the supplier documents. Existing NIST curve
packs remain separate and retain their original coverage.

The NIST expansion input is `epsilon'(T)/(1+epsilon(T))`, where epsilon is
relative to the source's **293 K** length. It is an instantaneous coefficient,
not interval-mean expansion; the 293 K fit residual is retained. The complete
profile is valid only at 293.15 K with `source-pressure-known=0` (unknown),
not over the full individual NIST curve ranges. Together with the existing
2024-T3 NASA profile at 300 K, it supports two distinct reference-state plate
calculations through the existing material resolver. The temperatures differ,
so comparing those calculations is not a controlled same-temperature material
substitution, and neither point set authorizes a finite heating trajectory.

`stainless-316-nasa-tp216435` retains three literal rows at 273, 300, and 400 K
from NASA/TP-2010-216435 Table 23 for the report's generic “316 stainless
steel”: density, specific heat, and thermal conductivity. Linear interpolation
between the table knots is explicitly an engineering approximation. The NASA
record is public and permits public use, but the report does not identify UNS
designation, chemistry, product form, heat treatment, specimen, method,
pressure, or uncertainty. This source therefore cannot be merged automatically
with UNS S31600/316L records or the 316 20 °C engineering reference. It
deliberately excludes the table's separate CTE column, and supplies no Young
modulus, Poisson ratio, enthalpy, or complete heating profile. Discovery is
available in [`warm-316-thermal.json`](../../../examples/material-discovery/warm-316-thermal.json).

`dowtherm-sr1-eg50-40-60c` and `dowfrost-pg50-40-60c` retain density,
specific heat, conductivity and dynamic viscosity from the Dow-authored
2001 SR-1 and 2007 DOWFROST engineering guides. Each curve has five literal
5 °C source rows over 40–60 °C, with explicitly approximate linear
interpolation. The concentration is exactly 50% **glycol by volume**, not
50% product or 50% glycol by mass. These are distinct inhibited products;
the later technical data sheets' differing values are not pooled into them.
Only selected attributed numerical facts and original notes are distributed.

The guides do not state the volume-reference temperature or numerical
property-test pressure. Both remain explicitly unknown in the query flags;
the 20 °C water-density reference used for specific gravity does not establish
a concentration-reference temperature. The selected rows are below the
tables' atmospheric-boiling markings. The separate phase tables lack an
exact 50.0-volume-percent row, so freezing/boiling points and mass-to-volume
conversion are not invented for these profiles. Other concentrations, phase
changes and operation outside 40–60 °C remain unsupported.
The source-resolved tuple is exposed by `fs_material::liquid::resolve_liquid_state`;
the existing heat and flow owners consume its coefficients. Discovery requests
are [`dowtherm-sr1-eg50.json`](../../../examples/material-discovery/dowtherm-sr1-eg50.json)
and [`dowfrost-pg50.json`](../../../examples/material-discovery/dowfrost-pg50.json).

`air-dry-ussa1976` supplies five parameters for the existing calorically perfect
gas model: molar mass, heat-capacity ratio, reference viscosity, its absolute
reference temperature, and the Sutherland temperature interval. The static
reference viscosity is evaluated from the source equation at 273.15 K;
`fs_material::gas::resolve_sutherland_gas_state` recovers its coefficient and
uses the existing `GasState` equations. The caller explicitly selects
`ConductivityModel::Ussa1976AirFit`; the conductivity equation and historical
gas constant remain in that existing code, rather than a duplicate EOS.
The model then derives density, Cp/Cv, sound speed, viscosity, conductivity,
Prandtl number and impedance at the requested temperature and pressure.
Every selected parameter keeps its source receipt.

The declared application range is 273.15–313.15 K, 80–110 kPa, exactly zero
relative humidity, and the USSA-1976 reference composition. It is a bounded
engineering model, not measured coverage of every ambient condition. The
rounded composition list is retained without silently renormalizing it or
identifying it as today's universal atmosphere. The source's p. 3 gas constant
and p. 19 equations 51/53 take precedence over the errors in Table 2 and the
printed conductivity column; these errors are also identified in the appended
errata. Source use is permitted by the NTRS record. Discovery of the complete
parameter set uses [`dry-air.json`](../../../examples/material-discovery/dry-air.json).
Humidity transport, real-gas effects and experimental accuracy remain outside
this profile. It does not claim nine independent measured output properties.

`water-vapor-sutherland-ambient` supplies the five corresponding parameters
for a dilute water-vapor component. Its molar mass comes from NASA
TP-2002-211556, while its constant heat capacity uses NIST-JANAF H-064 at
298.15 K. The distinct NASA heat-capacity value is not silently substituted.
The COMSOL/White steam viscosity reference is retained literally at 350 K;
using its Sutherland relation at 273.15–313.15 K is an explicitly unbounded
engineering extrapolation. Eucken conductivity is likewise a model estimate.
The 80–110 kPa axis denotes total mixture pressure, not stable pure steam.
Discovery uses [`water-vapor-component.json`](../../../examples/material-discovery/water-vapor-component.json).

`fs_material::moist_air::resolve_moist_air_state` keeps both component cards
and ten receipts, takes relative humidity separately, and uses the existing
gas/acoustic path. It mixes ideal molar heat capacities; viscosity and
conductivity follow the [IDAES Wilke](https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/transport_properties/viscosity_wilke.html)
and [Wassiljewa–Mason–Saxena](https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/transport_properties/thermal_conductivity_wms.html)
formulations. The existing Buck-1996 liquid-water saturation equation is
documented by [CIRES's formula comparison](https://cires1.colorado.edu/~voemel/vp.html),
which attributes it to the Buck Research CR-1A manual, revised July 1996;
that manufacturer manual was not retrieved in this acquisition.
Pure-water saturation pressure omits the moist-air enhancement factor.
Supersaturation and unsupported component conditions refuse; condensation,
real-gas behavior and experimental accuracy are not qualified by this model.

`silicon-cubic-25c-nasa-rp1057` uses manifest v6 to carry all 36 explicitly
addressed engineering stiffness entries in crystal axes `[100],[010],[001]`,
order `xx,yy,zz,yz,zx,xy`, at exactly 298.15 K. They are a cubic-symmetry
expansion of **three** handbook constants, not 36 independent measurements.
The Mbar interpretation of the unlabeled stiffness table follows the chapter's
pressure convention and is an explicit unit inference. Density comes from a
separate pure-crystal reference at the same temperature. This combination is
an engineering reference, with unknown pressure and stiffness specimen/doping;
it does not qualify a particular wafer or supply thermal/electrical data.
Sources and the NTRS `GOV_PUBLIC_USE_PERMITTED` decisions are in its manifest.
Frame and tensor IDs are BLAKE3 hashes of these exact UTF-8 strings, without a
trailing newline:

```text
frankensim silicon crystal axes 1=[100],2=[010],3=[001]; engineering xx,yy,zz,yz,zx,xy; NASA-CR-275 App A
frankensim silicon cubic stiffness NASA RP-1057 p101 at25C; C11=165773MPa C12=63924MPa C44=79619MPa; explicit cubic expansion
```

The focused source-to-consumer suite is executable with
`scripts/e2e/material_sources.sh --run` (required-remote RCH), or with an
explicit `MATDB_PACK_TEST_BIN` from an existing native test build. `--list`
lists exact cases and `--check` checks prerequisites without building.
Each run retains command, stdout/stderr, physics diagnostics and exit status
in the unique receipt directory printed at startup. This suite checks source
transport and bounded numerical responses, not experimental qualification.

`water-liquid-iapws-sr6-08` supplies five liquid-water curves (density, isobaric
heat capacity, specific enthalpy, conductivity and dynamic viscosity), each
with 15 knots from 283.15 to 353.15 K at exactly 100000 Pa. The source is
[IAPWS SR6-08(2011)](https://iapws.org/documents/release/LiquidWater.download),
which permits republication with attribution. The enthalpy convention is
inherited from its Gibbs-energy equation, with an explicit query flag; liquid
phase and pressure are also required. The reusable discovery request is
[`liquid-water.json`](../../../examples/material-discovery/liquid-water.json).
Linear interpolation adds sampled relative discrepancies of approximately
0.00393%, 0.00788%, 0.00796%, 0.01888% and 0.38335%, respectively, against
the source equations. These are numerical sampling results, separate from
source uncertainty and not certified bounds. Independently interpolated
enthalpy and heat capacity differ by at most 1.100 J/kg per 5 K interval
under the integral consistency check, and 3.519 J/(kg K) under the derivative
check. This bounded dataset does not supply boiling, steam, ice or a fluid solve.

| Source family | Material records | What the records support |
| --- | --- | --- |
| NASA-CR-71699 | Annealed pure aluminum and copper; pure iron, nickel and titanium; 7075-T6 aluminum; X-750; 304A and 347 stainless; titanium A-110AT | Selected temperature/conductivity tables with explicit linear engineering interpolation, excluding source-marked extrapolation and unrepresented phase discontinuities |
| NIST cryogenic material pages | 304 and 316 stainless, C26000 brass, 1100 aluminum | Stainless conductivity, heat-capacity, relative-expansion and model-derived instantaneous-expansion interpolation over 77–300 K; modulus over 77–293 K (304) or 77–294 K (316). Other materials retain exact temperatures. Complete coefficients, derivations and source ranges remain in the observations |
| NASA/TP-3287/REV1 | Aluminum, copper, lead, tin and zinc crystalline, liquid and fusion records | Phase-specific heat capacity and enthalpy, plus fusion energy from paired rows at the same transition temperature and reference state |
| NASA/TP-2010-216435, Table 23 | Generic “316 stainless steel” | Literal 273/300/400 K density, conductivity, and specific-heat rows with explicit linear engineering interpolation; no UNS/form/heat/pressure identity or complete thermoelastic profile |
| NASA-CR-144016 | Elemental liquid lead | Compiled temperature-dependent density, conductivity, viscosity and surface-tension rows, retaining source pressure/purity/surface gaps and the explicit conductivity calorie conversion |
| NBS RP500, Table 2 | Cast-lead sample 1144, 99.9% Pb, sand mold | Measured density at 25 °C and four mean linear-expansion coefficients from 20 °C to 60/100/200/300 °C on the first heating run; individual sample retained, without pooling heating/cooling histories |
| NBS RP668, Table 5 | Lead conductivity reference standard L.S., cast in a bottom-feed cast-iron mold and machined | Published calibrated conductivity estimate over 0–300 °C, explicitly acknowledging the assumed absolute calibration; no numerical purity or pressure inferred |
| NIST CFAST thermal library | Normal/lightweight concrete, cement mortar, clay/common brick, calcium-silicate board, cellulose/glass-fiber/urethane insulation, ordinary/Type X gypsum board | Named textbook/model inputs for conductivity, specific heat, density and layer/observer parameters; not specimen measurements or grade allowables |
| NISTIR 6475 | Normal-strength concrete Mixture IV, cylinder ST-IV-25-1 | Measured room-temperature strength, dynamic and static modulus, plus explicitly derived apparent cylinder density; the mix and conditioning remain attached |
| USDA FPL-GTR-282 | Southern Pine and Douglas-fir plywood; Southern Pine mill 1 and Aspen mill 7 OSB | Selected product bending, shear and internal-bond values with direction and source identity; unknown grade, moisture, layup and temperature remain explicit |
| USDA FPL-GTR-282, Chapter 5 | Waterborne-preservative-treated lumber | Historical processing thresholds and named missing inputs; no treatment-specific strength or thermal-property correction is inferred |
| NTRS 19860021558 | Duran borosilicate glass | Density, elasticity, heat capacity, conductivity, range-mean expansion and temperature ratings; conflicting tensile-strength columns are withheld |
| NIST SRM 1826b | Named soda-lime glass density standard | Whole-specimen density at 20 °C with the certificate's expanded uncertainty; no microhomogeneity or generic-glass certification claim |
| NIST PVC cryogenic page | Low-density PVC foam at three stated densities | Exact 300 K fit evaluations; air-filled and CO2-blown conductivity claims remain distinct from each other and from the specific-heat specimen |

The NASA table curves are piecewise-linear approximations of published values,
not reproductions of the authors' complete fitted laws. No interpolation error
bound or statistical confidence is invented. The 304/316 stainless thermal
curves use 11 fit evaluations over 77–300 K with linear interpolation between
them. The knots were evaluated with 60-digit decimal arithmetic, then rounded
to 15 significant digits. This corrects tiny cancellation errors in the earlier
77/293 K floating-point evaluations (largest change about `4.43e-8 J/kg/K`).
These curves approximate the published polynomial; they do not execute it.
An independent 0.01 K sampling comparison found a largest relative deviation
of about 1.374%, for 304 specific heat near 87.59 K. That is a sampled numerical
discrepancy, not a global bound or additional physical uncertainty. The source
fit-error notes remain separate from the `Unstated` statistical uncertainty.
Stainless modulus now has separate linear approximations over 77–293 K (304)
and 77–294 K (316); relative expansion covers 77–300 K for both. Their knots
use the same high-precision evaluation method. Numerical comparisons at
non-knot temperatures check modulus against its published polynomial within
0.04% and expansion within `1e-5` absolute strain. These are sampled
approximation checks, not global bounds or measurement uncertainty.
The source's small nonzero expansion fit residual at 293 K is preserved.
Relative expansion remains dimensionless. A separate curve now supplies the
model-derived instantaneous coefficient `alpha(T) = epsilon'(T)/(1+epsilon(T))`
in `K^-1`, calculated from the source polynomial and represented by 11 knots
over 77–300 K for each alloy. It is neither a relabeling of relative strain nor
a range-mean coefficient. A 0.01 K sampling comparison found maximum relative
interpolation discrepancy about 0.6981%. This is not a global bound; the source's
5% relative-length fit error provides no derivative-error or statistical
uncertainty bound. Both alloys publish the same expansion coefficients. The
derived curves retain pressure-unknown context and the fitted 293 K residual;
they do not supply missing density/Poisson ratio or qualify a thermal-strain,
loss or heating model.
Other NIST cryogenic point records still refuse temperatures between their
retained points.
The NASA iron conductivity curve stops before the source's duplicated
phase-transition temperature rather than joining two different phases.
The elemental thermodynamic packs likewise keep their solid and liquid tables
separate. Tin's crystalline table specifically describes the tetragonal state;
the table-phase flag cannot automatically select between solid allotropes.

The new pure-metal conductivity tables describe annealed aluminum of
99.996+% purity over 200–930 K and annealed copper of 99.999+% purity over
273–1350 K, with their stated residual resistivities. Their explicit solid
endpoints stop before the separate liquid rows. Pressure remains unstated;
these records cannot automatically join the 1-bar NASA thermodynamic tables
into a condition-matched specimen. NASA-CR-71699 contains no lead conductivity
table, so this tranche does not fill that gap.

`lead-cast-expansion-nbs-rp500` adds five measurements for one named solid
specimen: sample 1144, sand-cast 99.9% lead. Table 2 reports its density as
`11.310 g/cm3` at 25 °C and its first-heating mean expansion coefficients as
`28.3`, `28.6`, `29.5`, and `31.2` × `1e-6 K^-1` from 20 °C to 60, 100,
200, and 300 °C, respectively. The query requires the actual sample, heating
direction and both interval endpoints; it refuses cooling, another reference
temperature, and unmeasured intervals. Density remains an exact-temperature
record with no asserted before/after-heating condition. These are neither
instantaneous expansion coefficients nor a continuous density law. Pressure
remains unstated, and this specimen cannot automatically join the NASA
standard-state thermodynamic records into a complete heating bundle.

`lead-solid-conductivity-nbs-rp668` supplies the separate L.S. reference
standard's published conductivity estimate over 273.15–573.15 K: 35.2, 33.2,
31.2, and 29.2 W/(m K) at 0, 100, 200, and 300 °C with linear interpolation.
The experiment measured temperature dependence against an **assumed absolute
anchor** of 0.352 W/(cm K) at 0 °C. The source's illustrated measurement run
only reaches down to 58 °C; the 0 °C entry is its adopted calibration value.
Queries must acknowledge `source-calibration-assumed=1`. The reported roughly
2% comparison consistency and estimated roughly 3% anchor accuracy remain
source caveats, without statistical confidence semantics. The curve stops
below melting and does not identify the RP500 specimen.

After compiling the two stainless manifests to `304.fsmatpk` and `316.fsmatpk`,
compare their thermal coverage with:

```bash
frankensim --json discover \
  examples/material-discovery/stainless-thermal.json \
  /path/to/304.fsmatpk /path/to/316.fsmatpk
```

Both source bundles supply conductivity and specific heat throughout the
request's 90–297 K interval, with explicit acknowledgement of unknown source
pressure. Missing density, specimen/process identity and executable evolving
models still prevent a complete heating trajectory claim.

Use `examples/material-discovery/stainless-thermomechanical.json` with the same
command to request all four properties over 90–293 K. Both packs satisfy that
declared subset. Extending the upper temperature to 294 K leaves only 316
complete because 304's modulus fit stops at 293 K; at 295 K neither modulus
is supported. The longer thermal range does not extend mechanical coverage.

The CFAST records require explicit acknowledgement of missing source conditions.
Their nominal layer thicknesses describe CFAST layers, and their emissivities
describe model inputs; neither becomes an intrinsic, fully characterized bulk
material property. The entire required [CFAST notice](CFAST_NOTICE.md) is also
carried in the affected manifest license fields and compiled pack provenance.
The original CSV is pinned to commit
`d4a6e17c292851b80329581d0fb6018317ca3269` of `firemodels/cfast`.

For the measured concrete, approximately six months of water curing and the
final 24-hour drying condition are distinct from the 25 °C mechanical test
state. Density is calculated from the reported mass and cylinder dimensions;
its measurement temperature is not supplied. A single specimen does not define
mixture-wide statistics, a moisture law, or a structural design allowable.

The plywood table's proportional-limit SI header conflicts with its paired
imperial values. These records disclose the conflict and use the paired psi
values, rounded to the table's MPa precision. OSB parallel and perpendicular
properties remain separate; neither panel family supplies a full elastic
tensor. Specific gravity is retained without inventing a moisture-dependent
density conversion. The treated-lumber records describe the 2021 handbook's
retention and post-treatment redrying discussion, not current code compliance
or a universal correction for ACQ, CCA, copper azole or other treatments.

Duran's 250 °C thermal-shock rating is a **250 K temperature difference**.
Its 20–300 °C expansion entry is a range mean. The source's tensile-strength
columns disagree: 7.8 MPa is paired with 11,400 psi (about 78.6 MPa), and
3.9 MPa with 5,700 psi (about 39.3 MPa), under a safety factor of five.
Both conflicting pairs remain in the observation, but neither numerical
strength claim is admitted pending corroboration. The SRM 1826b density
interval is an expanded half-width at
approximately 95% confidence, without an inferred probability distribution,
coverage factor or degrees of freedom. Its temperature-correction coefficient
has no supplied range and therefore cannot establish heating coverage.

Lead's NASA records use pure crystalline/liquid standard states at 1 bar and
the source's `207.2 g/mol` molar mass. The crystalline table is retained over
`200..600.65 K` and the liquid table over `600.65..1300 K`. Both enthalpies use
zero for crystalline lead at `298.15 K`; their paired values at `600.65 K`
differ by `4.812 kJ/mol`, or about `23.224 kJ/kg`. Fusion is named
`latent-heat-fusion`, distinct from vaporization or a liquid formation enthalpy.
The phase flag selects a source table, not a mixture or melting-kinetics law.
These independent heat-capacity and enthalpy interpolants are not an exactly
thermodynamically consistent coupled constitutive model.

The additional `lead-liquid-transport-nasa-cr144016` bundle retains density
over 400–1000 °C, conductivity over 400–800 °C, viscosity over 441–844 °C,
and surface tension over 350–500 °C. These are compiled source values with
engineering linear interpolation. Pressure and specimen assay are unstated;
the surface environment is also unstated. The conductivity's historical
`Cal/(m sec °C)` unit has no identified calorie convention, so the pack
explicitly selects 4.184 J/cal and retains that ambiguity in query context.
These records do not supply solid-phase transport, across-melting coverage,
or a condition-matched 1-bar specimen bundle.

Compile and query the new liquid heat-capacity coverage with:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/lead-liquid-nasa-tp3287/manifest.tsv \
  --out /path/to/lead-liquid.fsmatpk
frankensim --json discover \
  examples/material-discovery/lead-liquid-heat-capacity.json \
  /path/to/lead-liquid.fsmatpk
```

This request concerns only heat capacity over `650..800 K`, at the source
pressure and liquid phase. Heating/melting still needs condition-matched
density, conductivity, evolving mechanical/liquid/interface properties and
executable models. The older C447 seed and its original gap demonstration
remain separate historical source records.

Primary source corpus:

- [NASA thermophysical tables, NASA-CR-71699](https://ntrs.nasa.gov/citations/19660014513)
- [Thermodynamic Data for Fifty Reference Elements, NASA/TP-3287/REV1](https://ntrs.nasa.gov/citations/20010021116)
- [Design of Refractory Metal Life Test Heat Pipe and Calorimeter, NASA/TP-2010-216435](https://ntrs.nasa.gov/citations/20100029805)
- [Diffusion in Liquid Metal Systems, NASA-CR-144016](https://ntrs.nasa.gov/citations/19760003214)
- [NIST cryogenic material-property index](https://trc.nist.gov/cryogenics/materials/materialproperties.htm)
- [NBS RP500, cast-lead expansion and density](https://nvlpubs.nist.gov/nistpubs/jres/9/jresv9n5p703_A2b.pdf)
- [NBS RP668, calibrated solid-lead conductivity](https://nvlpubs.nist.gov/nistpubs/jres/12/jresv12n4p429_A2b.pdf)
- [NIST CFAST source thermal library at the retained commit](https://github.com/firemodels/cfast/blob/d4a6e17c292851b80329581d0fb6018317ca3269/Utilities/for_bundle/Bin/thermal.csv)
- [NIST SP 1041, CFAST user's guide](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1041.pdf)
- [NISTIR 6475, concrete mechanical measurements](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6475.pdf)
- [USDA Wood Handbook, FPL-GTR-282](https://www.fpl.fs.usda.gov/documnts/fplgtr/fplgtr282/fpl_gtr282.pdf)
- [NASA contractor fuel-flowmeter report, Duran glass table](https://ntrs.nasa.gov/citations/19860021558)
- [NIST SRM 1826b density certificate](https://tsapps.nist.gov/srmext/certificates/archives/1826b.pdf)
- [NIST PVC cryogenic fits](https://trc.nist.gov/cryogenics/materials/PVC/PVCrev.htm)

Coverage remains incomplete. Missing grade/process, temperature, moisture,
frequency, surface preparation, uncertainty or optical-band information must
not be filled from a superficially similar material. Common-family coverage
and complete condition-specific use cases are different acquisition targets.
In particular, treated lumber still needs measured species/chemical/retention
bundles, and the added PVC data do not describe ordinary solid PVC pipe.
Broader alloy grades, glass optical spectra, plastic formulations and
continuous temperature/moisture/frequency coverage remain acquisition work.

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

## NIST SRM 1720 northern-continental-air tranche

`nist-srm-1720-northern-continental-air/` retains the information-value table
for NIST SRM 1720 Lot `1720-A-XX`. NOAA filled the lot from whole ambient air
at Niwot Ridge, Colorado, in 2010 and scrubbed it of moisture to represent
northern continental air. The table reports oxygen at `20.93 mol%`, argon at
`0.935 mol%`, carbon monoxide from `130` through `180 nmol/mol`, and nitrogen
as the balance gas.

NIST explicitly classifies those four entries as information values that
cannot establish metrological traceability. The carbon-monoxide endpoints are
therefore separate observed-range claims, not a confidence interval. The
generic archived certificate leaves the cylinder-specific certified carbon-
dioxide, methane, and nitrous-oxide values and expanded uncertainties as
`SAMPLE` placeholders; no runtime claims are invented for them. Residual
humidity, use-state temperature and pressure, remaining trace composition,
and covariance also remain unknown.

The certificate was issued on 2015-09-08 and expired on 2020-09-01. The pack
is a historical, named natural-air lot, not a current certificate, universal
dry-air composition, or standard-atmosphere model.

Northern-continental-air sources:

- <https://tsapps.nist.gov/srmext/certificates/archives/1720.pdf>
- <https://www.nist.gov/copyrights-disclaimers>

## NIST SRM 2728 auto-emission reference-gas tranche

`nist-srm-2728-auto-emission-reference-gas/` retains the named composition in
NIST's archived generic Certificate of Analysis for SRM 2728, *Carbon Dioxide,
Carbon Monoxide, and Propane in Nitrogen (Auto Emission Inspection/Maintenance
Gas)*. The three title values are `14 mol%` carbon dioxide, `8 mol%` carbon
monoxide, and `3000 micromol/mol` propane, with nitrogen named as the balance
gas. The pack normalizes those title values to dimensionless amount fractions
but leaves nitrogen as observation text: deriving a closed balance would be
dishonest because the certificate also reports an `8 micromol/mol` estimated
total-other-hydrocarbons information value.

The archived generic certificate leaves the per-cylinder certified-value and
95%-confidence-interval fields blank. The runtime claims therefore remain
explicitly nominal or informational with `Unstated` uncertainty; they are not
represented as certified values for a physical cylinder. Cylinder/sample
identity, use temperature and pressure, and oxygen, water, and NOx fractions
are also absent.

This is a versioned synthetic reference mixture intended to calibrate mobile-
source emissions measurements. It is not engine-generated exhaust, a
combustion or equilibrium state, or authority for thermophysical, reaction, or
transport properties. NIST's unmarked public-information reuse policy is
retained with attribution.

Reference-gas sources:

- <https://tsapps.nist.gov/srmext/certificates/archives/2728.pdf>
- <https://tsapps.nist.gov/srmext/certificates/archives/>
- <https://www.nist.gov/copyrights-disclaimers>

## Aluminum 6061-T6 cryogenic curves

`aluminum-6061-t6-cryogenic/` is the first committed bulk-material tranche.
It identifies the condition as Aluminum 6061-T6, UNS A96061, and retains
thermal conductivity, specific heat capacity, and Young's modulus as three
24-knot curves over `77–293 K`. Knots evaluate the polynomials printed on
NIST's Aluminum 6061-T6 material-property page with 60-digit decimal arithmetic,
rounded to 15 significant digits. The `77 K` and `293 K` anchors remain;
tiny specific-heat endpoint corrections remove prior floating-point cancellation
error. Between knots, the declared model is piecewise linear. The v3 manifest
retains absolute-temperature semantics and requires `source-pressure-known=0`:
the source did not specify pressure, which does not mean vacuum or atmosphere.

NIST reports curve-fit errors relative to the underlying data of `0.5%` for
thermal conductivity, `5%` for specific heat, and `1%` for Young's modulus.
Those figures do not state a confidence level or degrees of freedom, so the
source records retain them as observation caveats and explicitly encode the
runtime uncertainty as `Unstated`. They are not laundered into statistical
confidence intervals.

For the declared grid (`77`, every `10 K` from `80` through `290`, then `293 K`),
99 interior probes per segment measured maximum relative interpolation
discrepancies of `0.0858%` for conductivity, `0.1840%` for heat capacity and
`0.003558%` for modulus against the published equations. These sampled
discrepancies are separate from NIST's fit errors; they are not certified global
bounds or new experimental data. The curves refuse outside `77–293 K` and do
not provide the density, Poisson ratio or thermal expansion of a complete
heating specimen.

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

## OFHC Copper cryogenic curves

`ofhc-copper-rrr100/` retains NIST's combined OFHC Copper scope, UNS
C10100/C10200, without pretending that the two UNS designations are
interchangeable for properties the source does not cover. Thermal-conductivity
claims pin `RRR = 100`, a load-bearing condition because NIST publishes
distinct conductivity correlations for RRR 50, 100, 150, 300, and 500. The
specific-heat source publishes one combined OFHC correlation without an RRR
condition, so its observation names that omission explicitly instead of
silently borrowing the conductivity condition.

The tranche evaluates the published RRR-100 conductivity correlation and
specific-heat polynomial at the same 24 temperatures over `77–293 K` as the
aluminum curves, retaining both original endpoint temperatures and using
60-digit decimal evaluation rounded to 15 significant digits. Its declared
linear interpolant has sampled maximum discrepancies of `0.5555%` for
conductivity and `0.2766%` for heat capacity (99 interior probes per segment).
These are observed discrepancies from the equations, not global error bounds.
The typed query also requires `source-pressure-known=0`. NIST
reports `1%` curve-fit error for the RRR-100 conductivity correlation and `5%`
for specific heat at these temperatures, but does not give confidence levels
or degrees of freedom. The pack therefore retains the errors and complete
coefficient sets as observation caveats while encoding `Unstated` statistical
uncertainty.

NASA-CR-134806 independently reports typical room-temperature OFHC Copper
specific heat of `386 J/(kg K)` and thermal conductivity of `390 W/(m K)`.
Both are within `2%` of the NIST-derived `293 K` values. That coarse G3 check
is evidence only: it neither overwrites the NIST claims nor erases NASA's
distinct material/processing context. The NTRS record marks the NASA report
public and as a work of the U.S. Government whose public use is permitted.

The [cryogenic aluminum/copper discovery request](../../../examples/material-discovery/cryogenic-aluminum-copper.json)
asks for both conductivity and heat capacity over `90–290 K` with the source's
pressure-unstated context. Compile these two manifests and pass both packs to
`fsim --json discover <request.json> <aluminum.fsmatpk> <copper.fsmatpk>`.
Complete property coverage remains unbound to a named specimen; it does not
establish same-specimen compatibility or a transient thermal model.

OFHC Copper references:

- <https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm>
- <https://www.nist.gov/publications/properties-copper-and-copper-alloys-cryogenic-temperatures>
- <https://ntrs.nasa.gov/citations/19750021165>
- <https://ntrs.nasa.gov/api/citations/19750021165/downloads/19750021165.pdf>

## PTFE/Teflon cryogenic exact-point tranche

`ptfe-teflon-nist-cryogenic/` retains the material identity at exactly the
specificity published by NIST: the page is headed `Teflon`, while the index
places it among solid materials from cryogenic to room temperatures. The
source does not identify a resin grade, crystallinity, filler, processing
history, density, or supplier, so the observation names those omissions rather
than silently treating every PTFE formulation as interchangeable.

The tranche evaluates NIST's displayed thermal-conductivity and specific-heat
polynomials at exactly `77 K` and `293 K`. The resulting four scalars use
degenerate temperature validity intervals and therefore do not expose the
underlying polynomial as a continuous constitutive curve. NIST states a
`4..300 K` data and equation range and reports curve-fit errors relative to the
data of `5%` for thermal conductivity and `1.5%` for specific heat. Because the
page supplies no confidence level or degrees of freedom, those figures remain
observation caveats and every runtime uncertainty is explicitly `Unstated`.

The page also publishes Young's-modulus and linear-expansion correlations.
They are deliberately excluded from this first PTFE slice: the rendered page
does not expose the Young's-modulus coefficient table, and linear expansion is
a reference-length-relative quantity rather than the linear-expansion-
coefficient property needed by most thermal models. No value is reconstructed
from an unavailable plot or relabeled into a different physical quantity.

PTFE/Teflon references:

- <https://trc.nist.gov/cryogenics/materials/Teflon/Teflon_rev.htm>
- <https://trc.nist.gov/cryogenics/materials/materialproperties.htm>
- <https://www.nist.gov/copyrights-disclaimers>

## NASA LaRC PEEK THERMIC-plate tranche

`peek-nasa-thermic-plate/` pins NASA/TM-20210014330's instrumented PEEK test
plate rather than treating the polymer family name as a universal card. LaRC
removed three `7.9 mm`-diameter by `12.7 mm` plugs from the plate, instrumented
each with three Type K thermocouples through its thickness, and used the center
plug's transient temperatures with a one-dimensional Crank-Nicolson model and
Continuous Genetic Algorithm. The source does not report the PEEK grade,
supplier, crystallinity, filler, or processing history, so all observations
carry that missing identity explicitly.

The tranche evaluates NASA's final all-thermocouple thermal-conductivity cubic
and its LaRC DSC specific-heat cubic at exactly `300 K`, `400 K`, `500 K`, and
`525 K`. It also retains the `1264 kg/m3` commercial-laboratory density used by
the thermal model, with a fail-closed flag because the density measurement's
temperature is not reported. Equation 2 does not print a specific-heat unit
inline; `J/(kg K)` is retained because it is the SI unit required by Equation
1's dimensional balance with the paper's stated conductivity, density, length,
and time units. The source describes conductivity testing at atmospheric
pressure without assigning a numeric pressure, so conductivity claims carry
`source_pressure_atmospheric = 1` instead of inventing `101325 Pa`.

NASA's Abstract and Results state a `300..525 K` range, while the Concluding
Remarks say `300..550 K`; the pack deliberately adopts the narrower repeated
range. Four thermocouple-subset conductivity estimates differed by about `3%`,
and the chosen CGA fit was about `35%` below the prior laser-flash series. Those
comparisons are diagnostics, not confidence intervals, so all nine runtime
claims retain `Unstated` uncertainty. Exact-temperature validity prevents the
offline polynomial evaluations from becoming an implicit continuous curve.

PEEK THERMIC references:

- <https://ntrs.nasa.gov/citations/20210014330>
- <https://ntrs.nasa.gov/api/citations/20210014330/downloads/NASA-TM-20210014330.pdf>

## AISI 4140 Rockwell C33 low-temperature tranche

`aisi-4140-rc33/` is a condition-pinned structural-steel tranche from
NASA-TM-X-64791, not a generic 4140 card. The retained specimen is QQ-S-624
heat `137M186`, machined from nominal one-inch bar. NASA normalized the stock
at `899 degC`, hardened it at `851 degC`, oil quenched it, tempered it at
`566 degC`, and measured Rockwell C33. Those processing and hardness fields
are part of every observation identity.

At the report's exact `26.7 degC` and `-73 degC` test points, the pack retains
ultimate tensile strength, 0.2%-offset yield strength, two-inch tensile
elongation, reduction of area, MIL-STD-151 Charpy V-notch impact energy, and
ultimate and yield double-shear strengths. The NASA table reports averages
from five tensile specimen sets, four impact tests at each retained point, and
four shear specimens. It does not publish dispersion or confidence metadata,
so every runtime uncertainty is explicitly `Unstated`; sample counts are
observation caveats, not fabricated error bars.

The G3 pack battery checks NASA's redundant ksi/GPa columns for the retained
tensile values within their printed rounding precision. This is a source-unit
transcription check, not an independent-source agreement claim. No second
source with a sufficiently identical heat, geometry, process schedule,
hardness, and temperature condition is treated as interchangeable.

AISI 4140 references:

- <https://ntrs.nasa.gov/citations/19740002417>
- <https://ntrs.nasa.gov/api/citations/19740002417/downloads/19740002417.pdf>

## AISI 1045 cold-drawn tensile tranche

`aisi-1045-cold-drawn/` retains one experimental tensile series from Kang and
Lee's 2024 turning study, not a generic handbook card for AISI 1045. The source
identifies an approximately `0.45 wt% C` cold-drawn AISI 1045 bar with `37 mm`
outside diameter and `102 mm` starting length. Three ASTM E8 cylindrical
specimens used a `50 mm` gauge length, `12.5 mm` test-section diameter,
extensometry, and `10 mm/min` crosshead speed.

The pack retains the paper's reported means: `550.51 MPa` yield strength,
`695.31 MPa` ultimate tensile strength, and `14.1%` elongation over the
`50 mm` gauge length. The source also prints all three replicate values for
each property. The pack derives symmetric 95% Student-t half-widths for the
mean from those replicates using `t(0.975, df=2) = 4.302652729911275` under an
explicit iid-normal assumption: `21.7281028138069 MPa`,
`29.0129099774989 MPa`, and `0.943972330500626` percentage points,
respectively. These are small-sample confidence intervals for the reported
series, not grade tolerances or certified material allowables. Although the
three properties were measured on paired specimen rows, the pack does not
invent a joint covariance block.

The paper does not report tensile-test temperature. Every claim therefore
requires both the exact normalized crosshead-speed point and the explicit
dimensionless context `source_test_temperature_known = 0`. That second axis is
a fail-closed acknowledgement of missing source metadata, not a physical
temperature. No `temperature` validity interval is present, and a caller must
not interpret these values as temperature-independent. The paper is published
under Creative Commons Attribution 4.0 International; the manifest retains the
authors, title, journal, DOI, and license.

The paper's printed Vickers-hardness replicates round to a mean different from
the mean printed in the same table, so this tranche refuses to select either
number. No second source with the same bar, cold-drawn state, specimen geometry,
test rate, and unreported temperature is treated as interchangeable.

AISI 1045 reference:

- <https://www.mdpi.com/2227-9717/12/6/1171>

## AISI 52100 CVM bearing-steel tranche

`aisi-52100-cvm-hot-hardness/` pins the bearing steel to the single
consumable-vacuum-melted ingot studied in NASA-TN-D-6632. Table I reports the
actual composition as `0.96% C`, `0.22% Si`, `0.36% Mn`, `0.012% S`,
`0.007% P`, and `1.36% Cr`, balance iron. The analytical method, heat
identifier, and composition uncertainties are not reported, so those six
claims carry `Unstated` uncertainty and do not infer missing trace values.

Table II gives a common process spine: austenitize at `1116..1144 K` for
`30 min`, oil quench at `325 K`, then first-temper at `394 K` for `60 min`.
Five separately retained states use a `60 min` second temper at `505 K`,
`450 K`, `433 K`, or `394 K`, or no second temper. Their exact `294 K`
Rockwell C scale readings are `59.7`, `62.3`, `63.4`, `64.6`, and `65.1`.
NASA used a standard Rockwell C tester with a `150 kg` load and diamond
indenter in a low-oxygen furnace, taking at least two readings per material and
temperature. The table does not publish dispersion, so these measurements
remain `Unstated`; the report's separate `+/-1` Rockwell-point accuracy claim
for its predictive equation is not laundered into measurement uncertainty.

X-ray diffraction reports retained-austenite volume fractions of `12.8%`,
`15.6%`, `18.4%`, and `11.8%` for the `450 K`, `433 K`, `394 K`, and
no-second-temper states. The `505 K` state is printed only as `<2%`; that
censored result remains an observation caveat and never becomes an exact
scalar. Hardness and retained-austenite observations stay method-separated,
with no invented covariance. NASA reports ASTM grain size 12 for the 52100
material.

NASA's NTRS record marks the report public and a work of the U.S. Government
whose public use is permitted. The manifest retains NASA attribution.

AISI 52100 reference:

- <https://ntrs.nasa.gov/citations/19720007811>
- <https://ntrs.nasa.gov/api/citations/19720007811/downloads/19720007811.pdf>

## AISI 9310 CVM carburized gear-steel tranche

`aisi-9310-cvm-carburized/` pins one consumable-electrode-vacuum-melted AISI
9310 spur-gear lot from NASA-TM-104352. All test gears came from one heat. The
pack retains Table I's nominal composition by weight: `0.10% C`, `0.63% Mn`,
`0.27% Si`, `3.22% Ni`, `1.21% Cr`, `0.12% Mo`, `0.13% Cu`, `0.005% P`, and
`0.005% S`, balance iron. These are nominal grade values, not a reported
analysis of the tested heat, so the observation says so explicitly and every
composition uncertainty remains `Unstated`.

Table II supplies the process identity: carburize at `1172 K` for `8 h`, air
cool, copper plate, reheat at `922 K` for `2.5 h`, air cool, austenitize at
`1117 K` for `2.5 h`, oil quench, subzero treat at `180 K` for `3.5 h`, double
temper at `450 K` for `2 h` each, finish grind, and stress relieve at `450 K`
for `2 h`. The detailed Test Materials text reports a Rockwell C58 case,
Rockwell C40 core, and `0.97 mm` case depth.

The same report's abstract and summary instead describe the gears as hardened
to Rockwell C60. The pack preserves C58 and C60 as two separate
`case_rockwell_c_scale_reading` claims, each linked to the exact source
statement that supports it. It does not average, silently select, or fabricate
uncertainty around the discrepancy. Rockwell C values are named empirical
scale readings in dimensionless storage, not ratio quantities.

As G3 plausibility evidence only, NASA SP-410 reports a different VAR AISI 9310
gear lot with nominal Rockwell C62 case, Rockwell C45 core, and `1 mm`
effective case depth. Those nearby values can catch order-of-magnitude
transcription errors, but the different melt route, lot, and process prevent
fusion or substitution. NASA's NTRS record marks NASA-TM-104352 public and a
work of the U.S. Government whose public use is permitted.

AISI 9310 references:

- <https://ntrs.nasa.gov/citations/19910020285>
- <https://ntrs.nasa.gov/api/citations/19910020285/downloads/19910020285.pdf>
- <https://ntrs.nasa.gov/citations/19750018303>

## NASA/NAPC polyol-ester gear-oil tranches

`napc-pe-5-l-1274-gear-oil/` and
`napc-pe-5-l-1307-1553-gear-oil/` retain application-specific oil evidence
from NASA-TM-104352. PE-5-L-1274 is NASA identification A, the program's
reference polyol-ester synthetic gear oil. PE-5-L-1307 and PE-5-L-1553 are
NASA identification B: two tested batches of the same polyol-ester lubricant,
reported as meeting MIL-L-23699. The batch codes remain visible rather than
being collapsed into a single averaged row.

For PE-5-L-1274, Table IV supplies an exact `516 K` flash point, reported
specific gravity `0.998` at `298 K`, and total acid number
`0.07 mg KOH/g oil`. Its pour point is printed only as `<200 K`, so that
censored value remains an observation caveat rather than an exact scalar. For
the two B batches, the pack separately retains flash points of `539 K` and
pour points of `220 K` (PE-5-L-1307) and `213 K` (PE-5-L-1553). Values spanning
both B columns remain common-formulation claims: reported specific gravity
`1.000` at `289 K` and total acid number `0.03 mg KOH/g oil`.

The table also prints kinematic-viscosity numbers at `311 K` and `373 K`, but
does not state their unit. The source observations preserve those numbers and
temperatures verbatim, while the normalized packs deliberately contain no
`kinematic_viscosity` claim. Supplying a conventional unit from outside the
source would violate the Five Explicits. Measurement standards, replicates,
dispersion, and confidence metadata are likewise absent, so all admitted
claims carry `Unstated` uncertainty.

NASA says the additive and chemical data are proprietary to the respective
manufacturers. These packs therefore identify the oils only by the public
NASA/NAPC codes, basestock class, reported specification relation, and measured
Table IV properties. They do not claim a redistributable formulation.

NASA/NAPC gear-oil reference:

- <https://ntrs.nasa.gov/citations/19910020285>
- <https://ntrs.nasa.gov/api/citations/19910020285/downloads/19910020285.pdf>

## Rheolube 2000 Pennzane bearing-grease tranche

`rheolube-2000-pennzane-grease/` pins the named rolling-bearing grease in Paul
A. Bessette's NASA-CP-3350 paper. The source identifies Rheolube 2000 as the
then-available grease prepared from Pennzane SHF X-2000 multiply alkylated
cyclopentane base oil. Its sodium octadecylterephthalamate thickener is
described as approximately `20%` of the formulation and as suitable for
high-speed rolling-element-bearing applications.

Table 7 contributes three bounded claims: NLGI consistency grade `2`, density
`0.89 g/cm3` at exactly `25 degC`, and oil separation `3.3%` after exactly
`24 h` at `100 degC`. NLGI grade is a named empirical consistency scale, not a
ratio quantity. The table labels its values typical and reports no lot,
manufacturing or vacuum-hardening state, measurement standards, replicates,
dispersion, or confidence metadata, so all three claims carry `Unstated`
uncertainty.

The prose gives an approximate `260 degC` dropping point and the table prints
unworked/60-stroke worked penetrations without a unit or method. Those values
remain observation-only. The wear-scar, oxidation-pressure-drop, and
extrapolated vapor-pressure rows also stay outside the bulk card: the first two
are conditioned test-system outcomes, while the PDF text does not preserve an
unambiguous vapor-pressure exponent. Color, odor, and ultrafiltration status
remain descriptive provenance rather than quantitative claims.

The NASA proceedings record is public with public use permitted. The manifest
retains the paper, proceedings, table, and NTRS identifier.

Rheolube 2000 reference:

- <https://ntrs.nasa.gov/citations/19970021613>
- <https://ntrs.nasa.gov/api/citations/19970021613/downloads/19970021613.pdf>

## Pennzane SHF X-2000 aerospace bearing-oil tranche

`pennzane-shf-x-2000-bearing-oil/` pins the multiply alkylated cyclopentane
fluid behind the Rheolube grease to its public Pennzane SHF X-2000 identity.
Bessette identifies it chemically as Tris(2-octyldodecyl) cyclopentane with an
approximate molecular weight of `910 g/mol` and discusses the fluid as an
aerospace lubricant compatible with antioxidant and boundary-lubricant
additives.

Table 6 contributes three exact-temperature kinematic-viscosity points:
`80,500 mm2/s` at `-40 degC`, `107 mm2/s` at `40 degC`, and `14.3 mm2/s` at
`100 degC`. It also supplies viscosity-index scale reading `137`, flash point
`300 degC`, pour point `-55 degC`, and density `0.84 g/mL` at `25 degC`.
Viscosity index is retained as a named empirical scale reading. The table calls
these typical properties and provides no lot, additive state, measurement
standards, replicates, dispersion, or confidence metadata, so all seven claims
carry `Unstated` uncertainty and no continuous viscosity curve is inferred.

The printed `8e-4 cc/cc/degC` expansion coefficient remains observation-only
because the source-unit grammar has no degree-Celsius interval token; silently
substituting affine `degC` would be dimensionally wrong. The conditioned
`0.34 mm` wear scar is a tribometer-system outcome, not a bulk property. Table
5 vapor pressures remain excluded because they are lot-specific and the
retained PDF text extraction does not preserve their superscript exponents
unambiguously.

Pennzane SHF X-2000 reference:

- <https://ntrs.nasa.gov/citations/19970021613>
- <https://ntrs.nasa.gov/api/citations/19970021613/downloads/19970021613.pdf>

## S2-S pearlitic gray-cast-iron tranche

`gray-cast-iron-s2-s/` retains one engine-relevant experimental ingot from
Wang et al., not a universal gray-iron card. Sample S2-S was melted in a
`500 kg` medium-frequency induction furnace from `70 wt%` steel scrap and
`30 wt%` pig iron, superheated to `1530 degC`, transferred to a ladle holding
`0.4 wt%` Sr-FeSi inoculant (`Ino_2`, itself `2.0 wt% Sr`), and poured into an
EN-1561 Type II mould. The resulting matrix was fully pearlitic with type-A
graphite.

Table 1 pins the sample to `3.54% C`, `1.62% Si`, `0.51% Mn`, `0.025% P`,
`0.028% S`, `0.35% Mo`, `0.58% Cu`, and `0.060% Sn`, with reported carbon
equivalent `4.05%`. The compiler battery independently checks the paper's
`CE = C + 0.31 Si + 0.33 P` formula within printed rounding. No analytical
uncertainty, balance-iron scalar, or unreported trace value is invented.

Table 2 contributes four microstructure means: `9.0%` graphite area,
`273 um` maximum graphite-flake length, `15.6%` primary-dendrite area, and
`371 cm^-2` eutectic-colony density. The paper labels their error bars as one
standard deviation and reports eight cross-section fields, but provides no
confidence level or field-level observations. The SD values remain explicit
observation caveats while runtime uncertainty stays `Unstated`.

Figure 8 is the only published numerical presentation of strength and thermal
conductivity for S2-S. The retained centers are deliberately low precision:
`326 MPa` ultimate tensile strength and `58.8 W/(m K)` room-state thermal
conductivity, digitized from the authors' `2141 x 1490` PNG to the nearest
`1 MPa` and `0.1 W/(m K)`. The displayed one-SD bars are approximately
`8 MPa` and `0.3 W/(m K)`, but are not confidence intervals and therefore do
not become runtime half-widths. The tensile series averaged three
GB/T T228.1-2010 specimens; the thermal path used NETZSCH LFA 457 results and
Archimedes density. Neither path reports an exact test temperature, so both
claims require `source_test_temperature_known = 0` rather than implying a
temperature validity interval.

As broad G3 plausibility evidence only, ORNL/TM-2012/506 Appendix C reports a
generic gray-cast-iron thermal-conductivity range of `42..62 W/(m K)`. The
S2-S digitization falls inside that range, but the ORNL row is not
condition-matched and does not overwrite or fuse with the primary claim. The
article is CC-BY-4.0 and the manifest retains full attribution.

Gray-cast-iron references:

- <https://www.mdpi.com/1996-1944/11/10/1876>
- <https://www.osti.gov/servlets/purl/1148409>

## NASA-CR-115153 inhibited water/ethylene-glycol coolant tranche

`nasa-cr-115153-water-ethylene-glycol/` retains the exact inhibited coolant
specified in Table 5 of *A Fundamental Study of Sublimation Through a Porous
Surface*, not a generic water/glycol mixture. The source identifies a
water/ethylene-glycol solution containing `0.10..0.25 wt%` sodium nitrite,
`1.33..1.57 wt%` sodium benzoate, and `36..38.5 wt%` water. Each printed
formulation endpoint is a separate bound claim. No midpoint, exact batch
composition, or ethylene-glycol balance is inferred.

Table 5 also reports density `67.5 lbm/ft^3`, thermal conductivity
`0.22 BTU/(hr ft degF)`, and the approximate relation
`cp = 0.67 + 0.008 T_degF BTU/(lbm degF)` over `0..100 degF`. The pack stores
deterministically converted SI density and conductivity plus heat-capacity
evaluations at exactly `0`, `50`, and `100 degF`; it does not expose a
continuous law. The source omits density/conductivity test temperature,
pressure, methods, dispersion, confidence metadata, and the BTU convention.
Those gaps remain fail-closed validity flags, every runtime uncertainty is
`Unstated`, and the extra SI digits preserve conversion reproducibility rather
than claiming source precision. The declared normalization uses the exact
international pound and foot plus `Btu_IT = 1055.05585262 J`; callers still
must acknowledge that the report itself does not distinguish International
Table from thermochemical BTU.

As coarse G3 comparison evidence only, NASA/TM-2019-220019 Table VIII lists a
separate, composition-basis-unspecified `50-50 water ethylene glycol` fluid at
`1082 kg/m^3` and `0.402 W/(m K)`. Those figures are within `0.1%` and `6%` of
the retained density and conductivity transcriptions, respectively. They are
not condition-matched observations and do not overwrite, fuse with, or assign
uncertainty to the NASA-CR-115153 claims. Both NTRS records are public and
permit public use; the manifest retains the primary report attribution.

Water/ethylene-glycol references:

- <https://ntrs.nasa.gov/citations/19710026875>
- <https://ntrs.nasa.gov/api/citations/19710026875/downloads/19710026875.pdf>
- <https://ntrs.nasa.gov/citations/20190001449>
- <https://ntrs.nasa.gov/api/citations/20190001449/downloads/20190001449.pdf>

## N0602-001 nitrile O-ring JP-8 compatibility tranche

`n0602-001-nitrile-jp8-compatibility/` retains the named N0602-001 nitrile
rubber O-ring from Graham et al.'s NASA Seal/Secondary Air System Workshop
study. It does not generalize the source's results to every nitrile or Buna-N
compound. The source omits the compound formulation, supplier, cure, hardness,
lot, and O-ring dimensions; those identity gaps remain explicit in every
observation.

The study used thermogravimetric analysis to estimate dry-material
semi-volatiles, direct thermal-desorption GC-MS to measure absorbed fuel and
fuel/polymer partitioning, and optical dilatometry to measure volume swell.
Its matrix covered nine JP-8 fuels, one Fischer-Tropsch fuel, `25`, `37.5`,
`50`, and `75%v/v` FT blends, and small `1..5 mg` specimens in `1..10 mL` of
fuel. The pack retains `10.1 wt%` TGA semi-volatiles; `8.7` and `27.9%v/v`
absorbed fuel at `0` and `25%v/v` aromatic content; alkane and aromatic
partition coefficients `0.120` and `0.412`; their printed ratio `3.4`; and the
volume-swell regression's `R^2 = 0.948`.

The source is internally inconsistent about the volume-swell slope: its
summary row prints `0.451`, while the adjacent plot prints
`y = 0.463 x - 1.167`. Both slopes remain separate, provenance-linked claims;
the intercept is explicitly a regression parameter, neither slope is averaged
or preferred, and no raw-point refit is invented. The source's `~57%` overlap
between 90% prediction intervals remains observation-only because it is
printed as approximate. Exposure temperature, duration, allocation, raw
points, coefficient uncertainty, and prediction limits are absent, so the
runtime claims use `Unstated` uncertainty and fail-closed missing-condition
axes.

The NTRS record marks the conference paper public with public use permitted.
No condition-matched independent source exposes the same N0602-001 formulation
and fuel matrix, so this tranche records no synthetic agreement band.

N0602-001 reference:

- <https://ntrs.nasa.gov/citations/20080003822>
- <https://ntrs.nasa.gov/api/citations/20080003822/downloads/20080003822.pdf>

## NASA-TN-D-8184 M-19 material-deck tranche

`nasa-tn-d-8184-m19-material-deck/` retains the public M-19 silicon-steel
example that Bollenbacher used to document NASA Lewis's induction-motor
analysis program. Figure 10 prints `14` magnetization points. The pack converts
magnetizing force from ampere-turns per inch to `A/m` and magnetic flux density
from kilolines per square inch to `T`, then stores those exact printed points as
a tabulated-only curve. The report's program reconstructed between points with
straight lines on semilog paper; the current runtime has no semilog policy, so
the pack deliberately does not substitute ordinary linear interpolation.

Figure 11 and the printed `FELOSS` card supply a separate fixed-condition
frequency-loss approximation: `9.4 W/lb` at `400 Hz`,
`64.5 kilolines/in^2`, and `0.014 in` lamination thickness, with frequency
exponent `1.47`. Exact SI normalization yields
`20.723452645378494 W/kg`, `0.9997519995039992 T`, and `0.0003556 m`.
The report defines only
`WFE = WCORE * (frequency/FCORE)^SLOPE` at that fixed flux density and
thickness. It supplies neither a flux-density exponent nor a separable
Steinmetz coefficient, so the pack does not label this as a complete
Steinmetz law.

The source identifies no manufacturer, chemistry, fully/semi-processed state,
anneal, magnetic test method, specimen preparation, waveform, direction,
temperature, repeats, dispersion, or measurement uncertainty. It also does not
bind Figure 10's magnetization curve to a lamination thickness. Those gaps are
explicit fail-closed axes and all runtime uncertainties are `Unstated`. This is
therefore a source-faithful public M-19 material-deck tranche, not the fully
pinned M-19 design card required by bead `1sxe`; process/test-method-matched
curve and loss evidence remains open.

NASA-TN-D-8184 M-19 references:

- <https://ntrs.nasa.gov/citations/19760013282>
- <https://ntrs.nasa.gov/api/citations/19760013282/downloads/19760013282.pdf>

## NASA-CR-4538 Tempel 24N208 M19 rating tranche

`nasa-cr-4538-tempel-24n208-m19/` retains the separate M19 state selected for
the laminated return path of Longsworth's split-Stirling cryocooler compressor.
Unlike the anonymous NASA-TN-D-8184 material deck, this report identifies the
manufacturer as Tempel Steel Company and the product as `24N208`, and describes
it as `0.025 in` thick, non-oriented, annealed, nominal `3%` silicon steel of
AISI M19 grade. Exact normalization stores the thickness as `0.000635 m` and
the nominal silicon fraction as `0.03`.

Section 3.1.3 prints one manufacturer-data rating: `2.08 W/lbm` at `15 kG`,
`60 Hz`, and `50% w/ the grain`. Exact SI normalization gives
`4.585615053445453 W/kg` at `1.5 T`. The report does not define the
`50% w/ the grain` phrase, so the pack retains it as a dimensionless source
condition rather than generalizing it to a direction-independent property.

The heading and surrounding text discuss hysteresis loss, but the report does
not state whether the manufacturer rating is hysteresis-only or total core
loss. The property is therefore deliberately named
`specific_hysteresis_loss_rating`, not `specific_core_loss`. Its bound
semantics, test method, excitation waveform, test temperature, repeats,
dispersion, and uncertainty are absent. The manufacturer sheet/revision,
product lot, heat chemistry, anneal schedule, and surface insulation state are
also missing. The report authors' later `8 W` core and `10 W` shell assembly
estimates are model outputs and are not admitted as material claims.

This pack remains separate from NASA-TN-D-8184. It improves manufacturer,
product, orientation class, nominal composition, annealed-state, and thickness
identity, but supplies no B-H curve, anneal schedule, magnetic test method, or
complete Steinmetz model. The fully pinned M19 card required by bead `1sxe`
therefore remains open.

NASA-CR-4538 Tempel 24N208 references:

- <https://ntrs.nasa.gov/citations/19940009424>
- <https://ntrs.nasa.gov/api/citations/19940009424/downloads/19940009424.pdf>

## Torrent 2018 M19/M290-50A reported loss-fit inputs

`torrent-2018-m19-steinmetz-inputs/` retains the numeric inputs printed for two
M-19 loss equations in Torrent, Perat, and Jiménez's motor study. Section 3.2
maps AISI M-19 to EN 10106 M290-50A and states that the manufacturer loss data
span `50..1000 Hz` and `0.1..1.5 T` under sinusoidal excitation. The Equation 5
state uses a `0.5 mm` sheet.

Equation 4 prints `k_h = 4.8`, frequency exponent `a = 1.2`, induction
exponent `n = 2`, and an output scale of `10^-2`. Equation 5 prints `k_f = 60`,
frequency exponent `x = 2.05`, sheet-thickness exponent `y = 2`, induction
exponent `z = 2.05`, and an output scale of `10^2`. The pack stores these as
literal source-equation numeric fields. In particular, the `k_h` and `k_f`
numbers are not asserted to be dimensionless physical coefficients: their
compound units and semantics depend on the printed equations. No executable
Steinmetz or loss model is synthesized from them.

Reference 13 names Cogent Electrical Steel as the manufacturer-data source, but
the article does not identify a Cogent product designation or revision,
production site, heat, chemistry, processing, anneal, coating, rolling
direction, magnetic sheet-loss test method, test temperature, repeats,
dispersion, fit covariance, or B-H curve. It also reports larger
analytical-versus-experimental differences at high flux density. These inputs
therefore remain separate from both the anonymous NASA-TN-D-8184 material deck
and the Tempel `24N208` state. They do not close the fully pinned M19 card
required by bead `1sxe`.

Torrent 2018 M19/M290-50A references:

- <https://doi.org/10.3390/en11061549>
- <https://creativecommons.org/licenses/by/4.0/>

## NGYC N42 sintered NdFeB magnet tranche

`ngyc-n42-sintered-nickel-coated/` retains the N42 sintered NdFeB magnets
identified by Telfah et al. for their multilayer Halbach array. The source pins
the supplier to Ningbo Permanent Magnetic Materials Ltd. (NGYC), Yinxian,
Ningbo; a nickel coating; and cube side lengths of `20`, `10`, `6`, and `3 mm`
with `+/-0.05 mm` dimensional tolerance. It does not provide a production lot,
chemistry, sintering schedule, coating stack or thickness, or allocation of
each reported magnetic value to a particular cube size.

The article reports remanence `1350 mT`, coercivity `923 kA/m`, and maximum
energy product in two mutually inconsistent forms: `318.3 kJ/m^3` and
`42 MGOe`. Exact unit normalization with `1 G = 1e-4 T` and
`1 Oe = 1000/(4*pi) A/m` maps `42 MGOe` to
`334.2253804929802 kJ/m^3`, not `318.3 kJ/m^3`. The pack therefore retains
the printed SI value and normalized CGS value as separate claims with separate
observations. It does not average them, silently correct either representation,
or reinterpret the discrepancy as uncertainty.

The source does not say whether the three values were measured by the authors
or copied from supplier nominal data, and it gives no magnetic test method or
test temperature. Every claim therefore carries `Unstated` uncertainty plus
fail-closed `source_magnetic_test_temperature_known = 0` and
`source_magnetic_test_method_known = 0` axes. No intrinsic coercivity, recoil
permeability, demagnetization curve, temperature coefficient, irreversible
loss boundary, or maximum operating temperature is admitted. Those omissions
mean this is a supplier-pinned temperature-state-unknown seed tranche, not the
temperature-dependent N42 design card still required by bead `1sxe`.

The article is distributed under CC-BY-4.0; the manifest retains its authors,
title, DOI, and license identifier.

NGYC N42 reference:

- <https://doi.org/10.1038/s41598-023-47689-2>
- <https://www.nature.com/articles/s41598-023-47689-2>
- <https://creativecommons.org/licenses/by/4.0/>

## Jinshan N42 pristine temperature-endpoint tranche

`jinshan-n42-pristine-temperature/` retains the pristine commercial N42
sintered NdFeB state tested by Chang et al. The source pins the supplier to
Jinshan Magnetic Materials Co., identifies wire-cut `10 mm x 10 mm x 6 mm`
pieces, and names the NIM 6500C permanent-magnet measuring instrument. The
paper calls those unequal-sided pieces cubes; the pack preserves both the
printed wording and the dimensions rather than silently regularizing them.

Table 1 supplies intrinsic coercivity `Hcj`, remanence `Br`, and maximum energy
product at `25 degC` and `120 degC`. The six CGS values are normalized with
`1 G = 1e-4 T` and `1 Oe = 1000/(4*pi) A/m`; their exact temperature states
remain claim validity rather than being interpolated. Table 2 additionally
prints interval-average remanence and coercivity coefficients of
`-0.129 %/degC` and `-0.602 %/degC`, normalized to `-0.00129 K^-1` and
`-0.00602 K^-1`. The rounded Table 1 endpoints do not reproduce both printed
coefficients exactly, so the coefficients remain separate source claims and
are not recomputed from the rounded endpoint values.

The pristine `N42` specimen is distinct from the authors' `N42-HT` state: it
did not receive their `800 degC` for `5 h` plus `550 degC` for `2 h` treatment.
The source omits the supplier lot, composition, original sinter schedule,
coating, magnetization protocol, temperature-control details, field sweep,
repeats, dispersion, and uncertainty. Figure 10 depicts `25`, `80`, and
`120 degC` demagnetization curves, but the pack does not digitize graph pixels.
It admits no `80 degC` scalar, recoil permeability, irreversible-loss boundary,
continuous constitutive law, or extrapolation beyond `25-120 degC`. This card
therefore adds supplier-, state-, method-, and temperature-pinned endpoints;
it does not pretend to complete the bead's still-missing recoil-limit evidence.

The article is published under CC-BY-4.0; the manifest retains its authors,
title, DOI, publisher notice, and license identifier.

Jinshan N42 reference:

- <https://doi.org/10.1016/j.jmrt.2024.12.235>
- <https://www.sciencedirect.com/science/article/pii/S2238785424030382>
- <https://creativecommons.org/licenses/by/4.0/>

## Y30 catalog-midpoint and model-input tranche

`sjolund-2020-y30-catalog-model-inputs/` retains the Y30 inputs that Sjölund,
Leijon, and Eriksson used for their linear-generator finite-element study.
Table II reports the grade data as midpoint plus or minus half-range:
remanence `385 +/- 15 mT`, coercive field strength
`192.5 +/- 17.5 kA/m`, intrinsic coercive field strength
`200 +/- 20 kA/m`, and maximum energy product `28 +/- 2 kJ/m^3`.
The authors state that the average value of every magnitude is used in their
model, so the pack admits the four midpoints as the exact source-used values.

The printed half-ranges came from an online ferrite-grade catalog and carry no
statistical confidence semantics. The current runtime uncertainty types model
confidence-bearing half-widths, so the pack does not mislabel these catalog
bounds as a confidence interval. Instead, it preserves all four ranges in the
shared observation caveat and marks each runtime uncertainty `Unstated`.

The paper's Equation 2 derives a constant model relative permeability from the
remanence and energy-product midpoints,
`mu_r = Br^2 / (4 mu_0 BHmax)`. Exact evaluation gives
`1.0531580999909012`. That value is retained under the deliberately narrow
property name `model_relative_permeability`: it depends on an assumed linear
second-quadrant branch, and the authors assume isotropy because transverse data
are unavailable. It is not relabeled as measured recoil permeability or as a
minor-loop law. The source explicitly says that no adequate Y30
demagnetization curve was found to validate the modeled knee point.

The simulation is held at `20 degC`, but the catalog property's measurement
temperature is not reported. The pack therefore records `20 degC` only as
`source_simulation_temperature` while keeping
`source_magnetic_test_temperature_known = 0`. The catalog does not identify a
purchased product, product supplier, manufacturer, composition, production
lot, sinter schedule, magnetization state, or magnetic test method. Those gaps
remain fail-closed axes. This tranche supplies redistributable grade-level
model inputs; it does not complete bead `1sxe`'s supplier/process/temperature-
pinned Y30 material card or its temperature-dependent recoil-limit evidence.

The article is distributed under CC-BY-4.0; the manifest retains its authors,
title, DOI, journal citation, and license identifier.

Y30 catalog/model-input references:

- <https://doi.org/10.1063/1.5129303>
- <https://research.chalmers.se/publication/516199>
- <https://research.chalmers.se/publication/516199/file/516199_Fulltext.pdf>
- <https://creativecommons.org/licenses/by/4.0/>

## Y30 optimized-motor application-demagnetization tranche

`kim-baek-2026-y30-afcp-demagnetization/` retains two exact
temperature-conditioned outputs from Kim and Baek's optimized axial-flux
consequent-pole pump-motor model. At five times rated current and unchanged
`6000 rpm`, the modeled stator field reaches a spatial maximum of
`256.49 kA/m`. The reported maximum demagnetization coefficient is `1.654%`
at `20 degC` and `22.396%` at `-40 degC`.

The source defines that coefficient as `1 - Br1/Br0`, where `Br1` is the recoil
flux density recovered after the external field. Each claim is bound to the
optimized model's `750 W` output, `60.16 cm^3` magnet volume, five-times-rated
current, `100 Hz` rotational frequency, exact temperature, and maximum modeled
field. The two values are application-level 3-D FEA spatial maxima, not generic
Y30 material limits or a continuous temperature law.

Figure 13 plots the Y30 input B-H curves at both temperatures, but no points are
tabulated or digitized. The source does not identify the magnet manufacturer,
grade specification, chemistry, process, lot, magnetization state, or curve
provenance. It also omits FEA software, mesh-convergence evidence, complete
boundary-condition detail, numerical uncertainty, repeats, dispersion, and
prototype validation. The pack therefore admits no remanence or coercivity
endpoint, recoil permeability, reversible recoil line, standalone irreversible
loss boundary, extrapolation, or supplier-qualified Y30 card. This tranche is
useful overload-demagnetization evidence for one fully named motor model; it
does not complete bead `1sxe`'s supplier/process-pinned Y30 requirement.

The article is distributed under CC-BY-4.0; the manifest retains its authors,
title, journal, DOI, publisher notice, and license identifier.

Y30 optimized-motor application-demagnetization references:

- <https://doi.org/10.3390/app16021094>
- <https://www.mdpi.com/2076-3417/16/2/1094>
- <https://creativecommons.org/licenses/by/4.0/>

## NACA TN 2680 iso-octane flame-speed tranche

`naca-tn-2680-isooctane-flame-speed/` retains the
2,2,4-trimethylpentane fuel and Bunsen-flame observations reported by Dugger
and Graab. The source identifies the chemical name and a supplier-claimed
minimum purity of `99.6 mol%`, but not the supplier, lot, analytical method,
measured assay, or impurity composition. The pack therefore stores the purity
statement as a lower-bound claim with `Unstated` uncertainty; it does not turn
`99.6 mol%` into an exact pure-fluid composition.

The flame-speed observation transcribes all `15` experimental rows from Table
I. Each scalar retains the initial mixture temperature, oxygen mole fraction
in the oxygen-nitrogen stream, stream-flow Reynolds number, burner inside
diameter, and the experiment-wide average atmospheric pressure. Two repeated
nominal conditions at `422 K`, `0.210` oxygen fraction, `Re = 1000`, and
`1.256 cm` diameter remain distinct (`56.1` and `59.0 cm/s`); the two repeated
`422 K`, `0.294`, `Re = 900`, `0.617 cm` rows likewise remain distinct
(`108.0` and `102.1 cm/s`). No averaging or synthetic dispersion is invented.

The report locates each maximum only between equivalence ratios `1.0` and
`1.1`, so every flame-speed claim exposes an exact-ratio-missing flag and the
two source bounds. Its average `744 mm` mercury pressure is retained as a
source-precision `99.2 kPa` condition and explicitly marked as not known per
row. The `0.210` oxygen-fraction runs used service air with approximately
`0.4 vol%` water; that approximate value remains observation-only. The other
oxygen-nitrogen mixtures retain the supplier's `+/-0.1 vol%` analysis width.

The report also presents least-squares lines and an empirical maximum-speed
equation. Those are deliberately absent: the tranche records apparatus-bound
mixture observations, not a continuous constitutive law or a geometry-free
bulk-fluid flame-speed constant. NASA marks NTRS record `19930083861` public
and states that public use of the US Government work is permitted.

NACA TN 2680 references:

- <https://ntrs.nasa.gov/citations/19930083861>
- <https://ntrs.nasa.gov/api/citations/19930083861/downloads/19930083861.pdf>

## FACE G CDTRF-G 2023 v1 gasoline-surrogate tranche

`face-g-cdtrf-g-2023-v1/` retains Rao, Zheng, and Yang's five-component
`CDTRF-G` surrogate formulation targeting FACE G gasoline. The repository's
`2023-v1` suffix versions this exact source transcription; it does not imply
that the authors named a separate v1 formulation.

Table 2 states the composition on a volume-fraction basis: `23.75%`
isooctane, `19%` n-heptane, `42.75%` toluene, `9.5%` DIB
(diisobutylene), and `5%` cyclohexane. The five retained claims sum exactly to
`100%`. The source's surrounding FACE G hydrocarbon-class targets use a molar
basis, so the pack does not conflate those targets with the surrogate's volume
fractions. It also does not add supplier, lot, purity, preparation order,
mixing temperature or pressure, volume-contraction treatment, or a broader
commercial-DIB specification.

The article is internally inconsistent at its displayed RON precision. Table
2 prints `94` for CDTRF-G, while Table 7 prints `93.9` for Fuel 2 with the
identical five fractions. Both are described as calculated, not CFR-engine
measurements. The pack retains the two prints as separate `Unstated` claims
with separate observations; it does not average them, prefer one, or treat
their difference as statistical uncertainty.

This is the bead's first explicit versioned multicomponent gasoline-surrogate
composition, but not a general gasoline material card or reaction mechanism.
The article is distributed under CC-BY-4.0, and the manifest retains the
authors, title, DOI, journal citation, and license identifier.

FACE G CDTRF-G references:

- <https://doi.org/10.3390/molecules28114273>
- <https://www.mdpi.com/1420-3049/28/11/4273>
- <https://creativecommons.org/licenses/by/4.0/>

## WO 2018/125520 Formulation 8 5W-30 reference-oil tranche

`wo2018-125520-formulation-8-5w30/` retains the exact source identity of
Lewis, Fitch, and Han's Table IX comparator Formulation 8. Unlike a card named
only by viscosity grade, this one records all four as-added, source-era
commercial components on the table's stated mass basis: `60 wt%` SpectraSyn 4,
`10 wt%` Synesstic 5, `18 wt%` SpectraSyn Elite 150, and `12 wt%` Infineum
P6003. The fractions sum exactly to `100 wt%`, and the pack ID pins the patent
family, table, formulation number, and grade.

The same Table IX row reports kinematic viscosities of `61.49 mm^2/s` at
`40 degC` and `10.62 mm^2/s` at `100 degC`, viscosity-index scale reading
`164`, pour point `-66 degC`, CCS `4886 mPa*s` at `-30 degC`, MRV
`10782 mPa*s` at `-35 degC`, HTHS `3.395 mPa*s` at `150 degC`, and Noack
mass loss `9.2%` at `250 degC`. Dashes for CCS at `-25 degC` and MRV at
`-30 degC` remain absences; they are not converted to zeros.

This is a digital reference for that disclosed formulation, not a statement
that current products with the same trademarks are fungible. The source gives
no component lots, manufacturing dates, detailed package chemistry, final
blend protocol, sample storage history, test-method editions, repeats,
dispersion, or confidence metadata. Every runtime uncertainty therefore stays
`Unstated`, and the missing axes remain machine-readable. The patent source is
also not a license to practice patent claims or use trademarks. USPTO's Terms
of Use says patent text and drawings are typically not subject to copyright
restrictions; the manifest retains the inventors, title, family publication
identifiers, and that narrow redistribution decision.

WO 2018/125520 Formulation 8 references:

- <https://patents.google.com/patent/WO2018125520A1/en>
- <https://patents.google.com/patent/US20180179462A1/en>
- <https://www.uspto.gov/terms-use-uspto-websites#patent-information>

## NASA UAM motor-winding insulation constituent tranche

Three packs retain the constituent identity and actual process state from
Tallerico et al.'s NASA electric-aircraft winding-lifetime campaign. The source
is one coupled application record, but each bulk fact stays on the constituent
to which it belongs:

- `nasa-uam-mw16c-polyimide-magnet-wire/` identifies the campaign's `16 AWG`
  heavy-build `MW-16C` polyimide-insulated magnet wire under NEMA MW 1000-2003.
  Its `240 degC` and `20000 h` thermal-endurance statements are cross-bound as
  one ASTM D2307-2013 classification basis. Neither number is admitted as a
  free-standing lifetime law.
- `nasa-uam-nomex-410-slot-liner/` identifies DuPont Nomex 410 aramid paper and
  records the `0.08 mm` thickness selected between the M15 core and winding in
  the in-house motorettes. It is a stack selection, not a product-wide measured
  thickness distribution.
- `nasa-uam-cooltherm-ep2000-180c-cure/` identifies Parker LORD CoolTherm
  EP-2000 and distinguishes the highest completed `180 degC` post-cure step
  from the manufacturer-recommended `210 degC` final step that NASA explicitly
  omitted because of the nylon end-board temperature limit.

NTRS record `20240007451` marks the paper public, non-export-controlled, and
government public use permitted. The source does not publish wire vendor or
lot, constituent uncertainty, the wire-classification raw data, Nomex moisture
conditioning or thickness metrology, epoxy lot, complete cure times/ramps, or
degree of cure. Those gaps are retained as validity axes. Campaign PDIV values,
motorette hot-spot estimates, and modeled insulation stresses are deliberately
not bulk claims: they depend on the assembled geometry, winding, cooling,
electrical procedure, and thermal-cycle history.

NASA UAM winding-insulation reference:

- <https://ntrs.nasa.gov/citations/20240007451>
- <https://ntrs.nasa.gov/api/citations/20240007451/downloads/Winding%20Lifetime%20Testing%20V4.pdf>

## NASA-CR-195445 OMC rotary PS-200 coating-system tranche

`nasa-cr-195445-omc-ps200-rotary-coating/` is a configuration card for the
air-cooled OMC rotary-engine experiments in NASA-CR-195445, not a generic
aluminum, Wankel-housing, zirconia, or wear-coating card. The source explicitly
connects its PS-200 coating to US Patent 4,728,448. The pack therefore retains
the patent's preferred PS200 feedstock composition as three separately linked
mass fractions: `80 wt%` nickel-bonded chromium carbide, `10 wt%` silver, and
`10 wt%` barium-fluoride/calcium-fluoride eutectic. These are powder-feed
component fractions, not an assay of the as-sprayed coating.

The engine report applied PS-200 over a plasma-sprayed thermal-barrier
underlayer that it calls zirconia in the narrative and SX-331 elsewhere. The
pack preserves both names without asserting that they are chemically
identical. For Test 3's coated aluminum-alloy side-housing context, it retains
the reported RMS finish change from `21 microinch` (`0.5334 um`) to
`7 microinch` (`0.1778 um`) around the narrative's `2.5 h` run. For Test 6's
aluminum-alloy trochoid context, it retains the change from `24 microinch`
(`0.6096 um`) to `17 microinch` (`0.4318 um`) around the foreshortened `1.5 h`
run. Stage indices keep before and after claims distinct; they are not combined
into a wear rate.

The failure evidence is part of the identity. Test 3 developed a drive-side
TBC crack attributed to a probable spray defect and that housing was scrapped.
The Test 6 trochoid had uneven spray and honing, local PS-200 breakthrough,
leakage paths, and hot spots; the planned two-hour run stopped after `1.5 h`
when the exhaust-port temperature exceeded `500 degF`, near the report's stated
aluminum-alloy limit. Aluminum grade/temper/casting route, powder and substrate
lots, SX-331 chemistry, coating thickness, as-tested profile, spray parameters,
surface-finish method, metrology settings, repeats, and uncertainty remain
explicitly unavailable.

Rotary-coating references:

- <https://ntrs.nasa.gov/citations/19950016805>
- <https://ntrs.nasa.gov/api/citations/19950016805/downloads/19950016805.pdf>
- <https://patents.google.com/patent/US4728448A/en>

## NASA AISI 52100 dry and GXL-320A interface tranches

`nasa-52100-dry-air-interface/` is an ordered rider-on-disk system card for the
SAE 52100 atmospheric endpoint printed in Figure 11(a) of Buckley, Swikert, and
Johnson's NASA-authored friction study. It retains the source's Table 1
composition range and Rockwell C60..63 condition, separate 3/16-inch-radius
rider and 2-inch-diameter disk roles, finish-ground 4..8-microinch surfaces,
specified cleaning sequence, fresh specimens, one-hour duration, `75 degF`,
`760 mm Hg`, `390 ft/min`, and `1000 gram-force`. The sole claim is the printed
kinetic-friction coefficient `0.45`. SI validity points record the conventional
conversions explicitly; the source does not name the roughness metric, repeat
count, dispersion, humidity, oxide state, wear rate, or transferable friction
law, and no such claim is inferred.

`nasa-52100-gxl320a-vacuum-interface/` is a separate ordered four-ball system
card from NASA/TM-1999-209064. It binds one rotating and three stationary
grade-25, 9.5-mm AISI 52100 balls to the reported cleaning sequence, named
GXL-320A base grease, vacuum strictly below `6.7e-4 Pa`, `200 N`, `28.8 mm/s`,
initial `3.5 GPa` Hertz mean stress, approximately `23 degC`, four-hour test,
and four-run aggregate. Table 2's `0.11` mean and observed `0.09..0.23` range
are three distinct claims with `Unstated` uncertainty; the observed extrema
are not represented as a confidence interval. Supplier lot, thickener, full
formulation, steel heat/hardness/roughness, grease quantity, run-level values,
wear rate, friction law, service lifetime, and transfer to another grease or
environment remain outside the card.

NTRS marks both records public, U.S.-Government works, and Public Use
Permitted. The 1999 PDF's preprint page separately requests author permission
to cite or reproduce. Its manifest preserves that conflict and the pack uses
only attributed factual rows; it does not reproduce a table or figure.

Interface references:

- <https://ntrs.nasa.gov/citations/20150020881>
- <https://ntrs.nasa.gov/api/citations/20150020881/downloads/20150020881.pdf>
- <https://ntrs.nasa.gov/citations/19990031946>
- <https://ntrs.nasa.gov/api/citations/19990031946/downloads/19990031946.pdf>

## NASA SAE 4340/high-lead-bronze hexane journal tranche

`nasa-tn-d-2223-4340-high-lead-bronze-journal/` binds NASA-TN-D-2223's
ordered Rockwell C35 SAE 4340 journal and high-lead bearing bronze
(`70 wt% Cu`, `26 wt% Pb`, `4 wt% Sn`) to the report's pressure-fed hexane
screening rig. The journal was finished to `5..10 microinch rms`; the nominal
`1.5 in`-ID by `1.5 in`-long bearing had `0.0015 in` room-temperature
diametral clearance and the source's Figure 2(a) supply groove. The pump inlet
was `10 psig` through a `0.125 in` hole. Hexane viscosity is retained as
`0.043e-6 reyn` at the source's `75 degF` reference point.

The sole runtime claim is Table I's maximum demonstrated unit bearing load of
`220 psi` (`1516846.6044970395 Pa`) for the run at `15000 rpm` and `107 ft/s`.
The source says the maximum speed/load was held for `30 min`, the total test
time was `90 min`, the test objective was attained, and the bearing surface was
excellent with light polishing in the load zone. This is a source- and
apparatus-bounded demonstrated result, not a design allowable. Exact run
temperature, surrounding atmosphere, hexane grade/purity, bearing roughness,
material heats, repeats, dispersion, endpoint torque, quantitative wear, and
confidence remain unstated. The pack therefore supplies no friction
coefficient, wear rate, lifetime, or transferable journal-bearing law.

Journal-interface references:

- <https://ntrs.nasa.gov/citations/19650001032>
- <https://ntrs.nasa.gov/api/citations/19650001032/downloads/19650001032.pdf>

## Carbon-filled PTFE/chrome-coated piston-rod seal tranche

`zhang-2021-carbon-ptfe-cr-piston-rod/` binds one GMSS Sterling stepped seal
to the chrome-coated counterface and named lubricant tested by Zhang et al.
The seal's wear ring is PTFE filled with approximately `15%` carbon fiber,
mixed at high speed, molded at ambient temperature, sintered, and turned; its
supporting D-ring is nitrile rubber. The counterface is a stainless-steel
piston rod carrying approximately `50 micrometers` of electroplated chromium.
Kunlun 15 aviation hydraulic oil is the retained medium.

The four runtime claims are the paper's printed single-seal instroke and
outstroke friction forces at its observed extrema: `97.8955 N` and `75.0238 N`
at `10 MPa`, `0.5 m/s`, and `516.9906 N` and `404.4382 N` at `35 MPa`,
`0.1 m/s`. All points retain the `20 degC` experiment temperature, `300`
reciprocations, `2.4 kHz` sensor sampling, two seals per experiment, and the
reported two-experiment average. The paper does not publish the two replicate
values or pointwise dispersion, so all four uncertainties remain `Unstated`.

This is a coated piston-rod card, not coated-bore evidence. It supplies no
coated-bore force, friction coefficient, leakage rate, quantitative wear rate,
seal life, or transferable friction law. The source observed seal-lip wear and
Cr/Fe material transfer; the card therefore does not turn the printed endpoint
forces into durability authority. Exact resin, fiber, rubber, steel, coating,
surface-finish, and lubricant lots; rod roughness; stroke; preload/contact
pressure; and full oil formulation remain unavailable.

Seal-interface reference:

- <https://doi.org/10.3390/app11041370>
- <https://creativecommons.org/licenses/by/4.0/>

## PTFE CF10/steel-cylinder-bore/ISO VG 32 interface tranche

`deaconescu-2020-ptfe-cf10-steel-bore-interface/` binds the source's coaxial
hydraulic-piston seal to the cylinder's inner steel surface. Surface A is
`90 wt%` PTFE plus `10 wt%` carbon fiber, backed by a `70 Shore A` NBR O-ring
mounted at `15%` initial radial deformation. Surface B is the source-described
steel cylinder bore. The medium is the source-linked CITGO Anti-Wear Hydraulic
Oil Grade 32 reference, while its production lot and additive formulation
remain unknown.

The seven runtime claims transcribe the PTFE CF10 rows of Table 2. They retain
the exact pressure, oil temperature, piston speed, displayed voltage drop, and
Equation 11 fluid-film thickness for two velocity-sweep rows, four
pressure-sweep rows, and the sole populated temperature-sweep row. Thickness
ranges from `2.5 micrometers` at `200 bar`, `333 K`, `0.2 m/s` to
`10.8 micrometers` at `100 bar`, `333 K`, `0.5 m/s`. Table 2 prints
`8.3 micrometers` twice and `8.1 micrometers` once at the same nominal
`100 bar`, `333 K`, `0.2 m/s` condition. All three source rows remain separate;
they are not averaged or converted into an invented uncertainty interval.

The film thicknesses are derived by the paper's resistive method from the
displayed voltage using `Ra = 10 megaohm` and
`CR = 63.735 megaohm/micrometer`. The paper's `3.5%` maximum error compares
experimental thickness with its theoretical Equation 6 result; it is not a
measurement-confidence statement, so runtime uncertainty remains `Unstated`.
The source omits repeats, raw traces, dispersion, calibration uncertainty,
exact PTFE constituent grades, supplier/lot/process, bore steel grade, measured
finish, coating state, oil lot/formulation, opposed-chamber state, leakage,
friction force, wear, and life.

This is genuine piston-seal/cylinder-bore evidence, but it is not a coated-bore
card: the source never identifies a bore coating. Recommended catalog
roughness elsewhere in the paper is not metrology of the tested bore. The pack
therefore refuses coated-bore film thickness, friction coefficient, leakage,
wear, service life, and a transferable lubrication law; bead `1sxe` still owns
the exact coated-bore gap.

PTFE CF10/steel-bore references:

- <https://doi.org/10.3390/polym12010155>
- <https://www.mdpi.com/2073-4360/12/1/155>
- <https://creativecommons.org/licenses/by/4.0/>

## MAHLE PTFE/nickel-SiC-coated-bore patent tranche

`mahle-wo2019072721-ptfe-nickel-sic-coated-bore-interface/` records an actual
coated-cylinder/PTFE-ring pairing disclosed by MAHLE. Surface A is the patent's
PTFE piston ring. Surface B is the cylinder running surface with a liner layer
whose main component is nickel and which contains silicon carbide. Claims 9 and
10 specify that surface as honed to an approximate arithmetic mean roughness of
`0.1 micrometer`.

The runtime claims are design-disclosure facts, not test results. They retain
the claim-defined approximate SiC window of `2..6 wt%`, the strict
`1 < particle diameter < 4 micrometers` window and preferred
`2.5 micrometer` particle diameter, and the strict
`40 < liner thickness < 100 micrometers` window and preferred
`50 micrometer` thickness. The disclosed embodiment separately says
`2..4 wt%` SiC. That narrower statement and claim 5's wider `2..6 wt%`
statement remain distinct rather than being fused into a fabricated assay.

This card establishes a named PTFE/coated-bore design interface, but its
boundary is deliberately narrow. The patent gives no manufactured specimen,
PTFE grade or filler, cylinder substrate, nickel-alloy chemistry, coating
application process, as-built metrology, working-gas state, operating point,
uncertainty, or quantitative friction, wear, leakage, or lifetime result. It is
also explicitly oil-free: no lubricant is present or identified. Consequently
the card does not complete bead `1sxe`'s lubricant-qualified coated-bore
requirement and authorizes no transfer of the patent's qualitative low-friction
or durability language into a constitutive law.

The German Patent and Trade Mark Office states that officially published patent
documents are excluded from copyright protection under Section 5(2) UrhG,
subject to source acknowledgement and the prohibition on altering the patent
document. The manifest retains inventor, assignee, title, and the DE and WO
family identifiers. That narrow redistribution basis is not a license to
practice patent claims, use trademarks, qualify a product, or claim performance.

MAHLE PTFE/nickel-SiC-coated-bore references:

- <https://patents.google.com/patent/WO2019072721A1/en>
- <https://patents.google.com/patent/DE102017218156A1/en>
- <https://www.dpma.de/english/services/customer_care_services/faq/index.html>
- <https://www.gesetze-im-internet.de/urhg/__5.html>

## A2017/Seiken long-life-coolant/air wetting tranches

`yilmaz-2026-a2017-seiken-llc-ra005-wetting/` and
`yilmaz-2026-a2017-seiken-llc-ra3-wetting/` are separate
solid-liquid-gas system cards for the two texture states printed by Yilmaz et
al. Surface A is the source-composition A2017 aluminum specimen; surface B is
the source's Seiken Chemical Industry long-life-coolant free surface; the same
liquid is the intervening medium and ambient air is the environment. The
mirror-finish `Ra = 0.05 micrometer` card retains `62.70 deg`
(`1.0943214410004447 rad`). The EDM `Ra = 3 micrometer` card retains
`102.89 deg` (`1.7957692673769656 rad`). Separate packs prevent roughness from
being demoted from interface identity to an optional query field.

Both cards pin the pre-boiling B100W contact-angle-meter measurement, theta/2
method, approximately `1 microliter` droplet, five-measurement average,
isopropyl-alcohol/nonwoven cleaning, `24.5 degC`, and `60%` relative humidity.
Table 2 reports calculated contact-angle measurement uncertainty of
`+/-0.01 deg` but gives no confidence semantics, and the five measurements and
their dispersion are absent. The runtime uncertainty therefore remains
`Unstated`, not a fabricated confidence band.

The source names the coolant manufacturer but omits product code,
concentration, formulation/additive package, lot, and age; it also omits exact
A2017 temper/lot, native-oxide state beyond the cleaning procedure, air
pressure/composition, and cleaning-to-measurement elapsed time. The cards claim
only pre-boiling static angle. They supply no advancing/receding angle,
hysteresis, surface energy, post-boiling/corroded-state angle, temperature law,
or transferable wetting law.

The article is CC BY-NC-ND 4.0. Its manifest retains the noncommercial and
no-derivatives limits verbatim. These packs contain attributed factual scalar
rows and reproduce no article figure; commercial or adapted reuse still needs
independent rights review.

Coolant-wetting references:

- <https://doi.org/10.1080/08916152.2026.2615768>
- <https://creativecommons.org/licenses/by-nc-nd/4.0/>

To compile the sources into canonical runtime packs:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/methane/manifest.tsv \
  --out /path/to/CH4.fsspcpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nist-srm-1720-northern-continental-air/manifest.tsv \
  --out /path/to/nist-srm-1720-northern-continental-air.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nist-srm-2728-auto-emission-reference-gas/manifest.tsv \
  --out /path/to/nist-srm-2728-auto-emission-reference-gas.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aluminum-6061-t6-cryogenic/manifest.tsv \
  --out /path/to/aluminum-6061-t6-cryogenic.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/ofhc-copper-rrr100/manifest.tsv \
  --out /path/to/ofhc-copper-cryogenic.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/ptfe-teflon-nist-cryogenic/manifest.tsv \
  --out /path/to/ptfe-teflon-nist-cryogenic.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/peek-nasa-thermic-plate/manifest.tsv \
  --out /path/to/peek-nasa-thermic-plate-2021.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aisi-4140-rc33/manifest.tsv \
  --out /path/to/aisi-4140-qq-s-624-rc33.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aisi-1045-cold-drawn/manifest.tsv \
  --out /path/to/aisi-1045-cold-drawn-tensile.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aisi-52100-cvm-hot-hardness/manifest.tsv \
  --out /path/to/aisi-52100-cvm-nasa-tn-d-6632.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aisi-9310-cvm-carburized/manifest.tsv \
  --out /path/to/aisi-9310-cvm-carburized-nasa-tm-104352.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/napc-pe-5-l-1274-gear-oil/manifest.tsv \
  --out /path/to/napc-pe-5-l-1274-polyol-ester-gear-oil.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/napc-pe-5-l-1307-1553-gear-oil/manifest.tsv \
  --out /path/to/napc-pe-5-l-1307-1553-mil-l-23699-gear-oil.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/rheolube-2000-pennzane-grease/manifest.tsv \
  --out /path/to/rheolube-2000-pennzane-shf-x-2000-bearing-grease.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/pennzane-shf-x-2000-bearing-oil/manifest.tsv \
  --out /path/to/pennzane-shf-x-2000-mac-aerospace-bearing-oil.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/gray-cast-iron-s2-s/manifest.tsv \
  --out /path/to/pearlitic-gray-cast-iron-s2-s-sr-fesi.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-cr-115153-water-ethylene-glycol/manifest.tsv \
  --out /path/to/nasa-cr-115153-inhibited-water-ethylene-glycol-coolant.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/n0602-001-nitrile-jp8-compatibility/manifest.tsv \
  --out /path/to/n0602-001-nitrile-o-ring-jp8-compatibility.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/torrent-2018-m19-steinmetz-inputs/manifest.tsv \
  --out /path/to/torrent-2018-m19-steinmetz-inputs.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/ngyc-n42-sintered-nickel-coated/manifest.tsv \
  --out /path/to/ngyc-n42-sintered-ndfeb-nickel-coated-cubes.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/naca-tn-2680-isooctane-flame-speed/manifest.tsv \
  --out /path/to/naca-tn-2680-2-2-4-trimethylpentane-flame-speed.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/face-g-cdtrf-g-2023-v1/manifest.tsv \
  --out /path/to/face-g-cdtrf-g-2023-v1.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/wo2018-125520-formulation-8-5w30/manifest.tsv \
  --out /path/to/wo2018-125520-formulation-8-5w30.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-uam-mw16c-polyimide-magnet-wire/manifest.tsv \
  --out /path/to/nasa-uam-mw16c-polyimide-magnet-wire.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-uam-nomex-410-slot-liner/manifest.tsv \
  --out /path/to/nasa-uam-nomex-410-slot-liner.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-uam-cooltherm-ep2000-180c-cure/manifest.tsv \
  --out /path/to/nasa-uam-cooltherm-ep2000-180c-cure.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-cr-195445-omc-ps200-rotary-coating/manifest.tsv \
  --out /path/to/nasa-cr-195445-omc-ps200-rotary-coating-system.fsmatpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-52100-dry-air-interface/manifest.tsv \
  --out /path/to/nasa-52100-dry-air-interface.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-52100-gxl320a-vacuum-interface/manifest.tsv \
  --out /path/to/nasa-52100-gxl320a-vacuum-interface.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nasa-tn-d-2223-4340-high-lead-bronze-journal/manifest.tsv \
  --out /path/to/nasa-tn-d-2223-4340-high-lead-bronze-journal.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/zhang-2021-carbon-ptfe-cr-piston-rod/manifest.tsv \
  --out /path/to/zhang-2021-carbon-ptfe-cr-piston-rod.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra005-wetting/manifest.tsv \
  --out /path/to/yilmaz-2026-a2017-seiken-llc-ra005-wetting.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/yilmaz-2026-a2017-seiken-llc-ra3-wetting/manifest.tsv \
  --out /path/to/yilmaz-2026-a2017-seiken-llc-ra3-wetting.fsintpk
```

## No-claim boundary

The species-association gas tranches are identity/standard-state associations,
not complete material cards. By themselves they do not specify an air or
exhaust mixture, composition basis, humidity level, combustion state, or
validity domain. The separate NIST SRM 2728 pack does define a synthetic auto-
emission calibration mixture at its three nominal title values, but it is not
engine-generated exhaust and does not fill the certificate's blank cylinder-
specific certified values. In particular, the H2O association does not select
a wet-air composition and the CO association does not claim incomplete
combustion. None of these packs supplies heat-capacity coefficients, an
equation evaluator, thermodynamic-property uncertainty, a reaction mechanism,
an equilibrium result, or transport data. A decimal agreement check is not an
uncertainty estimate.

The Aluminum 6061-T6 tranche is not a general-purpose design card. It contains
six polynomial-derived scalar points for one named temper, not density,
composition tolerances, anisotropy, batch/process variation, plasticity,
strength, fatigue, fracture, corrosion, joining, or a continuous-temperature
constitutive model. NIST's curve-fit error is retained but not represented as a
confidence interval. The NASA comparison uses nearby temperatures and is only
a coarse independent agreement check. These missing claims and the remaining
named materials and interface systems keep bead
`frankensim-ext-matdb-seed-dataset-1sxe` open.

The two AISI 52100 interface tranches are apparatus- and history-bound system
cards, not bulk coefficients for bearing steel. The dry value applies only to
the retained rider/disk roles and atmospheric endpoint. The GXL-320A values
apply only to the retained four-ball vacuum system; its printed extrema are an
observed range, not uncertainty bounds. Neither pack supplies static friction,
wear, a constitutive friction law, lifetime, air/vacuum transfer, surface-state
evolution, or authority for another heat, finish, cleaner, or lubricant.

The SAE 4340/high-lead-bronze tranche records one successful NASA
pressure-fed-hexane screening endpoint. `maximum-demonstrated` does not mean
safe working load, fatigue capacity, hydrodynamic closure, or a material-pair
constant. No coefficient of friction, torque at that endpoint, wear rate,
thermal state, lifetime, or transfer across clearance, finish, groove, fluid,
load history, or material heat is claimed.

The OFHC Copper tranche likewise does not select between UNS C10100 and C10200,
define electrical resistivity, density, mechanical strength, temper, grain
state, purity beyond the source's OFHC class, or a continuous-temperature
model. Only its conductivity claims bind `RRR = 100`; the source leaves RRR
unstated for specific heat, and the stored observation preserves that gap.

The PTFE/Teflon tranche is not a generic PTFE design card or a seal,
tribology, dielectric, radiation, outgassing, chemical-compatibility, creep,
strength, or lifetime model. It binds only NIST's source label and four
polynomial-derived exact-temperature points; it does not select a resin grade,
crystallinity, filler, processing state, supplier, or density. The published
fit errors are not confidence intervals, and no interpolation, extrapolation,
Young's-modulus, or thermal-expansion-coefficient claim is admitted.

The PEEK THERMIC tranche is likewise not a generic PEEK card, seal design, or
service-temperature qualification. It binds one NASA test plate and its thermal
model inputs without supplying grade, supplier, crystallinity, filler,
processing history, thermal-expansion behavior, mechanical properties,
tribology, permeability, chemical compatibility, radiation response, creep,
aging, or lifetime authority. The fit diagnostics are not statistical
uncertainty, and the density's unknown test temperature must be acknowledged.

The AISI 4140 tranche applies only to the exact NASA heat/process/hardness/bar
condition above and the two retained temperatures. It does not provide a
continuous-temperature law, elastic modulus or Poisson ratio, a constitutive
plasticity model, fatigue or fracture curves, general cryogenic qualification,
or permission to substitute another 4140 heat treatment. The report's Rockwell
C44 branch remains a separate conflicting condition and is not fused here.

The AISI 1045 tranche applies only to the source's cold-drawn bar and three
ASTM E8 specimens at the retained test rate. The source omits test temperature,
supplier, heat number, and any heat-treatment history beyond `cold-drawn`;
callers must explicitly acknowledge the missing temperature metadata through
the fail-closed validity flag. The derived intervals assume three iid normal
replicates and are not population scatter, minimum design values, process
capability, or permission to substitute a different AISI 1045 condition. The
tranche makes no hardness claim and admits no joint covariance.

The AISI 52100 tranche is not a generic bearing-steel allowable. It applies to
one CVM ingot, the reported chemistry, the common austenitize/quench/first-
temper spine, and the five separately keyed second-temper states. Rockwell C is
stored as a named empirical scale reading with dimensionless storage, not as a
ratio quantity. The tranche does not provide a hot-hardness curve, long-time
tempering stability, rolling-contact fatigue life, wear, elastic/plastic
constitutive behavior, cleanliness-based life adjustment, or permission to
apply the short-term NASA measurements to another heat. The `<2%` retained-
austenite result remains censored, and all table values retain `Unstated`
measurement uncertainty.

The AISI 9310 tranche applies only to the single CVM heat, spur-gear lot,
nominal chemistry, and complete carburize/quench/subzero/temper/grind/stress-
relief schedule reported by NASA-TM-104352. It does not provide a generic
grade allowable, actual heat chemistry, hardness dispersion, a case-hardness
profile, residual stress, microstructure, elastic/plastic constitutive law,
fracture or fatigue properties, or permission to transfer the gear-system
results to another geometry or process. The report's lubricant-dependent
surface-fatigue lives are system properties and are deliberately excluded from
the bulk-material pack. Its internally conflicting C58 and C60 case statements
remain separate `Unstated` claims; neither is promoted as the resolved value.

The NASA/NAPC oil tranches are not complete lubricant formulations or generic
MIL-L-23699 cards. Proprietary additive and chemical identities remain
unknown. The packs provide no viscosity claim or curve because Table IV omits
the viscosity unit; no density, pressure-viscosity law, oxidation stability,
wear coefficient, friction coefficient, compatibility, aging, or service-life
model is inferred. The report's EHD film thickness, lambda ratio, and
surface-pitting lives belong to the tested oil/gear/roughness/load/temperature
system and are deliberately excluded from these bulk-oil packs. The two B
batch columns remain distinct where the source reports distinct values.

The Rheolube 2000 tranche is not a universal bearing-grease card or a complete
formulation. It binds the public product name, Pennzane SHF X-2000 base-oil
association, thickener identity, and the three admitted typical values, but
does not claim exact thickener fraction, dropping point, penetration,
volatility, wear, oxidation life, torque, friction, EHD film, bearing life,
vacuum-hardening state, contamination state, compatibility, or service
temperature range. NLGI grade is retained as an ordinal named-scale reading.
Application admission still requires the actual bearing, load, speed,
temperature, atmosphere, preparation, and life-test evidence.

The Pennzane SHF X-2000 tranche is not a complete formulated bearing oil or a
generic aerospace-lubricant allowable. It binds one named MAC base fluid and
seven typical source values, but supplies no lot identity, additive package,
continuous viscosity law, pressure-viscosity behavior, volatility claim,
thermal-expansion claim, wear/friction law, oxidation life, material
compatibility, contamination state, or bearing-life authority. Flash and pour
points are characteristic test temperatures, not a certified operating range.
Application admission still requires the actual bearing system and its
environmental/life evidence.

The S2-S gray-iron tranche applies only to the reported charge, composition,
Sr-FeSi inoculation, mould, and fully pearlitic/type-A-graphite state. Its
graph-digitized strength and conductivity are source-presentation estimates,
not design minima, grade limits, or certified allowables. It supplies no
elastic/plastic constitutive model, compression or fatigue law, elevated-
temperature transport curve, thermal-expansion law, casting-section-size
transfer, wear law, or universal engine-block/housing identity. Exact test
temperature, statistical confidence, raw replicates, and joint covariance are
not claimed.

The NASA-CR-115153 coolant tranche is not a generic `50/50` coolant card,
freeze/boil envelope, corrosion-inhibitor qualification, wetting system,
electrical-conductivity model, compatibility result, degradation law, or
service-life authority. It binds only Table 5's named water/ethylene-glycol,
sodium-nitrite, and sodium-benzoate specification. The exact within-range
composition, inferred glycol balance, supplier, grades, preparation history,
test pressure, and most property temperatures remain unknown. Its approximate
heat-capacity equation is represented only by three exact-point transcriptions;
no interpolation or extrapolation is authorized. The comparison report has a
different, composition-basis-unspecified fluid and remains evidence only.

The N0602-001 tranche is not a generic nitrile seal card, formulation,
hardness specification, constitutive law, permeability model, compression-set
law, ozone/aging qualification, fuel-system compatibility approval, or service
life. Its TGA, absorbed-fuel, partitioning, and swell-regression claims bind the
source's exact O-ring code and JP-8/FT test matrix. Missing exposure conditions
must be acknowledged. The two printed swelling slopes remain conflicts, the
approximate prediction-interval overlap is not promoted to an exact claim, and
the regression intercept does not certify shrinkage in a particular fuel.

The NGYC N42 tranche is not a generic N42 allowable, full magnet identity,
demagnetization model, recoil law, thermal derating curve, irreversible-loss
boundary, or service-temperature qualification. It binds only the supplier,
sintered grade, nickel coating, cube family, and four conflict-preserving claims
reported in the CC-BY-4.0 article. The source omits the property test method and
temperature, production lot, chemistry, exact process, intrinsic coercivity,
recoil permeability, temperature coefficients, and second-quadrant curves.
The `318.3 kJ/m^3` and `42 MGOe` energy-product representations remain separate;
neither authorizes selecting a resolved value for motor design.

The NACA TN 2680 tranche is not an exact pure-iso-octane fluid card, gasoline
surrogate specification, reaction mechanism, ignition model, or continuous
flame-speed law. It binds one supplier-minimum-purity statement and fifteen
atmospheric Bunsen-flame observations to their reported mixture and apparatus
conditions. Supplier, lot, exact assay, impurity composition, row-level
pressure, exact maximizing equivalence ratio, raw images, dispersion, and
confidence metadata remain absent. It supplies no density, viscosity, surface
tension, heat capacity, latent heat, vapor pressure, octane rating, transport
curve, or permission to transfer the measured maxima to another burner,
pressure, oxidizer, humidity state, or practical gasoline. The empirical fit
printed by the report remains excluded as model-only evidence.

The FACE G CDTRF-G 2023 v1 tranche is not an assay of FACE G, a fungible
gasoline recipe, a set of pure-component cards, or a combustion mechanism. It
binds five published volume fractions and two conflict-preserving calculated
RON prints to one named surrogate. Component supplier, lot, purity, detailed
DIB component specification, preparation history, mixing conditions,
contraction behavior, measurement uncertainty, CFR-engine RON, and target-fuel
batch identity remain absent. It supplies no density, viscosity, vapor pressure,
distillation curve, heat capacity, latent heat, flame speed, ignition delay,
emissions behavior, storage stability, material compatibility, or authority
outside the source's reported modeling context.

The NASA UAM insulation tranches do not constitute a qualified winding system,
generic constituent allowables, or a lifetime model. The MW-16C temperature
and duration are one thermal-class basis and may not be separated or
extrapolated. The Nomex claim is only the source-selected liner thickness. The
EP-2000 claims preserve an intentionally incomplete cure state; the omitted
`210 degC` step is not evidence of a completed cure. No pack claims dielectric
strength, thermal conductivity, moisture response, fatigue, adhesion, partial
discharge inception voltage, Arrhenius activation energy, degree of cure, or
service life. The paper's PDIV measurements and FEA stress/temperature results
belong to the assembled motorette and its geometry, cooling, excitation, and
aging history, so they are excluded from these bulk constituent cards.

The NASA-CR-195445 coating tranche is not a transferable PS-200 property card,
an aluminum-alloy allowable, or a qualified Wankel coating stack. Its
composition values describe named feedstock components, not post-spray phases.
Its four roughness values belong only to the reported OMC engine configurations
and their imperfect surfaces. No friction coefficient, wear rate, adhesion,
porosity, hardness, thermal conductivity, coating-thickness tolerance,
temperature envelope, fatigue life, specific-fuel-consumption benefit, or
service life is admitted. The cracked Test 3 housing and locally breached,
leaking Test 6 trochoid prevent any universal durability claim. Redistribution
of source text does not grant patent-practice, trademark, or application-
qualification rights.

## Lip and perioral tranche, October 2026 (pexfp.194, M8-W2d, lane P5)

`lip-resonance-human-video-newton-thesis2009`, `lip-lower-aspiration-luboz-abme2014` and
`perioral-stiffness-chu-barlow-jslhr2009` hold the only measured mechanical numbers this
program has found for the human lip, and every one of them is an **apparatus response**, not a
constitutive property of a tissue layer. Property names end in `_apparatus_effective` and
observation ids start with `apparatus-`, so none resolves as an intrinsic modulus.
Newton's thesis (Sec. 4.4, Table 4.2) gives the resonance frequencies and Q factors of two
trombonists' whole embouchures, measured by frontal video under an acoustic sweep: an outward-
striking resonance at 32 to 60 Hz with Q 1.5 to 2.4 and a second resonance, for target notes
Bb1 to F3; the printed second-resonance range in Table 4.2 (70 to 100 Hz) disagrees with the
text of the same section (140 to 200 Hz) and both are kept. Luboz et al. give an initial
Young's modulus of the lower lip, 33.7 +- 7.3 kPa, from an aspiration device and a
Neo-Hookean inverse analysis on 16 young subjects (abstract values; the statistic behind the
+- is not named there). Chu, Barlow and colleagues give the lateral stiffness of relaxed mouth
corners, 0.073 N/mm on average and a quadratic regression on a 24 mm stretch, with the unit of
the regression variable inferred, not printed. The Takemoto upper-lip pack is the
teeth-backed indentation it always was; its observation id and property names now carry the same
apparatus form (`apparatus-` prefix, `_apparatus_effective` suffix), and the cheek values with them. No layered lip constitutive law, no loss tangent, no
orbicularis oris active stress and no lip mass or density has a measured value here: they are
UNKNOWN, and the ledger rows of this tranche record what was read and what was not (one lip-tissue
source sits behind a paywall and is a `pending-permission` row).

## Human vocal-fold tranche, October 2026 (pexfp.194, M8-W2d, lane P1)

Four packs hold the measured mechanical numbers this program has found for the human vocal fold, typed
from publisher-copyright papers read as text (numeric facts with attribution, licence token
`Published-Scientific-Facts-Attribution-Retained`, no table layout reproduced).
`vocalfold-human-layers-kelleher-bmmb-pmc3745772` has the cover and the contralateral ligament of six
subjects measured separately (longitudinal modulus from cyclic tension with optical strain; transverse
modulus and longitudinal shear modulus from an inverse fit of an anisotropic contact model, which are
strongly coupled: a 10 percent change of one moves the other by up to 30 percent). `vocalfold-human-elastic-set-alipour-jvoice-pmc3508138`
has true and false fold moduli, a Poisson ratio and the first measured human fold density (Archimedes,
N = 10). `vocalfold-human-cover-shear-chan-rodriguez-jasa-pmc2561277` has the shear rheology of seven
covers at 1 to 250 Hz as power laws fitted to the mean, with the G'/G'' crossover near 50 Hz and a loss
factor G''/G' of 0.74 at 100 to 250 Hz. `vocalfold-human-indentation-layers-chhetri-jvoice-pmc2891297`
has indentation moduli of the intact hemilarynx, the isolated cover and the isolated thyroarytenoid
muscle; it is an **apparatus response** (properties end `_apparatus_effective`, observations start
`apparatus-`), each row carrying its indenter regime. Layer thickness is geometry and is not stored in
any of them. Every number comes from 1 Hz cyclic or quasi-static loading except the cover rheology; donors are mostly
geriatric and cadaveric; test temperature is unstated in three of the four. No constitutive law family,
no active stress of the muscle, no tongue and no fingertip value is admitted here: they are UNKNOWN. The
fold model's per-mass damping ratio of 0.1 in `glottis.rs` has no counterpart in these packs (the cover
loss factor is a different quantity).

## Clear-wood thermal completion, October 2026

The existing `pine-loblolly`, `pine-longleaf`, `pine-slash`, `pine-ponderosa`,
`oak-northern-red`, `oak-white`, `cedar-northern-white`, and `fir-subalpine`
FPL-GTR-282 bundles now also carry `specific-heat-capacity` and
`thermal-conductivity-transverse`. These are the same handbook's engineering
correlations from [Chapter 4, equations (4-17)–(4-20)](https://research.fs.usda.gov/download/treesearch/62243.pdf),
evaluated at fixed 12% dry-basis moisture, with explicit
`engineering-reference-model=1`, `clear-specimen=1` and
`source-pressure-known=0` conditions. Cp is affine in absolute temperature
over 280–420 K. Conductivity uses each bundle's own specific gravity and is
admitted only at 297.15 K across the grain; it is neither an isotropic value
nor a temperature-dependent conductivity law. The reported possible 20%
conductivity variation is retained as a model caveat, not fabricated confidence
metadata. These additions do not extend the original mechanical data's unknown
test temperatures, supply a thermal-density law, or qualify graded lumber.

## Clear-wood moisture response

All 33 existing FPL clear-wood bundles carry radial, tangential and independent
volumetric green-to-ovendry shrinkage, plus the explicit 30% fiber-saturation
approximation from [FPL-GTR-282, Chapter 4, Tables 4-3/4-4 and Equation (4-9)](https://research.fs.usda.gov/download/treesearch/62243.pdf).
Baldcypress shrinkage instead comes from [FPL-GTR-102, p. 115](https://research.fs.usda.gov/download/treesearch/5648.pdf);
imported-wood records retain the source's warning that literature values may not
represent true species averages. `fs_material::wood_moisture::ResolvedWoodMoisture`
consumes the exact source claims and 12% density reference to calculate radial
and tangential dimension ratios, volume ratios and density at prescribed uniform
dry-basis moisture within the retained model domain. The five imported woods
use explicitly derived reference densities: `rho12 = 1120*Ggreen/(1-0.6*Sv)`,
where `Sv` is fractional volumetric shrinkage, following Equations (4-9), (4-11)
and (4-14). These references combine green-volume gravity with the separate
volume model; they are not new density measurements. The consumer preserves
source receipts and ovendry mass while accounting for the prescribed water mass.
Independent shrinkage averages do not supply longitudinal strain, drying rates
or a temperature-dependent density law.

## Source-backed relaxation, creep, recovery and rate-dependent cycles

Two attributed CC BY 4.0 bundles add 24 claims, including 15 discrete measured
points. [Deng et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12348320/)
supplies the five-arm Hongda PU65 foam spectrum and density.
[Golovkin et al. (2026)](https://www.techscience.com/cmc/v86n3/65492/html)
supplies the two-arm printed REC TPU creep-fit spectrum, original viscosities,
first-cycle work, and tangent-modulus means and reported specimen standard
deviations. Printed TPU at 100% infill is distinct from flexible foam.

`fs_material::visco::source::ResolvedPronyReference` resolves explicit axial
Maxwell coefficients from original material claims and retains their receipts.
Its specimen shares one persistent history across strain and nominal-force
control. A strain hold relaxes stress; a force hold creeps; zero force allows
recovery. The response separates external model work, recoverable storage and
nonnegative integrated dissipation. Work may decrease during unloading.

These are small incremental, fixed-temperature models about the source operating
point. Strain, rate and duration bounds are authored model-use limits. Fitted
coefficients and these calculations do not establish an independently validated
finite-strain foam law, densification, permanent set or damage. Original length,
area and volume define conjugate engineering strain, nominal stress and energy.
Optional source-geometry bindings reject a mismatched specimen; missing density
remains `unstated` rather than preventing quasistatic mechanical evaluation.

For the Hongda PU65 source, the 25% preload does not replace the original 80 mm
height with 60 mm. The spectrum was fitted from 1–20 Hz measurements at that
preload. Its transient predictions are explicitly selected with
`--time-domain-reference`; the source retains `source-time-domain-validation-known=0`.

Compile to a new output filename, then run a 1 N creep and zero-load recovery:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/hongda-pu65-deng2025-prony-reference/manifest.tsv \
  --out /tmp/hongda-pu65-history.fsmatpk
cargo run -p fs-material --example material_viscoelastic_reference -- \
  /tmp/hongda-pu65-history.fsmatpk 0.08 0.0064 -0.1 0.1 1 2000 \
  --time-domain-reference \
  --source-length-axis source-original-height \
  --source-area-axis source-original-area \
  --axis source-sample-hongda-pu65-deng2025 1 \
  --axis source-test-temperature-known 0 \
  --axis source-moisture-condition-known 0 \
  --axis source-pressure-known 0 \
  --axis source-fitted-spectrum-reference 1 \
  --axis source-time-domain-validation-known 0 \
  --axis source-precompression 0.25 \
  --axis source-original-height 0.08 \
  --axis source-original-area 0.0064 \
  --axis source-frequency 1 \
  --axis source-acceleration-peak-to-peak 0.25 \
  --axis source-preconditioning-section2-2 1 \
  --force 1 1 --force 1 120 --force 0 1 --force 0 300
```

This reference calculation gives strains 0.00496924, 0.02197976, 0.01708397
and 0.00486161 at the four endpoints. Replacing the four controls with
`--strain 0.001 1 --strain 0.001 30 --strain 0 1 --force 0 60` demonstrates
relaxation and mixed control: the hold decreases stress from 30.84127 Pa to
9.63770 Pa. These values are model predictions, not added measurements.

The REC TPU EasyFlex pack preserves a two-arm creep-fit spectrum and measured
first-cycle loop work for a 20 mm × 20 mm, 100% infill printed cylinder. Its
density is absent. An explicit triangular reference at 12 mm/min is runnable as:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/rec-tpu-easyflex-golovkin2026-creep-reference/manifest.tsv \
  --out /tmp/rec-tpu-history.fsmatpk
cargo run -p fs-material --example material_viscoelastic_reference -- \
  /tmp/rec-tpu-history.fsmatpk 0.02 0.0003141592653589793 -0.05 0.2 1 6000 \
  --time-domain-reference --source-length-axis source-original-height \
  --axis source-sample-rec-tpu-easyflex-golovkin2026 1 \
  --axis source-density-known 0 --axis source-infill-fraction 1 \
  --axis source-original-height 0.02 --axis source-original-diameter 0.02 \
  --axis source-moisture-condition-known 0 --axis source-pressure-known 0 \
  --axis source-fitted-spectrum-reference 1 --axis source-fit-method-creep 1 \
  --axis source-canonical-parallel-maxwell-interpretation 1 \
  --axis source-time-domain-validation-known 0 \
  --axis source-time-domain-experiments-reported 1 \
  --axis source-loading-level-definition-ambiguous 1 \
  --axis source-test-temperature-known 0 --axis source-test-duration 6000 \
  --strain 0.15 15 --strain 0 15
```

The source's measured loops and this ideal triangular reference are different
histories. The paper's accuracy statement used reconstructed experimental
stress histories and must not be transferred to this prescribed triangle.
Changing both ramp durations gives the following explicitly reported comparison:

| Speed, mm/min | Each ramp, s | Measured loop work, J/m³ | Model loop work, J/m³ | Model dissipation, J/m³ | Remaining model storage, J/m³ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.12 | 1500 | 119000 | 33494.882 | 32353.493 | 1141.389 |
| 1.2 | 150 | 134000 | 67255.760 | 54974.683 | 12281.078 |
| 12 | 15 | 159000 | 124389.558 | 116327.531 | 8062.028 |
| 60 | 3 | 194000 | 146415.741 | 112166.839 | 34248.901 |
| 120 | 1.5 | 218000 | 102422.473 | 81137.174 | 21285.298 |

A closed strain path can retain internal spring energy, so loop work is not
automatically all dissipative heat. The source measurements remain discrete
`TabulatedOnly` claims. The example emits every source coefficient and receipt,
authored limits, per-step work/storage/dissipation and deterministic identities.
It admits at most 64 axes and 256 steps and validates the complete history before
printing a successful result; a late failure emits no partial history.
