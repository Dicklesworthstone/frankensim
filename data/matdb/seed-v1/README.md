# FrankenSim material seed data v1

This directory retains the human-reviewable source inputs compiled by
`xtask matdb-pack`. Generated binary packs are deliberately not committed:
the compiler must reproduce them from the pinned manifest and source record.

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

## OFHC Copper cryogenic exact-point tranche

`ofhc-copper-rrr100/` retains NIST's combined OFHC Copper scope, UNS
C10100/C10200, without pretending that the two UNS designations are
interchangeable for properties the source does not cover. Thermal-conductivity
claims pin `RRR = 100`, a load-bearing condition because NIST publishes
distinct conductivity correlations for RRR 50, 100, 150, 300, and 500. The
specific-heat source publishes one combined OFHC correlation without an RRR
condition, so its observation names that omission explicitly instead of
silently borrowing the conductivity condition.

The tranche evaluates the published RRR-100 conductivity correlation and
specific-heat polynomial at exactly `77 K` and `293 K`. As with the Aluminum
tranche, degenerate validity intervals prevent interpolation claims. NIST
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
