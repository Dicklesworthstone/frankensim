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
