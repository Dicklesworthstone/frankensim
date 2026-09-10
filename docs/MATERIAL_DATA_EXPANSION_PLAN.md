# Material database expansion plan

Planning snapshot: **2026-09-08**. Based on the current shared working tree, the [material inventory](MATERIAL_DATA_INVENTORY.md), the seed manifests and selected source records, and primary-source research. Proposed work below has not been implemented by writing this document.

The largest everyday gaps at planning time were **liquid water, commodity plastics, ordinary structural steel grades, and usable rubber compounds**. The first implementation now adds bounded liquid-water source data and demonstrates its use in steady conduction. The largest opportunity to reuse existing work is to complete **condition-compatible metal, wood, glass, and construction-material datasets** and connect them to actual consumers. Adding another isolated melting point or hardness reading usually delivers less than supplying the missing density, heat-capacity curve, or modulus that makes an existing simulation usable.

The catalog now has 180 source bundles, of which 162 are bulk-material bundles. Those bulk bundles have a median of four populated property names. Multiple bundles can describe the same material under different conditions; these are not 162 complete materials. The [inventory](MATERIAL_DATA_INVENTORY.md) gives exact counts, sizes, paths, and every populated property. “Missing” below means missing from the sourced seed catalog, not absent from every Rust constant, model, example, or test fixture.

P01 has a first named HDPE acquisition in `ensinger-tecafine-pe300-natural-2017`:
four producer comparison facts for density, tensile modulus, yield strength and
yield strain. The source's internal 2017 AA revision is preserved. Unknown test
temperature, process details and rate remain explicit; CLTE lacks its interval
and convention and is observation-only. These records do not complete P01 or
P11a. The missing positive delivery is a condition-compatible elastic or thermal
profile connected to its actual consumer; service-temperature limits cannot fill
that gap. A second reviewed producer candidate, INEOS HD6070FA cast film, states
density at 23 C but omits temperatures for Cp and conductivity; it is not merged
with Ensinger stock shapes to manufacture completeness. The next acquisition,
`roechling-polystone-g-natural-reference`, retains three same-source nominal
thermal inputs (density, Cp, conductivity) under an explicit frozen 25–26 C
reference opt-in. Its producer sheet does not state the property test
temperatures or selected product form/process. The existing PVC heat regression
now also runs this HDPE reference at identical geometry and boundary, with
source/discovery/store/resolver, independent response/energy and refusal checks.
Runtime verification is pending. P01/P11a remain unfinished for their original
condition-qualified delivery; the reference does not establish measured
temperature coverage, creep or product qualification.

P02 now has six producer comparison facts in
`mcam-proteus-homopolymer-pp-natural-2023`. Its conductivity is explicitly
at 23 °C and feeds the new source-card slab resistance consumer through
compilation, persistent store/reopen, discovery and resolution. The focused
integration passed remotely on 2026-09-09 after the sibling pager fix and
source-unit spelling correction. Cp, Poisson ratio and complete processing
qualification remain missing, so P02 and the complete PP consumer remain open.
No generic PP heat capacity is merged into this named grade.

P03 now includes `iplex-pvc-u-pipe-engineering-reference`: five nominal
producer-family inputs with extrusion retained and unknown exact formulation,
schedule and property test conditions. The 25–26 C frozen model requires
explicit opt-in; source service limits do not define property coverage.
Specific gravity uses a declared nominal conversion, and CLTE remains
observation-only. The source/store/reopen/discovery/resolver/lumped-heat test
passed remotely on 2026-09-09, including independent response and energy,
replay and unsupported-state/reference refusals. Exact compound qualification,
condition-dependent properties and pressure/creep consumers remain open.

P05 adds exact Makrolon 2405 PC and uncolored Terluran GP-22 ABS bundles,
each with density, tensile modulus and conductivity at the source's 23 C
table condition. PC additionally retains 50% RH for modulus/conductivity,
1 mm/min for modulus and through-plane conductivity. ABS conditioning,
rate and direction remain unknown. The existing source-card slab consumer
is exercised by a new import/store/reopen/resolver regression for both grades;
the PC discovery request is an executable example. The focused integration
passed remotely on 2026-09-09: both source packs survived store/reopen, PC
discovery completed, and 2 mm walls of area 0.01 m² gave 1 K/W for PC and
1.1764705882352942 K/W for ABS, with thickness scaling and unsupported-state
refusals. Neither producer source supplies Cp, and CLTE remains observation-only
where interval/convention is incomplete. Full P05 thermal/elastic profiles,
creep and housing qualification remain open.

Priorities are engineering judgment based on breadth of applications, FrankenSim's existing consumers, the size of the current gap, and likely acquisition effort. They are not a measured global consumption ranking. **A** means the next delivery tranche; **B** means the following expansion; **C** means a targeted application should pull it forward. Effort estimates concern the first useful, bounded dataset, not full physical qualification.

B01 now adds two explicitly opted-in reference bundles:
`pilkington-float-glass-reference` and `schott-borofloat33-reference`.
Each retains six manufacturer inputs, including interval-mean expansion under
a separate property name. The 25–26 C frozen-coefficient model domain is an
author-selected calculation boundary, not measured temperature coverage or an
accuracy guarantee. Source temperatures and omissions remain in each property's
observation; callers must acknowledge unmatched source temperatures. Pilkington
is a generic producer float-glass bulletin, not a qualified named pane product.
The shared heat-comparison regression exercises these reference inputs through
the existing compiler, store, discovery, resolver and lumped heat owner. This
focused test passed remotely on 2026-09-09 for both glasses, including independent
heat-response and energy checks, replay and unsupported-context refusals.
B01's fully condition-qualified profiles and the optical
and fracture consumer tasks remain unfinished.

Execution is tracked under **`frankensim-7sga6`**. The 2026-09-08 review preserved all 56 original tasks and added nine execution steps: separate HDPE/PP/PVC consumers, separate building structural/heat consumers, separate glass optics/fracture consumers, a sourced soft-sheet law adapter and the shared material-source E2E script. The graph has 71 records: six epics and 65 task records, including four tasks that now aggregate children. Every original material and consumer requirement remains; each task has explicit unit, source, domain, physical-reference and logged E2E acceptance. Nine material tasks now relate to MR10 instead of waiting for tensor infrastructure that already exists. M11's copper/aluminum conductors no longer wait for M02's potentially different metal pair. Existing MR14–MR18 obligations remain attached to the root.

Use `br show frankensim-7sga6`, `br ready` and `bv --robot-triage --label material-expansion` for current state. M00, M01, F01 and F02 retain their focused runtime closures. After the dependency correction, `bv` ranks B01 glass first; the bounded AGC/Pilkington/SCHOTT search has useful reference data but unresolved common-temperature/treatment and instantaneous-versus-mean expansion context, recorded in the still-open Bead. M02 now delivers the two reference-state metal profiles described below; the broader MR14–MR18 qualifications remain unfinished. MR10 is closed on the sourced silicon tensor reaching the actual oriented solid operator, with independent force/energy checks; this does not complete the alumina/silicon electronics task. The shared runner [scripts/e2e/material_sources.sh](../scripts/e2e/material_sources.sh), tracked by `.6`, passed eight exact source/consumer cases remotely on 2026-09-09 UTC, including dry air and both metal profiles. It retains numerical diagnostics and failing exit status; full DSR and experimental qualification remain outstanding.

M02 reuses the NASA 2024-T3 set at 300 K and adds
`stainless-316-20c-engineering-reference` at 293.15 K. Six source requirements
(rho, E, nu, Cp, k and instantaneous alpha) pass through compilation,
sealed/reopened storage, discovery and typed resolution into the existing
thermoelastic plate and radiator. Independent mass/modal/Zener references,
six receipts, unsupported-state/pin refusals and deterministic 128-step pressure
evolution pass. The same 0.20 × 0.15 × 0.0016 m geometry gives masses of
0.12984 kg and 0.384 kg. The 316 set explicitly joins supplier density/nu with
NIST equation evaluations at 20 °C; it is a nominal isotropic engineering
reference, not a common-coupon or product-form qualification. Different source
temperatures remain explicit. M03 warm-temperature acquisition, M13 optics and
MR15's broader acoustic/trajectory obligations remain open.

F03 adds the `air-dry-ussa1976` parameter pack and a small material-card adapter
into the existing `GasSpec`/`GasState` model. Its actual compiler, persistent
store/reopen, discovery, typed resolution and cylinder acoustic-loss path pass
at four temperature/pressure states. The five parameters derive the full gas
tuple over the declared 273.15–313.15 K / 80–110 kPa dry application range;
they are not independent measurements of every output property. Source-table
errata, the 1976 composition and the model's approximation limits remain
explicit. F04 now adds the water-vapor parameter source and reuses the existing
humidity path with ideal molar heat-capacity mixing, Wilke viscosity and WMS
conductivity. The ten-case source runner passes, including five humid states
through actual compilation, store/reopen, component discovery, resolution,
cylinder drag and classical sound absorption. Independent source-equation
references, exact dry recovery, humidity-basis conversions, receipts, replay
and unsupported-state refusals pass. Constant vapor heat capacity, ambient
Sutherland extrapolation and Eucken conductivity have no experimental error
bound; condensation and real-gas behavior remain outside this delivery.
Fourteen focused gas unit tests and three existing acoustic assembly tests
also pass remotely. The actual duct waveform shifts by 17.605 cents versus
the independent mixture prediction of 17.248 cents, within the unchanged
3-cent tolerance. This is numerical regression evidence, not experimental
qualification of the mixture or duct model.

M11 adds standard annealed 100% IACS copper (10–30 °C, NBS HB100) and
EC-H19 aluminum (0–30 °C, NBS HB109) volume-resistivity curves. These are
explicit handbook reference conditions, not inferred C11000 or 1350 grades.
The source compiler, persistent store/reopen, discovery and typed conductor
adapter now reach the existing circuit DAE. The nine-case source runner passed
remotely on 2026-09-09 UTC, including independent resistance, voltage and
Joule/supplied-energy checks at eight knots and two interpolated states,
geometry scaling, receipts, replay and unsupported-input refusals. This closes
the missing source-to-DC-consumer connection; temperature feedback, ampacity,
AC/contact effects and broader electrical qualification remain separate work.

F05 now supplies two bounded coolant profiles: DOWTHERM SR-1 ethylene glycol
and DOWFROST propylene glycol, each at exactly 50% glycol by volume and 40–60 °C.
Four five-knot curves per product resolve through the new liquid adapter into
existing lumped heat and LBM channel calculations. The eleven-case source runner
passed remotely on 2026-09-09 UTC, covering literal source values, interpolation,
persistent store/reopen, discovery, receipts, independent heat/flow references,
replay and unsupported-state refusals. F05 remains **in progress**: the guides
do not state the volume-reference temperature or exact 50-vol% phase endpoints.
Neighboring phase-table rows cannot silently supply those missing claims. This
delivery does not qualify mass/volume conversion, other concentrations, evolving
transport coefficients, boiling/freezing or experimental accuracy.

F06 delivers a bounded mineral-oil engineering profile from Shell's May 2011
Heat Transfer Oil S2 data sheet: density, heat capacity, conductivity and
kinematic viscosity at four common temperatures over 0–200 °C. An explicit
kinematic adapter derives dynamic viscosity using density at the same resolved
point and preserves the original four receipts. The twelve-case source runner
passed remotely on 2026-09-09 UTC, including real compilation, persistent
store/reopen, discovery, gated heat and LBM channel flow at four knots and one
interior interpolated state. Independent formulas, replay and domain/basis
refusals pass; three liquid unit tests also pass. The source's printed Prandtl
numbers disagree with its other tables, and that unresolved inconsistency is
retained explicitly. Sparse interpolation and derived transport values are
engineering approximations, not measured validation. Pressure dependence,
aging, friction, hydraulic-grade substitution and experimental qualification
remain outside this delivery.

## Materials to add or complete

| Priority / effort | Family and first targets | What exists now | Most valuable missing coverage and resulting use |
| --- | --- | --- | --- |
| A / low–medium | **Pure liquid water; dry and humid air** | New IAPWS source pack: five 15-knot liquid-water curves over 10–80 °C at exactly 0.1 MPa, with declared enthalpy reference and interpolation/consistency limits. Its compiled/stored conductivity now drives two actual steady-conduction solves with state-dependent flux and source-domain refusals. Water-vapor species metadata, air reference-gas composition and one water/ethylene-glycol formulation remain separate. Humidity already affects an acoustic gas-state path. | Extend water to the selected transient-heat or flow consumer; complete air transport and thermodynamics matched to composition, pressure, temperature and humidity. Enables ordinary cooling and fluid experiments and sourced acoustic inputs. Extend to boiling/steam/ice only with the appropriate phase models. |
| A / medium | **6061-T6 aluminum, C11000 copper, annealed 304/316 stainless; existing 2024-T3 reference** | 6061 has three 24-knot curves over 77–293 K; OFHC has two. 304/316 retain five cryogenic curves each, including derived instantaneous expansion, with NIST product form/heat treatment unspecified. The separate 316 engineering reference now adds a complete six-property set at 293.15 K using attributed supplier density/nu and NIST equations; both it and the 2024-T3 300 K profile drive the actual thermoelastic plate consumer. The complete profiles do not inherit the wider individual curve coverage. | Extend compatible full profiles to finite heating domains and additional everyday grades/forms. Add electrical resistivity for conductor consumers and sourced loss models for broader acoustics. Preserve 304/316, copper purity/temper and cross-source qualification distinctions. |
| A / medium | **A36 plate, A572 Grade 50, A500 Grade C tube, 1018 mild steel; A615 Grade 60 rebar** | 1045/4130/4140 and specialized bearing/carburizing steels; no dedicated bundles for these proposed ordinary grades. | Begin with one product form and grade: density, elastic properties, thermal properties, then stress–strain/hardening where plastic deformation is requested. Expands buildings, frames, brackets and machinery. Rebar remains a separate reinforcement constituent; it does not turn concrete into one isotropic material. |
| A / medium | **HDPE, LDPE/LLDPE, polypropylene, rigid PVC and plasticized PVC** | No dedicated PE/PP bulk bundles. The PVC source is low-density foam, not pipe-grade solid PVC or flexible cable insulation. | First one named HDPE grade, one PP grade and rigid PVC formulation, with density, elastic response, heat capacity, conductivity and expansion at stated conditions. Then creep/viscoelasticity and processing dependence. Enables containers, pipes, housings, films and insulation without substituting foam data. |
| A–B / medium | **PET, ABS, polycarbonate, PMMA, POM, PA6/PA66; PLA and PETG for printed objects** | Makrolon 2405 PC and Terluran GP-22 ABS now have three-property room-state bundles; Cp and complete elastic/thermal profiles remain missing. Other listed everyday resins still lack dedicated bulk bundles. | Complete the ABS/PC profiles, then PET packaging and nylon/POM parts. Retain molding direction, fillers, moisture conditioning and rate. Printed PLA/PETG require build direction and process-specific data; injection-molded values cannot qualify printed parts. |
| A–B / low–medium | **Architectural soda-lime float glass and borosilicate glass** | Soda-lime SRM has density and expansion only. Duran borosilicate has eight metrics, with ambiguous strength entries deliberately withheld. | Complete mechanics and thermal properties for one named annealed float-glass product. Add wavelength-dependent optical data for rendering. Treat heat-strengthened/tempered glass, surface damage and fracture statistics separately; a tabulated breaking stress is not a universal material constant. |
| B / medium–high | **EPDM, silicone, NBR, natural rubber/SBR, neoprene and TPU** | One NBR compound's fuel-sorption/swelling data; specialized tissue/contact sources. No broad rubber mechanical/thermal bundles. | Start with one documented EPDM or silicone compound for a seal or compliant pad: density, compression/stress–strain response, relaxation and thermal behavior. Add matched friction only when a counterface is selected. Hardness or a chemical-compatibility chart alone cannot supply a nonlinear constitutive law. |
| B / medium–high | **Concrete, cement paste/mortar, brick and gypsum** | CFAST thermal model inputs and one concrete specimen with strength, density and static/dynamic modulus. | Complete one documented normal-weight concrete mix and one wallboard/mortar product: moisture, cure age, density, thermal curves, expansion and mechanical response as needed. Fire, dehydration, cracking and creep require their own supported models. Cement powder, cured paste, mortar and concrete remain different identities. |
| B / medium | **Douglas-fir and Southern Pine lumber; plywood, OSB and pressure-treated lumber** | Many wood species and directional ratios; a few plywood/OSB source populations. The treated-lumber pack contains treatment/process thresholds, not a complete treated-material law. | Complete one species/grade/moisture condition and one panel layup with compatible stiffness, density, thermal and moisture expansion data. Add damping for acoustic use and time-dependent response for sustained loading. Treated stock needs preservative, retention, incising, redrying and moisture state. |
| B / medium | **MDF, particleboard, LVL, glulam and CLT** | No dedicated grade/layup bundles for these products. | Start with MDF or particleboard for furniture/panels, then engineered structural products. Preserve adhesive, density, layer axes, layup and moisture; do not relabel clear-wood data as a finished panel. |
| B / medium | **EPS/XPS, mineral wool and flexible polyurethane foam** | Glass-fiber, cellulose, rigid urethane and other CFAST insulation inputs; some PVC foam values. | Extend temperature/moisture/density dependence for actual insulation products, then compression/viscoelastic data for cushioning. Open/closed cells, fill gas and aging matter. Thermal insulation data do not establish porous acoustic absorption. |
| B / medium–high | **Epoxy adhesives, silicone sealants, FR-4 PCB laminate, alumina and silicon** | A specialized epoxy/insulation source and magnet-wire data; no general PCB laminate, silicon or alumina bundles. | Choose one cured adhesive and one exact FR-4 construction first. Add anisotropic thermal expansion/conductivity, dielectric properties versus frequency, and mechanics. Adhesive joints need bondline and interface data; an epoxy resin card cannot stand in for glass-fiber laminate. |
| B / medium | **Hydraulic/mineral oils and water–glycol coolants** | Shell Heat Transfer Oil S2 now has four typical-design curves at 0–200 °C, with its printed Prandtl inconsistency retained. Both Dow glycol families have four curves at 40–60 °C and exactly 50% glycol by volume. | Extend to a named hydraulic-oil grade and additional coolant concentrations. Preserve sparse-interpolation limits, pressure/additive/aging context, and missing coolant phase endpoints and concentration-reference temperatures. Heat-transfer mineral oil does not qualify a hydraulic grade. |
| B–C / high | **Asphalt/bitumen, sand, gravel, clay/soil, granite and limestone** | No dedicated bulk bundles for these proposed classes. | Start with a selected pavement mix or rock application, then soils with density/porosity/saturation and stress history. A generic “soil” modulus would be misleading; these need stateful granular/porous constitutive support as well as data. |
| B–C / medium–high | **GFRP and CFRP laminates** | No dedicated general engineering laminate bundles; isolated fiber-filled interface materials are not substitutes. | Name fiber, resin, fraction, cure and layup. Supply orthotropic stiffness and thermal response, then failure and interlaminar properties for the requested analysis. Reuse MR10's tensor/frame machinery. |
| B–C / medium | **Paper/cardboard, cotton, polyester textiles, felt and leather** | Specialized cellulose insulation, Nomex and biological source records, not complete everyday sheet/fabric bundles. | Packaging and soft objects need direction, thickness or areal density, moisture, bending response and loss. Woven textile effective behavior differs from a bulk fiber's modulus. |
| B–C / medium | **5052 sheet, 6063 extrusion, A356 cast aluminum; 430/316L stainless; ductile iron and bearing bronze C93200** | Related metals and a specific gray-cast-iron microstructure exist, but not these named product conditions. | Expand common metal product forms after the first compatible substitution set works. Keep cast/rolled/extruded states and gray/ductile graphite morphology separate. Brass and bronze already have source fragments worth completing before adding many more names. |
| C / medium–high | **SAC305 and Sn63Pb37 solder; titanium Ti-6Al-4V, magnesium AZ31, tungsten, silver and gold** | Pure tin/lead and titanium-related packs cover only portions of the required behavior; no dedicated bundles for these proposed alloy/product states. | Solder deserves early promotion for electronics: phase interval, enthalpy, creep, conductivity and joint conditions. The other materials follow aerospace, tooling, contact or optical demand. Pure tin/lead data cannot simply be mixed into a qualified solder alloy. |
| Existing flagship / high | **Complete the lead heating/melting specimen; electrical steel and magnets** | Multiple lead phase/transport/expansion sources; several M19, N42 and Y30 source-specific records. | Continue MR14's lead bundle because it blocks an explicit project demonstration. For machines, complete one supplier/process-specific magnetic dataset including temperature, field/history and loss domain. More names do not resolve these existing compatibility and model gaps. |

The table intentionally includes both absences and partial entries. Aluminum, copper, tin, iron, wood and concrete are already represented; calling them wholly missing would obscure the actual work.

## Define completeness by the calculation

Use the existing material requirements and gap-query path. A material need not have every property to support a useful calculation. These are acquisition profiles, not a new database schema or certification system.

| Intended calculation | Minimum data to assemble for the chosen model | Additional data only when the calculation needs it |
| --- | --- | --- |
| Linear elastic solid or undamped modal frequencies | Density and elastic constants, with condition and direction; for a general isotropic solid, E and Poisson ratio | Plasticity, fatigue and damage for failure; frequency/rate-dependent loss for damping and decay; geometry-specific reductions can require fewer constants |
| Transient heat conduction | Density, heat capacity/enthalpy relation and conductivity over the requested trajectory | Surface radiation/convection inputs for those boundary conditions; latent heat and phase law for phase change |
| Thermally induced deformation | Elastic properties plus thermal free strain or expansion with a reference state, joined to the heat profile | Restraint/prestress belongs to the specimen; yield, creep and annealing when the trajectory reaches those regimes |
| Liquid heat transport | Density/EOS as required, viscosity, heat capacity/enthalpy and conductivity | Compressibility/sound speed for acoustics; surface tension, wetting and phase relations for free surfaces or boiling |
| Rubber, polymer or wood time response | Appropriate nonlinear/directional stiffness and measured relaxation, creep or dynamic loss over the declared rate/frequency range | Temperature/moisture shift models, aging and cyclic history only with applicable evidence |
| Electrical, magnetic or optical response | Conductivity/resistivity, dielectric/magnetic response or spectral optical properties required by the selected observer | Frequency, temperature, field, surface finish, coating and hysteresis dependence where relevant |
| Contact or joint | Constituent properties plus the selected ordered interface law and its load/rate/environment domain | Wear, thermal/electrical contact conductance, adhesion and wetting for those actual consumers |

For example, six suitable mechanical/thermal properties can make a metal useful for elastic thermal deformation, while a dozen composition or hardness metrics may not. Heat simulation does not need an optical spectrum unless that observer or boundary model requests it.

## Execution order

### 1. Turn the best existing source sets into usable selections

- [x] Inspect the existing six-property thermoelastic requirement set and reuse the 2024-T3 300 K input set. Its typed conductivity already works in the resolver; no new adapter was needed. The 6061/C11000/304 wider complete-profile gaps remain separate acquisition work.
- [x] Select the first two materially distinct bundles: 2024-T3 at 300 K and an explicitly cross-source 316 engineering reference at 293.15 K. Both satisfy the same six-property point profile; this does not certify arbitrary stock or a same-temperature substitution.
- [x] Extend the 6061 and OFHC source-fit point samples into supported curves where the published equations and existing compiler can do so. M01 delivered five 24-knot curves over 77–293 K; source, compiler/store/discovery, deterministic rebuild and original endpoint checks passed. Fit/data limits and sampled interpolation discrepancies remain explicit.
- [x] Add 316 density and Poisson ratio from explicit 20 °C supplier facts, with documented compatibility limits, and directly evaluate four NIST inputs at that temperature. Keep alpha=epsilon'/(1+epsilon), source reference length and unknown uncertainty distinct from interval-mean expansion. Existing curve packs and semantic mappings remain intact.
- [x] Run both selected sets through the actual compiler → store/reopen → discovery → resolver → thermoelastic plate/radiator. The new case and all seven earlier source/consumer cases passed remotely; independent equations and evolved pressure sensitivity provide the physical positive result.

**Temperature targets must be stated honestly.** A useful subsequent heating target is 20–100 °C, but it is a proposed acquisition domain, not current coverage. Cryogenic curves ending near 300 K cannot meet it. Even the existing stainless modulus endpoints differ: 304 ends at 293 K and 316 at 294 K. A room-temperature point does not support a finite heating trajectory. Use the source-supported interval for the first delivery, then obtain compatible warmer data explicitly.

### 2. Add water and the first commodity plastics

- [x] Implement a bounded pure-liquid-water dataset using the IAPWS liquid-water release: F01 delivered five 15-knot curves over **10–80 °C at 0.1 MPa**, with source equations, units, Cp/enthalpy approximation consistency, actual compiler/store/discovery and unsupported-state checks passing.
- [x] Bind water to an existing heat/fluid consumer. F02's real source/compiler/store/resolver/conduction path passed manufactured temperature/flux and energy checks at **20–25 °C and 60–65 °C**, exactly 0.1 MPa, liquid phase. The warmer interval gives 8.5% greater mean flux for the same 5 K/m gradient. Unsupported temperature/pressure/phase queries refuse. This establishes steady conduction, not fluid motion, transient heating or boiling.
- [ ] Acquire one named HDPE grade, one PP grade and one rigid PVC formulation. Follow with ABS and PC. Select a bounded elastic or thermal requirement set first; obtain the missing quantities instead of declaring a product complete from its datasheet title.
- [ ] Add one ordinary mild-steel product condition and one annealed float-glass product. Reuse the same consumer-level checks where applicable.
- [x] Complete a bounded dry-air profile around the existing gas model: F03's five sourced parameters resolve through the real compiler/store and gas adapter into acoustic cylinder loss, with independent property/loss references, source sensitivity, replay and input refusals. No parallel EOS was added.
- [x] Complete bounded humidity-dependent source coverage and transport for the existing acoustic consumers (F04): the ten-case source runner passes the five-state humid-air mixture, cylinder-drag and sound-absorption checks. Separate full airflow/heat evolution and experimental qualification remain unclaimed.

### 3. Make buildings, furniture and compliant parts useful

- [ ] Complete one lumber species/grade/moisture condition and one plywood or OSB layup; add MDF or particleboard next.
- [ ] Complete one named concrete mix and one wallboard/insulation product for a chosen thermal or mechanical calculation. Extend existing CFAST input sets without promoting their model inputs into experimental measurements.
- [ ] Acquire one EPDM or silicone compound with usable deformation and relaxation data; pair it with a named counterface only when a contact simulation needs it.
- [ ] Add the next commodity polymers and coolant/oil formulations from the priority table. Attach directional, moisture and temperature support to the first consumer rather than treating these as future metadata cleanup.

### 4. Extend by demonstrated application demand

- [ ] Add FR-4 and solder for electronics; engineered laminates for structural composites; asphalt/soil for civil simulation; paper/textile systems for packaging and soft objects.
- [ ] Continue the existing lead heating/melting acquisition alongside the bounded deliveries above. Its phase, surface and specimen compatibility gaps remain open until the MR14 requirements are met.
- [ ] Add further metal grades, magnetic products and optical spectra when an existing demonstration names the need. Do not hold useful thermal or elastic releases for unrelated missing optical or failure data.

## Primary sources and what to take from them

These are verified source leads, not permission to copy whole publications or a claim that every required value is available. Apply the existing source/redistribution policy to each selected table or dataset. Manufacturer typical values, government model inputs, measurements and derived fits retain their distinct meanings.

| Source | Immediate use | Boundary to preserve |
| --- | --- | --- |
| [IAPWS SR6-08(2011), liquid water](https://iapws.org/technical-guidance/release/LiquidWater) | Compact thermodynamic, viscosity and thermal-conductivity correlations near 0.1 MPa | The release includes metastable regions and property-specific lower limits. Select a stable-liquid subset first; it is not a universal water/steam or saltwater model. |
| [NIST dry-air formulation](https://www.nist.gov/publications/thermodynamic-properties-air-and-mixtures-nitrogen-argon-and-oxygen-60-2000-k-pressures) | Reference thermodynamics for dry air and nitrogen/argon/oxygen mixtures | Humidity and transport need their own applicable sources/models. Use published equations through the existing pure-Rust model architecture, not a foreign runtime dependency. |
| [NIST structural-steel physical properties](https://www.nist.gov/publications/physical-properties-structural-steels) and [mechanical properties](https://nvlpubs.nist.gov/nistpubs/Legacy/NCSTAR/ncstar1-3d.pdf) | Thermal/elastic/mechanical source investigations | Retain the actual steel/specimen and test direction; a similar specified yield strength does not establish A36/1018 interchangeability. |
| [Borealis BE52 polypropylene](https://www.borealisgroup.com/products/product-catalogue/be52) | Example of a named commodity resin and downloadable grade datasheet | This is a blow-molding/thermoforming homopolymer, not every PP formulation; missing thermal or dynamic properties still need acquisition. |
| [Covestro Makrolon 2405 polycarbonate](https://solutions.covestro.com/en/products/makrolon/makrolon-2405_000000000000945088) | Named-grade mechanical, creep, thermal and electrical leads | Preserve test speed, conditioning, direction and property-specific temperature coverage. Heat-deflection temperature is not a melting point or a full creep law. |
| [USDA FPL Wood Handbook, revised 2021](https://research.fs.usda.gov/fpl/wood-handbook) | Extend sources already in the repository into moisture, thermal and engineered-product coverage | Clear wood, graded lumber and panels have different populations and allowable-property meanings. |
| [Pilkington ATS-129 float-glass bulletin](https://www.pilkington.com/-/media/pilkington/site-content/usa/window-manufacturers/technical-bulletins/ats129propertiesofglass20130114.pdf) | Mechanical and thermal starting set for ordinary soda-lime float glass | Its fracture figures specify treatment, loading duration, surface condition and breakage probability; preserve these rather than inventing one universal strength. |
| [Parker O-Ring Handbook](https://www.parker.com/content/dam/Parker-com/Literature/O-Ring-Division-Literature/ORD-5700.pdf) | Identify suitable compounds, application conditions and further data needs | Handbook family/compatibility information is not sufficient identification of hyperelastic or viscoelastic parameters. |
| [Isola FR408HR](https://www.isola-group.com/pcb-laminates-prepreg/fr408hr-laminate-and-prepreg/) | Exact laminate/construction candidates with directional and frequency-dependent data | Keep glass style, resin content, test method and direction; it is one FR-4 product system, not universal FR-4. |

Combining compatible sources is legitimate when their material population and conditions support it. Requiring every property to come from the same physical coupon would unnecessarily obstruct useful engineering models; combining unrelated grades or heat treatments is equally unsound. Document the actual compatibility basis in the existing material selection. If only a family-level estimate is supportable, expose it as an explicit engineering estimate through the existing research/override policy, with its assumptions visible. Do not silently fill a missing measurement or relax normal admission.

## What counts as delivery

For each selected material and use case:

1. Preserve the source identity, condition, units, validity domain and uncertainty meaning in the existing static files under `data/matdb/seed-v1/`. Add curves, tensors or model parameters only where supported. Different temperature/frequency knots do not count as additional material properties.
2. Compile and ingest with the existing pack/store tooling. Exercise discovery and material resolution at representative interior points and boundaries, including an unsupported-state refusal. Reuse current tests; build no new general harness.
3. Run the selected existing physical consumer and show a material-dependent result with a suitable independent reference or bounded numerical check. If consumer integration is still blocked, report the dataset as acquired/queryable and keep the physical-delivery item open.
4. Report coverage by material **and use case and operating domain**, with measured, fitted and authored inputs distinguished. Update the existing inventory. Track completed consumer requirement sets as the main progress measure; material-name count and megabytes are secondary.

Reuse the current implementation ownership: **MR14** for the lead specimen, **MR15** for ordinary metal substitutions, **MR16** for its oriented-solid/compliant-contact targets, **MR17** for interfaces and **MR18** for corpus-to-consumer integration. The [existing plan](MATERIAL_REALITY_IMPLEMENTATION_PLAN.md) retains their acceptance criteria and exact Beads IDs. New families outside those scopes should receive small implementation tasks when their first consumer and source tranche are selected, not another speculative all-materials epic.

The immediate recommendation is to finish a small compatible metal set, add liquid water and the first commodity plastics, and deliver their actual simulations. Then deepen construction, wood and rubber coverage. This sequence addresses common objects while continuing the project's lead-melting flagship, and makes existing data useful before accumulating another long list of partial entries.
