# Coupled finite-rate compartment example

Run the complete source-owned path from a heated, receding solid into reacting,
ventilated gas and finite radiating walls:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- all 0.1 0.001 12 kinetic be
```

Arguments are `[all|closed|vented] [end_s=0.1] [max_dt_s=0.001]
[solid_cells=12] [kinetic|frozen] [be|sdirk2]
[hydrogen|methane|charring-methane]`. The default fuel is hydrogen and the
default chemistry method is backward Euler (`be`). CSV state histories go to stdout; work counts
and maximum conservation residuals go to stderr. The timestep is an upper bound:
the original network and chemistry owners can retry shorter intervals. The
example refuses an exhausted source domain or work budget instead of clipping
temperatures, supplying seed species, or repairing physical balances.

Select the optional two-stage SDIRK2 homogeneous reactor method with the sixth
argument:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- all 0.1 0.001 12 kinetic sdirk2
```

Both methods use the same original source, constrained Newton solver and
physical acceptance checks. The CSV `method` column and summary log distinguish
them. Frozen-chemistry mode does not advance the reactor, regardless of the
selected method. SDIRK2 is a second-order reactor formula on smooth states; it
does not upgrade the first-order surface/gas/network composition to global
second order.

The seventh argument selects the full methane mechanism and, optionally,
retained carbon in the condensed fuel:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.001 12 kinetic sdirk2 methane
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- all 0.1 0.001 12 kinetic sdirk2 charring-methane
```

These use the same surface, vent, gas and audit code. Species are resolved from
the original canonical source names once, and every source species remains in
the network and chemistry inventories. The CSV includes selected CO, CO2 and
NO masses, fuel and hydrogen masses, and retained char. Zero columns for
species absent from the hydrogen mechanism do not imply a chemical prediction
for those species.

## Physical owners

- `fs-thermochem::kinetics::hydrogen::cantera_h2o2` provides the complete pinned
  Cantera 3.2.0 hydrogen/oxygen submechanism: ten species and 29 reversible
  reactions, including its third-body, explicit-collider, duplicate and Troe
  channels. The example explicitly selects `ContinuousAt300K`, retaining the
  original heat capacities and low-temperature reference while using the source
  owner's bounded upper-region enthalpy/entropy continuity treatment. This is
  the hydrogen subset.
- `fs-thermochem::kinetics::gri30::cantera_gri30` supplies the complete 53-species,
  325-reaction GRI-Mech 3.0 source for methane. Its 16 originally irreversible
  channels retain their direction, with signed chemical entropy diagnostics;
  no reverse rates are invented. The same explicit continuity policy retains
  original Cp and bounds high-region enthalpy/entropy adjustments. Its common
  source domain is 300–3000 K. Loading this source does not establish detailed
  balance, a general entropy theorem, or experimentally accurate emissions.
- `ZoneKinetics` advances the actual transported gas inventory in each network
  pressure trial. Formation energies, reverse rates and pressure dependence
  come from the same original mechanism. There is no prescribed heat-release
  rate, additional heat of combustion, global reaction or radical seed.
- `RecedingSlab` supplies implicit finite-volume conduction, finite-rate
  decomposition, volatile enthalpy and conservative remapping of the shrinking
  solid. `with_char` additionally retains the separate carbon mass, calorics
  and phase volumes through remapping and composition-dependent conduction.
  `VolatileMixture` transfers actual integrated pure-H2 or pure-CH4 release
  using the same NASA gas enthalpy and energy references.
- `Wall1d` supplies a separate finite wall's stored heat and conduction.
  `radiation::reactive` solves both surface temperatures with the existing
  participating enclosure. The original material-card owner supplies the
  explicitly authored surface emissivities.
- `MolecularPlanck` uses actual local H2O, CO2, CO and CH4 partial pressures and
  temperature in the published TNF/RADCAL Planck-mean fits. Other species are explicitly treated
  as transparent. Using this emission mean as gray absorption across both
  surface spectra and the finite path is an approximation, not spectroscopic
  transfer or a certified optical error bound.
- `ZoneNetwork` and the original stratified opening law solve ventilation
  against the chemistry-inclusive pressure. The returned solid and wall
  candidates are installed only with their matching accepted gas step.

Mechanism source:
[Cantera 3.2.0 h2o2.yaml at its pinned release commit](https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/h2o2.yaml).
Methane source:
[Cantera 3.2.0 gri30.yaml at the same pinned commit](https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/gri30.yaml).
Optical source: [TNF radiation reference](https://tnfworkshop.org/radiation/),
corrected 7 June 2000; its fit is historical and deliberately simple.

## Authored numerical fixture

**The condensed fuel is a manufactured H2- or CH4-releasing solid. It is not PMMA,
wood, metal-hydride data, or an experimentally validated material.** Its
composition is a bookkeeping choice that permits an elemental balance with
the actual source mechanism. The charring variant is defined by mass as
80% releasable CH4 and 20% retained condensed carbon, so its initial carbon
and hydrogen inventory is explicit. Its kinetic and thermal values are explicit
inputs, not a calibration or a claim that this material exists.

| Input | Value |
| --- | --- |
| Fixed gas volume | 0.02 m³: 0.1 m² floor area and 0.2 m height |
| Initial gas | Hydrogen: 1000 K; methane: 1500 K. 101325 Pa, 10% O2 / 90% N2 by mole |
| Gas layering | One active upper layer; the lower layer remains absent |
| Ambient reservoir | Same authored diluted oxidizer at 300 K and 101325 Pa |
| Vented opening | 0.2 m high, 0.001 m wide, discharge coefficient 0.7 |
| Solid | Area 0.01 m², thickness 0.006 m, density 850 kg/m³ |
| Solid thermal law | 1500 J/(kg K), conductivity 0.2 W/(m K), initially 500 K |
| Solid decomposition | First-order Arrhenius: 10⁸ s⁻¹, activation energy 100 kJ/mol |
| Solid formation reference | Yield-weighted original volatile and char enthalpy at 300 K minus 1 MJ/kg virgin |
| Optional char | 20% mass yield, 300 kg/m³, 1000 J/(kg K), zero authored h at 300 K, 0.08 W/(m K) |
| Char geometry/conduction | Additive phase volumes, volume-fraction parallel conductivity, fixed back |
| Solid temperature domain | 300–1200 K |
| External absorbed irradiation | 450 W = 50000 W/m² × 0.01 m² × 0.9 |
| Finite wall | Area 0.01 m², thickness 0.02 m, density 500 kg/m³, four cells |
| Wall thermal law | 1000 J/(kg K), conductivity 1 W/(m K), initially 300 K |
| Solid / wall emissivity | 0.9 / 0.7, constant diffuse-gray authored cards |
| Gas contact coefficients | 10 / 8 W/(m² K) for solid / wall |
| Optical geometry | Reciprocal normal chord between two 0.01 m² patches; 0.2 m gas path |
| Molecular fit domain | 300–2500 K; no extrapolation |
| Chemistry local controls | Relative tolerance 10⁻⁵; absolute amount tolerance 10⁻¹² mol |

The patch geometry is not a whole-room enclosure. Unmodeled room surfaces are
adiabatic. Solid recession does not change vessel volume or introduce pressure
work. Both condensed back faces are insulated. The hot gas initialization
provides a definite chemical initial state without pretending to model an
ignition device.

## Balance and comparison

Every accepted step checks clocks, nonnegative gas inventories, species
increments, complete gas-plus-solid molecular mass and elemental inventories,
and gas-plus-solid-plus-wall physical energy. The gas owner's aggregate source
report contains internal surface heat, volatile enthalpy and numerical
chemical-reference changes. The example separates these to recover the
external vent streams, then adds the absorbed heater exactly once to the
combined storage balance. The same separation is applied to vent species.

`net_chemical_fuel_consumption_kg` is the integrated **net** selected-fuel species change
from all source reactions, including reversible dissociation. It is not a
single reaction extent or an instantaneous complete-burning estimate.
`gas_radiation_heat_w` is positive into gas. `endpoint_kappa_m_inv` is evaluated
at the actual accepted gas endpoint; the coupled surface calculation used the
previous accepted gas state over that interval.

Useful comparisons keep the source data fixed:

```bash
# Finite-rate chemistry versus frozen chemistry; radiation still follows composition.
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.001 12 frozen

# One timestep and one spatial refinement of the same closed case.
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.0005 12 kinetic
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.001 24 kinetic
```

Gas temperature, opacity and optical paths are frozen over each surface-source
interval. Chemistry uses the selected original implicit method inside the
pressure solve. These compositions do not establish a globally second-order
scheme or a fully discrete entropy theorem. The model is homogeneous ideal
gas and authored condensed material, with no resolved flame structure,
diffusive boundary layer, soot-production model, sprays, char oxidation, pore flow,
real-gas phase or experimental fire validation claim.

## Executed backward-Euler baseline

The following five historical hydrogen runs used the backward-Euler reactor and reached
**0.1 s** with the optimized Rust binary, using the
original production fire and chemistry modules through the focused native
manifest. This is executed numerical integration evidence, not a full-workspace
DSR gate, G2 validation, or comparison against an experimental fire. Each run
completed its mass, elemental, energy, source-domain and shared-clock gates.

| Case | Maximum step (s) | Solid cells | Gas T (K) | Pressure (Pa) | Solid face T (K) | Gas H2O (kg) | Opacity (m⁻¹) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Closed, kinetic | 0.001 | 12 | 1331.047977 | 137365.460154 | 575.581792 | 1.620425103e-4 | 0.153803437 |
| Vented, kinetic | 0.001 | 12 | 1295.158522 | 114614.158554 | 573.930881 | 1.481052751e-4 | 0.144886061 |
| Closed, frozen chemistry | 0.001 | 12 | 989.028059 | 103913.655032 | 563.412549 | 0 | 0 |
| Closed, kinetic, smaller step | 0.0005 | 12 | 1330.604570 | 137322.649575 | 575.647817 | 1.619122982e-4 | 0.153736972 |
| Closed, kinetic, finer solid | 0.001 | 24 | 1331.611157 | 137428.397829 | 548.250225 | 1.623552811e-4 | 0.154027588 |

The closed kinetic baseline emitted 1.816062521e-5 kg H2 and consumed a net
1.815734218e-5 kg through the source reactions. Its remaining gas H2 was
3.283029734e-9 kg. The vented case transported a net 9.957055456e-4 kg out of
the gas volume, with the actual donor composition and enthalpy. At the final
baseline intervals, radiation removed 102.880729 W from the closed gas and
86.920753 W from the vented gas. In the frozen comparison no water forms, so
the selected gas model stays transparent while surface-to-surface radiation
and gas convection remain active.

| Run | Maximum combined mass residual (kg) | Maximum combined energy residual (J) | Maximum elemental residual (mol) | Accepted chemistry substeps |
| --- | ---: | ---: | ---: | ---: |
| Closed, kinetic | 5.636e-17 | 8.293e-11 | 5.593e-14 | 79392 |
| Vented, kinetic | 1.522e-17 | 7.617e-11 | 1.548e-14 | 73235 |
| Closed, frozen | 3.123e-17 | 7.348e-11 | 3.098e-14 | 0 |
| Smaller step | 2.288e-17 | 2.517e-10 | 2.237e-14 | 94417 |
| Finer solid | 4.834e-17 | 1.637e-10 | 4.775e-14 | 79505 |

All runs used 100 accepted network intervals, except the smaller-step run's
200. None needed a network halving. The maximum vent pressure residual was
9.509e-6 Pa; the maximum surface-power residual over all five runs was
4.318e-7 W. Chemistry counters describe the payloads of accepted network
intervals; discarded pressure trials perform additional work. This
backward-Euler stiff-chemistry baseline is expensive; use the optimized build.

Halving the network step changes final gas temperature by 0.443407 K. That
small bulk change does **not** establish convergence of trace species: final
residual H2 changes from 3.2830e-9 to 9.9499e-9 kg, and OH from 3.9220e-7 to
5.2215e-7 kg. Increasing solid resolution from 12 to 24 cells changes gas
temperature by 0.563181 K and emitted mass by about 0.193%, but changes the
solid face temperature by **27.331567 K**. These two coarse meshes do not
resolve the early surface thermal transient well enough for an accurate
surface-temperature or pyrolysis prediction. The runs expose that sensitivity;
they do not claim a convergence order, a mesh-independent flame prediction,
or that small balance residuals imply small modeling/discretization errors.

## Executed optional SDIRK2 comparison

The same hydrogen source data, geometry, solid/wall owners and local chemistry controls
were rerun with `kinetic sdirk2`. Both the closed and vented 0.1 s cases passed,
followed by the single closed timestep refinement. All use 12 solid cells.
Backward Euler remains the default; these are explicit method selections.

| Case | Maximum step (s) | Gas T (K) | Pressure (Pa) | Gas H2O (kg) | Opacity (m⁻¹) | Accepted chemistry substeps | Rejected chemistry steps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Closed, SDIRK2 | 0.001 | 1331.048718 | 137365.531842 | 1.620427281e-4 | 0.153803548 | 7310 | 106 |
| Vented, SDIRK2 | 0.001 | 1295.158292 | 114614.581605 | 1.481052000e-4 | 0.144886016 | 6909 | 88 |
| Closed, SDIRK2, smaller step | 0.0005 | 1330.605417 | 137322.731028 | 1.619125959e-4 | 0.153737145 | 9629 | 158 |

At the original timestep, final temperature differs from the corresponding
backward-Euler run by +0.000742 K closed and −0.000230 K vented. Accepted
chemistry substeps drop from 79392 to 7310 closed and from 73235 to 6909
vented. A SDIRK2 step contains more implicit stages, so a substep ratio is not
a runtime-speedup measurement. The optimized two-case invocation took
40.5295 s wall time in one observed run; the closed smaller-step invocation
took 19.4146 s, excluding compilation. There was no paired wall-time benchmark
of backward Euler, and these times are not a general performance guarantee.

No network halvings were needed. Across these three SDIRK2 runs, maximum
combined residuals were 4.670e-17 kg mass, 2.655e-10 J energy and 4.647e-14 mol
elements. The vented pressure residual remained below 9.990e-6 Pa, and the
largest surface-power residual was 2.191e-7 W. The same positivity, shared-clock,
source-domain and physical-accounting checks remained active.

The outer timestep sensitivity remains: halving it changes closed gas
temperature by 0.443301 K, residual H2 from 3.2733e-9 to 9.9238e-9 kg, and OH
from 3.9195e-7 to 5.2204e-7 kg. The higher-order reactor method does not resolve the
surface/gas splitting or coarse solid mesh. No new global time-order,
unconditional positivity, discrete entropy theorem or experimental-validation
claim follows from selecting SDIRK2.

## Methane and retained-char modes

The full methane source uses the same prepared constrained reactor, including
accepted-step local reaction coordinates and a checked direct Newton fallback
for cancellation-sensitive trace channels. This changes neither the physical
acceptance tolerances nor the example's first-order coupling. Summary counters
include successful direct recoveries from chemistry calls in accepted network
payloads; they do not count chemistry in discarded pressure-root trials.

The original hydrogen consumer was rerun after this extension at 0.1 s,
0.001 s maximum step and 12 solid cells, using SDIRK2. It reproduced
1331.048718306 K and 137365.531842355 Pa, with 7310 accepted chemical steps,
106 chemical rejections and no direct fallback. Its maximum combined residuals
were 2.10e-17 kg mass, 9.39e-11 J energy and 1.73e-14 mol elements.

The frozen charring-methane control was executed with the same duration,
timestep and mesh in both closed and vented configurations:

| Case | Gas temperature (K) | Pressure (Pa) | Retained char (kg) | Emitted CH4 (kg) | Molecular opacity (m^-1) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Closed, frozen chemistry | 1483.269144014 | 100761.760097715 | 3.687486464256e-6 | 1.474994585702e-5 | 1.166101486895e-2 |
| Vented, frozen chemistry | 1468.963507895 | 101312.977035876 | 3.687413601255e-6 | 1.474965440502e-5 | 1.172167248934e-2 |

Both controls use 100 accepted outer intervals, 100 conservative remaps and no
network halvings. They retain exactly one quarter as much char as emitted
methane, as required by the 20% char yield. With homogeneous chemistry
disabled, no CO, CO2, H2O or NO is generated; original methane still absorbs
and emits in the declared molecular-opacity approximation. Across both
controls the largest residuals were 6.43e-17 kg mass, 3.26e-7 J energy and
1.41e-14 mol elements. The vented pressure residual stayed below 8.25e-6 Pa.
These are manufactured physical-accounting checks, not measured-material or
experimental fire validation.

The reacting vented charring case was separately executed through 0.020 s
with the same 0.001 s maximum step and 12 cells:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- vented 0.02 0.001 12 kinetic sdirk2 charring-methane
```

It completed 20 outer intervals and 20 remaps, with no network halving,
1814 accepted chemical steps, 15 chemical rejections and two successful direct
Newton recoveries. Its gas endpoint is 1522.866245311 K at
102199.510551066 Pa. Retained char is 7.307958338754e-7 kg; emitted methane
is 2.923183335501e-6 kg and net chemical methane consumption is
2.922890229529e-6 kg. Vent transport carries the remaining difference in
addition to the final retained methane.

The retained gas contains 6.385709128275e-6 kg H2O,
7.828288815611e-7 kg CO, 6.749418836361e-6 kg CO2 and
2.062289609032e-9 kg NO. These are finite-time outputs of the supplied source,
not equilibrium replacement or validated emissions. Endpoint molecular opacity
is 1.667583073060e-2 m^-1 and signed net vent mass is
-3.603331934282e-5 kg. Maximum residuals are 4.14e-17 kg mass,
6.75e-8 J energy, 8.65e-15 mol elements, 8.51e-6 Pa pressure and
6.02e-8 W surface power. This check's 20 ms duration is separate from the
100 ms frozen controls and closed reacting cases.

### Closed reacting methane and char endpoints

Both closed reacting cases completed 0.100 s using 0.001 s maximum outer
steps, 12 solid cells and SDIRK2. Each completed 100 accepted outer intervals,
100 conservative remaps and zero network halvings:

```bash
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.001 12 kinetic sdirk2 methane
cargo run --release -p fs-airflow --features compartment-fire \
  --example coupled_compartment_fire -- closed 0.1 0.001 12 kinetic sdirk2 charring-methane
```

| Quantity at 0.100 s | Methane | Charring methane |
| --- | ---: | ---: |
| Gas temperature (K) | 1683.198402902 | 1646.380168650 |
| Pressure (Pa) | 114521.127420776 | 111856.846191578 |
| Solid face temperature (K) | 586.790256637 | 583.112693413 |
| Wall face temperature (K) | 354.223466992 | 349.135884797 |
| Retained char (kg) | 0 | 3.691197223245e-6 |
| Emitted methane (kg) | 1.845524418064e-5 | 1.476478889298e-5 |
| Gas methane (kg) | 7.098127257085e-15 | 2.781055449400e-13 |
| Gas H2O (kg) | 4.102870996967e-5 | 3.282265397031e-5 |
| Gas CO (kg) | 2.429458799602e-7 | 2.579431951179e-7 |
| Gas CO2 (kg) | 5.024451160383e-5 | 4.009733030090e-5 |
| Gas NO (kg) | 2.527904623452e-8 | 2.015138822055e-8 |
| Molecular opacity (m^-1) | 1.010265908974e-1 | 8.359188712794e-2 |
| Solid surface recession (m) | 2.171205197717e-6 | 9.408934098501e-7 |

| Work / maximum residual | Methane | Charring methane |
| --- | ---: | ---: |
| Accepted chemical steps | 12140 | 11384 |
| Rejected chemical steps | 278 | 159 |
| Successful direct Newton fallbacks | 2321 | 98 |
| Combined mass residual (kg) | 5.284130338151e-17 | 7.156750777941e-17 |
| Combined energy residual (J) | 1.979501007554e-10 | 3.259047289816e-7 |
| Element residual (mol) | 8.646729166006e-15 | 1.193966800428e-14 |
| Surface power residual (W) | 7.151322023446e-8 | 5.622507615044e-8 |

The char case retains one quarter of the emitted methane mass as condensed
carbon; the corresponding original atomic inventories close independently.
These single-resolution coupled checks do not establish temporal or spatial
convergence, a flame-profile prediction, or calibrated methane/wood emissions.
The separate 60 s charring-slab mesh/time study is retained in
`data/reference/compartment_fire/pyrolysis/README.md`; it does not substitute
for coupled-gas refinement.
