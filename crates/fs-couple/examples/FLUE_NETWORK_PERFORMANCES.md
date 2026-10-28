# Branched flue resonators and regional atmospheres from data

The `music_render wind` command accepts `flue-network.performance` through the
same `frankensim-flue-performance-v1` source format as the single-pipe example:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-network.performance flue-network.wav
```

This is a synthetic geometry and performance, not a characterized flute. One
jet drives the complete connected acoustic graph. Branches are not separate
voices. Every section derives travel time, characteristic impedance and any
selected viscothermal loss from its own coherent material gas state. Every
selected opening contributes its actual terminal volume flow to the existing
multirate Rayleigh receiver. The source never substitutes internal pressure for
exterior microphone pressure.

## Geometry and material gas

### Molecular gas response inside DG bores and along exterior paths

`flue-dg-molecular.performance` uses the existing DG flue renderer with explicit
terrestrial-air species. The same shared material relaxation law acts inside
the cylindrical bore and on every surface-to-microphone path. Its jet and short
tube remain exploration fixtures, not a characterized flute.

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-dg-molecular.performance molecular.wav --decimate
```

After the complete `dg_tube` or `dg_bore` block (including its radiation and exterior rows),
before any `dg_wall`, select:

```text
dg_molecular iso9613-terrestrial-air
```

The material-modal DG reed path accepts the same row for its bore and optional
upstream tract. Every finite side chimney evaluates the law against its own gas
record. The complete state allowance includes all thermal/molecular states
before propagator or radiation preparation. `dg_bore` retains the original
spatial wall quadrature for cones, flares and other admitted geometry with a
uniform gas. Each axial sample uses the local radius and its own radial thermal
diffusion rates. Axial gas fields still refuse: substituting their mean or inlet
state would discard the prescribed temperature/composition field.

`flue-dg-conical-molecular.performance` exercises this path with a declared
1-to-2 mm radius cone, twelve radial poles and four axial quadrature points.
It retains the same independent internal and exterior gas declarations:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-dg-conical-molecular.performance conical.wav --decimate
```

This is an exploration of the shared local-cylinder Webster model, not a
measured flute or oboe. Axial quadrature, radial-pole and acoustic-mesh errors
remain separate convergence axes; low-frequency storage is exact only to the
admitted spatial quadrature's accuracy. Molecular state counts follow the
actual thermal coordinates: `elements * quadrature_points * poles * reservoirs`
for `dg_bore`, rather than the nodal count used by `dg_tube`.

The exterior choice stays independent:

```text
observer_gas molecular 4096 1000000 50000000 iso9613-terrestrial-air
```

The three integers bound distinct paths, retained source-filter taps and kernel
work. The exterior law binds the original exterior composition. Native plate
observers accept this same declaration against the string performance's gas.
ISO rates require admitted terrestrial-air composition; adding helium or
enriching oxygen/CO2 does not silently reuse those rates. Legacy pseudocomponent
air has no explicit molecular species and must be replaced by a `mixture` or
source-bound `thermochemical` record for this option.

For other gases, replace `iso9613-terrestrial-air` by `declared N`, followed on
the same line by exactly N species records `ID RESERVOIRS`. Each reservoir has
one of the following definitions and one rate definition:

| Definition | Fields |
| --- | --- |
| Vibrational reservoir | `vibration MODES THETA_K DEGENERACY ...` |
| Rotational reservoir | `rotation DOF` (2 or 3) |
| Partner-independent rate | `power P_TAU_REF_PA_S T_REF_K EXPONENT T_MIN_K T_MAX_K` |
| Collision-partner rates | `partners COUNT ID LOG_P_TAU_PA_S SLOPE_K_1_3 T_MIN_K T_MAX_K ...` |

These are the original `fs-material::relaxation` parameters, not an instrument
loss table. Every present species and required collision partner must have data;
`ID 0` explicitly declares no reservoirs. Supplied coefficients remain authored
constitutive data and gain no measured-source authority from parsing. Missing
data, excess counts, invalid temperatures and nonphysical heat capacities refuse.

Inside the bore, molecular and radial thermal states share one thermodynamic
pressure equation and the original DG exponential. They are not independent
loss filters: that would count low-frequency compressibility twice. Shear and
compliant solid-wall mechanics retain their existing owners. The equilibrium
Bessel error admission is not a new accuracy certificate for the combined law.
No extra bulk Stokes–Kirchhoff term, axial heat conduction, finite wall thermal
feedback, evolving chemistry or live reservoir-transfer rule is added. Exterior
propagation keeps its original shared causal material filter and travel phase.

Replace the original `tube LENGTH RADIUS LENGTH_ERROR REFLECTION` row with:

```text
network NODE_COUNT SECTION_COUNT
duct_node KIND ...
... one duct_node per node ...
duct_section NODE_A NODE_B LENGTH_M RADIUS_M LENGTH_ERROR_M [gas GAS_RECORD]
... one duct_section per section ...
```

Node indices and section indices are zero-based in source order. The shared
network owner admits one degree-one inlet, connected sections, valid junction
and load degrees, and passive loads. The file limits are 64 nodes and 128
sections. All dimensions and numerical error allowances are explicit.

| Node declaration | Physical meaning |
| --- | --- |
| `duct_node inlet` | The unique port driven by the jet |
| `duct_node junction` | Common-pressure junction with no added storage |
| `duct_node driven` | Degree-one chimney with an explicitly declared material moving boundary |
| `duct_node reflection R` | Explicit memoryless pressure reflection in [-1,1] |
| `duct_node impedance R L C` | Degree-one series resistance, inertance and compliance |
| `duct_node series R L C` | Two-section inline pressure-drop load |
| `duct_node shunt R L C` | Load from an interior common-pressure junction to its reference |

For R-L-C records, units are Pa s/m³, Pa s²/m³ and m³/Pa. Zero `C` explicitly
omits compliance. Positive compliance and nonnegative resistance/inertance are
checked by the existing passive impedance owner; the file parser contains no
replacement impedance or propagation equations. Supplied coefficients are not
promoted to measured radiation, tissue or structural material properties.

A section without `gas` uses the top-level gas record. A supplied record uses
the shared material gas parser, including explicit species mixtures and loaded
thermochemical sources. All regions must have the same static pressure. The
jet and a selected windway use the **actual inlet section's** density,
viscosity, sound speed and other gas properties, even when that section is not
first in the file. `exterior same` uses the top-level default gas; an explicit
exterior record describes the receiver environment independently.

The example includes a warmer N₂/O₂/CO₂ inlet, an ambient main tube and a
separate branch temperature. These are prescribed frozen states. They do not
claim to solve breath advection, thermal exchange, gas mixing or condensation.

## Material walls in the original graph

After the complete duct records, an optional wall block assigns source material
to an unused degree-two junction between equal-radius sections:

```text
walls COUNT MIN_HZ MAX_HZ
wall NODE isotropic PACK T_K OMEGA_RAD_S THICKNESS_M AXIAL LOSS
wall NODE orthotropic PACK T_K OMEGA_RAD_S THICKNESS_M AXIAL HOOP_AXIS AXIAL_AXIS LOSS
```

Supply exactly `COUNT` wall rows. `AXIAL` is `stress-free` or `strain-fixed`.
Isotropic loss accepts `elastic`, `kelvin-voigt`, or the existing sourced Maxwell
spectrum. Orthotropic loss accepts `elastic`,
`proportional-kelvin-voigt TIME_PROPERTY`, or
`normal-kelvin-voigt MAX_DT_OVER_TAU MAX_ANGULAR_STEP HH HA AA`.
The proportional property is a nonnegative time in seconds declaring
`C_viscous = tau C_elastic` in the selected membrane plane.
Its distinct principal axes are numbered 1 through 3 and follow the **local
hoop and axial directions**. They do not describe a fixed Cartesian grain angle.
The original directional moduli, major Poisson ratio, density and time claim
retain their source identities; no equivalent isotropic material is fabricated.

For independent normal viscosity, `HH`, `HA`, and `AA` are the exact names of
Pa·s properties in the local hoop/axial plane. Their symmetric matrix must be
positive semidefinite; the cross term may be negative. Fixed axial strain selects
the hoop response. Free axial stress eliminates the complete `C + sD` law,
retaining its resulting Maxwell memory together with residual Kelvin resistance
in the original pressure/flow load. Separately condensing elasticity and viscosity
would lose that physical memory. The declared `MAX_DT_OVER_TAU` is in `(0, 2]`
and `MAX_ANGULAR_STEP` in `(0, 1]`; the graph's actual mechanics clock must satisfy
both limits when memory is present. Native provenance retains the tensor source
names and values, material axes, reduced memory, and residual wall resistance.

This is the same local constitutive reduction as the DG and material-reed paths.
It changes the existing pressure/volume-flow load, wall storage and physical
loss while retaining the jet, regional gas, propagation and microphone owners.
Each wall covers half of each adjacent physical section; unselected strips stay
rigid. Source claims must cover the wall band, which covers the declared acoustic
band. The existing thin-wall and topology limits remain: off-axis grain, full
three-dimensional viscoelastic condensation, general directional relaxation,
axial shell waves and exterior shell radiation are not
supplied by this local membrane image. Omission or `walls none` retains rigid walls.

## Distributed loss and exterior openings

Before `mouth_radiation`, optional `mouth_aperture rectangle WIDTH_M LENGTH_M`
sets independent acoustic aperture dimensions. Both the receiver and a selected
mouth reaction use this rectangle; the jet's width, height, cut-up, source
coefficient and convection remain unchanged. Omission uses flue width times
cut-up, as before. A supplied aperture requires `mouth_radiation`; nonpositive
or nonfinite dimensions refuse. This row is supported by the network route only;
the DG extension remains separate. Render metadata reports the resolved dimensions.

The flute cards use 10 by 12 mm from Boehm's mouth-hole dimensions in
`data/instruments/geometry/flute-boehm-silver-bfoot-boehm1871/dimensions.tsv`.
This is a rectangular approximation to the published rounded hole, with unknown
corner radius and lip coverage. The roughly 90 mm² PF-661 measurement belongs
to another specimen; the rectangle's 120 mm² area is not that measurement.
These rows correct the former use of the much smaller jet span as the radiating
hole. They do not enable mouth reaction or change the authored linear loss.

A network uses `viscothermal none`, or explicitly addresses its original
physical sections:

```text
viscothermal sections 2
loss_section 0 CELLS MINIMUM_FREQUENCY_HZ MAXIMUM_FREQUENCY_HZ
loss_section 2 CELLS MINIMUM_FREQUENCY_HZ MAXIMUM_FREQUENCY_HZ
```

Each addressed section is split into 4 x CELLS quarters, and every quarter must
span at least one acoustic transit sample: `length / (4 x CELLS) >= c / f_mechanics`
(3.6 mm for air at 96 kHz, 7.3 mm at 48 kHz). A finer split refuses with
"fractional section needs finite geometry and at least one full transit sample".
Each cell must also stay within 0.5 rad of phase at MAXIMUM_FREQUENCY_HZ. A short
section at a low clock can admit no cell count for a wide band (0.06 m at 48 kHz
with a 1 kHz band needs 3 cells for phase but allows only 2 for transit), so lower
the loss band or raise the mechanics clock.
The original section indices survive the shared loss subdivision. No gas is
averaged across branches. Existing gas, wide-tube, fractional phase, memory and
frequency-domain admission checks remain in their physics owners. A plain
single tube retains its existing short viscothermal record.

The material moving-boundary compiler is shared with the
[two-sided material-modal wind source](material_modal_valve/README.md#material-moving-keys-on-either-acoustic-side).
The same `boundaries`, `boundary_plate`, contact, optional viscous slit and
Maxwell records compile one original cantilever from its material pack and
actual chimney gas/load. The flute keeps its existing typed `closing_force`
gesture paths and source-clock ramp compiler. The two-sided reed format uses
that same physical source grammar with typed source/key gestures or explicit
mechanical-sample force and exterior-reservoir programs on either side.

After `receiver`, an explicit network must select its openings:

```text
outlets PRIMARY_NODE ADDITIONAL_COUNT
outlet OTHER_NODE CENTER_X_M CENTER_Y_M
... one outlet per additional opening ...
```

The primary opening defines the origin of the common baffle plane. Additional
centers and receiver coordinates use that same frame. The shared receiver
checks terminal topology, opening radii, overlaps, common pressure and its
quadrature limits. The example's two outlets contribute to one pressure
history, preserving propagation timing and interference. A single-tube source
may omit this block and retain its original sole terminal receiver.

Prescribed terminal loads retain their one-way infinite-baffle receiver. The
explicit MIMO selection below instead supplies shared exterior feedback. These
remain plane-wave cylindrical sections and ideal junctions; the source file
does not establish measured flute fidelity, a resolved jet or a validated
three-dimensional radiation model.

## Explicit mutual radiation from a supplied artifact

`flue-mimo.performance` and its adjacent `flue-mutual.mimo-radiation` exercise
the existing common terminal operator through the native loader:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-mimo.performance flue-mimo.wav
```

The artifact is an explicitly **synthetic passive two-port example**. Its
coefficients demonstrate the complete source-to-feedback-to-microphone route;
they are not a measured flute or a validated exterior-radiation bake.

After optional wall and moving-boundary records, select:

```text
radiation mimo flue-mutual.mimo-radiation 32 1048576 257
```

The fields are the relative artifact path, maximum rational states, maximum
runtime memory bytes and sampled-passivity grid points. The common MIMO decoder
also limits artifact input to 8 MiB. The original runtime enforces the declared
state/memory limits and checks both continuous and actual sampled responses.
The sampled result remains `GridOnly`; it is not a discrete storage certificate.

Artifact ports follow the exact `outlets`/`outlet` order. Their circle centers
and radii must match those openings. The artifact separately declares its
exterior gas using the same gas-record grammar as the performance, and its
frequency band must cover the entire mechanical Nyquist band, even when the
receiver emits a narrower audio band. Existing reactive terminal loads cannot
be replaced. No matrix is inferred from the instrument name, fitted at load
time or rescaled to conceal a stale gas declaration.

Each accepted terminal flow drives the one shared rational state, whose full
matrix pressure returns to the same duct and jet solve. The same openings and
exterior gas drive the multirate microphone. Changing off-diagonal coefficients
therefore changes the jet's returning wave as well as emitted pressure. Material
walls, regional section gases, windway dynamics, mouth contraction, moving-key
mechanics and the gesture clock remain attached to their original owners.

An atmosphere change that affects the baked operator requires a matching
artifact. Changing only the performance gas while retaining an incompatible
artifact refuses before rendering. The library entry point
`FlueFilePerformance::from_bytes_with_all_sources` accepts the existing material
loader plus a radiation-byte loader; the CLI resolves both relative to the
performance directory. `mimo_radiation_record()` retains the exact artifact
hash, declared identities, gas record, geometry, band and limits. The existing
constructors and `prescribed`/`compact-baffled` selections retain their behavior.

To derive that matrix from the actual aperture geometry and exterior instead,
use `radiation computed-baffled` with the fourteen explicit numerical controls
in [the shared source guide](MIMO_RADIATION_ARTIFACTS.md).
`flue-computed-radiation.performance` is a complete 48 kHz example. Changes to
the gas or openings rerun the existing physical array and fitting owners;
no radiation artifact or edited coefficient table is needed. The same terminal
host supplies reaction and observation, and `computed_radiation_record()`
retains the actual geometry/gas, attempted orders and accuracy/passivity results.
The finite-band baffled approximation must pass every requested gate; it does
not imply measured flute accuracy or complete delayed-radiation qualification.

## Pressure and jet direction in one performance

`flue-embouchure.performance` adds a prescribed jet-angle phrase, the
existing separated-flow mouth loss, and physical mouth pressure to the
branched, regional-gas example:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-embouchure.performance flue-embouchure.wav
```

The source is a manufactured fixed-window demonstration. Its geometry, player
controls and loss coefficient are authored inputs; none was inferred by
matching an instrument recording.

The canonical gesture suffix requires one `blowing-pressure` track and accepts
one additional `jet-angle` track, whose values are `angle-rad`. Track IDs are
arbitrary and distinct. Both tracks use the original gesture sampler, including
interrupted ramps. A tick applies at
`ceil(tick * mechanics_rate / control_rate)`, before that mechanical sample.
Changing callback size does not alter the physical control clock or discard
partial microphone frames. The `controls` limit bounds the number of tracks
times control ticks before held-value coalescing.

The direction enters the existing geometric cut,
`labium_offset + cut_up * tan(jet_angle)`. It changes the jet volume entering
the mouth, and the same source-pressure derivative and implicit mouth solve
propagate that change through every duct and selected outlet. Angles must be
finite and strictly inside `(-pi/2, pi/2)`, with a representable transverse
offset. The complete phrase is checked before any source sample.

This is the existing instantaneous tilted-cut approximation at fixed mouth
area, cut-up, exit width and exit height. It does not resolve transient jet
turning or lip movement, change lip coverage or remesh convection. The
windway velocity field, jet feedback history, wave propagation, material-wall
memory and radiation history persist. A pressure release still leaves stored
fluid and acoustic motion intact. Dynamic cut-up or mouth-area changes need a
moving-domain and boundary-work model and are not exposed as simple controls.

## Mouth and outlet sound in the same exterior

An optional record after the outlet records and before optional
`control_source` and `controls` adds the mouth to the physical microphone:

```text
mouth_radiation -0.3 0 4 4
```

The first two numbers place its center in the common baffle plane; the last
two select quadrature counts along width and cut-up. The rectangle's actual
dimensions come from `flue`, so changing the source width or cut-up changes
its acoustic area and directivity automatically. Overlap with an observed
terminal opening refuses construction.

The receiver consumes the actual signed outward inlet-node volume flow,
which is the negative of flow into the bore. That rectangular Rayleigh source
and the terminal sources share one exterior gas and physical clock. Their
pressures add with signs, travel delays and every accepted mechanical sample
preserved, allowing cancellation as well as reinforcement. The summed signal
uses the declared pascal scale; it is never normalized to a target loudness.

This is explicitly a coplanar infinite-baffle description, not a free-field
flute model. Mouth observation adds no second jet force, turbulent noise,
mutual radiation or new reaction load. Existing mouth resistance and duct
termination laws still determine the simulated motion. Omitting the record
retains the original terminal-only output.

## Compact reaction from the actual mouth and exterior gas

`flue-mouth-reaction.performance` connects the rectangular mouth's radiation
pressure to the source mechanics. It retains the original network, regional
gases, unsteady windway, jet excitation, separation and receiver. The mouth's
physical self reaction replaces the equivalent-circle end correction and
authored linear mouth loss. After the mouth observation record, select:

```text
flue 0.010 0.001 0.004 0 0 0
... existing windway, jet, network and receiver records ...
mouth_radiation -0.3 0 8 4
mouth_reaction compact-baffled
```

The exact static Rayleigh integral of the actual width by cut-up determines
near-field inertance `M`. The exterior gas determines both `M` and the exact
compact resistance coefficient `C = density / (2*pi*sound_speed)`. Those two
leading physical terms determine the passive impedance `R*s/(s+p)`, with
`p=M/C` and `R=M*p`; no measured note, equivalent mouth radius, fitted gain or
authored damping multiplier enters that load. The existing `fs-vfit` midpoint
relaxation owner supplies its state, near-field storage and escaped energy.
Its current-flow pressure and retained history enter the original simultaneous
jet/mouth/duct solve, including sub-sample jet feedback and optional separation.
No extra source is added to the output.

Changing just `exterior` to a supported gas mixture changes the physical mouth
motion and returning bore waves, as well as outgoing propagation. For example:

```text
exterior mixture 293.15 101325 mole 2 N2 0.5 CO2 0.5
```

Reaction and observation retain the same rectangular cells. Cold construction
allows at most 1024 cells, requires `k*mouth_diagonal <= 0.5` and
`dt*fmax <= 0.1`, and compares the finite-cell static inertia with the exact
rectangular integral. Coarse quadrature exceeding 5% refuses and requests
refinement. At 32 frequencies across the declared band, both analogue and
actual bilinear responses must also agree with the shared Rayleigh surface
within 5%, independently for complex impedance and weak radiation resistance.
Those finite checks retain finite-frequency quadrature error and do not
establish broadband transient accuracy or a measured flute claim. The band
is not an output low-pass filter. Mouth-to-terminal mutual loading, free-field
body scattering and resolved jet turbulence remain outside this selection.
The existing outlet MIMO option still covers its declared terminal ports only.

`mouth_reaction none` and omission preserve the original mouth mechanics.
Enabling reaction requires the actual `mouth_radiation` receiver, zero legacy
`mouth_loss_rel`, one cold attachment and matching static pressure. An explicit
`mouth_loss contraction ...` remains a separate physical separated-flow loss.
The block ledger records returned/stored near-field energy separately from
nonnegative radiation loss; the replaced backward-Euler mouth inertance no
longer contributes numerical dissipation. CLI metadata records the actual gas,
passive coefficients, checked band and both response discrepancies.

## Mouth separation from the actual geometry and gas

After the optional `windway` row and before `jet_card`, select:

```text
mouth_loss contraction 0.85 0.3 authored-fixed-window-demonstrator
```

The fields are contracted jet area / physical opening area, maximum contracted
Mach number, and one explicit source token. The existing `ContractionSpec`
owner validates them and binds the loss to the physical mouth window
`flue_width * cut_up` and the **actual inlet section's complete gas state**.
There is no independently entered density or sound speed. The example's source
label explicitly identifies a manufactured coefficient, not measured edge data.

The mouth solve includes the signed pressure drop
`rho * Q * abs(Q) / (2 * (Cc * area)^2)` at the accepted flow, together with
mouth inertance, the supplied linear mouth resistance and the same acoustic
port. Its power is recorded as separation loss; it is distinct from radiation
and windway viscosity. The solved contracted Mach limit remains active during
rendering. `Cc = 1` still supplies this loss law. Omit the row or use
`mouth_loss none` for the original linear mouth.

The CLI sidecar includes the supplied coefficient/source, bound effective area
and actual gas density/sound speed. This remains a quasi-steady separation
closure with no vortex memory, resolved shedding or measured flute accuracy.
It does not qualify the full Boehm assembly or the plan's register experiment.

## Causal DG bore through the original jet and wind renderer

`flue-dg.performance` selects the explicit
`frankensim-flue-dg-performance-v1` format. The original jet, mouth equation,
material-gas windway and gesture compiler then use the same prepared DG port
as moving reeds. The exact propagator exposes its current step-mean pressure
law; the original mouth solve supplies the same accepted flow to that owner.
Every acoustic and viscothermal memory state survives callback boundaries.

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-dg.performance flue-dg.wav --decimate
```

The fixture declares a 96 kHz mechanical clock. The command's existing
`--decimate` path produces 48 kHz PCM with its reported filter delay. Its
short 1 mm-radius bore is an integration fixture, not a flute specimen.

Source records, in order, are `clock RATE`,
`output SAMPLES FULL_SCALE_PA MAX_ABS_PRESSURE_PA`, `gas GAS_RECORD`,
the shared `dg_tube` or `dg_bore` declaration and its geometry/gas-field rows,
`flue WIDTH HEIGHT CUT_UP OFFSET ANGLE MOUTH_LOSS_REL`,
`windway none|slit LENGTH|unsteady LENGTH TRANSVERSE_CELLS CONVECTION_CELLS`,
`jet_card tonal-interim|json CANONICAL_CARD`,
`jet DISPLACEMENT_GAIN SEED_AMPLITUDE VOICE_SEED VOICE_INDEX`, an optional
`jet_response forced-bickley SPACING WINDOW`,
`jet_limits MIN_SPEED MAX_MACH MAX_HISTORY`, `controls MAX_CONTROLS MAX_WORK`,
and `observation inlet-pressure` or `observation baffled X Y Z RINGS ANGLES FMAX_HZ`,
followed by the existing canonical `gesture` suffix. Pressure and optional jet-angle tracks retain their existing
physical units and interpolation. The actual DG inlet supplies mouth bore area,
density, viscosity and sound speed; there is no second authored inlet medium.

`DISPLACEMENT_GAIN` is the default jet chain's AUTHORED gain. The cards' 90 is
about 7-18x the forced jet's response at each note's speaking mode (mail 7491).
It stands until the forced jet passes fg_003, fg_004, fj_001 and gas_direction_amount
(ruling 7518, mail 7859, bead pexfp.135.1). The optional `jet_response forced-bickley SPACING
WINDOW` row, in this file and in the network files after the `jet` row, selects the
forced Bickley jet's causal response by name in place of that chain. It is
integrated at grid SPACING [b] over WINDOW [b/U0] at the Reynolds number of the
gesture's highest blowing pressure, and register physics is graded on 0.1 90. The
gain is then unused. In network files it excludes `jet_transfer nolle-price-2015`.

During a pressure gesture the response uses the accumulated convective clock
`tau = integral U0(t)/b dt`. Each accepted physical-time velocity interval keeps
its own speed and age; changing the pressure no longer retimes the entire past
using the newest speed. The current velocity remains an affine unknown in the
same mouth solve. Speed is held and transverse velocity is linear within each
mechanical interval. A constant positive-speed trajectory retains the original
hat-weight arithmetic. Zero speed disables forcing/output and freezes this
response clock without consuming history capacity.

This is an explicit **fixed-Reynolds response approximation**,
`dz/dt = (U0(t)/b) A z + B v`, using the already integrated operator. It does not
resolve changing Reynolds number, molecular decay at rest, jet spreading or the
labium's hydrodynamic reaction. The finite response window and existing history
budget remain binding; no register or acceptance gate is relaxed.

`MIN_SPEED` in `jet_limits` (this file and the network files) is the jet's declared
validity edge U_min, not a stability limit (ruling 8132, bead pexfp.135.2). It sizes
the retained convection history, so a smooth gesture may pass below it, at any control
rate. In the onset, a lookup older than the first motion reads the rest state, which is
exactly zero. Any other lookup below U_min that leaves the history (a release, or a dip)
gives no jet deflection, and the windway's steady flow splits at zero deflection. The
forced response takes the same state, its clock frozen as at rest, once its retained
intervals no longer reach its window. The receipt names the edge.

The bore declaration is shared with plate valves; see
[PLATE_VALVE_PERFORMANCES.md](PLATE_VALVE_PERFORMANCES.md) for its numerical
limits and geometry/gas-field syntax. Static gas/material substitutions rebuild
the source and causal bore together. They do not mutate an accepted history.

The default observation emits **internal step-mean inlet pressure**. The explicit
computed-baffled selection below instead emits exterior pressure from its solved
terminal flow. Neither selection adds solid-wall motion or a room.
The retained tonal jet card and authored receptivity remain exploration-only;
`--acceptance` refuses this source. This connection does not implement the
resolved compressible jet or zonal interface required by the instrument plan.
Dense DG dissipation remains an inferred energy deficit where the original
ledger reports `balance_checked = false`. Numerical spatial/time convergence
and measured instrument fidelity require separate evidence.

### Computed DG radiation and exterior sound

`flue-dg-baffled.performance` adds a causal radiation reaction and microphone
to the same jet, windway, bore, material gas and gesture owners. Select `baffled`
in the shared DG header. After every geometry and gas-station row, provide:

```text
dg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257
dg_exterior mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2
```

The fourteen controls use the shared source grammar documented in
[PLATE_VALVE_PERFORMANCES.md](PLATE_VALVE_PERFORMANCES.md#computed-dg-outlet-radiation-and-exterior-pressure).
The actual outlet area supplies the circle; the complete exterior gas record
supplies its radiation medium and retains thermochemical source receipts.
Bore and exterior static pressure must agree. Temperature and composition may
differ. A material change rebuilds the physical source and causal load.

At the existing observation position, select:

```text
observation baffled 0 0 0.2 8 32 2000
observer_gas stokes-kirchhoff 4096 1000000 50000000
```

These are microphone X/Y/Z in metres, radial/angular quadrature counts and
maximum declared frequency in Hz. Every accepted outward outlet-flow sample
enters the original differentiating retarded Rayleigh observer once. It uses
the retained source's exact opening and exterior, without adding a second
reaction. Inlet observation remains available with the same reaction active.

The optional `observer_gas` record follows `observation` and precedes any
`mouth_radiation` row. It selects `stokes-kirchhoff` or an explicit `lossless`
comparison, followed by positive limits for distinct exact path lengths,
retained source-filter taps and kernel construction terms. Their shared
ceilings are 16,384, 4,000,000 and 500,000,000. Omission preserves the previous
lossless observation. An internal-pressure observation cannot select an
exterior propagation law.

Classical propagation uses the complete bound `dg_exterior` state: density,
sound speed, viscosity, conductivity and heat capacity. A change of gas
composition, pressure or temperature therefore rebuilds the physical transfer
of every original surface path. It applies to all declared circular terminals
and to the optional rectangular flute mouth before their signed pressures are
summed. Physical flight time and gas phase remain separate from one common,
reported numerical latency; no post-mix equalizer or extra reaction is added.
Metadata retains the applied law, actual gas, path work and source-filter tail
certificate bounds; the library report keeps every original path receipt.
Molecular relaxation is explicitly reported as omitted by this
classical law; the native record does not invent species relaxation spectra or
apply a terrestrial-air law to an arbitrary gas mixture. Library callers may
select the existing validated ISO or sourced molecular law through the same
`GasPropagation` owner. This selection does not add bulk molecular relaxation
inside the bore.

This observer runs on the mechanical clock within its original 8–192 kHz rate
domain. The existing `--decimate` path remains available for explicit 96-to-48
kHz conversion. Metadata includes actual receiver gas, propagation support,
uncompensated filter delay and original source fitting/passivity reports.

The example requests a 10% fit budget for integration coverage; its numerical
admission has not been executed in the authoring environment. Original fitting,
held-out, causality, storage and passivity failures remain errors. Reaction
covers mechanical Nyquist independently of the microphone band. The observer
band must also fit the wall-loss admission and original quadrature limits.

This is one circular opening in an infinite rigid baffle. It does not add
unflanged-bell diffraction, mouth-to-outlet mutual radiation, body scattering
or resolved jet turbulence. The declared observer band is not an additional
output low-pass filter. No measured flute qualification is implied.
