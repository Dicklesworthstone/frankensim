# Mesh and material driven valve performances

A complete supplied plate, spatial lay, material history, air column and pressure
phrase can now produce a WAV without editing Rust or replacing the model with the
primitive reed fixture. The source header selects its physical owner explicitly.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve.performance /tmp/plate-valve.wav \
  --decimate --block 37
```

The committed file is an **illustrative numeric specimen, not measured cane**:
21 original nodes, 24 triangles, a supported root, asymmetric 10–390 micrometre
clearances, area-weighted lay contact, two supplied Maxwell arms and a pressure
phrase with interrupted ramps, releases and a second attack. Material, contact
and traveling-wave history survive the complete phrase. No attack envelope,
pitch correction, guessed mass/stiffness, output normalization or new integrator
is added. Pressure remains an **internal tube observation, not an exterior
microphone or a validated instrument recording**.

## Complete source records and SI units

The example supplies every record in this order. Header records are whitespace
separated; blank/unknown/trailing records and extra fields refuse. Counts are
checked before allocation. Node, triangle and section indices are zero based in
source order. No coordinate welding, geometry repair or region inference occurs.

```text
frankensim-plate-valve-performance-v1
audio MECHANICS_HZ SOURCE_SAMPLES FULL_SCALE_PA
ambient TEMPERATURE_K PRESSURE_PA RELATIVE_HUMIDITY
tube LENGTH_M RADIUS_M TERMINAL_REFLECTION LENGTH_ERROR_M MAX_WAVE_BYTES
observation inlet|terminal
plate clamped|simply-supported PRETENSION_N_PER_M DAMPING_RATIO LOW_HZ HIGH_HZ MODE_INDEX MAX_ANGULAR_STEP
aperture REST_COORDINATE_M MAX_SLIT_VARIATION MAX_SLOPE
# Optional: slit_length STREAMWISE_LENGTH_M
initial OPENING_COORDINATE_M VELOCITY_M_PER_S
contact K_PA_PER_M_ALPHA ALPHA CHI_S_PER_M MAX_PENETRATION_M
source SINGLE_TOKEN_SOURCE_LABEL
sections COUNT
section isotropic THICKNESS_M DENSITY_KG_PER_M3 YOUNG_PA POISSON_RATIO fixture LABEL
# Alternatively: section orthotropic THICKNESS_M DENSITY_KG_PER_M3 E1_PA E2_PA NU12 G12_PA ANGLE_RAD fixture LABEL
# Or pack-backed: section material isotropic THICKNESS_M PACK STATE [AXIS=VALUE ...]
#             or: section material orthotropic12|orthotropic21 THICKNESS_M ANGLE_RAD PACK STATE [AXIS=VALUE ...]
nodes COUNT
node X_M Y_M REST_GAP_M SUPPORT_FLAG
triangles COUNT
triangle NODE_A NODE_B NODE_C SECTION_INDEX LAY_FLAG
slit_edges COUNT
edge NODE_A NODE_B
relaxation REGION_COUNT
# Material records below only when REGION_COUNT > 0.
time_limits MAX_DT_OVER_TAU MAX_MATERIAL_ANGULAR_STEP
memory_initial relaxed|unrelaxed
# Alternatively: memory_initial explicit COUNT VISCOUS_DISPLACEMENT_M ...
region SECTION_INDEX E_INF_PA POISSON_RATIO BAND_LOW_HZ BAND_HIGH_HZ BRANCH_COUNT
branch ADDITIONAL_MODULUS_PA TAU_S
compile_limits MAX_SAMPLE_VISITS MAX_CONTROL_TICKS
schedule
```

Repeat section/node/triangle/edge/region/branch records according to their counts.
The `#` lines above explain alternatives; they are **not literal file records**.
Flags are exactly `0` or `1`. A supported node uses the selected plate support law;
a lay-covered triangle contributes its real area to existing contact quadrature.
The supplied gaps determine both local contact and positive-gap slit flow. The
opening coordinate is the retained structural coordinate, not necessarily the
physical mean clearance of an asymmetric slit. Contact provenance retains the
source label as supplied numeric data; it does not mint material-identification
or measured-data authority.

Every section names where its material comes from (bead pexfp.32). Bare numbers
must end in an explicit `fixture LABEL`; without one the file refuses. The
committed examples' blade (E 4 GPa, rho 900 kg/m^3, nu 0.3) is labelled
`fixture authored-illustrative-blade-not-measured`: authored, not cane. A
`section material` record instead names a compiled material-card pack, relative to
this file, the state its card was compiled at (`MaterialStateId` process: a cane
reed must say `water-soaked` or `dry`), and the query point in the card's own
normalized units. Its density, Young's modulus and Poisson ratio (or
`young_modulus_1`, `young_modulus_2`, `shear_modulus_12` and `poisson_ratio_12`)
resolve through `fs_material::state_point`. A card compiled at another state, or a
missing property, refuses by name. The render sidecar's `section_sources` lists
each section's pack, hash, state and point, or its fixture label.

To include viscous flow through the reed channel, insert `slit_length 0.001`
immediately after the `aperture` record, replacing that illustrative one-millimetre
length with the specimen's measured streamwise channel length. The solver uses
the actual local positive gap along each slit edge. The gap resistance therefore
changes as the blade moves or contacts its lay. It solves viscous pressure loss
and separated outflow together with the same-sample bore feedback.

Viscosity and density come from the same inlet gas. A regional gas declaration
on the inlet section takes precedence over the ambient default, and a separate
exterior gas changes radiation without changing the channel medium. Temperature
and mixture changes consequently update channel flow as well as acoustic
propagation. The optional record also works with `symmetric-opposed-pair`; the
declared gap is the total opening between the blades and is counted once.

The WAV sidecar records the selected `slit_flow` law, channel length and resolved
viscosity. Omitting the record retains the existing inviscid channel. This is a
quasi-steady Newtonian slit with separated outflow, not resolved unsteady channel
flow, pressure recovery or validation of an oboe/clarinet specimen. The supplied
solid sections and their constitutive memory remain independent material inputs.

Immediately after `schedule`, append the existing exact canonical bytes of one
`GestureSchedule` blowing-pressure track. The example uses 700 Hz controls on
96 kHz mechanics. The existing compiler applies tick k at
`ceil(k * mechanics_rate / control_rate)` and preserves interrupted ramps.
Controls commit only with successful physical samples. A command starting after
the last observed tick refuses; a ramp may continue beyond the finite window.
The trace never invents the unobserved portion. Body-flow excitation is zero.

Every section must be used. With relaxation enabled, exactly one region names
each section; its equilibrium isotropic bending law must match the actual
covered elements. Positive arms are retained, never truncated. Separate modal
damping must be zero. Orthotropic elastic sections are supported, but a scalar
isotropic Maxwell law is not silently converted into an anisotropic memory law.
The full source plate and all regional inputs remain available through the
runtime. `relaxation 0` selects the ordinary elastic/damping path and omits the
material-only records. Explicit initial viscous displacements use positive-arm
order across regions, in metres of the retained opening coordinate.

The selected mode, rather than the search-window upper endpoint, must satisfy
`dt * sqrt(K/M) <= MAX_ANGULAR_STEP`, with the allowance in (0,1]. Material
attachment also checks its instantaneous frequency, supplied use bands and
`dt/tau` limits. Those are model-resolution checks, not proof that all nonlinear
transient content lies below an observation bandwidth.

The input cap is 4 MiB, 512 nodes, 2048 triangles, 64 sections, 64 total supplied
Maxwell terms, 16384 gesture events, 262144 sampled control ticks and 16777216
sampling visits. Mechanics must be positive and no faster than 192 kHz; duration
is positive, at most 600 seconds and at most 28800000 mechanical samples. Tube
memory has an explicit ceiling no greater than 64 MiB. The existing plate,
contact, material and propagation owners can still refuse within these caps.

## Observation, mixing and retained state

The file's gas state determines fluid density, sound speed and tube impedance;
it does not change the separately supplied solid material constants. The legacy `tube`
interior is lossless characteristic propagation with a supplied passive
memoryless terminal. Its integer transit must satisfy the explicit length-error
allowance; requested and represented lengths are recorded in the WAV sidecar.
This is not the primitive reed fixture's Bessel/TMM bore image, nor a viscothermal
or geometry-derived exterior-radiation claim.

At 96 kHz the command requires explicit `--decimate` for 48 kHz WAV output. Only
observed pressure is filtered: mechanics, contact, memory and pressure controls
stay at 96 kHz. Filter delay and the complete declared window are retained; there
is no tail padding or peak normalization. Native 48 kHz files need no conversion.
Noninteger ratios and incomplete output intervals refuse before output creation.
Existing output and sidecar files are never overwritten.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/valve-and-plate.wav --full-scale-pa 20 --decimate \
  --valve crates/fs-couple/examples/plate-valve.performance \
  --plate crates/fs-couple/examples/plate-mesh.performance --block 37
```

This demonstrates pressure-stream composition, **not physical equivalence of an
internal valve pressure and a plate microphone**. The caller must establish a
compatible observation and common time origin. All parts need equal physical
duration. Source PCM scales are recorded but never become per-part gains.

Library consumers can construct `AperturePerformance` directly around an
`ApertureTube`, `ApertureNetwork` or `ApertureDg`, including already-attached passive loads and
material-bound specimens. The observation explicitly selects inlet, tube terminal
or a network node. `PlateValvePerformance::into_renderer()` returns this complete
finite object, suitable for the existing `DecimatedRenderer`, `PressureEnsemble`
and incremental PCM interfaces. Keep that whole chain during cancellation/resume.
A refused individual mechanical step retains its previous state; a failed output
callback poisons the pressure adapter because an earlier prefix may have advanced.
The accepted prefix and pending controls remain inspectable, never reported as a
complete successful callback.

This remains one linear structural mode with uniform face pressure, a fixed lay
and quasisteady gap flow. Supplied constitutive data is not wet-cane
identification, mesh conformity is the caller's responsibility, and there is no
new continuum-convergence, calibrated-microphone or real-time-throughput claim.

## Exact DG bore with causal viscothermal memory

`plate-valve-dg.performance` selects the original `fs-duct::DgDuct` propagation
owner inside this same file-driven valve and renderer. After `ambient`, replace
the legacy tube/network record with:

```text
dg_tube LENGTH_M RADIUS_M OUTLET ELEMENTS ORDER POLES F_MIN_HZ F_MAX_HZ MAX_REL_ERROR EXP_TOL MAX_DOFS
observation inlet
```

`OUTLET` is `open`, `closed` or `matched`. These are ideal acoustic boundary
conditions; none implies an exterior radiation model. Geometry is one uniform
circular cylinder in the declared ambient gas. Mesh order, element count, pole
count, checked frequency band, relative Bessel error, exponential tolerance and
complete state budget are explicit source data. The state budget includes both
acoustic fields and every retained viscous/thermal memory coordinate, and is
checked before constructing the dense exponential.

The gas record is resolved once through the existing material owner. Its
density, sound speed, viscosity, thermal conductivity and heat capacity drive
the DG column and its Bessel-derived memories. The same inlet gas binds the
original aperture density and optional slit viscosity. Pressure or composition
changes therefore reach propagation and the nonlinear valve through data.
The original authored plate, spatial contact, Maxwell material history,
pressure gestures and physical force events stay in the existing runtime.

At each step the actual DG exponential supplies its mean-pressure relation.
The retained strict aperture solver resolves the corresponding total bore
flow, including swept volume, and both candidate states are accepted together.
No fitted reflection, midpoint relaxation update or second nonlinear pressure
solver is introduced. These files require the existing ExactZoh blade path.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-dg.performance \
  /tmp/plate-valve-dg.wav --decimate --block 37
```

The same input is accepted by `ensemble --valve`. Its output is **internal
step-mean inlet pressure**, not a radiated microphone. Terminal/exterior
observation, independent exterior gas and two-sided material-modal network
conversion refuse rather than silently substituting another model. The WAV
sidecar reports the actual operator, mesh, complete DOFs, pole counts, band,
checked Bessel errors and exact mean impedance.

The example is illustrative source data, not a measured reed or instrument.
Finite-band Bessel admission does not establish DG spatial convergence or
nonlinear time convergence; those refinements remain separate. The dense DG
ledger infers wall/terminal loss from the energy deficit and does not provide
an independent heat integral. This bounded offline route retains a single structural plate/reed mode.
A paired modal valve and exterior radiation need their own coupled paths. Its solid sections and Maxwell
spectra remain authored numeric data, not identified source material packs.


### Authored tapered bores and primitive gas fields

`plate-valve-dg-bore.performance` uses the same runtime with an explicit
cylinder/cone bore. Replace the single `dg_tube` row with:

```text
dg_bore SEGMENTS OUTLET ELEMENTS ORDER POLES F_MIN_HZ F_MAX_HZ MAX_REL_ERROR EXP_TOL MAX_DOFS WALL_QUADRATURE_POINTS
dg_segment cylinder RADIUS_M LENGTH_M
dg_segment cone INLET_RADIUS_M OUTLET_RADIUS_M LENGTH_M
```

Supply exactly SEGMENTS geometry rows. The original duct partitions at geometry
boundaries. WALL_QUADRATURE_POINTS must be at least ORDER+1 and at most 64.
The complete declared state budget is checked before dense construction:
2 * ELEMENTS * (ORDER+1) acoustic coordinates and
2 * ELEMENTS * WALL_QUADRATURE_POINTS * POLES wall coordinates. The existing
owner enforces its own mesh and global state limits.

The spatial viscothermal owner evaluates local radius and material gas at its
positive quadrature points, projects through the original consistent mass
matrices and retains all causal diffusion states in the same exponential.
On a taper this is the **local-cylinder Webster approximation**, not a full
three-dimensional conical thermoviscous solution. Changing segment data changes
the same source-bound physical operator; it does not select an instrument preset.

An optional header suffix supplies a frozen primitive gas field:

```text
dg_bore 2 open 2 2 12 50 2000 0.01 1e-11 224 4 gas_field 3 mole 2 A B
dg_segment cylinder 0.001 0.025
dg_segment cone 0.001 0.002 0.025
gas_station 0     300 1   0
gas_station 0.025 320 0.5 0.5
gas_station 0.05  310 0   1
observation inlet
```

A and B are explicit IDs in the caller's thermochemical source registry, each
bound to its normalized species, caloric and transport packs. They are not
built-in materials. The same species-column reader also serves prescribed and
thermally generated reed fields. Exact absent-species fractions remain zero.
The original material owner evaluates temperature/composition between stations;
derived density, sound speed, viscosity and conductivity are never averaged.
The ambient record supplies static pressure. Stations must cover the exact bore
length, and the original mesh must resolve every geometry and field boundary.

The **actual inlet field state** supplies aperture density and slit viscosity.
The rest of the field drives propagation and viscothermal memory. Thus changing
primitive gas data affects both the excitation and the resonator without editing
their equations. The same shared DG source/load also supports other nonlinear
sources; no second gas law, bore fit or pressure observer is introduced.

The output remains internal step-mean inlet pressure. The sidecar retains the
actual inlet/outlet gas and area, spatial wall rule, canonical field stations and
source-pack receipts. Gas-field data is a frozen background: thermal evolution,
mean-flow acoustic corrections, moving solid walls and an exterior microphone
are separate physical couplings.
