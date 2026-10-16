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
initial OPENING_COORDINATE_M VELOCITY_M_PER_S
contact K_PA_PER_M_ALPHA ALPHA CHI_S_PER_M MAX_PENETRATION_M
source SINGLE_TOKEN_SOURCE_LABEL
sections COUNT
section isotropic THICKNESS_M DENSITY_KG_PER_M3 YOUNG_PA POISSON_RATIO
# Alternatively: section orthotropic THICKNESS_M DENSITY_KG_PER_M3 E1_PA E2_PA NU12 G12_PA ANGLE_RAD
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
it does not change the separately supplied solid material constants. The tube
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

Library consumers can construct `AperturePerformance` directly around either an
`ApertureTube` or `ApertureNetwork`, including already-attached passive loads and
material-bound specimens. The observation explicitly selects inlet, tube terminal
or a network node. `PlateValvePerformance::into_renderer()` returns this complete
finite object, suitable for the existing `DecimatedRenderer`, `PressureEnsemble`
and incremental PCM interfaces. Keep that whole chain during cancellation/resume.
A refused individual mechanical step retains its previous state; a failed output
callback poisons the pressure adapter because an earlier prefix may have advanced.
The accepted prefix and pending controls remain inspectable, never reported as a
complete successful callback.

This remains one linear structural mode with uniform face pressure, a fixed lay
and quasisteady Bernoulli flow. Supplied constitutive data is not wet-cane
identification, mesh conformity is the caller's responsibility, and there is no
new continuum-convergence, calibrated-microphone or real-time-throughput claim.
