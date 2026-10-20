# Material modal valves through the wind renderer

This input path makes the existing multimode material-wall solver available to
`music_render wind`. One or two independently declared plates load both a tract
and a bore. Original regional material cards determine density, stiffness and
the complete selected viscous tensor. Every requested signed mode, slit trace,
pressure traction and swept-volume projection stays in the same coupled system.
A shared exterior declaration drives both compact radiation reaction and the
microphone; pressure is never replaced by an internal bore signal.

The supplied geometry and material values are **manufactured integration data**.
They are not wet-cane measurements or a qualified clarinet/oboe specimen.

## Run

With the repository's admitted Rust toolchain and dependencies:

```bash
modal_demo_dir="$(mktemp -d)"
cp crates/fs-couple/examples/material_modal_valve/*.tsv "$modal_demo_dir/"
cp crates/fs-couple/examples/material_modal_valve/*.performance "$modal_demo_dir/"
cargo run -p xtask -- matdb-pack \
  --manifest "$modal_demo_dir/blade-manifest.tsv" --out "$modal_demo_dir/blade.matpack" \
  --material-state synthetic solid manufactured 0
cargo run --release -p fs-couple --bin music_render -- wind \
  "$modal_demo_dir/paired.performance" "$modal_demo_dir/paired.wav" --block 37
```

The example declares six physical mode coordinates, two thickness regions on
each original mesh, different regional gas temperatures, distributed bore loss,
fractional travel, gap viscosity and a physical exterior receiver. Mechanical
steps run at 768 kHz and pressure frames at 48 kHz. The renderer consumes every
mechanical step; do not add `--decimate`. WAV and sidecar paths must be new.
The example explicitly selects `propagator exact-zoh`. The clock and retained
mode bounds stay authored data; exact free structural propagation is not a
reason to lower the clock without measuring nonlinear time error.

To change one blade or region, normalize another material source and change its
`section` pack path. Original geometry and numerical choices remain explicit.
No pitch normalization or output equalizer compensates for that substitution.
Solid state queries are independent of the declared gas state. Material and
gas changes construct a new complete source; they are not live evolving fields.

## Source grammar

Use LF lines, fixed record order, and no empty header lines. The complete example
is the grammar fixture. Gas records and duct records use the existing
[plate-valve grammar](../PLATE_VALVE_PERFORMANCES.md); numerical loads,
junctions, cavity/wall nodes, regional gas mixtures and optional viscothermal
section suffixes retain their original owners.

```text
frankensim-material-modal-valve-v1
audio MECHANICS_HZ OUTPUT_HZ OUTPUT_FRAMES FULL_SCALE_PA
propagator midpoint|exact-zoh
ambient GAS_RECORD
exterior GAS_RECORD
upstream
TUBE_OR_NETWORK_RECORDS
downstream
TUBE_OR_NETWORK_RECORDS
propagation fractional MAX_HZ SUMMED_PHASE_ERROR_RAD
supply none
slit_length LENGTH_M
receiver X Y Z RADIAL_RINGS ANGULAR_POINTS MAX_HZ
openings COUNT
opening upstream|downstream NODE X Y
valve TOTAL_REST_OPENING_M
limits MAX_DISPLACEMENT_M MAX_VELOCITY_M_S MAX_ENERGY_J EVALUATIONS VELOCITY_TOLERANCE_M_S RELATIVE_TOLERANCE
coupling MAX_ENERGY_J ABS_ENERGY_J REL_ENERGY ABS_FLOW_M3_S REL_FLOW
contact EFFECTIVE_STIFFNESS_N_PER_M_ALPHA EXPONENT INTERNAL_LOSS_S_PER_M
loss tensor V11 V22 V66 V12 V16 V26 REFERENCE_RESISTANCE_N_S_M PROJECTION_TERMS
blades COUNT
BLADE_RECORDS
compile_limits WORK CONTROLS
schedule
CANONICAL_GESTURE_BYTES
```

Both inlets must contain the same complete gas state. Remote gas states and the
exterior must share static pressure. Gas records accept the original named
mixture and thermochemical-source selections. Radiating openings must name cold
memoryless endpoints; the common exterior installs each load once. Do not place
`baffled-low-ka` loads inside these input ducts. Other explicit loads retain
their declared mechanical role. Outlet radii come from the original geometry.

`propagation integer` selects the original integer delay image.
The separate `propagator` row is optional and appears immediately after `audio`.
Omitting it preserves the original `midpoint` method and its complete loaded
one-cent clock admission. `exact-zoh` uses the original mass, stiffness and full
material viscosity matrix in the shared exact structural propagator. Pressure,
contact and external forces are held over each mechanical sample; swept volume
and force work use actual displacement increments, and viscous dissipation is
integrated along the same trajectory. The wind sidecar reports the selection
and its physical scope. Both methods retain the original two ducts, supply
passage, gas states, contact, materials and matched exterior microphone.

With `exact-zoh`, each blade's `MAX_ANGULAR_STEP` is an explicit positive bound
strictly below pi. Every retained mode must satisfy that bound. Midpoint keeps
its original cap and stricter loaded phase gate. Raising the exact bound admits
a coarser numerical rung; it does not establish convergence of the nonlinear
valve, contact, channel flow or acoustic network.

`slit_length none` selects the original inviscid gap law.
`supply slit LENGTH_M WIDTH_M HEIGHT_M TRANSVERSE_CELLS` selects the shared
unsteady no-slip fixed passage; it is additional geometry, separate from the
moving gap. The selected structural propagator and the original channel,
network, nonlinear work and flow admissions all remain binding. This ingress
neither changes the passage nor substitutes a static resistance for its state.

Each blade repeats:

```text
blade MODES MIN_HZ MAX_HZ MAX_ANGULAR_STEP MAX_TRACE_SLOPE MAX_SLIT_DEVIATION_M PROJECTION_TERMS
support clamped|simply-supported PRETENSION
sections COUNT
section PACK isotropic|orthotropic12|orthotropic21 SOLID_T_K SOURCE_OMEGA_RAD_S THICKNESS_M MATERIAL_ANGLE_RAD
nodes COUNT
node X Y SUPPORT_0_OR_1 UNIT_FORCE_WEIGHT
triangles COUNT
triangle A B C SECTION
slit_edges COUNT
edge A B
initial MODE_COUNT
state DISPLACEMENT_M VELOCITY_M_S
```

Sections cover every triangle exactly once and every declared section must be
used. Shared nodes bond the regions. Fixed thickness means mass changes with
density. Ordered orthotropic planes and material angle are resolved by the
existing regional plate adapter. Boundary edges define the physical slit width;
both blades must have the same width. The force weights define each blade's
existing force column; the file schedules blowing pressure, and the retained
library performer also exposes the original held blade-force ports.

`loss elastic` explicitly selects zero material dissipation. The tensor form
requires all six named plane-stress viscosity properties in every region,
including explicitly sourced zeros. Diagonal damping ratios are not synthesized.
The shared source adapter admits passivity, fixed-state applicability and actual
retained frequency support, then projects the full matrix in the same basis as
pressure and displacement. No damping matrix is copied from a replaced material.

The shared limit is 32 total modes across one or two blades, with at most 512
nodes and 2048 triangles per blade. Each normalized pack is capped at 8 MiB and
the complete distinct-source set at 64 MiB. Mechanical/output ratio is an integer
1–16. Whole-program work and control budgets are declared once.

## Physical scope and verification

This is a source/renderer integration over flat linear DKT plates, a finite
declared mode basis, an effective uniform slit/contact, plane-wave duct graphs,
and compact independent infinite-baffle radiation. It does not supply distributed
reed-surface collision, shell/solid wet-cane geometry, resolved moving-channel
flow, nonlinear jet separation, unbaffled bells, mutual radiation, or live
material/gas evolution. It does not qualify instrument recordings.

Focused regressions are in
`bernoulli_aperture::performance::file::material_modal::tests`: original source
and complete coordinate retention, physical pressure changes from material and
atmospheric-pressure substitutions, callback/PCM identity, and missing-source,
clock, count, gas-pressure and double-radiation refusals. Additional selector
regressions preserve omission/explicit-midpoint parity, exercise full-tensor and
elastic exact material sources at a coarser clock, and retain material changes
in the actual microphone signal. Native compilation and
execution remain required; committed tests are not evidence of a passed run.
