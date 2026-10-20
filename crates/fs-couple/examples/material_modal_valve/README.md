# Material modal valves through the wind renderer

This input path makes the existing multimode material-wall solver available to
`music_render wind`. One or two independently declared plates load both a tract
and a bore. Original regional material cards determine density, stiffness and
the complete selected viscous tensor. Every requested signed mode, slit trace,
pressure traction and swept-volume projection stays in the same coupled system.
A shared exterior declaration drives both selected radiation reaction and the
microphone; pressure is never replaced by an internal bore signal.

The supplied geometry and material values are **manufactured integration data**.
They are not wet-cane measurements or a qualified clarinet/oboe specimen.

## Run

With the repository's admitted Rust toolchain and dependencies:

```bash
modal_demo_dir="$(mktemp -d)"
cp crates/fs-couple/examples/material_modal_valve/*.tsv "$modal_demo_dir/"
cp crates/fs-couple/examples/material_modal_valve/*.performance "$modal_demo_dir/"
cp crates/fs-couple/examples/material_modal_valve/*.radiation "$modal_demo_dir/"
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

`paired-keyed.performance` adds one material cantilever key on an explicitly
declared downstream chimney. Supply `key.matpack` through the same normalizer,
or copy the manufactured `blade.matpack` to that name for the integration
fixture, then run the unchanged renderer command with a new WAV path. Its
initial 0.0005 N closing force is explicit data. Replace only `key.matpack` to
hear the key material change through the retained reed and bore. Typed force
ramps use the optional control declarations described below.

`paired-keyed-sampled.performance` adds material keys on both sides and explicit
mechanical-sample force and exterior-pressure assignments. It reuses
`blade.matpack` for both keys and the original six reed coordinates. The two
chimneys retain their own gas temperatures; the downstream chimney has
distributed viscothermal loss. This fixture uses the alternative
`boundary_groups` declaration described below and needs no additional pack.

`paired-forced.performance` adds independent prescribed loads on both original
reed plates while retaining those material keys, their reservoir programs, the
no-slip supply and all six reed modes. After normalizing `blade.matpack` above,
render this complete example with:

```bash
cargo run --release -p fs-couple --bin music_render -- wind \
  "$modal_demo_dir/paired-forced.performance" "$modal_demo_dir/paired-forced.wav" --block 37
```

`paired-mutual.performance` selects an externally supplied two-port operator
in `mutual.radiation`. Run it through the same command with a new output path.
It adds a second explicitly declared downstream branch and uses the shared
MIMO terminal host, so off-diagonal reaction returns to the original material
valve before the same exterior microphone emits pressure. The supplied matrix
is synthetic integration data, not a measured or geometry-validated radiation
bake. The original `paired.performance` keeps compact baffled reaction.

`paired-walls.performance` adds independently declared material walls in both
the upstream tract and downstream bore. Normalize its additional source and
render it through the same ordinary entry point:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$modal_demo_dir/wall-manifest.tsv" --out "$modal_demo_dir/wall.matpack" \
  --material-state synthetic solid manufactured-membrane 0
cargo run --release -p fs-couple --bin music_render -- wind \
  "$modal_demo_dir/paired-walls.performance" "$modal_demo_dir/paired-walls.wav" --block 37
```

This is a manufactured compliant-cylinder demonstration, not acquired player
anatomy. Changing its wall pack changes actual wall motion and returning duct
pressure in the existing valve solve, and therefore the emitted sound.

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
radiation mimo GROUP_COUNT
radiation_group upstream|downstream ARTIFACT MAX_STATES MAX_MEMORY_BYTES PASSIVITY_POINTS
walls COUNT MIN_HZ MAX_HZ
wall upstream|downstream NODE isotropic PACK SOLID_T_K SOURCE_OMEGA_RAD_S THICKNESS_M stress-free|strain-fixed elastic|kelvin-voigt
valve TOTAL_REST_OPENING_M
limits MAX_DISPLACEMENT_M MAX_VELOCITY_M_S MAX_ENERGY_J EVALUATIONS VELOCITY_TOLERANCE_M_S RELATIVE_TOLERANCE
coupling MAX_ENERGY_J ABS_ENERGY_J REL_ENERGY ABS_FLOW_M3_S REL_FLOW
contact EFFECTIVE_STIFFNESS_N_PER_M_ALPHA EXPONENT INTERNAL_LOSS_S_PER_M
loss tensor V11 V22 V66 V12 V16 V26 REFERENCE_RESISTANCE_N_S_M PROJECTION_TERMS
blades COUNT
BLADE_RECORDS
upstream_boundaries
BOUNDARY_RECORDS
downstream_boundaries
BOUNDARY_RECORDS
control_source ASSEMBLY PART
source_forces COUNT
source_force MECHANICAL_SAMPLE BLADE_INDEX FORCE_N
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

The radiation rows are optional. Omission or `radiation compact-baffled` retains
the original compact reaction. MIMO mode requires one group per observed side,
at most two groups, each in that side's original opening order. Unobserved
sides and unselected outlets retain their authored loads and histories;
independent upstream/downstream groups do not couple through the exterior.
The operator's independently declared gas must exactly match the receiver,
and its geometry and band must match the actual ports and mechanics clock.
See [the shared artifact grammar](../MIMO_RADIATION_ARTIFACTS.md), including
thermochemical source-pack and complete gas-state binding. The importer never
relabels old coefficients with a changed current atmosphere.

Alternatively, `radiation computed-baffled GROUP_COUNT` accepts the same side
groups with fourteen numerical controls in place of the artifact path/limits.
It derives the matrix from that side's actual opening radii/centers and common
exterior, so data-only geometry or gas changes rebuild the physical operator.
The [shared source guide](../MIMO_RADIATION_ARTIFACTS.md) documents the fields,
finite-band approximation, three-order admission and separate memory budgets.
Both routes attach through the same original terminal host and preserve the
material reed, source-force projections, supply and complete duct histories.

The `walls` block is optional and follows any radiation records. Omission or
`walls none` preserves the original duct loads and arithmetic. Each wall binds
an unused degree-two junction between equal-radius cylindrical sections. Its
physical strip covers half of each adjacent section; unselected strips remain
rigid. Walls never replace terminal radiation, cavities, preexisting wall loads
or branching junctions. There are at most 64 selected wall nodes on either side.

Every wall resolves density, Young modulus and Poisson ratio from its declared
pack at the explicit solid temperature and angular frequency. `kelvin-voigt`
also requires the dimensioned `kelvin_voigt_young_viscosity` property; `elastic`
explicitly chooses zero internal loss. Source claims must be constant at that
solid temperature and valid across the whole wall band. That band must cover
the declared propagation and microphone bands. Unsupported material data,
thickness/radius above 0.05, missing coefficients or invalid graph sites refuse.

The wall shares the plate source cache and total material-byte budget, while
keeping its independent state query. The original regional gas stays in each
duct section; neither solid temperature nor its material name selects a new
gas. The same material-wall owner participates in pressure/flow work and storage
with compact or MIMO terminal radiation, and with or without a supply channel.
It is a fixed-state local isotropic membrane approximation: no full shell,
orthotropic tissue, axial wall coupling or exterior shell radiation is implied.

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

### Physical force on each original source plate

The optional `source_forces` block follows every blade and any key/control
declarations, immediately before `compile_limits`. Omission or `source_forces
none` retains the original pressure-only source drive. Each `source_force` row
names a zero-based mechanical sample, an existing zero-based blade index and a
finite signed force in newtons. Positive force acts toward opening; a negative
force closes. Each blade begins at zero force and holds its last assignment.
Zero releases only the actuator while the original mechanics keep evolving.

```text
source_forces 4
source_force 97 0 -0.0001
source_force 225 1 -0.00015
source_force 769 0 0
source_force 2049 1 0
```

Each blade's `node` records already declare the physical `FORCE_WEIGHT`
footprint. The load multiplies those original nodal weights; the original plate
reduction projects it into **every retained mode**. The same signed projection
maps modal velocity back to actuator velocity, so mechanical work remains
force times displacement at that physical port. A force is never substituted
for blowing pressure or applied to an independently animated reed. Material,
thickness, support and force-site changes therefore alter the physical response
through the existing source basis.

Assignments may be unsorted. At one sample, input order is retained and the
last assignment to a blade wins; intermediate assignments do not advance
physics. Both blades share control storage and auxiliary compilation work with
all key programs and the existing pressure schedule. Events execute before
their mechanical step, including within one output audio frame, and commit
only with accepted physics. The input and sidecar retain the node weights,
modal force projections and mechanical-time force program. These are prescribed
actuator loads; they do not supply compliant lip or tongue contact, muscle
dynamics, or a new source of pressure work.

### Material moving keys on either acoustic side

The optional `upstream_boundaries` and `downstream_boundaries` blocks follow all
blade records, in that order, before `control_source` and `compile_limits`.
Each selected side uses the **same material-boundary grammar and physical
builder as the flute input**. Add a `duct_node driven` endpoint connected by its
own chimney section; every driven endpoint must have exactly one key in
ascending node order. An omitted block is valid only if that side has no driven
endpoints. A driven key is a coupled source boundary and is not a passive
radiating `opening`.

```text
downstream_boundaries
boundaries 1
boundary 3 instrument key-1 0 0
boundary_plate key.matpack 293.15 1000 0.025 0.01 0.0003 6 2 0.0002 0.0002 0
boundary_mode 0 600 0 0.25 0.1
boundary_contact 10000000000 2 0.5
boundary_slit 0.0003
boundary_loss elastic
control_source instrument reed
```

The `boundary` fields are graph node, assembly name, part name, initial closing
force in newtons and exterior reservoir pressure in pascals. The plate fields
are original material pack, solid temperature, source angular frequency,
length/width/thickness in metres, mesh intervals, rest gap, initial gap and
opening velocity. The material card supplies density, Young's modulus and
Poisson ratio. `boundary_loss kelvin-voigt` also selects the card's
`kelvin_voigt_bending_viscosity` [Pa s] at the same solid temperature and query
frequency. The existing plate owner projects that loss on the retained mode;
all selected properties must cover the mode search band. Changing only this
material property changes key motion and coupled pressure without retuning
mass, stiffness or authored damping. Zero viscosity remains an explicit
source-bound loss selection. `boundary_loss maxwell COUNT MIN_HZ MAX_HZ
MAX_DT_OVER_TAU MAX_ANGULAR_STEP` plus `boundary_term MODULUS_PROPERTY
TIME_PROPERTY` records select the same material Maxwell histories as flute
keys. These solid-loss selections are mutually exclusive; `boundary_loss elastic`
retains lossless bending. `boundary_slit LENGTH_M` selects the actual chimney gas's viscosity;
`boundary_slit none` or omission keeps inviscid gap flow.

All blades and keys use one cached material loader, with one 64 MiB aggregate
allowance and an 8 MiB limit per normalized pack. Reusing a path reuses the
same immutable artifact. Changing a key's material, thickness or gas changes
the physical motion and pressure returned to the original reed/bore solve.
The source retains every selected reed mode, both acoustic networks, fractional
and distributed-loss histories, and compact or supplied mutual outlet reaction.

`control_source` selects an all-typed canonical v2 gesture suffix. Use
`two_sided::observer::schedule::boundaries::source_pressure_path(assembly, part)`
for the declared source's `PressurePa` track and the shared
`material_boundary::boundary_closing_force_path(assembly, part)` for each key's
`ForceN` track. Positive force closes the material slit. A force ramp is sampled
by the existing gesture owner, then applied before the corresponding mechanical
sample, including events within a single audio frame. Zero force releases the
actuator while mechanics, contact, material memory and trapped pressure continue.
Untracked keys retain their explicit initial drive. Without `control_source`,
the original one-track pressure suffix remains valid and all keys retain those
initial drives. Unknown consumers and mixed legacy/typed tracks refuse.

For explicit force and exterior-reservoir assignments at mechanical sample
indices, replace the optional side blocks and `control_source` with one
`boundary_groups` block, still **after every blade and before `compile_limits`**:

```text
boundary_groups 1
boundary_side downstream
boundaries 1
boundary 3 instrument bore-key 0 0
boundary_plate blade.matpack 293.15 100 0.025 0.01 0.0003 6 2 0.0002 0.0002 0
boundary_mode 0 600 0 0.25 0.1
boundary_contact 10000000000 2 0.5
boundary_slit 0.0003
boundary_loss elastic
boundary_controls 2
boundary_control 321 3 0.003 0.3
boundary_control 1505 3 0 0
compile_limits 100000 8192
```

There may be one or two groups, each naming a different side. Each group uses
the unchanged shared `boundaries` grammar and must declare at least one key.
An omitted side must have no driven endpoints. The `boundary_control` fields
are zero-based mechanical sample, graph node, nonnegative closing force in
newtons and finite exterior reservoir pressure in pascals. Each row replaces
both drive fields for that key before the named mechanical step; other keys
retain their previous drives. Initial drives come from the `boundary` records,
and every drive remains held until its next assignment. The existing compiler
sorts assignments by sample and node. Simultaneous changes on different keys
are allowed; repeated assignments to the same key at one sample and samples
outside the source horizon refuse. Both sides share control storage with the
source pressure and source-force assignments. Boundary and source-force
compilation share the existing auxiliary work allowance; the original pressure
compiler keeps its independent allowance. Sampled groups retain the original
one-track pressure gesture suffix. Combining them with `control_source`,
`upstream_boundaries` or `downstream_boundaries` refuses, so one scheduler owns
each key's complete force and reservoir drive.

Both declaration routes attach keys to the final lowered acoustic graph before
the original microphone is admitted. Each key uses that chimney's actual gas,
represented impedance, mechanics clock and horizon, including distributed-loss
lowering. Original source cards share the blade cache, and physical key paths
must be unique across the two sides. Compact, imported MIMO and computed MIMO
radiation keep their original selected exterior bindings and histories. The
sidecar preserves the shared material-key source records and adds complete
sampled drives with their mechanical clock; it does not reinterpret reservoir
pressure as a static atmosphere or change the gas mixture during a run.

The shared key model is a material-derived one-mode cantilever slit with explicit
contact, not a resolved circular pad, linkage or tongue. Only declared passive
outlets feed the existing exterior microphone; direct key/body radiation is not
added. Source pressure and independent source-mesh forces remain distinct from
key actuation. New integration regressions require native execution before any
numerical or measured-instrument qualification claim.

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
and either compact independent or source-bound mutual infinite-baffle radiation.
It does not supply distributed
reed-surface collision, shell/solid wet-cane geometry, resolved moving-channel
flow, nonlinear jet separation, unbaffled bells, delayed broadband radiation, or live
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
Native MIMO regressions additionally retain source byte hashes and declared
bake identities separately, preserve all material/gas/loss/supply owners,
exercise off-diagonal feedback and PCM partition parity, reject stale
atmospheres and source mismatches, and compare the actual CLI against the
original library performer before any WAV/sidecar publication.
The sampled-key regressions additionally compare both-side subframe drives
against the original direct physics and receiver, retain Maxwell/supply/loss
and mutual-radiation memories, and exercise material, local-gas, force and
reservoir substitutions in the same pressure and PCM path. The CLI regression
retains the upstream typed-key example and separately renders the sampled
example; malformed key programs refuse before either output is created.
