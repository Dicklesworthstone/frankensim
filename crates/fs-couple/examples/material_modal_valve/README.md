# Material modal valves through the wind renderer

This input path makes the existing multimode material-wall solver available to
`music_render wind`. One or two independently declared plates load both a tract
and a bore. Original regional material cards determine density, stiffness and
the selected viscous tensor, regional hereditary spectrum or complete 3-D
tensor-Maxwell law. Every requested signed mode, slit trace,
pressure traction and swept-volume projection stays in the same coupled system.
A shared exterior declaration drives both selected radiation reaction and the
microphone; pressure is never replaced by an internal bore signal.

The supplied geometry and material values are **manufactured integration data**.
They are not wet-cane measurements or a qualified clarinet/oboe specimen.

The DG variant also accepts `dg_molecular LAW` immediately after each complete
`dg_tube` or `dg_bore` block and before its optional `dg_wall`, independently for the bore and
upstream tract. It uses the same evaluated species/reservoir law as the flute,
finite side chimneys and exterior observers. `observer_gas molecular MAX_PATHS
MAX_FILTER_TAPS MAX_KERNEL_TERMS LAW` binds the exterior's own composition.
See [the shared gas declaration and its physical limits](../FLUE_NETWORK_PERFORMANCES.md#molecular-gas-response-inside-dg-bores-and-along-exterior-paths).
The native allowance includes every molecular state at the actual nodal or
spatial wall quadrature coordinates. Tapered geometry and prescribed axial
gas fields use the original local-cylinder thermal diffusion owner. A field
evaluates its molecular law at the actual thermochemical material points,
retains the integrated frozen pressure mass and preserves pressure/flow work.
The main bore and upstream tract may bind independently declared field laws.
Absent species data, unsupported local temperatures and ineligible ISO mixtures
refuse; the source never substitutes a generic air-loss coefficient. Fields
remain frozen during a performance; live molecular-history transfer is absent.

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

`paired-spatial.performance` uses the same material packs and renderer with
actual nodal clearances and distributed blade-to-blade contact. Every retained
mode changes local opening, partial edge closure and the gas-dependent channel
flow. Its interface is deliberately manufactured, with a 20–60 micrometre
rest-gap profile; it is not a measured oboe reed. Render it by selecting that
filename in the command above.

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

`paired-pressure.performance` supplies independent distributed pressure loads
on both source blades through the same renderer. The first trace changes sign
from root to tip and has zero net force, yet it bends the blade because the
retained displacement modes are nonuniform. The second trace rises from zero
at the root to its reference amplitude at the tip. These are manufactured
external pressure fields, not pressures inferred from the reed gap or acquired
player measurements. Selecting a different material pack rebuilds their modal
force/velocity overlaps on the new specimen and changes actual exterior sound.
Both fields retain separate signed work accounts; uniform blowing pressure and
the existing physical blade swept flow remain coupled exactly once.

`paired-wetted.performance` exposes only the tip half of the first material
blade to the coupled upstream/downstream pressure port. Original triangle
addresses select the physical upstream face. Its entire material mass,
stiffness, viscosity, force footprint and contact mechanics remain present;
the selected geometry changes pressure traction and reciprocal swept flow
together. This manufactured fixture uses the same `blade.matpack` and ordinary
renderer command. Change a facet's `positive` side to `negative` to move the
upstream fluid to the opposite physical side, or replace the material pack to
rebuild all retained modal properties and geometric pressure overlaps.

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
OPTIONAL_WETTED_SURFACE_RECORDS
OPTIONAL_PRESSURE_SHAPE_RECORDS
OPTIONAL_SPATIAL_CLOSURE_RECORDS
upstream_boundaries
BOUNDARY_RECORDS
downstream_boundaries
BOUNDARY_RECORDS
control_source ASSEMBLY PART
source_forces COUNT
source_force MECHANICAL_SAMPLE BLADE_INDEX FORCE_N
source_pressures COUNT
source_pressure MECHANICAL_SAMPLE BLADE_INDEX PRESSURE_PA
compile_limits WORK CONTROLS
schedule
CANONICAL_GESTURE_BYTES
```

Optional physical wetted surfaces follow all blade records, before prescribed
pressure traces, spatial closure or boundary declarations:

```text
wetted_surfaces COUNT
wetted_surface BLADE_INDEX PROJECTION_TERMS SOURCE_LABEL
wetted_facets COUNT
wetted_facet ORIGINAL_TRIANGLE_INDEX positive|negative
```

Repeat the facet row for each exposed original triangle and the surface block
for each selected blade. `positive` and `negative` identify which physical side
of the triangle's original winding faces **upstream**. The opposite side faces
downstream. No signed effective pressure area is supplied. The shared oriented
surface owner integrates each retained displacement mode over the actual
selected faces. Those same geometric overlaps drive both pressure traction
and swept-volume feedback in the existing nonlinear acoustic solve. Swapping
the upstream face changes pressure forcing relative to opening motion; selecting
only the exposed part of a blade changes its loading without removing any mass,
stiffness, material loss, contact or physical force footprint.

At most one surface can be attached to each blade. Original specimen identity,
triangle addresses, source label and projection budget are retained. Repeated
or unknown triangles, foreign specimen bindings and late attachment refuse.
Omission keeps the previous complete positive-side pressure surface and its
arithmetic. A `wetted_surfaces 0` block also makes no change. A selected blade
with `wetted_facets 0` is explicitly sheltered from this coupled pressure port;
it can still move under external force/contact and change the flow opening.

The same attachment works with full tensor viscosity, exact structural
propagation, proportional or regional Maxwell memory, spatial contact and
independent prescribed-pressure ports. Only the cached acoustic response is
rebound; the material operators and initial internal strains are retained.
The approximation remains a fixed P1 normal trace of flat DKT modes. It does
not model moving wetting boundaries, local channel-pressure recovery, a curved
shell blade or a measured oscillation threshold.

Optional additional pressure traces follow the blade and wetted-surface records, before spatial
closure or boundary declarations:

```text
pressure_shapes COUNT
pressure_shape BLADE_INDEX PROJECTION_TERMS SOURCE_LABEL
pressure_nodes ORIGINAL_NODE_COUNT
pressure_node DIMENSIONLESS_SHAPE_VALUE
```

Repeat the node row for every original node, and the shape block for each
selected blade. At most one trace per blade is admitted. All values are signed
and finite; they are **not normalized**. An amplitude `P` in pascals applies
the local pressure `shape(x) * P`. This permits balanced and zero-net-force
loads without dividing by their resultant force. The shared consistent P1
triangle quadrature projects each trace onto the blade's already retained
material modes. Complete source identity, all original node addresses and the
explicit projection-work budget are retained. Missing nodes, duplicate blade
bindings, ineffective or nonfinite projections, and late attachment refuse.

The optional `source_pressures` block follows `source_forces` (when present)
and precedes `compile_limits`. Count zero explicitly supplies no assignments;
omission preserves the original pressure-only or force-driven render. Pressure
events use mechanical sample addresses, not output frames. Unsorted records
are stably ordered; the last assignment at an equal sample wins. Zero releases
the source while keeping vibration, contact and material histories. The same
combined control/work budgets include these records, ordinary blowing pressure,
mechanical forces and moving keys.

Spatial pressure work is `P * integral(shape * blade_velocity dA) * dt`, with
the displacement increment from the original exact or midpoint solver. It is
reported separately per blade and included once in total mechanical work.
That generalized volume flow is not added to acoustic through-flow. These
explicit external pressure inputs do not solve a moving-gap fluid field and
must not repeat the uniform bore-pressure difference already owned by the
valve. Original elastic, tensor-loss, Maxwell and spatial-contact paths all
use the same traction input and pressure observer.

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

Every isotropic wall resolves density, Young modulus and Poisson ratio from its declared
pack at the explicit solid temperature and angular frequency. `kelvin-voigt`
also requires the dimensioned `kelvin_voigt_young_viscosity` property; `elastic`
explicitly chooses zero internal loss. Source claims must be constant at that
solid temperature and valid across the whole wall band. That band must cover
the declared propagation and microphone bands. Unsupported material data,
thickness/radius above 0.05, missing coefficients or invalid graph sites refuse.

The same original graph sites also accept a principal-aligned orthotropic wall:

```text
wall upstream|downstream NODE orthotropic PACK SOLID_T_K SOURCE_OMEGA_RAD_S THICKNESS_M AXIAL HOOP_AXIS AXIAL_AXIS LOSS
```

`AXIAL` retains `stress-free` or `strain-fixed`. The distinct principal axes are
one-based integers in 1..=3, ordered hoop then axial, and follow the local
cylindrical directions. `LOSS` is `elastic`,
`proportional-kelvin-voigt TIME_PROPERTY`, or
`normal-kelvin-voigt MAX_DT_OVER_TAU MAX_ANGULAR_STEP HH HA AA`.
The original nonnegative proportional time claim
declares `C_viscous = tau C_elastic` on this membrane plane. Density, the two
directional Young moduli and their major Poisson ratio are resolved directly.
The same constitutive owner serves graph jets and DG ducts, including reciprocal
Poisson reduction under fixed axial strain. No equivalent isotropic card or
material-name rule is introduced.

The independent normal-viscosity selection names the three Pa·s entries of a
symmetric positive-semidefinite viscosity matrix in the local hoop/axial plane;
`HA` may be negative. Fixed axial strain keeps the original hoop entry. Free
axial stress eliminates the complete `C + sD` law and retains the resulting
Maxwell arm and residual Kelvin resistance in the original duct load. Changing
only a source cross-viscosity value can therefore change reed feedback and
external pressure even when the base wall mass, stiffness and resistance stay
the same. Its explicit resolution limits satisfy `0 < MAX_DT_OVER_TAU <= 2` and
`0 < MAX_ANGULAR_STEP <= 1` and are checked against the actual graph clock when
memory is present. Native provenance preserves the source tensor names and
values, material axes, reduced memory and residual resistance. Full 3-D
viscoelastic condensation and unsupported directional spectra still refuse.

The wall shares the plate source cache and total material-byte budget, while
keeping its independent state query. The original regional gas stays in each
duct section; neither solid temperature nor its material name selects a new
gas. The same material-wall owner participates in pressure/flow work and storage
with compact or MIMO terminal radiation, and with or without a supply channel.
It is a fixed-state local membrane approximation: no full shell, off-axis grain,
axial wall coupling or exterior shell radiation is implied.

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

### Spatial reed/lay or blade/blade closure

After all blade records, before any boundary or source-force records, select:

```text
closure spatial PRESSURE_STIFFNESS_PA_PER_M_ALPHA EXPONENT INTERNAL_LOSS_S_PER_M MAX_PENETRATION_M SOURCE_TOKEN
clearances NODE_COUNT
gap TOTAL_REST_GAP_M
...
lay_triangles TRIANGLE_COUNT
lay_triangle ORIGINAL_TRIANGLE_INDEX
...
```

There must be one rest gap per original node, in original node order. Triangle
IDs select the physical contact patch; its triangle areas are the quadrature
weights. `SOURCE_TOKEN` identifies the explicitly supplied contact-pair law.
The two-blade form requires identical planar node/triangle/slit indexing, while
preserving each blade's independent material, support, thickness and mode basis.
One blade uses a stationary lay; two blades share one contact interface.

This selection replaces the earlier effective `contact` law in both dynamics
and stored energy. It does not add a second contact force. Each local gap is the
rest clearance plus every retained signed displacement contribution. Contact
uses the existing pointwise nonadhesive discrete force, then projects reaction
into every original mode, including modes with zero mean opening. Local
Hunt–Crossley loss is applied before those reactions are summed. The same-time
solve includes all modes, the original pressure/swept-volume faces and either
full material viscosity or the selected exact Maxwell history.

Slit flow integrates the positive part of the physical P1 edge gap, including
partially closed edges. With `slit_length`, the existing laminar/Bernoulli law
uses the actual local height and source gas viscosity; it does not replace the
profile by a mean-height channel. Final states retain the original trace-slope
and new local penetration bounds. The effective-slit deviation limit applies
only to the uniform interface after attachment; original constructors still
admit the supplied initial material state before this optional attachment.
Omission or `closure none` retains the existing uniform path.

This is a finite-mode, matching planar interface. Contact does not change the
basis, mask the pressure face, introduce curved/nonmatching surfaces or resolve
moving-wall fluid inertia. The supplied nonlinear evaluation budget includes
all modal Jacobian probes and backtracks; failure publishes no new modal or
material-memory state.

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
triangle A B C SECTION [thickness THICKNESS_M]
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


A triangle may override only its geometric thickness, independently of the
authored material section. For example, in the **manufactured** paired fixture:

```text
triangle 0 1 8 0 thickness 0.00015
triangle 0 8 7 0
```

The first triangle has thickness 0.15 mm; the second retains section 0's
declared 0.30 mm. Both use the same original material pack, state, grain axes
and constitutive selection. This lets a measured vamp or another plate
thickness map be refined without inventing a separate material for every
geometric value. All thicknesses must be finite and strictly positive.
The original 64-section material limit, 512-node/2048-triangle geometry limits,
pack-byte limits and declared projection budgets still apply. At most one
occupied chart region is needed per original triangle; equal thicknesses in
one material section share a region.

The adapter subdivides the original assignment through
`PlateRegionMaterial` and `compile_plate_material_chart`. The same regional
geometry therefore changes integrated mass, anisotropic bending stiffness,
full tensor viscosity or Maxwell memory, signed modes, pressure traction,
swept volume and source-mesh contact. No independent thickness correction,
modal frequency fit or output filter is applied. Material replacements retain
the authored thickness map. The native source records distinguish `section`
(the original material declaration) from `region` (the resolved thickness
subdivision), and every subdivision shares the original source artifact.

Omitting triangle thickness, or explicitly repeating each section's existing
value, retains the original region layout and arithmetic. Both the two-sided
and DG material-modal consumers use this same reader. This is a piecewise
constant finite-element thickness field; physical accuracy still requires
geometry/mode/time refinement and a complete measured specimen. The shipped
regressions use manufactured materials and geometry and do not assert a
qualified clarinet or double-reed reconstruction.

`loss elastic` explicitly selects zero material dissipation. The tensor form
requires all six named plane-stress viscosity properties in every region,
including explicitly sourced zeros. Diagonal damping ratios are not synthesized.
The shared source adapter admits passivity, fixed-state applicability and actual
retained frequency support, then projects the full matrix in the same basis as
pressure and displacement. No damping matrix is copied from a replaced material.

`loss proportional-maxwell relaxed COUNT MODULUS TIME ...` selects a sourced
generalized Maxwell solid with one to eight modulus/time pairs. Use
`unrelaxed` instead of `relaxed` to start each viscous strain at zero; `relaxed`
equilibrates the arms at the explicitly supplied initial displacement. The
ordinary Young modulus is the equilibrium modulus. Every additional modulus
is in pascals and every relaxation time is in seconds. For example, replace
the `loss tensor ...` row in `paired.performance` with:

```text
loss proportional-maxwell relaxed 1 maxwell_modulus maxwell_time
```

The supplied `blade.tsv` includes these manufactured coefficients. Normalize
it with the same command above, then render the changed file through the same
wind command. Changing only the modulus/time claims changes every retained
mode's restoring history and loss, the returning duct pressure, and the actual
exterior microphone. No modal damping ratios or output filters are fitted.

This selector requires `propagator exact-zoh` and zero plate pretension. Each
blade's regions must share the same normalized spectrum `(E_j/E_inf, tau_j)`;
its regions may still have different geometry, density, thickness, equilibrium
moduli and orientation. For orthotropic plates, the declared axis-1 Young
modulus sets the ratios that scale the whole original elastic tensor. This
proportional law makes branch stiffness diagonal in the original plate basis;
nonproportional spectra require full matrix projection and are refused. Every
mode retains its own internal strains, storage and independently integrated
dissipation in the same pressure/contact solve. Full tensor viscosity and
Maxwell relaxation are alternative explicit material laws.

The highest instantaneous elastic frequency must satisfy the original blade's
angular-step cap and every source coefficient's frequency validity. Positive
branches also require `dt/tau <= 2`. Neither exact free propagation nor the
source law establishes convergence of nonlinear contact or acoustic coupling.
These are fixed-temperature manufactured inputs, not wet-cane measurements.

`loss regional-maxwell relaxed|unrelaxed COUNT MODULUS TIME ...` allows each
original material region to have different modulus ratios and relaxation times.
For example:

```text
loss regional-maxwell relaxed 1 maxwell_modulus maxwell_time
```

Keep the original geometry and reference different normalized material packs
in the `section` records. The properties are read separately from each region's
original material state. The complete regional bending tensors are projected
using the same retained DKT displacement and rotation eigenvectors as the
pressure, gap and force ports. All cross-mode material forces remain in the
same exact structural solve, so changing one region's data reaches the duct
pressure and rendered PCM. No independent replacement reed is synthesized.

This selector requires `propagator exact-zoh`; it permits equilibrium
pretension, which is not included in the relaxing bending branches. The local
plane-stress tensor is scaled by each region's own modulus ratios. This is not
arbitrary three-dimensional tensor Prony condensation. The complete exact loss
block must satisfy `2*(3*modes+channels) <= 256`, in addition to the original
projection, angular-step, material-frequency and `dt/tau <= 2` limits. Singular
nonzero regional forms refuse instead of receiving artificial stiffness or a
rank cutoff. The proportional selector keeps its existing diagonal path.

The shared limit is 32 total modes across one or two blades, with at most 512
nodes and 2048 triangles per blade. Each normalized pack is capped at 8 MiB and
the complete distinct-source set at 64 MiB. Mechanical/output ratio is an integer
1–16. Whole-program work and control budgets are declared once.

### Full 3-D tensor material memory

Both native material-modal selectors accept the same source-bound law:

```text
loss tensor-maxwell relaxed|unrelaxed STRAIN_LIMIT_PROPERTY COUNT MAX_POINTS MAX_MEMORY_STATES MAX_PROJECTION_TERMS
tensor_branch TIME_PROPERTY C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66
```

One `tensor_branch` row follows for each branch. Every entry names a source
property; the 21 stiffness entries are the symmetric upper triangle in Mandel
order `(11,22,33,sqrt(2)12,sqrt(2)23,sqrt(2)31)` and have units Pa. Times are
positive seconds. The strain limit is a positive dimensionless Frobenius bound.
Zero couplings are explicit source claims and can reuse a declared zero property.
The canonical equilibrium elastic properties remain the same ones used to build
the plate: isotropic `young_modulus`/`poisson_ratio`, or the full three-dimensional
orthotropic engineering properties. No second equilibrium stiffness is added.
All claims retain their source temperature and intersected angular-frequency
domains. Missing values, active tensors or unsupported states refuse.

The original section orientation rotates the entire constitutive law. At every
original DKT area point the shared plane-stress owner eliminates the complete
normal traction constraint, including its coupled histories. An exact change
of history coordinates feeds the existing structural propagator and pressure/
contact solver; it does not fit scalar modal damping or condense each Maxwell
branch separately. Original six-component strains are reconstructed at the
physical outer faces before an accepted sample is published. Homogeneous linear
bending integrates thickness analytically. Nonzero 3-D Kelvin normal-rate terms
and evolving material temperature/moisture are outside this attachment.

`paired-tensor-dg.performance` supplies a complete manufactured example with
two two-triangle blades, two modes per blade, 12 material points and all 72
material-history coordinates of one full tensor branch. Its computed baffled
terminal and Stokes–Kirchhoff exterior paths produce physical microphone
pressure on the original 96 kHz mechanical clock. The intentionally coarse
geometry and authored tensor demonstrate integration; they are not measured
cane or an accuracy-qualified instrument. Build the additional pack and render:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$modal_demo_dir/tensor-blade-manifest.tsv" \
  --out "$modal_demo_dir/tensor-blade.matpack" \
  --material-state synthetic solid manufactured 0
cargo run --release -p fs-couple --bin music_render -- wind \
  "$modal_demo_dir/paired-tensor-dg.performance" "$modal_demo_dir/paired-tensor-dg.wav" --block 37
```

The copy commands in the run recipe include these source files. Changing only
`branch_33` or `relaxation_time` in `tensor-blade.tsv`, then rebuilding the pack,
changes the normal-strain relaxation, actual reed motion, bore feedback and
external pressure. Geometry and equilibrium elasticity remain declared
independently. The original dense exact propagator still bounds the complete
state: `3 * total_modes + material_histories <= 128`; each source branch retains
six histories at each of the three area points per triangle. Exceeding a budget
refuses the source rather than discarding material directions or blade modes.

## Physical scope and verification

This is a source/renderer integration over flat linear DKT plates, a finite
declared mode basis, either an effective uniform slit/contact or the selected
source-mesh spatial interface, plane-wave duct graphs, and either compact
independent or source-bound mutual infinite-baffle radiation. It does not supply
curved/nonmatching reed-surface collision, shell/solid wet-cane geometry, resolved moving-channel
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
Hereditary regressions separately check modal memory storage and work, positive
material loss, unchanged accepted histories after a refused step, the exact
elastic limit, data-only spectral changes at the microphone, callback identity,
and refusal of unresolved source laws. The regional selector additionally has
cross-mode storage/work tests, agreement with a uniform proportional spectrum,
and data-only regional changes through the actual CLI PCM path. These tests are
authored; native compilation and execution remain required.
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

## Use the same material plates with the DG bore

`paired-dg.performance` selects `frankensim-material-modal-dg-valve-v1`.
It retains both original regional plates and all six requested modes, with
the same material viscosity, slit flow and contact owners. The actual DG
inlet gas and mean pressure/flow impedance enter that retained mechanical
system. Change the material pack or gas record in the source data to rebuild
the coupled physical properties.

After the normalization steps above, render the supplied finite source with:

```bash
cargo run --release -p fs-couple --bin music_render -- wind \
  "$modal_demo_dir/paired-dg.performance" "$modal_demo_dir/paired-dg.wav" \
  --decimate --block 37
```

This fixture declares a 96 kHz mechanical/output source clock and observes
**internal mean inlet pressure**. `--decimate` performs the existing explicit
conversion to 48 kHz. The synthetic bore, plate and material data demonstrate
integration; spatial, nonlinear-time and measured-instrument validation remain
separate work.

The DG schema starts with `audio RATE RATE SAMPLES FULL_SCALE_PA`, requires
`propagator exact-zoh`, and uses `ambient` followed by the shared `dg_tube` or
`dg_bore` records. After `slit_length`, select `observation inlet-pressure` or
`observation baffled-outlet X Y Z RINGS ANGLES FMAX_HZ`. The latter requires
the shared computed DG radiation declaration and uses that reaction source's
original outlet geometry and exterior gas. It does not infer an exterior
from the reed's inlet gas.

Immediately after a baffled observation, the optional shared record
`observer_gas stokes-kirchhoff 4096 1000000 50000000` enables physical classical
gas absorption on each original terminal-to-microphone quadrature path. The
three positive budgets bound distinct exact path lengths, retained source
filter taps and kernel construction terms; their ceilings are 16,384,
4,000,000 and 500,000,000. `lossless` is an explicit comparison, and omission
retains the previous lossless observer. Inlet-pressure output rejects this
exterior-only record.

`paired-dg-chimney.performance` uses this selection for both its bore outlet
and finite chimney. The complete bound exterior gas supplies density, sound
speed, viscosity, thermal conductivity and heat capacity. Changing its source
data rebuilds the original physical path filters while retaining signed
terminal flow, geometry and relative phase. The shared owner also serves the
native DG flute, including its optional mouth. It does not change the reed,
bore or radiation reaction, and reports physical gas, path work, numerical
latency and the classical law's omitted molecular relaxation in the native
sidecar. No species relaxation data is inferred from an instrument or material
name; the existing programmatic molecular law requires its sourced state.

`valve`, `limits`, `contact`, `loss`, `blades`, the complete plate blocks,
optional wetted surfaces, pressure shapes, spatial closure, source forces,
source pressures, and the canonical schedule keep the same record meanings
shown above. This DG route has one explicit upstream pressure supply; it
does not declare the separate upstream tract or moving graph boundaries of
the two-sided schema. Every force/pressure event uses the mechanical clock.

### Finite DG terminal passages and closures

Every source using the shared DG bore grammar can now place an explicit fixed
slit or circular channel, or a rigid closure, between a terminal and the complete exterior
radiation load. This includes the flute, scalar aperture and material-modal
reed paths. After any `dg_side_openings` block and before `dg_radiation`, use:

```text
dg_terminal_channels 2 16
dg_terminal_closed 1
dg_terminal_slit 7 0.004 0.00008 0.001 16 0.1
```

The header is `COUNT MAX_TOTAL_TRANSVERSE_CELLS`. A slit row contains `PORT`,
width, height and length in metres, transverse finite-volume cell count, and
the maximum admitted local Mach number. A closure row needs only its existing
port ID. Port 1 is the axial outlet; additional IDs come from the declared side
openings or finite chimneys. Omitted ports retain their direct connection.
The separate cell allowance sums every slit and circular field; closures require no fluid
states and may use an allowance of zero.

For a circular no-slip passage, replace the slit row with:

```text
dg_terminal_circular 7 0.0001 0.002 16 0.1
```

The fields are `PORT RADIUS_M LENGTH_M RADIAL_CELLS MAX_MACH`. Annular finite
volumes resolve the radial velocity profile, including axis symmetry and the
outer no-slip wall, through the same momentum and energy owner as the slit.
The radius must not exceed one tenth of the length; at least two radial cells
and the same Mach limits apply. Both rows describe **additional passage
geometry** between the duct's pressure plane and its exterior load. Do not
count a chimney's already resolved length or viscous loss again here.

Each passage uses the actual terminal's original gas, including a finite chimney's
own gas record. Density and viscosity therefore change the channel inertia,
velocity profile and dissipation when the source data change. Its original
unsteady transverse momentum solver participates in the same interval solve
as every DG duct and the full radiation matrix. A closure imposes exactly zero
flow while preserving the exterior matrix and its histories. Source feedback,
stored energy, viscous loss and the microphone all use the accepted flow.
The original channel's numerical loss is reported separately from viscous loss.

The admitted slit has fixed parallel no-slip walls, at least two transverse
cells, height no more than one tenth of either width or length, and a declared
Mach ceiling in `(0, 0.3]`. Its area cannot exceed a known attached duct area.
The declared exterior circular surface remains unchanged; each passage is a finite
series channel, not a replacement radiation surface. Moving-pad volume
work, lateral-wall corrections, entrance/exit shedding and turbulent losses
are outside this reduction. Changing a pad gap during playback needs the
corresponding moving-boundary operator and is not exposed by these records.

`../flue-dg-slit-chimney.performance` combines this shared channel with a finite
DG chimney, three distinct source gas records and mutual mouth/terminal
radiation. To close that chimney, replace its slit row with
`dg_terminal_closed 7` and set the header's cell allowance to zero. The source
remains a manufactured integration fixture; no measured flute or reed accuracy
is asserted.

## Ensemble playback

The existing `music_render ensemble --valve` input accepts both material-modal
schemas as well as the original plate-valve schema. Each input keeps its own
physical renderer, material histories, controls and microphone. Supply the
ensemble's explicit `--full-scale-pa`; individual source PCM scales do not
become mix gains. Sources must declare the same physical duration.

Rate conversion follows the source's pressure output clock. The original
two-sided source already produces 48 kHz frames from its faster mechanics;
`--decimate` is required when a selected producer, such as the 96 kHz DG
fixture, still needs conversion to the ensemble's 48 kHz output. Numerical
receiver-filter and outer-decimator delays are aligned together. Physical
distance/sound-speed travel remains in each original receiver, and fractional
final-output alignment delays refuse rather than being rounded.

Actual exterior microphones must agree on the complete bound listener gas.
A selected internal DG pressure trace is identified as internal in the output
metadata; combining it with other traces does not make it an exterior
microphone. The ensemble superposes independent physical pressure sources and
does not introduce mechanical coupling between instruments.

# Three-dimensional solid material source

`solid-hereditary-dg.performance` exposes the existing complete 3D
material-valve owner with its actual tensor relaxation histories through the
same native DG source and Lab runner. `solid-dg.performance` selects the
explicit frozen-frequency image of the same source as a comparator.
These are manufactured 1 cm
tetrahedron with two free displacement coordinates, positive equilibrium
elasticity and a nonproportional tensor-Maxwell branch. It is an integration
fixture, **not measured cane, tissue or a qualified instrument specimen**.
The remaining structural family is chosen by data: replace the volume mesh,
supports, material assignments, axes, source claims and frequency window.

Normalize its source data using the same material normalizer as the plate
examples, selecting `solid-valve-manifest.tsv` and writing
`solid-valve.matpack` alongside the performance. The source clock is 96 kHz;
use the ordinary native mechanical render and explicit decimation for playback.
The computed baffled terminal and exterior gas observation retain their
existing physical owners. A density, elastic modulus or tensor-branch change
rebuilds the actual 3D mass/stiffness/full viscosity and changes the coupled
aperture, bore and exterior pressure.

The `solid` block replaces `loss` and `blades`, immediately after `contact`:

```text
solid tet4 SOURCE_TOKEN
solid_window LOW_HZ HIGH_HZ MAX_SELECTED_MODES MIN_RETAINED_COMPLIANCE
solid_approximation REFERENCE_OMEGA MAX_STORAGE_ERROR MAX_LOSS_ERROR MAX_MODAL_LOSS MAX_DECAY_ERROR MAX_FREQ_CENTS MAX_CLUSTER
solid_limits MAX_NODES MAX_TETS MAX_FREE_DOFS MIN_SCALED_JACOBIAN MAX_INTERNAL_STRAINS
solid_materials COUNT
solid_material PACK isotropic|orthotropic T_K QUERY_OMEGA STRAIN_LIMIT_PROPERTY
solid_axes R11 R12 R13 R21 R22 R23 R31 R32 R33
solid_branches COUNT
tensor_branch TIME_PROPERTY C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66
solid_nodes COUNT
solid_node X_M Y_M Z_M FIX_X FIX_Y FIX_Z FORCE_X FORCE_Y FORCE_Z
solid_elements COUNT
solid_tet N0 N1 N2 N3 MATERIAL_SECTION
solid_orifice OPEN_X OPEN_Y OPEN_Z WIDTH_X WIDTH_Y WIDTH_Z OPPOSING_PLANE_M
solid_face_nodes COUNT
solid_face_node NODE
solid_pressure_faces COUNT
solid_pressure_face N0 N1 N2 SIGN
initial COUNT
state DISPLACEMENT_M VELOCITY_M_S
```

For complete material memory, replace the `solid_approximation` row with:

```text
solid_memory relaxed|unrelaxed MAX_ANGULAR_STEP MAX_BASIS_SCALARS MAX_POINTS MAX_PROJECTION_TERMS MAX_FACTOR_TERMS MAX_BASIS_ERROR
```

This selects the original equilibrium modes and every source tensor branch.
The shared solid projection uses its original strain matrices, integration
volumes and consistent mass. Source terms with exactly equal relaxation time
share the existing orthogonal energy-preserving channel compression; no
material rate is fitted or averaged. The same exact aperture owner evolves
the resulting histories with pressure, swept volume, contact and physical
force. Maxwell energy and heat appear in its existing material ledger, and no
frequency-frozen viscosity is added for those same branches.

`solid_window MAX_SELECTED_MODES` also caps the hereditary projection
coordinates. The additional limits bound actual free-DOF basis entries,
material quadrature points, projection work and equal-time factor work.
`MAX_BASIS_ERROR` checks mass/equilibrium-stiffness projection identities,
not omitted-mode accuracy. The source clock must resolve both the declared
instantaneous angular-frequency ceiling and `dt/tau <= 2`. Relaxed and
unrelaxed initial memories are distinct physical initial conditions; they
never reset when the output callback size changes. The source frequency band
and complete three-dimensional strain domain remain in force.

Repeat material/axes/branch records per section, and the node, element, face and
state rows according to their counts. Tensor entries and times are **property
names in the pack**, not coefficients authored in the performance. Full
orthotropic equilibrium requires all three Young moduli, Poisson ratios and
shear moduli; no missing out-of-plane data are inferred. All tensors use the
original Mandel convention and the explicit material-to-mesh rotation. Each
declared material section must own at least one original tetrahedron. The same
query resolves density, equilibrium elasticity, strain limit, branch entries
and times; their original validity domains remain in force.

The existing `fs-solid::loads3` owner integrates each selected exterior
triangle. `SIGN=+1` is unit inward pressure; `SIGN=-1` is unit outward
pressure relative to the actual solid. Connectivity winding cannot change
this convention. This is the conjugate reference load: the aperture applies
its **negative** for positive upstream-minus-bore pressure, so an opening-positive
reference map gives a closing blowing-pressure force. The resulting signed nodal pressure map drives **both**
structural traction and swept volume. An interior, duplicate or absent
triangle refuses. Opening-axis/width geometry supplies the mean slit and its
rest gap; `valve REST_GAP` must agree with that geometry. The three force
columns per node are an explicit dimensionless nodal pattern per applied
newton; optional `source_forces` rows address source zero and drive the same
complete mode basis without dividing by mean opening.

Both choices retain every selected mode, including zero-opening modes, and the
entire applicable coupled viscosity matrix. The `solid_approximation` choice admits a **frozen-frequency** K/C image
using the original storage/loss error samples, equilibrium compliance check,
modal bandwidth and exact propagation limits. It does not evolve the complete
hereditary volume history, establish an error bound between frequency samples,
or certify omitted modes. The `solid_memory` choice instead evolves the full
projected hereditary source, while modal/spatial refinement, sampled nonlinear
coupling and within-sample strain extrema still require separate resolution.
It does not reconstruct per-element heat deposition. Native volume geometry currently uses Tet4/P1
surface pressure with mesh refinement; Tet10 input refuses instead of
replacing its curved faces with linear loads. The gap contact is the existing
effective mean-slit law; plate spatial closure, plate pressure traces,
follower pressure and large-strain mechanics are not implied. Material
strain-domain checks operate on the original solid strain matrices.

