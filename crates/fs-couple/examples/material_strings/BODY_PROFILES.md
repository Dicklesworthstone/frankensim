# Material plate and curved shell bodies in the shared string performer

`material_strings PERFORMANCE NEW_OUTPUT.wav --body-profile PROFILE` selects
an explicit rectangular material chart (v1/v2), a curved shear-deformable
shell (v3), or an assembly of separately sourced shell components (v4).
Omitting the flag retains the original construction and material profile.
The performance's `body` record still supplies its pack, dimensions, thickness,
mode count and observer distance. String layers, support prescriptions, torsion,
loads, plucks, bow tracks, shared gas and optional cavity are unchanged.

The rectangular profiles normally use `from_viscous_material_plate_chart`; an
explicit regional-Maxwell selection uses the original elastic chart plus the
shared coupled material-memory owner described below. The curved
profile uses `from_material_shear_shell` and the shared MITC3+ shell owner.
The assembled profile uses `from_material_shell_assembly` and the original
work-conjugate mortar ties. All project material tensors on the SAME eigenbasis as bridge displacement,
force and radiation. Full modal viscous cross terms participate in the shared
string/body/friction solve; the body is not a bank of independent decay filters.
No extra scalar body damping or an elastic-modulus-derived viscosity is added.

## Run the synthetic fixture

From the repository root, with the existing toolchain and pinned dependencies:

```bash
demo="$(mktemp -d)"
cp crates/fs-couple/examples/material_strings/*.tsv "$demo/"
cp crates/fs-couple/examples/material_strings/pluck.material-strings "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wire.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/directional-body-manifest.tsv" --out "$demo/body.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/pluck.material-strings" "$demo/directional.wav" \
  --body-profile crates/fs-couple/examples/material_strings/orthotropic.body
```

The coefficients are synthetic, not measured spruce, maple or any other material.
Change only the body pack or its grain angle to compare results; there is no pitch
normalization, per-mode gain fitting or replacement of the physical string input.
Complete source and model admission still precedes `create_new` output.

## Profile grammar

The fixed-order UTF-8/LF records, with no ignored fields, are:

```text
frankensim-string-body-v1
mesh NX NY
elastic orthotropic12
grain ANGLE_RAD
support simply-supported
pretension FORCE_PER_LENGTH_N_M
window LAMBDA_LOW LAMBDA_HIGH
bridge IX IY
```

`NX`/`NY` are interval counts in `2..=12`. Node `(IX,IY)` is at
`(IX*length/NX,IY*width/NY)` and must be interior. This is a physical nodal point
attachment, not a supplied modal vector; all strings currently use this one
footprint. The separate string records still specify each bridge's stiffness,
damping and station on its complete span including afterlength.

Elastic models are `isotropic`, `orthotropic12` or `orthotropic21`; the ordered
plane and grain angle are passed unchanged to the original material chart.
Support is `simply-supported` or `clamped`. Pretension is nonnegative conservative
plate force per unit length and does NOT acquire material bending loss. The
spectral window contains squared angular frequencies in rad²/s², not note names
or frequencies in Hz. Existing modal/Nyquist/source-band admission still applies.
The finite spatial basis and coupled dashpot time discretization require their
usual mesh, mode-count and time-step refinement.

Every selected body pack supplies `density` and either `young_modulus` plus
`poisson_ratio`, or `young_modulus_1`, `young_modulus_2`, `poisson_ratio_12` and
`shear_modulus_12`. The default tensor loss also requires ALL six reduced in-plane viscosities in Pa s:
`plane_stress_viscosity_11`, `_22`, `_66`, `_12`, `_16`, `_26`. The engineering
shear convention is the original shared material owner's convention. Signed
couplings and explicit zero entries are preserved; the tensor must be positive
definite or explicitly all zero. Missing entries never mean zero. The legacy
`kelvin_voigt_bending_viscosity` cannot substitute for this complete tensor.
All original state, units, source-band, coefficient and matrix checks apply.

### Source-driven material memory in a rectangular body

After the last `bridge` row of either a v1 or v2 profile, an optional loss row
selects the body's constitutive law. An absent row or `loss tensor` keeps the
existing full Kelvin–Voigt tensor. A Maxwell selection is:

```text
loss regional-maxwell relaxed 32 1000000 1 maxwell_modulus maxwell_time
```

The fields are initial history (`relaxed` or `unrelaxed`), maximum internal
material states (`0..=256`), maximum cold projection work, branch count (`1..=8`),
then exactly one excess-modulus property and relaxation-time property per branch.
Every region resolves those keys from its own original pack, in Pa and s. This
selection replaces the tensor-viscosity requirements; unselected viscosity data
are not applied in parallel. A zero excess modulus is explicit lossless material
data. A missing property, incompatible units, unresolved relaxation time,
insufficient budget or unsupported source band refuses construction.

The actual region bending forms, including cross-mode terms, use the same DKT
eigenbasis as all bridge and pressure maps. Within each region each excess
modulus scales that region's equilibrium tensor. Regions may have different
spectra and time constants; this is not a general three-dimensional tensor Prony
reduction. Original conservative prestress does not relax. Every material state
contributes storage, loss and reaction through the common string/body solver.
Gas memory, longitudinal waves, torsion and the cavity are constructed before
this attachment, and the original bow/force/PCM pipeline retains it afterward.

Use `regional-relaxing.body` with the same commands and directional pack above
to exercise this path. Its two added `maxwell_modulus` and `maxwell_time` source
properties are deliberately synthetic, just like the other fixture coefficients.
Changing a region's source values changes the coupled response without changing
physics code. The complete source band must cover the coupled small-signal
frequency range, including its instantaneous elastic bound. A full general
memory model for v3/v4 shells and the mechanical-only `prepared-linear` backend
is not supplied by this selector; an explicitly selected prepared backend
refuses before reading material files. Existing exact-state capacity and spatial,
modal and time refinement requirements still apply.

Native regression sources cover regional material substitution reaching the
string/bridge/PCM output, partition-independent histories, zero-law parity and
composition with bonded strings, both transverse planes, torsion, gas memory,
cavity and axial waves. Syntax checks passed in the implementation environment;
these new native regressions still require execution.

Profile input is bounded to 64 KiB; the original source-pack and complete-network
limits remain. Viscous projection is capped at 16,777,216 operations and 16 body
modes. The dashpot factor uses a 1 N s/m coordinate scale; this changes neither
the projected damping matrix nor material coefficients. The exact profile bytes,
original packs, regional material states, tensors and projected damping remain
attached to the performer. No live source-state mutation is implied.

The v1/v2 model is a fixed-state, flat, linear rectangular DKT plate. It does
not resolve arching, a soundpost, bass bar or glue interface. Pressure uses the
inherited compact baffled body/neck observation; it does not include full
exterior radiation or radiation back-loading. The v3 shell model below resolves
an explicitly authored curved midsurface.

Native regression target:

```bash
cargo test -p fs-couple --example material_strings body::tests -- --test-threads=1
```

## Regional packs and separate string attachments

A v2 profile assigns material data to nonoverlapping grid-cell rectangles on ONE
shared plate mesh. It also supplies one physical bridge node per material string,
in the performance's string order. It does not create independent body voices or
recombine their PCM output. Source values are never averaged across regions.

```text
frankensim-string-body-v2
mesh 4 4
support simply-supported
pretension 0
window 1 100000000
regions 2
region left @body orthotropic12 1 0.17 0 2 0 4
region right @body orthotropic12 1 -0.23 2 4 0 4
bridges 1
bridge 2 2
```

This is the supplied `regional.body` fixture. It is runnable with the same
performance and directional material pack as the v1 example: replace the final
profile argument with `crates/fs-couple/examples/material_strings/regional.body`
and choose a new output name.

Each region record is:

```text
region NAME PACK ELASTIC_MODEL THICKNESS_RATIO GRAIN_RAD X0 X1 Y0 Y1
```

The half-open CELL rectangle `[X0,X1) x [Y0,Y1)` owns both original triangles in
each named cell. Every cell must belong to exactly one region; gaps, overlaps,
empty/out-of-grid rectangles and repeated names refuse. There are at most 16
regions. `THICKNESS_RATIO` multiplies the performance's physical body thickness;
no material stiffness is adjusted to disguise a thickness change. Density,
elasticity and all six viscosity coefficients are resolved from that region's
original pack, at the unchanged performance material query point.

`@body` references the pack in the original performance `body` record. At least
one region must use it, so selecting a profile cannot silently discard that
authored material. Other regions may select another normalized pack, for example
`right-body.fsmatpk`; paths are relative to the PERFORMANCE directory, not the
profile directory, and retain the same no-traversal policy. Repeated paths load
once and share the exact decoded pack; different source selections keep their
original material identities. A missing tensor entry in even one region is a
refusal, not permission to borrow the adjacent region's viscosity.

The v2 `bridges` count must equal the actual string count (1..=8), including
sympathetic/unplayed strings. Each `bridge IX IY` must be an interior node.
Different nodes project onto the same actual body eigenbasis and can excite
different signed mode combinations. Repeated nodes are allowed and retain
separate connections. For a two-string performance, for example, use:

```text
bridges 2
bridge 1 2
bridge 3 2
```

The existing spatial-port constructor accounts for these forces, their body
reaction, cavity volume response, full material damping and pressure together.
The first bridge is only the reducer's reference drive; it does not add a duplicate
physical load. Each string retains the stiffness/damping/station of its original
bridge declaration, and surface torsion and preload/release remain unchanged.

Regions are perfectly connected parts of one linear plate field. This is not a
model of a glue joint, separated component, three-dimensional bridge contact or interlayer slip.
The same mesh, modal, time-step and compact-radiation limitations still apply.

Five additional native regressions cover same-material subdivision, changing
only one region's loss, different spatial bridge maps and sympathetic motion,
callback equivalence, original geometry indexing and source/coverage refusals.
They are included in the existing `--example material_strings` test command.

## Reciprocal rotational bridges

A performance can additionally connect sourced string twist to scalar body
rocking without supplying a modal rotation vector. After its last complete
string/section record and before `cavity`, include:

```text
rotational-bridges 1
rotational-bridge 0 0.8 0 1 0.001 1e-7 0
```

This selects string 0 at station 0.8, the ordered pair of body bridge footprints
belonging to strings 0 and 1, angular stiffness 0.001 N m/rad, damping
1e-7 N m s/rad, and stress-free relative angle zero. The selected string must
declare positive twist-mode count and supply its original shear material data.
The body profile must be v2 or v3 with distinct physical nodes for the chosen
pair. V3 additionally requires identical authored translation directions at the
two nodes and derives its rotational lever from their three-dimensional chord.

If the retained displacements along the common force direction are
`w_positive` and `w_negative`, body rotation is their difference divided by
the physical rotational lever. For the flat v2 plate this is the Euclidean
node separation. V3 uses `|chord × force_direction|`, with the full physical
3D chord; dividing by the full chord length would mis-scale an oblique
attachment. Parallel chords refuse because they measure extension rather than
rotation. The lever is recalculated when geometry changes. Swapping the
ordered pair reverses the scalar coordinate; its axis follows the chord/force
cross product. The string joint still declares the matching scalar orientation,
without an inferred complete 3D string pose. This coordinate is exact for rigid
rocking and a secant approximation for a deforming body, not a resolved local
slope claim.

The body angular row is derived from its original spatial bridge projections.
The string angular row uses its sourced polar inertia and its actual outer
radius. One angular spring/dashpot then supplies opposite string/body torques in
the existing joint mechanical solve. Its stored energy and dissipated work enter
the same ledger as transverse bridge reactions, all material body damping cross
terms and the cavity. Twist can reach body pressure while bows are released;
no direct torsional pressure gain or separate audio voice is added.

Using the `$demo` directory and directional body pack built above:

```bash
cp crates/fs-couple/examples/material_strings/torsional-bridge.material-strings "$demo/"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/torsional-bridge.material-strings" "$demo/torsional-bridge.wav" \
  --body-profile crates/fs-couple/examples/material_strings/torsional-bridge.body
```

The supplied fixture has two strings, one shared regional body, a cavity and a
surface preload released at sample zero. All bow loads stay zero. Its 0.2 m
footprint baseline and explicitly synthetic joint law demonstrate the input
path; they are not fitted violin/cello bridge mechanics. The normal render
report includes each declared angular law and its derived physical separation,
while the original performance/profile bytes and material sources stay attached
to the performer.

The optional block works with v1, v2 and v3 performance files. At most 16 joints
are accepted, within the existing complete-network connection/work limits.
Omitting it preserves the old construction. Missing profiles, v1 single-footprint
profiles, absent string indices, coincident or unrepresentable baselines,
unsourced twist and nonpassive laws refuse before a WAV is created. Native tests
exercise file admission, version compatibility, geometry refusals before pack
loading, and actual released-body pressure/PCM from an explicitly preloaded
rotational spring while retaining the original full material damping.

## Gas memory with a regional body and physical contacts

The performance's optional `cylinder-gas-loading X_MIN X_MAX INTERVALS` record
also works with these uniform and regional body profiles. Place it after the
rotational-bridge block, when present, and before `cavity`. The shared string
owner derives gas loading from each actual outer section and the existing gas.
Fluid coordinates follow every physical string/body/neck component, so original
spatial footprints, angular joints, complete material viscosity and pressure
projections retain their physical indices and sources.

The same construction accepts `--stop-profile` before static preload. The native
composition regression uses two layered two-plane strings, a regional viscous
body, one torque joint, a cavity, one receiver-bound stop and six fluid branches
per transverse mode: 33 complete coordinates. It also changes one region's
viscosity without changing the fluid law or the other materials.

See [README.md](README.md#causal-gas-loading-of-the-moving-strings) for resolution
semantics, the runnable source fixture, continuum assumptions and the current
42-mode connected-island limit. The small regression is a coupling check; body
mesh, modal, fluid-spectral and time refinement remain independent requirements.

## Curved bodies with material transverse shear

A `frankensim-string-body-v3` profile retains an explicit 3D midsurface,
unit nodal directors, homogeneous supports, material regions and one physical
force direction per string. This connects the shared six-DOF MITC3+ shell to
the same reciprocal string/body/cavity solve. Curvature, membrane stiffness,
bending stiffness, transverse shear, rotary inertia and full material viscosity
are derived from that geometry and its original material packs.

The complete input is [arched.body](arched.body). Its fixed-order record layout
is shown below; the repeated node, triangle, support and bridge records are
abbreviated here:

```text
frankensim-string-body-v3
director-min-dot 0.9
mass consistent
drilling 0.001
window 1 100000000
projection-work 67108864
regions 1
region body @body orthotropic12 1 0.17 0.8333333333333334 1 0 0
nodes 25
node X_FRACTION Y_FRACTION Z_M DIRECTOR_X DIRECTOR_Y DIRECTOR_Z
triangles 32
triangle NODE_I NODE_J NODE_K REGION_NAME
supports 16
support NODE 111111
bridges 1
bridge NODE DIRECTION_X DIRECTION_Y DIRECTION_Z
```

Each region has nine fields after `region`:
`NAME PACK MODEL THICKNESS_RATIO GRAIN_RAD SHEAR_CORRECTION AXIS_X AXIS_Y AXIS_Z`.
Pack paths follow the v2 rules, including retaining at least one `@body`
region. Every triangle belongs to exactly one nonempty named region.
Thickness ratio multiplies the performance's body thickness.

The explicit unit global grain axis is projected onto each physical facet;
the additional grain angle rotates within that plane. Grain therefore does not
follow arbitrary triangle-edge numbering. A near-normal grain projection refuses.
The axis names the first selected section axis: material axis 1 for
`orthotropic12`, material axis 2 for `orthotropic21`. The original ordered
material interpretation is retained. Orthotropic packs additionally supply `shear_modulus_13` and
`shear_modulus_23` in Pa. The isotropic model derives its elastic shear from
the same Young's modulus and Poisson ratio.

In addition to all six plane-stress viscosity entries, every selected pack
supplies `transverse_shear_viscosity_13`, `transverse_shear_viscosity_23` and
`transverse_shear_viscosity_13_23` in Pa s. Zero cross coupling must be explicit.
The shared owner validates both viscous tensors, rotates them with their original
material axes, and projects them using the SAME elastic bubble condensation and
physical eigenvectors as stiffness, force and volume response. It retains the
full modal damping matrix. Missing shear or viscosity data never means borrowing
another material's coefficients.

Node x/y fractions lie in `[0,1]` and scale with the performance's body length
and width; z is in metres. Directors and force directions must already be unit
vectors; the parser does not normalize or smooth them. The shared midsurface and
director owners check actual geometry, topology and facet alignment. Support
masks use six global DOFs `(u,v,w,rx,ry,rz)`, with 1 fixed and 0 free.
Each bridge is a nodal translation along its authored direction and must retain
an unconstrained translational component. The bridge count equals the actual
string count, including unplayed strings. Source bytes, geometry, packs,
physical eigenvectors and force maps remain attached to the performer.

The file admits at most 2048 nodes and 4096 triangles per component, 16 regions, eight string
attachments and 128 retained shell modes, within 1 MiB (flat plate profiles keep 64 KiB and 16 modes). `mass` explicitly selects
`consistent` or `lumped` inertia. `drilling` supplies positive numerical
drilling stabilization; `director-min-dot` declares the incident-facet alignment
floor in `(0,1]`. The spectral window contains squared angular frequencies.
`projection-work` declares a separate cold shell-projection budget, at most
34,359,738,368 terms (the count is `triangles * (10^6 + 324 modes^2) + 6 nodes modes` plus
force and pressure projections, about 2.6e10 for a cap-sized component). It does not enlarge the existing complete-network limits.
These discretization choices require mesh, mode, stabilization and time-step
refinement for the intended specimen.

The synthetic arched fixture uses the same directional body pack and a layered,
two-plane, torsion-bearing string with a cavity and causal gas loading:

```bash
demo="$(mktemp -d)"
cp crates/fs-couple/examples/material_strings/*.tsv "$demo/"
cp crates/fs-couple/examples/material_strings/cylinder-flow.material-strings "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wire.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/sleeve-manifest.tsv" --out "$demo/sleeve.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/directional-body-manifest.tsv" --out "$demo/body.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/cylinder-flow.material-strings" "$demo/arched.wav" \
  --body-profile crates/fs-couple/examples/material_strings/arched.body
```

The body can be bowed, plucked, preloaded, connected through rotational bridge
chords and used with the existing stop and gas profiles. Its observer integrates
normal velocity over the actual shell facets before the shared compact baffled
pressure projection.

This model is an unloaded, fixed-state, linear midsurface. It does not yet
provide a complete assembled violin/cello/bass corpus, prestressed shell
equilibrium, separately connected ribs, resolved soundpost and bass-bar solids,
distributed glue mortars or full exterior acoustic back-loading. The fixture's coefficients are synthetic;
it establishes the data-to-pressure path, not measured whole-instrument fidelity.

## Separate exterior and cavity pressure surfaces

A v3 profile can declare which original triangles participate in exterior
radiation and which change the enclosed air volume. Place this optional block
after the supports and any interfaces, immediately before `bridges`:

```text
pressure-surfaces 2 2
exterior-facet 10 positive
exterior-facet 11 positive
cavity-facet 10 positive
cavity-facet 11 positive
```

The header gives the exterior count followed by the cavity count. Each count
is in `0..=TRIANGLE_COUNT`; the corresponding rows follow in that order.
Triangle indices are zero based and refer to the original mesh, independently
of material regions. `positive` uses its original winding normal and `negative`
reverses that normal. For the cavity row, the supplied normal points **outward
from the cavity**, so positive displacement increases cavity volume. A wall
bounding both fluid domains normally appears once in each list with the same
outward-from-cavity sign. Repeated triangles within one role refuse.

Only the listed facets enter each role. Internal mechanical members can be
excluded from exterior monopole radiation and net enclosed volume while keeping
their stiffness, mass, shear deformation, full material damping, bridge reactions
and material interfaces. Both independent signed volume rows use the same
physical modal basis and shared P1 pressure operator. The cavity row drives the
existing collective air spring and reciprocal body force; the exterior row
drives the existing compact pressure observer. Changing receiver participation
alone preserves mechanical history.

An explicit zero count means zero participation, including when the other
role remains active. `pressure-surfaces 0 0` retains a mechanical body with no
selected compact acoustic coupling. Omitting the entire block preserves the
original whole-midsurface selection for both roles.

These are oriented **midsurface** pressure declarations. They do not extrude
physical skins, infer cavity volume, certify an enclosed/watertight air domain,
or establish that an omitted component has zero dipole or higher radiation.
The performance continues to supply its explicit cavity geometry and gas.
The existing compact observer and lumped cavity retain their approximation
limits; separately resolved interior/exterior fluid meshes require those owners.

## Compliant holding and bonded interface footprints

A v3 body can now include an optional material interface block after its last
support and before its bridges. [held-arched.body](held-arched.body) adds an
explicit holding footprint to the same curved body, layered string, gas and
cavity example above:

```text
interfaces 1 2 1048576
interface holding holding.fsintpk 293.15 1 normal_stiffness_per_area normal_viscosity_per_area 0.0001 0 0 1
side region body 2
point 11 0.5 0.01 0 0
point 12 0.5 0.01 0 0
side fixed
bridges 1
bridge 12 0 0 1
```

The header is `interfaces COUNT MAX_TERMS MAX_WORK`. Each interface record is:

```text
interface NAME PACK TEMPERATURE_K OMEGA_RAD_S STIFFNESS_PROPERTY VISCOSITY_PROPERTY AREA_M2 DX DY DZ
```

It is followed by exactly two ordered sides. A moving side is
`side region NAME POINT_COUNT` followed by `point NODE WEIGHT OX OY OZ` rows.
A fixed side is `side fixed`; it explicitly sets zero counterpart motion and
retains that ordered surface's material identity from the interface card. At
least one side must move. Reversing the surface roles requires corresponding
source data, not an automatically reordered card.

Each footprint is a weighted physical application point. Weights are positive,
sum to one and name distinct nodes belonging to that source region. Offsets
are global metres from the midsurface nodes. Translation and rotation both
participate: at a point with offset r, displacement is u + theta cross r.
The signed difference of the two weighted displacements along the authored
unit global direction is the interface extension. The same row returns the
equal-and-opposite forces and moments. Two moving sides must represent the
same physical application point, so common rigid motion does not strain a
bond. Independent component meshes and continuum mortar interpolation are
not inferred from region names.

The original normalized interface pack supplies stiffness per area in N/m3
(`Pa/m`) and viscous resistance per area in N s/m3 (`Pa*s/m`) under the
property names in the record. Multiplying by area produces one N/m spring
and one N s/m dashpot. Both coefficients must be explicitly sourced, finite
and nonnegative; at least one is positive. A pure dashpot with sourced zero
stiffness is supported. These are fixed-state constant coefficients with an
explicit frequency validity interval; all retained shell modes must fit the
intersection of the original bulk and interface source bands.

The shared MITC owner adds the complete spring outer product before the body
eigensolve. The existing reciprocal solver applies interface damping on that
same physical eigenbasis alongside full bulk material viscosity, bridge
reactions and cavity coupling. Spring stiffness is not added again at runtime.
The interface's stored energy and dissipated work participate in the shared
mechanical accounting. Changing only the interface pack's stiffness therefore
changes the body modes; changing only its viscosity changes loss and pressure
without retuning the conservative body. No output-stage effect is involved.

Using the temporary directory and three bulk packs from the arched example:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/holding-manifest.tsv" --out "$demo/holding.fsintpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/cylinder-flow.material-strings" "$demo/held-arched.wav" \
  --body-profile crates/fs-couple/examples/material_strings/held-arched.body
```

The original wildcard copy of the TSV files includes
[holding-manifest.tsv](holding-manifest.tsv) and [holding.tsv](holding.tsv).
The source explicitly names the directional body pack's current material state:
`source-pack:synthetic-directional-string-body`, phase
`declared solid model`, process `frozen source-pack realization`, revision 0.
A replacement body with another source identity requires its corresponding
ordered interface source identity as data. The second surface is the explicitly
immobile synthetic support. These manufactured coefficients demonstrate the
complete source-to-pressure path; they are not measured hand, glue or corpus data.

Up to 64 interfaces and 4096 total nodal footprint points are admitted.
`MAX_TERMS` counts points across both sides of every interface, not scalar
degrees of freedom. `MAX_WORK` bounds each shared cold binding, assembly and
projection operation and may not exceed 268,435,456. Every declared interface
reserves a possible damping channel under the existing complete-network limits
before pack I/O. Interface and bulk artifacts share the original 8 MiB per-pack
and 64 MiB total limits; repeated interface paths load once. Exact profile
bytes, original normalized packs, ordered interface states and projected rows
remain attached to the performer.

`supports 0` is accepted when explicit interfaces are present. The shared
solver still requires an admissible eigenproblem; a single scalar support
cannot restrain every rigid mode of a free three-dimensional body. The supplied
fixture retains its original outer supports. These bilateral linear resultants
resolve compliant holding or a finite bonded footprint at the reference state.
They do not supply unilateral separation, frictional slip, a resolved soundpost
solid, prestressed corpus equilibrium or a distributed glue mortar.

## Assembled material shell bodies (v4)

A v4 profile describes independently oriented structural parts and their actual
connections. Each component retains the v3 regional material, mesh, director,
support and optional material-interface inputs. The original MITC owner
assembles every complete component, and the original shell mortar owner joins
their stiffness and mass matrices **before one body eigensolve**. All strings,
full material loss, cavity pressure and exterior observation use those same
assembled modes. There is no per-component modal truncation or mixing of
independently rendered audio.

Use [assembled.body](assembled.body) with the existing `pluck.material-strings`
performance and directional body pack:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/pluck.material-strings" "$demo/assembled.wav" \
  --body-profile crates/fs-couple/examples/material_strings/assembled.body
```

This synthetic fixture has two coplanar material pieces connected along their
common seam. It demonstrates component assembly; it is not a violin or cello
geometry. The first piece is clamped at its outer edge, and the second receives
the string force. Changing only the second component's pack or grain changes
the coupled mechanics without retuning the string or changing physics code.

The record order is:

```text
frankensim-string-body-v4
window LAMBDA_LOW LAMBDA_HIGH
projection-work MAX_WORK
components COUNT
component NAME
director-min-dot MIN_DOT
mass consistent
drilling ALPHA
regions COUNT
region ...
nodes COUNT
node ...
triangles COUNT
triangle ...
supports COUNT
support ...
interfaces ...
pressure-surfaces N_EXTERIOR N_CAVITY
exterior-facet TRIANGLE positive
cavity-facet TRIANGLE positive
end-component
tie-budget MAX_SEGMENTS MAX_SLAVE_DOFS GAP_TOLERANCE_M ANGLE_TOLERANCE_RAD
edge-ties COUNT
edge-tie MASTER_COMPONENT SLAVE_COMPONENT OFFSET ROTATIONS N_MASTER N_SLAVE
master-node NODE
slave-node NODE
face-ties COUNT
face-tie MASTER_COMPONENT SLAVE_COMPONENT OFFSET ROTATIONS N_MASTER N_SLAVE
master-triangle TRIANGLE
slave-triangle TRIANGLE
bridges COUNT
bridge COMPONENT NODE DIRECTION_X DIRECTION_Y DIRECTION_Z
```

Repeat each component and its records in declaration order. The `region`,
`node`, `triangle`, `support`, optional `interfaces`, and pressure records have
their unchanged v3 meanings. Components inherit the assembly's window and
projection-work allowance; they do not declare their own window or bridges.
Individual components may have `supports 0` without an interface. The complete
assembly must still pass the original mechanics and spectral checks.

Every component must explicitly declare both fluid roles. `pressure-surfaces
0 0` makes an internal structural member mechanically active without inventing
an exterior or cavity pressure surface. Exterior and cavity rows are summed
independently over the complete body. They are not inferred from part names.

An edge tie supplies all ordered master nodes followed by all ordered slave
nodes. A face tie similarly supplies all master triangles followed by all slave
triangles. Indices are local to their named components. `OFFSET` is
`midsurface`, `positive-skin` or `negative-skin`; skin attachment arms come from
the master's actual resolved thickness. `ROTATIONS` is `full` or
`translations`. These records select the original rigid mortar law. Nonmatching
seams and overlapping patches remain subject to that owner's geometric and
work limits. Prescribed slave coordinates that the tie would eliminate refuse;
the profile cannot silently discard a support.

Each final `bridge` is one physical nodal force footprint, in actual string
order, on any component. Rotational bridges can use footprints on different
components; they retain the same physical chord and perpendicular-lever rule
as v3. The other string controls, material bow interfaces, afterlength, torsion,
preloads and shared gas continue through the original performer.

Profiles take the 1 MiB shell limit and admit at most eight components, 8192 total
nodes, 16384 total triangles, 64 ties and 128 shared body modes. Per-component v3
mesh limits still apply. Material/interface source bytes and every damping
channel remain under the existing complete-performance budgets. At least one
region across the assembly must retain `@body`; other regions can select
different source packs. Repeated bulk paths share the original material cache.

This is unloaded linear component mechanics with rigid ties. It does not add
prestressed assembly equilibrium, compliant glue between components, solid
soundpost contact, a certified closed cavity skin, or resolved interior and
exterior acoustic fields. The inherited compact observer and Helmholtz cavity
remain explicit approximations. The new native regression sources check tied
motion, complete material-to-pressure/PCM propagation, material-only changes
and callback equivalence; they still require execution in the pinned build.
