# Material plate and curved shell bodies in the shared string performer

`material_strings PERFORMANCE NEW_OUTPUT.wav --body-profile PROFILE` selects
an explicit rectangular material chart (v1/v2), a curved shear-deformable
shell (v3), or an assembly of separately sourced shell components (v4).
Omitting the flag retains the original construction and material profile.
The performance's `body` record still supplies its pack, dimensions, thickness,
mode count and observer distance. String layers, support prescriptions, torsion,
loads, plucks, bow tracks, shared gas and optional cavity are unchanged.

V3/v4 profiles can also supply an explicit tetrahedral interior with the
`interior-fe-v1` trailer below. Its gas is the performance's original resolved
atmosphere; its standing waves and uniform compression load the same material
body before the existing force, bow, fluid and radiation attachments.

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

`NX`/`NY` are interval counts in `2..=12` for perimeter support, or `1..=12`
for `clamped-x-min`. Node `(IX,IY)` is at
`(IX*length/NX,IY*width/NY)` and must be an unconstrained transverse node. This is a physical nodal point
attachment, not a supplied modal vector; all strings currently use this one
footprint. The separate string records still specify each bridge's stiffness,
damping and station on its complete span including afterlength.

Elastic models are `isotropic`, `orthotropic12` or `orthotropic21`; the ordered
plane and grain angle are passed unchanged to the original material chart.
Support is `simply-supported` or `clamped` on the perimeter, or
`clamped-x-min` on only the x=0 edge; the other three edges are free and may
carry the declared force footprint. Pretension is nonnegative conservative
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

### Complete three-dimensional tensor memory

The same source-bound tensor law used by material reeds can also drive the
shared string body. All six normal/shear histories per branch and DKT area
point survive the shared plane-stress elimination. This includes effects of
out-of-plane coefficients on in-plane hereditary bending; condensing each
branch separately would lose those effects. The native selection is:

```text
loss tensor-maxwell relaxed STRAIN_PROPERTY COUNT MAX_POINTS MAX_STATES MAX_TERMS
tensor_branch TIME C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66
```

Every token after `tensor_branch` is a property name. Coefficients use the
same Mandel order and Pa units as the canonical `TensorPronyProperties`
owner, and time is in seconds. A shared explicit zero claim may supply several
zero entries. `STRAIN_PROPERTY` supplies the full tensor-shear Frobenius strain
limit. `COUNT` is 1..=8, `MAX_POINTS` is 1..=4096, `MAX_STATES` is 0..=256 and
`MAX_TERMS` is 1..=67,108,864. The original exact propagator can impose a
stricter complete-state limit. No modes, points, branches or zero-coupling
histories are silently removed to meet a budget.

The original chart owns equilibrium stiffness and mass. Tensor memory adds
only the source's declared excess stiffness; selected or authored viscosity
cannot also be counted. Isotropic equilibrium uses E and Poisson ratio;
orthotropic equilibrium requires all three Young moduli, all three Poisson
ratios and all three shear moduli. A two-axis plate source cannot invent the
missing third direction. Each original region's axes and grain angle rotate
the entire law, not just its equilibrium section. The same complete coupled
frequency band and mechanical clock are checked against source support.

The shared coupled owner checks the full history-dependent strain at the
physical outer faces before publishing each accepted step or preload.
Refused strain, energy, source-domain and cancellation trials preserve the
accepted displacement, internal history and clock. Cold restarts and appended
internal coordinates retain that same source admission. The force, pressure,
friction, cavity and nonlinear owners continue to use the original trajectory.

`tensor-body.material-strings` with `tensor-body.body` uses a declared coarse
cantilever body: two original triangles, two body modes, six physical area
points and all 36 history coordinates. It reuses the existing manufactured
tensor material source from the reed fixture; there is no body-specific
material constant or alternate loss law. After building the wire pack above:

```bash
cp crates/fs-couple/examples/material_strings/tensor-body.material-strings "$demo/"
cp crates/fs-couple/examples/material_modal_valve/tensor-blade*.tsv "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/tensor-blade-manifest.tsv" --out "$demo/tensor-blade.fsmatpk" \
  --material-state synthetic solid manufactured 0
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/tensor-body.material-strings" "$demo/tensor-body.wav" \
  --body-profile crates/fs-couple/examples/material_strings/tensor-body.body
```

Changing only `branch_33` or `relaxation_time` and rebuilding that pack changes
the body, reciprocal string reactions and physical pressure. The material and
geometry are manufactured integration data, not measured wood or a qualified
violin/cello/bass. This remains fixed-state homogeneous DKT bending; a complete
tensor shell and large-state propagation remain separate capabilities. Even
the smallest 2×2 perimeter-supported mesh requires 144 histories for one full
branch, so the complete-coordinate selection may exceed the current dense
exact owner and must refuse.

### Exact shared material-history coordinates

The optional `tensor-maxwell-shared` selection carries the same original
point laws and histories in exact equal-time coordinates:

```text
loss tensor-maxwell-shared relaxed STRAIN_PROPERTY COUNT MAX_POINTS MAX_STATES MAX_TERMS MAX_FACTOR_TERMS
```

Branch rows are unchanged. `MAX_FACTOR_TERMS` is a separate positive work
budget, at most 67,108,864. `MAX_STATES` still caps evolved coordinates and
all original exact-network limits remain. `tensor-body-shared.body`, paired
with the same `tensor-body.material-strings` and material packs above, refines
the original cantilever from two to eight triangles. It retains 24 original
material points and reconstructs all 144 original histories while keeping
the same declared 36-state ceiling:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/tensor-body.material-strings" "$demo/tensor-body-shared.wav" \
  --body-profile crates/fs-couple/examples/material_strings/tensor-body-shared.body
```

Only bit-identical constrained relaxation times share coordinates. For one
time family, the original factor satisfies `A = Q R` and `Q^T Q = I`;
the source histories obey `z = Q y`. Relaxed `z = A q` and unrelaxed `z = 0`
initial conditions stay in this invariant subspace, so forces, storage and
heat are unchanged. Every QR direction survives, including rank-deficient
directions, with explicit original-factor and orthogonality residual checks.
No nearby times are rounded together and no mechanical modes are removed.
Every source strain and normal traction still uses the reconstructed
original point history. Independently authored original histories outside
this subspace are not accepted by this selection. This is a coordinate
change within the declared plate model, not additional physical validation.

The complete-coordinate route remains available for direct comparisons.
The focused native regressions compare actual pressure, original point
histories, normal traction, energy and callback partition, then exercise the
refined mesh at the unchanged state cap and its material/work/domain refusals.

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

A pack may also state `constant_loss_reference_omega` in rad/s: the angular
frequency at which its viscosities are its loss factors (viscosity = loss
factor x modulus / omega_ref). That region's viscosity is then projected on its
own elements and read at every retained mode's own frequency, so each mode's
damping carries the pack's loss factors instead of growing with frequency
(`fs_couple::shear_shell::constant_loss`). The anchor's claim domain and the
law's 50 Hz-10 kHz band clip the region's frequency support. A pack without it
is Kelvin-Voigt exactly as before.

Node x/y fractions lie in `[0,1]` and scale with the performance's body length
and width; z is in metres. Directors and force directions must already be unit
vectors; the parser does not normalize or smooth them. The shared midsurface and
director owners check actual geometry, topology and facet alignment. Support
masks use six global DOFs `(u,v,w,rx,ry,rz)`, with 1 fixed and 0 free.
Each bridge is a nodal translation along its authored direction and must retain
an unconstrained translational component. The bridge count equals the actual
string count, including unplayed strings. Source bytes, geometry, packs,
physical eigenvectors and force maps remain attached to the performer.

The file admits at most 3072 nodes and 6144 triangles per component, 16 regions, eight string
attachments and 256 retained shell modes, within 4 MiB (flat plate profiles keep 64 KiB and 16 modes). `mass` explicitly selects
`consistent` or `lumped` inertia. `drilling` supplies positive numerical
drilling stabilization; `director-min-dot` declares the incident-facet alignment
floor in `(0,1]`. The spectral window contains squared angular frequencies.
`projection-work` declares a separate cold shell-projection budget, at most
2,199,023,255,552 terms (the count is `triangles * (10^6 + 324 modes^2) + 6 nodes modes` plus
force and pressure projections, about 1.4e11 for a cap-sized component; a v4 assembly adds the owner's
conservative dense mortar bound `(6 total nodes)^2 (modes + bridges)`, about 2.0e12 in all at the v4 caps).
It does not enlarge the existing complete-network limits.
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

### Source-bound varnish and other measured films

V3 shells and every V4 component accept an optional block **after the last
`region` and before `nodes`**:

```text
coatings 1
coating body varnish-oil-copal-lammlein2021-s1.fsmatpk positive 0.000020 1000 0.30 0.001 0.35 0.25
```

The nine fields after `coating` are `REGION PACK SIDE THICKNESS_M REFERENCE_HZ
POISSON_RATIO STRAIN_LIMIT CONDITIONING_RH MAX_STORAGE_RELATIVE_DRIFT`.
`SIDE` is `positive` or `negative` relative to the original shell directors.
Repeated rows on one region/side stack outward in declaration order, so a
ground and outer varnish can have distinct sources and thicknesses. An absent
block or `coatings 0` retains the uncoated path. Each component admits at most
64 films; the compiler also bounds all retained film/facet records and charges
each distinct film's loss projection to the original projection budget.

The measured-film schema supplies density and `varnish_e0`, `varnish_einf`,
`varnish_alpha`, `varnish_beta`, `varnish_tau`. The source query retains the
performance temperature as typed absolute kelvin, reference frequency as typed
cyclic Hz, the authored conditioning RH, `source_fitted_response=1`, and
`source_density_sample_matched=0`. The published group density is not claimed
to be paired with an individual DMA sample. Historical stratigraphy alone
cannot satisfy these mechanical requirements. See
[`VIOLIN_VARNISHES.md`](../../../../data/matdb/seed-v1/VIOLIN_VARNISHES.md).

The isotropic, real, constant Poisson ratio is an explicit closure of the
measured uniaxial film spectrum. `STRAIN_LIMIT` declares the constitutive
small-strain envelope; this linear shell route does not yet monitor coating
strains during a nonlinear render. The example's values are assumptions,
not measured historical properties.

The shared laminate and MITC owners integrate complete `A/B/D` stiffness,
`I0/I1/I2` inertia and the film's offset viscous coupling. The reference remains
the wood midsurface: the original joint attachment coordinates and cavity
geometry are retained. Film thickness adds physical mass and does not consume
a fixed-mass substrate allowance. Perfect bonding is assumed. Penetrating oil,
curing stress, moisture transport, damage and a displaced exterior acoustic
boundary are not modeled by an added surface layer.

By default, storage is frozen at `REFERENCE_HZ`. Every retained body mode must satisfy
`abs(Eprime(reference)-Eprime(mode))/Eprime(mode) <= MAX_STORAGE_RELATIVE_DRIFT`,
with an authored tolerance in `[0,1)`, and all original material frequency
bands. The reference itself must be in both film and wood source bands. Each
film's modal viscosity instead uses its full source spectrum `Edoubleprime(omega)/omega`.
After physical projection, positive diagonal congruence scales its full matrix:
diagonal modal losses read the measured law at their own frequencies and
off-diagonal terms use the geometric mean of those factors. Passivity and
cross-mode loss survive. Substrate loss keeps its own existing law, so wood's
constant-loss rescaling cannot accidentally alter varnish loss.

This is a bounded approximation of dispersive stiffness plus source-dependent
modal loss. It is not the full hereditary, frequency-dependent shell
eigenproblem. The material owner separately offers passive tensor-Prony
lowering with independent storage and loss error bounds; that memory is not
silently substituted into this shell renderer. The reproducible section study
in [`VIOLIN_VARNISH_PHYSICS.md`](../../../../docs/VIOLIN_VARNISH_PHYSICS.md)
evaluates the exact source spectrum and quantifies the mass/stiffness/loss
competition.

For table conversion, the options file accepts
`coating COMPONENT PACK SIDE THICKNESS_M REFERENCE_HZ POISSON_RATIO STRAIN_LIMIT
CONDITIONING_RH MAX_STORAGE_RELATIVE_DRIFT`. A literal component name selects
all its regions; `@exterior` selects only components explicitly marked exterior
in the source table. A selector matching no component refuses. The committed
`cannone-h8-varnished.assembly` applies a 20 µm modern oil/copal surrogate to
the il Cannone top, back and ribs (45 components, 60 regions), excluding the
bass bar, soundpost and bridge. It does not establish the original Guarneri
finish or transfer Stradivari observations to that instrument. Compile the
varnish pack and use the violin recipe's `ASSEMBLY_OPTIONS` selection described
in [`violin/README.md`](../../../../data/instruments/descriptions/violin/README.md).

### Hereditary varnish and a structural frequency response

An explicit record **after `projection-work`**, before `regions` in V3 or
`components` in V4, selects the original passive material-memory owner:

```text
coating-memory prony-continuation 16 0.001 0.001 512
```

The fields are `prony-continuation TERMS STORAGE_RELATIVE_TOLERANCE
LOSS_RELATIVE_TOLERANCE MAX_MEMORY_COORDINATES`. Assembly options files accept
the same record. It applies to all declared films, with one complete assembly
history cap. The shell basis uses equilibrium film stiffness; projected
positive Maxwell branches supply frequency-dependent storage and loss through
the existing coupled renderer. Wood and interface losses keep their own
operators. Equal relaxation times are combined across regions before factoring
their complete modal matrices, preserving the cross terms without duplicating
one history bank for every small component. Body histories survive cold gas
and connection rebuilds and remain in the shared pressure/work/heat state.

`prony-continuation` explicitly accepts the passive fitted model outside the
measured fit band during a transient. The original source band is retained;
it is not expanded into a claim of measured DC or high-frequency behavior.
The two fit tolerances bound sampled numerical approximation errors separately,
not measurement uncertainty. Every retained body mode still needs source
support. A declaration above 256 histories requires `prepared-linear`, whose
actual combined mechanical/history dimension, memory and work limits remain
binding. The ordinary source-string and plate memory guards are unchanged.
New histories start relaxed at the current body displacement.

`cannone-h8-varnish-memory.assembly` selects this path for the same 20 µm O1
film and geometry as `cannone-h8-varnished.assembly`. Use it through the
existing `ASSEMBLY_OPTIONS` selector. Its 512-history declaration is an
independent cap, not a promise that every case at that cap fits the prepared
resources. The current O1 fit retains 13 positive scalar times and therefore
at most 390 histories for the complete 30-body-mode violin. With its 445
mechanical coordinates and at most 36 connection channels, the original
prepared-owner bound is 3,954,407,664 bytes and 145,345,993,394,700 multiply-add
terms, inside the existing 4 GiB / 1e15-work declaration. A full 512-history
case would exceed that byte budget and refuse. These are analytical resource
bounds for the declared model, not evidence of a completed whole-violin render.

The native renderer also accepts a read-only response export:

```text
--body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX
```

After staging a run with the launcher’s `--prepare-only`, for example:

```bash
target/release/material_strings \
  /path/to/new-memory-run/G3/performance.material-strings \
  /path/to/new-memory-run/G3/render.wav \
  --body-profile /path/to/new-memory-run/G3/body.profile \
  --body-response /path/to/new-memory-run/G3/body-response.tsv 101 2000 128 0
```

This writes logarithmically spaced complex bridge mobilities, signed exterior
and cavity volume velocities, wood/glue, interface and per-film dissipation,
supplied power, and power/solve residuals. Source identities travel in the TSV.
The transfer is normalized to a 1 N peak bridge force; this normalization does
not certify the source’s finite-amplitude strain envelope. The response uses
the exact source film law at each common drive frequency and refuses outside
the original wood/film/interface bands. It supports bare and hereditary shell
images; frozen-film images explicitly refuse this exact-source analysis.

The TSV describes the unloaded retained structural body and fixed reduced
wood/glue viscosity. The ordinary WAV calculation separately retains the
coupled strings, bow, cavity and chosen acoustic observation. For a controlled
comparison, stage bare and hereditary runs with equal mesh, mode count,
strings, loads, atmosphere and PCM scale. The response adds no sound-pressure
gain, amplitude normalization of the WAV or new time stepper. Reports admit
2–1024 frequencies within the complete work allowance and an 8 MiB output
limit; all frequency solves pass before a new report file is created.

For structural analysis alone, omit the WAV positional argument and select
`--body-response-only`. This resolves and projects the same body geometry and
material sources but constructs no string bank, cavity dynamics, bow performer
or prepared time propagator:

```bash
target/release/material_strings \
  /path/to/new-memory-run/G3/performance.material-strings \
  --body-response-only /path/to/new-memory-run/G3/body-response.tsv 101 2000 128 0 \
  --body-profile /path/to/new-memory-run/G3/body.profile
```

The performance still supplies the body dimensions, material state and mode
count; its declared input limits remain validated. Only the body/interface
source packs are loaded. Existing restricted-source options remain available.
The TSV explicitly excludes strings, bow/contact, cavity pressure feedback,
exterior reaction and configured microphone filters. Playback-only CLI
options such as contact heat, bow traces, stops and decimation refuse in this
mode. No WAV path is opened.

### Loaded bridge response and receiver pressure

`--loaded-body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX`, used with
the normal performance and new-WAV positional arguments, instead interrogates
the actual coupled performer before the WAV render. It uses the original
continuous harmonic owner with the retained strings, bridge connections,
selected cavity, complete material memory and admitted exterior radiation
reaction. Bow/contact traction is released for this linear transfer; it is
not a harmonic approximation to nonlinear bowed playing. Unsupported nonlinear
or adhesive attachments refuse through the original owner.

The report includes all bridge mobilities, structural exterior/cavity volume
velocities, the native and selected receiver pressures in Pa/N, supplied power,
component/connection loss, coating and other material loss, radiation loss,
mean stored energy and solve/power residuals. Coating heat comes from the
actual retained Maxwell branches. The unloaded structural report evaluates
the original film law; this loaded report evaluates the actual fitted
constitutive realization used for playback. Both retain the source band.

The original `prepared-linear` multiply-work declaration is the complete loaded
sweep allowance, divided across its frequencies. A small nonprepared run uses
`MAX_WORK`; no limit is silently increased. The 2–1024 frequency and 8 MiB report
limits and create-new semantics remain binding. Every solve passes before a
report file is created, and the query does not advance accepted state or
material histories. The normal WAV then continues through the performer.

An existing `shell_observer` declaration selects the original reciprocal
P1/P0 finite-distance observation owner, tied to the same complete closed-shell
BEM load, gas, basis and body/mouth traces. An existing `one_way_observer`
declaration selects the original ruled Burton–Miller producer on the profile's
acoustic skin, including every retained body mode and actual aperture-flow
share. Its report distinguishes loaded body modes from modes observed one-way;
partial or absent acoustic feedback is not certified negligible.

These receivers evaluate the original acoustic producer at the forcing
frequency. They preserve raw propagation phase without building a playback
FIR, adding its latency or taper, or charging another acoustic loss. The v2
TSV keeps `native_pressure` separate from `receiver_pressure` and adds the
three complex receiver pressure gradients when an exterior receiver is
selected. With no receiver declaration, both pressure columns contain the
native observer's pressure. `plate_observer` remains unsupported by this
shell-specific export. No exterior-wave convergence or finite-amplitude
strain certificate is implied. Structural and loaded response options remain
mutually exclusive on one invocation.

Set `BRIDGE_INDEX` to `all` to export the complete bridge transfer matrix.
Each frequency has one row per declared bridge, in physical bridge order;
`drive_bridge` identifies the column driven by a 1 N peak force. All bridge
velocities and both pressures retain complex phase. The original mechanics
owner factors once per frequency and solves each right-hand side against that
factor; the selected receiver also shares one BEM frequency sample across the
columns. The library's `material_shell_bridge_forces_harmonic_response` and
prepared receiver's `bridge_forces_harmonic_response` accept coherent complex
forces at all bridges and compute power and heat from their combined response.
Individual unit-column losses cannot be added to obtain coherent-drive heat.

For the same connected response without constructing a time propagator or
rendering a WAV, select `--loaded-body-response-only` and omit the WAV
positional argument:

```bash
target/release/material_strings \
  /path/to/new-memory-run/G3/performance.material-strings \
  --loaded-body-response-only /path/to/new-memory-run/G3/loaded-response.tsv 101 2000 128 0 \
  --body-profile /path/to/new-memory-run/G3/body.profile
```

This route uses the same source resolution, shell projection, string/body
connections, gas loading, material memory and selected acoustic loading as the
playback bank. The original coupled owner is explicitly constructed for
harmonic analysis: it retains the physical coefficients and the same
constant-loss damping columns, without allocating coupled time transition or
heat integral matrices or solving bow step compliance. Shell reduction and the
original physical components' initialization still run. Subsequent source/connection
rebuilds preserve that construction purpose. The frontend refuses to schedule
this bank, and attempts to advance its physical clock refuse.

An explicit `prepared-linear` record is required for this loaded analysis; its
mode, byte and multiply-work declarations are retained as resource limits,
including the complete sweep work allowance. This selects no new default
budget and does not prepare a time exponential. All original source and
physical admission checks remain active. A declared shell BEM or skin load
still binds its geometry-derived passive air reaction; a resource record alone
does not select the prepared playback backend. Unsupported contact/nonlinear
attachments and plate microphone selections still refuse; no connected
physics is silently removed to make the analysis fit. The original
`--loaded-body-response` option continues to prepare playback and render a WAV.

The complete sweep allowance includes every requested drive column and, when
selected, the receiver's boundary solve, quadrature and complex projection.
Receiver work is reserved before mechanics. Its declared frequency-sample
allowance and pass band also remain binding; use the intersection with the
body/material source bands. Matrix output remains subject to the same 8 MiB
report bound, so a many-bridge run can require fewer frequency points.

### Local coating strain and sufficient drive scale

Add `--body-response-strain MAX_BYTES MAX_WORK` to either structural or loaded
response command to observe the complete retained harmonic coating field.
The two positive integers explicitly admit the additional observer's storage
and its preparation plus every requested query. For `all`, queries equal the
number of frequencies times the number of bridge columns. These allowances
do not replace or increase the mechanics/receiver allowance.

```bash
target/release/material_strings \
  /path/to/new-memory-run/G3/performance.material-strings \
  --loaded-body-response-only /path/to/new-memory-run/G3/strain-response.tsv 101 2000 32 all \
  --body-profile /path/to/new-memory-run/G3/body.profile \
  --body-response-strain 536870912 1000000000000
```

The observer uses each original physical mode, exact director field, elastic
bubble condensation, material region and film placement. It prepares a fixed
48-by-18 Bernstein observation map once per coated triangle; each response
gathers coherent local nodal motion, applies this small map, and shares the
result across stacked films. Both film faces and all triangle interiors are
covered by the polynomial convex-hull bound. It includes plane-stress normal
strain and converts engineering shear to tensor shear before evaluating the
peak full-tensor Frobenius norm. No section shear-energy correction is applied
to the recovered strain.

The report records each film's source, component, signed faces, Poisson ratio
and authored small-strain envelope. It appends `12 * film_count + 3` columns:
whole-film Frobenius bounds, six tensor-component amplitude bounds, envelope
ratios, sufficient peak-force scales and ruling element/control/face indices,
plus the common ratio/scale/status. Component order is `[xx, yy, zz, xy, yz, zx]`
in facet axes; independent component bounds need not share a phase or location.
The controlling Bernstein coefficient is not generally a physical maximum
location. The 8 MiB output limit still applies; a many-film matrix may need
fewer frequency points.

For each 1 N row, `coating_strain_sufficient_peak_force_N` is a sufficient common
scale for these film bounds at this frequency. A ratio above one reports that
compliance is unresolved by the bound; it does not reject a normalized transfer
query or establish actual overstrain or damage. `inf` means no finite cap from
these film bounds because all observed film strains are zero or there are no
films. It is not an instrument strength claim. The authored envelope is a
constitutive modeling assumption, not a measured failure threshold.

Strain-enabled structural reports are v2 and loaded reports v3. Without this
option the existing v1/v2 schemas and preparation paths remain unchanged. All
queries finish before create-new output. Bounds apply to the retained linear
harmonic increment, with ordinary floating-point/condensation roundoff; they
exclude unmodeled preload/cure strains, transients, omitted modes, continuum
error and three-dimensional interlaminar recovery. They do not enforce local
strain during a subsequent nonlinear WAV render.

### Local coating strain during playback

Add `--coating-strain-guard MAX_BYTES MAX_WORK` to a WAV command to check the
initial retained shell displacement and every accepted mechanical endpoint.
It works with curved and assembled shell profiles, including hereditary film
stacks, contact-thermal playback, bow tracing and pressure decimation. The
positive allowances cover preparation and the entire scheduled horizon;
inadequate allowances refuse before the WAV is created. Response-only commands
use `--body-response-strain` instead.

The original cached MITC map observes real body coordinates without per-sample
allocation, source cloning, pressure normalization or a second dynamics solve.
The native log records per-film maximum Frobenius and component bounds, the
authored-limit ratio and mechanical sample, retaining source identity and the
ruling element/control/face. An unresolved envelope stops output and poisons
the performance. Its last accepted mechanics are not rolled back, and any
partial WAV remains explicitly incomplete under the existing writer policy.

The check includes static offsets present in the retained body coordinates.
It does not infer absent prestress or cure strain, bound motion between sample
endpoints, predict damage, or certify omitted modes/continuum behavior. Bounds
are deterministic floating-point spatial enclosures of the represented field,
not outward-rounded interval certificates. Bind the library monitor after
state-changing setup; a later preload cannot replace its checked initial state.

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
glue PACK MASTER_REGION SLAVE_REGION WIDTH_M
face-ties COUNT
face-tie MASTER_COMPONENT SLAVE_COMPONENT OFFSET ROTATIONS N_MASTER N_SLAVE
master-triangle TRIANGLE
slave-triangle TRIANGLE
glue PACK MASTER_REGION SLAVE_REGION WIDTH_M
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
`translations`. Without a `glue` record the tie uses the original rigid mortar
law. The optional record follows the complete node or triangle lists for that
one tie and selects the original compliant mortar. Nonmatching seams and
overlapping patches remain subject to that owner's geometric and work limits.
Prescribed slave coordinates that a rigid tie would eliminate refuse; compliant
ties retain both sides' supports.

### Material adhesive data

`glue PACK MASTER_REGION SLAVE_REGION WIDTH_M` reads an original normalized
interface pack relative to the performance directory. The named material
regions must own every edge segment or face triangle on their respective sides.
The card's ordered surface A/B material identities must match the actual
resolved master/slave regions, and the card must name its adhesive third body
and history. A bulk material swap therefore also needs a valid interface card
for the new pair. For an edge tie `WIDTH_M` supplies its actual bonded area
`width × seam length`; face ties use their actual overlapping facet area. The
shared glue law requires a positive width for either record.

All six area-specific source properties are required:

| Property | SI unit |
| --- | --- |
| `glue-normal-stiffness` | N/m³ |
| `glue-shear-stiffness-1` | N/m³ |
| `glue-shear-stiffness-2` | N/m³ |
| `glue-normal-viscosity` | N·s/m³ |
| `glue-shear-viscosity-1` | N·s/m³ |
| `glue-shear-viscosity-2` | N·s/m³ |

Stiffnesses must be positive. Each viscosity must be nonnegative; zero loss
must be explicitly sourced. Coefficients resolve at the performance's solid
`state T OMEGA`, and every selected claim must state its frequency applicability.
The complete retained eigenvalue intervals must fit the intersection of those
bands before the performer is published. Fixed-temperature source curves are
admitted; frequency-dependent coefficients and constant-loss anchors require a
different constitutive realization and refuse. Changing the declared solid
state requires reloading the original interface sources.

The profile retains original packs, ordered surface descriptions, adhesive and
history, claim receipts and uncertainties through `ResolvedGlueLine`. The joint
stiffness participates in the common body eigensolve, and its complete projected
viscosity tensor joins regional material loss on that same basis. No modal loss
or timbre coefficient is inferred from an adhesive name. These are fixed-state,
linear Kelvin–Voigt bonds; creep memory, debonding, slip and evolving cure or
temperature are outside this profile.

Each final `bridge` is one physical nodal force footprint, in actual string
order, on any component. Rotational bridges can use footprints on different
components; they retain the same physical chord and perpendicular-lever rule
as v3. The other string controls, material bow interfaces, afterlength, torsion,
preloads and shared gas continue through the original performer.

Profiles take the 4 MiB shell limit and admit at most 64 components, 12288 total
nodes, 24576 total triangles, 157 ties and 256 shared body modes. Per-component v3
mesh limits still apply. Material/interface source bytes and every damping
channel remain under the existing complete-performance budgets. At least one
region across the assembly must retain `@body`; other regions can select
different source packs. Repeated bulk paths share the original material cache.

This is unloaded linear component mechanics with rigid or source-bound compliant
ties. It does not add prestressed assembly equilibrium, adhesive failure, solid
soundpost contact, a certified closed cavity skin, or resolved interior and
exterior acoustic fields. The inherited compact observer and Helmholtz cavity
remain explicit approximations. The new native regression sources check tied
motion, complete material-to-pressure/PCM propagation, material-only changes
and callback equivalence; they still require execution in the pinned build.

A `frankensim.shell-assembly.v1` table (the Esteso guitar body,
`data/instruments/geometry/guitar-esteso-1964-zimmer2024/shell-assembly-h12.tsv`)
is written as this v4 profile by `body/shell/assembly_table.rs`, choice for
choice as `examples/guitar_body_modes.rs` assembles it: orthotropic regions on
their own axes at the part's thickness, a free body, the part flags' exterior
and cavity facets, full-rotation rigid edge ties on the master's skin or
midsurface. What the table does not carry is explicit: the window, projection
work, which wood is `@body`, the reference thickness, and one bridge per string
at the nearest node of a declared point (the snap distance is reported), each
pushing along one declared direction. The table's cavity row supplies the
performance's cavity record. The il Cannone violin body
(`data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h8.tsv`,
48 components and 119 ties; its 6 mm mesh `shell-assembly-h6.tsv`, 64 components
and 157 ties, sets the two caps above) converts the same way;
its frame has x across the body, y vertical and z along the strings, so its
strings push their bridge nodes along x. That x is also the bow's tangent, so
the il Cannone rungs give every string `polarizations two 1 0 1 0`: the bow
drives the plane the bridge holds. The string record carries one bridge
direction, so the bow-normal plane has no bridge, and nothing drives it. With
the violin-sized rung's `0.6 0.8 1 0`, 0.8 of the bow drove that free plane, and
the G string sounded its full 0.383 m span (167 Hz) instead of its 0.328 m
speaking length.

With the two violin shell packs, both built by `xtask matdb-pack` from
`data/matdb/seed-v1` (`packs/violin-spruce-norway-viala2020-shell.fsmatpk` as
`@body` for the top, bass bar and soundpost,
`packs/violin-maple-sycamore-sonderegger2013-shell.fsmatpk` for the back, ribs
and bridge), the window from 15 Hz to 2 kHz holds 30
assembled modes, 473.5 Hz to 1967.7 Hz. A request for more refuses once the
window is solved and names the count it holds. Both packs state
`constant_loss_reference_omega` (2 pi 200 rad/s), so the body modes' Q is
68.3-80.3 across the window; without it the same packs gave a Kelvin-Voigt Q
falling from 31.0 at 473.5 Hz to 7.2 at 1915.7 Hz. The ignored rungs in
`body/shell/tests.rs` bow the violin-sized G3 rung's synthetic steel strings
(not violin strings) on this body: a quarter-second cost rung printing every
body mode with its Q, the construction and render times and the pressure bits,
and a 2 s render for listening. On trj the body constructs in 315-350 s and the
coupled network (245 coordinates, 36 connections) renders at 100-125 s per
audio second.

### A table from the command line

`material_strings assemble-shell TABLE OPTIONS NEW_PROFILE` writes a table as
this v4 profile and admits it before it writes it. OPTIONS states the choices
the table does not carry:

```text
frankensim-shell-assembly-options-v1
window-hz LO HI
projection-work N
reference-thickness-m T
pack TABLE_PACK PROFILE_PACK          (one per table pack; @body is the performance body pack)
bridge-component NAME
bridge-direction X Y Z
bridge X Y                            (one per string, in string order, table coordinates [m])
```

Each record except `pack` and `bridge` appears exactly once, and a line that
starts with `#` is a comment. The command prints:
- the body length and width [m] and the reference thickness, for the
  performance's body record;
- each bridge's snap distance;
- the table's cavity row (volume, neck radius, neck length), for the cavity
  record.

`cannone-h8.assembly` holds il Cannone's choices. The rungs in
`body/shell/tests.rs` read the same file. `cannone-h6.assembly` makes the same
choices for the 6 mm table `shell-assembly-h6.tsv` beside it, with the
projection work its 11769 nodes need, except the window: it reaches 5 kHz, the
narrower pack's validity. The spruce pack declares its source's band (Viala's
model updating, 100-5000 Hz) and the maple pack Laemmlein's 150 Hz to 10 kHz
sweep; no source reaches 5-6 kHz. The cavity row and the bridge snaps differ
slightly from h8's, so `assemble-shell` prints them for the performance.

`cannone-dominant-g3.material-strings` bows il Cannone's G for 2 s at 144 kHz,
on the four Dominant string packs (README.md). It uses the rung's 0.383 m span
and 0.328 m speaking length; the contacts table puts il Cannone's own G at
323.1 mm bridge to nut. The G keeps its bending and twist modes to 30 kHz; D, A
and E keep theirs to 10 kHz, which is 474 coordinates with the body and cavity.
The bow presses 1 N and moves at 0.2 m/s, a tenth of the speaking length from
the bridge. The WAV's full scale is 1 Pa. The peak at the microphone is
0.34 Pa on this body and 0.37 Pa on h6 at 5 kHz, both in the first 50 ms.

Two inputs of this note are authored, not measured:
- the bridge joint (200000 N/m, 0.1 N s/m);
- the Stribeck friction law (0.8, 0.4, 0.1 m/s).

`cannone-dominant-g3-crown-compliant.material-strings` and
`cannone-dominant-g3-crown-stiff.material-strings` replace the authored joint.
They are the same note, but every string's joint comes from its interface pack
`packs/violin-cannone-crown-dominant-{e,a,d,g}-interface.fsintpk`.

`gen_crown_interface_packs.py` writes both notes and the four sources from
committed data only:
- the maple shell pack;
- the wrap-metal packs;
- the CT break angles and bridge-top thickness;
- this note's string radii and tensions.

Each joint is the string's lateral stick on the maple crown: Hertz crossed
cylinders and Mindlin's no-slip tangential stiffness, k_t = 8 G* a. The joint
loss is the maple's shear retardation time times its share of the compliance.
That dashpot is written at the maple pack's constant-loss anchor (2 pi 200
rad/s), and each source states the anchor as `constant_loss_reference_omega`.
The prepared network then gives every coupled mode the joint's loss factor
(share x eta_LR, about 0.016) times the joint's share of that mode's strain
energy, at the mode's own frequency, instead of the Kelvin-Voigt reading that
grows with frequency. The note needs the prepared owner, which it declares.

The stiff end's coupled stiffness bound exceeds the network's 0.9 pi f_s guard
at 144 kHz, so its note renders at 384 kHz over the same 2 s, with 2.67 times
as many steps.

Every claim covers 0 to 2 pi 200 kHz: the contact is quasi-static there, and the
network asks the law at the cavity neck's free coordinate and at its coupled
stiffness bound. The loss outside the maple's measured 150 Hz-10 kHz is a STATED
UNKNOWN; there it is the five-branch constant-loss shape, extrapolated.

Three inputs are STATED UNKNOWN. Each is taken at both ends of its range:
- the maple's orthotropy;
- the crown edge radius, from 0.1 mm to half the CT's 1.8 mm bridge top;
- the wrap metal Pickering does not name: aluminium to silver, and aluminium to
  steel for the E.

| string | compliant k_t, D | stiff k_t, D |
|---|---|---|
| G | 1.08e5 N/m, 1.37 N s/m | 1.30e6 N/m, 15.9 N s/m |
| E | 1.05e5 N/m, 1.33 N s/m | 1.30e6 N/m, 16.4 N s/m |

The D and A strings lie in between.

The elastic contact's mean pressure (170 MPa to 2.2 GPa) is beyond what end-grain
maple bears, so a played crown beds into a conforming groove. That only enlarges
the contact: the compliant end is a lower bound, and the stiff end can be
exceeded. The winding's turn-by-turn contact is not resolved.

Gather the packs and the performance in one directory, then run from the root:

```bash
demo=$(mktemp -d)
m=crates/fs-couple/examples/material_strings
cp $m/packs/violin-*.fsmatpk $m/packs/violin-cannone-crown-*.fsintpk $m/cannone-dominant-g3*.material-strings "$demo/"
cargo run --release -p fs-couple --example material_strings -- assemble-shell \
  data/instruments/geometry/violin-cannone-guarneri-1743-esrf2026/shell-assembly-h8.tsv \
  $m/cannone-h8.assembly "$demo/cannone-h8.body"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/cannone-dominant-g3.material-strings" "$demo/g3.wav" --body-profile "$demo/cannone-h8.body"
```

For a crown end, run the last command on
`cannone-dominant-g3-crown-compliant.material-strings` or
`cannone-dominant-g3-crown-stiff.material-strings`, with a new output.

## Authored distributed interior air on the material shell

An optional `interior-fe-v1` trailer follows every ordinary record of a v3
shell or v4 assembled body profile. It gives the original Helmholtz finite
element owner the actual volume mesh. No material-specific cavity preset,
assigned resonance, frequency tuning or independent pressure filter is used.

```text
interior-fe-v1 ORDER MAX_HZ MAX_PRESSURE_MODES MAX_DOFS MAX_BOUNDARY_FACES MAX_PROJECTION_TERMS RITZ_TOL
air-nodes N
air-node X_FRACTION Y_FRACTION Z_METRES
... exactly N rows ...
air-tets T
air-tet A B C D
... exactly T rows ...
```

Every node index is zero based. X and Y use exactly the shell's existing
fractional convention: multiply by the performance body's length and width.
Z is already in metres. Reuse the exact original boundary coordinate values;
the adapter neither welds nearby points nor closes gaps. Every tetrahedron
needs four distinct valid vertices and nonzero physical volume, every volume
vertex must be used, and repeated tetrahedra are refused. The volume must be
connected and conforming, and every boundary triangle must match one explicitly
selected original `cavity-facet` or a declared compact opening below. Solid
pressure facets are oriented outward from the gas. A profile
without an explicit `pressure-surfaces` block cannot bind a distributed
interior. Exterior pressure roles remain independent.

The order is 1–4. The positive frequency ceiling must stay below 0.45 times
the mechanical sample rate. The original inertia-counted eigensolver retains
the **entire** Neumann eigenwindow, including uniform pressure; exceeding
`MAX_PRESSURE_MODES` refuses instead of dropping modes. `RITZ_TOL` is the
relative eigen residual tolerance and must lie strictly between zero and one.
The source admits at most 256 pressure modes, 8192 pressure degrees of freedom,
8192 boundary faces, 32768 tetrahedra and one billion projection terms. It also
retains the complete profile's 4 MiB byte ceiling. Higher order contributes edge,
face and interior degrees of freedom; `MAX_DOFS` bounds the complete FE
space, not only its vertices.

The network reserves up to `min(MAX_PRESSURE_MODES, actual FE DOFs)-1`
additional standing-wave coordinates and one compression connection per
pressure mode before material-pack reads. Each declared opening also reserves
one fluid-inertia coordinate and one possible internal-damping connection,
including when its declared resistance is zero. The exact uniform pressure mode contributes
compression storage without an invented moving air mass. These reservations
share the existing complete mode, connection, setup-work and optional
`prepared-linear` budgets with string polarizations, torsion, interfaces and
cylinder-fluid memory. The original exact owner still applies its connected
block limit. Structural modes, the pressure frequency window, element order
and mesh resolution remain explicit convergence axes.

The performance must select `cavity none`: a simultaneous lumped neck cavity
would duplicate the same gas storage. Its existing `gas` record supplies the
same density and sound speed to strings and the original pressure FE, including
source-backed thermochemical records when supplied. A change to pressure,
composition or temperature rebuilds air storage and standing-wave frequencies.
A replacement solid material rebuilds the original shell basis and integrates
the same pressure geometry against its new displacement fields. The existing
work-conjugate coupling supplies both air reaction on the shell and shell work
on the air inside the complete friction solve.

### Run the closed-shell interior fixture

`closed-shell-cavity.body` extends the existing synthetic 42-node,
80-triangle shell with all 80 explicit pressure facets and a center-star
tetrahedral volume. It retains the identical boundary coordinates and material
pack. Its order-two FE has 205 pressure degrees of freedom; the authored
1100 Hz window has a 24-mode allowance. One constrained-shell structural
coordinate and the original four-mode string keep the integration fixture
small. The fixture is a manufactured source-to-performance example, not a
measured violin, cello or bass.

From the repository root:

```bash
demo="$(mktemp -d)"
cp crates/fs-couple/examples/material_strings/*.tsv "$demo/"
cp crates/fs-couple/examples/material_strings/closed-shell-cavity.material-strings "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wire.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/closed-shell-manifest.tsv" --out "$demo/closed-shell.fsmatpk"
cargo run --release -p fs-couple --bin material_strings -- \
  "$demo/closed-shell-cavity.material-strings" "$demo/closed-shell-cavity.wav" \
  --body-profile crates/fs-couple/examples/material_strings/closed-shell-cavity.body
```

The sample retains the existing compact body-pressure observer, so its WAV is
an integration diagnostic rather than a converged full-band external sound
field. It does not reuse the older example's 120 Hz BEM band for a distributed
air window extending to 1100 Hz. An independently admitted closed-shell BEM
load and finite-distance shell observer can still be selected through their
existing performance records; their frequency bands, spatial resolution and
error budgets need to cover the desired observation.

Without an explicit thermal declaration, this interior path represents frozen,
homogeneous gas with lossless linear acoustics on straight-sided tetrahedra and
the original linear shell kinematics. The thermal extension below adds gas
heat diffusion. Viscous wall shear and nonconforming volume/shell interfaces
remain separate physical operators. Missing boundary assignments are refused;
no fictitious solid closure surface is added.

### Thermal diffusion in the same cavity and mechanical solve

After the last `air-tet`, or after the complete optional `air-openings` block,
select the original gas heat equation with:

```text
air-thermal isothermal-walls MAX_DECAY_PER_S MAX_THERMAL_MODES MAX_THERMAL_DOFS MAX_PROJECTION_TERMS MAX_UNRESOLVED_FRACTION
```

The solid wall temperature perturbation is fixed to zero. Actual opening faces
have zero diffusive heat flux, so the declaration does not invent a heat
reservoir across an aperture. The retained gas supplies conductivity, density,
heat capacity and gamma; these determine the diffusivity `kappa/(rho*cv)`.
Heat modes come from the original volume FE stiffness and consistent mass
matrices, with the actual wall trace constrained. Every positive heat mode in
the declared decay-rate window is retained, subject to the explicit mode and
work allowances.

Those modes become passive compression-relaxation states in the existing
string/body/cavity propagator. They change both pressure reaction and decay,
including the uniform compression mode. The instantaneous compression remains
the original adiabatic law; the captured thermal field relaxes toward its
isothermal equilibrium. There is no independent heat timestepper or assigned
modal damping. Original pressure and neck observers include accepted thermal
history. A bare coordinate-only pressure readout refuses when it would omit
that history.

The part of the compression field outside the retained thermal space stays
adiabatic. `MAX_UNRESOLVED_FRACTION`, strictly between zero and one, bounds the
infinity norm of the missing normalized compression projector. It measures
the FE/window representation, not overall instrument accuracy. Refining the
volume and retaining more heat modes are the way to reduce it; the source
refuses an unresolved declaration. Each retained heat mode reserves a possible
first-order state in the complete material-memory budget before pack I/O.
The `prepared-linear` backend is refused because it has no such states.

`closed-shell-thermal.body` reuses the original manufactured closed-shell
geometry with one pressure mode and this declaration:

```text
air-thermal isothermal-walls 1.0 43 256 50000000 0.201
```

Use it with `closed-shell-cavity.material-strings` and the same run command
above, selecting the thermal body profile and a new output path. An independent
exact polynomial FE calculation for the authored 0.4 by 0.4 m body gives 205
volume DOFs, 43 free heat DOFs and a uniform-compression projector remainder of
0.20. That explicit coarse remainder is why the fixture declares 0.201. For
293.15 K dry air the retained decay rates are about 0.0084–0.1804 per second.
This source demonstrates physical coupling and retained history; its coarse
mesh does not resolve audio-frequency thermal boundary layers.

The model covers linear gas thermal conduction to fixed-temperature walls.
Viscous boundary shear, finite wall heat capacity and effusivity, mean flow and
aperture enthalpy transport are not supplied by this selector. Material-card
and gas substitutions rebuild their original physical owners; this fixture
does not establish measured violin, cello or bass fidelity.

### Explicit compact openings through the original volume and exterior

After the last `air-tet`, an optional opening section names original volume
faces and internal neck data:

```text
air-openings COUNT
air-opening INTERNAL_LENGTH_M INTERNAL_RESISTANCE_PA_S_M3 FACE_COUNT
air-opening-face A B C
... exactly FACE_COUNT rows, then the next opening ...
```

There may be one to eight openings. Indices are the original `air-node` IDs;
they are never inferred from the shell's triangle row numbers. Faces must be
distinct across all openings. The shared physical owner checks that each face
is on the volume boundary, has no remaining solid shell triangle over it, and
belongs to one connected planar patch with a single boundary loop. Every rim
translation must be fixed in every retained shell mode. Remove the actual
solid triangle and its exterior/cavity pressure-role references when declaring
a hole; retaining a solid triangle with zero pressure weight does not make it
an opening.

`INTERNAL_LENGTH_M` is positive and contains only the internal duct length and
any explicitly included interior correction. The nonnegative resistance is
only internal loss. The original shell BEM supplies exterior reaction once,
including its cross-coupling between the structural surface and opening flow.
No exterior end correction, radiation resistance, assigned resonance or
instrument-specific sound calibration is added by the file adapter.

Opening-bearing profiles require `body-radiation bem ...` before material-pack
reads. That band's upper endpoint must cover the selected interior pressure
window, and the original joint owner additionally checks the reference neck
frequencies and its compactness bound. The cavity binds before cylinder-fluid,
torsional, nonlinear axial and material-memory attachments. The complete
performance is returned only after the original exterior owner has bound the
joint structural and opening load. Existing sealed files keep their previous
grammar and construction path.

The existing `shell_observer` record selects the actual external pressure
output. Its shared BEM transfer takes original body modal velocities and each
accepted physical opening volume flow. It preserves their interference and
the observer's explicit numerical latency; no separately driven opening tone
or output normalization is introduced.

`open-shell-cavity.body` is a manufactured variant of the closed-shell fixture.
It removes solid triangle `[0,12,14]`, whose three vertices are fully constrained,
and retains that exact face in the original 80-tetrahedron volume. The explicit
opening has a 0.04 m internal length and zero internal resistance. Its order-two
interior selects the complete 0–180 Hz window with one pressure-mode allowance;
the same eigensolver refuses if that allowance cannot retain the full window.
The paired performance declares BEM loading over 0.4–180 Hz and an 8 kHz
microphone with a 60 Hz pass-band endpoint and a 180 Hz taper stop. The nominal
uniform-neck resonance lies within that declared taper, so this example is a
finite-band source/observer integration, not full-passband instrument audio or
a measured violin, cello or bass. Frequency, panel and FIR error gates remain
those of the original owners.

The joint load now uses the existing symmetric P1/P0 Galerkin source; the
receiver evaluates that same source's P1 pressure and P0 auxiliary density.
Independent body and mouth columns retain the unchanged `1e-12` reciprocity,
passivity and causal-fit gates. No collocation matrix is repaired. The owner's
slab resonance exclusion and angular-resolution chart limit the declared band;
spatial, angular, receiver and FIR refinement remain independent obligations.

The 80-panel, 42-vertex, two-source load requires 162,083,000 original work units
per frequency and 10,535,395,000 over its 65 frequencies. Its file therefore
declares a bounded 12-billion-unit source allowance. The microphone has an
80-billion-unit boundary-sweep allowance and separately declares at most two
million adaptive receiver-kernel evaluations per frequency. These are numerical
work ceilings, not time estimates or relaxed error tolerances.
`shell_angular_refinement 0` and `shell_quadrature 1e-7 14 2000000` expose the
existing owners' refinement controls as input data.

The integration inputs and command remain pending native execution. Their test
retains all construction, matrix, fit, field and sample-level assertions as hard
gates; this implementation is not measured-instrument qualification.

After generating the two material packs using the commands above, the pending
integration command is:

```bash
cp crates/fs-couple/examples/material_strings/open-shell-cavity.material-strings "$demo/"
cargo run --release -p fs-couple --bin material_strings -- \
  "$demo/open-shell-cavity.material-strings" "$demo/open-shell-cavity.wav" \
  --body-profile crates/fs-couple/examples/material_strings/open-shell-cavity.body
```


