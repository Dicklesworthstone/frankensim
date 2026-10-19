# Directional body materials in the shared string performer

`material_strings PERFORMANCE NEW_OUTPUT.wav --body-profile PROFILE` selects
an explicit rectangular material chart instead of the legacy isotropic body
image. Omitting the flag retains the original construction and material profile.
The performance's `body` record still supplies its pack, dimensions, thickness,
mode count and observer distance. String layers, support prescriptions, torsion,
loads, plucks, bow tracks, shared gas and optional cavity are unchanged.

This is the existing `from_viscous_material_plate_chart` path, not a new solver.
The material tensor is projected on the SAME eigenbasis as bridge displacement,
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
`shear_modulus_12`. It also supplies ALL six reduced in-plane viscosities in Pa s:
`plane_stress_viscosity_11`, `_22`, `_66`, `_12`, `_16`, `_26`. The engineering
shear convention is the original shared material owner's convention. Signed
couplings and explicit zero entries are preserved; the tensor must be positive
definite or explicitly all zero. Missing entries never mean zero. The legacy
`kelvin_voigt_bending_viscosity` cannot substitute for this complete tensor.
All original state, units, source-band, coefficient and matrix checks apply.

Profile input is bounded to 64 KiB; the original source-pack and complete-network
limits remain. Viscous projection is capped at 16,777,216 operations and 16 body
modes. The dashpot factor uses a 1 N s/m coordinate scale; this changes neither
the projected damping matrix nor material coefficients. The exact profile bytes,
original packs, regional material states, tensors and projected damping remain
attached to the performer. No live source-state mutation is implied.

The current model is a fixed-state, flat, linear rectangular DKT plate, not a
carved shell, soundpost, bass bar, glue interface or empirical instrument model. Pressure remains the inherited compact baffled
body/neck observation, not full exterior radiation or radiation back-loading.

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
The profile must be v2 with distinct physical nodes for the chosen pair.

If the retained nodal displacements are `w_positive` and `w_negative`, body
rotation is `(w_positive - w_negative) / separation_m`. The separation is the
Euclidean distance of those SAME nodes in the performance's declared body
length and width; it is recalculated when geometry changes. Swapping the ordered
pair reverses this scalar coordinate. Its sign and direction declare the
rocking axis assumed by the string joint; no 3D string/body orientation is
inferred or certified. This chord rotation is exact for rigid rocking and a
secant approximation for a deforming plate, not a resolved local slope claim.

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
