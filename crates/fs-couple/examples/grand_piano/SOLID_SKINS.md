# Real soundboard thickness and crown from paired OBJ skins

`piano_solid_import` turns explicitly selected upper and lower soundboard skins
into the **existing stiffened CST/DKT shell**, with crown and element thickness
computed from the supplied geometry. This closes the previous requirement to
hand-author a midsurface and populate every thickness from a nominal value.
There is no new vibration solver, sample bank, EQ, or mesh in the audio loop.

```sh
cargo run --release -p fs-couple --example piano_board_import -- inspect piano.obj
cargo run --release -p fs-couple --example piano_solid_import -- \
  piano.obj materials.fspi skins.fsps soundboard.fss
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry soundboard.fss \
  --midi performance.mid --midi-half-pedal --render piano.wav
```

All paths above name supplied inputs and fresh outputs, not bundled assets.
The existing Model D scale, felt/shank cards, pedals, strings and receiver
controls compose with the imported board. All required key bridge stations
must still be supplied. Downbearing is a separate operation on an explicitly
unloaded reference; do not settle an already loaded scan a second time.

## Geometry and material input

Use the existing `frankensim-obj-board-v1` FSPI (see `MESH_IMPORT.md`). Its `part`
selects the **upper** skin. Materials, grain, unit scale, orthonormal frame,
bridge stations, support selection and bonded rib/bridge sections remain
explicit. Support and beam vertex references identify upper OBJ vertices.
The shell supports clamped boundaries and zero authored membrane pretension.

A separate FSPS declares the lower skin and correspondence, independently of
OBJ vertex order, face winding, shader names or unrelated cabinet objects:

```text
frankensim-board-skins-v1
lower,soundboard_underside
thickness,geometry
thickness-range,0.002,0.020
pairing,projected,1e-7
```

For registered skins, `pairing,projected,TOLERANCE` finds each upper vertex's
unique lower counterpart in the declared frame's XY plane. The tolerance is in
**metres after the source-unit transform**, bounded to 1e-12..1e-4 m. This mode
handles arbitrary OBJ vertex ordering, not differently tessellated surfaces.
Matching uses a bounded spatial grid, not all-to-all search. Every vertex must
have exactly one candidate: ambiguous matches, missing matches, reused lower
vertices, unresolvable grids and excessive candidate work refuse. No nearest
candidate is silently preferred and no nodes are welded.

For skins whose paired vertices are offset in XY, omit the `pairing` row and
supply `pair,UPPER_INDEX,LOWER_INDEX` for every vertex instead. Indices are
one-based source OBJ vertices. Automatic and explicit modes cannot be mixed.
For example, `pair,1,101` maps upper vertex 1 to lower vertex 101.

The upper and lower sets must be
disjoint and their paired triangle connectivity identical. Side walls are
excluded explicitly with object/group labels. Missing faces, duplicate faces,
missing/duplicate pairs, overlapping selections and incompatible topology
refuse, rather than guessing a repair. The correspondence is a supplied
geometric hypothesis, not proof that a render asset has manufacturing accuracy.

`thickness,geometry` explicitly overrides the **nominal thickness** in each
FSPI material row. Density, E_L, E_R, nu_LR, G_LR and grain are retained; no
physical property comes from an MTL shader. Each pair's midpoint becomes a
shell node. For a facet with unit normal n, its section thickness is
`mean(dot(upper_i - lower_i, n))` in metres. Consequently curved or tilted
skins use normal separation, not the longer point-to-point gap. Element
thickness enters both the existing mass and stiffness assembly.

This is a single homogeneous, piecewise-constant section approximation; it
is **not** exact integration of arbitrary taper, a laminate layup or a general
solid-to-shell mesher. Refine supplied triangles to resolve rapid thickness
changes. The output source description reports the minimum/maximum vertex
normal thickness and maximum within-element thickness spread. Beam dimensions
and eccentricities are not inferred or changed: provide them relative to the
new midsurface. Supplied glue/bond behavior is the existing perfect-bond model.

Midpoint heights must be within 50 mm of the declared plane. Midsurface normals
must have chart z component >= 0.95. Skin facets and pairing columns must align
with the corresponding midsurface normal to cosine >= 0.95; every vertex's
normal separation must lie within the explicit SI thickness range. The range
itself must be within 10 micrometres..100 mm. Steep, folded, inverted or highly
sheared correspondences refuse. Native plate and shell admission runs before
creating the output. General self-intersection and specimen fidelity remain
unverified. Inputs are bounded and outputs never overwrite existing files.

## Steinway asset and material research (checked 2026-09-24)

Steinway's Model D specifications identify a 2.74 x 1.56 m envelope, Sitka
spruce soundboard tapering from 9 to 6 mm, sugar-pine ribs, maple-capped
hardwood bridges, steel/copper strings and wool felt. These support species
and dimensional choices, not specimen elastic tensors, damping or calibration:
https://www.steinway.com/pianos/steinway/grand/model-d

The seavenois **Steinway D274** model is page-declared CC0, Blender 2.6x:
https://blendswap.com/blend/7279 . Its download page requires sign-in:
https://blendswap.com/blend/7279/download . No mesh bytes were obtained or
redistributed for this implementation, and usable skin topology is unverified.

The CC0 **BWV846Prelude_piano_animation** derives from that model and supplies
procedural wood/metal appearance. Its author explicitly warns that assigned
material names can differ from actual piano-part materials; a cast-iron plate
is labelled AnodizedMetal. Never turn these visual names into constitutive laws:
https://blendswap.com/blend/28847 . This source was inspected, not downloaded.

## Focused native regressions

```sh
cargo test --release -p fs-couple --example piano_solid_import mesh_import::crowned::solid
cargo test --release -p fs-couple --example piano_solid_import cli_tests
```

Tests exercise supplied crown, geometry-controlled mass, independence from
nominal FSPI thickness, thickness-driven shell modes, automatic pairing under
vertex permutation and rotated millimetre frames, and invalid-input refusal.
End-to-end regressions run imported geometry through the existing source Model
D strings, felt and shanks to receiver pressure and PCM WAV, checking physical
response changes, energy accounting and exact block-splitting determinism.
Another regression preserves all 88 Model D bridge stations, supports and
bonded rib paths. The existing piano workflow's `piano_board_import` test target
also includes these shared tests. They are authored native tests; the editing environment
has no Rust compiler, so this change does not claim a native test pass or an
audibly validated Steinway digital twin.
