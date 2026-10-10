# Real geometry and physical material assignments

`piano_board_import` connects editable OBJ geometry to the existing
`grand_piano` soundboard solver. It does not introduce a plate solver, an
oscillator bank, a hammer law, or a mesh in the audio callback.

## Use the current Model D geometry as an editable asset

These commands generate the source-derived Model D reconstruction, export its
physical sections and topology, and bring the edited mesh back to the solver.
Use fresh output paths (the importer refuses overwrites).

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dump-geometry model-d.fsb
cargo run --release -p fs-couple --example piano_board_import -- \
  export model-d.fsb model-d.obj model-d.fspi
cargo run --release -p fs-couple --example piano_board_import -- \
  inspect model-d.obj
cargo run --release -p fs-couple --example piano_board_import -- \
  import model-d.obj model-d.fspi imported.fsb
cargo run --release -p fs-couple --example piano_board_import -- \
  render-steinway imported.fsb piano.wav 6
```

`render-steinway` combines the imported board with all 88 source-derived Model D
courses, per-key wool/Prony hammer cards and published shank geometry. Silent
keys remain in the coupled resonator. It preserves the source string tensions
rather than substituting MIDI frequencies or applying automatic retuning. The
default demonstration strikes key 69 at 2 m/s; append a Standard MIDI File path
after the duration to use the existing MIDI gesture importer (channel 1,
velocity 127 -> 4.5 m/s, switch sustain; an explicit uncalibrated mapping).

This offline path retains 24 partials per string, board modes through 400 Hz,
and four mechanics substeps per 48 kHz frame. Its receiver is at
[0.675, 1.0, 1.0] metres in the imported board chart. It uses the existing
infinite-baffle physical-pressure observer and PCM encoder at 2 Pa full scale,
without normalization; clipping and the mechanics energy balance are reported.
These defaults are not a full-audible-band or measured-SPL claim. Every one of
the 88 bridge stations must be present; missing stations refuse.

For other string/hammer cards, sample rates, modal budgets and microphone
positions, use `grand_piano --board-geometry imported.fsb --scale strings.csv
--hammers felt.fsh --render piano.wav`. The existing `grand_piano --preset
steinway-d --dump-scale strings.csv` command exports the source-derived scale.
Selecting a board in that general CLI alone does not enable preset voicing.

The exported OBJ is the **physical midsurface**, not a full cabinet rendering.
Every distinct native element section gets a `section_N` material assignment;
its thickness, density, longitudinal/transverse moduli, Poisson ratio, shear
modulus and grain angle travel in the sidecar. Explicit rib/bridge beam paths,
supports and all key bridge coordinates travel with it. No physical properties
are read from an optical MTL file. When editing or re-exporting an OBJ, preserve
vertex order or update every `fixed` and `stiffener` vertex reference. Moving a
bridge requires editing its physical location, not just moving a visual object.

## Import a separately sourced mesh

Run `inspect` first. Pick an exact object/group name for one conforming panel
midsurface, not the entire piano. Supply this UTF-8, comma-delimited sidecar:

```text
frankensim-obj-board-v1
source,estimated,Illustrative parameters only; replace with your documented cards
part,soundboard
units,1
frame,0,0,0,1,0,0,0,1,0
flatness,1e-8
material,spruce,0.008,450,1e10,8e8,0.3,6e8,0
support,clamped
boundary,all
damping,0.01
pretension,0
bridge,69,0.5,0.5
```

This small example is **not a Steinway calibration**. The bridge must lie in the
selected panel, and a complete scale requires a bridge row for every key.

`units` is metres per source coordinate unit (millimetres use `0.001`). `frame`
is origin x,y,z in source units, then unit U and V axes in source coordinates.
They must be orthonormal; U cross V defines the normal. Chart x,y and bridge
positions are in metres; grain angle is measured from chart U, not world X.
`flatness` is the admitted normal projection distance in metres (maximum 1 mm).
The actual maximum is reported. Crown beyond that band refuses: a curved shell
is not silently replaced with a flat plate.

`material,name,h,rho,E_L,E_R,nu_LR,G_LR,grain_rad` maps an exact OBJ `usemtl`
name to a physical orthotropic section. All dimensions are SI. `material,*,...`
is an **explicit** fallback section; without it any unmapped label refuses.
Different regions can carry different wood, taper and grain. MTL libraries are
retained as unresolved names only and are never opened or downloaded.

Use `fixed,OBJ_vertex` for individual supports, or `boundary,all` to explicitly
support every topological boundary, including holes. `support` selects clamped
or simply_supported. Positive OBJ vertex indices are one-based. They are
remapped after removing unselected geometry, never welded by guesswork.

`stiffener,E,G,A,I,J,eccentricity,rho,v0,v1,...` uses the existing FSB beam law
and source OBJ vertex indices. The path must be present in the selected panel;
a labelled rib render mesh is not automatically a measured beam section.

`bridge,key,x_m,y_m` locates a station barycentrically in the actual triangles.
Outside stations refuse instead of snapping to a nearby vertex. Degenerate or
duplicate faces, nonmanifold edges, coincident projected vertices, nonfinite
inputs, duplicate controls and over-budget inputs also refuse. Native plate
admission runs before writing an FSB. It checks local incidence and plate
quality, **not general triangle intersections or specimen accuracy**.

### Consistent rib and bridge inertia

Add this optional row to a **flat** FSB or its OBJ import sidecar to integrate
the mass of the bending ribs and bridges along their cubic Hermite displacement:

```text
stiffener-mass,consistent-hermite
```

Without this row, or with `stiffener-mass,lumped`, the original endpoint mass
law is retained. The consistent option uses each existing beam's density,
cross-section area and length. It replaces the two endpoint masses with the
exact integral of the beam's transverse velocity field, including its nodal
tangential slopes. Total physical beam mass is unchanged; no material parameter,
mode capacity or structural geometry is changed. This gives rib and bridge
bending inertia the same interpolation as their existing bending stiffness.

The row survives native FSB → OBJ/sidecar → FSB round trips and reaches ordinary
rendering, modal export, full-vector motion preparation and exterior harmonic
solves through the same board assembler. For example, after adding it to an
exported `model-d.fsb`:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry model-d.fsb \
  --equilibrate-board-mass --edge-cubic-board-mass \
  --board-band-hz 1200 --note 84 --duration 2 --render c6-hermite.wav
```

Panel inertia is selected independently: the usual lumped panel, exact P1
`--consistent-board-mass`, and `--edge-cubic-board-mass` each compose with this
beam option. Complete modal slices still must fit the unchanged 128-mode limit.
The example keeps the existing 1.2 kHz band; higher bands require their own
mesh and bridge-response convergence checks.

This option supplies **translational** Euler–Bernoulli inertia. It does not add
axial, eccentric rotary or torsional beam inertia; the existing offset stiffness
`EI + EAe²` is unchanged. Crowned shells have their own six-DOF beam model and
refuse this flat-board row, including during crown import. The numerical option
does not establish measured Model D mobility or improve a recording by itself.

## External Steinway asset research (checked 2026-09-22)

- **seavenois, Steinway D274**, BlendSwap 7279, page-declared CC0, Blender 2.6x,
  16.2 MB: https://blendswap.com/blend/7279 . The author describes a simplified
  model. The download page requires sign-in. No file bytes were obtained or
  redistributed in this change; topology, scale and physical region names
  therefore remain unverified. This is a visual asset candidate, not a measured
  structural mesh. https://blendswap.com/blend/7279/download
- **sandy2, BWV846Prelude_piano_animation**, BlendSwap 28847, page-declared CC0,
  Blender 2.9x/Cycles, 48.5 MB: https://blendswap.com/blend/28847 . This derives
  from the above D274 and credits CC0 procedural wood and metal materials.
  Its author explicitly warns that visual material names can differ from the
  actual piano part materials. Not downloaded or redistributed here.
- **Steinway Model D specifications**: https://www.steinway.com/pianos/steinway/grand/model-d .
  Published envelope 2.74 by 1.56 m; Sitka spruce panel with 9-to-6 mm taper,
  sugar-pine ribs, maple bridge construction, steel/copper strings and wool
  hammer felt. Species and dimensions do not determine an individual
  instrument's elastic tensors, damping or force-compression curves.
- **Boutillon, Ege and Paulello (2012)**: https://arxiv.org/abs/1210.3948 .
  The existing `steinway_d.rs` reconstructs the drawing's geometric facts;
  its source comments enumerate approximations. Export/import preserves those
  estimates and does not promote them to a measured digital twin.

This implementation connects admitted meshes and physical cards to the
existing modal preparation and pressure-rendering components. Full cabinet/lid scattering, crown/downbearing,
measured felt coupons, measured soundboard FRFs and human listening validation
are separate requirements; this adapter does not claim them.

## Focused regressions

```sh
cargo test -p fs-io obj::tests
cargo test -p fs-couple --example piano_board_import mesh_import::tests
cargo test -p fs-couple --example piano_board_import cli_tests
cargo test -p fs-couple --example piano_board_import mesh_render::tests
```

The tests exercise actual plate admission and modal mass, multi-material mass,
rigid-frame/unit covariance, geometry/material refusal, remapped stiffeners,
all 88 source-derived Model D bridge stations, output overwrite protection,
and an imported-panel -> source felt/shank -> physical-pressure -> WAV chain.
They were authored but not executed in the editing environment (no Rust
compiler/Cargo available); no test-pass or acoustic-validation claim is made.

The added example-only `fs-io` dependency needs a Cargo.lock refresh. The current
editing environment could not run Cargo to regenerate it; an ordinary Cargo
invocation updates the local path-dependency edge, whereas `--locked` refuses
until that refresh is performed. No external package version is intentionally
changed by the new dependency.
