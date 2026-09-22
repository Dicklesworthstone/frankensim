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
cargo run --release -p fs-couple --example grand_piano -- \
  --board-geometry imported.fsb --scale strings.csv --render piano.wav
```

The last command requires a supplied string scale. For a source-derived Model D
scale, the existing `grand_piano --preset steinway-d --dump-scale strings.csv`
command exports one. Explicit `--hammers` cards can be supplied to the render;
selecting an imported board alone does not implicitly enable preset voicing.

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

The realistic-geometry route is now executable once an admitted mesh and
physical cards are supplied. Full cabinet/lid scattering, crown/downbearing,
measured felt coupons, measured soundboard FRFs and human listening validation
are separate requirements; this adapter does not claim them.

## Focused regressions

```sh
cargo test -p fs-io obj::tests
cargo test -p fs-couple --example piano_board_import mesh_import::tests
cargo test -p fs-couple --example piano_board_import cli_tests
```

The tests exercise actual plate admission and modal mass, multi-material mass,
rigid-frame/unit covariance, geometry/material refusal, remapped stiffeners,
all 88 source-derived Model D bridge stations, and output overwrite protection.
They were authored but not executed in the editing environment (no Rust
compiler/Cargo available); no test-pass or acoustic-validation claim is made.
