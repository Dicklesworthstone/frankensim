# Imported lid and cabinet poses in the actual acoustic solve

`--rigid-assembly parts.fspr` adds explicitly selected, positioned OBJ parts to
`piano_exterior response`, `admittance`, `render`, or `render-loaded`. It works
with either a supplied complete BODY OBJ or the native `board-skin` and
`board-skin-continuous` selections. In particular, the latter keep the generated
board's exact source-facet and rotation mapping; no second, hand-modelled board
surface or nearest-facet remapping is required to add an imported lid.

All added surfaces are **acoustically rigid** and have zero prescribed velocity.
They participate in the same BEM solve as the moving board. Their pose changes
both the microphone transfer and the radiation impedance seen by loaded
playback. No extra source, empirical cabinet IR, equalizer, or synthetic stereo
spread stands in for scattering. In ordinary one-way `render`, rigid geometry
changes the observation, not the mechanical trajectory. In `render-loaded`,
its contribution to the admitted passive acoustic load also reacts on the piano.

## Supply the geometry and its placement

An assembly file contains one source row and a `part`/`pose` pair per alias:

```text
frankensim-piano-rigid-assembly-v1
source,estimated,Illustrative placement only; replace with supplied asset attribution
part,lid,assets/piano.obj,Lid,0.001,0,0,0
pose,lid,0,0,0.1,1,0,0,30,0,0,0
```

This example does **not** supply the OBJ, identify a real Steinway hinge, or
make its numerical pose a factory measurement. The row formats are:

```text
part,alias,obj_path,exact_object_or_group_label,metres_per_obj_unit,origin_x,origin_y,origin_z
pose,alias,pivot_x,pivot_y,pivot_z,axis_x,axis_y,axis_z,angle_degrees,translate_x,translate_y,translate_z
```

The OBJ origin is in **source units**. Pivot and translation are in **metres in
the board frame**, and the axis must be a supplied unit vector. Positive angles
follow the right-hand rule. The transform, with no implicit axis guess, is:

```text
point = translation + pivot
      + rotation(axis, angle) * (scale * (obj_point - obj_origin) - pivot)
```

For a fixed part, specify a valid axis and zero angle rather than omitting its
pose. Each alias has exactly one selection and one pose. Object/group labels
are exact and case-sensitive; spaces are retained, CSV commas are not supported
inside paths or labels. The chosen faces must form complete, outward closed
components. Omitted artist-model faces are counted and reported, not silently
classified as physical parts. Distinct poses can instantiate the same source
part. Coincident duplicate geometry is rejected, not automatically welded.

Relative OBJ paths resolve against the assembly file's directory. Each unique
OBJ path is read once; its admitted coordinates are retained before structural
modal preparation. MTL references, texture paths and other nested files are
never opened. Visual wood/metal/lacquer parameters are not constitutive laws.
Use only assets you have permission to use; this command does not download,
license, authenticate, or redistribute a third-party model.

## Inspect and play

```sh
# Export the posed rigid parts in SI coordinates, without a board or BEM solve.
cargo run --release -p fs-couple --example piano_exterior -- \
  export-rigid parts.fspr posed-parts.obj

# Preserve native panel geometry and motion while adding the supplied assembly.
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-section-skin.fspe \
  piano-with-lid.wav 2 --rigid-assembly parts.fspr --note 69 --velocity 1

# The same geometry choice in a harmonic bridge-force experiment.
cargo run --release -p fs-couple --example piano_exterior -- \
  admittance model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-section-skin.fspe \
  69 bridge-with-lid.csv --rigid-assembly parts.fspr
```

The native skin's acoustic specification still declares identity OBJ units and
`moving,soundboard_skin` only. Added parts use their own assembly selections;
do not add unsatisfied `rigid` rules for them to the board-only specification.
A supplied BODY OBJ retains its original explicit moving/rigid mapping. Avoid
including the same physical part in both BODY and the extra assembly.

The existing complete material maps, finite hammer faces, nonlinear strings,
MIDI/CSV force performances, pedals, modal budgets and mono/stereo playback
remain available. The original scale must match any preloaded board equilibrium.
Reports/CSV include the source paths, exact labels, selected/excluded counts,
units and poses actually used. Output paths must be new. A rejected asset,
fit or scene budget never selects a cabinet-free fallback or partial WAV.

## Boundaries

The **combined board plus rigid parts** must fit the existing 2,048-panel BEM
limit; source OBJ text is bounded to 32 MiB total and assemblies to 64 parts.
A divisions-4 native continuous Model D skin consumes 1,788 panels, leaving
260 for additional components. A large artwork mesh is not automatically
coarsened. Prepare an appropriate acoustic-resolution asset explicitly; do not
expect a high-poly rendering model to fit the dense acoustic budget unchanged.

This is a static, disjoint-component exterior scattering model. A lid can be
posed before preparation, not opened during a performance. It has no flexible
vibration, absorption, mechanical hinge inertia, added soundboard mass or glue
connection. Geometric validity checks cover faces, winding, exact seams and
closed volume, **not global intersections, overlap or cavity accessibility**.
Posed parts must not intersect or seal the board into an unmodelled interior
problem. No collision repair or CAD Boolean union is performed.

The enclosing geometry and receiver flight admission are recomputed after
placement. Near receivers may consequently refuse. Existing BEM wavelength,
passive-load fitting, receiver fitting, time-step, energy and pressure-scale
limits remain unchanged. A successful export is geometry inspection, not proof
of acoustic convergence or a successful played render. This increment adds no
new third-party Steinway mesh bytes or factory material measurements.
