# A finite acoustic body from the actual soundboard

`piano_exterior` can now build its moving acoustic surface from the **same
admitted structural geometry and section thicknesses** used to prepare the
piano. This removes the need to hand-author a separate two-sided board OBJ.
It does not invent a cabinet, lid, frame, or measured material data.

Use `board-skin-continuous` in the BODY position for a smooth panel whose
thickness is represented by cellwise structural sections. This is an explicit
volume-preserving reconstruction, not silent smoothing or a failure fallback:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --mesh-divisions 4 --dump-geometry model-d.fsb

cargo run --release -p fs-couple --example piano_exterior -- \
  export-skin model-d.fsb steinway-d \
  crates/fs-couple/examples/grand_piano/model-d-section-skin.fspe \
  model-d-skin.obj --continuous-thickness

cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-section-skin.fspe \
  model-d.wav 2 --note 69 --velocity 1 --modes 128 --substeps 8
```

These are calls to the physical preparation/render path, not a promise that
all requested BEM discretizations or passive fits will meet their existing
error bounds. Export runs structural admission, equilibrium when declared, and
modal preparation, but **no BEM or receiver fit**. A successful OBJ export is
not an acoustic-convergence or successful-render certificate. Use fresh output
paths. With a settled board, reuse the exact scale that supplied its equilibrium.

The generated geometry remains the source-derived Model D approximation in
`steinway_d.rs`: approximate transcription of Boutillon, Ege and Paulello,
Acoustics 2012, Figure 2 / Table 1, and an estimated taper using published
Steinway thickness endpoints. It is not a newly acquired factory scan or full
piano CAD. Structural ribs and bridges still affect the modes, but their exposed
acoustic solids are not added by this panel-skin operation.

## Two deliberately different thickness images

`board-skin` retains the piecewise-constant sections as columns. For a facet
normal with vertical component `n_z`, its vertical half-height is `h/(2*n_z)`.
This gives normal separation `h` and volume `facet_area*h`. Exposed thickness
steps receive walls; internal coincident faces are omitted. The rim and holes
are closed along their actual boundaries. Every incident height splits shared
vertical edges, including third-facet junctions. A pinched, nonmanifold step
refuses rather than being welded, smoothed, or capped. There is no implicit
switch to the continuous image.

`board-skin-continuous` explicitly reconstructs a continuous piecewise-linear
vertical half-height at each original node:

```
z_node = sum_incident(facet_area * section_thickness)
         / (2 * sum_incident(projected_XY_area))
```

This is a positive mass-lumped projection, not a global average. Before the
stated height rounding, the resulting closed panel's total volume is exactly
`sum(facet_area * section_thickness)`. It preserves local taper variation and
holes, but changes individual acoustic facets' mean thickness. The **maximum
facet mean-thickness discrepancy and both volumes are reported**. Structural
thicknesses, mass, stiffness, grain axes, damping and eigenvectors never change.
A reconstructed smooth surface is not a unique recovery of manufactured shape.

This distinction matters for the Model D preset: interpreting its smooth-taper
samples as literal steps introduces pinched edges and many tiny walls. An
independent transcription of the divisions-4 preset gives 1,788 continuous-skin
panels versus 9,012 step panels, with a maximum facet mean-thickness change of
about 0.824 mm. That numerical check is not native Rust or acoustic evidence.
Divisions 6 and 8 give 2,340 and 2,892 continuous panels: **they exceed the
existing 2,048-panel BEM wrapper limit**. The exporter allows up to 250,000
panels for inspection, but playback does not silently decimate them.

Both constructions use a declared **1 nm vertical-height grid** to avoid
sub-roundoff seams. No source file or structural section is modified. Output
closure and volume are checked; unresolved tiny panels, missing sections,
topology mismatch, folded graphs, excessive offsets and budget overruns refuse.
The operation assumes the admitted midsurface is a non-self-overlapping shallow
XY graph. It is not a general CAD Boolean or self-intersection certificate.

## Units, labels, and motion

For either generated BODY keyword, the acoustic specification must contain:

```text
obj-scale-m,1
obj-origin,0,0,0
moving,soundboard_skin
```

No other moving/rigid labels are allowed, because there is no separate component
to satisfy them. The rest of the existing explicit medium, receiver, band,
fit-order, offset and PCM-scale rows still apply. The included `.fspe` is an
estimated study setup, not calibrated air or microphone measurements.

Generation occurs **after** the existing static equilibrium preparation, so
an admitted loaded crown produces a loaded skin rather than a flattened or
unloaded substitute. Each acoustic panel retains its known source facet,
barycentric coordinates and through-thickness arm. Positive degree-two
quadrature integrates translation plus physical rotation at those sites.
There is no nearest-node assignment or extrapolation beyond the panel.

The existing string-mass-loaded transformation is applied once. `response`,
`admittance`, `render`, and `render-loaded` then use their original BEM, contact,
acoustic-feedback, score, stereo and PCM owners. All prior material/finite-hammer,
nonlinear-string and pedal options remain available on their supported commands.
An invalid OBJ FILE cannot select generated geometry merely by containing a
keyword; the BODY keyword is interpreted only as an explicit command argument.

`export-skin` defaults to step geometry; add `--continuous-thickness` to select
continuous reconstruction. The OBJ records its geometry choice, attribution,
height-grid and thickness/volume discrepancy. Its mesh can be inspected or used
in external asset work. An OBJ alone does not retain the native facet embeddings:
reimport uses the ordinary projection contract and may refuse crowned rim
points or differ at facet boundaries. Use the BODY keyword for the native
known-site motion path.

No fitted acoustic error limit, mechanical energy check, mode cap or pressure
scale is relaxed. Rigid cabinet/lid scattering still requires an explicitly
supplied complete acoustic mesh; this operation does not fill that missing
asset. Full-band realism, measured Model D fidelity and real-time operation
remain unclaimed.
