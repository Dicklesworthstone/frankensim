# Spatial fields in Observatory

The same local workspace now reads **actual mesh and field exports**, rather than
limiting visualization to the hard-coded concentric section. Start as before:

```sh
node apps/observatory/server.mjs --open
```

## Open an existing native export

Drop an ASCII `.vtu` file emitted by `fs-viz::vtu::VtuWriter::write_ascii` into the
workspace, or use **Import a run**. Declare the original coordinate unit in the
import dialog: metres, centimetres or millimetres. Only coordinates are converted;
field `Unit` attributes are retained without rescaling. A field without a unit is
shown as `unspecified`. No material properties or physical units are inferred from
an attractive color, filename or model label.

The reader supports the existing writer's single-piece `UnstructuredGrid` subset:
triangles, planar convex quads, tetrahedra and convex planar-faced hexahedra.
Point and cell fields remain separate. Scalar and three-component vector arrays
retain their source entity indices and full-precision values. Volume meshes are
shown through their **exterior faces**, with original cell connectivity retained;
interior cell values remain available by entering their original probe index.
This is not a cutting-plane or volume renderer. Mixed volume/surface cells,
nonplanar quads, binary/appended/compressed VTU, VTR, PVD and XDMF are not supported.
Unsupported files are refused, not partially imported or silently repaired.

## Inspect a field

Select a field and, for vectors, magnitude or a signed X/Y/Z component. Click a
surface element to probe an original entity: cell fields select its owning volume
or surface cell; point fields select the nearest barycentric vertex. The readout
shows the exact source sample, not an interpolated color value. Enter an original
point/cell index to inspect interior values as well.

Orbit by dragging, pan with Shift-drag, zoom using the wheel or buttons, or focus
the viewport and use the arrow keys, +/- and 0. The camera fits the reference mesh
without distorting its aspect ratio. Element edges can be toggled. The fixed color
scale spans **all recorded frames**, so heating or cooling cannot be hidden by
per-frame auto-normalization. Point colors use triangle interpolation only for
display; cell fields remain discontinuous across their interfaces.

The WebGL renderer uploads normalized floating-point coordinates for display while
retaining the source doubles for probes and export. A bounded, depth-tested
software renderer uses the same vertex data, camera and color ramp on hosts where
WebGL is unavailable. Its raster is capped at 400,000 pixels and 24 million pixel
visits per draw. Extreme overlap produces an explicit refusal, not dropped faces;
source data and probes remain usable. It is a fallback, not a GPU-speed claim.
Neither renderer is a simulation kernel or a mesh validity certificate.

For deformation, explicitly choose a **three-component point field in metres** in
the displacement selector. Set its display multiplier; zero displays reference
geometry. Imported JSON may predeclare `mesh.displacement_field`, but a VTU field
name is never guessed to mean displacement. A multiplier changes only rendering,
not the stored geometry or result. Excessive displacement falls back to reference
geometry with a warning; it never loses the accepted field history.

## Replay a recorded sequence

Convert already-exported native snapshots without running a replacement solver:

```sh
node apps/observatory/tools/vtu-to-observatory.mjs \
  --units m --output field-sequence.json --times 0,0.1,0.25 \
  frame-0.vtu frame-1.vtu frame-2.vtu
```

Then import `field-sequence.json`. The output path must be new. Times are supplied
explicitly, never guessed from filenames. Geometry, original connectivity,
field names/associations/components/units and source description must match.
Changing-mesh sequences must be imported as separate runs; no remeshing,
resampling, entity matching or new thermal/fluid state is invented.

Play/pause, previous/next, speed and the frame slider inspect **recorded time**.
Between records playback holds the last recorded frame. It does not integrate or
interpolate physics. Playback pauses at the end, when switching runs and when the
page is hidden. The existing signal panel becomes an exact probe history; clicking
it or using its left/right keys selects recorded frames. One static export never
produces an invented time history. Static time zero is only a viewer reference.

Pin a run to compare the same field/entity on exactly matching source geometry.
Cell comparisons require original source connectivity; a matching visible exterior
alone is not enough to establish interior cell identity. Incompatible comparisons
are labeled rather than forced onto the same plot. Probe charts connect samples
for visual reading; playback and probe readouts always select source records.

## Keep and share results

Mesh runs use the existing run library. JSON export retains every admitted source
coordinate, original cell, frame and field value. CSV exports the selected
entity/component's complete probe history. PNG includes the source classification,
field, color limits, frame and displacement multiplier. Restored/imported results
remain labeled as such; source assertions are not authenticated certificates.

Maximum import: **32 MiB**, 50,000 points, 50,000 original cells, 100,000 exterior
triangles, 512 frames, 24 fields and 2 million total numeric values. These are
explicit browser working-set limits, not a claim of unbounded production-scale
postprocessing. Signal CSV and signal JSON retain their original 8 MiB limit.
Large runs stay available in memory even if they do not fit the optional 3 MiB
browser storage budget; smaller runs are not evicted merely by importing a larger
one. Export JSON for durable retention.

## Focused tests

```sh
node --test apps/observatory/tests/mesh*.test.mjs
python apps/observatory/tests/mesh_browser.py --url http://127.0.0.1:7070
```

The Node tests include real loopback asset serving and execution of the actual
converter. The browser test exercises the real UI with explicitly synthetic field
fixtures, not native Rust solver output. `--in-memory` is available on test hosts
that prohibit browser navigation; it rewrites only module/stylesheet URLs and
cannot establish served-asset/CSP integration or durable browser persistence.
`--screenshots DIR` records the actual viewer with its synthetic-fixture label.
GPU and software execution are reported separately. No browser security policy is
disabled. Native Rust execution, large production meshes and other browsers remain
separate validations.
