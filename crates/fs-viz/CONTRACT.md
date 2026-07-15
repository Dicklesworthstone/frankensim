# CONTRACT: fs-viz

Scientific visualization: the verifiable topological-summary primitives.

## Purpose and layer

Layer L5 (LUMEN). No dependencies — pure Rust.

## Public types and semantics

- `Vec2` — a 2-D point/vector.
- `streamline(field, seed, dt, steps) -> Vec<Vec2>` — RK4 integration of a flow
  field into an ordered polyline (seed included).
- `CriticalKind` (`Minimum`/`Saddle`/`Maximum`/`Degenerate`), `CriticalPoint {
  kind, morse_index }`, `classify_hessian(hessian, tol)` — the Morse type + index
  (number of negative Hessian eigenvalues) of a critical point.
- `Grid2::from_fn(nx, ny, lo, hi, f)` — a scalar field on a regular grid; `at`,
  `point`, `isocontour_crossings(iso)` (linearly-interpolated edge crossings).
- `Grid3::from_fn(dimensions, lower, upper, node_limit, field)` and
  `Grid3::from_values(...)` — owned finite scalar samples with x-fastest
  addressing, strict world bounds, and an explicit node budget.
- `Grid3::isosurface(iso, triangle_limit) -> Result<IsoMesh3, IsoSurfaceError>`
  — deterministic six-tetrahedra-per-cell extraction with canonical global
  edge/node vertex sharing, outward winding from `< iso` toward `>= iso`, and
  an explicit all-or-error triangle budget. `IsoMesh3::into_parts` yields the
  renderer-ready indexed arrays; `surface_area` measures the PL surface.
- `ScalarField3` / `ScalarFieldSemantics` / `ScalarLayout3` — a versioned,
  dependency-free scalar-field artifact codec with explicit quantity, domain
  unit, value unit, node-vs-cell centering, x-fastest dimensions, origin, and
  spacing. `SCALAR_FIELD3_ARTIFACT_KIND` and
  `SCALAR_FIELD3_SCHEMA_VERSION` let an L6 caller validate the ledger envelope
  before bounded decoding. Node-centered fields convert to `Grid3`;
  cell-centered fields retain one-cell-thick LBM slabs without inventing fake
  nodes.

## Invariants

- STREAMLINES honor the flow: a rotation field `(-y, x)` conserves the radius (a
  circle); a saddle field `(x, -y)` conserves `xy` (a hyperbola).
- `classify_hessian` recovers the known Morse type: `x²+y²`→min (index 0),
  `x²−y²`/`xy`→saddle (index 1), `−(x²+y²)`→max (index 2); a zero eigenvalue is
  degenerate.
- `isocontour_crossings` of a circle SDF all lie on the circle (to grid
  resolution); a level set outside the field's range yields no crossings.
- A planar `Grid3` level set has exact area and increasing-field winding.
  Sphere area error decreases under refinement, and gyroid extraction is
  indexed, centrally symmetric, and exactly replay-deterministic.
- Scalar-field schema-v1 encoding and decoding are byte-exact and
  replay-deterministic. Byte/sample budgets are checked before allocation;
  dimensions, layout, world geometry, quantity/units, byte length, and sample
  finiteness are validated before any downstream visualization claim.
- All primitives are deterministic.

## Error model

`Grid2::from_fn` panics only on a degenerate grid (`< 2` points per axis).
`Grid3` construction is fallible and refuses degenerate/overflowing dimensions,
invalid or non-finite bounds/samples, length mismatch, node-budget excess, and
allocation refusal. Isosurface extraction refuses non-finite levels, a zero or
exceeded triangle budget, `u32` index exhaustion, and non-finite geometry. It
never returns a silently truncated mesh.
`ScalarField3Error` distinguishes sample/byte budget refusal, malformed or
unsupported schema bytes, ambiguous semantics, invalid geometry, non-finite
values, allocation refusal, and node/cell layout mismatch.

## Determinism class

Fully deterministic: RK4, classification, and contouring are pure functions.
Grid3 sampling is z/y/x with x-fastest storage; isosurface traversal is
z/y/x/cube-tetrahedron order and uses an ordered edge cache.
Scalar-field artifacts use fixed little-endian IEEE-754 f64 bits and fixed
length-prefixed UTF-8 semantics; their bytes are cross-ISA stable.

## Cancellation behavior

None here; the production viz shares LUMEN tiling + progressive streaming with
`Cx` cancellation.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/viz.rs`: a rotation field streams along a circle; a saddle
field conserves the hyperbola invariant; Hessian classification recovers the
Morse type; a circle-SDF isocontour lies on the circle (+ empty out-of-range);
2-D grid sampling/addressing; determinism; exact oriented plane extraction;
sphere-area refinement; gyroid symmetry/indexing/replay; and fail-before-work
Grid3 budget/non-finite admission. The scalar-field artifact battery covers
node-centered byte replay and Grid3 conversion, an honest one-cell-thick
cell-centered LBM shape, byte/sample admission, schema-version refusal,
truncation, semantic validation, and non-finite payload rejection.

## No-claim boundaries

- v0 is the ANALYTICALLY-VERIFIABLE core: RK4 streamlines, Morse critical-point
  classification, 2-D edge crossings, and regular-grid marching tetrahedra.
  DUAL contouring with sharp-feature QEF placement, DIRECT VOLUME RENDERING
  with preintegrated transfer functions, LINE-INTEGRAL CONVOLUTION,
  tensor/stress-ellipsoid glyphs, and full MORSE–SMALE complexes / Reeb graphs
  with persistence thresholding are staged.
- `IsoMesh3` is the piecewise-linear isosurface of trilinearly sampled node
  data under a fixed tetrahedralization. It does not claim topology recovery
  below grid resolution, sharp-feature preservation, Hermite normals,
  watertightness when the surface intersects the domain boundary, or an error
  certificate against an unsampled continuum field.
- The API does not read ledgers or depend on L6. A higher-layer orchestrator
  must compare the ledger artifact kind, bounded-read the bytes by content
  hash, decode this versioned schema, call the appropriate L5 primitive, and
  bind source/output content hashes into lineage. The codec carries quantity
  and units but makes no claim that a ledger hash exists or that the stated
  semantics are physically appropriate.
- Adaptive/embedded-pair integration (fs-time steppers) and Qty-labeled
  perceptually-uniform colormaps are staged with the rendering integration.
