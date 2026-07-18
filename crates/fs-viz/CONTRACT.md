# CONTRACT: fs-viz

Scientific visualization: the verifiable topological-summary primitives.

## Purpose and layer

Layer L5 (LUMEN). Safe Rust; the scoped streamline and contour paths consume L0
`fs-exec::Cx` and the dependency-neutral `fs-blake3` content-identity owner.

## Public types and semantics

- `Vec2` — a 2-D point/vector.
- `streamline(field, seed, dt, steps) -> Vec<Vec2>` — explicitly no-authority
  RK4 compatibility wrapper: valid finite work returns the ordered seed-plus-step
  polyline, while any invalidity/refusal/contained callback unwind returns empty.
- `StreamlineRequest` / `StreamlineBudget` / `StreamlinePlan` — explicit
  dimensionless units, RK4 method/version, finite signed step, optional closed
  domain, boundary/stagnation policies, polling stride, and checked limits for
  steps, field calls, points, bytes, diagnostics, identity, polls, and work.
- `streamline_with_cx(cx, field, request, budget) -> Result<StreamlineOutput,
  StreamlineRunError>` — caller-owned cancellation and ambient deadline/poll/cost
  enforcement, contained callback unwinds, private staging, structured terminal
  evidence, declared early termination, and a domain-separated BLAKE3 identity.
- `CriticalKind` (`Minimum`/`Saddle`/`Maximum`/`Degenerate`), `CriticalPoint {
  kind, morse_index }`, `HessianSymmetryPolicy`, and
  `classify_hessian_with_policy(hessian, eigenvalue_tolerance, symmetry_policy)`
  — the Morse type + index (number of negative eigenvalues) of a nominal
  Hessian admitted under caller-explicit symmetry authority.
  `classify_hessian(hessian, eigenvalue_tolerance)` is the compatibility wrapper
  fixed to `HessianSymmetryPolicy::Exact`. For a fail-closed `Degenerate`
  result caused by invalid input, index 0 is a no-index-claimed sentinel rather
  than evidence of positive inertia.
- `Grid2::from_fn(nx, ny, lo, hi, node_limit, field) -> Result<Grid2,
  Grid2Error>` — an admitted finite scalar field on a genuinely 2-D regular
  grid; `at` and `point` expose its x-fastest nodes.
- `Grid2::isocontour_crossings(iso, crossing_limit) -> Result<Vec<Vec2>,
  IsoContourError>` — bounded, deterministic unique point intersections of the
  finite level with piecewise-linear grid edges; this is the explicitly
  no-cancellation compatibility entry point.
- `IsoContourBudget` / `IsoContourPlan` — the complete caller envelope and its
  checked conservative plan for cells, nodes, edge visits, exact ownership,
  interpolation, output, scratch, diagnostics, live bytes, identity bytes,
  polls, and deterministic work units.
- `Grid2::isocontour_crossings_with_cx(cx, iso, budget) ->
  Result<IsoContourOutput, IsoContourRunError>` — caller-owned cancellation,
  ambient deadline/poll/cost enforcement, private staging, final-checkpoint
  publication, a fixed `IsoContourReport` retaining the requested envelope,
  checked plan, actual/peak counters, disposition, and a domain-separated
  BLAKE3 artifact identity. Actual counters distinguish canonical exact
  crossings, strictly interpolated crossings, and exact candidates discarded
  by static first-edge ownership.
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
- Scoped streamline admission validates the first non-finite seed component,
  finite nonzero signed `dt`, positive poll stride, finite increasing optional
  domain, and in-domain seed before allocation or callback work. Negative `dt`
  is reverse time; zero steps publishes exactly the finite seed; RK4 performs
  exactly four callback attempts per completed step. Every field result and
  constructed stage/state is finite or the operation refuses atomically.
- The complete worst-case plan uses checked arithmetic for `steps + 1`, four
  calls per step, output and peak live bytes, fixed scratch/diagnostics,
  identity bytes, polls, and work. Output capacity is reserved once and checked
  against the allocator-reported capacity before callbacks. `StopBeforeExit`
  and `StopBeforeRepeat` publish only the valid prefix with typed termination;
  `RefuseExit` publishes nothing. Retained repeated points are explicit policy.
- Scoped streamline work polls before allocation, at deterministic complete-step
  chunks, at deterministic identity-point chunks, and immediately before
  publication. Its identity binds method/version, units, seed/dt/step bits,
  domain/termination policies, actual termination/counters, and every point bit
  under `org.frankensim.fs-viz.streamline-rk4.v1`.
- Hessian classification recovers the known Morse type:
  `x²+y²`→min (index 0), `x²−y²`/`xy`→saddle (index 1),
  `−(x²+y²)`→max (index 2); a zero eigenvalue is degenerate. The exact policy
  admits only finite, numerically equal off-diagonal entries (`+0.0` and `-0.0`
  are equal), performs no implicit symmetrization, and fails closed with
  `Degenerate`/index 0 on asymmetry. `eigenvalue_tolerance` applies only to
  zero-eigenvalue classification and never widens symmetry admission.
- A `Grid2` has at least two nodes per axis, an admitted product and allocation,
  finite strictly increasing bounds with finite extent, and finite samples in
  deterministic row-major/x-fastest order. Every generated coordinate is
  strictly increasing; an over-resolved axis whose adjacent logical nodes round
  onto the same `f64` is refused before the field callback runs.
- `isocontour_crossings` of a circle SDF all lie on the circle (to grid
  resolution); a finite level outside the field's range succeeds with no
  crossings. Strict sign changes use overflow-resistant scaled interpolation,
  but the real parameter alone is insufficient evidence: every constant output
  coordinate must retain the endpoint bits and every varying coordinate must be
  strictly between the endpoint coordinates in IEEE-754 total order. A result
  that rounds onto an endpoint or off the closed edge is refused before it can
  consume output-count budget as `UnrepresentableIntersection`, never relabeled
  as exact.
  A single truly exact endpoint is emitted once at its original coordinate bits
  even when incident edges share it. Ownership is the statically selected first
  incident edge in declared traversal order: from the row below, otherwise from
  the left on the first row, otherwise the origin's positive-x edge. This is an
  O(1) decision per edge endpoint, needs no search or marker scratch, and makes
  total extraction work O(nodes + crossings). A wholly exact edge is refused as
  a coincident segment that a point-only result cannot represent.
- Scoped contour planning uses checked arithmetic for the exact node/cell/edge
  geometry, conservative ownership/interpolation maxima, output payload,
  streaming-hasher scratch, fixed diagnostics, simultaneous live bytes,
  identity preimage, poll count, and work units. Every explicit limit is
  admitted before the output allocation; allocator-reported capacity is then
  checked again before edge work. The private vector cannot escape on any
  refusal, cancellation, allocation fault, or unwind.
- For every terminal contour report, `crossings == exact_crossings +
  interpolated_crossings`. `unowned_exact_candidates` counts only exact
  endpoints rejected by the static first-incident-edge rule, while
  `exact_ownership_checks` also includes canonical owners and therefore never
  falls below it.
- Cancellation is polled before allocation, at deterministic edge chunks, at
  deterministic identity chunks, and immediately before publication. At most
  `items_per_poll` edge/identity-point items occur between checkpoints. The
  final artifact identity binds schema version, dimensions, world-bound bits,
  iso bits, output count, and every point bit in traversal order under
  `org.frankensim.fs-viz.isocontour-crossings.v1`.
- A planar `Grid3` level set has exact area and increasing-field winding.
  Sphere area error decreases under refinement, and gyroid extraction is
  indexed, centrally symmetric, and exactly replay-deterministic.
- Scalar-field schema-v1 encoding and decoding are byte-exact and
  replay-deterministic. Byte/sample budgets are checked before allocation;
  dimensions, layout, world geometry, quantity/units, byte length, and sample
  finiteness are validated before any downstream visualization claim.
- All primitives are deterministic when caller callbacks are deterministic.

## Error model

`Grid2Error` distinguishes dimensions below two, node-count overflow, explicit
node-budget refusal, invalid/non-finite/non-increasing bounds, the first
non-representable axis coordinate, the first non-finite x-fastest sample, and
allocation refusal. All layout and budget checks precede callback evaluation;
a panic raised by the caller's field closure remains a caller panic rather than
a `Grid2Error`.
`StreamlineError` distinguishes non-finite seed, invalid signed step, invalid
method version, domain/seed placement, invalid polling, plan overflow, every one-short operation
resource, ambient Cx refusal, allocation refusal, contained callback unwind,
the first non-finite field/stage component, and refused domain exit.
`StreamlineRunError` retains the exact request, caller envelope, checked plan,
completed steps, callback attempts, staged points, requested/peak bytes, polls,
work, termination/refusal disposition, and no-publication state. Callback panic
payloads are not retained; user callback side effects cannot be rolled back.
`Grid2::at` and `Grid2::point` require admitted in-range node indices and panic
on caller indexing errors; they cannot expose extrapolated coordinates.
`IsoContourError` distinguishes non-finite levels, zero/exceeded crossing
budgets, invalid poll strides, checked plan overflow, per-resource operation
budget refusal, ambient Cx cancellation/deadline/poll/cost refusal, exact-level
coincident edges, strict real intersections with no representably interior
point produced by the binary64 interpolation, allocation refusal, and
non-finite interpolation geometry. The representability refusal
retains bounded endpoint
indices; endpoint-coordinate, sample, and iso bits; scaled interpolation
distances and parameter; computed point bits; and the first collapsed axis.
Extraction is all-or-error: it never returns a partial crossing vector, and
malformed evidence never becomes the successful empty result reserved for a
finite absent level.
`IsoContourRunError` retains the typed root error plus a terminal report. A
pre-plan refusal has `plan: None`; every post-plan refusal retains the checked
requirements and exact completed counters. Only `Completed` carries
`published=true` and a nonempty artifact identity.
`Grid3` construction is fallible and refuses degenerate/overflowing dimensions,
invalid or non-finite bounds/samples, length mismatch, node-budget excess, and
allocation refusal. Isosurface extraction refuses non-finite levels, a zero or
exceeded triangle budget, `u32` index exhaustion, and non-finite geometry. It
never returns a silently truncated mesh.
`ScalarField3Error` distinguishes sample/byte budget refusal, malformed or
unsupported schema bytes, ambiguous semantics, invalid geometry, non-finite
values, allocation refusal, and node/cell layout mismatch.
Hessian classification is total and fail-closed: non-finite matrix entries,
invalid eigenvalue tolerances, and off-diagonal disagreement under the selected
policy return `Degenerate` with index 0 rather than a confident Morse claim.

## Determinism class

Fully deterministic for deterministic callbacks: RK4, classification, and
contouring use fixed operation/traversal order. Streamline reports and identities
are bit-stable for the same request and callback results; negative-time and
declared early-termination semantics are identity-bound.
Grid2 sampling is row-major/x-fastest; crossing traversal is row-major with the
positive-x edge before positive-y, and shared exact nodes use the statically
derived first-incident edge, retaining first-encounter order without mutable
deduplication state. Scoped poll sites, reports, and artifact preimages are
derived from logical traversal counts rather than scheduler timing; retry under
a fresh Cx is byte-identical.
Grid3 sampling is z/y/x with x-fastest storage; isosurface traversal is
z/y/x/cube-tetrahedron order and uses an ordered edge cache.
Scalar-field artifacts use fixed little-endian IEEE-754 f64 bits and fixed
length-prefixed UTF-8 semantics; their bytes are cross-ISA stable.

## Cancellation behavior

`streamline_with_cx` observes the caller-owned `Cx` before allocation, between
bounded groups of complete RK4 steps, between bounded identity-point groups,
and at final publication. Each RK4 step is indivisible and contains four field
calls. Completed chunks are charged to the ambient admitted cost plan; any
cancellation/budget refusal drops the private vector and returns terminal
no-publication evidence. The compatibility wrapper has no cancellation authority.

`isocontour_crossings_with_cx` observes the caller-owned `Cx` before allocation,
at bounded edge/identity chunks, and at the final publication cutoff. It admits
the checked worst-case work against the ambient `fs-exec::AdmittedBudget`,
charges completed chunks, and never publishes on cancellation or any other
budget refusal. The small `isocontour_crossings` compatibility entry point
derives the same complete resource envelope but deliberately has no
cancellation authority.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/viz.rs`: a rotation field streams along a circle; a saddle
field conserves the hyperbola invariant; scoped streamline zero-step/exact-plan
receipts, request overflow/invalidity, every one-short operation resource,
negative-time, coordinate-scaling, RK4-refinement, and polling-chunk
metamorphisms, boundary/stagnation policy,
pre-requested and injected cancellation, ambient deadline/cost/poll refusal, callback
panic/non-finite results, callback-requested cancellation observed at the next
indivisible complete-step boundary, arithmetic overflow, allocation fault, and identical
retry; Hessian classification recovers the Morse type, retains inertia under
large finite positive scaling, and exercises exact/asymmetric/near-mismatch,
transpose, power-of-two scale, signed-zero, and non-finite symmetry admission;
a circle-SDF isocontour lies on the circle (+ empty out-of-range);
2-D grid sampling/addressing; fail-before-sampling dimension, budget, bounds,
overflow, coordinate-collapse, and non-finite-value admission; non-finite-level
and crossing-budget refusal; exact-node ownership and coincident-edge refusal;
the adversarial exact/non-exact checkerboard under exact and one-short output
budgets, both axis shapes, and byte-identical replay;
complete scoped-plan accounting, exact/interpolated/static-owner telemetry,
and artifact identity; exact and one-short
refusal for every resource before polling/edge work; pre-requested and
mid-traversal cancellation; ambient deadline-without-clock, cost-plan, and poll
exhaustion refusals; injected allocation and
checkpoint-panic faults; atomic retry equivalence; geometric endpoint reversal
and affine-translation metamorphisms;
extreme finite and subnormal interpolation; exact next-up/down strict-crossing
collapse refusal before output budgeting with endpoint/value/iso/t/point-bit
evidence; G3 axis, endpoint-sign, signed-zero, and power-of-two neighbor refusal
metamorphisms plus representable sign-inversion and positive-power-of-two scaling
equivalence; determinism; exact oriented plane
extraction;
sphere-area refinement; gyroid symmetry/indexing/replay; and fail-before-work
Grid3 budget/non-finite admission. The scalar-field artifact battery covers
node-centered byte replay and Grid3 conversion, an honest one-cell-thick
cell-centered LBM shape, byte/sample admission, schema-version refusal,
truncation, semantic validation, and non-finite payload rejection.

## No-claim boundaries

- v0 is the ANALYTICALLY-VERIFIABLE core: RK4 streamlines, Morse critical-point
  classification, 2-D edge crossings, and regular-grid marching tetrahedra.
  Hessian classification accepts exact numeric symmetry only; it does not infer
  equality of mixed partials before binary64 rounding, average disagreeing
  entries, attach a symmetrization receipt, or provide absolute, relative, ULP,
  or backward-error symmetry admission. It does not prove that the gradient
  vanishes, that the supplied matrix is an analytic Hessian, or that a continuum
  critical point exists. `Degenerate` includes refused/unclassifiable numeric
  input and is not itself evidence that a valid critical point is mathematically
  degenerate.
  Grid2 crossings are points only: they do not provide marching-squares
  connectivity, resolve saddle cells, represent plateau segments/regions, or
  certify the unsampled continuum field or its topology.
  DUAL contouring with sharp-feature QEF placement, DIRECT VOLUME RENDERING
  with preintegrated transfer functions, LINE-INTEGRAL CONVOLUTION,
  tensor/stress-ellipsoid glyphs, and full MORSE–SMALE complexes / Reeb graphs
  with persistence thresholding are staged.
- Contour byte budgets cover logical Rust payloads, the borrowed `Grid2`
  representation and retained sample-vector capacity, fixed reports, and
  allocator-reported output-vector capacity. They do not claim knowledge of
  allocator bookkeeping, virtual-memory page
  granularity, or the caller's unrelated live heap. Scoped contour extraction
  is sequentially tiled; it does not yet claim parallel speedup or the
  reference-hardware 200-microsecond wall-clock cancellation target without a
  retained measurement.
- Fixed-step streamline RK4 v1 has no embedded local-error estimator (the
  structured `error_estimate` is explicitly `None`), adaptive step control,
  dense output, event root localization, stiffness claim, physical unit
  conversion, or callback-purity enforcement. Only dimensionless units are
  admitted in v1. Cancellation latency excludes time spent inside one caller
  callback; callback side effects and panic-hook output are caller concerns.
  Callback panic containment requires Rust unwinding and cannot intercept a
  process configured to abort on panic.
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
