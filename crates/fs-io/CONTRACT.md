# fs-io — CONTRACT

Import/export with QUARANTINE (plan patch Rev J): the world boundary.
Dirty geometry comes in, useful artifacts go out — and no imported
artifact becomes a trusted value without a certification receipt.

Ambition tags: STL/OBJ/PLY + quarantine + catalogs + 3MF/GLB/VTK [S];
bounded STEP Part-21 syntax, strict native triangular faceted-resource
decoding/re-emission, and estimated-SDF handoff [S]; broader CAD/EXPRESS
interpretation, surface tessellation, and SDF-to-NURBS/B-rep re-fit export
explicitly STAGED (no-claim below).

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, fs-rep-mesh (repair +
half-edge validity), fs-rep-sdf, fs-exec, fs-evidence, fs-geom, fs-obs,
fs-math. PNG/EXR export is fs-img's (L5). Ledger `imports` rows are
written HELM-side from the receipt JSON this crate emits — L2 never calls
L6. Consumers: the P4 frame flagship (AISC catalogs), fs-fab.

## Public types and semantics

- **Imports** (`stl`, `obj`, `ply` modules): binary STL (auto-detected
  by exact sizing, so binary files beginning with "solid" still parse)
  + ASCII STL; OBJ subset (`v`/`f`, fan triangulation, negative indices,
  `v/vt/vn` forms with vt/vn ignored); PLY ASCII + binary_little_endian
  (vertex x/y/z any scalar type, face index lists; other
  elements/properties skipped with correct stride accounting). Every
  parser: element-capped (`MAX_ELEMENTS`), length-checked, non-finite
  coordinates refused, structured `IoError` — never a panic.
- **Quarantine** (`quarantine` module): `import_mesh` → `Quarantined
  { raw, source_receipt, defects }`. The census detects degenerate
  faces, duplicate faces, unreferenced vertices, and non-manifold-or-
  open surfaces — the latter by DIRECT EDGE COUNTING (every undirected
  edge of a watertight 2-manifold appears exactly twice) because the
  half-edge builder alone legally accepts open boundaries (a real trust
  gap the conformance suite caught during development). `promote` runs
  the fs-rep-mesh repair suite, re-censuses, and either yields
  `Evidence<Soup>` (exact numerics, receipt-chained provenance) plus the
  `trust: promoted` receipt JSON, or a `PromotionRefusal` with blocking
  defects, ACTIONABLE fixes, and a `trust: refused` receipt.
- **Exports** (`export` + format modules): binary STL / OBJ / ASCII PLY
  (deterministic; OBJ and PLY carry f64 shortest-round-trip text, exact
  on re-import); 3MF (minimal OPC ZIP with STORED entries, fixed
  timestamps for byte determinism); GLB (glTF 2.0 binary container,
  f32 positions + u32 indices, chunk-accounted); legacy-ASCII VTK
  unstructured grid with optional scalar point field.
- **Catalogs** (`catalog` module): CSV (RFC-4180 subset with quoted
  fields and `""` escapes) and JSON (minimal array-of-flat-objects
  reader) validated against a `Schema` of `ColumnSpec`s (Text /
  bounded Number, required flags). Violations name the 1-based data
  row, column, offending text, and the expectation; missing header
  columns list what WAS found.
- **STEP structure** (`step` module): bounded, ASCII-only parsing of the
  ISO-10303-21 clear-text envelope, mandatory `FILE_DESCRIPTION`,
  `FILE_NAME`, and `FILE_SCHEMA` header records, simple and complex DATA
  instances, aggregates, typed parameters, strings, enumerations,
  numeric tokens, and forward references. Parsing rejects duplicate or
  dangling instance IDs after the whole DATA section is known. Canonical
  writing sorts instances by numeric ID, preserves parameter/component
  order, doubles string apostrophes, and revalidates caller-constructed
  documents before emitting bytes. `require_declared_schema` supplies an
  exact, case-insensitive declaration gate without treating a schema label
  as conformance evidence. The sealed `ParsedStep` keeps its immutable
  receipt from becoming stale; `StepStructureReceipt` records syntax/crate
  versions, exact admission limits, non-cryptographic source/canonical-layout
  FNV fingerprints, schemas, graph counts, and a strictly non-authoritative
  AP203/AP214 label hint. HELM must replace fingerprints with its
  collision-resistant artifact identity before authority-bearing use.
- **STEP tessellation handoff** (`step_import` module): accepts a materialized
  triangle soup only alongside the sealed `ParsedStep`, an explicit
  adapter name/version/configuration fingerprint, a declared
  tessellation-deviation certificate, one shared length-unit ID for coordinates,
  deviation, and sampling spacing, a positive sampling spacing, and `Cx`.
  It removes duplicate/degenerate faces and unreferenced vertices, may unify
  orientation, and then refuses every residual boundary, non-manifold edge,
  orientation conflict, or disconnected closed vertex link with a bounded
  deterministic defect prefix.
  Publication yields a sealed `StepImportOutcome`: `Evidence<TiledSdf>` plus
  a source-bound receipt that separately retains tessellation deviation,
  mesh-to-SDF numerical evidence, their outward-rounded combined estimate,
  repairs, quality counters, and adapter identity.
- **Strict native STEP faceted decoding** (`step_faceted` module): materializes
  one caller-selected, root-reachable `FACETED_BREP -> CLOSED_SHELL -> (FACE |
  FACE_SURFACE) -> FACE_OUTER_BOUND -> POLY_LOOP -> CARTESIAN_POINT` closure.
  Plane-backed faces must resolve through `PLANE -> AXIS2_PLACEMENT_3D` to a
  3-D location and optional valid directions; the decoder checks triangle
  coplanarity and winding against `same_sense`. It admits exactly one triangular
  outer loop per face, preserves loop order except for explicit `.F.` bound
  reversal, canonicalizes EXPRESS `SET` face traversal by numeric instance ID,
  and passes the resulting soup into the existing topology/SDF handoff. A
  separate decoder receipt retains the exact admitted schema label, root/shell
  IDs, syntax and semantic fingerprints, face-profile counts, resource limits,
  decimal-to-f64 conversion, plane-consistency, and their combined estimated
  spatial deviation. This is bounded resource-entity decoding, not AP203/AP214
  conformance; callers supply the length unit because the admitted closure
  deliberately excludes representation context.
- **Strict native STEP faceted re-emission** (`step_faceted_export` module):
  accepts only an immutable decoded resource made entirely of bare `FACE`
  entities, assigns canonical dense IDs, emits the pinned triangular
  `FACETED_BREP` resource closure under its retained AP203/AP214 declaration,
  and binds caller-explicit `FILE_NAME` metadata. Coordinates use deterministic
  binary64 round-trip decimal spelling. Before publication the emitted bytes
  must pass the bounded syntax parser and native semantic decoder, with every
  coordinate bit and triangle index equal to the sealed source. The receipt
  retains caller metadata, both limit sets, the complete source decoder
  receipt, output syntax receipt, and output decoder receipt. This is a
  resource-level replay rung, not a NURBS re-fit or supplier-ready AP artifact.

## Invariants

1. **Round-trip fidelity per format**: OBJ and PLY re-import bitwise-
   identical f64 positions; STL agrees to f32 precision (documented
   lossy: positions only, welded by exact coordinate match, normals
   recomputed).
2. **No import is trusted without promotion**: the census runs on every
   import; promotion refuses while blocking defects remain, and both
   outcomes emit ledger-ready receipt JSON with the source hash, parser
   version, defect census, and trust status.
3. **Hostile input never panics**: 13.5k byte-mutants, all truncation
   prefixes, and pure junk across all three formats produce structured
   results (CI-checked fuzz lane).
4. **Deterministic exports**: identical soups produce identical bytes
   (fixed ZIP timestamps, fixed chunk layout).
5. **Schema errors teach**: row + column + offender + expectation.
6. **Part-21 graph integrity**: instance IDs are positive and unique;
   forward references are permitted but every reference must resolve by
   end of DATA; mandatory header records occur exactly once and in the
   supported order.
7. **Part-21 resource bounds**: input/output bytes, tokens, instances,
   values, nesting, encoded strings, number tokens, identifiers,
   complex-instance components, and schema-count each have an explicit
   nonzero cap. Recursive nesting also has an implementation hard ceiling
   independent of caller configuration. Cap violations are `ResourceBound`,
   not partial parses.
8. **Canonical syntax, not canonical CAD**: Part-21 output has fixed
   whitespace/keyword casing and numeric-ID instance order. It never
   reorders parameters or complex components, whose schema meaning is
   unknown at this layer. Numeric lexical spelling remains identity-bearing:
   this is layout canonicalization, not schema-aware numeric normalization.
9. **No topology laundering at the STEP handoff**: repair is always invoked
   with a zero hole-fill budget. Residual leaks, non-manifoldness, orientation
   conflicts, vertex-link failures, and non-outward aggregate orientation
   refuse publication; localized diagnostics are bounded to 256 records and
   state when truncated.
10. **No deviation laundering**: declared deviation must be a finite, ordered,
    non-negative `Exact`, `Enclosure`, or `Estimate` band. It remains separate
    in the receipt, and its upper bound is added with outward rounding to the
    mesh-to-SDF upper bound. The combined result is always `Estimate`, never a
    stronger authority grade.
11. **Every semantic input moves provenance**: exact soup position bits and
    triangle indices are FNV-fingerprinted before and after repair. Output
    provenance also binds the Part-21 source/layout fingerprints, adapter
    identity, shared length-unit ID, target-spacing bits, complete deviation
    certificate, deterministic execution mode, repair result, and underlying
    mesh-to-SDF provenance. These 64-bit fingerprints are replay aids, not
    collision-resistant authority.
12. **STEP tessellation preprocessing is separately bounded**: one million
    vertices, one million triangles, a conservative 512 MiB auxiliary-memory
    admission estimate, and at most 256 retained localized defects. The
    receipt records these limits, crate versions, STEP-import semantics label,
    and tessellation-fingerprint domain.
13. **Native faceted traversal is explicit and closed**: callers select a
    positive `FACETED_BREP` root. Only its fixed-depth pinned entity closure is
    interpreted; every reachable instance must be simple, exact-arity, and the
    expected entity type. `FACE_SURFACE` geometry is restricted to `PLANE` with
    `AXIS2_PLACEMENT_3D`, a 3-D `CARTESIAN_POINT` location, and omitted or 3-D
    finite nonzero `DIRECTION` values. Unknown unrelated instances remain
    outside the claim.
14. **No implicit triangulation or welding**: every `POLY_LOOP` has exactly
    three unique point references. Shared point IDs become shared soup vertices;
    distinct IDs with equal coordinates remain distinct. Holes, extra bounds,
    non-triangular loops, reused bounds/loops, and complex reachable instances
    refuse instead of being guessed or repaired by the decoder.
15. **Canonical semantic materialization**: shell face references, point
    positions, and triangles are emitted in numeric instance-ID order. Shell
    `SET` permutation therefore preserves the soup and semantic fingerprint;
    source spelling remains separately fingerprinted. `.T.` preserves the
    `POLY_LOOP` order and `.F.` reverses it. Plane, placement, location,
    direction, and `same_sense` semantics also move the semantic fingerprint
    even when two closures materialize the same soup.
16. **Schema labels gate but do not certify**: the decoder admits one exact
    declaration, either `CONFIG_CONTROL_DESIGN` or `AUTOMOTIVE_DESIGN`. The
    declaration is recorded as provenance, never promoted into EXPRESS or
    application-protocol authority. Finite coordinate conversion and accepted
    point-to-plane residuals carry conservative `Estimate` bands. Plane-backed
    faces refuse non-coplanar vertices, numerically degenerate triangles,
    direction drift, and winding inconsistent with `same_sense`. The existing
    zero-hole-fill handoff remains the sole owner of its bounded edge-use, local
    vertex-link, and aggregate-orientation admission; neither receipt claims
    global shell connectedness, component nesting, or self-intersection
    certification.
17. **Decoder memory admission is portable and explicit**: the auxiliary cap
    covers checked logical element payloads for every simultaneously live
    decoder vector. Platform allocator rounding and container headers are not
    misrepresented as measured bytes; `try_reserve_exact` failure still returns
    a structured resource refusal.
18. **Faceted re-emission is sealed and self-replaying**: version 1 accepts only
    an existing native decoder result whose faces are all bare `FACE`; it
    refuses plane-backed inputs rather than discarding surface semantics.
    Dense point/loop/bound/face/shell/root IDs are checked for overflow, input
    coordinates and indices are re-admitted, and output cannot publish until
    bounded write-parse-decode replay returns bit-identical positions and exact
    triangle order. Source and output semantic fingerprints remain distinct
    because instance-ID normalization is an explicit transformation.

## Error model

`IoError`: `Malformed { at, what }`, `Unsupported`, `ResourceBound`,
`Schema { row, column, what }`. `PromotionRefusal` carries blocking
defects + fixes + the refused receipt. The STEP syntax kernel uses
`Malformed` for grammar/graph failures, `Unsupported` for staged encoded
characters and binary literals, and `ResourceBound` for every declared
limit. `StepImportRefusal` separates raw admission, localized mesh integrity,
preprocessing resource admission, SDF build/cancellation, and evidence-
composition failures; each variant keeps the source fingerprint and later-
stage variants keep repair receipts. `StepFacetedRefusal` separately reports
schema-gate, root-reachable entity, decoder-resource, and cancellation refusals
with the source fingerprint plus exact instance relationship or decoder stage.
`StepFacetedImportRefusal` preserves whether refusal happened during native
materialization or the downstream topology/SDF handoff; a downstream refusal
retains the successful decoder receipt and selected-root provenance.
`StepFacetedExportRefusal` separately names source admission, checked
construction/allocation bounds, cancellation, nested syntax write/parse,
nested semantic replay, and first exact-geometry replay mismatch.

## Determinism class

**D0**: fixed parse/emit orders, deterministic welds/topology sorts, no ambient
state. Native faceted decoding sorts schema-defined `SET` members and materializes
points/faces by numeric instance ID. The STEP tessellation handoff rejects
`ExecMode::Fast`; its receipt and provenance explicitly bind deterministic mode.
Strict faceted re-emission uses dense IDs and source decoder order, fixed header
fields plus exact caller metadata, and deterministic round-trip decimal
spelling. Identical sealed source, metadata, and limits produce identical bytes
and nested receipts on one target.

## Cancellation behavior

Legacy mesh/catalog parsers are single-pass and element-capped. The STEP
kernel is deliberately multi-pass (parse, shape/graph validation,
canonical-layout serialization) and cap-bounded, but it has no `Cx` and
makes no cancellation-latency claim. Native faceted decoding polls at entry,
publication, duplicate/deduplication scans, and every 4096 indexed instances,
faces, and points. Sorting is a deterministic sequence of at-most-4096-element
local sorts followed by a cancellable k-way merge, so no million-record
standard-library sort becomes an unpolled region. The identified tessellation
handoff polls `Cx`
at entry, around cap-bounded library calls, and every 4096 records in its owned
validation, fingerprint, vertex-compaction, edge-localization, and vertex-link
passes before forwarding the same `Cx` to mesh-to-SDF sampling. Cancellation is
reported as `StepImportRefusal::SdfBuild`. The existing `repair`, topology-sort,
and `MeshChart` construction calls have no internal poll, so this subset makes
no sub-call latency claim for those separately bounded stages.
Strict faceted re-emission polls at entry/publication, every 4096 points and
faces during owned graph construction, and around the bounded syntax writer,
parser, and semantic decoder. The syntax subcalls have no internal `Cx`, so the
re-emitter makes no sub-call latency claim inside those separately capped
operations.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (canonical fs-obs `ConformanceCase` aggregate outcomes,
suite `fs-io/conformance`): io-001 STL/OBJ/PLY round trips (exact where the
format allows) +
deterministic bytes + ASCII STL fixture; io-002 the defect zoo
(duplicate/degenerate/hole/unreferenced) censused, repaired, promoted
with receipts — and an over-budget hole REFUSED with actionable fixes
and a refused receipt; io-003 13.5k mutants + truncations + junk with
zero panics; io-004 PLY face-list integer validation; io-005 AISC-flavored
CSV + JSON catalogs, quoting, and the teaching-error battery; io-006 3MF ZIP
structure (EOCD, entry count, model XML), GLB chunk accounting, VTK section
counts. Reaching an aggregate outcome means its preceding checks passed;
pre-verdict assertions and parser `expect` failures remain ordinary Rust test
diagnostics and therefore do not emit a failed aggregate record. `io-003`
records its exact input seed (`0x10_0003`) for mutation-stream replay; the other
five outcomes use deterministic seed zero. The suite has no concurrent
aggregate case, so these records make no scheduler-replay claim. Existing
promotion-receipt and fuzz-measurement data use validated fs-obs `Custom`
companions, not canonical aggregate outcomes; the fuzz companion also retains
the mutation-stream seed.

`tests/step.rs` (G0/G3): forward-reference and complex-entity parsing,
canonical permutation-invariant DATA ordering, doubled-apostrophe string
round trip, AP-family hint/receipt binding, duplicate and dangling
reference refusal, malformed/truncated envelope/comment/string/value
refusal, mandatory-header shape checks, strict uppercase keywords, exact
typed-parameter arity, explicit resource/hard-depth-cap refusal, and
writer-side revalidation of caller-constructed invalid graphs.

`tests/step_import.rs` (G0/G3/G4): sealed source/adapter/error receipt
composition on a closed fixture; deterministic leak and non-manifold
localization; no-hole-fill repair behavior; hostile admission cases; retained
duplicate/degenerate/unreferenced repair receipts; repair-exhausted audit
retention; closed disconnected vertex-link refusal; pre-requested cancellation;
outward-rounding overflow refusal; and fast-mode refusal. Differential fixtures
require changed soup bits or deviation claims to move output provenance.

`tests/step_faceted.rs` (G0/G3/G4): unsorted tetrahedron closure and canonical
soup materialization; bound-orientation reversal; shell-`SET` permutation
invariance; exact supported and refused schema declarations; plane-backed
`FACE_SURFACE` equivalence, default-axis handling, `same_sense` reversal, and
plane-provenance binding; non-coplanar, misoriented, parallel-direction,
short-direction, non-triangular, duplicate-point, non-finite-coordinate,
vertex-cap, and auxiliary-memory refusals plus the independent triangle cap;
pre-requested cancellation; and proof that the native bridge reaches the
existing topology quarantine rather than laundering an open shell.

`tests/step_faceted_export.rs` (G0/G3/G4/G5): exact finite-coordinate and
triangle replay through dense-ID export/write/parse/decode; negative-zero and
nontrivial binary64 decimal preservation; nested receipt and AP203/AP214
profile binding; byte-identical repeated export; canonical apostrophe escaping;
syntax-instance and semantic-vertex bound refusals; invalid metadata refusal;
plane-semantic-loss refusal; and pre-requested cancellation.

## PLY element order (bead wqd.25.1)

Element order is the header's to define: faces may legally precede
vertices. Parsing collects triangulated faces as pending records
(structural checks and the 1024-item list cap and triangle cap apply
immediately); index RANGE validation runs once, after every element is
consumed, against the final vertex count — with the exact offending
triangle ordinal in the diagnostic. Vertex-first and face-first files
import identically in both ASCII and binary (conformance-tested).

## No-claim boundaries

- **Full native STEP CAD semantics remain STAGED**: the syntax kernel does not
  load an EXPRESS schema or authorize AP203/AP214 conformance.
  `StepProfileHint` is label recognition only. The strict faceted decoder derives
  a triangle soup from one bounded resource closure, but does not
  interpret products, assemblies, shape-representation linkage, units/context,
  AP global rules, non-planar surfaces, voids, or general B-rep topology.
  Plane-backed `FACE_SURFACE` support proves only the pinned `PLANE`, placement,
  direction, coplanarity, and winding relationships. It does not parse or
  certify `FACETED_BREP_SHAPE_REPRESENTATION`, product/context correspondence,
  or the application protocol's global rules. External handoff adapters remain
  responsible for any semantics outside this native closure.
- **STEP-derived SDF authority is Estimate only**: the handoff does not certify
  component nesting, self-intersection freedom, generalized-winding sign, or
  full semantic correspondence between arbitrary Part-21 records and a
  tessellation. The native decoder claims correspondence only for its selected
  admitted closure and records decimal-to-f64 conversion plus accepted
  plane-consistency residual as an estimate.
  The strict re-emitter writes only its already-decoded bare triangular
  resource closure; it does not fit NURBS, establish an AP-conformant
  topological solid/product, compose a source SDF deviation bound, or establish
  manufacturing predicates.
- **Part-21 encoded characters and binary literals are refused** in this
  first subset. Source bytes must be ASCII; encoded-character directives
  and binary payloads need their own bounded conformance fixtures before
  admission.
- **Keywords/enumerations are strict uppercase Part-21 tokens**. Schema
  declaration admission may compare ASCII case-insensitively only because
  it operates on string payloads, not grammar keywords.
- **IGES and IFC are STAGED, not promised**; their quarantine paths have
  not shipped.
- **OBJ vt/vn and materials are dropped** (documented lossy subset);
  PLY color/normal properties are skipped, not preserved.
- **PLY binary_big_endian is refused** (structured `Unsupported`).
- **3MF/GLB are WRITE-ONLY** (import of container formats is follow-up);
  the 3MF is the minimal core-spec package, no extensions.
- **VTK export is legacy-ASCII**, one optional scalar field; XML VTK and
  vector/tensor fields land with fs-viz interop needs.
- **The census's manifoldness check is combinatorial** (edge counts +
  half-edge build); geometric self-intersection certification belongs to
  the validity-certificates machinery (wqd.23) and can be layered onto
  promotion by callers.
- **Receipts hash with FNV-1a**; HELM upgrades to the BLAKE3-class
  content address when writing the `imports` row (same field, stated in
  the receipt schema).
