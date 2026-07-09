# fs-io — CONTRACT

Import/export with QUARANTINE (plan patch Rev J): the world boundary.
Dirty geometry comes in, useful artifacts go out — and no imported
artifact becomes a trusted value without a certification receipt.

Ambition tags: STL/OBJ/PLY + quarantine + catalogs + 3MF/GLB/VTK [S];
STEP/IGES/IFC explicitly STAGED (no-claim below).

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, fs-rep-mesh (repair +
half-edge validity), fs-evidence, fs-geom, fs-obs, fs-math. PNG/EXR
export is fs-img's (L5). Ledger `imports` rows are written HELM-side
from the receipt JSON this crate emits — L2 never calls L6. Consumers:
the P4 frame flagship (AISC catalogs), fs-fab.

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

## Error model

`IoError`: `Malformed { at, what }`, `Unsupported`, `ResourceBound`,
`Schema { row, column, what }`. `PromotionRefusal` carries blocking
defects + fixes + the refused receipt.

## Determinism class

**D0**: fixed parse/emit orders, BTreeMap welds, no ambient state.

## Cancellation behavior

Parsers are single-pass and element-capped; export size is input-bounded.
P7 by boundedness.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-io/conformance`):
io-001 STL/OBJ/PLY round trips (exact where the format allows) +
deterministic bytes + ASCII STL fixture; io-002 the defect zoo
(duplicate/degenerate/hole/unreferenced) censused, repaired, promoted
with receipts — and an over-budget hole REFUSED with actionable fixes
and a refused receipt; io-003 13.5k mutants + truncations + junk with
zero panics; io-004 PLY face-list integer validation; io-005 AISC-flavored
CSV + JSON catalogs, quoting, and the teaching-error battery; io-006 3MF ZIP
structure (EOCD, entry count, model XML), GLB chunk accounting, VTK section
counts.

## No-claim boundaries

- **STEP/IGES and IFC are STAGED, not promised** (per the bead text):
  no subset ships here; the quarantine pipeline is where they will land.
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
