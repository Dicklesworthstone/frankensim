# fs-rep-voxel — CONTRACT

Occupancy/multi-material voxel charts, point clouds with normals, and
explicit lattice/strut graphs (plan §7.2) — the discrete-representation
side of MORPH, sharing the sparse VDB substrate with fs-rep-sdf.

Ambition tags: fields/morphology/DT/clouds/lattices [S].

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, fs-rep-sdf (`VdbGrid`),
fs-geom (`Chart`), fs-exec (`Cx`), fs-evidence, fs-math, and the
constellation's fnx-classes/fnx-runtime (FrankenNetworkx). Consumers:
validity/topology certificates (wqd.23), lattice/infill homogenization
(7tv.14), SIMP consumers.

## Public types and semantics

- `field`: `OccupancyField` (bool active set + voxel size + origin) with
  active-set booleans (union/intersect/subtract) and morphology
  (dilate/erode/open/close, 6-connected — the manufacturing-constraint
  primitives); `MaterialField` (u16 ids, 0 reserved for empty, per-id
  occupancy extraction); `DensityField` (SIMP fractions, out-of-range
  REFUSED not clamped). Manual `Clone` rebuilds through the active set
  (the substrate carries no derives).
- `dt`: `euclidean_dt` — EXACT Euclidean distance transform
  (Felzenszwalb–Huttenlocher separable lower envelopes) over the active
  bounding box; squared distances are integer-valued in voxel units and
  conformance checks EQUALITY against the O(n²) reference. Infinite
  (seedless) lines are handled by envelope construction over finite
  parabolas only.
- `cloud`: `PointCloud` with grid-hash radius/kNN queries (brute-force
  verified, deterministic tie order), PCA normal estimation (smallest
  covariance eigenvector via cyclic Jacobi), and orientation propagation
  that SURVIVES DISCONNECTED k-NN graphs (scan-line clouds): BFS per
  component with nearest-visited-anchor alignment on restart.
- `lattice`: `LatticeGraph` (nodes with junction blend radii, struts
  with cross-section radii). Coincident nodes, zero-length struts,
  missing indices, and non-positive radii are structured refusals.
  `to_fnx`/`from_fnx` round-trip through FrankenNetworkx with attributes
  preserved exactly. `sdf(p)` realizes the solid as a smooth-min of
  capsule SDFs — a continuous level set, watertight by construction;
  per-strut realization receipts are JSON log lines.
- `chart`: `OccupancyChart` implements fs-geom `Chart`: inside/outside
  from occupancy, distance magnitude from the exact DT on both
  polarities (complement DT inside), an exact active-center scan
  fallback outside the DT box, and an HONEST error certificate — an
  Enclosure of ± half a voxel diagonal, never "exact".

## Invariants

1. **Morphology matches brute force**; algebra laws hold on fixtures
   (opening ⊆ id removes sub-resolution spurs, closing ⊇ id,
   `(A∪B)\B ⊆ A`, `A∩B ⊆ B`).
2. **The DT is exact**, not approximate: equality with the O(n²)
   reference on scattered+slab fixtures; 1-Lipschitz in the voxel metric.
3. **Cloud queries match brute force** (radius and kNN); sphere normals
   are >97% outward-aligned after propagation — including on ring-sampled
   (kNN-disconnected) clouds.
4. **Graph round-trip is lossless**: node positions/radii bitwise, strut
   radii exact, undirected edges deduplicated deterministically.
5. **Realization behaves like a closed solid**: strut midpoints inside,
   far field outside, straight probes cross the boundary an even number
   of times.
6. **Resolution error is declared**: every chart sample carries the
   ±half-diagonal enclosure.

## Error model

`VoxelError`: `Parameters`, `Lattice` (offending element named),
`Cloud`, `Graph`. Nothing silently clamps or skips.

## Determinism class

**D0 on-target**: BTreeMap substrate ⇒ sorted iteration everywhere;
fixed DT pass order; deterministic query tie-breaks (distance, index);
BFS restarts at the lowest unvisited index.

## Cancellation behavior

Operations are bounded by active-set/box size; the chart polls nothing
itself (single-voxel lookups) — P7 by boundedness. Long DT boxes are the
caller's budget decision (see no-claims).

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-rep-voxel/conformance`):
rv-001 morphology vs brute force + algebra + spur removal; rv-002 DT
exactness vs O(n²) + 1-Lipschitz; rv-003 radius/kNN vs brute force +
sphere normals (the ring fixture DELIBERATELY breaks kNN connectivity —
it caught the propagation gap during development); rv-004 fnx
round-trip + degenerate refusals + level-set probe parity + realization
receipts; rv-005 the chart contract (inside/outside, DT-backed distance
near analytic, declared resolution error, out-of-box fallback).

## No-claim boundaries

- **The DT densifies the ACTIVE BOUNDING BOX**: memory is box volume,
  not active count — far-flung sparse sets pay for their box (tiled
  narrow-band DT is follow-up work with fs-rep-sdf's band machinery).
- **Out-of-box chart queries scan active centers**: exact but O(active);
  bulk far-field queries should go through the SDF representations.
- **Morphology is 6-connected**: 26-connected structuring elements and
  true Euclidean openings (via the DT) are follow-ups.
- **Normal orientation is heuristic** (BFS + nearest-anchor): globally
  correct on orientable, reasonably sampled surfaces; adversarial
  topology can still fool it — fitting beads consume the normals WITH
  their cloud provenance.
- **Realization emits an SDF**, not a mesh: meshing goes through the
  Rep Router to fs-rep-mesh; junction blending is smooth-min (bounded
  overshoot), not exact fillet geometry.
- **No FrankenNumpy views yet** (the fnp-* bridge lands with the SoA
  bead's consumers).
