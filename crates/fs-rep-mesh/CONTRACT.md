# CONTRACT: fs-rep-mesh

## Purpose and layer
Mesh charts (plan §7.2): half-edge surfaces with an edit core, oriented
tet complexes with exact-arithmetic incidence (the pre-FEEC δδ = 0
sanity), generalized-winding-number robustness for polygon soup, a
BVH-accelerated signed-distance chart with watertight raycasts, and the
repair suite with structured receipts. Layer: L2. Depends on fs-geom,
fs-exec, fs-evidence, fs-alloc, fs-obs.

## Public types and semantics
- `HalfEdgeMesh` — `from_triangles` (teaching refusals: non-manifold
  edges/pinched boundary vertices — the SOUP path is the repair suite),
  boundary loops as NO_FACE half-edges, `flip_edge` (Botsch relink;
  boundary/duplicate-creating flips are non-events), `check_invariants`
  (twin involution, 3-cycles, twin-origin coherence — the property
  battery's oracle), `euler_characteristic`.
- `TetComplex` — canonical sorted edge/face tables (BTreeMap order,
  deterministic) and signed `Incidence` operators d0/d1/d2 over INTEGER
  cochains; `HexComplex` is storage-only (no-claim below).
- `Soup` / `triangle_winding` (van Oosterom–Strackee) / `winding_exact` /
  `WindingOctree` — the Barill-style dipole hierarchy (β accuracy knob,
  area-weighted normal moments, exact leaves); octree stores indices
  only, callers pass the soup (no self-referential borrows).
- `point_triangle_distance` (Ericson region test, exact),
  `ray_triangle_watertight` (Woop permutation+shear; shared edges NEVER
  leak; exact-edge hits may double-count — documented for parity users),
  `Bvh` (median-split, deterministic tie-break by index; closest-point
  branch-and-bound + nearest-hit raycast).
- `MeshChart` — sd = BVH closest distance × winding sign; Lipschitz = 1
  is a RIGOROUS claim (distance-to-set); declared error covers fp slack
  only; `raycast` watertight.
- `repair(soup, max_hole_edges)` — dedupe → degenerate removal →
  orientation unification (flood fill + centroid winding vote) →
  fan-fill of small boundary loops, each action a `RepairReceipt`
  (defect/location/action, the fs-io quarantine format).
- `shapes` (PUBLIC fixture vocabulary): `cube`, `icosphere`,
  deterministic `corrupt` (dups/degens/flips/hole).

## Invariants
1. Half-edge invariants survive 2k random flips with the Euler
   characteristic preserved (rmesh-001).
2. `point_triangle_distance` never exceeds any sampled surface distance
   and matches 1830-point brute force within sampling gap; `MeshChart`
   tracks the analytic sphere within the icosphere chord band, satisfies
   inside ⇔ sd < 0, and honors its 1-Lipschitz claim (rmesh-002).
3. Winding classification is > 99% correct away from defects on the
   nightmare corpus (duplicates + degenerates + flipped patch + punched
   hole) — robustness on BROKEN input is the point (rmesh-003).
4. The dipole octree tracks exact winding within a measured, ledgered
   error off-surface (rmesh-004; β = 2).
5. The repair pipeline heals the corrupted corpus back to manifold,
   restores the face count and center winding, and emits receipts
   covering every defect class (rmesh-005).
6. d1∘d0 = 0 and d2∘d1 = 0 EXACTLY (integer cochains) on fixture
   complexes; axis rays through shared cube edges never leak; chart
   raycasts hit at analytic parameters (rmesh-006).

## Error model
Structured teaching errors (`MeshBuildError`); total functions elsewhere
(degenerate triangles yield well-defined distances; empty-soup handling
is the caller's constructor discipline). No panics across the boundary.

## Determinism class
Deterministic: BTreeMap orders, index-tie-broken BVH sorts, seeded
batteries; no clocks, no addresses in results.

## Cancellation behavior
`MeshChart::eval` is bounded per query (BVH descent); batch consumers
poll between queries per the fs-geom discipline. Repair/build are
bounded preprocessing passes (cancellation hooks join with the fs-io
quarantine bead where soups get large).

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. `[S]` solid-tier.

## Conformance tests
tests/conformance.rs, cases rmesh-001..rmesh-006 (JSON-line verdicts;
seeded cases carry seeds) covering invariants 1–6 with fs-obs-validated
evidence events (dipole error, repair receipts).

## No-claim boundaries
- Self-intersection FLAGGING is deferred to the validity-certificates
  bead (wqd.7), which owns certified broad/narrow phases via fs-ivl.
- `HexComplex` is storage only; its incidence operators land with
  fs-feec's tensor-product families.
- Curvature via discrete operators (cotan/normal-cycle) is deferred to
  its first consumer, with convergence-class documentation there.
- NO throughput claims (million-triangle dipole performance is the perf
  harness's); the octree's asymptotics are structural.
- The winding sign convention (+1 inside outward-oriented closed
  surfaces) is pinned by fixtures, not by a formal proof of the flux
  form; exact-edge raycast double-counting is documented.
- Attribute channels (normals/UVs/per-face data) beyond positions are
  deferred until a consumer (LUMEN materials) defines their semantics.
