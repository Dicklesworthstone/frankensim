# fs-mesh — CONTRACT

## Purpose and layer

L2 (MORPH). Body-fitted tet meshing (plan §7.5) for when a mesh is
WANTED — final verification, shells, export — remembering CutFEM-on-SDF
exists precisely so meshing stays optional inside optimization loops.
v1 is the Delaunay KERNEL — BRIO-ordered incremental Bowyer–Watson on
fs-ivl's exact predicates, with ghost tets carrying the hull, plus
radius-edge quality refinement — and SURFACE REMESHING: the
Botsch–Kobbelt split/collapse/flip/smooth loop measured in a Riemannian
metric (isotropic = identity metric), chart-projected, feature-locked.
Everything the crate claims about its output, it re-checks (`audit`,
half-edge round-trips, closed-manifold audits).

## Public types and semantics

- `delaunay(&[Point3], cx) -> Result<Tetrahedralization, MeshError>`:
  BRIO order (deterministic LCG shuffle → doubling rounds → Morton sort
  within rounds), visibility-walk location with locality hints,
  Bowyer–Watson cavity insertion. Bitwise-duplicate points are skipped
  WITH a stats receipt. Conflict rules are exact and canonical: real
  tets by strict `insphere` (cospherical `Zero` = NOT in conflict — the
  deterministic weak-Delaunay choice); ghost tets by `orient3d`, with
  exactly-coplanar cases delegated to an in-plane exact `incircle`
  (the halfspace-closed-by-the-disk rule). SoS appears ONLY in the
  walk; conflict regions are SoS-free, which is what makes the cavity
  star-shape argument (real boundary facets strictly visible) hold —
  the growth-repair path is a counted safety net, 0 on the whole zoo.
- `Tetrahedralization`: `tets()` (positively oriented, canonically
  ordered), `points()`, `hull()` (outward-oriented `Soup`),
  `complex()` (fs-rep-mesh `TetComplex`, δδ = 0), `stats()`,
  `audit(full_insphere)`.
- `audit`: exact self-audit — positive orientation, mutual adjacency,
  LOCAL Delaunay on every internal facet (the Delaunay lemma lifts
  local to global), Euler characteristic = 1, hull closed, hull
  EXACTLY convex; `full_insphere` adds the O(n·t) global
  empty-circumsphere check for fixture-scale belt-and-braces.
- `refine(&mut Tetrahedralization, RefineOptions, cx)`: worst-first
  radius-edge refinement by circumcenter insertion through the same
  kernel; offenders whose circumcenters escape the hull are SKIPPED
  AND COUNTED (`unrefinable_remaining`) — the honest v1 policy until
  constrained boundary handling lands. Steiner points append after
  `steiner_from`.
- `GHOST`: the at-infinity sentinel (slot 3 of hull tets), exposed for
  audit tooling.
- `remesh(&Soup, Option<&dyn Chart>, &dyn MetricField, RemeshOptions,
  cx)`: unit-METRIC-length remeshing — split above 4/3, collapse below
  4/5 (link condition, no-new-long-edge and normal-flip guards), flips
  toward valence 6 (fold-over guarded), Jacobi tangential smoothing —
  with Newton projection onto the chart for every placed or smoothed
  vertex. Dihedral creases, boundaries, and non-manifold fins are
  LOCKED (never flipped/collapsed; endpoints never smooth); split
  midpoints always project, which is a no-op on straight creases.
  Passes are FUNCTIONAL (connectivity rebuilt in `BTreeMap`s, ops in
  canonical order): auditable and P2-deterministic over raw
  throughput, until the perf lane profiles it.
- `MetricField` / `UniformMetric`: the SPD tensor input — isotropic
  remeshing IS `UniformMetric`; anisotropic fields (ultimately FLUX's
  DWR error metric) reuse the identical op set.

## Invariants

1. On general-position clouds the FULL exact audit is clean: global
   empty circumsphere, local Delaunay, orientation, adjacency, Euler,
   exact hull convexity (tmesh-001).
2. The degeneracy battery completes CORRECTLY on exact predicates:
   integer grids (massively cospherical/coplanar), exactly cospherical
   shells, collinear runs — all clean under the full audit; bitwise
   duplicates are skipped with receipts; all-coplanar input refuses
   with a teaching error (tmesh-002).
3. Determinism (P2/G5): identical input gives BITWISE-identical
   meshes; relabeled input gives the identical geometric tet set;
   dyadic translations preserve connectivity exactly with
   exactly-shifted coordinates (G3) (tmesh-003).
4. The hull soup is closed, 2-manifold, outward-oriented (winding +1
   inside), and the oriented complex satisfies δδ = 0 exactly
   (tmesh-004).
5. Refinement leaves NO interior-refinable offender above the
   radius-edge bound, keeps the full exact audit clean through every
   Steiner insertion, and is deterministic; hull-escaping survivors
   are counted, not hidden (tmesh-005).
6. Scale: 10k-point clouds build with clean O(t) audits and BRIO
   locality (order-10 walk steps per insertion, no exhaustive
   fallbacks) (tmesh-006).
7. Isotropic remeshing concentrates edges at unit metric length (>85%
   in [0.7, 1.4]), keeps every vertex ON the chart to fp precision,
   bounds centroid sag by the chord sagitta, stays closed/manifold/
   outward, is BITWISE deterministic, and is translation-equivariant
   in QUALITY PROFILE (threshold-driven ops legitimately flip borderline
   decisions under shifted fp arithmetic — the honest G3 statement)
   (tmesh-007).
8. Randomized remesh storms keep half-edge invariants, closed-manifold
   status, and Euler = 2 after EVERY round (tmesh-008).
9. Remeshing a cube keeps all 8 corners BITWISE, keeps every
   crease-grade output edge on a cube edge line, and stays on the box
   chart (tmesh-009).
10. The boundary-layer metric is realized: metric-unit conformity,
    physically stretched equator-aligned layer elements, and a MEASURED
    interpolation-residual win over isotropic at comparable element
    count — the adaptivity loop's value, demonstrated (tmesh-010).

## Error model

`MeshError` teaching errors: `TooFewPoints`, `DegenerateInput` (exact
all-coplanar detection, says to triangulate in 2D instead),
`Cancelled`. Kernel internals hold invariants by construction (no
flat tets: every created tet's apex is strictly visible); the audit
exists so any regression is LOUD rather than silently non-Delaunay.

## Determinism class

Fully deterministic and sequential in v1: fixed-seed BRIO shuffle,
exact predicate signs, canonical conflict rules, index-ordered
tie-breaks, `BTreeMap`/`BTreeSet` only. Identical input bytes →
identical output bytes (tmesh-003 is the trip-wire). The bead's
"same mesh at any thread count" criterion is met NON-trivially by
`delaunay_colored` (uee3 item 4): read-parallel conflict regions
(cavity + growth repair + one-ring, mirroring the insert transaction)
across scoped threads, FLIP-SAFE coloring (k = 1 + the largest
overlapping color — same-color members pairwise disjoint AND every
order-flipped cross-color pair disjoint, so cospherical TIE groups
keep their original order), canonical application. Thread count can
change only the wall clock; tmesh-013 gates raw thread-count
invariance, canonical kernel merge on general-position AND degenerate
fixtures, exact audits, adversarial within-color commutativity
(reversed application), and the width ledger. Two designs were
REJECTED on measurement: first-fit coloring (flipped tied pairs —
diverged on the 6×6×6 grid) and stop-at-first-clash prefix batching
(raw-order-preserving but BRIO locality collapsed width to ~3). Batch
width is STRUCTURAL (~6 at window 256: Hilbert-ordered windows form
mutually-overlapping chains, one color per chain element); strided
sampling would widen batches but reorders ties — rejected; the read
phase parallelizes independently of width.

## Cancellation behavior

`delaunay` polls `cx.checkpoint()` every 256 insertions; `refine`
polls per round; `remesh` polls per iteration. Cancellation returns `MeshError::Cancelled` between
insertions (request → drain → finalize; no torn mesh states escape
since the error consumes the builder).

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases tmesh-001..tmesh-010 — JSON-line
verdicts, seeded LCG randomness, fs-obs Custom events for kernel
stats, refinement stats, and scale/walk-locality stats. Any
reimplementation must pass the suite unchanged.

### Addendum (bead uee3, partial): policy floor, hull-split evidence, exudation

- `RefineOptions` gains `min_edge_factor` (the SMALL-INPUT-ANGLE
  POLICY: a minimum-new-edge floor from the input's closest-pair
  spacing; insertions below it YIELD and are counted as
  `protected_by_policy`) and `split_hull_facets` (default OFF): hull-facet
  splitting now runs under DIAMETRAL ENCROACHMENT PROTECTION
  (`facet_diametral_ball`) — the classical Ruppert rule, split a facet IFF a
  circumcenter lands in its minimum-enclosing sphere; an escaping circumcenter
  encroaching nothing is skipped (an unfixable boundary sliver). The in-plane
  split point is blended strictly into the facet interior (a point exactly on a
  hull edge is collinear-degenerate: the audit went red before the blend). It is
  exact-audit-clean and deterministic and MEASURABLY shrinks the convex-hull
  regression (~2.8e18 → ~3.5e17, ~8×, gated in tmesh-011 at `worst_after < 1e18`),
  but does NOT eliminate it: residual slivers come from near-boundary INTERIOR
  vertices, so true full-Ruppert quality stays coupled to constrained boundary
  recovery, exactly as the classical termination theory requires.
- `exude` / `ExudeOptions` / `ExudeStats`: sliver removal by
  deterministic Steiner PERTURBATION — offending Steiner vertices
  nudged by seeded deterministic offsets, full rebuild through the
  exact kernel, rounds kept only when the sliver census strictly
  drops AND the exact audit stays clean; input points are never
  touched (bitwise-checked). The weighted-Delaunay exudation pump
  needs a weighted exact predicate — recorded no-claim below.

## No-claim boundaries

- Weighted exact insphere predicate (the Edelsbrunner weight-pump
  exudation variant; the perturbation flavor ships).
- INTERIOR FACET recovery now ships in CONFORMING form for CONVEX
  planar facets (`recover_facets`, tmesh-015): batched longest-edge
  midpoint bisection of the fan triangulation (one-split-per-round was
  MEASURED to starve at the rounds cap; batching finished the fixture
  in 7 rounds), twin adoption via the shared coordinate-bits index,
  a facet correspondence table re-verified against the finished mesh,
  and honest starved-budget counters. Remaining no-claims:
  NON-CONVEX facets (fan triangulation needs convexity) and
  general-position planes (f64 midpoints stay exactly coplanar only
  when the plane is axis-aligned — the battery gates the bitwise case
  and the residual is measured, not assumed). Full-Ruppert QUALITY:
  the diametral encroachment machinery cut the hull-split regression
  ~8× (tmesh-011); the residual is coupled to boundary-layer
  refinement — still successor scope.

- SEGMENT recovery now ships in CONFORMING form (`recover_segments`,
  tmesh-014): recursive midpoint Steiner bisection with twin-vertex
  ADOPTION at shared midpoints (the four body diagonals of a box meet
  at its center — abandoning bitwise-duplicate midpoints was measured
  to strand 3 of 4 segments before adoption landed), a boundary
  CORRESPONDENCE table mapping every sub-edge to its parent segment
  (built by construction, re-verified against the finished mesh), and
  honest `unrecovered` counters at depth/budget caps. Convex
  hull-facet conformity is gated test-side; interior/non-convex FACET
  recovery (true constrained DT) remains the successor.
- Refinement is radius-edge with the minimum-new-edge policy floor;
  full local-feature-size Ruppert guarantees remain successor scope
  (sliver exudation ships in `exude`).
- Parallel domain coloring SHIPS (`delaunay_colored`, tmesh-013) —
  see Determinism; v1 is sequential.
- The 10⁷-point perf lane RAN (2026-07-09, ts1: Threadripper PRO
  5975WX, Linux x86_64, release): 10⁷ points in 100.0 s = 100,034
  points/s, 67.6M tets, throughput near-FLAT across the ladder
  (10⁴: 116k, 10⁵: 115k, 10⁶: 105k, 10⁷: 100k pts/s — BRIO locality
  holding over three decades), exact structural audit clean at every
  rung (full insphere certification at 10⁴). Ledger rows in
  `tests/perf_lane.rs` output and bead uee3. tmesh-006 pins 10⁴-scale
  behavior in CI; the nightly perf-CI cadence belongs to fz2.4.
- Remeshing no-claims: curved creases round under midpoint projection
  (straight creases are exact); boundary loops are locked, not
  remeshed; metric gradation control, log-Euclidean metric
  interpolation/intersection, and DWR-supplied discrete metric fields
  join with FLUX's estimator bead; the functional-pass architecture
  trades throughput for auditability until the perf lane says
  otherwise.
- `orient3d_sos` is a projection cascade, not the full 3D
  Edelsbrunner–Mücke ladder (fs-ivl's documented no-claim); it is used
  only for walk routing here, never for conflict decisions.
