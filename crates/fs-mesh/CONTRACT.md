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
The goal-oriented adaptivity seam additionally publishes bounded,
declaration-only accounting receipts over opaque retained L3/L6 identities;
it does not introduce an upward crate dependency or recreate the authority of
the estimator or Machine-IR lineage record it names.
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
- `AdaptivityReceipt::admit`: dependency-clean accounting seam for dynamic
  h/p and anisotropic mesh evolution. It binds a typed action and
  contact/wear/fracture/moving-mesh/goal-oriented trigger; explicit declared
  connectivity, physical-topology, and gradient-discontinuity effects; opaque
  source and target mesh-state plus lineage-record identities; state-bound
  before/after QoI-evidence identities; separate estimator and
  representation-conversion upper bounds; a retained remap-invariant identity
  (quantity, units, and balance convention), remap-evidence identity, signed
  balance defect, declared tolerance, and projection error; and a strict
  `Decreased`/`Unchanged`/`Increased` QoI-bound result. Effects are retained
  rather than inferred from a coarse action name. The only constructible
  authority is `DeclarationOnly`: `fs-mesh` validates complete accounting but
  cannot certify caller-supplied DWR, conversion, conservation, or lineage
  claims. Canonical JSON is available for an owning ledger to hash.
- `conservative_cell_remap(source_values, target_count, contributions,
  source_coverage_tolerance, balance_tolerance, cx)`: bounded sequential
  piecewise-constant remap for one cellwise EXTENSIVE scalar. Contributions
  are source-volume fractions and must arrive in strict `(source, target)`
  order, with every source row partitioning unity and every target receiving
  data. The kernel rejects duplicate/out-of-order pairs, gaps, invalid indices,
  zero/negative/non-finite/>1 fractions, non-finite values, arithmetic
  overflow, excessive local coverage defect, and excessive global
  target-minus-source balance before publishing any target vector. It polls
  cancellation through source, overlap, allocation-initialization, target
  coverage, and publication scans; static cell/contribution and 256 MiB
  auxiliary-storage ceilings precede work. The canonical report is explicitly
  `measured-f64`; `report.accounting(...)` requires the caller's own invariant,
  evidence, and projection-error declaration before feeding the receipt seam.

- `hexdom` module (plan §7.5, bead wqd.18; [F], behind
  `frontier-hexmesh`, OFF the critical path): hex-dominant meshing via
  octahedral frame fields. SH9 realized as a FIXED SPHERICAL SAMPLING
  of the degree-4 octahedral polynomial (a linear image of the SH9
  coefficients — exact Wigner-D machinery is the growth path); MBO =
  graph diffusion + deterministic seeded projection to the variety
  (energy decreases monotonically, boundary frames pinned); the
  24-element cube group drives matchings; SINGULARITIES are loop
  holonomies of matchings around lattice edges (winding, not local
  twist — a 45° isolated cell is NOT singular, measured);
  `extract_hex_dominant` routes frame-field / polycube-fallback /
  refusal by DOCUMENTED criteria, and refusals name IGA/CutFEM (the
  honest-alternatives doctrine); `accuracy_per_dof` reports both
  element classes whichever way it falls.

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
11. An adaptivity receipt names one declared QoI on both sides, retains exact
    source/target/evidence/lineage digest bytes, rejects non-finite or negative
    bound components and non-finite composed totals, and reports a strict
    QoI-bound trend without treating unchanged or regressed steps as
    successful. An error-free two-sum rounds outward only when the computed
    sum rounded down, so moving an identical error total between estimator and
    conversion-ledger components cannot change the trend. Signed zero is
    canonicalized for replay-stable JSON (G0/G3).
12. Before/after QoI snapshots must name the lineage source/target states.
    Declared effects refuse physical-topology change without connectivity,
    cannot suppress the gradient-discontinuity flag without a future evidence
    path, and must match the fixed semantics of h, p, and untangle actions.
13. A successful conservative cell remap has exactly one canonical overlap
    entry per source/target pair, covers every source and target, retains every
    source-row unity defect below the caller tolerance and the static `1e-6`
    ceiling, and retains its measured global extensive defect below the caller
    balance tolerance. Non-negative inputs plus admitted non-negative fractions
    cannot publish a negative target value. Signed-zero outputs and report
    totals canonicalize to positive zero.

## Error model

`MeshError` teaching errors: `TooFewPoints`, `DegenerateInput` (exact
all-coplanar detection, says to triangulate in 2D instead),
`Cancelled`. Kernel internals hold invariants by construction (no
flat tets: every created tet's apex is strictly visible); the audit
exists so any regression is LOUD rather than silently non-Delaunay.
`AdaptivityError` refuses a source state reused as its own target, before/after
QoI or source/target-state mismatch, contradictory action/effect declarations,
an unbacked continuous-gradient claim, non-finite or negative accounting, and
non-finite bound composition before any receipt is published. Digest parsing
deliberately adds no authority and accepts every 32-byte value, matching the
upstream identity adapter boundary.
`ConservativeRemapError` separately names empty/oversized requests, memory
admission, invalid tolerances/values/fractions/indices/order, missing local
coverage, excessive source-row or global balance defects, arithmetic overflow,
allocation refusal, and cancellation. A refusal owns no partial target vector.

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
Adaptivity accounting is a pure fixed-order operation. It canonicalizes signed
zero, preserves exact retained ID bytes, uses a strict total-bound comparison,
and serializes fields in one schema-fixed order (G3).

## Cancellation behavior

`delaunay` polls `cx.checkpoint()` every 256 insertions; `refine`
polls per round; `remesh` polls per iteration. Cancellation returns `MeshError::Cancelled` between
insertions (request → drain → finalize; no torn mesh states escape
since the error consumes the builder).
Adaptivity receipt admission is fixed-size synchronous metadata work and does
not accept a `Cx`; it publishes only after every input check completes.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

- `frontier-hexmesh` [F] (default OFF) — hex-dominant meshing (bead
  wqd.18); gates the `hexdom` integration target.

## Conformance tests

`tests/conformance.rs` defines 40 schema-validated fs-obs `ConformanceCase`
aggregates in a green run. The exact identities are the one-row cases `tmesh-001`,
`tmesh-002`, `tmesh-002b`, `tmesh-003` through `tmesh-012`, and `tmesh-017`;
`tmesh-013-threads-{1,2,4,8}`, `tmesh-013-audit-{1,2,4,8}`,
`tmesh-013-batch-width`, `tmesh-013-width-scaling`,
`tmesh-013-commutativity`, and `tmesh-013-degenerate-grid`;
`tmesh-014-{recovered,correspondence,audit-and-hull,replay}`;
`tmesh-015-zero-depth-existing-face` and
`tmesh-015-{recovered,coplanar,audit,replay,honest-caps}`; and
`tmesh-016-{recovered,tiles-L,audit,replay}`. Any reimplementation must pass
the suite unchanged.

Input-seed attribution follows the fixture, not the executor. `tmesh-001`,
`tmesh-002`, `tmesh-002b`, `tmesh-003`, `tmesh-004`, `tmesh-005`,
`tmesh-006`, and `tmesh-008` retain respectively
`0x1001_2026_0706_0021`, `0x1001_2026_0706_0022`,
`0x7115_ED00_C0B1_A11E`, `0x1001_2026_0706_0023`,
`0x1001_2026_0706_0024`, `0x1001_2026_0706_0025`,
`0x1001_2026_0706_0026`, and `0x1001_2026_0706_0028`. The
`tmesh-002` and `tmesh-002b` details name those roots because each aggregate
also contains fixed adversarial fixtures; each has only one stochastic input
root. `tmesh-011` and `tmesh-017` intentionally share
`0x1001_2026_0708_0011`; `tmesh-012` uses
`0x1001_2026_0708_0012`; every non-grid `tmesh-013-*` row uses
`0x1001_2026_0708_0013`; every `tmesh-014-*` row uses
`0x1001_2026_0708_0014`; the five non-zero-depth `tmesh-015-*` rows use
`0x1001_2026_0709_0015`; and every `tmesh-016-*` row uses
`0x160B_2026_0709_0016`. Fixed fixtures `tmesh-007`, `tmesh-009`,
`tmesh-010`, `tmesh-013-degenerate-grid`, and
`tmesh-015-zero-depth-existing-face` record input seed zero. Each listed LCG
root owns one sequentially derived fixture stream; no subordinate input stream
is silently replaced by the shared `Cx` provenance
(`0x7E7`, kernel 1, tile 0, iteration 0).

Eight object-shaped `Custom` companions retain the structured forensics:
`tmesh-001/measurement` (`mesh-delaunay-stats`),
`tmesh-005/measurement` (`mesh-refine-stats`),
`tmesh-006/measurement` (`mesh-scale-stats`),
`tmesh-007/measurement` (`mesh-remesh-iso`),
`tmesh-010/measurement` (`mesh-remesh-aniso`),
`tmesh-011/measurement` (`mesh-hull-split-evidence`),
`tmesh-012/measurement` (`mesh-sliver-exudation-evidence`), and
`tmesh-017/measurement` (`mesh-boundary-layer-pipeline-evidence`). Each has a
scope distinct from its aggregate, so fresh emitters cannot duplicate a
`(session, scope, sequence=0)` identity. Both companions and aggregates are
failure-linted, serialized, wire-validated, then printed. Ordinary aggregates
select Info/Error from their pass bit and perform the pre-existing terminal
assertion only after printing. For tmesh-011, tmesh-012, and tmesh-017, the
diagnostic remains at the old pre-assert boundary; their canonical passing
aggregate is deliberately withheld until all old gates have passed. Custom
measurements are diagnostic, not aggregate verdict authority, and this
observability migration does not enlarge any geometric, quality, determinism,
or performance claim. The primary target is default-feature compatible and
does not exercise `frontier-hexmesh`.
`tests/adaptivity.rs` adds G0/G3 admission, QoI-regression visibility,
state/lineage/accounting retention, effect and trigger coverage, a byte-pinned
schema-v1 JSON fixture, exact replay tests for the adaptivity receipt seam, and
piecewise-constant remap refinement/signed-cancellation, hostile overlap,
coverage-vs-balance, arithmetic-overflow, cap, and cancellation fixtures.
`tests/hexdom.rs`, cases hd-001..hd-005 behind `frontier-hexmesh`, emits
schema-validated fs-obs `ConformanceCase` verdicts and object-shaped,
wire-validated `Custom` measurements. hd-001 retains the actual MBO input-seed
family rooted at `0x5eed`; hd-002..hd-005 are fixed fixtures recorded with input
seed zero, and no execution/Cx seed is invented. Custom measurements are
diagnostic and do not substitute for or aggregate verdict authority. Existing
case-internal assertions may still abort before the verdict boundary; whenever
the verdict emitter is reached, it records an Info/Error event with a lintable
failure record, validates and prints the wire row, and only then performs its
terminal assertion. The target's central proof must explicitly enable
`frontier-hexmesh`; a default-feature workspace pass does not exercise it.

`tests/perf_lane.rs` remains an explicitly ignored ladder intended for the
documented release invocation. Each executed rung emits one failure-linted,
wire-validated `Custom` event under session `fs-mesh/perf-lane` at the dynamic
scope `mesh-perf/n-{n}/measurement`, so the 10⁴, 10⁵, 10⁶, and optional
10⁷ rows cannot reuse a sequence-zero identity. A compound row retains the
wall seconds, points/s, tet count, audit mode, the exact `debug_assertions`
mode, an escaped machine-configuration label, and its FNV-1a-64 label
fingerprint. The label uses OS, architecture, a best-effort CPU model, and
logical-CPU count; it identifies the recorded configuration and is not a
unique physical-host identity. It records cloud input seed `0xBEAD5EED`
separately from the Cx execution provenance (seed 41, kernel 77, tile 0,
iteration 0). The wall-clock values are finite-safe and explicitly
non-replayable, scoped only to that run and recorded configuration; `Custom`
is intentional because `BenchmarkResult` cannot bind the compound audit and
provenance receipt in one row. Without `FS_MESH_PERF_FULL=1`, the 10⁷ cloud, execution,
timing, and audit do not run; a Warn event instead records normalized status
`skipped`, the human instruction, null input/execution provenance, and the
configured seed at the distinct `mesh-perf/n-10000000/skip` identity. Event
emission remains at the old output boundaries, before the same throughput
assertions.

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
  rung (full insphere certification at 10⁴). The historical ledger rows
  remain in bead uee3; current runs emit canonical fs-obs evidence from
  `tests/perf_lane.rs`. tmesh-006 pins 10⁴-scale
  behavior in CI; the nightly perf-CI cadence belongs to fz2.4.
- Remeshing no-claims: curved creases round under midpoint projection
  (straight creases are exact); boundary loops are locked, not
  remeshed; metric gradation control, log-Euclidean metric
  interpolation/intersection, and DWR-supplied discrete metric fields
  join with FLUX's estimator bead; the functional-pass architecture
  trades throughput for auditability until the perf lane says
  otherwise.
- Adaptivity receipts are accounting artifacts, not DWR, conversion-error,
  topology, or gradient certificates. The remap kernel covers only one
  piecewise-constant cellwise extensive scalar and MEASURES algebraic f64
  balance; the caller still owns overlap geometry, units, geometric/projection
  error, admissibility of internal variables, higher-order accuracy,
  monotonicity beyond non-negative scalar input, vector/tensor frame semantics,
  and continuum conservation evidence. Marked-cell refinement, mesh
  untangling, and dependent invalidation remain unimplemented. Those algorithms
  must retain their own evidence and feed this seam rather than infer authority
  from it. Effect booleans are explicitly caller-declared and do not prove the
  named outcome. Opaque 32-byte IDs are one-way adapter boundaries for
  lower/higher-layer identity types, not a new competing identity scheme.
- `orient3d_sos` is a projection cascade, not the full 3D
  Edelsbrunner–Mücke ladder (fs-ivl's documented no-claim); it is used
  only for walk routing here, never for conflict decisions.

## No-claim boundaries (hexdom)

- v1 lattice tier: frame fields, singularity graphs, and extraction
  live on box-lattice domains; general tet-domain frame fields, true
  CubeCover parameterization, and curved hex extraction are the
  research core this tier deliberately does not claim.
- The SH9 sampling is a faithful linear image, not the coefficient
  basis; exact band-4 Wigner-D rotation is the named growth path.
- Scaled Jacobians are exact (1.0) only for the axis-aligned tier;
  warped-element quality arrives with the parameterization.
- IGA and CutFEM cover most hex use-cases at higher accuracy — this
  module's own refusal says so, by design.

## Boundary-layer quality: the measured decision (bead iw3l)

- tmesh-017 runs the refine(split_hull_facets) → exude pipeline on the
  ledgered convex-cloud fixture: exudation cuts the sliver census 27%
  (183 → 134) and lifts the worst dihedral off exact zero, and the
  longest/2·shortest EDGE-aspect of the final mesh is ~19 — the
  ledgered 3.5e17 "radius-edge" is CONFINED to near-coplanar hull
  slivers whose circumradius explodes while their edges stay tame.
- A hull-EDGE diametral protection tier was implemented and MEASURED
  COUNTERPRODUCTIVE (worst 3.5e17 → 4.3e18, reverted): convex-hull
  edges of a point cloud are not PLC features; the classical segment
  rule protects INPUT segments only.
- `split_hull_facets` therefore stays default-OFF: removing the
  residual near-coplanar hull sliver class needs WEIGHTED exudation on
  an exact weighted insphere predicate (the recorded fs-ivl no-claim)
  — the honest continuation, demanded-driven.
