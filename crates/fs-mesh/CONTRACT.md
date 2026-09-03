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
- `volumetricize(UnverifiedPlc, VolumetricPolicy, cx)`: constrained
  multi-region PLC volumetricization (bead s93ej.1). Type-state is
  `UnverifiedPlc` → `AdmittedPlc` → `ConstraintRecoveredPlc` →
  `LabeledTetComplex` → `AuditedLabeledTetComplex`. Each region is a
  closed oriented 2-manifold with a strictly interior seed. Segment and
  facet recovery reuse the existing conforming kernels; recovered faces
  are walls; seed-flood assigns chambers; exterior leftover that can
  reach a ghost without crossing a wall is carved; cavities are
  discarded; only the independently audited type is a geometry
  authority. The auditor reclassifies every retained tet by
  `winding_exact` of each region's own surface, requires positive
  orientation, and checks tet volume against the closed-surface
  triple-product identity with a second, distinct tet-volume formula.
- `rounded_cylinder_tet_mesh`: deterministic conforming P1 tetrahedra for a
  solid cylinder with equal circular outer-rim fillets. Geometry, fillet,
  radial/azimuthal/axial resolution, and count budgets are explicit. The
  cap/fillet tangent is a retained radial ring. Its derived boundary triangles
  carry outward normals, centroids, and areas for rendering/BEM consumers.
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
  throughput, until the perf lane profiles it. Scalar policy admission is
  allocation-free and precedes Soup cloning, cosine evaluation, polling, or
  metric work: `smoothing` is finite in inclusive `[0, 1]`, and
  `crease_angle` is finite in inclusive `[0, π]` radians so periodic cosine
  aliases cannot silently select a different threshold. Both signed-zero
  encodings are accepted and canonicalized to positive zero.
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
14. `volumetricize` on a closed axis-aligned multi-region PLC retains only
    declared solid tets, assigns exactly one region to each, carves hull
    fill that is not a declared solid, discards seeded cavities, and
    matches the closed-surface volume identity under two independent
    tet-volume formulas. Open, non-manifold, unlabeled-enclosed, boundary-seed,
    and duplicate-region inputs refuse without an audited mesh (s93ej.1
    corpus). A sheared (non-axis-aligned, exactly planar) pair of
    adjacent solids keeps both unit volumes and a labeled interface
    through the production `volumetricize` path. Parent conduction E2E
    remains a no-claim.
    Facet Steiner midpoints snap to the parent supporting plane, so a
    planar facet stays planar under bisection even when it is not
    axis-aligned.
    A facet counts as recovered when SOME set of coplanar mesh faces
    tiles it — its own bisected sub-triangles when they are all faces
    (the historical rows, byte-identical), otherwise any tiling whose
    only free edges lie on the facet boundary (vertex membership by
    provenance and the `1e-12` chord tolerance). The old sub-triangle
    test was a sufficient condition that co-circular ties defeat: an
    axis-aligned rectangle's corners share a circle, the kernel's tie
    break picks one diagonal, and bisection only manufactures smaller
    squares with the same tie. Facet recovery runs in PASSES until a
    whole sweep changes nothing, because a later insertion can flip
    away faces an earlier facet had conformed to; only the verification
    against the finished mesh decides, and the independent winding
    audit is what caught the single-sweep leak. A facet Steiner point
    within the chord tolerance of an existing vertex adopts it (no
    ulp-twin vertices). After carve-and-label and BEFORE the audit,
    `repair_flat_tets` removes zero-volume tets — coplanar co-circular
    quadruples the kernel's symbolic perturbation admits between two
    triangulations of one planar quad — by edge removal: the flat tet's
    in-plane diagonal is shared by a ring of tets on one side; the ring
    polygon is re-triangulated (a fan from a ring vertex, the first whose
    tets are all positively oriented, none flat, volumes conserved) and
    paired with both edge endpoints. Boundary faces, region labels and the
    union are unchanged by construction; the audit re-checks orientation,
    volume and winding, and the census after the pass is what consumers
    read. MEASURED on the comb: 3/220 and 9/722 flat tets, all removed;
    smallest dihedral 6.4° / 6.1° afterwards (`tests/comb_quality.rs`). Corpus: `tests/comb_prism.rs` — a finned
    heatsink comb prism (one, two and four fins; 36/60/108 facets)
    volumetricizes under the default policy with the analytic volume,
    where every earlier build refused with unrecovered facets at any
    budget. CONSEQUENCE, stated plainly: which coplanar tiling a facet
    ends up with follows the kernel's index-ordered tie-break, so the
    tet complex (and every downstream face numbering) depends on the
    input vertex ORDER as well as the point set; identical bytes still
    replay bit-for-bit, but a permutation of the input is a different
    mesh, not a relabelling of the same one. Consumers that need a
    permutation-invariant statement must state it on geometry (region
    identity, coordinates), never on face or slot indices.
15. Every public remesh call validates its two floating-point policy controls
    before geometry-dependent work. Exact endpoints admit; the adjacent
    representable value outside either interval refuses with stable field,
    rejected bits, and exact inclusive bound bits. Geometry translation or
    rescaling cannot change that scalar admission result (G0/G3).
16. `LabeledTetComplex::refine_uniform` is the h-ladder CONTROL, not a
    quality improver: one uniform 1→8 split (Bey/Zhang red refinement) of
    every retained tet by its edge midpoints — four corner copies in the
    parent's vertex order and four interior tets on the central
    octahedron's SHORTEST diagonal (canonical tie-break), each interior
    child given the parent's exact orientation sign by the predicate.
    Conforming by construction (a shared face splits identically on both
    sides); recovered walls (source faces) split four ways in their own
    plane with the parent facet inherited, so boundary classification by
    parent facet survives every rung; region labels replicate; per-region
    volume is preserved to the rounding of the midpoints; recovery and
    flat-repair evidence are carried unchanged as the base's provenance.
    Liu–Joe (1996) bounds every descendant's quality below by a constant
    times the parent's. MEASURED 2026-09-02 (`tests/uniform_refine.rs`,
    `src/uniform.rs`): two-fin comb 425 → 3,400 → 27,200 tets, smallest
    dihedral 6.86° → 4.738° → 4.738° (the interior class appears once and
    then persists; a regular tet and a skewed tet show the same over five
    generations), max radius-edge 45.0 unchanged (scale-invariant, so no
    uniform rung can improve it — that is constrained refinement's job).
    CONSEQUENCE for consumers: a dihedral floor applied to refined rungs
    must allow the first-generation class drop (here ×0.69); a 5° floor on
    the base would pass and then refuse rung 1 of this very comb.
17. Constrained refinement (`VolumetricPolicy::refinement`, off when `None`)
    runs between recovery and carving: rounds of worst-first circumcenter
    insertion restricted to the seeded chambers, every recovered wall
    (correspondence row) protected by its EQUATORIAL sphere — the smallest
    sphere through the three vertices, the only one whose emptiness keeps a
    face Delaunay (the hull code's minimum ENCLOSING sphere under-protects
    obtuse faces: guarding walls with it broke six comb facets in one
    round) — each round followed by segment and facet re-recovery and the
    exact audit, with the facet driver seeded INCREMENTALLY from the previous
    correspondence (`recover_facets_with_points`) so bisection repairs only
    what an insertion destroyed; re-recovery from a previous tiling with no
    insertions is a proven no-op (`tests/constrained_refine.rs`). A vertex
    that lies in a facet's plane strictly inside its loop is classified
    interior by geometry as well as by provenance. Evidence
    (`RecoveryEvidence::refinement`) discloses rounds, insertions, worst
    radius-edge before/after, offenders remaining, encroach-skips, wall
    splits and the stop reason. MEASURED 2026-09-02 on the two-fin comb
    with `split_walls: false` (the default): all 171 offenders' circumcenters
    encroach a wall, nothing is inserted, the mesh is untouched and the
    evidence says so — on thin bodies wall splitting is the whole game.
    `split_walls: true` (split the encroached wall at its in-plane point,
    handed to the driver as a known interior point) is opt-in and NOT yet
    claimed: 23 splits left 13 of 60 facets unrecovered within the recovery
    budget; the re-tiled facets were not recognised by the coplanar-tiling
    classifier (`tile:none` under `FS_MESH_TRACE_RECOVERY=1`), the open
    question for the next increment.
18. Delaunay kernel, coplanar ghost rule. A point exactly coplanar with a hull
    facet conflicts with that facet's ghost iff it lies strictly inside the
    facet's circumcircle, tested with the exact `insphere` against an apex
    lifted off the facet plane. The apex is `a` offset along the coordinate
    axis most aligned with the facet normal by the facet's longest edge —
    NOT `a + n`: on a hull sliver with nearly collinear vertices the normal is
    tiny (MEASURED 2026-09-02 on the rotated two-fin comb: |n| = 7e-20 against
    coordinates of 0.1), `a + n` rounds back to `a`, both predicates return
    Zero, `Zero == Zero` read as a conflict, the cavity's growth repair then
    absorbed a real tet that was not in conflict, the mesh ended with 227
    local-Delaunay violations and insertion of a later point swallowed hull
    vertex 9. Every sphere through a, b, c meets the plane in the same
    circumcircle, so the decision is apex-independent: a differential test
    over 2.4 million coplanar integer configurations shows zero disagreements
    with the previous rule wherever it was sound (`delaunay::ghost_rule_probe`).
    An exactly collinear ghost facet keeps the historical convention (it
    conflicts with everything and never survives the next insertion; the
    collinear-run battery mints such ghosts). The exact audit now also reports
    an input vertex absent from every live tet (bitwise duplicates exempt:
    whichever twin BRIO order met first is the present one), and
    `AdmittedPlc::recover` audits BEFORE counting unrecovered constraints so a
    kernel defect is named as one instead of as a recovery budget problem.
    Corpus: `tests/body_corpus.rs` — the rotated comb's point set keeps all 32
    vertices under the full audit (116 tets); the rotated comb volumetricizes
    with the exact volume (247 tets, 60/60 facets, 43 Steiner points, min
    dihedral 2.5°, no flat tets); a plate with a rectangular through-hole
    volumetricizes with the exact volume (180 tets).
19. Flat-tet repair, boundary cases. `repair_flat_tets` no longer refuses every
    flat tet that touches a wall; wall safety is decided per removal edge by
    the ring walk (a wall face on the ring refuses). A flat tet that sits IN
    the boundary — two wall faces, its other two faces shared with same-region
    tets, every vertex attached elsewhere — is DROPPED: rotation rounding
    leaves a segment Steiner point a hair off its edge, the Delaunay keeps the
    original triangle and mints a zero-volume tet between it and the
    degenerate sliver face (rotated comb: [0, 2, 3, 56], volume 1.9e-22), and
    neither diagonal is removable (one is a wall edge, the other's fan mints
    another zero-volume tet). Dropping removes zero volume, moves the two
    walls onto the tet's interior faces with the parent facet inherited, and
    the independent winding audit re-checks the result. Any other flat tet
    that resists removal stays disclosed as `unrepaired`.
20. Facet tilings close against their neighbours. A facet is recovered only
    by a tiling whose free edges are exactly the CHAIN of every mesh vertex
    on each of its edges (`recovery::chain_on_chord`: the endpoints and every
    vertex within the segment tolerance of the chord, in parameter order,
    evaluated with the endpoints in index order so both facets sharing a
    segment reach the same verdict). A face whose edge lies along a facet
    edge without being a chain sub-edge is never a tile, and a facet's own
    sub-triangle that is a mesh face but skips a chain vertex counts as
    missing, so longest-edge bisection adopts the vertex when its edge comes
    up (raw midpoints of dyadic chains coincide bitwise; `adopt_near` covers
    the ulp). After recovery, `AdmittedPlc::recover` and the constrained
    refinement re-recovery refuse with `Audit { reason: "recovered facet
    tiles do not close the region surface" }` unless every edge of each
    region's tile set is used exactly twice — the property the seed flood
    depends on and per-facet tiling cannot prove. MEASURED 2026-09-02 on the
    rotated four-fin heatsink shell: the crease edge (10, 13) between fin 1's
    top and side was tiled across the original edge by one facet and through
    the midpoint 114 minted on that edge by the other; both tilings passed
    the old test; the flood walked through the one sliver-shaped hole into
    the exterior, every zero-volume tet under the bottom plane was retained,
    and the winding audit tripped on the first of them far from the hole.
    The edge survived its own midpoint because the fin tops are a nearly
    coplanar HULL layer whose flat tets have circumspheres so large that the
    exact f64 midpoint, a rounding hair off the plane, falls outside them:
    Bowyer–Watson did not consume the whole star of the edge and a needle
    tet kept the edge alive beside its midpoint (14 such midpoints on that
    shell; none with eight far bounding-box points, which is not what ships).
    Dead end, measured: adopting the chain vertices by fanning the facet's
    triangle over them made interior splits run away (3393 of 3703 Steiner
    points on the four-fin comb that needs 136); balanced bisection does not.
    `FS_MESH_TRACE_RECOVERY` names the open edges, every midpoint that left
    its edge alive, and each unrecovered facet's anatomy; `FS_MESH_DUMP_MESH`
    writes the finished mesh as text for offline analysis.
21. Coplanar tiling is a one-sided sheet. Four coplanar co-circular points —
    every rectangle of a facet grid and every square bisection makes of one —
    are tetrahedralized with a zero-volume tet whose faces are BOTH
    triangulations of the quad (item 19), so once bisection has moved a facet
    away from its own sub-triangles the faces in its plane are a stack of
    double covers and no edge count can accept a tiling: on the rotated
    four-fin shell one base-top strip reached 4614 sub-triangles, twenty
    times its local feature size. `coplanar_tiling` therefore takes, when the
    in-plane faces are not already one clean tiling (that answer is kept,
    byte-identical), the sheet on one side of the facet's plane: the
    in-plane faces with exactly one incident tet whose apex lies strictly on
    that side — the boundary of that side's tets, manifold by construction;
    both sides are tried because a hull layer has real tets on one side only
    (`recovery::sheet_on_side`, unit-tested on a synthetic double cover). The
    flat tet between the two triangulations then lies either outside the
    walls (carved) or inside with two wall faces (item 19 drops it). Corpus:
    the rotated four-fin heatsink shell volumetricizes with the exact volume
    — 378 tets, 154 vertices, 98 Steiner points, 108/108 facets, two interior
    zero-volume tets disclosed — and the axis-aligned bodies are unchanged
    (four-fin comb prism 136 Steiner points).
22. Sliver repair by dihedral. `repair_flat_tets` treats a tet as removable
    when it is flat by volume (≤ 1e-9 of the largest, item 19) OR when its
    smallest dihedral angle is below 1° (`REPAIR_SLIVER_DIHEDRAL_DEG`, the
    conduction stage's mesh-quality floor): an STL imported at f32 precision
    (fs-io's weld) leaves the two triangulations of each facet a rounding
    hair apart, and the tets between them have a volume of 1e-6 of the
    largest — never "flat" — with dihedrals of 1e-6 degrees. Edge removal
    conserves the volume of the ring INCLUDING the removed tet (exact for a
    flat, its own volume for a sliver), and a fan that would mint a flat tet
    or a sliver is rejected, so every accepted flip strictly improves the
    census. Boundary drops stay reserved for true zero-volume flats: dropping
    a sliver would shave its volume off the region. MEASURED 2026-09-03 on
    the rotated four-fin shell at f32: 28 slivers → 2 (701 tets, exact
    volume); at f64 the two interior zero-volume tets of item 21 remain. The
    survivors are coplanar-cluster slivers whose ring apexes all lie in the
    same near-plane, so every re-triangulation of the ring mints another thin
    tet (`FS_MESH_TRACE_REPAIR` names the ring, the apex distances and the
    reason each fan was refused); item 23 clears them.
23. Needles and Steiner perturbation. Two more repairs finish the job. (a) A
    needle — a tet with one vertex ON the chord of one of its own edges (an
    original edge the kernel kept alive beside its midpoint, item 20, or two
    collinear input segments whose f32 rounding bent the line) — is repaired
    by removing THAT edge, whose star re-tiles through the on-chord vertex;
    its planar-quad diagonals are chord sub-edges and were rightly refused.
    (b) `volumetricize` then spends up to three rounds of Steiner
    perturbation on whatever survives: ONE Steiner vertex per surviving tet
    (moving every corner of a co-circular rectangle toward its centre only
    scales it — measured) is moved 5 % of the tet's shortest edge toward its
    centroid ALONG its own constraint — its input segment's chord, its wall
    tile's plane, or freely — the mesh is rebuilt from all points through
    the same kernel, recovery, carve and repair, and the result is kept only
    if the survivor count strictly drops; input vertices are never moved and
    a cancelled rebuild propagates. For the rebuild, `recover_segments`
    seeds each segment's chain with every mesh vertex already on its chord
    (a fresh PLC has none, so its chains are unchanged): bisection beside a
    moved chord vertex would otherwise mint midpoints until the depth cap
    (measured: 6 segments unrecovered). MEASURED 2026-09-03: rotated four-fin
    shell at f64 379 tets, smallest dihedral 1.64°, no flat tet, exact volume
    (survivors 2 → 1 → 0); the same shell at f32 (the CLI's PLC) 704 tets,
    smallest dihedral 1.09°, no flat tet, exact volume (2 → 0 in one round).
    `FS_MESH_TRACE_REPAIR` prints each round. Evidence fields for the rounds
    are deferred to the next conduction-receipt schema bump; the census
    (`flat_tets`, `min_dihedral_deg`) already discloses the outcome.

## Error model

`MeshError` teaching errors: `TooFewPoints`, `DegenerateInput` (exact
all-coplanar detection, says to triangulate in 2D instead), `InvalidFinite`
(preserves the rejected bits and refuses non-finite `crease_angle` or
`smoothing` before remeshing work), `InvalidControlRange` (preserves stable
field, rejected bits, and exact inclusive bound bits for finite policy
violations), and `Cancelled`. Kernel internals hold invariants by construction
(no flat tets: every created tet's apex is strictly visible); the audit exists
so any regression is LOUD rather than silently non-Delaunay.
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
Remesh control admission is bit-stable and geometry-independent. The admitted
intervals prevent scalar extrapolation and periodic policy aliasing; they do
not certify non-inversion, quality monotonicity, convergence, Newton-projection
stability, exact threshold robustness, or cross-ISA bitwise cosine
classification.

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

- Multi-region volumetricization does not claim quality/conformity
  beyond the recovered-face walls, the parent conduction E2E runner, or
  a cross-ISA mesh-byte identity. Facet Steiner midpoints are snapped
  onto the parent supporting plane; a point already on the plane is
  left bitwise unchanged. Constraint-edge (crease) midpoints stay raw
  so two non-coplanar parent facets adopt the same Steiner vertex
  (tmesh-019). Self-intersection of
  distinct region surfaces is refused only when both windings claim the
  same retained tet, not by a complete triangle-triangle predicate
  sweep. The discarded cavity/exterior volume fields on the witness are
  producer diagnostics and are not an independent certificate.
- The rounded-cylinder primitive is a piecewise-planar approximation of the
  exact circular meridian and circular revolution. Its output reports the
  meridian and azimuthal chord-error estimates; neither the volume mesh nor its
  boundary is an exact curved geometry. The v1 constructor refuses fillets that
  collapse the cylindrical side or planar cap instead of emitting degenerate
  elements. General line/arc profiles, annular topology, anisotropic sizing,
  and curved high-order volume elements remain outside this constructor.

- Weighted exact insphere predicate (the Edelsbrunner weight-pump
  exudation variant; the perturbation flavor ships).
- INTERIOR FACET recovery now ships in CONFORMING form for CONVEX
  planar facets (`recover_facets`, tmesh-015): batched longest-edge
  midpoint bisection of the fan triangulation (one-split-per-round was
  MEASURED to starve at the rounds cap; batching finished the fixture
  in 7 rounds), twin adoption via the shared coordinate-bits index,
  a facet correspondence table re-verified against the finished mesh,
  and honest starved-budget counters. Steiner midpoints snap to the
  parent supporting plane (tmesh-018), so a planar facet stays planar
  when it is not axis-aligned. Constraint-edge midpoints stay raw so a
  crease shared by two parent planes twins (tmesh-019). Remaining no-claim:
  Full-Ruppert QUALITY:
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
  honest `unrecovered` counters at depth/budget caps. Steiner midpoints
  snap onto the original parent chord (tmesh-020) so a general-position
  segment stays on the line. Crossing diagonals that no longer share a
  bitwise midpoint adopt the existing centre if it lies on the chord.
  Convex hull-facet conformity is gated test-side.
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
