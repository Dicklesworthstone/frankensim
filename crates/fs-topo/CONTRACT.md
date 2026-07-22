# fs-topo — CONTRACT

## Purpose and layer

L2 (MORPH). Validity and topology certificates (plan §7.8): three
certificate families, none of them sampling heuristics — manifoldness
with defect LOCALIZATION, self-intersection freedom as an EXACT PROOF,
and cubical topology (exact Betti numbers, true 0-dimensional
persistence) as the topology oracle ASCENT's castability/routing
constraints will consume.

## Public types and semantics

- `manifold_certificate(soup, interior_probe)` → `ManifoldReport`:
  combinatorial checks (edge-use census with direction bookkeeping,
  half-edge round-trip with the builder's vertex-link teaching
  errors) plus geometric red flags (degenerate faces, fold-overs) and
  an OPTIONAL outward-orientation winding probe.
  - The geometry is ADMITTED before it is read. Out-of-range vertex
    indices and non-finite coordinates are recorded as
    `VertexIndexOutOfRange` / `NonFiniteVertex` defects and every
    position-dependent check (zero-area faces, fold-overs, outwardness) is
    SKIPPED. Those checks are threshold comparisons that a NaN quietly
    falsifies, so without admission a NaN soup certified.
  - `manifold`, `closed` and `consistently_oriented` are purely
    COMBINATORIAL (index-only) and therefore stay meaningful on a refused
    soup. `consistently_oriented` means "every shared edge is traversed in
    both directions"; it makes NO claim about which side the normals face
    — a uniformly reversed mesh is consistently oriented.
  - Outwardness is a separate three-state `outward: Option<bool>`.
    `None` is a NO-CLAIM (no probe supplied, combinatorics already failed,
    geometry refused, or the winding was indeterminate), never a pass.
  - Combinatorial and geometric defects are typed and LOCALIZED to their
    face/edge. The GLOBAL inward-orientation verdict is
    `InwardOrientation { probe, winding }` — it is not disguised as a
    `MisorientedEdge` on a sentinel edge, because the edge-use census has
    already proved no edge is misoriented when that branch runs.
  - `certified()` ⟺ zero defects AND `outward == Some(true)`: a probe must
    have CONFIRMED outwardness. `combinatorially_certified()` is the
    weaker probe-free claim (closed, manifold, consistently oriented, no
    defect, nothing said about facing).
- `self_intersection_certificate(soup)` → `SelfIntersectReport`:
  sweep-and-prune broad phase; EXACT narrow phase — plane-separation
  early exits, then exact edge-vs-triangle tests (four `orient3d`
  signs each; complete for non-coplanar pairs because every
  intersection-segment endpoint lies on some edge), exact `orient2d`
  for the coplanar case.
  - ADMISSION comes first and its refusals are recorded, never absorbed:
    out-of-range vertex indices, non-finite coordinates (a NaN is DROPPED
    by the broad phase's `min`/`max` into an inverted box that culls every
    pair containing it, so `orient3d`'s own finiteness assertion is
    unreachable), and EXACTLY-degenerate faces (a zero-area face has no
    supporting plane, so every plane-separation argument is vacuous on
    it). `proven_free()` is false whenever any refusal is present, and
    `admitted()` exposes that directly.
  - Faces sharing a vertex or an edge are NOT skipped; they are decided
    against their SHARED FEATURE, exactly. Two shared vertices with
    distinct planes meet in exactly the shared edge (both planes contain
    its line, and a non-degenerate triangle meets the line through one of
    its own edges in exactly that edge); coplanar hinge pairs reduce to
    one `orient2d` per apex. One shared vertex `P` with distinct planes
    reduces to two segment-triangle tests on the edges OPPOSITE `P` — the
    chord `Tₖ ∩ (π₁ ∩ π₂)` has `P` as an endpoint, so the intersection
    exceeds `{P}` iff the far end of the shorter chord lies in the other
    triangle; coplanar shared-vertex pairs reduce to four exact
    corner-cone containment tests. Identical corner sets report
    `Touching`. `shared_feature_pairs_tested` counts the pairs so decided.
  - A PASS on an ADMITTED soup is a PROOF that no two faces meet outside
    their shared vertices/edges (exact arithmetic — no false PASS from
    rounding, and none from an excused pair either). Exact contact reports
    the conservative `Touching` kind (the bounded, LISTED false-FAIL class
    the acceptance contract allows); so does a coplanar overlap and an
    unresolvable degenerate corner.
  - `tri_tri_intersect` requires finite coordinates (`orient3d` asserts
    it) and returns the conservative `Touching` for an exactly degenerate
    triangle rather than claiming disjointness it cannot prove.
- `cubical::voxelize` / `voxelize_clipped` / `betti` / `persistence0` /
  `count_persistent` / `verify_topology` / `verify_topology_clipped`:
  - chart voxelization admits a finite positive-volume sampling domain and a
    nonzero longest-axis resolution before evaluation, derives every dimension
    with checked arithmetic, and refuses fields above the deterministic
    1,000,000-cell cap; ratio-first per-axis center placement stays inside the
    exact admitted spans even for extreme aspect ratios, and each derived
    dimension is incremented when nearest quotient arithmetic would make its
    realized width exceed the reported maximum `h`; non-finite source samples
    refuse; the paired clipped APIs sample the geometric
    intersection `chart ∩ clip`, not merely a replacement extent;
  - `betti`: exact Betti triple of the voxel solid — `b0` by
    26-connected union-find, `b2` as 6-connected bounded complement components
    against a virtual outside, `b1 = b0 + b2 − χ` with χ counted
    EXACTLY on the closed cubical complex (a k-cell is present iff an
    incident voxel is filled);
  - `persistence0`: true 0-dimensional persistence of the sublevel
    filtration — elder rule at every merge, essential classes kept,
    deterministic (voxels sorted by value then index); every `Bar` retains the
    exact `birth_index` whose activation created that component, so disconnected
    components with equal scalar endpoints remain distinguishable;
  - `verify_topology`: chart-level Betti at resolution `n`, HONESTLY
    framed — exact for the voxel solid, Estimate-grade for the chart
    (sub-cell features can be missed).

- `penalty` module (plan §9.5/§7.8, bead 7tv.15; [M], behind
  `moonshot-topo-persistence`): Betti targets as diagram penalties.
  `TopoSpec` (components/tunnels/enclosed-voids targets + τ, the
  feature-size floor), `evaluate` → `TopoPenalty` (graded total, zero
  iff the diagram matches the target up to τ; per-violation
  `Attribution` with voxel sets and fill/carve directions; the H₀ bars
  and Betti evidence ride along), `enclosed_voids` (the DUALITY route:
  an enclosed void is an empty-phase component that never reaches the
  domain boundary — localizes H₂ with union-find, cross-checked
  against `betti`'s b₂), `apply_attribution_step` (the descent
  primitive), `heuristic_cc_penalty` (the fallback the [M] gate
  compares against). Excess-component attribution floods from each retained
  H₀ birth representative below that bar's death/target-level cap; it never
  re-derives identity from a possibly shared scalar birth value.

## Invariants

1. Clean fixtures certify WITH an interior probe; punched holes localize
   exactly their boundary edges; duplicated faces read use-count 3;
   flipped patches read misoriented edges; degenerate faces are named
   (topo-001). A uniformly reversed surface stays
   `consistently_oriented` and reports the GLOBAL `InwardOrientation`,
   not a sentinel edge; without a probe `outward` is `None` and
   `certified()` is false while `combinatorially_certified()` holds
   (ci-topo-016/017). Non-finite coordinates and out-of-range indices are
   named as defects and suppress every position-dependent check
   (ci-topo-015).
2. Clean ADMITTED surfaces are PROVEN free of any face-pair contact
   outside the pairs' shared vertices/edges; planted piercings read
   `Crossing` with pair localization — including a piercing by a face
   that SHARES a vertex with its victim (ci-topo-014); exactly coincident
   patches read `Touching`; near-tangent surfaces at 1e-4 separation do
   NOT false-FAIL, and neither do the tetrahedron's six shared-edge pairs
   or the icosphere's shared-vertex fans (topo-002, ci-topo-014).
   Non-finite coordinates, out-of-range indices and exactly-degenerate
   faces REFUSE — `proven_free()` is false and the refusal is listed
   (ci-topo-015). Historical note: this certificate's first run caught a
   real latent bug in `shapes::icosphere` (off-origin midpoint
   projection) — the zoo's fixtures are honest witnesses.
3. Betti triples read exactly on the fixture zoo: ball (1,0,0), solid
   torus (1,1,0) — the tunnel via Euler duality — hollow ball (1,0,1),
   two balls (2,0,0) (topo-003).
4. Two planted wells stay exactly two persistent features under noise
   (1 essential + the shallow well with its analytic birth), against
   dozens of short noise bars (topo-004).
5. Stability: an ε-perturbation moves every surviving bar's endpoints
   by ≤ ε (the bottleneck stability theorem as a property test)
   (topo-005).
6. Persistence is BITWISE reproducible, and the ~10⁶-voxel scale run
   is ledgered with timings (topo-006).
7. Edge- and corner-touching closed voxels form one 26-connected
   contractible component and do not manufacture phantom tunnels
   (topo-007).
8. Sampling admission (topo-008): default voxelization and topology
   verification reject unresolved extended support before evaluation; paired
   clipped APIs sample the actual geometric intersection; zero resolution and
   checked voxel-cap refusals precede evaluation; translated chart+clip pairs
   preserve occupancy and Betti numbers (G3).
9. Equal-minimum disconnected H₀ components retain distinct birth
   representatives and disjoint excess-component attributions; one guided
   descent step repairs both islands without collapsing either attribution onto
   the first scalar match (tp-001b).

## Error model

Certificates are total: an input they cannot admit yields a typed,
enumerated REFUSAL rather than a verdict. The mesh certificates admit
in-range vertex indices, finite coordinates, and (for self-intersection)
exactly-nonzero face area; refusals ride in `SelfIntersectReport::refusals`
(`SelfIntersectRefusal`) and in `ManifoldReport::defects`
(`NonFiniteVertex`, `VertexIndexOutOfRange`, `IndeterminateWinding`), and
they force `proven_free()` / `certified()` false. `VoxelizeError` wraps
`SamplingDomainError` and names zero resolution, nonrepresentable cell size,
checked count overflow, deterministic work-cap refusal, non-finite chart
samples, and cancellation with completed-voxel progress.
Conservative flags (`Touching`) are typed and listed rather than silently
merged with strict crossings.

## Determinism class

Fully deterministic: BTree censuses, sorted sweep orders with index
tie-breaks, smaller-root union-find, value-then-index filtration
order. Identical inputs give identical reports bitwise (topo-006).

## Cancellation behavior

`voxelize` (and therefore `verify_topology`) polls `cx.checkpoint()` at most
every 256 completed voxels and once before publication. Mesh certificates are
single-pass and non-blocking.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

- `moonshot-topo-persistence` [M] (default OFF, bead 7tv.15) —
  persistence-based topology penalties; gates the `penalty` integration
  target.

## Conformance tests

`tests/conformance.rs`, cases topo-001..topo-008 — canonical fs-obs
verdicts, seeded LCG randomness, and a distinct fs-obs scale companion at
`topo-006/scale/measurement`. Any reimplementation must pass the suite
unchanged.

`tests/claim_integrity.rs` (cases ci-topo-014..ci-topo-017, named for the
beads `frankensim-extreal-program-f85xj.2.14`..`.2.17`) pins the
epistemic boundaries above with the filed minimal repros: a
vertex-sharing piercing must FAIL while the tetrahedron's six shared-edge
pairs and the icosphere's shared-vertex fans still PASS; NaN coordinates,
out-of-range indices, a non-finite probe and exactly-degenerate faces must
REFUSE both certificates; a uniformly reversed surface must report a
global `InwardOrientation` and never a `MisorientedEdge { edge: [0, 0] }`;
and an unprobed report must leave `outward` at `None` with `certified()`
false. Plain `#[test]`s, no fs-obs verdicts — they gate claims, not
conformance.

With `moonshot-topo-persistence`,
`tests/penalty.rs`, cases tp-001..tp-005, emits canonical fs-obs aggregate
verdicts. Fixed-input fixtures carry seed 0; tp-003 carries its actual jitter
input seed `0x1234`. The m-gate's graded trace is a finite-safe Custom event at
the collision-free `tp-004/measurement` identity, distinct from the tp-004
aggregate. `tp_001b_equal_minimum_islands_keep_distinct_birth_representatives`
locks the retained-representative, disjoint-attribution, and one-step repair
semantics for equal-valued disconnected islands.

## No-claim boundaries

- Persistence PAIRING is 0-dimensional; 1/2-dimensional persistence
  pairs (full cubical matrix reduction) and persistence-diagram
  matching penalties for ASCENT are the follow-up ([F/M] scope).
- Sequential reduction only; the chunked-PARALLEL reduction for
  10⁸⁺-voxel topology-optimization fields is routed to the perf lane
  with the deterministic-merge requirement attached.
- `verify_topology` is Estimate-grade at its resolution; interval-
  certified topology (no missed sub-cell features) is the sheaf
  certificates bead's.
- Geometric manifoldness is red-flag level (degeneracies, fold-overs);
  full local-injectivity proofs join the interval machinery.
- The self-intersection certificate REFUSES exactly-degenerate faces
  rather than analysing them: it has no answer for a face with no plane,
  and every pair containing one is excluded from the narrow phase. A
  refused soup carries NO non-intersection claim, not a weaker one.
- `ManifoldReport::outward` answers a SINGLE probe point. It certifies
  the winding sign there, not that the probe is interior, and not that
  the surface bounds a solid.
- A `ManifoldReport` with `outward == None` says nothing at all about
  facing. Do not read it as "not inward".
- Voxel connectivity is the closed-cube-consistent 26/6 pair
  (solid/complement); alternative connectivity conventions are not exposed.

## No-claim boundaries (penalty)

- TUNNELS (H₁) are counted at the level, not localized: the excess/
  deficit penalty is count-graded with no attribution map — the full
  cubical boundary-matrix PH is the growth path.
- Component-DEFICIT violations carry no localization either (there is
  nowhere to point when material must appear ex nihilo).
- The void "persistence" is its depth (min margin over the component) —
  a sublevel-filtration quantity on the dual phase, inheriting PH
  stability empirically (tested), not by the classical theorem verbatim.
- Resolution caveats inherit from `cubical`: features thinner than a
  voxel are invisible; τ must be chosen above the voxel scale.
- The tp-004 m-gate companion is deterministic fixture evidence for the
  stated attribution-vs-count-only comparison. It is not a wall-clock
  performance result or a general convergence claim for arbitrary fields.
