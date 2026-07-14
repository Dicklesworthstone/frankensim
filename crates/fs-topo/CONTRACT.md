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
  the outward-orientation winding probe. Every defect is a typed,
  LOCALIZED `ManifoldDefect`; `certified()` ⟺ zero defects.
- `self_intersection_certificate(soup)` → `SelfIntersectReport`:
  sweep-and-prune broad phase; EXACT narrow phase — plane-separation
  early exits, then exact edge-vs-triangle tests (four `orient3d`
  signs each; complete for non-coplanar pairs because every
  intersection-segment endpoint lies on some edge), exact `orient2d`
  for the coplanar case. Faces sharing a vertex are adjacency-excused.
  A PASS is a PROOF (exact arithmetic — false PASS impossible); exact
  contact reports the conservative `Touching` kind (the bounded,
  LISTED false-FAIL class the acceptance contract allows).
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
    deterministic (voxels sorted by value then index);
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
  compares against).

## Invariants

1. Clean fixtures certify; punched holes localize exactly their
   boundary edges; duplicated faces read use-count 3; flipped patches
   read misoriented edges; degenerate faces are named (topo-001).
2. Clean surfaces are PROVEN free; planted piercings read `Crossing`
   with pair localization; exactly coincident patches read `Touching`;
   near-tangent surfaces at 1e-4 separation do NOT false-FAIL
   (topo-002). Historical note: this certificate's first run caught a
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

## Error model

Certificates are total over admitted inputs. `VoxelizeError` wraps
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

`tests/conformance.rs`, cases topo-001..topo-008 — JSON-line
verdicts, seeded LCG randomness, the fs-obs scale ledger. Any
reimplementation must pass the suite unchanged.

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
