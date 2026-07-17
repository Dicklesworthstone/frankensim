# CONTRACT: fs-rep-mesh

## Purpose and layer
Mesh charts (plan §7.2): half-edge surfaces with an edit core, genuine
oriented 2-D triangle complexes and oriented tet complexes with
exact-arithmetic incidence (the pre-FEEC δδ = 0
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
- `TriComplex2` — an admitted oriented 2-D triangle cell complex, distinct
  from a surface half-edge mesh. It retains explicit topological dimension 2,
  embedding dimension 2, oriented face rows, canonical edges, exact integer
  d0/d1, caller-owned stable vertex keys, typed vertex/edge/face `EntityId`s,
  and prevalidated cell measures. Construction rejects malformed indices,
  payloads without a two-cell, repeated keys/cells, non-manifold or
  same-direction shared edges, non-finite coordinates, negative axisymmetric
  radii, degenerate faces, and non-representable measures. Its coordinate,
  stable-key, face, and canonical-edge tables are sealed behind read-only slice
  accessors because incidence rows, typed feature IDs, and measures are cached
  derived state; topology or embedding edits require constructing and admitting
  a new complex.
- `TriComplex2LineageId` / `TriFeatureId` — schema- and role-typed entity
  identities from `fs-blake3`. A lineage binds an exact caller-owned namespace;
  a feature binds that lineage, its topological dimension, and the canonical
  set of stable vertex keys. Storage order, coordinates, face order, and face
  orientation are deliberately not feature-identity inputs. Parsing or digest
  equality supplies no authority or provenance claim.
- `TraceMap2` / `TraceEdge2` — deterministic selected-side trace extraction.
  Exact d1 coefficients cancel selected-selected edges; selected-unselected
  edges become interface traces, and the complete face selection yields the
  outer boundary. Trace vertices and edges carry explicit parent maps and an
  exact trace-local d0.
- `Metric2::Planar` / `Metric2::Axisymmetric` — explicit measure metadata.
  Planar faces integrate area times a finite positive thickness. Axisymmetric
  coordinates are `(radius, z)` and integrate a finite positive sweep of at
  most one turn using exact linear-radius simplex quadrature; this weighting
  does not change the stored complex's topological or embedding dimension.
- `Soup` / `triangle_winding` (van Oosterom–Strackee) / `winding_exact` /
  `WindingOctree` — the Barill-style dipole hierarchy (β accuracy knob,
  area-weighted normal moments, exact leaves); octree stores indices
  only, callers pass the soup (no self-referential borrows).
- `point_triangle_distance` (Ericson region test, exact),
  `ray_triangle_watertight` (Woop permutation+shear; shared edges NEVER
  leak; exact-edge hits may double-count — documented for parity users),
  `Bvh` (median-split, deterministic tie-break by index; closest-point
  branch-and-bound + nearest-hit raycast).
- `MeshChart` — sd = BVH closest distance × winding sign. The unsigned
  distance-to-set magnitude is 1-Lipschitz, but a generic soup's winding sign
  is not topology-certified and can jump on open/self-intersecting/inconsistently
  oriented input. Raw `MeshChart` therefore retains the trait-default
  `TraceStepClaim::NoClaim`, returns `lipschitz: None`, and grades finite
  nominal samples as `Estimate` (non-finite samples as `NoClaim`) even for a
  clean-looking closed soup; `raycast` is watertight.
- `repair(soup, max_hole_edges)` — dedupe → degenerate removal →
  orientation unification (flood fill + centroid winding vote) →
  fan-fill of small boundary loops, each action a `RepairReceipt`
  (defect/location/action, the fs-io quarantine format).
- `shapes` (PUBLIC fixture vocabulary): `cube`, `icosphere`,
  deterministic `corrupt` (dups/degens/flips/hole).
- `dual_contour` / `dual_contour_clipped` / `DcOptions` /
  `bracket_certificate` — the SDF→mesh
  converter (plan §7.3 edge 2): uniform-grid dual contouring with
  explicit finite sampling-domain admission, finite positive cell spacing,
  checked corner-lattice products, and a 256-cells/axis cap before chart
  evaluation or allocation. Lattice coordinates use normalized interpolation
  over the admitted span and hit both admitted endpoints exactly, so a
  ceil-derived cell count cannot step past the admitted maximum. Cell spacing,
  QEF regularization, every nominal field sample, chart/fallback gradients,
  and all derived crossing/QEF coordinates are finite-validated. Edge secants
  scale opposite field magnitudes before division, Hermite means use convex
  online averaging, and QEF solves in cell-local coordinates, so finite
  extreme inputs either remain representable or return a structured refusal
  rather than leaking NaN/overflow into a mesh. Unbounded charts require `dual_contour_clipped`,
  which contours the geometric intersection `chart ∩ clip`, not merely a
  replacement grid extent. The converter uses
  Hermite data (secant crossings + chart gradients), regularized 3×3 QEF
  vertex placement (Schaefer mass-point regularization; `MassPoint` is
  the feature-blurring baseline), axis-uniform quad winding (cyclic
  (u,v) rings, circulation normal +d, flipped by crossing sign), and THE
  certificate: v1 accepts only the GLOBAL
  `TraceStepClaim::ExactDistance` theorem (therefore `L = 1`) and requires a
  finite rigorous `Exact`/`Enclosure` from `trace_value_enclosure` containing
  every nominal centroid evaluation. Recursive 4-way subdivision proves
  `max_abs(enclosure(centroid)) + radius ≤ tol` over the surface; local
  `ChartSample::lipschitz` values, `Estimate`, `NoClaim`, malformed evidence,
  invalid tolerances/indices, and non-representable geometry cannot produce a
  pass or fail. Completed rigorous failures localize to triangles with
  margins. `BracketCertificateError` distinguishes every refusal and
  cancellation with progress; `NoLipschitz` remains a compatibility type
  alias.
- `mesh_to_sdf` / `assess_quality` / `IncrementalMeshSdf` — the authority-graded
  mesh→SDF converter (plan §7.3 edge 1): point-triangle distance + winding sign
  sampled onto fs-rep-sdf dense grids. Edge-use, pairwise edge orientation, and
  aggregate signed-volume checks are diagnostics only: they do not certify
  per-component orientation/nesting, vertex-link manifoldness, or freedom from
  self-intersection. A disconnected soup with a large outward component and a
  smaller inward component can pass every diagnostic. Consequently every
  generic mesh payload and receipt is capped at Estimate and NAMES the
  `winding-sign-heuristic`, even when the basic screen passes. Defective inputs
  additionally record the screen's defect counts; out-of-range public soup
  indices are counted and force the screen false without being dereferenced or
  panicking. A private converter-only
  adapter supplies the unsigned distance-to-set's unit slope solely as a
  nominal reconstruction heuristic; it explicitly retains
  `TraceStepClaim::NoClaim`, forwards Estimate/NoClaim sample authority, and
  cannot authorize sphere tracing or an enclosure. Weak/NoClaim payloads remain
  NoClaim and can never be promoted by mesh quality. The QoI is the total
  abstract-distance estimate bound when available, otherwise the finite nominal
  reconstruction bound paired with a NoClaim certificate. The incremental path
  transactionally re-samples only a dirty box at exactly the original positions
  (bit-identical to full rebuild — G5); cancellation or refusal preserves the
  prior chart, field, authority, and refreshed-sample evidence as one consistent
  state. Initial and successfully edited `IncrementalMeshSdf` values apply the
  same Estimate cap before exposing their raw TiledSdf; a prior downgrade is
  never strengthened by an incremental edit.

## Invariants
1. Half-edge invariants survive 2k random flips with the Euler
   characteristic preserved (rmesh-001).
2. `point_triangle_distance` never exceeds any sampled surface distance
   and matches 1830-point brute force within sampling gap; `MeshChart`
   tracks the analytic sphere within the icosphere chord band, satisfies
   inside ⇔ sd < 0, and exhibits 1-Lipschitz behavior on that fixture
   (rmesh-002); this measurement is not a generic signed-field theorem.
   rmesh-002c locks the raw authority boundary on a clean-looking closed cube:
   default `NoClaim`, no Lipschitz bound, finite `Estimate`, non-finite
   `NoClaim`, and bracket-certificate refusal.
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
7. The converter matches analytic SDF fixtures within its recorded estimate, is
   translation-equivariant (G3), refreshes incrementally bit-identically
   (G5), keeps generic clean and defective soup at most Estimate, rejects
   mixed-component authority laundering, and never promotes weak sampled payload
   authority (rmesh-007 and focused converter unit tests).
8. Dual-contoured fixtures are manifold, closed, outward-oriented
   (winding +1 at the center — the ring-orientation law), vertex-accurate,
   translation-equivariant, bracket-certified, and the certificate DETECTS
   a fixed adversarial bad triangle with localized margins; QEF resolves the box
   corner at least twice as sharply as the mass-point baseline
   (rmesh-008).
9. Contouring admission (rmesh-009): default contouring rejects unresolved
   extended support before evaluation; the paired clipped API contours the
   actual geometric intersection; invalid spacing and checked grid-cap
   refusals also precede evaluation; every admitted positive axis has at least
   one cell even when `span / spacing` underflows, and near-integer quotients
   are incremented when their realized width would exceed the requested `h`;
   translated chart+clip pairs produce
   translated meshes with identical connectivity (G3).
10. Bracket authority (rmesh-010): `Some(1)` on a local sample cannot promote
    a `NoClaim` theorem, `Estimate`/`NoClaim` evidence, or malformed `Exact`
    evidence; a chart-requested cancellation is observed directly after the
    chart evaluation and returned with completed triangle/evaluation progress.
    Existing sphere/box exact-distance certificate coverage remains in
    rmesh-008.
11. The 2-D complex battery (rmesh-011) refuses face-less payloads rather than
    labeling a lower-dimensional object as topological dimension 2, checks
    d1∘d0 = 0 exactly over seeded admissible fans, and pins the complete square
    incidence tables. It refuses a one-face flip while a coherent global
    reversal negates d1 and boundary orientation, and matches hand-computed
    whole-boundary and selected-face trace tables. It preserves all typed IDs
    across cyclic/global orientation changes, and preserves existing vertex plus
    surviving boundary-edge IDs across append-only refinement. It reproduces
    planar and full-turn axisymmetric face/edge closed forms including a
    legitimate zero-measure axis edge. Compile-fail seals and detached-copy
    mutations prove safe callers cannot desynchronize admitted tables from
    cached incidence, identities, measures, or traces.

## Error model
Structured teaching errors (`MeshBuildError`, `TriComplex2Error`,
`Metric2Error`, `ContourError`, `BracketCertificateError`, and wrapped
`SamplingDomainError`). TriComplex2 errors retain a missing two-cell condition,
the offending cell, edge, face pair, coordinate bits, trace selection, or
canonical-identity refusal; no partial complex or feature table is published.
Contour errors
name invalid spacing/regularization, excessive resolution, checked
grid/coordinate overflow, non-finite chart samples or gradients,
non-representable derived arithmetic, cancellation, and empty zero sets.
Bracket refusals name invalid tolerance,
unsupported trace theorem, invalid indices/non-finite vertices,
non-representable geometry, malformed/non-rigorous evidence, and cancellation
progress. Total functions elsewhere (degenerate triangles yield well-defined
distances; empty-soup handling is the caller's constructor discipline). No
panics across the boundary.

## Determinism class
Deterministic: BTreeMap/BTreeSet orders, index-tie-broken BVH sorts, seeded
batteries; no clocks, no addresses in results. `TriComplex2` incidence and
trace tables are integer and storage-deterministic. Typed feature IDs are
bit-stable functions of the exact lineage ID, topological dimension, and
canonical stable-key set; admissible refinements that retain those inputs do
not move surviving IDs. Floating metric measures are deterministic exact-order
f64 evaluations on one implementation/ISA, not cross-ISA bit claims.

## Cancellation behavior
`MeshChart::eval` is bounded per query (BVH descent); batch consumers
poll between queries per the fs-geom discipline. Dual contouring polls directly
through its lattice-sampling, Hermite, vertex-placement, and stitching phases.
It also polls immediately after every chart evaluation, so a producer that
requests cancellation without polling cannot have its sample consumed.
Bracket certification polls `Cx` directly before and after each chart
evaluation (and after trace-enclosure retrieval), including recursive
subdivision, and refuses with completed triangle/evaluation progress instead
of publishing a partial verdict.
Incremental mesh-to-SDF refreshes stage all samples before committing any
chart/field state. Repair/build are bounded preprocessing passes (cancellation
hooks join with the fs-io quarantine bead where soups get large).

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. `[S]` solid-tier.

## Conformance tests
`tests/conformance.rs` has 13 aggregate case IDs (`rmesh-001`, `rmesh-002`,
`rmesh-002b`, `rmesh-002c`, and `rmesh-003` through `rmesh-011`) covering
invariants 1–11. Every reached aggregate result is an fs-obs
`ConformanceCase` with Info/Error severity, passes `lint_failure_record`, is
serialized and wire-validated, and is printed before a failing aggregate
assertion. Assertions, expectations, and panics reached before that point are
ordinary Rust test diagnostics: they stop the case before an aggregate verdict
exists and are not laundered into synthetic evidence.

Randomized verdicts carry their literal campaign-root input seed: rmesh-001
`0x1001_2026_0706_0001`, rmesh-002 `0x1002_2026_0706_D157`, rmesh-002b
`0x1002_B026_0706_DE9E`, rmesh-003 `0x1003_2026_0706_50FA`, rmesh-004
`0x1004_2026_0706_D1B0`, rmesh-006 `0x1006_2026_0706_DD00`, rmesh-007
`0x1007_2026_0706_C0DE`, and rmesh-011 `0x1011_2026_0714_DD00`. The
rmesh-002 chart-law substream derives as campaign root xor `0xC047`, while its
canonical verdict retains the campaign root. Fixed rmesh-002c, rmesh-005,
rmesh-008, rmesh-009, and rmesh-010 use input seed zero. The fixed Cx stream
seed `0x9E54` is separate execution provenance, recorded in verdict detail and
relevant Custom companions rather than presented as input randomness.

The dipole-error, repair-receipt, mesh-conversion, and dual-contouring Custom
companions remain object-shaped, fs-obs-validated, and printed. Randomized
companions carry a standalone numeric `input_seed`; Cx-backed companions also
record the standalone numeric `execution_seed` when applicable. The
dual-contouring companion encodes an unavailable non-finite certificate margin
as JSON `null`, never as a non-JSON numeric token.

If the fixed rmesh-008 bracket certificate completes with localized failures,
the suite emits a distinct `rmesh-008/bracket` Error `ConformanceCase`
diagnostic before continuing to the aggregate verdict. That diagnostic is
failure-record linted and wire-validated, uses input seed zero for the fixed
chart, and records the separate fixed Cx execution seed in its detail. It does
not assert on the bracket outcome or replace the rmesh-008 aggregate gate;
schema lint and wire validation remain fail-closed.

## No-claim boundaries
- Self-intersection FLAGGING is deferred to the validity-certificates
  bead (wqd.7), which owns certified broad/narrow phases via fs-ivl.
- `HexComplex` is storage only; its incidence operators land with
  fs-feec's tensor-product families.
- `TriComplex2` v1 is straight-sided and simplicial. It does not certify a
  curved-cell geometry, embedding injectivity, self-intersection freedom,
  domain topology, mesh quality, material/formulation correctness, or a 2-D
  FEEC mass matrix. The weighted-operator consumer owns those later claims.
- Stable vertex keys and lineage namespaces are caller authority. Typed feature
  hashes prevent role/schema confusion but do not prove that a key was minted
  by machine IR, that two namespaces denote the same entity, or that a
  refinement crosswalk is honest. Reindexing the convenience
  `from_indexed_triangles` constructor moves keys; durable graphs must provide
  explicit stable keys.
- Axisymmetric v1 assumes the first coordinate is nonnegative radius and uses
  the declared angular sweep. It does not certify an axis crossing, curved
  meridian interpolation, orientation of a generated 3-D body, or the physical
  applicability of an axisymmetric formulation.
- Curvature via discrete operators (cotan/normal-cycle) is deferred to
  its first consumer, with convergence-class documentation there.
- NO throughput claims (million-triangle dipole performance is the perf
  harness's); the octree's asymptotics are structural.
- The winding sign convention (+1 inside outward-oriented closed surfaces) is
  pinned by fixtures, not by a formal proof of the flux form. Generic soup has
  no component/nesting/self-intersection certificate and therefore no global
  ExactDistance promotion. Exact-edge raycast double-counting is documented.
- Attribute channels (normals/UVs/per-face data) beyond positions are
  deferred until a consumer (LUMEN materials) defines their semantics.
- The converter's incremental mode trusts the CALLER's dirty box (it
  records refreshed-sample counts for audits); automatic change-support
  inference joins with the edit-tracking half-edge attributes.
- Sparse VDB output and per-tile interval audit certificates for the
  converter join with the fs-ivl integration pass (the dense path's
  fp-envelope analysis is documented in fs-rep-sdf).
- A nominal sampled-field reconstruction bound is not abstract-region signed
  distance authority. Basic mesh quality diagnostics cannot strengthen a
  TiledSdf payload; generic mesh receipts are at best Estimate, and NoClaim
  payloads remain NoClaim and non-certifiable.
- Dual contouring is UNIFORM-GRID v1: adaptive octree contouring with
  crack-free stitching, the guaranteed-manifold MDC variant for
  ambiguous topologies, and quality post-passes (min-angle smoothing
  with certificate re-verification) are follow-ups; the fixture zoo's
  manifoldness is verified, not guaranteed for adversarial topology.
- The exact-distance centroid/radius bracket is conservative; v1 deliberately
  refuses `LipschitzImplicit` because that theorem preserves sign/zero set but
  does not provide a global Euclidean proximity upper bound. More specialized
  global theorems or fs-ivl triangle evaluation may tighten the same
  certificate surface later.
