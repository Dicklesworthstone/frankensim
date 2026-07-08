# fs-cutfem CONTRACT

## Purpose and layer

Layer: L3 (FLUX). CutFEM on SDFs (plan §8.1 frontend 2, bead tfz.8):
level-set geometry simulated with FEM-grade accuracy and zero meshing.
Quadtree background grids (the 2D restriction of the octree design,
FrankenVDB-dyadic-aligned), certified cut classification, cut
quadrature with error control, ghost-penalty stabilization,
aggregated-element fallback, Nitsche embedded Dirichlet conditions,
and estimate-driven h-refinement with hanging-node constraints.

## Public types and semantics

- `CutSdf` (trait): level-set function φ (negative inside) with
  `value`, `gradient`, and `enclose` — a CERTIFIED fs-ivl interval
  containing φ over any axis-aligned box. `Circle`, `HalfPlane` are
  the reference fixtures.
- `Quadtree` / `CellKey` / `NodeKey`: 2:1-balanced dyadic tree over
  the unit square; cells `(level, i, j)`, nodes as lattice coordinates
  at `max_level` (dyadic, exactly representable). `refine_where` is
  the adaptivity hook (dwr-adaptivity drives it); the specialized
  `refine_toward_interface` establishes the uniform interface band by
  refining every cell whose SELF-INFLATED box encloses zero;
  `refined_once` is the every-leaf-split copy fs-dwr's enriched
  adjoint solves on (bead tfz.23 addendum).
- `cut_cell_rules` → `CutRules`: bulk points (weights sum to inside
  area) and interface points (weights sum to interface length, with
  outward unit normals) for one cut cell; `depth` is the error-control
  knob.
- `Space::build(grid, sdf, params)`: classify → quadrature → hanging /
  aggregation / strong-outer constraints → free-DOF numbering → ghost
  faces. `assemble` produces the SPD stiffness + load for −Δu = f,
  u = g (Nitsche on Γ, strong on the outer box when enabled);
  `solve` runs Jacobi-preconditioned fs-solver CG; `l2_h1_error`
  measures against an exact solution with one-deeper cut quadrature.
- `FemParams`: `nitsche_beta` (applied as β/h), `ghost_gamma` (0
  disables), `quad_depth`, `agg: Option<AggPolicy>`, `strong_outer`,
  solver knobs.
- `AggPolicy`: `small_fraction` / `good_fraction` thresholds of the
  documented aggregation policy.
- `condition_estimate` → `CondReport`: full-spectrum (Jacobi
  rotations) conditioning of fixture-sized systems, gated at n ≤ 4096.

## Invariants

1. CLASSIFICATION IS CERTIFIED: a cell is Inside/Outside only when the
   fs-ivl enclosure excludes zero; every tangency lands in Cut. No
   misclassification, ever (cut-001, adversarial tangent battery).
2. Cut quadrature converges quadratically in `depth` on curved
   interfaces and is EXACT (to roundoff) for quadratic moments under
   linear level sets (cut-002). Sub-resolution features that never
   flip a corner sign can be missed by quadrature (bounded by one
   sub-cell area) — never by classification.
3. Ghost penalty restores conditioning independent of the cut
   fraction: eigenvalue-verified flat conditioning curves across cut
   fractions down to 1e-8, versus >1e3 blowup without (cut-003).
4. MMS orders hold across RANDOMIZED cut configurations with slivers
   deliberately included: median L2 order ~2, H1 order ~1 (cut-004).
5. Nitsche embedded-BC accuracy is within 3× of a body-fitted Q1
   reference at matched resolution on the shared strip fixture
   (cut-005).
6. A moving level set re-solves on one FIXED background grid — zero
   re-meshing — with stable accuracy (cut-006).
7. Aggregation (ghost OFF) restores conditioning by ≥100× on the
   pathological sliver and holds accuracy within 2× of ghost-ON
   (cut-007); every aggregated node is policy-logged.
8. Hanging-node constraints reproduce linears exactly on graded trees
   (patch test to solver tolerance, cut-008); the assembled system is
   symmetric positive definite on free DOFs.
9. Determinism: BTree-ordered traversal, straight-line IEEE
   arithmetic, deterministic CG — bit-identical across runs.

## Error model

`CutFemError` teaching errors: `EmptyDomain` (level set never enters
the grid), `CutBandNotUniform` (ghost faces need an equal-level
interface band; names the offending pair and the repair),
`AggregationNoAnchor` (no well-cut anchor within 4 BFS rings),
`ConstraintCycle` (corrupted constraint graph), `SolveNotConverged`
(residual gate missed; carries iterations and residual). Quadrature
and classification never silently degrade: conservative-Cut and
dropped-cell counts are surfaced in `BuildStats`.

## Determinism class

Bit-deterministic across runs on a fixed platform (BTree iteration
order, no threading, no ambient rounding-mode state). Cross-ISA
determinism inherits fs-math/fs-ivl discipline through fs-ivl
enclosures and fs-solver CG; no golden hash is recorded yet for this
crate (recorded follow-up).

## Cancellation behavior

Bounded synchronous loops (classification, quadrature, assembly, CG
with an iteration cap). Chunking to tile quanta with Cx poll points
between chunks is the fs-exec driver's job (the fs-feec/fs-la
discipline). No internal threading.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no unsafe capsules.

## Feature flags

None. The plan marks CutFEM-on-SDF [F]; per the crate-granular form
of the Ambition-Tag gating rule (the fs-feec precedent), the frontier
surface ships as this standalone crate and consumers opt in by
depending on it.

## Conformance tests

`tests/battery.rs`: cut-001 certified classification (six adversarial
fixtures, zero misclassification + enclosure containment law);
cut-002 quadrature exactness (linear-interface quadratic moments to
1e-15; curved-interface area/length depth-order ≈ 2); cut-003
ghost-penalty conditioning independence (eigenvalue-verified curves,
cut fractions 0.5…1e-8); cut-004 MMS orders across 7 randomized cut
configs incl. a 1e-9 sliver (median L2 order ≥ 1.8, H1 ≥ 0.85);
cut-005 embedded-BC vs body-fitted Q1 (ratio < 3, order ≥ 1.7);
cut-006 11-step moving interface on one fixed grid (stability ratio
< 2.5); cut-007 aggregation fallback (conditioning restored ≥100×,
accuracy within 2×, policy logged); cut-008 hanging-node patch test
(exact to 1e-8) + estimate-driven refinement efficiency (within 2× of
uniform-fine error at < 0.55× its DOFs). Grid unit tests: 2:1 balance,
leaf partition of unity.

## No-claim boundaries

- 3D octree instantiation (the design is dimension-generic — dyadic
  keys, face constraints, tensor rules — but only 2D ships here).
- Higher-order elements (Qk, k ≥ 2) and the matching higher-order
  ghost penalties (derivative jumps beyond first order).
- Moment-fitting quadrature (tessellation-based scheme ships; the
  moment-fitting alternative is future work under the same `CutRules`
  surface).
- DWR-driven refinement loops: this crate exposes the `refine_where`
  hook and proves the hanging-node machinery; the dual-weighted
  estimator itself is dwr-adaptivity's bead.
- Coupling to fs-scenario BC descriptors and FrankenVDB tile storage
  (the dyadic alignment is designed in; the wiring is a consumer
  bead).
- Time-dependent problems, non-Poisson operators, and the speculative
  Delaunay-vs-CutFEM race (the executor's bead).
