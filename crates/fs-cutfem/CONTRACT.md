# fs-cutfem CONTRACT

## Purpose and layer

Layer: L3 (FLUX). CutFEM on SDFs (plan §8.1 frontend 2, bead tfz.8):
level-set geometry simulated with FEM-grade accuracy and zero meshing.
Quadtree background grids (the 2D restriction of the octree design,
FrankenVDB-dyadic-aligned), certified cut classification, cut
quadrature with error control, ghost-penalty stabilization,
aggregated-element fallback, Nitsche embedded Dirichlet conditions,
and estimate-driven h-refinement with hanging-node constraints. This
crate also owns a vector Q1, plane-strain elasticity frontend on the
same certified cuts; its constitutive parameters come from
`fs-material` rather than another raw-parameter formula in this crate.

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
- `CutElasticity`: vector Q1 small-strain elasticity on a uniform active
  quadtree level. `IsotropicElastic` supplies the plane-strain Lamé
  parameters. Symmetric Nitsche imposes displacement data on the SDF
  interface with the cut-independent penalty
  `nitsche_beta * mu / h`; the componentwise ghost penalty scales as
  `ghost_gamma * mu * h` and controls degenerating cut fractions. The
  certified v1 material regime requires
  `(lambda + 2*mu) / mu <= 4`; larger ratios refuse rather than
  extrapolating the compressible-regime coercivity evidence toward
  incompressibility. An
  explicit natural,
  traction-free interface mode and zero design-box clamps support
  topology-optimization frontends. Loaded design-box edges must be
  certified wholly inside or wholly outside each boundary cell; an
  SDF-cut loaded edge refuses until certified 1-D clipping exists.
- `CutElasticityOperator`: the assembled symmetric CSR operator, load,
  public deterministic node/block map, clamp mask, dropped-cut count,
  and `fs_solver::LinearOp` implementation. Its `apply_vec` /
  `apply_transpose_vec` conveniences expose the exact-symmetric seam.
  `CutElasticitySolution` carries coefficients, nodal values, active
  cells, dropped-cut count, CG iterations, and final relative residual.
- With `adjoint-vjp`, `register_elasticity_apply_vjp` registers the
  exact matrix under a BLAKE3 content-addressed key beneath the stable
  `fs-cutfem.elasticity-apply.v1` prefix, and returns that key for the
  caller's tape. Different matrices cannot overwrite one another in a
  shared registry. The VJP is the explicit symmetric transpose apply;
  differentiation is through fixed-matrix `Kx`, not through CG,
  material/geometry construction, or a `K^-1 b` solve.
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
   (cut-007); every aggregated node is policy-logged. A graded
   aggregation-only interface band builds and reproduces a linear patch
   while recording zero ghost faces (cut-007b); equal-level refusal is
   specific to enabled ghost stabilization.
8. Hanging-node constraints reproduce linears exactly on graded trees
   (patch test to solver tolerance, cut-008); the assembled system is
   symmetric positive definite on free DOFs.
9. Determinism: BTree-ordered traversal, straight-line IEEE
   arithmetic, deterministic CG — bit-identical across runs.
10. VECTOR PATCH LAW (cte-001, G0): for a basis spanning vector-valued
    affine displacement fields, the exact coefficient vector satisfies
    the assembled system to relative infinity residual < 2e-12 on three
    closed piecewise-linear polygons. Deterministic CG separately gates
    relative residual <= 1.1e-13 and absolute L2/H1 forward errors below
    2e-8/2e-7. Polygon vertices lie on background nodes while edges make
    arbitrary cell cuts, so the reused rule represents both interface
    segment and normal exactly.
11. VECTOR MMS ORDER (cte-002, G1): a manufactured displacement on a
    curved SDF domain gates both successive Q1 refinement slopes within
    0.2 of the theoretical L2 = 2 and H1 = 1 orders, with monotone error
    decrease.
12. CUT-INDEPENDENT COERCIVITY (cte-003): the Nitsche penalty never
    scales with cut fraction. A cut-fraction sweep instead verifies
    that the vector ghost penalty bounds condition-number growth while
    the unstabilized reference deteriorates and is at least 100x worse
    on the most degenerate cut.
13. OPERATOR ADJOINT (cte-004): local matrices are canonically mirrored
    before COO accumulation, making stored CSR bit-symmetric. The
    content-addressed registered apply VJP is bit-identical to explicit
    transpose apply, remains correct with two operators in one registry,
    and passes deterministic central finite-difference directions.

## Error model

`CutFemError` teaching errors: `EmptyDomain` (level set never enters
the grid), `CutBandNotUniform` (enabled ghost faces need an equal-level
interface band; names the offending pair and the repair; ghost-free
aggregation does not impose this precondition),
`InvalidFemInput` (a scalar stabilization parameter is non-finite or
outside its documented range),
`ElasticityGridNotUniform` (the vector frontend refuses graded active
cells until componentwise hanging constraints exist),
`InvalidElasticityInput` (non-finite/non-coercive or unsupported
near-incompressible material, parameter, load, boundary data, or an
uncertified SDF-cut loaded box edge),
`AggregationNoAnchor` (no well-cut anchor within 4 BFS rings),
`ConstraintCycle` (corrupted constraint graph), `SolveNotConverged`
(configured residual gate missed; carries iterations and residual).
Quadrature and classification never silently degrade: scalar
conservative-Cut and dropped-cell counts are surfaced in `BuildStats`;
the vector operator and solution surface their dropped-cut count.

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

The plan marks CutFEM-on-SDF [F]; per the crate-granular form of the
Ambition-Tag gating rule (the fs-feec precedent), the frontier surface
ships as this standalone crate and consumers opt in by depending on it.
`adjoint-vjp` additionally enables fs-adjoint's frontier
`ledger-transpose` surface and the elasticity-apply VJP registration.
The core vector operator and solver do not require that feature.

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
accuracy within 2×, policy logged); cut-007b graded interface-band
aggregation with ghost OFF (successful build/linear patch, zero ghost
faces); cut-008 hanging-node patch test
(exact to 1e-8) + estimate-driven refinement efficiency (within 2× of
uniform-fine error at < 0.55× its DOFs). Grid unit tests: 2:1 balance,
leaf partition of unity.

`tests/elasticity.rs`: cte-001 affine-basis vector patch law over three
exactly represented piecewise-linear closed cuts, with a roundoff-scale
algebraic residual gate distinct from explicit solver forward-error
tolerances; cte-002 G1 curved-circle, nonzero-divergence manufactured
solution and both successive convergence slopes; cte-003
condition-number sweep over cut fractions 0.5…1e-8
with fixed Nitsche scaling and ghost ON/OFF; fail-closed invalid-input,
material, and cut-loaded-edge coverage. `tests/elasticity_adjoint.rs`
is a declared `adjoint-vjp`-required target: cte-004 exact registered
transpose, two-operator key isolation, and independent central-FD gate.
The acceptance tests emit deterministic JSON verdict rows suitable for
capture by the proof lane; stdout alone is not claimed as retained
ledger evidence.

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
- Vector elasticity on graded active trees: the scalar hanging-node
  constraint machinery is not silently reused componentwise. The
  vector frontend returns `ElasticityGridNotUniform` until that lift is
  implemented and proven.
- Vector constitutive scope is two-dimensional plane-strain isotropic
  small-strain elasticity. Plane stress, orthotropy, nonlinear
  material state updates, finite strain, and higher-order vector
  elements are not claimed.
- Nearly incompressible plane strain is not claimed. The v1 frontend
  refuses `(lambda + 2*mu) / mu > 4`; extending the `mu`-scaled
  stabilization beyond that compressible regime requires a separate
  material-robust coercivity analysis and locking/conditioning battery.
- The vector operator currently assembles canonical CSR. It exposes a
  matrix-free-shaped `apply` / `apply_transpose` seam, but a genuinely
  matrix-free cut-cell kernel is not yet claimed.
- The apply VJP covers a fixed assembled matrix and state vector only.
  A solve VJP, material/SDF/design derivatives, and a topopt compliance
  adjoint remain consumer/integration work.
- Certified clipping quadrature for nonzero traction on an SDF-cut
  design-box edge is not implemented; that configuration refuses rather
  than sample-masking a partial edge.
- `fs-solid::cutfront` and its `fs-topols` / `fs-lattice` consumers still
  use a legacy independent CutFEM elasticity facade. Migrating that
  facade to delegate to this operator is separate follow-up work; this
  bead does not claim workspace-wide de-duplication or consumer migration.
- Coupling to fs-scenario BC descriptors and FrankenVDB tile storage
  (the dyadic alignment is designed in; the wiring is a consumer
  bead).
- Time-dependent problems and the speculative Delaunay-vs-CutFEM race
  (the executor's bead).
