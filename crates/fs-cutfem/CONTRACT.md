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
  Read-only `grid`, `active_cells`, and `cut_rules` accessors let
  goal/error estimators consume the exact topology and certified rules
  retained by the built space instead of reconstructing a second view.
- `FemParams`: `nitsche_beta` (applied as β/h), `ghost_gamma` (0
  disables), `quad_depth`, `agg: Option<AggPolicy>`, `strong_outer`,
  solver knobs.
- `CutElasticity`: vector Q1 small-strain elasticity on uniform or 2:1
  graded active quadtrees. Hanging nodes are reduced componentwise through
  deterministic terminal expansions: bulk and embedded-interface Nitsche
  terms use `T^T K T` / `T^T f`, equal-level ghost terms use the same
  componentwise transform, and solution reconstruction restores every active
  mesh node. Design-box traction is assembled on outer-edge terminal nodes;
  that transform is necessarily identity for a valid leaf partition and the
  assembler refuses if the topology invariant is broken. A literal
  no-constraint path retains the original uniform-grid numbering, COO/RHS
  insertion order, clamp mask, and topology bits. `IsotropicElastic` supplies
  the plane-strain Lamé parameters. Symmetric Nitsche imposes displacement
  data on the SDF interface with the cut-independent penalty
  `nitsche_beta * mu / h`; the componentwise ghost penalty scales as
  `ghost_gamma * mu * h` and controls degenerating cut fractions. The
  certified v1 material regime requires
  `(lambda + 2*mu) / mu <= 4`; larger ratios refuse rather than
  extrapolating the compressible-regime coercivity evidence toward
  incompressibility. An
  explicit natural,
  traction-free interface mode and zero design-box clamps support
  topology-optimization frontends. The legacy `boundary_traction` field is an
  uncertified full-edge callback, so every active design-box edge must be
  certified wholly inside or wholly outside. The explicit
  `assemble_with_boundary_traction` / `solve_with_boundary_traction` methods
  require that field to be `None` and accept either the same
  `BoundaryTraction::Uncertified` semantics or typed support.
- `DesignBoxEdge` / `EdgeBand` / `BoundaryTraction`: a checked closed
  normalized support interval on one named unit-box edge and its traction
  callback. `EdgeBand::new` refuses non-finite, out-of-range, or reversed
  endpoints; its fields are private. The applied traction is definitionally
  zero on all other edges and outside the band. Assembly intersects the band
  with each boundary-cell edge exactly, retains endpoint contact as
  potentially nonzero support, classifies that supported subsegment with the
  SDF enclosure, skips certified-outside support, refuses a straddling
  subsegment, and applies two-point Gauss only to wholly-inside support. It
  never uses callback samples or Gauss-node zeros to infer support and never
  performs generic SDF clipping.
- `CutElasticityOperator`: the assembled symmetric CSR operator, load,
  public deterministic terminal-node/block map, clamp mask, canonical active-cell /
  cut-rule / ghost-face topology, dropped-cut count, exact algebraic
  compliance `b^T x`, and `fs_solver::LinearOp` implementation. Its
  `apply_vec` / `apply_transpose_vec` conveniences expose the
  exact-symmetric seam. `CutElasticitySolution` carries coefficients,
  nodal values, active cells, the retained cut rules and canonical ghost
  faces, algebraic compliance, dropped-cut count, CG iterations, and final
  relative residual. On a graded tree, `node_ids()` contains only algebraic
  terminal blocks while `nodal_values()` and solved fields reconstruct all
  active hanging nodes. A requested hanging-node clamp refuses unless every
  terminal in its trace is clamped; clamping a terminal alone is valid. Its
  Q1 value/gradient evaluator refuses inactive cells, absent corners, and
  non-finite state rather than substituting zero.
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
    curved SDF domain gates the fixed three-level log-log Q1 convergence
    fit within 0.2 of the theoretical L2 = 2 and H1 = 1 orders, with
    strict monotone error decrease and both cut-position-sensitive
    adjacent slopes retained in the evidence row.
12. CUT-INDEPENDENT COERCIVITY (cte-003): the Nitsche penalty never
    scales with cut fraction. The vector Q1 acceptance family uses the
    conservative fixed full-element trace constant `beta = 100`. A
    cut-fraction sweep instead verifies that the vector ghost penalty
    bounds condition-number growth while the unstabilized reference
    deteriorates and is at least 100x worse on the most degenerate cut.
13. OPERATOR ADJOINT (cte-004): local matrices are canonically mirrored
    before COO accumulation, making stored CSR bit-symmetric. The
    content-addressed registered apply VJP is bit-identical to explicit
    transpose apply, remains correct with two operators in one registry,
    and passes deterministic central finite-difference directions.
14. GRADED VECTOR REDUCTION (cte-005, G0/G3): a closed-polygon affine
    two-component patch on a 2:1 tree passes the exact algebraic residual,
    CG residual, bit-symmetry, field-error, explicit midpoint reconstruction,
    deterministic replay, and a nonzero equal-level ghost-energy oracle. A
    separate private equal-level-face oracle derives a real hanging midpoint
    from a 2:1 tree and verifies every terminal basis against `T^T K T`, with a
    nonzero constrained-corner contribution.
    A ghost-off mixed-level interface independently compares the reduced
    Nitsche matrix energy and load pairing against the reconstructed physical
    field, with a nonzero Nitsche contribution from cut cells carrying
    constrained corners.
    Separate all-inside graded fixtures exercise nonidentity body-load
    reduction and graded design-box traction conservation against independent
    analytic loads; every outer traction node is asserted terminal. Synthetic
    corrupted graphs refuse cycles, empty/non-finite/non-unit transforms, and
    partially clamped hanging traces. The uniform pre-graded operator is frozen
    in both mostly-clamped (cte-000) and fully-unclamped (cte-000b) fixtures by
    portable whole-operator FNV-1a goldens over CSR, RHS, terminal map, clamp,
    active/cut topology, certified rule bits, ghost faces, and dropped count.
15. TYPED BOX-EDGE SUPPORT (G0/G3): invalid normalized bands refuse; a right
    edge band disjoint from an SDF crossing integrates its exact constant-load
    length; a crossing through the supported band refuses; an uncertified zero
    callback still refuses; and an all-inside band aligned to cell edges is
    bit-identical to the legacy callback operator evidence.

## Error model

`CutFemError` teaching errors: `EmptyDomain` (level set never enters
the grid), `CutBandNotUniform` (enabled ghost faces need an equal-level
interface band; names the offending pair and the repair; ghost-free
aggregation does not impose this precondition),
`InvalidFemInput` (a scalar stabilization parameter is non-finite or
outside its documented range, or a scalar field/goal evaluation has
missing, non-finite, inactive, or topology-inconsistent evidence),
`InvalidElasticityInput` (non-finite/non-coercive or unsupported
near-incompressible material, parameter, load, boundary data, or an
invalid edge band or SDF-cut supported box-edge segment; malformed/non-finite coefficient,
Q1 field-evaluation, terminal transform, or incompatible hanging-clamp
requests also refuse),
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
enclosures and fs-solver CG. The cte-000/000b portable goldens freeze the
uniform operator's exact public evidence surface; cte-005 separately
requires bit-identical graded replay.

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
solution, three-level fitted convergence slopes, and adjacent-slope
evidence; cte-003
condition-number sweep over cut fractions 0.5…1e-8
with fixed Nitsche scaling and ghost ON/OFF; fail-closed invalid-input,
material, and cut-loaded-edge coverage, including typed edge-band validation,
exact supported-subsegment classification, load conservation, legacy
full-edge refusal, and aligned legacy/typed operator-bit equivalence; cte-005
graded vector reduction, including an independent mixed-level Nitsche
matrix/RHS oracle, body-load reduction, terminal outer-traction conservation,
solution reconstruction, deterministic replay, and malformed-constraint
refusals.
`tests/elasticity_adjoint.rs`
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
- Mixed-level vector ghost faces are not claimed. Componentwise hanging-node
  reduction permits graded active trees, but enabled ghost stabilization still
  requires the cut/interface band itself to have equal-level face neighbors;
  `CutBandNotUniform` remains the structured refusal and child work owns the
  mixed-level face quadrature/topology redesign.
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
- `b^T x` is the exact algebraic assembled-load compliance. With nonzero
  embedded displacement data, the Nitsche data terms are part of `b`; this
  value alone is not claimed to equal physical external work.
- Certified clipping quadrature for a supported traction subsegment that the
  SDF itself cuts is not implemented; that configuration refuses rather than
  sample-masking or clipping a partial segment. A named band may safely avoid
  a different cut portion of the same design-box edge because zero support
  outside the band is part of its checked input semantics.
- `fs-cutfem` is the canonical owner of certified CutFEM elasticity and typed
  design-box traction support. `fs-solid::cutfront` owns only consumer-facing
  delegation for `fs-topols` / `fs-lattice`; it must not grow an independent
  clipping or boundary-support implementation. Completing every consumer
  migration remains separately tracked integration work.
- Coupling to fs-scenario BC descriptors and FrankenVDB tile storage
  (the dyadic alignment is designed in; the wiring is a consumer
  bead).
- Time-dependent problems and the speculative Delaunay-vs-CutFEM race
  (the executor's bead).
