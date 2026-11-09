# fs-cutfem CONTRACT

## Complete 3-D strong-support force balance

Cartesian and adaptive elasticity expose `equilibrium_load` and
`equilibrium_response`. Complete body/surface loads and optional embedded
motion lifting retain supported entries; the original solve still masks its
strong rows. Recovered residuals include original scaled cell, B-bar,
Nitsche and ghost actions. Adaptive forces use the full hanging-node
transpose. Strong-support force/moment and the free residual norm remain
separate from embedded reactions.

Original full-field material admission and cancellation apply before
publication. A solved field is required to interpret supported residuals as
equilibrium reactions. These are discrete variational forces, not sampled
continuum traction. Existing masked/prepared load stencils cannot recover
discarded supported entries; use the complete load path. The six
`equilibrium3` regressions cover extension, supported loading, cut/graded
balance, mixed supports, pressure and refusal/retry behavior.

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
  refining every cell whose SELF-INFLATED box encloses zero, then repeatedly
  splitting the coarser side of any remaining cut-cell face mismatch. This
  closure also repairs finer neighbors introduced by independent adaptive work;
  `refined_once` is the every-leaf-split copy fs-dwr's enriched
  adjoint solves on (bead tfz.23 addendum).
- `cut_cell_rules` → `CutRules`: bulk points (weights sum to inside
  area) and interface points (weights sum to interface length, with
  outward unit normals) for one cut cell; `depth` is the error-control
  knob.
- `CutSdf3` / `HexCell` / `isolate_certified_height_root`: the narrow
  3-D Saye-style height-line primitive. It deterministically selects a
  derivative-certified monotone height direction and returns either one
  interval-enclosed root, a proven no-intersection result, or `Ambiguous`.
  It is not yet a 3-D volume/surface quadrature rule.
- `Space::build(grid, sdf, params)`: classify → quadrature → hanging /
  aggregation / strong-outer constraints → free-DOF numbering → ghost
  faces. `assemble` produces the SPD stiffness + load for −Δu = f,
  u = g (Nitsche on Γ, strong on the outer box when enabled);
  `solve` runs Jacobi-preconditioned fs-solver CG; `l2_h1_error`
  measures against an exact solution with one-deeper cut quadrature.
  Read-only `grid`, `active_cells`, and `cut_rules` accessors let
  goal/error estimators consume the exact topology and certified rules
  retained by the built space instead of reconstructing a second view.
- `Space::sample_scalar(nodal, p)` → `ScalarSample`: the one supported
  pointwise scalar read over a nodal field (bead ay40). `Active(v)` is
  a bilinear value backed by four present, finite corner values on the
  containing active leaf; `CertifiedOutside` is a build-time
  classification certificate (the leaf's SDF enclosure proved φ > 0),
  never a value default — callers own its physical meaning explicitly.
  Out-of-box or non-finite points, missing corner evidence, and
  non-finite corner evidence refuse as `InvalidFemInput`; zero is never
  fabricated. The background box is half-open, so rim probes must be
  clamped by the caller deliberately.
- `FemParams`: `nitsche_beta` (applied as β/h), `ghost_gamma` (0
  disables), `quad_depth`, `agg: Option<AggPolicy>`, `strong_outer`,
  solver knobs.
- `CutElasticity`: vector Q1 small-strain elasticity on uniform or 2:1
  graded active quadtrees. Hanging nodes are reduced componentwise through
  deterministic terminal expansions: bulk and embedded-interface Nitsche
  terms use `T^T K T` / `T^T f`, and ghost terms integrate exact dyadic
  equal- or mixed-level `SharedFacePatch` overlaps with positive-axis
  orientation and `h_F = min(h_a, h_b)` before using the same componentwise
  transform. Solution reconstruction restores every active mesh node.
  Design-box traction is assembled on outer-edge terminal nodes;
  that transform is necessarily identity for a valid leaf partition and the
  assembler refuses if the topology invariant is broken. A literal
  no-constraint path retains the original uniform-grid numbering, COO/RHS
  insertion order, clamp mask, and topology bits. `IsotropicElastic` supplies
  the plane-strain Lamé parameters. Symmetric Nitsche imposes displacement
  data on the SDF interface with the cut-independent penalty
  `nitsche_beta * modulus / h`; the componentwise ghost penalty scales as
  `ghost_gamma * modulus * h` and controls degenerating cut fractions. The
  effective modulus is configured via `CutStabilizationScaling`:
  `MuScaled` (canonical default, modulus = $\mu$) or `LongitudinalModulus`
  (physical modulus = $\lambda + 2\mu$, evaluated directly without dimensionless intermediate ratios). The
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
  residual as a typed `ResidualClaim::TrueEuclidean`, recomputed by an explicit
  `b - Ax` operator application rather than copied from CG's recursive estimate.
  If that true residual misses the requested gate after the recursive estimate
  stops, deterministic CG correction passes solve the latest true residual;
  every pass shares the configured aggregate iteration cap, and only a fresh
  `b - Ax` application can establish success. The legacy public magnitude and
  the fail-closed `euclidean_rel_residual()` accessor are derived from that
  same claim. On a graded tree, `node_ids()` contains only algebraic terminal
  blocks while `nodal_values()` and solved fields reconstruct all active
  hanging nodes. A requested hanging-node clamp refuses unless every terminal
  in its trace is clamped; clamping a terminal alone is valid. Its
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

## Native material fields in three-dimensional elasticity

`elastic3::ElasticMaterialField3` binds the existing admitted
`fs_material::elastic::acoustic::ElasticAcousticState` to Cartesian or adaptive
Q1 cut cells, uniformly or by physical background-cell bounds. The full global
Mandel tensor, density, source identity and original strain-domain shape remain
attached to each cell. Native isotropic/orthotropic adapters retain their
Frobenius ball; explicitly declared Mandel-component boxes retain that box. Bulk stiffness, Nitsche tractions and lifting, reactions,
stress observations and their pullbacks consume that same reference law.
Uniform fields share a prepared engineering tensor; no new material resolver,
constitutive admission or isotropic proxy is introduced.

Native tensor penalties use the retained maximum absolute Mandel row sum, a
conservative spectral upper bound. A face averages each side's reference
modulus times its own current stiffness scale; both scale derivatives retain
those distinct coefficients. The legacy isotropic constructor keeps its Lamé
bulk/traction formulas and shear-modulus ghost coefficient. These choices do
not certify coercivity, freedom from locking or contrast-independent conditioning.

`reference_cell_masses_kg` integrates original density over retained numerical
cut volume, including supported material. It is independent of the stiffness
interpolation and has no invented unit-density fallback for legacy inputs.
Adaptive transfer requires each child's complete native reference state to
match its coarse ancestor; child masses are re-integrated from their own
quadrature. Strain validity is checked on actual physical fields and complete
affine uncertainty families, never by treating algebraic basis strains or
arbitrary adjoint vectors as physical states. Geometry derivatives,
material evolution and continuum certification are outside this frozen-linear
operator's scope.

## Explicit isotropic mean-dilatation elasticity

`ElasticityOptions3::formulation = ElasticityFormulation3::BBar` selects
small-strain mean-dilatation elasticity on Cartesian and balanced-octree Q1
cells, including embedded displacement supports. The default `Displacement`
formulation retains its existing full integration and material restrictions.
B-bar admits finitely compressible isotropic laws, including explicitly
isotropic native material states; a generic or anisotropic Mandel tensor is
never silently replaced by an isotropic proxy.

The energy is `mu*dev(epsilon):dev(epsilon) + K*mean(div u)^2/2`.
The mean integrates the original PHYSICAL cut-volume rule, not the full
background box. Its gradients are cached once; assembly integrates the
volumetric contribution as one rank-one block and leaves the deviatoric
quadrature intact. No new pressure coordinates, global solve, per-apply
allocation or geometry/material reevaluation is introduced. The selected
formulation is retained through density updates, and inter-grid transfer
refuses a formulation mismatch. B-bar ghost jumps use the original shear
modulus, including native isotropic states, so the extension penalty does not
reintroduce a bulk-modulus constraint.

Bulk stresses use the projected volumetric strain. Nitsche stiffness, lifting
and reactions share its traction; physical/reference stress VJPs and all
scale contractions differentiate the same discrete formulation. Native source
admission checks BOTH actual Q1 strain and projected constitutive strain,
including complete affine uncertainty families. Projection cannot conceal a
material-domain exit. Relative displacement recovery removes rigid
translations before stress/domain contractions. Original density, material
identity, consistent inertia and physical-occupancy derivatives are unchanged.

`elastic3_bbar` specifies physical-cut averaging, bending energy refinement,
native inventory, isotropy admission, scale derivatives and adaptive transfer.
`elastic3_bbar_stress` specifies independent stress/work, state/scale VJPs,
raw and projected domain exits, rigid translations, embedded displacement and
weighted reaction checks. These are numerical formulation tests. B-bar relieves
volumetric locking; it does not remove shear locking, support exactly
incompressible materials, or establish an inf-sup/coercivity theorem for every
cut, pressure mode, penalty, contrast or refinement pattern. Conditioning and
physical convergence remain properties of the actual resolved problem.

## Prescribed eigenstrain in two-dimensional elasticity

`CutElasticityOperator::with_eigenstrain` adds cellwise bilinear
`Q1Eigenstrain` to the existing plane-stress or plane-strain operator. The
caller supplies global tensor-shear strains at the four Q1 corners, including
source-integrated isotropic expansion when applicable. The retained field is
immutable; duplicate attachment is refused. Geometry, material, constraints,
stiffness, Nitsche support and solver ordering remain with the original owner.

The added load is the physical-volume `B^T C epsilon_free` integral with the
negative eigenstress traction on embedded displacement supports. Plane stress
uses the original reduced constitutive law, including out-of-plane relaxation;
plane strain retains constrained out-of-plane stress and elastic storage.
`strain_stress` distinguishes total, free and elastic strains and reports
physical stress and recoverable elastic energy. No-eigenstrain observations
retain their existing arithmetic and material checks.

The caller supplies a separate total Frobenius-strain budget in `(0,1)`;
the original material calibration applies to elastic strain. Both domains are
checked at all Q1 corners, bounding the entire retained bilinear tensor field.
This is a small-strain constitutive model, not finite-deformation or thermal
back-reaction. The usual controlled solve enforces these two admissions before
publishing a solution. A new temperature/free-strain field requires attaching
it to an unattached operator, which callers may retain by cloning beforehand.

`elastic_eigenstrain` regressions cover stress-free expansion, restrained
compression, plane reductions near incompressibility, spatial compatibility,
natural and Nitsche boundaries, hanging nodes, source retention, cancellation
and distinct elastic/kinematic envelopes.

## Prescribed eigenstrain in three-dimensional elasticity

Cartesian and adaptive `prepare_eigenstrain` retain prescribed global-frame
tensor strains on the original physical cut rules. Thermal callers supply
source-integrated expansion from their actual stress-free reference; the same
interface also supports authored transformation or manufacturing strain. The
original material tensor supplies `C:(epsilon-epsilon_free)`. No replacement
material, invented expansion coefficient, density update or geometry rebuild is
introduced. The preparation is bound to the original immutable operator owner.

Reference cell columns contain `integral B^T C epsilon_free dV` and the
**negative** eigenstress traction on selected embedded Nitsche supports.
Current cell scales multiply both these loads and the original stiffness.
B-bar projects both total and free dilatation over physical cut volume while
retaining local deviatoric strain. Adaptive forces and observations use the
existing transpose/reconstruction maps, including homogeneous master clamps.

`eigenstrain_response` reports actual total, prescribed free, actual elastic
and constitutive elastic strains separately. Physical stresses and recoverable
elastic energy exclude ghost/Nitsche stabilization and thermal internal
energy. The original material envelope is applied to actual and projected
**elastic** strain; a separately supplied total Frobenius-strain budget in
`(0,1)` controls the declared small-deformation regime. This numerical budget
does not certify finite rotation, orientation, a linearization error bound or
an arbitrary callback between retained samples. Use the algebraic
`solve_controlled` followed by this physical response; the existing zero-free-
strain physical admission and stress queries do not describe thermal stress.

Prepared load actions and exact scale transposes never reevaluate the source
field. Stress VJPs differentiate the physical offset stress. Equilibrium
gradients include `z^T db_free/dscale - z^T dK/dscale u`; complete-load
compliance has **two** load-work terms. Native density and mass are independent
of this stiffness scaling. Geometry, free-strain and source-temperature
derivatives, thermal back-reaction and continuum certification remain outside
these fixed-field APIs. Point allowances and cancellation refuse unpublished
preparations/observations without changing the original operator.

Focused G0/G1/G3 regressions in `eigenstrain3`, `eigenstrain3_sensitivity` and
`adaptive_eigenstrain3` cover compatible nonuniform expansion, restrained
compression and storage, physical-cut B-bar averages, full coupled tensors,
Nitsche boundary signs, unchanged material inventory, distinct strain domains,
hanging-node reconstruction and independent finite-difference adjoints.

### Thermal embedded-support forces and moments

Cartesian and adaptive `eigenstrain_embedded_reaction` extend the original
Nitsche reaction by subtracting the retained `C epsilon_free n` traction on
exactly the selected embedded support points. Constant modes select forces;
rotational modes select moments about the caller's origin. This includes the
original scaled penalty and prescribed-motion terms. Additional strong box
supports remain outside this embedded reaction functional.

The exact displacement derivative is the existing reaction derivative; the
thermal correction changes the value and direct material-scale derivative.
Adaptive restriction is applied exactly once. Equilibrium gradients must also
include thermal-load and prescribed-motion derivatives before subtracting the
stiffness contraction. Physical thermal strain admission belongs to
`eigenstrain_response`; this reaction API evaluates a numerical functional.

Preparation retains surface eigenstress alongside bulk strains and load
columns, within the same original point budget. Repeated reaction queries
consume that cache without reevaluating the prescribed free-strain field.
Ownership, point alignment, finite arithmetic and point-level cancellation are
checked. `eigenstrain3_reaction` covers free-expansion force/moment cancellation,
analytic confined axial compression, adaptive partial and equilibrium
derivatives, source retention, invalid input and cancellation.

### Complete thermal strong-support force balance

Cartesian and adaptive `eigenstrain_equilibrium_load` retain the original
supported entries of both external loads and the prepared free-strain columns,
including any embedded prescribed-motion lifting. Reassemble after changing
material scales. Use the same complete vector for the algebraic solve and
`eigenstrain_equilibrium_response`; a masked thermal load cannot recover
the forces already discarded at its strong supports.

Thermal force recovery evaluates the original uneliminated cell, B-bar,
Nitsche and ghost actions, subtracts the complete RHS and reports strong-support
force and moment plus the free residual norm. Adaptive forces use the original
unmasked hanging-node transpose and master-node moment arms. Embedded support
reactions remain a separate functional, avoiding double counting.

Physical admission shares the thermal response owner: actual/projected total
strain obeys the explicit kinematic ceiling and actual/projected elastic
strain obeys the original material calibration at retained bulk points.
The admission-only path allocates no array of stress observations; adaptive
reconstruction remains cell-local. A solved state is needed to interpret
supported residuals as equilibrium reactions. The free norm measures algebraic
imbalance and does not certify continuum accuracy.

`eigenstrain3_equilibrium` covers compatible expansion beyond narrow elastic
calibration, restrained thermal face forces, supported source-density body
loads, graded B-bar force/moment balance and simultaneous strong/embedded
support reactions, together with ownership, budget and cancellation refusal.

## Consistent native material inertia in three dimensions

`CutElasticity3::prepare_material_mass` and the adaptive equivalent prepare
`MaterialMass3` from the SAME retained cut-volume quadrature, Q1 basis and
native source density as the elasticity/body-load owner. Each cell contributes
`occupancy[cell] * integral rho0 N^T N dV`. Stiffness scales, ghost penalties
and the SIMP ersatz modulus do not contribute inertia. Legacy elasticity
without source density refuses preparation.

The prepared cache stores symmetric scalar node pairs once and shares them
across three components. Adaptive preparation applies the existing hanging-node
map as `T^T M T`, with homogeneous strong clamps removed. Allocation-free
`apply_into`, the consistent diagonal and per-cell bilinear occupancy pullbacks
share this cache. An explicit scalar-entry cap bounds preparation; cancellation
is polled during quadrature, constraint expansion, application and pullback.
Errors during an in-place action require discarding the partial output.

Geometry/material ownership remains bound across stiffness updates and moves;
a separately rebuilt grid must prepare its own cache. Source states remain
available per cell. Coefficients are kilograms for SI geometry and density.
The matrix is the numerical consistent mass of the retained quadrature,
without row-sum lumping, stabilization inertia or artificial support masses.
Zero occupancy has zero mass. A positive diagonal does not prove positive
definiteness; singular cut or void configurations remain the consumer's
admission responsibility. This API alone claims neither a modal solution nor
a continuum mass/error certificate.

## Three-dimensional bulk stress observations

`elastic3::stress` evaluates the retained Cartesian or adaptive Q1 bulk
quadrature points with the same retained reference tensor as stiffness assembly.
Each `BulkStressPoint3` names its active cell, position, positive reference
volume weight, reference-material stress and current scaled physical stress.
Voigt order is `[xx, yy, zz, xy, yz, xz]` with physical shear components.
Ghost/Nitsche stabilization influences the solved displacement but is not
reported as material stress.

`bulk_stress_pullback` accepts one six-component cotangent per original point,
without inserting quadrature weights, and returns the exact partial derivatives
with respect to displacement and current cell scales. The separate
`reference_bulk_stress_pullback` differentiates the unscaled reference-material
stress. An equilibrium derivative additionally needs the solve adjoint and
`scale_bilinear_forms`, whose canonical bulk/Nitsche/ghost contraction is shared
by both backgrounds. Adaptive observations use the original hanging-node map
and its transpose; there is no nodal stress recovery, remeshing or new assembly.

### Three-dimensional strain, stress and elastic storage

`CutElasticity3::bulk_strain_stress` and the adaptive counterpart return
`BulkStrainStressPoint3` on the retained physical bulk quadrature. Each sample
contains the original stress response, actual displacement strain,
constitutive strain and elastic storage per reference volume. Strains use
tensor shear in `[xx, yy, zz, xy, yz, zx]` order. For B-bar the actual strain
retains local dilatation; only the constitutive strain uses the same
physical-cell mean dilatation as the operator. Storage is one half of their
work-conjugate stress/constitutive-strain contraction, including both shear
entries and current material scaling. It excludes stabilization and boundary
work. Both strain families obey the original material domain: a Frobenius
ball for native isotropic/orthotropic laws, including their source-state
adapters, or an explicitly declared Mandel-component box.

Every Cartesian stress/strain observation and its stress transpose removes
a common nodal translation before differentiation. Adaptive strain recovery
reconstructs only the current cell through the original hanging-node map and
strong-clamp mask; no global displacement copy or new assembly is required.
The point allowance precedes output allocation, with finite-data admission
and cancellation at input, cell, point and publication boundaries.

The `stress3` and `elastic3_bbar_stress` regressions check affine tensor shear,
material scaling, hanging fields, distinction of actual/assumed B-bar strain,
and integrated constitutive storage. These are small-strain numerical samples;
they do not bound unsampled extrema or certify continuum error.

`affine_bulk_stress` samples one offset field and a family of uncertain
directions on those same points. It shares native Mandel-strain extraction and
constitutive action with ordinary stress observations, and stores reference
stress columns in original field/point order. Coordinate columns are not
independent physical load cases. For an original Mandel-component box it checks
`abs(epsilon0) + ||epsilon_basis||_2` against each retained component ceiling,
which admits the complete coefficient unit ball in real arithmetic. The scaled
norm avoids squaring large or tiny components; computed admission is numerical,
not an interval enclosure. Native isotropic and orthotropic constructors retain
the supplied tensor Frobenius strain ceiling, including both symmetric shear
entries, even after conversion through `ElasticAcousticState`. Explicit bulk/shear,
cubic and generic Mandel constructors keep their declared component boxes.
Rigid translation is subtracted before the strain check; stiffness interpolation
never widens either domain. For a tensor-norm affine family, the owner normalizes
by that original ceiling
and accumulates a six-dimensional Gram matrix and offset coupling with constant
storage. A rotation-invariant spectral-moment bound admits the complete
coefficient ball without an eigensolve. Admission is exact in real arithmetic
for a single direction; multiple-direction affine families may conservatively
refuse. No material limit is applied to algebraic unit basis or adjoint vectors.
The checked `point_count * (basis_count+1)` allowance precedes output allocation
and adaptive field reconstruction. The adaptive implementation preserves the
same hanging-node map and strong-clamp masking for the entire family.

`admit_material_domain`, `admit_affine_material_domain` and
`admit_load_space_material_domain` expose that same material-strain admission
independently of stress output. They check shape, finite coefficients and scales
before any legacy/source distinction, then stream the eight corners of every
cell with constant auxiliary storage and no heap allocation. Q1
strain components, including B-bar's projected strains, are multiaffine in the
cell coordinates. Their coefficient-ball support, under either the component
or tensor norm, is convex in each coordinate separately, so its maximum lies
at a vertex. Both actual and projected B-bar strains are checked.
The bound covers the complete discrete Q1 field, including the fictitious
extension of cut cells; that extension can make refusal conservative. The
zero-offset load-space method uses an implicit local zero rather than allocating
an all-zero displacement vector. Adaptive fields reconstruct one cell at a time through the
identical scalar row owner used by full physical displacement export. This lets
compliance, response and optimization consumers enforce native validity without
allocating point-by-field stress arrays or performing constitutive stress work.
The whole-ball check retains bounded cancellation and leaves operator state
unchanged on success, refusal or cancellation. Arithmetic remains numerical,
not an interval enclosure or a continuum-error certificate. Stress observation
methods still return values only at their retained integration points.

`CutElasticity3::solve_physical_controlled` and
`AdaptiveElasticity3::solve_physical_controlled` run the existing algebraic
solver and then enforce the complete material domain before returning a
physical displacement. A material refusal or cancellation returns no iterate;
admission checkpoints report the final solver iteration count. The original
`solve_controlled` stays available for algebraic solves, adjoints and response
directions whose amplitudes do not represent physical loading. A true residual
alone is not proof of small-strain material validity.

Output point allowances are checked before allocating the stress array. State,
cotangents and arithmetic must be finite, and observation/pullback work polls
at each point and before publication. Errors return no partial result and leave
operator scales unchanged. The `stress3` tests cover analytical affine and
hydrostatic fields, rigid rotation, state/material directional derivatives,
constrained transpose identities, finite budgets and cancellation. These are
discrete stress samples and derivatives, not continuum peak-stress bounds,
material-failure guarantees or geometry derivatives.

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
16. GRADED VECTOR ACCURACY (cte-006/007, G0/G1/G3): a fixed half-domain
    2:1 refinement pattern retains mixed active levels across a three-level
    manufactured-solution ladder and gates fitted Q1 orders near L2 = 2 and
    H1 = 1. A separate affine dead-load fixture compares uniform and graded
    physical fields and algebraic compliance against the same analytic
    external work, then requires bit-identical graded solve replay.
17. SHARED GHOST TOPOLOGY: exact dyadic shared-face patches cover every
    balanced coarse/fine face without gaps or overlaps; reverse queries return
    the same positive-axis patch, and vector assembly/DWR consume that one
    geometry source. Mixed-level affine jumps vanish, while non-affine ghost
    energy remains symmetric and non-negative.

## Error model

`CutFemError` teaching errors: `EmptyDomain` (level set never enters
the grid), `CutBandNotUniform` (the scalar frontend's enabled ghost faces need
an equal-level interface band; names the offending pair and the repair;
vector elasticity instead integrates exact balanced 2:1 patches and
ghost-free aggregation does not impose this precondition),
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
(configured residual gate missed; carries iterations and the recomputed
Euclidean residual).
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
uniform-fine error at < 0.55× its DOFs); cut-009 canonical scalar
sampler fails closed (linear field reproduced on active leaves,
certified-outside classification, missing/non-finite corner evidence
and out-of-box/rim/non-finite points all refuse as `InvalidFemInput`).
Grid unit tests: 2:1 balance, leaf partition of unity.

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
refusals. The scalar graded-aggregation and vector graded-patch cases also pin
typed true-Euclidean residual provenance; the scalar case independently
recomputes `b - Ax` and requires bit-identical public evidence.
`tests/graded_elasticity.rs`: cte-006 fixed-pattern mixed-level MMS with
three-level fitted L2/H1 slopes; cte-007 uniform-versus-graded affine physical
field and compliance equivalence plus bit-identical graded replay.
`tests/elasticity_adjoint.rs`
is a declared `adjoint-vjp`-required target: cte-004 exact registered
transpose, two-operator key isolation, and independent central-FD gate.
The acceptance tests emit deterministic JSON verdict rows suitable for
capture by the proof lane; stdout alone is not claimed as retained
ledger evidence.

## No-claim boundaries

- 3D octree instantiation and 3-D bulk/surface quadrature: the narrow
  height-line root primitive ships, but no octree, tensor base rule,
  root-count subdivision, volume/surface enclosure, or 3-D ghost penalty
  path exists here yet.
- Higher-order elements (Qk, k ≥ 2) and the matching higher-order
  ghost penalties (derivative jumps beyond first order).
- Moment-fitting quadrature (tessellation-based scheme ships; the
  moment-fitting alternative is future work under the same `CutRules`
  surface).
- DWR estimation and refinement policy remain consumer responsibilities:
  fs-dwr supplies compliance indicators and fs-topopt owns
  `refine_dwr_cut_band`. Their feature-gated integration test composes two
  authentic estimate/refine/graded-re-solve cycles; this crate does not grow a
  second orchestration policy.
- Shared vector ghost patches are limited to balanced 2:1 Q1 quadtrees in 2D.
  Level jumps greater than one, incomplete/overlapping face coverage, 3D face
  geometry, and higher-order normal-derivative jumps refuse or remain outside
  the claimed surface.
- Vector constitutive scope is two-dimensional isotropic small-strain
  elasticity. Ordinary entries use plane strain; explicit plane-stress
  entries condense transverse stress with the original E and nu. Orthotropy, nonlinear
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

## Traction-free-thickness plane-stress elasticity

`CutElasticity::{assemble,solve}_plane_stress` and their checked
`_with_boundary_traction` entries use the same bulk, symmetric Nitsche, ghost,
2:1 hanging-node elimination, CSR operator, adjoint apply and residual-checked
solve. They analytically eliminate transverse stress before evaluating the
3-D first Lamé coefficient: `lambda_plane = E*nu/(1-nu^2)`, with the original
shear modulus. Both boundary consistency terms and stabilization use this
reduced law. Operator and solution retain the explicit `ElasticPlane` choice.
The original plane-strain arithmetic and material gate remain unchanged.

For thin in-plane bodies, stiffness stays finite as the original 3-D nu
approaches 1/2. This does not extend the admissible nearly incompressible
plane-strain model. All forces and energies use unit reference thickness;
there is no bending, variable-thickness, anisotropic or finite-strain claim.
Reference: A. F. Bower, *Applied Mechanics of Solids*, section 3.2.3,
https://solidmechanics.org/Text/Chapter3_2/Chapter3_2.php.

`tests/plane_stress.rs` checks independent affine strain energy on uniform
and graded cut meshes, auxetic and near-incompressible thin materials,
three-rung curved-domain manufactured L2/H1 displacement orders, checked
traction resultants, controlled-solve parity, and identical original operator
bits at nu=0. These are numerical model checks, not experimental qualification.

### Original-material strain and stress recovery

`elastic::strain::ElasticStrainStress` observes the original Q1 displacement
field through both `CutElasticitySolution::strain_stress` and
`CutElasticityOperator::strain_stress`. The latter also accepts the original
controlled solver's coefficients and reconstructs only the four required nodes,
including hanging-node expansions, without a heap allocation. Operator and
solution retain the original material, plane reduction and node lattice.
Observation refuses a different node lattice instead of silently relabeling
the retained displacement field. Gradients subtract one nodal translation
before differentiation so rigid motion does not generate cancellation strain.

Strain and Cauchy stress use fs-material's tensor-shear order
`[xx, yy, zz, xy, yz, zx]`. Plane strain retains transverse stress; plane stress
retains thickness contraction `epsilon_zz=-nu*(epsilon_xx+epsilon_yy)/(1-nu)`
and exactly zero transverse stress. Stress reuses the assembly's reduced Lamé
coefficients, avoiding cancellation of divergent 3D bulk terms near nu=1/2.
Physical bulk storage is `sigma:epsilon/2`, without artificial ghost/Nitsche
terms or prescribed-boundary work. The complete Frobenius strain norm, including
thickness and both symmetric shear entries, must fit the original material's
declared strain ceiling; algebraic convergence alone is not that admission.

`CutElasticityOperator::admit_strain_domain` checks all four corners of each
active cell. The Q1 gradient and the reconstructed thickness strain are affine,
so convexity of the tensor norm bounds the entire cell, including strain peaks
between quadrature sites. This is a floating-point bound on the discrete field,
not an interval or continuum-error certificate. It includes the fictitious
extension of cut cells and can therefore refuse conservatively. Point queries
do not classify physical SDF membership. Cancellation is checked per cell and
before returning a complete bound. All observation methods keep assembly and
solver arithmetic unchanged apart from translation-stable gradient recovery.

Ordinary and controlled physical solves apply that same whole-field admission
after the true residual gate and before publishing coefficients or observations.
The controlled path polls cancellation per cell during this final check. It
uses the operator's retained cell keys and node lattice, without rebuilding a
tree or allocating a second field. Plane-stress thickness contraction and the
cut-cell extension remain part of the bound; a domain exit is a physical-model
refusal even when the linear equations have converged. Algebraic operator and
adjoint actions remain available without treating their vectors as physical
displacements.

`tests/strain_recovery.rs` checks analytic extension/shear and energy against
the original stiffness on uniform/graded grids, auxetic and nearly incompressible
plane-stress materials, rigid translations/infinitesimal rotations, missed
quadrature-site peaks, thickness/shear strain limits, controlled/ordinary solve
parity, original lattice ownership and transactional cancellation. These are
G0/G1/G3 checks of the stated small-strain model, not finite-rotation, plasticity
or experimental validation.


## Signed thermal initial-stress stiffness

`CutElasticity3::prepare_eigenstrain_prestress` and its adaptive counterpart freeze
`G_ai,bj = delta_ij integral grad(N_a)^T sigma_0 grad(N_b)` on the original cut rules.
The stress is recovered from actual minus prescribed free strain before applying
the retained tensor, with the original B-bar projection and current scales.
Selected weak supports include both symmetric geometric Nitsche flux terms. The
original elastic penalty and ghost coefficients stay frozen in K_e; no extra
geometric penalty, ghost stress or artificial inertia is invented.

`PreparedPrestress3` condenses strong/hanging constraints once and shares scalar
upper-triangular node pairs across xyz. Actions, signed diagonals, bilinear work
and read-only sparse entries reuse that cache without new source evaluations.
It binds exact geometry/material/support coordinates and stiffness-scale bits;
changed scales or a rebuilt operator require preparation again. Bulk-plus-surface
point work and retained/current-cell pair storage are bounded and cancellable.
Both bulk observations and prestress now share one physical response kernel;
selected surface free strain is retained to avoid subtracting large stresses.

`eigenstrain3_prestress` adds seven regression sources for free expansion, tension
and compression, coupled tensors, nonuniform/B-bar correlation, adaptive virtual
work, geometric Nitsche flux, surface admission and refusal/rollback boundaries.
This is conservative reference-geometry initial stress at an admitted supplied
state. The caller must separately establish equilibrium; finite-strain dP/dF,
follower loads, coercivity and complete buckling/stability certification are not
inferred. In particular, negative geometric stiffness is retained, not clipped.


### Sparse handoff of original consistent inertia

`MaterialMass3::reference_scalar_entries` exposes an allocation-free read-only
stream of `(cell, upper node pair, reference mass)` from the original prepared
owner. Adaptive condensation and strong-row removal are already included. Sparse
consumers multiply physical occupancy, replicate xyz and sum/mirror pairs; they
must retain native source validation. No stiffness scale or mass lumping enters
this handoff. The signed CutFEM modal adapter uses it to factor only physical M.


## Retained-geometry native material refresh

Cartesian and adaptive 3-D elasticity can replace their complete native material
field using `refresh_material_field`. Bulk matrices, original B-bar modes,
selected symmetric Nitsche consistency/penalty terms, and both ghost sides are
staged using retained quadrature, support selections and face factors. Geometry
and adaptive constraint maps are borrowed; no SDF, patch or clamp is resampled.
Reference-density bits and physical reference inventory remain fixed, as do
current stiffness-design scales. Legacy states without native density refuse.

Explicit cell, retained-point and fixed native staging-payload limits precede
material callbacks. Dynamic source receipts, query allocations, allocator
overhead and kernel stack are outside the staging-payload bound. Returned errors
and cancellation preserve the original operator and its preparations. Success
rotates the common owner token: mass, body-load, eigenstrain and prestress caches
must be prepared again. This operation commits an operator refresh only; a later
equilibrium or thermal step has its own failure/commit semantics.


## Native isothermal thermal-source fields

`ElasticMaterialField3::ThermoelasticUniform` and `ThermoelasticByCell` retain
one complete native source response per active cell. The original material owner
supplies C(T), integrated free strain, reference density, identity, domains and
thermal derivatives. Temperature and constitutive data are cell-constant;
quadrature-level temperature variation within a cell is not implied. Explicit
isotropy is preserved independently of thermal-expansion orientation.

`prepare_thermoelastic` uses the existing free-strain bulk and negative Nitsche
traction loads with that same C(T). Its result feeds the existing equilibrium,
dynamics, design, stress and prestress APIs. Successful material refresh requires
new preparation. B-bar averages trace differences from a retained sample, so a
constant free dilatation survives normalization without an artificial offset.
Shared physical recovery enforces the thermal source's TOTAL and ELASTIC Mandel
component domains in addition to the caller's total Frobenius strain budget.

`thermoelastic_response` evaluates the original native potential at the actual
constitutive total strain and physical reference-volume weight; cell stiffness
scales multiply extensive elastic free energy, entropy, internal energy and
capacity exactly once. Capacity remains signed. The source's local mixed stress
partial is returned separately from any solved-equilibrium derivative. These
are physical bulk quantities: they exclude caloric charts, transport, numerical
ghost/support storage and whole-system heat-step balance. Adaptive recovery uses
the original homogeneous hanging reconstruction cell by cell. Allocation-free
admission follows the same physical checks. Counts bound retained bulk output;
cancellation returns no candidate observations.


## Prescribed-temperature residual transpose

Cartesian and adaptive `thermal_residual_temperature_pullback` return each
active cell's `z^T[(dK/dT)u - db_free/dT]` at fixed geometry, reference density,
stiffness scales and external load. Native C', integrated expansion slopes,
B-bar K'/G', both symmetric Nitsche consistency terms, the penalty modulus and
both ghost-side coefficients are included. The physical displacement passes
native total/elastic admission at retained bulk and selected surface points;
the finite algebraic dual has no physical strain ceiling.

Bulk, selected support and ghost point counts precede O(cells) scratch/output
allocation. The original adaptive reconstruction is applied locally to both
fields; no global expansion, derivative matrix, material fit, source callback
or finite-difference production kernel is introduced. Relevant nonsmooth
absolute-value/max-row stabilization branches refuse a derivative. The caller owns the
equilibrium-adjoint sign and explicit objective partials.

The Cartesian and adaptive `thermal_residual_temperature_pullback_with_motion`
entries additionally differentiate the original embedded Nitsche motion load.
The same prescribed displacement used by `prescribed_displacement_load` enters
both the constitutive consistency and penalty derivatives. An optional
displacement slope in metres per kelvin adds the full lifting product rule
with respect to the owning cell's temperature. Geometry and strong zero box
clamps remain fixed; the motion must be compatible with those clamps. Missing
embedded support or nonfinite callback data refuse. The original point budget
precedes callbacks, cancellation is checked immediately after each callback,
and no global lifting vector or derivative matrix is allocated. Existing
preparations and the operator remain unchanged on success or failure.
Temperature-dependent external forces and a heat-equation adjoint retain
their separate owners and terms.


Explicit-isotropic free-strain and signed temperature-stress actions share the
material owner's compensated bulk/deviator split. The same original native K/G
and their derivatives now drive the source thermal potential, including
independently oriented anisotropic expansion. This avoids recovering a weak
mode through cancellation between rounded dense tensor entries or a rounded
mean. Stored C(T) matrices are retained, generic orthotropy is not fitted, and
source domains remain mandatory. This improves numerical mode recovery; it does
not establish conditioning or convergence for arbitrary cut-grid solves.

Explicit isotropy also follows the native mode action in displacement-formulation
free-strain recovery, cold bulk observations and affine stress families. Stress
transposes apply the original constitutive modes to the physical-stress cotangent
before scattering nodal forces. B-bar retains separate local deviatoric and mean
pressure work; physical shear normalization and the rigid-translation transpose
remain unchanged. Generic tensors retain their complete source matrix and are
never fitted to isotropy. Extreme independent K/G regressions cover both weak
bulk and weak shear response and adjoint work, without asserting that the dense
assembled stiffness can resolve arbitrarily ill-conditioned global solves.


## Complete thermoelastic objective adjoints

Cartesian and adaptive `solve_thermoelastic_adjoint` solve the native source-owned
thermal equilibrium, admit its physical strain domain, and evaluate one supplied
objective of displacement, scaled bulk stress and cell temperature. Its declared
partials feed the original reference-stress transpose and the SAME native elastic
solver for the adjoint. The total temperature derivative is the explicit objective
partial plus local stress-temperature work minus the native residual-temperature
transpose. Original C(T), integrated expansion, B-bar, both ghost-side coefficients,
Nitsche consistency/penalties and optional temperature-dependent support motion
therefore enter once. No Krylov iteration is differentiated.

Both solves must pass their original true-residual gates within one aggregate
iteration allowance. The complete temperature point budget is admitted before
preparation, solving or objective callbacks. Adaptive reconstruction and restriction
remain with their original owners. The finite algebraic dual has no physical
strain-amplitude ceiling. Failed, cancelled and successful calls leave the original
operator and scales unchanged; no partial sensitivity is returned.

Force derivatives hold the authored objective and temperatures fixed. Geometry,
reference density, stiffness scales, strong box clamps and source selection are
fixed. Explicit force dependence, temperature-dependent external forces, heat-flow
adjoints and derivatives of the caller's constitutive/source-selection policy remain
separate. Stress cotangents include caller-specified quadrature weights and pair
with physical tensor shear components. Sampled stress functionals do not certify
continuum maxima; the original supported problem must be well posed. G0/G1/G3
regressions extend `fs-solid/tests/cutfem_thermal_sensitivity.rs` with nonlinear
stress/displacement objectives, physical source refresh, original hanging maps,
moving supports, shared iteration bounds and objective-callback cancellation.


## Native thermoelastic stress objective partials

Cartesian and adaptive `thermoelastic_stress_pullback` return independent-coordinate
displacement and per-cell fixed-state stiffness-scale and temperature partials.
One physical `[xx, yy, zz, xy, yz, zx]` cotangent per retained bulk point supplies
its own integration weight and any objective shear multiplicity. The native
potential supplies C'(T) and integrated, independently oriented expansion
derivatives. Actual/projected total and elastic domains remain enforced;
algebraic cotangents have no physical strain ceiling.

Streaming recovery retains only point cotangents and cell partials. The existing
reference stress and constraint transposes retain native material modes, B-bar,
hanging nodes and homogeneous clamps. For K^T z = dJ/du, solved temperature
gradients subtract the complete residual temperature VJP from the direct
temperature partials. Geometry, reference density and external-load derivatives
retain their separate owners. The focused `cutfem_thermal_stress_sensitivity`
fixtures cover Cartesian/adaptive fixed-state and solved stress gradients,
rotated orthotropic stiffness, separately oriented expansion and refusal
boundaries; no continuum maximum-stress or heat-step gradient is inferred.


## Prepared thermal residual Jacobian

`prepare_thermal_residual_jacobian` and its `_with_motion` entry freeze the
mixed block d(K(T)u-b_free(T)-b_g(T)-f)/dT at one admitted native physical
displacement. Temperatures are independent per retained cell. Geometry,
reference mass, cell scales and external forces remain fixed; optional pure
embedded-motion callbacks supply g and its local temperature slope. Strong
box supports remain homogeneous.

Preparation integrates one direct 24-entry force column per cell and retains
each original ghost jump and its two signed temperature weights. Bulk C'(elastic
strain)-C alpha, B-bar projection, both symmetric Nitsche consistency terms,
penalty derivatives and the complete motion product rule use the original
source/material owners. Native explicit isotropy uses compensated K'/G' for
both displacement and B-bar; only B-bar changes the ghost modulus to G'.

`ThermalJacobianLimits3` admits bulk/selected-surface/ghost point counts, peak
cached/scratch scalar storage and original stencil work before allocating or
calling boundary functions. Pointwise physical strain and smooth stabilization
branches are admitted during preparation. A failed preparation publishes no
state. Later forward and transpose `_into` actions allocate no scratch, query
no sources, integrate no quadrature, and share exactly the cached coefficients.
The original adaptive T/T^T and homogeneous strong mask remain authoritative.

Binding identity, coordinate space and exact scale bits reject stale material,
geometry or scale use before output writes. Changed physical displacement or
boundary data requires explicit new preparation. Directions and duals are
algebraic; no physical amplitude ceiling applies. Mid-action cancellation or
arithmetic overflow can leave partial output that the caller must discard.
This block supports coupled linearization and implicit adjoints; it does not
itself solve equilibrium, advance heat, or establish whole-system passivity.
