# fs-topopt CONTRACT

## Purpose and layer

Layer: **L4 ASCENT** (deps: fs-adjoint/fs-solver/fs-feec L3,
fs-ascent L4, fs-material/fs-solid L3, fs-la/fs-sparse L1, fs-math L0,
fs-rep-mesh L2). Density-based topology optimization (plan §9.5 [S]):
SIMP with the modern hygiene stack — Helmholtz PDE filtering,
Heaviside projection with continuation, exact chain-rule
sensitivities, and the classical optimality-criteria driver.

NAMING: the plan's atlas used "fs-topo" for this stack; that crate
name carries the L2 topology-CERTIFICATE machinery (persistence,
cubical homology), so the optimization stack lives here.

## Public types and semantics

- `DensityElasticity` — matrix-free K(ρ̄) = Σ_c E_c·K_c with per-cell
  reference-material 12×12 blocks (native solid tensors or the legacy
  fs-material isotropic tangent × fs-feec barycentric-gradient B-matrices) kept separate for the
  exact chain rule; Dirichlet dofs handled by identity-on-fixed
  masking (SPD on the full vector space); `cell_energies` =
  uᵀK_c u per cell (the compliance sensitivity kernel).
- `DensityFilter` — Helmholtz filter: volume-weighted cell→vertex
  scatter, (M + r²K)⁻¹M solve on the FULL vertex space (natural BCs —
  the correct filter behavior: no boundary droop), vertex→cell
  gather. Linear; `apply_transpose` is the exact chain-rule pullback
  (adjointness ⟨Fx, w⟩ = ⟨x, Fᵀw⟩ G0-tested; constants preserved to
  solver tolerance). One assembled operator, built once.
- `heaviside`/`heaviside_derivative` — tanh projection with β/η,
  exact endpoints, monotone, closed-form slope; tanh through the
  strict exp kernel (no platform libm in the pipeline).
- `DesignPipeline` — ρ → filter → projection → SIMP
  (E_min + (1−E_min)·ρ̄^p) with `pullback` reversing the chain
  exactly, and `compliance_and_gradient` exploiting self-adjointness
  (λ = u ⇒ dc/dE_c = −u_cᵀK_cu_c: ZERO extra solves — stated and
  FD-verified).
- `optimality_criteria` — the classical OC driver (documented choice
  for compliance/volume; fs-ascent's augmented Lagrangian is the
  general path): multiplicative update with move limits, volume
  multiplier by fixed 80-step bisection on the PROJECTED volume —
  fully deterministic, whole runs replay bitwise.
- `eigenfreq` — generalized eigenproblems K(ρ̄)φ = λM(ρ̄)φ (per-cell
  consistent P1 mass blocks; Cholesky reduction + Jacobi at fixture
  scale), exact eigenvalue design gradients through the full chain,
  the MASS-INTERPOLATION trap handled (linear mass above ρ = 0.1,
  continuously-matched ρ⁶ below — spurious void modes gated by
  regression), smooth-min aggregation with exact weighted gradients
  for CLUSTERED eigenvalues (FD-verified near a designed crossing).
- `stress` — relaxed von Mises measures with qp-RELAXATION (exponent
  q < penal; void cells cannot drive the constraint — share and
  floor-stability gated), p-norm aggregation with the ADAPTIVE
  normalization c = σ_max/PN reported per evaluation, and the
  NON-self-adjoint design gradient (one extra adjoint solve via
  `stress_pullback`, cross terms by the polarization identity on
  cell energies).
- `robust::{RobustPipeline, robust_optimality_criteria}` — the
  erode/dilate three-field formulation: one filter, three projections
  (η ± δ), POINTWISE-ORDERED realizations (tested); minimize ERODED
  compliance s.t. volume on the DILATED field whose target is adapted
  (DAMPED, 0.3 blend) so the NOMINAL design carries the budget — the
  undamped adaptation measured a period-2 limit cycle on cold starts
  (kept as a regression probe); reports carry the erosion-retention
  ratio vol(eroded)/vol(nominal), the measured minimum-length-scale
  signal.

- `marquee` module (bead b7d0; [F], behind `cutfem-marquee`): the
  CutFEM-quadtree topology marquee, the 2D analogue of the planned octree
  lane. `DensityDesign` (nodal densities;
  the SOLID region {ρ > ½} IS the CutFEM domain via the bilinear
  CutSdf with an exact-containment enclosure, plus a deterministic
  non-cryptographic `state_witness` for trace correlation), `run_marquee` (the
  volume-to-point heat fixture: interface-flux redistribution with a
  BAND-LOCAL volume projection, DWR-gated band refinement, an explicit final
  re-solve on the returned design, and zero rebuilds structurally),
  `refine_dwr_cut_band` (an estimator-agnostic,
  one-step cut-band planning policy over `CellKey` indicators, returning
  a `DwrBandDecision` reason, versioned mass evidence, analysis headroom,
  target leaf delta, and the ACTUAL halo/balance split count; invalid levels,
  non-leaf keys, non-finite indicators/masses,
  and non-finite SDF enclosures return `InvalidFemInput` before
  caller-visible mutation), `DWR_CUT_BAND_POLICY_VERSION` and
  `DWR_CUT_BAND_MASS_GATE` (the explicit policy identity and strict 0.15
  gate). Since policy version 2 the gate measures mass on the region it
  refines, `in_dwr_refinement_band` (cut cells plus their one-cell halo), and
  the band must also hold strictly more than its area share of the mass.
  Version 1 counted only cut cells, but the error it had to find sits one
  cell off the interface. On the seeded fixture at level 4, cut cells held
  11% of the mass on 10% of the cells, and the marquee never refined. The
  band held 48% on 26%. `void_components`
  (topology witness), `min_feature_cells` (the medial-axis-class thickness oracle).
  The `run_marquee` refinement argument is an enable flag: each enabled
  iteration may advance the whole band by at most one level; it is not a
  requested split budget.
  HARD-WON INVARIANTS from development, all conformance-guarded:
  fs-cutfem's ghost penalty demands the cut band AND ITS ONE-CELL HALO
  at a uniform level (`halo_cut`), conformed after EVERY design update
  because the interface moves; the next `Space::build` verifies rather than
  silently repairs that invariant; interface membership is a NEIGHBOR SIGN
  CHANGE, not a |φ| threshold (φ is a density gap, not a distance);
  DWR mass is classified against the EXACT solved/estimated design
  snapshot, then carried by an exact-grid-generation receipt that may
  authorize one level for the post-update design's band — old indicators
  are never reclassified against a new geometry; every post-update boundary
  is then conformed to the already-authorized level, including the final
  iteration when no later pre-solve repair exists;
  flux probes project THROUGH the interface so void-side nodes read
  real flux; and the volume projection lives ON THE BAND — a global
  shift silently fills voids from the inside.

## Shared physical materials and mass design

`DensityElasticity::with_material_field` accepts the existing
`fs_solid::linear3::TetMaterialField`, owning a uniform state or an ordered
heterogeneous field. Full Mandel stiffness tensors, physical densities and
material-state identities remain attached to one operator. The shared element
assembly converts both axes of the Mandel tensor to engineering-shear coordinates,
so compliance, stress recovery and their sensitivities use the same anisotropic
reference stiffness. Native material admission is repeated because callers can
edit a source state's public density or identity before construction.

Dense and matrix-free frequency objectives use each assigned reference density
in both consistent mass and its design derivative. Low-level mass actions still
accept absolute physical density weights over geometric unit-density blocks.
The matrix-free objective's legacy scalar `reference_density` must remain 1.0
when a material field is assigned; it cannot replace or multiply that field.

`DesignPipeline::mass_and_gradient` and its controlled counterpart return
kilograms and the full filter/projection gradient. Reference cell masses are
original tetrahedral volume times material density, including supported cells.
Physical inventory uses projected occupancy linearly; SIMP stiffness and the
modal low-density sixth-power interpolation do not define mass inventory.
`volume_and_gradient` retains its dimensionless volume-fraction semantics.

`SimpComplianceGoal` supplies mass gradients and matrix-free Hessian products
through the same weighted-projection helpers as volume. With `newton-design`,
`SimpNewtonModel::mass` exposes an SI-kilogram physics node for a minimum-mass
objective or mass constraint. It shares the exact-point compliance/volume cache
and lifetime work account; material-bearing preparation reserves one additional
transpose filter. Robust goals use dilated mass with eroded compliance. Legacy
stiffness-only models refuse the kg view. Material assignment is prescribed:
these derivatives vary occupancy, not the selected material card or its state.

`material_elasticity`, `material_modal` and `material_mass` cover tensor work,
physical density scaling, full-chain derivatives, robust mass curvature and
live-IR mass/compliance composition. The raw-SDF path below uses the same
physical-inventory convention with its own retained cut-cell quadrature.

## Physical material validity on raw-SDF designs

Physical 3-D equilibrium fields must satisfy each retained source material's
absolute Mandel-strain ceiling on the native bulk quadrature, independently
of whether the caller requests stress measurements or strength constraints.
The common native admission path covers ordinary and robust independent-load
compliance, displacement/reaction fitting and limits, and freshly resolved
binary endpoints. Every returned experiment is admitted, including a case
with zero objective weight. Source-free legacy elastic models receive no
invented strain ceiling.

The zero-offset spectral load-space evaluator also admits the COMPLETE
coefficient unit ball before its spectrum. Direct binary endpoints and their
threshold scanner share the same complete-family admission. Gray and binary refinement paths
independently admit both coarse and fine combined families against the actual
native operators. The fine family is checked before the shared DWR owner releases
its material and hanging maps; no reconstruction of a substitute operator or
extra equilibrium solve is needed. Independent-scenario APIs retain their
separate-experiment semantics. The native zero-offset entry point avoids even
allocating an all-zero displacement field.

Response and compliance refinement re-admit retained coarse physical fields
and independently admit fresh fine physical fields before producing goal
changes. The shared control adapter preserves cancellation and budget stops;
a material-domain refusal cannot publish a partial evaluation. Existing scale
rollback contracts remain in force. Admission adds no equilibrium solve and
does not allocate a retained stress-point family.

Adjoint vectors, virtual modes, residuals and DWR difference fields are
algebraic quantities, not physical experiments, and do not inherit physical
strain ceilings. The same native tensor and hanging-node/strong-clamp maps
still govern their actions and pullbacks. These checks establish numerical
validity at retained quadrature points, not continuum strain bounds,
constitutive extrapolation, yielding or a nonlinear material model.

## Physical mass on raw-SDF designs

All three `Sdf3Elasticity` backends expose the original native reference cell
masses from `fs-cutfem`: Cartesian, adaptive and explicitly preconditioned
adaptive operators. `CutDensityStudy3::mass_and_gradient` integrates
`rho_material * cut_volume * projected_occupancy` in kilograms, with geometry in
metres, and pulls its exact discrete derivative through the existing projection
and graph-filter transpose. Physical solid regions contribute inventory with
zero local derivative; physical void contributes none. The normalized geometric
weights defining the filter and volume fraction remain unchanged. No unit-density
fallback, SIMP ersatz mass or elasticity solve is added by this query.

`MinimumMaterialOptions3::objective` selects `VolumeFraction` or
`MassKg { scale_kg }`. The existing response-constrained AL owner uses the chosen
value and derivative, and ranks feasible accepted designs by actual physical
inventory. The mass value and pullback use the same prepared design as the
response solves; only the necessary transpose-filter work is added and charged.
History reports geometric volume, actual mass and the normalized AL objective
separately. Normalization that overflows or erases a nonzero inventory/gradient
is refused. Material tensors and densities stay fixed during a design session.

Refits preserve the objective and require native mass availability. Restart
binds the objective kind/scale and the ordered reference-cell mass vector, then
re-solves retained endpoints before installing them. The minimum-material byte
codec is version 3; older bytes are refused explicitly. Source/executable and
callback identity retain the existing external binding requirement. Adaptive
minimum-material refinement carries the selected objective and normalization
into every probe and refit. Mass mode measures actual coarse/fine kg inventory
on the inherited physical density field and adds each original parent's
absolute kg correction, divided by `scale_kg`, to Dörfler marking. Its
`mass_tolerance` checks that normalized mass change in addition to the existing
response and geometric-volume tolerances. Separate kg, volume and objective
reports remain visible. These are numerical checks on one enriched grid, not
continuum bounds; inventory adds no filter, equilibrium or adjoint solve.

Stress, response, joint response/stress and minimum-material optimizer
checkpoints bind the optional prepared material body load by its cached native
fingerprint. Restoration rejects a changed presence, sparse occupancy-to-force
map, source-state identity or density before reserving endpoint evaluations or
running physics. The fingerprint uses the existing domain-separated hash owner
over fixed-width layout and exact coefficient bits. It does not replace the
external binding of geometry, stiffness, other callbacks or executable identity;
fresh endpoint and derivative checks still precede resumed optimization.

`material_sdf3` covers heterogeneous inventory and finite differences, protected
regions, cancellation, density-independent geometric filtering, adaptive
inventory and real response-constrained mass minimization. The
`material_mass_sdf3` example exposes the path with independently tunable density
and stiffness in two explicitly synthetic material regions. These additions do
not claim a continuum certificate or completion of the full marquee milestone.

## Native 3-D vibration and eigenfrequency sensitivities

`sdf3::modal::controlled_sdf3_eigenpairs` connects every sealed
`Sdf3Elasticity` backend (Cartesian, adaptive and reusable adaptive solve space)
to the existing matrix-free `fs-spectral::generalized` owner. The same native
elasticity tensor, retained cut quadrature, strong/hanging constraints and
prepared `MaterialMass3` define `K phi = lambda M phi`. Fixed identity rows
are excluded. One existing backend preconditioner is prepared per evaluation,
and every inverse action uses the shared true-residual linear-solve gate and
the caller's original cumulative cancellation/iteration control. The adapter
reuses scratch vectors for stiffness and allocation-free mass applications.

`CutDensityStudy3::evaluate_eigenfrequency` returns source-material modes, an
unnormalized smooth minimum of the selected low eigenvalues, physical volume
and kg inventory, and all raw-design derivatives at the SAME design. The
gradient contracts both native stiffness and native consistent-mass cell
derivatives before the shared physical-occupancy/SIMP/projection/filter
pullback. Physical kg inventory remains linear in occupancy. Modal inertia
uses the existing shared low-density interpolation; it is not substituted for
physical inventory.

An extra computed mode checks the outer gap of a truncated selected cluster.
An unresolved boundary or a variable cell exactly at the piecewise mass-law
switch withholds the derivative. Protected cells retain their original zero
local design slopes. Mode normalization is not a physical excitation amplitude,
so source strain ceilings for physical responses do not constrain eigenvectors.
Both successful evaluations and returned failures restore incoming stiffness
scales; work already consumed remains charged.

The lane models undamped, infinitesimal vibrations about the undeformed,
zero-prestress reference state. Gravity and other external loads do not create
geometric stiffness in this formulation. Native source density is required.
Zero-mass free rows refuse; a positive mass diagonal alone is not an SPD
certificate. The shared eigensolver's residual-qualified Ritz estimates and
the numerical cluster check confer neither spectral completeness nor a
continuum frequency/error guarantee. The `cutfem-marquee` feature gate remains.

`sdf3::modal::design::Sdf3EigenfrequencyDesign3` connects this same evaluator
to the existing projected augmented-Lagrangian owner. Its two goals maximize
the smooth minimum under a physical volume cap, or minimize physical volume/kg
under a squared-frequency floor. The fixed eigenvalue reference is a numerical
normalization in the former and the actual required floor in the latter.
A source density cannot be overridden by a separate modal density scalar.

Both objective and constraint pass the existing full-chain numerical gradient
gate before optimization. Known cluster-boundary or piecewise mass-law
nonsmooth trials may be rejected by the bounded line search; cancellation,
linear budgets and other spectral errors retain their original causes.
Run segments preserve the same optimizer state. The feasible archive checks
original physical volume or both the smooth minimum and lowest returned
eigenvalue independently of AL tolerance, and a fresh recheck resolves the
entire selected modal group. Incoming stiffness scales remain unchanged.
The example `material_modal_sdf3` exercises the complete minimum-kg journey
with explicitly authored heterogeneous material inputs.

## Native material body loads on raw-SDF designs

Prepare an acceleration load with the native Cartesian or adaptive operator's
`prepare_material_body_load`, then attach it to that same operator's study using
`CutDensityStudy3::with_material_body_load`. Each sparse cell column is
`B_c = integral(rho_material,c * N^T * acceleration) dV`; it uses the same native
material state, cut quadrature, clamps and hanging-node map as the mechanics.
Geometry is in metres, acceleration in m/s², and assembled force in newtons.
The stencil is bound to the exact immutable operator owner, not just matching
dimensions. Preparation is bounded and cancellable; subsequent actions reuse
the columns without recalling the acceleration law or rebuilding quadrature.

Every independent experiment receives its declared dead load plus
`B * projected_occupancy`. Physical solid/void overrides and eroded or binary
occupancy apply to both material inventory and body force. SIMP stiffness and
its ersatz floor never create body-force mass. Acceleration remains fixed in
the reference configuration; deformation-dependent follower loads are outside
this model.

Fixed-grid compliance, displacement/reaction fitting, response-constrained
mass design, sampled stress and finite-case compliance envelopes use this
force in their actual equilibria. Compliance includes the load derivative
`2 B^T u`; other solved responses include `B^T adjoint`, composed with their
existing direct terms and the shared projection/filter pullback. Binary
response and finite-case compliance checks apply the actual binary occupancy
to the stencil. Actual external work remains separate from numerical support
lifting and is not relabeled as actuator work.

Explicit response and minimum-material refits accept a candidate study with
its own freshly prepared body load. Both studies must have matching body-load
presence, and the caller must keep the reference acceleration law unchanged.
Existing DWR and automatic refinement factories refuse attached material body
loads without an acceleration recipe. For adaptive work, configure the study
with `with_material_acceleration(law, max_terms, checkpoint)` instead. It retains
one pure, fixed acceleration law and prepares a fresh native stencil on each
probe or candidate operator, using the same explicit term allowance. Preparation
does not recur during density iterations. Replacing the stencil with
`with_material_body_load` discards the old recipe rather than silently reusing
a potentially different acceleration on the next grid.

Compliance DWR accepts the actual coarse projected field through
`estimate_reference_compliance_enrichment_with_occupancy` (and its coarse-level
variant); legacy field-only calls cannot infer physical occupancy from SIMP.
Response/reaction and response-limit estimators already retain that field.
Probes inherit it directly through the checked parent map. Native `B*rho`
enters both equilibrium and the localized weak residual, while total-force
compliance retains its load-dependent goal transfer and reactions retain their
affine offsets. Joint response/stress assessment reuses those same fine primals.

Automatic response, minimum-mass, robust-compliance, finite compliance-cap and
continuation refitters rebuild the candidate body load before optimization and
promote it together with the complete study. Candidate raw controls still pass
through the rebuilt filter/projection; they are distinct from inherited probe
occupancy. Failures retain the source study. `material_body_refinement_sdf3`
adds regressions for actual fine-grid force balance, reaction/displacement DWR,
protected regions, kg inventory and selective adaptive promotion.

Continuous zero-centered spectral load-space envelopes still refuse a material
body-load offset. The affine uncertainty owner below includes that offset in
its actual equilibria. Finite lists of independent load cases are also
supported. No follower-load, continuum-error or globally optimal design claim
is added by adaptive body-load support.

Run:

```bash
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example material_mass_sdf3 -- 24 250000 4 1.5 250 1 0.0002
```

The example minimizes kg under an axial reaction floor and a reference-volume-average downward
displacement ceiling. Its synthetic two-region material field supplies both
stiffness and mass; downward acceleration defaults to **1 m/s²**, not standard
Earth gravity. Output distinguishes actual kg, geometric volume, each original
response limit, body-force external work in joules, and the optimizer's actual
stop. A completed update allowance does not establish convergence, a local
strain bound, experimental validity or continuum accuracy.

## 3-D affine load uncertainty and constrained material design

`sdf3::robust::envelope::affine_load_space` evaluates the continuous physical
load set `f = b_dead + B*rho_eroded + F*a`, `||a||_2 <= 1`, for one to sixteen
fixed uncertainty directions. Columns of `F` encode physical amplitudes and
correlations. One positive compliance reference normalizes the complete set.
The shared forward filter prepares eroded, nominal and dilated realizations
once. Eroded occupancy controls both stiffness and native material self-weight;
dilated occupancy controls the reported volume and optional kg inventory.
Erosion is the prescribed manufacturing realization, not an assertion that it
maximizes compliance over manufacturing thresholds when self-weight is present.

One native stiffness preparation serves a direct baseline solve and each
direct basis solve. Reciprocity checks precede construction of
`G = F^T K^-1 F / reference`, `h = F^T K^-1 b / reference` and
`c = b^T K^-1 b / reference`. The shared
`fs_ascent::quadratic_ball::maximize_quadratic_ball` then solves the small
quadratic maximization without additional elasticity solves. The returned
coefficients reconstruct actual physical force, displacement and work; one
operator application checks the combined equilibrium residual. The report
retains the branch, numerical primal/dual gap, reciprocity and reconstruction
evidence. Incoming scales are restored on success and failure.

Every affine grid first admits the COMPLETE physical displacement ball against
the retained source materials' Mandel strain ceilings. This native check runs
before the compliance oracle, response adjoints and gradient publication,
independently of whether sampled strength is requested. It uses the same
material domain and strong/hanging reconstruction as native stress evaluation
without stress output allocation, constitutive stress evaluation or additional
linear solves. A source-valid baseline and valid coordinate endpoints do not
excuse an invalid correlated direction. Coarse and enriched grids, optimizer
trials, explicit refits and automatic campaigns all use this common grid owner.
A failed domain check or cancellation returns no report and restores incoming
scales. Unsourced legacy laws gain no invented ceiling. Retained-point numerical
strain admission is not a continuum maximum or an experimental validity claim.

A numerically separated regular maximizer supplies the discrete derivative
`(-u^T K' u + 2 B^T u) / reference`, composed with the existing physical-region,
SIMP, projection and filter-transpose operations. Hard, flat and unresolved
near-hard cases retain their worst-load witness but expose no smooth gradient.
All normalization and finite-range failures are explicit refusals.

`affine_load_space::design::AffineLoadSpaceDesign3` connects this evaluation to
the existing projected augmented-Lagrangian state. It can minimize worst
compliance under a dilated-volume cap, or minimize dilated volume/physical kg
under a worst-compliance cap. The latter uses the policy reference as the
physical cap; kg uses the native material inventory and an explicit kg scale.
Initialization checks the actual optimizer objective and complete constraint
with the shared directional audit. Single-cap checks retain the unshifted
quantity during differencing to avoid cancellation against a constant one.
Nondifferentiable baseline designs refuse; such line-search trials are
rejected by the existing bounded optimizer. Accepted state, history, work and
the best physically feasible point are retained across run segments. Physical
feasibility uses the actual cap, independently of AL stopping tolerances, and
`recheck_best` recomputes the complete envelope without advancing the optimizer.
No monotonic feasibility restoration is assumed for design-dependent gravity.

`affine_load_space::response` applies the existing signed displacement and
embedded support-force/moment limit types to the same continuous affine ball.
For an original scalar response `y(a)=y0+v^T a`, its upper/lower signed residual
has discrete worst value `[sign*(y0-bound)+||v||_2]/scale`. Every observation has
its own maximizing coefficient vector. The compliance adversary is never used
as a substitute. Each report retains the baseline, basis responses, radius,
coefficients, original bounds and a native remeasurement of the combined field.
Both the support formula and actual witness must pass the physical bound.

The joint evaluator shares the original filter, physical realizations, native
stiffness preparation and offset/basis primal fields. Response gradients use
the native observation pullback, embedded-reaction direct material term,
stiffness contractions and `B^T z`, then the same physical density pullback.
Few limits use one adjoint per observation; larger families factor the aggregate
by the affine primal fields and use at most `dimension+1` adjoints. Only the
offset field receives the material body-load derivative in that factorization.
The report identifies the retained adjoint layout. Measurement-only assessments
reuse the same support/witness calculations without adjoints or extra primals.

`AffineLoadSpaceDesign3::new_with_response_limits` retains these inequalities
through every gradient probe, optimizer trial, feasible archive update and
fresh recheck. The original compliance or volume constraint and response
aggregate are combined with the existing unnormalized log-sum-exp upper
aggregate. Its weights supply the exact discrete chain rule; no averaged
constraint can hide an original violation. The existing optimizer supports
both design goals, including dilated native kilograms. A zero response radius
retains valid discrete values but exposes no smooth norm derivative, so smooth
baseline admission refuses and such trials are rejected. Native material
scales are restored on all returned outcomes, including late component stops.

`affine_load_space::constraints` binds compliance/volume, original signed
responses and optional sampled strength to ONE constraint policy. The joint
evaluator shares the existing filter, native stiffness preparation and affine
primal fields. The optimizer's `new_with_constraints` and the initial/refitted
reference-grid `*_with_constraints` methods preserve that policy throughout
audits, trials, feasible archives and fresh physical rechecks. Response-only
entry points delegate to this same owner. The joint upper aggregate uses the
existing unnormalized log-sum-exp and exact component weights.

`affine_load_space::stress` applies the existing `SampledPeakBound` model to
the entire load ball. Squared von Mises stress at EACH retained quadrature point
is its own affine quadratic; the native bounded quadratic-ball solver finds
that point's adverse coefficients with the original oracle policy. The
unweighted p-norm of point envelopes bounds both qp-relaxed reference stress
and actual SIMP-scaled physical stress. This is conservative relative to the
maximum of a single global p-norm with one common load direction. Sample count
affects the problem; the result is neither a continuum bound nor an interval
certificate. The cap is explicit physical stress data, never inferred from E.

The stress adjoint factors by baseline and uncertainty fields, needing at most
`dimension+1` adjoints irrespective of point count. It reuses the native tensor
stress pullback, the existing paired-stress coefficient and density derivative,
stiffness contractions and baseline-only `B^T z`, followed by the common
physical-region/projection/filter pullback. Positive hard/flat/near-hard point
contributions retain values but withhold a smooth gradient. No stress basis
direction is admitted as a stand-alone physical state: the native sampler
checks the ENTIRE affine Mandel-strain ball against each retained source ceiling.
Budgets count all baseline/basis point fields before native output allocation.

Strength-bearing optimization uses the existing qp-turnover admission on the
ERODED projected density floor and checks the actual trial fields as well.
Immutable prescribed voids retain the existing exception; free densities cannot
exploit the ersatz foldback branch. Pure evaluation still exposes the declared
stress measure across its full admitted density range. Both original relaxed
and physical stresses and the conservative aggregate must pass the physical cap;
AL stopping tolerances never excuse excess.

`affine_load_space::refit` binds pure reference dead/basis laws to an initial
or explicitly refined native grid and invokes the same affine design owner.
`refit_reference_affine_load_space` immutably borrows the source, checks genuine
enrichment and complete native material compatibility, copies raw controls and
protected physical-region labels through the parent map, preserves the filter
and nominal SIMP/projection policy, and rebuilds the common acceleration map.
The candidate performs its own filter, gradient audit and constrained solves.
Its baseline may be infeasible; adding material is not assumed to restore the
constraint when that material also changes self-weight.

The owned refit result retains initial, latest accepted, best feasible and
freshly rechecked points separately, with the original AL outcome and spent
work. Numerical initialization failure retains the candidate and its failure;
later stops retain the completed optimization prefix. A completed bounded
allowance and fresh full-envelope physical-constraint check are required by
`completed_and_rechecked`. The source is never replaced by an explicit refit,
and objective improvement is compared within the candidate grid, without a
cross-grid descent claim. Both existing affine design goal modes are supported.

`response::reference::ReferenceAffineResponseLimits3` binds the existing pure
displacement-observation laws and live reaction modes to this same owner.
`solve_reference_affine_load_space_with_responses` and
`refit_reference_affine_load_space_with_responses` reintegrate observations in
each grid's native DOF ordering, then retain every original inequality through
the candidate audit, optimization and fresh recheck. The refit completion gate
requires the audited response policy and the initial ordered bounds to survive
the recheck. No coarse observation vector, frozen reaction row or compliance
maximizer substitutes for the original response family.

`estimate_affine_load_space_enrichment` independently recomputes the full
affine load model and maximizer on an actual enriched adaptive grid. The same
grid evaluator serves fixed-grid optimization and enrichment; fine stiffness
scales and physical occupancy inherit the coarse eroded fields through the
checked parent map. Pure dead/basis reference laws are reintegrated on each
grid, and the retained acceleration recipe prepares the fine native body map.
The fine operator's existing solve wrapper owns its selected preconditioner
and setup accounting. No refiltering, nodal-force prolongation, or second
elasticity/uncertainty implementation is introduced.

The assessment reconstructs both grids' fields at the coarse maximizing
coefficients and again at the fine maximizing coefficients. Both fixed-load
comparisons use the existing self-weight-aware native DWR estimator with total
force as the compliance goal and zero affine goal offset. The change caused
by selecting a different adverse load is reported separately and checked in
the two-grid compliance identity. Hard-case witnesses remain valid here;
smooth design gradients are unnecessary for this assessment. Cell markers
describe these frozen loads and do not become a full-envelope error bound.

The whole-set change observation includes every affine term:
`||G_f-G_c||_F + 2*||h_f-h_c||_2 + |c_f-c_c|`. In exact arithmetic this bounds
the normalized quadratic change for every coefficient in the unit ball; its
floating-point evaluation is numerical evidence, not an interval enclosure.
It detects a changed adverse direction even when grid maxima happen to agree.
All returned outcomes preserve the source study, and incomplete grid families
are never published as a usable assessment.
For the basic unbound assessment, `stable_and_feasible` explicitly interprets
the reference as a compliance cap: both grids must satisfy it in normalized
and physical units, and the whole model change must meet its separate tolerance.
A goal-bound assessment instead enforces its original design formulation.
Change tolerances never waive physical cap excess.

`estimate_affine_load_space_enrichment_with_inventory` additionally evaluates
the dilated objective on both grids. It shares the original forward filter,
copies the resulting physical dilated occupancy through the same native parent
map, and integrates each grid's own volume or reference-material mass weights.
It reports physical kg, normalized objective change and geometric volume change
separately. The ordinary mass evaluator, response-limit refinement and affine
refinement share the native weighted-occupancy integration; no SIMP stiffness
floor is interpreted as physical mass. The additional inventory needs no
elasticity solve, fine refilter or inventory adjoint.

`estimate_affine_load_space_enrichment_with_responses` adds the same reference
response family to this complete compliance/inventory experiment. Each grid
reintegrates its own observations and measures its already-solved offset and
basis fields. For each original inequality it retains both native witnesses
and the whole response-model change
`(|y0_f-y0_c| + ||v_f-v_c||_2)/scale`. This measures every direction in the load
ball, including changes hidden by equal worst values. The maximum of these
changes has its own tolerance. Both grids must pass every original bound and
the same joint unnormalized upper constraint used by the optimizer. These
checks need no extra primal or adjoint solve; they are numerical two-grid
evidence, with no continuum or interval-certificate claim.

`estimate_affine_load_space_enrichment_with_constraints` additionally carries
the original sampled-strength policy. Both actual native quadrature sets are
independently remaximized over the entire load ball. To compare complete stress
fields at matching sites, the owner natively prolongates every coarse affine
displacement field and samples that family and the actual fine family with the
SAME fine material and retained points. Each physical family separately passes
native source-domain admission; a difference field is never misrepresented as
a physical state. Occupancy and scales must be literal parent copies.

At each matched fine site the numerical bound
`VM(delta_sigma0) + sqrt(sum_j VM(delta_sigma_j)^2)` bounds the complete-ball
reference stress change by triangle/Cauchy inequalities. The report retains
reference, relaxed and physical forms normalized by the original stress cap.
Acceptance resolves the actually constrained relaxed/physical forms and the
change of the actual p-norm aggregate. Unscaled reference change remains a
diagnostic. The p-norm's sample-count dependence can keep a campaign unresolved
even for a nearly unchanged field; no adaptive-convergence or continuum-strength
guarantee follows. These checks require no extra primal or adjoint solve.

`estimate_affine_load_space_enrichment_for_goal` binds either original design
goal to that same assessment owner. Minimum-material mode retains the physical
compliance cap and selected dilated volume/kg objective. Worst-compliance mode
uses the compliance reference ONLY for normalization, evaluates native dilated
volume on both grids and enforces the original strict volume limit. Its joint
upper constraint uses `volume/volume_limit - 1`, while minimum-material mode
uses `maximum - 1`; all original response and strength bounds remain active.
The report retains the goal and inventory binding, and `objective_change()`
uses the normalized worst-compliance change or the selected inventory change
as appropriate. Neither mode adds a second solver, filter or optimizer.

`affine_load_space::campaign::AffineLoadSpaceCampaign3` automatically repeats
full assessment and refitting for BOTH design goals. Its original physical cap,
objective, filter/SIMP parameters, pure dead/basis laws and native acceleration
recipe remain fixed. Every probe refines all active source cells and
independently remaximizes the complete affine quadratic. A pass requires both
grids' original physical bounds, the full-model change gate, dilated-volume
change and the actual normalized-objective change. Grid histories report the
objective used by the optimizer. Uniform candidate refinement is the default;
optional frozen-witness marking remains a diagnostic heuristic with uniform
fallback for unresolved volume/material inventory or empty marks. Compliance
objective change can use that marking; selective refinement never substitutes
for the full final probe.

Candidates inherit raw controls and protected labels, rebuild native self-weight,
and reoptimize on their own filtered grid. A completed allowance and fresh
full-envelope cap check precede atomic promotion of the tree, study and fields.
Rejected candidate prefixes remain separate from the accepted design. Cumulative
refit and attempt caps span all run segments, and failed attempts are charged.
Completed probes and pending promotions are reused; original work-budget and
cancellation reasons survive native callbacks. `run(0)` assesses the current
design without attempting a candidate. Each successful promotion receives a
fresh full probe before a normal pause or count limit. The caller's builder
must preserve physical identity and enforce its cumulative geometry budget.

`AffineLoadSpaceCampaign3::new_with_response_limits` preserves the full response
family in every probe and refit. Admission checks the original audit policy and
ordered bound metadata before invoking the grid builder; the original `new`
constructor refuses a response-bearing design instead of dropping its limits.
The response model has an independent `response_change_tolerance`; unresolved
response change or physical bounds trigger uniform candidate refinement when
frozen-compliance-witness marking was requested. Promotion requires a fresh
full response/compliance recheck, and the final probe checks the joint upper
constraint, all original bounds, full load/response models and native inventory.
Reference callbacks remain pure model inputs whose physical identity the caller
must preserve across grids.

`AffineLoadSpaceCampaign3::new_with_constraints` retains compliance, responses
and sampled strength as one original family. Both older constructors delegate
to that owner with their explicit subsets and refuse a starting design whose
audited constraints would be omitted. The campaign always uses the same joint
refit/assessment dispatch, with independent `stress_change_tolerance`. Unresolved
strength fields or aggregate trigger uniform candidate refinement when frozen
compliance-witness marking was requested. Promotion and the final probe require
all original bounds and the same joint upper constraint, without spending
optimizer tolerance as physical slack. Histories retain stress aggregate,
constraint and resolution changes alongside compliance, responses and kilograms.

`material_affine_load_space_sdf3` exercises native force balance and load-space
sampling, full-chain compliance/volume derivatives, hard-case refusal and an
actual minimum-kg constrained run with fresh feasible-point reassessment.
Its adaptive regression introduces a coarse-unresolved spatial load and checks
the changed fine adversary, inherited material/body force, both DWR identities,
signed switching, source rollback and the direct-solve work bound. It also
assesses hard-case witnesses without a design gradient. Explicit and automatic
refit regressions check candidate-grid kg descent, a fresh final full-model and
inventory probe, native body-load replay, and preservation of an accepted source
when a completed candidate cannot meet the physical cap.
The response regressions cover independent adverse directions, upper/lower
limits, both adjoint layouts, material-body derivatives, native reaction direct
terms, cancellation rollback and nonsmooth-radius admission. The adaptive
response regression independently reintegrates observations on the promoted
and final probe grids and checks the complete response-model change and native
kg objective; a separate admission check prevents silent response omission.
Full-ball strength regressions compare point envelopes against native tensor
samples and angular loads, audit complete material/self-weight derivatives,
and exercise strict strength archives, budgets and rollback. Native sampler
regressions check the whole material strain domain and adaptive constraints.
Joint adaptive regressions compare matched stress fields independently, retain
all constraint families and native kilograms, check assessment solve counts,
and refuse campaign admission that would omit optimized strength requirements.
Worst-compliance campaign regressions retain original volume, response and
strength limits across Probe/Refit/Probe, including a normalized compliance
above one that is valid because the reference is a scale rather than a cap.
These APIs cover fixed geometry, homogeneous supports, fixed dead loads and
reference acceleration. Numerical load optimality evidence is not an interval
enclosure, a continuum guarantee, a global design optimum, or a serialized
restart contract for this design owner. Explicit refits and automatic adaptive
campaigns support both original design goals. Optimization descent is assessed
within each candidate's own discretization; cross-grid descent is not claimed.

## 3-D cut-cell continuation

`sdf3::continuation` provides `controlled_sdf3_continuation` and
`controlled_gradient_checked_sdf3_continuation` on both Cartesian and adaptive
`CutDensityStudy3` backends. Each finite SIMP/projection stage uses the existing
volume restoration and independently loaded OC solver; no geometry, quadrature,
support, load or filter rebuild occurs during the schedule. Sharpening a
projection first measures the incoming design under the new model and restores
its projected-volume feasibility before establishing a fresh equilibrium
baseline. Compliance descent is measured within each stage, never across
different material models. `MultiLoadContinuationReport` retains each model's
actual history, restoration scale, gradient evidence and stopping reason.

The optional required gradient gate and standalone
`controlled_sdf3_gradient_check` reuse the existing second-order one-sided
heterogeneous-direction experiment. They independently solve every load at up
to five designs and compare both compliance and physical-volume derivatives
through the complete graph-filter/projection/SIMP map. Checks restore incoming
stiffness scales on success and returned failure. A rejected gate or interrupted
stage restores both parameters and scales to the latest completely solved
accepted state. Restoration, probe, setup and rejected-trial work all consumes
the same `SolveControl` budget; no work is refunded by rollback.

Run the existing raw-implicit example with
`cargo run -p fs-topopt --features cutfem-marquee --example elastic_sdf3 -- --continuation 5 250000`.
It reports per-stage compliance, volume, derivative discrepancies and actual
stops before exporting the retained densities. The unchanged command without
`--continuation` still executes one fixed model. The schedule completion label
does not certify stationarity, global optimality, continuum accuracy,
manufacturing thickness or full elasticity-marquee acceptance.

`sdf3::adaptive_continuation::controlled_adaptive_sdf3_continuation` now joins
the existing real compliance DWR estimator to this gradient-checked schedule.
Before each later stage it enriches the accepted background, inherits physical
stiffness and re-solves the reference load family. Weighted absolute local
goal residuals select original leaves through deterministic Dörfler marking.
The marked candidate inherits raw densities, rebuilds the same physical-radius
filter, restores projected-volume feasibility and passes the compliance AND
volume gradient gates before optimization. Reference body/surface laws are
reintegrated on every background. Coarse nodal loads and solved fields never
initialize the new stage; the enriched stiffness experiment is distinct from
the new filtered density model and cannot establish cross-model descent.

Each refined proposal remains separate until its entire bounded stage and final
publication checkpoint complete. Geometry/level/leaf/transfer limits, numerical
failures, gradient-gate failures and cancellation leave the previous accepted
tree, model, scales, design and fields usable. All rejected numerical work is
still charged. Stage zero retains the fixed-grid driver's solved-prefix behavior.
The operator builder owns unchanged implicit geometry, material, support and
quadrature policy, plus cumulative geometry budgets; the driver checks active
tree membership and the existing transfer's geometric/constraint compatibility.
Octree refinement changes the background; it never meshes or moves the implicit
boundary. DWR is a numerical two-grid correction with consistency terms, not a
continuum enclosure or thickness/topology certification.

Run the real adaptive schedule with
`cargo run -p fs-topopt --features cutfem-marquee --example elastic_sdf3 -- --adaptive-continuation 3 250000`.
It shares a cumulative quadrature allowance and reports marked cells, actual
refinement acceptance, gradient discrepancies, restoration, stop causes and the
last accepted raw/projected design. Focused G1/G4 tests exercise three stages,
two actual DWR refinement cycles, each new design map's gradient/volume gates,
second-proposal cancellation preserving the first refined solution exactly,
and an enrichment leaf-cap stop. This does not close the full marquee bead's
moving-boundary, witness-oracle, continuum-bound or both-ISA golden requirements.

### Load-space optimizer restoration

`LoadSpaceDesignStudy3` and `LoadSpaceLimitedDesign3` share the existing projected
augmented-Lagrangian state and expose `checkpoint` / `restore`. Restoration
preserves accepted history, the original gradient experiment and the distinct
feasible incumbent. It freshly solves the accepted and incumbent endpoints and
compares physical responses and all objective/constraint gradients before
continuing. These one or two evaluations count against the original optimizer
allowance; the supplied control charges every solve. Earlier search steps and
the original gradient experiment are not replayed. The consumer must bind the
checkpoint to unchanged geometry, loads, material, physical regions and policy;
an endpoint comparison alone cannot establish global source equivalence.
G3/G4/G5 regressions cover changed responses, lost incumbents, cancellation,
exhausted allowances and split/restored trajectory equivalence. The checkpoint
conveys numerical state, not a continuum certificate or a new resource grant.

## Invariants

- Every stage of the density chain has an exact derivative; the
  composed sensitivity is FD-verified at MULTIPLE continuation
  stages (p = 1 → 3, β = 1 → 8) per the acceptance.
- The filter preserves constants and is symmetric in the
  volume-weighted pairing (mesh-independent length scale r).
- OC keeps designs in [1e−3, 1] with move limits; the volume
  constraint tracks the projected design.

## Error model

Structured panics on solver failures and invalid materials
(modeling errors). Filter construction likewise refuses negative, non-finite,
or non-representable radii, empty complexes, malformed/non-finite geometry, and
invalid cell/vertex measures; its assembled radius is immutable. Filter
application refuses wrong-length or non-finite vectors. Projection helpers refuse non-finite
densities, negative/non-finite sharpness, and thresholds outside `[0, 1]`, while
the well-defined zero-sharpness limit is the identity map. The design pipelines
also refuse nonphysical SIMP parameters, raw densities outside `[0, 1]`, malformed
physics sensitivities or loads, mismatched elasticity/design shapes, and invalid
robust threshold offsets before mutating the elasticity operator. Nominal and
robust OC refuse mismatched vector shapes, non-finite state, nonpositive cell
volumes, and volume targets or move limits outside `(0, 1]`; no `zip` truncation
can silently change an optimization problem. Optimization outcomes are reported
traces (compliance, volume, final change), never silent.
The public DWR band-planning helper is fail-closed: it validates its
level and every entry in the supplied indicator map, then plans recursive
halo refinement on a
clone, so a structured input refusal leaves both the caller's grid and
band level unchanged. The marquee's split-phase path additionally binds the
advance receipt to the exact ordered leaf set, maximum level, prior band
level, and policy version; a stale-grid receipt refuses before mutation.
Obsolete policy-version receipts likewise refuse before mutation.
Every iteration row names its analysis and post-update target generations,
with non-cryptographic state witnesses and separately labeled metrics. The
report re-solves compliance on the final conformed grid, so
`final_compliance` describes the returned `design`, not the prior analysis
snapshot.
The marquee flux probes read the solved field only through the canonical
fail-closed `fs_cutfem::Space::sample_scalar` (bead ay40): missing or
non-finite active nodal evidence propagates as
`CutFemError::InvalidFemInput` out of `run_marquee` instead of reading
as a plausible zero flux; a certified-Outside leaf reads the
homogeneous Dirichlet exterior u = 0 as an explicit mapping, not a
fallback.

## Determinism class

Bit-deterministic: fixed bisection schedules, deterministic solves
throughout; a WHOLE topology-optimization run replays bitwise
(G5-tested). Golden FNV-64 over pipeline stages, compliance
gradient, and a short OC run: `0x772a_2f8c_a720_dd64`; robust
three-field golden `0x519a_41e3_466e_4b7d`. Recorded on
Apple M4 Pro, verified on Threadripper (x86_64).

## Cancellation behavior

Iteration-granular through the resumable fs-solver states; OC
iterations are bounded and the driver can stop between them. Cx
wiring is driver scope.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

- `cutfem-marquee` [F] (default OFF, bead b7d0) — the CutFEM-quadtree
  marquee topology lane (`dep:fs-cutfem`, `dep:fs-dwr`, `dep:fs-ivl`);
  gates the `marquee` integration target.

## Conformance tests

`tests/topopt_battery.rs`: filter G0 laws (linearity ≤
1e−9, transpose adjointness ≤ 1e−9, constants preserved); projection
G0 (exact endpoints, monotone on a 100-point sweep, slope vs FD ≤
1e−8); FULL-CHAIN sensitivity vs FD at three continuation stages
(rel ≤ 2e−4 through solve+SIMP+projection+filter); OC cantilever
(kuhn(3), fixed face + edge load): compliance reduced ≥ 20%, volume
within 0.03 of the 0.4 target, design range > 0.5 (not gray), and
the ENTIRE run replaying bitwise; three-field pointwise ordering on
random designs + eroded-compliance sensitivity FD gate (rel ≤ 2e−4);
robust OC vs the non-robust baseline AUDITED WITH THE SAME
three-field probe — eroded compliance descends, volumes ordered,
nominal volume on budget, and erosion retention at least matching
the baseline (the min-length-scale claim, measured); EIGENFREQUENCY —
aggregate gradient vs FD through the whole chain (rel 1.9e−10),
clustered aggregate FD-verified near the symmetric-bending
near-crossing (3.0e−8), spurious-void-mode gate (λ_min stays at
9.3e−2 on a mostly-void design), +67% λ_agg ascent demo at fixed
volume; STRESS — aggregate gradient FD at two continuation stages
(≤1.3e−8; non-self-adjoint path), singularity-trap regression (void
share 0.000% of σ_max AND the gradient FD-verifiable AT the void
floor), −22% max-stress descent demo at fixed volume; FOUR cross-ISA
golden hashes (pipeline, robust, eigenfrequency
`0xbb7e_5ad3_851a_2bf1`, stress `0xc539_ad97_34d8_1b66`). Plus
`tests/probe_robust.rs`: the limit-cycle regression. Feature-gated
`tests/marquee.rs` additionally covers the heat solve/update/refine loop,
zero rebuilds, cut-band concentration, estimator-time gate evidence,
monotone one-level headroom, exact `leaf_delta = 3 * splits` accounting,
and deterministic synthetic
`CellKey` indicators driving `refine_dwr_cut_band` for exactly one
planning step. G3 checks show that indicator negation, exact power-of-two
scaling, reverse insertion, and within-partition mass-preserving
redistribution preserve authorization and target topology.
An area-proportional indicator map is refused even though its band share
exceeds 0.15, and one extra unit on a band leaf authorizes the advance.
The iteration log emits policy version, cut/band/total mass, band area
fraction, cut fraction, band
levels, exhaustive decision reason, analysis/target design witnesses and
metrics, target leaf counts, policy splits, post-update motion-conformance
splits, the fully conformed target-grid leaf count, and their reconstructed
total. Leaf-count chaining, together with the split-only implementation and
absence of an intervening grid mutation, accounts for continuity into each
next analysis and the final re-solve; cardinality alone is not claimed as a
topology-identity certificate. Non-cryptographic design witnesses correlate
each target with the next analysis snapshot and the returned design. The
returned design receives an explicit finite final solve.
The battery separately pins disabled, zero-mass, strict-threshold,
and exhausted-headroom no-op reasons; invalid-level,
non-finite-indicator/accumulated-mass/recursive-enclosure, and non-leaf-key
cases assert structured refusal with no grid or band-level mutation. The
literal vector integration uses `estimate_elasticity_compliance` on an
off-grid embedded disk large enough to retain coarse interior active leaves,
and passes its UNMODIFIED real vector indicator map into the same helper for
two consecutive cycles. Each cycle asserts authentic cut-mass admission,
positive actual halo/balance splits, a mixed-active-level vector re-solve, and
bit-identical structure, coefficients, reconstructed nodes, and compliance on
replay; coarse/enriched compliance and estimator metadata are emitted per
cycle.

## No-claim boundaries

- Scope: compliance/volume, robust three-field, eigenfrequency, and
  stress-aggregate objectives on FIXED kuhn meshes. The medial-axis
  thickness oracle (geometry-layer audit) and the CutFEM-quadtree
  heat marquee (zero remeshing + DWR-gated cut-band adaptivity on its
  recorded volume-to-point fixture) are the feature-gated extension;
  elasticity benchmark envelopes remain outside this contract.
- Eigen solves are dense (fixture scale); LOBPCG-scale pencils join
  via fs-solid stability's machinery when the consumer needs them.
- Stress constraints ship as the aggregate + gradient; the
  constrained DRIVER (AL with adaptive c_k updates per iteration)
  composes from fs-ascent's augmented Lagrangian at the consumer.
- OC is the compliance/volume driver; MMA is not implemented
  (fs-ascent AL is the general constrained path — documented
  choice).
- No multi-load/worst-case formulations, no continuation SCHEDULER
  (drivers own β/p ramps; the primitives take fixed parameters).

## No-claim boundaries (marquee)

- Heat-conduction (volume-to-point) benchmark class: the recorded
  compliance envelope is THIS fixture's golden band; no MBB/cantilever
  ELASTICITY compliance envelope is claimed here.
- The shared refinement helper accepts indicators from scalar or vector
  estimators. Its real-vector integration proves exactly two deterministic
  estimator-to-policy-to-graded-re-solve cycles on the recorded embedded-disk
  fixture. This is not an elasticity optimization benchmark, a general
  adaptive stopping policy, or a certified DWR bound. The logged
  `eta_signed / (J_h2 - J_h)` value remains an enriched-delta diagnostic, not
  an independent-reference effectivity claim.
- Per-iteration wall times are DEBUG-build measurements, labeled; the
  interactive-cadence targets are perf-CI's gates.
- The flux redistribution is a smeared shape-derivative heuristic with
  measured descent, not a certified gradient (the certified route is
  the marquee-demo crate's exact bilinear identity).
- Split concentration is bounded by the halo contract at ~2/3, not the
  0.8 a halo-free marker could reach — the ceiling is the solver's
  ghost-penalty stencil. The 0.60 final-snapshot threshold is a recorded
  fixture golden across its short moving-boundary trajectory; the stronger
  advancing and non-advancing motion tests derive the target halo refinement
  independently, require exact leaf-set equality, and build the next scalar
  space on the conformed result.
- The heat marquee gates refinement on the fraction of absolute DWR mass
  carried by zero-straddling cells of the SAME snapshot that produced those
  indicators. This is a global one-level authorization, not Dörfler/top-k
  cell marking. Its target is the post-update cut band and halo used by the
  next solve; no claim is made that old per-cell residuals localize error on
  the moved geometry. It does not multiply indicators by `|grad rho|`, and
  no such weighting is claimed. Absolute split counts are fixture outcomes,
  not policy invariants: the corrected independent coarse-adjoint estimator
  can legitimately change their value while the versioned gate, structural
  accounting, and boundary-concentration invariants remain intact.
- `state_witness` is deterministic forensic correlation metadata only. It is
  not collision-resistant, not a content address, and not authority for
  scientific identity; the full `DensityDesign` remains the artifact of
  record.
