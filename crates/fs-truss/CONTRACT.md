# fs-truss CONTRACT

## Unified prescribed-load catalog design

`material::MaterialLayout::solve_and_size` returns the same complete
`MaterialBodyLoadDesign` used by the existing self-weight solve, with the
original prescribed load case, zero body acceleration and zero body-load
iterations. Both public paths call the same convergence and sizing acceptance
owner. Original member ordering, pruned zero areas, material mass and source
assignments remain consistent. This reports a numerically accepted volume
layout and catalog design; it does not promote diagnostics to an optimum proof.

## Purpose and layer

Layer: L4 (ASCENT). Ground-structure truss layout optimization (plan
§9.5 [S/F], bead 7tv.13): candidate members under fabrication rules →
plastic-design LP solved by an in-house first-order primal-dual
iteration with explicit convergence diagnostics → Euler/code sizing
with catalog snapping → fs-solid rod re-analysis. The
steel-and-concrete flagship's engine (§15.2).

## Displacement-compatible elastic reanalysis

`elastic::solve` and `MaterialLayout::elastic_response` analyze the actual
selected areas with each original member's Young's modulus. The native checked
CG owner applies `K = B diag(EA/L) B^T` without assembling a global matrix and
uses Jacobi preconditioning. Memory and each action are linear in node/member
count. Member/recovery loops poll cancellation every 256 entries; native CG
retains its own iteration checkpoints. The total iteration budget also covers
at most two defect corrections. No new numerical recurrence is used.

`ElasticImposedState` adds optional complete original-order member stress-free
engineering strains and nodal prescribed translations. Omission means zero
without zero-array allocation. Prescribed nonzero values are admitted only on
supported coordinates; free-coordinate values must be zero. All supplied
entries must be finite, including removed members. A caller supplies thermal
strain from its own expansion law and reference temperature, or explicit
assembly misfit. No expansion coefficient, stress or modified modulus is
inferred. Geometry remains the original infinitesimal reference geometry.

`elastic::solve_with_imposed` and
`MaterialLayout::elastic_response_with_imposed` share the same native solve.
With elastic extension `e_el = B^T u - L epsilon_0`, the member force is
`N = (EA/L) e_el`. Known supported motion is lifted exactly into the free RHS:
`K_ff u_f = f_f + g_f`, where
`g_f = B_f diag(EA/L) (L epsilon_0 - B_p^T u_p)`.
The material consumer retains each original source-resolved modulus at its
admitted material state. Its separate stress-free reference is caller-owned.

The returned original-order displacements generate total strain
`member_strain`, elastic strain `member_elastic_strain`, signed forces and
stresses. Supported displacements equal the prescribed values exactly.
Support reactions are recovered internal forces minus actual nodal loads,
including loads applied directly at moving supports. Physical stored energy is
`U = sum(1/2 N e_el)`; external-load compliance is `f dot u`. With eigenstrain
or support motion compliance need not equal `2 U`. Equilibrium instead gives
`2 U = f dot u + R dot u_p - sum(N L epsilon_0)`, where R is the support reaction.

Acceptance recomputes physical equilibrium from recovered forces on every
original free coordinate. The normalized Euclidean residual uses
`max(||f_f||, ||g_f||, ||f_f + g_f||)` as its reference; this avoids a singular
relative measure when external loads vanish or imposed and external loads
cancel. If that reference vanishes, the residual is absolute in newtons.
A recursive CG residual alone cannot publish an accepted response. Zero-area
members contribute no stiffness, strain, force or energy. Unloaded
zero-stiffness free coordinates receive zero displacement; loaded ones refuse
without regularizing springs.

This explicitly supplements the existing **plastic-equilibrium** layout/sizing
audit: a redundant truss's statically admissible LP/refitted force distribution
need not be elastically compatible after snapping. Callers may pass the paired
`areas_m2` and `load_case` of a `MaterialBodyLoadDesign`, retaining its actual
self-weight loads. The original audit's semantics are unchanged. Reanalysis
does not certify strength, displacement limits, uniqueness, stability, buckling,
finite deformation, bending, geometric prestress stiffness or inelastic
redistribution. Numerical equilibrium alone does not exclude an unloaded
mechanism. Explicit stress-free strain supplies axial prestress, but does not
establish the stability of that prestressed configuration.

G1/G3 regressions exercise serial members, a redundant three-member analytic
solution, stiffness-dependent load redistribution, supported loads/reactions,
energy, rotation covariance, load reversal and pruning. G4 regressions target
loaded mechanisms, invalid material arrays, exhausted work and cancellation.
Imposed-state regressions cover free thermal expansion, fully restrained
thermal stress, heterogeneous serial expansion, support settlement with
external loads, the full energy-work identity, zero-state compatibility,
pruned-member semantics, malformed inputs and cancellation. Material consumer
regressions carry source-resolved moduli through both paths.

## Compatible elastic design sensitivities

`elastic::sensitivity::displacement_functional` computes
`J = sum(g dot u)` and its implicit design derivatives using the same
matrix-free elastic owner for the primal and one adjoint. The adjoint has
homogeneous prescribed motion and solves `K_ff lambda_f = g_f`; supported
weights remain in the case so their direct support-motion derivative is
retained. The total native iteration budget is shared between both solves.
Both recovered-force residuals must pass before any derivative is returned;
no differentiation through Krylov iterations or second recurrence is used.

For each retained original member, let `e_el = B^T u - L epsilon_0` and
`z = B_f^T lambda_f`. The returned partial derivatives are:

- Area: `dJ/dA = -(E/L) e_el z`.
- Young's modulus: `dJ/dE = -(A/L) e_el z`.
- Stress-free strain: `dJ/d(epsilon_0) = EA z`.
- Actual nodal loads: `dJ/df = lambda`, with zero supported entries.
- Prescribed motion: `dJ/du_p = g_p - K_pf lambda_f`, with zero free entries.

Using elastic extension includes the derivative of the eigenstrain/support
right-hand side as stiffness changes. Substituting the usual homogeneous
compliance formula can give the wrong sign under imposed motion. The original
node/member ordering is retained; zero-area members return `None`, explicitly
excluding derivatives that activate removed topology. Compatible mechanisms
remain possible: reported load derivatives apply only to perturbations that
retain solvable equilibrium on the same topology. An unsupported adjoint
objective refuses even when the original load case alone is solvable.

`MaterialLayout::elastic_displacement_sensitivity` retains each original
source-resolved E. Geometry, support flags, material source states, weights and
external loads are frozen for the member partials. Passing `case.loads()` as
weights gives fixed-load compliance. If a design parameter changes loads,
the caller adds `load_gradient dot df/dtheta`; changing objective weights
also adds `u dot dg/dtheta`. Both terms are needed for load-dependent
compliance, including self-weight: when `g = f` tracks the loads, its total
load derivative is `lambda + u`. This reduces to `2u` only for zero eigenstrain
and homogeneous supports. No material-state, density, catalog,
geometry, geometric-stiffness or fixed-point derivative is inferred.

G1/G3 tests check a settlement/thermal compliance case with a positive area
derivative, redundant rotated trusses against centered finite differences of
area, E, eigenstrain, support motion and actual loads, and source-bound
heterogeneous members. G0/G4 tests preserve removed-member indexing and refuse
malformed weights, unsupported adjoints, exhausted shared work and cancellation.

## Global truss buckling with a fixed prestress baseline

`elastic::buckling::solve` and `MaterialLayout::elastic_buckling` use the
existing compatible elastic owner twice: once for fixed nodal loads, fixed
member eigenstrain and prescribed support motion, and once for the reference
case's nodal loads with homogeneous support motion and zero eigenstrain.
`BucklingLoadFamily` makes that split explicit. Omitted fixed loads mean zero;
computed self-weight can be supplied as fixed loads when gravity must remain
constant. Every reference-case load is multiplied by the returned factor.
The two solves share one `ElasticControls` iteration allowance.

For a member's original axial direction n and transverse direction t, the
initial-stress contribution on relative endpoint translations is
`K_g = (N/L) t t^T`. Tension stiffens and compression softens. Using the actual
recovered fixed and reference member forces gives the pencil

`(K_e + K_g(N_fixed) + lambda K_g(N_reference)) phi = 0`.

Thus eigenstrain, support settlement and fixed dead loads alter the baseline
tangent; they are not silently multiplied with the reference load. All
original free coordinates are retained in this pencil, including unloaded
zero-stiffness coordinates that equilibrium alone can leave at zero motion.
The native positive-factor gate refuses a singular or indefinite baseline
without adding springs. A numerical/factorization refusal is not itself a
classification of the physical instability.

The domain assembles canonical sparse matrices with `fs-sparse`; it calls
`fs-solid::stability::buckling_loads_cancellable` for scaling, positive-base
admission, eigenpairs and physical-pencil residuals. Large pencils retain the
existing bounded sparse factorization and generic `fs-spectral` dispatch;
small systems retain the native bounded dense path. No new eigensolver is
introduced. Member assembly polls cancellation every 256 members. Combined
staged entries have an explicit caller allowance and a hard 8-million cap,
bounding native COO sorting between checkpoints. Native default DOF, stored
matrix, factor-memory, fill and operator-work limits also apply.

`TrussBuckling` retains both admitted equilibrium states and the native
`BucklingResult`, including the complete original-coordinate map and
recomputed pencil residuals. At multiplier lambda, displacements and member
forces superpose as fixed plus lambda times reference. Energies are quadratic.
The source-material wrapper retains each original member's resolved E.
Zero-area members supply no stiffness or prestress; their free coordinates
still participate in stability admission.

This is a small-displacement initial-stress analysis on frozen original
geometry and areas. It does not follow geometric nonlinearity, follower
forces, plastic redistribution, limit points, postbuckling, out-of-plane modes
or local member flexural/torsional buckling. Local Euler/section/strength checks remain
independent. Returned positive factors are residual-qualified numerical
predictions, not proofs of missed-mode exclusion, mesh convergence or actual
collapse. An empty factor list says no positive root was returned for the
queried reference-load directions.

G1/G3 tests use an analytically braced unequal-length two-bar system. They
exercise thermal compression/tension, support settlement, independently fixed
dead load, reference-force scaling, rotation and mode direction. Tensile and
zero reference loads produce no positive root. G4 tests ensure unloaded
mechanisms remain visible and refuse unstable baselines, malformed data,
exhausted shared work/storage and cancellation. A material consumer regression
checks the original source E assignments through the same path.

## Prestressed vibration with original reference mass

`elastic::modal::assemble` and `MaterialLayout::elastic_modal_pencil`
bind one accepted compatible equilibrium to its immutable sparse tangent,
consistent reference mass and complete original free-coordinate map.
The material route retains each original assigned source state's modulus and
reference density. It does not infer a dynamic material law or change density
when an explicit eigenstrain or support settlement changes.

Each member contributes `M_e = rho A L / 6 [2I I; I 2I]` in the two Cartesian
translations. This retains transverse kinetic energy and endpoint coupling;
mass is not projected onto the bar axis or silently lumped. Optional explicit
nonnegative joint masses add isotropic translational inertia. Zero-area
members contribute neither stiffness nor mass. The reported physical inventory
is `sum(rho A L) + sum(joint_mass)`, including supported material and joints;
it is not the reduced mass trace. No thermal or prestress mass is created.

The same private member kernel used by global buckling supplies
`K_t = K_e + K_g(N_actual)`, with `K_g = (N/L) t t^T`.
The actual nodal load case, stress-free strains and prescribed translations
establish the base state using the existing checked-CG equilibrium owner.
Modal perturbations keep those loads/imposed states fixed and have zero
displacement on supported coordinates. Original reference geometry is frozen:
this is small-displacement initial-stress vibration, not a finite-deformation
path, follower-load tangent, damping model or local member bending model.

Every original free coordinate remains in the pencil, including zero-stiffness
coordinates and mechanisms. The pencil's `slice_windows` method calls the
existing `fs-modal::slice_windows` owner with caller-supplied windows, explicit
factor/payload budgets and the same cancellation authority. Its mass-SPD gate,
endpoint inertia counts, sparse shift-invert solve, mass normalization and
original-pencil residual checks remain in force. Singular/massless subsystems
are not repaired with invented mass or static condensation. Assembly alone
does not claim positive mass, positive stiffness or complete stability.

Returned eigenvalues are signed `omega^2` in inverse square seconds, never
clipped. Only an interval resolved strictly above zero permits an oscillation
frequency; a strictly negative interval indicates linearized exponential
growth. The native relative interval-width gate can refuse unresolved
near-zero modes, which remain present in the physical pencil. Successful
windows do not establish stability outside the requested spectrum. Repeated
window requests reuse the original equilibrium and physical assembly.

Equilibrium retains its existing budget. Assembly admits at most 100000 free
coordinates and 800000 original coordinates; combined tangent/mass triplets
are caller-bounded with an 8000000 hard ceiling. Loops poll every 256 entries;
the bounded native COO sort is checked immediately before and after.
Density, areas, moduli, optional joint masses and all assembled sums must remain
admissible and finite. G1/G3 regressions cover consistent-mass bar/serial
frequencies, source material assignments, joint/support mass inventory,
rigid/transverse inertia, rotation covariance and thermal/settlement prestress
softening into negative modes. G4 covers malformed inertia, massless free
coordinates, work limits and cancellation.

## Public types and semantics

- `GroundRules` is immutable and admitted only through `try_new` or a finite
  default. `GroundLimits` applies caller limits beneath crate hard ceilings.
  `GroundStructure::try_grid` and `try_from_parts` return immutable nodes,
  canonical members, length-consistent vectors, and FrankenNetworkx graph state
  only after complete validation, bounded pair/triplet/member/retained-byte work,
  fallible vector reservations, and a final cancellation checkpoint. Imported
  within-tolerance lengths are replaced with recomputed canonical `hypot` bits;
  graph compatibility timestamps are cleared before publication.
  Generation is reproducible and `stats()` is the ledger row (counts + FNV
  hash).
- `LayoutCase` immutably admits exact-shape support flags and finite nodal
  loads. `LayoutLimits` bounds free DOFs, split variables, staged sparse
  triplets, and retained LP storage beneath hard ceilings.
- `LayoutLp::try_assemble`: the member-force LP — split tension/compression
  variables `q⁺, q⁻ ≥ 0`, volume objective `Σ l(q⁺+q⁻)/σ_y`, nodal
  equilibrium on free DOFs. Assembly validates canonical ground identity,
  surviving connected load DOFs, exact sparse dimensions/nnz, finite positive
  load norm at least `1e-30` with finite squared norm, costs and norm state,
  and cancellation before publishing immutable LP state. Sparse A/A-transpose
  construction uses fallibly reserved canonical CSR buffers directly and polls
  throughout counting, fill, validation, transpose, and norm multiplies.
  `solve` = PDHG (Chambolle–Pock) with
  power-iteration step sizing, sparse matvecs, warm starts across load
  cases, deterministic iterations. `solve` is fallible, caps direct solves at
  one million iterations, and validates controls plus warm-start shape/domain
  before work. `diagnostics` returns
  (relative primal/dual objective separation, equilibrium residual, volume):
  under this saddle the nominal dual objective is `−bᵀy` with feasibility
  `c + Aᵀy ≥ 0`, approximately restored by a floating uniform shrink of y —
  the battery pinned the OPPOSITE textbook
  convention (`+bᵀy`, `Aᵀy ≤ c`) reporting gap = 2 on exactly-solved
  instances.
- `LayoutLp::certify_optimum` is the separate cold certificate path. It turns
  split iterates into signed member forces, selects a deterministic square
  member basis after dropping only identically-zero/zero-load rows, and proves
  basis invertibility plus an exact equilibrium correction by the outward
  Neumann condition `rho = ||I - H M||_inf < 1`. The retained member-force box
  contains one correlated exactly equilibrated force; its positive/negative
  split is therefore nonnegative by construction. Independent residual boxes
  are sanity checks only: zero containment by itself is not the existence
  proof. The primal upper endpoint uses each sign's own positive cost and
  repaired split-force enclosure. Equal costs retain the tighter original
  absolute-force enclosure. Unequal tensile/compressive limits therefore pass
  through the same outward proof rather than borrowing the tensile coefficient
  for compression.
- The dual witness uses this crate's `c + A^T y >= 0`, `-b^T y` convention.
  A representable uniform shrink is re-evaluated directly from authoritative
  `A` with `fs-ivl`; every outward slack lower endpoint must be nonnegative.
  The finite lower endpoint is combined with the primal endpoint only when
  they are ordered. `LayoutOptimalityCertificate` fields are private and bind
  `A`, `b`, `c`, both input iterates, solver settings, proof budgets, correction
  and arithmetic versions, all witness endpoints, and the final receipt with
  domain-separated `fs-blake3` identities. `verify_optimum_certificate`
  recomputes the proof. `certify_optimum_for_report` additionally requires the
  private identity of the exact `solve` output before retaining bounds in
  `PdhgReport`; unrelated reports cannot acquire them. Domain-separated solve
  snapshots bind every public diagnostic and the complete trace. Public
  diagnostic mutation makes all verified accessors return `None` and omits
  verified fields from `to_json`; the trace itself is private and exposed only
  as an immutable slice. Certificate admission rejects a retained trace-length
  mismatch in O(1), then re-hashes the complete trace with bounded polling.
- `sizing::size_and_snap` → `CatalogAudit`: areas from yield, EULER
  floors for compression members (solid square `A ≥ √(12|q|l²/π²E)`),
  joint parsimony pruning with MANDATORY least-squares equilibrium
  re-verification on survivors (CG on the normal equations), catalog
  UP-snapping (feasibility preserved by construction), member-by-
  member post-snap re-checks as fs-constraint `Code` rows.
- `rodcheck::rod_buckling_check`: the critical compression member as
  an fs-solid Cosserat rod with a seeded bow, loaded to factor×design
  — stable/bow-ratio outcome (the tfz.14/tfz.15 spot check). The initial
  bow perturbs the Newton iterate, not the straight stress-free reference.
  Acceptance requires equilibrium, a positive free-coordinate conservative
  tangent and bounded bow. A converged unstable straight branch is refused,
  including zero/tiny initial bows. The clamped-free boundary has one quarter
  of the pinned-pinned Euler load used by LP sizing. This is a local discrete
  endpoint screen, not a stability proof for the whole loading path or a
  certified continuum buckling load.

## Invariants

Applied nodal loads use the global geometry axes. With a member direction from
endpoint a to b, positive tension maps to external loads `[-direction, +direction]`.
LP assembly and post-prune force refitting share that exact column owner, so a
push toward a support activates compression buckling and a pull activates tension.
The axial-direction regression checks both signs along all four coordinate directions.

### Source-bound material assignments

`material::TrussMaterial::resolve` resolves density, axial Young's modulus and
an explicitly admitted symmetric yield stress. `resolve_axial_limits` instead
requires separately named positive tensile and compressive magnitude limits
through `fs-material::state_point`, retaining the original card, physical query,
property roles and selected claims. Equal property names resolve once. Neither
constructor infers a missing strength, freezes a dynamic strength law, converts
an ultimate strength into yield, or applies an unstated safety factor.
`MaterialLayout` binds one
material row per ground member and reuses the ordinary equilibrium assembly,
PDHG, force refit, Euler sizing and rod spot-check owners. Each sign's own
strength changes the actual split-force objective; local modulus changes compression
buckling, and density changes the mass of the actual snapped member geometry.
No material names, default steel properties or guessed strengths select physics.

The LP costs use one common positive reference-strength scaling for numerical
conditioning. `try_new` selects volume; `try_new_with_objective` can select
`MaterialLayoutObjective::Mass`, multiplying each signed member cost by its
source density divided by one common reference density. Density therefore
changes the actual optimum load path for unlike materials. Uniform density
recovers the original volume operator exactly. The equilibrium matrix, source
strengths, Euler checks and catalog sizing stay on their original owners.
`lp()` and its optional certificates describe that normalized problem;
`MaterialLayoutSolution` reports physical approximate `volume_m3` and `mass_kg`
with the selected objective. For mass, divide a normalized bound by the
reference strength and multiply by the reference density to obtain kg.
`strength_scale_pa()` is the common divisor back to physical volume. Sizing
chooses the tensile or compressive limit from the force AFTER equilibrium
refitting, including when refitting reverses an input iterate's sign. The
existing `SizedMember::area_yield` field is the area from that declared limit;
for `resolve` and the uniform legacy path it is still the yield-based area.
The source-bound catalog audit labels axial limits and effective-length Euler as
physical model checks, not building-code certification. Absent catalog sections
leave mass/volume totals unavailable. Immutable geometry/material assignments
prevent using a different ordering or temperature state during sizing.

`MaterialLayout::with_effective_buckling_length_factors` optionally binds one
finite positive dimensionless `K` per original member; the default is pinned
`K = 1`. The caller owns the independent restraint/bracing justification.
The common compression sizing and post-snap audit both use `P_cr = pi² E I/(K L)²`,
and each sized member retains `K` and `K L`. The immutable assignment persists
through pruning and self-weight rebuilds. Physical length still determines
equilibrium, volume and mass; tensile limits and the axial LP/certificate problem
are unaffected. This does not infer rotational fixity from nodal translation
supports, establish global frame stability, or model inelastic/postbuckling
resistance. The separate Cosserat spot check remains its documented clamped-free
model. Tests in `material.rs` and `material_body_load.rs` check analytic Euler
scaling, refitted compression signs, unchanged tensile sizing/physical mass,
self-weight persistence, invalid assignments and cancellation.

The Euler effective-length relation follows FHWA's
[*Engineering for Structural Stability in Bridge Construction Reference Manual*](https://www.fhwa.dot.gov/bridge/pubs/nhi15044.pdf),
NHI Course 130102, Equation 4-2 and Figure 4-5, printed page 4.10. That reference
treats `K L` as the length of an equivalent pin-ended column and ties `K` to
the actual structural restraints; it does not justify inferring `K` from the
translation-only supports used by this axial layout model.

This is a frozen-state pin-jointed force-equilibrium sizing model with homothetic
sections (solid squares by default),
not elastic-compatible force distribution, rotating anisotropic failure,
post-failure redistribution, creep or thermally restrained strain. The continuous
axial mass objective does not imply globally optimal Euler/catalog sections or
a globally optimal self-weight design. Body-load iterations preserve the chosen
objective when rebuilding the LP at each actual load state. A named ultimate strength is a caller-admitted axial limit,
not evidence of a plastic mechanism or a validated brittle-failure prediction.
Body acceleration remains an explicit physical input. The existing synchronous
PDHG/refit/rod kernels retain their cancellation limitations; material admission
and result publication poll `Cx`. `tests/material.rs` checks source-temperature
substitution, heterogeneous optimal load-path changes, per-member Euler floors,
snapped mass, unsupported states, malformed inputs and cancellation. Its
asymmetric regressions exercise load reversal, refitted sign selection,
independent hand-computed volume/certificate bounds, equal-yield equivalence,
and persistence through self-weight rebuilding. Frame material layouts consume
this same path.

`MaterialLayout::body_loads` derives uniform member mass and equivalent endpoint
forces from those same source densities and actual sections. `solve_with_body_loads`
closes the existing LP/refit/yield/Euler/catalog path against those sections' own
body loads, always rebuilding from the immutable external load case. Acceptance
requires unchanged catalog sections and topology plus a successful common sizing
audit, so loads, mass, sections and equilibrium describe one state. The bounded
fixed-point procedure can refuse nonconvergence or insufficient catalogs and does
not establish a global optimum. Endpoint lumping preserves total force and first
moment but does not resolve distributed bending or within-member axial stress
gradients. Three `material_body_load` integration tests cover heterogeneous force
and moment, density-driven section changes, and numerical/catalog/budget/cancellation
refusals. Native execution remains unverified in the memory-limited authoring worker.

1. Ground rules hold member-by-member and generation is bitwise
   reproducible. Exact-cap admission succeeds; cap-plus-one, malformed parallel
   vectors, noncanonical member identity, inconsistent lengths, unsafe numeric
   state, allocation pressure, and cancellation return structured refusal with
   no partial `GroundStructure` (truss-001 and admission battery).
2. PDHG reaches hand-provable optima (aligned tie `PL/σ`; symmetric
   two-bar `2PL/σ`) to 1e-4 with objective separation < 1e-5,
   equilibrium residual < 1e-5, complementary slackness and observed dual
   feasibility violation < 1e-4 (truss-002). These are numerical oracle checks,
   not by themselves an outward-rounded certificate.
3. A finite optimum interval exists only after both independent certificate
   sides pass. A Neumann contraction proves existence of an exactly feasible
   nonnegative primal split; outward dual slacks prove the lower bound. Tiny
   basis enumeration independently checks the aligned-tie and symmetric
   two-bar optima. Approximate infeasible iterates contribute no raw upper
   bound; rank deficiency, poor conditioning, non-finite arithmetic, budget
   excess, identity mismatch, and cancellation fail closed (truss-002h-k).
4. Densifying the ground structure does not worsen the returned-iterate volume
   beyond the declared diagnostic tolerance
   (truss-003); the Michell closed-form catalogue comparison is a
   LEDGERED PENDING row until its vetted constants land via the
   fs-fab oracle spec — stated, never silently skipped.
5. PDHG cost per (iteration × nnz) is flat across problem sizes
   (spread < 3×) and warm starts reduce iterations on perturbed load
   cases (truss-004; the 10⁶-member wall-clock target is perf-lane
   scope, ledgered).
6. Sizing: post-prune equilibrium re-verified < 1e-6; Euler floors
   active on compression members; post-snap member-by-member audit
   all-pass (truss-005).
7. The rod spot check has teeth: catalog area stable at 1.3× design,
   an under-sized member fails or bows an order harder (truss-006).

## Error model

`GroundRules`, `GroundLimits`, `GroundStructure`, `LayoutCase`, `LayoutLimits`,
and `LayoutLp::try_assemble` return `TrussConstructionError` for invalid input,
parallel-vector shape mismatch, work/retained-byte excess, failed vector
reservation, or observed cancellation. No public ground or layout constructor
panics on caller input or publishes partial state. `LayoutLp::solve` separately
returns `PdhgError` for zero iteration/check intervals, invalid tolerance,
malformed warm-start shape, or non-finite/out-of-domain warm state. The
objective-separation and KKT numbers remain diagnostics. Certificate malformed
state, invalid limits, allocation failure, report substitution, and
cancellation are hard `LayoutCertificateError`s. Sound numerical inability is
`LayoutCertificateStatus::Unavailable` with a typed resource, rank,
conditioning, residual, non-finite, or endpoint refusal; it never fabricates a
wide finite bound. `NaN` catalog area marks an un-satisfiable member in the
audit rather than silently clamping.

## Determinism class

Bit-deterministic across runs on a platform (BTree generation, fixed
iteration/basis/pivot order, deterministic solvers and canonical receipt
hashing). Cross-ISA goldens not yet recorded.

## Cancellation behavior

Ground construction and LP assembly poll an explicit `Cx` at deterministic
bounded strides and immediately before publication. Cancellation returns a
structured `Cancelled { stage }` refusal; no partially built authoritative
value escapes. Certificate construction and admission validation, matrix
visits, elimination, interval work, hashing, verification, and final
publication likewise poll at bounded strides. Report admission re-hashes the
complete solver trace with the same deterministic stride before retaining
bounds; the cheap promotion verifier rejects shape or retained-operation-cap
excess before traversing caller arrays. A cancelled attempt clears report
bounds and publishes no proof. The
PDHG solver remains iteration-bounded but does not yet poll `Cx` inside its
fixed synchronous solve loop; its initial private trace snapshot is part of
that same unpolled solve. Only the cancellation-polled admission re-hash can
authorize later certificate retention.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/battery.rs`: truss-001 rules + determinism plus adversarial construction
admission, exact-cap/cap-plus-one, malformed-part, numerical, load, and
cancellation refusals; truss-002 provable oracles with numerical diagnostics
and malformed-solver-input refusal; truss-002h-k exact bound bracketing against
an independent tiny LP oracle, deterministic replay, private report binding,
post-certification diagnostic mutation plus immutable trace access,
rank/resource/overflow/near-zero fail-closed behavior, tamper identity, and
cancellation-before-publication;
`src/certificate.rs` unit tests: private receipt tamper detection and pre-hash
verification-budget admission;
`src/lp.rs` unit tests: same-length non-tail report-trace admission rejection,
pre-traversal length rejection, and bounded-stride interruption;
truss-003 refinement monotonicity within declared tolerances;
truss-004 scale trend + warm starts; truss-005 sizing/snap audit;
truss-006 rod spot check.

## No-claim boundaries

- SOCP extensions (elastic-compatible layout, stress constraints
  beyond plastic design) — the LP ships; SOCP is the recorded
  successor under the same PDHG surface.
- The vetted Michell closed-form catalogue (0.08-tolerance
  comparisons land with the fs-fab `:oracle (michell …)` spec
  constants).
- 10⁶⁺-member wall-clock budgets (perf lanes; the trend is ledgered
  here).
- 3D ground structures; frame (moment-carrying) layout; connection
  families beyond angle sets; discrete member-count MILP.
- Multi-load-case simultaneous layout (warm starts ship; the
  worst-case envelope LP is follow-up).
- Rank-deficient nonzero equilibrium systems and bases whose Neumann
  contraction cannot be proved remain `Unavailable`, even when some such
  systems are feasible. Generic rectangular exact-feasibility certificates are
  follow-up scope; no residual-containing-zero shortcut is claimed.
- Active-set membership and tropical load-path weights are not implied by an
  optimum interval. `fs-truss-e2e` owns their separate proof: it consumes this
  certificate's repaired signed/split-force intervals and identities, then
  independently requires threshold, member-product, orientation, endpoint,
  unique-path, and unique-bottleneck separation before promotion. This crate
  alone makes no load-path claim.


## Homothetic compression-section geometry

`sizing::HomotheticSection` admits fixed section families through
`SectionFamily`; `MaterialLayout::with_sections` borrows one admitted shape
per **original** ground member. Omission preserves the solid-square default.
Pruning, catalog snapping and every self-weight rebuild retain that assignment
and the independent effective-length factor K.

For actual material area A, both the continuous Euler floor and post-snap
resistance use the same weakest principal relation I_min = c A²:

| Family | Admitted shape parameter | c = I_min / A² |
| --- | --- | --- |
| Solid square | None | 1/12 |
| Solid circle | None | 1/(4 pi) |
| Solid rectangle | Major/minor ratio r >= 1 | 1/(12 r) |
| Circular tube | Inner/outer diameter ratio 0 <= r < 1 | (1+r²)/(4 pi (1-r²)) |

The tube coefficient factors 1-r² as (1-r)(1+r), avoiding thin-wall subtraction.
Nonfinite/out-of-domain ratios refuse. The checked
`minimum_second_moment_m4` query refuses nonpositive/nonfinite areas or
unrepresentable inertia. Nonsquare compression audits require positive finite
inertia and Euler resistance. Unit rectangles and zero-bore tubes canonicalize
to square and solid-circle families. Square sizing/audit arithmetic and default
audit rows remain unchanged.

Catalog values are net material areas, including wall material only for tubes.
Density*A*physical_length determines mass; neither enclosed area nor K changes
that identity. sqrt(A) is an equivalent-square side for nonsquare families,
not an actual diameter, side or wall thickness. Each `SizedMember` retains its
admitted shape; nonsquare audit rows include its family and weak-axis coefficient.
Signed axial strengths, tensile areas and the LP/certificate problem are unchanged.

These are elastic member Euler checks. No local-wall, distortional, torsional,
inelastic, imperfection-sensitive or code buckling capacity is inferred.
Fabrication and wall-thickness limits remain explicit caller responsibilities.
The existing square Cosserat spot check refuses a nonsquare **requested member**;
other members' section choices do not disable it.

Focused material tests use independent geometric moments for all four families,
combined K/shape scaling, refitted compression signs, catalog resistance,
original indexing after pruning, unchanged tension/mass semantics, invalid
geometry and cancellation. Body-load regressions check shape-dependent converged
sections and independent original K through the actual self-weight rebuild.

## Prestress-aware isolated eigenvalue design adjoints

`elastic::modal::sensitivity::PreparedModalDesign` binds immutable borrowed
ground, load case, areas, moduli and reference inertia to one original modal
pencil. Repeated derivative queries reuse its accepted compatible equilibrium,
original free-coordinate map, consistent mass and sparse tangent. Each query
uses the native window owner and exactly one original matrix-free elastic
adjoint; no eigensolver/Krylov iteration is differentiated.

One requested window must contain exactly one inertia-counted eigenvalue, with
its original-pencil residual interval strictly inside the window. Empty windows,
multiple modes, repeated eigenvalues and unresolved endpoint overlap refuse.
Signed lambda = omega² remains signed. This provides local derivatives of an
isolated eigenvalue, not a unique derivative at a crossing or a full-spectrum
stability assertion.

For mass-normalized phi, let `a = B^T phi`, `c = (T^T phi)^2/L`,
`K_e,ff z_f = [B diag(EA/L) c]_f`, and `d = B^T z`.
With elastic engineering strain e, direct stiffness and actual compatible
preload contribute `E a²/L + E e(c-d)` to the area derivative and
`A a²/L + A e(c-d)` to the modulus derivative. Stress-free strain contributes
`-EA(c-d)`. The actual reference mass contributes `-lambda phi^T M_theta phi`,
including original off-diagonal consistent member mass, every Cartesian
translation and optional isotropic joint mass. Nodal-load and imposed-motion
covectors are `z_f` and `g_c-K_cf z_f`; they include preload redistribution.

Outputs retain original node/member order, with None for removed-member
activation derivatives. Free/support covectors retain their physical zero
components. Load variations must retain solvable compatible equilibrium when
unloaded axial mechanisms exist. The caller composes source-temperature,
self-weight/body-load, density/area and geometry-dependent load or support-motion
chain rules. Material states, global force/motion components, support flags and
topology are held fixed in the reported partials. These are derivatives of the
declared initial-stress approximation, not finite deformation, follower forces,
catalog snapping or nonlinear collapse.

A separate explicit sensitivity payload allowance covers gradient arrays and
the original compatible-adjoint workspace/history; native WindowPlan admission
still governs spectral factors/vectors. Checkpoints bound source, sparse and
recovery loops. Neither failed isolation nor failed/cancelled adjoint work
publishes a derivative. The returned native mode interval and adjoint residual
are numerical diagnostics, not error enclosures on the gradients.

G1/G3 regressions in the existing elastic_modal target cover analytic member/
joint-mass derivatives, rotated thermally stressed and settled redundant members
against independent input differences, off-diagonal consistent mass, actual
original-pencil residuals and signed unstable eigenvalues. G4 cases refuse
nonisolated windows, insufficient payload/iterations and cancellation.


### Original node-coordinate derivatives

The same `ModalSensitivity` now returns `node_coordinate_gradient` in original
node/x-y order [s^-2/m]. Member directions and reference lengths vary together.
Areas, moduli, free engineering strains, reference densities, explicit joint
masses and GLOBAL load/prescribed-displacement components stay fixed. A support
node's reference position remains a geometry variable; its derivative is not
masked by its displacement constraints. Only variations preserving connectivity,
axial-equilibrium rank and solvability are admissible. No topology activation,
support-frame rotation, deformed-geometry update or self-weight chain is inferred.

For an original member, write `r=x_b-x_a=L*n`, perpendicular unit direction `t`,
`a=n.dot(phi_b-phi_a)`, `b=t.dot(phi_b-phi_a)`, `c=b²/L`,
`d=n.dot(z_b-z_a)`, `k=EA/L`, actual force `N`, and total strain `e_total`.
With `u_t=t.dot(u_b-u_a)`, `z_t=t.dot(z_b-z_a)`, the derivative with respect
to r is `n*radial + t*angular/L`, where

- `radial=-k*a²/L-(c-d)*k*e_total-N*c/L-lambda*rho*A*inertia_work/3`;
- `angular=2*(k-N/L)*a*b+(c-d)*k*u_t-N*z_t`;
- `inertia_work=sum_axis(phi_a²+phi_a*phi_b+phi_b²)`.

Opposite endpoint contributions assemble the original node covectors in one
member pass. This differentiates axial stiffness, both direction and magnitude
of geometric stiffness, compatible force redistribution, and full consistent
reference mass. It reuses the existing single elastic adjoint and forms no
geometry-derivative matrix. The additional payload admission includes the
16-byte gradient per original node.

The existing `elastic_modal` target adds G1/G3 coordinate finite differences on
a redundant thermally strained and settled truss, rigid translation and rotation
covariance, full homothety chain rules, and a serial-member exact 2-by-2 pencil
with off-diagonal mass. Existing analytic bar tests also check joint-mass length
derivatives and the signed compressed-mode branch. These are local discrete
gradient checks, not physical calibration or global shape-optimality claims.


## Total critical-load design adjoints

`elastic::buckling::sensitivity::PreparedBucklingDesign` binds the original
ground, reference case, selected areas and Young's moduli to the existing
fixed/reference buckling analysis. It reuses both native compatible equilibria,
the original support map and the existing positive-base spectral owner.
`controls.modes` names the derivative range; one additional root is requested
as an upper guard under the same native factor/iteration admission. Guard roots
cannot be requested as derivative results.

`solve_mode` admits a positive returned root only when both available adjacent
returned positive roots have relative separation greater than
`max(1e-6, 10*(relative_residual_i+relative_residual_neighbor))`.
The residuals are the native dimensionless original-pencil residuals. This is
the existing native solid sensitivity owner's numerical branch policy with an
explicit residual margin. It does not certify eigenvalue enclosures, exclude
missed modes, or turn fewer returned roots into spectrum completeness. Repeated
or unresolved adjacent roots refuse rather than supply an individual gradient.

At the selected load lambda, use `u*=u_fixed+lambda*u_reference` and
`N*=N_fixed+lambda*N_reference`. For original mode phi, the normalization is
recomputed from original member work:
`s=-sum(N_reference*(t.dot(phi_b-phi_a))²/L)>0`.
Nonpositive or nonfinite s refuses. One original elastic adjoint solves
`K_e,ff*z_f=[B diag(EA/L)c]_f`, with `c=(T^T phi)²/L`.
The direct stiffness plus combined-state geometric/preload numerator is divided
by s. No eigensolver iteration, dense derivative matrix or third physical solve
is introduced.

Outputs retain original indexing:

- Member A, E and fixed stress-free engineering-strain partials are the original
  elastic/preload contractions divided by s; removed members return None.
- Fixed-load and reference-load covectors are respectively z/s and lambda*z/s,
  zero on supported coordinates.
- Fixed prescribed-displacement covectors include geometric-source work on
  supports and are zero on free coordinates.
- Node-coordinate derivatives include length, direction and both compatible
  preload contributions. The private original `MemberTangent` contraction is
  shared with modal design sensitivity. Supported reference nodes remain
  geometry variables. Static buckling includes no mass derivative.

All reported partials hold source properties, global force/motion components,
support flags and topology fixed except for the named input. Geometry variations
must preserve axial-equilibrium rank and compatible solvability. Geometry-
dependent loads/motion, thermal source changes, self-weight and other parameter
maps need their explicit caller-owned chain rules. This is the frozen-reference
initial-stress approximation; it is not a local Euler-capacity, finite-strain,
follower-force, postbuckling, physical-calibration or nonlinear-collapse claim.

The additional sensitivity allowance bounds original-order gradients and native
adjoint payload/history separately from the prepared analysis. Its node term is
544 bytes, member term 192 bytes, iteration-history term 16 bytes per allowed
iteration, plus 4096 bytes of fixed numerical overhead. Member/node loops poll
at bounded strides. Exhausted work, isolation or adjoint refusal, nonfinite
arithmetic and cancellation return no partial derivative.

G1/G3/G4 additions to `tests/elastic_buckling.rs` cover an analytic biased
critical factor and its distinct fixed/reference-load derivatives, all member,
load, settlement and geometry finite differences, reference-load scaling,
translation/rotation/homothety identities, both sides of a repeated-root packet,
guard exclusion, absent positive roots, work and cancellation refusal. The
existing modal tests continue to exercise the shared geometry contraction and
its separate consistent-inertia contribution.

## Objective network elasticity and nonlinear continuation

`geometric::GeometricTruss` analyzes the complete selected pin-jointed network
using current chord lengths and directions. `MaterialLayout::geometric_truss`
retains each original member's source-resolved Young modulus through the same
owner. Geometry, connectivity, original areas, material state, support flags
and authored stress-free strains remain immutable during the solve/path.

For original length L, current length l, reference area A and fixed stress-free
engineering strain e0, the admitted axial model is

`e = (l-L)/L`, `e_el = e-e0`, `N = EA e_el`, `U = EA L e_el²/2`.

Its exact current-chord tangent on relative endpoint translations is
`H = (EA/L) n nᵀ + (N/l) (I - n nᵀ)`. Tension stiffens and compression softens
the transverse response. This uses actual deformed directions and l, rather
than superposing the existing frozen-reference buckling pencil. Rigid finite
rotations and common translations preserve the scalar material strain and
energy to floating-point accuracy. Original coordinates are not overwritten.

Large rotations do not establish a finite-strain material law. The caller
must author a finite `maximum_absolute_elastic_strain` in (0,1); every trial
and reported state must satisfy it. This bound is not source calibration,
a yield/failure limit, or proof of validity of the linear axial model. e0 uses
the original reference length, and `L (1+e0)` must remain positive. No thermal
expansion coefficient, current area, Poisson contraction or evolving material
state is inferred. Reported axial stress is **nominal** N/A, not Cauchy stress.
Zero area removes that member's energy/force/tangent while retaining its
original output index; removed members have zero strain/force/stress outputs.
All original free translations remain, including mechanisms.

`GeometricLoading` separates fixed physical dead forces, fixed thermal/misfit
and fixed support displacement from reference-case dead forces and support
motion per factor. Thus `f = f_fixed + lambda f_reference` and
`u_support = u_fixed + lambda g`. Both support-load components are retained
in reactions. `GeometricScales` defines `u_free = L_scale q` and
`R = R_physical / F_scale`; the exact tangent is `K_ff L_scale/F_scale`, and
the parameter column is `-dR/dlambda = (f_reference,f - K_fc g)/F_scale`.
The deformed K_fc term is essential for displacement-controlled paths.

`equilibrate` uses native `fs-solver::nonlinear::NewtonKrylovState` with its
original globalization, cancellation and iteration budget. Right Jacobi
preconditioning changes numerical coordinates only; zero diagonal rows gain
no spring. Accepted responses independently re-evaluate the physical member
force equation against the requested dimensionless residual target.
`start_from` admits an equilibrated complete physical displacement, including
exact current support values. A fixed thermal/gravity baseline can first use
`equilibrate`. `advance` uses the existing `fs-solid::continuation` Keller
bordered corrector and step control to traverse simple load-factor folds.
It stages a clone: any requested-step failure or cancellation returns no
candidate and leaves the supplied checkpoint unchanged. Resume with the same
immutable problem and scaling; public native path values do not authenticate
problem identity. Starting/resuming rechecks the original force equation.

The graph's canonical CSR pattern and each member's row scatter offsets are
prepared once with `fs-sparse`. Warm force/tangent assembly is O(nodes + members
+ stored entries), without repeated COO sorting. Native residual, tangent and
parameter callbacks at an identical state share one physical assembly. No
new nonlinear solver, Krylov recurrence or continuation algorithm is introduced.
The preparation admits at most 800000 original coordinates, 100000 free
coordinates and 8000000 staged contributions under explicit caller limits.
The conservative payload allowance is `512*original_dofs + 256*members +
128*staged + 2048*free_dofs + 65536` bytes, plus `32*(free_dofs+1)²` when the
native border uses its bounded dense branch. Solver/path histories reserve
512 bytes per admitted entry. Native Newton restart/cycles are bounded by
64/320; outer attempts by 10000. Path advances admit at most 1000000 steps
and 128 corrector iterations each. Malformed checkpoint vector lengths refuse
before cloning. Node/member/operator loops poll cancellation at bounded strides;
the one bounded native COO sort is checked immediately before and after.

Physical responses retain original displacements, signed member forces and
nominal stresses, elastic energy, supported reactions and **all** free force
residuals. A raw response does not assert equilibrium. The reported potential
`U - f dot u` holds the parameter fixed; it is not accumulated work along
moving-support/loading histories. Numerical equilibrium and path events do not
prove stability, missed branches, collapse, local Euler capacity, out-of-plane
behavior, buckling of an individual member, contact, plasticity or bending.

G1/G3/G4 tests in `tests/geometric_truss.rs` cover analytic loaded thermal
extension with prescribed support motion and forces on supports; independent
finite differences of the full coupled tangent, energy and parameter column;
finite rigid rotations with free thermal expansion; the closed-form shallow
two-bar force curve through both folds; bit-identical split path resumption;
collapsed/domain-invalid members, removed mechanisms, exhausted work/payload,
unbalanced seeds, malformed checkpoints and cancellation.
