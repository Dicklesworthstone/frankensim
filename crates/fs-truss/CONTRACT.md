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
