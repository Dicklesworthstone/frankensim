# fs-solid CONTRACT

## Purpose and layer

Layer: L3 (FLUX). The elasticity core (plan §8.2, bead tfz.13):
small-strain elasticity through finite-strain hyperelasticity, with a
MEASURED locking-free formulation, on both the body-fitted and the
CutFEM-on-SDF frontends. Constitutive laws live in fs-material; this
crate owns body-fitted kinematics, weak forms, boundary operators, and
the globalized Newton loop. Vector CutFEM assembly and solution are
owned canonically by fs-cutfem; fs-solid retains only a legacy
constructor facade.

The `linear3` module is the generic body-fitted three-dimensional
small-strain rung. It assembles a sparse stiffness/consistent-mass pencil from
P1 tetrahedra plus uniform or per-element `fs-material` state points. Its
constitutive seam is a finite, exactly symmetric, positive-definite 6x6 tensor
in an explicitly documented Mandel basis; isotropic and oriented-orthotropic
constructors lower into that shared representation. Object names and material
presets are deliberately absent: downstream vibration and vibroacoustic
consumers receive the same `(K,M)` artifact for any admitted
geometry/material composition. Anisotropy is therefore supported at the
operator level, but material-card resolution of engineering/crystal constants
and texture/orientation fields remains a separate evidence-bearing admission
step; absence of that data is not permission to assume isotropy.

`TetThermalStrainState` and `assemble_thermal_load` add the matching generic
small-strain thermomechanical load boundary. An upstream material law supplies
the already-integrated, possibly anisotropic Mandel thermal strain between an
explicit positive reference and current temperature. Each state is bound to
the exact elastic material-state identity and a nonzero thermal-law identity;
the isotropic binder consumes fs-material's exact signed
`integral(alpha(T) dT)` only when its current endpoint and tangent elasticity
share the exact card and complete state point. Otherwise it refuses before
assembly.
the operator assembles `integral B^T C epsilon_free dV` on the same tetrahedra
and publishes both full nodal forces and the constrained free-DOF projection.
`solve_thermal_displacement` advances a Jacobi-PCG solve in cancellation-bounded
chunks and accepts only an independently recomputed true Euclidean residual;
it reconstructs the full zero-Dirichlet nodal displacement field for an
explicit later geometry-update transaction. That transaction is
`update_geometry_from_displacement`: it identity-checks the reference mesh,
enforces a caller displacement ceiling, preserves connectivity, and refuses
every degenerate or orientation-flipped tetrahedron before publishing the
updated fixed-topology chart.
No material name, constant-alpha assumption, or Euler-disc geometry appears in
this layer.

## Public types and semantics

- `Mesh2` / `Patch`: structured body-fitted meshes — P1 triangles, Q1
  quads, mapped Q1 panels (`cooks_membrane`), with the four sides as
  named patches. `shapes_at`/`quad_points` expose the isoparametric
  machinery (chain rule with the transposed inverse Jacobian — the
  conformance battery pinned this exact slip on triangles and mapped
  panels).
- `LinearProblem` (plane strain/stress via `PlaneKind`, `lame`):
  small-strain elasticity with strong Dirichlet and dead-load traction
  patches; `Formulation::Standard` vs `Formulation::BBar` (Hughes'
  dilatation projection — element-average volumetric gradient,
  equivalent to condensed constant-pressure mixed elements).
- `HyperProblem` + `NewtonSettings`/`NewtonReport`: plane-strain
  finite strain; the 2D displacement gradient embeds into the 3D
  deformation gradient (F₃₃ = 1); residual `∫ P : ∇δu` and consistent
  tangent `∫ ∇δu : A : ∇Δu` come from fs-material's exact AD Piola
  stress and 9×9 tangent. Newton globalizes by Armijo backtracking on
  the potential energy (material refusals — det F ≤ 0 — halve the
  step) under linear load stepping; systems solve with MINRES
  (symmetric, indefinite-safe near buckling). `residual_and_tangent` /
  `potential_energy` are the public consistency probes.
- `CutElasticity`/`CutSolution`: source-compatible legacy facade over
  `fs_cutfem::CutElasticity` and its canonical solution type. The facade
  constructs `IsotropicElastic::new(E, ν, 1.0)`; the fixed strain limit
  supplies compatibility metadata only and is not evidence that a solved
  field remains inside that domain. Historical stabilization knobs retain
  their physical coefficients by translating
  `β_canonical = β_legacy(λ+2μ)/μ` and
  `γ_canonical = γ_legacy(λ+2μ)/μ` before delegation. Canonical assembly,
  graded-tree constraints, Nitsche and ghost terms, deterministic CG,
  retained topology, and error integration are not reimplemented here.
  The facade fixes the legacy solver controls at tolerance `1e-12` and
  `60_000` iterations. The removed implementation's post-hoc
  `1e-8` acceptance after targeting `1e-12` is not reproduced: a miss
  of the canonical `1e-12` gate fails closed. `BoundaryTraction`,
  `DesignBoxEdge`, and opaque checked `EdgeBand` are direct fs-cutfem
  re-exports. `solve_with_boundary_traction` delegates the typed support
  descriptor without inspecting or reimplementing it; the legacy callback
  field must be `None`, and an SDF crossing through supported load still
  refuses canonically.
  Physical-coefficient equivalence is claimed only when each active term's
  dimensionless translated coefficient is finite. Canonical validation refuses
  an overflowing Nitsche translation when `traction_free_interface` is false
  and any overflowing ghost translation. When `traction_free_interface` is
  true, the unused translated Nitsche value is deliberately ignored and
  carries no physical-coefficient claim. A future physical-coefficient input
  surface is required to represent active terms without a non-representable
  intermediate.
- `select_formulation(RegimeIndicators)`: the documented
  element-selection policy — B-bar when ν ≥ 0.45 or slenderness ≥ 5
  (fs-regime's indicators feed this; the locking battery measures the
  boundary).
- `stability` (bead tfz.15): `reduced_pencil` (free-DOF K and
  geometric stiffness K_G(σ₀) from the linear prebuckling state),
  `buckling_loads` — the pencil `(K + λK_G)φ = 0` solved MATRIX-FREE
  by fs-la LOBPCG on the Cholesky-reduced operator L⁻¹(−K_G)L⁻ᵀ
  (dense reduction fixture-gated at n ≤ 4096), modes −K_G-normalized;
  `lambda_indicator` (Richardson) + `evidence_row`;
  `eigenvalue_derivative` (direct pencil derivative at frozen
  prebuckling stress), `group_stiffness` (per-group Young's-scale
  lever), and the clustered-eigenvalue policy: `ks_aggregate` /
  `ks_aggregate_derivative` — the conservative smooth KS lower
  envelope with softmax-weighted derivatives, because individual
  eigenvalues are not differentiable where branches cross (measured
  in the battery).
- `continuation` (bead tfz.15): `PathResidual` (generic equilibrium
  residual R(u,λ) = R_int − λF_ext; `HyperProblem` implements it for
  homogeneous clamps + ramped dead tractions), `advance` — Keller
  bordered pseudo-arclength with MINRES/dense-LU tangent solves
  (indefinite-safe through limit points), step halving/growth,
  limit-point events (λ̇ sign flips) and branch-point candidates
  (Cholesky definiteness flips without a λ̇ flip, gated n ≤ 1024);
  `PathState` — the checkpointable trajectory (clone to checkpoint,
  hand back to resume; split runs are BITWISE identical to straight
  runs); `switch_branch` — null-direction predictor pinning
  (`pending_switch`): the arc constraint forces the first step along
  the buckling mode so the corrector lands on the bifurcated branch
  (a state-only perturbation relaxes back to the fundamental branch;
  a basin-scale jump inverts elements — both measured failure modes).
- `rod` (bead tfz.14): geometrically exact Cosserat rods — nodal
  positions + unit quaternions updated MULTIPLICATIVELY through
  fs-time's exponential map; strains Γ = Rᵀr′ − e₁ and
  κ = 2log(qᵢ⁻¹qᵢ₊₁)/L₀ from RELATIVE quantities (rigid motions give
  exactly zero strain — energy-invariance battery). Statics: FD
  residual/tangent Newton with residual backtracking and load
  stepping, fixture-gated; fixture slenderness matters — a
  penalty-stiff EA puts a quadratic spurious-stretch wall around
  every Newton iterate (measured, documented in the battery).
  Analytic tangents and SE(3) DYNAMICS under fs-time symplectic
  integrators are recorded successor scope.
- `fiber` (bead tfz.14): fiber sections over fs-material `Uniaxial`
  cards (Menegotto–Pinto steel; Mander concrete with the
  compression-positive convention flipped at the section boundary);
  `respond`/`commit` split (pure response vs state commitment);
  `update_sections_batched` — the §15.2 hot loop through fs-la's
  batched Cholesky with a per-outlier direct fallback, throughput
  measured and ledgered; `rc_section` fixture (confined core +
  unconfined cover + steel layers).
- `beamcol` (bead tfz.14): force-based (Spacone) cantilever element —
  EXACT equilibrium interpolation, five-point Gauss–Lobatto section
  integration (a point AT the hinge), per-section Newton with
  residual backtracking, virtual-work tip deflection; pushover and
  cyclic histories with committed, cloneable state (G4).
- `koiter` [F], feature `koiter-asymptotics` (off by default): FD
  energy expansion along the buckling mode at the critical state →
  a/b coefficients and `Bifurcation` classification, with the
  SAMPLED-CONTINUATION fallback oracle (the imperfect-geometry path)
  cross-checking imperfection tolerance.
- `contact` [F], feature `contact` (off by default): SDF-native
  log-barrier contact around `HyperProblem`, with `Barrier`,
  optional lagged smoothed-Coulomb `Friction`, `ContactProblem`,
  `ContactSolution`, and the frictionless `translation_gradient`
  adjoint probe. The solver requires a finite, intersection-free
  initial state; bad parameters, non-finite SDF values/gradients,
  invalid pins, and unsupported contact-adjoint inputs return
  `SolidError::InvalidInput`.
- `accept_scenario_bc`: fs-scenario integration — `Dirichlet` and
  `Traction` under `Physics::Elasticity` accepted; dimensioned-value
  resolution stays with the scenario consumer.
- `shell` (bead frankensim-b8bxd.9.1): a safe-Rust, deterministic,
  flat, single-patch triangular **estimate-only surrogate** for the Euler-disc
  base. It combines lumped inertia, edge-spring in-plane resistance, and
  normal-slope/rotation penalties; it is not a CST plane-stress membrane or a
  Kirchhoff--Love/isotropic continuum discretisation.
  `ShellPlate::assemble` returns type-separated dense
  `ShellMatrix<Mass>`, `ShellMatrix<Stiffness>`, and optional Rayleigh
  `ShellMatrix<Damping>` operators in node-major `(x,y,z,rx,ry,rz)` order.
  Each triangle contributes edge-spring, normal-slope, and rotation-gradient
  surrogate energy; translational and rotary inertia are lumped from the
  material card. `ShellSupport` applies three distinct, normal-aligned point
  pinned supports by reducing all translations at those nodes while retaining
  rotations; its normal identifies the plate side and it does not claim a
  frictional or rotational support law. Inputs carry model,
  source, and state IDs and fixed SI/right-handed-Cartesian interpretation.
  Assembly bounds node, element, dense-memory, deterministic
  validation-and-assembly work, and diagnostic sizes before allocation; it
  returns no partial operator. Validation charges finite-node scans,
  scale-relative coincident-node pairs, canonical face insertion, oriented
  edge incidence/adjacency, connectivity traversal, used-node marking, and
  triangle assembly against the declared work budget. Small assemblies report
  raw algebraic mass-array definiteness, stiffness-array nullity, and symmetry
  residue via a deterministic Jacobi eigensweep. Its raw mixed-unit spectral
  spread is diagnostic only, not a physical condition number or continuum
  certificate; larger valid assemblies explicitly report that those bounded
  algebraic diagnostics were not computed.
  Scale-relative coincident nodes, duplicate index faces, geometrically
  duplicated faces through coincident nodes, non-manifold edges (more than two
  incident faces), same-orientation shared edges, and components disconnected
  by an edge are refused before assembly; callers requiring several components
  must model separate plates.

- `reduce` module (music bead 3ez8g.3.3): the offline lip/reed reduce
  lab. `reduce_valve(request, cx) -> ValveCard` runs
  `TetLinearElasticProblem` -> `fs_modal::slice_window` -> per-mode
  re-normalization from the solver's `phi^T M phi = 1` convention to
  UNIT MEAN TRANSLATION of the orifice face (mass normalization would
  destroy the units), minting physical `m_eff = phi^T M phi` [kg],
  `k_eff = lambda m_eff` [N/m], and `c = eta sqrt(k m)` from the card's
  AUTHORED loss factor (the lab does not invent Q values). The
  retained-compliance law: the modal compliance
  `sum (phi^T r)^2 / lambda` of the retained modes — `phi^T r` is the
  FORCE-pattern projection; the first build used the base-excitation
  `phi^T M r` and retained 1e-11 of the compliance, an executed
  falsification — is compared against the EXACT static compliance
  `r^T K^{-1} r` (sparse LDLT), and a reduction below the authored
  floor REFUSES (`RetainedComplianceTooLow`). Orifice geometry is
  MEASURED from the mesh: the face is the surface strip facing the
  opposing plane, gap = closest approach, width/channel-thickness =
  face extents. The card carries the provenance chain (mesh digest,
  material state-point receipt identity, source id) under a
  domain-separated schema-versioned content identity with
  canonical-bytes round-trip and tamper refusal; `beating_reed()`
  mints the `fs_scenario::BeatingReed` with the runtime's own face
  convention `P_c = k H / (w * 0.025)` so the massless and massive
  branches agree.

## Invariants

1. Patch tests exact: linear displacement fields are reproduced to
   solver tolerance on P1, Q1, B-bar, and the CutFEM frontend
   (all-embedded fixture) — sol-001.
2. MMS optimal orders per implemented family × frontend: L2 ≈ 2,
   H1 ≈ 1 for P1/Q1/CutFEM-Q1 (sol-002).
3. Objectivity through the fs-material interface: W(RF) = W(F) and
   P(RF) = R·P(F) for both cards (sol-003).
4. Consistency merge gate: the assembled tangent equals the FD
   directional derivative of the residual, and the residual is the
   exact gradient of the potential energy (rel ≤ 1e-6) — sol-004.
5. Locking is MEASURED, not assumed: on thin near-incompressible
   bending, the standard element's error degrades by ≥10× from
   ν = 0.3 to ν = 0.4999 while B-bar stays within 2× (sol-005).
6. Cook's membrane envelope: converged tip at the literature reference
   point (48, 52) inside [23.5, 24.5] (plane stress, ν = 1/3), ≤2%
   self-convergence deviation; near-incompressible B-bar
   self-converges while the standard element's shortfall is logged
   (sol-006).
7. Newton robustness: large-strain bending and a buckling-adjacent
   compressed strip converge under load stepping + line search with
   fast terminal contraction; histories and backtracks are evidence
   (sol-007).
8. Determinism: BTree/insertion-ordered assembly, deterministic
   Krylov solvers; bit-identical across runs on a platform —
   including continuation trajectories (checkpoint/resume proven
   bitwise-equal in stab-002).
9. Buckling: the Euler strut's EXTRAPOLATED critical load lands
   within 3% of the analytic value with the Richardson indicator
   covering the raw fine-mesh gap (Q1 parasitic-shear inflation,
   reported, stab-001); pencil symmetry to 1e-10 and mode
   K-orthogonality to 1e-6.
10. Continuation: validated against the CLOSED-FORM von Mises truss
    (limit points within 5–10% of analytic, manifold deviation
    ≤ 1e-7, load control jumps where arclength traces, stab-002);
    the continuum tent traces the full snap-through Z-curve on one
    path — snap and snap-back limit points, load recovery beyond the
    load-control failure point (stab-003).
11. Branch handling: bifurcation detected within 15% of the pencil
    prediction on the compressed column; null-direction switching
    lands on the bent branch (transverse growth ≥ 20× measured,
    stab-004).
12. Flat plate-surrogate energy consistency: each assembled stiffness contribution is a
    positive weighted outer product, so the free operator is symmetric
    positive-semidefinite while positive lumped translational/rotary inertia
    makes the raw mass array positive-definite. A single free triangular panel has the six
    infinitesimal rigid motions in its stiffness null space; the focused
    `shell_operator` battery verifies this and independently reconstructs the
    quadratic energy from `Kx`.

## Error model

`SolidError` teaching errors: `SolveFailed` (linear gate missed),
`NewtonStalled` (carries the residual history; repair = more load
steps), `MaterialRefused` (det F ≤ 0 escaped globalization; repair =
smaller steps), `UnknownPatch`, `UnsupportedBc` (fs-scenario mapping
outside Dirichlet/Traction × Elasticity), and `InvalidInput`
(non-finite/inconsistent caller data such as invalid contact
parameters, infeasible starts, or malformed adjoint seeds).
`InternalInvariant` distinguishes a corrupted lower-level solver graph from
repairable caller input.

The CutFEM facade maps only `CutFemError::SolveNotConverged` to
`SolidError::SolveFailed`, preserving its iteration count and relative
residual. A canonical `ConstraintCycle` maps to `InternalInvariant`; every
material, parameter, geometry, callback, topology-admission, and
unsupported-regime refusal, including an invalid typed edge band or a cut
supported segment, maps to `SolidError::InvalidInput` with the canonical
diagnostic. No canonical refusal is converted into a solution.

## Determinism class

Bit-deterministic across runs on a fixed platform (ordered assembly,
fs-solver CG/MINRES, no threading, no ambient state). CutFEM determinism
and topology ordering are inherited directly from fs-cutfem rather than
duplicated in this crate. Cross-ISA golden hashes are not yet recorded
(follow-up).

## Cancellation behavior

Bounded synchronous loops (assembly, capped Krylov iterations, capped
Newton with capped backtracking). Chunked Cx polling belongs to the
fs-exec driver (the L3 discipline; fs-feec/fs-cutfem precedent).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`koiter-asymptotics` (off by default): the plan-flagged [F] Koiter
post-buckling module and its battery, per the Ambition-Tag gating
rule.

`contact` (off by default): the plan-flagged [F] SDF-native barrier
contact prototype and its battery, per the Ambition-Tag gating rule.

### Structural invariants (tfz.14)

12. Rod objectivity: superposed rigid motion leaves strain energy
    invariant to 1e-12 (str-001).
13. Elastica: large-deflection tip matches an in-test SHOOTING oracle
    (RK4 + bisection on the elastica BVP) within 0.015·L at
    PL²/EI = 2 (str-002).
14. Helical family: pure bending gives the EI/M circle (2%), pure
    torsion stays straight (1e-6) at rate M/GJ (2%), combined moments
    give the constant helix strain state (3%) — str-003.
15. RC hysteresis: per-cycle dissipation positive and growing with
    amplitude, peak moment inside the hand capacity band, section
    states bitwise-deterministic across runs (str-004).
16. Batched section updates are scalar-consistent to 1e-8 with
    ledgered throughput (str-005).
17. Force-based pushover softens past yield (secant < 0.5× initial),
    dissipates under reversal, and resumes from a checkpoint
    BITWISE (G4) — str-006.
18. Contact barrier calculus and solve hygiene (feature `contact`):
    barrier derivatives match finite differences, invalid/penetrating
    contact inputs return structured `InvalidInput` instead of
    panicking, accepted iterates keep positive gap, reactions balance
    the applied normal load, tight-clearance and curved-obstacle
    fixtures retain positive gap, lagged friction separates stick/slip
    cases, and the frictionless translation adjoint matches FD on a
    plane — cnt-001..007.

11. Reduce-lab gates (rl-001..rl-004): the thick cantilever-slab
    oracle (f above the slender Euler-Bernoulli line by the plate-strip
    factor, measured +11.4% inside an authored 0-35% band;
    strip-normalized effective-mass ratio 0.355 in (0.30, 0.42);
    `k = lambda m` at machine precision; the
    `k_ratio = f_ratio^2 * mass_ratio/0.25` self-consistency identity;
    retained compliance 0.98 > 0.95), the minted card driving the
    existing fs-phs valve island to its quasi-static equilibrium at
    1e-3 with supply defect < 1e-8, refusals (empty/out-of-range face,
    the retained-compliance floor via a window that excludes the
    compliant fundamental, a face behind the opposing plane), and the
    canonical-bytes round-trip + tamper refusal.

## Conformance tests

`tests/battery.rs`: sol-001 patch tests (body-fitted max nodal error
≤ 1e-9; CutFEM at the CG-tolerance floor ≤ 1e-7); sol-002 MMS orders
(P1 2.00/1.00, Q1 2.00/1.00, CutFEM 2.19/1.04); sol-003 objectivity
(energy rel ≤ 1e-10, stress rotation ≤ 1e-8, 8 random states × 2
cards); sol-004 tangent/energy consistency (≤ 1e-6); sol-005 locking
battery (standard degrades ≥10×, B-bar ≤2×); sol-006 Cook's envelope;
sol-007 Newton robustness (5-step large bending, 8-step compression,
terminal contraction ≥ 100×); the selection-policy probe. Unit tests:
mapped-mesh areas, selection thresholds. The CutFEM ownership probes
compare facade and direct canonical coefficients, nodal values,
residuals, iteration counts, active cells, and error norms bitwise after
the documented stabilization translation. A typed right-edge band likewise
compares coefficients, nodal values, and assembled-load compliance bitwise and
proves that an occupied legacy callback slot refuses. The probes also require
invalid material data, the uncertified near-incompressible regime, and an
empty domain to fail closed as `InvalidInput`.

`tests/stability.rs` (bead tfz.15): stab-001 Euler pencil (analytic
G2 + Richardson indicator); stab-002 von Mises truss closed-form
continuation oracle (+ bitwise resume); stab-003 continuum tent
snap-through with the self-calibrating load-control-failure probe;
stab-004 branch detection + null-direction switching; stab-005
eigenvalue-derivative gradient gate (rel ≤ 1e-3 vs frozen-K_G FD) and
the clustered-eigenvalue trap (min() kink ≈ 2 measured, KS aggregate
smooth with derivative matching FD to 1e-6, conservative).
`tests/koiter.rs` (feature-gated): stab-006 symmetric-stable
classification + imperfect-geometry sampled oracle.

`tests/contact.rs` (feature-gated): cnt-001 barrier calculus and
force divergence; cnt-002/003 intersection-free hard press and
reaction balance; cnt-004 tight-clearance and curved-obstacle gap
audits; cnt-005 lagged friction cone separation; cnt-006
frictionless contact translation adjoint versus FD; cnt-007
structured invalid-input returns.

`tests/structural.rs` (bead tfz.14): str-001 objectivity; str-002
elastica vs shooting oracle; str-003 circle/twist/helix; str-004 RC
moment-curvature hysteresis + determinism; str-005 batched
consistency + throughput ledger; str-006 force-based pushover,
reversal dissipation, G4 bitwise resume.

## No-claim boundaries

- `linear3` claims infinitesimal-strain elasticity and work-conjugate thermal
  eigenstrain loading on a fixed, conforming, non-degenerate tetrahedral mesh.
  The thermal operator does not infer expansion coefficients, solve a transient
  response, remesh, establish mesh/stress convergence, or claim finite-strain
  validity. Its coordinate update is fixed-topology and inversion-gated, not a
  finite-strain constitutive update. Its static solve establishes only its reported linear
  residual on the supplied discretization. It does not claim plasticity,
  physical damping, phase change, fracture, moving interfaces, or topology
  evolution. Callers must escalate or refuse when the resolved material state
  leaves this solid rung.
- Its free-body stiffness correctly contains rigid modes. Removing, fixing,
  or slicing around those modes is a boundary-condition/eigenproblem decision;
  the assembler never adds hidden grounding stiffness.

- The facade inherits fs-cutfem's no-claim for clipping a supported traction
  segment cut by the SDF and for quadrature error bounds on arbitrary traction
  callbacks. It only delegates checked support; it does not certify callback
  smoothness or grow a second boundary assembler.
- Extreme `(E, beta, gamma)` combinations whose active physical stabilization
  products are finite but whose fs-cutfem dimensionless coefficients are not
  representable are not claimed. Nitsche-active and ghost-active terms fail
  closed pending a canonical physical-coefficient API; an unused Nitsche value
  in traction-free mode is ignored rather than validated.
- TDNNS-proper / weakly-imposed-symmetry FEEC stress elements: awaits
  the simplicial H(curl)/H(div) family bead (dcng); B-bar is the
  shipped, measured locking-free path with the same acceptance metric.
- 3D and shells; higher-order families (Qk/Pk, k ≥ 2); IGA frontend
  (bead tfz.9).
- Hyperelasticity on the CutFEM frontend (linear ships there; the
  finite-strain cut path composes this crate's Newton loop with
  fs-cutfem quadrature in a successor bead).
- The legacy CutFEM facade makes no independent coercivity,
  near-incompressible, quadrature, graded-tree, convergence, or
  cancellation claim. In particular, material ratios beyond
  fs-cutfem's certified compressible regime are refused rather than
  routed to fs-solid's body-fitted B-bar formulation.
- fs-opdsl-generated residual/adjoint paths (the DSL coverage bead);
  the hand path here passes the same consistency gates the DSL output
  must pass.
- Contact no-claims: mesh-mesh point-triangle fallback, rigorous
  interval CCD over curved SDFs, adaptive barrier-stiffness policy,
  frictional adjoints/VJPs, dynamics/time stepping, and fs-scenario
  contact-law coefficient wiring. Plasticity flow remains successor
  scope (fs-material's plastic module wires in later).
- Ogden energies (fs-material's recorded no-claim propagates).
- Rod SE(3) DYNAMICS (fs-time symplectic wiring) and analytic rod
  tangents; 3D frames/lattices assembled from rods (topo-truss's
  consumer pass); J2 continuum element wiring (the fs-material card
  ships with algorithmic moduli; the element battery is successor
  scope); published El Centro fiber-model plot digitization (the
  shipped envelope is the capacity band + self-consistency).
- fs-scenario dimensioned-value resolution (units plumbing stays with
  the scenario consumer; kind/physics validation ships here).
- Production-scale pencil solves (shift-invert/sparse factorizations
  beyond the fixture-gated dense reduction) and the
  prebuckling-adjoint chain dλ/ds through σ₀(s) (fs-ad/ASCENT
  integration successor); LOBPCG preconditioning hooks.
- Plate/shell buckling references (needs bending elements — the
  shells bead); the 2D canonical here is the strut.
- CST plane-stress patch accuracy, Kirchhoff--Love/isotropic-continuum claims,
  physical condition-number authority, curved shells, multi-patch/mortar joins,
  NURBS geometry, transverse-shear
  (thick/Reissner--Mindlin) behavior, nonlinear/finite-rotation response,
  loads, time integration, constraint multipliers for oblique support normals,
  shell error enclosures, and NAFEMS qualification. The `shell` operator is a
  flat, linear spatial base only; fs-iga's dedicated surface-shell frontend
  remains the required successor for those claims.
- Koiter coefficients versus LITERATURE tables (v1 validates
  classification + the sampled oracle; quantitative a/b benchmarking
  awaits the shell fixtures where the classic tables live).
- The reduce lab is Estimate-authority end to end until measured
  tissue/cane packs exist: the material card in the battery is
  synthetic-authored (no soft-tissue elastic pack exists in
  data/matdb — the recorded upgrade path), the loss factor is an
  authored card field, and the "lip-shaped" fixture is an authored
  tapered slab. The lab claims reduction correctness, never tissue
  truth; claim-side promotion happens in the per-track gates beads.
