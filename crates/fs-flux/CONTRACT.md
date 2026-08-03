# fs-flux CONTRACT

## Purpose and layer

Layer: L3 (ORGANS). FEEC-native incompressible Navier–Stokes (plan
§8.3 [F], bead tfz.17): H(div)-conforming BDM1 velocities against P0
pressures, giving EXACTLY divergence-free discrete velocities and
therefore PRESSURE-ROBUSTNESS — velocity errors independent of the
pressure, the correctness property most production CFD lacks. 2D
triangle-mesh instantiation: interior-penalty viscosity, upwinded DG
convection, Picard steady solves, IMEX BDF1 transients, discrete
adjoints.

## Public types and semantics

- `TriMesh` (`trimesh`): triangles + a GLOBALLY ORIENTED edge table
  (`Edge { verts a<b, tris, len, normal, mid }`, normal owned by the
  lower-index adjacent triangle; `tri_edges[t][k] = (edge, ±1)` with
  local edge k OPPOSITE vertex k). Built from `fs_solid::Mesh2`
  triangles (CCW). RT0 helpers retained as utilities.
- `audit_affine_triangle_gcl2` (`ale`): a fixed-connectivity 2-D ALE
  geometry audit for triangles whose vertices move linearly through one
  explicit time step measured in coherent SI seconds. It validates the
  public `TriMesh` incidence tables, rejects endpoint and mid-step
  collapse, integrates owner-oriented swept edge areas, and returns
  per-cell plus global GCL defects without rejecting a finite mismatch
  or certifying itself.
- `reduced_aero`: a deterministic exterior-flow screening seam for a moving
  circular disc, suitable for a spinning tilted disc without importing the
  Euler-disc domain. `ReducedAeroInput` names one world frame and uses coherent
  SI fields: metres, seconds, kilograms, newtons, and newton-metres. A
  `ReducedAeroModel` separately exposes quadratic translational/form drag,
  face rotational skin friction, exterior-rim edge flow, and disc-axis-rate
  damping. `CandidateWrench` is force/moment **on the body** about the supplied
  reference point; it retains every component rather than hiding a thin-gap or
  fitted term in a total coefficient. Gas, roughness, and caller-supplied
  coefficient-uncertainty card identities flow into the `EstimateOnly` receipt.
  `WorkWindow` refuses duplicate caller keys,
  providing one domain-level exactly-once work receipt before the coupling
  driver supplies its chosen closed accounting chart.
- `bdm` module: `cell_basis` builds the per-cell BDM1 basis (all of
  P1², 6 dofs = mean + signed-arclength normal moments per edge,
  BOTH against the global edge normal and orientation — orientation
  consistency for free) by inverting a 6×6 dof-Vandermonde (fs-la
  LU); constant per-cell gradients and divergences cached. RT0 was
  measured rank-deficient for the viscous form (per-cell gradient
  c·I) and rejected — recorded in the bead trail.
- `FluxSystem` (`ns`): assembly + solvers. Dirichlet BCs split by
  component: u·n STRONG (boundary-edge BDM dofs ARE the normal
  moments — identity rows), tangential weak via symmetric interior
  penalty; every remaining basis function has v·n ≡ 0 on the whole
  boundary by dof duality, so the ∮p(v·n) consistency term vanishes
  identically. Pressure pinned by REPLACING cell 0's continuity row
  (its div-freeness is implied by the other cells plus the prescribed
  net flux — never an additive penalty, which would corrupt it).
  Convection: upwinded DG on the SINGLE-VALUED face flux w·n (H(div)
  conformity makes upwinding well-posed), volume term in
  divergence form −∫(u⊗w):∇v. `picard` for steady NS, `bdf1_step`
  for IMEX transients (implicit Stokes, lagged convection),
  `solve_adjoint` for the transposed system. Dense LU under n ≤ 1500
  (fixture path), GMRES restarts above.

## Invariants

1. BDM1 basis reproduces its defining dofs (Kronecker to 3.7e-15),
   satisfies the divergence theorem per basis function (1.4e-15), and
   has single-valued normal traces across interior edges (3.6e-15)
   (flux-001).
2. Stokes MMS: velocity L2 slope 1.61 on h = 1/8 → 1/12 (theoretical
   2, preasymptotic march 1.14/1.40/1.61 measured over h = 1/4..1/12;
   asymptotic confirmation is perf-lane scope, LEDGERED); pressure L2
   slope 0.72 toward 1; per-cell divergence ≤ 1e-12 on EVERY mesh
   (flux-002).
3. Pressure-robustness is a discrete IDENTITY, not an estimate:
   adding a 1e4-amplitude gradient forcing A·∇φ changes the velocity
   by 1.8e-12 relative and shifts the pressure by exactly
   A·(Π₀φ − Π₀φ|cell0) to 1.4e-12 against a scale of 547 (flux-003).
   The identity requires the load quadrature to integrate ∇φ·v_h
   EXACTLY (degree-4 φ here; a degree-8 φ measurably leaked 1.7e-4).
4. Lid-driven cavity Re=100: Picard contracts (13 iterations to
   7.5e-10), the NONLINEAR solution is div-free to 1.2e-13 per cell,
   and the vertical-centerline u_x matches Ghia–Ghia–Shin at three
   stations to worst deviation 0.029 on an 8×8 mesh (HONESTY BAND
   ±0.15; fine-mesh table comparison ledgered to perf lanes)
   (flux-004).
5. Taylor–Green IMEX BDF1: Richardson temporal ratio 1.93 (order 1),
   kinetic energy tracks exp(−4π²νt) to 0.4%, div-free to 9e-14 at
   every step (flux-005).
6. Discrete adjoint dJ/dν matches central FD to 1.0e-4 relative
   (linear-in-ν operator, exact matrix derivatives); a GRADIENT
   weight functional is annihilated to 1.1e-10 — pressure-robustness
   read backwards: ∫∇χ·u_h = 0 for exactly div-free u_h with zero
   boundary flux (flux-006).
7. Affine-triangle ALE geometry uses one deterministic edge receipt per
   global edge. Each admitted interior neighbor traverses its shared edge
   in reverse, so independently evaluated directed sweeps cancel, while
   `tri_edges` signs bind the owner sweep to each cell and each cell measures
   `Δarea - Σ signed_edge_sweep`; a quadratic-in-time signed-area check
   refuses a collapse between two otherwise valid endpoint meshes
   (flux-007, G0/G3).
8. `reduced_aero` is a passive exterior model in relative motion: for a
   stationary ambient, `F dot v + M dot omega <= 0`, so its retained relative
   dissipation is non-negative. It uses translational Reynolds (diameter),
   rotational Reynolds (`omega * radius^2`), relative roughness, and tip-Mach
   gates as explicit inclusive correlation domains. Moving ambient flow retains
   the identity `body_power = relative_power + F dot ambient_velocity`, but it
   requires a total gas-plus-boundary accounting window rather than asserting
   a sign for body power. Zero density/motion are exact zero limits; forbidden
   thin-gap-pressure and target-fitted families fail admission (aero-001
   through aero-009, G0/G3).

## Error model

Structured panics on programmer contracts (non-triangle elements,
degenerate cells, singular saddle factorizations) with teaching
messages. Solver quality is reported, not asserted: `FluxSolution`
carries the relative residual and iteration count. The ALE geometry audit
instead returns `AleGclError2` for malformed topology, non-finite input,
collapse, arithmetic overflow, and receipt allocation failure. Its finite
roundoff defect is evidence, not an admission failure.

## Determinism class

Bit-deterministic across runs on a platform: BTree edge ordering,
fixed assembly and quadrature order, deterministic dense LU / GMRES.
The affine-triangle audit also has exact same-instance replay evidence.
Cross-ISA goldens and a full G5 audit are not recorded.

The reduced aerodynamic wrench has no iteration, random sampling, or
arrival-order reduction. Alternative candidates are sorted by canonical
correlation identity, so equivalent input order replays bit-for-bit on the
same platform. Cross-ISA equality is not claimed.

## Cancellation behavior

Bounded synchronous loops (Picard caps, Krylov caps, fixed step
counts); chunked Cx polling belongs to the fs-exec driver layer.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/battery.rs`: flux-001 basis exactness + conformity; flux-002
Stokes MMS orders + exact divergence; flux-003 pressure-robustness
identity; flux-004 cavity Picard + Ghia coarse band; flux-005
Taylor–Green BDF1; flux-006 adjoint gradient + annihilation.

`tests/ale_gcl.rs`: flux-007 fixed mesh, rigid translation, affine
expansion/shear, two-cell interior cancellation, replay, translation
covariance, continuous-trajectory minima, and fail-closed endpoint,
topology, stale-cache, orientation, and finiteness cases.

`tests/reduced_aero.rs`: aero-001 zero-density/static limits; aero-002
reversal and stationary-ambient dissipation; aero-003 quadratic scaling and
source receipt mutation; aero-004 moving-ambient accounting identity;
aero-005 frame objectivity; aero-006 correlation-domain and missing-property
refusal; aero-007 alternative-correlation disagreement; aero-008 exactly-once
work; aero-009 hostile thin-gap/target-fit refusal.

## No-claim boundaries

- Turbulence: NO LES/RANS closure ships here and nothing pretends
  otherwise; the DG dissipation is a numerical property, not a
  subgrid model.
- 3D (tets, BDM on faces), BDM2+/higher-order pressures, and curved
  boundaries — recorded successors on the same dof-Vandermonde
  surface (bead dcng's high-order FEEC would unblock them).
- Projection/splitting time integrators and BDF2+ — the IMEX BDF1
  ships; the segregated variants are documented successors.
- Fine-mesh Ghia tables, Re=1000 cavity, cylinder benchmarks, and
  asymptotic-regime convergence confirmation — perf-lane scope,
  ledgered in invariants 2 and 4.
- Preconditioned iterative saddle solvers at scale (the dense-LU
  fixture path is the shipped truth; block preconditioners are
  successor work).
- Outflow/do-nothing and traction boundary conditions (Dirichlet
  ships).
- The affine-triangle GCL receipt does not remap fields, assemble or
  advance ALE Navier–Stokes, generate mesh motion, support curved faces
  or 3-D cells, couple FSI, or certify the continuum solver. It carries
  no ledger identity and the synchronous audit has no `Cx`
  cancellation lane; those remain integration work.
- `reduced_aero` is an Estimate-only, caller-coefficient screening model. It
  neither resolves Navier--Stokes/CFD nor models turbulence, wakes, thin-gap
  pressure, squeeze films, contact, ambient vorticity, or a physical device's
  target outcomes. It makes no claim that air dominates a device. Its compact
  domain gates are correlation applicability checks, not independent physical
  validation. `WorkWindow` prevents duplicate application inside one caller
  window but is not an `fs-couple` closed first-law/entropy audit; a moving
  ambient/boundary specifically needs that separate total accounting closure.
