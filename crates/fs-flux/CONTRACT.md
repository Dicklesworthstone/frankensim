# fs-flux CONTRACT

## Purpose and layer

Layer: L3 (ORGANS). FEEC-native incompressible Navier–Stokes (plan
§8.3 [F], bead tfz.17): H(div)-conforming BDM1 velocities against P0
pressures, giving EXACTLY divergence-free discrete velocities and
therefore PRESSURE-ROBUSTNESS — velocity errors independent of the
pressure, the correctness property most production CFD lacks. 2D
triangle-mesh instantiation: interior-penalty viscosity, upwinded DG
convection, Picard steady solves, IMEX and fully implicit rheological BDF1
transients, discrete
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
  coefficient-uncertainty card identities flow into the `EstimateOnly` receipt;
  a private evaluation snapshot makes `WorkWindow` refuse any later public
  correlation, source-identity, or uncertainty-width relabelling.
  `WorkWindow` refuses duplicate caller keys,
  providing one domain-level exactly-once work receipt before the coupling
  driver supplies its chosen closed accounting chart.
- `lc::step_quarter_wave`: a lossless distributed-line coordinate transform.
  `QuarterWaveParams` admits disturbance frequency [Hz], propagation speed
  along the circuit [m/s], and developed conductor length [m]. The result
  retains wavelength, quarter-wave target, signed length/phase errors,
  electrical length, length ratio, and `abs(sin(beta*l))` as a normalized
  standing-wave profile only. It does not promote that profile into volts.
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
  (fixture path, with two residual-correction solves reusing the factors),
  sparse block-preconditioned FGMRES above.
- `FluxSystem::solve_sparse` also exposes the sparse path below that cutoff.
  A velocity ILU(0) and pressure approximate Schur ILU(0),
  `C - B diag(A_u)^-1 G`, are built from the actual assembled blocks, including
  local viscosity, transient mass, convection, boundary identity rows and the
  pressure pin. Block LDU applications reuse their scratch buffers; no dense
  system or dense Schur matrix is formed. Rows and loads are equilibrated by
  their largest absolute matrix entry so large viscous momentum loads do not
  conceal constraint errors. Every restart independently recomputes both
  residuals against the original matrix. `FluxSolution::rel_residual` remains
  the original Euclidean relative residual; `scaled_rel_residual` exposes the
  row-equilibrated relative residual. `converged(tol)` requires both to pass,
  and the implicit rheology and physical material wrappers use that guard.
  Dense solves report both diagnostics after residual correction as well.
  `max_iters` is a strict Krylov iteration budget, including the last short
  restart. An approximate factorization breakdown falls back to identity
  preconditioning without changing the equations or asserting convergence.
- `FluxSystem::solve_adjoint_checked` solves the actual transposed assembled
  operator, including strong boundary rows, the pressure datum and nonsymmetric
  convection. It uses sparse counting-sort transposition, refined dense LU
  below the 1500-unknown cutoff and the same block-preconditioned FGMRES above
  it. Both independently recomputed transposed residuals must pass before a
  `FluxAdjointSolution` is returned; `FluxAdjointError` reports invalid inputs,
  detected singularity and exhausted sparse iterations. The legacy
  `solve_adjoint` uses default controls through this path and no longer has a
  dense-only size restriction. For A(theta)x=b(theta), sensitivities contract
  the returned multipliers with b_theta - A_theta*x and add any direct objective
  derivative. Parameter-dependent boundary loads belong in that derivative.
  The supplied matrix must be the intended linearization: a frozen Picard
  matrix omits constitutive and convecting-velocity derivatives and cannot
  supply a complete nonlinear rheology/Navier–Stokes gradient. Algebraic
  convergence does not establish a condition-number or gradient-error bound.
- `FluxSystem::assemble_variable_viscosity` assembles the full Newtonian
  stress `-div(2 nu D(u))`, where `D=(grad u + grad u^T)/2`. Its callback
  receives the cell index and quadrature point and returns positive finite
  kinematic viscosity. This admits spatial/temperature variation and
  discontinuous one-sided viscosity traces. The volume term is `2 nu D(u):D(v)`;
  boundary and interior consistency terms use the matching symmetric
  traction. Interior weights use the opposite side's `nu/h_normal`, and the
  penalty uses its harmonic mean, with `h_normal=2*cell_area/edge_length`.
  Material queries occur once per quadrature point, outside basis-pair loops.
  `assemble_viscous_transient` adds the existing BDF1 mass operator and load.
  The legacy `assemble`/`assemble_transient` remain the homogeneous
  vector-Laplacian discretization.
- `FluxSystem::shear_rates` evaluates the objective `sqrt(2 D:D)` invariant
  from each cell's affine BDM1 velocity. `picard_rheology` updates this rate,
  apparent viscosity and upwind convection together. A callback of
  `(cell, point, shear_rate)` can capture the temperature/material field and
  supply a Carreau or regularized power law; its finite zero-shear value
  initializes the solve. The returned iteration count and relative velocity
  update disclose budget exhaustion rather than claiming convergence.
  These APIs retain constant reference density and kinematic pressure; a
  dynamic material viscosity is divided by that SAME density throughout.
  They do not solve variable-density momentum, advect material interfaces,
  supply a thermal field, or add viscoelastic memory. For a smooth supplied
  constitutive tangent, `linearize_rheology` below supplies the state Jacobian
  needed by the checked transpose solve.
- `FluxSystem::bdf1_rheology_step` advances generalized-Newtonian flow with
  viscosity, symmetric viscous traction, and convection at the new time level.
  Its BDF1 mass load always uses the immutable previous velocity. The law,
  loads, and boundary callbacks describe the endpoint; the law can capture
  caller-supplied endpoint temperature and material fields. The previous
  velocity supplies the initial iterate. Acceptance requires both the relative
  velocity update and both original and row-equilibrated independently
  reassembled nonlinear residuals to meet the supplied tolerance. The returned
  scalar nonlinear diagnostic is the larger of those two residuals. This extra
  residual assembly is deferred until the
  update is small or the iteration budget is exhausted.
  `FluxRheologyStep` returns the accepted solution, iteration count, update, and
  nonlinear residual. `FluxRheologyError` refuses invalid controls/state,
  unconverged linear or nonlinear solves, and incompatible prescribed boundary
  flux; a refused step never modifies the caller's previous state. The flux
  check uses the same edge quadrature and outward normals as the constrained
  BDM moments, with compensated summation and a roundoff-scaled tolerance.
  It prevents pressure pinning from hiding the omitted continuity equation.
  Picard convergence is conditional on the constitutive response, flow state,
  and timestep; exhausted budgets return an error rather than an accepted
  lagged-viscosity result. This is a fixed-mesh, constant-density BDF1 step,
  without a coupled thermal solve or nonlinear constitutive adjoint.
- `MaterialFluxSystem` (`material`) binds that same FEEC operator to the shared
  `fs-material::fluid_transport::FluidTransport` owner. Gas, mixture and liquid
  evaluations supply one frozen `FluidTransportState`; `nu=mu/rho`, force
  density [N/m³] divided by rho, physical mean-zero pressure [Pa] and kinetic
  energy per unit out-of-plane depth [J/m] are derived together. Numerical
  controls cannot independently override viscosity. Stokes assembly/solve,
  Picard iteration and IMEX BDF1 all use the existing operator; original source
  receipts remain available through `material()`. Invalid inputs, nonfinite
  outputs and failed convergence return `MaterialFluxError` without publishing
  a physical solution. The primitive singular-mesh panic contract still applies.
- `MaterialFluxSystem::bdf1_rheology_step` resolves a retained `FluidTransport`
  owner per triangle at the current objective shear rate. The caller receives
  a `MaterialFluxRheologyPoint` containing the cell, centroid in metres, and
  shear rate in inverse seconds, plus the reference transport. Its resolver
  can capture prescribed endpoint temperature, composition and material fields;
  those fields are fixed during the implicit solve. Viscosity is piecewise
  constant per cell and the same evaluated `mu/rho` supplies the symmetric
  volume stress and one-sided interface traction. One response is cached per
  cell and exact rate, so repeated quadrature queries reuse the material owner.
  The final residual assembly refreshes that cache at the accepted velocity.
  `MaterialFluxRheologyStep` returns the physical solution, convergence
  diagnostics and `MaterialFluxRheologySample` owners from that exact endpoint,
  without a misleading homogeneous-viscosity accessor. All samples must retain
  the reference density used for force, inertia, pressure and kinetic energy.
  A changed density or original material-source refusal stops before another
  linear solve; failure leaves the caller's previous state unchanged and
  publishes no fallback solution. This method supplies no temperature evolution,
  interface transport or viscoelastic memory.

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
   the identity `body_power = relative_power + F dot ambient_velocity`; a
   nonzero `F dot ambient_velocity` is flagged and requires a total
   gas-plus-boundary accounting window rather than asserting a sign for body
   power. `WorkWindow` reconstructs every public power from private kinematics
   and refuses forged public power or receipt-identity fields. Zero
   density/motion are exact zero limits; forbidden thin-gap-pressure and
   target-fitted families fail admission (aero-001
   through aero-009, G0/G3).

## Error model

Structured panics on programmer contracts (non-triangle elements,
degenerate cells, singular saddle factorizations) with teaching
messages. Solver quality is reported, not asserted: `FluxSolution`
carries the relative residual and iteration count. The ALE geometry audit
instead returns `AleGclError2` for malformed topology, non-finite input,
collapse, arithmetic overflow, and receipt allocation failure. Its finite
roundoff defect is evidence, not an admission failure.

The fully implicit rheological step instead refuses failed linear/nonlinear
convergence through `FluxRheologyError` and validates global boundary mass
compatibility before replacing a continuity equation with the pressure gauge.
Positive-finite constitutive samples and nonsingular mesh/factorization remain
primitive programmer contracts; callbacks must return finite physical data.

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

- `continuum` (default): the CFD, ALE, gas-film, and reduced-aerodynamics
  surfaces with their numerical dependencies.
- `quarter-wave` (default): the dependency-free distributed-line coordinate
  transform. Browser boundaries that only need this law may disable default
  features and select `quarter-wave`, avoiding unrelated continuum solvers and
  runtime infrastructure without copying the equation.

## Conformance tests

The `lc` unit tests reproduce Tesla's printed 925 Hz, 185,000 mi/s,
200-mile-wavelength, 50-mile-quarter-wave example after exact mile-to-metre
conversion. A second test establishes the inverse frequency/quarter-length
relation without adding a device-specific coil model.

`tests/battery.rs`: flux-001 basis exactness + conformity; flux-002
Stokes MMS orders + exact divergence; flux-003 pressure-robustness
identity; flux-004 cavity Picard + Ghia coarse band; flux-005
Taylor–Green BDF1; flux-006 adjoint gradient + annihilation.
The symmetric-stress lane adds flux-008 variable-viscosity rigid rotation
(zero strain/stress), flux-009 exact nonuniform-viscosity shear with its
reciprocal force and BDF1 fixed point, flux-010 layered Couette flow with
continuous traction through viscosity contrasts up to 10,000:1, flux-011
the objective constitutive-rate update, and flux-012 a nonuniform Carreau
manufactured flow with near-second-order refinement and an independently
reassembled nonlinear residual. The fully implicit lane adds flux-013 nonlinear
Carreau startup against its independently differentiated stress and endpoint
equations, flux-014 temporal refinement of a spatially exact shear field,
flux-015 unforced viscous energy decay, and flux-016 refusal of exhausted
nonlinear iterations, invalid timesteps, and incompatible boundary mass flux.
The sparse lane adds flux-017 agreement with dense factorization under
discontinuous viscosity, mass and convection; flux-018 a spatially varying
viscosity flow above the 1500-dof dense cutoff; flux-019 strict iteration budgets
and original-residual reporting; and flux-020 a 10,000:1 layered Couette flow
above the cutoff with independent all-cell divergence and boundary-moment gates.
These remain 2-D finite-element checks,
not a high-Re, large-scale, or cross-ISA qualification.

`tests/adjoint.rs` checks the actual transpose against independently scattered
rows of a nonsymmetric constrained flow matrix; viscosity sensitivities above
1500 unknowns against central finite differences with heterogeneous viscosity,
transient mass and frozen convection; an exact viscosity-dependent boundary
load derivative; and finite input, singularity and strict iteration-budget
refusals.

`tests/material_flow.rs` (G1/G3): hydrostatic pressure from physical force
density at different material densities, inverse-dynamic-viscosity Stokes flow
under material swaps, density-weighted kinetic energy, shared-gas pressure
changes in a manufactured accelerating BDF1 flow, the Picard creeping-flow
limit, and invalid/stale material state, boundary/load and timestep refusal.
The generalized-Newtonian adapter adds a two-temperature, two-material shear
flow with an independently solved continuous physical traction, pressure and
energy checks at two densities, endpoint rate/source ownership, cached query
counts, and a frozen-reference-viscosity control. Immediate and later material
source errors, density changes and nonlinear-budget exhaustion all refuse
publication without changing the previous velocity.

`tests/ale_gcl.rs`: flux-007 fixed mesh, rigid translation, affine
expansion/shear, two-cell interior cancellation, replay, translation
covariance, continuous-trajectory minima, and fail-closed endpoint,
topology, stale-cache, orientation, and finiteness cases.

`tests/reduced_aero.rs`: aero-001 zero-density/static limits; aero-002
reversal and stationary-ambient dissipation; aero-003 quadratic scaling and
source receipt mutation; aero-004 moving-ambient accounting identity;
aero-005 frame objectivity; aero-006 correlation-domain and missing-property
refusal; aero-007 alternative-correlation disagreement; aero-008 exactly-once
work; aero-009 hostile thin-gap/target-fit refusal; aero-014 forged wrench,
power, and receipt-identity refusal.

## No-claim boundaries

The base `MaterialFluxSystem` Stokes/Picard/IMEX methods bind one homogeneous
Newtonian state; its generalized-Newtonian BDF1 method admits cellwise
material/temperature/shear-dependent viscosity at constant reference density.
Neither path implements spatially varying density, compressible gas dynamics,
thermal evolution, phase change or material-interface advection. Gas
thermodynamic pressure selects material properties; the flow pressure is a
mean-zero gauge and must not be used as an EOS absolute pressure. Both BDF1
time levels use the same density; a density change during evolution needs an
explicit conservative state transfer. The binding neither adds material source
authority to raw transport tuples nor turns numerical convergence into physical
validation. Difficult constitutive laws may exhaust the Picard budget and
return a typed refusal.

The `lc` quarter-wave result is not an absolute-voltage, current, energy,
coupling, quality-factor, loss, load, air-breakdown, corona, or streamer-length
prediction. Those require independently sourced boundary data and owning
physics. `remote_terminal_profile_fraction` is dimensionless standing-wave
geometry and must never be labelled as volts or discharge reach.

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
- Mesh-independent iteration counts, multilevel saddle preconditioning and
  production-scale throughput are not established by the sparse finite-element
  fixtures. `linearize_rheology` supplies constitutive/convection state
  Jacobians for prescribed-velocity boundaries. Mixed-boundary nonlinear
  Jacobians and time-dependent adjoint checkpointing remain outside that API.
- Nonlinear convection/backflow and transient open-boundary conditions.
  Steady generalized-Newtonian Stokes velocity/total-traction boundaries
  are provided by the APIs below.
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
  window but is not an `fs-couple` closed first-law/entropy audit; nonzero
  far-field boundary power specifically needs that separate total accounting
  closure.

### Bounded gas-film foundation (bead frankensim-b8bxd.10)

`gas_film` adds `solve_isothermal_gas_film_1d`, a deterministic, safe-Rust
structured one-dimensional finite-volume foundation for the **isothermal ideal
gas** Reynolds mass balance, at unit out-of-plane width.  Its fields use SI:
absolute pressure [Pa], gap [m], temperature [K], viscosity [Pa s], gas mass
per width [kg m^-1], mass flux per width [kg m^-1 s^-1], and wall power/heat
per width [W m^-1].  It retains canonical case/model/species/EOS/viscosity/
thermal/frame identities, a caller or synthetic input authority, a deterministic
seed, and a bounded iteration budget.  The input checks a caller-declared
density against the same ideal-gas EOS used in storage; density is not silently
treated as an independent incompressible property.
Caller-declared viscosity/gap/pressure uncertainty bounds are retained in the
receipt but are deliberately not propagated, calibrated, or certified.  A
restart binds a deterministic invariant-configuration fingerprint over model/
source/gas cards, grid length/count and fixed mask/topology kind, slip/
roughness/applicability, uncertainty, gauge reference, timestep, budget, and
authority. Dynamic current gaps and wall motion remain per-step inputs: restart
checks each new active gap against accepted gap plus declared `dh/dt * dt`
within a named numerical tolerance, while allowing changed tangential speed.
`dh/dt` may be a uniform fallback or an exact one-value-per-geometric-cell
profile; the latter is used consistently for fresh-state reconstruction,
restart evolution, storage, and normal pressure-work accounting so changing
wedge slopes are not collapsed to a mean rate.
The checkpoint also retains the full accepted-request fingerprint for
provenance, but never mistakes it for a static restart constraint.

The discretization retains pressure-driven Poiseuille mass flux plus moving-wall
Couette mass flux, implicit storage, uniform squeeze through exact old/new gap
values, upper-wall Couette shear/relative dissipation heat, and signed normal
moving-gap power `-integral(p dh/dt dx)`. It also carries
signed boundary/vent enthalpy transport from a caller-declared specific
enthalpy; this transport receipt is not a thermochemical calculation. Gauge
load is relative to a separately declared absolute pressure reference, never
implicitly the initial condition. Per-cell receipts name both gas-on-upper-wall
and upper-wall-on-gas normal/tangential tractions with opposite signs; the old
`upper_wall_shear_pa` remains only as the upper-wall-on-gas alias. `Sealed`,
two-reservoir `Open`, and one named-cell `Vented` topology are admitted.  The
receipt retains gas storage, every boundary/vent mass flux, their global closure
defect, pressure load, wall shear, wall work/heat, residual, iteration count,
and continuum-envelope maxima. A checkpoint binds active pressure and previous
gap values to its model/case/authority identity only after a converged step; no
partial state is published.

`ContactExclusionMask` is an explicit contact handoff, not a gap floor: this
first slice admits only a contiguous active prefix followed by excluded cells.
An interior hole, topology change, too-few active cells, non-positive gap or
absolute pressure, EOS-density mismatch, excessive slope/Mach/Knudsen number,
or requested slip/rarefied/roughness model is refused with `GasFilmError`.
Smooth-wall admission requires a named roughness bound strictly below every
active gap. There is no
anonymous `h_min`, pressure clamp, liquid cavitation substitution, or
Elrod--Adams semantics.

The implementation is bit-deterministic on one ISA by fixed cell/face order and
bounded Jacobi/Picard sweeps. It has no RNG, no `Cx` cancellation hook, no
entropy model, and no thermochemical or material-data admission; the
retained source identities and `GasFilmInputAuthority` state those limits rather
than minting material/EOS authority.

`tests/gas_film.rs` supplies G0/G1/G3 coverage for uniform sealed equilibrium,
analytic planar Couette shear/heat, uniform squeeze mass preservation,
open Poiseuille direction/closure, vent closure, spatial/time refinement of
analytic limits, reversal and pressure-density scaling, static contact
exclusion, rarefaction/negative/slip/roughness refusal,
iteration-budget refusal, and deterministic checkpoint replay.  These are
software/numerical checks of declared equations only.

No gas-film result establishes thermochemical validity, lubrication-regime
selection, roughness or slip closure, continuum validity outside the stated
gates, liquid/mixed/EHL/cavitation behavior, moving-domain topology changes,
2-D/axisymmetric geometry, energy/entropy closure, coupled contact, or any
device outcome.  The blocked fs-matdb, fs-thermochem, fs-tribo lubrication, and
fs-couple port-schema prerequisites remain required for those claims.

## Mixed steady-Stokes velocity and total-traction boundaries

With the continuum feature, FluxSystem::assemble_stokes_mixed accepts exactly
one StokesBoundary owner for every mesh boundary edge. Its constant-viscosity
entry delegates to assemble_stokes_mixed_variable_viscosity, which admits a
pure (cell, point) -> nu field at constant reference density. Both use the
existing BDM1/P0 symmetric-stress discretization:
-div(2 nu D(u)) + grad(p) = f and div(u) = 0. Velocity edges constrain normal
moments strongly and receive tangential SIP terms. Natural edges receive
t = (2 nu D(u) - p I)n and retain unknown normal velocity moments. They have
no velocity penalty or extra pressure-boundary matrix term. The variable field
uses the same one-sided face traces and harmonic viscosity/normal-height SIP
weights as the closed-boundary operator. The caller supplies total traction
with its physical viscosity already included; the assembler does not multiply
that boundary stress by viscosity again. The field entry ignores params.nu,
while the legacy constant entry preserves its positive finite nu admission.

One connected fluid domain and at least one velocity edge are required.
Duplicate, missing, interior and out-of-range ownership is refused before
callbacks run. Boundary quadrature samples are retained once for consistent
strong moments and SIP. Non-finite loads and assembly, incompatible net flux
for all-prescribed velocity, singular factorization and failed explicit
original or row-equilibrated residual admission return StokesBoundaryError
without a solution. An invalid viscosity sample retains the offending cell,
point and value in InvalidViscosity, stops subsequent viscosity queries and
discards the private assembly before any solve or publication. No hidden
grounding is added for pure-traction rigid modes.

All-velocity problems retain the original pressure pin and report MeanZero
pressure. With any total-traction edge, all continuity equations are retained,
no pressure pin is added and StokesSolution::pressure reports TractionDatum
without subtracting its mean. The pressure and traction are kinematic
(stress divided by one constant reference density); TractionDatum does not
mean thermodynamic/EOS absolute pressure. A pressure shift dp requires changing
the prescribed total traction by -dp*n. A pressure-only load is -p*n, but a
fully developed channel's end traction also contains tangential viscous stress.

The legacy assembly entry points keep their existing forms. The new checked
solve shares the same LU/sparse-FGMRES kernels and additionally recomputes both
original and row-equilibrated residuals before publishing a result. Its
read-only matrix is available to the existing discrete transpose solver.

tests/stokes_boundary.rs supplies the pressure-driven channel analytic
refinement/flux oracle, exact affine shear and pressure-datum shift patch,
all-velocity operator/solution bit comparison, and boundary-admission cases.
The spatial-viscosity cases independently check a two-layer pressure-driven
channel's velocity, continuous shear stress, pressure datum, flux and mesh
refinement; an exact heterogeneous affine patch; constant-law operator bit
identity; and checked coefficient refusal without publishing a fallback.
examples/pressure_channel.rs is runnable with
cargo run -p fs-flux --example pressure_channel.

The frozen-field assembly above covers steady linear, constant-density Stokes with frozen
positive spatially varying Newtonian viscosity on the existing straight-sided
2D triangulation. It does not add
nonlinear convection/backflow, open-boundary time integration, pure-traction
rigid-mode quotients, disconnected domains, material-bound SI pressure
conversion or continuum error bounds. Existing synchronous cancellation and
dense/iterative solver limitations remain in force.

## Horizontal capillary meniscus transport

The existing fully filled Hagen–Poiseuille API remains a distinct model.
`step_newtonian_circular_meniscus` advances a liquid column in a horizontal
rigid circular tube with constant Newtonian viscosity, no slip, negligible
inertia and a vented gas of prescribed pressure. Its SI input explicitly names
radius, tube length, current wet length, surface tension, cosine of the
through-liquid equilibrium angle, applied liquid-minus-gas pressure and time.
The caller binds that angle to the actual solid/liquid/gas system; material
names never select a coefficient.

The meniscus pressure is \(p_g-p_l=2\sigma\cos\theta/R\).
The exact constant-input advance satisfies
\(d(x^2)/dt=R^2(\Delta p+2\sigma\cos\theta/R)/(4\mu)\).
This handles wetting, nonwetting resistance, pressure-driven reversal and
explicit zero surface tension without dividing by the initial wet length.
An inlet/outlet event returns its actual elapsed time and endpoint; the
caller must handle any remaining requested time with the next boundary model.
There is no silent throughflow after the meniscus leaves this domain.

The returned volume is signed into the tube. Imposed-pressure work, the
Young wall-surface energy change and viscous loss close
\(W_p-\Delta E_{wall}=D\geq0\). A translating meniscus has constant area in this
model. Initial meniscus formation and final end-topology energies are excluded.
All inputs and final transfers are checked before returning a state;
nonfinite or out-of-domain values produce the named capillary error. The
fixed-size scalar computation allocates nothing and adds no dependencies.
It has no internal cancellation point or cross-ISA bit-equality claim.

Seven focused tests in the existing capillary unit module exercise analytic
penetration and wall energy, parameter scaling and timestep composition,
nonwetting drainage/reversal, event timing, zero-motion work, the integrated
Poiseuille limit and input refusal. The three existing fully filled tests
remain. These are G0/G1/G3 checks of this reduced law, not experimental proof
of dynamic contact-line or general free-surface behavior. Inertial startup,
gravity, gas drag, entrance losses, hysteresis, evaporation and evolving
material coefficients require other formulations.

Reference: E. W. Washburn, Physical Review 17, 273 (1921),
https://doi.org/10.1103/PhysRev.17.273.

## Consistent generalized-Newtonian residual linearization

`FluxSystem::linearize_rheology` supplies the sparse state Jacobian and actual
residual of steady or BDF1 generalized-Newtonian flow on the existing 2D BDM1/P0
mesh with prescribed velocity on every boundary edge. The pure constitutive
callback returns kinematic viscosity and its analytic derivative with respect
to `gamma = sqrt(2 D:D)` at the same quadrature point. Reference density,
previous velocity, force, boundary data, temperature and composition stay fixed.
The callback may be queried more than once per point; it must not mutate its law.

The implementation adds local sparse derivatives to the actual frozen primal
operator: cell stress, both one-sided interface viscosities, harmonic traction
weights and penalties, and viscosity-dependent weak boundary loads. Optional
convection differentiates both transported and advecting velocity. Strong normal
constraints, pressure pin, incompressibility and BDF1 history retain the primal
rows. No dense matrix or one-global-assembly-per-column differencing is used.
For interface delta_s=nu_s/h_s, q=delta_0 delta_1/(delta_0+delta_1), the traction
coefficient is h_s q and dq/dnu_s=weight_s^2/h_s; this also differentiates the
penalty 2 sigma q without multiplying disparate viscosities together.

The returned Jacobian can be passed to `solve_adjoint_checked`:
`dJ/dtheta = J_theta - lambda^T R_theta`. Parameter-dependent external fields
and loads still belong in R_theta. The reported `nonsmooth_upwind_points` counts
zero normal-flux sites with unequal traces. At these sites the >=0 active branch
is a semismooth linearization, not a unique classical derivative; smooth-gradient
claims must exclude switches or establish that admissible directions cannot
cross them. Invalid state/tangent/load data and nonfinite assembly refuse before
publication. A smooth radial viscosity must supply zero derivative at zero strain.

Four G1/G3 regressions in `tests/rheology_linearization.rs` compare Jacobian
products against independently perturbed primal residuals for shear-thinning
and shear-thickening fields with a 100:1 coefficient contrast, with and without
convection. They check a nonlinear BDF1 viscosity adjoint against independent
nonlinear re-solves, expose the frozen-matrix gradient defect, preserve the
Newtonian operator bits, and exercise zero-strain, upwind-switch and refusal
semantics. This supplies a fixed-mesh state derivative; mesh/shape derivatives,
material-temperature evolution, mixed open-boundary Jacobians and time-history
checkpointing require their own operators and evidence.

## Gravity-aware capillary rise and drainage

`step_newtonian_inclined_circular_meniscus` extends the horizontal
Newtonian/no-slip model with constant liquid density and signed gravity
projected along the tube. Its input embeds `CircularMeniscusInput` and names
`density_kg_m3` and `opposing_gravity_m_s2`; positive projection means uphill.
The inlet is at the reservoir elevation and the vented gas pressure is
prescribed uniformly. Zero projected gravity delegates exactly to the
horizontal path.

The motion satisfies `x dx/dt = K (A - b x)`, with
`K = R²/(8 mu)`, `A = delta_p + 2 sigma cos(theta)/R`, and
`b = rho g_parallel`. The separated analytic travel-time map is inverted
within a monotone bracket. Stable log1p endpoint weights retain the small
gravity limit; a bounded Newton/bisection solve approaches the finite Jurin
height without crossing it. The gravity-only branch is linear. Inlet/outlet
events retain actual travel time; an equilibrium is approached asymptotically
and does not authorize throughflow.

The result embeds the existing meniscus result and adds gravitational energy
change and accepted hydrostatic pressure. The viscous account now closes
`W_pressure - delta E_wall - delta E_gravity = D >= 0`.
The affine driving force's displacement average evaluates dissipation without
subtracting nearly equal large work terms. The computation allocates nothing
and adds no dependency. It retains the horizontal model's no-claim boundaries
for inertia/startup, entrance loss, gas drag, dynamic wetting and changing end
topology. It has no internal cancellation point or cross-ISA bit-equality claim.

The focused G0/G1/G3 tests compare with the dimensionless Washburn/Jurin
analytic solution, rise and drainage around equilibrium, inclination scaling,
nonwetting downhill entry thresholds, the gravity-only and horizontal limits,
time composition, end events, signed energy closure and input refusal.
Reference: Fries & Dreyer (2008), doi:10.1016/j.jcis.2008.01.009.

## Conservative triangular capillary surface forces

`capillary::surface` supplies allocation-free triangle force kernels and a
linear-time patch assembler over existing coordinate/index slices. Constant
or prescribed per-face tension multiplies geometric area. The instantaneous
force is the negative area-energy gradient; open patches retain their boundary
line forces. Shared vertices assemble the actual surface-stress resultant.

`triangle_surface_tension_step` uses endpoint Gram invariants and midpoint
edges to construct an area discrete gradient. Its forces have zero resultant
and zero torque about midpoint coordinates, and their displacement work equals
minus the geometric surface-energy change. Compensated product differences
avoid false area for collinear edges. A symmetric invariant correction restores
small work differences lost to floating cancellation on skinny faces without
adding a rotational force. The response retains both endpoint areas/energies,
the energy change, actual force work and prescribed-coefficient free-energy
supply. For changing tension, `W + delta E = A_mean * delta sigma` per face.
That supply is explicitly external free energy, not inferred thermal heat.

Patch assembly publishes only complete results and uses compensated scalar
reductions in fixed face order. Invalid coordinates, zero-area faces,
negative/nonfinite tension, inconsistent field sizes, invalid connectivity
indices and unrepresentable outputs refuse. No intermediate deformation path,
manifoldness, topology change, surfactant transport, thermal surface entropy,
fluid advection or time integration is claimed by the force operator. Triangle operations are bounded;
patch callers own their tiling/cancellation boundary.

Nine G0/G1/G3 tests check analytic triangle line forces, independent energy
gradients, finite-deformation work and midpoint torque, rigid-motion objectivity,
unit scaling, coefficient supply, variable-tension patch scatter, skinny-face
work, invalid geometry/data and analytic endpoint tangents under nonuniform
deformation and finite rotation. The discrete-gradient principle follows
Gonzalez (1996), doi:10.1007/BF02440162; the Gram formulas are derived locally.

## Conservative inertial capillary material surfaces

With `continuum`, `capillary::surface::dynamics::CapillarySurfaceDynamics`
binds fixed triangle connectivity and strictly positive reference nodal masses
to the existing `fs-time::galpha::SecondOrderState`. It advances actual xyz
positions and velocities through the original matrix-free structural
generalized-alpha owner at `rho_inf=1`. No other time integrator or liquid
solver is introduced. Initial acceleration is physical force divided by mass;
subsequent acceleration retains the native Newmark auxiliary-state semantics.

The midpoint resistance is the negative original endpoint capillary force.
Its analytic tangent differentiates the Gram formula and its symmetric
skinny-face work correction, then composes the endpoint/midpoint chain rule.
The original Newton--Krylov driver owns globalization and linear iteration;
positive lumped mass supplies its right preconditioner. Face callbacks allocate
nothing and never form a dense mesh Jacobian. Explicit topology/restart limits
and native iteration budgets bound work. Cancellation is checked between native
callbacks and at most every 64 faces within surface evaluation.

Every accepted endpoint independently satisfies the caller's SI tolerances for
`M delta v - dt (F_surface + F_external)`,
`delta x - dt v_mean`, and
`delta kinetic + delta surface - external force work - A_mean delta sigma`.
The nonlinear solver's convergence report does not substitute for these gates.
Kinetic change uses mean velocity times momentum change to avoid subtraction
of large aggregate energies. Rejection leaves positions, velocities,
accelerations, physical time, accepted-step count and history unchanged.
The trial does not clone the accumulated trajectory history.

Prescribed per-face tension change supplies external free energy, never
inferred heat. This rung represents a tension-only material surface with fixed
reference mass, including explicitly authored ideal inertial films. It does
not replace an elastic-solid material law, solve bulk incompressibility or
assign physical mass from current area. The caller identifies every interface
(including both sides of a film) and owns fixed-topology applicability.
No contact/support constraints, bending, viscosity, intermediate-path
nonintersection, remeshing, surface entropy, thermal coupling or dynamics
adjoint is claimed. Temporal and spatial refinement remain necessary.

Five focused G1/G3/G4 tests exercise the analytic midpoint map for an equilateral
breathing patch and second-order temporal convergence; actual rotating,
nonuniformly deforming surface evolution with energy and linear/angular momentum
conservation; explicit coefficient supply, held-force equilibrium and the
analytically accelerated zero-tension limit; finite
rotation/translation objectivity and checkpoint replay; and transactional
invalid-data, cancellation, nonlinear-budget and independent-balance refusals.


## Damped Newton BDF1 for stiff rheology

`FluxSystem::bdf1_rheology_newton_step` advances the same fully implicit
constant-density generalized-Newtonian equations using `linearize_rheology` and
the existing checked dense/sparse linear solver. The caller supplies the exact
viscosity tangent, immutable previous velocity and endpoint loads/boundaries.
It admits all-velocity boundary ownership and net flux through the existing
mixed-boundary owner before the pressure pin can omit a continuity equation.

`FluxNewtonControls` bounds Newton correction solves and halvings per line
search (0 through 52), and supplies the nonlinear residual tolerance. The search
holds Jacobian row scales fixed and requires Armijo decrease (coefficient 1e-4)
of the actual weighted nonlinear residual, or an already converged endpoint.
Trial assemblies build no Jacobian or factorization. Invalid trial constitutive
states are rejected; neither the material law nor physical timestep is clipped,
regularized or silently substepped. Upwind switches use the existing semismooth
active derivative. This is conditional convergence, not a global Newton theorem.

Both original and row-equilibrated residuals must converge at the returned
endpoint. A small update alone is insufficient. Accepted fields retain these
nonlinear residuals, actual Newton counts and rejected-trial counts. Invalid
inputs, incompatible boundary flux, failed linear work, exhausted searches or
unmet nonlinear budgets return typed errors without changing the input velocity
or publishing a partial flow. The callback must remain pure during trial work.

Three G1/G3 regressions in `tests/rheology_newton.rs` check stiff regularized
shear thickening that exhausts the existing Picard budget, agreement with an
independently iterated Carreau solution, separately assembled momentum residual
and divergence, active backtracking, zero-state acceptance, incompatible flux
and strict iteration/search budgets. The four constitutive-Jacobian regressions
remain green. These establish the numerical fixed-mesh BDF1 capability, not
experimental fluid calibration, mixed/open-boundary Newton flow, evolving
thermal fields, variable density or a time-history adjoint.

## Generalized-Newtonian mixed velocity/traction Stokes flow

`FluxSystem::solve_stokes_mixed_rheology` closes steady Stokes viscosity on the
current cell strain rate `gamma = sqrt(2 D:D)`, with the same prescribed-velocity
and total-traction boundary ownership as the frozen-field assembly. Boundary
samples are admitted and frozen once. The initial velocity supplies the first
rates; the solver does not probe the material at an invented zero-shear state.
The first solved flow establishes the boundary and continuity constraints.
Later corrections may be under-relaxed by the caller in (0, 1], while the
undamped update is measured so that tiny relaxation cannot manufacture
convergence. Acceptance requires a small update and independent reassembly of
the original and row-scaled nonlinear residuals at the accepted velocity.
The returned `FluxSolution` residuals describe those nonlinear equations and
are checked against the nonlinear tolerance, separately from the inner linear
solve tolerance. Exhausted budgets or failed admission publish no solution and
leave caller input unchanged.

`MaterialFluxSystem::solve_stokes_mixed_rheology` owns the SI conversion:
dynamic viscosity, force density in N/m^3 and total traction in Pa are divided
by the same immutable reference density. Returned pressure is in Pa with its
traction datum retained, and kinetic energy is in J/m. The exact original
material response and query at every accepted endpoint cell are retained.
A shared cell/rate cache also serves the existing material BDF1 entry; unchanged
cell/rate queries are reused instead of re-querying sources in quadrature
loops. Original material errors take precedence over generic solver errors,
and changed density refuses.

Four focused native tests on the pinned compiler cover an independently
integrated open shear-thinning channel, a source-owned two-layer thickening
flow, Newtonian reduction and refusal boundaries. Channel velocity L2 error
falls from 0.0119366 to 0.00328519 under refinement (order 1.8613), with 2.1039%
fine-grid outlet-flow error, divergence 6.30e-16 and nonlinear residual
1.55e-11. The layered flow preserves physical velocity, stress and pressure at
reference densities 2 and 7, scales kinetic energy coherently, and retains
the exact source objects with one query per cell/state. Its velocity L2 error
is about 3.50e-10, stress error 1.22e-10 and physical pressure error 9.68e-9 Pa.
Both existing material BDF1 ownership/refusal regressions pass after the cache
extraction.

This entry covers fixed-geometry, steady generalized-Newtonian Stokes with
fixed external temperature/material context and one reference density.
It adds no nonlinear convection/backflow policy, open-boundary time integration,
thermal evolution, viscoelastic memory or mixed-boundary constitutive tangent.
Convergence is conditional and bounded by the caller's explicit solve budget;
these finite fixtures do not establish continuum error bounds or production
throughput.

## Experimental multidimensional shock transport

Feature `experimental-compressible` exposes `compressible::Euler2d` on an
admitted `CartesianGrid`. This is the stationary, frozen-composition first
multi-D rung of extension-plan §5.8 Tier 4. `fs-flux` owns spatial evolution;
`fs-gas` owns the original NASA-9 thermodynamics, normal-face HLLC/HLL closure
and exact one-sided reflecting-wall shock/rarefaction solve. No constant-gamma
model replaces the retained variable-Cp source.

The state evolves density, both momentum components and full referenced total
energy, including both kinetic components. Each face is solved once and shared
by its neighbouring cells. Unsplit MUSCL reconstructs density, two velocities
and pressure with MC slopes and a common source-admissibility reduction; SSP
RK2 integrates both directions together. The accepted face CFL is
`dt*(max(a_x_left,a_x_right)/dx + max(a_y_lower,a_y_upper)/dy)` in every cell
and both stages, including any HLL wave-speed expansion. Its caller-selected
bound is at most 0.45. A CFL refusal halves the step without dropping order;
an actual state-admissibility refusal retries at first order with a smaller
step. No conservative state or source temperature is clipped.

MUSCL supplies reconstructed primitives directly to the checked `fs-gas`
primitive-face entries. That owner performs the same caloric and full-energy
precision admission, reusing its just-computed thermodynamic properties instead
of immediately inverting the new conservative energy. First-order faces retain
the original conservative decoding path. This removes redundant source work
without changing the Riemann, wall, fallback, CFL or timestep equations.

Axes may be independently periodic, transmissive or stationary reflecting.
Periodic pairs reuse identical flux bits. Transmissive means zero-gradient
continuation, not a nonreflecting subsonic outlet. Walls impose exact zero
mass, tangential-momentum and energy flux and use nonlinear source-owned wall
pressure. Grid depth converts areas, volumes, impulses and energies to coherent
SI extensives; constructor admission also refuses overflow of extensive totals.

Every accepted step retains compensated mass, vector momentum and energy
balances plus the time-integrated outward exchange on EACH boundary face.
The two RK stages supply these impulses. At a wall they are the force impulse
on its support, allowing `fs-couple` to preserve a spatially varying shock load
and its moment instead of sampling an endpoint or averaging the whole wall.
Failed or cancelled steps preserve all cells, time, last report and accepted
boundary impulses. Workspace is fallibly allocated at construction under the
shared 1,048,576-cell cap and reused; accepted steps allocate no numerical
workspace. Fixed-order cell/face loops poll cancellation at bounded boundaries.
Parallel execution, SIMD throughput and cross-ISA determinism are unmeasured.

Seven focused native regressions cover periodic uniform flow, source/geometry
and extensive admission, planar Sod against the existing independently checked
1-D owner, x/y rotation, nonlinear wall reflection and support impulse, a
circular pressure blast, smooth oblique contact refinement and atomic refusal/
cancellation/restart. Oblique-contact density L1 errors on 12²/24²/48² grids are
3.449212181169e-3 / 8.636199681318e-4 / 2.117512986278e-4, giving observed orders
1.998 and 2.028. The focused no-default-feature `fs-flux` library configuration
passes 47 tests, including the molecular-transport extension below, on the
pinned nightly. These are manufactured/conformance
results, not experimental validation or a production performance benchmark.

This lane remains Estimated engineering shock capturing. It does not supply
the separate fully discrete entropy/invariant-domain theorem, vacuum states,
reacting/composition-changing 2-D flow, viscous/thermal shock structure,
internal obstacles, curved/cut cells, moving meshes, reciprocal moving-wall
FSI, turbulence or a shock adjoint. Unsupported shocked temperatures refuse
against the original source interval instead of extrapolating it.

`cargo run --release -p fs-flux --no-default-features --features
experimental-compressible --example shock_box -- 32` runs an off-centre circular
nitrogen pressure pulse through more than one box acoustic crossing and wall
reflection. It uses the retained 200–1000 K NASA-9 coefficients also used by
`fs-gas`'s `shock_tube` example, with the same explicit example-level material
association. CSV contains position, density, both velocities, pressure and
temperature; stderr reports accepted steps, maximum CFL, mass/energy residuals
and integrated wall-support impulse. The 32² run advances 456 accepted steps
to 2.578101380129e-3 s with zero mass residual, 2.27e-13 J referenced-energy
residual and less than 1e-15 kg m/s momentum/support mismatch. This is an
executable engineering example, not a nitrogen shock experiment or benchmark.

A bounded before/after diagnostic of the primitive-face reuse used the same
48² example, 687 accepted steps, pinned nightly-2026-07-06, optimized release
profile, and shared Linux x86_64 container reporting AMD EPYC 9V74. One warm
binary run fell from 43.5289 s to 26.2137 s (1.6605× observed wall-time ratio).
Both runs retained zero mass and referenced-energy residuals. Across all 2304
exported cells, pressure CSV values agreed; maximum differences were 1.00e-13
kg/m³ in density, 1.00e-11 m/s in velocity and 1.00e-10 K in temperature.
This is a single-run diagnostic under shared machine load, not a portable
throughput guarantee or a repeated statistical benchmark. The direct numerical
regressions remain the correctness gate.

### Source-owned molecular transport in the same Cartesian gas

`compressible::molecular::MolecularTransport2d` binds the original
`fs-gas::multicomponent::dissipation::TransportLaw` to that same `FrozenGas`,
grid and boundary conditions. Admission compares the complete original NASA-9
models, not just species labels. `KineticTransport` therefore carries its
original collision records into the spatial solver; changing those data changes
actual momentum and temperature. `ConstantTransport` remains an explicitly
declared approximation. Every active face evaluates viscosity, bulk viscosity
and conductivity at its local midpoint temperature, pressure and original
frozen mole composition. Negative, nonfinite or source-refused values stop the
transaction without substituting air constants.

The molecular flux is `(0, -tau_xn, -tau_yn, -u_face dot tau_n - k grad_n T)`.
The complete planar velocity gradient includes both cross derivatives and the
three-dimensional Newtonian trace `tau = mu*(grad u + grad u^T - 2/3 div u I)
+ bulk*div u I`, with zero out-of-plane velocity. Shared face fluxes conserve
mass, both momenta and total energy. Dissipated kinetic energy remains in that
same total-energy state and raises the original caloric temperature; there is
no second viscous-heating source. Fourier diffusivity uses the source's
volumetric **Cv** because the transported variable is internal energy at fixed
cell mass and volume. There is no independently authored heat-capacity field.

Periodic pairs reuse the same flux. Stationary reflecting walls enforce zero
heat and work and zero tangential viscous traction, while retaining their normal
viscous support impulse. Transmissive boundaries have zero molecular flux.
These are adiabatic slip and zero-gradient conditions, not no-slip or finite
thermal reservoirs. Molecular impulses and balances include physical grid depth.

`advance` integrates a complete requested molecular interval with SSP RK2 and
bounded retries, then atomically publishes the field and per-face impulses. It
does not advance the Euler clock or replace its prior advection metadata.
The global parabolic estimate includes normal and cross stress coefficients,
both cell widths, the smallest current density/Cv and doubled wall gradients;
both RK stages check it and the original EOS admits every trial state. This is
a conservative engineering step-size estimate, not a proved nonlinear
invariant-domain or entropy bound. Source or work-limit refusal and cancellation
leave the caller's accepted field and operator impulses unchanged.

`MolecularEuler2d` composes molecular half step, original Euler step and molecular
half step as one atomic interval. If post-heating advection shortens its CFL
interval, the entire split restarts from the original field. Its public report
and per-face impulses combine all three phases; the borrowed Euler owner's
metadata still describes only the middle advection phase. The symmetric
composition permits second-order smooth evolution when the constituent paths
retain their order, while limiting and first-order retries reduce that order.
Workspace is bounded and allocated at preparation, then reused. Runtime material
law implementations may mutate only their evaluation workspaces.

Five focused G0/G1/G2/G3/G4 regressions exercise a 2-D Fourier mode, oblique
longitudinal and transverse Newtonian modes, collision-diameter-only pressure
and shear changes, wall impulse/energy closure, late cancellation/replay,
zero-transport recovery and atomic work refusal. Synthetic source and collision
data test composition and conservation; they are not measurements of a real gas.
This extension keeps spatially uniform frozen composition. It does not supply
species diffusion, chemistry, Soret/Dufour effects, no-slip boundaries, turbulence,
rarefied gas, physical shock-thickness validation or mutual fluid/solid thermal
evolution. Those capabilities require their original owners and additional
admitted source data, rather than inferred closures.

## Three-dimensional conservative shock propagation

`compressible::spatial::Euler3d` supplies a volumetric Euler field behind the
same opt-in `experimental-compressible` feature. The original `fs-gas` source
and its shared planar/spatial HLLC/HLL kernel own every caloric and Riemann
decision. All three momentum densities and their kinetic energy evolve in the
same five-component flux; the two transverse components are never separate
passive scalars appended after the energy calculation.

Uniform Cartesian cells retain explicit dx/dy/dz, actual face areas and cell
volumes. Cell order is x fastest, then y, then z. Boundary faces use the other
axes in increasing order: y/z on x sides, x/z on y sides, x/y on z sides. Six
physical sides support stationary reflecting or zero-gradient transmissive
ends; an axis may instead be periodic. Transmissive does not mean a transparent
subsonic acoustic boundary. Periodic endpoints reuse identical face flux bits.
Exact thermally perfect wall reflection has zero mass, tangential momentum and
energy flux and retains the nonlinear shock/rarefaction pressure response.

The spatial step uses unsplit MC-limited primitive reconstruction and SSP RK2.
Its CFL sums x/y/z contributions, using each cell's maximum adjacent numerical
signal speed in every direction, and checks both stages including HLL fallback
speeds. A source-domain limiter contracts all six reconstructed traces from
the admitted cell average. Source refusals retry at first order with halved
time steps; CFL-only retries retain the requested reconstruction. No cell
average or temperature is clipped. All cell, face, stage and boundary buffers
are fallibly allocated at construction under the existing cell-count budget.
Accepted steps reuse them. Cancellation is polled through cells, faces and
reductions; failed trials retain field, clock, last report and accepted loads.

Every physical boundary face publishes its outward RK-integrated mass,
three-axis momentum and referenced-energy transfer. Extensive totals use the
same compensated reduction as the planar owner. Closed-wall fluid/support
momentum therefore closes through the actual update, including nonuniform
pressure footprints; endpoint pressure times duration is not substituted for
the numerical impulse. Structural consumers use these original per-face loads.

Seven focused tests cover noncubic physical geometry and uniform flow; smooth
three-axis contact refinement; Sod evolution matching the original 1-D source
on each axis; extrusion matching the original variable-caloric 2-D solver;
spherical blast symmetry/conservation; all-wall impulse and late-cancel replay;
and geometry/source/extensive admission. Contact density L1 errors at
12³/24³/36³ are 1.08724113e-3, 2.86537347e-4 and 1.03126986e-4, with observed
orders 1.924 and 2.520 over this finite resolution range. These are numerical
fixtures, not measured gas validation or an entropy/invariant-domain theorem.

`cargo run --release -p fs-flux --no-default-features --features
experimental-compressible --example shock_box_3d -- 12` runs the original
200–1000 K nitrogen NASA9 source in a reflecting 1 m³ box with an off-centre
800 K sphere in a 400 K gas. It writes a full x/y/z/rho/u/v/w/p/T CSV and reports
balances. The retained 12³ run advanced 164 steps through 0.7 ambient acoustic
crossings, maximum combined CFL 0.316812; mass residual was -1.110223e-16 kg,
referenced energy residual 0 J and fluid/support momentum residual components
at most 1.95e-16 kg m/s. The observed 15.70 s wall time is one container run,
not a general performance claim.

This fixed-box Euler rung does not itself provide evolving composition,
reactions, internal obstacles, cut cells,
moving geometry, reciprocal moving-wall FSI, a shock adjoint, physical shock
thickness or experimental validation. The full EOS-specific fully discrete
entropy/invariant-domain and moving-domain programs remain separate work.

### Three-dimensional Newtonian stress, Fourier heat and the coupled field

`spatial::molecular::MolecularTransport3d` binds the original gas's exact species
models to the existing `TransportLaw`. It evaluates shear viscosity, bulk
viscosity and conductivity at each evolving face temperature and pressure.
The full tensor is `tau_ij = mu*(du_i/dx_j+du_j/dx_i) +
(bulk-2*mu/3)*div(u)*delta_ij`. Both tangential stresses, all cross derivatives
and all three diagonal contributions participate. Its momentum flux is
`-tau*n`; total energy carries `-u dot tau*n - k*grad(T) dot n`. Internal
heating follows from conservative total energy, with no duplicated viscous
source. One stored face flux transfers equal and opposite amounts to neighbors.

The operator reuses the volumetric Euler update and boundary capture, original
caloric admission, shared molecular controls and named source/refusal errors.
It preallocates cell gradients and face buffers; each cell gradient is computed
once per stage and reused on adjacent faces. The first prepared face pass
provides both the parabolic estimate and first SSP stage, avoiding duplicate
source/gradient work on successful molecular substeps. The bound includes all
three cell widths and full normal/cross-stress coefficients, while both stages
must pass source and parabolic admission. This is an engineering time-step
estimate, not a nonlinear entropy or invariant-domain theorem.

`MolecularEuler3d` composes molecular half step, original Euler step and molecular
half step over one accepted physical interval. If heating changes the Euler CFL
and shortens the interval, the entire trial restarts from its original field.
The complete report and per-face loads combine all phases. The borrowed Euler
owner's report describes only the middle advective stage; consumers needing
physical total wall traction must use the molecular owner's report and loads.
All accepted field, clock, history and impulses remain unchanged on failed work,
source refusal or cancellation, including the final publication checkpoint.

Stationary reflecting faces are adiabatic slip walls: normal viscous traction
acts on the support, and tangential stress, heat and wall work vanish. Open ends
have zero molecular flux. These are neither no-slip nor fixed-temperature
walls. The combined normal wall traction is signed `p-tau_nn`; it may be tensile
and cannot be treated as a positive absolute pressure or clipped to zero.

Five focused tests check all-three-axis heat diffusion, oblique longitudinal
and shear modes against independent discrete Fourier symbols, exact planar
extrusion, zero-law Euler recovery and source/work refusal, and six-wall total
impulse with late cancellation/replay. Thermal L2 errors for 6³/12³/18³ cells
on an anisotropic physical box are 4.45971174e-4, 1.13831840e-4 and 5.07037363e-5,
giving observed orders 1.97 and 2.00. These retain source-owned Cv and all
velocity components. They do not establish gas-specific molecular calibration,
physical shock thickness, turbulence or experimental validity. Composition
remains spatially uniform and frozen; 3-D species diffusion, reactions and
radial/no-slip boundary-layer physics require their own admitted evolution.

### Evolving composition on the shared Cartesian discretization

`compressible::mixture::MixtureEuler2d` advances every canonical species mass
density, x/y momentum and full species-referenced total energy on the existing
stationary `CartesianGrid`. Density is derived from the species equations.
The original `MixtureGas` owns composition-dependent calorics, formation-energy
references and planar HLLC/HLL/wall closures. No separate equation of state,
passive-scalar wave solver, trace seed or species normalization is introduced.

Frozen and evolving-species consumers share the private Cartesian MC/common-
theta reconstruction, face topology, combined CFL, SSP RK2, retry, conservative
divergence and boundary-integration core. Concrete views are monomorphized.
Mixture traces reconstruct density, every original mass fraction, velocities
and pressure. The source-owned common composition limiter preserves uniform
inert fractions and linear composition constraints through density waves;
independent partial-density limiting is not used. One source-admission
contraction scales the complete trace. The additional density slope on each
axis is included in the declared numerical workspace.
Flat buffers are prepared once under a checked total byte budget, including
caller state-vector excess capacity, and a 1 GiB hard ceiling. Successful steps
allocate no per-cell or per-step vectors. Original source storage and bounded
stack scratch are separate from the numerical-buffer budget.

Accepted reports retain per-species mass, vector momentum, energy and per-face
outward transfers. Physical-energy residuals use the original species-reference
transform without subtracting two large absolute formation-energy totals.
Failed source, budget or cancellation checks publish no cells, clock, report or
boundary impulse. Six added native tests exercise oblique material-contact
refinement, evolving local calorics, source-data-dependent shocks, original 1-D
and frozen-2-D compatibility, axis rotation, closed-wall balances, complete
workspace admission, late cancellation/replay and shear-heating source refusal.
All 18 compressible tests, including the existing molecular path, pass.

This field is Estimated stationary inviscid evolving-composition transport.
The reactive and molecular adapters below compose chemistry and species/heat/
momentum diffusion. Moving domains, exact pressure preservation at moving unequal-caloric contacts,
fully discrete entropy inequalities and experimental validation remain outside
this lane.

### Finite-rate chemistry in the same two-dimensional species field

`compressible::mixture::reactive::ReactingEuler2d` composes the original
`fs-gas::multicomponent::reactive::CellChemistry` with `MixtureEuler2d`.
One original source owner converts conservative species densities to amounts,
advances the prepared fixed-volume reactor and restores the conservative row.
Both momentum components enter kinetic energy and remain fixed during local
chemistry. The same species energy references convert chemical release into
sensible energy exactly once; there is no additional heat-release source,
species renormalization, trace seed or Cartesian reaction-rate implementation.
Source identity and the complete admitted caloric domain must match exactly.

Each accepted interval is chemistry(dt/2), shared Cartesian advection(dt),
chemistry(dt/2). The post-reaction face CFL still owns the admissible acoustic
step. If it shortens the interval, the whole trial restarts from the original
field before either chemistry half step is accepted. The two existing full
Cartesian workspaces, original prepared reactor and combined boundary buffers
are reused. Their combined numerical heap capacity, including caller state-vector excess
capacity, is checked against an explicit byte budget and the existing 1 GiB
ceiling. Immutable sources/prepared caloric coefficients and bounded stack
scratch are excluded. Accepted steps allocate no new numerical arrays.

The coupled report includes both chemistry halves, the original transport
report, per-species source changes, both momentum balances, referenced energy,
total mass and physical chemical/thermal/kinetic energy. Physical energy is
formed from conservative differences and original reference changes, avoiding
subtraction of large absolute totals. Chemistry has no boundary flux; without
molecular transport the original advective impulses are the complete boundary
exchange. With it, the reactive owner integrates both molecular half intervals
and advection into its combined loads. The borrowed field's Euler metadata describes only the middle
advection phase; the reactive owner's `last_step` describes the full interval.
Failed source, chemistry-work, coupling-retry or cancellation checks preserve
all accepted cells, clock, both reports and face impulses.

Five focused G0/G1/G3/G4 regression groups cover a moving exothermic field
against the original homogeneous reactor, heat-release-induced CFL retries,
formation-energy-data changes reaching pressure and sound speed, oblique
reaction/advection refinement against an analytic reversible-isomer solution,
the original reacting 1-D owner and x/y rotation, closed-wall impulses,
final-publication cancellation/replay, and total-workspace/source/work refusal.
The manufactured sources establish numerical composition and accounting, not
a calibrated chemical mechanism or measured combustion behavior.

The original reactor's selected backward Euler or SDIRK2 method is retained.
Symmetric composition alone does not establish coupled temporal order;
limiting and first-order retries may also reduce observed order. No fully
discrete entropy theorem follows from splitting. This rung remains stationary;
the explicit molecular extension below supplies viscous/thermal/species
transport. Flame-speed/thickness calibration, turbulence, moving boundaries
and experimental validation remain separate work.

### Evolving-species molecular transport and reacting diffusion in two dimensions

`compressible::mixture::molecular::MolecularTransport2d` prepares the original
`fs-gas::multicomponent::dissipation::PlanarTransport` on the exact field source
and caloric interval. The source owner supplies mass-averaged Maxwell--Stefan
species diffusion, original partial enthalpies relative to the same species
energy references, Fourier heat and the complete planar Newtonian stress.
Normal and tangential momentum and full viscous work share the conservative
energy flux. Neither chemistry nor viscous dissipation adds a second heat source.
Every actual face temperature, pressure and changing mole composition reaches
the original constant, kinetic or material transport law.

The existing Cartesian face traversal, reflected/periodic neighbor lookup,
conservative divergence and RK boundary integration are reused. Primitive
inversion and four velocity derivatives are cached once per cell/stage. The
source kernel re-admits face primitives without repeating energy inversion.
One stored face flux transfers equal and opposite amounts to its neighbors;
periodic duplicates copy the identical bits. The first prepared face pass also
sets the initial parabolic interval and supplies the first SSP stage, avoiding
a duplicate pass on successful substeps. Both stages enforce the original
source domain and a conservative engineering diffusion-scale estimate. This
is not an unconditional positivity, entropy or stability theorem.

Periodic axes and stationary adiabatic impermeable slip walls are admitted.
At a wall, both normal velocity and tangential divergence contribute to its
normal viscous traction; species, tangential traction, heat and wall work are
zero. A transmissive Euler end without the explicit source-owned reservoir
below still refuses molecular evolution: zero-gradient wave continuation
supplies no molecular composition/temperature boundary law. Normal support
traction includes signed `p-tau_nn` when composed with advection.

The operator advances a complete requested duration through bounded SSP RK2
substeps without changing the Euler clock or advection metadata. Its field,
report and integrated face loads publish only after the final cancellation
check. Flat state/stage/face/gradient/MS/boundary buffers are prepared under an
explicit combined field/operator byte budget and the existing 1 GiB ceiling.
Caller state-vector excess capacity is counted and follows its actual owner
during O(1) buffer exchanges. Immutable source/prepared coefficients, bounded
stack scratch and supplied transport-law storage/evaluation allocations are
outside this numerical-buffer budget. No whole-process/provider memory or
generic material-law allocation-free claim is made.

`ReactingEuler2d::step_with_dissipation` composes D(dt/2), C(dt/2), A(dt),
C(dt/2), D(dt/2). Heat-driven acoustic shortening restarts all five phases at
one common interval. Final molecular heating must retain the chemistry source
temperature interval. One explicit combined budget covers both owners.
Complete physical loads and source/energy balances include all phases; the
molecular owner's own last report describes the second half interval.
Rejected complete trials preserve both owners' accepted reports and loads,
including cancellation after the final molecular half interval.

Seven added G0/G1/G3/G4 regression groups exercise oblique species-diffusion
refinement with unequal heat capacities and formation references, data-only
diffusivity changes, independent oblique vector stress symbols and viscous
heating, the original 1-D molecular owner and x/y rotation, stationary wall
impulses, complete memory/excess-capacity admission, transmissive refusal,
requested-duration work/cancellation/replay, uniform reacting-field recovery
with heat-driven CFL restarts, and atomic full reacting/molecular wall loads.
These manufactured cases do not establish a measured chemical mechanism,
flame speed or thickness, Soret/Dufour transport, no-slip boundary layers,
rarefaction, turbulence or experimental validation.

### Held source reservoirs for evolving-species planar gas

`MixtureEuler2d::new_with_reservoirs` binds one immutable `ReservoirBoundary`
to each selected transmissive side. A boundary supplies the field's canonical
partial mass densities, Cartesian velocity and absolute pressure. The original
`MixtureGas` admits temperature, caloric precision, sound speed and species
energy references. Duplicate sides, periodic/reflecting sides, foreign species
dimensions and inadmissible original-source states refuse before publication.
Unbound sides retain their existing Euler condition; every open side must have
an explicit admitted bath before molecular evolution is allowed.

The same held state enters the original HLLC exterior Riemann trace and the
original Maxwell--Stefan/Fourier/Newtonian molecular face. The wave solve
determines actual inward/outward transport and backflow from both traces;
neither a prescribed species flux nor an inflow-only switch replaces it.
Exterior wave speeds also participate in the initial interval estimate.
The bath lies at the physical boundary face, so the molecular normal distance
is half a cell width. Cell velocity derivatives use reflected derivative
neighbors about that same boundary velocity. A uniform bath has zero
tangential derivatives; cross derivatives use the shared two-trace average.

The original two-trace midpoint approximation supplies constitutive state,
Maxwell--Stefan coefficients/gradients, traction and Fourier conduction.
`PlanarTransport::reservoir_face_from_primitives` supplies boundary energy
effort from the admitted bath itself: stress work uses its Cartesian velocity
and species energy uses its original referenced partial enthalpies at the
bath temperature. The adapter adds no second enthalpy table, work correction
or viscous heat source. A stationary bath therefore supplies zero traction
work, and prescribed linear Couette shear supplies the correct mechanical
power independent of the number of cells. Boundary-layer refinement remains
necessary for general nonuniform states and constitutive laws.

Signed species, both momentum components and referenced/absolute energy reach
the existing integrated face loads and conservation reports. The reacting
owner retains the exact held source metadata in its independent full-step
trial and includes both molecular halves plus advection in complete boundary
transfers. Reservoir numeric rows and primitive caches count toward field,
reacting-copy and combined molecular byte admission. Cancellation or refusal
preserves both owners' accepted cells, clock, reports, loads and bath states.

Five regression groups exercise supersonic composition inflow/outflow in both
directions on both axes, refinement to an analytic prescribed-composition
diffusion solution with unequal calorics and formation references, exact
one-cell Fourier heating and Couette boundary work, original-source/topology/
byte admission, late cancellation and replay, and coupled reacting boundary
species/momentum/energy balance. These are externally maintained infinite
reservoirs. Finite bath storage and feedback, time-varying/nonuniform boundary
data, characteristic-pressure-only admission and nonreflecting acoustics are
not implemented by this boundary model.

### Volumetric evolving-species shock fields

`compressible::spatial::mixture::MixtureEuler3d` advances the original canonical
species partial densities, three Cartesian momenta and full referenced total
energy on the physical 3-D grid. Density is their compensated species sum;
local gas constant, heat capacity, pressure and sound speed follow the evolving
original composition. Both transverse momenta participate in the same
source-owned HLLC/HLL energy and admissibility checks as the normal momentum.

The frozen and evolving-species 3-D fields now share the original volumetric
topology, MC/common-theta six-face reconstruction, face traversal, combined
numerical CFL, conservative divergence and RK boundary integration. The
existing molecular field consumes that same conservative update. No second
spatial shock solver or separate passive-scalar species transport is introduced.
All axes use physical face areas and cell volumes; periodic endpoints reuse
the exact stored flux. Source refusal contracts traces or retries the whole
step, without clipping cell averages or renormalizing species.

Density and all mass fractions use the same original compatible composition
reconstruction as the 1-D and 2-D owners. One common MC factor scales the
centered fraction vector, retaining linear composition constraints rather
than independently limiting species or assigning a remainder species. A
relative-only roundoff guard keeps constant-fraction division noise from
flattening resolved trace-species gradients. The common source-admission
contraction includes density, every fraction, all velocities and pressure.

Flat accepted, stage, slope, face, primitive and boundary arrays are prepared
once under an explicit checked byte budget and the shared 1 GiB ceiling.
Excess capacity in the supplied state vector is included. Source-model storage
and bounded stack scratch remain separate. Successful steps allocate no
per-cell vectors. A final cancellation check precedes publication of field,
clock, report and loads. Reports include species masses, three-axis momentum,
energy and original per-face transfers; the physical-energy residual applies
the source reference transform to differences rather than subtracting two
large absolute formation-energy totals.

Seven focused G0/G1/G2/G3/G4 tests cover original 1-D shock equivalence on each
axis, full two-dimensional extrusion, oblique species contact refinement,
closed-wall species/energy/support impulse, complete workspace and last-poll
cancellation/replay, and combined two-shear source refusal. The smooth contact
L1 errors at 12³ and 24³ are 2.48975629e-3 and 5.64993110e-4 (observed order
2.140). All twelve pre-existing volumetric Euler and molecular tests retain
their results, including the previously reported contact and heat errors.

An exact periodic advection case varies density and two species while holding
a third mass fraction uniform. Independent partial-density limiting previously
produced a 1.30e-3 local inert-fraction error after five steps on 8³ cells,
despite global conservation. Compatible reconstruction reduces that error to
5.56e-17; the corresponding 2-D case gives 1.12e-16. The checks also retain
pressure, velocity, original calorics and global species/physical energy.

The `mixture_shock_box_3d` example uses the original pinned nitrogen and oxygen
NASA9 sources, a hot nitrogen-rich pocket in an oxygen-rich surrounding gas,
and six reflecting walls. It writes all three velocities, temperature,
pressure and both species fractions, and reports actual support impulse and
species/physical-energy balances. A 12³ run completed 169 accepted steps over
1.713309878178e-3 s, maximum combined CFL 0.222233 and 1,078,272 numerical
workspace bytes. Both species masses and physical energy had zero printed
residual; fluid/support momentum residuals were at most 3.89e-16 kg m/s.
The exported nitrogen fraction remains within its original 0.25--0.95 range;
this observation is not a fully discrete maximum-principle certificate.
This stationary inviscid species owner does
not itself integrate chemistry or molecular species diffusion, preserve
pressure exactly at moving unequal-caloric contacts, move the geometry, or
establish an entropy/invariant-domain or experimental combustion claim.

### Finite-rate reactions coupled to three-dimensional compressive waves

`spatial::mixture::reactive::ReactingEuler3d` composes the original volumetric
species field with `fs-gas::multicomponent::reactive::CellChemistry`. The
source-owned fixed-volume reactor retains every species, original mechanism,
NASA caloric law and formation reference. All three momenta remain fixed in
each local chemical update. Changes in chemical reference energy enter the
stored sensible gauge exactly once; temperature, pressure and sound speed then
feed the next original conservative wave solve.

The original `step` composes chemistry(dt/2), unsplit waves(dt), chemistry(dt/2).
`step_with_dissipation` adds the original molecular operator before and after
that split: D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2). Actual evolving species,
temperature and pressure reach the original molecular law. If heating changes
the acoustic interval, all phases restart from the accepted field at the same
shorter duration. Source or work refusal retries only within explicit bounds.
The final molecular temperature must also remain inside the chemistry source
interval. A final poll precedes field/report/load publication; cancellation
preserves both flow and molecular owners' accepted reports and face transfers.
Chemistry is local. Complete boundary exchanges sum both molecular halves and
the wave stage; neither endpoint pressure nor the middle wave load substitutes
for the full interval. The molecular owner's own report describes the second
half, while the outer report retains both halves and their physical balances.

Two prepared complete field workspaces and one original prepared reactor are
reused. The total admitted numerical byte budget includes both fields, excess
input-vector capacity, combined boundary buffers, actual kinetics/Newton/
chemistry buffers and the molecular operator when enabled, under the shared
1 GiB ceiling. Actual capacity follows O(1) state-buffer exchanges between all
three owners. Immutable source data, prepared caloric
coefficients, fixed-size owners and bounded stack scratch are excluded.
`with_chemistry` consumes an already prepared source without rebuilding it.
Accepted source reductions are compensated in fixed cell/species order.
The complete report retains original transport and chemical diagnostics,
per-species reaction changes, three-axis momentum, mass, referenced energy and
physical chemical-plus-kinetic energy residuals formed from differences.

Four focused numerical tests compare a moving reacting field to the original
homogeneous reactor, force heat-release-induced acoustic restarts, refine a
three-axis reaction/advection solution against analytic reversible decay,
launch compressive waves from an exothermic pocket while closing chemical
energy and wall support, and exercise complete workspace/source/work refusal
with last-publication cancellation and replay. Manufactured mechanisms test
the conservative numerical composition and do not calibrate real combustion.
The analytic reaction/advection L1 errors at 12³/24³ are 1.14633243e-3 and
2.79239902e-4, giving observed order 2.037; the test requires order above 1.8
in this resolved range. The earlier 6³ sample is outside that refinement claim.

The `reacting_shock_box_3d` example retains the original pinned ten-species,
29-channel hydrogen/oxygen mechanism and nitrogen diluent, starting with
exactly absent radicals. An initially hot region reacts and drives the same
3-D gas field. Explicit resolution, chemistry, time, step and memory budgets
bound the run. Its CSV includes every species and all solved thermodynamic
and velocity fields; reported balances independently reconstruct original
element inventories, physical energy and fluid/support momentum.

Before the volumetric molecular extension, the default 4³/100 µs run completed
33 accepted intervals and 4,224 cell solves
with maximum combined CFL 0.339843 and 221,800 total numerical workspace bytes.
All 64 accepted cells were finite and species nonnegative; final temperatures
were 1200.633–2438.868 K, pressures 167869–294049 Pa and maximum speed
137.135 m/s. Printed mass, physical-energy and all fluid/support momentum
residuals were zero; the largest elemental residual was 5.21e-18 mol.
The example explicitly reports local inert-nitrogen preservation: maximum
fraction drift fell from 6.90e-4 before compatible reconstruction to 4.45e-16
afterward. All 89 focused flux-library tests passed in that earlier planar
molecular/reacting baseline. These are coarse numerical source/solver results,
not experimentally validated flame or detonation predictions.

Appending `molecular` to the example arguments enables the original 3-D
Maxwell--Stefan, Fourier and full viscous-stress/work owners with state-dependent
kinetic coefficients. All ten Lennard-Jones records come from the same pinned
Cantera mechanism revision. The nonpolar Chapman--Enskog/Wilke/Eucken estimate
explicitly omits polar H2O corrections, rotational relaxation and Soret/Dufour
effects, and declares zero bulk viscosity. The source-owned explicit
`2e-8 J/(mol K)` Cp-fit deficit allowance admits the retained atomic-H fit
without clipping its original Cp or changing its conductivity calculation.

The native molecular 4³/100 µs run completed 33 coupled intervals, 66 molecular
substeps, 4,224 chemistry cell solves and 40,639 accepted source steps without
coupling or molecular retries. Its complete numerical workspace was 345,688
bytes; maximum wave CFL and parabolic number were 0.339861 and 0.005563510.
All 64 exported cells were finite with nonnegative species. Mass and physical
energy residuals were -5.421011e-20 kg and 5.684342e-14 J; each fluid/support
momentum residual was at most 8.68e-19 kg m/s and each element residual at most
3.47e-18 mol. Final temperatures were 1201.141--2438.854 K. Local nitrogen
fraction change was 4.574416e-7: differential diffusion can change an inert
species' local fraction. This coarse run exercises the complete physical
consumer and its balances; it does not resolve flame or shock thickness.

The selected backward Euler or SDIRK2 source method remains explicit.
Symmetric ordering alone is not a coupled temporal-order or entropy proof;
shock limiting and first-order retries may reduce observed order. This
stationary reacting rung does not supply moving-domain GCL/remap, wall heat
transfer, physical flame/detonation thickness,
turbulence, reciprocal FSI or experimental ignition/detonation validation.

Three additional reacting molecular regression groups recover the original
uniform chemical solution with heat-driven CFL restarts; refine a complete
three-axis reaction/advection/diffusion solution against analytic cell averages;
and retain every phase's support load, complete source/energy balances, excess
capacity, work refusal and final-cancellation replay. The analytic L1 errors
on 12³/24³ cells are 8.09463491e-4 and 1.04199187e-4. This observed refinement
does not imply a general third-order method. The integrated source suite passed
207 gas tests and all 98 focused flux tests (20 spatial-mixture plus 78 other
cases). These runs preceded the opt-in Cp-fit allowance; all nine kinetic
transport tests and the complete molecular example passed after that change.
The repository-level DSR gate remains a separate unavailable check.

### Three-dimensional evolving-species molecular transport

`spatial::mixture::molecular::MolecularTransport3d` integrates the original
`SpatialTransport` source over the existing physical Cartesian grid. Every
species diffusive flux, all three viscous momentum fluxes, species-referenced
partial enthalpy, Fourier heat and full stress work enter one conservative
divergence. Actual local composition, temperature and pressure reach the
original transport law; pressure diffusion remains in the source's original
Maxwell--Stefan driving force. No frozen mixture coefficient or additional
viscous heat source replaces these terms.

Nine velocity derivatives are cached once per cell/stage. A face consumes both
transverse derivative directions with the complete three-dimensional Newtonian
trace and shear terms. Topology, physical areas/volumes, conservative updates
and RK boundary quadrature reuse the original spatial owner. One stored face
flux serves both adjacent cells; periodic copies retain identical bits. The
initial gradient/source pass supplies both the first SSP stage and its step
estimate, avoiding a duplicate successful pass.

Planar and volumetric molecular fields share one requested-duration controller
for timestep selection, exact final-remainder accounting, bounded retries and
accepted-stage extrema. Conservative stage updates and private face-load
buffers remain with the original field owners. A failed attempt reconstructs
the initial source faces before retrying; cancellation or exhausted work
leaves accepted fields and published loads/reports unchanged.

The prepared operator advances an exact requested duration using bounded SSP
RK2 substeps, original source admission and an engineering parabolic estimate
covering all axes, cross-stress stencils and wall gradients. This is not a
spectral, unconditional-positivity or entropy certificate. Field values,
integrated boundary transfers and reports publish after the final poll; Euler
time and its advective metadata remain owned by the composing flow integrator.
Actual state capacity follows O(1) exchanges between owners. All numerical
state/stage/face/gradient/boundary/MS buffers are counted under the explicit
combined field/operator budget and existing 1 GiB ceiling; immutable sources,
bounded stack scratch and supplied-law storage/allocations remain excluded.

Periodic axes and stationary adiabatic impermeable slip walls are admitted.
Wall normal stress includes both transverse dilatations; species, both shear
tractions, heat and stationary-wall work are exactly zero. Transmissive Euler
ends refuse because they do not define a molecular reservoir law. Soret/Dufour
effects, no-slip layers, moving geometry, dense-gas transport and calibrated
shock/flame thickness remain separate capabilities.

Focused checks use independent oblique three-dimensional species-diffusion,
vector-stress and Fourier symbols, physical energy converted from damped
kinetic motion, the original planar owner under extrusion/axis rotation,
closed-wall momentum balances, exact absent species, source/budget/work refusal
and last-publication cancellation with deterministic replay.
The oblique diffusion L1 errors on 8³/16³ cells are 6.45856775e-4 and
1.69469136e-4, respectively, for an observed order of 1.930. Full vector stress
is checked against independent discrete longitudinal/transverse symbols, and
the integrated thermal gain equals lost kinetic energy without a heat source.

## Two-fluid capillary displacement

`capillary::step_newtonian_two_fluid_circular_meniscus` adds the displaced
phase's Poiseuille resistance to the same horizontal rigid-tube model:
`8 [mu_l x + mu_d (L-x)] dx/dt = R^2 [delta_p + 2 sigma cos(theta)/R]`.
The affine resistance integrates exactly. A normalized, rationalized quadratic
inverse avoids subtracting nearly equal roots; inlet/outlet event times are
checked before evaluating that inverse. The kernel is allocation-free, with
bounded scalar work and no time-stepping iteration. Zero displaced viscosity
delegates exactly to the existing Washburn owner; equal viscosities give
constant speed. Positive displaced viscosity gives finite initial speed,
without claiming to resolve inertial startup.

The actual signed volume determines pressure work and the change in wetted-wall
free energy. Positive analytic resistance weights integrate each phase's
viscous loss over the entire displacement, including drainage. Their sum is
the primitive's total loss, `pressure_work - wall_surface_energy_change`.
The translating meniscus has constant area. Formation/exit topology, gravity
in this two-fluid rung, wetting films, inertia, gas compression and dynamic
contact angles remain outside its domain. The displaced phase must be
incompressible or explicitly approximated by uniform frozen gas properties;
it vents to the specified outlet pressure and is never a trapped gas pocket.

Seven G0/G1/G3 unit tests cover the integrated Poiseuille relation, finite
startup, equal-viscosity and original-law limits, signed split-step composition,
end events and independent spatial quadrature of the phase losses. These are
numerical checks of the reduced model, not experimental validation. The
two-fluid resistance follows Budaraju et al. (2016),
<https://doi.org/10.1021/acs.langmuir.6b02788>.

## Material-bound circular capillary transport

`capillary::material::MaterialCapillaryTube` binds the original
`ResolvedLiquidSurfaceState` and `ResolvedLiquidWettingState` to the existing
horizontal/inclined meniscus operators. Density, dynamic viscosity, tension
and equilibrium angle come from those same immutable owners; callers provide
only geometry, current wet length, pressure, gravity and requested time.
Admission checks the ordered interface's liquid identity and complete query
point, refuses neglected capillarity and unowned rheology evolution, and
retains the exact source owners in every accepted result. A kinematic-viscosity
source uses its own density to derive dynamic viscosity. Signed liquid mass
uses that same density; wall energy, gravity work and end-event time remain
the primitive operator's actual transfers.

The horizontal two-fluid method consumes one checked `FluidTransport`
evaluation and retains its original owner, both signed mass transfers and the
two separately integrated viscous losses. The displaced phase has frozen,
homogeneous properties and vents to a prescribed outlet pressure. Passing a
gas owner does not activate compressibility, sealed-gas storage or inertia.
The caller must assign the displaced phase to the fluid pair admitted by the
tension/wetting sources; this generic transport view cannot establish that
equivalence. There is no inferred source data, dynamic contact-line law,
thermal evolution or flow beyond an inlet/outlet event. A raw transport
fixture remains unsourced. These bounded scalar bindings allocate no memory
during an accepted step and contain no internal cancellation point.

Focused G0/G1/G3 tests in `tests/material_capillary.rs` exercise manufactured
Washburn/Jurin motion, physical material swaps, kinematic-viscosity conversion,
signed mass/energy and end-event accounting, refusal of source mismatches and
unevolved rheology, and finite displaced-fluid resistance with original-owner
retention. They establish numerical/source integration, not measured wetting
accuracy or arbitrary two-phase applicability.

## Finite-compartment osmosis and selective membrane networks

`osmosis` (the existing `continuum` feature) evolves the original
`fs-material::osmosis::IdealDiluteSolution` inventories across prepared
`MembraneTransport` links. Solvent volume and every independent solute amount
are extensive state; concentrations and osmotic pressure are recalculated at
each trial. Species share one explicit ordering, temperature and incompressible
solvent molar volume. Mechanical pressure has an independent linear volume
response, admitting prescribed pressure, hydrostatic columns and compliant
chambers without counting osmotic pressure twice.

At a proposed water transfer, the original exponentially fitted membrane
coefficients reduce all backward-Euler solute equations analytically to
positive, conservative endpoint inventories. A bounded scalar solve closes
the water transfer. Near-complete donor depletion preserves the small donor
remainder using positive implicit weights. The bracket grows from zero along
the initial driving direction; a distant inadmissible volume-floor state does
not veto an admissible interior root. General leaking-solute residuals do not
carry a uniqueness claim.

Membrane sweeps follow their declared fixed order. The first-order Lie split
compares one full sweep with two half sweeps, retains the latter, and controls
both water and every species with separate SI absolute tolerances plus a
relative tolerance. Maximum substep, accepted steps, rejected trials and scalar
iterations bound work. This is a local temporal error estimate, not a certified
global error bound. A volume floor is an explicit model boundary; depletion is
refused rather than clipped. Each face exchange is applied conservatively to
both compartments, and reports reconstruct closed water and species balances.

Four complete state buffers are prepared under the caller's byte budget;
successful advances allocate no memory. Cancellation polls each membrane and
scalar evaluation and once before publication. A refused or cancelled full
interval leaves public compartments and clock unchanged; retries reconstruct
all scratch state from the original accepted endpoint. Determinism follows
fixed cell/species/link traversal and the original material's deterministic
elementary functions. Workspace accounting includes retained input capacities
and original boxed coefficient/inventory arrays; bounded scalar stack scratch
is additional.

Focused G1/G3/G4 checks cover analytic unequal-volume diffusion, compliant
pressure transients and hydrostatic equilibrium, retained-solute equilibrium,
nonselective convection, leaking-species free-energy decrease, nonlinear branch
selection, graph refinement, stiff donor remainders and atomic failure/replay.
`cargo run -p fs-flux --example osmotic_columns` exercises two finite columns:
water transfer changes dilution and hydrostatic pressure until these balance.
Its explicit numerical coefficients are illustrative, not a calibrated source.

This owner assumes well-mixed, isothermal dilute compartments with negligible
solute volume and fixed membrane area/coefficients. It supplies no electrolyte
activity law, Donnan potential, resolved boundary-layer concentration
polarization, nonisothermal heat/entropy balance, mechanical inertia, ionic
reaction chemistry or experimentally validated membrane selectivity. Ordinary
nonselective fluid links use zero reflection; gas solution-diffusion retains
its separate existing physical law and permeability units.

## Material-driven spherical cavitation

The `cavitation` feature (included by `continuum`) exposes
`cavitation::RayleighPlessetBubble`. It consumes the existing material owner's
`CavitationState` and an explicit fixed polytropic gas inventory. Absolute
liquid pressure, including metastable tension, drives spherical radius and wall
velocity through vapor, gas, capillary and viscous stresses. The state and
material definition are separate. No pressure floor, injected nucleus, gas EOS,
minimum-radius rebound, or cavitation sound preset is inferred.

The numerical owner uses a two-component Dormand--Prince 5(4) pair, fixed stack
storage and endpoint-derivative reuse. Radius and velocity have separate SI
absolute tolerances; maximum step, minimum step and trial-attempt count are
explicit. An out-of-domain intermediate stage can retry with a smaller step,
but cannot relax the admitted radius interval or wall-Mach limit (at most
0.1). Cancellation polls before every RHS and publication. Failed or cancelled
requests leave the caller's entire accepted state unchanged. Identical inputs,
callbacks and call boundaries replay deterministically on the same ISA.
Embedded local-error ratios are numerical indicators, not global enclosures;
the caller resolves the pressure waveform and splits known discontinuities.

Accepted observations expose geometric volume, its first and second time
derivatives, liquid kinetic and surface energy, and instantaneous nonnegative
viscous dissipation. Volume acceleration is a physical acoustic source input,
not an observer pressure. The spherical infinite-liquid incompressible model
excludes acoustic back-reaction, shocks, nonspherical collapse/erosion, nearby
walls, interacting bubbles, gas diffusion, and evolving interfacial temperature.
It uses the small vapor/liquid density-ratio approximation to interface
kinematics: supplied vapor density must satisfy the caller's explicit ratio
limit, which cannot exceed 0.1. Without optional vapor density this remains a
caller-declared negligible-mass assumption. Geometric bubble volume is not
silently interpreted as vapor mass or deposited latent heat. A
finite-compressibility or thermally limited
successor must supply those equations explicitly.

Focused G0/G1/G3/G4/G5 sources check exact equilibrium, independent Rayleigh
collapse quadrature and first-integral energy, the damped Minnaert limit,
tensile-pressure growth, mechanical power and acoustic-volume identities,
radius/Mach/work refusals, and cancellation with bit-identical retry. Physical
equations follow [Brennen, chapters 2--4](https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm);
the independent collapse oracle uses
[Obreschkow et al., equations 3--8](https://arxiv.org/abs/1205.4202).
These checks establish mathematical-model behavior, not experimental accuracy
or a resolved multiphase CFD/film-cavitation result.
