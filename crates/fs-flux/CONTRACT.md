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
  supply a thermal field, or add viscoelastic memory. The nonlinear adjoint
  of the constitutive update is not provided by the linear transpose solve.
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
  fixtures. Full nonlinear constitutive/convection Jacobians and time-dependent
  adjoint checkpointing are not supplied by the assembled linear adjoint API.
- Nonlinear/backflow and transient open-boundary conditions. Steady linear
  Stokes velocity/total-traction boundaries are provided by the API below.
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
one StokesBoundary owner for every mesh boundary edge. It uses the existing
BDM1/P0 symmetric-stress discretization with constant kinematic viscosity:
-div(2 nu D(u)) + grad(p) = f and div(u) = 0. Velocity edges constrain normal
moments strongly and receive tangential SIP terms. Natural edges receive
t = (2 nu D(u) - p I)n and retain unknown normal velocity moments. They have
no velocity penalty or extra pressure-boundary matrix term.

One connected fluid domain and at least one velocity edge are required.
Duplicate, missing, interior and out-of-range ownership is refused before
callbacks run. Boundary quadrature samples are retained once for consistent
strong moments and SIP. Non-finite loads and assembly, incompatible net flux
for all-prescribed velocity, singular factorization and failed explicit
original or row-equilibrated residual admission return StokesBoundaryError
without a solution. No hidden grounding is added for pure-traction rigid modes.

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
examples/pressure_channel.rs is runnable with
cargo run -p fs-flux --example pressure_channel.

This extension covers steady linear, constant-density, constant-viscosity
Stokes on the existing straight-sided 2D triangulation. It does not add
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
