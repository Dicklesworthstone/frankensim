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

## Laminar flow through external hydraulic resistances

`PreparedStokesPressureResponse::solve_resistance_ports` couples the resolved
native flow to `p_i = p_reservoir_i + R_i Q_outward_i` at each prepared port. Pressure
is an increment to the original total traction; tangential traction remains
unchanged. The full cached port-response matrix couples every port, so only a
small dense system is factored for new reservoir/resistance loads. It rechecks
the original continuum equations and each port law at the returned field.
Zero resistance prescribes reservoir pressure; reverse flow reverses the
pressure drop and preserves nonnegative external `R Q²` loss.

Homogeneous and frozen cell-material entries share this native owner for
steady Stokes, BDF1 startup and variable-step BDF2. Source density consistently
normalizes pressure, resistance, inertia and power; each material retains its
own dynamic viscosity. Physical resistance is Pa s/m² for volume flow per
depth, pressure is Pa and external dissipation is W/m. External heat is kept
separate from resolved bulk viscous heat and wall-slip losses. Material owners,
histories and prepared caches remain immutable on success or refusal.

G1/G3 tests in `stokes_boundary/flow_control.rs` and `material_flow/density.rs`
compare reversible transient plug momentum and pressure work, an independent
full-field Robin assembly with interacting ports, and exact layered Couette
stress/heat. Scope remains fixed-mesh linear Stokes; every traction port may
now be controlled, with pressure selected by the reservoir laws. The external
resistance law itself has no compliance,
inertia, turbulence or automatic thermal deposition; frozen nonlinear operators
remain frozen. This does not claim full hydraulic-network or FSI closure.

These regressions are authored but were not executed. The original single-worker
Cargo build was killed by signal 9 while compiling the pinned asupersync
dependency, before it reached the changed FEM/material modules; no source/type
diagnostic was emitted. Source parsing and whitespace checks passed before the
final concurrent-change merge, and bounded numerical/API review found no
concrete issue. Native integration compilation and execution remain unverified.

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
  The Mach envelope covers fluid-relative center translation and all surface
  rotation through the conservative bound
  `(|v_center-v_gas| + |omega| hypot(radius, thickness/2))/c`; a translational
  or advancing-rim violation cannot bypass admission because rotational drag
  is absent. Vacuum and fully inactive contributions retain their exact zero
  limit. This bound may refuse a state whose exact maximum surface speed is
  lower; it does not extend the incompressible correlations into compressible flow.
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
  identically. Pressure pinned by REPLACING the first continuity row of
  each edge-connected fluid component (its div-freeness is implied by the
  other cells and that component's compatible prescribed net flux — never
  an additive penalty, which would corrupt it). Component topology is
  prepared once and reused by assembly, output and boundary admission.
  Convection: upwinded DG on the SINGLE-VALUED face flux w·n (H(div)
  conformity makes upwinding well-posed), volume term in
  divergence form −∫(u⊗w):∇v. A BDM1 normal trace can reverse within ONE
  face. The shared primal/Jacobian quadrature splits at that linear zero,
  integrating each cubic upwind branch with two Gauss sites; a non-reversing
  face keeps the original two sites, with no per-face allocation. Incoming
  boundary velocity uses the fixed affine trace through its original two
  Gauss samples, consistent with the prescribed normal moments. This also
  keeps the discrete Jacobian correct for non-affine user boundary data:
  sample locations do not move with the upwind split, whose endpoint flux
  is zero. The new local regressions compare piecewise cubic antiderivatives,
  the exact DG jump-dissipation identity, and full residual differences;
  no arbitrary nonlinear boundary-quadrature accuracy is implied.
  `picard` for steady NS, `bdf1_step`
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
  row-equilibrated relative residual. Both dense and sparse solves remove a
  constant pressure origin independently per connected component, accumulate
  `b-A*origin` with double-double arithmetic, and restore the original physical
  pressure afterward. This affine translation leaves the operator, traction
  datum and sealed pressure pins unchanged. Its O(nnz+n) setup prevents a large
  ambient traction pressure from becoming the stopping scale for weak flow.
  The scaled diagnostic also retains the maximum of both translated relative
  defects, and the original equations are checked after restoring pressure.
  The actually returned pressures are translated back for a final residual
  check, so a gradient lost when adding the datum refuses convergence even
  when the internal translated solve met tolerance. This reuses solver scratch.
  Dense LU uses the translated load for prescribed-column elimination and the
  same two residual corrections. Zero origins preserve the original load bits,
  dense arithmetic and diagnostics. After factorization, the dense input buffer
  supplies returned-pressure scratch; the LU owns its independent factors.
  Boundary assembly and returned absolute pressures retain f64 precision;
  translation cannot recover input pressure differences already rounded away.
  Manufactured affine-shear and disconnected hydrostatic regressions in
  `tests/stokes_boundary.rs` compare flow and pressure at large pressure datums.
  A minimal pressure-gradient system separately checks representable and
  unrepresentable differences under the same large datum in both solver routes,
  with a zero-origin case retaining the exact dense pressure difference.
  `converged(tol)` requires both diagnostics to pass,
  and the implicit rheology and physical material wrappers use that guard.
  Dense solves report both diagnostics after residual correction as well.
  `max_iters` is a strict Krylov iteration budget, including the last short
  restart. An approximate factorization breakdown falls back to identity
  preconditioning without changing the equations or asserting convergence.
- `TriMesh::containing_cell` and `FluxSystem::try_velocity_at` locate the
  containing triangle using its oriented half-planes. Velocity samples use
  that cell's affine BDM field, never a nearest-centroid extrapolation across
  a viscosity interface, solid gap or exterior boundary. Shared edges and
  vertices select the lowest-index cell's one-sided trace with a local
  roundoff allowance; tangential jumps are not averaged. Lookup allocates
  nothing and stops at the first containing cell. The checked velocity entry
  returns `None` outside the mesh or for invalid local data; `velocity_at`
  panics in those cases. Neither entry performs a spatial-error estimate.
- `ns::kinematics`: `try_cell_kinematics` and `cell_kinematics` reconstruct
  native affine BDM velocity gradients, signed out-of-plane curl, divergence
  and objective shear rate from the cached basis. Rates are in 1/s for a
  metre mesh and m/s velocity; material density does not rescale them.
  Material-flow callers pass the accepted `velocity_dofs()` to the native
  `FluxSystem` originally bound at construction, retaining its SI geometry.
  H(div) normal continuity does not imply tangential continuity:
  `vorticity_sheet` retains the affine signed tangential jump at both endpoints
  of each interior edge, including zero-mean sheets. Its strength is in m/s,
  not a bulk 1/s vorticity or a square-integrable enstrophy density.
  `vorticity_balance` adds cell and sheet circulation and compares with the
  represented boundary trace. `weak_vorticity` integrates that distributional
  curl exactly against continuous P1 vertex tests. Both use deterministic
  compensated sums in linear, allocation-free traversals. Boundary orientation
  follows fluid outward normals, including holes and disconnected components;
  no outside-fluid extension, artificial wall sheet or prescribed-wall trace
  is substituted. Nonfinite/incomplete queries are refused and trailing pressure
  entries are ignored. Stokes balance is an accounting identity, not a spatial
  accuracy claim; native curl convergence is tested on solved Poiseuille flow.
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
  initializes the solve. Steady constant-viscosity and rheological Picard
  early exit require a small measured velocity update, independently assembled
  endpoint residuals in both original and row-equilibrated equations, and a
  converged last linear solve. A rejected endpoint operator is reused in the
  next iteration. Returned linear residuals and measured relative updates retain
  their original meanings; reaching the iteration budget does not certify
  nonlinear convergence, even when the update is small.
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
  Picard iteration and IMEX BDF1 all use the same symmetric Newtonian stress
  `2 nu D(u)` and normal-cell-height SIP stabilization as the rheological
  material path, including for a spatially constant viscosity. The primitive
  legacy vector-Laplacian APIs remain available separately. Original source
  receipts remain available through `material()`. Invalid inputs, nonfinite
  outputs and failed convergence return `MaterialFluxError` without publishing
  a physical solution. Steady Picard publication checks the update and both
  nonlinear endpoint residuals against its nonlinear tolerance, independently
  of the final linear-solve tolerance. The primitive singular-mesh panic
  contract still applies.
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

Material-bound homogeneous and generalized-Newtonian Stokes/Navier–Stokes
solutions also retain `MaterialViscousDissipation`: cell heat production
`2 mu D:D = mu gamma²` in W/m³ and its area integral in W/m. The same accepted
velocity gradients and retained dynamic material viscosity determine both;
nonlinear and rigid-ALE endpoints use their final cell responses without
source requery or reference-viscosity substitution. The fields are constant
on each BDM1/material cell. Compensated summation accumulates power in mesh order.
Rigid translation and rotation contribute no strain heat. Wall-slip, Darcy,
SIP, upwind and BDF1 losses remain separate. These are endpoint power fields;
time integration and thermal deposition need an explicit consumer, and the
result alone is not a full discrete energy-balance or thermal-evolution claim.
Analytic Couette/extensional, rigid-motion, mesh-refinement, temperature-layered
rheology and moving-chart tests cover the physical values and source ownership.

`FluxSystem::edge_volume_transport` retains both nonnegative directional volume
rates on each BDM1 face, including simultaneous counterflow at zero net flow.
It reuses the momentum operator's affine-trace sign split and retains the signed
net directly from the original mean moment. Rates are per unit depth in m²/s;
consumers supply density, depth, carried state and fixed-inventory continuity
admission. The operator neither repairs divergence nor infers ALE mesh velocity.

### Conservative density and specific-inventory transport

`FluxSystem::prepare_density_transport` binds a prescribed BDM velocity to the
original fixed triangle volumes and an explicit physical depth. It reuses the
original sign-split face arithmetic, including both streams on a reversing edge.
Every cell's incoming and outgoing volume rates must balance within a caller
allowance no larger than 1e-12 relative; no velocity, rate or inventory is repaired.

`PreparedDensityTransport::advance_into` advances positive P0 density and N
finite specific quantities by synchronous donor-cell Euler. A directional mass
transfer is the original donor density times that directional swept volume.
The SAME transfer carries each specific quantity: N=3 may carry physical
[u_x,u_y,h], evolving cell momentum and enthalpy inventories with mass.
Boundary counterflow uses both the original interior state and an explicit
exterior donor state; zero net volume does not imply zero mass or energy exchange.

An outgoing-volume Courant cap strictly below one leaves positive donor mass.
Density and specific-quantity bounds have the usual convex-mixture meaning
up to the admitted volume imbalance and floating-point roundoff. Incoming
specific differences preserve uniform specific fields exactly and retain small
changes beside a common datum. If that difference update suffers severe
cancellation, the same three native face transfers reconstruct the extensive
convex mixture, preserving a small donor's finite inventory at large density
contrast. Compensated product differences audit represented cell inventories;
reported mass/carried balance residuals are not clamped or normalized away.

Preparation and reusable workspace allocation have explicit logical byte caps
and fallible reservations. The warm step is O((cells+edges)*N), allocates no heap
storage, and polls cancellation at cell/edge boundaries. No refusal changes
output cells. Every attempt invalidates the workspace's accepted-transfer view;
`face_transfers` returns both actual directional mass/carried amounts only after
the complete step succeeds. A changed advecting velocity requires preparation
of a new plan; neither density nor specific quantities are frozen into the plan.

This is first-order prescribed-flow transport, with numerical mixing to be
assessed by refinement. It does not solve variable-density BDM momentum or
pressure, update material viscosity/EOS, supply heat conduction or pressure work,
move meshes/interfaces, or represent phase-volume sources. Cell momentum mixing
loss is not silently converted to thermal enthalpy. `tests/density_transport.rs`
contains G1/G3 analytical exchange/refinement, unequal-density counterflow,
open-boundary inventory closure, high-contrast retained-inventory, throughflow
constant-specific, depth/density/datum scaling and atomic-refusal regressions.

### Compatible mass flux for native BDM momentum

`PreparedDensityTransport::prepare_mass_flux` retains the same original
sign-split volume faces and positive donor densities as `advance_into`.
Its per-depth quadrature carries volume [m²/s] separately from donor-weighted
mass [kg/(m s)]; opposing streams are never collapsed into net volume flow.
The plan retains its original normal moments, so a second velocity or mesh
cannot silently supply a different reconstruction.

Each cell reconstructs a vector P2 H(div) mass flux with the original face
moments through degree two and constant divergence equal to net outward mass
rate divided by original area. The remaining zero-normal, divergence-free
bubble minimizes physical L2 distance to the cell's rho times original BDM1
velocity. The explicit reference-triangle construction uses the contravariant
Piola map and fixed positive degree-four quadrature, with no local matrix
factorization or global solve. Geometry, orientation and finite arithmetic
are admitted before publication. Uniform density reproduces the original
affine velocity flux up to floating-point roundoff.

These quadratic moments are necessary for integration by parts against every
BDM1 kinetic-energy trace, including reversing faces. The reconstructed normal
polynomial is a moment projection; original sign-split face samples still own
upwind transport. This does not promote the P0 donor-cell density update to
second-order spatial accuracy. Its constant divergence and the represented
Euler density increment agree up to the original transport roundoff.

Preparation and storage are O(cells+edges), with fallible bounded vector
reservations and per-cell/per-edge cancellation. Warm cell and face sampling
allocates nothing. `tests/density_transport/mass_flux.rs` covers all quadratic
integration-by-parts identities, single-valued normal traces, original BDM1
reproduction, accepted inventory continuity, large density contrast, reference
depth, local corner permutations, rigid frames and refusal boundaries. Native
momentum, pressure and thermodynamic closure remain consuming-owner work.

Pressure/divergence coupling uses the exact integrated BDM face moments:
signed edge length for mean moments and zero for linear moments. The momentum
transpose and retained pressure reactions use the same coefficients, preserving
shared-face cancellation without reconstructed-derivative roundoff. Dense LU
isolates exact prescribed velocity and pressure values before factorization;
its existing two correction solves evaluate residuals of the original CSR
equations with native double-double arithmetic. The original/scaled residual
checks, reaction equations and sparse solver remain authoritative. These changes
do not alter transported fluxes after solving or relax local continuity bounds.

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
- `experimental-sgs`: includes `experimental-compressible` and the shared
  `fs-lbm` local closures and original thermochemical source types. Only an
  explicit `with_sgs` constructor activates subgrid transport on frozen or
  evolving-species three-dimensional gas.

## Conformance tests

The `lc` unit tests reproduce Tesla's printed 925 Hz, 185,000 mi/s,
200-mile-wavelength, 50-mile-quarter-wave example after exact mile-to-metre
conversion. A second test establishes the inverse frequency/quarter-length
relation without adding a device-specific coil model.

`tests/battery.rs`: flux-001 basis exactness + conformity; flux-002
Stokes MMS orders + exact divergence; flux-003 pressure-robustness
identity; flux-004 cavity Picard + Ghia coarse band; flux-005
Taylor–Green BDF1; flux-006 adjoint gradient + annihilation.
The steady shear-thickening regression `flux_024` distinguishes a small velocity
update from a failed nonlinear endpoint balance on the original BDM field.
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
The point-sampling lane adds flux-021 a solved 1,000:1 layered Couette channel
on stretched triangles whose nearest centroid lies in the wrong liquid;
flux-022 deterministic one-sided tangential traces under rigid coordinate
changes; and flux-023 disconnected fluid domains, outside points and invalid
local sample refusal.
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
refusal, translation/advancing-rim/finite-thickness Mach and common velocity
boosts; aero-007 alternative-correlation disagreement; aero-008 exactly-once
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

- The incompressible BDM solver has no LES/RANS closure; its DG dissipation
  is a numerical property. The separate `experimental-sgs` Cartesian gas
  path supplies explicit local subgrid stress/heat, with the limits below.
- 3D (tets, BDM on faces), BDM2+/higher-order pressures, and curved
  boundaries — recorded successors on the same dof-Vandermonde
  surface (bead dcng's high-order FEEC would unblock them).
- Projection/splitting time integrators and BDF3+ remain successors.
  Equal-step mixed-Stokes and fixed-mesh variable-step mixed Navier–Stokes
  BDF2 are available below; IMEX convection remains BDF1.
- Fine-mesh Ghia tables, Re=1000 cavity, cylinder benchmarks, and
  asymptotic-regime convergence confirmation — perf-lane scope,
  ledgered in invariants 2 and 4.
- Mesh-independent iteration counts, multilevel saddle preconditioning and
  production-scale throughput are not established by the sparse finite-element
  fixtures. `linearize_rheology` supplies constitutive/convection state
  Jacobians for prescribed-velocity boundaries. Mixed-boundary nonlinear
  Jacobians and time-dependent adjoint checkpointing remain outside that API.
- Mixed-boundary Newton Jacobians, nonlinear adjoints and high-Re open-boundary
  qualification remain absent. Linear and nonlinear BDF1 velocity/traction
  flow, with explicitly selected natural or dissipative backflow, are provided
  by the APIs below.
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
Couette mass flux, implicit storage, and squeeze through exact old/new gap
values. With the default `NoSlipContinuum` policy, both walls receive the full
no-slip Couette--Poiseuille tangential
traction: wall-on-gas values are `mu*(U_upper-U_lower)/h + h*p_x/2` at the
upper wall and `-mu*(U_upper-U_lower)/h + h*p_x/2` at the lower wall. The
pressure contribution uses the accepted Reynolds face gradient and face gap,
averaged onto each cell. Gas-on-wall tractions have exactly opposite signs;
the old `upper_wall_shear_pa` remains the upper-wall-on-gas alias. Upper-wall
normal traction and gauge load retain their explicit sign/reference conventions.

For no-slip, viscous dissipation includes both `mu*(U_upper-U_lower)^2/h` and
`h^3*p_x^2/(12*mu)`. Pressure-driven heat is integrated on full interior face
intervals and half intervals at the boundaries, before averaging gradients,
so opposing gradients cannot erase dissipation. Closed exterior faces infer
`p_x = 12*mu*U_mean/h^2` from their zero total flux; moving walls do not
silently become zero-gradient boundaries. The wall-work receipt is signed work
by both tangentially moving walls plus normal moving-gap work
`-integral(p dh/dt dx)`; pressure-driven gas can deliver work to moving walls.
A separate `pressure_gradient_power_to_gas_w_per_m` records
`-integral(q*p_x dx)`, where `q` is the accepted mass flux divided by the
same face ideal-gas density. Tangential wall work plus this pressure-force work
equals viscous dissipation up to roundoff. Normal gap work is separate from
that identity. This force-work receipt is not boundary enthalpy transport,
thermostat heat, or a total-energy closure; no gas temperature equation is solved.

Signed boundary/vent enthalpy transport still comes from a caller-declared
specific enthalpy and is not a thermochemical calculation. Gauge load is relative
to a separately declared absolute pressure reference, never implicitly the
initial condition. `Sealed`, two-reservoir `Open`, independently sealed/open
`MixedEnds`, and one named-cell `Vented` topology are admitted. A mixed end
uses `None` for zero mass flux and `Some(p)` for an absolute pressure bath;
its sealed/open kind is fixed across restart while bath pressure may vary.
The receipt retains gas storage, every boundary/vent
mass flux, their global closure defect, pressure load, both-wall shear,
wall/pressure work and dissipation, residual, iteration count, and
continuum-envelope maxima. A checkpoint binds active pressure and previous gap
values to its model/case/authority identity only after a converged step; no
partial state is published.

`ContactExclusionMask` is an explicit contact handoff, not a gap floor: this
first slice admits only a contiguous active prefix followed by excluded cells.
An interior hole, topology change, too-few active cells, non-positive gap or
absolute pressure, EOS-density mismatch, excessive slope/Mach/Knudsen number,
or an unadmitted slip/rarefied/roughness model is refused with `GasFilmError`.
Smooth-wall admission requires a named roughness bound strictly below every
active gap. There is no
anonymous `h_min`, pressure clamp, liquid cavitation substitution, or
Elrod--Adams semantics.

The declared molecular mean free path belongs to the declared temperature and
initial absolute pressure. At fixed temperature the solver uses
`lambda(p) = lambda_ref * p_ref / p`, checking the actual starting/checkpoint
state, the accepted cells, and each imposed reservoir/vent at its local gap.
Reservoir pressure is checked directly, without averaging it with a denser
interior cell. The reported maximum Knudsen number covers the accepted cells
and boundary states. The reference pressure is fixed across restarts. In a
uniform sealed isothermal squeeze/expansion, `p*h` and therefore `lambda/h`
remain constant; a smaller gap alone does not imply greater rarefaction.

The Mach envelope includes pressure-driven gas speed, not only wall motion.
Each accepted face gradient and gap reconstruct the same selected-wall-law
Couette--Poiseuille velocity profile used for traction and dissipation. Its
two wall values and interior stationary point, if present, give the exact
maximum absolute tangential speed; an off-centre peak cannot be hidden by a
mean-flow or mid-gap check. The receipt reports `maximum_gas_mach_number`.
Exceeding either pressure-dependent applicability envelope refuses the step
before a checkpoint or receipt is published; no wall-law or compressible-inertia
correction is silently substituted.

#### Bounded tridiagonal pressure Newton solve

Both strip and radial geometry use the exact tridiagonal Jacobian of the same
implicit gas-mass residual. With fixed face gap, the compiled no-slip/Maxwell
mass conductance has the affine form `K(p)=a*p+c`. A face of measure `W` and
resistance distance `d` therefore has pressure derivatives
`dF/dp_left = W*(K(p_left)/d + h*U_mean/(2*R*T))` and
`dF/dp_right = W*(-K(p_right)/d + h*U_mean/(2*R*T))`.
Storage, prescribed pressure ends, ideal pressure-pinned cells, finite capillary
conductance and the finite reservoir's reduced Schur derivative all enter this
same system. The numerical change does not alter the transport law, geometry,
time discretization, or physical applicability bounds.

Factorization uses the existing native
`fs_la::tridiagonal::TridiagonalWorkspace` with adjacent partial pivoting;
gas-film code adds no parallel linear-algebra implementation. Pressure, trial,
flux, residual, Jacobian and linear-solve storage is allocated once per step
and reused through all corrections and line-search trials. Each assembly and
linear solve has O(N) work/storage. No mesh-independent iteration count or
wall-clock speedup is claimed for arbitrary admitted inputs.

`maximum_iterations` bounds Newton corrections, and `relaxation` remains a
maximum step fraction in `(0,1]`. Each correction is limited to 99% of the
distance to zero for decreasing pressures, then tested with at most 32
successive halvings. The accepted residual either meets the original caller
tolerance or satisfies a strict decrease with Armijo coefficient `1e-4`.
All trial and final checks use the same cell mass residual, including the
finite reservoir storage defect. No pressure floor, artificial diagonal,
tolerance relaxation or unconverged checkpoint is introduced. Linear-solve,
line-search and exhausted-iteration failures are typed refusals. An initially
converged state uses zero corrections. The existing API remains synchronous
without a cancellation callback.

A focused stiff-drainage regression uses a 1 s step on 32-, 128- and 256-cell
annuli, a fixed 12-correction cap and the unchanged mass tolerance. The starting
pressure is uniform; a prescribed smooth gap-motion profile makes the analytic
logarithmic pressure solution an exact backward-Euler endpoint. Tests compare
pressure and signed mass flow with that independent solution and refine its
integrated dissipation. This checks convergence of a stiff numerical step,
not transient accuracy from using a large timestep or the physical validity
of unresolved startup flow. Existing squeeze, Maxwell, moving-wall, vent,
finite-reservoir, admission and restart regressions remain required.

#### Explicit first-order Maxwell slip

`SlipPolicy::MaxwellFirstOrder` admits the classical isothermal Maxwell boundary
condition on two smooth walls with the **same declared** tangential momentum
accommodation coefficient `sigma` in `(0,1]`. Its source identity and coefficient
are bound into restart fingerprints. Existing no-slip fingerprints and wall-law
behavior are retained; an unspecified `RarefiedSlipRequested` still refuses.
No accommodation value is inferred or calibrated. The local slip length is
`b(p) = ((2-sigma)/sigma) * lambda_ref * p_ref/p`, using the same pressure-scaled
mean free path as the Knudsen checks. This explicit policy always refuses local
`Kn > 0.1`, even when the caller supplies a looser cap; a stricter caller cap
remains effective. It models slip flow, not transition or molecular flow.

The closure follows the Maxwell/Navier law described by
[Lockerby (2025), §2](https://doi.org/10.1017/flo.2024.34).
The bounded Knudsen regime is supported by the opposing-wall DSMC study of
[Gavasane et al. (2011)](https://doi.org/10.1260/1759-3093.2.2-3.109).
These references motivate the model and its scope; they do not certify these
caller-declared surfaces or establish a kinetic-theory error bound. There is no
temperature jump, thermal creep, Knudsen-layer correction, or fitted effective
viscosity in this isothermal first-order policy.

Writing `g = dp/dx`, `Delta_U = U_upper-U_lower`, and
`U_mean = (U_upper+U_lower)/2`, the shared wall law gives

```
q = h*U_mean - h^3/(12*mu) * (1+6*b/h) * g
tau_upper =  mu*Delta_U/(h+2*b) + h*g/2
tau_lower = -mu*Delta_U/(h+2*b) + h*g/2
```

The Reynolds mass flux is `(p/(R*T))*q`. The same slip-enhanced conductance enters
the bounded pressure iteration and closed-boundary backpressure reconstruction.
Since `b*p` is constant, the extra mass conductance is pressure independent and
is evaluated without a per-face allocation. The wall-law coefficient is compiled
once per step; the neighbor stencil uses fixed stack storage.

The gas-side wall velocities and interior extremum come from that same quadratic
profile, so the Mach check includes pressure-driven slip and the shifted peak
under unequal wall velocities. Bulk dissipation is
`mu*h*(Delta_U/(h+2*b))^2 + h^3*g^2/(12*mu)`. The additional nonnegative interface
dissipation is `2*mu*b*(Delta_U/(h+2*b))^2 + b*h^2*g^2/(2*mu)`.
`viscous_heat_w_per_m` includes both terms; the new
`wall_slip_dissipation_w_per_m` exposes the interface part separately and is
exactly zero for no-slip. Consequently tangential wall work plus pressure-force
work still equals total dissipation. This is mechanical accounting; it neither
assigns the interface heat to a wall/gas temperature equation nor solves one.

Focused tests cover analytic slip Couette traction and bulk/interface heating,
the exact compressible Poiseuille profile with `(p+6*b*p/h)^2` linear in `x`,
convergent total dissipation, off-centre gas-speed limits in both flow directions,
and moving wedges with sealed/open/vented boundaries and restart. Invalid
accommodation, out-of-regime pressure states, and changes of wall law refuse.

#### Axisymmetric disks and annuli

`solve_isothermal_axisymmetric_gas_film_1d(input, inner_radius_m, checkpoint)`
extends the same pressure iteration, ideal-gas storage, no-slip/Maxwell transport,
pressure-dependent applicability and restart machinery to a fixed radial domain.
`input.grid.length_m` is the radial span. A zero inner radius gives a disk;
a positive inner radius gives an annulus. The frame identity now declares the
symmetry axis and radial direction for this entry point. No second solver or
independent radial wall law is introduced.

The no-slip continuum equation is
`d(p*h)/dt = (1/r)*d[r*h^3*p*(dp/dr)/(12*mu)]/dr`, as in
[the isothermal disk model, equation (2.8)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/role-of-fluidstructure-coupling-in-the-generation-of-an-attractive-squeezefilm-force/22F5E6D5F77C61657F102EC733893EDE).
The code uses exact cell areas
`A_i = 2*pi*r_i*dr = pi*(r_outer^2-r_inner^2)` and circular face measures
`W_f = 2*pi*r_f`. Its implicit mass equation is
`(p_new*h_new-p_old*h_old)/(R*T*dt) + (M_right-M_left)/A_i = 0`.
The shared local face mass flux is multiplied by `W_f`; load and normal work
use `A_i`. Static contact exclusion still removes an outer suffix, and the
pressure/flux boundary lies at the last active circular face.

Each radial face uses resistance distance
`d_f = r_f*ln(r_right/r_left)` between its neighboring pressure samples (a
boundary radius and cell center at an open end). The logarithm is evaluated
with `ln_1p` and the annular area in factored form to avoid subtracting nearly
equal radii. Metrics are compiled once per step; pressure sweeps perform no
radial logarithms. This preserves the exact constant-gap steady radial
potential: with `B = b*p` and `s = 6*B/h`, `(p+s)^2` is affine in `ln(r)` and

```
M = pi*h^3*((p_inner+s)^2-(p_outer+s)^2)
    / (12*mu*R*T*ln(r_outer/r_inner)).
```

No-slip has `s=0`. The pressure-driven heat and pressure-force work use the
same face gradient and dual measure `W_f*d_f`; their discrete equality is
preserved. This dual measure is an energy quadrature, distinct from exact
cell areas. Face pressure/gap and cell radial shear are reconstructed averages;
the Mach gate bounds the reconstructed transverse velocity profiles, not an
unresolved continuum field between radial samples. Integrated dissipation
converges to `M*R*T*ln(p_inner/p_outer)` in the steady constant-gap limit.

The disk center has zero circular measure, regular zero radial flux, and no
evaluation of `ln(0)`. A pressure bath at `r=0` is refused. A drained disk uses

```rust
input.boundary = GasFilmBoundaryTopology::MixedEnds {
    left_absolute_pressure_pa: None,
    right_absolute_pressure_pa: Some(ambient_absolute_pressure_pa),
};
let step = solve_isothermal_axisymmetric_gas_film_1d(&input, 0.0, checkpoint)?;
```

Normal uniform or per-cell squeeze/expansion is supported. An ideal `Vented`
cell represents pressure pinning and total exchange over that retained
positive-area annular cell (the central cell is a disk), not a point port at
the symmetry axis. Per-width `CapillaryVented`/`CapillaryReservoir` connections
are refused in this radial API. Instead, `AxisymmetricCapillaryVented` and
`AxisymmetricCapillaryReservoir` take an explicit nonzero integer **total tube
count** and, for a finite reservoir, a **total volume [m³]**. The Cartesian entry
point refuses these total-unit variants. No arbitrary width or circumference
converts one input basis into the other. Nonzero in-plane wall speed is also refused:
radial wall motion, swirl and centrifugal pressure require additional physics.

The total capillary flow is `M = N*pi*r_tube^4*(p_cell²-p_reservoir²)/(16*mu*L*R*T)`.
The existing no-slip capillary owner, positive-root reservoir elimination and
Newton derivative are shared without an extra coupling iteration. Total `M`
enters the selected cell's residual as `M/A_cell`; it is never multiplied by
that cell's circumference. The reservoir uses `B=V_total/(R*T*dt)`, and its
storage defect is also divided by that cell's annular area for the existing
local mass-residual gate. Tube Reynolds/Mach numbers use the **single-tube**
flow `M/N`. Both film ends remain sealed, so a finite reservoir's exchange is
internal to the combined mass account.

The radial capillary/reservoir receipts expose total work/enthalpy/required
isothermal heat [W], mass [kg] and mass rate [kg/s]. The shared checkpoint stores
`capillary_reservoir_mass_kg` only; its Cartesian per-width field stays absent.
The inverse holds for Cartesian checkpoints. Tube count, reservoir volume and
initialization, attachment cell, and radial geometry bind restart compatibility;
missing or wrongly based retained mass refuses. Bath pressure can change per
step. As for an ideal radial vent, finite-tube exchange is **homogenized around
the selected positive-area annular cell**, including a central disk cell.
Individual holes, azimuthal pressure variations and entrance/junction regions
are unresolved; the total tube count does not claim to resolve discrete ports.
All existing capillary creeping-flow/Kn/Mach and quasi-steady/isothermal limits
continue to apply, independently of the film's Maxwell wall policy.

`AxisymmetricGasFilmReceipt` has distinct total units: mass [kg], mass rate
[kg/s], load [N], and power/enthalpy transport [W]. Its local radial traction
is gas-on-wall, equal on the two stationary in-plane walls; wall-on-gas has
the opposite sign. This radial scalar field does not imply a net Cartesian
radial force, which cancels on azimuthal integration. The pressure field is
the upper-wall normal traction, and total wall work is the reported normal
gap work. Geometry and inner radius bind both checkpoint fingerprints:
changing radius or replaying a strip checkpoint as a disk/annulus refuses.

Focused G0/G1/G3 tests check sealed disk/annulus squeeze mass, actual pressure
load and normal work, contact-excluded areas and radius-squared scaling;
analytic logarithmic annulus flow in both directions with no-slip/Maxwell
walls and heat refinement; and three-rung radial refinement of a manufactured
fully compressible open-rim disk squeeze/expansion pressure and load. The disk
reference comes from differentiating the continuum radial Reynolds equation,
not from the discrete flux. Mixed Cartesian ends additionally retain sealed
wall backpressure and mechanical-work closure. Existing strip, capillary and
reservoir regressions remain applicable to the shared engine.
Additional radial connection tests check total tube-count/radius-fourth-power
flow, each annulus's mass balance, independent closed-volume equilibrium,
three-volume ODE time refinement, small-capacity stiffness and the large-bath
limit, squeeze/expansion transfer reversal, and unit/restart refusals. Tiny
reservoir fixtures test numerical conditioning, not continuum validity at those
physical volumes.

These are axisymmetric lubrication results, with the existing isothermal,
thin-film, small-slope and viscous-dominance assumptions. They do not add
resolved radial inertia/acoustics, edge entrance flow, azimuthal variation,
moving radial boundaries, structural/contact coupling or a thermal equation.

#### Finite-conductance capillary vents

`GasFilmBoundaryTopology::CapillaryVented` connects one retained film cell to
an isothermal pressure bath through declared parallel circular capillaries.
`CapillaryVent` supplies radius `r`, length `L` and tube count per out-of-plane
width `N`. Unlike the existing ideal `Vented` boundary, it does not pin that
cell's pressure or replace its mass-storage equation. The signed outward flow
is `m_dot = N*pi*r^4*(p_cell^2-p_bath^2)/(16*mu*L*R*T)`, the ideal-gas
long-capillary limit in [Berg (2005), equation (5)](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=915339).
The geometry/gas coefficient is compiled once per step. Its exact positive
pressure derivative enters the existing implicit pressure iteration, and the
same signed flow enters the accepted local/global mass balances and boundary
enthalpy receipt. There is no separate explicit discharge timestep restriction.

The capillary walls are stationary and **no-slip independently of the film's
Maxwell accommodation**. The restrictive envelope requires `r/L <= 0.01`,
diameter-based mean-flow `Re <= 1`, radius-based `Kn <= 0.01`, and centreline
`Mach <= 0.1`; stricter caller Kn/Mach caps remain effective. The initial and
accepted tube flows are checked. Mean free path and peak parabolic speed use
the lower of the two absolute endpoint pressures, including reverse flow.
The tube's radius-based Kn is separate from the film's gap-based Kn. These
gates do not certify an error bound on neglected entrance/exit corrections.

The optional `capillary_vent` receipt reports these tube descriptors and its
nonnegative pressure-gradient work/viscous dissipation,
`m_dot*R*T*ln(p_cell/p_bath)`. This work is outside the film: it is neither
added to film viscous heat nor counted again as boundary enthalpy. No wall/gas
thermal partition is supplied. Bath pressure may change between steps while
connection cell and capillary geometry remain bound to restart identity.

This is a quasi-steady fully developed creeping-flow connection to an infinite
bath, not an orifice, a finite gas inventory, a stored tube-pressure field or
a choked-flow law. The cell-local branch does not resolve junction flow or
finite-aperture surface traction. The caller must justify isothermal behavior and forcing
slow compared with the tube's viscous/thermal adjustment; timestep size alone
does not establish those physical assumptions. The cited paper's slip, inertia,
nonideal-gas and thermal corrections and metrological accuracy are not claimed.
Focused regressions cover signed discharge/ingress, the radius-fourth-power
law and work accounting, pressure relaxation against independent two-cell ODE
integration, time refinement, squeeze pressure between sealed and ideal-vent
limits, separate tube applicability, and restart/reversed-bath behavior.

#### Finite reservoir pressure feedback

`GasFilmBoundaryTopology::CapillaryReservoir` replaces the capillary's infinite
bath with a rigid finite volume of the same isothermal ideal gas.
`GasFilmReservoir` declares volume per width and initial absolute pressure;
the accepted checkpoint retains reservoir **mass**. Filling raises its pressure
and reduces discharge, and draining lowers its pressure. Both film endpoints
remain sealed, so capillary exchange is internal to the combined system.

The existing pressure iteration includes this feedback at the same backward-
Euler time level. With `G=N*pi*r^4/(16*mu*L*R*T)` and `B=V/(R*T*dt)`, the
reservoir equation is `G*p_r^2+B*p_r = G*p_f^2+B*p_r_old`. Its positive root
is eliminated analytically from each film update. The reduced derivative is
`d(m_dot)/d(p_f)=2*G*p_f*B/(B+2*G*p_r)`, not the infinite-bath derivative.
Scaled coefficients and a series-conductance expression avoid cancellation
for small/large capacities. There is no extra nonlinear coupling loop,
per-iteration allocation, explicit discharge stability limit, or pressure floor.
The accepted mass receives exactly the same signed transfer used by the film;
the represented reservoir storage residual, divided by the connection cell
width, also participates in the existing mass-residual acceptance budget.

`finite_reservoir` reports accepted pressure/mass, storage, receiving enthalpy,
and combined film-plus-reservoir mass closure. Its pressure uses the original
`p=m*R*T/V` closure. The receiving enthalpy is the film's outward vent enthalpy,
with opposite outward orientation. At fixed temperature, `e=h-R*T`, so a rigid
reservoir requires heat `Q_dot=-m_dot*R*T`: filling rejects heat and draining
requires heat. This explicit thermal requirement is separate from capillary
viscous dissipation and supplies no invented thermal conductance. It does not
turn the isothermal film into a total-energy or temperature-evolution solver.

Reservoir volume/initialization, connection geometry, and gas identity are
bound to restart compatibility. Missing/nonpositive/nonfinite retained mass
refuses; no trial reservoir state is published after failed pressure admission.
The no-slip tube's initial and accepted Kn/Mach/Re checks use the actual
evolving reservoir pressure. The finite volume is spatially uniform and has no
external feed, moving wall, changing composition or resolved internal acoustics.
Its well-mixed/isothermal approximation and negligible capillary storage must
be justified on the physical forcing timescale; a small numerical step does
not supply that justification. Very-small-capacity regressions exercise
algebraic conditioning, not continuum or acoustic validity at those volumes.

Focused regressions check closed-mass equilibration from either pressure order,
independent ideal-gas equilibrium pressure, enthalpy and required heat,
three-storage ODE time refinement, small-capacity stiffness and the large-bath
limit, squeeze/expansion flow reversal, restart identity and failed-step recovery.

The implementation is bit-deterministic on one ISA by fixed cell/face order and
bounded pivoted Newton corrections/backtracking. It has no RNG, no `Cx` cancellation hook, no
entropy model, and no thermochemical or material-data admission; the
retained source identities and `GasFilmInputAuthority` state those limits rather
than minting material/EOS authority.

`tests/gas_film.rs` supplies G0/G1/G3 coverage for uniform sealed equilibrium,
analytic planar Couette shear/heat, uniform squeeze mass preservation,
open Poiseuille direction/closure and exact steady compressible flux/wall load,
convergent Poiseuille dissipation, sealed co-moving-wall backpressure,
pressure flow doing work on translating walls, opposing-gradient heat,
vent closure, spatial/time refinement of analytic limits, reversal and
pressure-density scaling, static contact exclusion,
rarefaction/negative/slip/roughness refusal,
pressure-scaled cell/reservoir/vent rarefaction, sealed compression/expansion
Knudsen invariance, and pressure-driven/off-centre-profile Mach limits,
iteration-budget refusal, and deterministic checkpoint replay.  These are
software/numerical checks of declared equations only.

No gas-film result establishes thermochemical validity, lubrication-regime
selection, roughness or accommodation calibration, continuum validity outside the stated
gates, liquid/mixed/EHL/cavitation behavior, moving-domain topology changes,
general 2-D/3-D geometry, energy/entropy closure, coupled contact, or any
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

Each connected fluid domain requires at least one velocity edge.
Duplicate, missing, interior and out-of-range ownership is refused before
callbacks run. Boundary quadrature samples are retained once for consistent
strong moments and SIP. Non-finite loads and assembly, incompatible net flux
for all-prescribed velocity, singular factorization and failed explicit
original or row-equilibrated residual admission return StokesBoundaryError
without a solution. An invalid viscosity sample retains the offending cell,
point and value in InvalidViscosity, stops subsequent viscosity queries and
discards the private assembly before any solve or publication. No hidden
grounding is added for pure-traction rigid modes.

All-velocity components retain one pressure pin each and independently report
MeanZero pressure. A component with a total-traction edge retains every
continuity equation, has no pressure pin and reports TractionDatum without
subtracting its mean. A mix of these component kinds reports Componentwise;
StokesAssembly::pressure_reference_at and StokesSolution's per-cell reference
metadata disclose which datum applies, also through the material wrapper.
Independent mean-zero components have no known relative pressure offset.
Boundary flux compatibility and Brinkman drag anchors are likewise local:
opposite leaks or a physical anchor in another component cannot cancel a
sealed component's incompatibility or rigid modes. The pressure and traction are kinematic
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

The disconnected-domain regressions in the same test file cover independent
affine pressure-gradient equilibria, dense/sparse agreement, a pressure-datum
shift confined to one open chamber, component-local flux refusal, and the
original discrete adjoint against independently re-solved pressure-difference
sensitivities. The existing Brinkman test also checks independently anchored
porous domains and refusal when only the other domain has drag.

The frozen-field assembly above covers steady linear, constant-density Stokes
with frozen positive spatially varying Newtonian viscosity on the existing
straight-sided 2D triangulation. It does not add nonlinear convection/backflow,
steady pure-traction rigid-mode quotients or continuum error bounds. Existing
synchronous cancellation and dense/iterative solver limitations remain in force.

## Transient mixed-Stokes pressure and total-traction boundaries

`FluxSystem::assemble_stokes_mixed_bdf1` and its
`assemble_stokes_mixed_bdf1_variable_viscosity` field entry assemble
`(u_new-u_previous)/dt - div(2 nu D(u_new)) + grad(p_new) = force_new`, with
`div(u_new)=0`. They add the original BDF1 mass operator and immutable previous
velocity load to the same checked symmetric-stress/mixed-boundary assembly.
All force, velocity and total-traction callbacks describe the new time level;
capturing its time supplies pressure pulses, ramps and changing wall motion.
The field entry admits frozen positive endpoint `nu(cell, point)`, retaining
the existing one-sided harmonic SIP stresses. The constant entry uses
`params.nu`. There is no convection or implicit constitutive update.

`bdf1_stokes_mixed_step` assembles and applies the existing dense LU/sparse
FGMRES selection, requiring independently recomputed original and row-scaled
residuals. A finite positive timestep and finite positive reciprocal, together
with exactly `n_u` finite history moments, are required before callbacks run.
Failed admission, assembly or convergence never modifies caller history or
publishes a partial field. The matrix is also available for an individual
step's existing transpose solve; no time-history adjoint is supplied.

The physical positive mass term removes rigid-velocity nullspaces at finite
dt. Therefore entirely traction-driven components are admitted by these
transient APIs without adding an artificial velocity pin or resistance;
their steady counterparts retain the existing unanchored-domain refusal.
Pressure datums still follow boundary ownership separately on every connected
component: total-traction components retain every continuity equation and the
imposed datum, while sealed components enforce their endpoint flux compatibility
and report their own mean-zero pressure. A datum in another chamber cannot
absorb an incompatible sealed-component flux.

`MaterialFluxSystem::bdf1_stokes_mixed_step` supplies the same linear law with
physical force density [N/m³] and `MaterialStokesBoundary` total traction [Pa].
The immutable material snapshot determines `mu/rho`; that same density scales
inertia, force, traction, reported pressure [Pa] and kinetic energy [J/m].
`MaterialStokesSolution` retains that snapshot, actual per-cell pressure datums
and the velocity moments for the next step. Both time levels require the same
reference density. Mechanical traction pressure is not an EOS absolute state.

`tests/stokes_boundary/transient.rs` exercises variable-step pressure pulses
and free coasting, the exact backward-Euler work/kinetic balance, dense/sparse
agreement, an independently integrated pressure ramp with first-order temporal
error, disconnected sealed/open pressure datums, original mass/stress operator
identity, an exact variable-viscosity affine unsteady solution, and invalid
history/timestep and incompatible-flux refusals. `tests/material_flow.rs`
checks actual shared-gas density in pressure-driven acceleration and physical
pressure, dynamic-viscosity traction, inertia and energy in an affine shear.
The homogeneous entries remain fixed-mesh, constant-density linear unsteady
Stokes, without nonlinear open-boundary/backflow conditions, thermal evolution,
variable-step/BDF3+ integration, continuum error bounds or new cancellation/throughput claims.

### Bounded conservative density subcycling

`ns::transport::density::subcycle` advances one caller-declared macro interval
through the original positive donor-cell Euler owner, with an explicit maximum
substep count and outgoing-volume Courant cap strictly between zero and one.
The prepared BDM volume advector, original mesh, depth, fixed cell volumes and
exterior donor states remain frozen throughout that interval. Equal dyadic
subintervals are admitted before updates; a duration that cannot split exactly
or work beyond the stated cap is refused. This decouples density CFL admission
from a caller's macro cadence; it does not upgrade frozen-advector temporal
accuracy or add compressibility, geometric motion, sources or interface repair.

Every density substep reuses the original mass and ALL N specific-inventory
updates, preserving each actual directional carried transfer. The complete
interval report audits represented initial/final inventories and cumulative
physical boundary amounts, retains macro Courant separately from maximum
substep Courant and count, and scales signed carried residuals with integrated
absolute boundary exchange even when the carried quantity changes sign.
Specific quantities are never inferred by multiplying averaged primitives.

A scale-normalized compensated mean of PRE-substep donor densities avoids
raw density-sum overflow and premature underflow from dividing each tiny
density by the substep count. Since each frozen directional volume rate and
the original compatible mass-flux reconstruction are linear in donor density,
`prepare_mass_flux` reconstructs one time-averaged H(div) mass flux whose
integrated directional rates agree with the accumulated substep mass amounts
and whose cell divergence gives the macro density increment, up to represented
state, summation and reconstruction roundoff. The original quadratic face
moments and nonnegative donor signs are retained. These averages are flux
integration weights, not accepted material/thermodynamic states. Changing
advection, inflow or geometry requires a new plan/interval.

The reusable workspace is bound to its original plan. Its byte cap admits the
original Euler scratch, two stages, compensated means and cumulative face
amounts together. Warm advances allocate no heap storage. Refusal/cancellation
never publishes partial output and invalidates accepted average/transfer views;
the workspace can be retried. Reconstructing a later mass-flux operator has its
own explicit storage allowance. The original single-step API and arithmetic
are unchanged.

G1/G3 regressions in `tests/density_transport/subcycle.rs` compare closed
counterflow and open throughflow against independent rational matrix
recurrences for density AND extensive carried inventory. They check cumulative
directional/boundary accounts, one reconstructed flux's macro divergence,
constant-specific preservation, exact one-substep parity, large/subnormal
positive density means, bounded work/storage and final-poll cancellation/retry.
Native compilation and execution remain required validation gates; these
regressions alone make no convergence or performance certification claim.


### Conservative density and native BDM momentum

`FluxSystem::assemble_density_momentum_bdf1` and its
`density_momentum_bdf1_step` convenience entry couple the positive P0 density
owner to the original affine BDM1 velocity and P0 pressure system. A prescribed,
explicitly lagged divergence-free advector transports density and any N passive
mass-specific quantities. Momentum then solves
`M(rho_new) u_new/dt + C(j) u_new + A(mu) u_new + G p = M(rho_old) u_old/dt + f`
with the original continuity equations, mixed stress/velocity/slip boundaries,
pressure datums, linear residual gates and retained support reactions.

The H(div) mass-flux reconstruction above uses the SAME original directional
donor densities and volume faces. Its constant cell divergence agrees with the
P0 density increment up to reported represented-state roundoff; its quadratic
face moments also match the square of every affine velocity. Full consistent
old and new density metrics remain distinct. Old momentum loads are integrated
directly, avoiding subtraction of a new-density history correction at large
contrast. There is no cell-average velocity projection or lost subcell mode.

Inputs are physical dynamic viscosity [Pa s], body-force density [N/m^3],
reference density [kg/m^3], depth [m], and the existing
`MaterialNavierStokesBoundary` owners. Force and viscosity callbacks receive the
transported candidate cell and must describe frozen endpoint laws. One fixed
positive reference normalizes all momentum, stress, traction, slip, pressure
and inertia terms. The accepted result retains physical cells, mechanical
`pressure_pa`, consistent kinetic energy [J], linear momentum [kg m/s], and the
original donor-cell inventory report. Its `flow` view retains native normalized
pressure/stress units. Boundary transport keeps the actual volume rate separate
from density-weighted momentum and kinetic rates; a volume flux never gains an
extra density multiplier.

For homogeneous closed velocity data and the shared discrete mass balance,
`K_new-K_old + 0.5*||u_new-u_old||_(rho_old)^2 + dt*D_upwind
 + dt*a(u_new,u_new) = dt*(f,u_new)` up to quadrature/solve/roundoff defects.
`D_upwind` uses both directional mass streams and is nonnegative. The original
viscous/SIP form a, boundary kinetic transport, boundary upwinding and explicitly
chosen backflow resistance keep their separate meanings. The result reports
time and interior-upwind losses; they are numerical dissipation, not physical
heating. The N carried slots remain passive (e.g. enthalpy/species): momentum
forces act on the native BDM velocity, not on a velocity-like passive slot.

This is a first-order fixed-mesh IMEX discretization: mass flux is explicitly
lagged while momentum, viscosity and pressure are implicit. Density positivity
retains the outgoing-volume CFL restriction. It does not implement a fully
implicit coupled density/velocity iteration, nonlinear constitutive feedback,
EOS/compressibility, phase-volume sources, moving interfaces with sharp
reconstruction, heat diffusion or pressure-work closure. The exposed native
matrix permits its existing fixed-operator transpose solve; it is not a full
transport-step Jacobian. Per-owner transport storage is bounded; the original
native linear assembly/solver keeps its own storage and synchronous cancellation
limitations. Cancellation is checked around that solve and during preparation
and observations. Caller histories never mutate and no partial coupled state
is returned on refusal.

`tests/density_momentum.rs` covers unequal-density free-stream transport through
physical reservoirs, large density contrast and reference-normalization
invariance, actual boundary momentum/kinetic exchange, closed counterflow
kinetic and wall-momentum balances, first-order coupled temporal convergence
against an independent closed-form solution over three timestep refinements,
physical density-weighted acceleration and traction pressure, recovery of the
original uniform-density Oseen matrix/load, and cancellation, CFL, invalid-source
and unrepresentable-history refusals.

### Implicit SGS stress with conservatively transported density

Under `experimental-sgs`,
`FluxSystem::density_momentum_bdf1_sgs_step` closes molecular rheology and
Smagorinsky, WALE or Vreman stress at the same native BDM1 velocity endpoint.
The molecular source receives the original cell/centroid, actual candidate
shear rate and this step's transported density and specific inventories. It
returns one positive dynamic viscosity per cell; the same value is used in
volume stress, harmonic SIP traces and physical boundary reactions. The shared
full-gradient kernels supply kinematic eddy viscosity, converted with the
actual transported local density: `mu_eff = mu_molecular + rho_new * nu_sgs`.
The arbitrary conditioning reference never replaces that density.

Density and carried inventories advance exactly once. Their original mass
flux, distinct old/new inertia, physical force and admitted boundary samples
remain frozen through the original bounded Picard solve. Only constitutive
stress changes with each candidate. Acceptance requires the unrelaxed velocity
correction and both original and row-scaled residuals reassembled at the
returned endpoint. Existing donor/limited density reconstruction, native
pressure datums, conservative implicit momentum transport and support
reactions retain their original owners. The fixed-viscosity entry and density
subcycle/SSPRK2 paths retain their original arithmetic.

The result keeps the exact accepted molecular source samples and modeled
coefficients, without a reporting-time source query. Cell-interior molecular
heating `mu_molecular * gamma^2` and modeled resolved-to-subgrid transfer
`rho_new * nu_sgs * gamma^2` are separate [W/m^3] fields and physical volume
integrals [W], including the authored depth. Neither is conflated with native
temporal, upwind or SIP losses, nor silently deposited in a passive heat slot.

This remains an Estimated planar turbulent-stress screen, with explicit
isotropic filter widths and no DNS/three-dimensional LES validation,
backscatter, SGS energy equation, turbulent density/species diffusion or wall
model. Momentum is BDF1 with an explicitly frozen mass advector, not a fully
implicit density/velocity iteration. Source laws are cellwise; no subcell
viscosity reconstruction is inferred. Cancellation surrounds the original
synchronous solve and polls source/observation loops. Refusal publishes no
partial state and never modifies either caller history.

`tests/density_momentum/sgs.rs` exercises all three closures against changing
unequal-density counterflow, the native kinetic and wall-momentum identities,
source/force/boundary sampling and reference-density invariance, an affine
physical-traction solution with equal strain but different rotation, the exact
zero-model limit, nonlinear refusal and cancellation.

### Bounded linear density with compatible native momentum

`FluxSystem::with_density_reconstruction(DensityReconstruction::BoundedLinear)`
selects an opt-in spatial density reconstruction for the existing Euler,
subcycled and SSPRK2 transport owners and all native conservative momentum
entry points that prepare transport from that system. `DonorCell` remains the
default, including its original transfer arithmetic and allocation-free warm
path. No alternate momentum solver or material inventory is introduced.

The original mesh prepares weighted least-squares density slopes once. Each
step limits the three reconstructed vertex densities against the minimum and
maximum densities of their incident cells and actual incoming boundary states.
Those vertex patches preserve affine interior profiles; face-neighbor bounds
alone could unnecessarily limit every affine triangle. The limiter contracts
a bound-touching slope towards its positive cell mean rather than clipping an
accepted density or repairing mass. Singular one-sided stencils and
unrepresentable slopes retain the original constant donor. No exterior ghost
density is invented at an outflow or impermeable boundary.

Each sign-separated BDM face stream uses the exact volume-flux-weighted
centroid of its affine normal velocity. Evaluating the limited affine density
there gives that stream's integrated reconstructed mass flux. Both directions
survive a reversing face. Specific inventories keep their original donor values
and use those SAME mass transfers, so constant specific states and convex
specific bounds are preserved. This raises density's smooth spatial accuracy;
it does not raise the spatial order of specific enthalpy/species transport or
the temporal order of native BDF1 momentum. The optional carried-inventory
reconstruction below addresses the spatial heat/species limitation.

For a triangle, the mean of the three affine vertex densities is the stored
cell density. Every nonnegative vertex value is therefore at most three cell
means. A sufficient positivity restriction is the original volume Courant
allowance divided by three, applied automatically by `maximum_timestep_s`,
Euler admission and dyadic stage selection. Reports retain ACTUAL volume
Courant numbers. Nonnegative retained vertex weights and incoming streams
preserve positive density and global density bounds, subject to the existing
volume-divergence and floating-point allowances. This is a conservative
sufficient bound, not a claim of the largest stable timestep.

The exact directional mass is shared with the original compatible BDM2 mass
flux. Each effective directional density multiplies the original affine normal
velocity, retaining its quadratic test moments and constant cell divergence.
This adopted momentum trace need not equal the pointwise product of the two
affine fields; no higher spatial-order momentum claim is made. SSPRK2 and
subcycling accumulate actual limited directional donors with their accepted
time weights. Reconstructing an averaged cell state would be incorrect because
the limiter is nonlinear. The native density-weighted inertia, pressure,
momentum flux and kinetic-energy identity consequently remain consistent with
the accepted density transfers.

Prepared stencils, one vertex-bound pass and reusable reconstruction/mean
buffers keep warm work linear in cells, vertices and edges, without new heap
allocations during a step. All added storage is included in existing per-owner
caps. No measured speedup or sharp-interface reconstruction is claimed.
`tests/density_transport/linear.rs` checks smooth translating density error,
affine counterflow integration, positivity and carried-inventory bounds across
a large density jump, and every accepted Euler/RK face transfer against its
momentum flux. `tests/density_momentum/ssprk2.rs` also checks reconstructed-density
free streams and the closed native kinetic-energy identity.

### Bounded linear heat and species on the accepted mass flux

`DensityReconstruction::BoundedLinearInventories` retains the same bounded
density, directional mass transfers, native momentum coupling and sufficient
Courant limit, while reconstructing the carried specific quantities. The
least-squares stencil locates each cell's specific mean at the **mass centroid**
of its reconstructed density. Its weights are rebuilt once per stage and shared
by all components, with variations centered to have zero **mass-weighted**
triangle mean. Affine physical specific fields supplied as mass means are
reproduced where the density reconstruction and inventory limiter retain their
affine slopes and the local stencil has full rank.
Each original sign interval uses its density-weighted flux centroid; this
exactly integrates the selected affine density/specific traces against the
affine normal velocity, including their cubic covariance term.

A single cell limiter bounds all components at the original vertex patches
and bounds the retained donor inventory. For retained mass `r=m-sum(M_out)`
and reconstructed outgoing correction `C=sum(M_out*delta_s)`, it admits only
`s-theta*C/r` within the patch bounds. The updated state is therefore a
positive mixture of bounded retained and incoming states. Before emission,
actual represented face corrections are checked against those inventory
margins; at most eight common slope contractions precede donor-cell fallback
if roundoff cannot resolve a bounded correction. Equal/opposite corrections
conserve each extensive inventory, and common limiting preserves complementary
species sums to roundoff. Constant specifics remain exact. No accepted-state
clipping, species renormalization or change to the mass flux occurs.

Euler, subcycling and SSPRK2 use the same reconstruction and transfer accounts.
Work is linear in component count and native mesh size, with bounded local
roundoff retries and no warm heap allocation; all extra scratch is budgeted.
`tests/density_transport/linear.rs` checks spatial thermal-pulse refinement at
matched time steps, reduced heat diffusion, affine mass-mean reproduction and
independent signed cubic face quadrature, large-density-contrast species bounds/sums, flux compatibility and
transactional cancellation. Limiter activation can reduce local accuracy;
this adds neither sharp interfaces nor diffusion, reactions, pressure work or
thermodynamic closure.

### Conservative SSPRK2 density stages and accepted flux weights

`ns::transport::density::ssprk2` adds bounded second-order TIME integration of
the original fixed-advector spatially discrete density/inventory problem.
Each physical subinterval evaluates the original positive donor-cell Euler
map twice, `U1=E_h(U0)` and `U2=E_h(U1)`, then accepts `0.5*(U0+U2)` in
EXTENSIVE coordinates `(mass, mass*specific[0..N])`. The fixed mesh, depth,
original divergence-free BDM volume advector and exterior donor states remain
unchanged throughout the macro interval. Default P0 donor-cell spatial transport
retains its original first-order mixing; the bounded linear density option is
described above. This owner neither sharpens interfaces nor raises
the temporal order of a caller's native BDF1 momentum step.

Both stages satisfy the original outgoing-volume Courant cap in `(0,1)`.
The accepted extensive convex combination preserves density positivity and
the original donor-specific bounds up to represented-state roundoff. Specific
values are not averaged, and small mass fractions are not formed before
inventory products: representable light-donor inventories survive even when
the corresponding mass fraction underflows. The blend uses compensated actual
physical mass-times-specific products and their FMA remainders, switching to
safe half-sums when unweighted sums overflow. Equal specific values copy
exactly. Original Euler arithmetic and its admission range are unchanged;
this does not claim to admit every state accepted by the isolated blend.

Accepted directional face amounts are half-sums of the two actual stage
transfers, accumulated over physical subintervals. The two stages do NOT count
as two full physical intervals of boundary exchange. Robust pair midpoints
retain two equal subnormal amounts and near-cancelling signed inventories.
Macro reports reuse the original endpoint and physical boundary audit, with
RK-weighted absolute boundary exchange as the signed-inventory residual scale.
They distinguish physical substep count, exactly twice as many Euler stage
evaluations, stage Courant and macro Courant.

The mean of the two PRE-Euler donor densities per physical subinterval has
exactly the same RK weights as the accepted face amounts. Its existing
scale-normalized compensated accumulator avoids raw density-sum overflow or
premature division underflow. One original compatible mass-flux reconstruction
then retains macro density-divergence compatibility and quadratic native
kinetic moments, up to represented-state, summation and reconstruction error.
The existing `mass_flux::momentum` action remains the momentum convection owner;
no separate convection discretization is introduced.

Equal dyadic physical subintervals, representable duration splits and a caller
maximum count bound work. One workspace allowance covers original Euler
scratch, all three state buffers, first-stage transfers, accepted transfer
sums and donor averages. Warm advances allocate no heap storage. Cancellation
or refusal leaves caller output unchanged, invalidates accepted views and
allows retry. A reconstructed mass-flux operator retains its separate explicit
storage allowance. No complete nonlinear transport adjoint is claimed.

G1/G3 regression sources compare counterflow and throughflow against independent
RK matrix polynomials, verify second-order temporal refinement on a FIXED
spatial mesh, check all extensive/directional/boundary accounts and one
compatible macro flux, exercise admitted high-contrast stages and subnormal
face weights, and cover work/storage/final-poll refusal and exact retry. Private
blend tests separately cover an extreme density ratio whose direct mass fraction
underflows. Native Rust execution remains a validation gate; source inspection
does not establish a runtime or performance result.


### Density subcycles with one native macro momentum solve

`ns::boundary::conservative::subcycle` adds
`assemble_density_momentum_subcycled_bdf1` and
`density_momentum_subcycled_bdf1_step`. The explicit
`DensityMomentumSubcycleControls` retains the original physical momentum
controls and adds a maximum density-substep count. Its outgoing-volume Courant
cap applies to each density substep; the accepted report retains the separate
macro Courant, substep Courant, dyadic count and interval.

The bounded density owner advances all actual mass/specific inventories under
one frozen original BDM volume advector, fixed geometry and fixed exterior
donors. Its compensated PRE-step density mean reconstructs ONE compatible
mass flux. The resulting final density, original initial density and complete
macro timestep enter the original consistent native inertia and conservative
momentum assembly. Both public time-integration paths call the same extracted
native assembly tail with unchanged matrix/RHS arithmetic. There is one native
momentum/pressure assembly and one original global solve per convenience-step
call, independent of the admitted density substep count. Force and viscosity
callbacks see the final transported endpoint cells; they are not sampled or
solved at density stages. No averaged specific or thermodynamic state is
invented.

The original pressure normalization, velocity/slip/traction/reservoir/backflow
boundaries, complete affine BDM velocity modes, reaction rows and transpose
surface remain in use. The one time-averaged mass flux has the macro density
increment as its divergence and retains original quadratic face moments up
to represented-state/reconstruction roundoff. Thus free-stream compatibility
and the native conservative kinetic identity use old-density temporal loss,
macro integrated upwind loss and the original viscous/SIP/external work terms.
Physical boundary inventory totals come from actual density substeps;
momentum/kinetic boundary observations use the same averaged native flux.
Numerical kinetic losses remain separate from transported thermal inventories.

Each original owner retains its explicit storage allowance; the density
subcycle workspace admits all stages and compensation buffers together and
reuses them without warm allocations. Cancellation/refusal returns no coupled
candidate and leaves immutable caller state unchanged. The original
synchronous linear solve is checked before/after, not advertised as
interruptible. Native matrix/factor storage remains separately owned. A larger
density-CFL-admitted interval does not improve the first-order frozen-advector
temporal approximation or introduce nonlinear density/velocity iteration,
compressibility, diffusion, phase sources or moving geometry.

G1/G3 regression sources in `tests/density_momentum/subcycle.rs` check exact
one-substep equality of the original native matrix, RHS, solution and physical
observations; high-contrast macro free streams with physical pressure and
boundary momentum/energy balance; and the closed macro kinetic/wall-momentum
identity against an independent one-mode native solution. Endpoint force-call
counts stay one original six-point cell quadrature regardless of substeps.
Bounded work and late cancellation return no partial candidate. Native Rust
execution remains a validation gate; no runtime, convergence-order improvement
or performance certification is inferred from source inspection.


### SSPRK2 density inventories with one original BDF1 momentum solve

`ns::boundary::conservative::ssprk2` adds
`assemble_density_momentum_ssprk2_bdf1` and
`density_momentum_ssprk2_bdf1_step`. Their names and reports distinguish the
second-order fixed-advector density/inventory stages from the original
FIRST-ORDER native BDF1 macro momentum update. `DensityMomentumSspRk2Controls`
retains original physical controls and a maximum physical density-subinterval
count; each subinterval has two original Euler evaluations. The density stage
report records both counts, stage Courant and the complete physical macro
inventory/boundary report.

Density, all passive extensive inventories and accepted physical face amounts
use the same SSPRK2 weights. The one compatible mass-flux reconstruction uses
the corresponding PRE-stage donor-density mean. Original initial density,
full-interval final density, macro timestep and that single flux enter the
existing shared conservative native assembly tail without modification.
There is one native momentum/pressure assembly and one original global solve
per convenience-step call, independent of density-stage count. Endpoint force
and viscosity laws see the final transported cells; no native solve or source
sampling is performed at density stages.

The original affine BDM velocity space, pressure normalization, boundary
policies, reaction rows and native transpose surface remain in use. Macro
mass-divergence compatibility and quadratic kinetic face moments retain the
original free-stream and conservative kinetic balances up to represented
state and reconstruction error. Old-density temporal loss and macro upwind
loss remain separate numerical observations; they are not deposited into
passive thermal inventories. Physical boundary mass/carried amounts retain
the actual RK weights and are not doubled with the stage count. The existing
`mass_flux::momentum` operator is unchanged.

The original per-owner storage allowances, immutable caller history and
no-candidate-on-refusal/cancellation behavior remain in effect. Complete
density-stage scratch is admitted together. Native matrix/factor memory
remains separately owned, and the synchronous linear solve is polled before
and after rather than claimed interruptible. P0 spatial mixing and the
frozen-advector BDF1 momentum coupling remain first order. No second-order
fully coupled flow claim, interface reconstruction, moving geometry, EOS,
phase source, diffusion or complete nonlinear adjoint is added.

G1/G3 regression sources in `tests/density_momentum/ssprk2.rs` check weighted
macro free streams across a large density contrast and multiple physical
reference densities, mechanical pressure and boundary momentum/energy
accounts, closed macro kinetic/wall balances and an independent one-mode
BDF1 solution. Uniform density recovers the existing Euler-subcycled native
macro matrix, RHS, solved field and kinetic observations exactly. Force-call
counts remain one original cell quadrature independent of RK stages. Work,
source and cancellation refusals publish no partial coupled candidate.
Native Rust compilation and execution remain validation gates; no runtime or
performance result is claimed from source inspection.


### Fixed stratified density in native Stokes inertia

`assemble_stokes_mixed_bdf1_variable_density` adds the spatially weighted
consistent mass to the same original mixed-boundary BDF1 solver:
`r(x) (u_new-u_previous)/dt - div(2 nu_reference D(u_new)) + grad(pi) = f`.
The caller chooses one positive constant `rho_reference`. The density callback
returns `r=rho(x)/rho_reference`; viscosity is `mu(x)/rho_reference`, force is
physical force density divided by `rho_reference`, and pressure, total traction
and Navier slip resistance use that same normalization. Dividing pressure or
viscosity by the local density would solve a different equation. The primitive
result retains its existing reference-normalized units and pressure datums.

Positive density is sampled once at each of the original six positive volume
quadrature points per cell; the cell index permits unsmeared one-sided material
jumps. Stored weights supply both the full signed consistent mass and the
immutable history load, including retained support-reaction equations. No mass
lumping or cell-mean substitution is used. BDM1 velocity is affine, so the
degree-four rule integrates density profiles through degree two per cell
exactly; general fields require spatial quadrature/refinement. Unit relative
density retains the original matrix and history arithmetic.

`StokesAssembly::kinetic_energy` observes that same retained metric without
re-querying density and returns energy per depth divided by `rho_reference`.
Multiplying by `rho_reference` gives J/m. It accepts either original velocity
moments or the complete solved field and checks representability. Invalid
density samples refuse with their cell, position and value before constitutive
or force assembly; failed steps never publish a partial state.

This density field is fixed at both time levels: linear perturbations about a
stationary stratification and shear parallel to density layers are supported.
The entry does not advect density, move material interfaces, introduce
nonlinear momentum convection, or model mixing, thermal feedback or gas
compressibility. Its Stokes boundary transport diagnostics remain zero because
convection is absent. `tests/stokes_boundary/transient.rs` adds independent
unequal-density pressure acceleration/coasting and BDF1 work checks, quadratic
stratified-density shear with analytic energy and supported stress/torque,
invariance under the reference-density choice, original unit-density matrix
identity, retained sampling and typed refusal regressions.

### Source-bound spatial Stokes material fields

`MaterialFluxSystem::bind_stokes_materials` prepares one original
`FluidTransport` owner per triangle for the density-weighted mixed-Stokes
entry. Owned sources and heterogeneous `&dyn FluidTransport` sources are
supported. Shape admission precedes source queries; each owner is queried once,
and its exact density/dynamic-viscosity snapshot and original source are retained.
Original material refusals propagate unchanged, with no query after the first
refusal and no partial field publication.

`MaterialStokesField` uses the existing binding's density only as ONE fixed
normalization. Every actual cell supplies `rho/rho_reference` to inertia and
`mu/rho_reference` to one-sided viscous stress. Force density [N/m³], total
traction [Pa] and physical Navier resistance [Pa s/m] share that reference;
pressure [Pa], reactions [N/m, N] and slip power [W/m] recover through the same
reference. The homogeneous reference viscosity is unused. Positive conversions
that overflow or underflow to zero are refused.

Repeated `bdf1_stokes_mixed_step` calls reuse the frozen source snapshots.
Consistent density-weighted mass, immutable history, support reactions and
kinetic energy use the existing native operator and its retained quadrature.
The resulting `MaterialStokesFieldSolution` reports energy [J/m] from that
weighted metric, and bulk heat from the actual cell dynamic viscosities.
Snapshots and original sources are shared into each result without cloning
owners or querying material laws again; physical output remains tied to the
admitted state even if a source can subsequently change.

The piecewise-constant material field remains fixed across time levels on the
original mesh. This supports stationary stratification and flow parallel to
material layers. It does not transport density, mix fluids, move interfaces,
add convection, evolve temperature or infer gas compressibility; changing
density requires a separately conservative state transfer. Existing rheological
paths still require constant reference density rather than silently changing
their inertial metric. Cancellation and spatial-order limitations remain those
of the native mixed Stokes solver.

`tests/material_flow/density.rs` specifies G1/G3 regressions using original
dry-air evaluations in two layers: unequal density and viscosity, continuous
interface shear stress, physical Navier slip, repeated-step momentum and
history, analytic weighted energy, actual-mu heating, physical pressure and
boundary force/moment, and invariance under two reference normalizations.
Focused admission tests retain original material errors, freeze source
snapshots across solves/reporting, and refuse unrepresentable normalization.

### Second-order time integration with frozen density

`assemble_stokes_mixed_bdf2`, its variable-viscosity and variable-density
variants, and `bdf2_stokes_mixed_step` add equal-step BDF2 to the original
linear mixed-Stokes operator. The homogeneous material binding and retained
`MaterialStokesField` expose the same physical-unit step. `previous_velocity`
is u_n; `older_velocity=Some(u_(n-1))` uses
`(3 u_new-4 u_n+u_(n-1))/(2 dt)`. Both supplied history intervals must equal
the requested timestep; these slice APIs do not carry timestamps that could
verify the spacing. `None` executes the original BDF1 startup exactly. For
smooth solutions, that one startup step preserves second-order convergence at
fixed positive final time.

The implementation assembles the existing operator with effective interval
`2 dt/3` and immutable history `u_n+(u_n-u_(n-1))/3`. A constant history is
preserved exactly in this combination. One retained density quadrature supplies
the full consistent mass, both historical contributions, support reactions and
kinetic energy. Source-material fields retain their original frozen density and
dynamic viscosity at all three time levels; pressure, force, traction and slip
use the same fixed reference normalization. Invalid histories and unrepresentable
effective timestep/history refuse before field or boundary callbacks.

For zero forcing, homogeneous boundary data and a dissipative original spatial
form a, equal-step BDF2 satisfies
`G_new-G_n + ||u_new-2u_n+u_(n-1)||_rho^2/4 + dt*a(u_new,u_new)=0`,
where `G_n=(||u_n||_rho^2+||2u_n-u_(n-1)||_rho^2)/4`.
Thus this two-state energy decreases when a is nonnegative; single-state
physical kinetic energy need not decrease at every step. The spatial term
includes the original viscous/SIP form, not just physical bulk heating.
No variable-step, arbitrary-boundary-work, noncoercive-penalty or nonlinear
stability claim is made. Density transport, changing mass, moving interfaces
and Navier–Stokes BDF2 remain outside this entry.

`tests/stokes_boundary/transient/bdf2.rs` specifies exact quadratic-in-time
graded-density shear, common-reference invariance, both history contributions,
physical stress force/moment, exact BDF1 startup, second-order time refinement,
constant-history free stream, and the weighted G-energy identities with and
without viscous/SIP loss. `tests/material_flow/density.rs` carries original
two-layer gas sources through BDF2 and checks physical pressure/energy, source
reuse, both immutable histories and the homogeneous source adapter.

### Local Coriolis coupling on native incompressible momentum

`FluxSystem::with_coriolis_parameter(f)` admits a finite signed planar
`f = 2*Omega_z` [s^-1] for frame-relative velocity. The common BDM momentum
owner adds `+f*J*u` on the left, with `J(u_x,u_y)=(-u_y,u_x)`; positive
f turns an unforced velocity clockwise. Velocity, loads and boundary data
must use the same local frame.

The cell operator integrates `r*f*v dot J*u` using the original physical
vector basis and the same current density ratio `r=rho/rho_reference`
quadrature as inertia, including fixed and conservatively transported density.
BDM edge-normal coefficients are never treated as Cartesian velocity pairs.
Each local pair is accumulated once and mirrored with its negative, with a
zero diagonal before constraints. Existing boundary-reaction recovery and
shared rheology, Newton, material and porous assemblies retain the coupling.
The added cell block uses fixed-size storage; zero f skips it entirely.

The skew source supplies no instantaneous kinetic-energy work or viscous heat.
Existing BDF1/BDF2 time integration is unchanged: for example, uniform BDF1
inertial oscillations retain amplitude `1/sqrt(1+(f*dt)^2)` per step. The
skew spatial operator therefore does not imply energy conservation of the
complete time discretization. This is a constant-f local planar model;
centrifugal, Euler and origin acceleration remain explicit body forces.

The translating and rigid ALE owners evolve WORLD/absolute velocity and
refuse nonzero relative-velocity Coriolis configuration; the checked rigid
history assembly also enforces this boundary. Regression definitions in
`tests/coriolis.rs` cover both signs of inertial oscillation, a density-weighted
analytic bilinear integral, zero source work, a geostrophic pressure/reaction
case, material-density scaling, sparse solution, zero equivalence and invalid
or incompatible inputs.

## Nonlinear pressure-driven flow with explicit backflow behavior

`solve_navier_stokes_mixed`, its frozen spatial-viscosity variant, and
`solve_navier_stokes_mixed_rheology` solve the true steady equations with
convection and symmetric viscous stress. Their initial velocity is a nonlinear
guess only: the original assembly receives zero mass coefficient and no
history. They reuse the mixed Picard endpoint checks, physical boundary
reactions, Navier-slip law, reservoir momentum and explicit backflow policies.
Steady boundary admission retains the original rigid-mode anchoring refusal;
no fictitious timestep or inertia makes an unanchored domain admissible.
Nonlinear controls and finite initial moments are checked before callbacks.

`bdf1_navier_stokes_mixed_step` advances constant-density incompressible
Navier–Stokes with backward Euler and endpoint convection. Its frozen spatial
viscosity counterpart admits `nu(cell, point)` without changing the physical
density. Both reuse the original BDM1/P0 mass, symmetric one-sided SIP stress,
conservative cell advection and upwind interior faces. Reversing affine face
fluxes use the existing split quadrature. The separate Oseen assembly entries
expose this same operator with an explicit frozen advecting velocity; solving
one such linear system is not claimed to solve the nonlinear endpoint.

Every `NavierStokesBoundary` edge owns prescribed velocity, open traction
with an explicit `OpenBoundaryPolicy`, or `OpenReservoir` with literal traction
and a separately supplied incoming velocity. On the original `Open` variant,
`NaturalTraction` imposes the supplied total stress literally and uses the
interior convected velocity on either flow sign. Backflow may bring kinetic energy into that domain. In contrast,
`DissipativeBackflow` imposes, in kinematic units,
`(2 nu D(u)-p I)n = t + min(u.n,0)u`. Its added inflow resistance is a declared
change to the boundary law, not the same physical traction with a hidden
numerical fix. In physical units that resistance is multiplied by the same
source density. It is relative to the specified fixed frame.
The conservative boundary momentum flux becomes `max(u.n,0)u`. For a
divergence-free advecting field, its open-boundary contribution to kinetic
energy loss is `|u.n||u|²/2`, in addition to interior DG jump dissipation and
separate imposed-traction work. Natural traction instead has the signed
`(u.n)|u|²/2` contribution. This is a statement about the declared convective
operator, not an unconditional whole-scheme theorem. For related modified
traction and backflow analysis see Dong and Shen (2015), JCP 291, 254–278,
<https://doi.org/10.1016/j.jcp.2015.03.012>; their pressure-correction algorithm
is not implemented here.

The existing bounded mixed Picard loop is shared with constitutive Stokes.
Boundary samples are admitted once. Transient calls retain the same old mass
history; steady calls have no mass/history term. Acceptance requires the
unrelaxed correction and independently
reassembled original and row-scaled endpoint residuals. Under-relaxation never
changes the timestep or physical loads. Exhausted work or any invalid input
returns no endpoint. Transient inertia anchors entirely open components.
Steady solves require the original physical rigid-mode anchors; closed
components retain independent pressure gauges and flux compatibility.
`NavierStokesStep` retains each open edge's selected policy with the nonlinear
work and residuals. The SI `MaterialFluxSystem` entry binds the same source
density to inertia, convection, `mu/rho`, force, traction, added resistance,
reported pressure and kinetic energy, retaining the material snapshot.

The focused regressions use an exact affine mixed flow with nonzero history,
variable viscosity and live convection, a pressure-datum shift, a
convection-omitted falsifier, independently integrated natural/backflow matrix
work including a reversing boundary face, exact original all-velocity operator
identity, invalid/exhausted nonlinear steps, and two actual gas-source states
with coherent SI stress and energy. They were added but not executed in the
resource-limited implementation session. These original entries use frozen
viscosity. The constitutive extension below closes viscosity at the endpoint;
all of these entries remain fixed-mesh 2-D incompressible flow, without
variable-density storage, thermal or turbulent closure, nonlinear adjoints,
whole-scheme stability, scalable throughput or synchronous cancellation claims.

The steady `MaterialFluxSystem::solve_navier_stokes_mixed` and
`solve_navier_stokes_mixed_rheology` entries preserve the same source density
for convection, viscosity, force, traction, slip, pressure and kinetic energy.
The source-aware steady solve uses the original cell/rate cache and retains
the exact accepted material owners; their dynamic viscosities also supply
physical bulk heating. Source failure or changed density refuses publication,
including when detected only during final endpoint assembly. Temperature and
composition remain prescribed source context. Initial guesses are immutable
and contribute no transient load. Existing material regressions now cover
steady natural/backflow cases, retained responses, SI stress and heat, and
the same endpoint from distinct initial guesses.

### Reservoir momentum through pressure openings

`NavierStokesBoundary::OpenReservoir` supplies total traction and the velocity
carried by incoming fluid. The existing pressure/continuity solve determines
the normal flow rate; the reservoir's normal velocity is carried momentum,
not a strong flow-rate constraint. For retained advector w, the conservative
flux is `(w.n)*u` on outflow and `(w.n)*u_reservoir` on inflow. The original
symmetric stress is the literal supplied traction, with no added resistance.

Ownership is validated before reservoir callbacks. Two original edge-Gauss
samples define a fixed affine reservoir trace, reused by every Picard assembly,
endpoint residual and reaction. Existing sign-split quadrature integrates each
polynomial branch exactly; it does not resample a nonlinear callback at changing
split positions. Invalid samples or exhausted native solve budgets publish no
endpoint and preserve the borrowed velocity history. There is no additional
solver, per-iteration source query, or new synchronous cancellation guarantee.
The separate reservoir sample vector is allocated only when this boundary is present.

`BoundaryReactionOwner::OpenReservoir` identifies the returned exchange.
`outward_momentum` uses the actual admitted incoming trace and retained advector;
mechanical stress stays separate and backflow resistance is zero. Reported
open-boundary stress policies remain NaturalTraction for these edges. Both
Newtonian and implicit rheological material routes retain the original density
for inertia, momentum flux and SI reaction units, and the actual accepted
viscosity/source owners. The translating-domain adapter evaluates reservoir
callbacks at world positions and subtracts its accepted mesh velocity before
the native solve; world flux recovery uses the existing retained mass moments.

The reservoir test modules under `tests/stokes_boundary` and
`tests/material_flow`, plus `tests/ale_gcl/translating.rs`, cover independent
polynomial force/torque integrals on a reversing face, live BDF1 momentum
injection and volume balances, callback/refusal behavior, actual gas-density
swaps, coupled endpoint rheology, and translating-frame covariance. These
tests do not establish a whole-scheme stability guarantee.
This prescribed infinite-reservoir boundary adds no reservoir inventory,
variable density, thermal/species transport or whole-scheme stability claim.

### Coupled rheology and pressure-driven transient inertia

`FluxSystem::bdf1_navier_stokes_mixed_rheology_step` advances
`(u-u_old)/dt + div(u tensor u) - div(2 nu(x,gamma) D(u)) + grad(p) = f`,
where `gamma = sqrt(2 D:D)` [1/s] is evaluated from each cell's actual affine
BDM1 velocity. The caller's law supplies positive finite kinematic viscosity
[m²/s] at the given cell, position and strain. Spatial viscosity gradients
and one-sided interface stresses remain in the original symmetric-stress
operator. There is no separate pressure approximation or Laplacian shortcut.
The actual old velocity initializes both constitutive strain and convection;
the new path does not probe an invented zero-shear state or add a floor.

Convection and viscosity use the SAME candidate at every Picard assembly and
at independent endpoint acceptance. The shared bounded iteration still
requires the unrelaxed correction and original/row-scaled nonlinear residuals.
Accepting a frozen constitutive response or frozen history advection is
insufficient. Boundary samples, component-local pressure datums, continuity,
history and the natural/dissipative backflow laws above retain their existing
owners. The frozen-viscosity path does not allocate a new strain field.

`MaterialFluxSystem::bdf1_navier_stokes_mixed_rheology_step` resolves the
original `FluidTransport` source at each cell centroid and actual candidate
strain. The existing `RheologyField` caches unchanged cell/rate queries and
reuses one source response throughout cell and face quadrature. Every source
must preserve the reference density: that SAME rho scales inertia, convection,
dynamic viscosity, body force, physical traction, added backflow resistance,
pressure [Pa] and kinetic energy [J/m]. Source failure, including at the final
reassembly, retains its original typed error and refuses before publication.
No changed-density response or fallback viscosity is an accepted endpoint.

`MaterialNavierStokesRheologyStep` returns the physical flow and explicit
boundary policies through `step`, and the exact accepted per-cell source
owners through `samples()`. The original homogeneous transport snapshot remains
the constant-density reference; `samples()` supplies the endpoint viscosities.
Reporting does not query the source again. Captured temperature/composition
stay fixed, and this generalized-Newtonian law adds no material fitting,
thermal evolution, viscoelastic memory or density evolution.

The native and material-bound `bdf2_navier_stokes_mixed_rheology_step`
entries use `(3 u_new - 4 u_n + u_(n-1))/(2 dt)` on the same fixed mesh and
constant reference density. Both history intervals must equal the supplied
timestep. `older_velocity=None` is the original BDF1 startup, bit-for-bit;
two finite histories are otherwise required. One per-step mass-load vector
`(4 u_n-u_(n-1))/3` and effective mass timestep `2 dt/3` reuse the original
quadrature and reaction rows. The actual u_n initializes strain and convection;
the load combination is never used as a physical constitutive query. Invalid
or unrepresentable histories refuse before boundary, force or source callbacks.
No history buffer is allocated in the startup path or inside nonlinear loops.

This extends temporal accuracy for smooth resolved, converged nonlinear flows
to second order, including the BDF1 startup at a fixed positive final time.
Endpoint material owners, physical pressure and instantaneous `mu gamma²`
heating retain their original meaning. The endpoint kinetic energy is not
BDF2's two-history energy, and endpoint heating multiplied by dt is not a
whole-step heat balance. This equal-step entry makes no ALE-history, thermal
evolution or general nonlinear stability/passivity claim.

The native `bdf2_navier_stokes_mixed_variable_step`,
`bdf2_navier_stokes_mixed_variable_viscosity_variable_step` and
`bdf2_navier_stokes_mixed_rheology_variable_step` accept an optional pair
of older velocity and its actual interval to the previous velocity. For current
interval h and previous interval k, r=h/k, the shared mass load is
u_n + r²/(1+2r) * (u_n-u_(n-1)), with effective timestep
h*(1+r)/(1+2r). Actual u_n still initializes the nonlinear solve. Equal
intervals retain the original arithmetic; no older pair uses original BDF1.

The explicit admission 0<r<=2 limits the parasitic zero-equation increment
factor r²/(1+2r) to 4/5. This is a conservative history-stability margin,
not a nonlinear flow energy guarantee or an automatic timestep controller.
Intervals, histories and representable coefficients are checked before source
or boundary callbacks. One history vector is formed per step, never inside
the nonlinear loop; the frozen-viscosity entries perform no strain queries.
The original mixed traction/slip/reservoir policies, endpoint residual checks
and immutable-history refusal semantics remain authoritative.

The material `bdf2_navier_stokes_mixed_variable_step` and
`bdf2_navier_stokes_mixed_rheology_variable_step` bind those same unequal
intervals to the original physical density, dynamic viscosity, force density
and traction units. Accepted source responses remain attached to the endpoint;
reporting preserves physical pressure, kinetic energy and instantaneous
viscous dissipation without requerying material. Changed density or an original
source failure refuses the complete step. The focused material regressions use
quadratic-time affine flow with increasing and decreasing intervals, distinct
densities, exact pressure/heat, original-history strain and final-source refusal.

The focused `navier_stokes_rheology.rs` regressions under `tests/stokes_boundary`
and `tests/material_flow` add an independently derived affine unsteady flow
with nonzero history, live convection and shear-dependent spatial viscosity;
both explicit backflow laws and traction pressure datums; independent endpoint
assembly with lagged-strain and lagged-advection falsifiers; physical material
swaps and analytic kinetic energy; retained final source identities and bounded
query counts; and final-reassembly source/density refusals, invalid coefficients
and exhausted nonlinear work without changed history. They claim no general
nonlinear convergence, empirical material validity or error envelope beyond
their explicit manufactured physics.

The G1/G3 BDF2 additions in those same files cover an independently differentiated
time-quadratic affine flow with spatial rheology and both backflow policies;
temporal refinement of exponential flow with real BDF1 startup; bit-identical
startup; first source queries at actual historical strain; invalid two-history
refusal; final source retention, physical pressure and instantaneous heating
at different material densities; and endpoint source refusal with both
histories preserved. These are numerical manufactured-flow checks, not
empirical validation of a named material or an integrated thermal model.

## Discrete pressure and viscous boundary reactions

Every checked mixed `StokesSolution` now retains `BoundaryReactions` at the
returned field. This includes steady/transient Stokes, Brinkman, frozen Oseen
and accepted nonlinear Navier–Stokes/rheological endpoints. The primitive
assembly keeps only the two discarded normal momentum equations per prescribed
velocity edge, in a rectangular sparse matrix. Its RHS is retained before
identity replacement. Pressure, symmetric stress, BOTH SIP consistency terms,
mass/history, actual convection and any porous drag remain in those equations.
Primal insertion order and strong-constraint arithmetic are unchanged; legacy
assembly without mixed-boundary results does not retain reaction rows.

This is a discrete reaction functional, not a pointwise P0 Cauchy-stress trace.
Let R0 and R1 be the original normal-equation residuals evaluated with the
RETURNED pressure datum. Their virtual work is `R0*m0(v)+R1*m1(v)`, where
`m0=integral(v.n)/edge_length` and `m1=integral(s*v.n)/edge_length`, with
`s` from -1/2 to 1/2 along the edge. Translation and rotation are exactly
represented by the BDM1 trace, including the first normal moment.

The original velocity-boundary kernel assembles full-vector Nitsche terms
before eliminating normal rows. Consequently the recovery restores the FULL
numerical viscous flux `2 nu D(u)n - eta*(u-g)`, not just its tangent, where
`eta=sigma*nu/(2*cell_area/edge_length)` is the actual one-sided penalty.
The adjoint-consistency contribution remains in R; it is not discarded or
replaced by a pressure-only trace. For a global rigid extraction velocity,
`D(v)=0` and its interior jump is zero. The constraint functional plus restored
numerical flux therefore reconstructs the original rigid weak momentum work.
The reported stress is its negative: force and moment exerted by fluid on the
exterior. General Nitsche weak-form reaction extraction is discussed by
D'Angella et al. (2022), Journal of Mechanics 38, 60–76,
<https://doi.org/10.1093/jom/ufac006>; their mesh formulation and convergence
results are not claimed for this BDM1/SIP implementation.

`BoundaryReactionOwner` distinguishes velocity constraints, open traction and
open traction with prescribed reservoir momentum.
A velocity condition may be an inlet/outlet and is never inferred to be an
impermeable mechanical wall. Open edges retain the supplied total stress and
any declared backflow resistance separately. `outward_momentum` uses the
original upwind trace and reversing-face split quadrature; it is advective
transport, not a wall force. For dissipative open boundaries the physical
outward momentum remains `(w.n)u`, while the separately reported stress adds
`min(w.n,0)u`. Their difference is the actual `max(w.n,0)u` operator. No pressure,
stress or momentum transfer is counted twice.

`BoundaryReactions::resultant` selects original boundary edges, an authored
moment origin, and component-indexed pressure offsets. Selected mean-zero
components REQUIRE an explicit offset, including an explicit zero if desired;
traction-datum components retain their imposed stress with None. The pressure
offset adds the corresponding normal-pressure force and moment. Actual applied
references are retained in the result; one chamber cannot supply another's
unknown reference. Empty, duplicate, nonboundary or incompatible selections
and nonfinite resultants refuse. No thermodynamic absolute pressure is inferred.

The native force/moment are divided by the constant reference density. Material
`boundary_loads` accessors return force per depth [N/m] and moment per depth
[(N m)/m = N], using the SAME retained density as inertia and reported pressure.
`rigid_power_w_per_m` evaluates `F.V + M*omega` [W/m] at an explicitly authored
rigid virtual motion. This is conjugate mechanical work, not advective kinetic
energy or a body/mesh time update. Endpoint viscosity, prescribed data and
advection are retained from the accepted assembly; reaction reporting makes no
new source query, nonlinear solve or full-operator assembly. The existing
material-rheology cache and pressure/solver refusal semantics remain intact.

`tests/stokes_boundary/reaction.rs` adds a static pressure box with a traction
datum; affine shear AND normal strain with exact forces, moments and work;
disconnected hydrostatic fields with explicit offsets and a naive cell-pressure
trace falsifier; and non-affine BDF rheological flow whose whole-domain linear
and angular momentum balances include actual SIP slip, nonzero history and
separate natural/dissipative open fluxes. The material rheology regression also
checks SI endpoint reactions and source-query reuse under density swaps. These
regressions were added but not executed in the resource-limited implementation
session. No pointwise traction error bound, moving-domain/FSI integration,
whole-scheme passivity, performance benchmark or new cancellation guarantee is
claimed.

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

## Viscous interfacial surface stress

`capillary::surface::viscosity` supplies the isotropic Newtonian
Boussinesq--Scriven law on the same material triangles. Surface shear viscosity
`mu_s` and dilatational viscosity `kappa_s` are independently prescribed,
nonnegative coefficients in N s/m. They are interfacial transport properties;
neither is inferred from bulk viscosity or numerical mesh thickness. With
tangential projector `P` and projected symmetric velocity gradient `D`,
`tau = 2 mu_s (D - tr(D) P/2) + kappa_s tr(D) P`. The weak nodal force is
`-A tau grad_s(N_i)`, including boundary tractions on open patches. Its power
is minus `A (2 mu_s D_dev:D_dev + kappa_s tr(D)^2)`, a separately reported
nonnegative mechanical loss, never stored surface energy or inferred heat.

Triangle force and simultaneous position/velocity directional-derivative
kernels allocate nothing. The derivative includes changing area, tangential
projector and shape gradients. Uniform and per-face assembly is linear in
mesh size, uses deterministic scalar loss accumulation, and rejects invalid
coefficients, geometry, kinematics and connectivity. The force balances linear
and angular momentum, is independent of winding, and dissipates no rigid
translation or rotation. Five tests check independent analytical shear and
dilation, superposed rigid spin, force tangents, nonuniform patch power balance,
winding and invalid input. No surfactant transport, nonlinear interface
rheology, surface elasticity, bending viscosity or bulk-fluid law is claimed.
Reference: Scriven (1960), doi:10.1016/0009-2509(60)87003-0.

## Inertial capillary material surfaces

With `continuum`, `capillary::surface::dynamics::CapillarySurfaceDynamics`
binds fixed triangle connectivity and strictly positive reference nodal masses
to the existing `fs-time::galpha::SecondOrderState`. It advances actual xyz
positions and velocities through the original matrix-free structural
generalized-alpha owner at `rho_inf=1`. No other time integrator or liquid
solver is introduced. Initial acceleration is physical force divided by mass;
subsequent acceleration retains the native Newmark auxiliary-state semantics.

The inviscid midpoint resistance is the negative original endpoint capillary force.
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

Optional `initial_state_viscous` and `advance_viscous` apply the same
Boussinesq--Scriven surface stress at midpoint geometry and velocity. Its exact
position and velocity partials enter the original matrix-free Newton--Krylov
solve. Acceptance includes the viscous force in momentum and adds
`dt * D_surface` to the energy balance; the report exposes this mechanical loss
as `surface_viscous_dissipation_j`. The existing inviscid entry points retain a
zero-viscosity fast path without new per-face kernels or geometry allocations.
An analytical damped equilateral breathing mode checks the actual midpoint
trajectory and viscous work; a nonuniform rotating surface checks accumulated
loss plus storage and conserved linear/angular momentum over repeated steps.

Prescribed per-face tension change supplies external free energy, never
inferred heat. This rung represents a capillary material surface with fixed
reference mass, including explicitly authored ideal inertial films. It does
not replace an elastic-solid material law, solve bulk incompressibility or
assign physical mass from current area. The caller identifies every interface
(including both sides of a film) and owns fixed-topology applicability.
No contact/support constraints, bending, bulk-fluid viscosity, intermediate-path
nonintersection, remeshing, surface entropy, thermal coupling or dynamics
adjoint is claimed. Temporal and spatial refinement remain necessary.

Five focused G1/G3/G4 tests exercise the analytic midpoint map for an equilateral
breathing patch and second-order temporal convergence; actual rotating,
nonuniformly deforming surface evolution with energy and linear/angular momentum
conservation; explicit coefficient supply, held-force equilibrium and the
analytically accelerated zero-tension limit; finite
rotation/translation objectivity and checkpoint replay; and transactional
invalid-data, cancellation, nonlinear-budget and independent-balance refusals.

## Conserved insoluble surfactant on material surfaces

With `surface-surfactant` (also enabled by `continuum`),
`capillary::surface::surfactant` binds a fixed amount in moles to each material
triangle. Its concentration is `Gamma = N/A`; material motion supplies exact
Lagrangian dilution without creating or discarding surfactant. The isothermal
ideal-dilute and finite-packing Langmuir laws own a surface Helmholtz density
`f(Gamma)`, and tension is `sigma = f - Gamma * df/dGamma`. The supplied thermal
energy is `R*T` in J/mol. Clean tension, concentration reference or saturation
coverage, and the actual interface environment remain explicit physical inputs.
The domain rejects nonfinite/negative amounts, saturation packing and negative
tension; it does not continue an unstable compressed film by clipping its law.

Each triangle stores `E(A,N) = A*f(N/A)`. At fixed amount, its finite-step force
uses the area secant of this SAME energy with the existing Gram area discrete
gradient. The analytic endpoint tangent includes the area dependence of that
secant and its coincident-endpoint limit. Triangle leaves allocate nothing;
uniform-law/per-face-law assembly retains boundary tractions, variable-tension
surface stress, zero resultant and midpoint torque, and exact mechanical work
to roundoff. Zero coating delegates directly to the original clean-tension
force and tangent. Deterministic logarithms and cancellation-resistant energy
increments keep the constitutive evaluation on the existing numerical path.

`initial_state_surfactant` and `advance_surfactant` use the original
`CapillarySurfaceDynamics` midpoint clock, exact Newton--Krylov callbacks and
physical acceptance gates. Fixed borrowed amounts cannot change during a
mechanical trial. The accepted account is `delta K + delta E_surface +
surface_viscous_loss = external_work`; dilution has no prescribed-coefficient
energy supply. Surface viscosity remains optional and uses the original
Boussinesq--Scriven implementation. Invalid coverage, cancelled trials and
failed solves publish no trajectory state. Native regressions compare an
ideal coated breathing triangle against an independently reduced scalar
nonlinear midpoint equation and check repeated rotating nonuniform Langmuir
motion, viscous loss, momentum, energy and invalid-state rollback.

This is an isothermal, fixed-topology insoluble coating model. It supplies no
adsorption/desorption, bulk solute exchange, micelles, surface phase transition,
temperature evolution, thermal heat deposition or remeshing transport.
The constitutive relation and ideal/Langmuir free energies follow Barrett,
Garcke and Nurnberg (2015), equations 2.12--2.14,
doi:10.1051/m2an/2014039. The finite-area force and its tangent are derived from
those energies through the original discrete-gradient surface operator.

`surfactant::diffusion::diffuse_surface_surfactant` evolves these SAME face
amounts on fixed stage areas and a common isothermal law. Each symmetric edge
supplies a nonnegative conductance [m²/s] for the concentration-driven molar
flux. An exact two-face exponential propagator is an area-weighted convex
averaging map: it preserves concentration bounds and moles and decreases the
same convex free energy. Palindromic edge sweeps provide second-order temporal
splitting; a maximum substep and explicit face/edge/work limits bound accuracy
and cost. Propagators are prepared once, with no allocation or exponential in
a sweep. There is no diffusion CFL restriction or global solve. Actual amount
defects and free-energy increments are checked before publication, and signed
admitted roundoff is reported without hiding it behind clipping.

`diffuse_partitioned_surface_surfactant` extends the same exchange to per-face
ideal references or Langmuir site capacities. The authored edge conductance
`k` [mol/s] transports `k*(q_j-q_i)`, with `q=Gamma/Gamma_ref` for ideal dilute
patches and `q=Gamma/Gamma_sat` for Langmuir patches. The conserved pair weights
are the actual capacities `A*Gamma_ref` or `A*Gamma_sat`. Equilibrium therefore
has equal chemical potential, including unequal physical concentrations across
different patches. Ideal references must encode physical standard-state
affinity differences on a common chemical-energy datum; independent arbitrary
reference shifts would change the physical partition equilibrium. A conducting
edge requires a common held thermal energy and constitutive family. Different
clean tensions remain admissible because fixed-area clean energy is constant.

Every pair remains an exact convex map of normalized concentration, so moles,
nonnegativity, the normalized concentration hull and Langmuir packing are
preserved without a CFL bound. The original per-face law is checked at each
exchange to refuse negative tension or any unrepresentable constitutive state,
and the existing original-energy and total-amount acceptance gates remain
atomic. Conductance is explicit interface data; setting `k=Gamma_ref*g` for a
uniform law recovers the original concentration-driven exchange. No geometry,
thermal diffusion or adsorption law is inferred. Analytic unequal-affinity and
stiff Langmuir pairs, an independent capacity-weighted graph spectrum,
homogeneous recovery, extreme normalized concentrations and incompatible-law
or cancelled-state rollback exercise this path. The partitioned operator acts
on fixed stage areas.

`surfactant::diffusion::adaptive` adds mixed absolute/relative concentration
error control for both homogeneous and partitioned held-surface transport. One
full palindromic sweep is compared with two half sweeps using the second-order
fine-path error estimate. Only the unextrapolated fine trajectory is accepted:
the original pair maps retain chemical partition equilibrium, nonnegativity,
packing, amount and free-energy admission. A whole-duration account checks
original-to-final amount and energy again without accumulating tolerances.
Trial and pair-exchange budgets charge coarse and rejected work, while physical
reports include accepted work only. A minimum-step, source, work, cancellation
or physical refusal preserves the entire caller inventory, including after
successful private prefixes. The suggested next step is bounded by the original
controls. The estimate concerns local splitting error, not global or spatial
accuracy. Independent three-face spectra, stiff fast/slow Langmuir modes,
concentration-scaling invariance, equilibrium and late rollback are tested in
`surfactant/diffusion/adaptive/tests.rs`.

`advance_surfactant_diffusive` composes half diffusion on the initial areas,
the original mechanical trial, and half diffusion on its solved final areas.
It publishes trajectory and inventory together only after all stages pass.
The complete account is `delta K + delta E_surface + surface_viscous_loss +
diffusion_free_energy_loss = external_work`. The trial reuses the original
mechanical candidate and never clones accumulated trajectory history. The
returned report separates substage accounts from complete initial/final
storage, amount defect and energy residual. The graph and its conductances
remain held physical inputs throughout the macro step; a changing transport
metric requires its actual owner. No arbitrary-mesh two-point consistency or
bulk-species exchange is inferred. Spatial and temporal refinement remain
necessary. Surfactant inertia is neglected; reference nodal masses remain fixed while
the transported coating contributes interfacial free energy and stress.
Six direct transport tests cover the exact unequal-area mode,
stiff Langmuir bounds, an independent three-face spectral convergence oracle,
stationary/vacuum limits, invariance and rollback. Native composed regressions
exercise moving/viscous surface storage, diffusion loss, momentum, total
amount, empty-graph trajectory recovery and late cancellation of both states.

`advance_partitioned_surfactant_diffusive` reuses this same composition with
`SurfacePartitionedSurfactantTransport`. Each original law stays attached to
its material face through diffusion, the exact mechanical force/tangent, and
both endpoint storage evaluations. Actual evolving areas determine the
capacities `A*Gamma_ref` or `A*Gamma_sat`; temperature and material parameters
remain held. One bounded per-face coating array stages the intermediate
mechanical inventory, without cloning trajectory history. Native moving
ideal/Langmuir regressions independently reconstruct original-law storage and
check complete energy, amount and momentum conservation; homogeneous recovery,
disconnected-law recovery and late cancellation preserve the original behavior
and atomic trajectory/inventory publication.


### Adaptive diffusion through moving material surfaces

`advance_surfactant_diffusive_adaptive` and
`advance_partitioned_surfactant_diffusive_adaptive` use the same original
mechanical composition with adaptive diffusion on its actual initial and solved
final areas. Original per-face chemical laws supply partition transport,
mechanical stress and both energy accounts. The before/after intervals share
one trial and pair-exchange budget, including every coarse and rejected
evaluation; phase-specific budget refusals include all work already spent.
Only accepted fine diffusion trajectories contribute physical loss. The
complete amount allowance and combined diffusion-energy allowance are checked
once more before either the mechanical state/history or coating amounts publish.
The result retains the original coupled physical report plus both adaptive
receipts and total work.

Refinement here resolves internal diffusion edge splitting without rerunning
the mechanical solve for each diffusion substep. It does not estimate the
mechanical macro-step error or diffusion/motion composition error; the caller
must refine that interval independently. Graph conductances remain authored
held inputs. Four native composed regressions define independent exact pair
flows on expanding areas, translating ideal/Langmuir partition equilibrium,
empty-graph mechanical recovery, and shared-budget/final-cancellation rollback.

## Finite aqueous exchange with material surfaces

The optional `surface-soluble` feature adds
`surfactant::soluble::advance_soluble_surface`. The original `AqueousSolution`
and original Langmuir face amounts exchange one explicitly selected nonionic
species through a single finite-bath balance. The existing scalar site owner
supplies backward-Euler adsorption/desorption; its rates obey detailed balance
with the original nonideal bulk activity and one authored surface-minus-bulk
standard chemical potential. Each face retains its own capacity and desorption
rate. The original convex virial model supplies co-solute effects; a zero
interaction matrix is its ideal limit. Missing bulk chemical-reference data
refuses instead of receiving an invented datum.

Existing partitioned surface-diffusion sweeps surround the simultaneous exchange
stage. Sealed faces may receive solute by lateral diffusion. All faces share
one actual bulk debit, while co-solutes and solvent volume remain unchanged.
The complete method is first order because the kinetic stage is backward Euler;
its maximum physical step and work bounds remain explicit. Geometry, topology,
temperature and coefficients are held. Bulk spatial transport, moving mechanics,
micellization, ionic adsorption and adsorption heat are outside this rung.

The physical account restores the changing bulk reference term omitted by the
source's reference-cancelled energy primitive, adds original surface energies
and the authored standard-potential contribution, and checks both stage and
original-to-final free-energy changes. Amount conservation, finite packing,
source support, tension, work and final cancellation pass before either original
inventory is replaced. Six actual-source regressions cover the independent
finite-bath quadratic and Riccati transient, nonideal activity and reference
invariance, stationary gross turnover, connected diffusion/exchange equilibrium,
and late physical/work/cancellation rollback.

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

## Mixed-boundary Newton flow for stiff viscosity

`FluxSystem::solve_stokes_mixed_rheology_newton` resolves steady
generalized-Newtonian Stokes flow under velocity, total-traction and finite or
perfect Navier-slip boundaries. It reuses the existing sparse constitutive
Jacobian and bounded Newton/Armijo owner. The supplied positive finite
`nu(gamma)` and analytic `dnu/dgamma` are evaluated at the same physical rate;
a smooth radial law has zero derivative at zero strain. No fitted viscosity,
rate floor, artificial timestep or altered boundary law is introduced.

Boundary data is sampled once. Natural faces retain their normal unknowns,
continuity equations and pressure datum. Prescribed total traction and fixed
Navier resistance have no bulk-viscosity derivative; only velocity-Dirichlet
faces add constitutive boundary-load terms. All cell and harmonic interface
terms retain the original derivatives. `linearize_stokes_mixed_rheology`
exposes this same steady Jacobian for the existing checked transpose solve.

The original, row-scaled and component-pressure-translated residuals must all
pass before publication. The translation reuses the linear solver's origin
arithmetic so a large common traction pressure cannot hide weak unconverged
flow. Final pressure interpretation and boundary reactions are rebuilt at the
accepted viscosity and velocity. Initial state remains unchanged on refusal.
Temperatures, composition and reference density stay fixed; this adds no
convective mixed Newton solve, thermal evolution, viscoelastic memory or
unregularized yield stress. The material-bound entry below supplies physical
units and retains the original source owners through this same Newton solve.

Focused tests in `mixed_rheology.rs` use exactly representable pressure-driven
extension and layered Couette flow to check nonlinear stress, pressure-datum
shifts, bulk dissipation, Navier-slip loss and undamped Picard failure.
`rheology_linearization.rs` compares mixed Jacobian actions with independent
primal differences at finite/perfect slip and spatial viscosity contrast.
Invalid tangent refusal preserves the input state. These new tests are authored
but **not executed**: the session's executor became unavailable before native
compilation. Their presence is not a convergence, performance or validation claim.

## Material-bound mixed Newton flow

`MaterialFluxSystem::solve_stokes_mixed_rheology_newton` connects the original
source-owned viscosity and analytic `dmu/dgamma` to the existing mixed-boundary
Newton owner. One immutable reference density normalizes dynamic viscosity,
its derivative, physical force density [N/m³], total traction [Pa] and Navier
resistance [Pa s/m]. No fitted duplicate rheology, scalar pressure shift or
alternate nonlinear solver is introduced. The material response and derivative
share the same cell, physical position, strain rate and fixed temperature /
composition. Density changes, unrepresentable derivatives, nonfinite tangents
and nonzero radial derivatives at zero strain refuse the solve.

`MaterialStokesRheologyNewtonSolution` retains the existing physical mixed
flow owner, the accepted dynamic-viscosity derivatives [Pa s²] and the actual
backtrack count. Its `flow` supplies velocity, pressure in the original datum,
kinetic energy, original endpoint sources and physical boundary reactions.
Bulk heat remains the original `mu * gamma²`; finite-slip wall heat remains
a separate boundary loss. The shared cell/rate cache reads each source once
and supplies the final stress and heat without a reporting requery.

Source-support failures are terminal even in a line-search trial or the final
reaction assembly. They cannot be hidden by taking a smaller step, accepting
fallback viscosity or returning a partial physical endpoint. The original
typed material error has precedence, and caller state is immutable. Scope is
steady creeping generalized-Newtonian flow with smooth finite viscosity and
fixed constitutive context; no thermal evolution or variable-density momentum
is claimed.

Focused tests in `material_flow/newton/mixed.rs` cover exact pressure-driven
extension at two densities and two pressure datums, failure of undamped Picard
for the same stiff law, retained original layer sources, finite Navier slip,
independent bulk/wall/work closure, and an actual thermal Carreau--Yasuda source
with its analytic derivative. Refusal tests exercise a genuine Newton trial,
density changes, invalid tangents and invalid initial state. A private checked
owner test isolates source loss in the final reaction assembly. All changed
Rust sources pass the pinned rustfmt parser. These new tests are authored but
**not executed** in this environment; parser acceptance is not compilation,
convergence, performance or physical validation evidence.

## Original material Newton flow with mixed transient transport

`MaterialFluxSystem::solve_navier_stokes_mixed_rheology_newton`,
`bdf1_navier_stokes_mixed_rheology_newton_step` and
`bdf2_navier_stokes_mixed_rheology_newton_variable_step` bind the original
material response and analytic `dmu/dgamma` to the existing mixed convective
Newton owner. They reuse the same sparse constitutive/advection Jacobian,
bounded Armijo search, history normalizer, boundary admission and reaction
assembly. No additional nonlinear solver or fitted viscosity is introduced.

One fixed reference density scales dynamic viscosity, its derivative, force
density, mass, convection, total traction, Navier resistance and recovered
pressure. Original reservoir inflow, literal traction and declared dissipative
backflow retain their distinct momentum laws. Finite BDF mass anchors all-open
components without an artificial velocity constraint. Unequal BDF2 retains
both original histories and the native step-ratio limit; only the mass load is
combined, while actual previous velocity initializes the constitutive solve.

`MaterialNavierStokesNewtonStep` retains physical flow, exact endpoint source
owners and derivatives, Newton/backtrack counts, pressure datums, bulk heating
and boundary reactions. Unchanged cell/rate queries use the existing cache;
reporting neither resolves nor reads original material transport again.
Bulk `mu*gamma²`, Navier-slip loss, advected kinetic energy and numerical
upwind dissipation remain separate. These are endpoint rates, not integrated
thermal deposition or an asserted BDF2 energy identity.

Source failures remain terminal during Newton trials and final reaction
assembly; original typed errors precede numerical refusals. Invalid histories
are rejected before source/load callbacks. No fallback material, changed
history, clipped timestep or partial endpoint is returned. Density and captured
temperature/composition stay fixed. Existing synchronous solver/cancellation
limits, semismooth upwind qualifications and conditional Newton convergence
remain; no thermal evolution, variable-density flow or viscoelastic memory is
added.

Four focused regressions in `material_flow/newton/mixed_transient.rs` cover
unequal BDF2 and BDF1 layered shear with live transverse convection, original
source/history identity, physical pressure/density/heat, finite slip, reservoir
transport, independent BDF1 work/kinetic balance, an actual thermal Carreau
all-open acceleration, steady natural/dissipative boundary policies and
refusal/retry. A private native test isolates source loss in final reactions.
The changed Rust sources pass the pinned rustfmt parser. Native compilation
and these regressions remain **unexecuted** because the current workspace lacks
its constellation dependencies and persistent storage capacity; parsing and
independent manufactured algebra are not native numerical validation.

## Original material SDIRK2 with mixed transient boundaries

`MaterialFluxSystem::sdirk2_navier_stokes_mixed_rheology_newton_step` extends
the source-owned mixed Newton flow to the existing two-stage, second-order,
stiffly accurate SDIRK2 method, with `gamma = 1 - 1/sqrt(2)`. It needs only
the actual previous velocity and current interval. The first stage and endpoint
take separate borrowed boundary slices, allowing original total traction,
reservoir velocity, normal wall velocity and finite Navier resistance to vary
at their actual stage times. Force and material callbacks receive that same
stage. No alternate stress, convection, boundary or nonlinear solver is added.

The shared time owner validates history and representable stage mass before
callbacks. The existing sparse mixed Newton owner evaluates stress, analytic
rate derivatives, convection and reactions with the same immutable reference
density. Its accepted first-stage velocity initializes the endpoint solve;
the SDIRK history combination enters only the inertial load. Each stage has
its own constitutive cache and retains its exact original sources and tangents.
Finite stage mass anchors all-open components without artificial velocity pins.

`MaterialNavierStokesSdirk2Step` returns both accepted stages only after the
complete step succeeds. `FluxSdirk2Error` preserves the failing stage and
original material/numerical error. Source refusals remain terminal during
Newton trials and final reaction assembly; a failed endpoint changes no caller
history and publishes no partial pair. Newton work limits apply separately
to both stages. There is no hidden clipping, substepping or boundary fallback.

Physical pressure, bulk molecular power, Navier loss, advected kinetic power
and upwind dissipation retain their existing separate stage owners. The
returned `dissipated_heat` also integrates original bulk and slip powers with
the SDIRK2 weights, as specified below. This is generated heat, not thermal
deposition or a full discrete energy balance.
Temperature/composition are fixed external fields within each stage and may
differ between stages; density remains fixed. Smooth temporal data and resolved
spatial/nonlinear errors are required for second order. Linear L-stability
does not establish nonlinear positivity, passivity or global flow uniqueness.
The existing synchronous solver and cancellation limits remain unchanged.

Three regressions in `material_flow/newton/mixed_sdirk2.rs` cover a layered
original thermal Carreau flow with live transverse convection, stage-specific
slip/traction/reservoir data, pressure/density/heat and source identity; all-open
cubic-in-time translation with second-order refinement; and terminal source,
input/boundary refusal and clean retry. Independent scalar quadrature gives
translation errors 2.426406871e-4, 6.066017178e-5 and 1.516504294e-5 for 4, 8
and 16 intervals. All changed Rust files pass pinned rustfmt parsing. Native
compilation and these new regressions remain **unexecuted** in the current
workspace because constellation dependencies and persistent storage capacity
are unavailable; the scalar oracle is not native PDE validation.

## Adaptive physical time for original mixed material SDIRK2

`MaterialFluxSystem::advance_navier_stokes_mixed_rheology_sdirk2` advances a
positive-duration forward interval using the SAME original mixed material
SDIRK2 and Newton owners. `material::adaptive` exposes its controls, physical
time boundary data, stage-aware errors and accepted result. Every comparison
runs one coarse interval and two half intervals from the same accepted history.
Only the fine endpoint advances history, without Richardson extrapolation.
The result retains the last accepted fine pair and its exact source owners;
coarse and rejected candidates never become accepted material states.

Force density, original material and borrowed velocity/traction/reservoir/slip
laws receive actual absolute stage times, including every retry. The original
boundary policies, source density, rate tangents, sparse Jacobian, line search
and physical reaction conversion remain authoritative. Topology stays fixed;
wall resistance and other physical data may vary with time. Callbacks must be
pure in their declared time, position and state inputs. Initial histories,
controller bounds and representable interior stage times are admitted before
stage callbacks. The exact requested final time is retained.

The error indicator is
`RMS((u_fine-u_coarse)/3)/(atol + rtol*max(RMS(u_old),RMS(u_fine)))`.
RMS uses the existing finite-element kinetic-energy quadrature divided by mesh
area, preserving physical velocity units independently of edge-moment count,
orientation and domain area. One difference vector is reused, and the accepted
history norm is cached across rejected trials. The shared `fs-time::PiController`
settings use order-three local-error gains; its RK45 integrator is not reused
as a different flow solver. Proposal and previous accepted error are returned
for continuation. Clipped final steps do not replace the untrimmed proposal.

The global attempt limit includes every accepted and rejected comparison,
with at most six native Newton stages per attempt. Each retains the original
Newton, backtracking and Krylov limits. Only numerical correction-solve,
line-search and nonlinear-convergence failures may shorten a step. Original
material, boundary, input and representability refusals are terminal, preserving
the comparison, failing stage and physical time. Exhausted attempts, a blocked
step floor or cancellation return no partial advancement. Cancellation is
checked before each native two-stage solve and before publication; synchronous
Newton is not interrupted inside that pair. Caller history remains immutable.

This controls local velocity time error under smooth data and resolved
spatial/nonlinear errors. It does not certify pressure or stress, global error,
adaptive derivatives, positivity, passivity or thermal evolution. Two strongly
damped endpoints can agree while underresolving fast-mode dissipation history.
Retained stage powers and transport remain available. Whole-interval
`dissipated_heat` accumulates only accepted fine-stage production as specified
below; velocity acceptance is not a heat-accuracy certificate.

Three focused regression groups cover absolute stage times and retained
thermal Carreau sources; tolerance refinement and physical RMS invariance
under changes in mesh and domain area; time-varying layered Carreau slip,
traction and reservoir data; and work/source/boundary/cancellation/time-
representability refusal with clean retry. All six changed Rust files pass
pinned rustfmt parsing. Native compilation and these regressions remain
**unexecuted** because the current workspace lacks constellation dependencies
and persistent storage capacity. An independent scalar cubic-translation
calculation gives accumulated errors 1.43913e-4 and 5.94509e-6 for absolute
tolerances 1e-4 and 1e-6; this is not native PDE validation.

## Accepted original-stage viscous and Navier-slip heat

The original mixed SDIRK2 owner now returns `MaterialDissipatedHeat` with every
successful fixed step. For original triangle i, its bulk energy is
`h*area_i*((1-gamma)*Phi_i(stage1) + gamma*Phi_i(stage2))` [J/m], where Phi
is the existing source-owned `mu*gamma_dot^2` power density [W/m³]. Each
boundary edge uses the SAME positive weights on its retained Navier-slip
power `integral beta*|(u-wall)_t|^2` [W/m]. The method supplies its own
original timestep and cell areas; callers cannot reconstruct or substitute
those metadata after a solve. No constitutive, strain or heating kernel is
duplicated, and no material or boundary callback is repeated.

`DiscreteBoundaryReaction::slip_dissipation` exposes the existing retained
kinematic scalar. `MaterialStokesSolution::boundary_slip_dissipation` provides
an allocation-free checked view in physical W/m using the same immutable
reference density. Every entry retains its original edge and connected
component. Heat needs neither a pressure datum nor a moment origin and does
not allocate/evaluate full boundary wrenches. Non-slip and perfect-slip edges
remain present with zero power.

The integrated result exposes elapsed seconds, original-cell bulk energies,
original-edge/component slip energies, and separate compensated bulk/slip
totals in J/m of out-of-plane span. Dividing a cell energy by its original
area gives integrated density in J/m³. Original fixed-step material owners
remain on the returned two stages. A shared compensated accumulator uses
O(cells + boundary edges) storage and refuses nonfinite, negative or
unrepresentable contributions and sums.

Adaptive advancement sums these integrals from BOTH accepted fine half steps
over the ENTIRE requested interval. Coarse comparison heat, rejected attempts
and failed stage pairs contribute nothing. Earlier accepted source objects
are not archived: spatial identity and cumulative physical heat are retained
in bounded storage, while the final fine stages retain their exact original
sources. Any later source, boundary, cancellation or representation refusal
publishes neither partial flow nor partial heat. An accumulated-heat refusal
is terminal and retains the physical endpoint time; it is not a smaller-step
retry signal.

Bulk heat and Navier-slip production remain distinct. Slip heat is explicitly
unassigned between fluid and wall; a thermal owner must choose that partition
and deposit energy. These outputs exclude mechanical wall/pressure work,
advected kinetic or thermal energy, numerical upwind/SIP loss, and caloric
state evolution. They are not a full discrete energy identity, positivity
theorem or passivity certificate. SDIRK2 heat quadrature is second order under
resolved smooth source/strain history and spatial/nonlinear error. The default
adaptive entry estimates velocity error only: agreement of strongly damped
velocity endpoints cannot certify stress or dissipation-history accuracy.
The optional heat-aware entry below controls additional physical indicators.

Focused regressions extend the layered original thermal-Carreau fixture to
check stage-weighted physical bulk/slip heat and original cell/edge identities.
A manufactured affine shear u=(1+t)y with mu=0.6 Pa s and two Navier walls of
beta=1.2 Pa s/m has exact bulk and total slip heat each 1.4 J/m over [0,1] s.
An independent scalar quadrature gives errors 7.582521472e-4, 1.895630368e-4
and 4.739075920e-5 for 4, 8 and 16 intervals. A separate cubic translation plus
steady shear fixture exercises rejected adaptive comparisons, whole-interval
heat, a late original-source refusal and clean retry. Changed Rust sources
pass pinned rustfmt parsing; native compilation and these PDE regressions
remain **unexecuted** because constellation dependencies and persistent storage
capacity are unavailable. The scalar quadrature is not native PDE validation.

### Heat-sensitive original-material time integration

`advance_navier_stokes_mixed_rheology_sdirk2_with_heat_control` adds
`MaterialSdirk2HeatControls` to the existing original-material adaptive driver.
Existing velocity-only controls and entry points retain their signatures.
The new result contains the same accepted flow and integrated original heat,
plus a separate `MaterialSdirk2HeatReport`. The public types are available
through `material::adaptive`.

Each comparison uses the same six native Newton stages already required by
step doubling. It compares the coarse cell bulk heat against the sum of
both fine halves, and independently compares original edge/component slip
heat. Local differences are divided by three for the second-order fine-path
estimate. Separate bulk/slip sums of absolute spatial discrepancies prevent
opposite local errors from cancelling in a total. Absolute heat tolerances
are in J/m; relative scales use only that comparison's physical production.
No extra constitutive callbacks, fitted viscosity or substitute heating
operator is introduced.

Each original fixed SDIRK2 step additionally retains its signed
`temporal_kinetic_defect_j_per_m`. For original stage velocities Y1, Y2 and
initial Un, define `d1=Y1-Un` and
`d2=Y2-Un-((1-gamma)/gamma)*d1`. The time-quadrature defect is

`D = rho*(K(d2)-K(d1)) = (rho/2)*<d2-d1, M*(d2+d1)>`.

Here K and M are the native BDM velocity mass metric and rho is the original
immutable material density. The direct signed pairing uses the original
triangle/basis quadrature with double-double accumulation; centered stage
increments avoid subtracting two large kinetic totals. The consumed
mass-history allocation is reused and one O(velocity DOFs) sum vector is
added. There is no new mass assembly or linear solve. The diagnostic concerns
the complete time discretization, including force and boundary histories;
it is not a viscous-only loss, a mesh/SIP/upwind loss or a full energy balance.

SDIRK2 permits BOTH signs of D. The heat-aware controller therefore charges
`|D_first_half|+|D_second_half|`, before any sign cancellation, and separately
reports the signed accumulated value. With absolute allowance A, relative
fraction r, total requested duration T and accepted original bulk-plus-slip
heat H, the whole-call allowance is `sum(abs(D)) <= A+r*H`.
A comparison receives `A*dt/T+r*H_comparison`; unused local credit is not
borrowed. This deliberate allocation prevents an early transient from
consuming later absolute credit, but a very short stiff transient inside a
long requested window can require many steps. It is a resolution policy,
not an optimal-cost or certified global-error theorem.

Local costs, double-double accepted histories and a separate whole-call
inequality keep both defect signs and small positive increments visible.
The common PI proposal sees the largest velocity, heat or defect indicator.
Rejected/coarse histories contribute no accepted heat or defect. Attempt,
Newton, source, boundary, representability and cancellation limits remain
authoritative, including a final cancellation check after the accounts are
complete. Refusal returns no partial flow, heat or time-defect result.
Neither sign of D is deposited as physical heat or used to correct the
original molecular/slip heat.

The existing native mass and material-flow test seams cover affine moment
pairing, exact-velocity/non-exact-heat refinement, a stiff unforced vortical
flow, source/work/cancellation refusal and clean retry. An independently
executed scalar decay calculation exposes the missed mode directly: with
lambda*h=1e6, the fine endpoint velocity estimate is only about 1.61e-6 and
the heat estimate about 5.02e-6, yet its two fine steps carry about 0.49997
absolute time defect from an initial kinetic energy of 0.5. At lambda*h=4,
the coarse quadrature instead overproduces heat, giving about 0.62252 and
negative defect about -0.13222. The signed diagnostic and absolute budget
distinguish both failures. Pinned rustfmt parsing and these independent
equation checks do not replace native Rust compilation or PDE execution;
those remain unexecuted under the current dependency/storage constraints.

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

### Hydrostatic density stratification in the original planar owner

`Euler2d::new_isothermal_atmosphere` initializes exponential volume averages
and retains uniform gravity along x or y. It reuses the original `fs-gas`
`IsothermalAtmosphere` integration and cached thermodynamic cell/face states;
only a compact axial profile survives initialization. At least two cells and
physical ends are required along gravity. The transverse axis retains its
periodic, reflecting or transmissive semantics. `bind_isothermal_atmosphere`
admits caller-authored disturbed averages before evolution, changes no state
or energy gauge, and refuses replacement of an existing gravity declaration.

Both first-order and MUSCL reconstruction use density and pressure relative
to that reference, with physical velocities. This keeps a positive reduced-
density atmosphere reconstructible instead of producing negative faces from
large additive deviations. Relative coordinates must be representable; an
inadmissible zero-slope physical face uses a local physical-cell trace without
clipping its conservative state. Every evolving face retains the original
source-owned Riemann or stationary-wall closure. Transmissive ends extrapolate
the disturbance while continuing the reference to the physical face.

Momentum cancels the source-resolved reference face pressure before time
scaling, then adds `g*(rho-rho_reference)`. Gravity work uses the same shared
mass fluxes as continuity. Consequently the SSP RK2 update conserves
`sum V*(E + rho*phi)` with `phi=-g*(x-x_reference)`, up to floating-point
arithmetic and boundary transport of `energy_flux + phi_face*mass_flux`.
`gravitational_potential_energy()` exposes this midpoint potential coordinate;
`StepReport::source_change` includes the actual gravitational impulse and work.
The original accepted state, reports and boundary loads remain transactional.
Additional authored forcing composes with this retained gravity, and the
coupled molecular owner retains the same reference in its trial field.

The isothermal reference is not a thermostat: disturbances evolve full gas
energy and the original variable heat capacity. Native regressions cover
exact resting states with more than 300:1 density contrast, both axes and both
reconstructions, a two-dimensional disturbance, coarse strongly reduced-density
states, closed gas-plus-potential conservation, open potential-energy transport,
additional acceleration and atomic cancellation/reference refusal. This
bounded reference method does not provide arbitrary hydrostatic profiles,
rotating equilibria, changing-composition well balancing, an invariant-domain
theorem or experimental validation.

`Euler2d::new_isentropic_atmosphere` and `bind_isentropic_atmosphere` use the
same pressure/gravity residual with the original `fs-gas` isentropic profile.
That owner integrates conserved density and energy under explicit
`AtmosphereQuadrature` controls. Temperature follows the source enthalpy
difference and pressure follows its entropy difference, retaining variable
heat capacity rather than freezing gamma or imposing a fixed lapse rate.
Only the reference has uniform entropy; disturbances retain the original
Euler energy evolution. `hydrostatic_reference()` reports either physical
family, while `atmosphere()` continues to identify only isothermal references.
The molecular/Euler/molecular owner retains the original family and quadrature
in its trial state. No profile integration or allocation is added to stepping.

The isentropic regression target independently integrates a linear-Cp source
for face temperature, pressure and total mass, then checks resting states,
two-dimensional pressure disturbances, combined gas/potential conservation,
source-window and preparation-budget refusal, immutable reference binding,
and cancellation without state publication. Molecular tests distinguish a
resting viscous isentrope without heat conduction from the physically evolving
temperature gradient that nonzero Fourier conduction would introduce.

### Hydrostatic pressure balance in the original volumetric owner

`compressible::spatial::Euler3d::new_hydrostatic_atmosphere` extends the same
source-owned isothermal, isentropic and prescribed thermal-lapse references
to any Cartesian axis. `bind_hydrostatic_reference` admits an authored initial
disturbance without changing its conserved averages. The gravity axis requires
at least two cells and physical ends; both transverse directions retain their
native boundary conditions. Preparation integrates only one axial profile;
stepping adds no profile integration or allocation.

Both native reconstruction orders use relative density and pressure with all
three physical velocities. The momentum residual removes the actual reference
Riemann traction before differencing, then adds the density-perturbation force.
Gravity work uses the shared continuity flux. Closed inviscid evolution thus
conserves gas plus midpoint potential energy, including variable source heat
capacity. `StepReport::source_change` includes gravity impulse and work, and
pressure-only outlets admit the reconstructed physical face rather than the
stratified volume average. Additional forcing composes with gravity; molecular
coupling retains the reference through accepted/trial swaps without adding
gravity during the separate diffusion substeps.

Native regressions cover all three axes, both gravity signs and reconstruction
orders, independently integrated mass, exact resting states, three-dimensional
disturbances, transverse momentum, pressure outlets at both physical ends,
combined gas/potential conservation, additional forcing, molecular trial reuse
and transactional cancellation/refusal. This bounded constant-gravity method
does not claim arbitrary or rotating equilibria, mixture well balancing or
suppression of physical heat flow through a nonisothermal atmosphere.

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

Periodic pairs reuse the same flux. Reflecting walls default to stationary
adiabatic slip, retaining their normal viscous support impulse with zero heat,
work and tangential traction. `MolecularTransport2d::with_wall` and
`MolecularEuler2d::with_wall` additionally admit `MolecularWall2d::NoSlip` with
finite tangential velocity and adiabatic or prescribed-temperature heat transfer.
Resolved shear and Fourier heat use the physical half-cell distance; their
conservative energy flux uses the actual wall velocity, so a stationary
adiabatic no-slip wall exerts drag with zero boundary work. The original law is
evaluated at the cell/wall mean temperature for an isothermal surface. No
separate viscous-heating term is added. Normal wall motion, nonreflecting sides
and out-of-domain wall temperatures refuse. Transmissive ends retain zero
molecular flux. Molecular impulses and balances include physical grid depth.

The authored wall velocity is constant along its side, so its tangential
derivatives at the physical wall are zero even in a nonuniform boundary layer.
No-slip velocity ghosts enter the adjacent interior cross derivatives. Four
cell-gradient entries are prepared once and reused by neighboring faces,
replacing repeated neighbor reconstruction. The shared `MolecularWall<D>`
condition and frozen-source shear/heat/work kernel also serve the existing
three-dimensional owner; `MolecularWall3d` retains its original public paths.
The evolving-species owner continues to use its original shared gas-source
molecular face kernel.

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
Changing an authored wall invalidates the complete coupled owner's prior report
so old exchange cannot be relabeled under a new wall condition. Identical
reconfiguration preserves accepted history; cancellation and source refusal
publish no partial field, time or per-face exchange.

### Stabilized frozen-composition heat and momentum diffusion

Both `compressible::molecular::MolecularTransport2d::advance_stabilized`
and `compressible::spatial::molecular::MolecularTransport3d::advance_stabilized`
apply the original gas-owned second-order RKL2 coefficients and the shared
multidimensional interval controller to these frozen-composition operators.
`MolecularControl` retains physical macrostep, CFL and retry limits;
`MolecularStabilization` separately declares maximum stages, cumulative actual
operator evaluations, and an explicit cap on EXTRA retained numerical storage.
The existing flow, source-law, molecular and optional SGS allocations are outside
this extra-storage cap. Required and actual retained extra bytes are queryable.

The three original state arrays are reused. Additional anchor-face and
boundary-stage arrays are allocated lazily once. Every stage and the final
endpoint use the original local transport data, caloric inversion, full
Newtonian stress/work and Fourier heat. Three-dimensional SGS, when explicitly
selected, is reevaluated by its original owner. One conservative face flux
serves neighboring cells; the SAME signed Legendre recurrence integrates
per-face boundary heat and mechanical impulse/work.

The controller selects the smallest sufficient stage count up to the declared
limit and bounds the negative-real linear parabolic number by
`CFL * (s*s+s-2)/4`. This is not a nonlinear positivity or entropy guarantee.
Physical refusals halve private trials without clipping; the caller must choose
an accuracy-appropriate maximum macrostep because no time-error estimate is
provided. Work includes failed stage prefixes, recomputed retry anchors and
endpoint admission. Source failures remain typed to the original owner.
Cancellation, exhausted work and insufficient extra storage cannot publish a
partial field or boundary load; successful pure molecular updates do not
advance Euler time.

The same frozen-composition `MolecularEuler2d` and `MolecularEuler3d` owners
expose `step_stabilized` and `step_stabilized_with_forcing`. Both route through
the ORIGINAL complete molecular/advection/molecular retry loop, with RKL2
replacing only the two molecular integrators. The Euler reconstruction, forced
impulse/work, Mach admission, reservoir CFL, wall conditions and selected SGS
operators retain their original implementations. Any shorter or refused middle
step restarts the complete split from its accepted state.

One actual-evaluation counter spans BOTH molecular halves and every discarded
whole-step trial. It is never reset by a retry, and exhaustion terminates the
call immediately. Returned stage/parabolic extrema describe the final accepted
halves; total work also includes discarded attempts. The stabilized report
wraps the original physical report, which remains available through
`last_step()`. Complete coupled refusal preserves the accepted field, clock,
all boundary loads and prior report. Regressions cover accelerated thermal
coupling, exact uniform-force impulse/work, shared half-step budgets, spent
work across CFL retries, final-duration replay, late cancellation and heated
wall energy across both halves.

Focused regressions compare planar and anisotropic volumetric thermal modes
with independently derived discrete Laplacian eigenvalues, check time
refinement and reduced operator work, recover viscous shear decay and conversion
of kinetic energy to internal energy, close heated-wall transfer using the
signed recurrence, and exercise byte/work/cancellation refusal and replay.
These are manufactured source fixtures, not experimental calibration.

Five focused G0/G1/G2/G3/G4 regressions exercise a 2-D Fourier mode, oblique
longitudinal and transverse Newtonian modes, collision-diameter-only pressure
and shear changes, wall impulse/energy closure, late cancellation/replay,
zero-transport recovery and atomic work refusal. Synthetic source and collision
data test composition and conservation; they are not measurements of a real gas.
This extension keeps spatially uniform frozen composition. It does not supply
species diffusion, chemistry, Soret/Dufour effects, turbulence,
rarefied gas, physical shock-thickness validation or mutual fluid/solid thermal
evolution. Those capabilities require their original owners and additional
admitted source data, rather than inferred closures.

Seven additional planar wall regressions exercise independent all-four-side
shear/heat/power, exact co-moving translation, native Couette momentum and wall
heating under a common tangential boost, absence of spurious wall-normal stress,
agreement with the actual three-dimensional extrusion, late cancellation/history
and source-temperature admission. These are resolved continuum stencil and
conservation checks, not a coarse-grid wall function or experimental validation.

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

### Translating and piston-driven volumetric gas

`Euler3d::step_moving` advances the original frozen-gas field on a translating
Cartesian box with at most one deforming axis. `MeshMotion` gives the held
laboratory velocity of each minimum/maximum plane; interior planes interpolate
those speeds. The accepted grid retains its origin and spacing, and physical
cell/face coordinates include that origin. A periodic axis admits translation
only. Positive endpoint widths exclude collapse throughout the affine motion.

The existing unsplit MUSCL/SSP RK2 owner advances extensive `V*U` using the
actual stage face areas and divides by the actual stage volumes. A single
deforming axis makes volume affine, so the same quadrature satisfies the
discrete geometric conservation law. Shared interior faces retain one ALE
HLLC/HLL flux with opposite neighbour signs; periodic endpoints copy its bits.
The original source owns relative advection, all three momenta, caloric
admission and laboratory energy. Moving slip walls use the original nonlinear
shock/expansion closure and reflect MUSCL ghosts about the actual wall speed.
Their mass transfer is zero and their energy transfer is wall speed times
normal pressure impulse. Accepted per-face loads use both actual RK stages.

Relative face-wave CFL checks, caloric refusals and bounded retries apply to
the complete field/geometry transaction. Existing preallocated buffers are
reused. Refusal or cancellation preserves cells, grid, clock, loads and report;
exactly zero motion delegates to the original stationary entry. Focused native
regressions cover deforming free streams/GCL on all axes, comoving periodic
shear and density contacts, piston mass/impulse/work, refinement toward exact
homogeneous adiabatic expansion, Galilean covariance and rollback/replay.

This entry admits moving reflecting walls and unbound zero-gradient openings.
Retained hydrostatics, prescribed exterior traces/profiles, pressure-only and
characteristic outlets, simultaneous multi-axis deformation, remapping,
molecular moving-volume composition and reciprocal FSI remain unsupported.
The frozen source and existing opt-in engineering authority are unchanged.

### Reciprocal finite-piston volumetric gas

`Euler3d::step_with_inertial_walls` evolves the original three-dimensional
gas together with both finite guided pistons on a selected Cartesian axis.
It reuses `fs-gas`'s simultaneous two-wall midpoint-speed solve, original
finite masses, bounded Newton/backtracking, and independent wall/combined
impulse and work acceptance. Each proposal advances the same accepted field
with the native ALE MUSCL/SSP RK2 operator at one fixed duration. All retries
reuse the field scratch; no trial clones the gas or publishes elapsed time.

The actual stage-integrated boundary pressure impulses accelerate the two
walls, and those same impulses paired with midpoint speeds supply the gas
work. Exterior pressure and additional centre-applied axial forces remain
explicit. The selected axis has two reflecting ends; transverse directions
are periodic or stationary reflecting walls. Combined axial momentum and
gas-plus-wall energy are checked against actual changes; transverse support
impulses remain in the original gas report. Nonuniform pressure also produces
per-piston angular impulse about each plane's geometric centre, with the
equal and opposite ideal-guide angular impulse reported explicitly. The
orientation guides have zero rotational velocity and do no rotational work.

Accepted cells, grid, time, wall velocities, per-face loads and the gas report
publish together after the final cancellation gate. Unrepresentable plane
travel refuses even if a small displacement disappears only at the upper
endpoint. Native regressions cover equilibrium/translation on every axis,
pressure-area/mass acceleration, repeated unequal-mass expansion and combined
energy/momentum, exact replay through prescribed native ALE motion, nonuniform
pressure torque, axial Galilean boosts, complete retries and atomic rollback.
The original frozen-composition inviscid authority remains Estimated; gravity,
leakage, open transverse faces, prescribed exterior profiles, molecular heat
or stress, piston rotation, surface deformation and contact stops are outside
this entry. The guides and stationary transverse walls are explicit supports.


### Conservative body forcing on the original Cartesian owners

The planar `Euler2d::step_with_forcing` now supplies the same physical source
law and symmetric source/advection/source composition on its original SSP RK2
transaction. `BodyForce2d` distinguishes uniform or per-cell acceleration
[m/s²] from force density [N/m³]; acceleration multiplies the density of each
current source stage. The existing 2-D face closure, reconstruction, boundary
loads and retry policy remain the transport owner. The 2-D and 3-D source
types and kick arithmetic share one dimension-generic implementation, including
the frozen, mixture and reactive 3-D consumers.

Planar `StepReport::source_change` exposes body impulse and signed mechanical
work separately from its per-face boundary exchange. Source-only work is the
represented momentum increment dotted with midpoint velocity, so braking
removes kinetic energy without becoming artificial internal heating. Source
calorics, Mach admission at accepted/kicked/RK cell states, both advective
face CFL checks and the final kicked-cell acoustic CFL all apply before atomic
publication. An exactly zero authored source follows the ordinary step after
validating source shape, controls and scratch.

`ForcingWorkspace2d` admits and reuses one 32-byte conservative state per cell,
with checked required bytes and actual capacity admission. Source loops and
retries allocate no additional arrays. Focused regressions cover unequal
density acceleration versus volumetric forcing, variable-Cp signed work,
spatial source fields, zero-source bit equality, force-driven Mach/CFL retries,
budget refusal and final-checkpoint cancellation/replay. The same held-source
boundaries below apply. The split alone is not hydrostatic well balancing;
the separately selected planar atmosphere reference above supplies that path.

`MolecularEuler2d::step_with_forcing` reuses this forced interval between its
original Fourier/Newtonian molecular halves. The complete split restarts from
the accepted field if the central source/advection/source interval shortens,
or if the final molecular state fails its Mach or acoustic CFL guard. Molecular
heat and viscous work remain in the original transport operator; the combined
report distinguishes body impulse/work from wall and advective exchange. The
caller reuses the same planar forcing scratch, with no added retry arrays.
Regressions cover braking with variable calorics and physical wall loads,
force-driven common-duration retries/replay, final-checkpoint rollback including
prior reports/loads, zero-source bit equality and zero-transport recovery of the
forced Euler owner.

`MixtureEuler2d::step_with_forcing` admits the same uniform or per-cell
`BodyForce2d` and the same conservative kick between its original SSP RK2
transaction boundaries. Density is the compensated sum of each current cell's
canonical species masses, including after advection; no independent density
or frozen-composition surrogate is introduced. Species and their chemical
energy references are unchanged by the kicks. The mixture report separates
source impulse/signed work from wall and reservoir exchange, and subtracts
both when checking referenced and absolute chemical-energy conservation.
The original composition-dependent calorics, HLLC/HLL closure, reconstruction
and atomic retry owner remain in use. Mach and acoustic CFL admission cover
accepted, kicked and RK cell states, with exactly zero forcing following the
native path. Its `mixture::ForcingWorkspace2d` adds one flat `(species+3)`-double
row per cell under both the explicit additional-byte budget and the existing
total mixture cap; source sweeps and retries allocate no arrays. Regressions
cover distinct composition/density, variable-Cp braking, formation-reference
invariance, current-stage density across a compatible material contact,
force-driven retries, zero-source identity and complete cancellation/replay.
This entry drives the inviscid mixture; planar reactive/molecular mixture
forcing and hydrostatic well balancing remain outside this extension.

`Euler3d::step_with_forcing` composes source half step, original SSP RK2
advection and source half step inside the existing retry/publication owner.
`BodyForce3d` distinguishes acceleration [m/s²] from force density [N/m³],
with uniform or borrowed x-fast per-cell fields. Acceleration uses the density
of each current source stage. The authored field is held over one attempted
interval; callers with time-varying forcing must provide each new accepted
interval's field explicitly.

Each kick changes only momentum and total energy. For represented momentum
increment `delta_p`, its work is
`delta_p dot ((p_old/rho + p_new/rho)/2)`. This preserves source-only internal
energy to floating-point accuracy and prevents a momentum increment lost to
rounding from becoming spurious heat. Mass and composition are untouched.
`ForceChange3d` reports extensive impulse and mechanical work separately from
the original per-face boundary exchange; conservation residuals subtract both.
There is no additional sensible-heat source or subtraction of formation energy.

An explicit positive finite `ForceControl3d::maximum_mach` admits the accepted,
kicked and native RK cell states. Original face CFL checks remain at both
advective stages; the final source endpoint also passes a cell acoustic CFL
check. The returned maximum CFL includes that endpoint. A failed kick,
caloric recovery, Mach/CFL check or cancellation retains the accepted field,
clock, loads and report. Every retry recopies the accepted state. An exactly
zero force follows the ordinary step after validating its authored inputs.

`MolecularEuler3d::step_with_forcing` places this same forced Euler interval
between its original molecular halves. A shortened central interval or refused
final molecular endpoint retries the whole composition from the accepted
state. Current molecular/SGS stress, heat, physical wall work and the final
Mach/acoustic checks keep their original owners and a common accepted duration.
Its combined report carries both actual boundary exchange and body work.

The added frozen-gas `ForcingWorkspace3d` stores one five-component conserved
state per cell (40 nominal bytes/cell), with checked required bytes and actual
capacity admission. It is reusable across retries; the added source sweeps
allocate no arrays and poll cancellation. This supplies driven compressible
flows and channel forcing. Gravity remains an explicit split source, without
a hydrostatic well-balanced or turbulent-channel accuracy claim. Focused
native regressions accompany the implementation. The five new Euler tests
and two new coupled tests pass in the native debug build, including variable-Cp
source work, force-driven stability retries, active dynamic channel SGS,
physical wall energy and final-checkpoint rollback.

### Local Coriolis sources on the shared compressible momentum operator

`BodyForce2d` and `BodyForce3d` also admit `Coriolis` (uniform effective
acceleration) and `CoriolisField` (borrowed cellwise effective acceleration).
Both carry explicit angular velocity [rad/s] in the velocity axes. The planar
owner refuses in-plane spin; spatial owners admit arbitrary axes. These are
local rotating-flow sources. Centrifugal, Euler and origin acceleration must
be supplied explicitly in the non-Coriolis field when applicable; the existing
pressure flux remains the sole owner of pressure gradients.

The common source kick integrates `dm/dt = rho*a - 2*omega cross m` with an
analytic skew midpoint inverse. It solves for the impulse directly, avoids
per-cell allocation/trigonometry/iteration, and uses each source stage's actual
density. Pure Coriolis changes momentum direction with exactly zero reported
source work and bit-preserved conservative energy. The speed/axial-velocity
identity holds to floating-point roundoff. With effective acceleration, the
represented kinetic change supplies signed mechanical work using the original
source convention. Neither case adds a sensible-heat or chemical source.

The existing planar/spatial gas, mixture, reactive and molecular forcing
consumers inherit this shared source. Native retries, cancellation, workspace,
Mach and advective CFL admission still govern the complete interval. Zero spin
uses the original acceleration arithmetic. The local midpoint phase is second
order rather than exact; resolve the inertial period. This is not global
spherical geometry or a well-balanced pressure/Coriolis discretization.

Focused regression cases accompany planar geostrophic balance and zero-spin
identity, a repeated arbitrary-axis periodic gas rotation with an independent
Rodrigues phase oracle, spatial drive/work and density coupling, mixture
species/caloric-reference preservation, and invalid source admission. Native
Rust execution of these added cases was unavailable in the authoring session.

### Spatial Coriolis coefficients and latitude-dependent local flow

`BodyForce::CoriolisSpatial` borrows angular-velocity and effective-acceleration
slices. Each independently contains one uniform value or one value per original
x-fast cell. Uniform gravity therefore requires no duplicate cell-sized array.
Both shapes and every coefficient are admitted before a candidate is advanced;
cancellation is polled while inspecting fields. The same source descriptor is
held across both half steps and every retry, with no coefficient allocation.

Each cell uses its own rotation vector and current stage density in the original
analytic midpoint source. In two dimensions only the normal component is legal.
For a traditional beta-plane source, callers sample
`Omega_z(y) = (f0 + beta*y)/2`, where f0 has units s^-1 and beta has units
m^-1 s^-1. Signed and zero local coefficients are supported. Three-dimensional
owners also admit changing local axes. This describes a Cartesian local source;
it neither supplies spherical metric terms nor infers a deforming frame.
Authored coefficient fields must respect the physical domain and its boundaries.

Pure rotation still supplies zero source work and preserves conservative energy
in the source substep. The inherited pressure/transport discretization, timestep
phase accuracy and Mach/CFL controls are unchanged. Uniform broadcast and
cellwise representations produce identical uniform-source arithmetic; the
zero-spin field uses the original acceleration path. Regression definitions
cover hemispheric sign/equatorial rest, heterogeneous density/axis/drive,
periodic conservative transport, allocation reuse, full-field admission and
atomic cancellation/retry.

### Exact held Coriolis source flow

`BodyForce::CoriolisExact` integrates the same physical source
`dm/dt = rho*a - 2*Omega cross m` with its exact held-coefficient affine flow.
Its borrowed angular-velocity and acceleration slices independently broadcast
one vector or supply one per original cell. It uses the current source-stage
density and the existing planar normal-axis/spatial arbitrary-axis admission.
All existing compressible forcing consumers can select it; the midpoint
variants retain their original arithmetic.

The source acts on the net initial impulse with the analytic rotation integral.
This preserves represented geostrophic rest, resolves the inertial phase
`2*|Omega|*dt` to roundoff, and integrates both parallel and transverse drive.
Small-angle series and half-angle trigonometry avoid subtracting nearly equal
cosines. Each cell needs at most one `sin_cos`, fixed stack storage, and no
source subcycling. Zero spin retains the original acceleration path. Pure
rotation reports zero work and retains conservative-energy bits; a driven
source adds only its represented kinetic-energy change, without a thermal or
chemical source. Unrepresentable states and coefficients are refused.

The complete transport/source composition remains second order and retains
native Mach/CFL checks, retries, cancellation and workspaces. Exact local
rotation does not remove time-refinement needs for transport or changing
fields, nor add a pressure/Coriolis well-balanced discretization. Regressions
cover physical phase over many turns, speed and axial momentum, analytic
constant-drive work, the affine flow composition law, geostrophic and zero-spin
identities, broadcast/planar admission, and the original Euler trajectory with
variable heat capacity and full impulse/energy accounts.

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

Reflecting faces default to stationary adiabatic slip walls. Both molecular
owners expose `with_wall(side, MolecularWall3d)` using the same wall definition
as evolving-species flow. No-slip sides admit finite tangential wall velocity
and either adiabatic or prescribed-temperature heat transfer. Their velocity
ghosts mirror all components about the wall, supplying the full SGS gradient.
Wall shear and Fourier heat use the physical half-cell distance; stress work
uses the prescribed wall velocity. The original transport law is evaluated at
the cell/wall mean temperature for an isothermal surface. Modeled SGS stress
and heat vanish at impermeable faces by default. An explicitly configured
smooth equilibrium SGS wall may supply adiabatic tangential stress as described
below; other walls require a resolved molecular layer. Normal wall motion, out-of-domain
temperature and assigning a wall to periodic/transmissive sides refuse.
Open ends retain zero molecular flux unless the complete-reservoir exchange
is explicitly enabled as described below. The combined normal wall traction is
signed `p-tau_nn`; it may be tensile and is not clipped to positive pressure.

Five wall regressions cover half-cell traction/heat/work on all six sides,
zero transverse wall derivatives despite a nonuniform interior velocity,
uniform co-moving flow through the coupled owner, actual heated/moving-wall
energy and momentum exchange, and active interior SGS with unchanged molecular
wall flux. These are conservation and stencil checks, not turbulent-channel
validation or a coarse-grid wall function.

Five focused tests check all-three-axis heat diffusion, oblique longitudinal
and shear modes against independent discrete Fourier symbols, exact planar
extrusion, zero-law Euler recovery and source/work refusal, and six-wall total
impulse with late cancellation/replay. Thermal L2 errors for 6³/12³/18³ cells
on an anisotropic physical box are 4.45971174e-4, 1.13831840e-4 and 5.07037363e-5,
giving observed orders 1.97 and 2.00. These retain source-owned Cv and all
velocity components. They do not establish gas-specific molecular calibration,
physical shock thickness, turbulence or experimental validity. Composition
remains spatially uniform and frozen; evolving-species diffusion and reactions
use their separate admitted owners. Boundary-layer accuracy requires adequate
normal resolution and appropriate source coefficients.

### Molecular transport on moving volumetric frozen gas

`spatial::molecular::MolecularEuler3d::step_moving` composes the original
Newtonian/Fourier SSP half step on the old volume, native extensive ALE
advection, and the second molecular half on the actual new volume. The same
source, cell/face buffers, conservative energy flux and boundary quadrature
remain the numerical owners. No second stress model or viscous-heating source
is added. Every complete trial retains one physical duration; shortened
advection or a later source refusal restarts both halves from accepted cells.

The compatible gradient and its stress adjoint now admit the actual normal
wall velocity as well as authored tangential velocity. No-slip, finite Navier
slip and free-slip walls retain their original traces, including normal
viscous recoil and cross-strain coupling. Their full laboratory traction work
enters energy with the same impulse. `MovingMolecularStepReport3d` separately
retains mechanical molecular wall work and nonmechanical wall heat on all six
sides. Fourier heat and the authored Navier friction-heat partition are
integrated directly at each SSP stage, without subtracting large work terms
from a combined energy account.

Cells, geometry, time, source binding, boundary loads and reports publish
atomically. The prepared operator binds to the new accepted geometry only
after complete acceptance; failed moving proposals cannot poison subsequent
stationary operations. Zero-motion comparisons retain the stationary
numerical result, including nonzero thermal and tangential wall transfer.

Six native moving-field regressions cover all-side stress/work/heat for three
slip laws, comoving full momentum, zero-motion parity, actual-volume balances,
cancellation/replay, complete-interval retry and affine viscous-heating
refinement. Sixteen compatible-stress regressions, including three new moving
normal-trace cases, independently close the quadratic dissipation and wall
power identity at corners, singleton axes and variable viscosity.

This entry admits closed reflecting or periodic domains, translation and
single-axis affine deformation. Reservoir profiles, hydrostatics, SGS,
simultaneous multi-axis deformation and moving chemistry remain outside this
composition. Original source-domain and parabolic/acoustic limits still apply;
these numerical checks do not establish general FSI stability or experimental
accuracy.

### Reciprocal molecular stress and heat with finite volumetric pistons

`spatial::molecular::MolecularEuler3d::step_with_inertial_walls` uses the
original simultaneous two-wall midpoint solve with the complete native
molecular/ALE/molecular trial. Every proposal starts from the same accepted
gas at one fixed physical duration. Both actual pressure and molecular
traction impulses accelerate the finite masses, and their conjugate wall
work uses those same midpoint speeds. No endpoint-force surrogate, full-field
clone per proposal, artificial damping or kinetic-energy repair is introduced.

`MolecularInertialWallsStep` retains the complete moving molecular report,
each wall's mechanical balances, and explicit stationary-side and ideal-guide
reactions. Fixed transverse walls can absorb axial momentum through viscosity;
their actual inward impulse is subtracted in the combined axial system gate.
Piston translation guides carry transverse force, and orientation guides carry
the pressure-plus-viscous torque from every accepted face about its plane
centre. Those guide motions are constrained and do no mechanical work.

Directly integrated Fourier/Navier heat remains external thermal exchange.
Both the source's physical and referenced gas-plus-piston energies must pass
the authored energy tolerance after accounting for that heat. Physical energy
uses differences of conserved mass and energy, preserving the original source
reference without subtracting large absolute totals. Each wall's impulse and
mechanical work still pass their independent tolerances. Exterior pressure
and additional centre-applied axial forces remain explicit.

The selected axis has two reflecting ends; transverse directions are periodic
or stationary reflecting walls. All authored tangential wall velocities must
be zero. No-slip, Navier slip, free slip and prescribed wall temperatures retain
their original laws. The thermal reservoir has no evolved finite heat capacity.
Gas, geometry, source binding, clock, loads, reports and both finite masses
publish together only after complete acceptance and final cancellation.

Six native regressions cover exact zero-transport inviscid agreement on every
axis, independently predicted initial longitudinal viscous recoil, isothermal
and Navier heat versus mechanical work, side support and full guide torque,
whole-interval retry/cancellation/replay, and Galilean work/heat covariance.
The original nine frozen molecular regressions also pass. Original moving
transport/source restrictions remain; this guided-plane composition does not
admit piston rotation, contact stops or surface deformation and does not
establish general FSI stability.

### Instantaneous molecular stress on moving volumes

Both frozen and evolving-mixture `MolecularEuler3d` owners expose
`stress_field_moving`. The caller supplies the current `MeshMotion`; the
mixture also receives its original `DissipationControl`. After a finite-piston
interval, current endpoint piston velocities supply this instantaneous trace,
while the interval's midpoint speeds remain appropriate for integrated work.

The view evaluates the same native compatible molecular stress with actual
normal wall motion in prepared buffers. It borrows the existing
`StressField3d` facade and provides original conserved momentum density,
source-decoded primitives and direct face traction. No normal velocity is
inferred from historical motion or stored as hidden state. Existing
`stress_field` keeps its authored stationary-normal-wall semantics.

All cells, geometry, time, accepted loads and reports stay unchanged.
Cancellation/refusal returns no view; reusable source caches may change.
Observation shares the moving transport's closed/periodic/source admission
without inventing a timestep or advancing any phase. Original viscosity,
bulk stress, cross-strain and no-slip/Navier/free-slip traces remain the
constitutive owners, with no independently reconstructed force model.

Eight native regressions cover independent all-side normal/Navier traction,
exact zero-motion snapshot identity, comoving zero stress for all three slip
laws, actual piston endpoint versus midpoint velocities on all axes, and
cancellation preserving complete accepted field and load histories.

### Exploratory subgrid stress and heat in the three-dimensional gas

With `experimental-sgs`, `MolecularTransport3d::with_sgs` and
`MolecularEuler3d::with_sgs` consume the same explicitly authored Smagorinsky,
WALE or Vreman closure as `fs-lbm`. Each SSP stage caches cell gradients for
transverse face traces. Local closures use the full physical face tensor:
normal derivatives are adjacent-cell differences over the face distance,
and transverse derivatives average the two cached cell traces. This matches
the conservative stress stencil and resolves alternating-cell shear invisible
to centered cell derivatives. Smagorinsky/WALE use the declared multiplier
times the cube root of cell volume; Vreman retains all three directional
widths. No additional numerical scratch is allocated. The zero-coefficient
control preserves molecular-only arithmetic and source calls.

Faces add `mu_t = mean(rho)*nu_t(face)` to deviatoric viscosity, retaining the
original molecular bulk viscosity. Heat adds `mu_t*cp(T)/Pr_t`, using the original
caloric source at the face temperature and the caller's positive turbulent
Prandtl number. The existing shared momentum and total-energy flux carries
both SGS stress work and heat, so internal heating is not added a second time.
The total coefficients enter the parabolic bound at every stage and retry.
A caller-authored maximum `mu_t/mu` is a typed refusal, never a silent clamp.
State, time and published boundary impulses remain atomic on failure.

The six focused SGS tests cover exact zero controls through the coupled Euler
trajectory, all three active closures, periodic momentum/total-energy balance,
resolved kinetic-to-internal energy conversion, an independent helical-mode
stress/heat formula, temperature-mode decay, refusal atomicity and bitwise
replay. These are numerical consistency tests, not DNS/experimental turbulence
validation. The local stress/heat model follows Nicoud's 2007 course equations
(24)-(26); no backscatter, SGS kinetic-energy equation,
pressure-dilatation or shock/turbulence closure or evolving-species mixing is
supplied. Reflecting sides retain the molecular wall treatment unless the
equilibrium SGS wall below is explicitly selected. Acceptance renders refuse. This opt-in low-Mach exploration lane does
not change the default solver, incompressible BDM, or reacting-mixture paths.

### Smooth equilibrium wall stress for coarse SGS gas meshes

Both frozen-gas and evolving-mixture `MolecularTransport3d` owners expose
`with_equilibrium_sgs_wall(side, EquilibriumWallConfig)` under `experimental-sgs`.
Configure an adiabatic no-slip physical wall and an explicit local or channel
dynamic SGS operator first. Frozen `MolecularEuler3d` exposes the same selector;
the mixture coupled/reactive owners consume the already configured operator.
No default wall or zero-interior-coefficient behavior changes without this
additional selection. The current original source density, molecular viscosity,
sound speed, wall velocity and physical half-cell distance are used at every
face evaluation, SSP stage and retry. The source law is not evaluated twice.

This integration reuses `fs-lbm::conjugate::turbulence::law_of_the_wall_ratios`
with its Spalding constants `kappa=0.41`, `E=8.6`. The profile satisfies
`y+ = u+ + [exp(kappa*u+) - 1 - kappa*u+ - (kappa*u+)^2/2 -
(kappa*u+)^3/6]/E` and `u+*y+ = rho*|u_t-u_wall|*y/mu`.
Its **secant** `r=y+/u+` supplies `mu_wall=mu*r`, giving the physical traction
`rho*u_tau^2` over the actual half-cell distance and a continuous molecular
limit at rest and within the viscous sublayer. Tangential traction opposes
the wall-relative velocity; bulk viscosity and normal stress retain their
original molecular source values. The existing wall-work flux uses the same
wall velocity and modified traction. No separate viscous heating, SGS kinetic
energy or chemical source is introduced; an adiabatic stationary wall exchanges
momentum without removing total gas energy, while a moving wall transfers its
native mechanical work. Species flux and the original zero wall heat flux
remain unchanged.

The nonlinear traction derivative enters the original parabolic bound.
For `s=dy+/du+`, the largest tangential Jacobian eigenvalue multiplied by the
wall distance is `mu*2*r*s/(r+s)`, rather than the smaller secant viscosity.
The frozen operator uses this shear coefficient in its full stress row-sum
estimate; the mixture owner adds its excess over the already accounted
molecular shear bound. Both retain the original half-cell and three-axis
inverse lengths. Six fixed-size wall controls add no per-face allocation or
numerical heap workspace, and dynamic heat/species fits remain independent.

The explicit controls admit `0 < maximum_wall_distance_plus <= 1000` and
`0 < maximum_mach_number <= 0.3`; current states outside either authored bound
refuse the complete trial without clipping wall stress. The Mach check uses
the full wall-relative velocity. Positive source density and molecular
viscosity are required. The momentum-only selector refuses slip or isothermal
walls; the separate thermal selector below admits isothermal heat. These numerical
checks do not establish physical applicability: the caller must select smooth,
attached, approximately equilibrium flow and an adjacent sample within its
inner layer. Roughness, separation/adverse-pressure-gradient corrections,
transition, nonequilibrium wall history and shock interactions are absent.
This is an engineering wall closure, not a turbulent-channel DNS qualification.

Focused regressions manufacture Spalding samples across the viscous, buffer
and logarithmic layers and differentiate traction independently to check the
stability coefficient. Native frozen and mixture faces cover all normal axes,
both wall signs, local and channel dynamic momentum/heat/species owners,
moving-wall work, species impermeability, no added workspace, enhanced
coarse-cell braking and total-energy conservation. Late wall-domain refusal
preserves accepted fields, loads and mixture reports. The coefficient check
also covers rest, the viscous limit and wall-relative Galilean invariance.

### Equilibrium turbulent heat transfer at prescribed-temperature walls

Both spatial frozen-gas and evolving-mixture `MolecularTransport3d` owners expose
`with_equilibrium_sgs_thermal_wall(side, EquilibriumThermalWallConfig)` under
`experimental-sgs`. Configure an isothermal no-slip wall and an explicit SGS
operator first. Frozen `MolecularEuler3d` exposes the same selector; mixture
coupled and reactive owners retain it through the supplied molecular operator.
The config combines the existing momentum applicability bounds with an explicit
`WallHeatConfig` (turbulent Prandtl number, quadrature tolerance and bounded
refinements). The original adiabatic selector and default wall behavior remain.

The shared `fs-lbm` thermal wall kernel integrates
`dT+/dy+ = Pr / (1 + Pr/Pr_t * nu_t/nu)` on the same continuous Spalding
profile. Its area-specific resistance `R_w` gives the effective half-cell
conductivity `k_wall = y/R_w`. This includes the molecular sublayer and buffer
layer; it does not multiply molecular conductivity by an arbitrary interior
SGS ratio. At rest it returns the exact original conductivity. At `Pr=Pr_t`
it recovers `k_wall/k = y+/u+`, independently of thermal quadrature.

The original source supplies molecular viscosity and conductivity once per
face at the arithmetic fluid/wall temperature. Source Cp is sampled at that
same temperature; evolving mixtures use the current partial densities and
original species molar masses and calorics, excluding absent species. Density
and wall-relative velocity use the same adjacent sample as momentum. No new
material law, frozen-composition surrogate, species leakage, or per-face heap
workspace is introduced. Properties are constant within the wall integral;
the closure assumes small property variations and does not model strong
heating, compressible wall corrections, or viscous heating within that layer.

The full conductance enters the existing signed total-energy wall flux along
with mechanical work at the prescribed wall velocity. The mixture owner also
updates its nonmechanical heat receipt without subtracting mechanical energy.
The explicit rate uses the modeled thermal conductance and the original
volumetric Cv, in addition to the nonlinear wall-traction bound and physical
half-cell/three-axis inverse lengths. No separate volume heating is added.
Zero molecular conductivity, invalid properties, applicability exit, exhausted
thermal quadrature or cancellation refuse the complete trial. The shared
callback adapter polls the enclosing trial during inversion and quadrature;
accepted states, clocks, wall impulses and reports remain transactional.

Authored regressions cover the independent Reynolds-analogy flux and timestep
on all axes/signs; moving, co-moving, heated and cooled walls; local and channel
dynamic SGS; changing mixture composition and temperature-dependent Cp;
complete-step energy and species inventories; and late quadrature refusal.
Runtime execution of these new regressions was blocked by exhausted build
storage and unavailable dependencies in the implementation session.

### Periodic dynamic gas transport

Both frozen-gas molecular owners also expose `with_dynamic_sgs`, which consumes
the original `fs-lbm::sgs3::dynamic::PeriodicDynamicSgs` Favre--Germano fit at
every SSP stage and retry. Its 27-cell test box has physical widths
`[3hx, 3hy, 3hz]`; domain or local box averaging operates on the numerator and
denominator separately. Negative
fits follow the core's explicit zero-backscatter policy; positive coefficient
and eddy/molecular-ratio bounds refuse. Density and momentum-derived physical
velocity feed the filter, and `rho*nu_t` feeds the same conservative stress,
stress-work, source-owned Cp heat and total parabolic estimate described above.
The turbulent Prandtl number remains authored rather than dynamically fitted.

Dynamic admission requires fully periodic, uniform Cartesian spacing and at
least three cells per axis. Independent `hx`, `hy` and `hz` enter the shared
core's physical directional derivatives, with the cube root of cell volume
as its equivalent filter length. Actual face geometry and directional
parabolic bounds remain those of the original finite-volume operator. Equal
spacings retain the previous cubic-grid arithmetic. The caller budgets all added array capacities,
including density/velocity staging, the filter and face-viscosity cache; the
nominal numerical payload is 176 bytes/cell, excluding the existing molecular
and Euler workspaces and fixed-size/allocator metadata. Arrays are reused
without per-stage allocation. Core cancellation uses the original transaction's
cancellation branch. Four focused regressions cover zero-fit trajectory
identity, variable-density stress/heat routing and energy balance, topology/
byte/coefficient refusals, and cancelled coupled-step replay. They do not
establish turbulent DNS agreement, shock interaction, wall behavior, dynamic
heat closure, or physical coefficient calibration.

`with_dynamic_thermal_sgs` on the same frozen-gas owners independently fits
momentum viscosity and temperature diffusivity using
`DynamicThermalSubgridConfig { momentum, heat }`. The shared scalar Germano
fit consumes the same-stage density, velocity and temperature; it does not
derive diffusivity from a prescribed Prandtl number or divide by the momentum
fit. Thus a zero eddy viscosity can coexist with nonzero modeled heat flux.
The scalar fit retains `C_phi` independently of centered cell strain. Interior
faces interpolate density and the nonnegative coefficient separately, then
evaluate `alpha_t = C_phi,face Delta^2 |S_face|` on the same full velocity
gradient as the conservative stress. Thus a resolved face shear still mixes
heat when both adjacent centered derivatives cancel. The existing face flux
adds `rho_face*alpha_t*Cp(T_face)` to conductivity and uses the total in the original
parabolic stability estimate and conservative total-energy update. The
authored heat-to-molecular diffusivity bound is checked as `k_t/k_m` using the
same face material state; active modeled heat with zero molecular conductivity
refuses this finite ratio bound. Zero modeled heat preserves original
molecular arithmetic. The complete additional SGS payload is nominally
256 bytes/cell, including the reusable scalar filter and thermal staging.
Existing local and fixed-Prandtl dynamic constructors retain their semantics.
Both dynamic constructors support the same rectangular cells without extra
per-cell storage. This remains a periodic low-Mach gradient-diffusion
exploration model; it does not close all filtered compressible energy terms
or qualify arbitrary cell aspect ratios or turbulent DNS agreement.

### Dynamic stress and heat in resolved planar gas channels

With `experimental-sgs`, both frozen-gas molecular owners expose
`with_channel_dynamic_thermal_sgs` and
`ChannelDynamicThermalSubgridConfig { momentum, heat }`. Admission requires a
fixed uniform Cartesian grid, at least three cells on every axis, two periodic
directions and one reflecting wall pair. The three constant cell widths may
differ. The caller supplies explicit lower/upper `MolecularWall3d::NoSlip`
conditions: finite tangential wall velocities and adiabatic or source-admitted
isothermal temperatures. Each preparation reads the current owner walls, so
later valid `with_wall` changes participate in the fit; changing an active
channel wall to slip refuses.

The shared channel Favre--Germano core test-filters only the homogeneous
periodic plane with a nine-cell box and averages each fit's numerator and
denominator over that plane. Physical directional derivatives retain all three
cell widths and the current half-cell wall conditions. Adiabatic temperature
uses the physical zero-normal-gradient condition without imposing that
condition on its variable-density Favre-filtered field. Momentum viscosity
and temperature diffusivity are fitted independently; zero eddy viscosity does
not disable heat. Negative fits follow the core's zero-backscatter policy;
positive coefficient and eddy/molecular bounds refuse instead of clipping.

Interior faces use the original conservative stress, stress work and heat
operator. The momentum and scalar coefficients remain independent of centered
cell strain; interpolated coefficients and density act on the full face strain
to supply `rho_face*C_s,face^2*Delta^2*|S_face|` and
`rho_face*C_phi,face*Delta^2*|S_face|*Cp(T_face)` from the original caloric source.
The original molecular bulk viscosity is retained. Modeled stress and heat
vanish at the physical wall by default; the explicit equilibrium wall selector
can replace adiabatic tangential traction while retaining the same wall-work
owner. Other walls keep molecular half-cell traction and Fourier heat. Fits and total
parabolic estimates are refreshed at every original SSP stage and retry,
including both molecular halves of `MolecularEuler3d`. The existing complete
interval transaction owns accepted cells, time, reports and boundary impulses;
no separate heating source or acceptance path is introduced.

The caller's additional SGS capacity budget covers all filter, primitive
staging and coefficient arrays: nominally 256 bytes/cell, excluding existing
molecular/Euler workspaces, fixed metadata and allocator overhead. These arrays
are reused without per-stage allocation. Periodic constructors keep their
existing geometry and numerical arithmetic. This is a low-Mach exploration
model requiring either adequate wall resolution or the explicitly admitted
equilibrium stress treatment; it supplies no curved/open channel treatment, SGS kinetic-energy
equation, compressible pressure-dilatation or shock closure. Acceptance renders
refuse, and no turbulent DNS or arbitrary-aspect-ratio accuracy is established.

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

### Mechanical forcing of reacting planar gas

`mixture::reactive::ReactingEuler2d` accepts `step_with_forcing` and
`step_with_dissipation_and_forcing` using the original `BodyForce2d`,
`ForceControl2d` and prepared `ForcingWorkspace2d`. The same conservative
K/2-A-K/2 map replaces only the central advection phase of C/2-A-C/2 or
D/2-C/2-A-C/2-D/2. Every kick uses its current compensated species mass and
represented momentum increment; signed mechanical work never becomes a
second chemical or heat source. Original calorics, reactions, molecular
species diffusion, viscous stress, reservoirs and face-load owners remain
authoritative.

The complete report exposes `source_change: ForceChange2d`. Momentum,
referenced-energy and physical-energy balances subtract that accepted source
separately from boundary transfer and chemical changes. The final physical
state after chemistry and molecular heating must also satisfy the authored
Mach bound and complete-interval acoustic CFL. Any shortening repeats every
phase from accepted state, and rejected trials contribute no source or loads.
Cells, clock, both owners' reports and boundary impulses publish together.

The supplied force workspace is reused without per-step allocation and is
included in the complete numerical byte cap even for an exactly zero field.
After input/workspace admission, zero forcing retains the unforced numerical
route. Four regression groups cover signed work with unchanged source-owned
chemistry, exothermic retries and zero parity, an accelerating oblique reacting
contact against its analytic phase on three grids, and molecular wall/source
impulse balance with late cancellation and deterministic replay. This remains
prescribed forcing on the stationary Cartesian discretization; no hydrostatic
well-balanced scheme, implicit velocity-dependent force law or general
coupled temporal-order theorem is asserted.

### Finite Navier slip in planar and volumetric molecular gas flow

`MolecularWall::NavierSlip` adds finite wall friction to the original frozen-gas
`MolecularTransport2d` and `MolecularTransport3d` operators and both evolving-
species owners: `mixture::molecular::MolecularTransport2d` and
`spatial::mixture::molecular::MolecularTransport3d`. Their coupled and reactive
steps consume the same wall law. The authored wall velocity is constant on a
flat reflecting side and has zero normal component. Slip length `b` is finite
and nonnegative; the friction-heat fraction `chi` is finite and in `[0,1]`.
The original gas transport source owns shear and bulk viscosity at each actual
cell state for compatible mechanics. Its original cell/wall temperature and
pressure sampling still owns Fourier conductivity. No accommodation coefficient,
mean free path or rarefied-gas applicability is inferred from `b`.

For physical cell-to-wall distance `h = dx/2`, outward sign `s`, and
centered cross derivative `C_t u_n`, the compatible tangential fluid trace
and traction on the fluid satisfy

```text
r   = u_wall-u_cell+s*h*C_t u_n
u_f = u_wall-b*r/(h+b)
t   = mu*r/(h+b).
```

The shared local trace solve couples adjacent off-diagonal gradients through
the complete symmetric Newtonian strain. Normal velocity remains impermeable;
normal traction retains the adjacent cell's full cross divergence required by
the compatible adjoint. The simpler series formula with
`r = u_wall-u_cell` applies when the cross derivative vanishes. Natural
free-slip uses the matching zero-traction trace instead of an even tangential
ghost. Fourier conduction uses `h`, independently of momentum slip.
`b = 0` uses the same compatible no-slip branch, and increasing `b`
approaches the complete natural free-slip operator. Per-stage evaluation adds
no allocations. Conservative face fluxes, boundary impulse integration,
admissibility retries and atomic publication remain with the existing owners.
Bath and explicitly configured equilibrium-SGS domains retain their complete
native wall discretization.

The massless interface dissipates
`D = -t·(u_f-u_wall) = mu*b*|r|²/(h+b)² >= 0` per unit area.
Actual wall mechanical power is `t·u_wall`, while fluid traction work is
`t·u_f = t·u_wall-D`. The caller explicitly sends fraction `chi` of this
friction heat to the fluid and the complement to the prescribed wall. The
inward fluid energy flux is therefore
`t·u_f + chi*D + q_heat = t·u_wall - (1-chi)*D + q_heat`.
The outward boundary-energy receipt includes the complementary wall heat.
Setting `chi = 1` retains all dissipated kinetic energy in the fluid without
adding a second volume heat source. Setting `chi = 0` delivers all interface
friction heat to the wall while resolved fluid viscous heating remains in the
original conservative energy update. The conductive thermal condition and
this authored interface heat partition are separate physical inputs; wall
temperature and conjugate wall heat storage are not solved here.

The evolving-species owners use the same compatible trace and traction/heat
kernel as frozen gas. Their original `PlanarTransport` and
`SpatialTransport` sources admit the impermeable species boundary and Fourier
heat at actual local mole composition and cell/wall temperature. Actual cell
source viscosity supplies compatible stress and its mechanical bound; the
already evaluated source-owned nonmechanical energy remains intact. This
preserves zero wall flux for every species and adds no chemical or viscous
volume heat source. Original conservative boundary integration carries the
authored interface heat partition through both molecular halves of a reactive
step. Zero slip uses the same compatible no-slip path and arithmetic.

The planar evolving-species operator binds four fixed-size side values through
`with_wall`; periodic/open sides and normal wall motion refuse. Its SSP RK2
and stabilized RKL2 paths use the identical source/ghost/face owner, so adaptive
and composed stepping retain the same physical boundary law. Complete held
reservoirs keep their original exterior traces and source-owned heat/species
efforts. Compatible free-slip gradients use the natural zero-traction trace;
the numerical wall trace changes when normal-velocity cross derivatives are
nonzero.

Channel filters and equilibrium SGS wall closures that require no slip still
refuse this condition explicitly, including when reconfiguring an existing
operator. Prescribed wall acceleration and finite-inertia wall bindings remain
limited to their supported no-slip interfaces. The law does not add slip to
structural wall-load couplers.

Focused native wall tests cover all six oriented spatial faces, series-
resistance traction and analytic interface energy, the finite-slip normal-
stress term, exact Couette velocity with physical viscous heating, actual
conservative energy changes under different heat partitions, bit-identical
zero-slip stepping and receipts, parameter refusals, and agreement between
the frozen planar owner and its three-dimensional extrusion. Evolving-species
tests additionally check actual temperature/mole-composition-dependent
transport on all four planar and six spatial sides, zero species leakage,
finite-slip normal stress, mixture Couette heating, bit-identical zero-slip
updates and receipts, reactive species conversion with the complete absolute-
energy balance under both heat partitions, and late cancellation/replay.
Planar regressions exercise the same exact Couette heat partition in SSP and
RKL2, and independent finite-volume Poiseuille curvature and wall tractions.
These manufactured tests do not calibrate a slip length or heat fraction,
establish rarefied-gas validity, or certify arbitrary boundary-layer resolution.

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

Periodic axes and stationary adiabatic impermeable free-slip walls are the
defaults. At a free-slip wall, normal velocity and tangential divergence
contribute to normal viscous traction; species, tangential traction, heat and
wall work are zero. `with_wall` also admits the shared no-slip and finite
Navier-slip conditions described above, with fixed tangential wall velocity,
adiabatic/isothermal conductive data and an explicit interfacial heat fraction.
Compatible wall mechanics uses actual cell source viscosities and the shared
complete stress/trace law, while original wall-temperature sampling owns
conduction. The same traction supplies mechanical work and physical interface
friction heat. Species remain impermeable in every wall condition.
A transmissive Euler end without
the explicit source-owned reservoir below still refuses molecular evolution:
zero-gradient wave continuation supplies no molecular composition/temperature
boundary law. Normal support traction includes signed `p-tau_nn` when composed
with advection.

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

`mixture::molecular::MolecularEuler2d` consumes the original mixture field and
its already-prepared molecular operator without requiring chemistry. It composes
D(dt/2), A(dt), D(dt/2) over one accepted physical duration. `step_with_forcing`
uses the existing held acceleration/force-density transaction in the middle;
the actual evolving species sum remains the source-stage density. Molecular
species enthalpy, Fourier heat and stress work retain their original owners.
There is no separate reaction, heat-release or viscous-heating source.

A shortened central interval or refused phase restarts both molecular halves
from the original accepted field. Forced paths also admit the complete endpoint
Mach and acoustic rate after the second molecular half. The published report
and face loads include all three phases, with body impulse/work separate from
boundary exchange. Absolute chemical-energy residuals apply the original species
references to conservative differences, so source work is subtracted once.
The borrowed Euler metadata describes the middle phase; the molecular operator's
published report describes the second half. Final cancellation or exhausted work
preserves both owners' accepted fields, reports, clock and face loads.

Construction admits the exact original gas, grid, boundary conditions and held
reservoir metadata before preparing the duplicate field. Only immutable owner
access is exposed. The stored total byte budget includes both fields, prepared
molecular buffers, combined boundary buffers, incoming excess capacity and
retained reservoir rows. Supplied forcing scratch counts even for an exactly
zero source; zero forcing otherwise preserves the ordinary path's bits. Buffers
are prepared once and reused across retries, with the same immutable-source,
bounded-stack and provider-allocation exclusions as the standalone operator.

Six focused G0/G1/G3/G4 tests exercise an accelerating advective/diffusive species
wave with analytic velocity, signed work and impulse; complete reservoir species,
momentum and chemical-energy balance; common-duration replay; late cancellation
and work refusal; zero-source/budget behavior; and molecular cooling that violates
the force Mach bound only after the admitted central transport. The manufactured
wave uses identical occupied-species caloric/formation laws to retain uniform
temperature, while the reservoir case uses unequal calorics and references.
These checks do not certify general temporal order, entropy stability or a
calibrated transport model.

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
The ordinary stationary inviscid species step does not itself integrate
chemistry or molecular species diffusion, preserve pressure exactly at moving
unequal-caloric contacts, move geometry, or establish an entropy/invariant-domain
or experimental combustion claim.

### Moving volumes with evolving species

`MixtureEuler3d::step_moving` advances the original canonical species, all three
momenta and referenced energy with the shared `MeshMotion` geometry. Translation
is admitted on every axis and expansion/contraction on one axis. The original
composition reconstruction, source-owned moving HLLC/HLL flux and nonlinear
wall closure remain coupled; no nominal frozen mixture or independent passive
species advection replaces them.

The same extensive SSP RK2 update uses actual stage volumes and face areas,
satisfying the supported affine-volume GCL. Periodic faces share flux bits;
moving walls transfer zero species mass and pair pressure impulse with wall
speed for work. MUSCL wall ghosts reflect about the actual plane speed.
Prepared cell, face and boundary buffers are reused. Cells, geometry, time,
loads and report publish together; cancellation/refusal retains all accepted
owners. Exactly zero motion follows the ordinary step and polling sequence.

Ten native G0/G1/G3/G4 regressions cover zero-motion parity, all-axis GCL,
comoving unequal-caloric composition/shear contacts, piston species/impulse/work,
adiabatic and relative-advection refinement, velocity boosts, source-reference
shifts, and rollback/replay. The smooth translating-grid advection fixture
reduces L1 error from 0.0046171329 on 8³ cells to 0.0012048216 on 16³ cells,
an observed order of 1.938. This is numerical evidence for that fixture.

This entry excludes moving-domain chemistry, molecular transport, contact-energy
projection, prescribed exterior states/profiles and pressure outlets. Existing
source precision, unequal-caloric advecting-contact and engineering-authority
limits remain; no general entropy, positivity or FSI-stability theorem follows.

### Reciprocal finite-piston evolving mixtures

`MixtureEuler3d::step_with_inertial_walls` advances every original species and
both finite guided piston planes on one Cartesian axis. The same native
`solve_coupled` mechanics and `SpatialInertialWallsStep<StepReport>` used by
the frozen field consume the mixture's actual stage-integrated pressure
impulses and work. Local composition, calorics, all momenta and source energy
references remain with the original mixture owner. Every proposal restarts
from the same accepted field at one duration and reuses prepared scratch.

Both referenced and physical gas-plus-wall energies must satisfy the authored
energy budget. Physical gas increments use referenced `delta E` plus the
compensated sum of each actual species-mass change times its original energy
reference, avoiding subtraction of large absolute formation-energy totals.
No species renormalization, nominal frozen mixture or caloric correction is
introduced. Combined axial momentum, each wall's work and actual represented
plane travel are admitted independently. Nonuniform pressure produces the
same explicit piston torque and opposite ideal-guide angular impulse.

Six native G1/G3/G4 regressions independently reconstruct accepted species,
axial momentum and physical gas-plus-wall energy, exercise both reconstructions
and every axis, and check translation, boosts with external work, source
reference shifts, pressure torque and cancellation/refusal/replay. The ten
moving-mixture regressions pass again after sharing the trial/report helpers.

This entry admits sealed, adiabatic inviscid flow with periodic or stationary
reflecting transverse sides. Reactions, molecular transport, leakage, wall
rotation, contact stops, prescribed exteriors and deformable piston surfaces
remain separate. Original source and unequal-caloric-contact limits remain;
no general FSI stability or experimental validation is asserted.

### Exploratory turbulent mixing on evolving-species three-dimensional gas

The mixture transport owner also exposes `with_dynamic_thermal_sgs` with
`MixtureDynamicThermalSubgridConfig`. Momentum and temperature diffusivity
are fitted separately from the same SSP-stage density, velocity and temperature;
species mixing retains its explicitly supplied common Schmidt number. The
source mixture Cp converts `rho_face*C_phi,face*Delta^2*|S_face|` to conductivity.
Fitted coefficients and density are interpolated independently; the full
conservative-face strain preserves mixing when centered cell derivatives cancel.
The physical strain invariant is evaluated once for all dynamic transports,
while each coefficient retains its scale-safe length multiplication. Heat
therefore remains active when modeled stress and species diffusivity vanish,
and source Cv is still refreshed for its contribution to the stability bound.
The original molecular face reports its evaluated conductivity alongside
viscosity, so heat-ratio admission requires no second transport-law query.
Total energy, referenced species enthalpies, Maxwell--Stefan diffusion and the
reactive splitter retain their existing owners. The combined field/operator
budget includes 264 nominal additional numerical bytes/cell. This independent
heat closure has the same periodic, uniform Cartesian and low-Mach
applicability as the frozen-gas dynamic thermal model. The fixed-Schmidt dynamic
constructors pass all three physical cell widths to the shared fit; the
original directional face distances and combined species/heat stability
estimate remain in force. Rectangular cells add no per-cell storage. This
does not establish accuracy at arbitrary aspect ratios, dynamic species or
chemistry closure, or complete compressible SGS-energy qualification.

`with_dynamic_species_sgs` and `MixtureDynamicSpeciesSubgridConfig` also fit
one common species diffusivity from every resolved mass fraction. The fit
minimizes the sum of all species' Germano residuals using a shared variation
scale; independent normalization would incorrectly increase trace-species
weight. Momentum, temperature and species coefficients are independent and
use one preparation of the same SSP-stage Favre fields. No turbulent Schmidt
number or marker species is chosen. The face adds
`j_i = -rho_face*D_t,face*grad(Y_i)`, where
`D_t,face = C_species,face*Delta^2*|S_face|`, with the existing present-species pivot,
exact absent species, referenced partial enthalpies and combined stability
estimate. Species transport and its enthalpy flux remain active at zero eddy
viscosity. Original Maxwell--Stefan differential diffusion is retained.

The species fit has an explicit coefficient bound and a physical admission
bound `maximum_species_diffusivity_m2_s`. Both the fitted cell diffusivities
and the actual conservative-face diffusivity must satisfy this bound; a
centered strain null cannot bypass face admission. The molecular owner is a coupled
Maxwell--Stefan operator, so its combined viscous/thermal/species stability
envelope is not misrepresented as a scalar molecular diffusivity for a ratio.
The complete field/operator capacity budget includes 352 nominal additional
bytes/cell over molecular-only flow, independent of species count; the fit
visits original mass fractions without allocating a transposed species field.
Refusals retain accepted fields, reports and boundary exchanges. Applicability
remains fully periodic, uniform Cartesian, low-Mach exploratory mixing; there
is no SGS reaction-rate, wall, shock-interaction or DNS-accuracy qualification.

With `experimental-sgs`, the existing spatial mixture molecular owner admits
`sgs::MixtureSubgridConfig { transport: SubgridConfig, turbulent_schmidt }`
through `MolecularTransport3d::with_sgs`. Every SSP stage/retry evaluates shared
Smagorinsky, WALE or Vreman from the complete Favre-velocity face gradient.
Normal derivatives use adjacent cells; transverse derivatives average the
cached cell traces, transposed from derivative-first to the shared
component-first convention. The local face coefficient is
`mu_t = mean(rho)*nu_t(face)`, so modeled stress, heat and species mixing see
the same resolved face shear, including alternating-cell modes. Vreman
retains directional filter widths. Both Prandtl and Schmidt numbers are
explicitly authored positive finite values.

Original Maxwell--Stefan transport, including its Soret/Dufour terms, remains
unchanged. Each interior face adds deviatoric stress/work, `mu_t*cp(T)/Pr_t`
heat, and `j_i = -(mu_t/Sc_t)*grad(Y_i)` species flux. A present-species pivot
closes modeled mass flux to zero while retaining exactly absent species.
Matching original referenced partial enthalpies enter the same total-energy
flux; their source energy differences avoid cancellation of large formation
offsets. No separate chemical or viscous heating source is added. The original
face receipt now exposes its actual molecular `dynamic_viscosity_pa_s`, so
the eddy/molecular ratio is checked without evaluating a stateful law twice.

Physical walls keep molecular transport with zero modeled wall flux by default. A
prescribed constant no-slip wall uses zero physical tangential derivatives
in its normal traction, even when adjacent-cell divergence is nonzero. This
requires a resolved wall layer unless the smooth equilibrium stress model
above is explicitly selected on an adiabatic side. The existing
`ReactingEuler3d::step_with_dissipation` composes this transport with its
unchanged chemistry; no turbulent reaction-rate or flame-wrinkling closure
is inferred.

The owner adds molecular and modeled parabolic scales, including species/
enthalpy coupling, and checks both SSP states. Two reusable scalar arrays add
16 bytes per cell to the numerical budget, exposed by
`required_workspace_bytes_with_sgs`. Hot loops allocate no numerical scratch;
zero model coefficients preserve molecular arithmetic and source calls.
Source, ratio, work, byte-budget and cancellation refusals retain the existing
atomic field/report/boundary publication. Eight focused tests cover independent
full stress/heat/species/energy formulas, gradient orientation, separate Pr/Sc,
zero controls, species/energy conservation and formation-reference invariance,
physical wall fluxes, small-Schmidt stability and refusal behavior, and actual
reactive split-step mixing. This is an algebraic exploration model, without
SGS kinetic energy, backscatter, pressure dilatation or shock/turbulence closure.

### Independent turbulent transport in resolved mixture channels

With `experimental-sgs`, the same mixture molecular owner exposes
`with_channel_dynamic_species_sgs` and
`MixtureChannelDynamicSpeciesSubgridConfig`. Admission requires a fixed uniform
Cartesian grid with at least three cells per axis, two periodic directions and
one reflecting wall pair. Rectangular cells retain their physical directional
widths. The caller supplies two explicit no-slip walls with finite tangential
velocities and adiabatic or source-admitted isothermal conditions. Every
preparation reads the current wall data; rebinding a channel wall to slip
refuses.

Momentum viscosity, temperature diffusivity and one common species diffusivity
are fitted independently from the same SSP-stage Favre fields. Only the
homogeneous tangential plane is test-filtered. All mass fractions contribute
their Germano contractions with one common variation scale before plane
averaging; no marker species or turbulent Schmidt number is selected.
Impermeability in a coupled Maxwell--Stefan/Soret law does not imply a zero
mass-fraction gradient. Both physical and Favre-filtered species derivatives
therefore use three interior fluid centres at each wall, without inventing a
composition boundary condition. Actual zero species wall flux remains enforced
by the original molecular wall owner.

Interior faces retain original Maxwell--Stefan diffusion, Soret/Dufour terms,
source calorics and molecular bulk viscosity. The modeled common species flux
`j_i = -rho_face*C_species,face*Delta^2*|S_face|*grad(Y_i)` uses the existing
present-species pivot, exactly
absent species and matching source-referenced partial enthalpies in conservative
total energy. Independent heat and species transport remain available when
the momentum fit vanishes or centered cell strain cancels. Density and all
three fitted coefficients are interpolated separately and act on one shared
physical face-strain invariant. The actual face species diffusivity also obeys
its configured physical admission bound. Modeled heat and species flux are zero on physical
walls; an explicit equilibrium wall may supply tangential stress using the
same prescribed-wall work. Other walls retain the molecular half-cell stress
and heat. Every SSP stage and retry refreshes the fits and the
combined viscous/thermal/species-enthalpy stability bound. The existing complete
interval and reactive split owners retain accepted state, reports and boundary
exchanges.

The complete field/operator capacity budget includes 352 nominal additional
numerical bytes/cell over molecular-only flow, independent of species count.
All arrays are reused without per-stage allocation or a transposed species
field. Positive coefficient, molecular stress/heat ratio and physical species
diffusivity bounds refuse rather than clip; negative fits use zero backscatter.
This is a low-Mach exploratory closure with resolved walls by default and an
optional smooth adiabatic equilibrium wall stress, without curved/open geometry, catalytic or turbulent reaction-rate
model, SGS kinetic-energy equation or shock/pressure-dilatation closure.
Acceptance renders refuse; no turbulent DNS or arbitrary-aspect-ratio accuracy
is established.

### Body forcing with evolving species and finite-rate reactions

`MixtureEuler3d::step_with_forcing` consumes the shared `BodyForce3d` and
`ForceControl3d` within its original advection, retry and publication owner.
The same midpoint kinetic-work map surrounds SSP RK2 transport. Source kicks
retain every species density, including exact zeros, and derive current total
density from the original compensated species sum. Acceleration and force
density therefore remain distinct in a heterogeneous mixture. The first
kicked state's caloric admission uses the native spatial primitive sweep;
there is no duplicate source-only EOS recovery. Native RK states and the
final kicked state retain Mach and combined acoustic admission.

`mixture::forcing::ForcingWorkspace3d` reuses one flat original-state row per
cell, requiring `8*(species_count+4)` nominal bytes per cell. The explicit
additional-byte cap, actual capacity, geometry/species layout and original
total mixture memory ceiling are checked. Added source sweeps reuse this
storage and accumulate impulse/work in fixed compensated order.

The existing reactive owner adds `step_with_forcing` and
`step_with_dissipation_and_forcing`. The latter keeps its original
D/2-C/2-A-C/2-D/2 composition with K/2-A-K/2 in the central transport interval.
Both molecular halves still prepare their current molecular and optional SGS
stress, heat and species fluxes. Source, chemistry, molecular and final
Mach/acoustic refusals retry the whole physical interval before publishing
any accepted field, load or report. Its explicit total numerical budget
includes the supplied force workspace.

Both Euler and reactive reports expose `source_change: ForceChange3d`.
Referenced-energy and original physical-energy residuals subtract body work
separately from boundary exchange and chemical changes. Mechanical forcing
adds no species source, formation-energy subtraction or duplicate heating.
An exactly zero force preserves the ordinary numerical route after input and
workspace admission. Six focused regressions cover all force representations,
source-only calorics/reference invariance, heterogeneous density, retries,
rollback and a reacting molecular/SGS channel. All six pass in the native
debug build alongside the shared/frozen forcing tests. This supplies a
conservative drive; it does not calibrate turbulent mixing or reaction rates.

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
ends require complete source-admitted reservoir traces; unbound and
pressure-only sides refuse molecular evolution. Soret/Dufour
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

### Stabilized resolved multicomponent molecular diffusion

`MolecularTransport2d::advance_stabilized` and
`MolecularTransport3d::advance_stabilized` select bounded second-order
Runge–Kutta–Legendre (RKL2) integration on the ORIGINAL resolved molecular
operators. Both use the same requested-duration controller and the original
`fs-gas::multicomponent::dissipation::Rkl2Scheme` coefficient owner shared with
axial and network transport.
The negative-real linear stability extension is `(s*s+s-2)/4`, for 2–256
stages. `StabilizedDissipationControl` separately limits physical macrostep
duration, attempts, halvings and actual spatial-operator evaluations.

Every intermediate state and final endpoint evaluates the original face law
at actual local temperature, pressure and composition. Maxwell–Stefan
counterdiffusion, source-defined Soret/Dufour coupling, species-referenced
enthalpy, Fourier heat, complete viscous stress and stress work stay in the
same conservative face flux. The original Cartesian divergence uses one
shared face for neighboring cells. Each signed boundary transfer follows the
SAME Legendre recurrence, including its negative weights; it is not estimated
from the final flux or from an unrelated positive-weight quadrature.

The existing three field arrays are reused. One anchor-flux array per axis and
one boundary-stage array are retained lazily and included in both the required
stabilized workspace estimate and actual owned capacities. There are no
per-stage field allocations. Source-provider allocations retain their original
scope. `required_stabilized_workspace_bytes` admits the combined field/operator
numerical buffers before allocation, subject to the existing 1 GiB ceiling.

All original caloric, species, source-domain and stage-dependent parabolic
checks remain mandatory. A refused physical trial can be halved within its
declared budget. Work counts include failed prefixes, retry anchor evaluations
and final endpoint qualification. Cancellation or exhausted work leaves the
accepted field, physical clock, boundary loads and published reports unchanged.
A successful requested duration publishes the original physical balances plus
actual operator work, largest accepted stage count and the parabolic number
scaled by the Legendre extension.

The planar path retains periodic axes, reflecting slip walls and the existing
source-owned held gas-reservoir faces. The volumetric path retains
periodic axes, complete held reservoir traces and static prescribed molecular
walls, including the original thermal and no-slip conditions. Its optional SGS law, when selected, is
recomputed through the existing spatial owner at every stage. Finite wall
inertia is explicitly incompatible with this stabilized path; prescribed
wall-acceleration histories continue through their original SSP stage-time and
mechanical-work owner. No accelerated-wall trajectory is silently held fixed.

This opt-in path is an Estimated numerical capability. Its linear stability
extension does not prove nonlinear positivity, entropy decay, unconditional
stability or time accuracy. The caller supplies an accuracy-appropriate maximum
macrostep. Focused regressions use independent multidimensional finite-volume
decay symbols and Legendre stability polynomials, signed boundary/energy
balances and whole-duration work/storage/cancellation refusals. These are
implementation checks, not experimental material validation.

### Coupled nonreacting evolving-species molecular flow in three dimensions

`spatial::mixture::molecular::MolecularEuler3d` owns a complete nonreacting
`D(h/2)-A(h)-D(h/2)` interval around the original `MixtureEuler3d`. It consumes
a prepared `MolecularTransport3d`, preserving its admitted wall conditions and
optional SGS configuration, and shares the planar `MolecularStepControl`.
Held body acceleration or force density enters the original conservative
source/advection owner between the molecular halves; acceleration uses the
current compensated species density at each source stage.

Every complete trial starts from the same accepted field. A shorter central
advective step, final forced Mach/acoustic refusal or recoverable molecular
refusal restarts both molecular halves at one common duration under the
whole-step retry bound. Final force admission includes the state after the
second molecular half, so wall
cooling cannot evade the velocity/sound-speed guard. The final cancellation
poll precedes publication of flow time, both owners' reports and all six faces'
integrated loads. Immutable flow/operator access keeps their source, grid and
boundary binding intact.

The complete report and boundary iterator include both molecular halves and
the advective interval. Held-source impulse and signed work remain separate
from wall momentum, wall work and heat, and enter the total conservation
residual exactly once. Species-referenced and absolute chemical-energy
residuals use the original mixture caloric model. Actual field capacities,
molecular/SGS scratch and complete boundary buffers enter the stored numerical
budget; forcing scratch is admitted even for an identically zero source.
Zero forcing retains the native unforced result bits.

Focused checks cover an accelerating three-axis molecular composition wave,
exactly absent species, negative source work, moving isothermal wall/source
balance with unequal species calorics, whole-duration replay, last-publication
cancellation, storage/work refusal and the final wall-cooling Mach guard.
These checks exercise the existing molecular closure; they do not establish
new SGS calibration, moving geometry or an entropy certificate.

### Molecular transport on moving evolving-species volumes

`spatial::mixture::molecular::MolecularEuler3d::step_moving` executes the
original Maxwell--Stefan/Fourier/Newtonian half step on the old volume,
source-owned extensive ALE advection, and the second molecular half on the
actual new volume. Every species, all three momenta, referenced partial
enthalpies and original calorics remain native. No frozen-composition
surrogate, independent passive-species transport or added viscous-heating
source is introduced.

Moving normal wall traces enter the compatible stress and its conjugate
laboratory work, while authored no-slip, Navier and free-slip laws remain
intact. The complete report retains original species and physical/reference
energy balances and separately integrates mechanical molecular work and
Fourier/Navier wall heat on all six sides at the same SSP stages. Wall species
flow remains zero. Prepared source, diffusion, state and face buffers are
reused under the original complete-owner workspace bound.

A shortened central ALE step or a refused later phase restarts the whole
split from the same accepted field. Geometry/source binding, cells, clock,
boundary loads and reports publish together after final acceptance; failed
proposals leave subsequent stationary and direct molecular evolution valid.
Zero normal mesh motion preserves the original stationary numerical result
and still reports actual prescribed tangential work and heat.

Five native regressions cover exact stationary parity, comoving unequal-caloric
Maxwell--Stefan refinement, Galilean normal-stress/work covariance for all
three slip laws, deformation and source rebinding on every axis, complete
boundary balances and cancellation/replay. Translation and single-axis
affine deformation are admitted on closed reflecting/periodic domains.
Moving chemistry, prescribed reservoirs, retained SGS and simultaneous finite
tangential-wall inertia are not composed by this entry. Source-domain,
species, energy-chart, parabolic and relative acoustic checks remain active.

### Reciprocal molecular mixture transport with finite volumetric pistons

`spatial::mixture::molecular::MolecularEuler3d::step_with_inertial_walls`
connects the full original evolving-species molecular/ALE transaction to the
same finite two-piston mechanics. Actual pressure-plus-viscous impulses and
conjugate work drive both masses. Maxwell--Stefan counterdiffusion, referenced
partial enthalpies, Fourier heat, all three momentum components and original
calorics remain active in every fixed-duration nonlinear proposal.

The mixture aliases of `MolecularInertialWallsStep` and
`MolecularInertialWallsError` retain the original mixture report/error types
by default. The shared mechanical report includes side-wall support impulse,
transverse translation-guide impulses and full pressure/viscous guide torque.
Stationary transverse walls can absorb axial viscous momentum and are counted
as external supports in the combined axial balance.

Directly integrated Fourier/Navier wall heat is separate from piston work.
Both species-referenced and physical gas-plus-wall energies must meet the
authored tolerance. Physical increments add each actual species-mass
difference times its original formation reference to referenced energy,
avoiding subtraction of large absolute chemical totals. Every species remains
impermeable at the piston; no normalization or caloric correction is applied.

All proposals reuse original field/source scratch under the complete-owner
workspace budget. A shortened/refused proposal restarts a complete interval,
and cells, geometry, source binding, time, loads, reports and both pistons
publish together. Closed reflecting/periodic directions, zero authored
tangential wall velocities, and the existing single-axis moving-volume
domain apply. Prescribed wall temperature is an external thermal reservoir,
without a finite wall thermal-capacity model.

Four native regressions independently reconstruct all momentum components,
side and guide reactions, heat and both energy charts on every axis; preserve
the piston response under shifted source formation references; recover the
original zero-transport finite-piston owner; and check refusal/cancellation
with no partial publication or workspace growth. The six original coupled
mixture regressions also pass, including source heat, forcing and reservoirs.
The existing source and unequal-caloric-contact limitations remain.

Finite-reservoir report storage is prepared lazily on the heap only for bound
sides. Empty native mixture fields no longer carry two maximum-species report
banks inline, and stage/retry/clone restoration reuses prepared report capacity.
Actual capacities remain in the original required/retained workspace account.
In the native 64-bit build, the combined molecular owner's inline size falls
from 141,992 to 40,616 bytes (99 KiB less), eliminating the reproduced default
test-stack overflow. Original finite-tank regressions and empty-owner/report
lifecycle guards cover source transfers, publication, cancellation and reuse.

### Accelerating no-slip walls in evolving-species three-dimensional flow

`spatial::mixture::molecular::WallAcceleration3d` authors a constant tangential
acceleration on an existing reflecting no-slip side for one accepted interval.
`MolecularTransport3d::advance_with_wall_acceleration` and
`MolecularEuler3d::step_with_wall_acceleration` use the current accepted wall
velocity as the initial value. The combined
`step_with_forcing_and_wall_acceleration` also retains the native body-source
impulse/work owner. Unlisted walls hold their velocity. Successful calls retain
the accelerated endpoint wall velocity, so successive calls can author a
piecewise-constant acceleration history without reconstructing numerical buffers.

The existing requested-duration controller supplies the exact accepted substep
start to the original SSP RK2 stages. Both wall-gradient ghosts and physical
half-cell shear/work use `u_wall(t)` at the first stage and `u_wall(t+dt)` at
the second. The complete molecular split uses wall-time offsets zero and
`h/2`; advection retains its fixed reflecting geometry. Internal parabolic/
source retries reuse the same substep start, and whole-interval shortening
replays both halves from the original wall history. Conservative total energy
includes the original stage wall work; no extra viscous heating or external
momentum kick is added.

A fixed six-side stack descriptor adds no numerical heap capacity or per-stage
allocation. Valid all-zero acceleration uses the held-wall arithmetic unchanged.
Duplicate sides, nonfinite acceleration, normal acceleration, or a side without
an existing reflecting no-slip condition refuse. Overflowing stage/end velocities
also refuse. Endpoint wall velocities publish at the same successful transaction
boundary as accepted cells, reports and face impulses; work exhaustion, source
failure and final cancellation preserve all accepted owners. Native transport
continues to leave the Euler clock unchanged; the coupled owner advances its
single physical clock.

Focused regressions compare native and coupled evolution with an independent
exact one-cell semidiscrete accelerating-Couette solution on three time rungs:
momentum, integrated wall work and internal heating converge at second order.
They also cover simultaneous three-axis body/wall transfer, shortened-step
replay, cancellation/work-limit rollback, and zero-acceleration state bits.
This time-integration check is not continuum shear-profile accuracy on a
one-cell grid. Geometry, wall temperature and transport law remain unchanged;
normal motion, moving-volume work, finite wall inertia and time-varying thermal
data are outside this addition.

### Accelerating walls coupled to finite-rate chemistry

`spatial::mixture::reactive::ReactingEuler3d` accepts the same
`WallAcceleration3d` history through
`step_with_dissipation_and_wall_acceleration` and
`step_with_dissipation_forcing_and_wall_acceleration`. These extend the
existing `D(h/2)-C(h/2)-A(h)-C(h/2)-D(h/2)` transaction, with held body
forcing inside its original central mechanical step. The original molecular
operator supplies every wall ghost, shear and stage wall-work flux; the
original cell reactor supplies evolving species, referenced chemical energy,
temperature and pressure. No separate wall impulse, reaction heating term or
frozen-composition substitute is introduced.

One immutable wall history is prepared from accepted wall velocities before
the whole-interval retry loop. Molecular SSP stages sample physical wall-time
ranges `[0,h/2]` and `[h/2,h]`; the two chemistry half durations are not
added to that wall clock. Heat-driven acoustic shortening, force admission or
molecular refusal restarts the complete composition and resamples the same
history at the new duration. The final molecular state must still satisfy the
original chemistry temperature interval, including when wall work heats it
after the last chemical solve.

Wall exchange, body impulse/work and chemical reference-energy change remain
separate in the existing complete report. The chemical-plus-kinetic absolute
energy residual includes boundary transfer and body work exactly once.
Endpoint wall velocities publish together with accepted fluid time, both
owners' reports and all face loads after the final cancellation checkpoint.
The added history is bounded stack data; existing numerical capacities and
the complete forcing-workspace budget are retained. Zero wall acceleration
recovers both original held-wall reacting paths without changing state bits.

Focused regressions exercise a reacting accelerating-Couette temporal ladder
against independent momentum, wall-work, heating and reversible-species
solutions; activated exothermic feedback from actual wall work; heat-driven
retry replay; late cancellation and work refusal; final chemical-temperature
refusal; and zero-drive/storage behavior. The Couette oracle is an exact
one-cell semidiscrete time check, not continuum spatial or experimental
validation. Normal wall motion and time-varying thermal boundary data retain their
existing exclusions; finite tangential wall inertia uses the separate owner
below.

### Finite tangential wall inertia in native molecular and reacting flow

`spatial::mixture::molecular::MolecularTransport3d::with_wall_inertia` binds
positive finite complete carrier masses to up to six already configured
reflecting no-slip sides. Each planar carrier has two uniform tangential
velocity components. Its normal coordinate remains fixed; the guide receives
normal molecular traction and the original Euler wall supplies pressure
traction. There is no normal carrier work or moving-volume work. Initial
velocity and the thermal condition come from the existing wall model; mass
and model cannot be rebound after finite inertia is attached.

Every original accepted SSP RK2 substep solves
`2*M*(v_mid-v_start) = J_wall(v_mid)` using the existing
`fs_solver::{NewtonKrylovState, NewtonKrylovConfig}` and fallible numerical
Jacobian actions. Each residual replays the same native two-stage molecular
trial with the finite wall held at its candidate midpoint. The original
source computes full stress, species diffusion, Fourier heat and wall work;
the candidate field retains its original positivity, parabolic and source
domain checks, including configured SGS wall-gradient dependence. No second
PDE, wall traction law or explicit viscous-heating source is introduced.

Only the final accepted physical replay contributes native face quadrature
and swaps the prepared candidate buffer. Rejected probes never advance the
current fluid or private wall endpoint. The endpoint is `2*v_mid-v_start`;
componentwise midpoint and represented-endpoint impulse gates and an
independent kinetic/mechanical-work gate precede acceptance. Each substep
therefore pairs received tangential impulse with actual carrier kinetic
storage through midpoint work. The interval report independently checks
aggregate momentum and kinetic work. Native outward boundary energy minus
that mechanical work is the separately reported thermal-reservoir exchange;
normal impulse is separately reported for the fixed guide. Thermal exchange
is never silently added to carrier kinetic energy.

Standalone `advance`, `MolecularEuler3d` and `ReactingEuler3d` all consume
this same molecular owner. The first molecular half's unpublished wall
endpoint seeds the second half. Advection/body-force work, isochoric
chemistry, referenced and absolute chemical energy, and final temperature
and source-domain checks retain their original owners. Whole-step shortening
restarts the wall and fluid together from their accepted state. Only complete
acceptance publishes the field, endpoint walls, face loads and reports.
Prescribed acceleration may coexist on other walls; nonzero prescribed
acceleration on a finite-inertia side refuses. A zero descriptor preserves
the undriven path exactly.

`WallInertiaControl3d` bounds the cumulative native spatial evaluations and
attempted SSP substeps across every nonlinear probe, both molecular halves,
and all outer retries of one public call. Nonlinear dimension is at most 12,
Newton iterations at most 32, restart at most 12 and linear cycles at most 8.
A conservative 64 KiB numerical-workspace reservation covers these bounded
solver vectors and telemetry; callers add it to their original complete
operator budget. Existing field buffers are reused without adding a full
field copy. Allocator metadata and arbitrary constitutive-provider
allocations retain their existing exclusions.

Midpoint carrier storage is energy balanced but is not an L-stable relaxation
method. In particular, a very light carrier can require a much smaller time
step than the fluid-only parabolic cap; `maximum_step_s` controls the actual
coupled substeps and must resolve that carrier relaxation. A coarse,
energy-balanced oscillation is not evidence of resolved physical recoil.

Validation uses the independent one-cell half-cell shear oracle
`delta(t)=delta(0)*exp[-(2*mu*A/H)*(1/m_fluid+1/m_wall)*t]`,
including heavy and light wall masses, distinct temporal refinement rungs,
reciprocal momentum and heat release. Additional cases exercise two finite
walls in the same solve, isothermal reservoir separation, body forcing,
actual reacting species/exothermic feedback and absolute energy, shortening
replay, last-checkpoint rollback and a budget exhausted during the second
half. These are semidiscrete time and conservation checks; they do not
establish continuum boundary-layer resolution, elastic wall modes, external
carrier actuators, normal motion or time-varying thermal models.

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

## Carreau capillary filling and drainage

`capillary::rheology::CarreauCapillaryTube` advances a meniscus through the
original `LaminarCarreauPipe` law, frozen temperature and explicit laminar
shear-rate support. At each liquid length the actual radial viscosity profile
supplies Q; bounded travel-time quadrature integrates A dx/Q, and safeguarded
Newton uses the actual endpoint derivative. The prepared material pipe's
`with_length` changes geometry without repeating its thermal admission. No
apparent viscosity or timestep-lagged constitutive state replaces the law.

`step_horizontal` retains the fixed-pressure horizontal owner.
`step_inclined` additionally consumes positive constant density and signed
gravity projected against increasing wet length. It resolves the combined
drive `p_applied + p_capillary - rho*g*x` through the original radial law at
every spatial sample. Both endpoints admit the whole declared wet interval:
the maximum magnitude of `p_total/x - rho*g` occurs at an endpoint. Signed
gravity supports rising, downward filling and drainage. A stable hydrostatic
equilibrium uses logarithmic remaining distance for the travel-time integral;
finite steps remain on their original side, and unrepresentable remaining
motion refuses rather than snapping to or crossing equilibrium. Accepted
viscous pressure is integrated in physical distance, keeping heat quadrature
well conditioned after many equilibrium relaxation times. Both integrations
share one explicit profile-evaluation budget and retain positive weights.

Reservoir pressure work, wall wetting-energy change, column gravitational
energy change, positive viscous dissipation, constitutive/quadrature work
defect and travel-time residual remain separate. Gravity energy uses the
actual signed volume transfer and mean column height. The inclined result
retains density, projected gravity, signed mass transfer and endpoint
hydrostatic pressure. Its common meniscus result retains the actual liquid
length and original material profile. Filling and drainage stop at the
declared geometric endpoint, returning only its arrival time; no post-event
throughflow is inferred. Zero projected gravity delegates to the horizontal
owner.

The caller supplies surface tension, equilibrium angle and a strictly positive
minimum wet length. Density and gravity are explicit physical inputs, not new
material-source authority. This is a fully developed, no-slip, quasisteady,
vented-gas reduction with no dry-entry singularity model, inertia, displaced
fluid resistance, variable density, dynamic wetting or thermal evolution.
Success adds no experimental authority. Two consecutive Simpson refinements
supply an error estimate, not a rigorous quadrature certificate. Numerical
root, refinement and total profile-evaluation budgets are explicit; refusal
publishes no partial state. Focused numerical tests cover analytic Washburn
and gravity-dependent rise/drainage, independently inverted nonlinear radial
and travel solutions, positive heat and signed mass/energy, temporal
subdivision, reversal, equilibrium approach from both sides, exact geometric
arrival times, zero-gravity delegation and support/work-budget refusal.

`capillary::rheology::two_fluid::CarreauTwoFluidCapillary` additionally resolves
horizontal and inclined displacement through two immiscible original thermal
Carreau fluids.
Both owners have the same physical tube radius and total length, but retain
independent constitutive laws, frozen temperatures and wall-rate supports.
Surface tension and angle must describe the actual ordered pair and substrate;
this numerical consumer supplies no interface calibration. Positive minimum
column lengths are declared fully developed operating limits. Reaching one
returns its arrival time and leaves the remaining requested duration unused.
There is no invented entrance response when either column disappears.

At every interface position, the two original prescribed-flow inverses resolve
the pressure required by a common volume flux. A bounded scalar Newton solve
matches their summed pressure to the imposed reservoir-plus-capillary drive.
Its positive derivative is the sum of the original differential resistances;
no apparent-viscosity blend or finite-difference material tangent is used.
Exact wall-rate profile queries prepare the common source-supported flow
ceiling. Pressure at that fixed common ceiling is affine in interface position,
so both interval endpoints admit the entire declared wet-length range.
Geometry changes retain each prepared original thermal state.

The moving front inverts its positive travel-time integral. Shared fixed-size
Simpson accumulators integrate time and BOTH actual profile dissipation powers
along the path, with compensated sums and two consecutive refinement agreements.
The resulting separate phase heats are not an endpoint resistance split.
The result retains both actual radial profiles, their flow residuals against
the common transport flux, the series-pressure residual, signed volume and wall
energy, reservoir work, and a separate constitutive/quadrature work defect.
Positive physical heat is never increased to absorb those numerical defects.
The invading inventory gains the signed transferred volume, the displaced
inventory loses it, and the same signed volume passes their respective open
reservoirs. No mass transfer is inferred without explicit density data.

The checked entry point polls before every bounded original pipe solve and
before publication; refusal/cancellation leaves the immutable model unchanged
for an identical retry. Numerical pressure, time, quadrature and aggregate
profile work limits remain explicit. The default entry point uses the same
implementation without an external check. The path adds no per-sample heap
allocation. This owner does not include transient radial inertia,
compressibility, mixing, deposited films, dynamic wetting or thermal evolution.
Using a gas requires an explicitly admitted uniform-property incompressible
approximation and vented pressure reservoirs. Numerical residual/refinement
estimates do not certify material uncertainty or the fully developed regime.

The inclined entry point requires both constant phase densities and signed axial
gravity. At wet length `x`, hydrostatic pressure is
`g*(rho_invading*x + rho_displaced*(L-x))`; reservoir pressures belong to the
actual tube ends. Both column resistances use their original radial profiles
under the remaining affine drive. The source-support check compares absolute
affine drive with the common-flow pressure ceiling at both interval endpoints.
This admits the entire interval because absolute affine drive is convex.

A stable density-contrast equilibrium uses the same logarithmic distance
coordinate as the single-fluid inclined owner. Time is integrated in that
coordinate while each actual phase heat is integrated in physical distance,
keeping both integrands resolved after many relaxation times. Finite steps do
not cross or snap to equilibrium; motion below representable position
resolution refuses. Unstable equilibria, equal densities and signed drainage
use the same bounded pressure and travel solves. Zero gravity preserves the
horizontal result exactly.

The inclined report retains both signed phase inventory mass changes, endpoint
hydrostatic pressure, stored gravitational potential change and net potential
energy transported outward to the reservoirs. Their sum is the gravitational
work; counting only stored energy would omit displaced fluid leaving at the
outlet elevation. Pressure work minus wall-energy change, gravitational work
and the two actual heats leaves the explicit constitutive/quadrature defect.
Density is authored separately and adds no source calibration authority.

Focused inclined regressions recover exponential equilibrium approach and
individual phase heats, independently integrated unequal-column Newtonian and
nonlinear Carreau travel, signed drainage and exact arrival, phase/axis
exchange, equal-density and zero-gravity limits, interval support, unresolved
equilibrium refusal and late cancellation with identical retry.

Focused regressions recover independent two-fluid Newtonian travel and phase
heats, identical nonlinear-fluid translation, independently partitioned mixed
Carreau/Newtonian pressure and travel, phase exchange, reversal, subdivision,
frozen source/temperature retention, exact declared-limit arrival, whole-range
source refusal, numerical work refusal and late cancellation with exact retry.
The original horizontal and inclined Carreau regressions cover the shared
quadrature implementation as well.

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
equivalence. Contact-line friction/hysteresis require the explicit sourced
`with_contact_line` opt-in described below. There is no inferred source data,
thermal evolution or flow beyond an inlet/outlet event. A raw transport
fixture remains unsourced. These bounded scalar bindings allocate no memory
during an accepted step and contain no internal cancellation point.

Focused G0/G1/G3 tests in `tests/material_capillary.rs` exercise manufactured
Washburn/Jurin motion, physical material swaps, kinematic-viscosity conversion,
signed mass/energy and end-event accounting, refusal of source mismatches and
unevolved rheology, and finite displaced-fluid resistance with original-owner
retention. They establish numerical/source integration, not measured wetting
accuracy or arbitrary two-phase applicability.

## Binary liquid coexistence in open capillaries

`capillary::material::miscibility::step_binary_liquid_horizontal` consumes the
two actual source-resolved phases of `ResolvedBinaryLiquidEquilibrium`. A
homogeneous miscible result refuses a liquid-liquid meniscus before any flow
or work evaluation. Both composition-dependent viscosities feed the existing
exact two-fluid Washburn operator. Density and actual component concentrations
provide signed inlet/outlet mass and species transfers, including reversal;
their difference is the net component inventory change inside the tube.
End events retain their actual elapsed time and exchanged volume.
Both phases' exact density and viscosity source requirements must satisfy the
shared static-at-fixed-external-conditions policy. Typed held temperature,
pressure and composition remain admissible; rate/frequency/field dependence or
applicability and ambiguous untyped external axes refuse before flow evaluation.
A sampled coefficient alone is not a Newtonian constitutive law.

Pair tension and the wall angle measured through the selected invading phase
are explicitly authored for the actual coexisting pair/substrate. They have no
material-source authority; no liquid/gas surface or pure-liquid wetting value
is substituted. Zero tension explicitly neglects capillary work. Phase/source
receipts remain on the borrowed equilibrium owner, whose finite inventory is
never mutated: reservoirs maintain frozen compositions and pressure during
the accepted step. Bulk coexistence is assumed unchanged by the capillary
pressure jump. No pressure-dependent phase diagram, interfacial transfer,
diffusion, films, nucleation, inertia or dynamic-angle law is inferred.
The bounded scalar adapter allocates no memory. Synthetic-source tests check
the independent integrated resistance law, work/loss balance, both invasion
choices and flow directions, component/mass closure, end events, unchanged
reservoir data and refusal of a meniscus for a miscible solution.

## Finite-compartment osmosis and selective membrane networks

`osmosis` (the existing `continuum` feature) evolves the original
`fs-material::osmosis::IdealDiluteSolution` inventories across prepared
`MembraneTransport` links. Solvent volume and every independent solute amount
are extensive state; concentrations and osmotic pressure are recalculated at
each trial. Species share one explicit ordering, temperature and incompressible
solvent molar volume. Mechanical pressure has an independent volume response,
admitting prescribed pressure, hydrostatic columns, linear compliance and the
original nonlinear convex bulk materials without counting osmotic pressure twice.

Ideal and nonideal networks, finite concentration films and fast-ion Donnan
swelling cancel a common pressure reference BEFORE adding the individual
chamber responses. Mechanical work uses the exact primitive of those same laws.
Consequently an arbitrary common gauge cannot erase a physically resolved
restoring force. Metamorphic regressions repeat actual compliant osmotic
motion under large common pressure offsets and require identical accepted
inventories, and identical pressure work where reported.

At a proposed water transfer, the original exponentially fitted membrane
coefficients reduce all backward-Euler solute equations analytically to
positive, conservative endpoint inventories. A bounded scalar solve closes
the water transfer. Near-complete donor depletion preserves the small donor
remainder using positive implicit weights. The bracket grows from zero along
the initial driving direction; a distant inadmissible volume-floor state does
not veto an admissible interior root. Scaled secant interpolation stays inside
that bracket; bisection follows an interpolation that failed to halve it.
Acceptance requires the water-residual tolerance or a sign-changing bracket
narrow enough to bound water-transfer error. The latter avoids false refusal
when steep osmotic feedback amplifies floating-point residual noise. General
leaking-solute residuals do not carry a uniqueness claim.

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
activity law, Donnan potential, transient boundary-layer storage or resolved
crossflow fields, nonisothermal heat/entropy balance, mechanical inertia, ionic
reaction chemistry or experimentally validated membrane selectivity. Ordinary
nonselective fluid links use zero reflection; gas solution-diffusion retains
its separate existing physical law and permeability units.

### Original nonlinear material resistance to swelling

`OsmoticPressureVolume::from_bulk_material` retains an original
`DiscreteBulkPoint` from `fs-material`, with its quadratic, log-squared or
Simo--Taylor volumetric energy, original prestress, held natural reference and
explicit convex volume domain. Given reference solvent volume Vw, retained
non-solvent volume Vs and the original total volume ratio J0, the mapping is
V0=(Vw+Vs)/J0 and J(V)=J0+(V-Vw)/V0. Background pressure is separate from the
original elastic pressure. The aqueous owner still supplies osmotic pressure
exactly once.

Pressure, positive endpoint tangent and finite storage changes use the SAME
original material law. The tangent includes both reference-volume chain rules;
it is not the finite-interval secant compliance used by the material's discrete
gradient solver. Evaluations are direct and allocation-free on success, with
no nested pressure inversion. Every trial retains the original material's
domain refusal instead of extrapolating a spring or clipping solvent volume.

All consumers of the common compartment law inherit this response: ideal and
nonideal water networks, retained solid/solution inventories, concentration
films, coupled solution-potential accounting and finite Donnan swelling.
Donnan Newton predictors use the actual trial stiffness. Donnan and transient
film reports integrate exact material storage plus background-pressure work;
a pressure trapezoid would be inaccurate for nonlinear chambers. The coupled
nonideal owner uses that same primitive in its total-potential acceptance gate.
Common pressure references cancel before work accounting, with actual solvent
balance checked separately.

Regressions compare all three original bulk laws against independent pressure,
tangent and energy formulas, including dry volume, original prestrain and held
natural reference. Finite Donnan and retaining-network trajectories distinguish
the nonlinear response from its initial tangent spring; finite films check
exact nonlinear storage work. A nonideal coupled network checks independent
mixing-plus-material storage and identical dynamics under a large common
pressure gauge. Domain refusal preserves inventories and clock.
These tests were added but could not be executed in the resource-exhausted
implementation environment. This homogeneous volumetric restraint supplies
no shell/shear geometry, poroelastic transport, evolving reference, polymer
mixing or experimentally calibrated biological-membrane claim.

## Material-driven spherical cavitation

The `cavitation` feature (included by `continuum`) exposes
`cavitation::RayleighPlessetBubble` and `KellerMiksisBubble`. Both consume the existing material owner's
`CavitationState` and an explicit fixed polytropic gas inventory. Absolute
liquid pressure, including metastable tension, drives spherical radius and wall
velocity through vapor, gas, capillary and viscous stresses. The state and
material definition are separate. No pressure floor, injected nucleus, gas EOS,
minimum-radius rebound, or cavitation sound preset is inferred.

Keller--Miksis adds first-order liquid compressibility and acoustic
back-reaction while retaining the same spherical gas/vapor state. Its pressure
input explicitly includes the analytic time derivative. The derivative of the
complete pressure jump includes gas, surface tension, imposed pressure and
Newtonian viscous stress; the latter contributes `4 mu/(rho c)` to the
acceleration denominator. No finite difference of an opaque drive or extra
authored damping constant is used. The wall-Mach cap remains at most 0.1,
and an additional pressure/bulk-modulus limit admits the constant-density,
constant-sound-speed approximation. This does not authorize shocks or sonic
collapse. The Newtonian specialization follows
[Warnez and Johnsen, equations 1--2](https://pmc.ncbi.nlm.nih.gov/articles/PMC4474959/).

Both equations implement `SphericalBubbleModel` and share one integration
kernel. Their native `BubbleLinearization` exposes the equilibrium stiffness,
viscous/radiative damping, viscous inertia and forcing-derivative coefficients
for fractional-radius response. Nonpositive restoring stiffness refuses a
stable harmonic response. Frequency-domain consumers reuse these coefficients
rather than reconstructing another bubble law.


`cavitation::inception` evaluates the original gas inventory's static curve
`p_eq(R) = p_v + p_g(R) - 2 sigma/R`, its restoring slope and the Blake
fold for positive gas, positive surface tension and `3 kappa > 1`.
The same operations apply to the native RP and Keller--Miksis models because
their static balances agree. A pressure below saturation alone does not imply
loss of a gas nucleus's stable equilibrium. Changing imposed pressure never
reconstructs gas inventory with a fresh equilibrium constructor.

Stable and unstable pressure continuation use a bounded, cancellable bisection
of the retained law, requiring both the radius bracket and pressure residual
targets. The admitted branch is intersected with the original radius window;
a usable root can therefore remain available even when the analytic fold is
outside that window. The public threshold query still refuses an out-of-domain
critical radius. Missing roots, nonfinite derived values and unmet tolerances
return errors without pressure floors, radius clamps or state mutation.

This is static mechanical inception under an explicitly selected polytropic
thermal limit, following Brennen, Cavitation and Bubble Dynamics,
[section 2.5, equations 2.41--2.46](https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm).
It supplies neither homogeneous nucleation nor a frequency-dependent acoustic
threshold, gas diffusion or evolving thermal equilibrium. The nine regression
cases include independent isothermal closed forms, stable/unstable native
response, below-threshold growth, restricted native windows, work/cancellation
bounds, numerical refinement and coordinated physical scaling. All nine passed
within a 54-test native cavitation run using the original numerical modules and
material coefficient definitions through focused dependency facades. This does
not replace full Cargo workspace or DSR validation.

The material-water cavitation example uses this API to obtain both equilibria
below saturation, then crosses the Blake threshold and advances the original
RP model with the same retained gas inventory. Its isothermal gas exponent is
explicitly prescribed; no thermal transfer is inferred. The complete example
executable has not yet been exercised against the full dependency stack.

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
derivatives, leading-order near-field liquid kinetic and surface energy, and instantaneous nonnegative
viscous dissipation. Volume acceleration is a physical acoustic source input,
not an observer pressure. The kinetic quantity is the incompressible reference
integral, not a complete compressible or radiated-energy account.
Rayleigh--Plesset excludes acoustic back-reaction; Keller--Miksis includes its
first-order correction. Both exclude shocks, nonspherical collapse/erosion,
nearby walls, interacting bubbles, gas diffusion, and evolving interfacial temperature.
It uses the small vapor/liquid density-ratio approximation to interface
kinematics: supplied vapor density must satisfy the caller's explicit ratio
limit, which cannot exceed 0.1. Without optional vapor density this remains a
caller-declared negligible-mass assumption. Geometric bubble volume is not
silently interpreted as vapor mass or deposited latent heat. A thermally
limited or strongly compressible successor must supply those equations explicitly.

Focused G0/G1/G3/G4/G5 sources check exact equilibrium, independent Rayleigh
collapse quadrature and first-integral energy, the damped Minnaert limit,
tensile-pressure growth, mechanical power and acoustic-volume identities,
radius/Mach/work refusals, and cancellation with bit-identical retry. Physical
equations follow [Brennen, chapters 2--4](https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm);
the independent collapse oracle uses
[Obreschkow et al., equations 3--8](https://arxiv.org/abs/1205.4202).
These checks establish mathematical-model behavior, not experimental accuracy
or a resolved multiphase CFD/film-cavitation result.

## Fixed gas inventory in volume coordinates

The original `BubbleGas` also exposes `from_reference_volume_m3`,
`reference_volume_m3` and `volume_response`. The same polytropic inventory
supplies pressure, analytic positive compression stiffness -dp_g/dV, and the
mechanical work potential -integral(V0..V) p_g dV. The kappa=1 limit is
logarithmic; exp_m1 preserves the near-isothermal limit. An authored zero
inventory gives exactly zero response, while unrepresentable positive
pressure, stiffness or nonzero work refuses. Four native regressions compare
radius and volume coordinates, independent pressure/work derivatives, the
isothermal limit and representability. This is barotropic mechanical work,
not a new gas caloric model, phase-change energy or thermal closure.

## Nonlinear polydisperse bubble evolution

`BubblePopulation` retains each original Keller--Miksis model and its absolute
number density per fixed reference mixture volume. `EvolvingBubblePopulation`
advances each class through the original radial equation and adaptive kernel
under one prescribed smooth pressure history. It preserves gas inventories,
number weights and one accepted clock. Actual native volumes and their
geometric derivatives supply compensated sums of void fraction, its first and
second time derivatives, and the number-radius moment.

The sum of authored per-class maximum-radius volumes must fit the declared
dilute fraction (at most 0.01) before integration. This whole-interval envelope
does not rely on asynchronous endpoint checks or silently tighten native
radius limits. One shared trial budget bounds all class work. Two preallocated
state vectors support proposal inspection and infallible commit without
per-advance state/model cloning or allocation. Failure, cancellation or dropping
the proposal preserves every accepted class and the common clock.

Ten new native regressions check independent geometric derivatives, original
per-class trajectories, split-bin invariance, proposal replay, shared-work
exhaustion, late-class refusal, cancellation, dilution/clock boundaries and
weighted-source representability.
This fixed-number owner supplies no transport, nucleation, breakup, coalescence,
gas diffusion, thermal phase kinetics or inter-bubble pressure feedback.
Dilution alone does not certify weak collective cloud interaction.

## Anisotropic Stokes–Brinkman permeability (`continuum`)

`FluxSystem::assemble_brinkman_mixed` adds `nu K^-1 (u-u_s)` to the existing
BDM1–P0 velocity operator. The steady, fully saturated rigid-frame model
is `-div(2 nu_eff D(u)) + nu K^-1 (u-u_s) + grad(p/rho) = f/rho`, `div(u)=0`.
At rest, velocity is superficial; geometry is metres, intrinsic permeability is m²,
and `nu=mu/rho` comes from the actual fluid. `BrinkmanCell` retains the
original `fs_material::porous::permeability::IntrinsicPermeability` and an
explicit `mu_eff/mu` shear closure; effective shear viscosity never replaces
molecular viscosity in Darcy drag. The reference for that distinction is
[Sandia Aria, Brinkman Momentum](https://www.sandia.gov/files/sierra/Aria_Users_5_28/theory/gov_eq/brinkman_mom.html).

`BrinkmanCell::with_skeleton_velocity` explicitly prescribes a uniform world
velocity per porous cell. Motion uses the one-domain/high-porosity Brinkman
approximation (phi=1 velocity convention) of [Sherwood (2020), equations
2.4–2.5](https://doi.org/10.1017/jfm.2020.265). This does not equate an arbitrary
finite-porosity fluid-only superficial flux with intrinsic fluid velocity:
that would require the phi factor in that paper's equation 2.3 and a different
flux/storage treatment. Geometry, permeability and prescribed motion are a
frozen steady snapshot; interfaces, porosity and solid dynamics are not evolved.
Nonzero skeleton motion is refused in free-fluid cells. Zero prescribed motion
keeps the original stationary matrix and load exactly.

All in-plane tensor couplings enter the same momentum matrix. Positive
principal rates and world-axis projections compile once per cell; symmetric
cell integrals use the original quadrature without a second global matrix.
The same spectral projections add `nu K^-1 u_s` to the original momentum and
constraint-reaction load rows. The matrix and stress law do not change with
prescribed velocity, and no force is injected after the constrained solve.
The original SIP interface law, strong normal velocity moments, total-traction
loads, exact cell continuity, pressure datum and checked residuals remain in
force. Positive pore drag anchors rigid velocity modes in its own connected
component, allowing entirely traction-driven porous components. Every other
component still needs its own prescribed-velocity boundary.
Free-fluid cells have exactly zero drag and preserve the original Stokes law.

`MaterialFluxSystem::solve_brinkman_mixed` uses the original immutable material
snapshot for viscosity, force/traction division and physical pressure. The
accepted result retains cell laws and reports nonnegative cell/total Darcy
dissipation `integral(mu (u-u_s).K^-1.(u-u_s))` in W/m of out-of-plane depth.
`cell_darcy_exchange` and `darcy_exchange` also retain force ON the skeleton
[N/m], signed fluid drag work and signed reaction work ON the skeleton [W/m].
The two powers sum to minus relative dissipation. A driven skeleton can do
positive work on the fluid while producing heat; its mechanical power is not
mislabelled as dissipation. The primitive `cell_darcy_exchange_per_density`
uses the same velocity reconstruction, tensor and quadrature. A common velocity
translation preserves force and loss, while changing the two subsystem powers
by equal and opposite reaction-force work. These are constitutive accounts,
not a complete energy balance or error certificate for the discretization. It does not
automatically deposit that loss into a thermal model or reinterpret superficial
velocity as pore-scale inertia. Raw pore tensors and viscosity ratios do not
gain measured source authority through this binding.

Admission refuses zero principal permeability (no artificial leak floor),
coupling to unresolved out-of-plane flow, unrepresentable drag/shear
coefficients and all original boundary/solve failures. There is no claim for
zero-shear pure Darcy discretization, pore storage, deformation, slip,
Forchheimer inertia, unsaturated/capillary coupling or thermal feedback. The
steady entry inherits the original synchronous Stokes execution and iteration
budget; no new cancellation or performance claim is made.

`cargo test -p fs-flux --test brinkman` checks G0/G1/G3: a pressure-driven
rotated Darcy slab without artificial velocity pins, affine tensor-flow and
dissipation identities, hyperbolic Brinkman-channel refinement, symmetric
positive drag power, mixed free/porous cells, original-source viscosity versus
density scaling, free-fluid equivalence and explicit refusal boundaries.
Moving-frame checks add the independently known translating anisotropic Darcy
slab, Galilean force/loss invariance, positive fluid work from a driven skeleton,
cellwise motion with exact affine moments and original boundary-reaction
closure, zero-relative-flow and exact stationary limits, and invalid/overflow
refusals. No arbitrary-porosity or moving-geometry claim follows from those tests.

### Native transient Brinkman momentum

`FluxSystem::assemble_brinkman_mixed_bdf1` adds the original fluid mass term
`(u-u_previous)/dt` to the same Stokes–Brinkman equation. The full anisotropic
Darcy resistance and moving-skeleton load are implicit at the endpoint; no
explicit drag substeps or second flow integrator are introduced. Coefficients
compile once per cell, outside quadrature, through the same admission as the
steady entry. Endpoint geometry, viscosity, permeability, skeleton velocities
and mixed boundary data remain fixed within the step. The caller retains the
clock and immutable accepted BDM1 velocity history.

Physical transient inertia anchors every connected velocity component,
including free-fluid components with only prescribed total traction. All
original cell continuity equations, per-component pressure datums, SIP terms,
normal constraints, reaction functionals and checked linear residuals remain
active. Free-fluid cells produce the same matrix and right-hand side as native
mixed BDF1 Stokes. Invalid/nonrepresentable timesteps and invalid history are
refused before load callbacks; no candidate updates caller-owned history.

`MaterialFluxSystem::bdf1_brinkman_mixed_step` uses the existing bound material
snapshot and the same physical result conversion as the steady solve. Density
is constant at both time levels: changing material density requires an explicit
conservative state transfer. Original molecular viscosity sets Darcy drag,
while the authored effective viscosity sets shear stress. Physical force/work,
pressure, cell laws and boundary reactions retain their source and units.
Returned drag powers are endpoint rates, not accumulated heat. Backward Euler
is first order in time and has a separate numerical kinetic-energy loss; that
loss is not added to physical Darcy dissipation or to a thermal model.

This is fixed-mesh linear unsteady Stokes–Brinkman with the same one-domain,
phi=1 velocity convention. It does not supply convective momentum, finite-
porosity storage, pore-scale added mass, moving interfaces or skeleton dynamics.
The native synchronous solver's existing execution limits still apply.
`tests/brinkman/transient.rs` exercises independently solved rotated relaxation
at a timestep beyond an explicit drag stability limit, endpoint skeleton
changes, the native mass energy identity with signed skeleton work, first-order
convergence to exponential relaxation, once-bound material ownership with
heterogeneous porous/free cells and constrained boundary reactions, exact
free-fluid operator equivalence, and history/callback refusal boundaries.

### Native nonlinear Brinkman–Forchheimer momentum

`ForchheimerCell` retains the complete original `BrinkmanCell` and adds an
explicitly authored isotropic beta [1/m]. The physical pore resistance is
`mu K^-1 w + rho beta |w| w`, with `w = u-u_s`; the two contributions are
independent, and there is no factor of 1/2, diagonal permeability substitute,
velocity floor or fitted effective permeability. Positive beta requires an
actual porous cell; beta=0 is admitted in both porous and free-fluid cells.
The coefficient convention matches the existing directional Forchheimer law
in `fs-material::porous::permeability`. That owner's specimen/gas/direction
calibration is not silently promoted into an isotropic two-dimensional material
law. Authored beta alone carries no measured source authority or applicability
support envelope; the caller must establish its physical validity.

`FluxSystem::solve_brinkman_forchheimer_mixed` solves the steady nonlinear
Stokes–Brinkman equation. `bdf1_brinkman_forchheimer_mixed_step` consumes immutable
previous BDM velocity and adds the original implicit fluid mass/history term.
Both use one bounded Newton/backtracking loop and the original dense/sparse
linear solve. Linear subproblems retain both native residual checks. The
consistent quadratic tangent is `beta (|w| I + w tensor w / |w|)` for nonzero
relative speed, and exactly zero at co-motion. The original anisotropic Darcy
terms compile once outside iteration; all-zero beta follows the unchanged
original steady/BDF1 Darcy solve and its diagnostic arithmetic.

The nonlinear owner reassembles `J(u)` and the affine load `J(u)u-F(u)` at every
candidate and accepts only the actual endpoint equations. Its convergence
normalization is fixed: row weights and physical reference loads come from the
zero-velocity operator, including skeleton forcing, original boundaries and
BDF1 history. A zero reference load norm uses the initial actual defect instead;
it is never replaced by the changing artificial Newton affine source. The full
Newton velocity correction uses the fixed larger norm of the projected initial
velocity and first Newton candidate. Both original and fixed-row-scaled defects
must satisfy the requested relative tolerance. The full velocity correction
must satisfy either its relative threshold or the explicitly authored
`ForchheimerControls::absolute_velocity_tolerance_m_s`. This absolute tolerance
is finite and nonnegative and defaults to zero, preserving strict relative-only
stopping. It measures the Euclidean norm of native BDM velocity-moment
corrections [m/s]; the caller selects it for the mesh's DOF normalization and
physical accuracy. It is independent of pressure datum and cannot replace either
force/continuity residual criterion. No implicit velocity scale or field clipping
is introduced. The native result's `velocity_update_norm_m_s` and matching
material accessor expose the actual FULL unrelaxed correction, including when
an explicitly accepted absolute correction exceeds the relative threshold.

The fixed relative reference admits homogeneous relaxation to rest without
dividing by a vanishing quadratic Newton load. An explicit absolute correction
allowance also permits pressure-balanced rest, where the first candidate's
velocity may consist only of roundoff. Such an allowance is never inferred from
large pressure loads; resolvable weak nonlinear flow retains the same explicit
accuracy requirement after a pressure-datum shift. The line search uses
the fixed-row-scaled absolute defect and bounded halvings. Exhaustion,
stagnation and unrepresentable states return refusal without a candidate flow.

Native boundary ownership selects the exact normal moments and component-local
pressure pins; their original RHS values are projected into the initial iterate
and every candidate. No numerical row-shape heuristic selects constraints.
SIP interfaces, total traction, Navier slip, continuity, pressure datums and
constrained momentum reaction rows remain with the original owner. The accepted
reaction report uses the endpoint constitutive force, rather than a previous
frozen tangent. Existing rigid-chart history and spatial-density inertia paths
are preserved by the internal cell-resistance hook; this new public closure
still uses one fixed reference density and fixed geometry.

`ForchheimerSolution` separates Darcy from Forchheimer cell exchanges. In each
channel the skeleton reaction is opposite the fluid force, signed fluid and
skeleton powers sum to minus nonnegative relative loss, and an imposed skeleton
can transfer positive work to the fluid. The inertial relative loss is
`integral(beta |w|^3)` per reference density. Force, tangent and inertial work use
the same original positive degree-four cell quadrature. General varying-slip
norms remain nonpolynomial: this is a discretized closure requiring mesh
convergence, not an exact continuum integration claim.

The matching `MaterialFluxSystem` steady and BDF1 entries preserve the original
once-bound transport snapshot. They apply the same density to inertial force,
pressure, boundary reactions and physical work; molecular viscosity still owns
Darcy drag and the authored effective-viscosity ratio only modifies shear
stress. `MaterialForchheimerSolution::flow()` retains the full accepted nonlinear
flow and original Darcy account; its separate Forchheimer exchanges report
physical force [N/m] and signed power/loss [W/m] per unit out-of-plane depth.
Neither physical loss is automatically deposited as heat or conflated with
skeleton work or backward-Euler numerical energy loss.

Scope remains saturated incompressible unsteady Stokes on frozen geometry with
the parent's one-domain/phi=1 velocity convention. This is not convective
Navier–Stokes, turbulent stress, finite-porosity storage, pore-scale added mass,
poromechanical evolution, moving-interface flow or a nonlinear adjoint. No new
async/cancellation execution guarantee is added. Native regressions in
`tests/brinkman/forchheimer.rs` exercise independent scalar roots, rotated
anisotropic/moving-skeleton force and work, homogeneous relaxation to rest,
exact boundary projection, heterogeneous physical source/reaction closure,
implicit relaxation and BDF1 energy balance, exact zero-beta arithmetic, and
failure without accepted state. The production constitutive leaf has direct
force-difference, symmetry, positive tangent, frame-work, scaling, zero-limit
and representability tests. The production convergence helper separately tests
actual correction norms, explicit absolute thresholds, weak-flow protection,
strict default behavior and representability. Native pressure-rest regressions
cover steady/BDF1 equilibrium with affine offset pressure and a large synthetic
beta, material propagation of the actual correction norm, and weak nonlinear
flow under a large pressure datum. These synthetic tests assert no empirical
coefficient validity.

## Inclined two-fluid capillary displacement

The circular meniscus owner now includes both phase columns' hydrostatic head
and Poiseuille resistance. `step_newtonian_inclined_two_fluid_circular_meniscus`
integrates the constant-property separated travel-time law, with a bounded
monotone inverse and actual inlet/outlet event times. Stable equilibria are not
crossed; a negative density contrast reverses the equilibrium's stability.
Equal densities retain the full reservoir-to-reservoir head. Zero gravity and
the massless/inviscid displaced-phase limit preserve their existing operators.

The result separately carries tube gravitational storage, potential energy
transported through the outlet, and each phase's integrated viscous loss.
Pressure work minus wall-energy change, storage and outlet transport closes
against their sum. Positive quadratic resistance weights avoid cancellation
when viscosities nearly match; no allocated quadrature workspace is used.
The material adapter's `step_two_fluid_inclined` samples the original displaced
transport owner once and uses both original densities for hydrostatics and
signed mass transfer. Its existing horizontal method retains the same result
view; `inclined_discrete` exposes the expanded gravitational account.

Focused G1/G3 tests cover analytic exponential motion, density-matched head,
stable/unstable equilibria, repeated equilibrium states, subdivision, end
events, and independent spatial integration of the two phase losses. This
remains a rigid, quasisteady, Newtonian, frozen-property two-reservoir model:
it does not supply trapped-gas storage, inertia, deposited films, thermal
evolution or a dynamic contact-angle law.

## Capillary filling against a closed thermochemical gas pocket

`capillary::material::sealed` connects the original material-bound tube to
`ClosedGasParcel`, retaining the shared NASA-9 caloric/transport model, species,
mass, source query context and liquid/wetting assignments. The tube geometry
explicitly includes any downstream dead volume and must agree with the actual
gas volume. Inlet pressure is absolute. Gas heat is a signed power integrated
over the requested positive interval; changing gas pressure and temperature
feed back into the meniscus motion.

The bounded implicit-midpoint solve uses the SAME gas midpoint pressure for
liquid motion and gas boundary work, and its original analytic volume tangent
for safeguarded Newton iteration. The receiving immutable candidate retains
the accepted gas state and separate inlet work, wall surface-energy change,
liquid gravity, bulk viscous dissipation, optional contact-line and pinning
losses, signed liquid mass, gas heat and actual gas storage increments. Whole-system energy and mechanical pressure/length
tolerances must pass independently. Viscous dissipation is not silently added
to supplied gas heat or duplicated as a second loss.

Cancellation, original source-domain failures, unrepresentable volume/mass/
heat, exhausted convergence, unsupported downward-gravity conditioning, and
oversized unheated equilibrium crossings refuse before publication. Tube-end
crossing requires a smaller step or an appropriate boundary/event model;
the owner never clamps motion and consumes a full interval's heat at an
invented partial event. Smooth branches use a second-order midpoint time
discretization; nonsmooth wetting transitions need refinement. This is not the
exact finite-time Washburn or exact finite isentropic gas operation.

Focused tests exercise independently manufactured perfect-gas midpoint states,
compression/drainage and supplied heat, shared work balance, source retention
and caloric substitutions. The scalar production solver has an 8/16/32-step
refinement check against an independently integrated continuous adiabatic law.
The liquid remains frozen, incompressible and Newtonian at its original source
point. The gas is spatially uniform, with its original frozen composition;
internal gas pressure gradients/viscous resistance, liquid/gas inertia,
liquid heating, dissolution and phase change are outside this rung. Optional
linear contact-line friction and static hysteresis follow the shared sourced
law below; no general dynamic-wetting closure is inferred. The caller explicitly assigns the gas to the declared
wetting environment; material-name strings do not establish that assignment.

### Polarized membrane faces in finite networks

`OsmoticMembrane::new_polarized` consumes the original shared
`PolarizedMembraneTransport`. At each scalar water trial it obtains effective
bulk donor coefficients, eliminates the implicit solute transfers conservatively,
and reconstructs actual membrane-face concentrations before evaluating osmotic
water flow. Each species fit is evaluated once and reused. The membrane and its
films carry the same solute flux; no hidden film inventory, artificial solute
floor or separate nested water root is introduced. Pure-solvent layouts and
explicitly absent films preserve the original limits. Preparation memory counts
the actual boxed membrane and film arrays under the original network budget.

G1/G3 checks compare finite diffusion with the series-resistance analytical
solution, tiny-step retained-solute water rates with independent film roots,
slow polarized relaxation with the unchanged finite-reservoir equilibrium,
and leaking-mixture conservation, orientation and ideal free-energy decrease.
`cargo run -p fs-flux --example osmotic_filtration` compares finite forward and
pressure-driven reverse osmosis with and without films. It reports actual
bulk inventories and conserved interval-mean water flux, using illustrative
explicit coefficients. Films prescribe quasi-steady mass transfer and negligible
mechanical pressure loss; transient film storage, concentration-dependent
transport, and empirically inferred crossflow coefficients remain separate.


## Finite soluble stocks in selective membrane networks

`OsmoticCompartment::with_solubility` binds the original material
`FiniteSolubilitySolution` to the native well-mixed membrane controller.
Each compartment retains its own solid stock. Only dissolved molecules enter
the existing membrane donor coefficients and osmotic force; total dissolved
plus solid inventory is conserved and used in network balance reports.
At each implicit water trial the original cached pair laws supply finite
capacities. A bounded piecewise-affine component equation resolves the common
solute transfer, including crystal exhaustion or precipitation within the
same step. It adds no iterative chemical solve, rate constant, independent
transport integrator or hot-path allocation. Positive donor remainders use
the native trace-preserving implicit recovery.

All four prepared integration buffers own independent mutable phase storage.
Adaptive error control includes total stocks as well as dissolved amounts:
a reservoir cannot drain invisibly while its saturated concentration stays
fixed. Cloning a compartment or requesting a `solubility_snapshot` prepares
an independent checkpoint, so a caller retaining a material clone cannot
invalidate later buffer rotations. Complete phase inventories and original
pair-identity allocations are included in the workspace budget. Accepted
states remain unchanged on refusal or cancellation.

The extension assumes independent ideal neutral solutes, instantaneous local
equilibrium, a common fixed temperature and incompressible solvent, negligible
solute/solid volume, and the caller-declared species ordering. Phase-enabled
cells must share identical ordered original solubility laws. Ordinary aqueous
cells may share the network under its existing solvent and species convention.
A membrane touching a phase-enabled cell must be well mixed and bare:
quasi-steady or transient films would need their own precipitation inventory.
Donnan ion owners refuse these neutral phase stocks; nonideal diffusion owners
continue to require their original shared virial model. This instantaneous-equilibrium option
does not evolve dissolution heat, electrolyte speciation or crystal kinetics;
the finite-rate option below retains the same original constitutive sources.

Focused numerical cases cover exact phase-boundary limits, strong-transfer
traces, orientation, all affine phase combinations, finite reservoir depletion,
water-driven precipitation and redissolution, independently manufactured
coupled water/solute endpoints, ordinary unsaturated parity, conserved
checkpoints, unsupported mode admission and atomic cancellation.

### Finite-rate dissolution and crystallization on the native timeline

`OsmoticCompartment::with_kinetic_solubility` instead binds the original
`KineticSolubilitySolution` and retains its actual dissolved and crystal stocks.
Original fixed-temperature, fixed-volume phase half-stages surround the existing
implicit membrane sweep. Only dissolved molecules are exchanged in that sweep;
crystals remain in place. The original material boundary then updates solvent
and dissolved inventories without repartitioning or erasing a trace solid.
Equilibrium, kinetic and ordinary aqueous compartments can share one network
when their original laws and species ordering are compatible. Conductances and
explicit growth-surface policies may differ by compartment. The original
opt-in self-similar particle law changes area with actual retained solid stock;
reference stock and film coefficient survive membrane sweeps and checkpoints.

The native whole-step/two-half-step comparison controls phase/transport coupling,
including separate dissolved, total and solid amount scales. The membrane
substage remains first-order implicit Euler, so symmetric phase placement does
not confer second-order temporal accuracy on the complete method. Every trial
uses prepared buffers; copied states retain full phase receipts. No additional
transport integrator or per-sweep heap allocation is introduced. Positive-time
phase evaluations, including discarded trials, have a separate report counter;
zero-time target refreshes during inventory replacement are additional work.

Complete phase changes enter the conservation report directly, preserving small
solid contributions against dominant dissolved stocks. Cancellation and work
refusal publish no partial state or clock. The independent kinetic checkpoint
exposes actual phases and the final half-stage receipt. Finite-storage films,
Donnan exchange, and nonideal/coupled molecular diffusion explicitly refuse
these phase owners until those pathways account for their evolution and energy.

This is still a well-mixed isothermal ideal-neutral model with fixed declared
film coefficients and either fixed area or the original particle-area response.
It supplies no phase heat, particle breakup/agglomeration, spontaneous nucleation
or reactive electrolyte transport. Focused
regressions cover an independent continuous dissolution/permeation oracle,
finite-rate concentration/dilution, seeding, temporal refinement, complete
rollback and retry, trace crystals, mixed equilibrium/kinetic cells and
unsupported-path refusal. Particle cases additionally compare the native coupled
timeline with an independent continuous dissolution/permeation ODE, verify
temporal refinement and decreasing conductance, retain a receiver's inert
crystals, and close the complete solid-plus-dissolved inventory.

## Nonideal finite-volume osmosis

`OsmoticCompartment` retains the shared `AqueousSolution`; original
`IdealDiluteSolution` inputs convert directly. Nonideal membrane endpoints use
that actual model's pressure at each implicit candidate volume and retained
inventory. The original hydraulic mobility then applies mechanical pressure
minus the full osmotic-pressure difference. A membrane touching a nonideal
compartment must be bare and fully retaining (all reflection coefficients one,
all solute permeances zero). Existing ideal leaky and concentration-film laws
are refused for such endpoints because they lack matching nonideal chemical
potential forces. Ideal network behavior and its original selective transport
remain available.

The existing adaptive pair exchange retains positive volumes and individual
solute amounts, evaluates declared concentration support before publication,
and can recover a supported nearby root after refusing a distant trial. Shared
model storage is conservatively counted in every prepared state buffer; trial
thermodynamics and successful network advances allocate no memory. Native
focused fixtures include an independent nonideal equilibrium, cross interactions
at equal ideal osmolarity, long-step conservation, concentration-boundary roots,
unsupported combined-law admission and atomic failure/cancellation.


## Neutral crystals coupled to Maxwell--Stefan diffusion

`OsmoticNetwork::with_kinetic_solution_diffusion` couples finite neutral crystal
stocks to the original stationary-solvent `SolutionMaxwellStefan` transport.
Original material kinetic half-stages surround each native molecular stage.
Only actual dissolved molecules cross faces; crystals remain local and retain
their original fixed-area or evolving-particle laws. The existing full/two-half
adaptive timeline controls the complete coupled temporal error. Existing
explicit, implicit binary and implicit multisolute engines remain available.

Every kinetic compartment must share identical original ordered solubility
laws and pure-solid reference, with the same solvent, temperature and molar
volume. Rate coefficients, particle reference stocks and growth policies may
differ. Ideal phase owners require exactly zero virial interactions. Nonzero
interactions require every original kinetic owner to have explicitly bound
that identical virial model through `with_nonideal_thermodynamics`, so its
conditional saturation and molecular forces share the same potential. An ideal
phase source cannot silently acquire nonideal saturation through its aqueous
view. The original concentration support remains binding. Ordinary
dissolved-only cells may join the same network; isolated reacting cells are
also supported. Original
material-backed molecular faces requalify their endpoint sources after phase
exchange. No saturation coefficient is substituted for molecular diffusivity.

Whole-window conservation sums dissolved AND solid changes. The free-energy
report uses the original material owner's calibrated solid-relative chemical
work, retaining the saturation reference and any explicitly bound virial excess
through phase conversion. An increase beyond the original
absolute physical budget refuses before publishing inventories or time. This
is an isothermal chemical-work account, not dissolution heat or a caloric
temperature trajectory. This fixed-volume constructor supplies no water exchange,
ionic/electrical constraint, instantaneous-equilibrium solid projection,
nucleation or crystal transport. Constant-source explicit and specialized binary
advances reuse prepared storage; implicit multisolute/source-query allocation
costs retain their existing separate accounting.

Focused cases compare coupled trajectories with an independent continuous
three-stock RK4 oracle, check phase conservation and actual osmotic particles,
exercise implicit multisolute transport, and check free-energy decrease and
atomic cancellation. Nonideal cases reuse explicit and implicit multicomponent
stages, with an independently evaluated calibrated virial potential.
`examples/crystal_diffusion.rs` connects independently
selected saturation/particle source models to two finite baths and reports
their concentration, crystal stock, evolving conductance and chemical work.
Its coefficients are explicitly authored demonstration data. `--nonideal`
selects an explicit virial binding while preserving the original pure-pair
saturation source.

### Water movement, molecular diffusion and local crystals

`OsmoticNetwork::with_kinetic_coupled_solution_transport` combines original
finite-rate neutral phase owners, bare fully retaining water membranes and
independent stationary-solvent Maxwell--Stefan links at one fixed temperature.
It preserves the original solvent, species ordering and thermodynamic binding
described above. Original phase half-stages surround the existing
water-then-molecular sweep, and the shared full/two-half controller measures the
complete splitting error. Water publishes through the compartment inventory
handoff, keeping the actual phase-owner volume synchronized before subsequent
dissolution or precipitation. Only dissolved molecules cross molecular links;
particle reference geometry and local crystals remain in their compartment.

The coupled account uses the same original material owner's calibrated
solid-relative chemical work, including nonideal excess and variable-volume
contributions, plus the exact pressure primitive and prescribed load work.
The pure-solid calibration is linear in component moles and adds no pressure
term at fixed amounts; water therefore uses the same original virial osmotic
pressure. Complete solvent and dissolved-plus-solid component balances remain
independent absolute gates. Pressure-driven water removal may increase chemical
free energy while mechanical potential falls; acceptance uses the complete sum.
Late source-domain, work or cancellation refusal publishes neither partial
phase inventories nor elapsed time. This isothermal potential account is not
a heat balance. Ionic precipitation, caloric coupling, nucleation, crystal
transport and solvent-drag solute transport require their own constitutive laws.

Focused regression sources compare water/molecular/particle trajectories with
an independent four-state continuous RK4 solution and temporal refinement,
check pressure-driven crystal growth and load work, retain concentration-reference
invariance, reuse implicit transport, and exercise complete rollback. A
two-solute nonideal case independently evaluates calibrated chemical plus
mechanical potential while water moves and local components dissolve and grow
in opposed directions. `examples/crystal_diffusion.rs --water` selects this
coupled path; combining it with `--nonideal` binds virial activities throughout.

## Nonideal solution diffusion in finite stationary-solvent networks

`osmosis::OsmoticNetwork::with_solution_diffusion` evolves finite solute
inventories using the original material `SolutionMaxwellStefan` owner and its
virial chemical-potential forces. Every compartment and connection retains one
original thermodynamic model, fixed temperature/solvent convention and the same
explicit solute ordering. Each link supplies complete solvent/solute and
solute/solute binary Maxwell–Stefan coefficients and original molar masses;
thermodynamics supplies no kinetic data. Canonical endpoint ordering makes
connection input order and orientation immaterial. Duplicate endpoint pairs
refuse. Exact absent solutes remain part of the original species layout. The
material owner's logarithmic-mean face admits finite transport into pure solvent
without concentration floors or seeded trace mass. A species absent at both
ends has zero transfer and zero withdrawal fraction, and stays exactly absent.
All proposed inventories remain nonnegative; a nonzero withdrawal from an empty
donor still refuses. Singular instantaneous endpoint chemical-potential entropy
is reported as unavailable by the constitutive owner, while the finite network
still checks its actual reference-cancelled free-energy change at every step.

This constructor admits only fixed-solvent diffusion connections, with no
membrane water path. Solvent volumes stay exactly unchanged. Original mechanical
pressure laws remain attached to compartments but do no pressure/volume work.
The material owner performs the sole friction solve and the explicit
molar-to-solvent frame conversion; neither gas EOS nor barycentric zero-mass
flux is substituted for the stationary-solvent constraint.

The existing prepared state buffers and complete-duration adaptive timeline
compare one canonical pair sweep with two half sweeps. Diffusion pair stages
use conservative explicit fluxes and remain first order in time; the local
step-doubling estimate is not a certified global bound or an implicit-stiffness
claim. Typed per-pair withdrawal/passivity refusals can halve the interval.
Original source/domain/friction refusals are hard failures. Nonzero transfers
that disappear in both endpoint inventories return `StepUnderflow`, while one
shared actual molar balance budget governs representable rounding. Every face
evaluation, including rejected probes, enters the reported work count.

Before state and clock publish, compensated actual original-to-final species
amount changes and the original source's reference-cancelled mixing-free-energy
change must fit explicit whole-window absolute budgets. Budgets do not grow
with accepted step count. The common linear concentration-reference term is
cancelled analytically; inventory rounding is checked separately. Only this
fixed-volume mode reports `diffusion_free_energy_change_j: Some(...)`; it does
not claim a caloric, electrochemical or mechanically open free-energy balance.
Constant-data numerical/source workspaces are admitted at preparation;
the explicit and specialized binary paths allocate nothing during successful
constant-data advances. The native multicomponent option below accounts for
its allocating nonlinear work separately. Late failure publishes no partial
inventory/time.

`SolutionDiffusionLink::from_material_sources` uses the same timeline with
`fs_material::osmosis::diffusion::source::MaterialSolutionMaxwellStefan`.
Every binary pair retains its original immutable card, exact property and
selected claim. Source resolution refreshes actual typed temperature and all
solvent/solute amount concentrations at both endpoints and the logarithmic-mean face;
the face value enters the existing shared friction law. An admissible face
does not excuse an unsupported endpoint or joint interpolation gap. No source
is switched, extrapolated or inferred from the thermodynamic model.

Each proposed source-backed pair update is qualified before either inventory
is changed, and all sourced links are qualified again after the complete sweep:
later links can change an earlier link's endpoint. Thus the final accepted
window remains inside every pinned source's actual support. These are hard
source refusals, not requests to relax a source boundary. Constant and sourced
connections can coexist under the same whole-window species/free-energy gates.
`solution_diffusion_links` exposes original cards and pinned IDs for explicit
terminal queries. Runtime source receipts are returned by the material owner
and discarded by the timeline; no last-query cache is claimed.

Source queries allocate through the original material resolver. Each sourced
face flux requires three resolutions per binary pair; candidate and end-of-sweep
qualification repeat those resolutions. This source work is additional to the
reported physical face-evaluation count. Prepared workspace admission includes
the full estimated retained card payload per pair, including unselected claims,
observations and models, plus model, mapping and numerical arrays. Shared cards
may be charged repeatedly. The charge excludes allocator-specific overhead and
transient queries/receipts; it is not an RSS or peak-memory guarantee. The
constant-data allocation-free numerical route remains unchanged.

Focused neutral zero-inventory regressions cover uptake into pure solvent against the
independent dilute two-volume exponential and temporal refinement, multi-edge
propagation with a permanently absent species, exact empty equilibrium, late
cancellation/work rollback and retry, and source-qualified zero concentrations.
The joint water/species owner uses this same molecular sweep: its pure-solvent
case retains changing solvent volumes, conservative species exchange and the
combined mixing/mechanical potential gate. No caloric closure or general
nonideal phase-boundary model is inferred by admitting zero solute inventory.

Focused native tests compare binary relaxation with an independent nonlinear
solvent-frame ODE, tighten temporal error, demonstrate finite uphill ternary
transport with decreasing free energy, permute branched links and concentration
references, and exercise actual source-boundary, tiny-step, work, cancellation,
sealed-face and layout/workspace refusals.
Source-backed tests additionally compare unequal-volume evolution to an
independent concentration-dependent mobility ODE, distinguish frozen mobility,
and cover later-edge support loss, late failure atomicity, unsupported initial
endpoints and retained source workspace admission.

### Optional implicit binary molecular exchange

`with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls, byte_budget)`
selects backward Euler for networks with exactly one original dissolved solute.
The existing constructors retain explicit molecular stages by default. The
option applies to both fixed-solvent and coupled water/species networks and
prepares two reusable original `AqueousSolution` trial buffers within the total
workspace budget. It introduces no alternative thermodynamic or friction law.

Each canonical pair solves the conserved amount equation
`q = dt * area * N(n_left - q, n_right + q)` between zero transfer and the
concentration-equalizing transfer at the actual stage volumes. A bounded linear
estimate starts a safeguarded bisection; every trial samples the original
constant or pinned material-source owner at its actual finite inventories.
Nonlinear concentration-dependent diffusivity is refreshed at each candidate.
Every accepted root satisfies BOTH original rounded endpoint amount equations
against an explicit absolute allowance plus a relative allowance scaled by
the concentration-equalizing transfer, not by a large common inventory.
Overflow-safe multiplication and finite-residual admission prevent a huge
duration from producing a false equilibrium or an `Inf <= Inf` success.

An explicit per-pair root-iteration limit and a separate complete-window
molecular-evaluation limit bound initial samples, nonlinear trials, rejected
steps, and full/two-half sweeps together. Source queries are additional work
under the original source contract. Root exhaustion/stagnation returns
`ImplicitDiffusionConvergence`; work exhaustion, source failure, cancellation,
or unrepresentable changes publish no state or clock. Exact empty equilibrium
remains exact, and nonzero changes erased by rounding refuse.

Implicit stability does not waive the original `maximum_species_fraction`,
nonnegative-inventory, source-support, amount-balance or finite free-energy
gates. Callers can deliberately relax the withdrawal allowance independently
of their temporal and nonlinear tolerances. The original full/two-half timeline
still controls first-order temporal and pair-splitting error. This binary
capability does not claim arbitrary multicomponent implicit transport, a global
error bound, or root uniqueness for arbitrary material diffusivity curves.
Focused tests cover an independent stiff unequal-volume BE limit, a complete
stiff window within a small step budget, nonlinear ODE refinement, original
source refresh, coupled water rollback/potential, exact zeros, late cancellation,
whole-window work, physical withdrawal and extreme-exposure refusal.

### Implicit multisolute Maxwell--Stefan exchange

`with_implicit_multicomponent_diffusion(ImplicitMulticomponentDiffusionControls,
byte_budget)` selects coupled backward-Euler species equations for two or more
declared solutes. It uses the existing `fs-solver` `NewtonKrylovState` and
`FallibleFiniteDifference` owners for iteration, FGMRES and globalization. Every
residual, numerical Jacobian probe and line-search candidate evaluates the
original solvent-frame Maxwell--Stefan owner with all binary friction and
nonideal cross forces; material coefficients refresh at the actual candidate.
Neither independent Fick equations nor frozen source coefficients replace them.

Unknowns are individual conserved transfers. Bounded logit coordinates span
their nonnegative, thermodynamic-concentration-supported inventory intervals.
The initial interior equal-concentration trial permits invasion from an exact
zero endpoint without modifying the original physical state or adding trace
mass. A species absent at both endpoints keeps exact zero amount and transfer.
The original selected material source must admit every probe; unsupported
endpoint/face states or joint interpolation gaps remain hard refusals.

Native residuals use FIXED absolute molar scaling. After each native Newton
checkpoint, a fresh physical replay checks BOTH actual rounded endpoint amount
equations in moles against absolute accuracy plus relative accuracy times the
actual transfer magnitude. Thus neither a large background inventory nor a
transformed convergence flag can hide physical error. Finite integration,
original withdrawal/balance/free-energy gates, whole-window atomic publication
and the existing full/two-half temporal controller remain unchanged. The same
option uses actual post-water volumes in the coupled network. Exhausted/stalled
native nonlinear work is a named hard refusal, not a successful equilibrium.

Every original flux evaluation, including initial states, derivative probes,
discarded trials and checkpoint replays, consumes one complete-window budget.
Source qualification queries retain their additional existing cost. Native
dimension/restart, linear-cycle and outer-attempt caps admit a conservative
512 KiB numerical payload allowance, plus two retained original solution
buffers. This estimate includes bounded staged clones/telemetry and is NOT an
allocated arena, allocator-overhead allowance or peak-process-memory guarantee.
This optional native solve allocates transient vectors even for constant data;
the explicit and specialized binary paths retain their original allocation
behavior. It does not assert arbitrary-source root uniqueness, universal
nonlinear convergence or a certified global temporal error bound.

Focused regressions cover exact finite-concentration multisolute counterdiffusion
at large stiffness, one-sided-zero invasion with an absent spectator, an adaptive
stiff window, nonideal uphill transport against independent 2x2 solvent-frame
equations and ODE refinement, all-pair source refresh, coupled water/potential,
whole-window work, cancellation/retry and extreme-exposure rollback.

### Finite-rate zero-current ionic solution networks

`with_zero_current_diffusion(ZeroCurrentDiffusion, byte_budget)` binds one common
explicit valence/charge-budget assignment to every prepared molecular link
before first advancement. `with_zero_current_solution_transport` prepares the
same network with optional fully retaining water paths and the ionic face from
its FIRST network source query; original joint material support need not admit
the off-neutral componentwise-logarithmic face of the neutral transport mode.
All compartments must be initially neutral within the absolute molar budget.

Positive-area molecular paths must form a forest. Their independently floating
zero-current voltages then admit a potential on each connected component.
Cyclic ionic paths, independent ionic membrane leaks and both existing
independent-transfer implicit choices refuse in this forest mode. General
closed graphs use the separately selected global electrical closure below.
No electrode model is inferred. Closed zero-area links constrain neither
current nor potential. A single shared ionic assignment is charged to prepared
storage, with a bounded fallibly allocated forest check at preparation.

The original explicit conserved-transfer sweep, full/two-half adaptive error
controller and all withdrawal/source/free-energy/balance/rollback gates remain.
Each face evaluates the original material's charge-neutral convex composition,
chemical forces, full binary friction and zero-current electrical reaction.
Both actually rounded updated inventories must retain neutrality before either
is published. Source queries retain their original allocation costs; constant
data allocate nothing during advance. A reported physical face evaluation uses
at most three bounded native friction solves. Coupled retained-water paths use
actual changed volumes and the original total mixing/mechanical potential.
There is no concentration seed, ionic activity or caloric model invention.

Focused regressions cover analytical ambipolar salt flux and diffusion voltage,
divalent nonideal uphill salt transport, independent Na/K/Cl salts on a neutral
intermediate face, original source-coefficient refresh, salt-ODE refinement,
coupled water/ion dynamics, exact empty states, singular ionic fronts,
late cancellation/retry, complete-window work and topology/implicit refusals.

### Implicit charge-neutral binary-salt exchange

`with_implicit_salt_diffusion(ImplicitSaltDiffusionControls, byte_budget)` selects
backward Euler on a prepared zero-current forest with exactly two declared
opposite ions. Their gcd-reduced positive stoichiometry defines one conserved
salt extent: `Delta n_i = nu_i*q`, with `sum z_i*nu_i=0`. Original species order
and valences are retained. Neither a counterion nor salt dissociation is inferred.
Neutral spectators, multiple independent salts and global-current graphs refuse
this particular integration choice.

Every nonlinear probe uses both actual rounded endpoint inventories and the
original zero-current Maxwell--Stefan law. All original source pairs refresh at
the trial's endpoints and neutral face. The existing bounded scalar bracket
owner is shared with neutral one-solute diffusion; it does not copy a second
constitutive or root-solving implementation. An original-law bracket check and
bounded interior expansion admit small declared charge-roundoff offsets without
sampling a fabricated concentration floor.

Convergence requires BOTH rounded backward-Euler amount equations for EACH ion,
each with absolute plus actual-transfer-relative tolerance. A charge projection,
large background inventory, infinite predicted exposure or one ion's small
residual cannot replace those equations. Two original solution buffers are
prepared within the caller's byte budget; constant-data stages allocate no
vectors. All probes, full/two-half comparisons and rejected steps share one
complete-window evaluation limit. The original withdrawal, source, charge,
passivity, time-accuracy and atomic publication gates remain mandatory.

The same stage supports fully retaining water paths at their actual changed
volumes. Focused regressions cover unequal-volume stiff ambipolar exchange,
nonideal mixed/reversed/reduced valences, stiff adaptive evolution, nonlinear
source refresh and source-domain refusal, moving-water coupling, work and late
cancellation rollback, exact empty states and singular ionic-front refusal.

### Globally electroneutral diffusion on closed graphs

`with_electroneutral_diffusion(ionic_policy, ElectroneutralDiffusionControls,
byte_budget)` selects a global floating potential for an unadvanced molecular
network, including cycles and disconnected components.
`with_electroneutral_solution_transport` prepares that closure with optional
fully retaining water paths and the correct ionic face from its FIRST source
query. The policy supplies explicit valences and absolute inventory-charge
accuracy; finite edge current is admitted. This does not add electrodes,
imposed current/voltage ports, Poisson charge storage or caloric heating.

Each original material owner supplies its affine current response
`q = q_chemical - G*(psi_right-psi_left)`, with `psi=F*phi/(RT)` and positive
conductance when carriers exist. Area-weighted charge continuity gives a sparse
graph Laplacian. One gauge per conducting component removes the nullspace;
independent component scaling and diagonal preconditioning feed the ORIGINAL
native checked-CG owner, including its recomputed residual and bounded defect
corrections. No dense graph matrix or replacement Krylov recurrence is introduced.
Exactly carrier-free edges contribute no electrical constraint and require no
seeded ions. Every component, including isolated nodes, retains charge admission.

After the potential solve, each edge performs one original combined-force replay
at its bound physical face and EXACT resolved source matrix. The response and
replay cannot see different callback coefficients. Actual replayed edge currents
must satisfy the explicit absolute charge-rate budget at EVERY node, including
the gauge nodes. Circulating currents are allowed; independently zeroing every
edge current would violate potential compatibility on a general loop.

All species transfers are assembled simultaneously from this common state.
Gross outgoing withdrawals are summed per donor species, both signs of each
transfer are retained, and actual rounded inventories satisfy charge and global
species conservation before publication. Pairwise chemical free energy may rise
while electrical work passes between links; the finite-step gate therefore
checks the GLOBAL original mixing free-energy change. Original material replay
independently checks combined chemical/electrical entropy against friction.
Post-transfer source qualification covers every original face before the
private trial state changes. Retained-water stages and whole-window total
mixing/mechanical potential remain on the original adaptive timeline.

Face preparation/replay counts and native potential iterations have separate
whole-window limits, including discarded adaptive probes. The step report gives
those counts and maximum checked node-current, edge-current and voltage-drop
magnitudes; it does not present an earlier trial potential as a final snapshot.
Prepared admission counts reusable arrays, transient bound responses/resolved
matrices and a conservative native-CG numerical-payload allowance. Native CG and
source queries allocate; original receipt strings and allocator overhead remain
outside this payload estimate. This is not an allocation-free or peak-RSS claim.
Independent-transfer and pairwise-salt implicit options are excluded because
their trial coordinates do not enforce the global node-current equations.

Focused regressions cover an independently reduced circulating-current circuit,
compatible loop potentials, per-species conservation and global energy decrease,
isolated/carrier-free components, complete-window work and late-cancellation
rollback, and explicit storage/integration-mode refusal.

### Accepted-state ionic transport observations

`OsmoticNetwork::ionic_diffusion_snapshot(byte_budget, cancelled)` solves the
original global electrical closure at the CURRENT accepted inventories and
time. It returns dimensionless node potentials, one reference node per
conducting component, actual outgoing node charge rates, and each active
connection's original species fluxes, current, voltage and native
chemical/electrical entropy and friction diagnostics. Original connection
indices and left-to-right orientations are retained; zero-area links are
omitted. Potential differences between disconnected components have no
physical meaning.

The readout performs fresh source qualification and original combined-force
replay, rather than exposing scratch from a prior coarse, fine or nonlinear
probe. Its own face and potential work obey the prepared global controls.
Accepted time, inventories, solvent volumes and advance work counters are
unchanged on success or failure. Cancellation or any material/current refusal
returns no partial observation. Additional retained output arrays have a
separate caller-admitted payload budget; the existing network workspace and
source-receipt/allocator exclusions retain their original scope.

Focused regressions compare the accepted-state readout against direct native
face replay, verify gauges and node-current closure, and cover storage/work
refusal, late cancellation and deterministic subsequent evolution.

### Globally implicit multisalt diffusion

`with_implicit_electroneutral_diffusion(ImplicitElectroneutralDiffusionControls,
byte_budget)` selects simultaneous backward Euler on an already prepared global
electrical network. It supports cycles, multiple independent salts, neutral
co-solutes and disconnected molecular components in the original species order.
The existing specialized binary-salt forest option remains available. This
global mode also composes with fully retaining water paths at their ACTUAL
updated solvent volumes on the original adaptive timeline.

A shared feasible-inventory coordinate owner eliminates one anchor cell per
molecular component and one present charged pivot species, preserving each
component's original species totals and each cell's original signed charge.
One component-wide rational radial factor bounds the full inventory direction
inside the original nonnegative concentration-cap box. This chart is piecewise
differentiable; it does not project fluxes, invent salt identities or impose
zero current on each edge. Interior original inventories are the initial
center. Boundary states use a conserved feasible reference solely as nonlinear
coordinates, never as an accepted physical state. Exactly absent species,
saturated pools, isolated cells and zero-dimensional components remain fixed;
unrepresentable interiors refuse without concentration floors.

The ORIGINAL `FallibleFiniteDifference` and `NewtonKrylovState` owners supply
Jacobian products, FGMRES, line search and nonlinear iteration. Every trial
decodes complete physical inventories and calls the original global electrical
owner, including actual source resolution, gauge-fixed checked CG and bound
combined-force face replay. At each checkpoint, EVERY original node/species
backward-Euler equation must satisfy its explicit absolute/relative molar
budget, including all eliminated anchors, counterions and absent spectators.
Convergence of only the reduced nonlinear system is insufficient.

The final original-law replay supplies gross outgoing withdrawals at the
converged candidate. Publication retains the same rounded inventory-charge,
species-balance, global mixing-free-energy and post-transfer source gates as
the explicit global mode. The enclosing full/two-half comparison still controls
local temporal error, and coupled water stages retain their total mechanical
and mixing-potential account. Source refusals, exhausted nonlinear work and
cancellation publish neither physical state nor time.

The nonlinear face budget counts actual preparations/replays across the complete
requested window, including finite differences, line searches, discarded time
probes and failed prefixes. Native electrical work has its own whole-window
bounds. `maximum_unknowns` explicitly limits independent coordinates. Prepared
storage admits reusable trial solutions, the shared chart and a
dimension-scaled native nonlinear payload allowance; it reserves topology's
possible rank because water can release an initially saturated constraint.
The estimate retains original source-receipt and allocator exclusions and is
not a peak-RSS bound. No dense global Jacobian or replacement numerical solver
is added. Implicit stability is not a universal Newton-convergence or temporal
accuracy guarantee.

Focused regressions cover an exact finite-concentration cyclic decay mode,
stiff adaptive evolution, heterogeneous circulating currents and complete
original-equation replay; coordinate tests cover mixed valences, absent and
neutral-only pools, disconnected/saturated states and extreme probes. Existing
retained-water and source-dependent cyclic fixtures exercise both global
integration modes. Work, storage, unknown-count and cancellation refusals retain
atomic publication. Accepted-state ionic observations use the same native
readout for either mode. This extension adds no electrode reactions, imposed
electrical ports, Poisson charge storage or double-layer model.

## Coupled nonideal water and molecular solution transport

`OsmoticNetwork::with_coupled_solution_transport` composes the existing fully
retaining nonideal water exchanges and original solvent-frame Maxwell–Stefan
links on ONE existing complete-duration adaptive timeline. All compartments
retain the same virial model, fixed temperature, solvent-volume convention and
original species ordering. Water membranes must be bare, with reflection one
and zero solute permeance. Molecular connections remain independent stationary-
solvent pathways; sharing endpoint compartments does not infer a solvent-drag
law or Onsager cross coefficient. Either path set can be empty for its physical
limiting case. Both sets have canonical endpoint order; duplicates within one
set and an oversized combined graph refuse.

Each trial first sweeps the original implicit water owner, then the original
molecular owner. The existing full/two-half comparison controls the combined
first-order Lie-splitting error in BOTH volumes and species inventories. No new
transport solver, EOS, pressure closure, hidden water inventory or ideal leaky
membrane law is substituted. Material-backed pair coefficients refresh at the
actual concentrations produced by water motion, including endpoints and face;
all source links are requalified after water and after the complete molecular
sweep. Pinned-source and physical-domain failures stay hard refusals. Only typed
water root/bracket failures and finite-step withdrawal/passivity failures may
shrink the step.
The optional neutral binary, multisolute or charge-neutral salt backward-Euler
choice changes only the molecular stage on this same timeline; its source samples use the actual post-water
volumes, and its nonlinear work/refusal limits retain whole-window rollback.

The coupled potential is `Phi = F_mix + sum(Psi(V))`, with `dPsi/dV=p(V)` from
the original linear chamber-pressure law. Thus `dPhi/dV=p-Pi`, exactly the force
already consumed by the retaining water law. Each water stage checks the
original inventory free-energy change plus the exact pressure-primitive
change; molecular stages retain their original fixed-volume passivity gate.
Before the complete interval publishes, compensated actual water/species
changes satisfy explicit absolute conservation budgets and actual total
`Delta Phi` satisfies ONE whole-window allowance. The allowance does not grow
with accepted/rejected stage count. Mixing free energy may increase during
pressure-driven separation when the declared pressure loads supply work.

Common linear concentration and pressure reference terms cancel analytically
before potential arithmetic, with inventory residuals checked independently.
Pressure differences also cancel the common reference before evaluating elastic
response. `coupled_potential_change` separates mixing and mechanical changes
and their compensated total; the negative mechanical change is pressure work
delivered to the solution. `diffusion_free_energy_change_j` remains exclusive
to the fixed-volume mode. This is an isothermal potential/work account, not a
caloric-energy, thermal-reservoir, inertia or cavitation model.

All existing publication, bounded-work, source-payload and allocation rules
remain: a failed/cancelled window publishes no volume, amount or clock. The
original explicit/binary constant-data paths advance without allocation; the
optional native multisolute solve has its declared numerical payload allowance,
and material queries allocate their original transient receipts. Focused tests use an independent two-state nonideal
water/species ODE with state-dependent mobility, temporal refinement,
pressure-driven separation, concentration/pressure datum and link-orientation
invariance, water-induced source-boundary failure, cancellation, work and
workspace refusal, and unsupported leaky membranes.

## Finite-storage transient concentration films

`osmosis::transient::TransientOsmoticPair` advances two finite bulk solutions
and explicit one-dimensional film cells. Film solvent volumes are fixed
`area * cell_width`; the same normal water flux crosses every interface while
bulk volumes change conservatively. Film thickness, cells, initial profiles and
per-species diffusivities are explicit independent inputs. Initial film water
and solute are additional physical inventory, never silently taken from bulks.

`TransientFilm::from_cells` additionally accepts nonuniform cell widths and
cell-major species diffusivities. Each interior face combines the two actual
half-cell diffusion resistances in series, preserving slow barrier layers and
the original exponential advection fit. Storage and balance use each cell's
own volume; `film_cell_volume_m3(side, cell)` exposes that volume. The uniform
constructor remains available. Coefficients are positive and fixed in time;
mobile concentration is continuous across internal faces. Optional explicit
linear or finite kinetic sorption below partitions total storage while preserving
the same mobile phase. Focused regressions cover analytical resistance, both flow directions,
high material contrast and nonuniform spatial refinement.

The `step` call accepts one backward-Euler step. Original exponentially fitted donor
coefficients form a positive implicit amount system for each species. Native
`fs_la::tridiagonal::TridiagonalWorkspace` solves the nonsorbing/linear case in
linear work per species and supplies each kinetic Newton correction below;
all storage is prepared once. Terminal half cells use the existing polarization
composition and actual reconstructed membrane-face concentrations drive water
flow. The finite scalar solve is bounded and safeguarded; stiff sign-changing
brackets may control transfer-position error, with the actual residual retained.

Reports include every bulk and film inventory, original linear-equation
residuals, membrane-face concentrations, water-root residual and exact boundary
pressure work for the declared linear chamber laws. Failure, cancellation,
nonphysical amounts, volume depletion or exhausted work publish no state or
clock. Only neutral ideal solutions are admitted: nonideal chemical transport,
ionic fields, resolved crossflow, film hydraulic losses and heat storage are
outside this owner. Root tolerances alone are not temporal error estimates;
`step` callers refine time independently, and all callers still refine the mesh.

The `transient_osmotic_filtration` example matches initial whole-system water
and total solute against a nonsorbing steady-film reference. It exposes startup
mobile/bound storage and face enrichment in layers with R=1,2,6. The reference
uses the true series mobile diffusion resistance but has different storage
physics; it is not the steady limit of the sorbing film. Focused regressions cover analytical diffusion,
multiple permeating species on unequal meshes, steady forward/reverse limits,
time/space refinement, steep roots, inventory checks and atomic refusal.


### Explicit liquid-film surface and interface resistance

`TransientFilm::with_interface_resistances` accepts finite nonnegative
mass-transfer resistance [s/m] in face-major species order, including both
bulk/film and film/membrane ends. Spatial left-to-right face order is retained
on either membrane side. At zero water flux the contact law is
Js=(c_left-c_right)/resistance. With solvent motion, the same original stagnant-
layer advection-diffusion equation is integrated through the added resistance.

Preparation composes contact and actual half-cell resistances before the
existing exponential transport coefficients are evaluated. Stable series
arithmetic handles extreme contrasts and overflowing k*r products. Zero
resistance retains the exact original coefficient and storage path. The
original tridiagonal solve, mobile-only transport, osmotic reconstruction,
linear/kinetic sorption and adaptive rollback consume these prepared faces.

This is an explicit passive boundary-layer limit with no contact storage,
partition coefficient, extra selectivity or hydraulic loss. The same solvent
flux crosses each contact. It does not infer a general moving phase boundary
or a measured interface coefficient. Focused regressions compare independent
series transfer and whole-film profiles in both flow directions, signed Peclet
limits, original zero resistance, and contact plus kinetic adaptive evolution.

### Equilibrium immobile sorption in transient films

`TransientFilm::with_linear_sorption` accepts an explicit cell/species total-to-
mobile ratio R>=1 for instantaneous, unsaturated, immobile linear sorption.
Conserved film amount is n_total=R*V_water*c_mobile. Initial profiles remain
mobile concentrations, so the added bound inventory is explicit; callers
prescribing total concentration supply c_total/R. Prepared R*V capacities
convert total moles to mobile concentration at every diffusion/advection face
and membrane reconstruction. D, actual water volume and hydraulic laws stay
unchanged, avoiding double-counted retardation.

Common mobile concentration and ideal chemical potential define the no-flow
equilibrium; total concentration may jump with R. This is not a change in the
mobile phase's chemical-potential standard state. Existing film/system amounts
and adaptive inventory norms include both pools; separate getters expose mobile
concentration, mobile moles and sorbed moles. R=1 retains the exact original
calculation without an additional capacity array. Prepared stepping remains
allocation-free and additional preparation storage is included in its budget.
The finite kinetic option below additionally represents saturation and delayed
release; swelling and sorption heat require other owners.

Focused authored tests check exact R=1 compatibility, stationary heterogeneous
partition without a false osmotic force, analytical finite-bath retardation,
reduced fixed-volume chemical-potential dissipation, and mobile-only advection
in both directions at high Peclet number. They are not measured sorption or
filtration validation.

### Finite-capacity kinetic sorption in transient films

`TransientFilm::with_langmuir_sorption` accepts an optional independent
`FilmLangmuirSite` for each cell/species: fixed capacity in moles, adsorption
rate in m3/(mol s), desorption rate in 1/s and explicit initial trapped moles.
Mobile concentration drives dn_trapped/dt =
k_ads*c_mobile*(capacity-n_trapped)-k_des*n_trapped. Site capacities do not
change during water-root trials. Aqueous, instantaneous linear and finite
kinetic storage can coexist. Only mobile concentration enters diffusion,
advection, membrane reconstruction and osmotic water transfer.

The shared scalar-site backward-Euler owner eliminates each local trapped
endpoint and supplies its analytic derivative. Bounded Newton/backtracking
then uses the original tridiagonal workspaces for the remaining mobile/linear
amount equations. Explicit sorption iteration and backtrack limits supplement
the existing water-root work bounds. Nonsorbing and linear-only films retain
their direct solve and allocate no kinetic arrays.

Accepted and trial mobile/linear amounts and finite trapped amounts remain
separate, avoiding subtraction of a large retained inventory to recover a tiny
mobile pool. Independent kinetic and coupled conservation residuals must fit
the caller's increment/flux-scaled budget; an unrepresentable uptake cannot
silently disappear beside a large preloaded site. Kinetic residual and total
site-capacity closure are reported separately: mobile/flux accuracy audits
kinetics, while the original capacity-scaled closure gate remains intact.
One-ULP capacity summation cannot force endless refinement of otherwise
representable tiny uptake. The shared scalar owner preserves exact empty and
balanced-site limits without erasing the uptake derivative. Getters expose trapped inventory and the original site coefficients.
Preparation accounts for all new buffers; prepared steps and adaptive windows
allocate no storage.

Adaptive snapshots, rollback and separate error norms include both stored
pools. Failed nonlinear solves, water-root exhaustion, source refusal or late
cancellation publish neither pools nor the clock. Authored tests compare an
independent finite-bath quadratic, saturation/release refinement, unchanged
osmotic equilibrium at unequal bound loading, small mobile pools, sub-ULP
uptake refusal and transactional adaptive failures. These ideal independent
sites do not infer competitive occupancy, swelling, sorption heat or measured
material applicability.

### Accuracy-controlled transient osmotic films

`TransientOsmoticPair::advance` adds opt-in prepared temporal control. It
compares one full backward-Euler step with two half steps from the same state,
retains the fine conservative state without extrapolation, and controls every
bulk/film solute amount, both bulk water volumes, and both reconstructed
membrane surface concentrations. Supplied algebraic controls remain explicit
caps and are tightened below temporal allowances. Root work includes coarse
and refused trials; accepted-step, rejection and root budgets plus a minimum
step bound the work.

The exact representable requested endpoint is reached or the entire window
refuses. A one-ULP terminal clock remainder is included before solving, and
rounded unequal halves receive their actual algebraic budget shares and
first-order error factor. Any refusal or cancellation, including after the last
internal step, restores the original clock and inventories. Reports reconcile
original/final inventories directly with compensated summation and accumulate
only accepted physical pressure work. Prepared advancement allocates no storage;
the optional snapshots require a separate admitted byte budget, leaving
single-step storage unchanged.

These are local temporal estimates, not spatial or accumulated global-error
certificates or guarantees of nonlinear-root uniqueness. Authored regressions
check analytical finite-reservoir relaxation, layered-film startup and surface
polarization, exact/unequal-half clocks, work limits and late cancellation.
The existing filtration example uses explicit nonuniform layers, matches its
steady comparison by true series resistance and actual initial inventories,
and reports adaptive work and balances through the physical transient.

## Finite Donnan swelling and mechanical water exchange

`osmosis::donnan::DonnanOsmoticPair` couples an original finite-bath ionic
partition to finite water motion and the original compartment pressure/volume
laws. Preparation canonicalizes mobile ions at the supplied volumes while
retaining each species' total. Fixed charge and independent trapped osmotic
particles remain explicit. Every water trial re-equilibrates the finite ions
and uses FULL mobile-plus-trapped osmotic pressure exactly once, without ionic
reflection factors or an infinite salt bath. Compartment solution getters hold
mobile ions only; the Donnan report holds full osmotic pressure including traps.

One step uses backward Euler. Analytical nonnegative ionic stiffness augments
the actual mechanical tangent at each trial in a safeguarded Newton solve.
The original convex equilibrated free energy makes the scalar residual monotone. Bounded bracketing handles
volume applicability and stiff roots, and a narrow sign-changing bracket may
control transfer-position error while retaining the actual balance residual.
All mutable storage is prepared once. Every ion total, water, both fixed charges,
mechanical boundary work and accepted clock remain explicit. Failure, exhausted
inner/outer work or cancellation publishes neither inventories nor diagnostics.
The ionic work diagnostic counts completed equilibria; refused trials remain
bounded by the explicit outer and inner controls.

Both constructors represent fast-ion limits, not transient ionic permeability,
Nernst--Planck transport, double-layer storage or imposed voltage/current.
The ideal constructor refuses nonideal inventories; `new_nonideal` requires
both compartments to retain the original shared virial model.
Root tolerances do not control first-order temporal error. Mechanical work is
the exact original pressure primitive for linear or nonlinear chamber laws,
not a thermal closure.

The `donnan_swelling` example follows finite charged-bath swelling against
compliance, compares an independent monovalent algebraic equilibrium, and
reports voltage, conserved inventories and ideal-plus-mechanical energy.
Focused regressions cover linear and nonlinear analytical solutions,
counterion-only limits, multivalent charge conservation, stiff hydraulics,
orientation, energy decay, volume depletion, workspace limits and atomic retry.


### Activity-consistent finite-bath Donnan swelling

`DonnanOsmoticPair::new_nonideal` retains a prepared
`fs_material::osmosis::donnan::VirialDonnan` and the same original
`VirialSolutionModel` in both accepted solution owners. Ion partition obeys
`ln(cL/cR)+B*(cL-cR)+z*psi=0`, every species total, and both fixed charge
constraints. Water responds to the pressure from that SAME original free
energy. Concentration references cancel from the ionic difference. Trapped
particles remain explicit, ideal, noninteracting osmoles; fixed charge does
not infer a particle count or a second osmotic contribution.

The material owner uses a prepared native LU workspace and a diagonally scaled
exact constrained Hessian, preserving representable trace and exact absent
inventories without a concentration floor. When the ideal partition exceeds
the original concentration caps, a bounded charge-feasible seed starts the
original equilibrium solve. Caps remain applicability boundaries; a capped
state never replaces the original electrochemical equations. Individual
chemical-potential and charge tolerances, all residual evaluations, original
support, representability and cancellation remain explicit acceptance gates.

Analytic ionic redistribution under solvent transfer supplies the relaxed
osmotic stiffness. It is evaluated as positive original-Hessian quadratic
forms, avoiding cancellation of large direct and relaxation terms. The shared
water solve uses it alongside the original nonlinear mechanical tangent.
The existing exact material work, conservative adaptive half-step endpoint,
window rollback and bounded work remain active. Adaptive workspace admission
also compares original thermodynamic models, preventing snapshot reuse from
silently substituting a different solution law. Successful repeated material
and water solves allocate no memory after preparation.

Material regressions cover independent monovalent roots, original activity
and pressure formulas, multivalent cross interactions, finite-volume
sensitivity, trace-ion unit scaling, exact counterion-only limits, supported
roots outside the ideal initial guess, and atomic refusal/retry. Flux
regressions cover zero-interaction recovery, independent nonideal swelling
against Simo--Taylor mechanics, total-potential decrease, common-gauge
invariance, late cancellation and incompatible-model rejection.
The `nonideal_donnan_swelling` example combines finite charged baths, retained
dry volume, adaptive water motion and nonlinear bulk resistance, with an
independent scalar equilibrium and energy account.

These regression sources were added; native execution was blocked by shared
storage/memory exhaustion in the implementation environment. The physical
model is isothermal, single-phase and within the original convex virial
support, with negligible solute partial volumes. It does not supply a universal
electrolyte calibration, binding/association, double-layer storage, transient
ionic permeability, heat transport or spatial poroelasticity.

### Adaptive finite-bath Donnan advancement

`DonnanOsmoticPair::advance` uses an optional `DonnanAdvanceWorkspace` with
three prepared state snapshots. Each trial compares one backward-Euler step
with two actual representable half intervals and retains the conservative fine
endpoint. Both solvent volumes and EVERY mobile ion amount on both sides enter
the local temporal-error norm, including initially absent species. Unequal
rounded halves use their first-order Richardson factor. The estimate does not
certify global temporal error or fast-ion model validity.

Explicit maximum/minimum steps, accepted/rejected work limits and existing
inner root/charge limits bound work. Algebraic tolerances only tighten the
caller's caps; absolute water/charge allowances and relative charge accuracy
are distributed over the requested window. Charge tolerances retain the
original equilibrium owner's charge-scale semantics. An unresolvable terminal
interval or accuracy below the admitted step floor refuses; the clock is never
advanced without integrating transport.

Success reports accepted/rejected steps, all attempted implicit/root work,
accepted-half-step mechanical work and original-to-final water/species balances.
Refused ionic solves retain the original limited work-reporting semantics.
Failure or late cancellation restores the complete original inventories,
equilibrium diagnostics and clock. Reused successful advancement allocates no
memory. Single-step callers incur no extra snapshot storage.

Regressions cover analytic compliant relaxation and temporal refinement,
amount-controlled ionic repartition with loose water tolerances, counterion-only
square-root swelling without a salt floor, whole-window rollback/retry,
workspace limits and temporal-floor refusal. The `donnan_swelling` example now
uses this adaptive owner while retaining independent finite-bath equilibrium
and energy oracles.


## Finite-rate ionic osmosis and charged swelling

`osmosis::electrodiffusion::IonicOsmoticPair` evolves an ordered mixture of
up to 64 mobile species, two finite solvent volumes, and the original
mechanical response. Charged and neutral permeants can have independent
channel permeances; zero-permeance and globally absent species retain their
exact inventories. Its constructor preserves the supplied nonequilibrium
partition. Each bath must satisfy the explicit absolute charge budget with
its own retained fixed charge; no salt, counterion, or equilibrium reset is
supplied. A fixed charge does not itself create an osmotic particle count.

The shared material `ElectrodiffusiveMembrane` supports the original ordered
ideal or virial solution model. Water-selective pores use the full original
osmotic-pressure difference, including separately declared ideal trapped
particles. Independent species channels use exponentially fitted Nernst--Planck
transport with their original permeances. For
`a_i = ln(gamma_L,i)-ln(gamma_R,i)+z_i*F*(phi_L-phi_R)/(RT)`,

```text
J_i = P_i * [Bernoulli(-a_i)*c_L,i - Bernoulli(a_i)*c_R,i]
Bernoulli(a) = a / (exp(a)-1)
J_v = L_p * (DeltaP - DeltaPi)
sum_i z_i*J_i = 0
```

The floating potential solves the last equation from the actual species
fluxes. Both concentration endpoints, original activities, and full solvent
pressure remain attached to the same free energy. Positive fitted donor
coefficients admit exact zero endpoints without artificial trace
concentrations. Carrier-free and counterion-only limits retain their
indeterminate or infinite potential semantics. A singular instantaneous
entropy at a zero endpoint is unavailable rather than a fabricated finite
heat rate.

This is an explicit effective-channel model: electrical and excess chemical
potentials are reconstructed linearly through each channel, permeances are
constant caller inputs, and water and ion channels have no cross-drag
coefficient. It does not solve a nonlinear membrane interior, Debye layers,
capacitance, active pumps, binding, or electrode reactions. The floating voltage
assumes electrical relaxation is fast; finite species relaxation is evolved.
The pair is isothermal and does not turn its free-energy dissipation into
thermal-state evolution. No calibration or experimental authority follows
from numerical coefficients.

One backward-Euler step solves the original finite-rate species equations
inside each trial water transfer. The two-opposite-ion layout retains its
scalar salt solver. Larger and neutral mixtures use independent transfer
coordinates with one charged pivot eliminated by the exact charge constraint.
A prepared dense native LU solve supplies damped Newton updates. Original
inventory and concentration bounds constrain the private numerical trials;
a charge-feasible box initializer does not impose Donnan equilibrium.

The shared material response differentiates the same original membrane law
analytically. It retains the fitted donor coefficients, channel slopes and
original virial interactions after one potential solve, and eliminates the
floating-potential derivative through the zero-current constraint. Jacobian
columns require no repeated electrical root solves or concentration floors.
Derivatives at singular electrical limits are restricted to their supported
directions. Derivative output is atomic on failure.

Every accepted endpoint must satisfy every original species equation on both
sides, the original water equation, each finite-inventory bound, original
solution/mechanical support, and the full mixing-plus-mechanical potential
budget. The eliminated charged row is checked explicitly. Root convergence
alone cannot publish a state that fails these physical equations. Nonlinear
work and cancellation remain bounded. The original mechanical pressure
difference cancels a common pressure gauge before adding the actual restoring
response; work uses the exact original pressure primitive. Per-species reports
use `IonicSpeciesValues`: indexing, `as_slice()` and iteration expose the actual
active length without allocating a new numerical vector.

Prepared adaptive advancement compares a full step with two actual half steps
on water and every species inventory, retains the fine conservative endpoint,
and checks the whole requested window before publication. This remains a
first-order method with a local temporal estimate, not a certified accumulated
error bound. Failure, exhaustion, and late cancellation restore inventories,
diagnostics and clock. Workspace compatibility checks species shape and the
original solution law. Mechanics, membrane, fixed charge and trapped particles
remain attached to the active pair; snapshot restoration cannot substitute a
different physical owner.

The `finite_rate_ionic_swelling` example compares blocked, slow, and fast
binary channels and then evolves separate sodium/potassium/chloride channels,
a neutral permeant, and a retained tracer. All coefficients and species names
are illustrative roles, not fitted electrolyte data. Focused regression
sources retain the binary independent oracle and add a nonideal three-ion
backward-Euler oracle with opposite sodium/potassium transfer, original
mechanical work and energy, species permutation, exact blocked/absent pools,
and adaptive rollback.

## Sourced contact-line friction and hysteresis

`MaterialCapillaryTube::with_contact_line` opts the original tube owner into
the shared `fs-material::liquid::surface::contact_line` binding. The original
ordered wetting bundle owns Young angle, advancing/receding thresholds and
independent line friction. Its checked borrowed view preserves source identity,
rejects unsupported evolving rate/history coefficients, and is also used by
the existing planar bridge through its compatible resolver re-export.
Explicit speed and apparent-angle departure limits belong to the consumer.

The open two-fluid operator adds `zeta*R/4` to
`mu_l*x + mu_d*(L-x)`, selecting advancing, receding or pinned force branches.
It reuses the exact affine-force travel integral, bounded inversion and positive
loss moments. At zero speed the apparent angle balances the static reaction;
a pinned meniscus creates no mechanical work or dissipation. The moving law
obeys `cos(theta_app) = cos(theta_branch) - zeta*u/sigma`.
Both interval endpoints bound speed and apparent angle throughout this
constant-property open-tube trajectory. The scalar core allocates nothing.

Young wall free energy stays separate from both bulk-fluid losses,
`D_line/dt = 2*pi*R*zeta*u^2`, and nonnegative threshold pinning work.
The existing nested viscous fields remain bulk-only; the wrapper's
`total_dissipation_j` supplies the complete mechanical balance. Source-owned
horizontal, inclined and two-fluid methods consume the same optional binding.

The sealed midpoint owner includes these same resistance and threshold terms
in its safeguarded scalar solve and its independent energy balance. Gas heat
still enters the original closed parcel exactly once. Initial/endpoint branch
checks reject hidden reversals, oversized stopping crossings and heated pinned
steps that end outside the static interval. Mean and physical endpoint speeds
and apparent angles must satisfy the declared operating limits. Adaptive
substeps sample the interior; these checks are not a continuous-trajectory
domain enclosure. Deposited films, changing meniscus end-cap storage, temperature-dependent
liquid evolution and a general Cox--Voinov or molecular-kinetic law remain
outside this bounded linear-friction closure. Column inertia is admitted
through the explicit opt-in described below.

## Adaptive sealed capillary intervals

`SealedGasCapillary::advance_candidate` integrates a complete constant-boundary
interval or refuses without changing the original owner. One coarse midpoint
step is compared with two half steps in position, gas temperature and pressure,
using independent absolute tolerances and a relative tolerance. Only the
actual half-step states and transfers can be accepted; the gas state is never
Richardson-extrapolated. Smooth-branch errors use the midpoint factor of three;
detected pinning/depinning/reversal changes use the full discrepancy, reset the
existing PI controller's history and limit immediate growth.

Rejected trials add no time, heat, work or dissipation. Compensated sums retain
accepted liquid mass/volume, pressure work, wall/gravity storage, all three loss
channels and original gas heat/work/storage/entropy increments. The complete
interval independently passes its energy gate. Step/attempt limits,
cancellation, unresolved minimum-step error and represented-time stagnation
refuse publication. A shorter final interval reaches the exact requested time;
unsupported tube-end events never become clipped partial successes. These are
local error estimates, not certified global solution-error bounds.

Source-backed regressions cover complete-window conservation, rejected-trial
accounting, continuous adiabatic refinement, heated depinning against an
independent onset split, explicit source-domain refusals, cancellation,
bounded work, minimum-step failure and final-tail accounting.

## Inertial capillary motion and reservoir exchange

`capillary::inertia` adds physical velocity to the circular-tube state.
Both incompressible columns share one mean velocity, with moving mass
`M(x) = A*(rho_l*x + rho_d*(L-x)) + added_mass`. Added mass and advancing
and receding reservoir coefficients are explicit inputs. The force
`M*u_dot + M_prime*u^2/2` and directional reservoir loss `C*abs(u)*u`
are separate. For a single liquid, `C=rho*A/2` gives the bidirectional
reservoir convention in Fricke et al., Physica D 455 (2023) 133895,
doi:10.1016/j.physd.2023.133895, corrected preprint
https://arxiv.org/html/2305.05939v3. The two-fluid form follows the same
energy account; no universal entrance coefficient is inferred.

The midpoint discrete gradient balances actual endpoint kinetic energy,
Young wall free energy, signed open-system gravity and pressure work against
separate liquid, displaced-fluid, contact-line, threshold and reservoir losses.
A moving meniscus can coast through the static hysteresis band. A resolved
zero-velocity event returns its actual elapsed time and exact zero endpoint;
the caller restarts the remaining interval to pin or reverse. Pure viscous
coasting is not replaced by a fabricated finite stopping time. Tube-end
crossings, unsupported angle/speed domains and non-monotone oversized trials
are typed refusals. Cancellation and bounded iterations publish no state.

This optional closure admits inertial startup and damped rise/oscillation.
It uses lumped mean-flow inertia and instantaneous laminar resistance; developing
radial profiles, turbulent or compressible columns, films and universal
entrance-flow corrections remain outside its scope. Dry startup requires an
explicit positive effective mass. Temporal accuracy requires refinement;
a small mechanical residual alone is not a trajectory-error estimate.
Focused executable checks cover independent startup, drainage and damped-rise
limits, both motion directions, hysteretic stopping, all energy channels,
timed downstream pressure and refusal/cancellation boundaries.

## Source-owned inertial capillary intervals

`MaterialCapillaryTube::with_inertia` and `with_inertia_two_fluid` retain
the original liquid and ordered wetting owners. The latter also borrows the
original displaced-fluid owner and evaluates its admitted transport once.
Both phases' densities determine moving mass; their viscosities determine
separate drag. Explicit added lengths, four directional reservoir coefficients
and the operating speed limit are declared in `CapillaryInertiaParameters`.
No material name, guessed entrance length or hidden regularizer supplies them.

`SourceInertialCapillary::step_candidate` consumes a complete interval
through bounded stopping events. It restarts at exact rest, allowing pinning
or reversal without reflecting velocity. `advance_candidate` adds
independent length, velocity and stopping-time error scales around the existing
PI controller. Early-stop trials preserve a post-stop remainder, avoiding
repeated shrinking toward an earlier coarse stopping estimate. The full
coarse/fine discrepancy is used across branch changes; smooth intervals use
the midpoint factor of three. Only accepted conservative half steps contribute.

Complete-interval mass, pressure work, wall energy and gravitational terms
are reconstructed from actual endpoints; kinetic change uses a stable product
identity. Accepted bulk, line, pinning and reservoir losses are summed with
compensation. Source restrictions, cancellation, event/attempt budgets,
time stagnation and unsupported tube ends refuse the candidate atomically.
Zero duration and a shorter final remainder preserve the exact requested clock.
Estimates are local temporal controls, not certified global error bounds.

Regressions exercise added-mass dry startup, directional two-fluid exchange,
coast/stop/pin behavior, purely viscous asymptotic coasting and refusal boundaries.
The numerical owner/driver paths have executed with explicit numeric source-view
doubles; native original-card integration tests are supplied for the normal
workspace build and must not be confused with that narrower execution evidence.

## Inertial sealed-gas capillary dynamics

`SealedGasCapillary::with_inertia` adds physical liquid velocity and explicit
added/end inertance to the existing owner. The declared geometry must match
the actual sealed tube; gas bulk inertia remains excluded from the spatially
uniform closed-parcel model. The same original wetting source supplies friction
and advancing/receding thresholds. `SealedCapillaryInertiaAccuracy` provides
independent velocity and stopping-time scales without changing the existing
overdamped controls or choosing physical coefficients implicitly.

Every trial obtains midpoint pressure, endpoint pressure and both pressure
derivatives from the original NASA-9 gas heat/volume candidate. The shared
inertial mechanics receives the analytic chain-rule derivatives
`dP/dx = -A*dP/dV_final` and `dP/dt = heat_power*dP/dQ`. The accepted gas
candidate and its authentic receipt move together. A stopping segment consumes
only its actual elapsed heat. The whole-system energy gate balances inlet
work plus external heat against gas internal-energy change, liquid/end-mass
kinetic change, Young wall and gravity storage, and independent bulk, line,
pinning and reservoir losses. Dissipation is never silently added to gas heat.

The existing adaptive driver now controls optional velocity and event time
alongside meniscus position, gas temperature and pressure. A stopped coarse
path and refined path reach a common post-stop time through their own authentic
gas operations. At most five bounded step evaluations are needed per attempt;
the normal overdamped path still uses three. Accepted fine receipts are summed
individually, including any post-stop remainder. This avoids fake aggregate
gas derivatives, heat from rejected attempts and repeated shrinking toward an
earlier coarse stop. A second unresolved stop in a tail requires refinement.

An initially resting state must honor the original pressure's depinning branch:
heating cannot average an initially active force into false pinning or immediate
opposite motion. Nonzero incoming velocity remains free to coast through the
static band. Source/angle limits, tube-end refusal, atomic cancellation and
complete-interval energy checks remain in force.

Source-bound regressions cover independently manufactured advancing/receding
motion, pressure sensitivities, actual stopping heat and direct gas-receipt
identity, adaptive stop/pin completion, initial thermal branch changes, and
unchanged owners on refusal. Focused numerical execution uses explicit
constant-heat-capacity gas and material-view doubles with the actual new
mechanics/owner/driver; full NASA-9 dependency integration remains a separate
native workspace test requirement.

## Nonlinear interactions among restrained cavitation bubbles

With the `cavitation` feature, `BubbleCluster` advances up to sixteen original
`RayleighPlessetBubble` owners at fixed centers on one common physical clock.
The leading incompressible monopole interaction solves for
`S_i = d(R_i^2 R_i')/dt` through the symmetric, radius-normalized liquid inertia.
Native `fs-la` factorization and retained workspace avoid a second linear
algebra owner. Each observation uses the solved joint acceleration and neighbor
pressure; no gas inventory, equilibrium radius or source property is invented.

The matching liquid kinetic energy includes the pair cross terms. Original
gas volume work, capillary storage and the reference-minus-vapor pressure term
complete the mechanical energy. External pressure work and native viscous loss
are integrated as coordinates of the same trajectory, then checked against
the endpoint energy change before publication. The original DP54 tableau and
controller are shared through a private fixed-array kernel; the isolated
Rayleigh--Plesset and Keller--Miksis entry points retain their two-coordinate
path and arithmetic. Cluster state, work coordinates and scratch storage are
bounded independently of the requested duration.

Admission enforces exact shared liquid properties and clocks, separated
maximum-radius spheres, radius/separation at most 0.1, normalized inertia row
sum at most 0.5, and pressure imbalance divided by liquid bulk stiffness at
most 0.01. The caller declares the actual nonlinear volume-source bandwidth;
its propagation phase across the cluster is limited to 0.1 radian. This is an
explicit approximation assumption, not a waveform certificate. Native radius,
wall-Mach, attempt, step and cancellation gates remain active. Failed solves,
energy checks, domain exits and cancellation leave the accepted owner unchanged.

Eight native regressions cover isolated-bubble bit identity, independent
two-bubble acceleration and normal-mode frequencies, Lagrangian power,
nonlinear work/loss balance, member permutation, atomic refusal and domain
limits. Thirty-seven existing cavitation regressions also pass with the shared
kernel. This model retains fixed-center radial interaction only: translation,
jets, near-contact deformation, finite propagation delay, radiation loss and
cloud shocks require additional physical owners.

## Explicit velocity correction at pressure-supported equilibria

`FluxSystem::bdf1_navier_stokes_mixed_step_with_velocity_tolerance` and the
same-named `MaterialFluxSystem` entry retain the original nonlinear BDM1/P0
owner while admitting an explicit finite nonnegative absolute correction
tolerance in m/s. The original entry points pass zero and preserve their
relative-only behavior and public control structs. No flow, material, pressure
or history is clipped or replaced when the absolute criterion is selected.

The absolute norm is max_edges(|delta_mean| + 6*|delta_linear_moment|), exactly
the largest endpoint change of the original affine normal-velocity traces.
It is measured before under-relaxation and is independent of edge count.
It controls a normal-trace correction, not a continuum velocity error bound.
Pressure coefficients never enter this norm. A correction must meet this
authored absolute bound OR the original relative-L2 bound; both freshly
reassembled original and row-scaled nonlinear equation residuals must still
meet the unchanged relative tolerance. The reported relative_update remains
the actual ratio and may be large for an accepted near-zero field.

This resolves the concrete failure in which pressure supports a conservative
force with zero or tiny relative flow: solver roundoff could keep the old
relative correction ratio large even with an independently assembled equation
residual near machine precision. Focused native regressions exercise zero and
small physical velocities on two mesh sizes, one-iteration absolute admission,
zero-absolute parity with the existing path, and refusal when a generous
correction allowance does not satisfy the nonlinear convection residual.

## Nonlinear pressure flow on exactly translating domains

The original material BDM1/P0 Navier–Stokes solver is now advanced by
ale::translating::TranslatingAleFlow on x = X(t) + y. The relative velocity
v = u_world - W obeys the same nonlinear convection and stress operators,
with the physical body force sampled at y + X1 and inertial source
-rho*(W1-W0)/dt. The original pressure remains physical mechanical pressure;
its mean-zero or imposed-traction datum is not converted into an EOS state.

Translation has identity deformation gradient and unit Jacobian. The domain
actually changes world position while its metric, topology, bases and material
snapshot are retained. BDF1 uses X1 = X0 + dt*requested_W1; the represented
W1 = (X1-X0)/dt must meet an explicit velocity tolerance and then owns the
force, relative boundary velocity and world reporting conversions. Independent
position/speed histories cannot enter. Each endpoint callback sees world
coordinates. Relative history, geometry, frame velocity and clock publish
together after successful native nonlinear admission, physical conversion and
the final cancellation poll. Native synchronous solvers are polled around,
not internally interrupted by this composing owner.

The reported world BDM velocity, linear momentum and kinetic energy use the
same source density. Physical stress and shifted moment origins remain
separate from outward momentum through the moving control surface:
rho*u_world*((u_world-W).n). DiscreteBoundaryReaction now retains signed volume
flux and its first spatial moment from the SAME original advector quadrature
as its momentum receipt. A native selected-edge accessor and a material
mass-flux accessor expose those moments without reassembly or source queries.
The world receipt adds W times those retained mass moments; it never mixes a
later endpoint normal trace with an earlier frozen-advector momentum flux.

Velocity, NaturalTraction and OpenReservoir boundaries are supported. Reservoir
velocities use the same world-to-relative conversion as prescribed velocities;
world momentum uses the upwind reservoir velocity during inflow. The original
DissipativeBackflow law is fixed-frame and refuses rather than silently
changing its reference velocity. Direct `step` calls prescribe translation.
`candidate` solves exactly the same endpoint while uniquely borrowing its
original owner. Dropping it preserves all accepted histories; `commit` publishes
the admitted candidate without another source query or flow solve. The unique
borrow prevents stale or wrong-owner acceptance. Coupled owners finish their
own mechanical and cancellation gates before this infallible publication.
The candidate explicitly accepts the original nonlinear owner's absolute
normal-velocity correction tolerance; direct step calls pass zero. Neither
criterion removes the independent nonlinear equation residual gates.
`fs-couple::stokes_motion::translating::GuidedTranslationMotion` uses these
candidates to solve finite guided body translation from the original discrete
fluid reaction. Neither entry supplies rotating/deforming ALE, topology
change, a nonlinear adjoint, or a whole-domain energy/passivity theorem.

Five authored G1/G3/G4 regression groups cover acceleration-induced pressure
and reaction; preservation of uniform world flow under changing mesh speed
with independent moving-surface momentum flux; Galilean covariance of a
nonlinear mixed-boundary affine flow; a retained advector deliberately distinct
from the endpoint velocity; and final cancellation, continuation and motion/law
refusal. Native execution and formatter checks were unavailable after the
session worker failed; these are authored cases, not claimed passing results.

## Dynamic channel heat transport in evolving mixtures

`MolecularTransport3d::with_dynamic_channel_thermal_sgs` exposes
`MixtureDynamicChannelThermalSubgridConfig` on the existing evolving-mixture
finite-volume transport owner. One reflecting axis has two explicit no-slip
walls; the two tangential axes are periodic. The shared channel core fits
momentum and temperature independently within homogeneous wall-distance planes
using current SSP-stage density, velocity, temperature, Cartesian spacings and
current wall velocity/thermal data. The original mixture face retains source
Cp, referenced partial enthalpies, Maxwell–Stefan diffusion, conservative
stress work, molecular wall heat exchange and the combined stability estimate.
Heat remains active when the momentum fit vanishes; species mixing retains
an explicitly supplied common Schmidt number. All normal orientations and
rectangular uniform cells are supported. Tangential directions require at
least three cells; two normal cells suffice for fixed-temperature walls,
while any adiabatic wall requires three. Complete storage adds 264 nominal
numerical bytes per fluid cell, with reused buffers and atomic field/report
publication. Applicability remains low-Mach exploratory LES, with resolved
walls by default and optional smooth adiabatic equilibrium wall stress;
no reaction-rate model or
DNS-accuracy qualification is claimed.

Four new native tests cover original source Cp/partial enthalpy and wall work,
all axes and rectangular geometry, independent heat with zero momentum fit,
formation-reference invariance, evolving walls, conservative advancement,
capacity admission, cancellation and retry. The broader 25-test molecular/SGS
group also passes with the actual current fs-lbm, fs-gas and fs-tribo dependencies.

## Boundary momentum and kinetic-energy transport

The original mixed-flow reaction assembly now retains its SAME sign-split
advector/upwind quadrature in four inline records: position, carried velocity
and weighted normal flux. `BoundaryReactions::transport` and the single-edge
`DiscreteBoundaryReaction::transport` expose volume flow, linear/angular
momentum, carried kinetic power and boundary upwind dissipation together.
Selections reduce in canonical edge order with fixed-size compensated sums.
No field, material or reservoir callback is queried during reporting.

For normal relative flow q, interior velocity u and the admitted upwind trace g,
carried kinetic flux is q|g|²/2. The numerical inflow loss is
max(-q,0)|u-g|²/2. Since g=u on outflow, their sum equals the boundary kinetic
term q(u.g-|u|²/2). The loss is separate from interior-face upwinding, viscosity
and any declared backflow resistance, and is not silently deposited as heat.
Both terms are cubic on an affine BDM1 branch and use its existing quadrature.

A uniform carried-velocity offset preserves the original control-surface-
relative flow. Each retained velocity is shifted before its momentum and
kinetic-power products, and each position is translated before its torque
product. No sum over different samples is formed before these conversions.
This avoids subtracting large absolute kinetic fluxes merely to recover small
world-frame exchange, and preserves tiny reversing-flow contributions even
when their volume contribution rounds away from the total. Counterflow can
carry nonzero momentum and energy while its net volume vanishes.
Finite-precision admission and existing solve tolerances still apply; no exact
state, large-coordinate accuracy or whole-scheme passivity theorem is asserted.

`MaterialStokesSolution::boundary_transport` uses the original material density
for mass, momentum and kinetic power in SI per unit depth. Newtonian and
coupled-rheology Navier-Stokes endpoints share that owner.
`TranslatingAleStep::boundary_transport_world` shifts carried velocities by the
accepted mesh rate and translates the requested moment origin. The existing
world boundary-load method now consumes the same momentum reconstruction.
Pressure/Cauchy-stress work remains separate; these kinetic quantities do not
supply enthalpy, internal energy, species, entropy or a complete StreamPort.

Four focused kernel regressions exercise large common velocity offsets,
zero-net-volume counterflow, energy/work identity, invalid frames, bounded
storage and tiny-flux preservation for both signs and sample orders. Native
reversing-face tests independently integrate cubic transport polynomials;
source-density, endpoint rheology and translating-domain tests retain their
original owners and verify kinetic transport, frame covariance, unchanged
callback counts and canonical selection. Stored transport is fixed-size per
boundary edge (176 nominal bytes on 64-bit targets), with a capacity shared
with the native quadrature generator and an explicit overflow refusal. Reporting does not resample fields and needs no nonlinear solve or per-face
heap allocation.

## Physical Navier-slip wall momentum

The original mixed BDM1/P0 Stokes and Navier-Stokes owners accept `NavierSlip`
alongside prescribed velocity and open total stress. Normal BDM moments remain
strong: u.n = wall.n. Tangential Cauchy traction obeys
t.(2 mu D(u)-p I)n = -beta*(u-wall).t with finite edgewise beta >= 0 [Pa s/m].
The source-material paths divide beta, viscosity, pressure and inertia by the
SAME retained density; nonzero physical resistance may not silently underflow
to perfect slip. Zero beta selects perfect slip, while the original Velocity
condition retains exact no-slip. For Newtonian viscosity mu, slip length is
mu/beta; beta itself is an authored physical interface coefficient.

The original sparse momentum operator adds the exact tangential Robin bilinear
form and wall-velocity load. It retains symmetric stress, interior SIP, mass,
history, nonlinear convection and constitutive endpoint ownership. No boundary
Dirichlet penalty is added to a slip edge. Free-slip faces add no Robin entries;
finite faces reuse six tangent basis traces per quadrature point. Steady
configurations that retain an unrestrained rigid mode retain the native singular
or convergence refusal; finite BDF1 inertia remains the transient anchor.

The SAME removed normal rows recover the normal constraint reaction; adding
the equal/opposite physical tangential Robin effort recovers fluid-on-exterior
force, moment and rigid virtual work. `BoundaryLoad::slip_dissipation` and
`MaterialBoundaryLoad::slip_dissipation_w_per_m` retain integral
beta*|(u-wall)_t|², separately from stress work, viscous volume loss, numerical
upwinding and transported kinetic energy. No thermal recipient is invented.

Slip does not prescribe an incoming reservoir velocity. Native convection and
the retained reaction/transport samples use the interior carried velocity on
both signs of prescribed normal flow. Nonzero normal data on a fixed mesh
therefore define an explicit transport boundary, not an impermeable stationary
wall. In translating ALE, physical wall velocity is transformed with the SAME
accepted mesh rate; an impermeable wall has zero relative normal flux. Slip
traction and loss are invariant under a common boost, while world momentum
transport retains the existing original quadrature owner.

Material Stokes, Brinkman, rheology, nonlinear flow and translating ALE all
retain the law. The existing linear guided/planar Stokes response keeps beta
when homogenizing boundary data. The finite two-axis body owner exposes
`PlanarTranslationMotion::new_with_slip`; its original constructor stays
no-slip. Native stress closes body momentum and work before body/fluid/mesh/time
commit, and physical interface slip loss is reported separately in joules.

Four native regression groups cover analytic Couette force, normal pressure,
torque, bulk-plus-slip power, zero-friction limit, nonlinear constitutive
endpoints, source-density scaling, accelerated translating frames, boundary
callback reuse, prescribed-normal momentum transport, and invalid coefficients
or normal flux. A coupled regression closes an actual translating body's
impulse/work and separate slip loss against an independent affine-flow oracle.
Pinned nightly rustfmt parsed all changed Rust in this session. Native Cargo
execution was unavailable; these added regressions are not claimed passing.

## Native nonlinear flow on rigid rotating and translating domains

`ale::rigid::RigidAleFlow` advances the SAME material BDM1/P0 Navier-Stokes
owner on x = X + R(y-c), with one fixed reference pivot c. Reference metrics,
basis functions, topology and physical area are retained under the isometry.
The velocity unknown a = R^T u_world is absolute fluid velocity in current mesh
axes. Its native advector is a - R^T V - omega J(y-c); the transported velocity
is a. World normal continuity remains the original Piola/BDM normal continuity.

The native BDF1 mass load integrates R_new^T R_old a_old from the ORIGINAL old
cell field at the original quadrature points. Rotating global edge moments
would mix a discontinuous tangential trace into a normal one and is not used.
Only six rotated load entries are accumulated per cell; the existing mass
matrix, symmetric stress, pressure/gauges, nonlinear Picard acceptance, removed
normal reaction rows and sign-split upwind quadrature keep their original
owners. One relative-advection workspace is reused throughout the nonlinear
solve and its independent endpoint residual reconstruction. The ordinary
stationary/translating entries retain their original history arithmetic.

`RigidAleMotion` supplies physical pivot velocity and angular rate.
`kinematics` previews the exact represented pivot, unit planar rotation,
signed rotation increment and rates used by `candidate`. Pivot displacement
uses backward Euler. Rotation uses exp(dt*omega); each increment is strictly
inside (-pi,pi), and both represented rates must satisfy authored absolute
tolerances. Larger turns require substeps. A coupled mechanical owner can
author its wall callback from this SAME admitted pivot/rate before the fluid
solve. Candidate ownership and commit mirror the original translating path:
geometry, flow history and time publish together only after complete admission.
Cancellation polls the coordinate/storage passes and surrounds the synchronous
native solver; it does not interrupt that solver's current internals.

Physical force density, velocity, traction points and normals all transform
through the admitted endpoint pose. One retained material density scales mass,
viscosity, force, pressure, boundary effort and storage. Navier-slip resistance
is unchanged; OpenReservoir retains its independently authored incoming
absolute momentum. The existing fixed-frame DissipativeBackflow law refuses.
Impermeable wall velocity is the admitted rigid mesh velocity at that point,
so native relative normal flux is zero. No artificial rotational inertia,
Coriolis force proxy, remapping solve or replacement flow discretization is
introduced.

The result exposes actual world mass, linear momentum, angular momentum about
the moving pivot or any authored world origin, and native physical kinetic
energy. World boundary queries rotate native force vectors and transform the
moment origin; planar torque, slip loss and carried kinetic power are rotation
invariant. The retained normal flux is relative to the moving surface while
its carried velocity is already absolute. Mesh momentum, stress and kinetic
transport are not counted twice. These are physical storage/flux queries, not
a complete StreamPort or whole-system energy/passivity certificate.

This is covariant first-order time integration on a prescribed global rigid
mesh. Finite-step angular conservation across changing geometry is NOT an
exact theorem of this backward-Euler map; its actual angular inventory and
boundary torque must be distinguished from temporal truncation. Six authored
G0/G1/G3/G4 regressions cover free-stream preservation and independent world
transport, discontinuous tangential history with linear-momentum closure,
parity with the original translating owner, rotating finite-slip Couette
stress/torque/source scaling, first-order physical rotation refinement, and
preview/drop/late-cancel ownership. Pinned nightly rustfmt parsed all changed
Rust. Native compilation and tests were unavailable and are not claimed passing.

Independent moving obstacles, mesh deformation/topology change, coupled body
dynamics, variable density, rotating-domain rheology, nonlinear adjoints and
automatic thermodynamic deposition remain outside this bounded owner.

## Source-owned nonlinear rheology on rigid moving domains

`ale::rigid::RigidAleFlow::{candidate_with_rheology,step_with_rheology}`
resolve generalized-Newtonian transport from each original absolute-velocity
candidate's objective strain invariant sqrt(2 D(u):D(u)). The original mixed
Picard owner rebuilds BOTH constitutive stress and relative mesh advection at
the returned endpoint. Constant-density inertia, covariant history quadrature,
physical source conversion, pressure, slip/reservoir laws, stress reactions and
transport retain the same native owners as the Newtonian rigid-domain path.

The existing material response cache now accepts an optional physical-position
map. Every resolver query and retained response uses the endpoint WORLD cell
centroid, while original cell identity and strain-keyed caching remain intact.
This adds no centroid vector or mesh remapping. The Newtonian entry performs
no new strain evaluation/allocation; both entries reuse the same relative
advection workspace and original sparse/native nonlinear solver.

`RigidAleStep<G = ()>` and `RigidAleCandidate<..., G = ()>` preserve ordinary
Newtonian call sites. The step's `rheology_samples()` returns exact original
endpoint material owners in cell order, including world position, physical
strain and validated transport. It is empty for ordinary Newtonian steps.
Typed source refusals and changed density refuse before publication. Geometry,
clock, history and all material samples remain one dropped/committed candidate;
reading returned responses does not query the source again.

Four focused native regressions cover changing translation/spin and
inhomogeneous thinning/thickening with an independent layered constant-stress
solution, world-position ownership and objective zero-strain free streams,
parity with the original stationary nonlinear material owner, and typed source/
density refusal plus dropped/finally-cancelled candidate retry. The layered
fixture changes both old shear rates and angular speed, checks nonzero physical
stress/torque and density-scaled momentum, and admits literal reservoir inflow.
Pinned Rust parsing passed. Native compilation and numerical tests remain
unrun because local storage is full and remote build tools are unavailable.

Reference density remains fixed. Captured thermal/composition fields are fixed
within each solve; thermal evolution, viscoelastic memory, constitutive fitting
and multiphase density transport are outside this owner. The existing rigid
BDF1 geometry/time truncation limits remain unchanged; nonlinear material
closure does not establish exact finite-step whole-fluid angular conservation.

## Original-source binary molecular diffusivity

`BinaryMiscibilityNetwork::bind_diffusivity_source` attaches an immutable
`fs-material::liquid::solution::MaterialBinaryDiffusivity` to the SAME accepted
binary field. The ordered species must match its original thermodynamic pair.
The returned candidate preserves both component inventories, compensated low
words, composition, geometry and elapsed time. Its source explicitly replaces
the authored constant face diffusivities; zero area still closes a geometric
face, while zero source diffusivity closes that stage's molecular kinetics.

At each canonical pair stage the source resolves the actual absolute temperature,
both current endpoint compositions and their arithmetic-mean face. The face
value supplies the nonnegative Maxwell--Stefan diffusivity in the original
`c D x_A x_B / (R T)` mobility. It is held during that scalar convex-split root;
later pairs, timesteps and adaptive coarse/fine trials resolve it again. Thus
the existing monotonic root, species conservation and free-energy inequality
are retained with first-order lagged nonlinear kinetics. This is not a fully
implicit composition-dependent D scheme, a temperature-gradient/Soret law, or
an inferred effective Fick coefficient in the spinodal region.

Prescribed temperature changes re-evaluate the pinned source without accumulated
rescaling, and the existing adiabatic and liquid-material adapters preserve it.
Published diffusion, advection and temperature endpoints qualify every original
cell and face, including disconnected cells. Source-domain refusal, negative
kinetics and cancellation publish no partial state. Explicit flow rebinding may
retain its donor/continuity caches across sourced-D changes only when the original
thermodynamics, capacities, geometry and gradient coefficients still agree.

`faces()` exposes binding/temperature snapshots; transport never treats them as
a last-successful-state cache. The numerical workspace retains its original
bound. Shared material-card storage and transient source-query/receipt allocations
are additional; O(cells+faces) bounded source queries per sweep occur outside the
scalar root iterations. `tests/miscibility_material/diffusivity.rs` checks constant
source equivalence, original-source temperature cycling, composition refresh,
adaptive stages, conserved inventories, dissipation and transactional refusal.

## Conservative spatial binary miscibility

The opt-in `experimental-miscibility` feature exposes `miscibility` without
requiring the continuum runtime feature. `BinaryMiscibilityNetwork` prepares a
closed finite-volume graph using the original
`fs_thermochem::solution::gradient::BinarySquareGradientV1`. Geometry and the
binary diffusivity on each face are explicit physical inputs; no coefficient is
inferred from a negative effective Fick diffusivity.

For cell capacities `N_i=c*V_i` and face weights
`K_e=kappa*A_e/L_e`, the retained energy is
`F=sum_i N_i*g_mix(x_i,T) + 0.5*sum_e K_e*(x_i-x_j)^2`.
Each face moves A and B in equal opposite molar amounts. This is a
molar/volume reference frame with equal fixed partial molar volumes; unequal
molar masses do not imply zero barycentric mass flux. Local density can change
with composition. Pressure and geometry remain externally fixed; uniform
temperature can be changed atomically through the temperature owner below.
The molecular lane does not advance heat, momentum, gravity, wetting or
nucleation rates. Prescribed volume-flow advection is a separate operator below.

One canonical face sweep solves an implicit convex-split scalar exchange on
each face. Logarithmic entropy, convex interactions and the entire graph's
gradient energy are implicit; positive-chi concavity uses the preceding pair
state. Every neighboring gradient contribution participates, including faces
outside the exchanging pair. A safeguarded Newton/bisection solve stays inside
the physically available paired inventories. Pure opposite components can
mix without a logarithmic floor or invented initial trace. Both component
amounts are represented independently, with retained addition residuals
preventing drift across long repeated exchanges.

Each immutable candidate audits both global species balances, every fixed
cell capacity, the direct actual free-energy increment and the convex-split
dissipation inequality. Nonzero transfers that disappear from represented
inventory arithmetic refuse. The free-energy allowance applies once to the
whole sweep, not once per face. Cancellation or failure publishes no partial
state.

`miscibility::adaptive` advances a requested interval using one full sweep
versus two half sweeps. The accepted two-half state remains first order in
time. The local estimator considers both species with explicit absolute and
relative mole-fraction tolerances; it is not a global error certificate.
Step-size bounds and a maximum attempt count bound work. The report includes
accepted/rejected intervals, completed-sweep work, integrated face transfers
and the suggested next step. A final independent species/capacity/energy
audit applies the original single allowance to the entire interval and sums
only accepted dissipation. Failure anywhere leaves the input unchanged.

Canonical faces and CSR incidence are prepared once. Retained arrays and
conservative transaction scratch, including adaptive trials, are admitted
under a caller-supplied payload budget of 512 bytes per cell plus 512 bytes per
face and the fixed owner size; allocator overhead and pair-model strings are
excluded. Graphs are bounded to 65,536 cells and 262,144 faces; scalar roots
are bounded to 256 iterations. Work scales with faces times root iterations
and with `sum(degree^2)` for cache updates, which is linear on bounded-valence
finite-volume grids. The area/distance formula requires an appropriate
orthogonal mesh or explicitly authored compartment connection.

Focused native tests cover full-graph energy derivatives, pure contact,
stationary interfaces, independent spinodal growth/short-wave damping rates,
gradient-cache consistency, first-order refinement, canonical symmetry,
long evolution, adaptive error reduction against an independent RK4 scalar
oracle, integrated flux/species/energy closure and atomic budget/cancellation
refusals. These feature-level tests do not replace the full default-feature
workspace gate.

### Prescribed flow and finite-rate mixing

`miscibility::advection` transports the original paired species inventories
with a closed incompressible volume flow. Each canonical face retains two
nonnegative directional rates: opposing streams mix species even when their
net volume flux is zero. Preparation checks both rates, face areas and every
cell's volume balance, with an explicit relative roundoff allowance no larger
than `1e-12`. It does not project the flow or normalize away capacity drift.

Synchronous donor-cell updates use the same preceding state on all faces.
Outgoing-volume CFL substeps preserve positivity without clipping or seeded
traces; a caller-supplied substep cap bounds work. Both components retain their
compensated inventory remainders. The final independent audit applies the
original species/capacity allowance once over the entire interval. Reports
retain signed integrated component transfers and the actual bulk plus gradient
free-energy change. That change may have either sign under prescribed stirring;
it is neither dissipated heat nor a closed mechanical-work account.

`miscibility::transport` composes advection and the original implicit molecular
exchange over ONE physical interval. Complete stage reports survive; the
combined species/capacity audit uses the smaller original allowance once,
and independent total-energy increments expose the stage-accounting residual.
The physical clock advances once. Failure or cancellation in either stage
publishes no candidate. This Lie split and the donor transport are first order:
CFL subcycling alone does not control splitting error or numerical interfacial
mixing. The caller must refine the macro interval and spatial discretization.

`miscibility::transport::adaptive` advances the COMPLETE prescribed stirring
and demixing operator with one full Lie step versus two half steps, accepting
the latter. Both independently represented species and their retained low
words participate in the local composition estimate. Every trial uses exactly
one donor advection step, so the macro step is restricted by the original
outgoing-volume CFL limit. This prevents a subcycled coarse path from becoming
identical to its two half paths and falsely reporting zero advection error.
The original time-step minimum may refuse that restriction; a final shorter
remainder must satisfy the same accuracy and conservation checks. The frozen
prepared flow and fixed-volume reference frame remain unchanged.

The adaptive result retains separate accepted advection and diffusion face
transfers, their combined species transfers, each operator's free-energy
increment, accepted diffusion dissipation, and completed trial work. Rejected
and coarse trials contribute only to work counts. A final independent species
and cell-capacity audit uses the stricter original amount allowance ONCE over
the requested interval. Accepted diffusion and the direct total energy change
minus accounted advection must each satisfy the original diffusion-energy
allowance once; total free energy may rise under imposed stirring. The reported
advection increment is neither heat nor mechanical work. No trajectory history
is retained; bounded attempts and scalar solves use O(cells+faces) payload
within the original graph workspace admission. Cancellation or any failed
final audit publishes no partial field or flux account.

Focused G0/G1/G3/G4 tests compare the coupled trajectory with independently
integrated three-cell chemical-potential equations, exercise the CFL-alias
failure against an exact counterflow solution, and cover tolerance refinement,
zero-flow diffusion, pure and trace species, unequal volumes, positive energy
under stirring, separate integrated flux closure, endpoint clocks and atomic
refusals. The estimate remains local and the accepted method remains first
order; it does not bound spatial donor mixing or establish a global error
certificate, evolving momentum, thermal feedback or measured material validity.

Prepared flow data and reused O(cells+faces) transaction scratch fit the original
graph payload budget. Work is linear per advection substep. These are fixed
geometry, isothermal, equal-partial-molar-volume fields. A prescribed volume
velocity is not automatically a barycentric velocity for unequal molar masses;
momentum, capillary feedback, thermal evolution and open boundaries remain
outside this operator.

Focused G0/G1/G3/G4 regressions cover exact donor exchanges, zero-net
counterflow, periodic translation refinement, unequal-volume inventories,
integrated face balances, homogeneous/pure/trace compositions, positive energy
change under forced nonconvex mixing, combined diffusion/advection limits and
time refinement, the single physical clock, and atomic source/work/cancellation
refusals. The numerical fixtures are not measured liquid-pair validation.

### Native flow and original material sources

With `continuum`, `prepare_advection_from_bdm` consumes the existing
`FluxSystem::edge_volume_transport` sign-split normal traces. It requires the
same triangle volumes, physical depth, interior adjacency and face areas as
the miscibility graph, while retaining the graph's authored diffusion distances.
It does not infer an orthogonal diffusion stencil from arbitrary triangles.
Both physical directions survive orientation changes. Every exterior stream
must vanish; zero-net exterior counterflow is still an open boundary and refuses.
Continuity is checked on original mean moments before subtracting large
opposing streams, and again on the final directional rates. Geometry and
native extraction are bounded before allocations; cancellation is polled
around extraction and through the mapping. The frozen prepared flow is reusable.

`experimental-miscibility-material` enables source-bound fields independently
of the continuum runtime. `resolve_liquid_field` binds every original cell
inventory to the SAME `BinaryLiquidSquareGradientLaw`; the full original pair,
temperature, concentration and gradient coefficient must match. The field
privately retains the exact original liquid source across subsequent transport
and temperature changes.
Density, viscosity, heat capacity and conductivity come from each cell's
actual homogeneous composition and original selected material claims. Unstable
cells remain homogeneous; independently flashing them would erase the spatial
model's chemical and interfacial energy. Source-derived volumes must pass the
existing equal-partial-molar-volume check at every accepted composition.

The material field's `transport_candidate` refreshes sources only after the
entire numerical advection/diffusion step. A late source-domain failure,
incompatible volume or cancellation discards the complete candidate. An
explicit maximum resolved-cell count bounds source queries and retained
results; source-receipt allocations are separate from the numerical graph's
payload budget. These results expose the existing `FluidTransport` and
`ThermalFluidTransport` consumers. They do not add momentum or heat feedback.

Integration fixtures exercise an actually solved closed Stokes field,
sign-reversing traces, orientation/depth scaling, original mean continuity,
source-backed density changes with unequal molar masses, nonlinear viscosity
and thermal-property refresh, retained claim identity, homogeneous spinodal
states, source-domain refusal after otherwise admissible diffusion, and atomic
cancellation. They retain the existing prescribed-flow and source-data limits.

### Uniform temperature changes and original calorics

`miscibility::temperature` changes the prescribed uniform temperature of an
existing spatial field without resetting composition, interfaces, compensated
inventory remainders or elapsed time. The original pair law determines the new
chi and the same validated square-gradient constructor rejects unsupported
temperatures or a zero-kappa quench into a spinodal regime. Cooling can therefore
activate finite-rate demixing; no cell is independently flashed or phase seeded.

The returned candidate contains the new network and its admitted state.
Unchanged cells, capacities and CSR/gradient geometry share immutable arrays.
Face mobilities are recomputed from original concentration, authored diffusivity
and geometry with the new RT, avoiding accumulated temperature-ratio error.
Diffusivity itself remains the caller's explicit coefficient; this API does not
infer a temperature-dependent diffusivity or thermal expansion. Prepared volume
flow is explicitly reusable through `rebind_advection`, which requires the same
pair, concentration, kappa and every physical cell/face and shares the admitted
rate caches. An old flow or state cannot silently enter a different thermal owner.

For the original `chi=A+B/T` law and fixed composition, mixing entropy is
temperature independent and gradient energy is unchanged. The direct external
field free-energy increment is `-delta(T)*S_mix`, retaining small changes beside
large interface energies. Either sign is admissible. This parameter-change
account is neither supplied heat nor dissipation.

`BinaryMiscibilityTemperatureCandidate::caloric_change` optionally uses the
same full pair and explicitly qualified ordered pure-liquid Cp sources.
Direct original Cp and Cp/T integrals, evaluated once per species, supply the
separate sensible heat and system entropy increments using both words of every
paired inventory. No bulk equilibrium solve, thermal root, reference-energy
subtraction or replacement of unstable spatial cells occurs. Pure-liquid
qualification and the fixed-volume/maintained-pressure approximation remain
explicit caller responsibilities. These APIs do not advance heat conduction,
adiabatic feedback or reservoir entropy.

Focused native regressions exercise quench-induced spinodal growth, ideal
diffusion's temperature cancellation, retained traces/history/clock, exact
thermal-owner round trips, direct small energy increments, flow rebinding,
analytic piecewise Cp/CpT calorics and atomic domain/cancellation refusals.

### Source-bound temperature cycles

The resolved material field's `temperature_candidate` applies the same numerical
temperature change and then refreshes every actual homogeneous cell through
the original liquid-law reference. `BinaryLiquidSquareGradientLaw::at_temperature`
retains the original pair, concentration, kappa, complete liquid source binding
and original volume-mismatch allowance. New wrappers are owned through Arc;
subsequent transport shares the active wrapper without cloning the pair law
or replacing source claims. `field.network()` exposes the active immutable
numerical owner for explicit flow rebinding.

Actual density, viscosity, heat capacity and conductivity are queried at the
new temperature and unchanged paired inventories. Temperature-dependent
density must still match the fixed-site model; thermal expansion is not
silently suppressed or accepted by loosening the original tolerance. Any
material-domain/volume refusal or cancellation discards the entire numerical
and material candidate. Cooling retains homogeneous spinodal cells until the
spatial solver itself evolves them. The external free-energy report remains
distinct from supplied heat.

Integration regressions use real temperature-dependent property curves and
two-axis density data to exercise heating, cooling, source identity, continued
transport, strict density/domain refusal and cancellation throughout refresh.

### Adiabatic finite-rate diffusion

`miscibility::thermal` binds the existing field to the exact original
`BinaryLiquidCaloricModelV1` and ordered pure-liquid Cp sources. The binding
retains the original conserved species totals, including their compensated
low words, across all accepted intervals. It does not reconstruct calorics
from mixture transport Cp or replace the evolving cells by bulk equilibrium.

`adiabatic_candidate` advances one molecular-diffusion interval at the
preceding temperature, then solves one uniform-temperature energy balance at
the completed spatial compositions. The energy equation is
`delta(H_mix) + delta(E_gradient) + sum_s N_s integral Cp_s dT = 0`.
The original fixed volumes and maintained pressure make the change of pV zero.
Constant kappa contributes interfacial internal energy and no entropy. This
is the limit of fast internal thermal equilibration, with no resolved local
heat conduction, thermal expansion, thermodiffusion, flow work or external heat.

Direct original mixing-enthalpy and mixing-entropy increments retain traces
and small near-critical changes without subtracting large reference energies.
Gradient energy also uses direct face increments. Each thermal root trial
evaluates only the two original pure-Cp integrals; cells are not repeatedly
flashed or re-resolved. An explicit source-supported temperature bracket,
temperature/energy tolerances and at most 256 bisections bound the solve.
Unsupported heat targets refuse instead of clamping to a source limit.

Every accepted candidate independently audits the actual whole-field energy
residual and `delta(S_mix) + sum_s N_s integral Cp_s/T dT`. The entropy allowance
applies once. A split step can lower isothermal free energy and close energy
while losing total entropy; it therefore refuses if the separate entropy gate
fails. Refining the physical interval is necessary for this first-order thermal
split, even when the scalar energy root is very accurate. The physical clock
advances once, and the accepted temperature updates the original chemical
potentials and mobilities through the same temperature owner above.

The complete diffusion, temperature, sensible/mixing/interfacial energy,
entropy and scalar-work reports survive publication. Any solver, domain, work,
energy, entropy or cancellation failure leaves the original field unchanged.
The six thermal native regressions cover demixing and interface heating,
exact ideal isothermality, conservation, entropy-loss refusal with successful
refinement, atomic boundaries and first-order composition/temperature convergence
against an independent energy-conserving scalar RK4 reference.

### Source-qualified adiabatic fields

`miscibility::material::thermal` binds an already resolved spatial field to its
original `BinaryLiquidCaloricLaw` through `field.bind_caloric`. The exact
liquid-law reference must match, so a numerically identical replacement cannot
substitute another material card, claim selection or pressure context. The
source law retains the originally qualified ordered pure-liquid Cp curves,
molar-mass conversion, common caloric datum and typed static-pressure binding.
Its read-only model and liquid accessors expose those existing owners.

The resolved thermal field stores one numerical state, the exact borrowed
caloric law, the active gradient-source wrapper and shared immutable material
receipts. Binding shares the existing receipt array through Arc. Accepted
adiabatic diffusion first completes the numerical energy/entropy gates, then
resolves actual final inventories at the accepted temperature through the same
original liquid sources. Unchanged temperatures reuse the active wrapper.
Final density/domain failure or cancellation discards the complete candidate.

Mixture transport Cp remains an independently sourced property and is not
substituted into the pure-liquid caloric energy balance. This adapter does not
certify consistency between independently authored transport Cp and the caloric
pair. All fixed-volume, uniform-temperature and first-order thermal assumptions
of the numerical owner remain in force.

Three source integration regressions cover original piecewise-linear Cp with
analytic sensible heat, source identities and shared receipts, continued
evolution, rejection of an equivalent but different liquid-law instance, and
late density/domain/cancellation failures after a valid numerical thermal step.

## Same-energy binary interface tension in capillary flow

`capillary::material::miscibility::step_binary_liquid_planar_horizontal`
accepts the privately bound `ResolvedBinaryLiquidPlanarInterface`. It uses
that original interface's derived surface tension with its actual coexisting
material phases, then delegates to the existing open two-fluid Washburn
operator. The returned `BinaryLiquidPlanarMeniscusStep` retains the original
interface, actual/allowed width-to-radius ratios and the complete existing
flow, end-event, work and signed species/mass-transfer result.

The wall contact-angle cosine remains explicit. Before any flow evaluation,
the caller supplies a finite positive maximum characteristic-width/radius
ratio in (0,1]; an overly broad or unrepresentable interface refuses. This is
an explicit applicability check for a planar-to-curved sharp-meniscus
approximation, not an independent bound on curvature or contact-line errors.
The central-slope width is not the full extent of the profile's infinite
tails. All existing open-reservoir, frozen-coefficient and
pressure-independent-coexistence assumptions remain in force.

The original supplied-sigma API remains available. Neither path borrows a
liquid/gas surface or turns the calculated tension into measured material
authority. Focused integration regressions compare the entire derived-tension
step with the existing supplied-tension result in both directions and phase
orders, retain original inventory/owner identity, and reject thick interfaces
or invalid wall/phase controls. These continuum/material tests require the
normal runtime build.

## Source-owned pressure outlets in Cartesian compressible flow

`BoundaryCondition::PressureOutlet { pressure_pa }` supplies nonnegative finite
absolute backpressure to frozen and evolving-species Euler fields in two and
three dimensions. The original `FrozenGas` or `MixtureGas` characteristic owner
resolves outgoing compression, expansion, sonic fan sampling and supersonic
departure. Interior composition, entropy and tangential velocities follow the
departing gas. Pressure alone cannot supply incoming thermodynamics: reversal
refuses without publishing state, time, loads or a successful report.

Exactly zero backpressure denotes exterior vacuum. The existing source closure
must still admit the positive-pressure, positive-density sonic face, or retain
the interior trace when all waves already depart. No zero-density/zero-pressure
EOS state is constructed, and an unsupported sonic source interval still
refuses. Pressure alone supplies no incoming vacuum state. This does not
resolve an internal vacuum, the exterior expansion or rarefied flow.
The same rule applies to the legacy spatial `PressureOutletBoundary` adapter.

First-order and MUSCL paths use the same source closure. Every species, all
momentum components, pressure traction and full referenced energy share the
existing face flux and RK boundary impulse. No face-sized storage or per-step
allocation is added. Initial outgoing-wave bounds participate in Euler and
complete split intervals, with the existing stage CFL gates retained. The
planar hydrostatic path uses its physical reference face for this bound; a
resting high-pressure face is not mistaken for reversed flow from its lower
cell-average pressure.

The original `Euler3d::new_with_open_boundaries` and `PressureOutletBoundary`
API remain supported. The enum uses that same immutable outlet owner. Explicit
reservoir or outlet bindings still require unique transmissive sides and cannot
override a pressure condition. Frozen molecular open-face heat/stress policy is
unchanged; evolving-species molecular owners refuse pressure-only sides because
no species/thermal bath is specified. Backpressure defines neither a
nonreflecting boundary nor finite exterior storage.

G0/G1 regressions cover all four/six side orientations, both reconstruction
modes, steady advection, discharge and choking, absent species, shear and
physical-energy conservation, hydrostatic balance, legacy spatial API parity,
body forcing and complete split intervals. G4 covers invalid pressure,
conflicting reservoir bindings, reversal and cancellation without partial
publication. These are manufactured source/physics regressions, not an
experimental accuracy certificate.

## Source-owned damped Newton rheology

`MaterialFluxSystem::bdf1_rheology_newton_step` carries original material
transport and analytic `dmu/dgamma` [Pa s²] into the existing sparse constitutive
Jacobian and bounded Armijo Newton solver. Viscosity and its derivative use the
same reference density as inertia, pressure and force. Original cell/rate
responses are cached across volume and one-sided face quadrature, and the
accepted endpoint retains those owners, derivatives and molecular heat
`mu*gamma²`. No numerical differentiation or frozen-reference heat is used.

The native checked entry stops after an already-started assembly if a source
refuses, before another linear solve or candidate acceptance. A source refusal
during line search remains terminal and keeps the original material error;
it is not permission to extrapolate a source's domain. Invalid numerical trial
states retain the native backtracking behavior. Failed calls leave caller
history untouched. The native raw-coefficient API keeps its original behavior.

G1/G3 regressions use the existing stiff shear-thickening case that exhausts
Picard at the same budget, compare two source densities against the native
Newton result, and check pressure, kinetic energy, constitutive tangents and
local/integrated heat. The canonical thermally shifted Carreau-Yasuda source
is also exercised and compared with converged Picard. G4 covers actual trial
source refusal, changed density, invalid tangents and invalid timesteps.

This is constant-density incompressible BDF1 with velocity boundaries and
external endpoint temperature/composition held fixed. It supplies no evolving
thermal field, viscoelastic history, mixed pressure-boundary Newton solver or
unconditional convergence guarantee. The native smooth radial-law admission
requires zero scalar viscosity derivative at zero strain.

## Pressure-fed evolving species in three-dimensional flow

`spatial::mixture::MixtureEuler3d::new_with_reservoirs` binds immutable complete
gas traces to selected transmissive sides. Each trace retains canonical partial
mass densities, Cartesian velocity and absolute pressure. The original
`MixtureGas` admits its thermodynamic state and supplies composition-dependent
temperature, heat capacity, sound speed and referenced energy. First-order and
MUSCL boundaries use the existing spatial HLLC solve; each face's species,
momentum, pressure work and enthalpy enter the same conservative RK account.

The exterior three-axis signal bound participates before selecting either the
Euler interval or the complete chemistry/waves/chemistry interval. Prescribed
body forcing retains the same traces and its separate mechanical-work account.
The reacting trial clones the original binding once, includes both copies in
its explicit numerical-byte budget and retains it across accepted intervals.
No face-sized reservoir storage or per-step allocation is added.

G0/G1 regressions cover uniform stationary/moving mixtures, pressure-driven
species and physical-energy exchange on all six sides with both reconstruction
modes, original caloric recovery, hot-reservoir CFL, body forcing and reactive
persistence without coupling retries. G4 covers source/topology/byte admission
and late cancellation preserving accepted state, time, loads and reports.

These are prescribed infinite reservoirs, not finite inventories, pressure-only
or nonreflecting boundaries. The optional molecular owner uses those same
complete traces for source-owned species, heat and stress exchange with the
exterior, as specified below. It adds no pressure-equilibrium
guarantee for unequal-caloric moving contacts or fully discrete entropy theorem.


## Opt-in mixture contact-pressure equilibrium on Cartesian grids

The standalone `MixtureEuler2d` and `MixtureEuler3d` owners provide
`step_with_contact_equilibrium` with the original `fs-gas::multicomponent`
`pressure_equilibrium::{ContactEquilibriumControl, ContactWorkspace}`.
Both first-order and MUSCL transport retain the native species/momentum fluxes,
source thermodynamics, CFL admission and RK2 discretization. The same prepared
contact workspace is reused across steps; each owner admits its combined
numerical workspace against the existing hard storage bound.

Each cell's original pressure and full Cartesian velocity are checked against
both actual traces at every incident face. Physical boundary cells, HLL
fallback stencils and variations exceeding the caller's tolerance retain
the conservative update. The second stage can only remove eligibility.
All velocity components participate, so a normal-velocity match does not
erase a shear interface. Velocity tolerance is scaled by sound speed.

Eligible cells reconstruct energy through the original mixture calorics
after both forward-Euler updates and the final SSP average. All represented
momenta and transported species remain unchanged; zero species stay zero.
The correction includes every kinetic component and uses the physical cell
volume and RK weights 1/2, 1/2 and 1. The native arithmetic is preserved for
default stepping and cells without correction.

The generic `ContactEquilibriumStepReport` retains the native transport
report and adds the signed numerical energy defect, absolute-defect sum,
unique corrected-cell count and corrected energy residual. The original
energy residual still exposes the defect. The caller bounds the sum of
absolute corrections per accepted step, so opposite signs cannot cancel
the budget. A budget retry halves time without downgrading MUSCL.
Cancellation or refusal publishes no field, clock or load change.

This is an Estimated, explicitly nonconservative energy option for admitted
equilibrium contacts. It is not an entropy or general shock-accuracy theorem.
The default transport's moving unequal-caloric contact limitation remains.
Hydrostatic/body-force, molecular, reacting, moving-wall, rotor and network
compositions continue through their existing methods and do not enable
this contact correction implicitly.


## Face-resolved held exterior gas in frozen three-dimensional flow

`spatial::Euler3d::prepare_reservoir_profiles` binds a fixed, unique set of
otherwise unbound transmissive sides. Each physical face carries its own
`ReservoirFaceState` density, three Cartesian velocities and absolute pressure.
The original `FrozenGas` source constructs conservative energy and decodes the
same calorics as a native cell. Existing MUSCL reconstruction, HLLC/HLL fluxes,
SSP retries and per-face mass/momentum/enthalpy accounts remain the owners of
transport; the exterior trace is piecewise constant on each physical face.
Every cached exterior velocity and sound speed enters the unsplit initial CFL.
Uniform reservoirs and pressure-only outlets retain their original semantics
and cannot share a profile side.

Preparation admits active and candidate cache capacities against the explicit
additional profile byte budget before sampling. `update_reservoir_profiles`
visits each prepared face once in canonical side/face order and publishes all
sides together after source admission and cancellation checks. Updates allocate
no numerical arrays. Ordinary held stepping invokes no profile callbacks;
the explicit stage-time entry below observes its actual RK times instead.
The caller owns any stochastic field and its temporal model.
A failed update preserves active traces, cells, clock, reports and historical
boundary impulses. A successful update invalidates the old report without
changing the accepted cells, clock or historical impulses.

`MolecularEuler3d` retains both complete fields' profiles when composing an
already profiled Euler field. Its explicit prepare budget and workspace-byte
getter cover both copies; its update stages and synchronizes both candidates
before either publishes. Profiles survive accepted/trial ownership swaps and
share the original complete-interval retry. Transmissive molecular flux remains
zero by default; the explicit reservoir-exchange selection below enables
viscous and thermal bath coupling. These remain prescribed infinite exterior
states, without finite storage, nonreflecting behavior, nonlinear flux averaging
or a calibrated turbulence model.

Five focused G0/G1/G4 tests were authored: uniform bit-equivalence and canonical
sampling; analytical supersonic mass, all momentum components and enthalpy on
all six orientations; budget/source/conflict/cancellation transactions; the
hottest individual face's unsplit CFL; and evolving boundary controls through
molecular ownership swaps with zero-transport equivalence. These tests have not
yet executed in this environment; formatting passed. No new DNS, turbulence
spectrum, temporal convergence or experimental accuracy claim is made.


### Face-resolved held inflow for evolving-species three-dimensional gas

`MixtureEuler3d::prepare_reservoir_profiles` admits one complete primitive trace
per physical face on selected transmissive sides, in `BoundarySide::ALL` and
`CartesianGrid::face_center` order. The callback fills every canonical species
partial density and returns Cartesian velocity and absolute pressure. The
original `MixtureGas` supplies calorics, sensible-energy precision admission,
referenced energy and sound speed. Uniform reservoirs and pressure-only, wall
or periodic conditions cannot simultaneously own the same side.

The existing spatial HLLC/SSP path consumes each face's own exterior state.
Species, shear momentum, pressure work and enthalpy remain in the original
conservative flux and boundary-impulse account. The largest original-source
three-axis acoustic/advection rate across all prepared faces participates in
the initial interval selection; it is cached during profile admission rather
than rescanned at every step. No species clipping, renormalization or second
thermodynamic model is introduced.

`update_reservoir_profiles` prepares every candidate before one publication,
using fixed active/trial rows and primitive caches. Invalid source data, a
provider refusal or cancellation leaves all active traces, cells, clock and
loads unchanged. Both retained payloads count against the explicit additional
profile budget and the existing total numerical ceiling. Updates and native
face evaluation add no per-face numerical allocations. Uniform bindings retain
their previous storage and arithmetic.

The reacting owner clones both profile payloads under its complete field budget
and exposes the same replacement operation through a shared two-owner
transaction. The provider is evaluated once per face; the candidate is copied
into the existing trial buffers before either field publishes. Traces are held
through the full chemistry/transport split and all retries.

Five focused G0/G1/G4 tests cover uniform-path equivalence and canonical face
order, all-six-side supersonic species/shear/enthalpy fluxes with original
variable calorics, a single hot face's initial CFL bound, budget/source/topology
refusal and last-poll rollback, and reacting clone/update/advance. These tests
require the normal runtime build and were not executed in the constrained
editing environment.

These are prescribed infinite exterior traces, not finite inventories or
nonreflecting boundaries. Replacement does not advance a stochastic source or
interpolate its stages; callers author the sampling cadence and hold one
coherent realization. A face-averaged primitive is not an exact average of its
nonlinear Riemann flux. No turbulent-inflow spectrum, wall/DNS accuracy, or
experimental qualification is implied by the boundary adapter alone.


## Molecular exchange through volumetric mixture reservoir profiles

`spatial::mixture::molecular::MolecularTransport3d` admits uniform and
face-resolved complete reservoirs already bound to the original Euler owner.
Each molecular evaluation borrows those same accepted traces. The original
`SpatialTransport::reservoir_face_from_primitives` supplies Maxwell–Stefan
species flux, partial enthalpy transport, Fourier heat and all Newtonian
stress/work components at the physical half-cell normal distance. Midpoint
constitutive coefficients retain their existing source law; the prescribed
bath temperature and velocity supply boundary enthalpy and mechanical work.
No additional heat source, species normalization or finite bath inventory is
introduced.

The native nine-derivative cache uses velocity ghosts mirrored about each
physical reservoir trace. Tangential profile derivatives use centered
stencils, periodic wrapping or quadratic one-sided end stencils; two faces
use a linear derivative and one face resolves no tangential variation.
Interior and exterior derivatives enter the source midpoint approximation.
Local SGS therefore sees consistent adjacent-cell gradients, while modeled
SGS flux at the externally maintained reservoir remains zero. This does not
extend periodic/channel dynamic-fit admission to arbitrary open topologies.

SSP stages, stabilized RKL stages and existing finite-wall trials share this
same spatial operator and original boundary quadrature. No additional
numerical buffers or per-face allocations are needed. All complete exterior
traces remain held during an accepted interval and its retries. Unbound
transmissive ends and pressure-only outlets still refuse molecular evolution.

`MolecularEuler3d` clones and budgets complete reservoir state in its trial
field. Its `update_reservoir_profiles` reuses the original paired transaction:
one callback evaluation per face, original caloric admission, copying into
existing trial storage, then joint publication. Both molecular halves and the
middle HLLC solve use the same controls. Failed control replacement or a
cancelled flow step preserves accepted profiles, cells, clock and reports.
Prepare profile topology on the Euler field before constructing this owner.

Four focused regression groups were added for analytic six-face heat,
viscous traction and bath work; steady three-axis species/enthalpy flux with
both integrators; periodic and one-sided profile derivatives; and coupled
update/cancellation/replay. Rustfmt and source/API review passed. These new
full-owner tests were not executed in the resource-constrained editing
environment; no native integration or turbulence calibration result is claimed.
The boundary is a prescribed infinite bath, not a nonreflecting inlet,
characteristic wave-injection method or validated turbulent entry length.


## Coherent atmospheric turbulence at compressible reservoir faces

Feature `experimental-mann-inflow` exposes
`spatial::mann_inflow::MannReservoirInflow3` for the original frozen and
evolving-composition gas owners. It observes the original full-tensor
`fs-atmo` Mann realization at native physical face centres, with an explicit
proper field/grid rotation and a grid-corner position in metres. A frozen
snapshot can be convected to the actual fluid time; a fixed-step OU source
must already be at that exact accepted time. Applying a trace draws no noise
and advances neither source nor fluid. Caller-owned baseline callbacks supply
density or every canonical species density, absolute pressure and mean velocity;
only velocity receives the shared fluctuation. Each resulting complete state
passes the original gas source's caloric admission before publication.

Point observation allocates no gains. Optional rectangular face averaging uses
the original exact `MannBoxFilter`; opposite faces share one filter per normal
axis, and the combined actual capacities obey the declared byte limit. The
filter also remains available for original full-tensor covariance queries.
There is no RMS restoration, independent per-face noise, species normalization,
mean-flux correction or replacement thermodynamic law. Native active/trial
profile buffers supply atomic preparation and refresh. Every original face
remains distinct in the existing Riemann solve, CFL and conservation accounting.
The immutable `MannReservoirTrace3` can also feed composed molecular owners'
profile callbacks without duplicating the atmosphere or its phase state.

All three cell widths admit the represented wavenumber cutoff against spatial
Nyquist; transverse periodic axes refuse because arbitrary retained Mann modes
are not periodic on the CFD section. A checked face-count times mode-count
budget bounds one complete application. The declared maximum hold duration
bounds advected spatial phase; it is not a temporal bandlimit for OU noise.
The caller of `apply_*` must refresh the control within that duration. Those
published traces remain held through SSP stages and retries, so that interface
uses a first-order hold for varying boundary data. The explicit `step_*`
interface below instead samples the actual stage times. Face averaging is a linear velocity observation, not an exact average
of nonlinear mass, momentum or energy flux. These complete-state exterior
Riemann traces provide no nonreflecting-inlet, wall-layer, turbulent mixing,
resolved-DNS, measured-spectrum or instrument-acceptance qualification.

Focused fixtures cover original-source face states and conservative native
stepping, rotated noncubic face averages against direct spatial quadrature,
canonical species preservation, transactional late refusal, source clock/model
mismatch and sampling budgets. Three exact production geometry-kernel tests
passed natively (frame/phase mapping, rectangular integration, overflow-safe
advective speed); the four full consumer fixtures remain unexecuted in this
environment, where DSR/RCH and the complete constellation build are unavailable.


## Stage-time prescribed mixture inflow in three dimensions

`spatial::mixture::MixtureEuler3d::step_with_reservoir_profiles` supplies complete
canonical-species exterior traces at the native SSP RK2 times. With boundary
trace `b(t)`, the original finite-volume residual advances
`U* = U^n + dt L(U^n,b(t_n))` and
`U^(n+1) = [U^n + U* + dt L(U*,b(t_n+dt))]/2`.
Physical face impulses retain `A dt (F_n + F_*)/2` for every species, all
momentum components and source-referenced total energy. The original mixture
EOS, sensible-energy admission, pressure/domain limiter and HLLC/HLL owner
remain responsible for thermodynamics and transport.

The initial prescribed trace bounds the interval; each stage's exterior
three-axis acoustic/advection rate and actual Riemann signals pass the native
CFL gates. A shortened retry resamples the start and its own endpoint. The
explicit `maximum_face_evaluations` counts every provider call across all
attempts, including the initial sample used for interval selection. Provider
refusal, exterior-source admission failure and exhausted sampling are terminal;
numerical trial failures retain the native halving/reconstruction policy.

The provider must be pure and replayable at absolute time in seconds. No source
clock, stochastic realization or exterior inventory is advanced. Existing
admitted active/trial profile capacities are reused between residuals, with no
new numerical arrays or per-face allocation. Endpoint profiles publish only
with accepted cells, clock, boundary impulses and report after the final
cancellation poll. Any refusal preserves all those accepted values. Ordinary
held stepping and composed molecular/reacting owners keep their existing
semantics; this entry does not silently enable forcing or contact-energy
correction.

Five focused G0/G1/G4 suites cover held-path parity, canonical absolute-time
sampling, six-side time-varying species/shear/pressure/enthalpy, original variable
calorics, analytic second-order incoming species integration, shortened CFL
retries and cumulative sampling, source refusal, and last-poll rollback/replay.
All five suites pass in the focused native-source build described below;
pinned rustfmt parsing and formatting also pass. Second-order temporal
behavior is scoped to smooth prescribed data; no
nonreflecting boundary, evolving OU interpolation, DNS, turbulent-spectrum or
experimental validation is asserted.


## Stage-time prescribed profiles for native frozen three-dimensional Euler

`spatial::Euler3d::step_with_reservoir_profiles` advances the original HLLC/HLL,
MUSCL and SSP RK2 owner with complete exterior data at physical times `t` and
`t + dt`. Profile geometry is prepared through the existing explicit byte
admission. The first residual consumes the source-admitted start trace; only
then is the same trial array reused for the endpoint trace and second residual.
The original frozen gas supplies all conserved energy, calorics and sound
speeds. The integrated per-face exchange retains both RK weights for mass,
all Cartesian momenta and referenced energy.

The actual endpoint exterior's unsplit acoustic/advection bound and both
native face-wave CFL bounds must admit the proposed duration. A CFL retry
halves that duration without dropping reconstruction order, then samples the
start and the newly shortened endpoint again. Native cell-admissibility
failures retain the existing first-order retry rule. Sampling/source refusals
and exhausted callback budgets are terminal, not numerical retry requests.

`maximum_face_evaluations` is one cumulative budget for the whole call,
including rejected attempts. Every complete scan must fit before it starts;
a callback which refuses still consumes its evaluation. Samples visit all
prepared faces in canonical side/face order. The callback must be a pure,
deterministic observation of a read-only replayable source. It must not advance
an RNG, OU clock or external solver. Source time validity and any interpolation
model remain explicit caller responsibilities.

The existing active/trial profiles are the only face-sized storage. No profile
clone, new numerical array, replacement EOS or stochastic source is created
inside stepping. Uniform reservoirs and pressure-only outlets retain their
bindings. Default held stepping retains its original stage arithmetic.
Accepted profiles remain unchanged during every trial; the accepted endpoint
profile is published only at the same final commit as cells, clock, report and
boundary impulses. Provider failure, exhausted retries, failed admission and
late cancellation preserve all of that accepted state. Subsequent held
stepping uses the newly accepted endpoint profile.

Four focused G0/G1/G4 regression groups exercise constant-source equality of
native field and boundary-account bits; second-order integrated supersonic
inflow momentum/energy against closed-form time-dependent flux integrals with
original variable calorics and both reconstructions; actual endpoint-CFL
halving with exact repeated sample times and a call-wide work budget; and
provider/source refusal plus final-poll cancellation preserving a previously
accepted complete transaction. All four regression groups pass in the focused native-source build described
below. Pinned rustfmt parsing and formatting also pass; a complete constellation
or DSR/RCH build was not run.

The temporal convergence target concerns smooth prescribed boundary data and
its physical Euler flux, not a calibrated turbulence spectrum, DNS resolution,
nonreflecting inlet or experimentally qualified flow. This native Euler entry
point does not implicitly retime forcing, molecular, reacting or evolving
stochastic-source compositions. Those owners continue using their documented
held-source methods until an explicit stage-aware composition is supplied.


## Frozen-gas molecular exchange with prescribed volumetric reservoirs

`spatial::molecular::MolecularTransport3d::with_reservoir_exchange` and the
corresponding `MolecularEuler3d` builder explicitly enable heat and traction
through complete uniform or face-resolved reservoirs already admitted by the
original Euler owner. The default remains zero molecular flux at open ends.
Unbound transmissive ends and pressure-only outlets retain that original
zero-flux behavior even when reservoir exchange is selected. No incomplete
exterior thermodynamic state is inferred.

Every molecular stage borrows the original cached source-admitted bath trace.
The same original transport law evaluates the cell/bath midpoint temperature,
pressure and frozen composition. The normal stress, both shear tractions and
Fourier heat use the physical half-cell normal distance. Total-energy exchange
uses the prescribed bath velocity for all mechanical work; no separate viscous
heating source or finite bath storage is introduced. Source coefficient maxima
from these faces enter the original conservative parabolic rate, whose factor
of four includes doubled half-cell boundary conductance.

Bath velocity ghosts mirror the boundary trace into the existing cell-gradient
cache. Tangential profile derivatives use centered differences, periodic
wrapping or quadratic one-sided end differences; two faces use a linear
derivative and a single face resolves no variation. The original midpoint
stress approximation combines these with adjacent-cell tangential derivatives.
Local SGS sees the consistent adjacent gradients. Modeled flux at the
externally maintained bath remains zero unless the separate local SGS bath
selection below is enabled. Periodic/channel dynamic-fit admission is
unchanged. No additional numerical buffers or per-face allocations are
introduced.

SSP RK2, native RKL2 (including retry anchors and endpoint admission), both
coupled molecular halves and forced/stabilized coupled steps share this same
spatial operator and boundary quadrature. The selection survives accepted/trial
owner swaps. Complete bath profiles are held through each requested molecular
interval and its retries; the timed native Euler inlet API does not make these
molecular split stages sample a varying bath. Source refusal, work limits and
cancellation publish no partial field or boundary impulse.

Seven authored regressions exercise all six face signs and anisotropic spacings,
uniform/profile agreement, exact source midpoint sampling, hot-bath transport
rate bounds, quadratic tangential profile stress, uniform equilibrium, original
unbound behavior, SSP/RKL conservation and late cancellation/replay, forced
owner swaps, and active interior SGS with unchanged molecular bath flux. These
are source/stencil/balance fixtures, not gas-specific turbulence, DNS, wall-law,
nonreflecting-boundary or experimental validation. Composition remains frozen.


### Opt-in local SGS exchange with frozen-gas reservoirs

`with_local_sgs_reservoir_exchange` on `MolecularTransport3d` and
`MolecularEuler3d` explicitly extends an authored local Smagorinsky, WALE or
Vreman closure to complete source-admitted uniform/profile reservoirs. It
also enables molecular reservoir exchange. Molecular-only bath exchange
remains the default when only `with_reservoir_exchange` is selected.

The full constitutive tensor uses the existing half-cell normal differences
and midpoint cell/profile tangential derivatives. Authored physical filter
widths are unchanged, including directional Vreman widths. Face density is
the cell/bath average; heat transport uses the original caloric source at
midpoint temperature and pressure. Bulk viscosity stays molecular and
total-energy stress work uses the prescribed bath velocity. The augmented
coefficients enter the original conservative face flux and diffusion bound;
excess eddy-to-molecular viscosity refuses rather than clips.

This selection adds no numerical arrays and follows the original held or
physical-stage reservoir sampling in standalone and coupled SSP/RKL paths.
Missing local authorship and dynamic/Lagrangian fits refuse: exterior filter
samples and histories are not inferred. Unbound sides and pressure-only
outlets retain zero bath diffusion. The selection adds neither inflow
fluctuations nor a DNS, entropy or shock-interaction accuracy claim.

Five consumer regression groups cover six profile faces with all three
closures, variable density/source heat capacity, zero-model trajectory
identity, exact augmented diffusion bounds, ratio refusal, closure admission,
and complete SSP/RKL conservation/cancellation/replay. These sources were
parsed and reviewed but not compiled or executed. Separately, three tests of
the exact production `reservoir_face_gradient` and original `tangents`
functions passed under the pinned Rust compiler with warnings denied:
full affine tensors on all normals and half-cell distances, rigid rotation,
and common velocity boosts/uniform comotion. That isolated mathematical
check does not establish a full native consumer build or integration gate.

## Stage-time coherent Mann forcing of native compressible flow

`MannReservoirInflow3::step_frozen` and `step_mixture` advance the original
Euler owners with a frozen, coherently advected Mann realization observed at
the actual SSP RK2 start and endpoint times. Shortened numerical retries
resample their own endpoint; the atmosphere remains immutable and draws no
noise. Caller baselines receive the same absolute time and own all density,
canonical composition, pressure and mean-velocity data. Existing native EOS,
Riemann fluxes, reconstruction, boundary quadrature and retry admission remain
in control. Accepted endpoint profiles publish with cells, time and loads.

The same prepared spatial filters serve every stage. One stack-held trace
caches frame, topology and phase admission for each sampled time; no numerical
buffer is allocated in either step wrapper. `maximum_mode_samples` bounds the
complete call, including both stages and rejected attempts. The report exposes
actual face evaluations and represented-mode work. The adapter clips the
requested maximum duration to the frozen spatial field's conservative advective
phase limit; the native face CFL may shorten it further. Fixed-step OU sources
refuse this interface because an accepted endpoint does not define a continuous
stochastic bridge at arbitrary adaptive stage times. Their held-trace interface
remains available. Source callbacks must be deterministic and replayable.

Four focused tests cover independent rotated/advected sampling against the
native stage-time owner, conservation and source immutability, late endpoint
rollback and explicit OU refusal, complete-step sampling budgets, and canonical
mixture profiles at the accepted time. All four tests pass in the focused native-source build; pinned rustfmt parsing
also passes. This change
removes the held boundary approximation for this frozen-source Euler path; it
does not claim general spatial/temporal convergence, retime molecular or
reactive split stages, or qualify a DNS, nonreflecting boundary or measured
turbulence spectrum.

### Focused native validation of the reservoir extensions

Twenty focused regression tests pass: four frozen stage-time, five mixture
stage-time, six molecular bath, four Mann stage-time, and one SGS bath test.
The original `fs-gas`, `fs-rand` and `fs-atmo` production sources and the
complete compressible source tree were compiled with the pinned nightly,
including a second build with `experimental-sgs`. The test dependency slice
uses original thermal source and solver numerical kernels with unused executor
adapters excluded. For SGS it exports the original model, Vreman and dynamic
modules unchanged. No substitute physics or mocked physical owners are used.
This is focused native-source evidence, not a complete Cargo/constellation or
DSR/RCH gate, performance benchmark or experimental turbulence calibration.

The native check also verified the held-profile bridge in forced molecular
composition after the internal stage-time argument was introduced. The fresh
exterior CFL fixture uses a tangential speed jump to isolate its unsplit bound;
it does not suppress stronger normal-shock wave speeds or require those shocks
to bypass the existing retry policy.

## Prescribed physical-time baths in frozen volumetric molecular transport

`spatial::molecular::MolecularTransport3d::advance_with_reservoir_profiles`
and `advance_stabilized_with_reservoir_profiles` observe prepared exterior
profiles at the original SSP or RKL2 physical stage times. The explicit
`physical_start: Time` supplies the origin because a standalone molecular
operator does not advance the Euler clock. Existing held entry points retain
their original arithmetic. Molecular bath exchange still requires the
existing `with_reservoir_exchange` selection.

Each SSP substep samples its start and endpoint; a rejected trial rebuilds
its start trace and samples the shortened endpoint. Native RKL2 reuses its
original coefficients, anchor flux and signed boundary recurrence. Its state
abscissae are `c_0=0`, `c_1=4/(3D)`, and `c_j=(j*j+j-2)/D` for `j>=2`, with
`D=s*s+s-2`; the residual forming stage `j` observes state `j-1`, and final
admission observes the exact endpoint. The dimension-independent original
`Rkl2Scheme::abscissa` is shared by both frozen and mixture owners. No
replacement temporal or spatial transport kernel is introduced.

`MolecularEuler3d::step_with_reservoir_profiles` and its stabilized counterpart
compose the first molecular half on `[t,t+dt/2]`, Euler on `[t,t+dt]`, and the
second molecular half on `[t+dt/2,t+dt]`. A shortened central Euler step retries
the complete split from accepted cells. Split clocks retain the original
origin and phase offsets: the final source observation is the canonical
`t+dt` used by the native flow clock, including nonzero origins where adding
two separately rounded half-step clocks would leave the admitted interval.

The original gas, original transport law and its actual midpoint source
admission own every heat/stress/work evaluation. All six bath normals,
half-cell conductances, tangential profile derivatives, cell-gradient ghosts,
wall conditions and selected interior SGS follow the same immutable active or
staged reservoir view. Bath faces retain molecular-only exchange unless the
explicit local SGS reservoir policy is enabled; that policy observes the same
stage trace. Density/composition/caloric ownership does not change.

The provider must be a replayable read-only observation: identical arguments
return identical physical data, including repeated or nonmonotone time
visits. One cumulative face-evaluation budget covers initial exterior CFL
admission, both molecular halves, central Euler, accepted substeps and all
discarded retries. Each full scan must fit before it starts. A provider,
exterior-source or face-budget error terminates immediately with its original
error; native numerical cell/CFL failures retain bounded retries. RKL's
operator-evaluation budget remains independent and cumulative across the
complete coupled call.

The existing active/trial profile arrays remain the only profile storage;
time-aware sampling adds no numerical buffers or per-stage allocation. Source admission writes private
staged data. A complete standalone molecular interval publishes cells, its
endpoint profile and molecular loads together, leaving the Euler clock and
advection metadata unchanged. A complete coupled interval publishes the full
accepted owner, time, profile, combined loads and report. Failure, exhausted
work or final-poll cancellation preserves all accepted data. Before a later
coupled call, the private owner restores accepted profiles, so held stepping
after a failed or successful timed call uses the accepted endpoint.

Focused regressions cover held-bit compatibility, analytic oscillating
six-face Fourier response, manufactured nonautonomous isothermal shear with
bath work, actual accelerated RKL abscissae, native temperature-dependent-law
retry sampling, complete split replay after central CFL shortening, shared
provider budgets, exact nonzero-epoch endpoints, late source refusal and
cancellation, and subsequent held owner swaps. Their time-refinement targets
concern the original discrete transport operators with smooth prescribed
boundary data. They do not qualify a measured turbulence spectrum, a DNS,
nonreflecting boundaries, wall laws or experimental material data. These APIs
do not implicitly retime reacting chemistry, body-force compositions or an
evolving stochastic-source owner; a separately admitted replayable source
adapter must supply those observations.



## Stage-time molecular reservoirs for native volumetric mixtures

The explicit `spatial::mixture::molecular::MolecularTransport3d` methods
`advance_with_reservoir_profiles` and
`advance_stabilized_with_reservoir_profiles` sample complete prescribed
reservoir histories inside the original molecular stages. `physical_start`
supplies absolute seconds because this standalone operator leaves the Euler
clock unchanged. SSP uses each actual substep start and endpoint, including
shortened retries. RKL2 uses the original Legendre recurrence abscissae and
qualifies the accepted endpoint with an additional source evaluation.

The original Maxwell--Stefan face owner, source calorics, every partial
enthalpy, Fourier heat, full Newtonian stress and prescribed-bath work remain
authoritative. One bank selection reaches both mirrored velocity ghosts and
tangential reservoir derivatives, so any admitted SGS model sees consistent
stage gradients. The existing externally maintained reservoir SGS flux remains
zero. Native conservative divergence and signed boundary quadrature are
unchanged; there is no additional viscous-heating source or species clipping.

The composing `MolecularEuler3d` methods `step_with_reservoir_profiles` and
`step_stabilized_with_reservoir_profiles` use physical intervals
`[t,t+h/2]`, `[t,t+h]`, and `[t+h/2,t+h]` for first molecular half, Euler, and
second molecular half respectively. Phase clocks retain the original origin
and offsets, with explicit canonical phase endpoints. Consequently the final
profile timestamp is the accepted `t+h` bit pattern even when separately
adding two rounded half-step origins would exceed that endpoint. No provider
query is clipped or widened to conceal an out-of-range stochastic sample.

One existing reservoir-stage context counts every supplied face across initial
interval selection, all molecular substeps and recurrence stages, both split
halves, Euler stages, and rejected complete intervals. Sampling budget, provider
and exterior-profile source-admission failures terminate immediately. Original
numerical trial failures retain their native retry policy. RKL2's independent
operator-evaluation budget continues to cover both halves and discarded probes.
Providers must be pure and replayable at absolute time; mutable closure access
is for observation scratch, not advancing external inventory or randomness.

Prepared active/trial profile banks are reused without new numerical arrays.
The ordinary opt-in RKL buffers retain their existing preparation budget.
Standalone success publishes endpoint profiles with cells, molecular report
and molecular loads; complete split success also publishes the clock, Euler
report and all composed loads. Refusal or final cancellation preserves these
accepted values. Each complete retry restores the trial active profiles from
the accepted bank, including a subsequent ordinary held step after a failed
time-dependent call. Prepare reservoir topology before constructing the
composed owner. Ordinary held methods retain their existing semantics.

Eight focused native tests cover held-operator parity with variable calorics;
second-order sinusoidal species and heat reservoirs against exact one-cell
solutions in SSP and RKL2; original diffusive enthalpy transfer and energy
conservation; multistage Legendre sample times; physical Strang timing and
canonical endpoint bits; cumulative work through outer retries; terminal late
provider, source and budget failures; and final cancellation with complete
rollback and replay. Finite wall inertia is refused before provider work.
These new APIs do not add time-dependent reactive, adaptive, finite-inertia,
body-forcing or wall-acceleration compositions. Smooth prescribed-data order
and native conservation do not establish nonlinear RKL stability, a turbulent
spectrum, nonreflecting inflow or experimental validation.

## Conditional Mann OU paths in native compressible stage consumers

`MannReservoirInflow3::step_frozen_bridge` and `step_mixture_bridge` observe
one prepared `fs-atmo::MannOuBridge` at the native Euler SSP stage times.
`step_molecular_frozen_bridge` and `step_molecular_mixture_bridge`, with their
`step_stabilized_molecular_*_bridge` counterparts, observe the same path
through the original molecular/Euler/molecular composition and RKL2 stages.
The original atmosphere owns the full tensor, temporal correlations, shared
conditional innovations and exact fixed-clock endpoints. The original gas
owners retain EOS admission, canonical species, Maxwell--Stefan transport,
partial enthalpies, Fourier heat, complete stress/work and conservative fluxes.
The adapter only adds the shared velocity fluctuation to a replayable authored
baseline; it does not synthesize density, species, pressure or temperature noise.

`MannReservoirSource3::BridgeFrame` admits one held observation for existing
profile preparation. The dedicated bridge steps then reuse a caller-prepared
frame. One update supplies every consecutive face at the same stage time;
nonmonotone split stages and retries deterministically reconstruct the same
path. Point/frame rotation and exact rectangular face filters use the same
observation implementation as existing held/frozen inlet consumers. No
successful adapter stage allocates numerical storage, and no per-face random
draws occur. Original native lazy RKL workspaces retain their separate budgets.

The complete call charges two independent turbulence work limits:
`maximum_mode_samples` counts every face times represented modes, and
`maximum_normal_draws` counts every conditional frame reconstruction, including
discarded trials. `MannReservoirBridgeReport3` reports actual face/mode work,
normal draws, frame updates and the admitted maximum duration beside the
original native conservation report. Preparing the right OU endpoint has its
own explicit atmosphere budget. A source or draw-budget refusal terminates the
native transaction without advancing accepted fluid state or the original
stochastic owner. The mutable frame is observation scratch and may retain its
last successful observation after flow failure; querying again replays the same
path. It is not an additional physical source clock.

Each step is capped by the retained right endpoint and by the advected spatial
phase bound. OU fluctuations have no temporal bandlimit; this cap is not a
stochastic accuracy estimate. Flow time must start inside the retained interval,
and another interval must be explicitly prepared when its right endpoint is
reached. Neither the adapter nor its callbacks advance the fixed-step source.
Canonical molecular phase clocks ensure that the last bath observation is the
same represented physical endpoint as the accepted native flow, including
nonzero epochs where adding two rounded half-step origins would overshoot.

Six native consumer regressions exercise direct original-source filtered
synthesis against frozen and mixture Euler results, canonical species and
enthalpy balances, coupled SSP/RKL bath conservation, exhausted-draw rollback
and replay, source immutability, interval clipping, and the represented
`6e-6 -> 7e-6` split endpoint for both gas families and both integrators.
The original held/frozen Mann tests and molecular/SGS bath regressions also
remain passing. Validation uses selected original-source native builds;
it does not substitute for Cargo/constellation gates, measured atmospheric
calibration, mesh-resolution studies or turbulent-flow validation. The exact
stochastic path is the admitted finite-mode Gaussian OU model. Existing
independently keyed interval-average observations are not jointly coupled to
this path. Timed reacting, body-force and finite-wall-inertia compositions
remain outside these entry points.

### Held forcing with stage-time Cartesian exterior profiles

`Euler3d::step_with_forcing_and_reservoir_profiles` composes the original conservative source half-kicks around native SSP Euler, observing a replayable prepared exterior at the actual accepted start and endpoint of each trial. Uniform/per-cell acceleration and force density retain their distinct units and original density-dependent impulse and represented midpoint kinetic work. Source impulse/work and boundary transfer remain separate terms in the original conservation balance. Timed RK input, predictor and accepted candidate states, the final force kick, and endpoint CFL retain authored native admission.

`MolecularEuler3d::step_with_forcing_and_reservoir_profiles` and `step_stabilized_with_forcing_and_reservoir_profiles` connect the same forcing to original SSP or RKL2 molecular halves. First-half, central Euler and second-half physical intervals remain `[t,t+dt/2]`, `[t,t+dt]` and `[t+dt/2,t+dt]`, using the original canonical origin/offset clock. A shortened central interval replays the complete molecular/force composition from accepted cells. Original Newtonian stress/work, Fourier heat, configured wall/bath and optional exploratory SGS owners remain active.

One provider-face budget covers initial exterior admission, all phases and discarded numerical trials; stabilized operator work remains separately cumulative. Providers must be replayable under repeated and nonmonotone observations. Provider/budget failure terminates without numerical retry; original bounded numerical retries remain intact. Cells, time, endpoint profiles, face loads and reports publish only after complete acceptance. Failed timed calls preserve subsequent held-force recovery, and a zero force delegates to the existing timed path exactly. Existing caller-budgeted forcing scratch and reservoir caches are reused without per-stage allocation. Body forcing remains authored and held over the accepted interval; this adds no stochastic volume forcing, well-balanced gravity scheme or turbulence-accuracy certificate.

Focused native coverage reconstructs force/advection/force independently from kinetic energy, reconstructs both original molecular halves and central forced Euler (including active RKL stage acceleration), checks conservative source/boundary work and exact zero-force parity, and exercises Mach-driven whole replay, cumulative provider exhaustion, late source/cancellation rollback and held recovery.

## Held body forcing with stage-time volumetric mixture reservoirs

`spatial::mixture::MixtureEuler3d::step_with_forcing_and_reservoir_profiles`
composes the existing held `BodyForce3d` with complete prescribed exterior gas
at the original SSP RK2 stage times. The original force half-kicks surround
those residuals. Uniform and per-cell acceleration use the compensated sum of
canonical partial densities at each actual kick; force density remains per
volume. The exact momentum impulse and kinetic work use the original fixed-order
compensated source ledger. No species change, formation-reference shift,
additional heating source or stochastic volume-force model is introduced.

The timed residuals now carry the original force control through every cell,
reconstructed face and updated RK stage. Final kick admission retains the
Mach and combined three-axis acoustic CFL bounds, alongside exterior and
Riemann-wave checks. Numerical failures restart both half-kicks from accepted
cells at the shortened interval. Complete exterior profiles are then observed
at the retry's actual stage times. A single face-provider budget includes
initial admission and every actual provider call in discarded attempts.
Provider refusal, exterior-source admission errors and budget exhaustion are
terminal; they never become additional numerical retries.

The corresponding `MolecularEuler3d` methods
`step_with_forcing_and_reservoir_profiles` and
`step_stabilized_with_forcing_and_reservoir_profiles` place the same forced
Euler transaction between the original molecular halves. Molecular phase
clocks remain `[t,t+h/2]` and `[t+h/2,t+h]`; central Euler stages remain
`[t,t+h]`. One shared profile context and its non-resetting budget cover all
three phases and complete-interval retries. Native SSP substeps, RKL2
abscissae, source calorics, Maxwell--Stefan transport, partial enthalpies,
Fourier heat, full stress and any already configured SGS operator remain
owned by their existing implementations. RKL2's separate operator budget
still counts discarded residuals. Final coupled force admission includes the
second molecular half, and mechanical source work remains separate from
boundary exchange in both native and complete-split balances.

A zero force retains the original timed arithmetic and reports. All supplied
forcing scratch remains subject to the original complete-owner byte budget,
even when the force is zero. The new compositions add no numerical buffers
beyond existing prepared profile banks, forcing scratch and opt-in RKL work.
Only a complete success publishes accepted cells, physical clock, endpoint
profiles, native and composed reports, and every boundary load. Late refusal
or cancellation preserves all accepted quantities; a subsequent call restores
and reuses the existing trial banks. Providers remain pure, replayable absolute
time observations; caller-side effects are outside the transaction.

Focused regressions compare the combined native step against independent
analytic force impulses and kinetic-work differences around the original
timed residuals, including spatial force-density and acceleration fields.
They exercise canonical absent species, variable-caloric energy conservation,
zero-force and constant-history parity, native and complete-split shortening
with deterministic replay, one cumulative face budget, SSP/RKL physical
sample times, and terminal late source/provider/budget or final-cancellation
rollback with subsequent successful reuse. These manufactured checks do not
establish turbulent spectra or experimental validation. Finite wall inertia,
wall acceleration and adaptive molecular compositions retain their existing
explicit limits; held-profile APIs retain their existing semantics.

### Time-resolved reacting spatial reservoir transport

`ReactingEuler3d` can now consume prepared prescribed reservoir histories through
its original chemistry/waves/chemistry transaction, with optional source-owned
molecular outer halves using SSP or RKL2. Each of those three paths also admits
the original held `BodyForce3d` and forcing workspace. The provider fills every
canonical species density and supplies velocity and absolute pressure; the
existing gas, kinetic mechanism, molecular law, energy references and physical
boundary operators remain authoritative.

One prescribed-face evaluation allowance covers initial exterior CFL admission,
both molecular halves, the central native Euler stages and every discarded
whole-chemistry trial. Chemistry-driven acoustic shortening replays the same
prescribed history at the shortened physical times. Provider errors, exterior
source admission errors and provider work exhaustion are terminal. Original
RKL2 work accounting separately spans both halves and all discarded trials.

The physical molecular intervals are `[t,t+dt/2]` and `[t+dt/2,t+dt]`; the wave
residuals sample `t` and `t+dt`. All split times retain the original outer clock
and relative offsets, so the accepted endpoint has exactly the original `t+dt`
bits even where adding two rounded half increments would exceed a bounded
prescribed source interval. Local chemistry retains its two original
fixed-volume half intervals and does not independently advance physical time.

Accepted cells, physical time, endpoint profiles, complete boundary loads,
body impulse/work, chemistry changes, formation/reference-energy balances and
both owners' reports publish together only after the complete split succeeds.
The private field restores accepted profile banks before every whole trial and
later held-profile calls, including recovery after late source failures or final
cancellation. Provider mutation of an external accepted source is outside this
transaction: prescribed callbacks must be replayable, and may use mutable access
only for observation scratch or work counters.

Finite wall inertia is refused before provider or chemistry work on these timed
entry points; wall acceleration is also unsupported. Existing held wall and
reactive APIs retain their original behavior. This extension does not synthesize
species/temperature fluctuations, change the kinetic or molecular model, or
claim fully coupled stochastic-source evolution, resolved detonation thickness,
or a nonlinear positivity/entropy guarantee from RKL2 stabilization.

Focused native tests cover zero-rate and nonzero-rate constant-history parity,
all six public combinations, changing composition/pressure with nonzero chemical
heat and body work, original RKL2 stages above two, heat-release whole-trial
replay and cumulative source work, physical split stages and rounded endpoints,
terminal late provider/admission/budget failures, final cancellation, complete
accepted-state rollback, later held recovery and early finite-inertia refusal.

### Evolving Mann boundaries with mechanical forcing and reacting transport

The experimental `MannReservoirInflow3` bridge adapters now reach the original
forced frozen and multicomponent Euler transactions, including their SSP and
RKL2 molecular compositions. They also reach original finite-rate reacting
Euler, reacting SSP molecular transport, and reacting RKL2 transport, each
with an optional held native mechanical force. `MannMixtureForce3` carries
the original force, authored Mach control and caller-owned reusable scratch;
it adds no source law, independent clock or numerical workspace.

All adapters reuse the original conditional OU bridge observer and native
time-resolved profile transactions. One retained stochastic realization is
sampled at actual native stage times, including nonmonotone split observations
and shortened complete-interval retries. The same frame is reused across all
faces at one consecutive time. Exact rectangular face filtering, the original
coordinate rotation, source identity checks and the advective phase bound
remain active. The accepted `FixedStepMannField` is never advanced by these
observations; the caller owns its eventual endpoint advance.

One cumulative face/mode allowance and one conditional-normal allowance span
the entire call, including discarded trials. Native forcing, chemistry and
molecular work retain their existing separate budgets. Source or bridge work
failure preserves every accepted physical field, clock, profile, report and
boundary load, including the independently supplied molecular operator. The
mutable bridge frame is reusable observation scratch, not an accepted physical
state, and may retain a failed trial's last successful observation. Bridge and
atmosphere errors retain their original values; reacting errors are carried by
`MannReservoirError3::Reactive` without replacing their native cause.

These paths add the Mann velocity fluctuations to the caller's complete mean
exterior state. Species, temperature/composition histories, kinetics, calorics,
Maxwell--Stefan diffusion, Fourier heat, full stress, force impulse and mechanical
work stay with their native owners. The force is held over the accepted interval;
this does not create stochastic volume forcing or a joint scalar turbulence
model. Existing timed wall-inertia and acceleration exclusions still apply.
Stabilized diffusion enlarges the explicit stability interval and does not
certify stochastic temporal accuracy or physical turbulence spectra.

Focused native tests compare the adapters to independent rotated/filtered
bridge synthesis supplied directly to the original forced and reacting owners,
with nonzero chemical heat and mechanical work. Both SSP and RKL2 compositions
exercise late draw-budget rollback, identical subsequent replay, both owners'
conservative loads and source immutability. Reacting molecular tests include
the nonzero-epoch rounded endpoint case, so split clocks cannot observe beyond
the retained source interval. These are manufactured numerical composition
checks, not an experimental validation or a full workspace acceptance gate.

### Density and inventory continuation across prepared meshes

With `continuum`, `ns::transport::density::remap` adds
`PreparedDensityTransport::remap_from` for a receiving native transport plan.
It reuses `fs_mesh::conservative_cell_remap` for physical volume, mass, and each
extensive `mass*specific` quantity through the same canonical source-volume
fractions. Actual source/target plan volumes determine cell mass and reconstructed
density; specific fields are reconstructed from received extensive inventory and
received mass. Source fractions must cover every source, and compensated local
checks require every target's overlap volume to match its native physical volume.
The caller supplies overlap geometry and the separately admitted target advector.

Volume [m³], mass [kg], and each carried inventory have distinct absolute
conservation allowances. Final residuals accumulate individual source and actual
reconstructed target amounts with compensation and carried-product FMA remainders,
preserving small defects that subtracting rounded endpoint totals can conceal.
Positive representable masses
and densities are required; signed inventories and `N=0` are supported. There is
no clipping, mass repair, species renormalization, or numerical loss converted to
heat. Phase-specific errors retain original scalar-owner refusals. Adapter vector
payload is admitted separately from the original remapper's existing per-call
caps; allocations are fallible and cancellation is polled through all scans.
No input state or plan changes before a complete candidate is returned.

The returned cells continue through the original target `advance_into` or
density-subcycling owner. This is prescribed-advector P0 material transport,
not a BDM velocity/pressure projection, coupled Navier–Stokes checkpoint,
angular/kinetic conservation claim, geometric intersection generator, or
constitutive-history transfer. Focused `density_transport/remap.rs` regressions
advance, refine, advance, coarsen and continue enthalpy/species-like inventories;
they also cover extreme density contrast, small signed defects, physical-volume
mismatch, cancellation and storage refusal.

## Gravity from conservatively transported density

`assemble_conservative_momentum_bdf1_with_gravity` couples physical Cartesian
acceleration [m/s²] to the density accepted by the original mass-transport
operation. The endpoint force is `rho_new * g`; the same positive quadrature
weights already used by endpoint inertia integrate it into the native BDM1
momentum load. Additional force density remains a separate [N/m³] callback.
There is no point-to-cell search, additional source query, separate density
prediction, or extra cell traversal.

Pressure retains its original total-traction datum and sealed-component gauge.
This is full variable-density gravity, with no Boussinesq linearization or
implicit subtraction of hydrostatic pressure. Open boundaries must supply their
physical hydrostatic traction. The original pressure-robust operator balances
resolved conservative gravity loads; lateral density contrasts can drive
buoyant circulation. Zero acceleration retains the original arithmetic path.

The original first-order IMEX, fixed-mesh, density-CFL, solver and cancellation
limits apply. This entry does not add an EOS, heat transport, variable gravity,
an energy-conserving gravitational-work discretization or support reactions.
Invalid acceleration is refused before callbacks; a failed/cancelled candidate
leaves both histories unchanged.

The focused G1/G3 regressions in `conservative_density_momentum/gravity.rs`
cover transported-density free fall and reference-density invariance, layered
hydrostatic rest and pressure, sinking/rising under lateral density contrast,
gravity reversal, exact zero-gravity compatibility, and invalid/cancelled input.

## Dynamic SGS on conservative gas faces

Periodic, independently fitted thermal and planar-channel dynamic SGS now
retain the nonnegative Germano coefficient independently of the cell strain.
Each conservative face evaluates `mu_t = rho_face C_s_face² Delta² |S_face|`
using the same normal velocity difference and transverse traces as its stress
flux. Density and coefficient are averaged separately. The existing cache is
reused, with no additional arrays or change to its storage budget. The face
viscosity feeds the original total-energy flux, fixed-Prandtl heat transport,
ratio admission and parabolic timestep bound. Independently fitted scalar heat
transport keeps its original model. This fixes loss of modeled face shear
when adjacent centered gradients vanish; it does not resolve unrepresented
scales or validate a turbulent flow. Regression sources include a positive
Germano-fit field with zero adjacent centered strains, Cartesian face-flux
references, channel references and conservative evolution.

## Native compressible vorticity budget

`compressible::spatial::vorticity` provides an explicitly byte-budgeted,
reusable `VorticityWorkspace3d` and `Euler3d::vorticity`. Every query decodes
each accepted conservative cell once through its original frozen-gas EOS.
The result separates physical curl, divergence, vortex stretching,
compressible dilation, and baroclinic production
`grad(rho) cross grad(p) / rho^2`. Their sum is the smooth inviscid material
rate D(omega)/Dt, not the Eulerian partial time derivative. Existing
conservative pressure fluxes already generate baroclinicity; querying never
adds another source or changes accepted cells, time, loads, or reports.

Centered derivatives use each physical axis width and periodic wrapping.
A singleton periodic axis denotes homogeneous extrusion; every other axis
requires at least three cells. Nonperiodic boundary cells use a second-order
one-sided reconstruction through interior cell centers, without inventing
wall or reservoir ghosts. Cancellation or numerical/source refusal exposes
no partial output, and repeated queries allocate no new field buffers.

The G0/G1/G4 tests cover affine expansion/rotation and Galilean boosts,
variable-caloric-source baroclinicity in all orientations, barotropic zero,
anisotropic quadratic boundary stencils, periodic refinement, native Euler
baroclinic onset, and workspace/cancellation refusals. Independent numerical
reference calculations were performed; native Rust execution was unavailable
when this extension was introduced. This Estimated output does not certify
shock derivatives, exact discrete circulation, viscous torque, applied-force
curl, or unresolved turbulent production.

### Evolving-composition vorticity

`MixtureEuler3d::vorticity` reuses this same workspace, physical derivative
and inviscid budget on the original accepted multicomponent cells. Each cell
is recovered once through its actual mixture EOS, retaining evolving species
densities, energy references and temperature-dependent heat capacities.
Composition affects pressure-density misalignment through that EOS; no
separate composition torque is added. Storage is independent of species
count, and repeated queries allocate no field buffers. Byte caps, physical
grid identity, unresolved-axis refusals and cancellation remain explicit.

The mixture regressions cover all three baroclinic orientations, rotating and
expanding affine flow, a barotropic mixture with varying composition, and
baroclinic onset through the native conservative time step under refinement.
They also cover transactional cancellation, workspace reuse and admission.
These observations retain the smooth inviscid limitations above; molecular,
applied-force and modeled stress production require their own contributions.

## Dynamic mixture SGS face stress

Under `experimental-sgs`, periodic and planar-channel Favre dynamic mixture
closures retain the fitted nonnegative coefficient separately from cell strain
and evaluate eddy viscosity on the actual conservative face gradient:
`mu_t,f = rho_f (C_s²)_f Delta² |S_f|`. The normal derivative uses the two
adjacent cell velocities; transverse derivatives retain the existing centered
reconstruction. A nonzero fitted coefficient remains active even when a
centered cell strain vanishes. The existing cache is reused without increasing
the prepared per-cell workspace.

The same face viscosity drives momentum stress and fixed-Prandtl/fixed-Schmidt
heat and species transport. Independently fitted thermal and species
diffusivities keep their original owners. Molecular transport, source-owned
thermodynamics, stress work, species enthalpy, conservative advancement,
transport-ratio refusal and timestep bounds remain authoritative. The
face-reference tests include a positive-Germano-fit velocity field whose two
adjacent centered strains vanish while its resolved face shear is nonzero.

## Native planar Smagorinsky stress

With `continuum` and `experimental-sgs`,
`ns::sgs::PlanarSmagorinsky` binds an explicit nonnegative coefficient,
positive cell filter widths and molecular kinematic viscosities to the original
BDM1 mesh. It adds `nu_t = (C_s Delta)^2 sqrt(2 D:D)` to the existing symmetric
stress operator, with strain and convection evaluated at the same nonlinear
endpoint in BDF1 or equal-step BDF2. Original mixed velocity/total-traction
boundaries, backflow policy, harmonic SIP, history, convergence and constitutive
source-refusal checks remain authoritative. A zero coefficient retains the
molecular constitutive law. Prepared mixing lengths are reused and the local
constitutive callback allocates no storage.

Accepted-velocity diagnostics report cell-interior specific SGS transfer in
m²/s³ and its triangle-area integral in m⁴/s³; multiplication by constant density
converts the latter to power per span. These exclude molecular, SIP and upwind
contributions and are neither accumulated heat nor a whole-step energy balance.
This is an Estimated planar screen, without a resolved three-dimensional
turbulent cascade, validated LES, wall model, transition model or backscatter.
The deferred three-dimensional turbulence epic remains open.

`tests/incompressible_sgs.rs` covers analytical shear, filter scaling, rigid
rotation and translation invariance, mixed-traction BDF1/BDF2 endpoint
consistency, zero-closure molecular equivalence, no-slip decay and refusal of
unrepresentable transport without mutating either history.

## Source-owned stress vorticity production

`compressible::spatial::molecular` adds `StressVorticityBudget3d` and
`StressVorticityWorkspace3d` to both `MolecularTransport3d` and its coupled
`MolecularEuler3d` owner. The native `vorticity` query evaluates the original
configured face stress, reuses the advancement kernel's conservative divergence,
divides by each accepted cell's density, and computes
`curl((div(tau))/rho)`. Temperature-dependent shear/bulk viscosity, density
gradients, full cross stresses, selected SGS and authored no-slip/Navier-slip
or enabled reservoir traction remain in that operator. Molecular and SGS
contributions are combined explicitly, not mislabeled as independently
separated measurements. No additional momentum or heat source is applied.

One original EOS recovery per cell supplies both the inviscid budget and
the stress operator. Source primitives, velocity gradients and face buffers
reuse prepared molecular storage. Additional byte-budgeted storage is one
three-component acceleration and one combined result per cell (152 numerical
bytes per cell), with no new field allocation during queries. The same
physical derivative handles acceleration and primitive fields, preserving
anisotropic spacing, periodic wrapping, singleton periodic extrusion and
one-sided interior derivatives. The output is Estimated for smooth resolved
fields, without a wall-sheet, shock or exact-circulation claim.

Queries use accepted held reservoir profiles and do not sample a timed provider.
They preserve accepted cells, physical time and published interval loads;
transport-law and closure caches remain reusable scratch. Cancellation or
source refusal exposes no partial output. The combined material-rate method
excludes imposed body-force curl; Fourier heat affects subsequent states
without a separate instantaneous momentum source.

The G1/G4 tests in `molecular/vorticity_tests.rs` cover Fourier vortex decay
under the actual molecular advance, quadratic convergence, irrotational
vorticity production from variable viscosity/density, no-slip and Navier-slip
wall generation, and capped/cancelled/reused coupled observations. Native
rustfmt and isolated Rust checks of the exact shared divergence/derivative
kernels against independent stress fixtures passed. Full crate compilation
and these native owner-level tests were unavailable when introduced.

## Full-gradient closures in native planar flow

`ns::sgs::PlanarSgs::with_model` admits explicit Smagorinsky, WALE or Vreman
models. WALE and Vreman use the complete original BDM1 cell gradient, embedded
with zero out-of-plane components, through the shared SGS kernels. One
effective viscosity per triangle is prepared at each implicit BDF1/BDF2
candidate; stress, convection and the independently reassembled acceptance
residual share that candidate. Scalar rheology and Smagorinsky retain the
original strain calculation through the same prepared-response solver seam.

Filter widths are explicitly isotropic; molecular viscosities remain frozen
native cell inputs. Zero coefficient retains the molecular constitutive law.
The existing per-mass transfer and area-integral units are unchanged.
This planar embedding supplies no general three-dimensional near-wall
asymptote, turbulent cascade or transition validation. The native SGS
regressions distinguish equal-strain gradients with different rotation, check
the actual planar wall limits, and extend mixed-traction BDF1/BDF2 endpoint
and disabled-closure comparisons across all three models.

### Physical-time Lagrangian dynamic gas SGS

The experimental periodic molecular/Euler composition now accepts a caller-authored
Lagrangian Favre--Germano stress model. The shared signed numerator/denominator
history is sampled from the accepted physical start state, with memory time in
seconds, and its fitted nonnegative coefficient is held across both molecular
halves, SSP/RKL stages and whole-split retries. Conservative faces still use their
actual stage density and face strain, fixed authored turbulent Prandtl heat
transport, and the existing viscosity-ratio, stability and conservation admission.
Only a completely accepted interval commits history and its public fit report;
failed or cancelled intervals preserve the accepted sample, and the next fit uses
the actual elapsed sample time after any shortened interval. Coupled stress-vorticity
observations evaluate an uncommitted endpoint trial and restore the same accepted
history/report on success or failure.

The additional numerical payload is 208 bytes/cell (including both history banks,
physical-field staging and held coefficients), with prepared storage reused during
stepping. This is first-order pathline/coefficient time integration on a fully
periodic Cartesian grid with at least three cells per axis. It does not supply
wall/inflow memory, dynamic scalar history, backscatter energy, isotropic SGS
energy, shock closure, adaptive memory calibration, or benchmark validation.

Validation in this change is limited to rustfmt parsing/formatting and focused
source review; the added Rust regressions have not been compiled or executed in
this environment.

### True steady native SGS flow

With continuum and experimental-sgs, PlanarSgs::solve_steady connects the
existing Smagorinsky, WALE and Vreman closures to the original steady mixed
Navier--Stokes owner. The initial velocity is only a nonlinear guess: the
operator has zero mass coefficient and no time history. One transport response
is prepared per triangle per assembly, and acceptance checks convection and
SGS viscosity at the same endpoint using the original correction and nonlinear
residual gates. Velocity/total-traction/Navier-slip boundaries, backflow policy,
pressure datums, harmonic SIP stress and physical reactions retain their native
semantics. Original scalar rheology uses the same private prepared-response seam.

Three focused regression groups cover full-gradient affine flow with mixed
traction/slip and different guesses, exact-zero molecular equivalence, and
source/anchor/nonlinear refusals. Rustfmt parsing/formatting and bounded source
review completed; these Rust tests were not compiled or executed here. This
remains an exploratory frozen-material planar algebraic screen, without a
three-dimensional turbulent cascade or Reynolds-averaged calibration.

### Material rheology composed with planar SGS

The material-bound native mixed flow exposes explicit SGS variants of steady
rheology, BDF1, equal-step BDF2 and variable-step BDF2. A `MaterialSgsConfig`
provides the existing Smagorinsky, WALE or Vreman model and one physical filter
width per original triangle. At each nonlinear candidate, a single native cell
gradient supplies both its objective source shear rate and the full-gradient
SGS response. The original `RheologyField` resolves and caches the molecular
owner; `nu_effective = mu_source/rho_reference + nu_t` supplies the original
symmetric stress and harmonic SIP assembly. Source rates are not inflated by SGS.

Convection, molecular viscosity and modeled stress converge at the same
accepted endpoint under the existing unrelaxed-update and nonlinear-residual
gates. The steady entry has no inertia. Transient entries retain the original
physical histories and the native BDF1/BDF2 coefficients and step-growth
admission; an artificial history combination never becomes a constitutive
sample. Original physical traction, slip resistance, pressure datum, force
density and reference-density conversions remain intact. Changed-density
responses, source refusal, malformed closures and failed solves publish no
partial result and change neither caller-owned velocity history.

`MaterialSgsStep` retains the exact accepted material samples and SGS diagnostics
using those molecular responses, with no reporting-time source requery. Its
material viscous-heating report contains only `mu_source gamma^2`; the separate
`modeled_transfer_w_per_m` is the area integral of `rho nu_t gamma^2`.
Boundary reactions belong to the full solved molecular-plus-SGS stress. Modeled
transfer is not automatically deposited as heat, and excludes SIP/upwind and
time-discretization losses. The reference density remains fixed; temperature
and composition captured by the provider are prescribed, not evolved.

Manufactured regressions combine nonlinear source viscosity with mixed physical
traction/slip boundaries and all four steady/time entries. Equal-strain,
different-rotation fields distinguish the full-gradient models while retaining
the same material rheology. Other checks cover exact zero-SGS/native equivalence,
physical pressure and separate powers, original source samples, and late
density refusal with unchanged histories. These new consumer tests were authored
and source-reviewed; full runtime compilation remains blocked by unavailable
dependency artifacts and build storage. This remains an exploratory planar
screen, without a resolved three-dimensional cascade or wall/transition model.

### Physical-time Lagrangian SGS in evolving mixtures

The experimental spatial mixture MolecularEuler3d now attaches an explicit
MixtureLagrangianSubgridConfig to a prepared molecular operator without an
existing SGS closure. The shared Favre--Germano history samples the accepted
physical start state once, holds its coefficient through both molecular halves,
SSP/RKL stages and complete retries, and commits only with complete acceptance.
Current face density, strain, composition, Cp/Cv, Maxwell--Stefan/Soret/Dufour
transport, zero-net-mass common-Schmidt species flux and referenced partial
enthalpy transport remain with the original conservative owners. No extra
viscous or chemical heat source is added. The next history fit uses actual
elapsed physical sample time after shortened intervals.

Additional SGS numerical payload is 216 bytes/cell, and detached history
capacity remains counted against the complete coupled budget during stepping.
Prepared storage is reused. This is first-order pathline/coefficient integration
on fully periodic Cartesian grids with at least three cells per axis. It binds
the nonreacting physical owner; raw molecular advancement outside that prepared
transaction refuses. Fixed authored Prandtl/Schmidt numbers provide no independent
scalar-memory, reacting-history, wall/inflow, SGS-energy or shock qualification.

Five focused regression groups cover conservative stress/species/heat/enthalpy,
physical histories and formation-energy gauge, forced retries and cancellation,
multistage RKL replay, and budget/admission/molecular-limit behavior. Rustfmt and
bounded source review completed; Rust tests were not compiled or executed here.

### Independent Lagrangian temperature-flux history

The experimental periodic gas composition also accepts
LagrangianThermalSubgridConfig through
MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs. Momentum and temperature
use independent shared Germano histories and authored memory times, with one
accepted start-state sample and one atomic history commit for the whole physical
split. The original step and vorticity transaction is shared with the fixed-Prandtl
path. Internal SSP/RKL stages, retries, failed steps and observations cannot
independently age either history or publish a partial pair of fit reports.

Faces reconstruct their full current gradient once and evaluate both held
coefficients against current face density and strain. Heat uses the independent
C_phi coefficient and original source Cp(T_face); it remains active when modeled
momentum viscosity is zero. Both actual stress/molecular and heat/molecular
ratios are admitted, with the original conservative total-energy flux and
stability estimate and no second heating source. The paired history/filter and
all physical staging/held fields cost 320 numerical bytes/cell. The fixed-Prandtl
path retains its 208-byte/cell payload, and neither path allocates history during
stepping. Both remain first-order in the pathline/coefficient time update.

Four focused regression groups cover shared-history agreement, independent
nonzero heat with zero modeled stress, actual conservative face fluxes,
observation and late-RKL cancellation/replay, budgets, thermal refusal and the
exact uniform molecular limit. Rustfmt and independent bounded source review
completed; these gas integration tests were not compiled or executed here.
Periodic-domain, smooth-field, exploration-only and unqualified turbulence/model
boundaries remain; no wall/inflow memory or compressible SGS-energy model is added.

### Lagrangian turbulence in the complete reacting mixture step

`ReactingEuler3d::attach_lagrangian_dynamic_sgs` and
`attach_lagrangian_dynamic_species_sgs` attach fresh periodic histories to
the original borrowed molecular operator. The latter fits momentum,
temperature and a common all-species Germano objective independently.
The complete diffusion/chemistry/advection/chemistry/diffusion transaction
samples the accepted physical START state once. Every chemistry solve,
molecular SSP/RKL stage and whole-interval retry uses those held fits;
only complete physical acceptance publishes histories and fit reports.
The next sample uses the actual elapsed accepted duration.

The original chemical mechanism, finite-rate source integration,
Maxwell--Stefan/Soret/Dufour fluxes, evolving heat capacities, mass-frame
closure and referenced species enthalpies remain with their existing
owners. This supplies subgrid transport alongside resolved finite-rate
chemistry; it does not model turbulent reaction rates, unresolved flame
structure, scalar variance, SGS energy or shock interaction.

Attachment validates gas/grid/boundary ownership and complete current
reactor/operator capacities before publishing anything. Detached history
remains included in the byte budget; the added payload remains 216
bytes/cell for fixed transport ratios or 448 bytes/cell for three
independent histories, without species-count-dependent storage.
Cancellation and refusal preserve the prior accepted physical state,
loads, sample time and all fit reports. Public reports live on the
molecular operator that owns the history. Standalone molecular steps
cannot advance it outside a prepared physical interval.

Focused regression sources cover active finite-rate reaction and
conservative balances, accepted-start sampling, shortened/rejected
intervals and cancellation replay, source/budget refusals, fixed-ratio
compatibility and the zero-fit molecular limit. Source review and pinned
Rust formatting/parsing completed; the full consumer test target was
not compiled or executed in the constrained environment.

### Independent Lagrangian heat and species mixing in evolving mixtures

`MolecularEuler3d::with_lagrangian_dynamic_species_sgs` accepts
`MixtureLagrangianSpeciesSubgridConfig` with independent momentum, temperature
and common-species memory controls. The original complete nonreacting split
owns one shared triple-history filter. It samples the accepted physical START
state once, reads every original mass fraction without a transpose, and holds
all three coefficients through both molecular halves, SSP/RKL stages and
whole-split retries. The next sample uses the actual elapsed accepted duration.

Conservative faces apply those held coefficients to current face density and
the full current strain. Heat uses the original source Cp; common species flux
retains the present-species mass closure, exactly absent species, original
Maxwell--Stefan/Soret/Dufour transport and matching referenced partial enthalpies
in total energy. The original combined stability estimate includes each
transport. Stress and heat retain their actual molecular-ratio gates, and
`maximum_species_diffusivity_m2_s` bounds both fitted cell and actual face
diffusivity. Independent scalar transport remains active at zero momentum fit.

Complete acceptance commits all three histories and fit reports together.
Failures, cancellation and rejected retries publish none of them. The
`lagrangian_subgrid_heat_report` and `lagrangian_subgrid_species_report` getters
share the existing accepted-start sample clock; they are not endpoint fits.
Standalone molecular advancement outside the prepared physical step refuses.

The additional numerical payload is 448 bytes/cell: 376 shared history/filter,
40 physical rho/u/T staging, 24 held coefficients and 8 original rho*Cv cache.
Storage is independent of species count and reused during stepping; detached
history remains budgeted. The fixed-Prandtl/Schmidt path retains 216 bytes/cell.
Both paths retain first-order pathline/coefficient integration, explicit
authored timescales and periodic Cartesian exploratory admission. No turbulent
reaction-rate, wall/inflow history, SGS-energy or shock closure is supplied.

Six consumer regression groups cover independent zero-momentum heat/species
fluxes, budgets/refusals, uniform molecular bits, accepted-start/source-gauge
conservation, forced shortening/cancellation/replay and RKL publication. These
consumer tests were parsed and reviewed but not compiled or executed here;
full consumer validation remains outstanding after upstream asupersync compiler
SIGKILLs in earlier full and metadata-only builds. The shared core passed 68
real-source tests. A separate executed core probe confirmed the consumer shear
fixture has exactly zero momentum fit and positive heat/species diffusivities.

### Physical-wall Lagrangian histories for frozen-composition gas

The original gas `MolecularEuler3d` now provides
`with_lagrangian_dynamic_channel_sgs` and
`with_lagrangian_dynamic_channel_thermal_sgs`. They select momentum with
authored turbulent Prandtl heat transport or independent momentum and
temperature memories, respectively. The wall normal is inferred from the
unique reflecting Euler pair; the other axes must be periodic and every
dimension must contain at least three cells. Both original physical walls
are explicitly no-slip, with their admitted tangential velocities and
adiabatic/isothermal thermal laws.

Accepted-start preparation reads the current molecular wall data through
the shared channel helper. Both original molecular halves, SSP/RKL stages
and complete retries retain the held coefficients; actual face density,
strain and source Cp continue to update at each stage. Original wall
heat/work remain intact, including explicitly selected equilibrium stress.
Subsequent wall replacements retain no-slip admission and feed the next
history sample.

The same wall-aware preparation is used for vorticity observations.
An endpoint diagnostic may evaluate a trial but never commits history,
sample time or fit reports; success and refusal restore the accepted
owner. Physical steps keep the original complete acceptance transaction.
Normal history departures obey the shared half-cell wall bound; no
opposite-wall wrapping or inflow history is introduced.

Additional numerical payload remains 208 bytes/cell for momentum with
fixed-ratio heat, or 320 bytes/cell for independent momentum/heat.
These caps retain the gas API's existing meaning: SGS buffers only,
excluding the original Euler/molecular storage. No new history or staging
allocation occurs during stepping. The physical normal axis is exposed
by `lagrangian_channel_normal_axis`.

Three regression groups cover all axes/current walls/original wall flux,
accepted-start reference fits, observation nonpublication, late RKL
cancellation/replay, domain/storage refusals and exact zero-SGS molecular
trajectories. All changed consumer sources passed pinned formatting/parsing
and bounded review, but were not compiled or executed here. The shared core
passed 74 real-source SGS tests. First-order authored-memory and exploratory
model limits remain; no DNS, SGS-energy or shock-interaction claim is added.

### Lagrangian turbulence beside reacting and nonreacting mixture walls

`MolecularEuler3d::with_lagrangian_dynamic_channel_sgs` and
`with_lagrangian_dynamic_channel_species_sgs` attach the shared wall-aware
history producer to the original complete molecular/advection split.
The first selects momentum with authored Prandtl/Schmidt ratios; the second
fits independent momentum, temperature and common all-species histories.
Exactly one Euler axis has reflecting ends, the other two are periodic,
and every axis has at least three fluid cells. The supplied physical
lower/upper walls must be explicitly no-slip. Their tangential motion and
isothermal/adiabatic laws remain with the original molecular owner.

The reacting owner provides corresponding
`attach_lagrangian_dynamic_channel_sgs` and
`attach_lagrangian_dynamic_channel_species_sgs` methods. Configure the
borrowed molecular operator's no-slip walls before attachment. These methods
preserve the original mechanism and complete D/C/A/C/D transaction.
Source-owned finite-rate chemistry and subgrid transport coexist; no
unresolved reaction-rate or flame-structure closure is inferred.

Each accepted-start fit reads current physical wall velocities and heat
conditions. Only tangential directions are test filtered or periodically
wrapped. History interpolation uses the original core's bounded half-cell
extension and refuses a trajectory that crosses an impermeable wall.
Common species use interior normal derivatives. All coefficients remain
held across internal phases and rejected whole intervals; conservative
faces still use current stage strain, density, source calorics, mass-frame
closure and referenced enthalpy. Original molecular wall heat/work and
impermeable species flux remain intact.

The physical normal axis is exposed by `lagrangian_channel_normal_axis`;
configuration, sample-time and fit-report getters are shared with periodic
histories. Numerical storage remains 216 or 448 bytes/cell beyond the
original solver, with no new per-step history fields. The same accepted-step
transaction commits or rolls back every history and report.

Three nonreacting regression groups cover all axes and current wall data,
original wall fluxes, accepted-start evolution/conservation, late
cancellation/replay, budgets and the exact zero-SGS molecular limit.
One reacting-channel regression composes real finite-rate heat release,
moving isothermal walls, three fit reports and transactional replay.
These consumer sources passed pinned formatting/parsing and bounded
production review; they were not compiled or executed here. The wall-aware
core passed all 74 real-source SGS tests. First-order authored-memory,
resolved-wall exploratory and unqualified-turbulence boundaries remain.

### Explicit diffusive closure for evolving-mixture pressure outlets

The spatial mixture `MolecularTransport3d` accepts an explicitly authored
`PressureOutletPolicy::ZeroDiffusiveFlux` through
`new_with_pressure_outlet_policy`. The existing `new` constructor and the
policy default retain `Refuse`: pressure alone does not define a molecular
bath. The selected closure sets the complete outlet face vector to zero,
including every species diffusion flux, viscous momentum flux, stress work,
Fourier/Dufour heat and referenced species-enthalpy transport. Interior original
Maxwell--Stefan, thermal and full-gradient viscous transport remain active.

This is an explicit zero-normal-diffusive-transfer truncation assumption.
The existing interior-velocity extrapolation supplies nearby gradient samples;
no incoming composition or temperature is inferred. It is not a general
nonreflecting or backflow-stable outflow treatment. Original Euler pressure,
choking, backflow refusal, acoustic bounds and conservative advective species,
both shear momenta and energy remain with the source-owned pressure outlet.
Unbound transmissive ends still refuse, complete reservoirs retain their bath
exchange, and periodic and reflecting boundaries retain their native laws.

The stored policy follows standalone SSP, stabilized RKL, adaptive RKL and
coupled molecular/advection stages, retries and boundary accounting. It adds no
numerical field storage. The coupled owner admits the prepared operator's
policy and retains its transactional publication and cancellation behavior.
Existing SGS factories continue to use their original strict constructor;
this extension does not add an outlet SGS calibration or boundary model.

Five focused regression groups cover strict/default admission and unchanged
wall/periodic trajectories, independent homogeneous-Neumann species/shear decay
and its vorticity, uniform stabilized/adaptive stages, all-six-face original
advective fluxes with a complete opposite bath, and nonuniform draining species
with conservative energy accounting and final cancellation/replay. Rustfmt
parsing/formatting and bounded source review completed; these Rust integration
tests were not compiled or executed in this constrained environment.

### Locally nonisothermal binary-liquid diffusion

`miscibility::nonisothermal` binds an accepted binary field to per-cell temperatures and the same original homogeneous caloric model. It retains both compensated component inventories and elapsed physical time. Equal fixed partial molar volumes, kappa=0, strictly positive common species support, and a convex original law over the complete authored temperature bracket are required. This lane does not localize square-gradient interface energy or perform an equilibrium flash.

Each canonical face uses two implicit scalar exchanges: an entropy-conjugate chemical/enthalpy exchange with frozen original kinetic mobility and transported partial enthalpy plus authored Soret heat, then Fourier heat exchange. Both recover local temperatures through original pure and mixing calorics. Conductivity is the zero-diffusive-flux value, and the same signed heat-of-transport coefficient supplies reciprocal Dufour heat. Original D(T,x) claims are sampled at actual endpoint temperatures/compositions and the mean face outside scalar roots. Completed stages and final fields retain original source support. Exact physical zero coefficients are admitted; positive unrepresentable mobilities refuse.

Every proposal checks whole-field component amounts, per-cell capacity, total enthalpy and actual entropy before publishing. The first-order split has adaptive step doubling with independent composition and temperature tolerances, accepts the two-half-step result, and applies the original conservation/entropy allowance again to the whole requested interval. Cancellation, exhausted work, source refusal or failed numerical admission leaves the input unchanged. Conductivity and heat of transport are explicitly authored fixed face data; this owner does not infer experimental coefficients or model pressure/volume work or mechanics. Prescribed advection is provided by the extension below.

### Enthalpy-carrying spatial miscibility transport

`miscibility::nonisothermal::advection` uses the SAME prepared directional
volume flow, original geometry and donor-volume CFL admission as isothermal
miscibility. Every stream carries both independently represented species and
the original donor parcel's total homogeneous enthalpy, including nonideal
mixing. Every debit has the same opposite credit. All faces sample the preceding
substep synchronously, and bounded inverses of the original caloric model
recover local temperatures. Temperature is not volume-averaged and no second
mixing-heat source is added. Opposing streams remain distinct even when their
net volume rate vanishes. Exact zero flow performs no caloric inverse.

`miscibility::nonisothermal::transport` composes this advection with the original
chemical/Fourier/Soret-Dufour sweep over ONE physical interval. It retains both
original stage reports and total per-face species/enthalpy transfers. Every
complete candidate independently checks species, cell capacity, total enthalpy
and actual entropy using the original whole-interval allowances. Caloric
inverse tolerances share the energy budget across cells/substeps. The original
material diffusivity is qualified at actual local temperatures/compositions
after every advection substep and throughout the original diffusion stages.
Failures, late cancellation and unsupported caloric states publish no field.

The scheme is first order in time and donor-cell in space. CFL bounds establish
positivity, not temporal/spatial accuracy. The existing kappa=0, convex-law,
positive-common-support and equal-fixed-partial-volume scope remains explicit;
there is no momentum feedback, inferred stirring power or viscous heating.
Focused regressions compare independent unequal-Cp parcel balances, both signs
of nonideal mixing heat, directed circulation, unequal-volume analytic donor
decay, complete species/energy receipts, zero-operator limits, macro-time
refinement, finite work and cancellation/replay.


### Adaptive stirred nonisothermal mixtures

`miscibility::nonisothermal::transport::adaptive` controls the COMPLETE
advection/chemical/Fourier operator with independent composition and Kelvin
tolerances. Each coarse or half trial has exactly one donor step, and its macro
interval is capped by the original volume-flow CFL bound. This prevents
identical internal subcycles from falsely reporting zero temporal error.
The two-half-step result is accepted; rejected/coarse trials contribute solver
work but never physical species/enthalpy receipts. Both independently retained
species, including compensated low words, participate in the shared estimator.
A zero trace tolerance cannot hide a nonzero amount discrepancy.

The owner reuses original calorics, sourced molecular diffusivity, graph and
prepared flow. It applies the original whole-window component/capacity,
enthalpy and entropy allowances after the complete requested interval.
Temperature-bracket and numerical solve failures may refine within explicit
attempt/minimum-step bounds; original source failures and cancellation remain
refusals. Only accepted face receipts accumulate, the original clock reaches
the requested endpoint, and a short final remainder retains all admission
checks. This is first-order local time-error control, not a global or spatial
accuracy certificate. Six focused regressions cover independent thermal and
composition modes, CFL aliasing, accepted-path balances, work/minimum-step
refusal, continuation and cancellation/replay.

The published numerical sources through commit 1fabe824 passed all 76
miscibility unit tests on the pinned `nightly-2026-07-06` Rust toolchain. This
includes 16 new regressions across enthalpy advection, coupled transport and
adaptive stepping. Verification used a focused Cargo workspace compiling the
original miscibility module and its original internal dependencies, with no
replacement physics or test stubs. It is a numerical-owner check, not a full
workspace or DSR validation, and does not qualify the material adapter below.


### Original materials for locally nonisothermal miscibility

`miscibility::material::nonisothermal` binds an existing source-resolved spatial
field to its exact original `BinaryLiquidCaloricLaw` and one temperature per
cell. Pointer identity of the original liquid law is checked before numerical
binding. The field retains the original liquid/caloric bindings, independently
represented component inventories, elapsed time, resolved-cell budget and the
source's existing equal-fixed-volume allowance.

All five evolution paths are exposed: one molecular/heat sweep, adaptive
diffusion, parcel advection, coupled transport and adaptive coupled transport.
Each completes its numerical transaction and then refreshes density, viscosity,
transport Cp, conductivity and original material receipts at the actual final
temperature/composition of EVERY cell. The existing homogeneous source resolver
and volume-consistency check are reused. Unsupported source temperatures,
incompatible density/volume, late source errors and cancellation publish neither
the numerical field nor a partial material report. Source-domain checks for
calorics and any bound molecular diffusivity remain within the numerical stages.

This additional material qualification applies to the complete endpoint;
independently sourced transport properties are not certified for every discarded
or intermediate numerical trial. Mixture transport Cp does not replace original
pure-liquid calorics, and local conductivity receipts do not silently supply a
face interpolation law. Conductivity and reciprocal Soret heat remain explicitly
authored face data. No cell is flashed to equilibrium. Actual cell temperatures
come from `field().temperatures()` or `cells()`; the network/source reference
temperature is not a local observation.

Three focused material regression groups cover temperature-dependent source
properties, retained card/pressure/caloric identity, independent parcel-enthalpy
temperatures, all five evolution paths, exact-law mismatch, late source and
density refusal, cancellation and replay. These regressions were source-reviewed
but remain uncompiled and unexecuted. The final material check used the pinned
Rust toolchain, a single Cargo job, offline dependencies and the original-source
material dependency closure. It stopped when the compiler emitting metadata for
the exact pinned `asupersync` dependency received `SIGKILL`, before reaching
`fs-material` or the changed integration tests. No diagnosis against this adapter
was produced; the interrupted dependency build is not evidence that the adapter
compiles. DSR and RCH were unavailable, and the full workspace was not validated.

### Locally nonisothermal phase separation with interfacial energy

`miscibility::nonisothermal::phase_field` admits the original nonconvex binary
liquid law with positive, constant energetic square-gradient coefficient.
It retains the same rigid volumes, equal fixed partial molar volumes, paired
inventories, original caloric sources and optional original molecular
diffusivity. Every cell has its own supported temperature; both species remain
strictly positive. The existing convex, zero-gradient owner keeps its admission.

Each canonical edge stores `kappa A/(2 d) (x_left-x_right)^2`. Half belongs to
each endpoint's local energy. Entropy is the original homogeneous entropy:
constant energetic kappa contributes no entropy. The molecular entropy force
therefore uses the FULL incident stencil weighted by
`kappa A/(2 d) (1/T_left+1/T_right)`, including temperatures of neighbors that
do not exchange species on the current face. This is an explicit discrete
energy localization, not an inferred temperature for an unresolved interface.

A canonical face freezes its original transport coefficients and temperatures
during a bounded scalar convex-split solve. Original ideal logarithms, negative
chi and the full weighted graph gradient are implicit; positive chi is explicit.
The deterministic logarithm matches the original thermodynamic kernel. Paired
energy transport carries the original divided partial enthalpy and authored
reciprocal heat of transport once. Every cell incident to a changed interface
then recovers its actual temperature from its original caloric inverse and its
own half-edge energy increment, including neighbors with unchanged inventories.
The same implicit Fourier kernel as the convex local-temperature owner supplies
heat conduction. Sparse incident patches reuse scratch across faces; no global
matrix or full-mesh scan is performed inside a composition root.

Before publication, the complete interval checks both component inventories,
cell capacities, homogeneous enthalpy PLUS interfacial energy, and original bulk
entropy. This frozen-temperature split is first order and is not unconditionally
entropy stable. A temperature-domain, energy or entropy refusal requires a
smaller interval; an isothermal free-energy decrease never overrides the actual
thermal account. The owner has no flow-advection operation because interfacial
transport and mechanical work require a consistent momentum coupling.

Seven focused regressions cover spinodal growth, a neighbor's thermal response
without species transfer, an independent entropy-gradient finite difference at
fixed local energies, Fourier convergence with excess heat capacity, unequal
volumes and temperature-dependent original calorics, zero-gradient refusal and
late cancellation/replay. These validate the stated discrete model; they do not
supply nucleation kinetics, variable-kappa calorics, unequal partial molar
volumes, mechanical backreaction or an experimental accuracy claim.

The complete focused numerical suite passed 83 unit tests on the pinned Rust
compiler using original sources. This excludes the material adapter and full
workspace/DSR validation.

### Adaptive local thermal phase separation

`phase_field::adaptive` resolves the complete composition/interface/Fourier
sweep with whole-step versus two-half-step control. It reuses the original
compensated composition and Kelvin estimator and adds an independent
interfacial-energy allowance. That norm sums absolute differences of INDIVIDUAL
edge energies, with one whole-graph absolute allowance and a relative scale
formed from the larger energy on each edge. Equal total interface energy cannot
hide errors that moved between edges.

Only the accepted two-half-step path contributes species and total local-energy
transfers. All completed coarse and rejected sweeps contribute reported solver
work. The original interval, local states and material bindings stay unchanged
on exhaustion, source refusal, unresolvable minimum step or cancellation. The
clock follows actual trial endpoints, including continued runs and final
remainders shorter than the nominal minimum step. Complete-window species,
capacity, homogeneous-plus-interface energy and entropy budgets apply again
before a candidate is returned.

Six additional regressions cover interface-driven refinement, spatial error
cancellation, the independent Fourier relaxation mode, accepted-only local
energy/species receipts, work/minimum-step/cancellation refusal and exact
continuation. All 89 numerical miscibility unit tests passed in the original-source
focused workspace. This first-order controller estimates local temporal error;
spatial consistency still depends on the original authored graph geometry and
resolution, and no global trajectory error bound is supplied.

### Original material bindings for local-temperature phase separation

`miscibility::material::phase_field` binds a previously source-resolved
spatial field to its exact original liquid and caloric objects, retaining the
original square-gradient model and density/volume allowance. Its one-step and
adaptive candidates preserve the numerical phase-field owner, interface-energy
account, paired inventories, elapsed time and original resolved-cell budget.
Actual endpoint temperatures and compositions refresh density, viscosity,
transport Cp, conductivity and material receipts through the SAME resolver
used by the convex nonisothermal material adapter.

The source/caloric identity check precedes initialization. Every returned
endpoint must satisfy the original material domain and fixed-volume allowance;
a late source or cancellation refusal publishes no field or partial material
report. This is endpoint property qualification. Independent material property
receipts do not qualify every discarded/internal trial, supply an automatic
face-conductivity law or replace original homogeneous calorics. Diffuse-interface
cells retain their original homogeneous branch and are never individually
flashed to bulk equilibrium.

Focused material regressions exercise original source identity, actual local
temperature/property refresh, preserved energy accounting, adaptive evolution,
late source refusal and cancellation/replay. They are source-reviewed and
formatted, but not compiler- or execution-validated in this environment.
A fresh single-job check of the original material dependency closure with one
compiler frontend thread and incremental compilation disabled still stopped
when pinned `asupersync` metadata compilation received `SIGKILL`, before
reaching `fs-material` or these integration tests. The 89 passing numerical
miscibility tests do not replace this missing material validation.

### Compatible molecular stress with variable viscosity

Frozen-composition planar and spatial molecular owners sample the original
transport law at cells and use a shared compatible Newtonian mechanical
stencil on periodic, zero-diffusive-flux and reflecting domains, including
natural free-slip, prescribed no-slip and finite Navier slip. Normal face
derivatives use arithmetic cell viscosities; each centered cross derivative
carries its own cell coefficient before face averaging. For nonnegative shear
and bulk viscosity, prescribed wall power minus the semidiscrete molecular
kinetic-energy rate equals the nonnegative bulk strain/residual quadratic plus
physical interfacial slip dissipation. Smooth periodic variable-material
manufactured fields retain second-order force convergence. This prevents a
demonstrated positive-viscosity stencil from creating unforced kinetic energy
and spurious vorticity while conserving total energy and actual boundary work.

`uses_compatible_molecular_stress` reports the complete-domain selection.
No-slip and Navier walls may have constant tangential velocity and their
original adiabatic or isothermal thermal condition; zero-length Navier slip
is the same Dirichlet case. Enabled molecular bath exchange in 3D or an
explicitly configured equilibrium SGS wall law retains the entire original
mechanical discretization. The dissipation identity is not claimed for that
fallback, arbitrary time steps, shocks or the combined SGS model.
Ordinary interior SGS face viscosity, ratio/source admission and heat remain
authored contributions. Modeled viscosity is retained separately so replacing
a large molecular midpoint sample cannot erase finite modeled stress.

Every SSP/RKL stage and stress-vorticity observation refreshes a preallocated
16-byte/cell source-viscosity cache. Mechanical parabolic bounds use the actual
cell coefficients and the global minimum density; heat keeps the native face
coefficient. Matching constant cell/face coefficients at interior faces preserve
original stress, work and step-bound arithmetic. Otherwise traction and total-energy work are
formed directly, without cancellation against a possibly enormous unused
midpoint stress. No extra viscous-heating source is introduced; conductive
heat remains in the same conservative energy flux. Molecular cancellation and
publication remain with the original owners.

Six tests of the exact shared Rust kernel passed in an isolated rustc build:
275 mixed supported-boundary energy identities, the positive-energy
counterexample, constant-coefficient agreement, affine stress, second-order
variable-material refinement, and an extreme midpoint-coefficient cancellation
case. Four owner regression groups cover both dimensions, actual source
sampling and parabolic bounds, retained heat/SGS at extreme contrast, explicit
boundary scope/constant-law trajectory bits, and late cancellation/replay.
Those owner tests were parsed and source-reviewed but not compiled or executed
in the resource-constrained environment; the isolated kernel tests do not
replace a full Cargo or DSR gate.

### Compatible molecular stress for evolving mixtures

The planar and spatial evolving-mixture molecular owners use the shared compatible Newtonian stress discretization on supported impermeable or zero-exchange domains. Planar domains admit periodic boundaries and all valid reflecting wall laws: natural free-slip, prescribed no-slip and finite Navier slip. Spatial domains additionally admit pressure outlets explicitly configured with `PressureOutletPolicy::ZeroDiffusiveFlux`; admission still rejects a pressure outlet without that policy. Baths and explicitly configured equilibrium SGS wall laws retain the original mechanical discretization for the whole domain. `uses_compatible_molecular_stress()` reports this selection.

At every SSP, retry, RKL, and adaptive spatial evaluation, the original transport source samples molecular shear and bulk viscosities at each current cell's temperature, pressure, and composition. Normal face derivatives use arithmetic cell viscosities; cross derivatives are weighted by their originating cell viscosity before averaging. The shared local slip solve supplies the same complete gradient to interior stress, wall traction and ordinary interior SGS. Prescribed-wall power minus the semidiscrete molecular kinetic-energy rate equals nonnegative bulk viscous dissipation plus physical Navier interface friction. This does not assert monotone enstrophy, arbitrary-time-step kinetic monotonicity, or a dissipation guarantee for the independently authored SGS model or fallback boundary domains.

`PlanarTransport::coefficients_from_primitive` and `SpatialTransport::coefficients_from_primitive` expose original-source cell coefficient admission without a second material law or a face diffusion solve. The source validates the actual primitive, coefficients, longitudinal viscosity, and binary diffusivity matrix using prepared scratch. `TransportFaceReport` additionally exposes `bulk_viscosity_pa_s`, `nonmechanical_energy_flux_w_m2`, and `nonmechanical_diffusivity_m2_s`. The latter two receipts retain original species-enthalpy, Fourier, Soret, and Dufour contributions before mechanical work and separate their admission bound from mechanical diffusivity.

Compatible face traction and its conservative mechanical work are formed directly. This avoids cancellation when an unused nonlinear midpoint viscosity is much larger than either cell viscosity. Species fluxes and the source-owned nonmechanical energy flux retain their original values. Wall mechanics also uses the shared physical slip dissipation and authored fluid/wall heat partition; it replaces the complete native mechanical contribution without retaining an old friction-heat remainder. No additional viscous-heating source is introduced. Interior faces whose two cell coefficient pairs exactly equal the original face coefficient pair retain original flux arithmetic. Compatible mechanical rate admission always uses the cell stress bound and the minimum current stage density, including matching-coefficient faces and mechanically active spatial singleton/all-wall domains. The original thermal-wall source continues to supply Fourier heat and its diffusion bound. Ordinary interior SGS remains separate and retains its source model and ratio semantics.

Each owner prepares one shear/bulk coefficient pair per cell, adding 16 bytes per cell to the required and actual base workspace accounting. Stabilized and adaptive workspace requirements inherit the increased base size. Source sampling and all face work remain inside the existing staged commit/cancellation boundary.

Regression coverage adds five planar and three spatial owner test groups for the energy-producing variable-composition vortex, SSP/RKL/adaptive stepping, original species/heat/cross-transport preservation, an extreme nonlinear midpoint-viscosity contrast, actual-cell rate admission, workspace budgets, whole-domain boundary selection, and late cancellation/replay. These owner tests were parsed and reviewed but could not be compiled or executed in the resource-constrained session. The six exact-production shared-stress kernel tests from the preceding frozen-gas correction passed; they establish the shared algebraic identities and convergence evidence, not a full mixture-owner build gate.

### Prescribed-wall vorticity and compatible mechanical power

The prescribed no-slip branch of the shared `stress::wall_stress` kernel uses
`stress::no_slip_traction`. It retains the full cell-centered cross stress
paired with the same odd velocity ghost used by the interior gradient.
Replacing that cross stress with derivatives of the exact constant wall trace
would break the discrete mechanical energy identity. Each face uses its own
prescribed tangential velocity for mechanical work, including corners and
singleton axes; no velocities from different planes are averaged. For no
slip, total-energy flux is the same traction's wall work plus the original
source's Fourier heat. Source coefficients at thermal wall states continue to
own conduction, while cell viscosities own compatible stress and its rate bound.

The existing sixteen-byte/cell viscosity cache, SSP/RKL stage machinery,
stress-vorticity observation and transactional loads remain shared. No new
numerical field allocation is introduced. Nine tests of the exact Rust stress
kernel passed, including six retained tests and three new groups: 850 mixed
stationary/moving wall configurations, explicit wall-driven kinetic gain with
the matching power balance, and variable-viscosity Couette refinement. The
Couette solution and steady wall load approached second order, while traction
evaluated on the analytic sampled profile approached first order at the wall;
there is no blanket second-order local wall-stress claim. Full flow-owner
regressions were source-reviewed and parser-checked but not compiled or run in
the resource-constrained environment.

The same prescribed-wall kernel is consumed by both evolving-mixture owners.
Momentum transfer, boundary loads and conservative energy use the same actual
traction; original species, enthalpy, Soret/Dufour and wall heat remain with the
transport source. No extra buffers or stage arguments are introduced. Spatial
inertial-wall midpoint probes retain their existing fresh stage evaluation and
actual load transaction. Focused owner regressions cover stationary and moving
wall work through SSP/RKL/adaptive steps, composition-dependent shear and bulk
viscosity under axis rotations, thermal transport with a very large unused wall
viscosity, zero/positive slip selection, and explicit equilibrium-SGS fallback.
These owner regressions passed source parsing/formatting and bounded review,
but were not compiled or executed in this environment. The shared kernel's
nine native tests establish the numerical wall identity and the separately
qualified Couette refinement evidence.

### Finite Navier slip and natural free-slip mechanics

Compatible frozen-gas and evolving-mixture wall mechanics solves the fluid trace from
the complete Navier traction, including the tangential derivative of normal
velocity. With half-cell distance `d`, slip length `b`, outward sign `s`
and centered cross derivative `C_t u_n`, the trace satisfies
`u_f = (b u + d U_wall - b d s C_t u_n)/(b+d)`. The implementation uses
relative velocities and stable geometric weights instead of this literal
potentially overflowing expression. Natural free-slip imposes zero tangential
traction with `u_f = u - d s C_t u_n`; its former even tangential ghost is
therefore replaced in compatible domains. This makes large finite slip
approach the same complete free-slip operator, including its neighboring
stress and vorticity production.

At corners and singleton axes, each off-diagonal gradient pair is a local
two-by-two solve. Both `b/(b+d)` and `d/(b+d)` are retained independently;
the determinant is formed from positive complements, without subtracting
nearly equal products. Unrepresentable finite geometric weights and nonzero
physical Navier traction refuse instead of silently erasing friction or
changing the boundary to free-slip. A purely free-slip
singular pair selects zero rotation only when both provisional derivatives
are exactly zero. Zero-length Navier slip follows the prescribed no-slip
branch exactly. Zero viscosity needs no division by viscosity. These changes
reuse the existing stage gradient and viscosity caches, with no new field
buffers, global iteration, or caller/workspace arguments.

For each wall, physical interface friction is `D = -t dot (u_f-U_wall) >= 0`,
where `t` is the outward-normal traction exerted on the fluid. If `chi`
is the authored fluid heat share, the inward total-energy flux is
`t dot U_wall - (1-chi)D + q_heat`. Thus wall power minus kinetic-energy
rate equals bulk viscous dissipation plus interface friction, while internal
fluid energy receives bulk dissipation plus `chi D`. Fourier heat retains
its original thermal-wall source sampling. The prescribed-wall mechanical
rate bound remains conservative under the local dissipative slip elimination.
Enabled baths and explicit equilibrium SGS wall laws keep the complete native
fallback; the shared identity is not asserted for that fallback or for SGS.

This extension's Rust regressions were authored and source-reviewed after the
execution environment became unavailable; they have not been compiled or run,
and the earlier native test results above apply to the preceding no-slip
implementation. An independent scalar-resistance JavaScript reference for
`mu(x)=1+x`, `b=0.3` Couette flow gives 32-to-64-cell error ratios 2.0130
for wall load and 2.0003 for the solved velocity profile. The finite-slip
interface uses the adjacent cell viscosity, so this case is first order;
the earlier no-slip convergence result is not extended to finite slip.

## Finite aqueous exchange on moving material surfaces

With `continuum` and `surface-soluble`,
`CapillarySurfaceDynamics::advance_soluble_surfactant` composes the original
finite-bath adsorption/desorption and lateral surface diffusion with the native
mechanical trial. Chemical half steps use the actual initial and solved face
areas. Mechanics carries the same intermediate face amounts, so dilution,
packing and capillary force follow the original surface free energy. One
original aqueous candidate supplies both halves and publishes with the complete
trajectory and face amounts after all checks.

The whole endpoint independently balances bulk plus surface amount and original
bulk/surface free-energy change. Geometric capillary work, chemical change,
adsorption and lateral-diffusion losses, surface viscosity and external work
stay separate. Substep, scalar-root and pair-application budgets span both
chemical halves. A stable joint area/amount energy increment uses only admitted
physical endpoints and exact existing secants; expansion/adsorption need not
fit the final amount on the old area. Equal-area calls retain the original
stationary fast path. Late work, source, balance or cancellation refusal leaves
all three original states unchanged.

The complete composition remains first order because adsorption uses backward
Euler. Temperature, well-mixed solvent volume, topology, coefficients and nodal
reference masses are held. This is an inertia-free dilute coating exchanging
with an independently held-volume bath; enclosed volume feedback, solvent
displacement, osmotic traction, solute momentum and adsorption heat need their
own physical owners. Isothermal free-energy loss is not automatic heat deposition.

## Adaptive finite-bath soluble surface transport

`surfactant::soluble::adaptive::advance_soluble_surface_adaptive` now completes
a held finite-bath duration with local error control on every surface amount,
the selected bulk-species amount and its original nonideal dimensionless
activity. Each comparison uses one native BE adsorption step against two half
steps. The first-order fine-error factor is one for equal halves; the existing
surface-only diffusion controller keeps its second-order factor of one third
through a shared private schedule. Accepted states are the unextrapolated
physical fine endpoints. The original source datum and co-solute interactions
remain intact, with a log-space activity comparison that avoids overflow under
admitted reference changes.

All coarse, fine and rejected probes spend the same global substep, root and
pair-exchange allowances. Only accepted fine-path gross adsorption/desorption
and chemical/lateral losses enter the physical report. Before publication,
original-to-final combined amount and free energy are independently checked
against the fixed whole-duration allowances. Cancellation or a late failure
preserves both the original aqueous inventory and every surface amount.

Areas, temperature, solvent volume, topology and coefficients remain fixed.
This controller estimates local temporal error; it does not certify global or
spatial accuracy or add a moving mechanical macro controller. Six new
regressions cover exact finite depletion/desorption, nonideal noncommuting graph
dynamics against an independent ODE, first-order rejection and fine-state
recovery, equilibrium gross turnover, whole-call work and cancellation rollback,
and activity-only control through an extreme admitted reference change.

A focused Rust build using the actual surface, aqueous, adsorption and mechanical
source closure passed all 74 tests, including the prior 68 and these six new
regressions. The existing stationary adaptive diffusion tests also passed after
schedule extraction. This is focused native validation, not a full Cargo
workspace build.

## Adaptive coupled motion and finite-bath surface chemistry

`CapillarySurfaceDynamics::advance_soluble_surfactant_adaptive_to` now completes
an exact absolute-time interval through the original coupled chemical/mechanical
endpoint. One complete proposal is compared with two original half proposals
using the first-order estimator required by BE adsorption. Position, velocity,
every face amount, finite bulk amount and original nonideal species activity
receive explicit error scales. Each chemical half executes one native substep,
so hidden subdivisions cannot mask the coupled comparison. The shared schedule
uses actual representable half intervals on the accepted physical clock.

All coarse and rejected work spends the same mechanical, chemical-substep,
root and pair allowances. Trial checkpoints copy endpoint arrays, not accumulated
trajectory history. Only accepted fine transfers, physical losses and telemetry
are staged. Before publication, original-to-final bulk/surface amount, actual
area/amount free energy, kinetic change and held-force work independently close
against the accepted geometric, chemical and viscous accounts. Original history
capacity is reserved before the existing publisher appends accepted telemetry
once. A late failure or cancellation preserves every original physical state.

The native fixed-topology, reference-mass, held-temperature, held-volume bath,
coefficient and external-force scope remains. Only temporal-error rejection
retries; physical/source/solver refusals propagate. Local error estimates do
not certify global or spatial accuracy.

A focused native Rust build of the actual surface, aqueous, adsorption and
mechanical source closure passed 80 tests with zero failures: the prior 74 plus
six new moving-adaptive regressions. These compare a continuous breathing/site
ODE, exact two-half native recovery, first-order rejection, global budget and
late cancellation rollback, motion-only refinement and unequal representable
halves on a large clock. Existing stationary first/second-order adaptive tests
also passed after the shared schedule extension. Full Cargo workspace
validation is not claimed.

### Conservative momentum with spatial body acceleration

`assemble_conservative_momentum_bdf1_with_acceleration` accepts one held physical
acceleration field in m/s² on the original cells. The existing mass transport
first computes the endpoint P0 density; its original six positive triangle
mass-quadrature weights then supply both inertia and the acceleration load.
The field is sampled once per point, in the existing cell traversal, without
a source array. Affine acceleration moments are exact up to roundoff. Additional
volumetric force retains its independent N/m³ meaning, and the uniform/zero
gravity paths retain their original arithmetic.

A nonfinite sample identifies its original cell and point. Unrepresentable
loads, transport or solver refusals, and cancellation publish no partial
density/velocity endpoint. This is the original fixed-mesh, first-order IMEX
scheme with transported P0 density and unchanged quadrature for general fields.
The callback is a pure frozen endpoint source; prescribed centrifugal
acceleration is admitted, but velocity-dependent Coriolis integration is outside
this API.

Four pinned-compiler probes of the actual mass loop, original BDM basis and
quadrature passed: affine loads against independent degree-two moments,
centrifugal source/pressure/traction balance, constant/zero exact compatibility,
and located refusal/cancellation/retry. Their maximum load and weak-balance
errors were 1.155e-14 and 1.919e-13. The three public coupled-solver regressions in
`tests/conservative_density_momentum/gravity.rs` also exercise transported
density, reference normalization and exact angular velocity. Those full Cargo
regressions were added but not executed because the uncached dependency build
exceeded the available storage; the probes do not replace transport or the
Stokes solver.

### Prescribed integral flow with material SGS

The `MaterialFluxSystem` entries
`solve_navier_stokes_mixed_rheology_sgs_flow_control`,
`bdf1_navier_stokes_mixed_rheology_sgs_flow_control` and
`bdf2_navier_stokes_mixed_rheology_sgs_variable_step_flow_control` resolve
the pressure required by authored integral port flows while molecular rheology,
Smagorinsky/WALE/Vreman stress and convection converge together. The existing
native constrained Picard owner prepares one effective constitutive field per
candidate and uses it for all pressure responses. Pressure offsets relax with
the complete field. Final acceptance checks actual port flows and both
independently reassembled, pressure-inclusive nonlinear equation residuals.

Targets are signed outward volume rates per out-of-plane depth [m^2/s];
physical depth converts them to m^3/s. Each disjoint port uses original
traction edges and leaves its resolved velocity profile unknown. One uniform
normal pressure increment is added to its original stress shape, retaining
tangential traction and the selected open-boundary policy. Every controlled
connected component retains an uncontrolled traction datum. Required pressure
offsets, source force and boundary loads use the original material density and
coherent SI conversions.

The steady entry has no artificial inertia or timestep. BDF1 and equal or
unequal-interval BDF2 retain the original physical histories, interval admission
and mass operator; constitutive stress and convection are evaluated at the
returned endpoint. Accepted material samples, physical reactions, molecular
heating and separate modeled SGS transfer retain their original owners without
a reporting-time source query. Refusal leaves both histories unchanged and
returns no partial pressure/flow result. These entries retain the existing
Estimated planar SGS applicability limits.

The three groups in `tests/incompressible_sgs/flow_control.rs` exercise all
three models against an affine physical-pressure solution for steady, BDF1,
equal and unequal BDF2 with natural and dissipative backflow; they check integral
units on a nonunit-height port, pressure-inclusive reactions and separate
heating. Additional cases check increased required channel pressure at fixed
flow, physical density scaling, the exact zero-model limit against the original
material BDF2 control owner, and temporal, source, datum and nonlinear refusals.
Source parsing and bounded review passed; native compilation and execution
remain pending in the constrained environment.

### Opt-in local SGS exchange with evolving-mixture reservoirs

With `experimental-sgs`,
`spatial::mixture::molecular::MolecularTransport3d::with_local_sgs_reservoir_exchange`
explicitly extends an already authored local Smagorinsky, WALE or Vreman
closure to complete source-owned reservoirs. The default keeps molecular
bath exchange. Dynamic and Lagrangian variants refuse this selection because
an exterior test-filter field or pathline history has not been authored.
Pressure-only outlets retain their declared zero diffusive flux.

The modeled face uses the native half-cell normal derivative, midpoint
interior/profile tangential derivatives, and directional Vreman cell widths.
Midpoint density supplies `mu_t = rho_bar nu_t`; midpoint mass fractions and
original species Cp at midpoint temperature supply `k_t = mu_t Cp/Pr_t`.
Common-Schmidt species flux `j_i = -(mu_t/Sc_t) d_n Y_i` closes on a present
species, preserving zero summed modeled mass flux and exact absent species.
The original molecular transport law is evaluated once per physical face;
its constitutive samples and Maxwell--Stefan/Soret/Dufour terms remain intact.

The shared conservative interior/bath kernel retains source-referenced
species enthalpies. At the bath, their temperature and mechanical work
velocity are the prescribed boundary values: energy exchange is conductive
heat plus `sum_i h_i(T_b) j_i - u_b dot tau_n`. No separate heating source
or reservoir inventory is introduced. Molecular and modeled diffusion scales
are added, including stress, heat and species/enthalpy coupling with both
trace densities and heat capacities, before the native three-axis half-cell
stage bound.

The same native spatial callback serves SSP, RKL, coupled advection and
stage-selected reservoir traces. Existing acceptance, boundary-load
publication, cancellation and retry transactions remain the owners of the
complete update. No additional numerical heap workspace is required; an
exactly zero model returns the original molecular face and receipt.

Five focused regression groups are authored for all six boundary
orientations and three local models, nonlinear calorics and prescribed bath
work, original source-call counts, zero-model native equivalence and
nonlocal-filter refusal, low-Schmidt SSP/RKL admission with atomic refusal,
and coupled species/energy conservation, formation-offset invariance and
late-cancellation replay. Pinned rustfmt parsing and bounded source review
were completed; runtime execution of these new regressions remains
unverified in the disk-constrained implementation environment. This closure
retains the existing exploratory SGS and engineering-bound qualifications;
it adds no experimental-validation, positivity or entropy theorem.


## Accepted-state stress observations for acoustic sources

The frozen and evolving-composition spatial `MolecularEuler3d` owners expose
`stress_field`, returning the shared borrowed `StressField3d` from
`compressible::spatial::molecular`. This evaluates the original accepted endpoint
once with its actual transport law, caloric/species owner, velocity gradients,
compatible stress discretization, explicit bulk viscosity, wall/reservoir
conditions and any explicitly selected SGS closure. The mixture entry passes
the caller's original `DissipationControl` to that operator. No constitutive law,
SGS completion or velocity gradient is inferred in the acoustic consumer.

The view reads original source-decoded primitives and exact conserved momentum
density. `face_traction_pa(axis, coordinates)` returns the actual native
`tau[:,axis]` in Pa, opposite to the viscous conservative momentum flux.
Its axis is the positive Cartesian direction; the physical outward sign is
applied by the surface consumer. It contains neither inviscid pressure nor
advective momentum, and thermal/species energy flux is never treated as force.
The same face buffers supply advancement and observation, including all three
global momentum components of a multicomponent face.

`cell_stress_pa` averages the opposing face columns to the cell center and
symmetrizes the off-diagonal pairs, returning xx, yy, zz, xy, xz, yz in Pa.
These are undoubled physical shear entries. Reconstruction is exact for affine
stress and second order in smooth resolved interiors; walls, material jumps
and modeled SGS retain their native limitations. Surface loading can use the
unaveraged face traction. Source consumers apply their original rotation and
physical volume or face-area quadrature, with the Lighthill sign `-tau_ij`
and FW-H loading sign `-tau_ij n_j`.

Observation reuses prepared native field storage and allocates no additional
field buffers; original provider-owned evaluation costs retain their contracts.
The view's borrow prevents a step or another evaluation while it is live.
Flow cells, physical clock, last accepted interval reports, boundary impulses
and accepted Lagrangian closure histories are unchanged. A Lagrangian endpoint
fit is uncommitted and the original filter is restored on success or failure.
A failed/cancelled query returns no view, though reusable constitutive caches
may change. The returned stress is instantaneous at the held accepted time,
not the preceding step's integrated load.

The accompanying `fs-aeroac` owner tests specify independent affine full-stress
and bulk-viscosity identities, all-axis face signs, rigid-motion and zero-law
limits, smooth periodic shear refinement, original multicomponent composition,
stationary tensor rotation, accepted-clock sampling and final-cancellation
publication. They establish the resolved extraction contract; unresolved sound
generation and completeness of the selected acoustic source region remain
separate physical modeling choices.

## Pressure-continuous cavitation populations

Both original radial owners expose `linearization_at_radius`: gas pressure,
restoring stiffness, viscous damping and the Keller--Miksis inertia/radiation/
drive-derivative terms are evaluated at the specified equilibrium of the SAME
retained gas inventory. The established reference `linearization()` arithmetic
is unchanged. New tangents refuse unstable or numerically unresolved marginal
stiffness, including a roundoff residual at the Blake fold.

Static continuation also handles the explicitly authored gas-only and vapor-only
limits. For `3 kappa = 1`, the restoring sign selects the monotone branch;
a flat marginal curve refuses. For `3 kappa < 1`, the finite turning point is
a maximum and the stable branch lies at larger radii. This does not relabel a
maximum as a Blake minimum.

`BubblePopulationClass::at_equilibrium_pressure` and
`BubblePopulation::at_equilibrium_pressure` prepare stable equilibria while
retaining original models, gas references, class order and absolute number
densities. Each actual radius supplies its new occupied volume and cached native
tangent. Root brackets, pressure residuals and iterations remain visible.
Work is bounded by the existing class and root budgets, with cancellation before
publication. Empty bins require no fictitious equilibrium. Failed roots,
overfilled volume, unstable tangents and cancellation return no new population.
The ambient shift from the ORIGINAL reference equilibrium is additionally bounded
by the native pressure/bulk-modulus ratio so repeated queries cannot walk frozen
liquid coefficients through an unbounded pressure excursion.

`PopulationEquilibrium::into_evolving` starts the original nonlinear models at
these prepared radii and zero wall velocity, then applies the native whole-radius
dilution-envelope admission. It does not reset an existing trajectory. This is
held-volume equilibrium preparation with frozen liquid properties and an authored
polytropic exponent, not closed-parcel compression, transport, nucleation,
thermal relaxation or diffusion. Frequency-domain consumers still admit their
own dilution, amplitude and interaction bounds.

Focused G0/G1/G3/G4 regressions compare native tangents with independently
differenced nonlinear accelerations, closed-form pressure equilibria, retained
gas/number inventories, histogram splitting, pressure-cycle repeatability,
stationary nonlinear initialization and transactional refusals. All 65 native
cavitation tests passed, including the original-reference pressure-excursion
guard, with the original numerical source and restricted dependency facades.
Full Cargo/DSR remains unavailable because the seven constellation siblings
and DSR/RCH are absent.


## Thermal diffusion in retained cavitation nuclei

`cavitation::thermal::IsothermalWallBubble` prepares the original native
Keller--Miksis bubble at its current stable equilibrium with a coherent
`GasState` evaluated at the same gas partial pressure and interface temperature.
The static inventory must be explicitly isothermal; the material's gamma
enters only the harmonic gas heat equation. The native transport owner supplies
D = k/(rho cp), rather than a fitted damping rate or an independent gas table.
`BubbleGas::pressure_at_radius` exposes the checked retained pressure law for
material preparation without reconstructing the original inventory.

For exp(-i omega t), spherical gas diffusion with a fixed-temperature interface
supplies kappa = gamma/[1 + 3(gamma-1)(z coth(z)-1)/z^2],
z^2 = -i omega R^2/D. This follows Commander and Prosperetti (1989), equations
19--28, with the Fourier convention conjugated. The low-frequency series
preserves small thermal loss, the middle branch uses scaled exponentials, and
the high-frequency branch retains the algebraic loss while omitting only a
sub-roundoff exponential correction. DC returns the original native compliance.
Frequency queries require no allocation, iteration, gas-source lookup or root.

The complete complex restoring pressure enters the SAME native KM operator,
including its pressure-derivative forcing and viscous acceleration contribution.
The impedance form avoids cancellation between low-frequency radiation terms.
Cycle-mean gas thermal absorption per squared fractional-radius amplitude is
-6 pi R^3 omega p_g Im(kappa); it excludes viscous loss and radiated/scattered
power. It is not a nonlinear thermal state or a complete heat-deposition account.

Preparation rejects stale or inconsistent gas snapshots, non-isothermal static
inventories, volatile water in the noncondensable gas, unstable tangents,
excessive gas density or vapor fraction, and an excessive equilibrium pressure
shift from the original frozen-liquid reference. Queries independently admit
gas omega R/c_g; liquid compactness cannot replace that condition. Ideal-gas
continuum conduction, small amplitude and an isothermal interface remain
explicit assumptions. Liquid thermal impedance, vapor mass/latent heat,
molecular relaxation and nonlinear thermal collapse remain outside this owner.

All 71 native cavitation tests passed in the focused original-source build,
including six new groups: independent spherical heat eigenmodes with a bounded
truncation tail, isothermal/adiabatic limits and branch continuity, retained gas
and DC compliance, the complete expanded KM operator, cycle-averaged gas work,
and material/physical-domain refusals. Native gas coherence and transport code
were included verbatim through restricted dependency facades. Full Cargo/DSR
remains unavailable without the seven constellation siblings and DSR/RCH.

## Nonlinear finite-gas caloric bubbles

The explicitly opt-in `caloric-cavitation` feature exposes
`cavitation::caloric::CaloricBubble`. It combines either original
Rayleigh--Plesset or Keller--Miksis liquid momentum kernel with one source-backed
`ClosedGasParcel`. Radius, radial velocity, gas internal-energy increment,
supplied heat and gas pressure work evolve together through the original
fixed-stack DP5(4) integrator. The original `GasCaloricBranch` integrates the
actual temperature-dependent Cv; no fixed polytropic exponent or surrogate gas
is substituted. Source/composition identity, mass and auxiliary query axes
survive every accepted endpoint.

The gas law is `p_g = m R_g T/V` and its energy equation is
`U_g' = Q_heat - p_g V'`. Keller--Miksis consumes its analytic derivative
`p_g' = -3 (1+R_g/Cv) p_g R'/R + R_g Q_heat/(Cv V)` at the same RK stage,
including both compression and signed heat input. The existing viscous
acceleration denominator, pressure-drive derivative and acoustic
back-reaction remain in the native shared equation. The original polytropic
owners use those same liquid kernels and retain their own gas laws.

Original-source energy reconstruction cancels common formation references
before arithmetic. The complete-window residual `Delta U_g + integral p_g dV
- integral Q_heat dt` and original-source/integrated-coordinate discrepancy
must both satisfy an explicit absolute joule budget. Inner caloric tolerance
must fit the declared local and complete-window energy scales. This is a gas
heat/work account, not a total compressible-liquid/radiation energy theorem.
The inherited incompressible liquid kinetic-energy output remains only a
leading-order diagnostic for Keller--Miksis.

Preparation amortizes continuous caloric-source admission across the complete
advance; radial stages allocate no source model, transport receipt or history.
Each caloric inversion retains its original bounded work/cancellation checks.
A final endpoint resolves the original material transport/source domain once.
`prepare_advance_controlled` returns a borrowing proposal for downstream
acoustic admission; a dropped proposal, source refusal, numerical failure,
physical boundary or cancellation leaves all accepted gas and radial state
unchanged. Heat and pressure callbacks must be deterministic and bounded and
must have their bandwidth resolved by the maximum native time step.

Both phases remain in explicit restricted models: frozen Newtonian liquid,
fixed-mass/composition ideal noncondensable gas, spatially uniform gas pressure
and lumped temperature. The source's declared water species is refused.
Original radius, liquid wall Mach, vapor/liquid density and selected KM
pressure/bulk-modulus limits apply, and every stage independently bounds actual
gas/liquid density ratio, gas wall Mach using its current Cv/T, and vapor
pressure fraction. The latter three ceilings are at most 0.1. These are
applicability restrictions, not experimental error bounds. Other gas phase
validity, incident compactness and the uniform-temperature approximation
remain caller conditions. The signed heat callback is an explicit external
power port, not a derived conduction closure. No vapor calorics/latent heat,
phase mass transfer, liquid thermal diffusion, chemistry, dissolution, shocks,
wall jets or erosion prediction is implied.

Focused G1/G3/G4 tests are authored to compare both radial equations with the
constant-Cv adiabatic polytropic limit, independently integrate the entropy and
energy of `Cp/R = 3.5 + b T`, shift the formation-energy datum, verify the KM
heat term, and exercise signed heat, preserved inventories, source-domain
refusal and transactional cancellation. They are numerical model checks, not
experimental validation of cavitation or heat-transfer coefficients.

The implementation remains opt-in. The focused original-source subset now
passes all 82 native cavitation tests, including the caloric groups and the
new conservative-coordinate regression. Numeric RK stages may sample a gas
energy smaller than one represented temperature increment while retaining that
energy coordinate; strict original-source endpoint reconstruction and both
energy audits remain unchanged. Full workspace/DSR execution remains outside
this focused result. The workspace selector is `cargo test -p fs-flux
--no-default-features --features caloric-cavitation --lib cavitation`.

## Finite surrounding-liquid conduction for cavitation nuclei

`cavitation::thermal::ConductingLiquidBubble` retains an admitted
`IsothermalWallBubble` and the SAME liquid's `ThermalFluidTransportState`.
Temperature, density and dynamic viscosity must match the retained liquid
snapshot exactly. Gas and liquid heat capacity/conductivity determine their
own diffusion scales; no fitted heat-transfer coefficient, replacement gas
inventory or second radial/acoustic equation is introduced. The underlying
`ThermalBubbleResponse` uses the existing native Keller--Miksis operator.

For exp(-i omega t), write `z_j^2 = -i omega R^2/D_j`,
`B = z_g coth(z_g)-1`, `F = 3B/z_g^2`,
`Y_g = 4 pi R k_g B`, and `Y_l = 4 pi R k_l (1+z_l)`.
Continuity of temperature and heat flux between the spherical gas interior
and unbounded liquid exterior gives `T_i/T_ad = Y_g/(Y_g+Y_l)` and
`kappa = gamma/[1+(gamma-1) F Y_l/(Y_g+Y_l)]`. These are the heat equations
and interface conditions of Prosperetti (1977), Sec. IV, in the uniform-gas-
pressure limit: <https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/JASA_Thermal_77.pdf>.
Infinite liquid admittance recovers the original fixed-interface-temperature
owner; vanishing liquid admittance recovers adiabatic gas. Finite positive
conductivity returns the original isothermal compliance at DC.

The response exposes complex interface temperature and signed oscillatory
outward heat rate per fractional-radius amplitude. Those first-order heat
fluxes are distinct from mean irreversible loss. Green's identity evaluates
gas loss as `Re(Y_g)|T_ad-T_i|^2/(2 T0)` and liquid loss as
`Re(Y_l)|T_i|^2/(2 T0)`, each per squared fractional-radius amplitude.
Their sum closes the original gas p--V loop
`-6 pi R^3 omega p_g Im(kappa)`. Neither viscous loss nor outgoing acoustic
power is included. The common response's thermal power therefore denotes
combined gas and liquid conduction for this boundary condition.

The existing gas low-frequency series supplies the real admittance without
subtracting nearly equal terms. Each interface share is divided directly,
and both nonnegative entropy powers are evaluated independently; liquid loss
is not subtracted from total work to recover a small gas loss. Preparation
is constant work and each harmonic query is allocation-free and constant
work, with the same gas compactness/refusal behavior as the base owner.

The infinite liquid has a fixed far-field temperature, frozen coefficients
and a thermally perfect interface. Nonlinear thermal histories, interfacial
resistance, overlapping thermal fields of different nuclei, thermally driven
liquid expansion, pressure-driven liquid temperature changes, vapor mass,
latent heat, dissolution and molecular relaxation remain outside this model.
The retained small-vapor bound alone does not prove negligible latent heat;
the caller must establish the noncondensable-gas approximation. This owner is
a harmonic constitutive response, not a nonlinear heat-transfer callback.

Focused tests cover coherent retained states and refusals, exact DC response,
isothermal-wall and insulating-liquid limits, weak-loss/high-contrast energy
closure and interface flux, independent radial integration of both
`k |grad T|^2/(2 T0)` entropy fields, and the original expanded native KM
operator and cycle-integrated gas work. They establish numerical/model
identities rather than experimental validation.

The focused original-source native build passed all 81 cavitation tests,
including the five finite-liquid conduction regression groups. This build
uses the actual changed sources and original numerical dependencies with a
restricted module/dependency closure; it is not a full workspace or DSR run.

## Implicit RCR storage and series inertance at resolved laminar ports

The native `bdf2_stokes_mixed_variable_step_windkessel_ports` advances the
original fixed-mesh Stokes momentum and passive lumped external loads at one
common endpoint. Outward flow Q is volume per depth; with compliance pressure x
relative to the same reservoir, the external laws are
`p=p_res+Rp Q+L D_t Q+x` and `C D_t x=Q-x/Rd`.
Rp and L are finite/nonnegative; Rd and C are finite/strictly positive.
The reservoir increment is held during each interval; external coefficients,
port ordering and the compliance reference remain fixed through both histories.

Fluid, compliance and series inertia share the original admitted BDF reduction:
BDF1 startup, equal-step BDF2 arithmetic, or variable-step BDF2 with ratio at
most two. A history pairs velocity moments with compliance pressures at one
time level. Previous/older port flows are integrated from those SAME native
velocity moments; independently supplied inertance-flow histories are not
accepted. The original positive fluid mass and viscosity remain in the
continuum solve. Eliminating storage produces a positive effective port
resistance, solved through the existing full coupled mobility response.
No endpoint flow is lagged and no profile is imposed. Original continuum
residuals and BOTH unreduced dynamic laws are checked on integrated endpoint
flows before any field or storage state is returned.

The report separates external C*x^2/2 and L*Q^2/2 storage, nonnegative
Rp*Q^2 and x^2/Rd resistor powers, reservoir and added-traction work, the
actual-interval storage change and a signed BDF energy defect. The effective
resistance's apparent loss is NOT physical heat. BDF2's defect can be negative;
neither ordinary single-state energy monotonicity nor arbitrary-step coupled
passivity is asserted. Original traction work and resolved-fluid heat remain
separate. The load is linear lumped storage, not a distributed wave, cavitation,
thermal update, moving interface or nonlinear external device.

Focused exact-source scalar tests cover the original two-state equations,
BDF1/unequal BDF2, flow reversal, zero series inertance/proximal resistance,
closed-port exponential relaxation and second-order convergence, original-law
residuals and independent energy accounting. A native FEM regression is authored
for the complete perfect-slip channel plus external storage/inertia, physical
pressure datum and coupled energy. Full native FEM execution remains pending:
the workspace dependency build was killed in pinned asupersync before fs-flux.

## Characteristic pressure outlets on Cartesian gas fields

`BoundaryCondition::CharacteristicOutlet { reference_pressure_pa,
reference_velocity_m_s }` connects the original acoustic characteristic law to
frozen and evolving-species Euler fields in two and three dimensions. Reference
velocity is signed along the positive Cartesian axis and must point outward
(or be zero) on its side. Reference absolute pressure must be positive.
The original source projects that pressure onto the local interior isentrope;
no remote temperature, constant gamma or replacement composition is inferred.

Both first-order and MUSCL traces retain departing entropy, all original
species, tangential velocities and complete kinetic energy. The same face flux
supplies conservative RK transfer and boundary impulse. Source-owned sonic
sampling and supersonic causality remain active; reverse contact motion refuses
before publication. Initial face wave bounds participate in the complete
unsplit interval, including transverse advection. Hydrostatic bounds distinguish
the configured pressure reference from the reconstructed physical face state.
Uniform 3-D bindings survive trial clones without per-step allocation.

Frozen molecular open-face behavior is unchanged. Planar evolving-species
molecular transport refuses incomplete baths. The volumetric mixture's explicit
`PressureOutletPolicy::ZeroDiffusiveFlux` may also close a characteristic outlet;
its default still refuses. This is an independently declared molecular closure,
not a thermodynamic reservoir implied by the acoustic condition.

The added G0/G1/G4 regressions cover all four/six orientations, both
reconstruction modes, variable heat capacities, complete shear/species energy,
absent species, outgoing pulse reflection versus fixed pressure, admission and
transactional refusal. The source-owned boundary kernels passed native tests
in a bounded source-linked build. The added Cartesian integration regressions
were not executed: the workspace's constellation dependencies were unavailable
and the execution environment subsequently went offline.

This locally one-dimensional, smooth-wave outlet does not establish exact
nonreflection for oblique waves, shocks or vortical boundary interactions. It
does not fix instantaneous outlet pressure, supply backflow thermodynamics,
resolve a finite exterior gas inventory or certify discrete entropy stability.

## Source-normalized physical RCR and inertance coupling

The homogeneous `MaterialFluxSystem` and frozen-cell `MaterialStokesField`
both expose the dynamic Windkessel BDF entry. The original material owns fluid
density and viscosity; external resistance, compliance and inertance are
explicit device data. No bulk modulus or gas compressibility is invented from
transport properties. With Q in m²/s per unit depth, physical resistances are
Pa s/m², series inertance is Pa s²/m² and compliance is m²/Pa. A 3-D network's
resistance/inertance multiply by modeled depth; its compliance divides by depth.

One fixed reference density normalizes the entire system: pressure, both
resistances and series inertance divide by rho_reference; compliance multiplies
by rho_reference. Native velocity and all Q histories are unchanged. Each
spatial material cell retains its original rho/rho_reference mass and
mu/rho_reference stress; reference viscosity never replaces those sources.
The original source snapshots, weighted kinetic energy and boundary reactions
survive without re-query. Positive quantities or signed pressure histories that
underflow away during conversion are refused before publication.

All reported external energies/powers and pressure-law residuals return to SI
with that same reference density; flow residuals remain m²/s. Compliance and
inertance storage, both resistor losses and the signed BDF defect remain
separate from resolved-fluid viscous heating. No automatic temperature, EOS
pressure, material-interface transport or cavitation update is implied.

Two focused native material regressions are authored for physical density swaps,
arbitrary reference-density invariance, BDF1/BDF2 histories, energy units,
source retention, normalization refusal, and layered Couette flow while its
external compliance charges. Their full FEM/material runtime execution remains
pending at the pinned-asupersync build gate; the six shared exact-source scalar
BDF/RCR tests passed. These checks support the stated numerical model, not
experimental qualification of an authored external compliance.

## SSPRK2 inventory transport with nonlinear native SGS momentum

With `experimental-sgs`,
`FluxSystem::density_momentum_ssprk2_bdf1_sgs_step` returns
`DensityMomentumSspRk2SgsStep`. It advances density and every carried inventory
through the original configured reconstruction and bounded dyadic SSPRK2
transport, then solves the original nonlinear molecular-plus-Smagorinsky, WALE
or Vreman momentum endpoint.

The actual RK-weighted compatible H(div) mass flux, final positive density and
specific inventories, original old inertia and macro interval remain fixed
through the constitutive iteration. The flux is never reconstructed from the
final density or an intermediate Picard velocity. Molecular source queries
receive that final transported material state and the current candidate strain;
force and boundary data are sampled once. Retained endpoint samples and
molecular heating, modeled SGS transfer and numerical dissipation remain
separate. Transport, source, nonlinear-budget or cancellation refusal publishes
no partial candidate and does not modify caller histories.

Only the frozen-advector density/inventory ODE has second-order time
integration. Momentum remains BDF1; the existing planar Estimated SGS scope
and model-form limitations remain. This composition supplies no new closure,
EOS, species diffusion or automatic heat deposition.

Three focused regression groups are authored for an independent analytic
two-cell RK exchange, final-state material/closure sampling, bounded
complementary inventories, native discrete kinetic accounting, the exact
zero-SGS limit, and transactional refusal/cancellation. Scoped pinned formatting
and whitespace checks passed. Native execution remains pending because the
pinned asupersync dependency exceeded this execution environment's available
build resources before the changed crate could compile.

### Endpoint-advector density/momentum coupling

`FluxSystem::density_momentum_coupled_bdf1_step` iterates the native endpoint
velocity as the conservative density advector. Each trial starts from the SAME
immutable old density, specific inventories and momentum history. It reuses the
original directional transport, compatible H(div) mass flux, distinct old/new
inertia, physical force/viscosity conversion and BDM pressure/reaction assembly.
Rejected Picard updates do not consume physical time or transported inventory.

The existing bounded Picard owner requires the unrelaxed velocity correction
and both original/row-scaled residuals reassembled at the returned endpoint.
The accepted density and ALL carried inventories are exactly the original
transport evaluated with that returned velocity; the retained pressure,
reactions and kinetic observations use the same endpoint mass flux. A failed
source, missing inflow donor after reversal, excessive trial CFL, exhausted
nonlinear work or cancellation publishes no partial coupled state.

This remains first-order endpoint-advector coupling, with old-state donor
densities and specific values. It is not fully implicit density integration
and retains the original CFL/divergence admission, fixed geometry and passive
inventory semantics. No phase-volume source, EOS, heat diffusion, full coupled
adjoint or automatic numerical-heating deposition is added. Callbacks must be
pure through all trials. Storage does not accumulate iteration history; native
assembly/solve allocation and synchronous cancellation limits remain unchanged.

G1/G3 tests in `tests/density_momentum/coupled.rs` cover acceleration from rest
with same-interval density/inventory transport, reference-density invariance,
closed unequal-density mass/carried/kinetic balances, and atomic nonlinear,
CFL, missing-donor and cancellation refusals.

## Multicomponent liquid miscibility with Maxwell–Stefan transport

`miscibility::multicomponent` connects the original complete regular-solution
model to the original molar-frame Maxwell–Stefan evaluator on a closed finite-volume
graph. Every component has an independent compensated inventory. Canonical face
transfers debit and credit the same amount; each returned candidate checks every
global species balance, every fixed cell capacity and the actual original
homogeneous Gibbs-energy change. Unequal molecular masses are retained: zero net
molar flux does not assert zero diffusive mass flux in the barycentric frame.

Each face uses the finite original chemical-potential difference, including every
nonideal interaction. The common chemical-potential gauge is removed through
antisymmetric Gibbs–Duhem pair accumulation at the arithmetic face composition.
The full coupled friction system is solved with its original equation/frame
residual checks. This admits cross diffusion, including an individual component
moving uphill while total chemical free energy decreases; independent Fick laws
are not substituted. The arithmetic face rule is a consistent two-point spatial
approximation, not exact finite-jump Fick diffusion.

Preparation canonicalizes faces, shares explicitly authored complete pair
diffusivity matrices, computes the fixed-temperature interaction matrix once,
and admits a sufficient global convexity bound using the original tangent
spectrum and its numerical uncertainty. Local stability at one composition and
pairwise binary miscibility alone do not establish this admission. The bound is
conservative and numerical, not an interval certificate or a general phase flash.

The explicit synchronous step has no hidden subcycles: positivity or actual
free-energy failure refuses atomically. Common strictly positive species support,
constant equal partial molar volumes and a uniform prescribed temperature are
required. Absent-species contacts, nonconvex multicomponent phase separation,
local thermal evolution, variable molecular-diffusivity sources and momentum
coupling remain outside this owner. Time/space refinement remains necessary.
The instantaneous constitutive dissipation rate is a diagnostic; it is not an
integrated explicit-step dissipation certificate.

Numerical storage is contiguous by cell/species, with one reused original
Maxwell–Stefan workspace per sweep and a shared matrix table. Admission charges
retained payload and conservative simultaneous transaction/adaptive scratch,
excluding allocator and species-string overhead. Geometry/species limits,
checked byte arithmetic and cancellation bound the work.

Validation: all 97 numerical miscibility unit tests passed on the pinned
nightly-2026-07-06 compiler in the focused actual-source dependency closure
with one compiler frontend thread. Eight independent new regressions cover
binary reduction, ternary uphill diffusion, unequal masses/volumes, canonical
faces, positive traces, convexity/source/workspace refusals, large-step energy
and positivity refusal, cancellation and deterministic replay. This is focused
numerical execution, not a full-workspace or physical-validation certificate.


## Nonlinear laminar flow with implicit external storage

The native convective Windkessel entries now solve RCR/compliance/series-inertia
loads inside the ORIGINAL bounded Navier–Stokes Picard iteration, including
candidate-dependent generalized-Newtonian viscosity. They reuse the same
BDF1/startup or variable-step BDF2 history owner and dynamic elimination as the
linear entry. Each current Oseen/constitutive operator supplies its own coupled
port mobility; an initial or lagged mobility is not used for final acceptance.

The existing nonlinear port owner admits either prescribed integral flow or
finite-resistance laws. It keeps pressure increments as additional unknowns and
relaxes them with the complete field. At an acceptance candidate, viscosity and
convection are reassembled at that velocity, the ACTUAL accepted pressure load
is added to the original RHS, and both original/scaled nonlinear residuals are
checked. The same integrated resistance-law evaluator checks the accepted field.
The dynamic wrapper then checks both unreduced capacitor and inertance laws
using the accepted BDM flow. Failed nonlinear budgets or source laws return no
partial fluid/compliance state.

Natural/dissipative backflow and explicit incoming reservoir momentum retain
their original meanings. External storage, resistor losses and signed BDF
defect do not include convected kinetic energy, fluid viscous heat, or added
backflow-policy work. No combined arbitrary-step nonlinear passivity, heat
deposition, turbulence or material-memory claim follows.

Focused native regressions are authored with nonzero convection, live viscosity
updates, forward/reverse port flow, BDF1/unequal BDF2, pressure relaxation,
independently reassembled momentum equations and nonlinear-budget refusal.
Full FEM execution remains blocked by the prior pinned-asupersync SIGKILL;
format/parsing and bounded static API review are distinct from runtime evidence.
The original six scalar BDF/RCR model tests were previously executed; their
unchanged leaf model does not establish new nonlinear integration test coverage.

## Consistent Newton integration of planar SGS stress

`PlanarSgs::{solve_steady_newton,bdf1_newton_step,bdf2_newton_variable_step}`
connect the authored Smagorinsky, WALE and Vreman closures to the original
damped Newton owner. Smagorinsky retains the scalar constitutive tangent
`d(nu_eff)/d(gamma) = (Cs*Delta)^2`; WALE and Vreman use the full velocity
gradient response described below. The native Jacobian differentiates cell
stress, harmonic SIP weights, penalties, boundary loads and convection.
Filter widths and molecular viscosity remain the prepared cell fields.
Neither an artificial steady mass nor a frozen eddy viscosity is introduced
to obtain convergence.

These entries prescribe velocity on every boundary. The BDF2 history carries
its actual previous interval, admits a new/previous step ratio in (0,2], and
uses the real previous velocity as the Newton initial iterate. Missing older
history is exactly BDF1 startup. Accepted SGS diagnostics belong to the same
velocity endpoint and remain instantaneous modeled transfer, not added heat.

At exactly zero strain, Smagorinsky uses the zero generalized derivative of
the strain norm. Viscosity-dependent DG penalties may still be nonsmooth when
a tangential jump remains. WALE and Vreman use their full-gradient tangents,
including the shared kernels' explicit singular-state selections; their
original Picard paths remain available. Original work bounds, Armijo
backtracking, endpoint residual checks, and read-only history apply.
No new cancellation mechanism or three-dimensional turbulence qualification
is implied.

The initial Smagorinsky regressions covered an independently reassembled
stiff nonlinear momentum residual, unequal-step quadratic shear, exact
zero-SGS native reduction, and tangent/history/work refusals. The original
NS, Jacobian, sparse, Krylov and dense numerical owners passed the focused
pinned-Rust validation. Current all-model validation is described under the
full-gradient extension below. These focused builds omit unrelated
material/thermal/density adapters; they are not full workspace or runtime
integration qualifications.

## Original material rheology in convective dynamic-port solves

Both fixed-viscosity material Navier–Stokes and the original source-rheology
resolver expose the nonlinear BDF Windkessel entry. The latter retains the
existing `RheologyField` query/cache owner: accepted cell materials belong to
the final actual shear, while source refusal takes precedence over provisional
linear/nonlinear results. It consumes the original Carreau-Yasuda/temperature
source without copying a constitutive fit or replacing it with frozen viscosity.

The linear and nonlinear material entries share one physical-load/history
normalizer and one result converter. Pressure, resistances and inertance divide
by the SAME reference density; compliance multiplies by it; velocity and Q
histories remain native SI moments. Every original source preserves that fixed
density. Fluid stress, inertia, convection, pressure/reactions and instantaneous
bulk heat use the original accepted material state. External resistor losses,
compliance/inertance storage and signed BDF defect remain separate; no coupled
thermal or thermodynamic-pressure update is inferred.

Focused native tests are authored with actual ThermalCarreauYasuda source owners,
its Newtonian limit, nonzero convection, candidate-dependent rates, BDF1 and
unequal BDF2, under-relaxation, physical density swaps, accepted source identity,
independent physical viscosity/heat/storage checks and source-domain refusal.
The new FEM/material tests remain UNEXECUTED at the pinned-asupersync build gate.
Rust formatting/parsing and bounded static signature/ownership checks passed;
the prior unchanged scalar BDF/RCR tests are separate evidence.

### Adaptive multicomponent physical intervals

`miscibility::multicomponent::adaptive` compares one full explicit step with
two half steps, accepts only the fine trajectory, and independently controls
every species' mole-fraction error. The error estimate includes the conserved
amount low words, so a trace or sub-ulp transfer cannot disappear because a
dominant fraction rounds to one. A relative allowance that underflows to zero
requires exact agreement. This is a local temporal estimate, not a global
space/time error certificate.

Positivity and free-energy refusals trigger bounded refinement. Invalid original
constitutive data, unsupported states, arithmetic failures and cancellation
remain explicit refusals. No partial interval is published on exhausted work or
accuracy. The actual representable clock endpoints are retained, including the
last remainder below the nominal minimum step.

Physical face receipts include only accepted fine steps. Work reports distinguish
completed numerical steps and completed original Maxwell–Stefan evaluations;
refused incomplete solves can perform additional work within the attempt bound.
The original single species/capacity/free-energy allowance is checked again
over the entire requested interval rather than accumulated per accepted step.

The complete initial rates are shared by the full and first-half trial and
retained across smaller-step retries at the same field. A normal accepted
attempt therefore needs two constitutive sweeps rather than three; the midpoint
still resolves the full original nonlinear friction system. Cached rates are
invalidated on every accepted field change. Thermodynamic/material states use
the published representable component amounts; the carried amount low words
additionally enter conservation and temporal-error accounting.

Validation: all 103 focused numerical miscibility unit tests passed, including
six adaptive regression groups. These cover independent small-perturbation
relaxation, accepted-only ternary receipts, positivity-driven refinement and
rate reuse, work/accuracy/cancellation refusals, independent positive-trace low
words, and exact short-remainder endpoints.



## Phase-changing saturated-water vapor bubbles

`cavitation::vapor::SaturatedVaporBubble` evolves geometric radius, liquid interface velocity and saturation temperature on the original fixed-stack adaptive radial integrator. It binds `SaturatedWaterCavitation` and the original `SaturatedWaterVapor` material jet; no ideal-gas vapor EOS or fitted evaporation coefficient replaces IF97. Positive interfacial mass flux j satisfies `Rdot = U_l + j/rho_l`. The common incompressible liquid momentum kernel retains both this Stefan slip and vapor recoil `j^2(1/rho_v-1/rho_l)`. The impermeable gas-bubble path remains its zero-slip specialization.

Vapor occupies the geometric sphere with `m_v=rho_v(T)V` and `U_v=m_v u_v(T)`. The explicit liquid-side heat port obeys the leading thermal interface law `dU_v=Q_l dt+h_l dm_v-p_v dV`. The original saturation mass and energy equations form a two-by-two scalar solve for temperature rate and mass flux; neither heat nor mass flux is numerically differentiated. Signed transfer mass, heat, liquid enthalpy, vapor pressure work and an omitted-interface-energy scale are integrated on the same adaptive stages. Independent original-source mass and reference-invariant energy checks admit each whole window before a borrowing proposal can commit. Original-datum accounts remain separately visible.

The physical acoustic source is displaced-liquid flow `Q_l=4 pi R^2 U_l` and its analytic derivative, not geometric bubble volume rate when phase exchange is present. Stage checks cover radius, actual phase density ratio, liquid/interface/vapor Mach, vapor pressure-relaxation scale, actual original liquid density/viscosity/sound-speed drift, pressure departures, omitted interface enthalpy scales and flat-interface saturation corrections. Surface tension remains explicitly frozen for the declared interface environment. These are applicability screens; the model does not claim an exact nonequilibrium-interface error bound.

This lane assumes a spatially uniform, locally saturated pure vapor and an incompressible frozen bulk liquid. It does not supply noncondensable gas, kinetic phase resistance, nucleation, resolved liquid conduction, finite liquid thermal storage, Kelvin equilibrium, shock/collapse jets or phase-changing Keller-Miksis dynamics. Source: Sullivan et al., JFM 948 (2022), Eq. 2.1, with the opposite mass-flux sign convention. All 87 focused native cavitation tests passed, including five new groups covering stationary equilibrium, independent differential balances and recoil, nonlinear heating/cooling, refinement/analytic liquid-flow derivative, and failed/cancelled proposal rollback.

### Original liquid materials for multicomponent spatial mixing

`miscibility::multicomponent::material` binds an accepted numerical field to
its original `MulticomponentLiquidSolutionLaw`. The complete thermodynamic
model, canonical species and all original molar masses must agree. Independently
reconstructed identical numerical laws are admitted; the exact borrowed
material law is then retained through every candidate.

Every returned cell is resolved through that original law at the evolved full
composition and prescribed temperature, preserving fixed pressure, chart axes,
selected property claims and source refusals. Source mass divided by original
density must agree with the rigid numerical volume within the explicit retained
allowance. A late source-domain, density/volume or cancellation refusal publishes
no candidate. The adapter qualifies completed endpoints, not every internal
adaptive trial, and performs no local phase flash.

The original numerical full-pair diffusivity matrices remain authored data.
Binary-only composition charts and stationary-solvent aqueous diffusion sources
are not reinterpreted as general multicomponent liquid diffusivities. Material
receipt allocation has its own explicit cell-count bound.

Material regressions were source-reviewed and formatted but not executed in this
environment: the actual-source material dependency closure remains blocked by a
SIGKILL while compiling the pinned `asupersync` dependency, including the prior
single-thread compiler attempt. This is not a passing material build claim.

The finite-face force also retains the original abundant-carrier complement
logarithm when a carrier fraction rounds to one. Range-safe binary scaling
preserves representable trace driving forces across extreme geometry ratios;
a truly unrepresentable nonzero force refuses explicitly. The ordinary force
path retains direct arithmetic. Four independent precision regression groups
bring the focused numerical suite to 107 passing tests; this does not change
the separate unexecuted material-test status above.

## Compact stagnation reservoirs on existing Cartesian gas bindings

`with_stagnation_reservoirs` selects the original one-sided compact-reservoir
law for uniform held reservoirs on `Euler3d`, `MixtureEuler2d` and
`MixtureEuler3d`. Each selected binding must have exactly zero Cartesian
velocity. Its original density, pressure and species define upstream stagnation
temperature and enthalpy through the same source; no independent gamma, caloric
fit or discharge coefficient is supplied. Unselected sides keep their existing
complete-state Riemann condition. Empty selection restores that condition.

The original gas law selects acceleration, choking, discharge or reversal.
Incoming gas carries the tank's canonical species and stagnation enthalpy once,
with zero shear. Outgoing gas retains interior entropy, composition, all shear
momenta and complete kinetic energy. Both first-order and MUSCL traces enter
the same source law and accepted RK boundary impulse. The sampled normal wave
and transverse advection bound the full unsplit interval, including the shared
initial bound used by reactive and force-split owners. Hydrostatic face traces
retain their original source and face-pressure admission.

Uniform source energy and selection metadata are cached once; clone paths
retain them, and face evaluation introduces no heap allocation. Selection
refuses duplicates, missing uniform bindings, moving tanks and cancellation.
Spatial or timed exterior profiles remain complete-state Riemann boundaries.
The reservoir is an externally maintained infinite store, not a coupled finite
inventory or a resolved nozzle/contraction geometry.

Stagnation data do not specify a molecular physical-face bath. Mixture molecular
transport and explicitly enabled frozen reservoir bath exchange refuse this
combination. The frozen field's existing default zero open-face molecular
exchange remains available. Molecular thermodynamics is not inferred from p0/T0.

Four native regression groups cover all relevant side orientations and both
reconstructions, source-gauged tank enthalpy, original donor/absent species,
outgoing shear, hot-reservoir CFL, selection, cloning and cancellation. A bounded
actual-source native inviscid build passed these four groups and all six prior
characteristic-outlet field groups (10 passed, zero failed), including the
outgoing pressure pulse. The build retained the original gas EOS/characteristics,
hydrostatic preparation, reconstruction, forcing, RK and conservation owners;
unrelated chemistry/molecular wrappers were excluded. Only the new planar
molecular-admission assertion was omitted from this inviscid build. The source
gas adapters separately passed 18 actual-source boundary tests. Full crate,
workspace, molecular/chemistry integration and DSR verification remain pending.

### Prescribed stirring of multicomponent liquids

`miscibility::multicomponent::advection` transports every original component
through a closed prescribed volume flow. Forward and reverse streams remain
independent, so zero net face volume does not suppress physical mixing.
Synchronous donor updates carry each amount and its signed low word; no
dependent species is reconstructed. Prepared volume continuity, outgoing-volume
CFL, bounded substeps, common owner identity and the original whole-interval
species/capacity/free-energy checks govern publication.

The original convex homogeneous law makes donor mixing dissipative within its
actual numerical allowance. This does not identify chemical free-energy loss
with heat or mechanical stirring work. Unequal molecular masses remain in the
molar/volume frame. No interface work, momentum feedback, local thermal evolution
or variable diffusivity source is implied. Donor transport is first order and
requires spatial refinement.

`multicomponent::transport` applies advection then the complete original
Maxwell–Stefan operator over ONE physical interval. The second stage consumes
the full first-stage composition. Separate stage receipts, combined species
transfers, whole-window balances and a stage-energy accounting residual remain
available. The stricter stage allowance applies again to the complete endpoint.

`multicomponent::transport::adaptive` controls the complete composed operator,
including splitting error. Every coarse/fine trial has exactly one donor step;
its CFL cap prevents internal subcycles from aliasing the error estimate.
Independent component amounts and low words enter the same temporal norm as
molecular-only evolution. Only accepted fine paths contribute to advection,
diffusion and total physical transfer receipts. Requested endpoints and short
final remainders are retained. The suggested next step remains CFL-bounded even
when a short final interval is below an otherwise incompatible minimum step.

Prepared arrays and substep scratch are reused, and donor mole fractions are
resolved through the original model once at the final endpoint. The workspace
admission additionally charges 512 fixed bytes per face for prepared flow and
simultaneously live nested report headers, independently of species count.

Validation status: ten new numerical stirring regression groups were authored
and source-reviewed, covering independent circulation/counterflow solutions,
unequal volumes/masses, same-interval coupling, complete receipts, CFL aliasing,
time/work/cancellation refusal, replay and final remainders. They have NOT
completed execution. Initial compile import errors were fixed, then the focused
test build was SIGKILLed; a reduced per-package code-generation-concurrency
attempt was followed by loss of the execution environment. Reviewed source was
recovered through GitHub without a new successful compile or test run. The last
completed numerical result remains 107 passing tests BEFORE this extension.

### Original materials for prescribed multicomponent stirring

`miscibility::multicomponent::material::transport` retains the exact borrowed
original liquid source across prescribed advection, one-interval stirring plus
molecular diffusion, and complete adaptive stirring/mixing windows. Each wrapper
delegates the numerical operation and refreshes every final cell through the
existing source resolver exactly once at the complete endpoint.

Original pressure, chart axes, species/molar-mass mapping, property claims,
resolved-cell budget and source-volume allowance persist. The complete
numerical report is preserved unchanged. A late source-domain, volume mismatch
or final cancellation refusal publishes no field or partial physical report.
This qualifies final endpoints, not internal adaptive trials, and never
replaces spatial compositions with phase-equilibrium flashes.

Two additional integration regression groups exercise all three wrapper paths,
independent final property resolution, source identity, volume/clock/report
retention, replay and late failures. Their source/API review and pre-outage
formatting completed, but compilation and execution remain unverified: the
earlier original-material dependency build was blocked by asupersync SIGKILL,
and the later numerical stirring build lost the execution environment. These
claims do not promote material execution or the ten new numerical stirring
tests to passing status.

## Laminar devices connected entirely to external hydraulic loads

`StokesAssembly::prepare_pressure_response` retains the checked original
operator, baseline field and full native pressure responses without inverting
the port mobility. Every natural-traction edge in a component may belong to a
controlled port. A common pressure shift leaves volume flow unchanged, so a
singular mobility is physically expected; a single port on an otherwise
impermeable component can have zero mobility.

Finite resistance and BDF-eliminated RCR/inertance loads factor only
`I - diag(R) M`. Their reservoir laws select the uniform pressure level while
the original traction shape, tangential traction and pressure reactions remain
intact. No pressure pin, regularization or dummy open edge is introduced.
The original continuum operator must still be solvable: the reduced method
does not stabilize a singular free rigid-motion mode. Actual singular coupled
systems and invalid geometry/loads still refuse without publishing a result.

Prescribed-flow control remains a separate typed cache with its original
uncontrolled-datum and invertible-mobility requirements. Its finite-resistance
method forwards to the same pressure-response owner for compatibility. Linear,
nonlinear convective/rheological and SI material load entries all use the
relaxed pressure-response preparation. Original endpoint momentum residuals
and laws evaluated from the accepted integrated native fluxes remain enforced.
External loss/storage and the signed BDF energy defect remain separate from
bulk fluid heat.

Four focused production-reduction tests cover rank-deficient and zero
mobility, common and disconnected pressure shifts, series resistance, reversal,
conservation/power and actual coupled singularity/refusal. Native FEM tests
cover both ends of an inertial channel, a single sealed-volume port and two
RCR states with unequal-step history. Material regressions cover reference
density invariance and both ports on the original convective Carreau source
path. These new tests were authored but not executed before the execution
environment became unavailable. The previous local core parsed with rustfmt;
the connector-reconstructed final sources have bounded static review only.
This is not a native FEM/material compilation or runtime validation claim.


## Consistent Newton SGS with mixed boundaries

`PlanarSgs::solve_steady_mixed_newton`, `bdf1_mixed_newton_step`, and
`bdf2_mixed_newton_variable_step` connect the Smagorinsky, WALE and Vreman
constitutive responses to the existing bounded Newton/line-search owner
with velocity, Navier-slip, open-traction and reservoir boundaries. The
native `FluxSystem` methods retain the same mixed-boundary Newton path for
general admitted scalar viscosity tangents.

The Jacobian differentiates both transported and advecting velocity using
the original face flux: natural traction retains interior momentum on
either flow direction; dissipative backflow differentiates
`max(u.n,0)u`; reservoirs supply their retained incoming trace;
prescribed-normal slip retains its original interior transport. Boundary
samples are admitted once and reused throughout Newton and endpoint
reporting. The same accepted velocity, viscosity, mass history and boundary
policy determine pressure datums, stress reactions, advective transport
and SGS diagnostics. Steady components retain their original velocity/slip
support requirement; physical transient mass can anchor all-open components.
Unequal-step BDF2 uses the original history owner and actual intervals,
admitting new/previous ratios in (0,2]. Neither history is modified.

The shared closures retain their documented generalized derivatives at
singular gradients; these do not establish classical differentiability.
The `Estimated` planar qualifications remain. No additional cancellation
context, wall model, generated turbulent fluctuations, heat deposition,
global uniqueness or whole-solver backflow stability is implied.

The original Smagorinsky extension passed 21 selected planar/steady/Newton
regressions in a focused original-source Rust build. Its five groups checked
the directional
Jacobian against the independent public primal assembly across twelve
boundary/viscosity combinations, including face-flow reversal; fifteen
analytic steady/BDF1/unequal-step BDF2 pressure/traction/virtual-work cases;
exact disabled-SGS agreement; pressure-driven open flow with independent
endpoint residual and frozen-SGS negative control; and all-open transient
acceleration with steady/history/unsupported-model refusals. In the
pressure-driven fixture, Newton reached residual 3.8253e-14 in seven
iterations with six backtracks; the original Picard path exhausted thirty
iterations at residual 0.9934. This is fixture-specific convergence evidence,
not a general runtime or turbulence-validation claim. The original
numerical kernels were compiled without numerical stubs; material/flow-control
adapter tests and full Cargo/workspace validation were not run.

### Transported-density endpoint with integral volume-flow ports

`FluxSystem::density_momentum_coupled_bdf1_flow_control_step` solves each port's unknown pressure together with the original endpoint-advector density/BDF1 momentum equations.
Every Picard trial starts from immutable old density, carried inventories and velocity; accepted transport and compatible mass flux are rebuilt at the returned endpoint, including configured density reconstruction.
The native pressure-response owner constrains total outward volume flow without prescribing its velocity profile; disjoint natural-traction ports require one uncontrolled traction datum per controlled component.
Targets/tolerances use m²/s per unit depth; multiply actual rates by declared depth for m³/s. Physical mass exchange uses actual transported donors, never the reference density.
Port offsets are returned in Pa and retained in the same pressure-adjusted reactions; original traction and Open/OpenReservoir momentum semantics remain intact. Reference density only conditions the common equations.
Every trial retains original donor admission, CFL/divergence/storage limits and cancellation. Refusal publishes no density, inventory, flow or controller state; frozen pressure responses are rebuilt when candidate coefficients change.
Focused native regressions cover heterogeneous accelerating transport, exact accepted endpoint accounting, physical pressure/work under reference-density and depth changes, actual density-dependent inertance, reversed donors, port admission and atomic cancellation.

## Finite pressure tanks on planar frozen gas fields

`Euler2d::with_finite_reservoirs` binds source-owned stationary `Plenum`
inventories directly to selected physical transmissive sides. All faces of
one side share its independent rigid, well-mixed tank. The original frozen
caloric law supplies stagnation inflow, backpressure outflow, choking and
reversal; outflow retains transverse kinetic energy and incoming energy
thermalizes in the tank. No fixed gamma or independent pressure law is added.

Tank and field advance together in the existing SSP RK2 transaction. The
first weighted boundary capture forms the full forward tank stage, and that
evolved inventory drives the second residual. The final tank receives the
same native half-weighted mass and complete total-energy impulses that leave
the field. The fixed tank support receives the Cartesian momentum impulse
without mechanical work. Flow enthalpy and pressure work are counted once.

Initial admission and both residuals include the compact aperture-wave rate
`0.5 sum(A lambda)/V` alongside the field's unsplit acoustic CFL. Deterministic
compensated face reductions and fixed four-side stage arrays require no warm
tank allocation. The source-owned plenum candidate admits both caloric
endpoints and refuses exhaustion or unrepresentable storage without clipping.
Accepted cells, tank inventories, clock, loads and reports publish together
after the final poll. Explicit body-force work stays separate. Molecular
split preparation and retries preserve accepted tank inventories; molecular
exchange at these open faces remains zero.

The public tank report retains original/accepted inventories, complete
transfer and mass/energy residuals in the source gauge. Independent tanks
may fill or drain. These compact volumes do not resolve internal waves,
neck inertia, mixing gradients, heat transfer or moving volume, and one tank
cannot yet span multiple sides. Retained hydrostatic references refuse in
both attachment orders until tank elevations and port-potential transport
are explicitly declared; ordinary body forcing remains supported.

Seven focused G0/G1/G4 regression groups cover all four orientations and both
reconstruction orders, mass and full-energy conservation with shear, exact
equilibrium and empty-binding parity, tiny-tank initial CFL, common-stage
time refinement, independent opposing tanks with native body work, late
cancellation/retry and source/gravity admission. All seven passed
in a bounded actual-source inviscid build. Time-refinement errors were
3.2199504e-4, 7.9003992e-5 and 1.9365930e-5 at 16/32/64 steps, against a
512-step reference. The native caloric, flux, RK, forcing and hydrostatic
sources were retained. Molecular split hooks were parsed and source-reviewed
but not executed; full Cargo/workspace validation remains separate.

## Finite compact tanks on volumetric gas boundaries

`Euler3d::with_finite_reservoirs` converts selected existing stationary uniform
reservoir bindings into independent finite rigid tanks. Each side shares one
well-mixed plenum across its faces. Initial density, pressure, composition and
energy reference come from the original field source; the original native gas
`Plenum::exchange_candidate` owns inventory updates and caloric admission.
No second EOS, fixed gamma or held-pressure approximation is introduced.

Tank and field advance in the same SSP RK2 transaction. The first native face
impulses form a full forward tank stage, whose evolved pressure drives the
second residual. The accepted tank receives the exact opposite weighted mass
and total-energy exchange of its field boundary. Incoming shear kinetic energy
thermalizes. The same Cartesian momentum transfer is reported at the fixed
tank support, which does no work; no extra enthalpy or pressure-work term is
added. `finite_reservoir` exposes accepted storage and
`finite_reservoir_report` exposes transfer and mass/energy roundoff residuals.

Initial and stage CFL admission include the original compact-volume rate
`0.5 * sum(area * signal_speed) / tank_volume`, so small tanks cannot evade
stability control through large field cells. Source/positivity refusal uses
the existing bounded retry owner without clipping or a new tolerance knob.
Side transfers use deterministic compensated sums. Fixed six-side stage arrays
and cached source traces add no warm stepping allocation. Empty bindings keep
the existing native field path. Finite bindings retain their compact law.

All accepted field cells, tank inventories, exterior caches, time, impulses and
reports publish together after the final cancellation poll. Existing held body
forcing and other-side timed profiles share this same transaction. Molecular
split retries restore the accepted finite inventories and exterior caches with
the accepted field; default zero open-face molecular exchange remains available.
Explicit molecular bath exchange still refuses stagnation data. Moving-grid
admission continues to refuse bound reservoirs.

Six focused G0/G1/G4 regression groups cover filling/draining with all six
orientations and both reconstruction orders, complete energy/support impulse,
exact equilibrium and empty-binding parity, tiny-tank initial CFL, common-stage
RK2 time refinement, final-poll rollback and native force/timed-profile work.
All six groups passed in a bounded actual-source inviscid build, using the
original gas/Plenum, reconstruction, hydrostatic preparation, moving-grid, RK,
forcing and timed-profile code. The five source Plenum tests also passed.
Molecular split retry integration was source-reviewed but not executed; full
Cargo/workspace and DSR gates were not run. This owner does not resolve tank internal waves, mixing gradients,
neck inertia, heat transfer, volume motion or one shared tank across multiple
sides. Evolving-species Cartesian finite tanks remain a separate integration;
the source mixture plenum candidate API is already available.

## Full-gradient WALE and Vreman Newton integration

`PlanarSgs` supports Smagorinsky, WALE and Vreman through the existing
steady, BDF1 and unequal-step BDF2 Newton entry points, including prescribed
velocity and mixed velocity, Navier-slip, open-traction and reservoir
boundaries. `FluxViscosityGradientTangent` represents a cellwise viscosity
and its analytic derivative with respect to the full native `du_i/dx_j`
tensor. The response is prepared once per triangle per assembly and
contracted with the original BDM basis gradients. The same derivative owner
includes cell stress, harmonic SIP weights and penalties, viscosity-dependent
boundary loads, and both transported and advecting velocities. Existing
scalar rheology APIs and their arithmetic remain available.

WALE and Vreman use the shared `fs-lbm` analytic derivatives and original
primal closure values, retaining rotational dependence and authored filter
widths. No production finite differences, scalar-strain fit, viscosity
clipping, separate nonlinear solver or additional physical timestep is
introduced. Original source-validity checks remain authoritative through
final reaction assembly. All histories remain immutable. Singular-gradient
selections follow the shared kernels' generalized-derivative contract.
The `Estimated` planar scope does not imply a three-dimensional cascade,
turbulence qualification, new wall model, heat deposition, cancellation
context or global convergence guarantee.

All 24 selected planar/steady/Newton integration regressions and both
current-main source-refusal unit regressions passed in a focused
original-source native Rust build. Directional Jacobians matched independent
original primal assemblies under general and pure-rotation perturbations,
DG jumps, and reversing open flux across sixteen closure/boundary/convection
combinations, with frozen-derivative negative controls. Further checks cover
one response per triangle, analytic pressure/traction/virtual work, time
histories, exact disabled-SGS reduction, and active modeled-stress endpoints.
In the active-stress fixture, WALE converged in ten iterations with eight
backtracks to residual 8.6354e-14; Vreman converged in six iterations with
three backtracks to 1.13075e-13. These are fixture-specific numerical results,
not general runtime or turbulence-validation claims. The original numerical
kernels were compiled without numerical stubs; material-adapter integration
and full Cargo/workspace validation were not run.

## Local-temperature multicomponent miscibility

`miscibility::multicomponent::nonisothermal` binds independent accepted species
inventories to the same original multicomponent caloric law, one temperature
per cell, and explicit canonical face conductivity/heat-of-transport data.
The full Maxwell–Stefan friction solve retains the original molar frame,
species identities, complete pair diffusivities, and shared dense workspace.
No collection of independent binary fluxes replaces that coupled solve.

Original finite-face composition efforts exclude the source-derived thermal
chemical-potential contribution. Without authored Soret data, uniform
composition has zero molecular flux even with unequal temperatures, pure Cp,
and B/T interactions. Authored heat-of-transport contrasts enter both the
species force and reciprocal Dufour energy flux. Original divided partial
enthalpy contrasts are formed before large common references are rounded.
The source's sufficient whole-simplex convexity bound applies throughout the
complete temperature bracket, not only at the graph's reference temperature.

Every explicit stage evaluates all face rates at one accepted state, then
applies equal-and-opposite species and enthalpy transfers synchronously.
Temperature recovery uses original Cp integrals and changed inventories;
an open cell has no independent entropy gate. Whole-window component/cell
capacity, original total enthalpy, and total entropy checks precede publication.
The nonnegative initial face entropy rate is a constitutive diagnostic, not
a discrete-step dissipation proof. Positive support is never clipped or seeded.

Adaptive evolution compares one full coupled step with two half steps using
every component, its compensated low word, and each local temperature.
Unchanged initial face rates are reused for full/first-half trials and interval
retries. Only accepted fine-path species/energy transfers enter receipts;
completed rejected work is counted separately. A single physical allowance
applies to the final complete interval, including short terminal remainders.

Thermodynamic sampling remains on representable component high words. A
conservative original-caloric bound rejects compensated inventories whose
omitted energy exceeds the cell's share of the energy budget; transported
reference energy cannot silently become heat. This is an explicit numerical
range refusal, not a substitute two-word thermodynamic model. Fixed equal
partial molar volumes, fixed pressure, positive support and convex homogeneous
liquids remain required. Multicomponent interfaces, thermal expansion, variable
pressure/momentum feedback and dynamically selected material diffusivities
are not supplied by this numerical owner.

All 134 actual-source numerical miscibility tests passed in the focused Rust
workspace: 117 earlier cases, eleven new independent spatial/precision cases,
and six coupled temporal cases. Oracles include independent ternary friction,
variable-Cp/B enthalpy, reciprocal Soret/Dufour transport, heat-only analytical
relaxation, finite-step entropy refusal, trace support, cancellation/replay,
work/transfer receipts and a large-reference artificial-heating regression.
Seven original-source caloric transport cases separately passed through the
actual public API. These are numerical verification, not experimental liquid
validation or full workspace/DSR/material-build evidence.

### Transported-density flow with finite external RCR storage

`density_momentum_coupled_bdf1_windkessel_step` advances density/inventories, native momentum, external compliance and optional series inertance at one BDF1 endpoint.
Existing SI loads/history/result carriers and unit conversions are reused; native reduction freezes only immutable previous fluid/compliance histories, while each nonlinear trial resolves the current density-dependent pressure response.
All traction ports may carry loads. Reservoir laws supply pressure datums; caller-owned inflow density/inventory and Open/OpenReservoir momentum semantics remain authoritative under reversal.
Final storage, resistor loss, pressure/reaction work and original dynamic-law residuals use measured accepted port flow. The kinetic report derives from actual heterogeneous-density energy, never a uniform reference-density surrogate.
External energies/powers remain J/m and W/m, while the nested density-step inventories already include declared depth. Multiply external accounts by that depth for total J/W.
Physical resistor loss and the reported BDF1 numerical energy defect remain separate from transported caloric inventory; no external EOS, material reservoir inventory or thermal deposition is implied.
Original transport/hydraulic error owners and late cancellation refuse the complete fluid/storage result atomically. Two-step physical-oracle regressions cover heterogeneous inertia, both compliance histories, normalization, per-port work, storage-driven reverse flow and immutable refusal.

## Original-source local-temperature multicomponent material fields

`ResolvedMulticomponentMiscibilityField::bind_nonisothermal` requires the exact
original liquid-law instance retained by the supplied multicomponent caloric
law. It preserves canonical species/molar masses, the original interaction
model and pure Cp sources, fixed pressure, geometry, cell-count budget and
source-volume allowance. Mixture transport Cp is never substituted for the
original pure-liquid calorics. Face conductivity and heats of transport remain
explicit numerical data.

Fixed and adaptive molecular/heat candidates first complete the numerical
transaction, then qualify every endpoint through the original material law
at its actual local temperature and independent component amounts. Source
identity, composition axes and selected receipts remain observable. A late
source-domain, density/volume or cancellation refusal publishes no endpoint.
Final source qualification is not a claim that every rejected trial remains
inside a material-card domain.

Two added original-source regression groups exercise variable pure Cp plus
B/T enthalpy, actual-temperature transport properties, exact liquid/caloric
identity, pressure and geometry, unchanged numerical receipts, cancellation
and replay, and late original composition/volume refusal. They are formatted
and source-reviewed. The selected material metadata check was killed while
compiling the prerequisite asupersync crate, before these adapters/tests were
type-checked; no material test execution or full workspace validation is
claimed. The separate numerical lane has 134 passing miscibility tests and
seven passing original-source caloric transport cases.

## Evolving-species finite tanks on planar Cartesian gas boundaries

`MixtureEuler2d::with_finite_reservoirs` couples selected already-bound uniform,
exactly stationary reservoir sides to independent finite, rigid, well-mixed
species inventories. The explicit total workspace budget includes both prepared
inventory buffers and interval reports. Initial canonical partial densities,
pressure and per-species energy references remain owned by the original source.
The shared Cartesian `mixture::finite` storage calls the original
`MixtureGas` plenum's `exchange_candidate_into`; it introduces no duplicate
caloric model or dimension-specific inventory update.

The full forward tank stage drives the second native SSP RK2 residual. The
accepted tank receives the same outward canonical species and total-energy
face impulses, with native RK weights, as the field boundary account. Complete
incoming kinetic energy thermalizes. Flow pressure work is already present in
the total-energy flux and is not added a second time. The entire Cartesian
momentum flux is reported to the stationary support, which performs no work.
Reports retain before/after pressure, temperature, volume, species inventories,
referenced energy, per-species residuals and original-gauge energy residuals.
Accepted exterior getters return the evolving tank composition and pressure.

Initial and both residual CFL bounds include each tank's summed aperture-wave
rate `0.5 * sum(face_area * signal_speed) / volume`; large field cells cannot
hide small finite storage. Side transfer and wave sums are compensated. Warm
stages, retries, inventory restoration and publication allocate no storage.
Cells, tanks, time, boundary impulses and reports publish together after the
last cancellation poll. Failed trial inventories remain private. Prepared
reactive/split clones retain and restore accepted tanks on every retry; tank
chemistry itself is not inferred. Finite sides retain the compact law when a
later held-reservoir selection omits them. The explicit molecular-bath path
continues to refuse a compact tank without physical face transport data.

Focused regressions cover variable-Cp filling and blowdown with different
species mixtures on all four sides and both reconstruction orders, exact-zero
species retention, full shear/support transfer, global species and energy
balance, RK2 pressure-feedback refinement, tiny-tank initial CFL, final-poll
rollback followed by deterministic retry, capacity admission and prepared
clone/restoration. This remains an estimated compact-reservoir engineering
model: no internal tank waves, finite mixing time, heat exchange or discrete
entropy theorem is claimed. Broad gas-control-volume/network beads remain open.

## Conservative thermal stirring of multicomponent liquids

The local-temperature multicomponent field now supports prescribed-volume
`advect_candidate`, same-interval `transport_candidate`, and complete coupled
`advance_transport_candidate`. The original prepared canonical flow retains
both directional streams, closed cell-volume continuity and the outgoing
donor CFL bound. Each synchronous donor substep transports every independent
component, its compensated low word, and the original homogeneous parcel
enthalpy. Original changed-inventory calorics recover each receiving
temperature; temperature averaging and an additional mixing-heat source are
not substituted.

The combined kernel composes donor stirring with full Maxwell–Stefan and
reciprocal thermal transport over one physical interval. Full-versus-two-half
adaptation measures both composition and temperature after the complete
kernel, including splitting error. Each coarse trial uses at most one donor
CFL step so independent subcycling cannot make coarse and fine trajectories
identical and hide error. Suggested next intervals also respect this cap
after a short terminal remainder. Separate accepted advection, molecular and
combined species/energy receipts remain observable, as do completed rejected
work counts. Original inventories, enthalpy and total entropy are audited
again over the complete accepted window.

A conservative parcel-energy arithmetic-resolution estimate rejects large
reference energies whose finite heat carriage cannot meet the requested
allowance. It includes absolute pure and original pair contributions,
directed products and accumulation; internally cancelling B pairs cannot
hide a large estimate. Accepted adaptive donor estimates share one
whole-window budget. This is a numerical admission estimate, not an interval
or source-uncertainty certificate. The numerical workspace now explicitly
charges the fixed-size original caloric state retained for every donor cell,
including at small actual species counts.

All 145 actual-source numerical miscibility tests passed: the earlier 134,
seven new thermal donor/coupled cases and four complete-kernel adaptive
cases. Independent checks include three-cell synchronous circulation,
nonideal variable-composition enthalpy, analytical counterflow relaxation,
trace/low-word conservation, coupled-clock/transfer receipts, thermal-only
temporal refinement, CFL alias prevention, final-step bounds, atomic
cancellation, and explicit large-reference parcel refusal. The long
40-substep counterflow explicitly offers 1e-7 J arithmetic allowance while
retaining independent 2e-9 K and 1e-8 J endpoint checks. No mechanical-work,
two-way momentum, thermal expansion or multicomponent interfacial owner is
introduced by these prescribed incompressible-volume operators.

## Original material sources through thermal multicomponent stirring

The source-bound local-temperature multicomponent field also wraps
`advect_candidate`, `transport_candidate` and
`advance_transport_candidate`. All five thermal evolution paths retain the
same exact original liquid/caloric law references and share one endpoint
refresh. Every accepted cell is requalified at its actual temperature and
independent composition, with original pressure, material receipts, cell cap
and source-derived volume allowance. A late source failure leaves the
caller's accepted field unchanged.

Donor, coupled-step and coupled-adaptive candidates expose their numerical
receipts unchanged, including original parcel energy, all species transfers,
energy/entropy audits, accepted work/error metrics and parcel-resolution
estimates. Source qualification covers the complete accepted endpoint; it
does not imply dynamic material-derived face coefficients or qualification
of every internal/rejected trial.

The two original-source thermal regression groups now also compare raw and
source-bound endpoints and receipts for all three stirring paths, independently
check variable-Cp plus B/T total enthalpy, refresh actual-temperature material
properties, and exercise late original-domain/volume refusal after a
numerically successful coupled step. These extensions are formatted and
source-reviewed; material compilation/test execution remains blocked by the
asupersync prerequisite compiler SIGKILL. The separate actual-source
numerical suite passes all 145 tests, including the eleven thermal stirring
and coupled-adaptation cases.


## Transported-density endpoint with source molecular rheology

`FluxSystem::density_momentum_coupled_bdf1_rheology_step` extends the existing
endpoint-advector BDF1 owner with a fallible **cellwise** dynamic-viscosity source.
Every native Picard assembly re-transports the same immutable old density and
specific inventories, evaluates each original triangle's source at its centroid
with that candidate's transported state and native `sqrt(2 D:D)` shear rate, and
uses the retained coefficient at every volume/face stress quadrature site for
that cell. The existing spatial-viscosity, integral-flow-port and RCR entry
signatures retain their original behavior through one private shared selector.

The original native Picard owner requires the unrelaxed velocity correction and
both reassembled endpoint momentum residuals to converge. Accepted density,
compatible mass flux, source samples, velocity, physical pressure and reactions
therefore belong to one endpoint. An original `MaterialError` is preserved in
`DensityMomentumError::Material { cell, error }`; unsupported sources, invalid
viscosity, transport/CFL failure and cancellation publish no partial state.

The result retains molecular viscosity, objective shear rate and `mu*gamma²`
[W/m³] per cell, plus its compensated volume integral at the authored depth [W].
Reporting never queries the source again. Molecular heating remains separate
from temporal/upwind/SIP numerical losses and is not deposited into a carried
inventory. Rigid rotation contributes no physical shear. Native material laws
such as Carreau–Yasuda can be used directly; this entry invents no constitutive
law, temperature feedback, yield regularization or compressible density model.
It retains the original first-order endpoint-advector/CFL limits. One reusable
source buffer is checked against the per-owner storage allowance; that allowance
is not a native linear-solver or total-process memory cap.



## Source rheology in transported-density flow ports and finite RCR storage

`density_momentum_coupled_bdf1_rheology_flow_control_step` resolves integral port
pressures with the same accepted transported density and molecular shear law.
`density_momentum_coupled_bdf1_rheology_windkessel_step` couples that endpoint to
the original finite compliance, resistance and series-inertance histories.
Both use the existing private viscosity selector and native Picard/pressure
responses. The original RCR reduction, original-equation reporting, physical
pressure conversion and one molecular-heating finalizer are shared; no second
nonlinear solver, constitutive law or pressure response is introduced.

Each trial rebuilds its pressure response with its actual transported density
and retained cellwise molecular viscosity. Accepted source samples, flow,
pressure, boundary reactions and load histories belong to the same endpoint.
Source failures remain original `MaterialError` values inside the momentum
error, including through RCR reporting and SI conversion. Rejected trials
advance neither the old fluid inventories nor the external storage.

The original traction-datum, port-disjointness, reversal/donor, BDF1, CFL,
cancellation and storage restrictions remain. Molecular heating is reported in
whole-domain W using the authored depth; external RCR powers remain W/m. Neither
physical resistor loss nor numerical BDF defect is combined with molecular
heating or silently deposited in a carried inventory.


## Finite evolving-mixture pressure storage on volumetric Cartesian fields

`MixtureEuler3d::with_finite_reservoirs` connects each selected stationary,
uniform physical side to the same source-owned finite-mixture storage used by
the planar field. All faces on one side share its rigid, perfectly mixed tank;
different sides retain independent inventories. Both native SSPRK2 residuals
use the corresponding field and tank stages. The complete first tank stage
drives the second pressure/composition trace, and the accepted inventory
receives the same two weighted outward face transfers as the field account.

Every canonical species and the original referenced total energy is conserved
across a port. Incoming kinetic energy thermalizes, pressure flow work is
counted once, and complete Cartesian impulse goes to the fixed support. Public
uniform and physical-face readback follows accepted tank composition and
pressure. Reports preserve source pressure/temperature, species inventories,
support impulse and both referenced and original-gauge energy residuals.

Initial and residual admission include compensated finite aperture/volume wave
rates alongside the unsplit three-axis field CFL. Prepared accepted/trial
inventories, explicit total capacity admission and joint publication preserve
cells, tanks, loads, clock and reports on cancellation or refused trials. The
existing forcing and timed-profile routes share this stage owner; other sides
may retain their original prescribed profiles. An explicitly requested contact
equilibrium correction remains a separate field energy receipt, not tank heat.

Reactive and molecular prepared-field cloning/budget/reset sites retain the
finite inventory and its accepted reports without warm allocations. This adds
no reaction or heat-transfer law inside the tank. The existing molecular
admission still refuses compact ports without a physical diffusive bath, and
nonzero moving-grid admission still refuses bound exteriors.

Five new volumetric native regression groups passed together with the five
planar groups: all six orientations, both reconstruction orders, filling and
blowdown with different variable-Cp compositions, exact-zero species, complete
support/shear transfer, body-force work, global species/energy balance, tiny-tank
initial CFL, workspace refusal and late cancellation followed by exact retry.
The measured pressure-feedback errors decreased by approximately four on
successive timestep halvings. The bounded build used the actual inviscid field,
reconstruction, RK, forcing, moving-grid/profile interfaces and original gas
EOS/port/Plenum sources. Full reactive/molecular orchestration and the complete
workspace were not executed; their narrow clone/reset edits were reviewed and
Rust-parsed. These compact tanks retain the original estimated-model scope:
no internal tank waves, finite mixing time or discrete entropy theorem.


## Coherent cell material forces in native laminar BDF1 flow

`MaterialFluxSystem::bdf1_navier_stokes_mixed_rheology_with_cell_forces`
resolves an original material and additional physical force density [N/m³]
together at each actual native objective-shear state. The original
`RheologyField` cache prepares all cell responses before the original volume
quadrature; no coordinate search, constitutive clone or fallback source is
published. The returned source force is constant on its original triangle.
A separately supplied held spatial force is added once, and the same retained
reference density scales inertia, both forces, viscosity, traction and pressure.

Cell identity now reaches the private original load quadrature directly.
Original held-force entries adapt their existing callback without changing
their public API or physical treatment. The existing mixed BDF1 owner still
admits histories/boundaries, iterates convection and rheology, and independently
reassembles its accepted nonlinear residual. That final reassembly prepares
both the accepted viscosity AND the accepted force; constrained boundary
reactions retain the complete original volume load. The returned tuple contains
the original physical flow/material endpoint and exact accepted cell forces in
triangle order. No force work is automatically deposited as material heat.

This entry adds no density evolution, EOS, material model, thermal law or new
timestep scheme. Original source, changed-density, nonfinite-force, boundary,
history and nonlinear refusals return no replacement history. Three authored
native regression groups cover independent hydrostatic cell pressures and
constraint resultants, nonzero cavity circulation against the original spatial
load entry with original Carreau rheology, and late source/nonfinite-load refusal.

Bounded source-linked compilation used genuine native BDM/material sources and
original numerical dependencies. It reached this new entry and identified two
concrete API/type errors, both corrected to the original error variants. A
complete native artifact or test execution was not obtained because unrelated
local module dependencies were missing/stale; these new tests remain unexecuted.
The executor then became unavailable; final formatter checks and corrected
native integration execution remain pending.

## Passive scalar and constant-Cp thermal transport with planar SGS

`FluxSystem::prepare_scalar_transport` prepares a conservative P1
discontinuous-Galerkin scalar equation on the original BDM mesh. Its volume
velocity and sign-split upwind faces retain both directions of an affine
reversing trace. Symmetric interior-penalty diffusion uses actual triangle
gradients and harmonic diffusivity/normal-height weights, without assuming
orthogonal centroid connectors. Both interface weights are formed directly
so a small finite contribution survives large diffusivity contrast.

`PlanarSgs::prepare_scalar_transport` binds the same supplied velocity to the
existing Smagorinsky, WALE or Vreman field and adds `nu_t/ratio` to explicitly
authored molecular scalar diffusivity. Temperature callers supply turbulent
Prandtl numbers; passive-species callers supply Schmidt numbers. The returned
native transport advances the scalar with BDF1 or the original unequal-step
BDF2 history owner, including exact BDF1 startup and the existing step-growth
bound. Prepared velocity, diffusivity and boundary traces belong to the
supplied endpoint and must be refreshed when they change. Sparse CSR,
controlled ILU(0), FGMRES and independently checked original/equilibrated
residuals remain the numerical owners.

Value boundaries use consistent weak diffusion and explicit incoming values;
flux boundaries require a donor wherever any incoming stream exists,
including zero-net counterflow. Each accepted result reports separate outward
advective/diffusive numerical boundary rates, the actual BDF storage
derivative, source rate and unrepaired global balance defect. Rates are per
out-of-plane depth; constant-density, constant-Cp temperature budgets become
energy after multiplication by rho Cp and physical depth. Neither molecular
nor modeled momentum dissipation is silently deposited as heat.

Five focused native integration groups and one private high-contrast
physical-interface regression passed. Coverage includes all-model accelerated
thermal-pulse mixing versus molecular-only transport, no spurious heating of
a uniform field, exact unequal-step histories, skew-mesh diffusion refinement,
closed inventory, independent reversing-face flux integrals, material
contrast, disabled-SGS bitwise identity, cancellation and immutable refusal.
The interface oracle retains both one-sided contributions at contrast 1e18.
These focused builds used the original numerical kernels without numerical
stubs; full Cargo/runtime integration was not run. This is passive transport
with fixed storage properties, without variable-Cp enthalpy, phase change,
or three-dimensional LES qualification. The opt-in joint Boussinesq owner
below supplies thermal feedback under a separate linear-density scope. The original
BDF entries do not enforce scalar bounds; the explicit bounded lane is described
below.

## Coupled planar SGS temperature and Boussinesq momentum

`PlanarSgs::bdf1_boussinesq_step` couples the original mixed BDM1/P0 momentum
owner to its P1-DG thermal transport. Each original Picard candidate recomputes
Smagorinsky, WALE or Vreman stress and `alpha_m + nu_t/Pr_t`, solves thermal
BDF1 from the same immutable old temperature, then applies the resulting
affine temperature at the original cell-indexed momentum quadrature. Final
native residual and reaction assembly repeats both preparations at the
accepted velocity; the returned temperature is that same endpoint. This is
implicit partitioned BDF1 with the original sparse solves and work caps,
without a second coupling loop or a monolithic Newton claim.

The acceleration is `-beta*(T-Tref)*g` for explicitly signed physical gravity.
Pressure and prescribed total traction use reference-hydrostatic-reduced
pressure `p_physical/rho - g.(x-x_ref)`; the caller retains its boundary datum
and must not add reference-density gravity again. Molecular viscosity,
thermal diffusivity, density and Cp stay fixed. The authored signed expansion
coefficient supplies only the linear momentum anomaly. Positive kelvin and
the caller's `|beta*(T-Tref)|` limit in `[0,0.1]` are required throughout every
candidate and the actual affine boundary-temperature traces, including their
endpoints. This declared small-departure scope is not an approximation-error
certificate or an EOS model.

Results retain independent original momentum and thermal residuals, thermal
storage/source/boundary rates, final SGS diagnostics and the accepted integral
of buoyancy acceleration dotted with velocity. The latter has units m4/s3 per
unit depth; multiplying by constant rho and depth gives watts. It is neither
a heat source nor a whole-step mechanical/potential-energy balance. No
molecular, modeled SGS, pressure or buoyancy work is deposited as heat.
Picard convergence remains conditional on coupling strength and timestep.
Thermal BDF1 is not bound-preserving: invalid temperature/departure or either
solver's refusal returns no joint endpoint. Cancellation at existing assembly
and thermal checkpoints preserves both input histories; the synchronous
momentum linear solve has no new preemption guarantee.

Four new native physical groups passed: all three closures with heat-driven
circulation and independent final-temperature momentum/thermal oracles;
uniform source heating with an exact endpoint force and thermal budget;
both gravity signs with the native pinned hydrostatic pressure balance; and
nonconvergence, invalid source/temperature/affine-boundary traces, late
cancellation and deterministic retry. All five existing scalar groups also
passed. The focused build compiled the original numerical NS, boundary,
reaction, SGS, scalar, LU and controlled-Krylov owners without numerical
stubs; full Cargo/runtime integration was not run.

## Resolved moving liquid heat around a phase-changing vapor bubble

`cavitation::vapor::liquid::SphericalLiquidField` retains up to 64 spherical
liquid shells in fixed-volume coordinates `r_i^3-R^3`. Arbitrary positive
face grading can resolve the interface diffusion layer. The original source
liquid enthalpy supplies every shell's finite storage and all signed ALE
mass transport. The native vapor mass transfer gives the same relative
outward mass flux, `-dm_v/dt`, through every face. Original source
conductivity and the existing spherical resistance geometry supply
interfacial, interior and outer conduction.

`cavitation::vapor::conducting::ConductingVaporBubble` advances this
retained field, native phase mass/temperature and radial motion on the
existing fixed-stack adaptive DP5 clock. Every stage computes interface
Fourier heat from that stage's liquid profile, then evaluates the original
phase balance and reciprocal liquid transport. Liquid temperature increments
are the adaptive coordinates within each physical window; an absolute Kelvin
datum cannot dilute their relative error scale. Shared phase derivatives and
final accounting also remain the owners for the prescribed-heat entry.

The same fixed enthalpy reference in both phases makes their interfacial
mass enthalpy cancel. Independent original-source endpoint changes and
quadratures admit the liquid first law and the combined thermal account
`Delta(U_v-h_ref*m_v)+Delta H_liquid+W_v+Q_outer+H_outer-Q_visc=0`.
Outer conductive heat and advected enthalpy are positive out of the finite
grid. Native radial viscous dissipation is deposited through its exact
spherical integral in each represented shell; the exterior share is
reported separately. Phase pressure work remains in its mechanical owner,
so this thermal account is not a claim to integrate all liquid kinetic and
surface energies.

The original matched saturation-path liquid calorics remain an explicitly
bounded incompressible approximation. Actual source density, heat capacity,
conductivity, viscosity and sound-speed drift, saturation-path capacity
correction, Bernoulli pressure-density departure and estimated omitted
reversible pressure heating are checked. Interior advection is donor upwind.
Exact boundary enthalpy traces require the authored outflow half-cell Peclet
bound; unsupported coarse outflow is refused. Finite outer temperature and
enthalpy traces connect to an explicit bath, not a sealed liquid vessel.
Local phase equilibrium, uniform pure vapor, fixed bulk-liquid density and
the original radial/Mach limitations remain. There is no finite-rate
nucleation, noncondensable gas, shock or wall-jet model in this owner.

A complete proposal retains all liquid cells with the phase/radial history.
Cancellation, physical admission, work limits, independent balance refusal
or a dropped proposal leave every accepted state unchanged. Spatial
resolution, exterior-domain extent and temporal refinement still need
physical consideration; small balance defects alone are insufficient.

Focused builds of the actual material, spatial, phase and stepper sources
passed all 99 native cavitation tests (87 existing, six liquid-operator
groups and six joint thermal groups). The new checks include spherical
Laplace heat transport, both signs of conservative phase exchange, exact
resolved/exterior viscous heating, independent cumulative original-source
mass/energy balances, complete rollback and source/approximation refusal.
The joint phase/temperature/heat response also converges under 16/32/64-cell
and 80/40/20 ns refinement. This is focused numerical verification; the
complete workspace and experimental validation were not run.

## Multicomponent spatial demixing with explicit interface energy

The experimental miscibility module now has a separate
`multicomponent::phase_field::MulticomponentPhaseFieldNetwork`
owner for nonconvex original regular-solution liquids. It retains the original
full-species thermodynamic model, complete positive pair diffusivities and
coupled Maxwell–Stefan molar-frame solver. Every original component remains an
independent compensated inventory. Common fixed partial molar volume preserves
rigid cell capacity; unequal molar masses do not make this a barycentric frame.

A finite positive constant energetic coefficient kappa [J/m] defines each edge
energy as kappa A/(4 d) times the sum over all species of squared mole-fraction
jumps. This species-symmetric normalization reduces exactly to the existing
binary kappa A/(2 d) convention. The cell chemical gradient is the derivative of
that same discrete energy, divided by the actual admitted cell molar inventory.
One canonical edge pass assembles all gradient potentials in O(edges * species)
work before the shared original force/friction sweep. The sufficient bulk
convexity bound remains diagnostic for this owner; public homogeneous network
construction still requires it. No raw nonconvex homogeneous network or field
escapes through the public phase-field API.

The synchronous explicit candidate applies all paired component transfers
before checking positivity. Publication requires independent species/capacity
balances and the actual change in bulk Gibbs plus interfacial energy within one
offered energy allowance. Bulk and interface increments are reported separately
and evaluated through stable differences. A valid smoothing step can increase
bulk Gibbs while decreasing total free energy. Initial chemical dissipation is
a rate diagnostic, not a substitute for the finite-step energy audit.

Common strictly positive species support remains required. No component is
seeded, clipped, reconstructed from a carrier or replaced by equilibrium phase
fractions. The scalar energetic gradient model does not assert independently
calibrated pair surface tensions, variable-kappa physics, nucleation, spatial
convergence, unconditional time stability or hydrodynamic capillary work.
Prescribed-flow advection is not exposed by this owner. Checked additional
workspace covers interfacial geometry and simultaneous potential/candidate
storage; canonical loops retain bounded cancellation checks.

Focused G0/G3 regression sources include discrete spinodal dispersion,
short-wave damping with increasing bulk Gibbs, independent unequal-pair ternary
friction and dissipation equations, binary normalization, component permutation,
all-species conservation through the supported species limit, scoped nonconvex
admission, workspace limits and transactional refusal. These numerical oracles
do not constitute experimental calibration or a G2 physical benchmark.


## Bounded density stages inside the coupled momentum endpoint

All six coupled density/momentum, molecular-rheology, integral-flow and RCR
entries accept `Into<CoupledDensityMomentumControls>`. Existing ordinary
`DensityMomentumControls` retain strict single-Euler-step behavior and macro-CFL
refusal. Passing the existing `DensityMomentumSubcycleControls` or
`DensityMomentumSspRk2Controls` explicitly selects their bounded dyadic density
work; no additional public solver-entry family is introduced.

Every Picard trial starts again from the SAME immutable old density, inventories
and momentum history. Its candidate endpoint velocity is held fixed only within
that trial's density stages. Native Euler subcycling or SSPRK2 retains the actual
accumulated/RK-weighted donor and face mass flux. That original flux and the
macro old/new density metrics enter the same native macro momentum assembly;
flux is never reconstructed from final density. The accepted `density_work`
variant retains its original native count, duration, CFL and complete transport
report, identical to the physical inventory report in `step`.

These work counts describe the accepted transport assembly, not discarded
nonlinear trials. Each assembly enforces the original substep and per-owner
storage allowances; the existing finite nonlinear budget bounds repeated trial
work. A changing dyadic count can still exhaust Picard convergence. Original
unrelaxed-correction and reassembled endpoint-residual checks remain decisive;
failure or cancellation returns no partial fluid or hydraulic history.

Momentum, endpoint-advector coupling and external RCR histories remain macro
BDF1. SSPRK2 has second-order integration only for the fixed-advector density
subsystem. Rejected iterations and density stages do not advance extra physical
time or repeatedly debit inventories. Molecular heat, native numerical kinetic
losses, external resistor losses and RCR BDF defects remain separate.

Four new native integration regressions cover macro CFL above one under
donor-cell, bounded-linear density and bounded-linear inventory reconstruction;
exact independent replay of accepted native transport; nonuniform-density,
nonzero-shear macro kinetic-work closure; pressure/rheology propagation; and
original macro RCR storage equations with atomic density-work refusal. All
34 focused tests pass: these four, 14 existing coupled/rheology/flow/RCR tests,
10 original frozen Euler-subcycle/SSPRK2 consumers and six original single-step
density/momentum consumers. The run uses the original native dense assembly,
transport, Picard, pressure/RCR and LU/GEMM owners. It is not full-workspace or
sparse-route validation and makes no performance benchmark claim.

## Adaptive multicomponent demixing

`MulticomponentPhaseField::advance_candidate` compares one explicit interval
with two half intervals and accepts only the fine trajectory. All independent
component high/low inventories control composition accuracy. A separate
whole-graph L1 norm of individual edge-energy differences controls interface
accuracy; errors at different edges cannot cancel in a total-energy difference.
The relative interface scale sums the larger coarse/fine energy at each edge.

The original complete Maxwell–Stefan rates at an unchanged initial state are
shared between the coarse and first-half trial and retained across retries.
Changed midpoints evaluate fresh original rates. Accepted face receipts exclude
all coarse/rejected transfers, while completed-step and constitutive-work
counters include them. Explicit attempt and step-size bounds limit work, and a
final sub-minimum remainder must satisfy the same acceptance rules.

The complete requested interval independently rechecks every species, rigid
capacity and total bulk-plus-interface energy using the original allowances.
Failure or cancellation publishes no partial field or elapsed clock. Adaptation
provides local temporal estimates; this explicit fourth-order evolution can
require small steps on fine meshes and does not claim unconditional stability
or a global trajectory error certificate.

## Original liquid sources through multicomponent phase separation

The phase-field material adapter retains the exact borrowed multicomponent
liquid law and resolves every returned cell at its actual original species
inventory and fixed temperature. It shares model/species/molar-mass admission
and cell source/volume validation with homogeneous miscibility. Selected claims,
pressure and the full original composition chart remain authoritative; cells
stay on the original homogeneous free-energy branch even where it is unstable.

Single-step and adaptive endpoint wrappers preserve the numerical demixing
receipts and qualify source density against each rigid cell volume before
publishing. A source-domain, material-volume, budget or cancellation refusal
leaves the original field and retained bindings unchanged. The wrapper exposes
only the interfacial numerical owner. Material-cell storage has its own explicit
count budget, separate from numerical workspace.

These source receipts qualify completed endpoints, not each internal/rejected
adaptive trial. Full pair diffusivities and kappa remain explicitly authored
numerical inputs; endpoint source refresh does not imply dynamically sourced
transport coefficients or two-way variable-density momentum coupling.

## Local-temperature multicomponent phase separation

The experimental multicomponent phase-field owner can bind the original
multicomponent liquid caloric law and evolve all independent species together
with local cell temperatures. The original A+B/T mixing law and pure-liquid
Cp curves supply source-conjugate chemical efforts, divided partial-enthalpy
carriage and bounded caloric inversion. Authored heat-of-transport contrasts
remain reciprocal in species and energy transport; no Soret coefficient is
inferred from the interaction law.

Constant energetic interfaces retain exactly the existing species-symmetric
edge law E_e=kappa*A/(4*d)*sum_s(delta x_s)^2 and carry no interface entropy.
Each cell stores U_i=H_i+half the sum of its incident edge energies. The
interface entropy-conjugate species potential is
sum_j K_ij*(1/T_i+1/T_j)*(x_i,s-x_j,s)/N_i, with K_ij=kappa*A/(4*d).
Actual admitted published molar totals supply N_i. A canonical graph pass
assembles these potentials before the original coupled Maxwell--Stefan sweep.
The inverse-temperature weight uses (1+T_min/T_max)/T_min rather than a
potentially overflowing reciprocal sum or temperature product.

All species and paired total-local-energy transfers are simultaneous. After
the independent inventories determine their original fractions, each changed
edge deducts half its energy increment from BOTH endpoint bulk-enthalpy targets.
The same original caloric inverse then recovers each affected temperature,
including spectators whose published composition did not change. Bulk enthalpy,
interface energy and total-energy residual remain separate receipt quantities.

The original thermal inventory and caloric-resolution checks are shared, not
reimplemented. Publication requires complete all-species/cell-capacity closure,
actual homogeneous enthalpy PLUS interface energy conservation, and nondecreasing
original homogeneous entropy within the offered numerical allowances. Positive
initial constitutive entropy production does not replace that finite-step audit.

Only the private positive-interface owner may retain the original caloric
closure's explicit nonconvex-bracket result. Both original endpoint domains and
curvature evaluations still run; source, stability, arithmetic and cancellation
failures propagate. Public homogeneous thermal binding keeps its full-bracket
convexity requirement. Neither a raw homogeneous field nor a nonconvex
zero-gradient network escapes from the new public owner.

The model retains rigid geometry, equal fixed partial molar volumes, static
pressure, common positive support and constant energetic kappa. Its explicit
time update has no unconditional stability claim. Phase-field advection and
hydrodynamic capillary work remain separate unimplemented couplings.

Focused regression sources compare variable-Cp/B/T spinodal heating, unchanged-
composition spectator heating and cooling, unequal-temperature/unequal-pair
friction with reciprocal heat, an independent total-entropy directional
derivative, binary normalization, and species/energy conservation through the
supported species limit. Additional cases retain ordinary convexity admission,
workspace boundaries and transactional cancellation. These are numerical
oracles, not material calibration or spatial-convergence certificates.

## Conservative bounds for planar SGS scalar transport

`PreparedScalarTransport::prepare_bounded_transport` adds an opt-in invariant
interval for passive concentration or constant-Cp temperature. It uses the
original P1-DG advection/SIP matrix, original BDM velocity and the same
molecular plus authored SGS diffusivity. Symmetric algebraic graph diffusion
and original lumped P1 storage form a monotone Euler predictor. Limited paired
antidiffusive fluxes restore consistent mass and high-order transport where
the available scalar inventory permits. Equal-and-opposite interior corrections
conserve inventory; no scalar clipping, redistribution or mass repair is applied.

Complete signed weak boundary coefficients are retained on the original
two-Gauss traces, including SIP and reversing-flow interpolation. The discrete
constant-mode identity is explicitly admitted before claiming bounds.
Reservoir admission includes the affine trace's physical endpoint extrema.
The bounded SSPRK2 driver subcycles the complete interval under the actual
low-order advection/diffusion CFL and an explicit work cap. Both Euler stages
and the returned P1 field retain authored bounds up to the declared
floating-point allowance. Failure or cancellation returns no partial interval
and leaves the supplied history unchanged.

Accepted boundary accounts integrate the original advective and SIP diffusive
fluxes on those same stages, plus any boundary graph diffusion not removed by
limiting. This extra numerical exchange is reported separately from molecular
or SGS mixing. Source-free means no volumetric source; value and advective
reservoir boundaries may still supply scalar or heat. Homogeneous diffusive
flux is supported. Arbitrary volumetric sources and nonzero imposed flux need
a separately admitted bound policy and remain on the original BDF path.
Frozen coefficients/data, constant storage properties and the planar Estimated
scope remain unchanged. Individual scalar calls do not enforce a coupled
species sum; the joint-species entry below does. No new three-dimensional
turbulence claims are made. Limiting may reduce local order near extrema.

Ten focused integration tests and both scalar numerical unit tests passed
using the original numerical source build. Five new integration groups
demonstrate a sharp-front failure of unrestricted DG and the inventory error
of clipping, bounded propagation, skew-grid thermal conservation, all-model
and disabled-SGS coupling, reversing reservoirs, smooth-grid refinement,
explicit source refusal and atomic cancellation. A new private oracle on
obtuse triangles verifies signed SIP reservoir terms: inactive limiting
matches original consistent-mass LU evolution, while bounded heating closes
the budget with nonzero limiter boundary exchange. The prior six scalar
checks also remain green. Full Cargo/runtime validation was not run.


### Reacting gas on translating and deforming three-dimensional volumes

`spatial::mixture::reactive::ReactingEuler3d::step_moving` composes the original
prepared cell chemistry as C(dt/2), conservative ALE(dt), C(dt/2). The chemistry
halves use the actual old and new cell volumes, retain all three momenta and
physical chemical-plus-kinetic energy, and report their original species and
referenced-energy source changes. `step_moving_with_dissipation` surrounds this
composition with the original molecular D(dt/2) intervals on the same old/new
geometry. Maxwell--Stefan species transport, partial enthalpies, Fourier heat,
full Newtonian stress and normal moving-wall work retain their original owners
and prepared numerical buffers.

Closed reflecting and periodic domains admit translation and the existing
single-axis deformation. All phases share one actual accepted interval; an ALE
shortening restarts the complete source/transport composition. The final state
must remain inside the chemistry temperature interval and the mesh-relative
acoustic bound after the last source, including inviscid heat release. No
per-trial whole-field allocation or alternative reaction/transport law is added.

`MovingReactiveStepReport3d` retains the complete original reacting report and
separate six-side molecular mechanical-work and wall-heat integrals. Species
and referenced-energy residuals subtract the actual chemical source, while an
independent physical-energy residual restores the species energy references
using conservative differences. Pressure work remains in the central ALE
report. Chemistry contributes no boundary species flux or external force.

Accepted cells, geometry, physical time, source binding, reports and complete
face transfers publish together after the final cancellation checkpoint. A
refused trial leaves the accepted owners unchanged; reusable trial geometry is
reset from the accepted field before either moving or stationary stepping.
Open ports, retained SGS histories, finite tangential wall inertia and
simultaneous deformation of multiple axes are not composed by these entries.
The selected reactor method retains its own accuracy; symmetric ordering does
not establish a fully discrete entropy law or a splitting convergence order.

Focused regressions cover stationary native-path parity, translating exothermic
flow against the original homogeneous reactor, whole-interval acoustic replay,
all three deformation axes with stress/heat/source balance, energy-reference
invariance, source rebinding and early/late cancellation. An independent
bimolecular logistic solution checks concentration-dependent conversion on both
actual volumes under compression and expansion on every axis.


### Finite guided pistons driven by reacting spatial mixtures

`spatial::mixture::reactive::ReactingEuler3d::step_with_inertial_walls`
composes the original cell reactor and conservative ALE field with the shared
simultaneous two-piston mechanical solve. Every held midpoint-speed proposal
restarts the same accepted gas and executes C(dt/2), ALE(dt), C(dt/2), with
each chemistry half using its actual old or new physical cell volume. The
reactor fixes all three cell momenta; chemical conversion changes the original
calorics and pressure that subsequently drive finite piston acceleration.
There is no separate combustion force, duplicate heat-release source, or new
mechanical solver.

`step_with_dissipation_and_inertial_walls` adds the original molecular half
intervals around that reacting split. Maxwell--Stefan diffusion, partial
enthalpies, Fourier heat and full compatible Newtonian stress retain their
original source owners. Actual pressure and viscous face impulses drive both
finite masses with their conjugate mechanical work. Direct Fourier/Navier heat
remains external thermal exchange and is not inferred by subtracting large
mechanical-work values from total energy.

`ReactiveInertialWallsStep` reuses the existing generic guided-piston receipt.
Fixed transverse walls can absorb axial viscous momentum; that support impulse
is explicitly removed from the combined axial mechanical residual. Translation
and orientation guides carry the actual transverse and angular impulses from
all accepted end-plane faces and do no work. Both the actual chemical-plus-
kinetic gas energy and its stored referenced-energy equation must pass the
authored mechanical tolerance. Physical energy uses differences of species
mass before multiplying their individual source references. The independent
referenced gate subtracts `chemistry.referenced_energy_change_j` exactly once.

A shorter ALE interval or a source/mechanical refusal restarts the complete
physical interval. The prepared reactor, field and molecular buffers are
reused across proposals. Accepted cells, geometry, time, combined reports and
loads, external molecular binding/report/loads, and both finite piston states
publish only after the final cancellation checkpoint. The molecular variant's
explicit byte budget covers both complete reacting and molecular owners.

The selected axis has two reflecting ends; other directions are periodic or
stationary reflecting supports. Authored tangential wall velocities are zero.
Original moving-source restrictions remain, including refusal of unsupported
exterior/finite reservoirs, SGS histories and finite tangential wall inertia.
This is a bounded engineering composition with the original chemistry and
shock-capturing limitations. It does not establish general FSI stability,
global entropy stability, resolved flame/detonation thickness, or experimental
combustion accuracy. Mechanical energy tolerance must accommodate the authored
cell-reactor accuracy over the whole accepted interval.

Six focused native regressions pass on the default thread stack. They cover
reaction-driven piston acceleration against thermoneutral gas; compression
and expansion on all three axes; identical physical dynamics under unequal
species reference-chart changes with large formation energies and nonzero
conversion; molecular heat, side supports and guide loads; cancellation and
root-exhaustion rollback of both accepted owners; and whole-interval replay
after reaction-induced wave shortening. These use manufactured mechanisms
and original source solvers; they establish numerical composition behavior,
not experimental reaction-rate or combustion validation.

## Adaptive local-temperature multicomponent phase separation

The thermal phase-field owner advances a complete requested interval with
explicit time and work limits. One full interval is compared with two half
intervals of the SAME coupled species, total-local-energy and caloric-recovery
kernel. Independent original all-species inventory and Kelvin norms are joined
by a graph L1 norm of per-edge interface-energy differences. Opposing changes
on different edges cannot cancel in the temporal error estimate. Relative
interface scaling sums max(coarse, fine) energy over each physical edge.

The unchanged initial coupled face rates are reused by the coarse and first
half stages and across rejected retries. Midpoint rates are recomputed, and
acceptance invalidates the cache. The original shared thermal admission,
caloric-bracket retry policy, precision failures and work counters retain their
semantics. Positivity, bounded caloric inversion, energy and entropy refusals
can request smaller intervals; cancellation, original-source and unresolved
inventory precision failures remain transactional refusals.

Only the two accepted half-step paths contribute physical species and total
local-energy transfers. Completed steps, constitutive evaluations and caloric
iterations separately report numerical work, including completed rejected or
coarse work. A final complete-window audit rechecks all independent inventories,
cell capacities, original caloric precision, total homogeneous enthalpy plus
interface energy, and original entropy before publishing the endpoint.

The controller exposes local time-discretization estimates, not an
unconditional stability, global trajectory-error or spatial-convergence
certificate. Constant energetic kappa and the fixed-geometry transport model
remain the same as the fixed-step owner.

## Original-source thermal multicomponent interfaces

The original-liquid phase-field adapter can bind the exact original liquid's
caloric law and return a source-qualified local-temperature phase field.
Identity is checked against the retained liquid instance, not merely equal
parameter values. The original pure-liquid Cp curves and A+B/T mixture
calorics remain the numerical chemical, enthalpy and entropy authority.

Fixed and adaptive candidates preserve the entire numerical interface/thermal
receipt. Each completed endpoint resolves every cell through the original
liquid law at its actual local temperature and all independent component
amounts. Original composition axes, pressure, selected material claims,
source-volume allowance and cell budget persist. The shared existing thermal
source-resolution helper handles cell budgets, original-source errors,
density/geometry consistency and cancellation for both homogeneous and
interfacial owners.

Source, volume or late cancellation refusal returns no partially refreshed
endpoint and leaves the caller's field and receipts unchanged. A spinodal
cell remains an original homogeneous branch of the evolving diffuse
interface; it is not independently flashed to another material state.
Resolved mixture transport Cp does not replace original caloric integration.
Molecular pair diffusivities, face conductivity, authored heat-of-transport
contrasts and energetic kappa remain explicit transport/interface inputs.

Material receipts qualify completed endpoints, not every rejected adaptive
trial. The wrapper exposes only the full energetic thermal field, and retains
the exact borrowed caloric/liquid pair across repeated fixed and adaptive
advances. Focused integration cases cover identity refusal, actual-temperature
property resolution, independent total enthalpy-plus-interface accounting,
unchanged numerical receipts, and transactional original-source/volume gates.

## Conserved gas inside the resolved thermal vapor bubble

The existing `vapor::conducting::ConductingVaporBubble` accepts an optional
original `ConservedBubbleGas` through `with_noncondensable_gas`, under
`caloric-cavitation`. The pure-vapor constructor remains available. Gas
mass and dry composition are immutable; the bubble's existing temperature
and radius are the sole evolving contents state. There is no separate gas
clock, fitted polytropic exponent or second radial/thermal integrator.

One shared loaded phase kernel evaluates water partial saturation pressure
and original ideal dry-gas partial pressure. Its common-temperature balance
is `[V*C_phase+m_g*Cv_g] Tdot + [rho_v*L+p_g] Vdot = Q_liquid`.
The original Stefan relation then supplies interface velocity and vapor mass
flux. Total pressure drives liquid momentum and contents pressure work.
The pressure correction to the liquid enthalpy and the flat-interface
chemical-potential admission include the dry partial pressure.

At the interface the noncondensable species has zero relative mass flux,
while vapor has relative speed `-j/rho_v`. Summing species momentum gives
the retained recoil `j^2(1/rho_v-1/rho_liquid)`; substituting total gas
density would discard the species diffusion momentum. Existing vapor
velocity/omitted-energy screens remain, and the new owner additionally
admits total interior density and frozen-mixture Mach/pressure relaxation.
The latter sound derivative follows the additive original real-vapor and
ideal-dry Helmholtz approximation. It is an admission quantity; emitted
sound still propagates through the original surrounding liquid.

The same fixed-stack adaptive solve now carries two additional gas
quadratures, `integral m_g Cv_g Tdot dt` and `integral p_g Vdot dt`.
The original source `m_g integral Cv_g dT` independently checks gas
storage. Vapor-only energy remains separately reported, while contents
energy, pressure work and the combined liquid/contents thermal balance
include the retained gas. Formation references never enter a gas increment.
The original gas source is qualified at the accepted candidate's actual
temperature and dry partial pressure before publication.

Uniform interior temperature/composition is a bounded approximation.
`GasVaporTransportBounds` requires explicit lower bounds on WET-mixture
conductivity and water/gas molecular diffusivity over declared temperature,
total-pressure and vapor-fraction intervals. Dry-gas conductivity is not
substituted. Each stage checks implied interior temperature variation from
`Q_liquid-L*mdot_v`, relative gradients of both species, and thermal and
diffusion times against actual state-change rates and an authored imposed
forcing-rate/bandwidth bound. These are applicability estimates, not
measured transport data or mathematical error certificates.

This owner does not resolve gas concentration/temperature boundary layers,
finite-rate interfacial accommodation, dissolution, cross-virial/fugacity
interactions, shocks or nonspherical collapse. It refuses declared mixing
limits instead of claiming those dynamics. The liquid grid and native
acoustic source remain shared, and a refused candidate or final cancellation
publishes no gas/phase/liquid history.

Fresh verification after the executor returned: the focused original-source material lane passed 33/33 tests, including all four conserved-gas groups, and the native `caloric-cavitation` lane passed 105/105 cavitation tests, including all six new mixed-gas groups and the retained pure-vapor, liquid, caloric, population and cluster regressions. Sources were scoped-formatted with Rust 2024 rustfmt. These local fallback lanes compile the actual product modules through the retained focused facades; DSR/RCH are unavailable, and no full-workspace or release-gate result is claimed.


## Thermal pressure feedback in finite Cartesian mixture tanks

Both `mixture::MixtureEuler2d` and `spatial::mixture::MixtureEuler3d` can
author `FiniteReservoirHeat` on already bound stationary finite tanks.
`with_finite_reservoir_heat` binds each original `WallThermal` and common
`WallControl` before stepping. The builder re-admits total numerical storage;
one prepared species buffer per heated side supplements the existing accepted
and private inventories. Large accepted/trial receipts remain in the original
prepared Vec storage rather than on every field owner's stack.

The original `fs-gas::multicomponent::network::Plenum::heat_candidate_into` supplies each heat
half. The first half heats or cools the accepted tank and wall before the first
native SSP residual. The full first port transfer supplies the second-stage
tank; the weighted final transfer then receives the second original heat half.
Actual mixed composition, Cp(T), energy references, choking, backflow and both
Cartesian shears stay under the original gas source. Heated face signals enter
the existing field and compact-volume CFL checks. A rejected attempt restarts
both gas and wall from accepted storage.

Finite-wall heat is internal to tank plus wall. Prescribed-temperature bath
heat is external. `FiniteReservoirHeatReport` retains original before/after
thermal contacts, the summed source receipt and the coupled gas/wall residual.
The native port transfer still carries full total energy and support impulse;
no duplicate enthalpy, mechanical work or reaction heat is added. Tank pressure,
field cells, wall temperature, receipts and accepted getters publish together
after the composing owner's last cancellation poll. Public observations use
accepted inventory explicitly. Prepared clone/reset copies retain wall storage,
controls and accepted receipts without successful warm-step allocations.

The thermal source is the original backward-Euler law and remains first order
in time. Symmetric composition with the second-order native port method does
not upgrade that wall-source order. The feature supports prescribed baths and
constant-capacity finite walls; a caloric wall without its immutable material
binding refuses through the original source. Tank chemistry, changing volume,
internal thermal gradients and finite mixing times remain outside this owner.

Six new actual-source native groups passed: finite-wall heating/cooling and
pressure-driven flow on all ten 2D/3D sides with both reconstructions; external
bath balances and adiabatic bit parity; last-poll rollback and deterministic
retry; prepared heated storage and split restore; budgets and source refusals;
and hot-wall CFL retry equivalence to a direct accepted interval. All twelve
existing shared/planar/spatial finite-mixture groups passed, including off-stack
storage, two-stage convergence, support work and small-tank admission. The same
bounded inviscid closure also passed the seven frozen-planar finite-tank groups.
It retains the actual field reconstruction, forcing, hydrostatics, native RK
and timed-stage implementations with original source kernels. Molecular and
reactive orchestration, material-wall wrappers and unrelated dependency owners
were not compiled or executed; prepared split-copy behavior is verified at its
shared storage seam. This is not a full workspace Cargo/DSR gate.


### Held body forcing on moving frozen-gas volumes

`spatial::Euler3d::step_moving_with_forcing` surrounds conservative ALE transport
with the original mechanical source half-kicks. Acceleration acts on actual
stage mass; force density acts on the actual old/new physical volume. The
native source law supplies all three impulses and represented kinetic work,
including its original Coriolis behavior, without a separate thermal source.
`StepReport.source_change` keeps this exchange separate from boundary transfer.

Forced accepted, kicked, RK and reconstructed physical states obey the authored
laboratory Mach bound. Acoustic transport uses actual mesh-relative waves and
admits the final kicked endpoint. The reconstructed trace consumes source-owned
sound speed through `FrozenGas::primitive_from_density_pressure`, without a
conservative energy construction and caloric inversion. Prepared force and RK
buffers are reused; refusals restart the complete K/2--ALE--K/2 proposal.
Cells, geometry, clock, reports and loads publish after the final checkpoint.
Zero-force and zero-motion cases delegate to their existing native entries.

Cellwise coefficients are held in original cell order through each attempted
interval, not resampled at moving coordinates. Existing one-axis affine
geometry and boundary admission remain authoritative. Five native tests cover
exact vector impulse/work, boost/relative-CFL behavior, actual-volume force
quadrature on all axes, zero identities, reconstructed/final-kick Mach refusal,
and middle/final cancellation with exact replay.

## Conservative joint planar species mixing

`scalar::bounded::species::PreparedBoundedSpeciesTransport` transports all
authored mass fractions together on the original P1-DG/BDM scalar owner.
Admission requires the same borrowed mesh, effective molecular-plus-SGS
diffusivity, spatial matrix and complete weak boundary coefficients for every
species. Original value/inflow reservoirs may carry different compositions;
all original Gauss samples and physical affine endpoint extrema must lie in
the unit simplex.

One correction graph is shared. Each interior and signed-reservoir
antidiffusive vector receives the smallest admissible factor across its
nonzero species components. Thus every species keeps its own paired
conservative inventory, while the zero-sum correction vector preserves the
composition sum. All species evolve independently; none is reconstructed as
a complement. No clipping, renormalization or inventory repair occurs.
The original scalar preparation/application, consistent mass, SSPRK2 stages,
dyadic CFL/work bounds and boundary receipts are reused. Each component
retains [0,1] to the original scalar roundoff allowance; unit sums are
independently admitted at every Euler/RK endpoint with absolute allowance
1024 machine epsilons. Cancellation or refusal returns no partial mixture.

Per-species receipts retain original advective/SIP exchanges and separately
report unremoved numerical boundary diffusion. Work scales with species times
sparse transport work; the graph is stored once and stage scratch is reused.
This is frozen, constant-density, source-free passive mixing with shared
effective diffusivity, compatible with Smagorinsky, WALE and Vreman. It does
not add differential/Maxwell-Stefan diffusion, reactions, species enthalpy
carriage or variable-density coupling; nonzero imposed diffusive fluxes remain
outside its admitted boundary policy.

Seventeen selected native checks pass: five new species groups, ten existing
scalar integration tests and both existing scalar numerical oracles. They
demonstrate loss of unit sum under independent scalar limiting, joint front
transport, closed skew-grid conservation of every species including an absent
component, species permutation symmetry, all-model reversing signed-reservoir
exchange, zero summed species diffusion, smooth inactive-limiter recovery and
transactional input/work/cancellation refusal. The original consistent-mass
LU and extreme-contrast diffusion oracles remain green. This is focused
original-numerical-source validation, not full Cargo/runtime or experimental
mixing qualification.


### Body forcing in moving frozen molecular flow and finite pistons

`spatial::molecular::MolecularEuler3d::step_moving_with_forcing` composes
the original molecular and mechanical sources as
`D/2–K/2–ALE–K/2–D/2`. Each phase uses its actual incoming state and
old/new geometry. The final molecular endpoint obeys laboratory Mach and
mesh-relative acoustic admission. A shortened central interval or a later
refusal restarts the complete split. Direct molecular wall heat and traction
work remain separate from the source-owned `ForceChange3d`.

Both frozen Euler and molecular Euler provide
`step_with_inertial_walls_and_forcing`. Their original simultaneous
`solve_coupled` mechanics consumes actual pressure and molecular face
impulses. Prescribed gas-body impulse and kinetic work remain explicit external
terms and are removed from combined gas/piston mechanical acceptance. Physical
and referenced energy are checked independently; fixed-side support impulse,
transverse guide reactions, wall torques and thermal exchange retain their
original accounts. Forces on piston bodies remain explicitly authored through
their mechanical loads.

Prepared source scratch is reused without full-field allocation during trials.
Zero forcing delegates to the original unforced transaction. Accepted gas,
geometry, clock, loads, reports, molecular binding and piston states publish
only after the final cancellation checkpoint. Cellwise force coefficients are
held in original cell order throughout motion and retries. Existing closed or
periodic moving-domain and single-axis affine-deformation restrictions remain
authoritative.


### Body impulse and work in moving three-dimensional mixtures

`spatial::mixture::MixtureEuler3d::step_moving_with_forcing` composes the
original held mechanical source half-kicks around conservative ALE transport.
The first kick uses the old physical cell volume; the second uses the actual
new volume. Acceleration acts on stage mass, while force density acts on stage
volume. Canonical species and their original energy references are preserved
by the kicks. The native source law supplies represented vector impulse and
kinetic work, including its midpoint or exact Coriolis behavior.

The mixture step report retains all body impulses and work separately from
moving-boundary transfer. Referenced and restored physical energy residuals
include that source exactly once. Every forced cell, source endpoint, RK stage
and reconstructed physical trace obeys the caller's laboratory Mach bound.
Moving acoustic CFL uses the actual velocities of both faces of each deforming
cell. The face check consumes the original mixture primitive constructor and
its sound speed without a nonlinear caloric inversion.

Prepared source scratch and native extensive RK buffers are reused. Every
refusal retries both source halves and ALE from accepted state; gas, geometry,
clock, boundary loads and reports publish only after the final checkpoint.
Zero force delegates to original moving transport; zero motion delegates to
stationary forced transport. Cellwise coefficients are held in original cell
order over each attempted interval. Existing domain restrictions remain in
force; this entry does not imply hydrostatic balancing, force-field resampling
or a moving-frame transformation.

### Held body forcing in moving and finite-piston reacting mixtures

`ReactingEuler3d::step_moving_with_forcing` composes the original chemical
half-solves around the native forced ALE interval:
`C/2–K/2–ALE–K/2–C/2`. The molecular variant
`step_moving_with_dissipation_and_forcing` encloses that sequence in the
original `D/2` phases. Each stage uses its actual incoming state and old/new
geometry. Source impulse and kinetic work are retained in the reacting report
once; chemical reference changes, molecular heat and boundary work retain
their original independent accounts. Physical Mach admission applies at
chemical/molecular intermediate endpoints as well as the native forced
transport stages. The final actual-geometry endpoint obeys relative acoustic
CFL before acceptance.

The inviscid and molecular reacting finite-piston forcing entries reuse the
same simultaneous finite-mass mechanics and complete reacting trial. Gas body
impulse and work are prescribed external inputs, removed explicitly from
combined mechanical residuals in both physical and referenced energy charts.
Actual pressure and molecular traction drive the pistons. Direct wall heat,
fixed-side support impulse, guide reactions and wall torques remain explicit.
Gas forcing does not infer forces on the piston bodies; those loads are still
authored separately.

All chemistry, source, transport, molecular and mechanical refusals retry the
complete interval from accepted owners. Prepared forcing scratch is included
in the aggregate workspace budget, including inactive forcing, and is reused
across trial solves. Zero forcing retains the original numerical transaction.
Gas state, geometry, clock, wall impulses, reports, molecular binding and both
piston states publish together after the final checkpoint. Force coefficients
are held in original cell order through motion and retries. Existing closed
or periodic domains and single-axis affine deformation remain authoritative.


### Body forcing in moving molecular mixtures and finite pistons

Mixture `MolecularEuler3d::step_moving_with_forcing` composes
`D/2–K/2–ALE–K/2–D/2` on actual old/new volumes. Evolving species,
source-owned mechanical work, and original Maxwell–Stefan diffusion, species
enthalpy, Fourier heat and full molecular stress remain in their existing
accounts. Final laboratory Mach and mesh-relative acoustic admission precede
publication. The complete owner budget includes prepared forcing scratch.

Both mixture Euler and molecular Euler expose
`step_with_inertial_walls_and_forcing`, reusing the original bounded
simultaneous finite-mass solve. Every root proposal replays one unpublished
physical interval. Actual gas-body impulse and work are explicit external
inputs, excluded from piston face loads and subtracted from combined momentum
and both referenced/physical energy gates. Molecular heat, fixed-side supports,
guide impulse and torque remain explicit. Each piston receives its actual
accepted boundary traction and work; additional body loads on pistons remain
separately authored.

Zero-source reductions preserve the original compensated summation order.
Cancellation or refusal preserves accepted cells, geometry, clock, molecular
binding, reports, face loads and both piston states. The original source and
transport workspaces are reused without per-proposal full-field allocation.
Coefficients remain held in original cell order, and the existing moving-domain
and single-axis deformation limits continue to apply.


## Generalized prescribed-boundary source work

The original mixed reaction owner now retains `BoundaryPower` at the accepted
native endpoint. Its sign is explicitly power supplied BY the exterior TO the
fluid, opposite to the existing fluid-on-exterior wrench. On a prescribed
velocity face, the work is `R0*m0(g)+R1*m1(g)` plus the FULL retained numerical
Nitsche flux `(2*nu*D(u)*n-eta*(u-g)).g`, integrated with the original admitted
samples. This is the discrete work conjugate to the prescribed data; all
adjoint-consistency, pressure, history, convection and drag terms stay in the
original discarded-row residual. Arbitrary nonuniform motion is not replaced
with a rigid wrench or average velocity. No primal operator or boundary law
changes.

Navier-slip faces instead contract the original physical Robin effort with
the actual wall velocity. The separate nonnegative `beta*|(u-g)_t|^2` receipt
continues to describe relative slip dissipation, not a second wall-work term.
Supplied open traction contracts the accepted fluid trace. A resolved original
pressure port updates that work with `-p*Q`, where Q is its accepted fluid
normal flux, never a possibly different Oseen advector. Dissipative-backflow
stress work is retained separately with its nonpositive sign, using the SAME
original reversing-face convection quadrature or conservative relative-mass
samples as the operator. Natural traction and supplied reservoir momentum do
not invent this resistance. Transported kinetic energy and numerical upwind
loss keep their original independent meanings.

Per-edge and canonical selected-edge queries reuse retained scalars, with
compensated sums and no source, material or boundary requery. Mean-zero
components require an explicit pressure offset, including zero; traction
datums require None. A pressure offset changes prescribed power by
`-delta_p*integral(g.n)`. This matters for a selected inlet/outlet even when
the complete sealed-component work is gauge independent. Invalid selections,
incompatible references and unrepresentable observations refuse.

Native powers are per-depth divided by the reference density, in m4/s3.
`MaterialBoundaryPower` and homogeneous, nonlinear-rheology and heterogeneous
field accessors use the SAME retained reference density for kinematic
pressure and W/m conversion. Heterogeneous cells do not apply local-density
rescaling to pressure or work. A consumer supplies physical extrusion depth
and duration when forming joules. The result is a source-work observation,
not a pointwise Cauchy-stress estimate, a full discrete energy/passivity
certificate, a moving-domain update or an automatic heat source.

Independent regression definitions cover nonrigid affine strain with a linear
pressure field (including the edge pressure/velocity covariance missed by a
midpoint-force proxy), signed work and explicit pressure shifts; moving finite
Navier-slip walls with separate bulk/slip losses; exact cubic backflow work on
a reversing original frozen advector; resolved pressure-port work under flow
reversal; and original-source retention plus homogeneous/heterogeneous SI
density swaps. All thirteen changed Rust files parsed with the pinned
formatter. The ten production/test owners outside the three narrow existing
registration/export parents passed rustfmt --check; those parents retain
unrelated pre-existing formatting differences. Native FEM/material test
execution remains unverified in the current restored environment.
