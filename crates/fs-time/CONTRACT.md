# fs-time CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-ad, fs-blake3, fs-exec, fs-ga, fs-ivl, fs-la,
fs-math, and the shared fs-solver operator/Newton-Krylov spine — all downward
or same-layer services).
Structure-preserving time integration (plan §8.5): integrators that
preserve what the physics preserves — symplectic (Störmer–Verlet with
its discrete-Lagrangian equivalence tested), Lie-group SO(3) and SE(3)
integration over the validated canonical `fs-ga::{So3, Se3, Twist}` owners
(exp-map lane
plus a free-body discrete Euler–Poincaré variational lane),
generalized-α with controllable dissipation, IMEX and exponential
integrators for stiffness, and embedded-pair adaptivity with a PI
controller. The two universal obligations (P7 + §8.7) — resumable
state machines and discrete adjoints of the stepper — are discharged
where claimed below.

## Public types and semantics

- `symplectic::verlet_step(q, p, h, force, scratch)` — one
  kick–drift–kick step for separable H = ½‖p‖² + V(q), unit mass
  (scale q/p externally otherwise). `force` writes F(q) = −∇V.
- `symplectic::verlet_adjoint(q0, p0, h, steps, force, force_jvp,
  (bar_q_N, bar_p_N)) -> (bar_q0, bar_p0)` — discrete adjoint OF THE
  STEPPER (transposed tangent map, reverse order), propagated with
  `fs_ad::revolve` checkpointing (O(log N) memory).
  **Precondition: symmetric ∂F/∂q** (conservative forces): the reverse
  pass evaluates Jᵀv through the user's `force_jvp` as Jᵀv = Jv.
- `symplectic::verlet_momentum_step(q, p, h, inverse_mass, force, scratch,
  cancelled)` retains physical canonical momentum for a constant mass operator.
  Kicks use generalized force; the drift uses the supplied `LinearOp` action
  `M^-1 p`, including coupled consistent mass. It assembles no matrix, allocates
  nothing, and uses three caller-owned state vectors of scratch. Nonfinite
  actions, invalid dimensions/step and cancellation leave both input states
  unchanged. Actions overwrite their complete outputs. The original unit-mass
  API and its arithmetic are unchanged.
  `verlet_momentum_adjoint` pulls back the same discrete trajectory through the
  supplied force VJP and **inverse-mass transpose**, with fallible binary
  checkpointing and immediate cancellation propagation. The mass/force actions
  must remain pure and unchanged through replay. A force JVP is insufficient
  unless its Jacobian is known symmetric. The derivative is with respect to
  initial position and momentum, not mass parameters. Symplectic/physical
  energy claims require caller-established constant SPD inverse mass and
  conservative force; arbitrary `LinearOp` instances are not certified. Variable
  mass, constraints and impacts are outside this separable lane. Polls bracket
  actions and occur every 256 vector components; actions own their internal
  latency and accepted forward publication uses two final slice copies.
- `verlet_momentum_step_with_displacement` performs the same map with four
  caller-owned scratch vectors and returns a borrowed slice containing the
  resolved drift `h M^-1 p_half`. Held-force work uses this increment even when
  adding it to an absolute position rounds it away. There is no additional
  force/mass evaluation or allocation. Displacement arithmetic and cancellation
  are checked before state publication; a nonrepresentable increment refuses
  even if the fused position update is finite. Numerical storage defects remain
  separate from supplied work and physical dissipation. Focused `time_battery`
  regressions cover coupled-mass ballistic work below the position ULP, reversal,
  unchanged nonlinear step bits, overflow and final-poll cancellation.
- `verlet_momentum_parameter_adjoint` extends the same checkpoint/reverse sweep
  to design-dependent constant physical mass and force laws. The
  `MomentumVerletParameterVjp` actions supply force-parameter derivatives and
  derivatives of `M(theta) v`, not an assembled inverse or solver-iteration
  tape. Drift contributions are `-h z^T (dM/dtheta) v`, with the actual retained
  primal drift velocity and the original inverse-mass transpose action. Both
  force kicks contribute; no extra mass solve or state-sized column per
  parameter is needed. `MomentumVerletGradient` returns initial position,
  initial canonical-momentum and parameter cotangents. Fixed initial velocity
  requires chaining `p0=M(theta)v0`; direct objective derivatives are separate.
  Matching mass/inverse/derivative actions remain caller obligations, and an
  approximate inverse carries its own derivative defect. Completely written
  finite parameter outputs, bounded cancellation and atomic refusal are checked.
  Zero-step gradients have zero parameter part. The focused `time_battery`
  tests check analytic ballistic sensitivity, independent dual/finite-difference
  gradients for coupled mass and nonlinear internal forces, actual design
  descent, original initial-state arithmetic, replay and derivative refusals.
- `fs_ga::{So3, Se3, Twist, Vec3}` are re-exported canonical state/tangent
  types; fs-time does not own a second quaternion, motor, pose, or twist
  representation. `lie::{so3_body_exp_step, so3_space_exp_step}` name the
  right/body update `R Exp(h ω_body)` and left/space update
  `Exp(h ω_space) R` explicitly. `lie::rigid_body_step(R, ω, I, h)`
  advances free rotation with a symmetric Hamiltonian splitting and returns
  `(R′, ω′)`. `rigid_body_momentum_step(R, L, I, h)` exposes body angular
  momentum directly, avoiding repeated velocity/momentum conversions.
  Each exact subflow pairs the attitude rotation with its inverse acting on
  body momentum, preserving spatial angular momentum to roundoff. Extracting
  the commuting spherical Hamiltonian at the median principal inertia makes
  spherical and axisymmetric tops exact; asymmetric tops have second-order,
  reversible symplectic evolution with a tested bounded energy envelope.
  Both APIs allocate nothing, iterate no nonlinear solve, and do not
  renormalize the stored attitude; invalid/non-finite states refuse.
- `lie::compensated::{rigid_body_momentum_step,
  rigid_body_momentum_step_in_frame}` retains `[Dd; 3]` body momentum through
  the same symmetric split. The latter accepts cached positive principal
  moments and a fixed principal-to-body `So3`; `InertiaTensor3` exposes it as
  `compensated_momentum_step`. Residual-flow phases, the spherical rotation
  vector and quaternion compositions retain impulse tails through the
  existing deterministic DD trigonometry, so small impulses affect attitude
  as well as spin. Projective quaternion action pairs each increment with
  its inverse momentum transport, including the full-tensor frame change;
  the spherical flow leaves momentum analytically unchanged. The incoming
  attitude participates in DD composition before the one final native `So3`
  admission. Fixed-size arithmetic adds no allocations, factorization,
  nonlinear solve, renormalization or momentum projection. Nonfinite or
  unrepresentable derived arithmetic refuses; this is neither an interval
  certificate nor an exact-energy claim for asymmetric finite steps.
- `lie::compensated::{rigid_body_momentum_step_fourth_order,
  rigid_body_momentum_step_fourth_order_in_frame}` advances smooth free rotation
  with the symmetric fourth-order Yoshida composition. The cached full tensor
  exposes this as `InertiaTensor3::compensated_momentum_step_fourth_order`.
  Adjacent equal-axis flows merge, and the commuting spherical flow runs once:
  asymmetric rotation uses eight elementary flows instead of twelve for three
  complete second-order steps. All phases, products and paired inverse momentum
  transports retain DD tails. Spherical/axisymmetric cases retain the original
  exact flow bitwise. `tests/lie_compensated.rs` verifies fourth-order convergence
  against an independent continuous Euler/quaternion RK4 oracle, phase accuracy
  at equal elementary-flow count, spatial momentum, retained momentum norm,
  signed-time reversibility, and exact-top parity. Negative internal substeps
  apply only to this conservative free flow. This entry neither increases the
  order of separately composed loads/contact/damping nor claims exact asymmetric
  energy, a discrete adjoint, or a wall-clock speedup.
- `lie::rigid_body_torque_step(R, omega_body, I, torque_body, h)` and
  `lie::rigid_body_spatial_torque_step(R, omega_body, I, torque_space, h)`
  apply held center-of-mass torques through symmetric momentum half-kicks
  around the same Hamiltonian free drift. Both are second order and time
  reversible for the declared held load, with no allocation or nonlinear
  solve. The spatial version rotates the torque at **both** endpoints and
  satisfies `R_next I omega_next - R I omega = h torque_space` to roundoff.
  The body version represents a body-attached load; its spatial impulse is
  the trapezoidal sum of the two rotated endpoint torques. These helpers do
  not claim an adjoint or second-order coupled accuracy when a changing
  torque is frozen once per step. Such consumers must evaluate their force
  law at the required stages. Loaded steps retain `I omega`, both impulse
  tails and the final inertia division in DD around the compensated free
  owner. The public state stays f64; a representable full impulse is not
  discarded merely because its two separately rounded half-kicks vanish.
  Exactly zero torque (including signed zeros) delegates to the unchanged
  scalar free owner, preserving its established bitwise parity and refusals.
- `lie::rigid_body_torque_function_step` and
  `lie::gyrostat::InertiaTensor3::torque_function_step` evaluate a fallible
  configuration/time torque law at both physical endpoint attitudes and
  times. Each `HeldTorque` result selects that sample's body or spatial frame;
  it does not freeze the force law for the step. Two symmetric momentum
  half-kicks enclose the existing compensated free rotation. Smooth attitude
  and time loads are second order; deterministic autonomous conservative
  loads retain symmetric Hamiltonian splitting. Spatial momentum changes
  by the trapezoidal sum of the two spatial endpoint torques to roundoff;
  energy is not exact. The dense-tensor owner reuses its cached principal
  basis without an eigensolve, allocation, or nonlinear force solve per step.
  Invalid time/duration/inertia/state is refused before force evaluation;
  callback and arithmetic failures return no partial physical state. The
  endpoint time must be finite and distinct from the initial time. These
  APIs do not cover velocity-dependent loads, rotor exchange, contact events,
  moving mass or callback adjoints; those require their dedicated owners.
- `lie::gyrostat::{momentum_step, carrier_momentum_step, driven_step,
  HeldTorque}` advances a carrier with internal rotor momentum `h`.
  `momentum_step` retains total body momentum `L = I omega + h`;
  `carrier_momentum_step` and the angular-velocity-driven path retain
  `k = I omega` separately so a larger stored rotor cannot round away
  carrier motion. `I` is the reduced principal inertia (rotor transverse
  inertia retained, axial spin inertia excluded); rotor momenta are absolute
  axial momenta, including the carrier's axial angular velocity. A symmetric
  composition of exact axis Hamiltonian flows preserves `R L` to roundoff
  and is second order and reversible for fixed `h`. Zero `h` delegates
  exactly to the rigid-body drift; a stationary carrier remains stationary
  for any stored rotor momentum. Spherical carriers and exactly axisymmetric
  carriers with rotor momentum on their symmetry axis use commuting exact
  Hamiltonian flows that preserve their continuous free
  trajectory and kinetic energy to roundoff. Only exact stored inertia/axis
  equalities select that specialization; near symmetry keeps the original
  general model. Pure axial motion avoids cancellation between large rotor
  phases, and exact carrier rest is preserved bitwise. General gyrostats need
  at most five elementary rotations, no allocation and no nonlinear iteration.
  `driven_step` composes internal
  rotor transfers with the drift and endpoint external torque half-kicks;
  motor reaction changes carrier momentum without creating external impulse.
  The separate-carrier drift evaluates each rotor transport as
  `(F^-1-Id)h` using the projective quaternion action without subtracting
  two large rotated momenta. Symmetric flows retain the small attitude phase
  separately from the rotor phase and rotate carrier momentum directly.
  Endpoint kicks combine external impulse and internal transfer before
  adding them to `k`, retaining the rotor endpoint subtraction in fixed-size
  double-double arithmetic as well. The public
  total-momentum drift and zero-rotor driven arithmetic are unchanged.
  This cannot recover carrier information already lost in a total-`L`
  input, nor repair a physically misaligned finite-precision input frame.
  Body-held versus spatially held torques are explicit. This map is second
  order for linear rotor momentum change and the declared held torque;
  changing load laws need a coupled stage. `tests/gyrostat.rs` checks analytic
  precession/refinement, rotor spin-up recoil, spatial impulse, frame
  covariance, reversibility, replay, rest, rigid reduction and admission.
  `tests/gyrostat_carrier.rs` additionally checks high-spin carrier motion,
  external acceleration and internal recoil, opposing impulses, exact signed
  symmetric trajectories including zero total momentum, retained small
  attitude phases, and principal/full-tensor coordinate consistency.
  General asymmetric or misaligned gyrostats retain the finite-step splitting
  and its measured kinetic defect; exact-energy or arbitrary-step trajectory
  accuracy is not claimed for that general case. Gimbaled rotors, flexibility,
  changing inertia and a discrete adjoint remain outside this extension.
- `lie::gyrostat::InertiaTensor3` admits a full symmetric positive-definite
  inertia in the caller's fixed body frame. `new(Mat3)` uses the existing
  native 3-by-3 Jacobi eigensolver once and checks the principal-frame
  reconstruction; `from_principal` accepts an already known principal frame
  without allocation. Cached momentum/velocity transforms and free/driven
  integration then have constant work and no allocation or factorization.
  The free map conjugates the native rotation increment and transports the
  original body momentum by its inverse, avoiding repeated full-state basis
  roundtrips. Its separate-carrier map conjugates the same attitude
  increment and applies the principal carrier change to the original body
  carrier; it never transports a large rotor by an approximate basis
  roundtrip. Body/spatial external torque conventions and internal rotor
  momentum remain explicit; strictly diagonal input retains the original
  map's arithmetic. `tests/gyrostat_tensor.rs` checks independent dense Euler
  dynamics (second-order refinement), 10,000-step spatial momentum balance,
  full-tensor reconstruction and inversion across scales, recoil and external
  impulses, world-frame covariance, reversal, replay and refusal. Admission
  and its reconstruction residual are numerical, not an interval eigenvalue
  certificate or geometric realizability proof. The tensor must remain fixed;
  this does not supply changing-mass, deforming-body or gimbaled-axis transport.
- `lie::gyrostat::rotors::RotorBank<N>` binds up to 64 balanced axisymmetric
  rotors with fixed body axes and complete reduced inertia. `from_locked`
  removes only `sum_i J_i a_i a_i^T` from the original complete locked tensor;
  transverse rotor inertia remains. Its canonical state stores total body
  momentum, absolute axial rotor momenta and unwrapped relative shaft angles.
  `apply_impulse` applies all internal shaft and external body impulses
  simultaneously, retaining full-tensor cross-axis carrier recoil. Internal
  transfers preserve total momentum. Complete kinetic energy is
  `T = omega.I_reduced.omega/2 + sum_i p_i^2/(2 J_i)`; shaft power uses the
  relative rate `p_i/J_i - a_i.omega`, including the carrier's axial motion.
  Actual rounded momentum changes and polarized quadratic increments give
  individual signed shaft/brake work, external work and a visible balance
  residual without subtracting large endpoint energies. An independent
  compensated transfer balance checks the actual carrier recoil before work
  accounting; its allowance uses carrier/transfer scale so a large stored
  rotor cannot hide a lost or materially distorted shaft transfer.
  `step` composes two held-torque half-kicks with one native gyrostat drift.
  Relative phase uses the same second-order relative-speed quadrature as
  shaft work; representable held torque times angle change agrees with that
  work to rounding. Body/space torque frames remain explicit. The native
  drift's kinetic defect is reported separately and is never labeled motor
  loss or heat. Stepping is O(N), allocation-free and has no new iterative
  solve. Invalid or completely unresolved nonzero impulses/phase increments
  refuse without partial state. Seven rotor-bank regressions cover independent
  locked energy, multiple-axis recoil, analytic spin-up/braking and 3D changing
  spin precession, relative phase/work, external impulse frames, reversal,
  generic refinement, empty-bank reduction and refusal. With nine principal
  gyrostat and seven full-tensor tests, all 23 focused native groups pass.
  Fixed axes, inertia and held torques are required: moving gimbals, electrical
  states and variable-speed bearing/bending dynamics are not supplied.
- `RotorBank::viscous_stepper(damping, duration)` prepares the physical
  relative-shaft law `tau_i = drive_i - damping_i * relative_speed_i` for
  nonnegative damping and a positive fixed interval. Two coupled implicit
  midpoint shaft stages surround the existing native gyrostat drift. Diagonal
  elimination caches one native 3x3 Cholesky factor; each step is O(N) and
  allocation-free, with full carrier recoil, nonorthogonal rotor axes and
  simultaneous external torque. Signed drive work, nonnegative viscous heat,
  external work and the native drift defect remain separate. The actual
  rounded midpoint impulse law and complete stage work are checked at local
  transfer/work scales; large background spin cannot hide a work discrepancy.
  `viscous_relaxation_rate_bound` exposes a PSD trace bound for the coupled
  carrier/rotor dissipative mobility. Preparation refuses a half-stage
  relaxation number above two to prevent unresolved midpoint sign reversal.
  Accuracy still requires refinement; this is not an L-stable or arbitrary
  stiff-step decay claim. All-zero damping delegates exactly to the existing
  held-torque motion and work path, including its explicit accepted-impulse
  rounding. Five focused native groups cover closed exponential/phase/heat
  refinement, independent full-tensor Euler dynamics, resolved stiff drag,
  signed braking and drive/heat balance, zero-drag equivalence, atomic
  refusals and a reproduced large-carrier work-loss counterexample. They and
  the seven existing rotor-bank groups pass; the full Cargo target remains
  unexecuted where the constellation runtime dependency closure is absent.
- `RotorBank::torsional_stepper(stiffness, rest_phase, damping, duration)`
  adds finite rotor/carrier shaft compliance on actual unwrapped relative
  phases: `V_i=K_i*(q_i-q_rest_i)^2/2` and
  `tau_i=drive_i-K_i*(q_i-q_rest_i)-d_i*s_i`. Original positive inertia,
  nonorthogonal fixed rotor axes, simultaneous external impulses and full
  carrier recoil are retained. Native gyrostat drift evolves attitude/total
  body momentum at fixed shaft phases and momenta. Two shaft midpoint stages
  evolve canonical phase/momentum; the second uses the first accepted phase.
  With half duration delta, midpoint phase elimination gives
  `drive_eff=drive-K*(q_before-q_rest)` and `d_eff=d+delta*K/2`.
  This reuses the existing affine rank-three solve and one cached native 3x3
  factor; stepping is O(N), allocation-free, without a new numerical solver.
  Trajectory completion consumes already-integrated phases rather than
  advancing them twice. Polarized actual phase increments determine spring
  storage changes. Drive work, viscous heat, elastic storage, external work
  and native drift defect remain separate. Actual stored phase, original
  spring/drag effort and complete stage work must resolve at local scales;
  large initial phase or potential energy cannot excuse lost increments.
  Zero stiffness delegates bit-for-bit to the original viscous owner.
  Preparation admits the combined elastic/viscous half-stage resolution
  number at most two and exposes a numerical PSD-trace undamped-frequency
  estimate. These are not arbitrary-step accuracy or interval proofs.
  Four new groups plus five viscous/seven rotor-bank groups pass on actual
  native sources. Harmonic-period, analytic driven/damped and independent
  oblique full-tensor Euler ladders all demonstrate second-order refinement;
  conservative harmonic total kinetic-plus-spring drift stays below
  3.34e-16 J in the tested period. Unresolved phase/stiffness and malformed
  parameters refuse without publishing a partial step. No global
  symplecticity, moving gimbals, variable stiffness/rest phase or parametric
  preload work is claimed. Full Cargo remains unavailable with missing
  constellation siblings.
- `RotorBank::torsional_network_stepper(couplings, duration)` adds finite
  elastic/viscous transmission BETWEEN distinct original rotor shafts.
  Each authored port has strain `x=q_first-ratio*q_second-rest`, potential
  `k*x^2/2`, relative rate `s_first-ratio*s_second`, and reciprocal shaft
  torques `(-effort, ratio*effort)`. Negative ratios, parallel couplings and
  dependent closed coupling loops retain their actual physical coefficients;
  zero ratio is the first-shaft/carrier limit. Full-tensor carrier recoil and
  simultaneous external torque remain owned by the existing rotor bank.
  Two implicit midpoint shaft stages surround the native gyrostat drift.
  The sparse two-entry incidence rows and rank-three carrier term assemble
  `H=B*(diag(1/J)+A^T I^-1 A)*B^T`; one cached native Cholesky factor solves
  `I+delta/2*sqrt(D_eff)*H*sqrt(D_eff)`, with `D_eff=D+delta*K/2`.
  At most 64 couplings are admitted. Preparation is O(E^3), stored factors are
  O(E^2), and accepted steps use O(N+E^2) work and fixed stack arrays without
  allocation. Empty/zero-law networks delegate exactly to held-torque stepping.
  Per-coupling polarized storage changes and nonnegative viscous heat remain
  distinct from drive/external work and the original gyrostat defect. Original
  actual impulse, recoil, phase, constitutive and work checks reject unresolved
  transfers without publishing a half-step. The combined half-stage PSD-trace
  resolution number is at most two; frequency estimates and second-order
  splitting still require refinement, not an interval or exact-trajectory claim.
  `gyrostat_torsional_network` covers analytic two-shaft oscillation and loaded
  damping, full-tensor oblique-loop Euler refinement, reciprocal work, limits
  and refusals. Fixed topology, ratios, rest phases and coefficients are
  required; rigid gearing, backlash, clutch friction, gear-tooth/material
  validity, thermal allocation and parameter-change work are not inferred.
- `RotorBank::lock_mode(mask)` prepares retained ideal rotor/carrier shaft
  locks with `I_mode=I_reduced+sum_selected J_i a_i a_i^T`.
  Engagement applies simultaneous original-owner impulses with complete
  carrier recoil, retaining attitude, total momentum and every shaft phase.
  Aggregate perfectly inelastic mechanical loss is nonnegative; no per-lock
  heat partition or finite thermal destination is inferred.
  Locked `apply_impulse` admits unlocked-shaft and external BODY ports while
  solving zero-work internal reactions. Already-admitted zero ports preserve
  the canonical state exactly, without a numerical reprojection impulse.
  Retained drift uses the existing native gyrostat on `I_mode` with only
  unlocked stored rotor momentum. Selected absolute momenta evolve with
  precession, their relative phases remain bit-for-bit fixed, and unlocked
  phases use endpoint trapezoidal relative rates. Actual canonical kinetic
  reconstruction and native drift defect are separately measured; neither
  becomes physical damping. Releasing a lock means using the accepted
  canonical state unchanged in ordinary rotor dynamics or another selection.
  Preparation uses the original tensor owner; impulses and drift are O(N)
  and allocation-free, with the original actual-transfer/recoil gates and
  local-scale work admission. No finite friction capacity, event-time
  prediction, moving axes, controller dynamics or exact-energy claim is made.
  Five new native groups and sixteen existing rotor/viscous/torsional groups
  pass on original numerical sources. Coaxial engagement loss/recoil,
  nonorthogonal simultaneous locks/permutation/ports, retained phase, release,
  reversal, empty selection and erased carrier transfer are covered.
  An independent full-tensor Euler/RK4 reference gives 20/40/80-step combined
  errors 8.00266e-5/2.00061e-5/5.00150e-6, demonstrating second-order
  refinement in that fixture. Pinned formatting passes; the complete Cargo
  target remains unavailable with missing constellation runtime siblings.
- `RotorLock::step` advances held unlocked-shaft drives and a body- or
  space-held external torque while selected ideal shaft locks remain engaged.
  Two native constrained half-impulses surround the original retained-mode
  gyrostat drift. Locked inertia stays in the cached full carrier tensor;
  spatial torque is transformed at both actual endpoint attitudes. Selected
  shaft drives must be zero, and the input must already satisfy the locks;
  engagement remains a separate dissipative event. Locked phases stay
  bit-for-bit fixed. Unlocked phase uses the same rate quadrature as drive
  work. `LockedRotorStep` separates drive, ideal-reaction and external work,
  native drift error and the complete energy residual. This symmetric
  fixed-inertia approximation is second order, O(N), allocation-free and
  subject to the original 64-rotor bound and numerical admission. No brake
  capacity or state-dependent drive law is inferred. G1/G3 regressions in
  `tests/gyrostat_locked.rs` specify analytic coaxial recoil/work, independent
  Euler-equation refinement, spatial covariance, empty-selection and zero-drive
  parity, reversal and atomic refusal.
- `RotorBank::dry_clutch_stepper(law, duration)` advances one finite-capacity
  dry clutch or shaft/carrier brake. A two-shaft port has slip
  `v=s_first-ratio*s_second` and torques `(-effort, ratio*effort)`; signed
  ratios retain full carrier recoil. Authored capacities obey
  `static >= kinetic >= 0`. With fixed attitude and held loads,
  `v'=a-H*effort`, where `H=b*(diag(1/J)+A^T I^-1 A)*b^T`.
  Its exact scalar flow resolves finite slip, synchronization, static hold
  within capacity, and overloaded reversal. Each half-stage has at most two
  constant-torque pieces; two such stages surround one original native drift.
  Phase and work integrate each piece, including the stopped part of a step.
  Every original rotor impulse/recoil gate remains authoritative. An explicit
  retained branch admits synchronization rounding without a capture-speed
  deadband; fresh nonzero slip slides, and precession-created native slip is
  re-evaluated before the second stage. Sticking phases use the one constrained
  angle without projecting canonical momenta. Signed drive, friction and
  external work, nonnegative friction heat, integrated effort, sliding angle,
  half-stage event times and native drift defect remain distinct in receipts.
  Preparation and stepping are O(N), allocation-free and bounded by 64 rotors.
  Zero capacities delegate exactly to the original held-torque bank. Smooth
  sliding is second order; general precessing transitions require refinement,
  with no claimed event-time certificate or transition convergence order.
  `tests/gyrostat_clutch.rs` covers analytic stop/hold/breakaway/reversal,
  signed-ratio synchronization, external-work closure, independent oblique
  Euler refinement, precession release and atomic refusals. This fixed-port,
  isothermal dry law supplies no empirical material validity, wear, thermal
  destination, actuator dynamics, changing ratio or simultaneous clutch network.
- `InertiaTensor3::viscous_stepper` prepares the full material-frame law
  `tau = drive - D omega` for a fixed positive interval. Symmetric PSD damping
  is admitted with compensated principal-minor checks, including zero-row
  coupling refusal; no projection changes the supplied constitutive tensor.
  One normalized native Cholesky factor of `I + dt D/4` is cached. Each
  allocation-free step uses two implicit midpoint torque stages around the
  existing compensated free rigid-body flow; spatial drive is transformed
  at both actual endpoint attitudes. Impulse tails remain in DD until public
  momentum rounding. Physical heat is the positive midpoint quadratic form;
  drive work, viscous reaction impulse, native free-flow energy defect and
  final-state numerical balance are separate outputs. Actual accepted kinetic
  increments and free-drift defects use the polarized fixed-inertia quadratic
  form with the cached inverse-inertia action, retaining transverse work beside
  a dominant unchanged spin without subtracting total energies. This also
  applies to the rotating-bath and finite-pair entries below. The method is second
  order for the held linear law. The numerical trace bound
  `tr(I^-1 D)` admits only half-stage relaxation numbers at most two, avoiding
  unresolved midpoint sign reversals; it is not an interval certificate or a
  replacement for decay/phase refinement. The basic step entry uses a
  stationary reference environment. Nonlinear drag, dry friction, varying
  tensors, internal rotors and discrete adjoints remain outside this API. G0/G1/G3 tests
  in `tests/viscous_rigid_body.rs` specify exponential spin/phase/heat
  refinement, passivity, material-frame covariance, the native zero-drag
  limit and invalid/unresolved-input refusals.
- `ViscousRigidBodyStepper::step_in_rotating_bath` holds a prescribed angular
  velocity in space, transforms it at both accepted endpoint attitudes, and
  reuses the same cached midpoint factor for
  `tau_drag = -D (omega - R^-1 Omega)`. `ViscousRotatingBathStep` retains the
  body receipt and adds signed bath source work
  `sum(tau_drag dot Omega dt)`; heat is the nonnegative relative-slip
  quadratic form. The bath receives the opposite drag impulse. Body kinetic
  change closes against drive plus bath work minus heat and the native drift
  defect. Exactly zero bath uses the original stationary arithmetic. Four
  physical groups cover independent coaxial entrainment/reversal work and heat,
  full Euler/quaternion RK4 refinement, joint material/spatial covariance, and
  co-rotation, work extraction, parity and refusal. This is a prescribed source:
  finite bath inertia, changing bath laws, translational drag and thermal
  allocation are not inferred.
- `lie::viscous::pair::ViscousPairStepper` couples two original finite
  full-inertia rigid bodies through a damping tensor fixed in body 0's
  material axes. At both actual endpoint configurations,
  `r = omega0 - R0^-1 R1 omega1`; coupling torques are `-D r` and the
  opposite spatial couple. Two implicit midpoint half-kicks surround both
  original compensated free flows. Independently held body/space drives
  transform at the actual endpoint attitudes. Both bodies recoil.
  Inertia-scaled relative mobility `K = I0^-1 + Q I1^-1 Q^T` reduces each
  kick to two native 3x3 Cholesky factors and one solve, with bounded native
  factor storage and no damping inverse, iterative solve or per-step
  eigensolve. The original compensated principal-minor owner admits PSD,
  including rank-deficient, damping. Each actual configuration refuses
  unresolved mobility or `s tr(K D) > 2` with a numerical margin.
  Local actual-midpoint constitutive and per-body work checks precede return.
  Equal/opposite coupling impulses share one public world receipt; each
  body separately reports drive work, signed coupling transfer, kinetic
  change, native free-flow energy defect and final public-state residual.
  Shared nonnegative relative-slip heat has no assigned thermal destination.
  The symmetric splitting is second order for resolved smooth linear drag.
  Fixed inertias and body-0-attached damping are required; joints, translation,
  changing laws, dry friction, internal rotors and adjoints are not supplied.
  Extreme ill-conditioned mobility refuses instead of projecting. Six
  G0/G1/G3/G4 groups in `tests/viscous_pair.rs` cover independent analytic
  recoil/phase/work/heat and Euler/quaternion RK4 refinement, independent
  material/world covariance, zero-drag parity, common-spin/passivity and
  invalid or unresolved input.
  The complete crate and first pair revision compiled with the pinned native
  compiler and original cached dependencies. Five of six pair groups and all
  five stationary-viscous regressions passed. The analytic driven group
  exposed a component-relative constitutive refusal; the final check also
  includes the local common-coupling vector scale for frame/solve arithmetic,
  with the same numerical multiplier and no background-momentum or absolute
  floor. Per-body work checks remain separate. Executor loss prevented
  recompiling/rerunning this final correction, so final native acceptance and
  a Cargo/DSR gate are not claimed.
- `galpha::GeneralizedAlpha::new(m, c, k, n, h, rho_inf)` +
  `galpha_step(ga, q, v, a, f_next)` — Chung–Hulbert: αm = (2ρ−1)/(ρ+1),
  αf = ρ/(ρ+1), γ = ½ − αm + αf, β = ¼(1 − αm + αf)²; one prefactored
  LU of the effective matrix per (M, C, K, h); `f_next` is the load at
  t + (1−αf)h. The dense lane solves for endpoint acceleration with
  `(1−αm)M + (1−αf)γh C + (1−αf)βh² K` and updates position/velocity
  by Newmark increments, retaining resolved vibration below absolute-position
  resolution without dividing position differences by h². Fixed velocity and
  acceleration history operators are preassembled; history assembly uses one
  scratch n-vector and one fused traversal, plus the LU permutation scratch.
  `tests/galpha_probe.rs` checks
  tiny-step vibration and origin-independent rigid acceleration with a coupled
  mass matrix, alongside order and forced static response. Position itself is
  still ordinary `f64`; sub-ULP position accumulation is not claimed.
- `galpha::{SecondOrderProblem, LinearSecondOrderSystem}` describe
  `M a + C v + r(q) = f`; the linear adapter binds the exact shared
  `fs_solver::LinearOp` interface for M/C/K, while nonlinear problems provide
  matched internal-force and tangent actions. `OperatorGeneralizedAlpha`
  solves the same Chung–Hulbert residual through `NewtonKrylovState` and
  commits a `SecondOrderState` only after convergence. Its Newton unknown is
  the Newmark correction `beta*h²*a_next`; stage and endpoint kinematics use
  displacement/velocity increments, preserving resolved acceleration even
  when an absolute displacement increment rounds away. This leaves the
  physical residual, Jacobian and preconditioner weights unchanged and avoids
  two temporary vectors per residual/Jacobian evaluation. Position remains an
  ordinary `f64`, so sub-ULP position accumulation is not claimed.
  `SecondOrderProblem::resistance_balance` also retains the separate stage
  increment when evaluating resistance minus applied load. Its default calls
  the complete nonlinear resistance at the same fused stage coordinate.
  Declared linear models can use `linear_resistance_balance`, also used by the
  linear adapter, to subtract `K q_base - f` before adding the incremental
  elastic force and full `C v`. Consistent acceleration uses the same balance
  at zero increment; Newton and adjoint partials remain the original physical
  derivatives. This preserves small damping and acceleration around balanced
  preloads without claiming accumulation below absolute-displacement resolution.
  The preload regressions in `operator_integrators` and
  `generalized_alpha_initialization` compare coupled nonproportional damping
  with independent equilibrium-relative dense dynamics (G3).
  If the initial explicit predictor fails the model's stage admission, one
  bounded retry seeds the original configuration by negating the Newmark
  predictor drift. The resulting stage velocity is admitted independently;
  initial-state checks, physical residuals, nonlinear budgets, timestep and
  endpoint admission are unchanged. Cancellation is checked after both setup
  attempts. This does not guarantee a feasible Newton seed for every envelope.
  `tests/structural_admission.rs` checks a stiff spring whose exact implicit
  endpoint is inside its geometric domain despite a refused predictor, the
  acceleration pullback, velocity-envelope refusal and atomic cancellation.
- `galpha::FirstOrderGeneralizedAlpha` + `first_order_galpha_step` are the
  prefactored small-system lane for `M udot + A u = f`.
  `OperatorFirstOrderGeneralizedAlpha` consumes `FirstOrderProblem` (or the
  `LinearFirstOrderSystem` adapter) and uses the first-order-system
  Jansen–Whiting–Hulbert parameters αm = (3−ρ)/(2(1+ρ)), αf = 1/(1+ρ),
  γ = 1/2 + αm − αf. This is not the structural second-order formulation
  relabeled: state and rate are first-class and the residual is enforced at
  `(t_n + αf h, u_n + αf Δu)`. The dense lane solves for endpoint rate;
  the operator lane solves for `γh rate_{n+1}` and constructs `Δu` from
  that correction and the old rate. Neither recovers a rate by subtracting
  two absolute states. Thus a resolved acceleration survives even when a
  velocity increment rounds away against a large velocity offset. Stored
  state remains `f64`; sub-ULP state accumulation is not claimed.
- Both implicit problem traits accept an optional `preconditioner_apply`
  action for their forward Newton/FGMRES solve. It receives the exact stage
  state (and first-order stage time), named effective-operator coefficients,
  and logical outer/inner iteration indices. Structural weights are
  `((1−αm)/(βh²), (1−αf)γ/(βh), 1−αf)` for mass, damping and tangent;
  first-order weights are `(αm/(γh), αf)` for mass and tangent. The default
  remains identity. A callback must completely overwrite its output, bound
  its own work, and reproduce its policy during checkpoint replay. Invalid
  output fails before publishing the time step. The explicit adjoint
  preconditioner is independent; derivatives remain those of the converged
  physical residual, without differentiating preconditioning or iteration.
  `tests/first_order_preconditioner.rs` and
  `tests/structural_preconditioner.rs` cover scale-separated nonsymmetric
  solves with one Krylov column, independent dense endpoints, exact stage
  inputs, malformed-output/cancellation atomicity, identity parity, and
  checkpointed parameter/initial-state gradients against dense differences.
- `galpha::adjoint::FirstOrderVjp` and
  `OperatorFirstOrderGeneralizedAlpha::step_vjp` differentiate the converged
  first-order residual through an explicit transposed effective-system solve.
  The model supplies mass/tangent transpose actions and the parameter VJP of
  `M(p) rate + r(t,u,p)`, including parameter-dependent mass. Returned
  cotangents cover initial state, initial rate, model parameters and forcing;
  neither mass nor tangent must be symmetric. Chaining steps requires both
  endpoint state and rate seeds. The caller propagates the initial-rate
  cotangent through any consistency initialization, and adds direct objective
  terms and forcing-parameter chain rules. `FirstOrderAdjointConfig` bounds the
  independent adjoint FGMRES solve with an explicit preconditioner, while
  `adjoint_workspace_components` gives a conservative scalar-storage ceiling.
  The pullback borrows existing input history and returns complete primal and
  adjoint reports. Failed solves, incomplete derivative writes, nonfinite
  results, insufficient workspace and cancellation return no partial gradient.
  Stage reconstruction uses the same rate-derived increment as the forward
  map, and the correction-variable pullback avoids cancelling `1/h` terms
  in the initial-state gradient for an endpoint-rate objective.
- `galpha::second_order_adjoint::SecondOrderVjp` and
  `OperatorGeneralizedAlpha::step_vjp` differentiate the structural
  Chung–Hulbert residual and Newmark correctors. Explicit transpose actions
  cover mass, damping and the nonlinear internal-force tangent; the parameter
  VJP includes `M(p) a + C(p) v + r(q,p)`. Returned cotangents cover all three
  initial vectors (displacement, velocity, acceleration), parameters and
  forcing. A trajectory reverse sweep must carry all three endpoint seeds.
  If initial acceleration comes from an equilibrium calculation, its cotangent
  must enter that calculation's chain rule. `SecondOrderAdjointConfig` bounds
  an independently preconditioned transposed FGMRES solve. Workspaces are
  conservatively capped, prior input history is borrowed, and failed solves,
  incomplete derivatives or cancellation return no partial result.
  The reverse pass differentiates the primal's Newmark correction directly;
  it never recovers displacement/velocity cotangents by subtracting nearly
  equal acceleration cotangents and dividing by `h²`. Stage displacement and
  velocity are reconstructed from the same increments used by the forward
  residual, retaining motion that rounds away in the absolute endpoint.
  `structural_generalized_alpha_adjoint::acceleration_observation_gradients_survive_vibration_step_refinement`
  checks all state, forcing and material-parameter cotangents against an
  independent forward-dual acceleration solve over `h=1e-4..1e-10`, with
  nonsymmetric coupled mass/damping, nonlinear stiffness and three spectral
  radii (G3, `2e-10*(1+|reference|)`).
- `galpha::initialization::{first_order_rate, second_order_acceleration}` solve
  the initial mass equations `M rate = f-r(t,u)` and `M a = f-Cv-r(q)` with
  bounded FGMRES and an explicit mass preconditioner. The corresponding
  `first_order_rate_vjp` and `second_order_acceleration_vjp` solve the transposed
  mass system and propagate the rate/acceleration seed into initial state,
  velocity, model parameters and forcing. Parameter-dependent mass and damping
  enter the same residual-parameter VJPs used by the step adjoints. These
  routines provide the consistency chain needed when an inverse problem's
  initial rate or acceleration changes with its parameters. They hold time
  fixed, require a nonsingular mass model, return both residual reports, bound
  scalar workspace, and publish no partial result on cancellation or refusal.
  The caller supplies the force at the initial physical time and applies any
  remaining state/velocity/forcing parameterization. They do not establish
  uniqueness of a singular-mass model or solve constrained DAE initialization.
- `galpha::impulse::second_order_impulse` applies a prescribed instantaneous
  generalized impulse through the original matrix-free mass: `M delta-v = J`.
  The stored velocity change is independently checked component by component
  against the requested impulse; large background momentum does not enlarge
  that allowance. Signed work and kinetic increments use the mean endpoint
  velocity and double-double reductions, without subtracting large energies.
  The same initialization owner recomputes `M a+ = f+ - R(q,v+)`, including
  nonlinear velocity resistance, before velocity and acceleration publish
  together. Configuration, physical time, accepted step count and prior history
  are retained. Two bounded mass solves and O(n) scratch reuse the native mass
  preconditioner, with explicit workspace and transfer/work allowances. There
  is no timestep-dependent force surrogate, mass lumping or new integrator.
  The caller supplies independent coordinates with a constant symmetric SPD
  physical mass; solve convergence is not a proof of that assumption. Supports
  are stationary and already eliminated by the model. Collision laws, moving
  supports, singular DAEs and impulse/event adjoints remain separate; smooth
  trajectory replay must be split at the event. G1/G3/G4 tests cover full mass,
  consistent damped restart and its analytic convergence, simultaneous cross
  work, braking, lost stored increments, model refusal, capacity and cancellation.
- `galpha::second_order_adjoint::trajectory::RecordedStructural` records and
  differentiates complete structural trajectories under a
  `StructuralTrajectoryModel` with time-dependent, parameterized loading.
  Loads and their parameter VJPs use `OperatorGeneralizedAlpha::forcing_time`,
  the actual `t_n + (1-alpha_f)h` load time of the production residual. The
  recorder retains initial/current q/v/a with empty solver history and O(N)
  compact endpoint fingerprints. Binary reverse replay parks O(log N) q/v/a
  checkpoints and enforces explicit checkpoint and forward-replay limits.
  Each replay verifies q/v/a, time and the absolute step counter before using
  its derivative. Leaf pullbacks reuse the checked endpoint and Newton report,
  avoiding a duplicate primal solve, and add the forcing parameter contribution
  once. Forward caps/cancellation preserve the accepted prefix; reverse failure
  returns no partial gradient and leaves the recording retryable. Initial
  consistency derivatives remain the caller's chain rule, supported by the
  initialization routines above. The models must remain pure and unchanged;
  these fingerprints diagnose replay consistency and do not certify derivatives.
- `RecordedStructural::pullback_samples` accumulates displacement, velocity and
  acceleration observations through `trajectory::samples::StructuralSampleObjective`.
  Ordered endpoint indices include the initial state and may repeat for several
  sensors. Each callback receives the actual q/v/a, clock and absolute counter,
  writes all three state partials plus direct parameter partials, and runs once
  in reverse declaration order after endpoint verification. All observations
  share one checkpoint sweep, including initial-only zero-step objectives.
  The sample cap is checked before scanning the timetable; `3*n+p` callback
  scratch is dropped before each leaf adjoint and fits the per-step workspace.
  No observed-state history is retained. Failed, incomplete or nonfinite
  observations and cancellation publish no partial objective or gradient.
  The returned initial q/v/a cotangents still require the caller's consistency
  chain. Timetable, step, clock and solver policy are fixed; this API does not
  interpolate measurements or differentiate their sampling times.
- `galpha::{ImplicitSolveConfig, ImplicitStepTelemetry}` retain the full
  Newton report per accepted step, including outer residual decisions and
  inner Krylov counts. `SecondOrderState` and `FirstOrderState` retain time,
  solution/rate variables, accepted-step count, and complete telemetry as
  cloneable plain data.
- `stiff::Imex2::new(l, n, h)` + `imex2_step` — ARS(2,2,2),
  γ = 1 − 1/√2, diagonally implicit on L (one LU reused for both
  stages), explicit N, stiffly accurate, R(∞) = 0. Second order in
  BOTH parts: the explicit weights are (δ, 1−δ) with δ = 1 − 1/(2γ);
  trapezoidal (½, ½) is only first order in N (caught during
  construction, locked by the convergence test).
- `stiff::OperatorImex2` applies that identical ARS(2,2,2) tableau to any
  `LinearOp`. Each `(I − γhL)` stage uses FGMRES with an injected
  `FlexiblePreconditioner`; both true-residual reports are recorded in
  `ImexStepTelemetry`, and `ImexState` changes only if both stages converge.
- `stiff::adjoint::{ImexVjp, ImexStepGradient}` and
  `OperatorImex2::step_vjp` differentiate the actual ARS(2,2,2) step for
  autonomous `u' = L(p)u + N(u,p)`. Two primal and two transposed shifted
  FGMRES solves use explicit primal/adjoint preconditioners. The pullback
  includes parameter derivatives of both implicit matrices, the explicit
  linear term and both nonlinear stage evaluations. It reuses the production
  stage calculation and applies the implicit function theorem, without
  differentiating Krylov iterations. The caller supplies consistent pure
  derivative actions, a finite endpoint seed and a scalar workspace ceiling;
  incomplete/nonfinite derivative outputs, insufficient workspace, failed
  solves and cancellation return no partially accumulated gradient.
  `adjoint_workspace_components` bounds live scalar storage conservatively;
  allocator metadata and callback-owned memory are outside that bound.
  `examples/imex_thermal_fit.rs` composes the pullbacks to fit conductance and
  heater power through a two-body thermal trajectory.
- `stiff::adjoint::trajectory::RecordedImex2` records/resumes the production
  fixed-step method with O(N) compact endpoint fingerprints and only the
  initial/current full states. Its reverse sweep follows the binary
  `fs-ad::revolve` schedule, parks O(log N) states, enforces explicit checkpoint
  and replay-work caps, and checks every replayed endpoint before using its
  derivatives. Leaf pullbacks reuse the checked forward stages. Cancellation
  leaves the accepted forward prefix or complete reverse recording retryable;
  no partial gradient is returned. Models and both preconditioners must remain
  pure and unchanged across recording, forks and replay. The fingerprints are
  replay diagnostics, not proof of derivative correctness. The thermal fitting
  example now uses this sweep with six parked checkpoints for forty steps.
  `IdentityPreconditioner` is the explicit unpreconditioned fixture lane, not
  the field-scale recommendation.
- `RecordedImex2::pullback_samples` accepts a capped, nondecreasing list of
  accepted endpoint indices and the existing RK45 `SampleObjective` interface.
  Index zero denotes the initial state; repeated indices support multiple
  sensors at one endpoint. Each callback receives the actual recorded time,
  supplies the scalar loss plus state/direct-parameter partials, and runs once
  in reverse declaration order after replay verification. All terms accumulate
  in one checkpointed sweep, with no stored observation states or extra solve
  per observation. The result contains the total objective, initial/parameter
  gradients, observation count, and replay/checkpoint usage. The caller chains
  any initial-condition parameter dependence. Samples do not alter the fixed
  IMEX mesh or interpolate off-grid times. Objective failure, unwritten output,
  insufficient budgets and cancellation return no partial value or gradient.
- `stiff::ExpEuler::new(a, n, h)` + `.step(u, nonlin)` — exponential
  Euler for u′ = Au + N(u), **symmetric A** via the fs-la Jacobi
  eigenbasis. `try_new` admits finite matrices and positive finite durations,
  scales the dense eigenproblem, and checks eigenpair residuals and
  orthogonality. Opposite entries within 64 eps relative to that pair use
  their symmetric mean; larger asymmetry refuses. Unresolved matrix/rate
  scaling and unrepresentable growth refuse rather than publish a false
  decay. The source coefficient uses h·φ₁(hλ) near zero and expm1(hλ)/λ
  for stiff rates, retaining finite forced equilibria even if hλ overflows
  negatively. Split exponential factors preserve finite homogeneous tails
  when the unscaled retention is subnormal or zero. The method is exact
  for a constant additive source up to eigenbasis/roundoff error, and remains
  first order for nonlinear N; it adds no positivity or energy certificate.
  `workspace` + `step_checked` reuse three scratch vectors without per-step
  allocation. Every callback output must be written and finite. Shape,
  input, source and arithmetic refusals leave the entire accepted state
  unchanged; callback effects remain caller-owned. Construction and stepping
  are dense synchronous kernels with no internal cancellation point.
  `tests/exponential_decay.rs` checks analytic forcing, timestep partition,
  zero and extreme decay rates, large coupled spectra, retained tails, and
  transactional refusal. Existing `new`/`step` convenience calls delegate to
  these checks and panic on refusal.
- `adaptive::{AdaptiveState, PiController, rk45_adaptive}` —
  Dormand–Prince 5(4), max-norm error against atol + rtol·|u|, PI
  update h ← h·safety·err^{−0.14}·err_prev^{0.08} clamped to
  [0.2, 5.0]; rejections shrink with the classical −1/5 exponent.
- `validated::{validated_step, ValidatedStep, ValidatedOptions}` supplies
  a smooth explicit-ODE foundation for extension-plan §3.4
  (`frankensim-ext-time-validated-step-ow2o`). Outward interval arithmetic
  checks strict Picard inclusion `X0 + [0,h] F([t0,t1],Y) ⊂ interior(Y)`.
  Accepted objects retain the initial uncertainty, whole-step enclosure,
  endpoint enclosure and derivative box; subwindow queries reuse that proof.
  The caller must supply an honest interval extension of one deterministic
  vector field that is continuous and locally Lipschitz in state on an open
  neighborhood of the tested box. These mathematical callback properties
  are assumptions, not conclusions inferred from finite samples.
  Dimension, width, iteration, callback and component-work limits are explicit;
  ambient execution budgets and cancellation are checked at bounded tiles and
  around callbacks. Callbacks own their internal work and cancellation. Invalid
  or unwritten derivatives, overflow, unsuccessful inclusion and exhausted
  budgets publish no step. Four state-sized interval vectors are retained;
  five are live during construction. The interval elementary-function
  assumptions remain those of `fs-ivl`. This first-order enclosure can wrap
  badly; a refusal does not prove that the ODE itself is ill-posed.
- `validated::trajectory::{validated_trajectory, TrajectoryTube}` composes
  those conditional enclosures over a finite forward horizon. Each accepted
  endpoint is the next initial set verbatim; no midpoint replacement discards
  uncertainty. Failed Picard inclusion or excessive width triggers deterministic
  midpoint subdivision, subject to the minimum step. Callback refusals,
  cancellation and arithmetic failures stop immediately. Accepted step sizes
  are retained without heuristic growth, and the final time is clamped exactly
  to the requested horizon. A zero-duration request checks its initial set
  without invoking the RHS.
  One ambient accountant and whole-trajectory attempt, RHS and component-work
  caps cover initial retention, all accepted steps and all refused attempts.
  Per-step limits remain additional bounds. The step count and dimension caps
  bound retained storage, with fallible allocation. Publication occurs only
  after each step's final checkpoint. `TrajectoryFailure` exposes the completed
  contiguous prefix and recorded total consumption, never a tube for the
  uncomputed suffix; an empty prefix covers only its initial time. The caller
  must keep the same deterministic RHS and its interval-extension assumptions
  throughout. G0/G1/G3/G4/G5 tests use analytic flows, subdivision, exact budget
  edges, callback/cancellation refusal and deterministic repeated execution.
- `validated::events::scan_guard` encloses a caller-supplied scalar guard on
  every closed time cell of a `TrajectoryTube`. A strict interval sign excludes
  the entire cell; it is never inferred from a finite set of samples. Cells
  whose guard contains zero are subdivided deterministically until the requested
  width or floating-point time resolution is reached. Exact `[0, 0]` ranges
  remain one unresolved cell without futile subdivision. Uncertain cells are
  retained as `Unresolved`, including tangencies, clustered roots and shared
  step boundaries. A zero-duration tube checks its initial point explicitly.
  An optional strictly signed enclosure of the total derivative
  `g_t + g_x f`, together with strict opposite endpoint signs, proves a unique
  crossing in that cell for each admitted trajectory, conditional on the
  caller's continuous differentiability and interval-extension obligations.
  Monotonicity and matching strict endpoint signs can instead exclude a cell.
  `no_events_proved()` is true only after the complete covered interval was
  excluded. Refusal, cancellation or exhausted cell, callback, frontier, output
  or work budgets returns no partial successful scan. One ambient accountant
  covers every range query and scan operation. Windows may share endpoints;
  their number is not a root count. G0/G1/G3/G4/G5 regressions cover analytic
  crossings, missed-by-sampling roots, uncertainty, budget edges and repetition.
- `se3` module (bead 3ol0): canonical `fs-ga::Twist` uses
  `[angular, linear]` ordering. `se3_exp_step(T, twist_body, h)` performs the
  right/body update `T Exp(h twist_body)` and
  `se3_space_exp_step(T, twist_space, h)` performs the left/space update
  `Exp(h twist_space) T`; both delegate group exp/composition and double-cover
  canonicalization to the sole `fs-ga::Se3` owner. `se3_exp_step_renorm` adds
  drift-controlled versor renormalization gated by `RenormPolicy` and
  returns a `RenormReceipt` (defect before, whether renormalized,
  reported drift) — drift is ledger fodder, never silently absorbed.
  `se3_rigid_body_step` uses the same free-rotation splitting and advances the
  center of mass with its constant spatial velocity. The pose origin must
  be the center of mass and inertia is diagonal in the body principal frame.
  Body linear velocity is transported through the accepted attitude, keeping
  linear and total angular momentum to roundoff instead of integrating a
  separate drifting transport ODE. `dep_free_step` — discrete
  Euler–Poincaré free-body step in body-momentum form (fixed-point
  solve controlled by `DepSolveParams`, per-step `DepStepReceipt`);
  spatial angular momentum is conserved EXACTLY by construction.
  `run_dep_free` returns a `BalanceReceipt` whose `Se3ClaimClass` is
  decided by `claim_for(declaration, all_solves_converged)`: the
  conservative variational theorem class only for declared smooth,
  conservative, regular-constraint, fixed-step fixtures with every
  solve converged — dissipation, adaptivity, or divergence demote to
  `MeasuredOnly` with measured drift receipts. `dep_momentum_adjoint`
  pulls a terminal cotangent through the transposed implicit-function
  tangent of each step's ACTUAL residual. Both 3-by-3 residual Jacobians are
  analytic, using `fs-ga`'s right SO(3) Jacobian and the same inverse momentum
  transport as the primal. One exponential replaces twelve perturbed residual
  evaluations per reverse step; no finite-difference step size is selected.
  `RattleProjection` is the
  constraint hook fs-mbd's constrained lanes plug into
  (`Unconstrained` is the trivial impl).

- `slabs` module (addendum Proposal 4, bead bk0o.7; [F], behind
  `time-slabs`): time slabs as CELLS. `SlabEntry` carries the TEMPORAL
  COCYCLE — the split step's defect against the monolithic residual
  over the slab — and `SlabLedger::attribute` is the budget pie
  pointed at time ("your error is in the coupling handoff at
  t ∈ [2.0, 2.5]"). `march_adaptive` doubles subcycles where the
  cocycle exceeds tolerance (cap 64). `activation_report` encodes the
  Proposal-4 SEQUENCING gate: splitting error under 20% of budget →
  INSTRUMENT-ONLY (measure it, don't control it). Slab constructors
  and marchers fail fast on zero substeps/slabs, invalid time spans,
  non-finite states/couplings, and non-positive total error budgets;
  ledger reporting rejects hand-built malformed entries before
  attribution or JSON emission.

- `hybrid` module (RE.Z1; `[M]` object semantics, no theorem promotion):
  `ZenoProblemIrV1` is the finite, versioned object language for a hybrid
  model and immutable model version, original-versus-regularized lineage,
  initial-state set, frame and unit conventions, physical time scale, finite
  modes and continuous ODE/inclusion/admitted-DAE dynamics, oriented guards,
  transverse/grazing/Unknown crossing state, deterministic/set-valued/
  terminal/Unknown reset relations, contact/relay laws, simultaneous-event
  groups and policy, dwell semantics, event language, hybrid-time domain,
  compactness, numerical accumulation candidates, continuation category, and
  semantic-analysis budget. All artifact roles use non-confusable typed IDs.
  Unsupported infinite-dimensional dynamics, zero DAE index, invalid time or
  budget data, missing mode references, malformed simultaneous priorities,
  and self-referential regularization refuse.
- `validate_zeno_problem_v1(ir, cx)` caps modes/events before any nested reset
  scan, then aggregates bounded reset-target counts with permutation-stable
  diagnostics before sorting or graph work. It polls cancellation at bounded
  strides, canonicalizes set order, and mints a domain-separated
  `ZenoProblemIdV1` plus independent canonical-preimage receipt only on
  success. Its conservative finite graph analysis records zero/Unknown-dwell
  cycles and local nonuniqueness. A unique
  continuation category refuses when differential inclusions, grazing or
  unresolved guards, set-valued or unresolved resets, unresolved simultaneous
  ordering make it an overclaim. A zero-time cycle is recorded separately: it
  enables an accumulation question but does not by itself select a unique or
  nonunique post-Zeno rule.
- `ZenoClassificationV1` keeps five states disjoint: positive finite-event
  separation, theorem-shaped certified-Zeno descriptor, numerical
  event-density warning, a distinct regularized-model lineage, and Unknown.
  `PostZenoStateV1` independently records not-applicable, unique, set-valued,
  terminal, or unresolved post-accumulation semantics.
  `validate_zeno_claim_descriptor_v1` cross-checks the exact problem and
  continuation category. Event caps and numerical-only traces cannot enter a
  theorem-shaped classification; finite separation requires a positive dwell
  bound on every admitted event; certified Zeno requires a zero/Unknown-dwell
  accumulation cycle in this finite automaton; regularized claims require the exact second
  validated problem and matching source lineage. Claim and post-state choices
  participate in `ZenoClaimIdV1`.
- A successful `ValidatedZenoClaimDescriptorV1` always returns
  `ScientificCorrectnessNotProven`. Opaque analytic/interval witness IDs are
  replay references, not authority tokens. RE.Z2 must independently check and
  promote a descriptor before any Zeno, finite-separation, or continuation
  theorem exists.

## Invariants

- Verlet is symplectic ⇒ **bounded** (non-secular) energy error;
  its KDK positions satisfy the discrete Euler–Lagrange recurrence
  q_{k+1} = 2q_k − q_{k−1} + h²F(q_k) (Marsden–West).
- Exp-map attitude updates keep ‖q‖ = 1 by construction — no
  renormalization in the SO(3) lane. The SE(3) lane exposes optional,
  threshold-controlled motor renormalization only through
  `se3_exp_step_renorm`, with every decision and drift returned in a receipt.
- Body/right and space/left steps use named entry points and canonical
  `fs-ga` group operations; no array multiplication order implicitly chooses
  a frame convention.
- Generalized-α: high-frequency per-step contraction → ρ∞; order 2
  across the whole ρ∞ range.
- Dense and operator-backed generalized-α/IMEX paths use the same parameter
  formulas, intermediate-time convention, correctors, and ARS weights.
  Equality across those arithmetic paths is to the declared nonlinear/Krylov
  tolerance; it is not claimed bitwise.
- Implicit operator steps are transaction-like: a failed nonlinear or shifted
  linear solve leaves the public trajectory state unchanged. Every accepted
  step appends exactly one telemetry row before returning.
- `AdaptiveState` is COMPLETE: checkpoint = `clone()`; split runs are
  **bitwise-equal** to straight runs, controller memory (`err_prev`)
  and counters included (P7). A step shortened to land exactly on
  t_end does not feed the controller (the clamped h would poison the
  h carried into a later resumed segment).
- Hybrid identities are invariant to input ordering of mode/event sets and
  sensitive to every model version, initial set, frame/unit/time convention,
  guard/reset/law, event order, regularization, budget, claim, and no-claim
  artifact. An execution event cap is byte-visible policy, never proof.

## Error model

Construction panics with a structured message when a required
factorization is singular (`GeneralizedAlpha::new`,
`FirstOrderGeneralizedAlpha::new`, `Imex2::new`):
these are modeling errors (h incompatible with the operators), not
runtime conditions. Operator-backed methods return typed dimension,
non-finite explicit-stage, Newton-setup, and not-converged refusals; their
reports retain the exhausted budget or breakdown diagnosis. `verlet_adjoint`
inherits `fs_ad::revolve`'s budget assertion. Legacy dense steppers themselves
are panic-free on finite input; NaN/Inf propagate as NaN/Inf
(garbage-in, garbage-out, never UB). Lie-group lanes return `LieStepError` or
`Se3Error`: canonical representation failures preserve the structured
`fs_ga::GaError`, while invalid step, inertia, damping, renormalization, and
variational-solver controls carry an explicit parameter-family context.
`DepSolveParams` requires finite positive tolerance and a nonzero iteration
budget even for an empty requested trajectory; per-step damping is a fraction
in `[0, 1]`.

## Determinism class

Bit-deterministic across ISAs BY CONSTRUCTION: no platform libm in any
solver path (workspace law) — transcendentals go through
`fs_math::det` (`sin`, `cos`, `sqrt`, `exp`, `expm1`, `pow`), including
the PI controller's `pow`, because the adaptive step SEQUENCE is part
of the contract. Test-side oracles may use std (disjoint-path rule).
Golden FNV-64 over Verlet, rigid-body, generalized-α, IMEX, ExpEuler
and RK45 trajectories (controller state and counters included):
`0x41c7_740c_1a69_00dd`, recorded on Apple M4 Pro (aarch64) after migrating
the rigid-body representative to canonical `fs-ga::So3`. The predecessor
golden was verified on Threadripper (x86_64); this new canonical-state golden
still needs retained x86_64 re-verification and makes no fresh empirical
cross-ISA claim. The operator-backed paths use
`fs-solver`'s deterministic reductions and logical-iteration preconditioner
contract. Their dense-vs-operator and split-run fixtures are deterministic on
the exercised build; they do not yet add a retained cross-ISA golden.
Hybrid problem and claim identities use bounded `fs-blake3` canonical schemas.
Canonical sets sort by typed identity, ordered priorities remain ordered, and
signed zero is normalized before hashing. Identical admitted input and schema
version replay byte-for-byte on every ISA; no numerical trajectory participates
in identity unless named by its exact trace ID.

## Cancellation behavior

All entry points are synchronous. Single-step owners run one bounded step;
the validated trajectory owner additionally enforces whole-run step, attempt,
RHS and work caps and preserves only completed steps on refusal.
Operator generalized-alpha bounds work by the configured Newton outer and
FGMRES cycle/restart budgets; operator IMEX bounds each stage by its configured
FGMRES budget. `OperatorImex2::step_controlled` and `step_vjp` additionally
poll around model callbacks, between FGMRES restart cycles and immediately
before publication. One cycle has at most the configured restart length of
Arnoldi iterations. Long callbacks must bound their own work; this does not
claim interruption within an operator or one Krylov cycle. Cancelled forward
attempts leave the entire `ImexState` unchanged and are retryable.
Both `OperatorFirstOrderGeneralizedAlpha::step_controlled` and structural
`OperatorGeneralizedAlpha::step_controlled` pass their existing mutable
cancellation predicate into the shared Newton solver. They poll before
setup, before/after each Newton residual, Jacobian and preconditioner
callback, inside Arnoldi and backtracking, and before publication. A true
predicate result interrupts the current Newton attempt without continuing
its remaining numerical callbacks. Numerical callback unwinds become
`TimeSolveError::SolverCallbackPanicked`; cancellation still returns
`TimeSolveError::Cancelled`. Both preserve the entire input state/history,
including an already accepted prefix, and successful retry follows the
same numerical trajectory. The first- and second-order preconditioner
batteries cancel from real Jacobian, preconditioner, and trial-residual
callbacks and verify callback counts, complete rollback, and bitwise retry.

An individual shared residual/Jacobian callback may compose several model
actions and must finish before its cancellation is observed; vector
reductions are also not preempted. The constructor's initial residual
remains bounded by the caller's setup checks, outside numerical callback
panic containment. Predicates must be bounded and nonpanicking. The
`step_vjp` methods retain their separate checks between adjoint restart
cycles and around final derivative callbacks; reverse FGMRES is not yet
interruptible inside a restart cycle through this driver.
Long trajectories are resumable by cloning `SecondOrderState`,
`FirstOrderState`, `ImexState`, or `AdaptiveState` between calls; split runs
continue bitwise when the same operators, forcing, preconditioner policy, and
configuration are supplied. Hybrid problem validation polls `Cx` before work,
every 32 modes/events, during zero-time graph analysis, and inside canonical
identity streaming; cancellation publishes no token. No async, no internal
threading, no I/O.

## Unsafe boundary

None. `unsafe_code = "deny"` via workspace lints; no capsules.

## Feature flags

- `time-slabs` — [F] time slabs as temporal cells and splitting-error
  ledger/controller instrumentation. No optional dependencies.

## Conformance tests

`tests/frankenscipy_ode_oracle_casebook.rs` adds a structured endpoint-oracle
subset for `rk45_adaptive`. Canonical frames bind the complete stationary,
decay, and oscillator fixtures; every production PI-controller option and work
cap; every FrankenScipy `solve_ivp` RK45 option; the authoritative constellation
pin; analytic-reference formulas; versions; and numerical ceilings. Complete
receipts retain every public `AdaptiveState` and `SolveIvpResult` field, both
implementations' endpoint bits, independent closed-form endpoint bits, signed
component errors, max-norm deltas, and remaining margins. Each measurement must
replay byte-for-byte within one executing build. A disclosed seeded bit mutation
of the exact stationary oracle reference must produce one stable red record
that `assert_green` refuses without accepting baseline drift. This finite,
smooth, forward, nonstiff endpoint subset is G0 comparison evidence only:
FrankenScipy is a pinned independent implementation, not ground truth, and the
test-side standard-library analytic functions are not a cross-ISA bit oracle.
The subset makes no claim for general tolerance calibration, controller-path or
counter equality, dense output, events, backward integration, stiffness,
checkpoint equivalence, adjoints, cancellation, performance, Python SciPy, or
fresh cross-ISA/full-G5 execution.

`tests/time_battery.rs` (16 cases, JSON logging): 10⁶-step harmonic
energy boundedness vs RK4's secular decay at the same h + e = 0.6
Kepler orbit (~16 revolutions); discrete Euler–Lagrange residual
≤ 1e−12 on a nonlinear potential; canonical SO(3) quaternion norm drift
≈ 1e−12 over
10⁵ steps with no renormalization; gyroscope battery (ω₃ constant to
1e−9, analytic precession phase Ω = (I₃−I₁)/I₁·ω₃ to 1e−3, energy and
spatial angular momentum to O(h²)); generalized-α spectral radius → ρ∞
at ωh = 10³ to 2% with ρ∞ = 0 annihilation and ρ∞ = 1 energy
preservation; IMEX hλ = −100 monotone contraction + measured order 2
on the logistic equation; ExpEuler exact (1e−13) on a linear system
against a disjoint-path oracle + measured ETD1 order 1; RK45 accuracy
tracking rtol, rejection recovery from an absurd h₀, and bitwise
split-run resumability at 4 cut points; Verlet adjoint gradcheck vs
central FD ≤ 1e−7 relative; the retained deterministic golden hash (fresh
x86_64 replay pending after the canonical-state migration).
`tests/galpha_probe.rs`: order-2 sweep across ρ∞ ∈ {0, .3, .5, .8, 1}
at two horizons — kept from the probe that diagnosed the period-point
metric blindness (at t = 2π, cos′ = 0: a q-only error measures phase
error quadratically and fakes order ≈ 4; the honest metric is
max(q, v) error).
`tests/operator_integrators.rs` (bead ra2q; G1/G3, JSON-line measurements):
second-order structural convergence on prefactored dense and diagonal
matrix-free `LinearOp`
paths; dense/operator structural agreement plus bitwise split replay with
Newton/Krylov telemetry; first-order-system order and dense/operator
agreement plus bitwise split replay; the normalized small-Biot lumped thermal
decay at one time constant bound directly to the canonical `fs-vvreg` Level-A
row; nonlinear first-order self-convergence;
constant acceleration under an unresolved increment of a large absolute
velocity in both first-order lanes;
and ARS(2,2,2)
dense/operator agreement, nonlinear logistic order, split replay, and both
stage iteration counts.
`tests/slabs.rs` under `time-slabs`: temporal cocycle consistency,
budget-pie localization, adaptive-vs-uniform cost, G3 repartition
envelope, activation gate, and fail-fast validation for invalid slab
counts/substeps/tolerances/budgets, malformed public ledger entries,
and non-finite couplings.
`tests/se3.rs` (bead 3ol0, printed measurements on numerical gates):
SO(3)-lane agreement for pure rotations (se3-001); constant-twist
one-parameter composition exactness (se3-002); double-cover
canonicalization determinism, `M`/`−M` bit-equality, and bitwise
replay across a scalar-zero crossing (se3-003); free-body DEP spatial
momentum to roundoff plus bounded energy over 10⁴ steps with the
theorem claim class earned (se3-004); adjoint-vs-central-FD gate on
the DEP momentum map, four directions, rel ≤ 1e−6 (se3-005); invariance under
inverse time/angular-velocity rescaling from 1e-12 to 1e12; the analytic
spherical-body identity pullback across inertia scales; the
honesty fixture — a damped run demotes to `MeasuredOnly` with nonzero
measured drift despite converged solves (se3-006); 10⁵-step
renormalization receipts bounding the final unit defect (se3-007);
SE(3) rigid-body agreement with the SO(3) lane and spatial
free-velocity conservation to roundoff (se3-008); explicit
body/right versus adjoint-related space/left equivalence for both groups plus
bitwise replay of all four stepping lanes (se3-009); invalid solver controls
and out-of-domain damping refuse even when a zero-step horizon could otherwise
bypass the step loop (se3-010). G0/G1/G3 splitting regressions cover a
20,000-step asymmetric top's spatial momentum and bounded energy, exact
axisymmetric precession for each principal axis, second-order refinement and
time reversal, and a translating/rotating body's linear momentum, total
angular momentum, and analytic center-of-mass trajectory.
Held-torque G0/G1/G3 regressions additionally cover exact constant angular
acceleration about each principal axis, the accumulated spatial angular
impulse over 2,000 steps, second-order refinement and time reversal for
both torque conventions, zero-load parity, and non-finite-load refusals.
`tests/lie_compensated.rs` additionally checks signed sub-ULP angular kicks
whose missing midpoint phase changes attitude by a full radian, a full
off-diagonal tensor, independent Euler/RK4 second-order refinement and
first-variation responses in both attitude and nonlinear spin, DD momentum
norm and spatial momentum preservation, time reversal, finite large phases
and subnormal durations, and the actual held-torque public owners' distinct
body/spatial endpoint impulse balances.
`tests/hybrid.rs` (RE.Z1; G0/G3/G4/G5): relay-chatter and positive-dwell
microcases; mode/event input-order and schema replay; model-version identity
sensitivity; event-cap and numerical-trace theorem falsifiers; finite
separation versus zero/Unknown-dwell cycles; simultaneous guards, grazing, and
set-valued-reset uniqueness refusals; exact regularization lineage; numerical
warning post-state honesty; unknown-mode, zero-index DAE, out-of-domain window,
unsupported-version, resource-cap, and pre-cancellation refusals.
Resource-cap coverage includes top-level mode/event/priority fail-fast behavior
before nested reset inspection, plus permutation-stable diagnostics for both
per-relation and aggregate reset-target caps.

## No-claim boundaries

- Canonical `fs-ga::Se3` states ship the exp-map lane and the FREE-BODY
  variational lane only: no forced/damped/constrained variational
  steps here (constrained lanes are fs-mbd's, through
  `RattleProjection`); the backward-error/modified-energy THEOREM is
  claimed only for the declared smooth conservative fixed-step class
  with converged solves — impacts, adaptivity, and dissipation get
  measured `BalanceReceipt`s, never the theorem. The
  `dep_momentum_adjoint` differentiates the converged implicit residual using
  analytic 3-by-3 tangents, not the finite fixed-point iteration history.
  Its agreement with the accepted forward map remains limited by the primal
  solve tolerance. Revolve checkpointing (O(N) memory today) remains follow-up
  work. The variational lane's group update solves the
  midpoint transported-momentum fixed point; no claim of equivalence
  to any OTHER published DEP variant. The DEP lane has no higher-order
  compositions or splitting beyond kick–drift–kick; the separate compensated
  free-rotation lane explicitly provides the fourth-order composition above.
- No Krylov φ-actions for large nonsymmetric A (needs Arnoldi);
  `ExpEuler` is dense-symmetric only. No Rosenbrock/SDIRK families,
  no BDF/multistep.
- RK45 dense output and event brackets describe the numerical interpolant.
  The separate validated explicit-ODE owner encloses the true flow under its
  callback assumptions; it does not upgrade RK45 brackets into true-flow
  certificates. Its guard scanner provides conditional whole-interval exclusion
  and isolated monotone crossings, retaining unresolved windows explicitly.
  No complete root count, grazing/simultaneous reset proof, constrained DAE
  tube, stiffness detection or Zeno certificate is provided by that owner.
- Adjoint support is specific to each admitted map: Verlet, free-body DEP,
  frozen-step/recorded RK45, autonomous operator IMEX ARS(2,2,2), and
  operator-backed first-order and structural second-order generalized-alpha.
  Exponential-integrator pullbacks remain absent.
  Generalized-alpha holds time, timestep, spectral radius and solver policy
  fixed. It differentiates the converged residual rather than Newton/Krylov
  iterations. The first-order tests check independent state/rate input seeds
  with nonlinear nonsymmetric, parameter-dependent mass and tangent
  against an independent dense endpoint-rate solve with forward duals and
  five-point differences, plus full trajectory accumulation and descent.
  Structural tests independently solve for endpoint acceleration with dense
  forward duals and five-point differences, check parameter-dependent mass,
  damping and internal force, and verify all displacement/velocity/acceleration
  and forcing cotangents through multiple steps and an inverse-loss descent.
  Primal and adjoint convergence residuals limit gradient accuracy; there is
  no interval gradient enclosure or adaptive-controller/event derivative.
  IMEX gradients hold the step size, solver policy and initial-condition
  parameterization fixed, require correct transposed model actions, and are
  limited by the returned primal/adjoint residuals. They do not differentiate
  adaptive controllers, events/resets or nonautonomous splittings and do not
  provide an interval enclosure of the gradient. IMEX binary checkpoint
  scheduling is provided, while initial-condition/direct-objective derivatives
  remain caller responsibilities. It does not provide disk spill or serialization
  of a recording, and reverse sweeps restart after cancellation.
  `tests/imex_adjoint.rs` checks nonsymmetric state/parameter
  gradients against independent dense-solve differences, reverse trajectory
  accumulation, descent, every cancellation boundary and derivative refusals.
  The recording regressions additionally compare checkpointed/full-storage
  gradients bit-for-bit, check forward cap/clone/resume behavior, reject changed
  replay physics, and exercise exact replay/checkpoint budgets.
  `tests/imex_samples.rs` compares sampled losses and all total parameter/state
  derivatives against dense five-point differences, checks exact full-storage
  reverse equality, repeated sensors and actual times, and exercises zero-step
  initial terms, all cancellation boundaries and observation/budget refusals.
- Operator-backed generalized-alpha and IMEX remove the dense storage/API
  ceiling, but no roofline or field-scale iteration-count claim is made.
  Both generalized-alpha formulations use model-provided flexible right
  preconditioning when supplied, including during forward checkpoint replay.
  `OperatorImex2` accepts an independently injected flexible
  preconditioner. Dense-vs-operator agreement is tolerance-based, not a claim
  that LU and Krylov execute identical floating-point reductions.
- Public trajectory checkpoints are accepted-step boundaries. Although the
  shared Newton/FGMRES engines are internally resumable, this fs-time API does
  not yet expose a mid-Newton or mid-IMEX-stage checkpoint artifact.
- First-order generalized-alpha requires the caller's initial `rate` to be
  consistent with its modeled residual and the supplied forcing. The crate
  records but does not independently certify that model consistency. The
  nonlinear `SecondOrderProblem`/`FirstOrderProblem` implementer likewise owns
  tangent consistency; this lane does not yet run an automatic JVP audit.
- The Level-A lumped thermal binding exercises only the normalized scalar ODE
  `d(theta/theta0)/d(t/tau) + theta/theta0 = 0` at `t/tau = 1` and checks the
  catalog's declared `Bi <= 0.1` applicability context. It does not model heat
  capacity, geometry, spatial gradients, or establish that an arbitrary body
  is in the lumped regime. Its numerical envelope is distinct from the
  catalog formula-reproduction tolerance, and no ladder or machine fingerprint
  is persisted into `fs-vvreg`.
- RE.Z1 validates only finite hybrid object semantics and content identity. It
  does not locate events, enclose accumulation times, prove transversality,
  prove compactness, establish finite separation, certify Zeno behavior, or
  choose physically correct post-Zeno dynamics. `CertifiedZeno`, analytic, and
  interval-validated are descriptor labels bound to opaque evidence references;
  every returned descriptor says `ScientificCorrectnessNotProven` until the
  independent RE.Z2 checker promotes it. A compliant/penalty regularization is
  a distinct model with an explicit no-equivalence artifact, never evidence
  about the discontinuous source model.

## No-claim boundaries (slabs)

- Visibility and control of SPLITTING error only — NO claim of
  coupling STABILITY: added-mass FSI instabilities and stiff
  time-scale pathologies are per-coupling analysis problems (the
  proposal's own honest scope, verbatim).
- The monolithic reference is a fine RK4 over the slab — a numerical
  reference, not an exact solution; defects below its own error floor
  are not resolvable.
- The fixture family is the linear two-field testbed; PDE couplings
  ride fs-couple's port-Hamiltonian interconnection when both land.
- Parallel-in-time coarse propagators (the BDDC-pattern extension) are
  explicitly deferred.

## Accepted dense paths for physical applicability and work observers

The original competing-event RK45 owner now offers
`rk45_until_any_event_observed`. It exposes each accepted dense step together
with the end of its accepted prefix before committing that step. On an event,
the discarded continuation is excluded. Observer refusal or subsequent
cancellation preserves the last accepted state; observers stage their own
candidate data until successful completion. The original entry point delegates
with an empty observer and retains the same trajectory and RHS work.

`DenseStep::component_range` queries a component over the accepted prefix using
the quartic's Bernstein coefficient hull, accumulated with double-double
arithmetic and including actual accepted endpoints. This includes internal
extrema without sampling, another ODE evaluation, or a new root solver. It is
a floating numerical hull, potentially wider than the interpolant's actual
range, and is not an interval certificate for arithmetic or exact physical
trajectories. It enables speed-domain checks for nonlinear friction while the
original integrator remains the numerical owner.

Four focused tests cover an internal extremum, prefix restriction and invalid
queries, exact constant ranges, identical original/observed trajectories and RHS
counts, late observer refusal, and cancellation after observation without
committing the in-flight step.
