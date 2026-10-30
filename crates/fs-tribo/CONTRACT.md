# CONTRACT: fs-tribo

> Status: `frankensim-ext-tribo-dry-baseline-tgbj` baseline. This is not an
> admitted material-data result or an Euler-disc validation result.

## Purpose and layer

`fs-tribo` is an L3 dry-contact constitutive leaf. It provides typed dry
friction, elastic Hertz reference forms, scalar rolling/contour loss laws,
explicit dissipated-work accounting, a bounded flash-temperature candidate,
and caller-owned Archard wear state. It does not depend on contact-solver
internals: the consuming solver owns contact complementarity, tangent
reactions, and residual embedding.

## Public types and semantics

- `InterfaceSystemRef` names ordered surfaces, history, medium, and a nonblank
  caller `source_id`. `InputAuthority` is a caller-declared ceiling
  (`CallerDeclared`, `SyntheticFixture`, or `Estimated`), never an admission or
  verification receipt. Every constitutive response returns this provenance.
- `DryInterfaceSystemCard` is a dependency-independent, immutable
  InterfaceSystemCard-style input. It combines that ordered interface/history
  identity with one `FrictionLaw` and closed SI temperature, nominal-pressure,
  and slip-speed ranges. Its `query` refuses all out-of-domain states and a
  mismatch between the declared slip-speed magnitude and the supplied tangent
  vector. It is not an alias for, or a substitute for, the pending seed-data
  `fs-matdb::InterfaceSystemCard` query/receipt path.
- `ContactFrame` normalizes a finite contact normal. `TangentialSlip` accepts
  only a finite velocity with no material normal component; it refuses rather
  than silently projects a closing/separating velocity into friction.
- `FrictionLaw` provides Coulomb, velocity-dependent, and Stribeck rungs. At
  zero slip `evaluate` reports static capacity only; it never invents a stick
  reaction. `kinetic_coefficient` exposes `μ(|v|)` for a strictly positive
  speed. `regularized_traction_1d` is the declared ODE embedding: a linear
  ramp through the origin to the kinetic value at scale `v_reg`, then the kinetic rung,
  with force opposing `v_rel`. That map is a solver regularization (the 1-D
  analogue of the fs-contact smooth adapter), not complementarity. A bow, a
  brake pad, and a fault are the same law.
- `TangentialLaw` is the stateful, solver-independent pHS constitutive seam.
  `FrictionLaw` and `DryInterfaceSystemCard` implement it with an implicit static
  reaction followed by deterministic sliding bisection, never an ODE ramp.
  Input motion is the released driver-minus-body increment and the explicit
  isotropic body compliance gives `du = du_free - C F`. Returned traction is
  effort entering the interface, conjugate to the **solved** increment, so
  `Delta H + D - F dot du = 0`. Dry laws have `H = 0`, export `D >= 0` through
  `DissipationStep`/`WorkLedger`, and retain accumulated slip plus the optional
  caller thermal boundary. They do not evolve temperature or a third-body layer.
  `TangentialSolveConfig` requires explicit force/residual/iteration budgets
  (1 through 128 bisections); unresolved or unbracketed roots refuse. Candidate
  state is immutable until the consuming solver publishes it. The card impl
  requires temperature and pressure without defaults, and checks the solved
  speed against its closed applicability interval. The matrix adapter below
  admits anisotropic body response; the mechanics owner supplies that matrix.
  Thermal transport and pHS network wrappers retain their respective owners.
- `layer_contact::ThirdBodyLayerContact` implements the same seam with retained
  elastic and plastic shear. The midpoint radial return solves the full
  isotropic body load line in two tangent axes. Reversible energy is
  `H = k |e|²/2`; physical dissipation is `Y |dp| + c |dp|²/dt`, with no
  numerical energy deficit relabeled as heat. `Sticking` means no plastic flow;
  elastic relative displacement may remain nonzero. Every step requires
  positive compressive load, the original frozen temperature, the caller's
  explicit elastic-slip bound and force budget. Opening and foreign-law
  checkpoints refuse without changing state. The numerical law has caller
  authority; the separate fs-material source owner binds actual layer geometry
  and the distinct plastic-shear-viscosity claim. Thermal softening, layer
  transfer/wear, recontact and measured rosin behavior are outside this rung.
- `HertzSpherePlane` and `HertzCylinderPlane` are G1 elastic closed forms,
  expressed in SI. They use caller-supplied radius/modulus and retain the
  caller ceiling; a mathematical closed form does not validate its inputs.
- `ConstantRollingMoment` uses `M = -sign(omega) N a` (N m) and
  `P = N a |omega|` (W). `ConstantContourForce` uses
  `F = -sign(v) F_c` (N) and `P = F_c |v|` (W). They are separate mechanisms;
  callers must not combine them without independent source support.
- `HeatPartition`, `DissipationStep`, `WorkLedger`, and `WearState` have
  private invariant fields and checked getters. `WorkLedger` retains total
  dissipated work and every explicit heat/other channel, and commits all four
  candidate totals only after their closure checks pass. `ArchardLaw` computes
  `dV = k F ds/H` in m3 and commits only a validated candidate state.
- `flash_temperature_candidate` is a typed model-form branch. With complete
  caller thermal conductivity/diffusivity data for every nonzero heat share,
  it reports the uniform-flux semi-infinite candidate
  `2 q'' sqrt(alpha * t / pi) / k`, with `t = traverse_length/slip_speed`.
  Missing properties return a typed `Unknown`, rather than guessed material
  data or a partial temperature result.
- `surface_excitation` admits named, uniformly sampled finite or periodic
  surface-height traces and filters them by the exact one-dimensional marginal
  of an elliptic Hertz pressure footprint. It analytically integrates the
  piecewise-linear trace and its path derivative, then applies only the
  explicitly bounded contact tangent `delta_F = k_n delta_h`.
  `SelfAffinePeriodicProfileSpectrum` can construct a deterministic,
  random-phase, band-limited one-dimensional profile whose spatial PSD follows
  `Phi(q) proportional to q^-(2 H + 1)` and whose Fourier coefficients are
  normalized to an explicit RMS height. Its seed selects spatial phase only;
  it creates no temporal/audio spectrum or apparatus-specific frequency.

### Directional static and kinetic Coulomb friction

`tangential::elliptic::EllipticCoulomb` implements the shared `TangentialLaw`
with explicit principal static/kinetic coefficients and a retained material
texture direction. A feasible zero-slip reaction is tested against the static
ellipse before solving sliding on the kinetic ellipse. Sliding follows the
associated maximum-dissipation law, so traction need not be parallel to slip;
independent axis clamps or a scalar coefficient sampled by direction are not
substitutes. Zero principal coefficients retain frictionless directions.

The scalar body load line reduces to a bounded monotone ellipse projection.
Circular ellipses and principal-axis motion use an exact noniterative path.
The remaining projection uses fixed-size arithmetic with no solver allocations, honors the original
force and iteration budgets, and exports actual `F dot du` as heat with zero
elastic storage. The existing full-matrix and simultaneous-contact adapters
consume it directly, including off-diagonal and cross-contact body recoil.
Objective history transport rotates both the material direction and accumulated
slip without time advancement, heat, or source replacement. Foreign native
checkpoints and unresolved roots refuse before publication.

Explicit coefficient continuation validates both new principal pairs while
preserving accepted slip, texture and the last thermal boundary. Its zero
elastic storage needs no parameter work; the source owner remains responsible
for authorizing coefficient changes. Physical dissipation is evaluated in the
principal frame, where its terms are nonnegative. Rotation roundoff in the
world-frame work is retained in the closure residual, so a valid zero-work
frictionless direction is not rejected as negative dissipation.

`advance_kinetic` continues an interval whose mechanical owner already admits
sliding, using the kinetic ellipse without spending the larger static capacity.
The original full law owns the returned checkpoint. A zero-slip kinetic
proposal requires a separate physical arrest check in a dynamics consumer;
the constitutive return map alone does not establish an event time.

The nine G1/G3 groups in `tests/elliptic_coulomb.rs` compare the response with
an independent constrained-energy minimization and an angular maximum-work
sweep, and exercise directional breakaway, unequal static/kinetic limits,
isotropic and frictionless limits, full matrix/multicontact consumers, frame
changes and atomic retry. Coefficients and the effective texture axes are
caller data. These tests establish the numerical associated elliptical model,
not measured directional behavior, elastic microslip, wear, or a rule for
combining two differently textured surfaces.

### Direct elliptical friction with full body compliance

`elliptic::SymmetricTangentialCompliance` and
`EllipticCoulomb::{advance_with_compliance, advance_kinetic_with_compliance}`
resolve the associated directional law against a finite, symmetric,
positive-semidefinite two-axis compliance in the accepted tangent frame.
The explicit matrix is the sole compliance owner; the input scalar compliance
must be zero. Exactly isotropic matrices delegate to the original scalar path.

The direct implementation complements the existing generic matrix adapter:
it solves `(D C D + lambda I) y = D du_free`, `|y| = 1`, with a fixed-size
scaled eigensystem and a scalar search bounded by `max_iterations`, without
heap allocation or repeated constitutive probes. It retains off-diagonal
response, original texture/history, physical traction-slip work and
world-frame compatibility residuals. Kinetic continuation preserves the
original static coefficients in the returned checkpoint. No measured
throughput improvement is claimed.

Positive-definite compliance admits the unique static reaction. Rank-one
compliance admits resolved sliding and static reactions determined by a
one-axis friction law. Nonzero feasible static demands with two active
friction axes conservatively refuse, including a boundary solution; this
refusal does not assert that every such reaction is nonunique. Invalid,
unrepresentable or unresolved cases leave accepted history unchanged.
The local law supplies no body trajectory, normal impulse or arrest
localization.

`tests/elliptic_coulomb/compliance.rs` exercises independently constructed
associated solutions, scalar equivalence, rotated frames, degenerate axes,
singular response and refusal boundaries.

### Continuous directional sliding and finite arrest

`EllipticCoulomb::relax_velocity` in `elliptic::relaxation` evolves the
original associated kinetic law with constant normal load, coefficients,
temperature, geometry and positive-definite two-axis body mobility. It retains
the original static coefficients, texture and slip checkpoint.

Writing the principal capacities as `a*d`, `max(d)=1`, the transformed
velocity `u=d*v` satisfies `du/dtau=-(d W d)u` and
`dt/dtau=norm(u)/a`. The existing normalized Gram rotation and symmetric
eigensystem provide closed velocity, integrated impulse and physical-work
expressions. Bounded deterministic quadrature resolves physical time, net
displacement and unsigned travel on a finite interval whose endpoint is
complete arrest. These are distinct accounts: mean force times net displacement
does not represent the work of a curved friction trajectory.

Isotropic/eigen-direction motion has an analytic branch. A zero kinetic axis
retains free motion, including cross-axis recoil after active-axis arrest.
A separate short-time chart, conservative clock bracket and `expm1` retain
small impulses/work without integrating the entire stopping history.

Time/travel tolerances, evaluation and root budgets are explicit. Quadrature
errors are adaptive estimates, not interval certificates; a target inside the
unresolved stopping-time uncertainty refuses. Cancellation, invalid mobility,
foreign history and numerical/work-limit failures publish no checkpoint.
The owner supplies no external tangential drive, pose/normal evolution,
material-source refresh or thermal store. Six new G1/G3/G4 groups cover exact
and curved references, split continuation, frame covariance, free-axis motion,
stiff short steps and atomic refusals.

### Full two-axis tangential body compliance

`tangential::matrix::advance_with_matrix_compliance` couples the original
`TangentialLaw` to `du = du_free - C F` with the full symmetric positive
semidefinite 2-by-2 body matrix. The mechanical owner derives C; no effective
isotropic mass replaces its off-diagonal or unequal-axis response. A zero
scalar compliance in the input avoids two independent body-response owners.
Exactly isotropic and zero matrices delegate to the native law in one call.

General matrices use bounded inverse-secant iterations around that same law's
scalar-compliance solve. The full matrix's zero-slip reaction is tried through
the original law first, including the minimum-norm reaction for rank-one
response. This preserves a feasible static reaction's priority over a separate
kinetic root, and a refused static force budget cannot silently become sliding.
Every probe starts from the original accepted state,
and scalar probe compliance changes only the numerical equation. The smallest
positive matrix eigenvalue preserves the native static/kinetic selection;
rank-one null motion bounds the probe through the explicit force budget.
Backtracking retains native source/domain refusals. The returned candidate must
close the original full-matrix displacement equation, original reversible
storage, native dissipated work and actual body work under explicit tolerances.
No residual is relabeled as heat, and no candidate is published on refusal or
cancellation. Arithmetic symmetry/PSD admission never rewrites physical matrix
entries. Original-law calls are explicitly bounded (1 through 4096), with a
cancellation poll before and after every call; the adapter allocates no solver
workspace beyond the original law's candidate state.

G1/G3 groups in `tests/tangential_matrix.rs` compare static and sliding Coulomb
against an independent matrix/KKT multiplier solve, and retained two-axis layer
elasticity/plasticity against an independent flow-multiplier solve. They cover
off-diagonal response, unequal static/kinetic limits, release of retained shear,
nonproportional paths, tangent rotations, singular/free limits, actual heat and
atomic refusal/retry. This is a numerical embedding of each original law, not
a new friction law or a nonproportional Hertz--Mindlin solution. Nonmonotone
laws retain their disclosed native branch selection; uniqueness, contact-frame
transport and simultaneous multiple-contact convergence are not implied.

### Simultaneous retained tangential contacts

`tangential::multiple::PreparedTangentialContacts` admits and caches the full
row-major `2N`-by-`2N` mechanical compliance, including every cross-contact
block. Complete-pivot positive-semidefinite admission checks the full matrix;
positive diagonal blocks alone are insufficient. A range factor and its Gram
factor supply a supported minimum-norm zero-slip reaction, so duplicate contact
rows do not arbitrarily assign all static force to the first site. These factors
are numerical admission/initialization only: physical matrix entries remain
unchanged and all accepted work uses the original full response.

Borrowed per-contact law instances, accepted histories and inputs permit
different source parameters while keeping one held interval and caller contact
order. Deterministic block Gauss--Seidel supplies each existing two-axis matrix
solver with the latest other-contact recoil. Every evaluation starts from the
same accepted history; sweeps never advance time, discard preload or deposit
heat. Each complete sweep independently reconstructs every contact's actual
displacement, reversible storage and physical heat through the full matrix,
checks local work, and checks their summed energy balance before returning all
candidate states together. Source, force, cancellation or convergence refusal
publishes none. Original static priority is tried through the global supported
reaction; nonunique cone feasibility and universal convergence are not claimed.

Contact count, matrix terms, complete sweeps, original-law calls, per-contact
compatibility/work and collective work all have explicit ceilings. The prepared
owner retains at most three dense arrays and O(N) metadata; preparation uses one
additional dense scratch array. Reuse avoids repeated O(N^3) PSD admission.
Each sweep costs O(N^2) matrix work plus the bounded original constitutive calls;
no numerical regularization, source replacement or extra storage is introduced.
Caller order discloses the block/root selection when the physical problem has
multiple solutions. Changing geometry/frames, normal contact and body stepping
remain with their mechanical owners.

Six G1/G3 groups in `tests/tangential_multiple.rs` compare complete elastic
equilibrium against independent Gaussian elimination and multisite plastic
history against an independent global Newton equation. They exercise full
cross-contact recoil, globally feasible static reactions with kinetic rivals,
mixed branches, separate tangent rotations/contact permutations, minimum-norm
redundant rows, incompatible drives, full PSD refusal and atomic bounded retry.

### Objective retained tangential history

`tangential::transport::FramedTangentialState` ties the original constitutive
checkpoint to its accepted `partial_slip::TangentFrame`. The optional
`ObjectiveTangentialLaw` capability reuses that frame's existing minimum-normal
rotation followed by caller-supplied common material spin. A tangent-gauge
change only re-expresses coordinates. The common spin angle is distinct from
relative torsional slip; neither a time interval nor a constitutive relaxation
is inferred. Finite rotations follow the supplied normal path, including its
geometric holonomy, so intermediate geometry and spin refinement belong to the
mechanical owner. Antipodal normal changes require resolved intermediate frames.

The native dry ladder and source card rotate accumulated slip and retain the
optional thermal boundary. `ThirdBodyLayerContact` rotates both elastic and
plastic vectors and keeps the original private law binding. No history reset,
parameter continuation, source replacement, new work or heat occurs. Candidates
must preserve finite nonnegative storage under explicit absolute and relative
roundoff tolerances. Their signed storage residual remains visible for the
coupled energy gate; it is never relabeled as heat. Refusal and bounded
cancellation leave the accepted state and frame unchanged. Scalar Mindlin and
other unsupported histories do not implement this capability.

`tests/tangential_transport.rs` covers loaded gauge rebases, minimum-normal
tilts, physical common spin, global rigid covariance, inverse and path
composition, elastic/plastic norm and storage invariance, the original dry-card
domain, post-transport nonlinear full-matrix continuation, and atomic refusal,
cancellation and retry. This orientation transfer supplies no contact geometry,
patch migration, wear, material advection or general three-dimensional plasticity.

### Exact Hertz sliding and spin wrench

`sliding_patch::HertzSlidingPatch` evaluates the actual surface integral of
local kinetic Coulomb traction over a circular Hertz pressure distribution.
The supplied two-axis centre slip and normal-axis spin share that one local
velocity field, so the resulting force and spin moment are coupled. Pure
translation recovers `mu P`; pure spin recovers `3 pi mu P a / 16`, using the
pressure first-moment radius instead of the RMS radius used by the separate
lumped partial-slip rung. Arbitrary changes in sliding direction, frame rotation
and spin sign require no authored friction ellipse or direction interpolation.

The exact piecewise Hertz integrals follow Antali (2024), *Nonlinear Dynamics*
112, 6785--6811, DOI `10.1007/s11071-024-09462-6`, equations 62--63. A convergent
32-term small-spin expansion evaluates the same function without the severe
cancellation of its large slip/spin ratio formula. Coefficients are computed
at compile time and evaluated with bounded Horner work; there is no allocation,
quadrature grid, mutable state or velocity regularizer in production. Physical
power is the integral of local friction times local slip speed and equals the
negative resultant wrench/velocity pairing. It is not yet allocated to bodies.

The original normal and interface owners must supply load, radius and the
uniform kinetic coefficient. This source-neutral response does not upgrade
their authority. Complete rest returns zero wrench as one admissible selection,
not a solved sticking reaction. It adds no elastic shear memory, nonproportional
Mindlin partial slip, adhesion, rolling-deformation loss, contact discovery,
normal evolution, mechanical stepper or thermal partition. Consumers must not
discard an existing elastic contact's storage when switching to this rung.

Seven native regression groups compare the exact force, spin moment and heat
against independent two-dimensional Hertz-pressure quadrature, with a live
three-rung radial and angular refinement. They also check analytic pure-motion
limits, coupled capacity reduction, arbitrary tangent rotation/reflection,
positive velocity scaling down to tiny speeds, small-spin asymptotics, branch
continuity, deterministic replay, rest/frictionless selections and refusals.

### Coupled Hertz static and sliding wrench

`sliding_patch::contact::HertzWrenchContact` solves a candidate held force and
normal-axis moment against the full three-coordinate body response. It uses
the original analytic Hertz support function, including its stable small-spin
series. Radius-scaled conjugates are `x = (du_x, du_y, a dtheta)` and
`w = (F_x, F_y, M/a)`; their pairing has units of work. Mechanical owners supply
the complete response, including translation/spin cross terms, and retain
responsibility for contact geometry, normal evolution and integration accuracy.

A feasible global zero-slip reaction has priority under the original local
static coefficient, which must be at least the kinetic coefficient. Static
membership uses the true axisymmetric joint Hertz surface: bounded inversion
maintains inner/outer force-moment brackets. This is numerical admission using
the original analytic integral, not a formal rounding certificate; unresolved
boundary membership refuses instead of choosing a kinetic branch. Separate
pure-force and pure-moment limits do not define the joint friction surface.
An admitted static reaction that exceeds the caller's force or moment budget
also refuses; kinetic friction cannot conceal that budget error.

Sliding minimizes the strictly convex body-response energy over the original
kinetic wrench set. A deterministic fully corrective support iteration retains
at most four feasible vertices generated from Hertz support atoms and solves
their constant-size affine faces. At the objective's rounding floor, bounded
convex combinations refine the native wrench residual in the body-compliance
metric; all candidates remain feasible. Explicit support/iteration limits,
native force/moment residuals, dual gap, original-matrix compatibility and actual
physical work gates govern acceptance. The original local sliding integral is
reported as unpartitioned heat. Work/duality/compatibility residuals remain
separate and are never thermalized. Allocation-free arithmetic admits positive
definite response above the declared pivot floor and exact zero response;
general singular mobility is refused without added physical compliance.

The law is immutable and candidate-only, with bounded cancellation polls and
no elastic history to mutate. It adds no partial-slip memory, rolling loss,
adhesion, thermal partition or source authority. An existing elastic contact's
stored energy cannot be discarded when choosing this separate rigid-patch rung.
`tests/hertz_wrench_contact.rs` checks static priority and joint-capacity
rejection, analytical translation/spin limits, independent anisotropic Coulomb
KKT, full cross-response/local work against spatial Hertz quadrature with a live
three-rung convergence check, frame covariance, dimensional similarity, exact
zero limits, and atomic cancellation/refusal/retry.

### Partial-slip return-map rung

`partial_slip` is a solver-independent, two-region Jenkins/Iwan
*return-map* rung (`two-region-iwan-frame-transport-v3`). `NormalPatchView`, `PartialSlipInterface`, and
`GeneralizedWorkOwnership` retain the caller's named normal-patch, ordered
interface/history, and work-owner inputs; their authority remains explicitly
caller-declared. `PartialSlipLaw::advance` has no `fs-contact` dependency and
cannot solve normal contact, evolve patch geometry, or admit a material card.
Two parallel regions retain independent elastic shear/twist and kinetic-branch
history. Their stiffness fractions are `1-h` and `h`, their Coulomb-capacity
fractions are `a(1-h)` and `1-a(1-h)`, and `0<a,h<1`. Holding the displacement
preserves preload without relaxation or heat; reversal first unloads elastically.
One yielding region means partial slip, both means gross slide. The mixed
translation/spin return satisfies maximum-dissipation normality in the elastic
compliance metric; pure translation or spin uses a closed-form return.
This finite two-region approximation is not a resolved Cattaneo--Mindlin
pressure/traction field or a measured slipping-area fraction.

Work is independently computed as `-F dot du - M dtheta`. The balance separates
reversible storage, physical plastic-slip heat, and the backward-Euler spring
loss `sum(k_i |Delta e_i|^2 / 2)`. Numerical loss is never thermalized. Its
refinement claim applies to smooth fixed-branch histories; unresolved
static-to-kinetic jumps require event resolution and have no such claim.
The endpoint wrench is unchanged by energy bookkeeping. The work is a discrete
endpoint contraction, not an exact continuous-time integral or midpoint work.
Purely elastic closed cycles generate no physical heat, while their numerical
loss decreases at first order as the interval is refined.
Retained shear is transported before applying the next relative displacement.
The declared geometric convention is the minimum rotation between successive
contact normals, followed by the supplied common material spin about the new
normal. The final tangent-gauge change only re-expresses coordinates. A single
normalized plane rotation acts on both regions, preserving spring lengths;
normal transport therefore introduces neither projection damping nor heat.
Relative spin alone changes stored twist. Common spin is a distinct mandatory
kinematic input, `(omega_a + omega_b) dot normal / 2` for two rigid surfaces;
it cannot be inferred from an arbitrary tangent gauge. Normal tilt and spin
are composed as an interval rotation, with the same time-resolution limitation
as the surrounding constitutive step. Nearly antipodal normal jumps refuse
because the minimum-rotation path is unresolved; supply intermediate geometry.
Checkpoints contain both regional strains, kinetic branches, and the last
accepted frame, plus all law, patch and interface identities. A nonempty
frame-free or aggregate-strain checkpoint is insufficient. Equality of those
inputs is not external physical validation.
Rolling-deformation loss remains a separate, zero channel in this rung.

### Fixed-load Hertz--Cattaneo--Mindlin contact

`mindlin::MindlinDeresiewiczPoint` implements `TangentialLaw` for a circular
Hertz point patch at fixed positive normal load and one unchanged tangent
direction. `MindlinPointParameters` supplies the patch radius, normal force,
both bodies' shear moduli and Poisson ratios, local Coulomb coefficient and
memory budget. Stiffness is derived as `k0 = 8 a / sum((2-nu_i)/G_i)`;
no tangential spring or slip-onset fraction is authored independently.
The virgin curve is `Q = mu P [1-(1-u/d)^(3/2)]`, where
`d = 3 mu P/(2 k0)`. Gross slip recovers Coulomb and zero coefficient is
exactly frictionless.

The Hertz continuum of elastic Coulomb elements is integrated analytically,
retaining piecewise affine strains over its yield thresholds. This reproduces
proportional Mindlin/Masing reversal loops and return-point memory without a
finite spring grid. Work separates reversible strain energy from physical
microslip heat. The returned effort is the path-average force, and the caller
body compliance is solved using that same effort. The endpoint force remains
available on `MindlinState`; substituting it into a held-effort work audit
would be incorrect. Root trials allocate nothing. Accepted state uses at most
the caller's 1--4096 analytic intervals; nested reversals exceeding that cap
refuse atomically, while erased history is coalesced.

The model follows the decoupled Cattaneo--Mindlin approximation, not a full
three-dimensional elastic-friction solution. Changing normal load/radius,
foreign-law state, or second-axis motion refuse; line contact, nonproportional
shear, torsion, normal/tangential elastic coupling, temperature-dependent
coefficients and local traction/wear reconstruction remain outside this rung.
No material authority is supplied by numerical parameters alone. References:
Popov, Hess and Willert (2019), sections 4.4 and 4.6.3,
<https://doi.org/10.1007/978-3-662-58709-6_4>; Lyashenko and Popov (2015),
Methods, <https://doi.org/10.1038/srep08479>.

`tests/mindlin.rs` supplies G1 virgin-force verification (relative `1e-10`),
initial stiffness (`1e-12`), Masing reversals and closed-form cyclic heat at
100/200/400 intervals per half-cycle (`2e-11`). The constitutive integrals are
exact within floating-point arithmetic; an independent endpoint-work observer
converges at second order. Tests also check the cubic small-amplitude loss,
gross slip, retained preload, compliant work/release, G0 domain/budget/memory
refusals, G3 step-subdivision invariance and G5 checkpoint replay. These are
synthetic numerical verifications; they do not close the varying-normal-load
or source/consumer requirements of `frankensim-irp-v5-instrument-realism-pexfp.147`.

### Hertz--Mindlin normal loading, unloading and recontact

`mindlin::loading::HertzMindlinContact` supplies the next explicit load-path
rung. Immutable curvature, both elastic bodies and a local Coulomb coefficient
determine the Hertz radius/indentation and the original fixed-load shear law.
`normal_step` represents one monotonic normal segment at fixed macroscopic
shear displacement. `TangentialLaw::advance` subsequently advances scalar shear
at the retained normal load; it refuses implicit normal-load changes. Callers
therefore retain responsibility for ordering/refining simultaneous motion.

The evolving patch retains the original analytic spring-history intervals.
Normal growth preserves the existing strain and adds stress-free material at
the edge. Unloading clips only strains that exceed the decreased local Coulomb
bound; newly detached material releases its entire remaining shear storage.
Positive analytic energy-product integrals resolve this microslip heat without
subtracting large nearly equal stored-energy values. Opening empties the history
only after accounting for its energy, and recontact begins stress-free even if
the bodies moved tangentially while apart. No normal dashpot or restitution loss
is inferred. Normal Hertz work/storage and tangential storage/heat are distinct
receipt fields, and the shared tangential work port excludes normal energy.

This is the classical decoupled Hertz/MDR approximation on one fixed tangent.
The positive-half-line stiffness density is `2 G_MDR`, with
`G_MDR = 4 / sum((2-nu_i)/G_i)`, so its fixed-load limit has exactly the existing
`8 a / sum((2-nu_i)/G_i)` stiffness. The original fixed-load shear integrator,
body-compliance root, heat calculation and memory cap are reused. Normal-path
transfer is linear in bounded retained history; it introduces no spring grid.
Foreign owners, invalid loads or exhausted memory refuse without modifying the
accepted state. Variable material coefficients, rolling, spin, nonproportional
shear and a simultaneous normal/tangential trajectory solver remain outside it.

`tests/mindlin_loading.rs` checks Hertz work, exact fixed-load delegation,
stress-free normal growth, the independent saturated-patch `a^5` unloading
energy, a small normal decrement, separation/recontact, normal subdivision,
load/shear path dependence, atomic memory refusal and spatial refinement of an
independent fixed-coordinate MDR spring quadrature under mixed loading/reversal.

### Rolling and contour-deformation loss rungs

`rolling_loss` provides three immutable, separately sourced coefficient candidates:
the `Leine-style Coulomb contour-force rung` (`F_c = -sign(v_c)c_cN`), a
viscous contour alternative (`F_v = -c_vv_c`, with `c_v` in N s/m), and a
finite-patch hysteretic rolling moment (`M_h = -sign(omega_r)N l_h eta_h`).
The last card declares its loss length `l_h` [m] and bounded loss factor
`eta_h`; it is a model-form candidate, not a resolved deformation field.
Every response retains a named patch/load/curvature receipt, ordered
interface/history, temperature/frequency range, source-card authority, and a
dedicated generalized-work key. Principal curvatures are retained directly;
an equivalent-radius substitute requires explicit caller approximation
authority. The rungs have no blend/preference operation, no numerical error
certificate, and no torsional microslip or gross-slide contribution. Pure spin
therefore produces no torque in this slice; a separately sourced torsional
material-hysteresis rung would need a disjoint work channel.

### Material-derived first-order viscoelastic sphere rolling

`rolling_loss::viscoelastic` supplies a separate `RollingLossLaw::ViscoelasticSphere`
candidate for a deforming Kelvin–Voigt sphere on a rigid plane. Explicit radius,
Young modulus, Poisson ratio, density and standard shear/bulk viscosities determine
the response; no rolling coefficient is fitted or defaulted. The corrected
first-order normal-stress relation of Brilliantov, Pimenova and Goldobin
(2014, <https://arxiv.org/abs/1410.3289>, equations 54–56) gives
`A = (1+nu)/(E(1-nu)) [4 eta_s (1-nu+nu^2)/3 + eta_b (1-2nu)^2]`.
Applying this relation to an advected Hertz pressure and integrating its first
moment gives `M = -A N R omega`. This rolling inference uses the pressure-moment
argument of Brilliantov and Pöschel (1998), not its older nonlinear-viscosity
coefficient. `A` is linear and monotone in both viscosities, regular at `nu=0`,
and equals the common material time for proportional Kelvin–Voigt viscosity.

The original `HertzSpherePlane` checks load/area/material consistency; principal
curvatures must identify the actual sphere and signed contour speed must equal
`R omega`. Caller ceilings bound indentation/radius, speed/shear-wave-speed and
maximum material time divided by indentation/speed. These are explicit
applicability choices, not error certificates. This normal-pressure-only,
small-contact perturbation does not solve the shifted contact boundary,
arbitrary-speed rolling, contact discovery, normal evolution, static traction,
sliding or torsional microslip. A released patch, zero viscosities or zero rolling
rate produces zero loss; spin contributes no torque. Normal Hertz storage remains
with its original owner. Held rolling work enters the existing rolling heat and
checkpoint path, with all generated bulk heat assigned to the declared sphere.
The existing caloric consumer requires that sphere's temperature and heat
destination; subsequent transport remains owned by its thermal network.

Eight physical regression groups compare the material coefficient with an
independent complex reduced-modulus limit, integrate the two-dimensional Hertz
pressure moment under refinement, and check viscosity monotonicity, proportional
viscosity, zero/reversal/scale limits, original patch geometry, source and
approximation refusals, retained checkpoints and interval subdivision.

## Invariants

- All public construction and evaluation paths refuse missing identity,
  non-dry medium, negative/non-finite input, and every non-finite derived
  candidate from finite extreme inputs.
- Dry tangential and scalar resistance oppose their declared relative rates.
  Elastic layers can return stored work; total contact work equals storage
  change plus finite, non-negative physical dissipation.
- Heat shares are finite, non-negative, and sum to one. A dissipation step's
  channels are finite, non-negative, and close to its total before a ledger can
  mutate. All ledger channels and wear candidate totals are checked before any
  assignment.
- Card state ranges are closed and finite; absolute temperature has a positive
  lower bound. A direct `FrictionLaw` remains a low-level constitutive rung,
  while state-domain enforcement occurs at `DryInterfaceSystemCard::query`.
- A flash result is either complete for both receiving surfaces or an explicit
  `Unknown`; it never silently assigns missing thermal properties.
- Partial-slip refuses nonfinite or nonpositive patch and law inputs, malformed
  tangent frames/kinematics/work ownership, and checkpoint mismatches. It
  computes only its declared lumped return-map state and never creates a
  contact-solver or material-admission receipt.
- Each rolling-loss card refuses non-dry interfaces, malformed patch/source
  receipts, out-of-range temperature/frequency, nonfinite derived work, and a
  wrong work channel. A rolling key conservatively refuses sharing the same
  patch/interval pair with partial-slip because that module does not expose its
  coordinate labels for a finer disjointness proof.
- Rolling/contour candidates are passive: endpoint body power is non-positive,
  heat is non-negative, reversible storage is zero, and the retained work
  identity reconstructs heat over the caller interval.
- No mutable material table, hidden state, or unordered material-pair lookup
  exists. Ordered interface and history identities remain explicit.
- Surface excitation refuses missing texture/source identity, uncovered finite
  footprints, excessive periodic repetition/work, non-finite kinematics, a
  height outside the caller's tangent-linearization fraction, and a
  perturbation that would open the unilateral contact. Both surface sources
  and authority ceilings remain separate in the receipt.
- `PeriodicHarmonicSurface` admits either a measured Fourier decomposition or
  an explicitly declared hypothetical periodic material track. Signed
  sine/cosine height coefficients are keyed by positive integer spatial cycle,
  canonicalized by cycle count, and sampled only with at least eight points per
  shortest period. Duplicate cycles and underresolved grids refuse. An empty
  spectrum is an exact smooth track. No material name, acoustic frequency,
  stochastic realization, seam correction, or authority upgrade is inferred.
- `SelfAffinePeriodicProfileSpectrum` refuses nonpositive/nonfinite RMS height,
  Hurst exponents outside `0 < H < 1`, inverted or zero cycle bands, and
  harmonic bands whose shortest wavelength cannot satisfy the eight-samples-
  per-period rule within the one-million-sample trace ceiling, as well as
  unbounded harmonic counts. Every retained cycle owns a counter-based phase
  stream keyed by the explicit seed and cycle, so extending the band preserves
  shared-cycle phases. Phase changes preserve the declared Fourier power, and
  the continuous periodic profile retains the requested RMS height to
  floating-point roundoff.

## Quasistatic Reynolds film

`resistive_film` owns the solver-independent incompressible thin-film pressure
law shared with the compatibility surface `fs-flux::resistive_film`. At each
configuration it solves the bounded pressure island `L(h) p = A B v` and applies
the reciprocal resisting force `B^T A p`. The analytic tangent includes both
gap-dependent conductance and independent mechanical effort changes. Cell areas,
port weights and Poiseuille conductances are never independently normalized.
This law carries no clock or pressure history; the mechanical caller owns
stepping, cancellation between calls and integrated energy accounting.

The fixed limits are 64 pressure cells, 256 channels and 256 mechanical ports.
Nonpositive channel apertures close exactly; collapsed volume cells, trapped
pressure pockets, pressure/gap limit violations and unresolved factorizations
refuse without modifying output slices. No leakage, gap floor, pressure clipping
or diagonal repair is introduced. Compressibility, thermal/slip effects, fluid
inertia, evolving topology, fluid radiation and elastohydrodynamic closure are
outside this law. The existing nine pressure-law tests run through the fs-flux
compatibility surface, while fs-couple tests exercise reciprocal mechanical
reaction, joint dissipation and atomic retry in its existing time integrator.

## Error model

`TriboError` is a total refusal surface for missing identity, non-dry media,
invalid/overflowing physical input, non-finite vectors, normal slip, malformed
or out-of-domain applicability ranges, invalid partitions, and forged/invalid
dissipation states. A refusal makes no partial state change to `WorkLedger` or
`WearState`. Missing flash thermal data is a typed `Unknown`, not an error or a
candidate.

## Determinism class

Pure scalar operations are deterministic for the same input and ISA. Every
solver-state transcendental is routed through `fs_math::det` (bead
3ez8g.7.3): norms use `det::hypot`, which is overflow/underflow-safe BY
CONSTRUCTION (max-scaling, so the ratio never exceeds 1) — the original
platform-`hypot` rationale is preserved, and the kernel is bit-identical on
every conforming target; the Stribeck decay uses `det::exp`, tangent-frame
rotation `det::sin`/`det::cos`, and the Hertz `x^{3/2}` is the IEEE-exact
`x*sqrt(x)`. CROSS-ISA BIT-STABILITY IS CLAIMED for these friction paths and
enforced by a single committed golden hash executed on both aarch64 and
x86-64 (`cross_isa_golden_hash_over_the_routed_friction_paths`), plus the
workspace `check-libm` lint (fs-tribo is a doctrine crate). The recorded
counter-argument: platform `hypot` is near-correctly-rounded while
`det::hypot` spends up to 2 ULP — a cost orders below Stribeck-parameter
uncertainty, accepted for reproducibility. Composed images (e.g. the bowed
string) remain capped by their OTHER operands until audited.

## Cancellation behavior

This bounded scalar leaf has no asynchronous work and no cancellation scope.
Callers scheduling many contact points own cancellation between calls.

## Unsafe boundary

None. The crate forbids unsafe code.

## Feature flags

`benchmark-laws` (default off) exposes `fs_tribo::benchmark`, the M4.2-W4
benchmark tangential laws. They are also compiled under `cfg(test)`. Only
dev-dependencies may enable the feature;
`tests/benchmark_laws.rs::g0_benchmark_laws_not_reachable_from_render_paths`
scans every workspace manifest and fails if a `[dependencies]` table enables
it.

### Benchmark tangential laws (`benchmark-laws`)

These laws exist for named comparisons, such as the X-Bw2 bow-law bakeoff.
They are never render laws.

`RateState` is Dieterich–Ruina rate-and-state friction.
- Friction: `mu = mu0 + a ln(V/V0) + b ln(V0 theta/Dc)`.
- `Evolution::Aging` is `d theta/dt = 1 - V theta/Dc`.
- `Evolution::Slip` is `d theta/dt = -(V theta/Dc) ln(V theta/Dc)`.
- Admission: `a > 0`, `b >= 0`, positive `V0`/`Dc`, and a declared `speed_range_mps`. Outside that range a step returns `OutsideApplicability`. A zero-speed step refuses because `ln V` is undefined.
- Over one held-effort interval, both state laws integrate exactly in closed form using `expm1`.
- Friction uses the state at the start of the interval. It is then strictly increasing in `V`, so the compliant root is unique.
- `critical_stiffness_n_per_m(N) = N (b - a)/Dc`.

`LuGre` (Canudas de Wit et al. 1995) is one-dimensional.
- Bristle law: `dz/dt = v - sigma0 |v| z/g(v)`, where `g(v) = Fc + (Fs - Fc) e^{-(v/vs)^2}`. It is integrated exactly at held `v`.
- Force: `F = sigma0 z1 + sigma1 (z1 - z0)/dt + sigma2 v`.
- `sigma0 z^2/2` is reversible storage, and the remainder of the port work is heat.
- A step whose dissipation would be negative refuses, rather than assuming the `sigma1` passivity condition.
- A two-dimensional slip increment refuses.

Both laws implement `TangentialLaw`:
- The compliant branch bisects with the caller's `TangentialSolveConfig` budget.
- A rate-and-state root that needs a speed outside the admitted range returns `OutsideApplicability`.
- Parameter caps are checked at `cap` and at the next representable value above it.

`SmithWoodhouse` is the plastic thermal model of rosin from Smith & Woodhouse
(2000), J. Mech. Phys. Solids 48:1633-1681. Equation numbers are the paper's.
- Friction (Eq. 21): `mu N = A k_y(T) sgn v`, set-valued in sticking, with `A = pi a^2` proportional to `N`.
- Temperature: the heat balance of Eq. (6) per unit area. Conduction into the stationary and the moving solids uses the discrete half-space Green's function of Eqs. (10)-(14). The moving solid sees only the history within slip `8a/(3 pi)` (Eq. 16 at constant speed), with the cut pulse weighted exactly.
- `k_y(T)` is calibrated from the steady-sliding fit (Eq. 2) through Eq. (20) (section 5.3). `SmithWoodhouse::new` certifies with an interval bound that the calibration temperature increases with speed, and refuses otherwise.
- Time stepping: backward Euler in the losses, with friction from the step-start temperature, over a declared uniform step. The paper uses a forward difference.
- Refusals: a different step length, zero load, a temperature above the calibrated range, a history beyond the declared capacity (at most `MAX_HISTORY_STEPS`), and non-positive parameters.

## Conformance tests

Inline and external G0/G1/G3 tests cover identity/dry/tangent refusals;
ordered-card applicability and slip-speed consistency; analytic rigid-block
equilibrium; independent numerical Hertz values and force-pressure-radius
cross-relations; reversal/scaling/provenance; resistance sign; partition and
ledger-channel closure; typed flash-data insufficiency and power scaling;
forged negative/non-finite work rejection; and rollback/deterministic replay
for work and wear. Surface-excitation tests additionally cover periodic-spectrum
order invariance, exact sample-grid realization, duplicate-cycle refusal, and
the eight-samples-per-shortest-period gate. They also cover deterministic
self-affine replay, exact RMS normalization, the declared PSD slope, seed-only
phase changes, and malformed statistic/band refusal. Test coefficients are
explicitly synthetic fixtures.
The separately owned partial-slip tests cover its scalar return-map admission,
reversal, checkpoint replay, and receipt-mutation refusals.
`tests/tangential_law.rs` adds G0 DEF-06 continuity at `v_reg*(1 +/- 1e-12)`
with a `1e-12*mu_s*N` band for all three rungs; 10,000 Philox-seeded steps
per rung with non-negative heat and per-step closure within
`1e-12*max(abs(F dot du), f64::EPSILON J)`; G1 incline acceleration within
`1e-12 m/s2`; closed-domain next-up refusal; force/iteration-budget refusal;
and G5 checkpoint replay. The actual modal owner is compared bitwise in
`fs-couple`'s `friction::scalar::trait_tests`, avoiding an upward dependency.
These are named verification targets; their execution is required before
claiming numerical qualification.
Rolling-loss tests cover zero load/rate/factor, force and moment reversal,
load scaling and SI dimensions, pure roll/spin/mixed motion, passivity with an
independent work integral, temperature/frequency refusal, work-key overlap,
source-card/patch replay mutation, and a no-target-field source scan.
`tests/benchmark_laws.rs` (feature `benchmark-laws`) covers:
- **G1 rate-and-state steady state:** the `mu_ss = mu0 + (a - b) ln(V/V0)` log law at relative 1e-10, for both evolution laws across 1e-8..1e-2 m/s.
- **G1 velocity step:** the 1/e transient decays over `Dc` of slip, within 1 % at 12.5/25/50 steps per `Dc`.
- **G1 Ruina spring-slider:** stiffness sweeps through the compliant `advance` path bracket `k_c = N (b - a)/Dc` between `0.98 k_c` and `1.02 k_c`. They compare oscillation amplitude early against late, over 300 `Dc` of slip.
- **G1 LuGre:** the steady-state Stribeck curve at relative 1e-10, and presliding stiffness `sigma0` at relative 1e-8.
- **G0 LuGre bristle map:** passivity over Philox-seeded histories (12,800 intervals).
- **G1 LuGre drive-spring:** stick-slip exists and is periodic.
- **G0:** parameter cap and cap+1.
- **G1 Smith & Woodhouse Green's function:** the weights match Eq. (13), its telescoped partial sums and Eq. (15).
- **G1 Smith & Woodhouse calibration:** Eq. (20) temperatures and Eq. (2) friction at relative 1e-12, over 1e-4..0.3 m/s for four rig geometries.
- **G1 Smith & Woodhouse steady sliding:** friction approaches Eq. (2) as `t^{-1/2}`, the stationary half-space tail. One Richardson step leaves under 10 % of the 4 s deviation.
- **G1 Smith & Woodhouse sinusoidal sliding (Fig. 18):** clockwise loops confined to low speed, which close at the reversal as `O(h)`. The three-rung step ladder has order in [0.8, 1.2]. The temperature levels are printed and not gated: they do not reproduce Fig. 18(a) or the 75-115 K stated for a 0.25 mm radius.
- **G0 Smith & Woodhouse:** caps, refusals and the set-valued stick/slip branch.
- **G0 reachability scan.**

## No-claim boundaries

The rate-and-state and LuGre benchmarks carry synthetic parameters. The Smith
& Woodhouse tests use that paper's published rosin, perspex and steady-sliding
inputs, for comparison only. No benchmark claims a material admission, bow
data, or equivalence to any render law. A benchmark result names a comparison and does
not admit a law.

This crate does not mint material admission, calibration, the upstream
`fs-matdb::InterfaceSystemCard` query receipt, a flash-temperature error bound
or thermal-port solution, a roughness measurement or statistical roughness
model, normal adhesion or adhesion-induced contact-area evolution, general three-dimensional plasticity, finite-patch partial
slip as a resolved contact field, general lubrication/EHL, contact geometry evolution,
wear geometry updates, stop time, Euler-disc ranking, one-millimetre optimum,
sound pressure, loudness, or experimental/video correspondence. The surface
trace leaf consumes a caller-supplied one-dimensional profile and a local
contact tangent; it does not replace a two-dimensional roughness map, nonlinear
re-solution of contact, structural response, or acoustic radiation. Its flash candidate is a declared
semi-infinite uniform-flux model-form estimate, not a temperature measurement
or validation result. A caller may record an estimated or synthetic result,
but must not promote it based on this crate alone.

The self-affine constructor is likewise a geometry realization, not evidence
that a specimen is self-affine or that its RMS height, Hurst exponent,
wavelength cutoffs, or phase realization are correct. Those inputs require
profilometry or an explicitly weak estimate. The one-dimensional profile does
not establish an isotropic two-dimensional areal PSD or resolve asperity-scale
partial contact.

The rolling and contour rungs likewise do not establish deformation-loss
dominance, coupon calibration validity, normal-contact admission, an error
bound, or any apparatus-specific selection between their rival model forms.

## Explicit semi-infinite heat allocation

`HeatPartition::semi_infinite_surfaces` allocates an interface-generated heat
pulse in proportion to `e = k / sqrt(alpha) = sqrt(k rho cp)`. It uses the
existing checked surface properties and heat partition. The normalized ratio
does not overflow when the sum of two finite effusivities would. This is the
equal-exposure, homogeneous isotropic half-space approximation: the two
incremental surface-temperature rises agree and the two fractions sum to one.
The caller owns contact resistance, moving-source effects, finite thickness,
material state and applicability. No third-body heat share is inferred.

G1/G3 coverage in `tests/dry_baseline.rs` checks an independent effusivity
ratio, equal incremental temperatures through the existing flash-temperature
formula, ordered-surface reversal and extreme finite normalization. The
semi-infinite boundary assumptions and finite-size limitations follow A. Jain,
*International Journal of Heat and Mass Transfer* 202 (2023), 123721,
<https://doi.org/10.1016/j.ijheatmasstransfer.2022.123721>. This is numerical
verification of a selected model, not experimental material qualification.

## Adhesive junction shear

`adhesive::AdhesiveFrictionLaw` retains shear resistance of an attached junction
at zero or tensile net normal force. Its static and sliding capacities are
`tau * A + mu * N_rep`, with independently supplied real contact area, junction
shear strengths, and the nonnegative repulsive contribution to normal force.
The signed net adhesive force must not be substituted for `N_rep`. Setting the
junction strengths to zero recovers the existing pressure law; setting its
coefficients to zero selects constant junction shear strength. Calibrations
must avoid counting the same junction resistance in both terms.

The leaf borrows the original interface provenance, allocates nothing during
evaluation, opposes sliding, and reports nonnegative frictional power. The
scalar regularization ramps the complete capacity continuously near zero slip.
`tests/adhesive.rs` checks finite zero-load friction, detachment, pressure-law
recovery, extensive scaling, vector power, reversal, regularization, input
refusals, and subnormal velocities. Contact geometry, partial slip, shear-induced
peeling, bond aging and measured strength admission remain the caller's concern.
No strength is inferred from adhesion energy or bulk material names. The
area-proportional term follows the measured single-asperity relation in Carpick
et al., *Langmuir* 12 (1996), <https://doi.org/10.1021/la9509007>; the optional
pressure term is an explicit uncoupled model choice, not a universal fit.

## Changing-coefficient partial-slip continuation

`PartialSlipLaw::advance_with_coefficients` continues the same contact/source model with changed Coulomb coefficients or translational/torsional stiffness. It retains both regions' elastic strains and kinetic branches and runs the original return map; regional stiffness/capacity fractions and the torsional capacity factor remain fixed. Its independent parameter-work port supplies `C = psi(new coefficients, old history) - psi(old coefficients, old history)` before motion. Returned storage change refers to the actual old law, so `delta(psi) + physical slip heat + numerical loss = wrench work + C`. Softening can export parameter work; reduced yield strength can relax stored strains even at zero imposed motion. Ordinary restore remains coefficient-exact.

The coupled owner must accept the coefficient driver and contact history together and balance the signed parameter port. This first-order parameter-then-motion split supplies no temperature-dependent material data or caloric closure: temperature-dependent spring free energy requires mechanical internal energy `U = psi - T d(psi)/dT`, so parameter work alone is not heat. Numerical dissipation also remains nonthermal. G0/G3 tests in `tests/partial_slip_coefficients.rs` cover signed stiffness work, independent mixed shear/spin energy, retained preload and normal-load changes, zero-motion relaxation, unchanged-law identity, and strict source/decomposition refusal.

## Adhesive held-effort tangential contact

`adhesive::AdhesiveTangentialContact` freezes caller-resolved real contact area
and nonnegative repulsive force and implements the existing `TangentialLaw`.
The shared two-axis held-effort solver selects actual static reactions inside
`tau_static A + mu_static N_rep`, or solves compliant sliding using
`tau_kinetic A + mu(v) N_rep`. It does not substitute signed net normal force
or an effective Coulomb load. Finite real area carries shear at zero repulsion;
a detached patch carries no traction or heat. The redundant input normal load
must equal the frozen repulsive load.

The original dry checkpoint and immutable candidate protocol are retained.
Reversible tangential storage is zero, and exported heat is the applied force
dotted with the solved relative increment, not the released predictor. The
pressure-only limit retains the original law and arithmetic. Coulomb and
linear velocity laws have checked affine solutions, including a pressure
coefficient clamped to zero while junction shear persists. Other laws use the
existing bounded root rule; no new multiroot history or branch certificate is
claimed. Normal evolution, partial slip, shear-induced peeling, ageing and
thermal material refresh require their own owners.


## Simultaneous original Hertz contact islands — v1

`sliding_patch::contact::multiple::PreparedHertzContacts` retains the full
row-major `3N × 3N` radius-scaled response for the original integrated Hertz
shear/spin laws. All off-diagonal patch blocks remain in the actual load line.
One preparation admits and factors that response for repeated coefficient
trials; no sweep advances constitutive time or accepts an individual patch.
The bounded pivoted positive-semidefinite response factor is shared with the
existing two-axis tangential owner. Its supported minimum-norm reaction is a
numerical seed only; physical compatibility always uses the original matrix.

Dependent coordinates can admit another static load allocation even when
that seed exceeds one patch's capacity. A bounded static-feasibility phase
therefore seeks redistribution inside every original static Hertz set. Its
temporary equal-coefficient projections supply no physical heat or published
mode. Kinetic alternatives require a directly checked separation inequality
against the original static supports, with an explicit work margin and an
arithmetic guard that includes cancellation in the full transposed response.
This is a numerical separation check, not an interval certificate. Unresolved
static feasibility explicitly refuses. Frictionless patches have exactly zero
wrench and are not additional zero-motion constraints. Small unsupported
motion cannot disappear merely because it falls below displacement tolerance.

Acceptance independently rechecks every patch's original kinetic support or
static membership at the final common motion. Per-patch native force, moment,
work and displacement residuals and the total physical heat/work balance must
meet the declared budgets. Intermediate block reactions do not consume
accepted-effort limits; those limits apply to the fully resolved physical
wrench. Numerical residuals and static-search work are never heat. Exact zero
response preserves prescribed slip, even below displacement tolerance.

Preparation is bounded by at most 128 contacts and 147456 matrix entries;
static and kinetic work share at most 4096 sweeps and 1048576 native support
evaluations. The original matrix, factor and temporary/range factor require
at most three dense arrays, plus linear per-contact state. Cancellation is
polled during preparation, solves and final admission. Failure leaves the
reusable owner and all caller state unchanged. The caller supplies a stable
physical contact order. Every local three-coordinate block must retain the
original positive-definite or exact-zero admission; no regularizing compliance
or diagonal shift is introduced.

Nonunique static allocations and different static/kinetic coefficients retain
explicit deterministic root selection without general convergence or uniqueness
claims. Normal evolution, geometry, thermal storage, rolling loss and retained
elastic partial-slip memory remain with their original owners. Focused G1/G3
tests in `tests/hertz_wrench_multiple.rs` check dense static reactions, mixed
stick/slide, unequal-radius spin, independent pressure quadrature, exact zero
limits, dependent-load redistribution, unsupported motion and bounded refusal.
The original `tests/tangential_multiple.rs` exercises the shared factor through
its unchanged retained-contact and history semantics.
