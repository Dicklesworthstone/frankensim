## Analytic affine Hertz stopping, static hold and reversal

The original rate-driven Hertz owner now recognizes invariant pure translation
and pure spin under the complete supplied mobility and parallel external drive.
It uses the original affine-law classification, including genuinely affine
limits of the other source families. Translation retains the Hertz normal-load
moment; spin retains the dry moment `3*pi/16` and viscous moment `2/5` in
radius-scaled coordinates. Each fixed-sign branch integrates velocity, traction,
travel and physical power analytically.

Analytic admission bounds accumulated transverse rate and travel over the whole
phase before consuming a constitutive evaluation. Insufficient collinearity
accuracy falls back to the original full-matrix IVP.

Moving states enter the static solve only at their resolved first zero.
The original joint static surface then holds the load or releases a reversed
branch. Clipped affine translation splits at its actual zero-traction boundary,
including coasting and subsequent re-entry. Small-speed capture, arbitrary
weakening-root selection and endpoint energy residuals are not used as heat.
Asymptotic viscous decay and unstable weakening equilibria retain their original
initial-value behavior; unrepresentable moving endpoints refuse.

`exact_invariant` distinguishes this formula path from adaptive estimates.
Its first stopping interval is a point at the computed analytic time, not a
validated floating-point enclosure. Original effort, work, evaluation/attempt
budgets and cancellation gates remain active. Partly clipped spinning disks,
non-affine curves and non-invariant motion retain the existing original-pressure
integrator and its explicit unresolved-event boundary. Load, geometry, source
parameters and drive remain held within each phase.

## Driven original local-speed Hertz trajectories

`HertzFrictionPatch::advance_driven` follows the original three-coordinate
initial-value problem under held normal load, radius, source curve, geometry,
external acceleration and full SPD mobility. Every positive RK4 stage evaluates
the original local-speed pressure field and integrates actual power, impulse
and travel. Weakening evolution is not replaced by an arbitrary algebraic root.
The shared RK4 step-doubling primitive supplies accumulated rate/impulse/travel/
heat estimates, separately from pressure-panel refinement and original-matrix
work gates. The existing passive arithmetic operation tree remains unchanged.

Rest is an explicit retained mode. Original joint static support either holds
the drive or supplies the zero-speed kinetic right-limit release direction.
An original separating support plane bounds the first moving interval; later
intervals exclude zero through the original mobility norm. The numerical finite
horizon envelope includes positive drive quadrature and admitted work defect;
every endpoint/stage must stay within its admitted mass-radius/local-speed
domain. The whole-path native wrench bound includes static reactions as well
as the original finite-speed kinetic coefficient bound.

Nonconstant driven zero-rate events lacking a proved continuation return
`UnresolvedDrivenEvent`, without speed capture, passive stopping tails or invented
terminal heat. Genuine Coulomb delegates its existing event owner and explicitly
marks adaptive error estimates unavailable. Unloaded nonconstant curves retain
finite free-path bounds on their original coefficients. General curved or non-affine
arrest/restart, evolving normal geometry and source changes within a phase remain
outside this owner. Focused tests cover analytic forced affine and full-SPD
viscous motion, unstable weakening, clipped re-entry, independent Stribeck
pressure refinement, retained mode, covariance/replay and bounded refusal.

## Driven Hertz shear and spin with retained physical mode

`HertzWrenchContact::advance_driven` advances one held-load/radius/coefficient/
external-drive phase in the original three radius-scaled shear/spin coordinates
and complete symmetric positive-definite mobility. Its explicit `Resting` or
`Sliding` checkpoint controls static admission; a small moving rate never grants
access to the larger static set. Resting-drive feasibility and breakaway reuse
the original joint Hertz static/kinetic support owner, including coupled force
and moment demand. Release from rest and admitted invariant rays have constant-
acceleration solutions, exact first-zero splitting, and subsequent static hold
or overdriven restart.

Changing-direction motion uses bounded first-order kinetic endpoint projections
with held traction and linear velocity. `maximum_step_s` is a required temporal
refinement axis; curved endpoint arrest is explicitly timestep-resolved. Such
capture also requires a strictly interior original kinetic-drive feasibility
probe. Both represented endpoint powers must be nonpositive. Physical heat pairs
acting traction with midpoint travel. Generalized drive work is `(W^-1*a) dot
travel`, independently gated against original-matrix kinetic change plus heat
and endpoint momentum compatibility.

Original support evaluations, retries, effort limits and cancellation remain
bounded. Unresolved direction, work, and rounded zero without a resolved event
refuse before publication. Geometry/load evolution, exact curved-arrest timing
and speed-dependent driven branch selection remain outside this owner. Focused
tests cover analytic pressure-moment events, coupled static separation, reversal
and hysteresis, fixed-grid ODE refinement, and bounded atomic refusal.

## Local-speed Hertz traction through finite-body relaxation

`HertzFrictionPatch::advance_transient` evolves the original local-speed Hertz
pressure integral under the complete held positive-definite mobility. A finite
speed domain follows from initial mass-metric energy plus the admitted work
defect. Every accepted endpoint and tail start is checked; positive RK4 stages
use a 1.25-radius expansion with explicit speed and wrench admission.

Conservative coefficient bounds on the finite speed interval supply the original
physical-power RK4 integration. A positive global lower bound admits the existing
finite-extinction tail, including hyperbolic and shifted-power curves with zero
infinite-speed asymptote. Without global coercivity, the same owner advances the
requested finite horizon and admits a tail only after the checked mass-radius
ball lies in a proved positive-coefficient domain. For clipped affine friction,
all local speeds below half the cutoff guarantee at least half the original
zero-speed coefficient. Pure viscosity retains the exact Hertz linear wrench,
with second-moment factor 2/5; its local linear rate bound keeps the safe timestep
tied to decay time rather than vanishing speed. It has no finite stopping tail.
Rounded zero endpoints and unproven zero traction refuse instead of inventing
arrest or coasting. Genuine Coulomb retains its original exact-ray path; proven
rest, unloaded motion and fully clipped zero traction retain exact stationary
or coasting behavior. No coefficient floor or substitute curve is used.

For an original affine or Stribeck viscous term, the safe stage bound uses
`|w(q)| <= b0 + b1|q|`, with `b0 = N*mu_zero*sqrt(1+k*k)` (mu_static for
Stribeck), `k = 3*pi/16`, and `b1 = N*slope_per_speed` for affine
strengthening or `N*viscous_per_speed` for Stribeck. The 1.25 mass-radius stage
bound contracts this upper bound as motion decays; every stage still evaluates
the original local-speed pressure law. The strict 12-second small-patch affine
spin regression closes against analytic stopping time, travel and heat in
42,406 evaluations within the unchanged 400,001-evaluation budget that the
global bound previously exhausted.

Time error estimates exclude spatial pressure-quadrature error. Callers refine
these two axes separately. Load, radius, mobility and source/environment remain
fixed; driven stick/slip, source evolution and normal motion retain their owners.
Twelve focused law tests cover analytic affine and zero-asymptote stopping/travel,
25 e-fold viscous decay, independently constructed full-mobility modal exponentials
and integrated heat, clipped translation entering its true finite-stop domain,
original clipped-pressure/time integration, live spatial refinement, Coulomb
parity, exact coasting, and bounded atomic refusal/replay across coercive and
noncoercive curves. A native regression additionally checks simultaneous viscous
sliding and spin against the actual bodies' analytic decay, travel, heat and
momentum, without a false finite stopping event.

# CONTRACT: fs-tribo

## Curved Hertz sliding and spin through finite extinction

`HertzWrenchContact::advance_transient` integrates changing-direction slip/spin
under the full positive-definite mobility and original Hertz pressure integral.
Exact invariant directions retain their analytic one-evaluation path. Curved
Hertz paths select positive-stage classical RK4 step doubling in the shared
circular owner. Two physical half-steps supply accepted impulse, travel and
original-power heat; a full step supplies the difference/15 error estimate.
No Richardson-extrapolated physical state or negative heat quadrature is used.
Every coarse/fine stage stays within h*beta of the original mass-metric point,
retaining the existing no-zero-crossing bound. Separate dimensional error
estimates, cancellation, effort limits and bounded attempts/evaluations remain.
Adhesive paths keep their original Bogacki–Shampine 3(2) arithmetic.

For C = mu*N, q = [ux, uy, a*omega] and k = 3*pi/16, symmetry/convexity and the
triangle inequality give C*k/sqrt(1+k*k)*|q| <= H(q) <= C*sqrt(1+k*k)*|q|.
These Hertz-specific bounds and the original matrix bound the stopping tail.
A tail closes only when its entire stop bracket precedes the requested horizon
and its unresolved time/travel fit the requested tolerances. Its eventual-rest
impulse and remaining relative kinetic heat are retained; unknown travel stays
an explicit bound. Independent prefix and total work checks never turn a
numerical residual into heat.

Prefix errors remain numerical estimates, not validated-flow enclosures.
Geometry, normal load, coefficients and mobility remain fixed without drive or
pose evolution. `hertz_transient` tests cover independent curved RK4, analytical
anisotropic-mobility stopping/travel, full-matrix extinction, stationary dwell,
splitting/refinement, frame/reversal, pressure bounds and budget/cancellation
replay. A strict small-radius fixture compares 1e-11 m/s rate and 1e-12 J heat
accuracy with converged independent fixed-step integration, within the original
work budgets. Existing adhesive transient tests retain their original
constitutive law and norm constants.

## Hertz slip/spin transition accuracy

The closed-form Hertz traction uses the same factored complement
`(1-z)*(1+z)` in its angle and square-root terms near unit slip/spin ratio.
This preserves the cancellation in the original pressure integrals without
snapping the branch or changing the friction law. The `sliding_patch` rim
regression checks both sides at 1–512 ulps and three velocity scales against
the exact rim force, moment and derivatives.

## Continuous Hertz sliding and spin through arrest

`sliding_patch::contact::relaxation` integrates the original pressure-resolved
Hertz force/spin wrench through exact finite arrest when the complete
radius-scaled mobility preserves the incoming slip/spin direction. One kinetic
law evaluation supplies the moving phase. Its duration ends at the actual stop;
stationary dwell adds neither impulse nor heat, and static capacity does not
capture moving motion early. The report retains actual integrated slip, spin,
impulses, physical loss and stopping time.

The allocation-free mechanics are shared with
`adhesive::spinning::relaxation`, preserving its existing public types, error
variants and arithmetic. Original full-matrix impulse work independently checks
the integrated physical loss. Geometry, normal load, coefficients and inertia
are held over this unforced trajectory; curved velocity paths and pose
integration require their surrounding owners.

`tests/hertz_relaxation.rs` covers pure sliding/spin, the closed-form combined
Hertz rim response with off-diagonal mobility, time splitting, stopped dwell,
zero capacity and bounded refusal. The existing `adhesive_relaxation` tests
retain the shared owner's prior behavior.

## Pressure-fed sliding Reynolds films

`resistive_film::shear::pressure::PressureDrivenShearedFilm` composes prescribed
reservoir pressures with the original sliding-wall film. A single retained
graph, wall mapping, viscosity profile and pressure factorization own squeeze,
Couette entrainment, pressure flow and complete work-conjugate wall traction.
Supply potentials are compiled once; evaluation and the full directional
tangent use bounded stack storage and allocate no heap memory.

The instantaneous work identity is
`F dot v - sum_boundary(p_b * Q_total,out) = Qp * delta_p + shear_heat`,
summed over channels. Supply power uses total flow, including entrainment.
Irreversible heat uses pressure-driven flow and relative-wall drag, preserving
both pressure-driven mechanical output and wall-driven pumping against supplies.
Mechanical power and supply power may each be negative; their sum is the
nonnegative viscous heat. This is not a time-integrated energy receipt.

Supply pressures retain the original mechanical and zero-gauge viscosity datum.
Barus viscosity uses the exact pressure-dependent shear average on that datum;
retained layers keep their local pressure/shear stress cross terms. The force
tangent differentiates gaps, wall rates and all prescribed boundary pressures
together. Viscosity continuation retains supplies, geometry and wall rows.
Admission failures leave caller outputs untouched.

The additive `evaluate_force_into` omits channel/layer publication for repeated
mechanical trial queries; the original report-bearing API remains available.
Count/range accessors preserve original channels and layers for thermal routes.
Both supplied-film owners use the original passive arithmetic for explicit
zero supplies, including tangents when supply directions are zero. Signed reservoir work uses compensated sums.
Sorted-endpoint Barus recovery retains finite drops when a one-sided logarithmic
ratio would round to its singular endpoint. Pressure-only and fixed-wall paths
recover the physical drop without requiring a separately representable viscosity
multiplier. Moving passages cache the drop and reuse their already required
mean viscosity for traction, heat and geometry derivatives.

The model adds no finite reservoir inventory, supply-line dynamics, thermal
clock, fluid inertia, cavitation or general high-pressure EHL constitutive law.
Regressions cover analytic supplied Couette/Poiseuille work, static Barus load
and drag, independent piecewise Stokes layer integrals, complete directional
derivatives, the zero-supply limit, viscosity refresh and atomic refusal.

## Externally pressurized Reynolds films

`resistive_film::pressure::PressureDrivenFilm` binds an explicit gauge pressure
to every original reservoir channel. The original geometry, conductances,
Cholesky elimination, force transpose and common Barus law remain the owners.
This admits pressure-fed lubrication, stationary load support and hydraulic
actuation alongside squeeze motion. All supplies share the original mechanical
and zero-gauge viscosity datum; their absolute values are not discarded.

The solve subtracts one common reservoir potential before elimination, keeping
a stationary equally pressurized film exact and separating pressure bias from
the hydraulic residual. Physical pressure and the original reference viscosity
are recovered in that chart. Boundary potentials compile once; evaluation and
the full gap/velocity/supply-pressure force tangent use bounded stack scratch
and the native factorization, without an added nonlinear pressure iteration.

Each channel reports signed volume flow and nonnegative physical `Q*delta_p`
heat. The distinct driven report closes `F*v - sum(p_b*Q_out) = heat` and
permits negative mechanical input when the supplies drive the body. Reservoir
power is an external supply, not extra heat or proof of passive mechanical
resistance. Force, pressure and channel outputs publish together after checks.
Wrong reservoir topology, source/pressure limits, trapped islands and failed
work balances refuse. The caller owns time integration, reservoir inventory,
cavitation, supply dynamics and thermal destinations; this adds no compressible,
sliding-wall or high-pressure EHL model.

Seven focused G0/G1/G3 native tests cover pressure-fed flow/load, hydraulic
mechanical output, independent Barus fluidity, the complete directional tangent,
zero-supply native parity, closed supply/refusal behavior, and the annular
hydrostatic-bearing analytic flow with second-order load refinement. These
tests use the original complete film implementation and real fs-math/fs-rand
dependencies; they are numerical verification, not measured bearing validation.

## Combined pressure and sliding work within individual layers

`ShearedFilm::evaluate_with_layers_into` reports each original channel's
complete lower-to-upper layer work, with one implicit layer for uniform
channels. It integrates the complete affine shear stress squared over local
viscosity. Normalized centered pressure/Couette amplitudes retain signed
cross terms and a positive within-layer variance; individual layers cannot
be assigned pressure-only fractions of total sliding heat.

Pressure-only and shear-only limits use their respective physical fractions.
A single full-gap layer exactly recovers the native channel work.
`layer_count`, `channel_layer_range` and `MAX_LAYER_REPORTS` expose the
bounded flattened layout without changing original channel or layer order.
All representability and per-channel work-closure checks finish before any
output publication. A second pure local arithmetic pass publishes identical
validated inputs, avoiding another pressure solve, heap allocation or a
maximum-sized internal layer buffer. Existing force/tangent entries retain
their no-layer evaluation path.

Five focused G1/G3 regressions cover independent stress quadrature, the
analytic 13/14 versus 1/14 work split under sliding reversal, exact native
limits, layer subdivision/velocity scaling and late output rollback. The
profiles remain prescribed within a solve; finite thermal stores, interfaces
and any source feedback retain their separate physical owners.

## Sliding-wall and layered Reynolds work

`resistive_film::shear::ShearedFilm` binds explicit lower/upper tangential
velocity rows to the same work-conjugate mechanical ports as native normal
film motion. Uniform-channel and retained-stratified constructors compile
the harmonic shear viscosity and two independently integrated wall-distance
entrainment weights. The existing graph and factorization solve
`L psi = A B v - div(Qc)`. Combined normal pressure and transposed wall
tractions satisfy `F dot v = sum(G*delta_psi*delta_p + S*(delta_U)^2)`;
total flow times pressure drop alone is not the sliding passage's heat.

The directional tangent differentiates aperture conductance, entrainment,
pressure wall tractions and inverse-gap shear drag, reusing the native
pressure factorization. Bounded stack scratch, atomic outputs and original
pressure/gap limits remain. Zero wall rows recover the native pressure
operator and exact passage closure; channels with nonzero wall rows require
positive aperture. Per-channel reports retain total flow, entrainment,
pressure work, shear work, both signed resisting wall tractions and signed
wall input power. Stationary walls retain pressure-driven shear; layered
viscosity shifts the pressure-shear neutral plane through its fluidity
centroid. `ShearChannelReport::shear_stress_pa` observes the affine physical
fluid stress anywhere inside the gap, continuous across viscosity jumps.
It refuses extrapolation outside the gap. Closed passages report zero shear.
Traction times original wall area and velocity reproduces wall input power,
including negative input when pressure drives the walls; this is distinct
from nonnegative irreversible heat. Reports reuse the accepted coupled solve
without another factorization or heap allocation.

A common Barus law `mu(z,p)=mu_0(z)*exp(alpha*p)` retains the same pressure
factorization and entrainment weights. The exact length-average multiplier
`M=delta_p/delta_psi` multiplies the harmonic shear viscosity, including
`M=exp(alpha*p)` when both endpoints have the same pressure. Stable divided
differences retain this limit and its two-endpoint derivative. The complete
tangent differentiates potential, physical pressure and this wall drag; it is
generally nonsymmetric. Zero alpha preserves the original arithmetic. Layer
work uses the same pressure and shear amplitudes, so its cross terms and sum
remain consistent with the mechanical power without another flow solve.
Reported Barus wall traction and cross-gap stress are exact channel-length
averages, not pointwise streamwise maxima.

Four focused G1/G3/G4 groups check independent local Stokes-work quadrature
along and across a layered passage, positive and suction pressure with unequal
wall speeds, multicell continuity and full directional tangents, equal-pressure
Couette drag, fixed-wall native identity and atomic pressure-barrier refusal.
This is the existing classical thin-gap Reynolds approximation within the
declared pressure range; it adds no cross-gap pressure correction, fluid
compressibility, cavitation, thermal evolution or general EHL accuracy claim.

Mean and relative wall velocity rows are compiled once and reused in flow,
force and tangent evaluation. A small relative sliding coordinate therefore
survives a much larger common translation. Original wall rows remain queryable;
unrepresentable derived rows refuse during construction.

### Navier-slip walls

`resistive_film::shear::slip::NavierSlipFilm` binds separate nonnegative lower
and upper slip lengths to the original uniform-viscosity channel graph. Its
single Navier/Stokes profile changes hydraulic mobility, Couette entrainment,
wall tractions and the full geometry/rate tangent together. Pressure elimination
retains the native graph, limits, datum and factorization. Prescribed pressure
supplies and common Barus viscosity remain supported; slip lengths stay fixed
when viscosity changes. This does not infer slip from material names or roughness.

Channel power closes as `Q*delta_p + wall_power = bulk_heat + lower_slip_heat
+ upper_slip_heat`. Each term on the right is nonnegative; the two interface
losses are explicit thermal-routing inputs, not automatically deposited into
the fluid. All quantities use the same accepted solve, bounded stack storage
and atomic outputs. Exact no-slip, equal/unequal-slip Couette and Poiseuille
limits, independent local velocity/heat integration, Barus averaging and full
force derivatives are tested. Stratified viscosity, tensor slip, rarefaction,
cavitation and interfacial state evolution retain their separate owners.

`ShearedFilm::with_layer_viscosities` continues the original uniform-channel or
layered film with one positive coefficient per original layer. It retains graph
geometry, wall rows, normalized interfaces, pressure limits and the common
Barus law. Hydraulic mobility, both entrainment weights, harmonic wall drag and
combined layer work are recompiled together; updating only a hydraulic mean is
insufficient. Unchanged coefficients preserve the exact original arithmetic.
This operation has no thermal state or clock and never reprices completed work.
Four focused regressions compare independent two-layer Couette integrals,
uniform-layer and unchanged-coefficient limits, updated Barus tangents, retained
geometry and transactional invalid-coefficient refusals.

Profiles and wall rows are prescribed. Pressure-only layer work fractions
do not partition combined Couette heat: layerwise pressure/shear cross
terms generally do not cancel. This module introduces no thermal clock,
advection, layer mass, slip, dry friction or evolving interfaces. Seven focused
regressions cover native reduction, free and blocked Couette flow, layered
quadrature, reciprocal work, analytic tangents, translation conditioning
and atomic refusal.

## Accept independently resolved directional slip history

`EllipticCoulomb::accept_resolved_slip` records an independently admitted two-component displacement in the existing native state's tangent frame. Original law identity, coefficients and material axis remain intact. It validates explicit temperature and componentwise representability, returning actual represented increments and residuals. Lost nonzero increments refuse independently of the caller's absolute length tolerance. The state has zero elastic storage and no clock; physical time, contact admissibility and objective frame transport remain with their original owners.

This operation performs no force, normal-load, compliance or impulse solve and makes no work, heat, wear or pressure claim. It is specific to the existing zero-storage elliptical Coulomb state and does not authorize arbitrary internal-variable updates in other material models. Focused regressions cover reversal, frame covariance, identity and numerical-resolution refusal. Source/format review passed; these new native fixtures have not been executed under the current disk/memory block.

> Status: `frankensim-ext-tribo-dry-baseline-tgbj` baseline. This is not an
> admitted material-data result or an Euler-disc validation result.

## Directional instantaneous friction impulses

`tangential::elliptic::impulse` resolves a maximum-dissipation tangent
impulse inside an explicitly declared ellipse at a supplied normal impulse.
It shares the original scalar/matrix projection with sustained friction through
dimension-neutral numerical controls and response data. The borrowed original
law/history supplies coefficient identity and retained texture; impact does
not advance accumulated physical slip, temperature or texture.

The ellipse requires equal static/kinetic principal coefficients. Full SPD
tangent mobility and zero friction axes are admitted. This is a conditional
tangent solve: normal restitution and any normal/tangent coupling remain with
their mechanical owners. No fictitious timestep or finite-slip heat is used.
The applied impulse has the actual body sign. Loss is negative midpoint
impulse work and remains unallocated mechanical loss, including at arrest.

Admission checks actual endpoint velocity, ellipse capacity, associated
final-slip stationarity and an independent inverse-mobility kinetic-energy
change. Impulse magnitude, velocity and energy tolerances are explicit.
The solve allocates no workspace, admits at most 128 root iterations, and
polls cancellation before/after that bounded projection and after admission.
Failure leaves the borrowed history unchanged.

Five impulse groups cover analytical isotropic/principal-axis response,
full-SPD off-axis maximum dissipation, retained-frame and unit covariance,
degenerate directions and limits, and history/cancellation/budget refusal
with deterministic retry. Existing finite-slip and driven/relaxation tests
retain their original physical semantics.

## Accepted force segments for moving-body consumers

`EllipticCoulomb::advance_driven_segments` exposes each accepted held-traction
interval to a fallible visitor without allocating a trajectory collection.
Duration, force, initial/final velocity and contact mode, net slip, unsigned
travel, physical heat and event exactness remain separate. Rejected numerical
trials never reach the visitor. Static holds are real segments with reaction
force and zero slip work.

This entry and `advance_driven` share one integration implementation.
The aggregate entry retains exact zero-drive relaxation. The segmented entry
uses bounded first-order changing-direction evolution even at zero drive so
the original pose owner can integrate the actual force history. Constant-
direction stop/hold/restart intervals retain their analytic timing. Applying
one averaged phase force would generally change positions and work.

Visitors must mutate private staged owners and publish only after success:
a later native or visitor refusal cannot undo previous visits. Cancellation
is checked around visits, including after the final hold. Four focused groups
check independent segment reconstruction, event boundaries, curved zero-drive
refinement against a separate continuum reference, aggregate compatibility,
refusal and deterministic retry.

## Driven directional friction and retained rest

`EllipticCoulomb::advance_driven` evolves the same associated static/kinetic
ellipses under explicit held tangential acceleration and the original full SPD
body mobility. The caller carries a separate accepted Resting/Sliding mode.
Resting validates an explicit actual-velocity residual and uses the complete
static ellipse to hold the load; sliding never becomes rest through a speed
cutoff. Superstatic drive breaks away onto the kinetic branch.

Zero drive reuses the original continuous relaxation. Constant-direction
segments resolve stop, hold and restart analytically. General changing-direction
motion uses bounded first-order endpoint-Coulomb substeps with the existing
matrix projection. Held traction opposes the whole accepted linear-velocity
segment. Its endpoint arrests are timestep-resolved, not exact continuous
events. Ambiguous capture outside the kinetic ellipse is conservatively refined
or refused; this is not a theorem excluding all anisotropic arrests.

The report separates net displacement, unsigned travel, traction impulse,
physical friction heat, supplied drive work and effective kinetic change.
Heat is traction work along the accepted displacement, and the independent
balance is `delta K + D - W_drive`. No numerical energy defect becomes heat.
Force, projection, substep and original quadrature budgets remain explicit;
failed or cancelled proposals publish no native history.

Six driven test groups exercise full-matrix static hold and breakaway, exact
stop/hold/reversal, changing-direction convergence against an independent RK4
reference, isotropic transverse hysteresis, frame covariance, zero-drive and
frictionless-axis limits, and atomic refusal. Curved motion requires refinement
of `maximum_step_s`; source/thermal refresh needs a separate outer refinement.
The leaf supplies no evolving pose, normal-load law or material measurements.

## Finite-arrest quadrature with ordinary anisotropy

The full-stop integration chart in `EllipticCoulomb::relax_velocity` uses
`z = x^8` and its Jacobian `8*x^7`. Nearly equal spectral decay rates can give
the original finite clock an integrable fractional-power endpoint whose
derivatives defeat recursive Simpson error allocation. The transformed clock
and travel integrands vanish sufficiently smoothly at the stop endpoint.
This is a change of variable in the same physical curve; coefficients, force
limits, time/travel tolerances, recursion and evaluation budgets are unchanged.
The short-duration exponential chart and closed-form impulse/heat stay intact.

A G1 regression first reproduces the original `UnresolvedQuadrature` refusal
for capacities [4,3], mobility diag(1,2), velocity [.4,.1]. Its pre-stop and
complete-stop results are checked against an independently inverted analytic
clock and separately integrated displacement/unsigned path. The 21-group native
ellipse target retains analytic limits, frame covariance, stiff short steps and
cancellation. Quadrature errors remain estimates, not certified enclosures.

## Shared occupied-branch dry friction

`stribeck_friction` owns the original scalar Gaussian Stribeck, hyperbolic,
shifted-power and affine load-line continuation formerly embedded in `fs-couple`.
The coupler re-exports that same owner; its modal consumers retain their existing
stable-root and adhesive-offset behavior. `HystereticFriction` makes the shared
owner a checkpoint-bearing `TangentialLaw` for scalar/isotropic held body response.
A sliding branch remains occupied while its stable root exists, including when a
static reaction also becomes admissible. Loss of the root, resolved reversal and
normal release follow the original rule; a force budget cannot choose another
physical branch. The separate memoryless `FrictionLaw` interface is unchanged.

The checkpoint retains the original law or optional ordered dry-interface card,
occupied branch, tangent axis, solved slip and last explicit thermal boundary.
Cards admit actual solved speed and retain their original temperature/pressure
limits and caller authority. Exported heat is traction contracted with actual
slip; no reversible or thermal storage is invented. Preparation, refused trials
and numerical iterations never advance accepted history. A complete candidate
alone may be committed by its surrounding mechanical transaction.

Gaussian critical points and nonlinear sliding roots use safeguarded Newton
proposals inside the original branch brackets. Shared analytic slopes reuse the
coefficient's exponential/power evaluation; unavailable or unhelpful slopes fall
back to the original value-only bisection. Every narrowed endpoint is checked,
and termination still requires adjacent floating-point values. Preparation probes
and selected-root trials share the explicit step iteration budget, enforced
before each update, including derivative-refusal retries. Exact matching load
and mobility reuse a prepared line. Optional bounded
`prepare_state` returns a zero-motion candidate and preparation count for difficult
near-fold or changing-load cases; ordinary first sliding steps need no separate
preparation call.

The one-dimensional hysteresis rule retains its scope: the first loaded motion
selects an axis, and continued sliding or collinear reversal uses that axis.
`ObjectiveTangentialLaw` transports the retained axis and slip with the existing
contact-frame convention. Untransported oblique continuation refuses; arrest and
normal release permit a new direction. This adds no two-axis hysteresis model,
rate/state aging, source acquisition or experimental validation. The original
empirical/benchmark-law boundaries remain. Native regressions exercise coexisting
roots, release/reversal/restart, original-card domains, work/heat, bounded cached
preparation and the existing isotropic matrix-compliance consumer.

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
- `FrictionLaw` provides Coulomb, velocity-dependent, Stribeck, hyperbolic, and
  shifted-power rungs. At
  zero slip `evaluate` reports static capacity only; it never invents a stick
  reaction. `kinetic_coefficient` exposes `μ(|v|)` for a strictly positive
  speed. `kinetic_coefficient_and_slope` returns the same kinetic value and
  its analytic speed derivative in s/m, including the right-hand sliding
  limit at zero and the right derivative of the clamped linear law. It does
  not differentiate the static cone. Value-only calls retain their original
  operation order and do not evaluate derivatives; combined calls reuse the
  single exponential or power decay. Binary-scaled products preserve
  representable slopes across extreme physical scales and after adding a
  weakening tail to the kinetic floor rounds that tail away. An unrepresentable
  slope, or a vanished decay whose amplified derivative remains unresolved,
  explicitly refuses. G0/G1/G3 checks cover physical limits, independent
  differences, unchanged coefficient bits, dimensional scaling, and refusal;
  the retained routed-friction golden remains unchanged.
  `regularized_traction_1d` is the declared ODE embedding: a linear
  ramp through the origin to the kinetic value at scale `v_reg`, then the kinetic rung,
  with force opposing `v_rel`. That map is a solver regularization (the 1-D
  analogue of the fs-contact smooth adapter), not complementarity. A bow, a
  brake pad, and a fault are the same law.
- `TangentialLaw` is the stateful, solver-independent pHS constitutive seam.
  `FrictionLaw` and `DryInterfaceSystemCard` implement it with an implicit static
  reaction followed by a checked direct solution or deterministic sliding
  bisection, never an ODE ramp. `direct_sliding_traction_n` supplies Coulomb
  and linear-velocity affine roots and the hyperbolic law's smaller quadratic
  root. The latter is attempted only after static infeasibility, where the
  admissible sliding root is unique. Its normalized discriminant is a sum of
  nonnegative terms; neither dimensional squares nor a difference of nearly
  equal roots is needed. Every candidate still faces the original force budget
  and constitutive residual; unresolved arithmetic falls back to bisection.
  The modal friction owner consumes the same direct solution.
  The dry tangential owner retains static priority, then selects the stable
  largest-slip Stribeck intersection when three sliding roots coexist.
  It isolates the high-speed load-line minimum analytically before applying
  the effort ceiling: a negative residual at an arbitrary force cap cannot
  hide a physically admissible smaller-force root. The Gaussian turning point
  uses bounded dimensionless Newton updates, charged to the same iteration
  budget as the traction solve. A constant adhesive junction shares this
  isolation because it shifts the load line without changing its derivative.
  This explicit memoryless branch policy does not model hysteretic branch
  continuation; force/iteration failures preserve the original checkpoint.
  Input motion is the released driver-minus-body increment and the explicit
  isotropic body compliance gives `du = du_free - C F`. Returned traction is
  effort entering the interface, conjugate to the **solved** increment, so
  `Delta H + D - F dot du = 0`. Dry laws have `H = 0`, export `D >= 0` through
  `DissipationStep`/`WorkLedger`, and retain accumulated slip plus the optional
  caller thermal boundary. They do not evolve temperature or a third-body layer.
  `TangentialSolveConfig` requires explicit force/residual/iteration budgets
  (1 through 128 root updates); unresolved or unbracketed roots refuse. Candidate
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
  `advance_local_hardness_pressure_weighted_sliding` integrates independently
  supplied `k_i I_i/H_i` samples with their own positive SI indentation hardness.
  It preserves the spatial/temporal correlation of pressure, slip, coefficient
  and hardness, rather than dividing total sliding by an averaged hardness.
  The original common-hardness quadrature delegates to the same compensated,
  single-pass kernel. Auxiliary numerical storage is constant; all samples
  must pass before the caller's volume changes. Source support, mapping bulk
  history to the wearing surface, quadrature error and geometry feedback remain
  with the consumer. Synthetic G0/G1/G3 dry-baseline tests cover independent
  two-region integrals, field interchange, refinement, uniform parity and late
  invalid-hardness rollback; this is not a calibrated life prediction.
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
one-axis friction law. With two active axes it also admits the unique support
point where the static load line is tangent to the ellipse. Boundary equality
uses the represented support without a force-tolerance capture band, and the
returned reaction must satisfy the original world-frame compliance to floating
arithmetic accuracy. Strictly interior static demands remain nonunique and
refuse; no null-space force or regularizing compliance is invented. Invalid,
unrepresentable or unresolved cases leave accepted history unchanged.
The local law supplies no body trajectory, normal impulse or arrest
localization.

`tests/elliptic_coulomb/compliance.rs` exercises independently constructed
associated solutions, scalar equivalence, rotated frames, degenerate axes,
singular response and refusal boundaries. Rank-one boundary regressions cover
oblique reactions, load reversal, signed frame changes, compliance rescaling,
zero static heat, kinetic continuation and adjacent representable demands.

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

For laws permitting numerical scalar-compliance probes, general matrices use
bounded inverse-secant iterations around that same law's
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

Occupied-branch laws such as `HystereticFriction` disclose
`requires_physical_scalar_compliance`: changing their load line can remove a
stable sliding branch even when a spurious static solution closes displacement
and work. These laws receive the original released motion and the exact physical
principal-direction compliance in one constitutive evaluation, without a
static-priority probe. Anisotropic, rotated-principal and singular/free-axis
responses preserve the scalar law's branch continuation, arrest and restart.
Genuinely oblique body response refuses under the existing one-dimensional
continuation contract; no two-axis branch model is inferred. The original full
matrix still gates displacement and actual work, and the selected physical
branch retains the native force budget and cancellation rules. The physical
regressions are in `tests/hysteretic_friction.rs`.

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

### Local-speed friction across a Hertz sliding and spinning patch

`sliding_patch::rate_dependent::HertzFrictionPatch` resolves the original shared
`FrictionLaw` at each local speed `|u + omega cross r|`. The full Hertz pressure
weights force, spin moment and heat together. In particular, zero centre speed
does not erase the friction and heat of a spinning patch, and neither a centre
coefficient nor an RMS speed replaces the original local law.

Coulomb and pure translation retain their exact responses. Unclipped affine
speed laws use the exact Hertz pressure moment `integral r^2 p dA = 2 P a^2/5`,
giving additional force `-slope P u`, spin moment
`-2 slope P a^2 omega/5`, and their corresponding local-slip heat. This includes
velocity strengthening, weakening before the zero clamp, and a purely viscous
curve with zero Coulomb coefficient. A fully clamped patch returns exactly zero.
Clipped affine pure spin integrates its analytically known active central radius
directly. A small positive frictional core therefore cannot vanish between
full-disk nodes or be lost by subtracting large Coulomb responses.
Mixed translation/spin does the same whenever its active disk lies within the
Hertz patch, including an off-centre stationary point. Integration follows the
actual cutoff radius; paired local-slip directions use a rationalized pressure
difference to retain the small net force. Local positive heat, translation force
and spin torque share one quadrature. Spin torque may assist rotation while
translation supplies that work; total interface power remains dissipative.
Independent small-disk asymptotics and full-patch integration check these limits.
Partially intersecting active disks use their actual circle-intersection domain,
including an exterior stationary point and a shrinking active crescent. Local-slip
rays follow the exact Hertz pressure roots, with radial and angular endpoint
transforms and a split at the two circles' intersection. Rationalized near roots
retain the rim limit; sorted-side triangle geometry and `atan2` retain narrow
exterior views whose cosine rounds to one. The same bounded positive quadrature
produces force, moment and heat without Coulomb subtraction. Independent rim,
vanishing-crescent and far-stationary-point asymptotics, full-surface integration,
branch continuity and refinement check these formerly unresolved geometries.
An unrepresentable positive region still refuses explicitly. Very large opposing
translation/spin powers can lose their small difference in floating arithmetic;
directly integrated heat and the separate power residual remain visible.

Other cases use explicit radial/angular composite four-point Gauss resolution.
The integrand subtracts the exact zero-speed Coulomb response, removing the
direction singularity at an isolated stationary point without regularizing
velocity. The radial transform `r=a sin(theta)` resolves the pressure edge;
reflection symmetry, cached radial geometry and compensated sums keep work
bounded and allocation-free. Each panel count is in `1..=256`, bounding original
curve evaluations by 1,048,576; pure spin requires only radial integration.
Panel counts are discretization controls, not a certified error tolerance.
Consumers must refine both axes for mixed motion, or the radial axis for pure
spin. Unrepresentable output or negative integrated heat refuses rather than
clamping. A missed active region with no positive sampled traction also refuses
rather than reporting a cancellation remainder as physical heat. The work
residual retains independently accumulated local heat.

Regressions retain the exact Coulomb/affine limits, compare all non-affine shared
curve families with independent full-disk traction/heat integration, exercise
three resolution rungs, and check objectivity, reversal, replay, rest, lift-off,
clamping and bounded-control refusal. This instantaneous prescribed-kinematics
owner adds no static-onset solve, retained shear, normal evolution, branch
selection or thermal partition. Interface source admission and applicability
remain the material adapter's responsibility; no measured authority is inferred.

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

### Instantaneous coupled Hertz patch impulses

`sliding_patch::contact::impulse::HertzImpulsePatch` solves simultaneous
tangential and normal-axis angular impulses for a prescribed compressive normal
impulse and fixed circular Hertz pressure-impulse distribution. One declared
impact coefficient supplies the integrated local Coulomb budget. Translation
and spin share that original support surface; separate force and torsional
clamps cannot replace it. This instantaneous law does not infer a static/kinetic
history, a timestep or a finite slip displacement.

The mechanical owner supplies the complete positive-definite mobility mapping
`[p_x, p_y, L/a]` to `[dv_x, dv_y, a*domega]`, including all translation/spin
cross terms. The bounded allocation-free Hertz projection is shared with the
sustained wrench owner. It returns either an admitted arrest impulse or the
maximum-dissipation sliding impulse. Physical linear/angular impulse limits,
the original mobility equation and support stationarity of the actual returned
impulse must all pass; tolerances never enlarge the friction set.

The report retains the actual post-impact velocity, including any floating
arrest residual. Dissipation is minus the physical impulse paired with the
midpoint velocity. An independent inverse-mobility kinetic-energy calculation
checks that loss; the endpoint support gap remains a separate optimality
quantity and is not substituted for impact energy. Invalid mobility,
unresolved boundaries, exhausted budgets and cancellation return no candidate.
The borrowed inputs and immutable law remain unchanged.

This is a conditional tangent/spin impact solve: normal restitution, coupled
normal/tangent impulse determination, evolving patch radius, elastic memory,
rolling deformation and heat allocation remain with their physical owners.
`tests/hertz_wrench_impulse.rs` covers analytic translation and Hertz spin
limits, full cross-coupled arrest, mixed motion against independent local
pressure integration, exact zero limits and atomic refusal.

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
coefficients and accumulated local sliding/wear remain outside this rung.
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
shear and body trajectory integration remain outside it.

`tests/mindlin_loading.rs` checks Hertz work, exact fixed-load delegation,
stress-free normal growth, the independent saturated-patch `a^5` unloading
energy, a small normal decrement, separation/recontact, normal subdivision,
load/shear path dependence, atomic memory refusal and spatial refinement of an
independent fixed-coordinate MDR spring quadrature under mixed loading/reversal.

### Simultaneous Hertz--Mindlin indentation and shear

`HertzMindlinContact::advance_prescribed_path` advances simultaneous affine
indentation and scalar macroscopic shear over the caller's interval. The
supplied normal load is the endpoint load; the Hertz force itself follows the
nonlinear indentation law. The caller supplies a prescribed path with zero
body compliance. `advance_coupled_path` instead accepts released tangential
motion and scalar body compliance and solves `du = du_free - C * Q_average`.
Normal indentation remains prescribed. Curved trajectories and normal/body
feedback still require their mechanical owner and trajectory refinement.

The original analytic spring history is retained. In a fixed reference patch,
each local Coulomb bound varies affinely in time. The method resolves spring
entry, yield and detachment events and integrates their force and positive
plastic-work products analytically. Newly contacting material starts relaxed;
detached memory is removed only after its complete work/storage/heat balance.
It returns time-average shear effort, actual endpoint force, normal work,
reversible shear storage change and nonnegative microslip heat separately.
No spring grid, numerical time substeps or normal dashpot is introduced.
Cost is linear in the bounded history; cancellation is checked at every
interval/event partition, and refusals publish no history. Constant-normal
paths delegate exactly to the original fixed-load implementation.

The mechanically coupled path uses that same average traction for both the
body response and physical contact work. At each local spring, the derivative
of time-average strain with respect to total shear is half the square of its
elastic duration. The existing analytic event moments therefore give a
nonnegative constitutive tangent without differencing or an extra spring grid.
The signed force residual is strictly monotone for nonnegative compliance;
Newton steps are safeguarded by the original Coulomb/force-budget bracket.
Normal geometry is cached across trials, trial responses allocate no history,
and only the accepted candidate materializes the bounded history. Convergence
also checks the returned constitutive effort on its own body load line.
The zero-compliance result is identical to the prescribed path. A mechanical
owner can supply `C = dt^2/(2*m_effective)` for a midpoint inertial response;
contact storage and heat then balance work with the body's kinetic-energy
change. No body integrator or normal/shear cross-compliance is inferred.

Five G1/G3/G4 groups in `tests/mindlin_loading.rs` check the independent
proportional Hertz solution, energy recovery during combined unloading,
retained-history subdivision invariance, a three-rung independent spring/time
quadrature, and atomic cancellation/memory/force/domain refusal. They use
synthetic inputs and do not claim material calibration or a full three-dimensional
elastic-friction solution. The model is the same decoupled MDR formulation in
Popov, Hess and Willert (2019), section 4.4, equations 4.21--4.25.

Four additional groups check the independent compliant proportional-growth
solution in elastic and sliding regimes, conversion of retained shear energy
into finite-body recoil during unloading without spurious heat, returned-force
compatibility through changing-load reversals and recontact, exact delegation,
and atomic force/iteration/memory/cancellation refusal. These remain synthetic
model-equation checks rather than material calibration.

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

## Source-line generalized-Newtonian journal shear

`journal::GeneralizedNewtonianJournal` solves concentric, no-slip, finite-gap
Couette flow from a bounded piecewise-affine dynamic-viscosity line in physical
shear rate. The complete line must start at zero shear with positive viscosity;
the stress derivative is positive at both endpoints of every segment. This
admits monotone shear thinning and thickening without inventing a power law.
The original angular-momentum relation `tau(r)=tau_i*(ri/r)^2` and no-slip
angular-speed integral determine the actual nonlinear torque. Stress inversion
and segment integrals are analytic; a capped, bracketed Newton solve uses the
exact angular-speed derivative and carries an explicit angular-rate residual.
Its fallback is bisection. Factored radial stress widths survive narrow
clearances and source-knot crossings; segment rate differences come from the
quadratic stress secant instead of subtracting nearly equal wall rates.
Thinning inverses and integrals use nonnegative discriminant and integral
terms. Each response allocates nothing. The constant
line retains the exact Newtonian torque path. Geometry and the source line are
immutable during a solve; the coupled source owner compiles a new thermal
endpoint and keeps its original receipts.

Signed torque and nonnegative power use the same conjugate angular velocity.
This instantaneous power is not itself a time-integrated heat receipt. The
caller owns Reynolds/quasisteady admission, mechanical integration, heat and
source support. This owner does not infer yield surfaces, viscoelastic memory,
pressure-supported loading, fluid inertia/storage, radial motion, end leakage
or Taylor instability. Independent tests recover finite-gap Newtonian and
power-law solutions through representable narrow-clearance limits, cross a
constitutive kink inside the gap, compare a piecewise thinning line against
separate radial quadrature, and retain the nearly constant limit. A bounded
solve regression spans broad source support down to the zero-shear plateau.
Nonmonotone stress, unsupported shear and exhausted solve budgets are refused.

`journal::slip::NavierSlipJournal` binds explicit nonnegative slip lengths to
the original source-line journal. The curved-wall boundary law uses the
physical shear `du_theta/dr-u_theta/r`: each interfacial angular-speed drop is
`b*gamma_wall/r_wall`. Both wall shear rates and their stress derivatives come
from the same source viscosity line, so a wall's mobility is `b/mu(gamma_wall)`
rather than a frozen apparent viscosity. The stress root retains the analytic
bulk integral, adds both interface drops and their exact one-sided tangents,
and preserves the original bracket and source-rate limits. Constant viscosity
has a direct finite-gap solution; zero slip retains the original response.
Evaluation allocates nothing.

The response separates bulk viscous dissipation from the inner and outer
interface losses, returns fluid-minus-wall slip velocities, and carries the
heat-sum defect against the same torque/angular-velocity work. Interface heat
partition and time integration belong to the coupled owner; the total torque
power must not also be deposited as bulk heat. Common solid-body rotation
does not generate shear. Independent tests recover Newtonian torque and all
three heat integrals, nonlinear power-law Couette flow, a thinning line with
source knots inside the gap against radial volume quadrature, and the
representable thin-gap limit. Unsupported rates and invalid slip data refuse.
This extension does not infer slip lengths, thermal wall models, yield stress,
pressure-supported load, inertia, cavitation or turbulent transitions.

## Quasistatic Reynolds film

`resistive_film` owns the solver-independent incompressible thin-film pressure
law shared with the compatibility surface `fs-flux::resistive_film`. At each
configuration it solves `L(h) psi = A B v` and applies `B^T A p`. Constant
viscosity gives `psi=p` and the reciprocal Newtonian resistance. The optional
`with_barus_viscosity(alpha)` uses `mu(p)=mu_0 exp(alpha*p)`, where `p` is
gauge pressure and `mu_0` is the reservoir viscosity. The exact fluidity
transform `psi=(1-exp(-alpha*p))/alpha` preserves the existing single Cholesky
factorization; `p=-log1p(-alpha*psi)/alpha`. Physical dissipation is
`sum Q*(p_i-p_j)`. The complete tangent differentiates gaps and the inverse
transform and is generally nonsymmetric; callers must retain that full action.
Cell areas, port weights and conductances are never independently normalized.
This law carries no clock or pressure history; the mechanical caller owns
stepping, cancellation between calls and integrated energy accounting.

`with_viscosity` refreshes `mu_0` while retaining `alpha`. Zero `alpha`
preserves the original arithmetic. Flow beyond `alpha*psi<1`, unrepresentable
pressure or viscosity, and existing domain limits refuse atomically. The
coefficient and its supported pressure interval are caller declarations for the
classical incompressible isothermal Reynolds reduction, with no general
high-pressure or elastohydrodynamic accuracy or sourced-material claim.
The gas-storage owner still refuses nonzero `alpha`; the cavitating owner uses
the pressure/content transform described below. Four original test groups
check integrated fluidity, positive and suction branches,
volume/force/work balance, the finite-difference nonsymmetric tangent,
coefficient continuation, zero-alpha identity, barriers and incompatible
storage owners.

`with_channel_viscosities` compiles one positive dynamic viscosity per original
channel. Pressure, full mechanical tangent and physical viscous work consume
the same per-channel conductances; the pressure solve gains no iterations or
evaluation allocations. Each value is uniform across its passage gap. For
fixed-aperture serial segments, the hydraulic viscosity is the length-weighted
arithmetic mean because resistances add. A gapwise temperature profile cannot
be replaced by a wall viscosity. A common Barus coefficient remains admissible,
with each passage retaining its own zero-gauge viscosity. The uniform accessor
`viscosity_pa_s()` now returns `Some(mu)` only for an exactly uniform field and
`None` otherwise; the ordered channel accessor never invents an effective mean.
`with_viscosity` explicitly restores a uniform field. Gas-storage and cavitating
continuations preserve their original masses, storage/state coordinates and
ownership while changing these same conductances. This is a frozen transport
field, not a thermal solver or an inferred material source. G1/G3 tests recover
independent unequal-channel network pressures, Poiseuille mass fluxes, pressure
work and directional tangents, conserve cavitation refill/squeeze mass, and
recover uniform coefficients bitwise. Invalid updates leave prior owners intact.

`resistive_film::stratified::{ViscosityLayer, GapViscosityProfile, StratifiedFilm}`
admits one to 64 contiguous positive-viscosity layers per original channel,
with strictly increasing normalized interfaces ending exactly at one. It retains
the authored layers and compiles their exact piecewise-constant, no-slip
pressure-flow mobility through scaled, centered weighted variance, including
each layer's thickness²/12 contribution. This avoids cancellation of raw second
moments. Equal-viscosity partitions recover the original coefficient exactly.
Cached neutral-plane coordinates and positive layer pressure-work fractions
come from the same momentum integration; multiplying those fractions by native
channel work locates dissipation without inventing a thermal destination.

The compiled native film remains the pressure, analytic-tangent and channel-work
owner. A common Barus multiplier is separable across all layers and preserves
the fractions. Unrepresentable layer weights, mobility contributions,
equivalents or fractions refuse without floors. Profiles remain prescribed and
fixed within each solve; this adds no Couette/slip law, thermal evolution, layer
advection/inventory or moving-interface conservation. Eight focused regressions
specify independent two/three-layer oracles, velocity/dissipation quadrature,
uniform/reflection/scaling limits, a narrow core at 1e40 viscosity contrast,
native Barus/flow/work integration and bounded refusal. The complete real film
module and all 42 unit tests pass Rust type/borrow checking to MIR; native
linking/execution of the new tests remains pending host storage availability.

The admitted cells and channels share immutable geometry storage across clones,
Barus continuations and viscosity recompilations. Updating a uniform or spatial
viscosity field rebuilds only its two channel coefficient arrays, without copying
cell/face closure rows or revalidating unchanged topology. Source objects remain
independent; evaluation still allocates no scratch and uses the original pressure,
tangent and work arithmetic. This removes graph-size copying from repeated
thermoviscous refreshes; it is not a measured real-time throughput guarantee.

`evaluate_with_channels_into` returns signed passage volume flows and their
nonnegative `Q*(p_from-p_to)` work in original channel order, together with the
same force and pressures from one solve. Physical pressure, not transformed
potential squared, defines Barus dissipation. The channel powers sum to the
original resisting-port power, allowing explicit local heat destinations.
Closed passages return exact zero. All outputs publish together after shape,
pressure and reciprocal-work checks. The ordinary force/Newton image retains
its arithmetic and statically excludes channel-report scratch/writes; neither
image allocates during evaluation. Focused native tests independently integrate
heterogeneous Poiseuille/Barus fluidity for both flow directions, close every
cell's volume and total work, check ordinary-result bit parity, and exercise
closed passages and atomic output refusal. Reports add no time, fluid inventory,
heat partition or thermal/advection evolution.

The fixed limits are 64 pressure cells, 256 channels and 256 mechanical ports.
Nonpositive channel apertures close exactly; collapsed volume cells, trapped
pressure pockets, pressure/gap limit violations and unresolved factorizations
refuse without modifying output slices. No leakage, gap floor, pressure clipping
or diagonal repair is introduced. Compressibility, thermal/slip effects, fluid
inertia, evolving topology, fluid radiation and elastohydrodynamic closure are
outside this law. The existing nine pressure-law tests run through the fs-flux
compatibility surface, while fs-couple tests exercise reciprocal mechanical
reaction, joint dissipation and atomic retry in its existing time integrator.

## Pressure-fed compressible gas films

`resistive_film::isothermal::reservoirs::PressureDrivenGasFilm` binds one positive
absolute pressure per original exterior passage, including shut passages;
internal passages require `None`. It retains the original `IsothermalFilm`
mass coordinates, geometry, ideal-gas storage, normal-pressure gradient and
Hessian. The original ambient remains the fixed mechanical back pressure and
free-energy reference. Changing supplies or uniform/spatial viscosity does not
reset gas inventory. Reservoir logarithms are compiled once; evaluation and
the complete geometry/mass/effort/supply-pressure tangent need no heap allocation.

The same logarithmic-mean mobility recovers compressible Poiseuille mass flow
`G*(p_i^2-p_j^2)/(2*R*T)` at the physical gradient. At arbitrary joint discrete
effort it conserves internal mass and reports, together, signed passage mass
flows, nonnegative isothermal free-energy dissipation, signed reservoir supply
and mass-storage power. Supply uses `-R*T*log(p_b/p_ref)*mass_flow_out`; storage
power plus dissipation equals supply. Mechanical pressure work remains in the
original storage gradient. Shape, topology, pressure, finite-arithmetic and
work refusals leave both output slices untouched. Explicit common ambient
recovers the existing passive flow and tangent bitwise.

This distinct active-port type cannot substitute for a passive-only resistance.
A surrounding integrator must advance mass and motion using the same effort
and integrate the supplied work; the existing passive impact attachment is not
silently rebound. Reservoir temperature/composition stay equal to the reference
gas, viscosity is pressure independent, and ideal-gas, laminar, no-slip, thin-gap
assumptions remain. This is a free-energy account, not isolated thermal energy
or enthalpy, with no finite supply storage, channel inertia or thermal owner.
Seventeen original native isothermal tests pass, including seven new groups
covering independent bidirectional Poiseuille/work oracles, a pressure-supported
steady throughflow, arbitrary discrete effort, full supply-sensitive tangent,
legacy parity, viscosity continuation, sealed/reopened inventory and refusals.

## Mass-conserving cavitating liquid films

`resistive_film::cavitation` extends the same `ResistiveFilm` cell/channel/gap
geometry with retained liquid content and absolute pressure. Its Elrod--Adams
complementarity is `p >= p_c`, `0 <= theta <= 1`,
`(p-p_c)*(1-theta)=0`; liquid mass is `rho*area*gap*theta`. A bounded implicit
finite-volume step uses the original Poiseuille conductances and signed,
upwinded Couette entrainment. Both sides of each internal face use exactly the
same flux. Reservoir pressure and incoming liquid content are explicit; partial
reservoir filling is admissible only at cavitation pressure. No material name,
pressure floor, artificial leakage or replacement gap supplies missing physics.

`CavitatingFilm::with_reservoirs` binds independent absolute pressure and incoming
liquid content to every original exterior channel. Internal channels explicitly
carry no reservoir; closed exterior passages retain their authored condition.
The same active-set solve now supports pressure-fed through-flow, starvation,
rupture and reformation with separate supplies. Original liquid inventory and
state identity survive rebinding. The common ambient datum remains the opposing
mechanical pressure and reference of the original viscosity law; it is not
silently replaced by whichever supply is encountered first.

Boundary Barus potentials are compiled once, shared by immutable clones and
used by both assembly and the independent mass-flux audit. Pressure/content
complementarity and finite-pressure admission apply to each supply. Uniform
explicit supplies reproduce the original trajectory exactly. Viscosity refresh
preserves the supplied conditions; endpoint material refresh retains their
absolute pressures and contents, rebuilds their potentials at the new rupture
pressure and refuses an incompatible supply without publishing a candidate.

Five focused regressions compare the saturated result with the native full-film
pressure solver and an independent Barus integral; exercise starvation/refill,
flow reversal, material-density transition and exact uniform compatibility;
and check topology, pressure/content domain, barrier and cancellation refusal.
These are prescribed reservoirs with mass accounting, not finite reservoir
storage or a claim of time-integrated thermal/mechanical energy closure.

`resistive_film::cavitation::CavitatingPowerLawFilm` composes this retained
inventory owner with the existing regularized-power-law constitutive owner.
Each original passage uses the exact combined Couette--Poiseuille pressure
response at its authored relative wall speed, plus upwinded mean entrainment:
`Q = Q_pressure(delta_p, U_upper-U_lower) + Q_entrained*theta_upwind`.
The independently authored mean volume flux and relative wall velocity use
the original channel orientation. The local stress profile, finite Newtonian
core and positive differential pressure conductance are reused without an
effective-viscosity approximation or a second constitutive implementation.

One bounded pressure/content coordinate per original cell enforces the
same `p >= p_c`, `0 <= theta <= 1` complementarity. The nonlinear Jacobian
selects either the true pressure differential or the original implicit
storage/upwind column. A full/rupture corner first selects the phase entered
by its direction and relinearizes without changing physical state; these pivots
consume the same bounded active-set budget. Damped corrections preserve the
physical domain, and the last accepted correction may use the complete budget.
The existing dense elimination still checks saturated connectivity before
publication, including already balanced states. Two reusable bounded candidate
workspaces hold nonlinear trials; no per-iteration allocation is added.

Pressure and entrainment contributions remain separate through compensated
incidence assembly and the independent original cell/global mass audit.
Gross divergence-free translation cannot erase a small pressure flux or
squeeze source, nor enlarge the nonlinear admission tolerance. Each internal
face still contributes the same amount with opposite signs. The original
pressure/volume limits, signed reservoir mass receipt, endpoint normal-force
projection, state ownership and final cancellation check retain authority.
All-Newtonian laws dispatch directly to the original cavitation step, retaining
its trajectory and reports exactly. Law and reservoir refresh preserve accepted
inventory and compatible state identity; a foreign film cannot reinterpret it.

This combination retains homogeneous liquid density and rupture pressure
within each step. It does not supply cavitated tangential traction, viscous heat
partition, latent heat, thermal/material time evolution, pressure-dependent
viscosity or compressibility of trapped full liquid. Relative wall motion is
refused through a closed passage. Original native normal force is the only
mechanical effort in the returned cavitation step.

Nine focused G1/G3/G4 regressions cover thinning/thickening rupture and refill
against independent momentum-profile integrals; fully filled recovery of the
native combined sliding/pressure solve; a reversing multicell reformation
front and separate face mass balances; branching-graph reformation from a
full/rupture corner; refill-time refinement; exact Newtonian
trajectories; inventory-preserving viscosity refresh; large-entrainment
conditioning; and sealed, exhausted-budget or late-cancellation atomic refusal.
All 132 original-source resistive-film tests pass through the established narrow
Rust facade with real `fs-math` and `fs-rand` dependencies. This focused run does
not claim a full constellation build or cross-ISA execution.

`CavitatingFilm::with_power_law` now binds the nonlinear law to the original
compatible liquid inventory, preserving geometry, explicit supplies and mass
tolerances. Native Barus coefficients refuse through the original nonlinear
law constructor rather than being discarded. The new
`CavitatingPowerLawFilm::advance_with_liquid` uses explicit endpoint liquid
data and channel laws in the same original backward-Euler mass equation:
`rho_next*(V_next*theta_next + dt*div Q_next) = rho_old*V_old*theta_old`.
The original accepted state is validated at its old material epoch, and the
new pressure/content solution uses the endpoint density and rupture pressure.
Newtonian and nonlinear continuation share the original coefficient/reservoir
preparation. Explicit supplies must remain admissible at the new rupture
pressure. A changed liquid receives a distinct inventory identity; a law-only
update retains state compatibility. All-n=1 endpoint laws use the native
Newtonian kernel exactly, even if the previous laws were nonlinear.

The returned `CavitatingPowerLawContinuation` pairs the refreshed owner and
candidate state for joint acceptance by the surrounding thermal/mechanical
owner. This first-order endpoint constitutive split infers no temperature,
caloric energy, latent heat or cavitated tangential work. Four new G1/G3/G4
groups exercise nonlinear thermal-expansion pressure against independent
momentum quadrature, sealed heating/cooling mass, exact Newtonian continuation,
retained versus changed liquid identities, explicit supply/Barus refusal and
late cancellation. All 136 original-source film tests pass in the same focused
native-source lane; full constellation and cross-ISA validation remain open.

`cavitation::CavitatingCarreauFilm` now binds original prepared
`LaminarCarreauSlit` passages to that same liquid-content owner. Every face
retains the source's thermal Carreau-Yasuda law, viscosity plateaus, transition,
temperature, supported shear-rate ceiling and bounded constitutive controls.
Width/length must match the original graph; the actual affine aperture is
rebound through the source's `with_gap`. Its simultaneous pressure/sliding
profile supplies `Q_pressure` and the positive pressure derivative, while the
native conservative flux remains
`Q_pressure(delta_p, U_upper-U_lower) + Q_entrained*theta_upwind`.

A private borrowed constitutive choice shares the existing nonlinear
pressure/content engine, compensated incidence assembly, bounded backtracking,
phase pivots and independent cell/global mass audit with the power-law owner.
There is no second pressure solver or apparent-viscosity fit. Source checks
remain active for constant Carreau laws as well: dispatch to an unrestricted
Newtonian flux must not bypass the original rate ceiling. The existing native
Newtonian limit is recovered numerically; bitwise trajectory equivalence is not
claimed for this bounded-source nonlinear path. Iteration workspaces remain
bounded stack storage without per-iteration allocation.

`CavitatingFilm::with_carreau` binds the source to compatible accepted liquid
inventory. Reservoir and passage refresh preserve that identity. Explicit
`advance_with_liquid` continues old mass through endpoint density, rupture and
passage changes using the native backward-Euler balance; changed liquid data
receive a new identity, and the returned `CavitatingCarreauContinuation` pairs
the refreshed owner with its matching candidate. Source temperatures change
only through explicit replacement. Original source, pressure, mass, geometry,
work-budget and cancellation refusals publish no candidate. An existing Barus
law and relative sliding in a closed channel refuse explicitly.

The returned quantities are liquid inventory, signed face flow and endpoint
normal pressure force. Full-gap stress profiles do not determine liquid
streamer geometry in a partly filled passage; their dissipation is therefore
not returned as cavitated heat or traction. Cavitated tangential work, caloric
or latent energy, thermal evolution, compressible trapped liquid, bubble
inertia and EHL remain outside this owner. `CarreauCavitationControls` bounds
original slit work and line search; the native `CavitationSolveConfig` retains
authority over pressure/content iterations and mass admission.

Seven new G1/G3/G4 groups check the constant-source native limit; rupture and
refill against an independent analytic stress inverse and physical-gap
momentum integral; manufactured combined sliding/pressure transport at actual
apertures; multi-cell reformation, reversal and shared-face mass; source and
endpoint-liquid continuation; source ceilings and invalid-domain refusals; and
late cancellation with exact retry. All 13 existing nonlinear power-law
cavitation groups also pass against the shared engine. These are focused
original-source numerical checks, not experimental or full-workspace evidence.

A common Barus law now uses `chi=(1-exp(-alpha*(p-p_c)))/alpha`, with channel
conductances multiplied by `exp(alpha*(p_ambient-p_c))`; supplied viscosities
remain referenced to the reservoir pressure. Zero potential is exactly rupture,
so the same bounded linear active-set solve retains mass and complementarity.
The inverse transform supplies actual absolute pressure for forces and domain
checks. Unrepresentable viscosity and the finite-pressure barrier refuse before
publication. Zero `alpha` preserves the original arithmetic, and uniform or
per-channel viscosity refreshes continue accepted liquid content without reset.
No extra nonlinear solve or iteration buffer is introduced.

`CavitatingFilm::advance_with_liquid` continues the ORIGINAL accepted inventory
through an explicit endpoint material change. The same active-set equation now
balances `rho_next * (V_next*theta_next + dt*div Q_next)` against
`rho_old * V_old*theta_old`. Original mass is reconstructed with the old density,
never by reinterpreting retained volume at the new one. Uniform endpoint density
also applies to reservoir inflow. The original geometric graph and common Barus
coefficient stay fixed; new per-channel viscosities are referenced to the new
reservoir pressure, and the rupture-normalized transform is rebuilt at the new
`p_c`. Heating can increase filling or expel liquid through existing passages;
cooling can reduce filling or draw liquid in. No inventory is reinitialized.

The returned `CavitationContinuation` contains the new native owner and matching
candidate state for joint acceptance. Changed liquid data receive a distinct
state identity; neither epoch can reinterpret the other's accepted mass.
Viscosity-only changes retain the established compatible state identity, and
unchanged inputs preserve the original trajectory. Invalid data, an unresolved
sealed full-liquid expansion or late cancellation publish nothing. This is a
first-order endpoint constitutive split with caller-owned temperature and time
refinement, not a caloric, latent-heat or thermal-advection solver.

Rupture, starvation and reformation change retained liquid content, preserving
mass under wall motion and refill. Each candidate reports cell and global liquid
mass residuals, signed boundary mass, absolute pressures, face fluxes and the
original work-conjugate endpoint normal force. A deterministic active-set solve
is capped at the caller's budget; local dense elimination is required because
upwind transport is nonsymmetric. Cancellation is polled at assembly, active-set
and elimination boundaries. Invalid state, singular saturated trapped pockets,
pressure/gap-domain violations, exhausted budget or failed mass closure leave
accepted state unchanged. Geometry remains capped at 64 cells, 256 channels and
256 ports. An opaque shared owner token prevents another liquid/geometry model
from reinterpreting retained content as a different mass; immutable model clones
can continue the same state. Iterations allocate no memory.

This is a pressure-independent liquid-density, negligible-void-mass lubrication
model, with explicit endpoint source changes available. Content is not vapor
mass; missing liquid volume is not a latent-heat
receipt. Endpoint forces do not certify time-integrated coupled work. Finite
vapor density, bubble inertia, dissolved gas, EHL, thermal evolution and sealed
full-liquid compression require additional owners. Focused G0/G1/G3/G5 tests
cover the original full-film pressure limit, analytical rupture/refill, sealed
liquid conservation, bidirectional starved Couette transport, positivity,
transactional refusals and deterministic replay. Barus regressions additionally
compare the native full-film solution, reconstruct heterogeneous face fluxes and
cell mass, exercise rupture/refill/reformation and viscosity continuation, and
check absolute-pressure shifts, the zero-alpha history, and atomic refusals.
Six endpoint-material regressions independently check drained thermal expansion,
sealed cavity heating/cooling, shared-face mass at changing density, refreshed
Barus refill, unchanged-history identity and transactional refusals. The original
before-mass receipt remains exactly the preceding accepted after-mass receipt.
No experimental accuracy or cross-ISA execution claim follows from these tests.

### Barotropic storage and viscosity in compressed liquid films

`cavitation::CompressibleCavitatingFilm` binds the same native graph, viscosity
field, common Barus law, prescribed supplies and solve budget to an authored
positive constant isothermal bulk modulus `B`. The native liquid density is now
explicitly its density at the ambient reference pressure. The constitutive
law is `rho(p)=rho_ref*exp((p-p_ref)/B)`, with rupture density `rho_c` evaluated
from that same law. This enables physically compressed sealed liquid pockets
without artificial leakage, pressure clipping or a softened numerical modulus.

The conserved coordinate `u=rho(p)*theta/rho_c` has fractional filling below
one and full compressed liquid above one, where `p=p_c+B*ln(u)`. Storage and
upwind mean-wall transport are linear in `u`. Pressure transport integrates
`rho(p)/mu(p)` over each original constant-gap passage using the native Barus
coefficient and channel viscosity; no averaged density replaces the mass flux.
For `x=u-1 >= 0`, the normalized primitive is
`B*expm1((1-alpha*B)*ln(1+x))/(1-alpha*B)`. The implementation evaluates its
endpoint difference directly with compensated `log1p`/`expm1` ratios, including
the removable `alpha*B=1` limit. With zero Barus coefficient this pressure
operator is exactly linear in `x`. Compression and partial filling remain
separate during iteration so a small `p/B` is not rounded away in `1+x`.

One deterministic semismooth Newton system reuses the native nonsymmetric
elimination and compensated face-incidence sums. Positive physical storage
provides the pressure response of sealed saturated islands. The original
iteration cap includes phase pivots; each correction admits at most 32 bounded
backtracks. Iterations allocate no memory. Candidate publication follows a
separate cell/global mass reconstruction and the final cancellation check.
The result exposes actual channel mass flows in kg/s and endpoint normal
pressure forces, with no ambiguous single channel volume-flow convention.

Construction creates a distinct inventory identity because incompressible
states cannot retain their mass under a different pressure-density law. Explicit
compressed initial pressures are available, subject to the original pressure/
content complementarity. Clones and viscosity/supply refresh preserve accepted
state identity. The five G1/G3/G4 regression groups check sealed compression,
release and reformation analytically; connected-cell cavity refill against
independent dimensional quadrature of `rho/mu`; bidirectional supplied flow and
viscosity refresh; the incompressible limit and starved transport; and domain,
budget, foreign-state and late-cancellation refusal. The primary constant-bulk
formulation is described in Brewe and Jacobson, NASA TM-88826, eqs. (1)--(5).

This is a Newtonian, isothermal, negligible-void-mass lubrication image with
caller-admitted constant bulk modulus and pressure-viscosity validity. It does
not infer arbitrary real-fluid compressibility, non-Newtonian rheology, thermal
evolution, latent heat, finite vapor mass, trapped gas, fluid inertia or EHL.
Compressive storage is included in MASS accounting; endpoint force is not a
caloric, viscous-heat or time-integrated coupled-work receipt.

`CompressibleCavitatingFilm::advance_with_liquid` now continues that same
pressure-dependent inventory through explicit endpoint material data. The old
pressure and fraction are decoded with the OLD reference density/reference
pressure/rupture pressure. Only the resulting original mass is normalized by
the endpoint rupture density for the solve. The independent final audit uses
each epoch's own density for its before/after totals. In a sealed full cell at
fixed volume, the pressure change is therefore
`p_new-p_old = p_ref,new-p_ref,old + B*ln(rho_ref,old/rho_ref,new)`.
Heating represented by reduced reference density can raise pressure without
creating mass; cooling can create a cavity before later reformation.

The original native `prepare_liquid` supplies coefficient and reservoir
admission. Explicit supplies retain their absolute pressure and incoming
filling, with refusal when the new rupture/reference domain is incompatible.
The independently authored `B`, geometry and Barus coefficient remain fixed.
`CompressibleCavitationContinuation` returns the new owner and its candidate
together. Changed liquid data receive a distinct inventory identity; no-change
and viscosity-only transitions retain compatibility and reproduce the direct
step exactly. This is a first-order prescribed endpoint constitutive split,
not a temperature/energy integrator or an inferred temperature law for `B`.
Four additional G1/G3/G4 groups check initially compressed sealed expansion,
cavitated/compressed material cycles, pressure-fed transition against independent
dimensional quadrature, exact unchanged behavior and late-cancellation refusal.

`CompressibleCavitatingFilm::advance_with_liquid_and_bulk_modulus` also admits an
explicitly supplied endpoint isothermal bulk modulus. Old compressed inventory
uses the accepted reference density, pressure and bulk modulus; endpoint storage,
reservoir density and pressure-driven mass transport use the proposed law.
Changing B creates a distinct inventory identity even when other material data
are unchanged. For a sealed full cell with unchanged volume,
`p_new = p_ref,new + B_new * (ln(rho_ref,old/rho_ref,new)
+ (p_old-p_ref,old)/B_old)`; the existing complementarity solver handles a
transition below rupture by changing filling, and reverses it on reformation.
No source sound speed is silently converted to an isothermal modulus.
The constant-modulus entry point delegates with the accepted B and retains its
exact unchanged-law behavior. Source/domain/cancellation failures leave the
accepted owner and state untouched. This prescribed endpoint split supplies no
caloric, latent-heat or compressive-work claim. Focused test cases cover sealed
compressed-mass preservation, repeated rupture/reformation, independently
integrated Barus reservoir transport, and owner/refusal/cancellation behavior.
All 13 native compressible-film cases, including the four new modulus cases,
passed with the original film operators, inventory solver and deterministic
math compiled in a focused dependency facade. Only unrelated Carreau module
alternatives were excluded; executed numerical functions and test bodies were
unchanged. The five material source-card coupling tests and the full Cargo
workspace remain unrun.

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

`RateState` and `LuGre` exist for named comparisons, such as the X-Bw2
bow-law bakeoff. They are never render laws. `SmithWoodhouse` is re-exported
here under its benchmark name, but it lives in the ungated
`fs_tribo::rosin_thermal` and is a render law (next section).

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

### Rosin thermal law (`rosin_thermal`, ungated)

Promoted out of `benchmark` on 2026-11-06. It is the bowed-string thermal rosin
model (Smith & Woodhouse 2000; Woodhouse 2003; Woodhouse & Galluzzo 2025), and
it is promoted for the fs-couple bowed-string rosin row; a card without that
row never reaches it. Every input stays the caller's: the module admits no material and
carries no bow geometry. The g0 reachability scan still guards
`fs_tribo::benchmark` (rate-and-state, LuGre).

`SmithWoodhouse` is the plastic thermal model of rosin from Smith & Woodhouse
(2000), J. Mech. Phys. Solids 48:1633-1681. Equation numbers are the paper's.
- Friction (Eq. 21): `mu N = A k_y(T) sgn v`, set-valued in sticking, with `A = pi a^2` proportional to `N`.
- Temperature: the heat balance of Eq. (6) per unit area. Conduction into the stationary and the moving solids uses the discrete half-space Green's function of Eqs. (10)-(14). The moving solid sees only the history within slip `8a/(3 pi)` (Eq. 16 at constant speed), with the cut pulse weighted exactly.
- `k_y(T)` is calibrated from a steady-sliding curve through Eq. (20) (section 5.3): the paper's Eq. (2) fit, or the shifted power law of a `FrictionLaw` (`SteadySlidingCurve::from_friction_law` takes `ShiftedPower`, and `Hyperbolic` as exponent 1; it refuses the other laws by name). `SmithWoodhouse::new` certifies with an interval bound that the calibration temperature increases with speed, and refuses otherwise.
- `RosinCalibration::woodhouse_galluzzo_2025_eq13(rig_solid)` names Woodhouse & Galluzzo (2025) eq. (13) with its rig. Their Fig. 1 data are this paper's steady-sliding data, so the rig is this paper's: 500 um at 3 N and a 10 um layer (section 5.3; Galluzzo, Woodhouse & Mansour 2017, Table 2). Its solids are the caller's `rig_solid`; the rig's rod and wedge were perspex (Table 1). WG2025's eq. (5) puts the rosin's `K_r`, `D_r` in that conduction term instead, which is not the measured rig.
- Time stepping: backward Euler in the losses, with friction from the step-start temperature, over a declared uniform step. The paper uses a forward difference.
- Refusals: a different step length, zero load, a temperature above the calibrated range, a history beyond the declared capacity (at most `MAX_HISTORY_STEPS`), and non-positive parameters.

Each checkpoint retains its full original law parameters, including the uniform
time grid, heat kernels, calibration and history capacity. Rebuilding the same
law permits exact replay; changing those parameters refuses the old checkpoint.
The binding copies only parameters, not the precomputed Green-function weights.
A new endpoint temperature outside the calibrated domain refuses before
returning a candidate thermal history. Consumers retain ownership of mechanical
and thermal step publication. Loaded zero-slip cooling remains supported;
zero-load off-contact cooling still refuses.

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

`tests/sw_calibration_law.rs` (no feature) covers the calibration curves: the eq. (13) curve equals `FrictionLaw::ShiftedPower` bit for bit, `Hyperbolic` maps to rounding, the eq. (13) calibration reproduces eq. (13) through Eq. (20) at relative 1e-12 and its steady sliding approaches eq. (13), an over-fast decay and the other laws refuse by name.

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

## Simultaneous Hertz surfaces and original linear rolling loss — v0

`sliding_patch::contact::multiple::rolling::PreparedHertzRollingContacts`
retains one original full `5N × 5N` positive-semidefinite mechanical response.
Per-patch conjugate coordinates are
`[du_x, du_y, a*dtheta_spin, R*dtheta_roll_x, R*dtheta_roll_y]` and
`[Fx, Fy, Mspin/a, Mroll_x/R, Mroll_y/R]`. The mechanical/material owner supplies
each original fixed nonnegative rolling coefficient `d` in N/m. This algebra
creates no material card or applicability claim; the native viscoelastic sphere
adapter supplies `d=A*N/(R*dt)` from its original source law.

All rolling coordinates are eliminated together through the symmetrically
equilibrated factor of `B=I+sqrt(D)*Crr*sqrt(D)`. The identity belongs to the
original viscous equation, not additional mechanical compliance. One reusable
transfer matrix supplies the full Schur response and released motion to the
original simultaneous Hertz owner, retaining every cross-patch response and
its dependent static-load redistribution. Recovery includes rolling induced
by the solved surface wrench even when initial rolling is zero. Zero damping
constrains effort only; that coordinate retains all native recoil.

Final admission reconstructs motion using the unchanged original full matrix.
It independently rechecks original Hertz kinetic support, force, spin moment
and surface work at that motion, plus the original rolling equation and each
separate heat/work balance. Static candidates retain the original feasible
static reaction; arithmetic residuals do not query a discontinuous kinetic law.
Final support evaluations are reserved within the caller's original evaluation
budget. Small compatibility/work residuals cannot excuse reversed tiny-slip
traction. Numerical Schur pivot failure or cancellation of a physically positive
diagonal refuses without clamping, regularizing compliance or invented heat.

Preparation is limited to 128 patches and 409600 full-response entries. At most
three full-size dense arrays coexist; trials reuse the factors with linear
scratch storage. Existing Hertz sweep, support-evaluation, physical-effort and
residual gates remain explicit. Cancellation is polled at bounded preparation,
solve and audit boundaries; failure leaves the prepared owner unchanged.
The caller retains canonical physical ordering, geometry, normal evolution,
material history and thermal storage. This bounded held-response solve claims
neither universal convergence nor uniqueness of the original Hertz root.

Seven focused G1/G3 groups in `tests/hertz_rolling_multiple.rs` compare against
independent full-equation Gaussian elimination, analytic shear/spin support,
manufactured static traction with nonzero bulk loss, dependent static load
redistribution, induced rolling, zero-damping/load and prescribed-motion limits,
canonical permutations, original-matrix numerical refusals and cancellation
retry. The original nine `hertz_wrench_multiple` groups also pass unchanged.


## Accepted Hertz moving phases — v0

`HertzWrenchPhase` lets the mechanical owner distinguish an already moving
contact from one eligible for the original static hull. `advance_phase` and
the simultaneous `advance_phases` entries retain each original law, but use
its kinetic capacity throughout a declared moving phase. Stationary patches
can still carry their own larger static reactions in the same shared solve.
The combined five-coordinate rolling reduction preserves that per-patch
selection and its original full-response work audits. Zero kinetic capacity
preserves free motion even when the original static capacity is positive.

This prevents a larger static hull from choosing zero midpoint displacement
before physical arrest. The original displacement-only entries are unchanged.
A feasible zero-midpoint kinetic reaction still does not establish an endpoint
stop; the dynamics owner must locate or admit that event and any subsequent
static phase. These entries provide no trajectory, source refresh or stopping
clock. Four focused G1/G3 groups in `tests/hertz_phases.rs` check exact moving
translation and spin, unequal shared moving/static capacities, zero kinetic
traction and combined rolling recovery; original source parameters remain
unchanged on success, refusal and retry.


## Integrated original viscoelastic rolling paths — v0

`rolling_loss::viscoelastic::path::ExponentialRatePath` evaluates the held
scalar equation `q'=drive-decay*q` and keeps endpoint rate, signed integral,
squared-rate integral and integrated displacement as separate quantities.
Weak decay uses bounded compensated analytic series with a positive Gram
form for squared rate. Strong decay uses endpoint interpolation and positive
variance, retaining finite transient work even when decay times duration
overflows. Numerical nonrepresentability or unresolved positivity refuses;
no negative heat is clamped. Cancellation is checked at bounded intervals.

`RollingLossLaw::advance_exponential_sphere_path` admits only the original
viscoelastic sphere card and integrates its actual `M=-A*N*R*omega` law.
It returns moment impulse `-A*N*R*integral(omega)` and independently physical
bulk heat `A*N*R*integral(omega^2)`. The original rolling receipt owns the
integrated work, cumulative heat and lineage-bound checkpoint; its wrench is
explicitly the endpoint wrench. A mean rate cannot replace squared-rate work,
and an RMS rate cannot replace signed impulse. Both monotone rate endpoints
must meet the original source, patch, temperature, frequency and speed gates.

`advance_orthogonal_exponential_sphere_paths` integrates two orthogonal,
undriven rolling components over one common duration. Their decay rates may
differ. The original source gates use the full vector speed, which decreases
monotonically; individually valid components cannot bypass that speed limit.
The same original constitutive owner sums their squared-rate heat integrals,
retains separate signed impulses, and advances one state/checkpoint. Equal
rates are invariant under orthogonal coordinate changes, and a zero component
retains the original scalar result. Driven or unequal-duration vector paths
refuse rather than inherit an unsupported speed-envelope assumption.

The caller owns the held geometry, original normal load, angular axis and
mechanical trajectory. This path utility adds no material calibration, normal
work, thermal storage or contact event detection. Seven G1/G3 groups in
`tests/rolling_viscoelastic_path.rs` compare analytic limits and independent
refined time integrals, weak/strong decay, signed reversal, split continuation,
zero loss and original source/cancellation refusals. A reversal with zero net
impulse retains its positive original bulk heat.
Two further G1/G3/G4 groups check independent component quadrature, combined
heat and one checkpoint, split continuation, equal/zero rates and components,
full-vector source admission, and transactional cancellation.

## Full tensor retained shear layers

`layer_contact::ThirdBodyLayerContact::new_tensor` shares the original retained
layer state, mechanical ports, midpoint integration and native consumers. Its
equilibrium and each Maxwell spring have independent symmetric two-axis
stiffnesses; their principal directions need not coincide. Positive-semidefinite
memory is retained through at most six original Gram rows per arm, including
rank-deficient tensors. The equilibrium tensor is strictly positive definite.
The original isotropic constructor and its closed-form return remain intact.
Scalar stiffness/arm getters return `None` for the tensor material; callers must
use the complete tensor getters instead of silently discarding cross terms.

Storage is `e^T K e/2 + sum(r_j^T K_j r_j/2)`. Maxwell heat is independently
`sum(dt r_mid^T K_j r_mid/tau_j)`. Plastic flow remains the explicitly declared
isotropic Bingham branch. Its tensor return solves one bounded scalar plastic
compliance against the full stiffness, preserving the original force/work/slip
sign convention and `dt <= 2*tau_j` admission. No energy residual becomes heat.
Coefficient continuation retains every strain and branch in original order;
changing stored tensors requires explicit coefficient work. The scalar thermal
derivative API refuses tensor materials because it cannot express their entropy.

The fixed-law objective transport capability refuses a nonidentity tensor
rotation. `transport_tensor_state` instead rotates each Gram factor and every
retained strain together, returning the transformed law/state and an explicit
zero-work storage residual. Rotating history alone is not an objective transfer.

The production `fs-couple::friction_heat::third_body::rigid::tensor` source
adapter projects the existing `ResolvedTensorMaterial` through physical
`sym(t_i tensor n)/h` strain columns and volume `A*h`. It retains original
material orientation, source receipts, Maxwell order and Gram rows. Existing
separated face forces/couples drive actual rigid momenta, and independently
computed plastic/Maxwell/Kelvin loss enters the existing finite layer caloric cell in
the original joint nonlinear solve. Geometry, material mass, caloric identity,
temperature and clock must agree. This fixed-geometry rung also admits
source-affine stored tensors through the original thermal potential. Equilibrium
moduli share its admitted common affine scale; each Maxwell cross term may have
its own affine derivative. Source-defined relaxation, Kelvin viscosity and
plastic kinetics may vary with temperature. Non-affine storage, thermal expansion,
phase change, normal/bending/torsional response and unowned geometry remain refused.

`TensorLayerThermalSlopes` retains complete symmetric equilibrium and Maxwell
derivatives in the same tangent frame. `tensor_thermodynamics` uses the original
physical spring strains to evaluate `S=-e^T K' e/2-sum(r_j^T K'_j r_j/2)` and
`U=Psi+T*S`, including signed cross terms and changing rank-deficient branches.
The shared rigid/finite-heat Newton owner deposits plastic, Maxwell and Kelvin
loss minus coefficient free-energy change and `Delta(T*S)`. Reversible exchange
may cool a cell and is never relabeled as frictional dissipation. Force/history,
temperature, rigid momenta and finite caloric cells are accepted together.

`with_kelvin_viscosity` retains the original positive-semidefinite shear-rate
Gram rows, including zero and rank-deficient tensors. Viscosity acts on the
recoverable shear rate in parallel with the elastic/Maxwell network and before
the series plastic slider. Its physical heat is `Delta_e^T C Delta_e / dt`;
it carries no stored energy and is separate from plastic and Maxwell heat.
Changing only this dissipative tensor preserves accepted strain/history with
zero coefficient work. Joint frame transport rotates the rate tensor as well
as the stored tensors. The original rigid/finite-heat candidate evaluates all
source coefficients at its same trial layer temperature and accounts for each
heat channel once.

`tests/tensor_layer.rs` covers independent cross-term equilibrium, the isotropic
limit, noncommuting and rank-deficient memory, second-order held-strain
relaxation, independent plastic power, rotation covariance and transactional
refusals. The production fixture is
`fs-couple/tests/rigid_layer_contact/tensor.rs`; its manufactured orthotropic
card and finite heat inventory are numerical verification, not measured rosin,
wood or other specimen validation.

## Coupled adhesive translation and spin

`adhesive::spinning::AdhesiveSlidingPatch` integrates uniform junction shear
over the actual circular solid-contact area. One local velocity field,
`u + omega cross r`, supplies force, spin moment and positive dissipated power.
No signed normal load or Coulomb coefficient is substituted for the independently
declared shear strength. The pure limits are `F=tau*pi*a^2` and
`M=2*tau*pi*a^3/3`; combined slip and spin share the same capacity set.

The uniform-disk closed forms from Antali (2024), DOI
`10.1007/s11071-024-09462-6`, use paired convergent series near pure motion and
a bounded complete-elliptic AGM evaluation elsewhere. The zero-slip point at
the rim uses its finite analytic limit. Evaluation is allocation-free, with
32 series terms or at most 16 AGM iterations. Six direct groups in
`tests/adhesive_spinning.rs` check independent positive surface integrals,
pure and detached limits, rim continuity, frame/rate covariance, the convex
dissipation support and invalid input. The instantaneous leaf chooses zero
at complete rest; its separate mechanical owner selects static reactions.
Normal evolution, peeling, elastic microslip and heat allocation remain outside
this frozen-area uniform-junction model.

`adhesive::spinning::contact::AdhesiveWrenchContact` uses this same integral to
solve static support and kinetic motion against the complete three-coordinate
body compliance. It shares the existing bounded circular-support projection
with the Hertz owner while retaining the uniform adhesive distribution and its
`2a/3` first-moment radius. Original static and kinetic strengths stay distinct;
an accepted moving phase cannot gain the larger static capacity. Force, torque,
motion, native dissipation, dual gap and body-work residuals must all pass their
declared budgets. No contact elasticity or numerical heat is introduced.
Five focused groups in `tests/adhesive_wrench.rs` check finite torsional recoil,
all translation/spin cross responses, the joint static capacity, moving-phase
admission, covariance, zero strength, cancellation and retry. The original
Hertz and multiple-contact regression targets exercise the shared projection.

`adhesive::spinning::contact::multiple::PreparedAdhesiveContacts` resolves all
attached uniform junctions on the original complete radius-scaled compliance.
Each patch supplies its own actual radius and static/kinetic shear strengths;
there is no normal-load surrogate. A crate-private generic circular-support
engine shares the existing PSD preparation, dependent static-load allocation,
bounded block solve and final physical audit with the Hertz owner. Every
support evaluation, torsional bound and static separating margin uses the
actual traction law, including the adhesive `2a/3` first-moment radius.

Dependent response can redistribute a statically feasible reaction even when
the minimum-norm seed exceeds one junction's capacity. An unresolved static
decision explicitly refuses before trying a weaker kinetic branch. Supplied
moving/static phases keep their respective original capacities. Final work
uses the unchanged full matrix and final common motion; block sweeps never
advance physical time or deposit heat. Zero-strength rows remain unconstrained
and force-free. Zero-area detached junctions are omitted by the normal owner.
Cancellation, budgets and late closure refusals publish no partial candidate.
The shared solver retains its bounded numerical decision and nonuniqueness
limits; it introduces no regularization, compliance, or contact storage.

Five G1/G3/G4 groups in `tests/adhesive_wrench_multiple.rs` cover an independent
three-junction nine-coordinate solution, uniform-disk torsional limits,
dependent static redistribution, mixed accepted phases, unsupported motion,
zero-strength rows, physical reorder/tangent covariance, full-matrix work,
cancellation and deterministic retry. Existing Hertz and rolling multiple-
contact targets remain the regression gate for the shared engine.

`AdhesiveWrenchContact::relax_velocity` resolves a continuous fixed-ray branch
through a finite stopping event. It uses the original uniform-disk wrench and
complete positive-definite body mobility in `[u_x,u_y,a*spin]`. Homogeneity
keeps the wrench constant precisely when `W*w = -q/t_stop`; ray admission is
limited to arithmetic roundoff and retains the original matrix, displacement,
work and effort gates. The moving phase integrates impulse, slip, spin and
native heat analytically. After the located zero, a no-drive static remainder
has zero reaction and adds no heat. Actual nonzero rates never acquire the
larger static strength before that event. Curved velocity paths explicitly
refuse even for tiny requested durations; there is no velocity clipping,
friction ellipse, smoothing or artificial dissipation.

This branch freezes real area, shear strengths and body response. It does not
advance poses or the normal law, admit external drives, or approximate curved
paths. Seven groups in `tests/adhesive_relaxation.rs` cover independent pure
translation/spin and full three-coordinate coupled rim solutions, stationary
dwell, time partition, tangent-gauge/reversal covariance, zero kinetic strength,
tiny moving states, unchanged dimensional gates and cancellation/atomic retry.

### Analytic local Hertz--Mindlin surface fields

`MindlinDeresiewiczPoint::surface_response` and
`HertzMindlinContact::surface_response` reconstruct radial Hertz pressure,
signed tangential traction and relative elastic surface displacement directly
from the original retained scalar history. The inverse Abel transform integrates
each piecewise-quadratic spring strain analytically. Evaluation allocates
nothing and takes O(retained intervals), without a new spatial mesh. Rationalized
square-root differences preserve thin slip-ring traction; positive angular
moments and the explicit center limit avoid cancellation in displacement.

The shear field integrates to the accepted endpoint force, not the last step's
path-average force. Virgin partial slip recovers the Cattaneo difference of two
Hertz profiles and gross sliding gives `tau = +/-mu*p`. Pressure and shear vanish
at and outside the current patch, while elastic surface displacement retains
its physical inverse-distance exterior tail. Changing-load observations use the
current patch and retained history, including open states and recontact. These
are relative surface fields in the existing decoupled scalar approximation;
nonproportional vector friction, per-body displacement, accumulated local
sliding and wear evolution require their own owners.

The physical reference is Popov, Hess and Willert (2019), section 4.4, equations
4.25, 4.34 and 4.58, <https://doi.org/10.1007/978-3-662-58709-6_4>.
Four focused groups verify independent virgin pressure/shear/displacement
formulas, tiny slip and center/edge/exterior limits, signed reversal
superposition and area-force convergence, changing normal/shear histories,
local Coulomb bounds and original ownership/domain/frictionless behavior.
All 30 surface, fixed-load and changing-load Mindlin tests pass.

### Curved adhesive slip/spin and finite extinction

`AdhesiveWrenchContact::advance_transient` extends the original constant-source
disk law to curved three-coordinate motion on its complete positive-definite
body mobility. The original exact-ray branch remains the fast path and must
also pass the transient endpoint-rate gate. Other paths use a bounded positive-
stage Bogacki--Shampine 3(2) update of actual force, spin moment, travel and
original physical sliding power. Independent dimensional rate, integrated-
wrench, travel and heat estimates govern refinement; accepted plus rejected
attempts and support evaluations are bounded and cancellable. Original effort
limits and the separate actual-matrix work gate remain in force.

The uniform disk's dissipation is a norm. Conservative original-matrix inverse
and factor residuals supply mass-metric radial decay bounds. Moving substeps
stay strictly before a possible zero crossing. An uncomputed extinction tail
is accepted only when its entire stopping interval lies before the requested
horizon, its interval width fits the time tolerance, and its travel bound plus
accumulated quadrature estimate fits the travel tolerance. Its integrated
impulse follows actual eventual rest and its heat is the remaining relative
kinetic energy. Prefix physical-power quadrature and its independent work
residual remain unchanged; numerical error is never added to physical heat.
The unknown tail travel is explicitly bounded, and the stopping time remains
an interval. These are analytical tail bounds from a numerically resolved
prefix, not validated enclosures of the complete exact flow.

Six groups in `tests/adhesive_transient.rs` check an independent RK4 curved-flow
reference, original full-matrix stopping impulse and energy, requested horizons
before arrest, time/accuracy refinement, pure-limit fast paths, the stricter
near-ray endpoint gate, tangent covariance, cancellation and atomic budget retry.
The branch holds geometry, strengths and mobility fixed without external drive;
it does not evolve normal adhesion, poses or temperature.

### Irreversible local Mindlin sliding and pressure-weighted wear input

`MindlinDeresiewiczPoint::advance_with_surface_sliding` performs the original
fixed-load scalar solve once and returns a private-path observer borrowing the
old history and owning the accepted candidate. It observes one monotone segment;
arbitrary endpoint pairs cannot hide intervening reversals. The original source,
body-compliance, force and history budgets remain. Local sliding is macroscopic
motion minus the change of relative elastic surface displacement. Its analytic
evaluation integrates the nonnegative spring plastic increments through the
positive Abel kernel, avoiding false wear from cancellation in a sticking zone.
Queries use O(old history intervals), stack-only temporary pieces and no new
history clone, solve or radial mesh.

Within the fixed patch, irreversible work per area is `mu*p*abs(delta_s)`.
The exact pressure-weighted sliding integral `I = integral p*abs(delta_s) dA`
is original dissipation divided by positive mu. Explicit zero friction instead
has `I = N*abs(actual_delta_u)` with zero friction work. The original fixed
pressure and coefficient, scalar direction and monotone accepted displacement
are essential; changing load, contact entry/exit or multiple reversals require
their path owners. Reference: Popov, Hess and Willert (2019), equations 4.41
and 4.42, <https://doi.org/10.1007/978-3-662-58709-6_4>.

`ArchardLaw::advance_pressure_weighted_sliding` accepts this spatial load-distance
measure directly and computes `delta_V = k*I/H` through the existing checked,
provenance-preserving cumulative volume update. The original uniform-sliding
entry forms `I=N*distance` and uses the same implementation. This does not infer
surface recession, debris transport, hardness conversion or an energy share.

Four new physical groups verify virgin sticking and microslip, signed reversal
and gross sliding, accepted compliant motion, the frictionless limit, monotone
subdivision and original ownership/domain refusal. Independent area integration
on 512/2048/8192 radii converges to original native dissipation over six shear
history segments, with finest relative error below 3e-5. All 34 fixed-load,
changing-load, surface-field and local-sliding Mindlin tests pass.

### Pressure-weighted sliding through simultaneous Hertz loading and shear

`HertzMindlinContact::advance_coupled_path_with_sliding` executes the original
analytic normal/shear path once and retains its global load-distance integral
in an opaque candidate. Positive constant mu gives `I=D/mu`, including
microslip caused by a shrinking normal load with zero macroscopic shear and
sliding before lift-off. Normal elastic work remains a separate signed port.
The explicit frictionless branch integrates the nonlinear Hertz load along
affine indentation and multiplies it by actual accepted shear travel. A scaled
positive polynomial gives the mean load without endpoint subtraction or
unscaled powers; first contact and lift-off give `mean(N)=2*N_peak/5`.

The observer adds constant work and no history clone, new solve or spatial mesh.
Original state ownership, scalar direction, force/history/iteration budgets,
compliance and cancellation remain in force. It supplies integrated sliding,
not local varying-load sliding fields, evolving wear geometry or a wear law.
The MDR work equivalence used here is Popov, Hess and Willert (2019), equations
4.41--4.42, <https://doi.org/10.1007/978-3-662-58709-6_4>.

Four focused G1/G3/G4 groups exercise exact frictionless Hertz means, temporal
subdivision, unload-induced fretting, separation/recontact, fixed-load reduction,
compliant recoil and transactional refusal. An independent pressure-times-plastic
slip return map refines in space and time across simultaneous load and shear
reversals, including contact entry/exit; it does not reuse native dissipated work.

## Coupled Carreau self-heating and cooling

`resistive_film::carreau::thermal::ThermalCarreauFilm` advances the actual
temperature of an explicitly uniform, recirculating film assembly. The fixed
thermal inventory includes the returning supply fluid; all pressure ports
return at the same body temperature. Original `TabulatedHeatStorage` owns the
authored total sensible C(T), its supported interval and exact piecewise-affine
caloric integral. Density, mass, heat capacity, supply-loop losses and material
provenance are not inferred from viscosity or film volume. This is not an open
film with silently omitted advected enthalpy.

The fixed-coordinate advance admits held tangential wall motion and refuses
gap-changing velocity. `advance_moving_into` and its checked variant instead
advance the original affine configuration as `q_next=q0+dt*v`. Cell and channel
gaps are admitted at both endpoints and the rounded midpoint. Original source
pressure, wall-rate and work admission is checked at the physical starting and
final configurations. Squeeze, tangential walls and held reservoir pressures
share the existing Reynolds graph at `q_mid=q0+dt*v/2`.

At every thermal trial the original Carreau slits are prepared at the mean
endpoint temperature and that complete pressure/squeeze/sliding graph is
solved again. The original caloric change satisfies
`E(T1)-E(T0) = dt*(D(q_mid,(T0+T1)/2) + Q - H*((T0+T1)/2-Tbath))`.
Only actual profile dissipation enters stored heat. Mechanical/supply work,
prescribed and bath heat, caloric residual and native constitutive/continuity
work residual remain distinct. Reported force, pressure and channel responses
are the accepted midpoint mechanics; the returned owner retains endpoint
configuration and endpoint material temperature together.

Resolved volume-cell storage is separate from the zero-storage hydraulic
passages. The affine cell-volume increment uses the same original area/closure
rows as the normal pressure work. Signed midpoint reservoir outflow closes a
reported volume-continuity residual with that increment. Liquid expelled from
the resolved film remains in the returning thermal inventory; neither its
caloric table nor its declared mass is rescaled by the changing gap. Normal
squeeze work and tangential work are already in the native pressure/sliding
power balance and are not deposited a second time.

A bounded local bracket/secant solve starts from the accepted temperature,
obeys an explicit temperature-excursion cap and shared original caloric and
material support, and reserves work for endpoint admission. The net-power sign
guard evaluates both thermal endpoints at the SAME midpoint geometry. It
refuses detected frozen-load equilibrium overshoots, including stiff midpoint
cooling, while the separate physical endpoint check admits the evolving gap.
It does not isolate every interior equilibrium of a nonmonotone load, establish
unconditional uniqueness, certify source support throughout an arbitrary
nonlinear trajectory, or diagnose physical thermal runaway. Unresolved
intervals require caller refinement. One candidate film is allocated per
interval; root iterations reuse original ordered passage/coefficient arrays
and bounded stack scratch. Configuration, sources, outputs and thermal clock
publish together after the final cancellation check. Failed geometry, work,
source or energy admission preserves accepted state and caller output slices.

This rung supplies midpoint geometry and temperature feedback with a fixed
recirculating thermal inventory. Spatial temperature profiles, open-stream
enthalpy, cavitation/latent heat, thermal expansion, fluid inertia and cooling
coefficient calibration remain outside its scope. Focused tests cover linear
cooling, nonlinear heating across caloric knots and both source shift modes,
pressure-fed viscous work, thermally coupled moving-gap squeeze/sliding,
geometric volume and normal-work closure, and second-order nonautonomous
thermal refinement against an independent Arrhenius/Reynolds reference.
Endpoint-only gap/pressure failures, unsupported initial squeeze, work limits,
stiff cooling and late cancellation preserve exact retry. Existing full-film
and cavitation regressions remain applicable to their mechanical/mass owners.

## Open Carreau films with changing liquid mass and donor enthalpy

`resistive_film::carreau::thermal::open::OpenThermalCarreauFilm` couples the
existing moving full-film pressure/temperature solve to actual supply and return
enthalpy. One original `ResolvedHeatCapacityCurve` supplies liquid density,
Cp(T), support and source receipts. The entire chart must be liquid at the same
positive density, consistent with incompressible Reynolds continuity. Retained
mass is that density times the original resolved cell volumes; hydraulic
passages retain no independent storage. The existing recirculating owner keeps
its authored fixed thermal inventory and behavior.

Every exterior channel has a held supply temperature in the original caloric
support. The accepted native flow determines its actual donor: inflow carries
reservoir enthalpy and outflow carries film-midpoint enthalpy. Reversal changes
that donor without changing topology or fabricating return conditions. Original
pressure, squeeze and sliding work already becomes native viscous heat; pressure
flow work is not added again to sensible advection.

Write `dh(a,b)=h(b)-h(a)`, `Tm=(T0+T1)/2`, `dm=m1-m0`, with positive
incoming/outgoing transported masses and total viscous/prescribed/bath heat Q.
The bounded shared root solves
`R=m1*dh(T0,T1)-dm*dh(T0,Tm)-Q-sum_in(min*dh(Tm,Tin))=0`.
For the reported mass residual `delta_m=dm+sum_out(mout)-sum_in(min)`,
`R+delta_m*dh(T0,Tm)` is the full per-step energy residual. A retained initial
reference temperature Tref provides a continuing storage/port ledger whose
residual is `R+delta_m*dh(Tref,Tm)`. Both complete balances and the reduced root
must pass the energy gate. All caloric work uses native original-law enthalpy
differences, avoiding subtraction of a large arbitrary absolute datum. A mass
defect is reported and gated; it never becomes physical heat.

The returned proposal retains endpoint geometry, liquid mass, temperature,
original source and clock together. Per-reservoir signed transferred mass,
sensible energy and actual donor temperature let a surrounding reservoir owner
propose its matching update. The retained reference survives subsequent steps.
Original midpoint mechanics, physical endpoint admission, bounded thermal work,
source limits and transactional output publication remain shared with the
recirculating owner. Root iterations reuse bounded scratch; per-port receipts
are allocated only for the accepted proposal.

This model assumes one full-film liquid, uniform mixed temperature, constant
density and the chart's fixed pressure/composition caloric convention. It does
not supply spatial heat transport, cavitated-streamer heat, phase change,
thermal expansion, compressibility or automatic finite-reservoir evolution.
Focused regressions compare hot and reversed through-flow, squeeze expulsion
and two-source refill with independent Arrhenius/Reynolds caloric balances;
they also cover nonlinear flow across original Cp knots, large-datum invariance,
continuing changing-mass storage, incompatible sources, invalid donor/gap/work
proposals and late cancellation with exact retry.

## Original thermal Carreau Reynolds films

`resistive_film::carreau::CarreauFilm` couples original
`fs_material::fluid_transport::rheology::slit::LaminarCarreauSlit` passages to
the live incompressible, full-film `ResistiveFilm` graph. Normal gap motion,
prescribed reservoir pressures and both tangential wall velocities share one
pressure solve. Every channel retains its original thermal Carreau-Yasuda law,
prepared temperature, physical shear-rate ceiling and source convergence
controls. A replacement material state is installed explicitly through
`with_passages`; evaluation does not advance a thermal clock or fit a scalar
effective viscosity.

The private, statically dispatched Reynolds engine is shared with the existing
regularized power-law owners. Their pressure Newton solve, positive Cholesky
Jacobian, line search, continuity and work gates remain the same. Source
zero-shear viscosity supplies the geometry seed, and a truly constant original
law retains the native Newtonian pressure path. Nonlinear pressure flow and
wall shear are resolved together by the original slit profile. Mean wall
translation is added as exact entrainment while the source resolves centered
relative motion; large common translation therefore does not erase a small
physical sliding rate. Source-supported zero-shear plateaus, high-rate
plateaus, transition curvature and thermal shifts are all retained.

Reports retain signed wall and pressure work, individual wall tractions,
source-integrated nonnegative viscous heat and the original constitutive work
defect. Requested mechanical work is not relabeled as accepted heat. The
complete graph balance includes pressure-supply work. Analytic pressure and
relative-wall mobilities from the source drive both Newton and the full
gap/velocity/reservoir directional tangent. In particular, with physical gap
`h`, relative wall velocity `r`, signed pressure drop `dp`, pressure-flow
component `Qp`, and source derivatives `g = dQp/ddp`, `Qr = dQp/dr`, the shared
geometry derivative is `dQp/dh = (2*Qp + dp*g - r*Qr)/h`. No finite differences
are used in production tangents.

Passage width and length must match the original graph. `with_gap` rebinds the
actual affine aperture without preparing temperature again. A closed
fixed-wall passage has exactly zero flow; an authored moving-wall passage must
remain open. Existing Barus pressure dependence is refused rather than
discarded. Successful evaluations and tangents use bounded stack scratch and
publish all outputs only after the source, pressure, continuity, work and final
cancellation checks succeed. Original source support or exhausted work
refuses the candidate without changing accepted caller output.

The scope is quasi-steady, incompressible thin-gap Reynolds flow at a frozen
source temperature: no cavitation, fluid inertia, thermal evolution, slip or
yield memory is inferred. The existing graph size caps and nonlinear solve
budgets apply. Focused tests exercise an independent nonlinear manufactured
pressure/flow/traction/heat solution, the Newtonian supplied sliding limit,
full directional tangents and reciprocity, multi-cell continuity and signed
reversal, temperature refresh, common translation, closed passages, source
refusal and late cancellation with exact retry. Existing regularized power-law
regressions remain applicable to the shared engine.

## Regularized power-law Reynolds passages

`resistive_film::power_law::PowerLawFilm` adds a stationary-wall,
incompressible, full-film generalized-Newtonian owner on the original
`ResistiveFilm` graph and its immutable geometry. Each original channel has
explicit SI consistency `K`, positive stress index `n`, and a positive
physical core shear rate `gamma0`, with
`mu(gamma) = K * max(gamma, gamma0)^(n - 1)`. This is a continuous,
strictly increasing stress law: it supports thinning and thickening without
inferring yield stress, shear history, material authority, or an arbitrary
numerical viscosity clamp. The `n = 1` image retains the native Newtonian
force, pressure, channel work, and tangent arithmetic.

The no-slip velocity profile is integrated analytically through the
Newtonian core and power-law region. Positive differential conductance and
the exact aperture derivative drive a bounded, damped Newton solve on the
existing cell-volume balance and equilibrated pressure factorization.
Geometry/laws are compiled once; call scratch is bounded and hot evaluation
does not allocate. The accepted Jacobian factorization is reused by the
analytic directional tangent; the Newtonian tangent shares the original
native implementation without a second pressure solve.

`PowerLawControls` declares iteration and line-search budgets and absolute
and relative cell-flow and work-defect tolerances. Domain handling shortens
the whole numerical pressure correction; it never clips individual physical
pressures. Sealed pockets, unrepresentable constitutive terms, exhausted work,
excessive pressure, failed volume/work closure, and caller cancellation
refuse before publishing caller outputs. Both force and tangent require the
same accepted flow, source-work, and actual projected-force-work checks.
The numerical work defect is diagnostic, never extra heat; irreversible
power remains `sum(Q * delta_p)`. Channel diagnostics retain original order,
including closed passages.

The original binding retains stationary walls; `ShearedPowerLawFilm` below
explicitly adds moving walls to the same constitutive owner. Stratification,
temperature evolution, Barus pressure dependence, cavitation, fluid inertia,
and yield memory remain outside this reduction. The exposed
`pressure_operator()` is the core-viscosity geometry view, not a substitute
for the nonlinear evaluator.

Seven focused regressions cover the exact Newtonian image, independent
no-slip profile integration through the core transition, mixed-law network
continuity, full directional derivatives and reciprocity, SI consistency
scaling, refusal/cancellation atomicity, and shared force/tangent admission
at a one-ULP projected-work difference. The original film module and all 116
of its current native tests compile and pass through a narrow Rust test
facade using the real `fs-math` and `fs-rand` dependencies. This includes all
13 passive and pressure-fed power-law tests. A full crate/workspace build and
cross-ISA execution are not claimed by this focused run.

### Pressure-fed non-Newtonian Reynolds films

`resistive_film::power_law::PressureDrivenPowerLawFilm` binds one explicit
gauge pressure to each original exterior passage and `None` to each internal
passage. The existing exact regularized-power-law channel flux and bounded
Newton/Cholesky pressure solve remain the numerical owners. Pressure relative
to the retained reservoir datum determines hydraulic drops; physical pressure
determines normal force and the original pressure-domain checks. No supply
pressure is silently treated as ambient, including on shut channels.

The separate active type reports mechanical input, signed reservoir input
`-sum(p_boundary * Q_out)`, nonnegative channel heat `sum(Q * delta_p)`, and
their work residual. A supply can drive the mechanical ports or receive work
from them. These instantaneous powers do not supply a time integrator, fluid
inventory or a thermal destination. Every cell must meet the original flow
budgets; supplied films scale the relative tolerance by the largest actual
passage flow or swept-volume rate. Work tolerance includes the magnitudes of
the mechanical and reservoir inputs, retaining accuracy when those ports
exchange power with opposite signs.

The accepted differential conductances and the same factorization give the
analytic force tangent with respect to geometry, normal velocity and every
original reservoir pressure. Evaluation and tangents allocate no heap scratch.
Zero supplies preserve the original passive path; `n=1` uses the original
pressure-fed Newtonian solve. Supply or constitutive-law rebinding retains the
geometry and numerical controls without introducing a thermal clock. All
refusals, including caller cancellation before publication, preserve outputs.
This pressure-only binding omits moving walls; the explicit sliding binding
below adds them. Pressure-dependent rheology, cavitation and fluid inertia
remain outside both combinations.

Six focused G1/G3/G4 groups recover independent through-gap momentum integrals
for thinning and thickening supplied flow, reversed throughflow, mixed-law
cell continuity, both signs of external power, the native Newtonian limit,
geometry/rate/supply finite-difference tangents and velocity reciprocity,
zero-supply parity, law refresh, and atomic domain/budget/cancellation refusal.

### Moving-wall non-Newtonian Couette--Poiseuille films

`resistive_film::power_law::ShearedPowerLawFilm` binds native
`WallVelocityPorts` to a `PressureDrivenPowerLawFilm`. Normal motion, both
tangential walls, prescribed pressure supplies and the original regularized
power-law liquid participate in one complete Reynolds solve. The flow is
computed from the simultaneous nonlinear stress field, not by superposing
a pressure-only non-Newtonian flux and an independently evaluated drag.

For a centered through-gap coordinate `x` in `[-1,1]`, stress is
`tau=c-delta_p*h*x/(2*length)`. A bounded scalar solve finds `c` from the
integral of the original inverse stress law and the authored relative wall
speed. Integration is analytic across the Newtonian-core transitions.
Near-coincident stress endpoints use convergent divided-difference series;
pressure-dominated sliding retains the small antisymmetric stress component.
The positive weighted variance of the inverse constitutive differential
gives passage pressure conductance without subtracting nearly equal second
moments. The native graph Newton solve and pressure factorization are reused.

Mean wall translation adds exactly `width*h*Umean` to flux. Its reciprocal
pressure traction contributes signed wall input and no heat. Lower/upper
resisting tractions, total flow, mean entrainment and signed wall power are
reported in original channel order. Irreversible power is
`(Q-Qentrained)*delta_p + width*length*c*(Uupper-Ulower)`; both terms use the
same combined profile. Full mechanical input plus signed reservoir input
equals this heat within the original declared work budget. A pressure supply
can drive the walls and the walls can pump fluid back into a supply.

Wall-transport incidence is summed with compensation before combining it
with the normal swept source, including in the original native Newtonian
pressure solve. Nonlinear residuals and tangents retain pressure flow and
entrainment separately until the latter's incidence cancels. Moving-film
flow tolerances use reduced pressure forcing and pressure-only flux, so
large divergence-free translation does not hide a small pressure defect.
Channel reports retain pressure-only flow separately from total transport.

The consistent tangent differentiates geometry, normal/sliding velocities
and reservoir pressures, including the implicit stress offset and complete
wall traction. Pure-Couette pressure conductance uses the law's differential
viscosity (`n*mu` outside the core). At the exact uniform core-transition
shear rate, the authored law has no unique full tangent: force remains
available and tangent evaluation refuses atomically. Fixed-wall evaluations
retain the original supplied arithmetic; `n=1` uses the native moving-wall
pressure solve. Constitutive and supply rebinding retain the original wall
ports without introducing a thermal clock.

Every call uses bounded stack storage and no heap allocation. Two reusable
Newton trial workspaces avoid copies of the complete pressure factorization
and channel state on each line-search candidate. Local stress corrections
use the original iteration ceiling and poll caller cancellation. Closed
moving passages, failed local velocity closure, failed graph continuity or
work, nonfinite outputs and exhausted budgets refuse before publication.
The complete film remains homogeneous, no-slip, incompressible and full-film;
this increment does not add cavitation, thermal/layer evolution or inertia.

Seven new native G1/G3/G4 test groups cover independent through-gap momentum
and viscous-heat integrals for thinning and thickening fluids, coupled flow
reversal, both directions of external work, Galilean heat invariance, pure
Couette and native Newtonian limits, pressure/slide scale separation, finite
normal loading under large common translation in both n=1 and n=0.5, full
configuration/velocity/supply finite-difference tangents, reciprocity,
positive incremental resistance, stationary parity, and atomic refusals.
All 123 original-source resistive-film tests pass through the existing narrow
Rust facade with real `fs-math` and `fs-rand` dependencies, including all 20
power-law tests. This focused run does not claim a full constellation build
or cross-ISA execution.

## Elliptical Hertz footprints with original orthotropic sliding traction

`sliding_patch::elliptic::EllipticHertzSlidingPatch` integrates the original
`EllipticCoulomb` kinetic support over a supplied elliptical Hertz pressure
field. The two positive semi-axes and footprint direction are independent of
the retained material texture. Local slip is the actual rigid field
`u + omega J r`; local acting force is the negative original associated
ellipse support. It is not a direction-dependent scalar coefficient or two
independent component clamps. The original ordered static/kinetic coefficients
and texture remain inspectable, although this prescribed-sliding owner uses
only the kinetic support.

With `r=A z`, the original normal load has measure
`3 N/(2 pi) sqrt(1-|z|^2) dz` on the unit disk. Unoriented lines through the
zero-slip point give exact first and second signed radial Hertz moments; force,
centre torque and local physical heat therefore require only one angular
integral. For a stationary point outside the disk, `p=sin(theta)` removes the
rim square-root loss of regularity and avoids forming a remote displaced
origin. Inside and outside formulas use the same boundary limit. Four-point
Gauss quadrature has 1..=256 caller-supplied panels per half-angle interval;
rank-one kinetic friction splits that interval at its actual directional
discontinuity. At most 4096 native support evaluations use bounded stack
storage and poll cancellation at each angular node.

Symmetric directions retain both the support half-mean and half-difference
using a rationalized norm difference. This preserves torque at very small
representable spin rather than subtracting nearly equal resultant forces.
When the two weighted norms differ by more than a factor of four, original
direct supports avoid amplifying roundoff near a one-sided zero. A nonzero
motion ratio that cannot be represented refuses rather than silently becoming
pure translation or central spin. Pure translation uses the original point
support exactly; an isotropic law on a circular footprint delegates to the
existing exact `HertzSlidingPatch`, independent of quadrature resolution.

The returned `HertzSlidingWrench` contains the pressure-integrated acting force,
spin moment, independently integrated nonnegative local friction power and
their generalized-power residual. Complete rest selects zero without solving
static friction. No material history, body momentum, normal evolution, slip
onset, rolling-deformation loss, heat partition or calibration is inferred.
Strong finite shape/texture anisotropy can require many angular panels;
declared resolution and a small work residual are not convergence certificates.

Nine native integration tests compare exact circular and original point-law
limits, an independent two-dimensional Hertz surface integral, a live
32/64/128-panel angular refinement on a sharply directional fixture, arbitrary
frame rotation and axis gauges, velocity reversal/scaling, rank-one pressure
first moments, torque down to spin `1e-280`, the stationary-point rim crossing,
and transactional cancellation/input refusal. The narrow fixture's maximum
force/moment/power errors against 256 panels decrease from `2.48996e-4` through
`9.51483e-6` to `4.72019e-8`. One additional native unit regression checks the
original rank-one support immediately on both sides of a one-sided zero.
All ten tests pass through the isolated original-source dry-tribology build
with actual `fs-math` and `fs-rand`; this is not a full constellation build or
experimental validation.

## Allocation-free changing-area Mindlin coupling probes

`mindlin::loading::path::MindlinPathTrial` exposes the same fixed-load or
changing-area path response and scalar body-recoil solve as the original
Hertz--Mindlin owner, without copying or constructing reversal memory. Its
analytic tangent differentiates average force against actual shear travel;
the supplied mechanical compliance is already used to resolve that travel.
The probe retains the original force, iteration, input and cancellation checks.
Only retained-history capacity is deferred: callers must obtain a successful
`advance_coupled_path` candidate before publishing mechanical or contact state.
Three G1/G3/G4 regressions cover exact probe/candidate agreement through loading,
unloading, reversal and separation, the actual-travel tangent, frictionless
limits, and atomic refusals. No new normal, adhesive or friction law is selected
by this source-neutral probe.

## Local source laws on elliptical Hertz pressure

`sliding_patch::elliptic::resolved::ResolvedEllipticHertzPatch` integrates a
caller-resolved original `EllipticCoulomb` law at every actual rigid slip vector
on the original normal footprint. The pressure geometry uses the same admission
as the uniform-law owner. Each pressure node owns its associated directional
support, centre moment and independently accumulated nonnegative friction heat;
the resultant exposes generalized-power closure and a separate integrated-load
residual. Pressure weights are never renormalized to conceal quadrature error.

Signed rays through the zero-slip point remove the direction singularity. On
each chord a sine substitution removes the Hertz rim square root; rays split at
zero slip when it lies inside the footprint. The exterior rule uses inverse
distance and finite local coordinates, without subtracting two remote positions.
Four-point Gauss rules have explicit radial/angular panel counts and a preflight
cap of 65,536 local law evaluations. The native integration uses fixed stack
storage, polls cancellation before and after each resolver invocation, and
returns no partial wrench. Source callbacks must stage their own side effects.

This is prescribed instantaneous gross sliding, with independently supplied
footprint and texture. It does not choose static reactions, evolve distributed
memory, allocate heat or certify source validity. The original uniform-law
analytic radial evaluator remains available when applicable. Both resolutions
require convergence checks for sharp laws or nearly rank-one ellipses; small
power/load residuals are not spatial error certificates.

Six focused G1/G3/G4 groups exercise independently derived linear-speed Hertz
second moments, the original uniform directional integral, convergence against
the existing circular Stribeck owner, zero slip inside/outside/on the rim,
arbitrary-frame covariance, and preflight budget/source/late-cancellation
refusals. They extend the original nine elliptical-patch integration tests.

## Exact signed-gap timing for evolving-area Mindlin paths

`HertzMindlinContact::{trial_coupled_gap_path,advance_coupled_gap_path}`
compose the original scalar Hertz--Mindlin owner over the exact positive-area
interval of an affine signed indentation path. The signed endpoints must match
the accepted and endpoint repulsive Hertz geometry. Negative indentation has
zero shear load; a separate adhesive normal owner supplies attraction and its
conservative work. There is no endpoint clamping or inferred tensile pressure.

For active duration fraction f, the original solve receives f times the released
total shear displacement and f² times the scalar body compliance. Its average
effort and actual-travel tangent return to whole-step units with factors f and
f². Force ceilings and absolute tolerances are transformed consistently; heat,
elastic storage change and Hertz normal work are physical interval changes and
are not multiplied by f. The active and gap durations are computed separately
to preserve short intervals even when normalized event endpoints round together.
Unrepresentable transformed physical quantities refuse.

Candidate construction advances zero-load gap drift and the original prescribed
active path in chronological order. Departure releases retained strain through
native microslip heat; reentry adds unstrained material after gap motion.
Accumulated macroscopic slip includes the whole actual path. Probes allocate no
history, complete proposals enforce the original history ceiling, and errors or
cancellation publish no state. Fully attached and fully inactive paths delegate
to the existing path methods. This adds no friction law, JKR partial-slip claim,
normal trajectory solver, finite-rotation transport or experimental calibration.

Seven G1/G3/G4 regression groups check explicit touch/gap decomposition, whole
compatibility and analytic tangent, force-budget scaling, release/flight/reentry,
contact-duration representation, original-path parity and atomic refusals.
All 26 changing-load/path tests pass in a focused original-source Rust build
with actual `fs-math`; this is not a full workspace or cross-ISA build.

## Adhesive rolling spring and irreversible yield

`rolling_loss::adhesive::AdhesiveRollingContact` retains fixed-axis elastic
rolling history on one frozen attached JKR patch. Dominik and Tielens (1995),
DOI 10.1080/01418619508243800, equation 25, supplies the spring moment
`M = -4 Fc (a/a0)^(3/2) xi`, where `Fc = 3 pi W R/2`, `xi = R theta`, and W
is work of adhesion. An independently supplied critical shift selects an
elastic/perfectly plastic plateau, explicitly distinct from a microscopic
contact-migration law. A caller ceiling bounds critical shift/contact radius;
the normal owner remains responsible for JKR applicability and attachment.

Each prescribed monotone scalar increment has an analytic elastic/yield split.
Endpoint and work-mean resisting moments, recoverable spring storage, factored
storage change, irreversible rolling loss and applied work remain separate.
Existing double-double arithmetic retains sub-ULP angle and elastic history.
Evaluation is constant-work and allocation-free, with no state mutation on
refusal. Original parameter and limit mismatches reject. Fixed inputs give
deterministic same-ISA results; no new cross-ISA accuracy certificate is claimed.

G0/G1/G3 tests cover elastic recovery, yield and reversal, exact path work,
subdivision, sub-ULP history, and nonzero resistance at tensile JKR pull-off.
The force-controlled pull-off stiffness is half its zero-load value. Arbitrary
vector rolling, spin, sliding, changing normal area, detachment release, finite
rotation, material calibration and heat allocation require their own owners.

## Changing an explicitly prescribed adhesive rolling radius

`AdhesiveRollingContact::refresh_contact_radius` retains the exact elastic
fraction and accumulated angle under the original material and critical-shift
calibration. Only the supplied attached radius changes; the original small-shift
screen is rechecked. The signed spring-energy increment is evaluated from the
stiffness coefficients and retained memory, preserving weak work that subtracting
rounded absolute energies could erase. This parameter work is separate from
irreversible rolling loss. The immutable proposal permits atomic source/history
publication and does not select a normal contact law or infer a heat partition.

## Bonded adhesive torsion with calibrated yield

`adhesive::torsion::AdhesiveTorsionContact` supplies circular full-stick
torsional compliance `K_phi = (16/3) G_bar a^3`, with
`1/G_bar = 1/G_a + 1/G_b`. Both original half-spaces contribute; the Mindlin
translational compliance with Poisson factors is not the torsional compliance.
The actual attached radius remains fixed. An explicitly supplied critical
twist and small-angle ceiling select an elastic/perfectly plastic moment
plateau. This is a lumped bonded-junction idealization, not local Coulomb
partial slip, an inferred atomistic threshold or a mixed-mode peeling model.

The private scalar increment kernel is shared with adhesive rolling, preserving
its allocation-free, constant-work elastic return, compensated weak increments,
mean resisting moment, recoverable storage and separate irreversible loss.
Neither original law's parameter/history identity is interchangeable. G0/G1/G3
tests check unequal half-space compliance, radius-cubed scaling, unloading,
yield reversal, hysteresis, monotone subdivision, tiny increments and refusal.
The normal source owns area/attachment applicability; changing area, normal
axis, material state or temperature needs an explicit work/history transfer.

## Bonded torsion under a prescribed radius change

`AdhesiveTorsionContact::refresh_contact_radius` retains the complete double-double elastic fraction and accumulated twist while the original normal owner prescribes a new circular radius. Both original shear moduli and the calibrated critical twist stay fixed; stiffness and yield moment retain their cubic radius dependence. Signed storage work is evaluated from the coefficient difference before multiplying retained history. Shrinking area returns reversible work rather than inventing plastic loss or heat. Invalid radius, foreign history and unrepresentable arithmetic refuse without mutation.

This is homogeneous bonded twist on the changed disk. It supplies no reciprocal normal effort, stress-free recruitment of annuli, mixed-mode equilibrium, separation or thermal destination. The two added G1/G3 unit groups cover cubic storage/work, loaded yield continuation, reversible radius cycles, sub-ULP memory and invalid/foreign trials. The complete production torsion leaf's seven unit tests passed with pinned nightly and the original moment kernel plus real fs-math; this is not a workspace or cross-crate integration gate.

## Loaded adhesive angular material continuation

Both `AdhesiveRollingContact::refresh_material_response` and
`AdhesiveTorsionContact::refresh_material_response` continue an existing spring
at held angular coordinate with an explicitly supplied attached radius and
material coefficients. Rolling uses the original work of adhesion and reduced
elastic modulus; torsion retains both half-space shear compliances. Geometry,
critical shift or twist, and declared applicability ceilings remain unchanged.
The original constructors recheck the new response before any proposal exists.

The exact double-double elastic fraction and accumulated angle survive the
update. Signed parameter work is formed from the change of the original energy
coefficient before multiplying retained history. It is separate from angular
motion, irreversible loss, and heat. Radius-only refresh delegates to the same
owner with unchanged material coefficients. Source and thermal owners remain
responsible for temperature queries, calibration validity and work supplies.

Four new focused tests check independent analytic material/radius scaling,
signed work and round trips, reverse-yield continuation, sub-ULP history and
atomic foreign-state/invalid-coefficient refusal. All 20 production angular
leaf tests passed using the actual moment/rolling/torsion modules and existing
real `fs-math` artifact. This isolated execution does not claim a full workspace
build, coupled thermal validation or experimental calibration.

## A sliding interface under prescribed slip on the composition contract

`fs_tribo::tangential::component` (feature `component`, off by default; it
brings in `fs-ports` and `fs-exec`) puts a `TangentialLaw` on the
composition contract of `fs-ports` as `SlidingInterface`. The law still owns
traction, solved slip, regime, reversible storage and the dissipation it
exports. The adapter drives it with a released slip increment of
`rate * dt` at zero compliance and with the interface temperature its heat
port receives, and hands the dissipation to that port. It adds no heat
fraction, heat capacity or flash temperature.

- Port `heat` receives the interface temperature [K] and returns the heat
  rate into the interface [W], minus the dissipation rate. Port `traction`
  is a signal output of two values [N].
- Friction that depends on temperature enters as the caller's source
  `law_at(T)`; a law bounded in temperature also receives `T` as its
  thermal boundary. The source is called at every evaluation of a step and
  must be pure.
- No mechanical port: a dry law is set-valued in the slip rate, and a kernel
  that differences smooth port maps cannot solve for the slip through one.
  The work the prescribed motion does against the traction is the
  component's own account (stored energy = the law's reversible energy less
  the work received), and the law's closure residual
  `storage + dissipation - work` is booked as numerical loss.
- Entropy generated: the dissipation over the interface temperature. Class
  `Smooth`; no tangent is offered. `cx` is not polled: one call is one law
  step within the law's root budget, and `TriboError` has no cancellation
  value.
- Refused: a temperature or step that is not positive and finite; the
  source's or the law's own refusal, unchanged; a step whose work balance is
  off by more than `work_tolerance_j`, or whose reported stored change
  disagrees with the law's stored energies by more than it (the admission of
  fs-couple's `tangential_heat_candidate`); a source whose law gives the
  held state another stored energy than it was left with.

`tests/sliding_interface.rs` (five tests, `cargo test -p fs-tribo --features
component --test sliding_interface`): with a Coulomb coefficient
`0.8 - 0.01 (T - 300 K)` the interface passes the contract's conformance
check at three slip rates (one oblique, one reversed), three temperatures
and three steps, with heat `-mu N |v|`, traction `mu N` along the slip, the
law's own slip and thermal boundary in the state, and the work as the
account; fifty steps book 2.775 J; at rest it makes no heat. A synthetic
spring (not a material) shows storage in the account (1.28 J in the spring,
1.28 J of work, nothing net), the refusal of a stiffness that follows the
temperature, and the refusal of heat no work pays for. The refusals are
each checked by the field they name.

Against the hand-written coupling: `crates/fs-couple/tests/friction_heat/
graph.rs` declares fs-couple's friction-heat module from this component, a
thermal node and a thermal timeline, and the module accepts every step of
it at its own tests' tolerances (the record is in
`crates/fs-system/CONTRACT.md`).

Not offered: a slip port (it needs a non-smooth solve tier in the kernel);
any law's analytic tangent in temperature; a thermoelastic account for a
reversible law whose stiffness depends on temperature.
